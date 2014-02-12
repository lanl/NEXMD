! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"
!---------------------------------------------------------
! Code for doing QM periodic boundaries
! Ewald summation.
!
! Written by Ross Walker (TSRI 2005)
!
! Theory is based on:
!        Nam, K., Gao, J., York, D.M., JCTC. 2005, 1, 2-13.
!----------------------------------------------------------

subroutine allocate_qmewald(natm)

  use qmmm_module, only : qmmm_nml, qmmm_struct, qmewald, qmmm_scratch
  implicit none

  integer, intent(in) :: natm

  integer :: ier=0

  allocate (qmewald%mmpot(qmmm_struct%nquant_nlink), stat=ier ) !mmpot = QM-MM
  REQUIRE(ier == 0) !Dealocated in deallocate qmmm
  allocate (qmewald%qmpot(qmmm_struct%nquant_nlink), stat=ier ) !qmpot = QM-QM
  REQUIRE(ier == 0) !Dealocated in deallocate qmmm
  allocate (qmewald%coulpot(qmmm_struct%nquant_nlink), stat=ier ) !coulpot = QM-MM direct Coulombic
  REQUIRE(ier == 0) !Dealocated in deallocate qmmm  
  
!And for the reciprocal force array that stores the forces on the MM atoms due to the K space
  if (.not. qmmm_nml%qm_pme) then
    allocate (qmewald%d_ewald_mm(3,natm), stat=ier )
    REQUIRE(ier == 0)
  else
    !We need an natom long scratch array
    allocate(qmmm_scratch%qm_pme_scratch(natm), stat=ier )
    REQUIRE(ier == 0)
  end if

  return

end subroutine allocate_qmewald

subroutine qm_ewald_setup(totkq,kappa, kmaxqx,kmaxqy,kmaxqz,ksqmaxq, natom,nquant, nlink)

!Initial setup for qm_ewald - should be called only once per run.
!This will calculate the total k space vectors and other related values.
!These values should be constant during a sander run.

!Also allocates memory

      use qmmm_module, only : qmewald, qmmm_nml, qmmm_mpi
      use nblist, only : volume
      
      implicit none

!Passed in
      integer, intent(inout) :: totkq
      integer, intent(in) :: kmaxqx, kmaxqy, kmaxqz, ksqmaxq
      integer, intent(in) :: natom !Total number of real atoms.
      integer, intent(in) :: nquant, nlink
      _REAL_, intent(inout) ::kappa

!Local
      integer :: kx,ky,kz,ksy,ksz,ksq
      integer :: mpi_division
      integer :: ier
#ifdef MPI
   include 'mpif.h'
      integer :: i, istartend(2)
      !PART OF WORKAROUND FOR BUGGY INTEL COMPILER
      integer, allocatable, dimension(:) :: gather_array
      
#endif

      ! Set up the kappa value (Ewald coefficient) if not set by the user: should be different from the one used in sander
  
      if (kappa<0) then
          kappa=(10.0d0)/ (volume**0.333333333d0)      
      endif 

      !Calculate the total number of kspace vectors
      totkq = 0

      do kx = 0, kmaxqx
         if (kx == 0) then
            ksy = 0
         else
            ksy = -kmaxqy
         end if
         do ky = ksy, kmaxqy
            if (kx == 0 .AND. ky == 0) then
               ksz = 1
            else
               ksz = -kmaxqz
            end if
            do kz = ksz, kmaxqz
               ksq = kx*kx + ky*ky + kz*kz
               if (ksq <= ksqmaxq .and. ksq /= 0) totkq = totkq + 1
            end do
         end do
      end do

!Now we need to allocate enough memory for the kvectors and the kvector exponential tables.
      if (totkq == 0) then
        call sander_bomb('qm_ewald_setup','INVALID NUMBER OF K VECTORS','Need totkq > 0')
      end if

      mpi_division = (totkq + (qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
      qmmm_mpi%kvec_end = min(mpi_division*(qmmm_mpi%mytaskid+1),totkq)
      qmmm_mpi%kvec_start = min(mpi_division*qmmm_mpi%mytaskid+1,totkq+1)
      qmmm_mpi%totkq_count = qmmm_mpi%kvec_end-qmmm_mpi%kvec_start+1

#ifdef MPI
      if (qmmm_mpi%commqmmm_master) then
        write (6,'(/a)') '|QMMM: KVector division among threads:'
        write (6,'(a)') '|QMMM:                  Start       End      Count'

!The FOLLOWING CODE SHOULD WORK FINE BUT THE INTEL COMPILER SEEMS TO MISCOMPILE
!THE i=1,iminus loop DUE TO A COMPILER BUG. SO I HAVE WRITTEN A WORK AROUND BELOW.
!        !Already know my own.
!        write(6,'(a,i8,a,i8,a,i8,a)') &
!              '|QMMM: Thread(   0): ',qmmm_mpi%kvec_start,'->',qmmm_mpi%kvec_end, &
!                                     '  (',qmmm_mpi%kvec_end-qmmm_mpi%kvec_start+1,')'
!        iminus = qmmm_mpi%numthreads-1
!        do i = 1, iminus
!          call mpi_recv(istartend,2,mpi_integer,i,i,qmmm_mpi%commqmmm,istatus,ier)
!          write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
!              '|QMMM: Thread(',i,'): ',istartend(1),'->',istartend(2), &
!                                     '  (',istartend(2)-istartend(1)+1,')'
!        end do
!      else
!        !Send a message to the master with our counts in.
!        istartend(1) = qmmm_mpi%kvec_start
!        istartend(2) = qmmm_mpi%kvec_end
!        call mpi_ssend(istartend,2,mpi_integer,0,qmmm_mpi%mytaskid,qmmm_mpi%commqmmm,ier)
!      end if
      end if
      if (qmmm_mpi%commqmmm_master) then
!WORKAROUND FOR BUGGY INTEL COMPILER
        allocate( gather_array(2*qmmm_mpi%numthreads), stat=ier)
        REQUIRE(ier == 0)
      end if
      istartend(1) = qmmm_mpi%kvec_start
      istartend(2) = qmmm_mpi%kvec_end

      call mpi_gather(istartend, 2, MPI_INTEGER, gather_array, 2, MPI_INTEGER, 0, qmmm_mpi%commqmmm, ier)
      
      if (qmmm_mpi%commqmmm_master) then  
        do i = 1, qmmm_mpi%numthreads
          write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
              '|QMMM: Thread(',i-1,'): ',gather_array(2*i-1),'->',gather_array(2*i), &
                                     '  (',gather_array(2*i)-gather_array(2*i-1)+1,')'
        end do
      end if
      if (qmmm_mpi%commqmmm_master) then
        deallocate( gather_array, stat=ier)
        REQUIRE(ier == 0)
      end if

#endif
!We only allocate these as the number of kvectors this cpu does since totkq_count is
!qmmm_mpi%kvec_end-qmmm_mpi%kvec_start+1
      allocate ( qmewald%kvec(qmmm_mpi%totkq_count),stat=ier )
      REQUIRE(ier == 0)
      allocate ( qmewald%dkvec(3,qmmm_mpi%totkq_count),stat=ier )
      REQUIRE(ier == 0)
      allocate ( qmewald%dmkv(3,qmmm_mpi%totkq_count),stat=ier )
      REQUIRE(ier == 0)
      if (.not. qmmm_nml%qm_pme) then
        allocate ( qmewald%ktable(6,natom,qmmm_mpi%totkq_count), stat=ier )
        REQUIRE(ier == 0)
      end if
      allocate ( qmewald%qmktable(6,nquant+nlink,qmmm_mpi%totkq_count), stat=ier )
      REQUIRE(ier == 0)

      return

end subroutine qm_ewald_setup

subroutine qm_ewald_calc_kvec(vectors,dvectors,dmkv, kmaxqx, kmaxqy, kmaxqz, ksqmaxq)

!This routine calculated the wave-vector arrays for K-space ewald summation.
!This routine fills kvec with the vectors.

  use nblist, only : volume, recip
  use constants, only : PI2, INVPI, TWOPI, FOURPI
  use qmmm_module, only : qmmm_mpi, qmewald
  implicit none

!For ew_coeff
!#include "ew_pme_recip.h"--replaced by qmewald%kappa--TL

!Passed in
  integer, intent(in) :: kmaxqx, kmaxqy, kmaxqz, ksqmaxq
  _REAL_, intent(out) :: vectors(qmmm_mpi%totkq_count) !Will contain the calculated vectors
  _REAL_, intent(out) :: dvectors(3,qmmm_mpi%totkq_count) !Will contain the calculated vectors
  _REAL_, intent(out) :: dmkv(3,qmmm_mpi%totkq_count) !Will contain the calculated vectors

!Local
  _REAL_ beta, vfact
  _REAL_ rkx(3), rky(3), rkz(3), rksq
  _REAL_ mkv(3), vector
  integer :: loop_count, kx, ky, kz, ksq
  integer :: ksy, ksz, kvec_count
 
  beta = PI2 / (qmewald%kappa*qmewald%kappa)
 
  vfact = INVPI/volume

  loop_count = 0
  kvec_count = 0

  do kx = 0, kmaxqx
     if (kx == 0) then
       ksy = 0
     else
       ksy = -kmaxqy
     end if
     rkx(1:3) = kx * recip(1:3,1)
     do ky = ksy, kmaxqy
       if (kx == 0 .and. ky == 0) then
         ksz = 1
       else
         ksz = -kmaxqz
       end if
       rky(1:3) = ky * recip(1:3,2)
       do kz = ksz, kmaxqz
         rkz(1:3) = kz * recip(1:3,3)
         ksq = kx*kx + ky*ky + kz*kz
         if (ksq <= ksqmaxq .and. ksq /= 0) then
           kvec_count = kvec_count+1
           if (kvec_count >= qmmm_mpi%kvec_start .and. kvec_count <= qmmm_mpi%kvec_end) then
             !only do this threads subset.
             loop_count = loop_count + 1
             mkv(1:3) = rkx(1:3) + rky(1:3) + rkz(1:3) 
             rksq = mkv(1)*mkv(1) + mkv(2)*mkv(2) + mkv(3)*mkv(3)
             dmkv(1:3,loop_count) = mkv(1:3) * TWOPI
             vector = vfact*exp(-beta*rksq)/rksq
             vectors(loop_count) = vector
             dvectors(1:3,loop_count) = FOURPI * mkv(1:3) * vector
           end if
         end if
       end do
     end do
  end do
!Sanity check to see if loop_count ends up equalling totkq. If it doesn't something has
!gone terribly wrong. We should never trigger this.
  if (kvec_count /= qmewald%totkq .or. loop_count /= qmmm_mpi%totkq_count) &
               call sander_bomb('qm_ewald_calc_kvec','INVALID NUMBER OF K VECTORS', &
                                'Need ipt = totkq, should never happen.')

  return

end subroutine qm_ewald_calc_kvec

subroutine qm_ewald_calc_ktable(natom, nquant, nlink, coords, dmkv)

!Written by Ross Walker, TSRI(2005)
!Updated for qm_pme by Ross Walker Nov 2005.

!When doing qm_pme only the QM part of the ktable is required.

  use nblist, only : recip
  use constants, only : TWOPI, HALFPI
  use qmmm_module, only : qmmm_struct, qmmm_nml, qmewald, qmmm_mpi
  implicit none

  !Passed in
  integer, intent(in) :: natom, nquant, nlink  !Total atoms (exc link atoms), real quantum atoms
  _REAL_, intent(in) :: coords(3,natom) !Complete unimaged coordinate array
  _REAL_, intent(in) :: dmkv(3,qmmm_mpi%totkq_count)

!Local
  _REAL_ :: mkv(3)
  _REAL_ :: xcoord, ycoord, zcoord
  integer :: loop_count !Kvector we are working on
  integer :: i, ii, qmid, im

  loop_count = 0  
! do loop_count = 1, totkq
  do ii = qmmm_mpi%kvec_start, qmmm_mpi%kvec_end
    loop_count = loop_count+1
    mkv(1:3) = dmkv(1:3,loop_count)
    if (.not. qmmm_nml%qm_pme) then
      !The coordinates of mm_link_pair atoms should have been replaced with
      !the link atom coordinates hence the below loop should calculate all non-link MM atoms
      !all real QM atoms and ALL link atoms. We will then extract out the link atoms. non-link MM
      !atoms have had their charge neutralised so not calculating their entry in the ktable should
      !not be a problem.
      do i = 1, natom
        !For X coordinates
        xcoord = mkv(1)*coords(1,i)
        ycoord = mkv(2)*coords(2,i)
        zcoord = mkv(3)*coords(3,i)
        !cache the values for doing a vectored cos
        qmewald%ktable(1,i,loop_count) = xcoord
        qmewald%ktable(2,i,loop_count) = xcoord - HALFPI !sin(x) = cos(x-(pi/2))
        !For Y coordinates
        qmewald%ktable(3,i,loop_count) = ycoord
        qmewald%ktable(4,i,loop_count) = ycoord - HALFPI
        !For Z coordinates
        qmewald%ktable(5,i,loop_count) = zcoord
        qmewald%ktable(6,i,loop_count) = zcoord - HALFPI
      end do !i=1,natom

      !Do a vectored cosine -> since we subtracted pi/2 from every other x,y,z value so we
      !will end up with cos,sin,cos,sin...
      call vdcos(6*natom,qmewald%ktable(1,1,loop_count),qmewald%ktable(1,1,loop_count))

      !Cache the qm atom values for use later so we can access them linearly in memory
      do i = 1, nquant
        qmid = qmmm_struct%iqmatoms(i)
        qmewald%qmktable(1:6,i,loop_count) = qmewald%ktable(1:6,qmid,loop_count)
      end do
      !We also need to cache the QM link atoms values.
      do i = 1,nlink
        im = nquant + i 
        qmid = qmmm_struct%link_pairs(1,i) !position of MM link pair atom in main amber array
        qmewald%qmktable(1:6,im,loop_count) = qmewald%ktable(1:6,qmid,loop_count)
      end do
    else
      !if we are using pme for the QM-QM interaction we only need the qm ktable
      !We have to use the original unimaged coordinates for this as that was what the
      !reciprocal lattice vectors were calculated for.
      do i = 1, nquant
        qmid = qmmm_struct%iqmatoms(i)
        xcoord = mkv(1)*coords(1,qmid)
        ycoord = mkv(2)*coords(2,qmid)
        zcoord = mkv(3)*coords(3,qmid)
        !cache the values for doing a vectored cos
        qmewald%qmktable(1,i,loop_count) = xcoord
        qmewald%qmktable(2,i,loop_count) = xcoord - HALFPI !sin(x) = cos(x-(pi/2))
        !For Y coordinates
        qmewald%qmktable(3,i,loop_count) = ycoord
        qmewald%qmktable(4,i,loop_count) = ycoord - HALFPI
        !For Z coordinates
        qmewald%qmktable(5,i,loop_count) = zcoord
        qmewald%qmktable(6,i,loop_count) = zcoord - HALFPI
      end do 
      !Do a vectored cosine -> since we subtracted pi/2 from every other x,y,z value so we
      !will end up with cos,sin,cos,sin...
      call vdcos(6*nquant,qmewald%qmktable(1,1,loop_count),qmewald%qmktable(1,1,loop_count))
      !Do the link atoms as well.
      do i = 1, nlink
        im = nquant + i
        qmid = qmmm_struct%iqmatoms(im)
        xcoord = mkv(1)*coords(1,qmid)
        ycoord = mkv(2)*coords(2,qmid)
        zcoord = mkv(3)*coords(3,qmid)
        !cache the values for doing a vectored cos
        qmewald%qmktable(1,im,loop_count) = xcoord
        qmewald%qmktable(2,im,loop_count) = xcoord - HALFPI !sin(x) = cos(x-(pi/2))
        !For Y coordinates
        qmewald%qmktable(3,im,loop_count) = ycoord
        qmewald%qmktable(4,im,loop_count) = ycoord - HALFPI
        !For Z coordinates
        qmewald%qmktable(5,im,loop_count) = zcoord
        qmewald%qmktable(6,im,loop_count) = zcoord - HALFPI
      end do
      if (nlink>0) call vdcos(6*nlink,qmewald%qmktable(1,nquant+1,loop_count), &
                              qmewald%qmktable(1,nquant+1,loop_count))
    end if
  end do

  return

end subroutine qm_ewald_calc_ktable

subroutine qm_ewald_mm_pot(qm_xcrd, npairs, qm_coords, natom,rij_incore, rij_data,mmcharge,kvec)

!Calculates the potential at the QM atom position due to the ewald sum of the MM atoms.
!Called once per MD step, before doing the SCF since the MM charges and positions don't change
!during the SCF.

!Written by Ross Walker (TSRI,2005)
!Updated for QM_PME by Ross Walker Nov 2005. When doing the QM-MM interaction using PME
!only the real space part is calculated here as the reciprocal space part should already
!have been calculated and filled into mmpot.

!Only REAL qm atoms see mm atoms.
  use constants, only : PI, INVSQRTPI, zero, one, two, AU_TO_EV
  use nblist, only : volume   
  use qmmm_module, only: qmmm_struct,qmewald, qmmm_mpi, qmmm_nml
  implicit none

#include "qm2_array_locations.h"

!Passed in
  integer, intent(in) :: npairs !Number of QM-MM pairs per QM atom. All QM atoms have the same list.
  integer, intent(in) :: natom !Total number of atoms in system MM+QM not including link atoms.
  _REAL_, intent(in) :: qm_coords(3,*) !Coordinates of the QM atoms
  _REAL_, intent(in) :: qm_xcrd(4,npairs)  !Wrapped qm coordinates, in same order as the pair list, laid out as
                                 ! x,y,z,charge,x,y,z,charge...
  _REAL_, intent(in) :: rij_data(QMMMNORIJ,*) !Null if not rij_incore
  _REAL_, intent(in) :: mmcharge(natom) !Amber's charge array but in electron units.
  _REAL_, intent(in) :: kvec(qmmm_mpi%totkq_count)
  logical, intent(in) :: rij_incore !Are RIJ values already stored?

!Local
  _REAL_ :: erfcx, mmpottmp
  _REAL_ :: vec(3), r2, oneRIJ, RIJ
  _REAL_ :: x_cos, y_cos, z_cos, x_sin, y_sin, z_sin
  _REAL_ :: ktgs(8), ksum, mmchg, qmi(3), kvectmp
  _REAL_ :: kappa, total_mm_charge
  integer :: i, j, loop_count

!For ew_coeff
#include "ew_pme_recip.h"
 
  if (.NOT. qmmm_nml%qm_pme) then
    !======================================
    !  Calculate K Space Potential
    !======================================
    !This is done between QM and MM atoms. Compute the potential at the QM atom
    !position for the K-space Ewald summation.
    qmewald%mmpot(1:qmmm_struct%nquant_nlink) = 0.0d0 !Zeros entire array
    loop_count = 0
!    do i = 1, totkq
    do i = qmmm_mpi%kvec_start, qmmm_mpi%kvec_end
      loop_count = loop_count + 1
      ktgs(1:8) = zero
      kvectmp = two*kvec(loop_count)
      do j = 1, natom
      !Do this over all MM atoms (skip QM)
!       if (.not. qm_atom_mask(j)) then !Quicker just to loop through all atoms, mmchg will be zero for QM atoms
        !it is an MM atom
        x_cos =qmewald%ktable(1,j,loop_count)
        x_sin =qmewald%ktable(2,j,loop_count)
        y_cos =qmewald%ktable(3,j,loop_count)
        y_sin =qmewald%ktable(4,j,loop_count)
        z_cos =qmewald%ktable(5,j,loop_count)
        z_sin =qmewald%ktable(6,j,loop_count)
  
        mmchg = mmcharge(j)
        ktgs(1) = ktgs(1) + mmchg * x_cos*y_cos*z_cos
        ktgs(2) = ktgs(2) + mmchg * x_sin*y_sin*z_sin
        ktgs(3) = ktgs(3) + mmchg * x_cos*y_cos*z_sin
        ktgs(4) = ktgs(4) + mmchg * x_cos*y_sin*z_cos
        ktgs(5) = ktgs(5) + mmchg * x_cos*y_sin*z_sin
        ktgs(6) = ktgs(6) + mmchg * x_sin*y_sin*z_cos
        ktgs(7) = ktgs(7) + mmchg * x_sin*y_cos*z_sin
        ktgs(8) = ktgs(8) + mmchg * x_sin*y_cos*z_cos
!       end if
      end do

      !Now loop over quantum atoms
      do j = 1, qmmm_struct%nquant_nlink
        ksum = zero
        x_cos =qmewald%qmktable(1,j, loop_count)
        x_sin =qmewald%qmktable(2,j, loop_count)
        y_cos =qmewald%qmktable(3,j, loop_count)
        y_sin =qmewald%qmktable(4,j, loop_count)
        z_cos =qmewald%qmktable(5,j, loop_count)
        z_sin =qmewald%qmktable(6,j, loop_count)
        ksum = ksum + x_cos*y_cos*z_cos*ktgs(1)
        ksum = ksum + x_sin*y_sin*z_sin*ktgs(2)
        ksum = ksum + x_cos*y_cos*z_sin*ktgs(3)
        ksum = ksum + x_cos*y_sin*z_cos*ktgs(4)
        ksum = ksum + x_cos*y_sin*z_sin*ktgs(5)
        ksum = ksum + x_sin*y_sin*z_cos*ktgs(6)
        ksum = ksum + x_sin*y_cos*z_sin*ktgs(7)
        ksum = ksum + x_sin*y_cos*z_cos*ktgs(8)
        qmewald%mmpot(j) = qmewald%mmpot(j) + ksum*kvectmp
      end do

    end do
     
    !======================================
    !  End Calculate K Space Potential
    !======================================
  end if !not qm_pme

  !==================================
  !  Calculate Real Space Potential
  !==================================
  !This is between QM and MM atoms, since the charge on the MM atoms does not change it
  !needs only be done once on each MD step before the SCF is done.

  !This is to correct for the overcounting in the reciprocal calculation.
  
  !TL:
  ! The usage of Ewald coefficient must be consistent:
  if (qmmm_nml%qm_pme) then
    kappa=ew_coeff  ! The Ewald coefficient from MM part should be used if MM PME contribution is used.
  else 
    kappa=qmewald%kappa  !  Otherwise the Ewald coefficient from QM part should be used.
  end if   
  
 
  ! TL:
  ! Calculate the direct Coulombic interactions between QM and MM atoms within qm_cut
  ! This information is useful to verify/compare Ewald calculations
  qmmm_struct%coulombic_eng=zero
 
 
  ! TL:
  ! The total MM charge is needed to calculate the "self" plasma term missed in the MM Ewald part
  total_mm_charge=zero
  do j = 1, natom  
     total_mm_charge=total_mm_charge+mmcharge(j)
  end do 
  total_mm_charge=-PI/(kappa*kappa)/volume*total_mm_charge
  
  if (rij_incore) then
    !Our QM-MM distances have already been calculated and stored
    loop_count = 0
!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
      mmpottmp = zero
      qmewald%coulpot(i)=zero
      do j = 1,npairs
        loop_count = loop_count+1
        call erfcfun(rij_data(QMMMRIJ,loop_count)*kappa,erfcx)
        mmpottmp = mmpottmp - qm_xcrd(4,j)*(one-erfcx)*rij_data(QMMMONERIJ,loop_count)
                             !qm_xcrd(4,j) = charge in electron units
        qmewald%coulpot(i)=qmewald%coulpot(i)+rij_data(QMMMONERIJ,loop_count)*qm_xcrd(4,j)
      end do !j=1,npairs
      qmewald%mmpot(i) = qmewald%mmpot(i)+mmpottmp + total_mm_charge
      
    end do !i=1, nquant
  else
    !We have to calculate the distance between QM-MM pairs on the fly
!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
      mmpottmp = zero
      qmewald%coulpot(i)=zero
      qmi(1:3) = qm_coords(1:3,i)
      do j = 1, npairs
        vec(1:3) = qmi(1:3) - qm_xcrd(1:3,j)
        r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
        oneRIJ = one/sqrt(r2) 
        RIJ = r2*oneRIJ !one/oneRIJ
        call erfcfun(rij*kappa,erfcx)
        mmpottmp = mmpottmp - qm_xcrd(4,j)*(one-erfcx)*oneRIJ
        qmewald%coulpot(i)=qmewald%coulpot(i)+oneRIJ*qm_xcrd(4,j)
      end do !j=1,npairs
      qmewald%mmpot(i) = qmewald%mmpot(i)+mmpottmp + total_mm_charge
    end do !i = 1, nquant
  end if
  
  qmmm_struct%coulombic_eng=qmmm_struct%coulombic_eng*AU_to_EV
  !======================================
  !  End Calculate Real Space Potential
  !======================================
  return
end subroutine qm_ewald_mm_pot

subroutine qm_ewald_qm_pot(nquant,nlink,scf_mchg,qm_coords,kvec)
!Author: Ross Walker (TSRI, 2005)

!Calculates QM potential in qmpot
!Requires the Mulliken Charges (or equivalent).

  use qmmm_module, only : qmmm_struct, qm2_params, qmewald, qm2_rij_eqns, qmmm_mpi
  use constants, only : PI, INVSQRTPI, zero, one, two
  use nblist, only : volume  
  implicit none

!Passed in
  integer, intent(in) :: nquant, nlink
  _REAL_, intent(in) :: scf_mchg(nquant+nlink) !Charges on the quantum atoms due to the current
                                         !density matrix
  _REAL_, intent(in) :: qm_coords(3,*) 
  _REAL_, intent(in) :: kvec(qmmm_mpi%totkq_count)

!Local
  integer :: i, j
  integer :: loop_count
  integer :: iminus, iplus
  _REAL_ :: esfact, wfact, qmi(1:3), vec(1:3), kvectmp
  _REAL_ :: r2, oneRIJ, RIJ, erfcx, chg_qm
  _REAL_ :: ktgs(8), x_cos, y_cos, z_cos, x_sin, y_sin, z_sin, ksum

#include "qm2_array_locations.h"

!For ew_coeff
!#include "ew_pme_recip.h"  -replaced by qmewald%kappa: TL

!We only need to do most of this work if we are changing the mulliken charges on the QM image atoms
!at every step (qm_ewald=1). If qm_ewald=2 then the charges on the QM image atoms are fixed at the SCF
!charges from the previous MD step and were included in mmpot. But if we are not updating the image charges
!we still need to if this is our first MD step.

  !===================================
  !    Calculate K Space Potential 
  !===================================
  !Calculate the K space potential between only QM atoms
  loop_count = 0
  qmewald%qmpot(1:qmmm_struct%nquant_nlink) = zero !Zero's entire array
!  do i = 1, totkq
  do i = qmmm_mpi%kvec_start, qmmm_mpi%kvec_end
    loop_count = loop_count+1
    ktgs(1:8) = zero
    kvectmp = two*kvec(loop_count)
    do j = 1, qmmm_struct%nquant_nlink
      x_cos =qmewald%qmktable(1,j,loop_count)
      x_sin =qmewald%qmktable(2,j,loop_count)
      y_cos =qmewald%qmktable(3,j,loop_count)
      y_sin =qmewald%qmktable(4,j,loop_count)
      z_cos =qmewald%qmktable(5,j,loop_count)
      z_sin =qmewald%qmktable(6,j,loop_count)
      chg_qm = scf_mchg(j)
      ktgs(1) = ktgs(1) + chg_qm*x_cos*y_cos*z_cos
      ktgs(2) = ktgs(2) + chg_qm*x_sin*y_sin*z_sin
      ktgs(3) = ktgs(3) + chg_qm*x_cos*y_cos*z_sin
      ktgs(4) = ktgs(4) + chg_qm*x_cos*y_sin*z_cos
      ktgs(5) = ktgs(5) + chg_qm*x_cos*y_sin*z_sin
      ktgs(6) = ktgs(6) + chg_qm*x_sin*y_sin*z_cos
      ktgs(7) = ktgs(7) + chg_qm*x_sin*y_cos*z_sin
      ktgs(8) = ktgs(8) + chg_qm*x_sin*y_cos*z_cos
    end do
    do j = 1, qmmm_struct%nquant_nlink
      x_cos =qmewald%qmktable(1,j,loop_count)
      x_sin =qmewald%qmktable(2,j,loop_count)
      y_cos =qmewald%qmktable(3,j,loop_count)
      y_sin =qmewald%qmktable(4,j,loop_count)
      z_cos =qmewald%qmktable(5,j,loop_count)
      z_sin =qmewald%qmktable(6,j,loop_count)
      ksum = zero
      ksum = ksum + x_cos*y_cos*z_cos*ktgs(1)
      ksum = ksum + x_sin*y_sin*z_sin*ktgs(2)
      ksum = ksum + x_cos*y_cos*z_sin*ktgs(3)
      ksum = ksum + x_cos*y_sin*z_cos*ktgs(4)
      ksum = ksum + x_cos*y_sin*z_sin*ktgs(5)
      ksum = ksum + x_sin*y_sin*z_cos*ktgs(6)
      ksum = ksum + x_sin*y_cos*z_sin*ktgs(7)
      ksum = ksum + x_sin*y_cos*z_cos*ktgs(8)
      qmewald%qmpot(j) = qmewald%qmpot(j) + ksum*kvectmp
    end do
  end do
  !===================================
  !  End Calculate K Space Potential 
  !===================================

  !======================================
  !    Calculate Real Space Potential
  !======================================

  !We have to calculate QM-QM distances on the fly for qmewald
! do i = 1, nquant
!Note: Ross Walker - Should really split this loop at some point
!                    but be careful to do balanced in parallel

  esfact = -two*qmewald%kappa*INVSQRTPI
  wfact = -PI/(qmewald%kappa*qmewald%kappa)/volume
  
!  do i =1,nquant_nlink
  do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
    chg_qm = scf_mchg(i)
    qmi(1:3) = qm_coords(1:3,i)
    iminus = i-1
    do j = 1,iminus
      vec(1:3) = qmi(1:3) - qm_coords(1:3,j)
      r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
      onerij = one / sqrt(r2)
      rij = r2*oneRIJ !one/oneRIJ
      call erfcfun(rij*qmewald%kappa, erfcx)
      qmewald%qmpot(j) = qmewald%qmpot(j) - chg_qm*(one - erfcx)*onerij + wfact*chg_qm
    end do
    !i=j
    qmewald%qmpot(i) = qmewald%qmpot(i) + esfact*chg_qm + wfact*chg_qm
    iplus = i+1
    do j = iplus, qmmm_struct%nquant_nlink
      vec(1:3) = qmi(1:3) - qm_coords(1:3,j)
      r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
      onerij = one / sqrt(r2)
      rij = r2*oneRIJ !one/oneRIJ
      call erfcfun(rij*qmewald%kappa, erfcx)
      qmewald%qmpot(j) = qmewald%qmpot(j) - chg_qm*(one - erfcx)*onerij + wfact*chg_qm
    end do
  end do

  !======================================
  !  End Calculate Real Space Potential
  !======================================

  return

end subroutine qm_ewald_qm_pot

subroutine qm_ewald_add_fock(fock_matrix, qmpot, mmpot)
!Author: Ross Walker, TSRI 2005

!Adds ewald potential info (mmpot and qmpot) to the diagonal elements of the fock matrix

  use qmmm_module, only : qmmm_struct, qmewald, qm2_params
  use constants, only : AU_TO_EV, BOHRS_TO_A
  implicit none

!Passed in
  _REAL_, intent(in) :: mmpot(*) !Potential at each QM atom due to MM images
  _REAL_, intent(in) :: qmpot(*) !Self energy of the QM atoms to avoid overcounting
  _REAL_, intent(inout) :: fock_matrix(*) !Fock matrix

!Local
  _REAL_ :: temp_pot
  integer :: i, ia, ib, i1, i2

  !Now add the MMPOT and QMPOT array contributions to the diagonal elements of the fock matrix
  do i = 1, qmmm_struct%nquant_nlink
    IA = qm2_params%orb_loc(1,I)
    IB = qm2_params%orb_loc(2,I)
    temp_pot = (mmpot(i)+qmpot(i))*AU_TO_EV*BOHRS_TO_A
    do I1 = IA,IB
       i2 = qm2_params%pascal_tri2(i1)
       fock_matrix(i2) = fock_matrix(i2) - temp_pot
                         !AU_TO_EV*BOHRS_TO_A converts (electrons/angstrom) to (eV/Bohr)
    end do
  end do 

  return

end subroutine qm_ewald_add_fock

subroutine qm_ewald_correct_ee(ee, mmpot, p)
!This routine should be called after a call to qm2_helect. Up to this point
!the energy for the Ewald sum only included half of the term from the MM atoms
!and half from the QM atoms. The QM atoms should contribute only half but the
!MM atoms should contribute full. This routine corrects EE for this.

  use constants, only : AU_TO_EV, BOHRS_TO_A, zero, half
  use qmmm_module, only : qmmm_struct, qm2_params, qmewald
  implicit none

!Passed in
  _REAL_, intent(inout) :: ee
  _REAL_, intent(in) :: mmpot(*)
  _REAL_, intent(in) :: p(*)

!Local
  _REAL_ :: etemp
  integer :: i, i1, ia, ib, i2

  etemp = zero
  do i = 1, qmmm_struct%nquant_nlink
    IA = qm2_params%orb_loc(1,I)
    IB = qm2_params%orb_loc(2,I)
    do I1 = IA,IB
       i2 = qm2_params%pascal_tri2(i1)
       etemp = etemp - mmpot(i)*p(i2)*AU_TO_EV*BOHRS_TO_A
    end do
  end do
  ee = ee + half*etemp

  return

end subroutine qm_ewald_correct_ee

subroutine qm_ewald_core(ewald_core,core_chg,mmpot,qmpot, scf_mchg, coulpot)
!Computes the interaction of the Ewald potenital with CORE in QM atoms.
!Ewald_core = electron volts
!ewald_core = Sum(Core(i)*V(i))

  use constants, only : zero, half, AU_TO_EV, BOHRS_TO_A, PI
  use qmmm_module, only : qmmm_struct, qmewald, qm2_struct
  use nblist, only : volume
  implicit none

!In parallel all threads do a full 1 to nquant loop but have a part
!of qmpot and mmpot so return a part of ewald_core that is then
!reduced later.

!Passed in
  _REAL_, intent(out) :: ewald_core
  _REAL_, intent(in) :: core_chg(*), mmpot(*), qmpot(*)
  _REAL_, intent(in) :: scf_mchg(*), coulpot(*) 

!Local
  _REAL_ :: ewdpot, coul_temp
  integer :: i

  ewald_core = zero
  qmmm_struct%coulombic_eng=zero

  do i = 1, qmmm_struct%nquant_nlink
    ewdpot = (mmpot(i)+half*qmpot(i))
                         !AU_TO_EV*BOHRS_TO_A converts (electrons/angstrom) to (eV/Bohr)
    ewald_core = ewald_core + ewdpot*core_chg(i)
    qmmm_struct%coulombic_eng=qmmm_struct%coulombic_eng+scf_mchg(i)* coulpot(i)
  end do
    
  qmmm_struct%coulombic_eng=qmmm_struct%coulombic_eng*AU_TO_EV* BOHRS_TO_A
  
  ewald_core = ewald_core*AU_TO_EV * BOHRS_TO_A

  return
   
end subroutine qm_ewald_core

subroutine qm_ewald_get_forces(qm_xcrd, qm_coords, natom, &
                               scf_mchg, qmmmrijincore, &
                               dxyzqm, dxyzcl, dvectors,mmcharge)
!Calculates the gradient at the QM atom for the Ewald summation by all atoms.

!Needs real scratch space of size ar least nquant.

!The gradient from the reciprocal sum is split from the real space contribution
!since the virial needs to be taken care of for reciprocal contribution.

  use constants, only : zero, half, one, two, INVSQRTPI, AU_TO_KCAL, BOHRS_TO_A, TWOPI
  use nblist, only : recip
  use qmmm_module, only : qmmm_struct, qmmm_mpi, qm2_rij_eqns, qmewald, qm2_params, qmmm_nml
  implicit none
#include "qm2_array_locations.h"

!Passed in
  integer, intent(in) :: natom
  _REAL_, intent(in) :: qm_coords(3,*) !Coordinates of the QM atoms
  _REAL_, intent(in) :: qm_xcrd(4,qmmm_struct%qm_mm_pairs)  !Wrapped qm coordinates, in same order as the pair list, laid out as
                                 ! x,y,z,charge,x,y,z,charge...
  _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink)
  _REAL_ , intent(inout) :: dxyzqm(3,qmmm_struct%nquant_nlink)
  _REAL_ , intent(inout) :: dxyzcl(*) !Warning only gets force for MM atoms in the QM list!
  _REAL_, intent(in) :: dvectors(3,qmmm_mpi%totkq_count)
  _REAL_, intent(in) :: mmcharge(natom) !Sander's MM charge array but in electron units.
  logical, intent(in) :: qmmmrijincore

!Local
  _REAL_ :: qmmulik_chg, df_qmmma, oneRIJ2, rij_coeff, r2, onerij, rij
  _REAL_ :: qmi(3), df_qmmm(3), erfcx, d_erfcx, fqm(3)
  _REAL_ :: vec(3)
  _REAL_ :: ccfk(3)
  _REAL_ :: ktgs(8), x_cos, y_cos, z_cos, x_sin, y_sin, z_sin
  _REAL_ :: x_cosj, y_cosj, z_cosj, x_sinj, y_sinj, z_sinj
  _REAL_ :: ktg1, ktg2, ktg3, df_sum(3)
  _REAL_ :: fda(8), fd(3), mmchg, kappa
  integer :: i, j, ii, loop_count, inner_loop_count
  integer :: jstart, jend

!For ew_coeff
#include "ew_pme_recip.h"

  !======================================
  !    Calculate Real Space Gradient
  !======================================
  !Step 1 do QM atoms with MM atoms in the cutoff list
  !This is to correct for the overcounting in the reciprocal calculation.
  
  if (qmmm_nml%qm_pme) then
    kappa=ew_coeff  ! The Ewald coefficient from MM part should be used if MM PME contribution is used.
  else 
    kappa=qmewald%kappa  !  Otherwise the Ewald coefficient from QM part should be used.
  end if  
    
  if (qmmmrijincore) then
    !Our QM-MM distances have already been calculated and stored
    loop_count = 0

!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
      inner_loop_count = 1
      qmmulik_chg=scf_mchg(i)*AU_TO_KCAL*BOHRS_TO_A
      qmi(1:3) = qm_coords(1:3,i)
      df_sum(1:3) = zero
      do j = 1,qmmm_struct%qm_mm_pairs
        loop_count = loop_count+1
        rij_coeff = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)*kappa
        call erfcfun(rij_coeff,erfcx) !(2/3)exp(rij_coeff^3)*INVSQRTPI
        d_erfcx = two*kappa*exp(-rij_coeff*rij_coeff)*INVSQRTPI
        oneRIJ = qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        oneRIJ2 = oneRIJ*oneRIJ
        vec(1:3) = qmi(1:3) - qm_xcrd(1:3,j)
        df_qmmma = qm_xcrd(4,j) * (-d_erfcx+(one-erfcx)*oneRIJ)*oneRIJ2
        df_qmmm(1:3) = vec(1:3)*df_qmmma*qmmulik_chg
        df_sum(1:3) = df_sum(1:3) + df_qmmm(1:3)
        dxyzcl(inner_loop_count) = dxyzcl(inner_loop_count) - df_qmmm(1)
        dxyzcl(inner_loop_count+1) = dxyzcl(inner_loop_count+1) - df_qmmm(2)
        dxyzcl(inner_loop_count+2) = dxyzcl(inner_loop_count+2) - df_qmmm(3)
        inner_loop_count = inner_loop_count + 3
      end do !j=1,qmmm_struct%qm_mm_pairs
      dxyzqm(1:3,i) = dxyzqm(1:3,i) + df_sum(1:3)
    end do !i=1, nquant
  else
    !We have to calculate the distance between QM-MM pairs on the fly
!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
      inner_loop_count = 1
      qmmulik_chg=scf_mchg(i)*AU_TO_KCAL*BOHRS_TO_A
      qmi(1:3) = qm_coords(1:3,i)
      df_sum(1:3) = zero
      do j = 1, qmmm_struct%qm_mm_pairs
        vec(1:3) = qmi(1:3) - qm_xcrd(1:3,j)
        r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
        oneRIJ = one/sqrt(r2)
        oneRIJ2 = oneRIJ*oneRIJ
        RIJ = r2*oneRIJ !one/oneRIJ
        rij_coeff = rij*kappa
        call erfcfun(rij_coeff,erfcx)
        d_erfcx = two*kappa*exp(-rij_coeff*rij_coeff)*INVSQRTPI
        df_qmmma = qmmulik_chg * qm_xcrd(4,j) * (-d_erfcx+(one-erfcx)*oneRIJ)*oneRIJ2
        df_qmmm(1:3) = vec(1:3)*df_qmmma
        df_sum(1:3) = df_sum(1:3) + df_qmmm(1:3)
        dxyzcl(inner_loop_count) = dxyzcl(inner_loop_count) - df_qmmm(1)
        dxyzcl(inner_loop_count+1) = dxyzcl(inner_loop_count+1) - df_qmmm(2)
        dxyzcl(inner_loop_count+2) = dxyzcl(inner_loop_count+2) - df_qmmm(3)
        inner_loop_count = inner_loop_count + 3
      end do !j=1,qmmm_struct%qm_mm_pairs
      dxyzqm(1:3,i) = dxyzqm(1:3,i) + df_sum(1:3)
    end do !i = 1, nquant
  end if

  !Step 2 - do all Real space QM atoms with QM atoms
  !We have to do this on the fly because of link atoms which we don't
  !include here but were included in original RIJ calculation.
  
  ! now kappa is the Ewald coefficient from QM part 
  kappa=qmewald%kappa
   
#ifdef MPI
  do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
    jstart =  qmmm_mpi%nquant_nlink_jrange(1,i)
    jend = qmmm_mpi%nquant_nlink_jrange(2,i)
#else
  do I=2,qmmm_struct%nquant_nlink
    jstart = 1
    jend = i-1
#endif
    qmi(1:3) = qm_coords(1:3,i)
    df_sum(1:3) = zero
    qmmulik_chg = scf_mchg(i)*AU_TO_KCAL*BOHRS_TO_A
    do j = jstart, jend
      vec(1:3) = qmi(1:3) - qm_coords(1:3,j)
      r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
      onerij = one / sqrt(r2)
      oneRIJ2 = oneRIJ*oneRIJ
      rij = r2*oneRIJ !one/oneRIJ
      rij_coeff = rij * kappa
      call erfcfun(rij_coeff, erfcx)
      d_erfcx = two*kappa*exp(-rij_coeff*rij_coeff)*INVSQRTPI
      !Reuse qmmm scalars for convenience
      df_qmmma = qmmulik_chg*scf_mchg(j)*(-d_erfcx+(one-erfcx)*oneRIJ)*oneRIJ2
      df_qmmm(1:3) = vec(1:3)*df_qmmma
      df_sum(1:3) = df_sum(1:3) + df_qmmm(1:3)
      dxyzqm(1:3,j) = dxyzqm(1:3,j) - df_qmmm(1:3)
    end do
    dxyzqm(1:3,i) = dxyzqm(1:3,i) + df_sum(1:3)
  end do
  !======================================
  !  End Calculate Real Space Gradient
  !======================================

  !If qm_pme we only want to do the QM-QM k space gradients here.

  !======================================
  !      Calculate K Space Gradient
  !======================================

  if (.not. qmmm_nml%qm_pme) then
    qmewald%d_ewald_mm = zero !Natom array of force on mm atoms due to QM ewald field
    loop_count = 0
!    do ii = 1, totkq
    do ii = qmmm_mpi%kvec_start, qmmm_mpi%kvec_end
      loop_count = loop_count+1
      ccfk(1:3) = dvectors(1:3,loop_count)
  !1) QM-MM Interaction
      ktgs(1:8) = zero
      do j = 1, natom
        !Do this over all MM atoms (skip QM)
!       if (.not. qm_atom_mask(j)) then - often quicker to just loop through all atoms, mmchg will be zero for qm atoms.
        !it is an MM atom
          x_cos =qmewald%ktable(1,j,loop_count)
          x_sin =qmewald%ktable(2,j,loop_count)
          y_cos =qmewald%ktable(3,j,loop_count)
          y_sin =qmewald%ktable(4,j,loop_count)
          z_cos =qmewald%ktable(5,j,loop_count)
          z_sin =qmewald%ktable(6,j,loop_count)
          mmchg = mmcharge(j)

          ktgs(1) = ktgs(1) + mmchg * x_cos*y_cos*z_cos
          ktgs(2) = ktgs(2) + mmchg * x_sin*y_cos*z_cos
          ktgs(3) = ktgs(3) + mmchg * x_cos*y_sin*z_cos
          ktgs(4) = ktgs(4) + mmchg * x_sin*y_sin*z_cos
          ktgs(5) = ktgs(5) + mmchg * x_cos*y_cos*z_sin
          ktgs(6) = ktgs(6) + mmchg * x_sin*y_cos*z_sin
          ktgs(7) = ktgs(7) + mmchg * x_cos*y_sin*z_sin
          ktgs(8) = ktgs(8) + mmchg * x_sin*y_sin*z_sin
!         end if
      end do
      do j = 1, qmmm_struct%nquant_nlink
        x_cos =qmewald%qmktable(1,j,loop_count)
        x_sin =qmewald%qmktable(2,j,loop_count)
        y_cos =qmewald%qmktable(3,j,loop_count)
        y_sin =qmewald%qmktable(4,j,loop_count)
        z_cos =qmewald%qmktable(5,j,loop_count)
        z_sin =qmewald%qmktable(6,j,loop_count)
  !Temporary arrays
        fda(1) = x_sin*y_cos*z_cos
        fda(2) = x_cos*y_cos*z_cos
        fda(3) = x_sin*y_sin*z_cos
        fda(4) = x_cos*y_sin*z_cos
        fda(5) = x_sin*y_cos*z_sin
        fda(6) = x_cos*y_cos*z_sin
        fda(7) = x_sin*y_sin*z_sin
        fda(8) = x_cos*y_sin*z_sin
  !Force on x-axis
        fd(1) = -ktgs(1)*fda(1)+ktgs(2)*fda(2) &
              -ktgs(3)*fda(3)+ktgs(4)*fda(4) &
              -ktgs(5)*fda(5)+ktgs(6)*fda(6) &
              -ktgs(7)*fda(7)+ktgs(8)*fda(8)
  !Force on y-axis
        fd(2) = -ktgs(1)*fda(4)-ktgs(2)*fda(3) &
              +ktgs(3)*fda(2)+ktgs(4)*fda(1) &
              -ktgs(5)*fda(8)-ktgs(6)*fda(7) &
              +ktgs(7)*fda(6)+ktgs(8)*fda(5)
  !Force on z-axis
        fd(3) = -ktgs(1)*fda(6)-ktgs(2)*fda(5) &
              -ktgs(3)*fda(8)-ktgs(4)*fda(7) &
              +ktgs(5)*fda(2)+ktgs(6)*fda(1) &
              +ktgs(7)*fda(4)+ktgs(8)*fda(3)
  
        qmmulik_chg=scf_mchg(j)*AU_TO_KCAL*BOHRS_TO_A
        dxyzqm(1:3,j) = dxyzqm(1:3,j) + ccfk(1:3)*fd(1:3)*qmmulik_chg
      end do
  !2) MM-QM interaction
      fda(1:8) = zero
      do j = 1, qmmm_struct%nquant_nlink
        x_cos =qmewald%qmktable(1,j,loop_count)
        x_sin =qmewald%qmktable(2,j,loop_count)
        y_cos =qmewald%qmktable(3,j,loop_count)
        y_sin =qmewald%qmktable(4,j,loop_count)
        z_cos =qmewald%qmktable(5,j,loop_count)
        z_sin =qmewald%qmktable(6,j,loop_count)

        qmmulik_chg=scf_mchg(j)*AU_TO_KCAL*BOHRS_TO_A  
        fda(1) = fda(1) + qmmulik_chg*x_sin*y_cos*z_cos
        fda(2) = fda(2) + qmmulik_chg*x_cos*y_cos*z_cos
        fda(3) = fda(3) + qmmulik_chg*x_sin*y_sin*z_cos
        fda(4) = fda(4) + qmmulik_chg*x_cos*y_sin*z_cos
        fda(5) = fda(5) + qmmulik_chg*x_sin*y_cos*z_sin
        fda(6) = fda(6) + qmmulik_chg*x_cos*y_cos*z_sin
        fda(7) = fda(7) + qmmulik_chg*x_sin*y_sin*z_sin
        fda(8) = fda(8) + qmmulik_chg*x_cos*y_sin*z_sin
      end do
      do j = 1, natom
        x_cos =qmewald%ktable(1,j,loop_count)
        x_sin =qmewald%ktable(2,j,loop_count)
        y_cos =qmewald%ktable(3,j,loop_count)
        y_sin =qmewald%ktable(4,j,loop_count)
        z_cos =qmewald%ktable(5,j,loop_count)
        z_sin =qmewald%ktable(6,j,loop_count)
              
        ktgs(1) = x_cos*y_cos*z_cos 
        ktgs(2) = x_sin*y_cos*z_cos 
        ktgs(3) = x_cos*y_sin*z_cos 
        ktgs(4) = x_sin*y_sin*z_cos 
        ktgs(5) = x_cos*y_cos*z_sin 
        ktgs(6) = x_sin*y_cos*z_sin
        ktgs(7) = x_cos*y_sin*z_sin
        ktgs(8) = x_sin*y_sin*z_sin
  !Force on x-axis
        fd(1) = -ktgs(1)*fda(1)+ktgs(2)*fda(2) &
              -ktgs(3)*fda(3)+ktgs(4)*fda(4) &
              -ktgs(5)*fda(5)+ktgs(6)*fda(6) &
              -ktgs(7)*fda(7)+ktgs(8)*fda(8)
  !Force on y-axis
        fd(2) = -ktgs(1)*fda(4)-ktgs(2)*fda(3) &
              +ktgs(3)*fda(2)+ktgs(4)*fda(1) &
              -ktgs(5)*fda(8)-ktgs(6)*fda(7) &
              +ktgs(7)*fda(6)+ktgs(8)*fda(5)
  !Force on z-axis
        fd(3) = -ktgs(1)*fda(6)-ktgs(2)*fda(5) &
              -ktgs(3)*fda(8)-ktgs(4)*fda(7) &
              +ktgs(5)*fda(2)+ktgs(6)*fda(1) &
              +ktgs(7)*fda(4)+ktgs(8)*fda(3)
  
  !qm muliken charge has been multiplied already
        mmchg = mmcharge(j)
        qmewald%d_ewald_mm(1:3,j) = qmewald%d_ewald_mm(1:3,j) - ccfk(1:3)*fd(1:3)*mmchg
      end do
    end do
  end if !(if .not. qmmm_nml%qm_pme)

!3) QM-QM interaction - Do this for both QM_PME and QM_EWALD
  loop_count = 0
  do ii = qmmm_mpi%kvec_start,qmmm_mpi%kvec_end
    loop_count=loop_count+1
    ccfk(1:3) = dvectors(1:3,loop_count)
    do i = 1, qmmm_struct%nquant_nlink
      x_cos =qmewald%qmktable(1,i,loop_count)
      x_sin =qmewald%qmktable(2,i,loop_count)
      y_cos =qmewald%qmktable(3,i,loop_count)
      y_sin =qmewald%qmktable(4,i,loop_count)
      z_cos =qmewald%qmktable(5,i,loop_count)
      z_sin =qmewald%qmktable(6,i,loop_count)
      fqm(1:3) = zero
      qmmulik_chg=scf_mchg(i)*AU_TO_KCAL*BOHRS_TO_A
      do j = 1,qmmm_struct%nquant_nlink
        x_cosj =qmewald%qmktable(1,j,loop_count)
        x_sinj =qmewald%qmktable(2,j,loop_count)
        y_cosj =qmewald%qmktable(3,j,loop_count)
        y_sinj =qmewald%qmktable(4,j,loop_count)
        z_cosj =qmewald%qmktable(5,j,loop_count)
        z_sinj =qmewald%qmktable(6,j,loop_count)
        ktg1 = x_cosj*x_cos + x_sinj*x_sin
        ktg2 = y_cosj*y_cos + y_sinj*y_sin
        ktg3 = z_cosj*z_cos + z_sinj*z_sin
        fd(1) = -x_cosj*x_sin + x_sinj*x_cos
        fd(2) = -y_cosj*y_sin + y_sinj*y_cos
        fd(3) = -z_cosj*z_sin + z_sinj*z_cos
!Force on x,y,z-axis - temporary use of FDA array and mmchg array
        mmchg=scf_mchg(j)*qmmulik_chg !scf mulliken charge on qm atom j * atom i (adjusted to kcals)
        fda(1) = mmchg*ccfk(1)*fd(1)*ktg2*ktg3
        fda(2) = mmchg*ccfk(2)*ktg1*fd(2)*ktg3
        fda(3) = mmchg*ccfk(3)*ktg1*ktg2*fd(3)
        if (i == j) then
          fqm(1:3) = fqm(1:3) + fda(1:3)
        else
          fqm(1:3) = fqm(1:3) + half*fda(1:3)
          dxyzqm(1:3,j) = dxyzqm(1:3,j) - half*fda(1:3)
        end if
      end do
      dxyzqm(1:3,i) = dxyzqm(1:3,i) + fqm(1:3)
    end do
  end do

  !======================================
  !    End Calculate K Space Gradient
  !======================================

  return

end subroutine qm_ewald_get_forces



