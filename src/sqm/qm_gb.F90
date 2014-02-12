! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"
!---------------------------------------------------------
! Code for doing QMMM Generalised Born Implicit Solvent
!
! Written by Ross Walker (TSRI 2005)
!----------------------------------------------------------

subroutine allocate_qmgb(nquant_nlink)

  use qmmm_module, only : qmmm_nml, qmmm_struct, qm_gb
  implicit none

  integer, intent(in) :: nquant_nlink

  integer :: ier=0

  allocate (qm_gb%gb_mmpot(nquant_nlink), stat=ier )
  REQUIRE(ier == 0) !Dealocated in deallocate qmmm
  allocate (qm_gb%gb_qmpot(nquant_nlink), stat=ier ) !qmpot = QM-QM
  REQUIRE(ier == 0) !Dealocated in deallocate qmmm

  allocate ( qm_gb%qmqm_onefij(nquant_nlink*nquant_nlink), stat = ier )
  REQUIRE(ier == 0) !Deallocated in deallocate qmmm
  if (qm_gb%saltcon_on) then
    allocate ( qm_gb%qmqm_kappafij(nquant_nlink*nquant_nlink), stat = ier )
    REQUIRE(ier == 0) !Deallocated in deallocate qmmm
  end if
  allocate ( qm_gb%qmqm_gb_list(1+nquant_nlink,nquant_nlink), stat = ier )
  !+1 because first index stores the number of loops.
  REQUIRE(ier == 0) !Deallocated in deallocate qmmm

  return

end subroutine allocate_qmgb

subroutine qmgb_calc_qmqm_onefij(nquant_nlink, qmqm_onefij, iqmatoms, born_radii, one_born_radii, qm_coords )

!Calculates GB Fij term for each QM-QM pair. This is required
!on every step of the SCF but only depends on the distance between
!pairs so can be calculated outside of the SCF and stored in memory.

!Assumes the effective Born radii have been calculated.

  use qmmm_module, only : qm_gb, qmmm_mpi
  use constants, only : fourth
  implicit none

#include "qm2_array_locations.h"

!Passed in
  integer, intent(in) :: nquant_nlink
  integer, intent(in) :: iqmatoms(nquant_nlink) !Index into born_radii - amber atom number of QM atom
  _REAL_, intent(out) :: qmqm_onefij(*)
  _REAL_, intent(in) :: born_radii(*), one_born_radii(*) !Indexed for entire natoms
  _REAL_, intent(in) :: qm_coords(3,nquant_nlink) 

!Local
  _REAL_ aij2, one4aij2, rij2, xqm, yqm, zqm, vec1, vec2, vec3, exptmp, alphai
  _REAL_ alphai_inv
  integer i, j, loop_count, inner_loop_count, qmi, qmj

  loop_count = 0
!  do i=1, nquant_nlink
  do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
    xqm = qm_coords(1,i)
    yqm = qm_coords(2,i)
    zqm = qm_coords(3,i)
    qmi = iqmatoms(i)
    alphai = born_radii(qmi)
    alphai_inv = one_born_radii(qmi)
    inner_loop_count=0
    do j=1, nquant_nlink
      vec1 = xqm-qm_coords(1,j)
      vec2 = yqm-qm_coords(2,j)
      vec3 = zqm-qm_coords(3,j)
      rij2 = vec1*vec1+vec2*vec2+vec3*vec3
      if ( rij2 <= qm_gb%mmcut2 ) then
        qmj = iqmatoms(j)
        one4aij2 = fourth*alphai_inv*one_born_radii(qmj)
        exptmp = exp(-rij2*one4aij2)

        aij2 = alphai * born_radii(qmj)
      
        loop_count = loop_count+1
        inner_loop_count = inner_loop_count+1
        qmqm_onefij(loop_count) = rij2+aij2*exptmp
        qm_gb%qmqm_gb_list(inner_loop_count+1,i) = j !Index of which quantum atoms interact with current quantum atom.
      end if
    end do
    qm_gb%qmqm_gb_list(1,i) = inner_loop_count !First index stores the total loops for QM i
  end do

  call vdinvsqrt( loop_count, qmqm_onefij, qmqm_onefij )

  if (qm_gb%alpb_on) then
    qmqm_onefij(1:loop_count) = qmqm_onefij(1:loop_count)+qm_gb%one_Arad_beta
  end if

  !Adjust for kappa, from saltcon
  if ( qm_gb%saltcon_on ) then
    !Calculate exp(-kappa*fij)
    call vdinv(loop_count,qmqm_onefij,qm_gb%qmqm_kappafij)
    qm_gb%qmqm_kappafij(1:loop_count) = -qm_gb%kappa*qm_gb%qmqm_kappafij(1:loop_count)
    call vdexp( loop_count, qm_gb%qmqm_kappafij, qm_gb%qmqm_kappafij )
  end if

  return

end subroutine qmgb_calc_qmqm_onefij

subroutine qmgb_calc_mm_pot(natom,gb_mmpot,qm_atom_mask,scaled_mm_charges, &
                            real_scratch1,real_scratch2, int_scratch1, &
                            qm_coords,mm_coords,born_radii,one_born_radii, iqmatoms)

!Calculates the GB potential at each QM atom due to all the MM atoms.
!                 1     1         q(mm[j])
! gb_mmpot(i) = (--- - ---)  sum(---------)
!                 ei    e0 1->MM[j] f(i,j)
!
! Units are q in electrons, f(i,j) in angstroms
!
! The values stored in gb_mmpot(i) are calculated outside of the SCF and then converted to Au/Bohrs and added
! to the diagonal elements of the Fock matrix of each QM atom (i).

! Needs two real scratch arrays of at least natom long.
  use qmmm_module, only : qm_gb, qmmm_mpi, qmmm_struct, qmmm_nml
  use constants, only : fourth
  implicit none

! Passed in
  integer, intent(in) :: natom
  _REAL_, intent(out) :: gb_mmpot(*) !nquant_nlink
  _REAL_, intent(in) :: scaled_mm_charges(natom)
  _REAL_, intent(out) :: real_scratch1(natom), real_scratch2(natom)
  integer, intent(out) :: int_scratch1(natom)
  _REAL_, intent(in) :: qm_coords(3,*), mm_coords(3,natom)
  _REAL_, intent(in) :: born_radii(natom), one_born_radii(natom)
  logical, intent(in) :: qm_atom_mask(natom)
  integer, intent(in) :: iqmatoms(*) !Index into born_radii - amber atom number of QM atom

! Local
  integer i, j, loop_count, atomnum, qmi
  _REAL_ xqm, yqm, zqm, alpha_qm, alpha_qm_inv, rij2, aij2, one4aij2
  _REAL_ vec1, vec2, vec3, exptmp, temp_pot, diel_fac

  diel_fac = qm_gb%intdieli - qm_gb%extdieli !Recalculated in inner loop if saltcon_on

!  do i = 1, nquant_nlink
  do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
    xqm = qm_coords(1,i)
    yqm = qm_coords(2,i)
    zqm = qm_coords(3,i)
    qmi = iqmatoms(i)
    alpha_qm = born_radii(qmi)
    alpha_qm_inv = one_born_radii(qmi)

    loop_count = 0 
    !we need to skip both QM atoms and MM link pair atoms
    do j = 1, natom
      !Do this over MM atoms only, skip QM and MM link pair atoms.
      if (.not. qm_atom_mask(j) .and. .not. qmmm_struct%mm_link_mask(j)) then
        !Calculate one_fij values
        vec1 = xqm-mm_coords(1,j) !We don't need to worry about adjusting coords here for periodic
        vec2 = yqm-mm_coords(2,j) !boundaries since GB does not work with periodic boundaries.
        vec3 = zqm-mm_coords(3,j) !we also don't need to worry about the mm_coords array containing
                                  !the MMlink pair coords rather than the link atom coords since
                                  !we skip MMlink pairs anyway.
        rij2 = vec1*vec1+vec2*vec2+vec3*vec3
        if (rij2 <= qm_gb%mmcut2) then
          one4aij2 = fourth*alpha_qm_inv*one_born_radii(j)
          exptmp = exp(-rij2*one4aij2)
          aij2 = alpha_qm * born_radii(j)
          loop_count = loop_count + 1
          real_scratch1(loop_count) = rij2+aij2*exptmp
          int_scratch1(loop_count) = j
        end if
      end if
    end do
    !Now replace real_scratch with 1/sqrt(real_scratch1)
    call vdinvsqrt(loop_count,real_scratch1,real_scratch1)

    !Calculate kappa dependency if required
    if ( qm_gb%saltcon_on ) then
      call vdinv(loop_count,real_scratch1,real_scratch2)
      real_scratch2(1:loop_count) = -qm_gb%kappa*real_scratch2(1:loop_count)
      call vdexp( loop_count, real_scratch2, real_scratch2 )
    end if
    
    !Real Scratch 1 contains 1.0d0/fij values
    !Real Scratch 2 contains exp(-kappa*fij)
    !int Scratch1 contains the j values (atom numbers) that we are including.
    if (qm_gb%alpb_on) then
      real_scratch1(1:loop_count) = real_scratch1(1:loop_count) + qm_gb%one_Arad_beta
    end if

    !Now we calculate sum(qj/fij)
    temp_pot = 0.0d0
    if (qm_gb%saltcon_on) then
      do j = 1, loop_count
        atomnum = int_scratch1(j)
        diel_fac = qm_gb%intdieli - qm_gb%extdieli * real_scratch2(j)
        temp_pot = temp_pot + diel_fac*scaled_mm_charges(atomnum) * real_scratch1(j)
      end do
      gb_mmpot(i) = temp_pot !Potential for this QM atom.
    else
      do j = 1, loop_count
        atomnum = int_scratch1(j)
        temp_pot = temp_pot + scaled_mm_charges(atomnum) * real_scratch1(j)
      end do
      gb_mmpot(i) = diel_fac*temp_pot !Potential for this QM atom.
    end if
  end do !i=1, nquant_nlink


  return

end subroutine qmgb_calc_mm_pot


subroutine qmgb_calc_qm_pot(gb_qmpot,qmqm_onefij, scf_mchg )
!Calculates the GB potential at each QM atom due to all the QM atoms.
!                 1     1        q(qm)
! gb_qmpot(i) = (--- - ---) sum(-------)
!                 ei    e0 1->QM f(i,j)
!
! Units are q in electrons, f(i,j) in angstroms
!

  use qmmm_module, only : qm_gb, qmmm_mpi
  implicit none

! Passed in
  _REAL_, intent(out) :: gb_qmpot(*) !nquant_nlink
  _REAL_, intent(in) :: qmqm_onefij(*)
  _REAL_, intent(in) :: scf_mchg(*) !current Mulliken charges on QM atoms

! Local
  _REAL_ temp_pot, diel_fac
  integer i, j, loop_count, inner_loop_end, atomnum

!Note, in parallel this routine expects the same memory layout
!between cpus as the calculation of fij values.
  loop_count = 0 
  if ( qm_gb%saltcon_on ) then
!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end !1 to nquant_nlink
      temp_pot = 0.0d0
      inner_loop_end = qm_gb%qmqm_gb_list(1,i)+1
      do j = 2, inner_loop_end !From 2 since element 1 of the qmqm_gb_list contains the no. interactions for i.
        atomnum = qm_gb%qmqm_gb_list(j,i)
        loop_count = loop_count + 1
        diel_fac = qm_gb%intdieli - qm_gb%extdieli * qm_gb%qmqm_kappafij(loop_count)
        temp_pot = temp_pot + diel_fac * scf_mchg(atomnum) * qmqm_onefij(loop_count)
      end do
      gb_qmpot(i) = temp_pot
    end do !i=1, nquant
  else
    diel_fac = qm_gb%intdieli - qm_gb%extdieli
!    do i = 1, nquant_nlink
    do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
      temp_pot = 0.0d0
      inner_loop_end = qm_gb%qmqm_gb_list(1,i)+1
      do j = 2, inner_loop_end !From 2 since element 1 of the qmqm_gb_list contains the no. interactions for i.
        atomnum = qm_gb%qmqm_gb_list(j,i)
        loop_count = loop_count + 1
        temp_pot = temp_pot + diel_fac * scf_mchg(atomnum) * qmqm_onefij(loop_count)
      end do
      gb_qmpot(i) = temp_pot
    end do !i=1, nquant
  end if

  return

end subroutine qmgb_calc_qm_pot

subroutine qmgb_add_fock(loop_extent,fock_matrix, gb_mmpot, gb_qmpot)
!Author: Ross Walker, TSRI 2005

!Adds GB potential info (mmpot and qmpot) to the diagonal elements of the fock matrix

  use constants, only : AU_TO_EV, BOHRS_TO_A
  use qmmm_module, only : qm2_params
  implicit none

!Passed in
  integer, intent(in) :: loop_extent
  _REAL_, intent(in) :: gb_mmpot(loop_extent) !Potential at each QM atom due to MM GB
  _REAL_, intent(in) :: gb_qmpot(loop_extent) !Potential at each QM atom due to QM GB
  _REAL_, intent(inout) :: fock_matrix(*) !Fock matrix

!Local
  integer :: i, ia, ib, i1, i2
  _REAL_ :: temp_pot

  !Now add the mmpot and qmpot array contributions to the diagonal elements of the fock matrix
  do i = 1, loop_extent
    IA = qm2_params%orb_loc(1,I)
    IB = qm2_params%orb_loc(2,I)
    temp_pot=(gb_mmpot(i)+gb_qmpot(i))*AU_TO_EV*BOHRS_TO_A
    do I1 = IA,IB
       i2 = qm2_params%pascal_tri2(i1)
       fock_matrix(i2) = fock_matrix(i2) + temp_pot
                         !AU_TO_EV*BOHRS_TO_A converts (electrons/angstrom) to (eV/Bohr)
    end do
  end do

  return

end subroutine qmgb_add_fock


