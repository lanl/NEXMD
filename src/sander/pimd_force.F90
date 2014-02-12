#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------------------------
subroutine pimd_part_spring_force( x, f, mass, Epot_spring, Epot_deriv, dvdl )
! computes the harmonic-spring energy and forces arising from discretized
! path integral action.
! computes the kinetic term of the virial estimator for the internal energy
! Computes the estimators for the thermonynamic integration w.r.t. mass.
!------------------------------------------------------------------------------
  use constants, only: hbar, kB
  use pimd_vars, only: Eimp_virial, x_centroid, nbead, dmdlm, itimass
#ifdef MPI
#  ifdef LES
  use evb_pimd,  only: natomCL, vel0_bead, bead_dcrypt
  use miller,    only: dlnCdd_dl, i_qi, gradRC, div_ndx
  use evb_parm, only: nbias
  use evb_data,  only: f_v, F_QI, G_QI
#  endif 
#endif

  implicit none

#include "md.h"
#include "les.h"
#include "memory.h"

  _REAL_, intent(in) :: x( 3, natom ) !! coordinates
  _REAL_, intent(in) :: mass(natom)
  _REAL_ :: f( 3, natom )             !! forces
  _REAL_ :: Epot_spring, Epot_deriv, beta, kT, dvdl
  integer :: idim, iatom, icopy, istart, iend, ierr
  _REAL_ :: coeff, tmp, tmp_x, tmp_f, tmp_e, tmp_ti, sum_e, sum_d, sum_v

  _REAL_ :: coeff_F, coeff_G1, coeff_G2
  _REAL_ :: F_sum1, F_sum2, F_sum3, F_sum4, G_sum, Ftot, Gtot
  _REAL_ :: ftmp, fv, gnorm

  integer :: m, n, mm, i3, mdx

  intrinsic :: sqrt, dot_product

#ifdef MPI
   include 'mpif.h'
  _REAL_ :: ener(4), temp_ener(4) ! JVAN
#  ifdef LES
  _REAL_ :: dx(3*natomCL)
#  endif
#endif

#include "parallel.h"

#ifdef MPI
  istart = iparpt(mytaskid) + 1
  iend   = iparpt(mytaskid+1) 
#else
  istart = 1
  iend   = natom
#endif

  kT = kB * temp0
  beta = 1.0d0 / kT

  coeff = dble( ncopy ) / ( hbar * beta )**2 
  sum_e = 0.d0
  sum_d = 0.d0
  sum_v = 0.d0
  dvdl = 0.d0 ! JVAN

  !For QI.
  coeff_F = 2.d0 / dble( ncopy )
  coeff_G1 = 6.d0 * dble( ncopy ) / beta / beta
  coeff_G2 = 4.d0 * coeff / beta

  do iatom = istart, iend 
     if(cnum(iatom).eq.0) then
        x_centroid(1:3,iatom) = x(1:3,iatom)
     else if( cnum(iatom).eq.1) then
        x_centroid(1:3,iatom) = 0.d0
        do icopy=1,ncopy
           x_centroid(1:3,iatom)=x_centroid(1:3,iatom)+x(1:3,iatom+icopy-1)
        end do
        x_centroid(1:3,iatom) = x_centroid(1:3,iatom)/ncopy
        do icopy=2,ncopy
           x_centroid(1:3,iatom+icopy-1)=x_centroid(1:3,iatom)
        end do
     end if
  end do 

!  +----------------------------------------------------------------------+
!  |  Virial estimator for kinetic energy                                 |
!  |                                                                      |
!  |  K = f/(2B) + 1/P sum(s=1,P)[ 1/2 (x_s - x_c) . dV(x_s)/dx_s ]       |
!  +----------------------------------------------------------------------+
        
  do iatom = istart,iend
     if(cnum(iatom).eq.0) then
        ! Contribution from classical atoms to the estimator of dV/dl 
        ! for TI w.r.t. mass. (Same for thermodynamic and virial estimators.)
        if (itimass > 0) dvdl = dvdl - 1.5d0 * kT * dmdlm(iatom)
        cycle
     else
     do idim    = 1, 3
        tmp = x(idim, iatom) - x_centroid(idim,iatom)
        sum_d = sum_d - tmp * f(idim,iatom)
        sum_v = sum_v - f(idim,iatom)*x(idim,iatom) 
        if ( cnum(iatom) == 1 ) then
           tmp_e = ( x(idim,iatom) - x( idim,iatom+1) )**2
           tmp_f = 2.0d0 * x( idim, iatom ) - x( idim, iatom+1 ) &
                                            - x( idim, iatom+ncopy-1)
        else if ( cnum(iatom) == ncopy ) then
           tmp_e = ( x(idim,iatom) - x(idim,iatom-ncopy+1))**2
           tmp_f = 2.0d0 * x(idim,iatom) - x(idim,iatom-ncopy+1 ) &
                                         - x(idim,iatom-1)
        else
           tmp_e = ( x(idim,iatom) - x(idim,iatom+1) )**2
           tmp_f = 2.0d0 * x(idim,iatom) - x(idim,iatom+1)&
                                         - x(idim,iatom-1)
        endif
        sum_e = sum_e + mass(iatom)*tmp_e

        !--------------------------------------------------------!
        ! Calculation of estimators of dV/dl for TI w.r.t. mass. !
        ! dvdl is computed from Eq. 3.16 or Eq. 3.7  in          !
        ! Vanicek & Miller, JCP, 2007                            !
        !--------------------------------------------------------!
        if (itimass > 0) then
          select case (itimass)
            case (1)   ! virial-like estimator, Eq. 3.16
              tmp_ti = kT/ncopy - f(idim,iatom) * tmp
            case (2)   ! thermodynamic-like estimator, Eq. 3.7
              tmp_ti = kT - coeff * mass(iatom) * tmp_e 
          end select
          dvdl = dvdl - dmdlm(iatom) * 0.5d0 * tmp_ti
        end if
        
        f(idim,iatom) = f(idim,iatom) - coeff*mass(iatom)*tmp_f
     enddo
     end if
  enddo


! KFW (d/dl) C_dd term

#ifdef MPI
#  ifdef LES

   if( i_qi > 0 .and. itimass > 0) then

      dlnCdd_dl = dvdl

      do n = 1, nbias

         gnorm = 0.0d0
         do m = 1, natomCL
            mdx = ( m - 1 ) * 3
            mm = bead_dcrypt( m, div_ndx(n) )
            gnorm = gnorm + ( gradRC(mdx+1,n)**2 + gradRC(mdx+2,n)**2 &
                            + gradRC(mdx+3,n)**2 ) / mass(mm)
         enddo

         do m = 1, natomCL
            mdx = ( m - 1 ) * 3
            mm = bead_dcrypt( m, div_ndx(n) )
            dlnCdd_dl = dlnCdd_dl + kT *  dmdlm(mm) * 0.50d0 * ( gradRC(mdx+1,n)**2 + gradRC(mdx+2,n)**2 &
                                  + gradRC(mdx+3,n)**2 ) / mass(mm) / gnorm
         enddo

      enddo

   endif

#  endif
#endif

! KFW (d/dl) C_dd term

  Epot_spring = 0.5d0 * coeff * sum_e
  Epot_deriv = 0.5d0 * sum_d
  Eimp_virial = sum_v

#ifdef MPI
  ener(1) = Epot_spring
  ener(2) = Epot_deriv
  ener(3) = Eimp_virial
  ener(4) = dvdl
  call mpi_allreduce( ener, temp_ener, 4, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr )
  Epot_spring = temp_ener(1)
  Epot_deriv  = temp_ener(2)
  Eimp_virial = temp_ener(3)
  dvdl = temp_ener(4)
#endif 

end subroutine



!------------------------------------------------------------------------------
subroutine pimd_full_spring_force ( x, f, mass, Epot_spring, Epot_deriv, dvdl )
! computes the harmonic-spring energy and forces arising from discretized
! path integral action.
! computes the kinetic term of the virial estimator for the internal energy
! Computes the estimators for the thermonynamic integration w.r.t. mass.
!------------------------------------------------------------------------------
  use pimd_vars, only: nbead,Eimp_virial,x_centroid, dmdlm, itimass
  use constants, only: hbar, kB   
  use file_io_dat
  use full_pimd_vars, only: xall,mybeadid
  implicit none
#include "md.h"
#include "memory.h"
  _REAL_, intent(in) :: x( 3, natom ) !! coordinates
  _REAL_, intent(in) :: mass(natom)
  _REAL_ :: f( 3, natom )             !! forces
  _REAL_ :: Epot_spring, Epot_deriv, beta, kT, dvdl
  integer :: idim, iatom, istart, iend, st, ierr,ibead,prev_bead,next_bead
  _REAL_ :: coeff, tmp, tmp_x, tmp_f, tmp_e, tmp_ti, sum_e, sum_d, sum_v

#ifdef MPI
   include 'mpif.h'
  _REAL_ :: ener(4), temp_ener(4)
#endif
#include "parallel.h"

  REQUIRE(ng_sequential)

#ifdef MPI

  REQUIRE( allocated(xall) )

  if( sanderrank.eq.0) then
      call mpi_allgather(x   ,3*natom,MPI_DOUBLE_PRECISION,&
                         xall,3*natom,MPI_DOUBLE_PRECISION,&
                         commmaster,ierr)
  endif

  call mpi_bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#endif


  istart = iparpt(mytaskid)+1
  iend   = iparpt(mytaskid+1)


  kT = kB * temp0
  beta = 1.0d0/kT
  coeff = dble( nbead ) / ( hbar * beta )**2 

  sum_e = 0.d0
  sum_d = 0.d0
  sum_v = 0.d0
  dvdl = 0.d0

  x_centroid(1:3,1:natom)=0.0
  do ibead = 1, nbead
  do iatom = istart,iend
     x_centroid(1:3,iatom)=x_centroid(1:3,iatom)+xall(1:3,iatom,ibead)
  end do
  end do

  x_centroid(1:3,1:natom)=x_centroid(1:3,1:natom)/nbead   

  prev_bead = mybeadid-1
  next_bead = mybeadid+1
  
  if(prev_bead.eq.0) prev_bead=nbead
  if(next_bead.eq.nbead+1) next_bead=1

  do iatom = istart,iend
  do idim    = 1, 3
     tmp   = ( x(idim,iatom) - x_centroid(idim,iatom) )
     tmp_e = ( x(idim,iatom) - xall(idim,iatom,prev_bead) )**2
     tmp_f = 2.0d0 * x(idim,iatom) - xall(idim,iatom,prev_bead) &
                                   - xall(idim,iatom,next_bead)
     sum_e = sum_e + mass(iatom) * tmp_e
     sum_d = sum_d - f(idim,iatom) * tmp
     sum_v = sum_v - f(idim,iatom) * x(idim,iatom)

     !--------------------------------------------------------!
     ! Calculation of estimators of dV/dl for TI w.r.t. mass. !
     ! dvdl is computed from Eq. 3.16 or Eq. 3.7  in          !
     ! Vanicek & Miller, JCP, 2007                            !
     !--------------------------------------------------------!
     if (itimass > 0) then
        select case (itimass)
          case (1)   ! virial-like estimator, Eq. 3.16
            tmp_ti = kT/nbead - f(idim,iatom) * tmp
          case (2)   ! thermodynamic-like estimator, Eq. 3.7
            tmp_ti = kT - coeff * mass(iatom) * tmp_e 
        end select
        dvdl = dvdl - dmdlm(iatom) * 0.5d0 * tmp_ti
     end if

     f(idim,iatom) = f(idim,iatom) - coeff * mass(iatom)*tmp_f
  enddo
  enddo

#ifdef MPI
  ener(1) = sum_e
  ener(2) = sum_d
  ener(3) = sum_v
  ener(4) = dvdl
  call mpi_allreduce(ener,temp_ener,4,MPI_DOUBLE_PRECISION,mpi_sum,commworld,ierr)
  sum_e = temp_ener(1)
  sum_d = temp_ener(2)
  sum_v = temp_ener(3)
  dvdl = temp_ener(4)
#endif

  Epot_spring = 0.5d0 * coeff * sum_e
  Epot_deriv = 0.5d0 * sum_d
  Eimp_virial = sum_v 


end subroutine pimd_full_spring_force


subroutine dispatch( ntime, begin_in, step_in, begin_out, step_out, input, output )
   implicit none
   integer ntime, begin_in, step_in, begin_out, step_out
   integer i, id_in, id_out
   _REAL_  input(*), output(*)

   do i=1,ntime
      id_in = begin_in  + (i-1) * step_in
      id_out= begin_out + (i-1) * step_out 
      output( id_out ) = input( id_in )
   end do
end subroutine

subroutine rotate( rotat, n, xc )
   implicit none
   integer i,i3,n
   _REAL_  rotat(3,3), xc(3*n), tmpx, tmpy, tmpz

   do i = 1, n
      i3 = 3*i-3
      tmpx = rotat(1,1)*xc(i3+1) + rotat(1,2)*xc(i3+2) + rotat(1,3)*xc(i3+3)
      tmpy = rotat(2,1)*xc(i3+1) + rotat(2,2)*xc(i3+2) + rotat(2,3)*xc(i3+3)
      tmpz = rotat(3,1)*xc(i3+1) + rotat(3,2)*xc(i3+2) + rotat(3,3)*xc(i3+3)
      xc(i3+1) = tmpx
      xc(i3+2) = tmpy
      xc(i3+3) = tmpz
   end do

end subroutine 

subroutine dispatch3( ntime, begin_in, step_in, begin_out, step_out, input, output )
   implicit none
   integer ntime, begin_in, step_in, begin_out, step_out
   _REAL_ input(*), output(*)

   call dispatch( ntime, 3*begin_in - 2, 3*step_in, 3*begin_out-2, 3*step_out, input, output )
   call dispatch( ntime, 3*begin_in - 1, 3*step_in, 3*begin_out-1, 3*step_out, input, output )
   call dispatch( ntime, 3*begin_in    , 3*step_in, 3*begin_out  , 3*step_out, input, output )
end subroutine dispatch3

subroutine normalize(vector,dim)
   !normalize a vector of components of a dim dimensional vector
   !Very small vectors need to be considered zero for simulated annealing.
   !These are trapped below and not normalized.
   implicit none
    _REAL_ vector(*),vlength
    integer dim, i
   
   vlength = 0.d0
   do i = 1, dim 
      vlength = vlength + vector(i)*vector(i)
   end do
   if (vlength>1.0d-6) then
      vlength = 1.0d0/sqrt(vlength)
      vector(1:dim) = vector(1:dim)*vlength 
   else 
      vector(1:dim) = 0.d0
   end if      

end subroutine normalize

#ifdef MPI

! ENTIRE NEB force calc (full_neb_forces) is embedded in a master-only conditional.
! this avoids all of the coordinate communication between non-masters 
!    as was done in older versions
! forces will be sent to non-master procs in the existing fdist call

subroutine full_neb_forces(mass,x,f,epot,fitgp,rmsgp)

! NEB extensive overhaul by carlos simmerling and christina bergonzo 12-2008
! described in Bergonzo, Campbell, Walker and Simmerling, IJQC 2009

! extensive changes to allow selection of subset of atoms for NEB forces
! allow explicit water simulations- coordinates in MD no longer rotated (neighbors are copied and rotated)
! also improve parallel performance
! NEB now requires MPI and groupfile. LES type NEB was removed.

   use full_pimd_vars, only: mybeadid
   use neb_vars, only: springforce, tangents, xprev, xnext, &
                       neb_force, last_neb_atom, neb_nrg_all, &
                       neb_nbead, next_node, prev_node
   implicit none

#  include "md.h" 
! timers
#  include "def_time.h" 
!need access to skmin, skmax, and tanmode
#  include "memory.h" 
!need access to natom, nattgtrms, etc
include 'mpif.h'
#  include "parallel.h"
  integer ierr
   _REAL_ mass(*)
   logical rmsok
   integer i, j, rep, iatom, iatm3, index, irep, frep,st,jatm3
   _REAL_ eplus,eminus,spring,spring2,ezero,emax
   _REAL_ dotproduct,dotproduct2
   _REAL_ rmsdvalue, rotat(3,3)
   _REAL_ dottmp(1)
 !external function
  ! _REAL_ ddot
   
  !passed in variables:
     !x is the coordinates
     !f is the forces
     !amass is the atomic mass array

   _REAL_ x(*),f(*),epot
   _REAL_ energy,tempdot
   integer rmsgp(*),fitgp(*)

   if (mod(neb_nbead,2)/=0) then
       write (6,*) "must have even number of beads",neb_nbead
       call mexit(6,1)
   endif

!**********************************************************
! note that beadid runs 1 to neb_nbead and MPI nodes run 0 to #nodes-1
!
! note also that next_node and prev_node are circular, 
!  i.e. first and last still have both neighbors defined.
!**********************************************************

   ! these ranks imply we use commmaster for the communicator
   next_node = masterrank+1
   if(next_node>=mastersize) next_node = 0
   prev_node = masterrank-1
   if(prev_node<0) prev_node = mastersize-1
 
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3,a,i3)') "NEBDBG: next_node=",next_node," prev_node=",prev_node
   write(50+worldrank,'(a,i8)') "NEBDBG: last_neb_atom=",last_neb_atom
   write(50+worldrank,'(a,i8)') "NEBDBG: nattgtrms=",nattgtrms
   write(50+worldrank,'(a,i3)') "NEBDBG: neb_nbead=",neb_nbead
#  endif

   ! carlos: only send/receive neighbor coordinates!
   ! carlos: beads 1 and neb_nbead should just exit! (no coordinate update, and no neighbor)

      if (mod(mybeadid,2)==0) then  !even, exchange with bead+1

         call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,101010, &
                    xnext,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,101010, &
                    commmaster,st,ierr)

      else   !odd, exchange with bead-1

         call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,101010, &
                    xprev,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,101010, &
                    commmaster,st,ierr)
      endif

      ! even send coords to rep-1 


      if (mod(mybeadid,2)==0) then  !even, exchange with bead-1

         call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,202020, &
                    xprev,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,202020, &
                    commmaster,st,ierr)

      else   !odd, exchange with bead+1 (except bead neb_nbead)

         call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,202020, &
                    xnext,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,202020, &
                    commmaster,st,ierr)

      endif

#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: epot=",epot
#  endif

      call mpi_allgather(epot,1,MPI_DOUBLE_PRECISION, &
                         neb_nrg_all,1, MPI_DOUBLE_PRECISION, &
                         commmaster,ierr)

! fit neighbor coords to self
! ok to fit both neighbors since next/prev were defined as circular

! note that we pass fitgp and rmsgp, which were arguments to this routine in the call from force()
! first argument is refc (which gets fit to the second argument, the actual MD coords)
! refc in this case is the neighbor image

        call rmsfit( xprev, x, mass, fitgp,rmsgp, rmsdvalue, nattgtrms, nattgtfit, rmsok) 
        call rmsfit( xnext, x, mass, fitgp,rmsgp, rmsdvalue, nattgtrms, nattgtfit, rmsok) 

   !ezero is the higher energy of the two end points
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3)') "NEBDBG: Getting max (ezero) between 1 and ",neb_nbead
#  endif
   ezero = max( neb_nrg_all(1), neb_nrg_all(neb_nbead) )
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: ezero=",ezero
   write(50+worldrank,'(a,f12.4)') "NEBDBG:   neb_nrg_all(1)=",neb_nrg_all(1)
#  endif
   
   !emax is the highest energy of all the points
   emax = neb_nrg_all(1)
   do rep = 2, neb_nbead
#  ifdef NEBDEBUG
      write(50+worldrank,'(a,i3,a,f12.4)') "NEBDBG:   neb_nrg_all(",rep,")=",neb_nrg_all(rep)
#  endif
      emax = max(emax, neb_nrg_all(rep))
   end do
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: emax=",emax
#  endif
   
   !spring2 is the spring constant of the second point in path
   if (skmin.EQ.skmax) then
      spring2 = skmax   
   else if (neb_nrg_all(2)>ezero) then
      spring2 = skmax - skmin*((emax-max(neb_nrg_all(1),neb_nrg_all(2)))/  &
          (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

   energy = 0.d0

   rep = mybeadid

   !nebrms = 0.d0

! don't bother with forces for first or last bead (forces are zeroed anyway)
! DAN ROE: Zero out the forces on all beads here since first and last beads
!          are about to exit. Some compilers do NOT automatically zero
!          arrays. Call in force.f only uses nattgtrms*3 spots
   do i=1,nattgtrms*3
     neb_force(i)=0.d0
   enddo

! CHANGE LATER IF WE DON'T DO (FAKE) MD ON ENDPOINTS!
   if( mybeadid==1.or.mybeadid==neb_nbead) then
        !write (6,*) "I am first or last bead, returning "
#  ifdef NEBDEBUG
      write(worldrank+50,'(a)') "-------------------------------------------------"
#  endif
        return
   endif
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3,a,i3,a,i3)') "NEBDBG: Rep=",rep," rep+1=",rep+1," rep-1=",rep-1
#  endif
   !calculate spring constant for rep and rep+1
 
   spring = spring2
   if (skmin.EQ.skmax) then
      spring2 = skmax
   else if (neb_nrg_all(rep+1)>ezero.AND.emax/=ezero) then
      spring2 = skmax - skmin*((emax-max(neb_nrg_all(rep),neb_nrg_all(rep-1)))/ &
                (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

   !tangents = 0.0 
! DAN NOTE: tangents is set to 0.0 here, also in the following do loop. Is this
! kind of assignment ok? 

! carlos: how are rep+1 for last rep and rep-1 for first handled as indices??
! for now it doesn't matter since first and last rep exited already

! only fill the tangent elements for atoms in the rmsmask. 
! this way the normalization will be ok, and we get to keep the 
! tangent pointer and atom coordinate pointers in sync to make the 
! code easier to follow as compared to a shorter vector for the 
! tangent for just the rmsgroup atoms

! init tangent array since only some elements will be filled

   do i=1,last_neb_atom*3
      tangents(i)=0.d0
   enddo

   if (tmode.EQ.1) then
      !calculate the tangents (for all images except first and last)

      if (neb_nrg_all(rep+1)>neb_nrg_all(rep).AND.neb_nrg_all(rep)>neb_nrg_all(rep-1)) then
         do j=1,nattgtrms
            iatom=rmsgp(j)
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = xnext(iatm3+1) - x(iatm3+1)
            tangents(iatm3+2) = xnext(iatm3+2) - x(iatm3+2)
            tangents(iatm3+3) = xnext(iatm3+3) - x(iatm3+3)
         end do
      else if (neb_nrg_all(rep+1)<neb_nrg_all(rep).AND.neb_nrg_all(rep)<neb_nrg_all(rep-1)) then
         do j=1,nattgtrms
            iatom=rmsgp(j)
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) =  x(iatm3+1) - xprev(iatm3+1)
            tangents(iatm3+2) =  x(iatm3+2) - xprev(iatm3+2)
            tangents(iatm3+3) =  x(iatm3+3) - xprev(iatm3+3)
         end do  
      else if (neb_nrg_all(rep+1)>neb_nrg_all(rep-1)) then
         eplus = max( abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
            abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)) )
         eminus = min(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
            abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)) )
         
         do j=1,nattgtrms
            iatom=rmsgp(j)
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = (xnext(iatm3+1) - x(iatm3+1))*eplus
            tangents(iatm3+2) = (xnext(iatm3+2) - x(iatm3+2))*eplus
            tangents(iatm3+3) = (xnext(iatm3+3) - x(iatm3+3))*eplus
            tangents(iatm3+1) = tangents(iatm3+1) + (x(iatm3+1) - xprev(iatm3+1))*eminus
            tangents(iatm3+2) = tangents(iatm3+2) + (x(iatm3+2) - xprev(iatm3+2))*eminus
            tangents(iatm3+3) = tangents(iatm3+3) + (x(iatm3+3) - xprev(iatm3+3))*eminus
         end do  
      else !neb_nrg_all(rep+1)<=neb_nrg_all(rep-1)
         eplus = max(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
            abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)))
         eminus = min(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
            abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)))
         do j=1,nattgtrms
            iatom=rmsgp(j)
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = (xnext(iatm3+1) - x(iatm3+1))*eminus
            tangents(iatm3+2) = (xnext(iatm3+2) - x(iatm3+2))*eminus
            tangents(iatm3+3) = (xnext(iatm3+3) - x(iatm3+3))*eminus
            tangents(iatm3+1) = tangents(iatm3+1) + (x(iatm3+1) - xprev(iatm3+1))*eplus
            tangents(iatm3+2) = tangents(iatm3+2) + (x(iatm3+2) - xprev(iatm3+2))*eplus
            tangents(iatm3+3) = tangents(iatm3+3) + (x(iatm3+3) - xprev(iatm3+3))*eplus
         end do  
      end if 
   else !tmode.NE.1 so use strict tangents definition
      do j=1,nattgtrms
         iatom=rmsgp(j)
         iatm3 = 3*(iatom - 1)
         tangents(iatm3+1) = xnext(iatm3+1) - x(iatm3+1)
         tangents(iatm3+2) = xnext(iatm3+2) - x(iatm3+2)
         tangents(iatm3+3) = xnext(iatm3+3) - x(iatm3+3)
      end do
   end if

#  ifdef NEBDEBUG
   do i=1,nattgtrms
      iatom=rmsgp(i)
      iatm3 = 3*(iatom - 1)
      write(worldrank+50,'(a,f12.4,a,f12.4,a,f12.4)') &
         "NEBDBG: t1=",tangents(iatm3+1)," t2=",tangents(iatm3+2)," t3=",tangents(iatm3+3)
   enddo
#  endif

   call normalize(tangents, 3 * last_neb_atom)


   dotproduct = 0.d0
   dotproduct2 = 0.d0
  
!carlos: dotproduct scalar depends on atoms, so loop only over NEB region atoms

    do i=1,nattgtrms
      iatom=rmsgp(i)
      !write (6,*) "i,iatom ",i,iatom

      iatm3 = 3*(iatom - 1)
      !write (6,*) " f is ",f(iatm3+1),f(iatm3+2),f(iatm3+3)

      ! this is "real" force along tangent
      dotproduct = dotproduct + f(iatm3+1)*tangents(iatm3+1)
      dotproduct = dotproduct + f(iatm3+2)*tangents(iatm3+2)
      dotproduct = dotproduct + f(iatm3+3)*tangents(iatm3+3)

      ! now we do the spring forces
      ! note use of two different spring constants

      springforce(iatm3+1) = (xnext(iatm3+1)-x(iatm3+1))*spring2 - (x(iatm3+1)-xprev(iatm3+1))*spring
      springforce(iatm3+2) = (xnext(iatm3+2)-x(iatm3+2))*spring2 - (x(iatm3+2)-xprev(iatm3+2))*spring
      springforce(iatm3+3) = (xnext(iatm3+3)-x(iatm3+3))*spring2 - (x(iatm3+3)-xprev(iatm3+3))*spring

      dotproduct2 = dotproduct2 + springforce(iatm3+1)*tangents(iatm3+1)
      dotproduct2 = dotproduct2 + springforce(iatm3+2)*tangents(iatm3+2)
      dotproduct2 = dotproduct2 + springforce(iatm3+3)*tangents(iatm3+3)

      !write (6,*) " springf is ",springforce(iatm3+1),springforce(iatm3+2),springforce(iatm3+3)

      energy = energy + 0.5*spring2*(xnext(iatm3+1)-x(iatm3+1))*(xnext(iatm3+1)-x(iatm3+1))
      energy = energy + 0.5*spring2*(xnext(iatm3+2)-x(iatm3+2))*(xnext(iatm3+2)-x(iatm3+2))
      energy = energy + 0.5*spring2*(xnext(iatm3+3)-x(iatm3+3))*(xnext(iatm3+3)-x(iatm3+3))

      energy = energy + 0.5*spring*(x(iatm3+1)-xprev(iatm3+1))*(x(iatm3+1)-xprev(iatm3+1))
      energy = energy + 0.5*spring*(x(iatm3+2)-xprev(iatm3+2))*(x(iatm3+2)-xprev(iatm3+2))
      energy = energy + 0.5*spring*(x(iatm3+3)-xprev(iatm3+3))*(x(iatm3+3)-xprev(iatm3+3))
   end do

    !write (6,*) "dotproduct is ",dotproduct   
    !write (6,*) "dotproduct2 is ",dotproduct2 
    !write (6,*) "tangents are "
    !write (6,*) tangents(1:last_neb_atom*3)

    tempdot=0.d0

    do i=1,nattgtrms
      iatom=rmsgp(i)


      iatm3 = 3*(iatom - 1)
      jatm3 = 3*(i - 1)
      
      neb_force(jatm3+1) = 0.d0
      neb_force(jatm3+2) = 0.d0
      neb_force(jatm3+3) = 0.d0

      ! leave modified forces packed in rmsgp, less to communicate
      ! unpack them in force() when they are added to normal forces
      ! note that tangents still use normal coordinate index iatm3, while neb_force uses rmsgp pointer

      neb_force(jatm3+1) = neb_force(jatm3+1) - dotproduct*tangents(iatm3+1)
      neb_force(jatm3+2) = neb_force(jatm3+2) - dotproduct*tangents(iatm3+2)         
      neb_force(jatm3+3) = neb_force(jatm3+3) - dotproduct*tangents(iatm3+3)

      neb_force(jatm3+1) = neb_force(jatm3+1) + dotproduct2*tangents(iatm3+1)
      neb_force(jatm3+2) = neb_force(jatm3+2) + dotproduct2*tangents(iatm3+2)         
      neb_force(jatm3+3) = neb_force(jatm3+3) + dotproduct2*tangents(iatm3+3)

   end do
!
   
! carlos: what is nebrms for? just info, or is it used anywhere?
! if info, why do we calc dotproduct2 here? seems like it is a local variable.
! don't support it since it's not clear how to change for pNEB if 
! we don't know what it is.

   !dotproduct = 0.d0
   !dotproduct2= 0.d0
   !do iatm3 = 1,3*natom
      !dotproduct = dotproduct + f(iatm3)*f(iatm3)
      !dotproduct2= dotproduct2+ v(iatm3)*v(iatm3)
   !enddo

   !nebrms = dotproduct

#  ifdef NEBDEBUG
   write(worldrank+50,'(a)') "-------------------------------------------------"
#  endif

    return
end subroutine full_neb_forces

#endif  /* mpi */


subroutine pimd_neb_energy_report(filenumber)
   use neb_vars, only: neb_nrg_all,neb_nbead
   implicit none
   integer filenumber,ix
   _REAL_ sum

   sum = 0
   do ix = 1, neb_nbead
      sum = sum + neb_nrg_all(ix)
   enddo

   write(filenumber,'(a)') "NEB replicate breakdown:"
   do ix = 1, neb_nbead
      write(filenumber,'(a,i3,a,f13.4)') "Energy for replicate ",ix," = ",neb_nrg_all(ix)
   enddo
   write(filenumber,'(a,f13.4)') "Total Energy of replicates = ",sum

end subroutine pimd_neb_energy_report


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  QI correlation functions (LES-PIMD)                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine qi_corrf_les ( x, mass )

#ifdef MPI
#  ifdef LES

   use constants, only: hbar, kB
   use evb_pimd,  only: natomCL, vel0_bead, bead_dcrypt
   use miller,    only: i_qi, gradRC, gradRC_norm, div_ndx
   use evb_parm,  only: nbias
   use evb_data,  only: f_v, F_QI, G_QI

   implicit none

#  include "md.h"
#  include "les.h"
#  include "memory.h"
   include 'mpif.h'
#  include "parallel.h"

  _REAL_, intent(in) :: x(3,natom) !! coordinates
  _REAL_, intent(in) :: mass(natom)

   !  ..........................................................................

  integer :: idim, icopy, nbead
  integer :: m, n, mm, nn, mdx

  _REAL_ :: dx(3*natomCL), gnorm
  _REAL_ :: beta, coeff
  _REAL_ :: coeff_F, coeff_G1, coeff_G2
  _REAL_ :: F_sum1, F_sum2, F_sum3, F_sum4, G_sum, Ftot, Gtot
  _REAL_ :: ftmp, fv
  _REAL_ , intrinsic :: sqrt, dot_product

!  +----------------------------------------------------------+
!  |  QI dynamical factors                                    |
!  +----------------------------------------------------------+

   if( i_qi == 2 ) then

! write(6,*) 'hbar = ', hbar
! write(6,*) 'beta = ', beta
! write(6,*) 'kB = ', kB
! write(6,*) 'temp0 = ', temp0

      beta = 1.0d0 / ( kB * temp0 )
      coeff = dble( ncopy ) / ( hbar * beta )**2

! 03162009      coeff_F  = 2.0d0 / dble( ncopy )
      coeff_F  = 2.0d0
!     coeff_F  = 2.0d0 / dble( ncopy )

! 03132009      coeff_G1 = 6.0d0 * dble( ncopy ) / beta / beta
      coeff_G1 = 2.0d0 * dble( 3*natomCL) * dble( ncopy ) / beta / beta
      coeff_G2 = 4.0d0 * coeff / beta

!  +---------------------------------------------------------------------------+
!  |  f_v [ Eq. 2.24, JCP 120, 3086 (2004) ]                                   |
!  |                                                                           |
!  |  f_v(x_1,x_2,...,x_P) = m(iP/2hbarB)^2 grad_RC(x_P) . [x_1 - x_(P-1)]     |
!  |                         grad_RC(x_(P/2)) . [x_(P/2+1) - x_(P/2-1)]        |
!  +---------------------------------------------------------------------------+

      ftmp = 1.0d0

      do n = 1, nbias

         gnorm = 0.0d0
         dx(:) = 0.0d0

         do m = 1, natomCL

            mdx = ( m - 1 ) * 3

            if( div_ndx(n) == ncopy ) then

               mm = bead_dcrypt( m, 1 )
               nn = bead_dcrypt( m, ncopy-1 )

               dx(mdx+1) = x(1,mm) - x(1,nn)
               dx(mdx+2) = x(2,mm) - x(2,nn)
               dx(mdx+3) = x(3,mm) - x(3,nn)

            else if( div_ndx(n) == ncopy/2 ) then

               mm = bead_dcrypt( m, ncopy/2+1 )
               nn = bead_dcrypt( m, ncopy/2-1 )

               dx(mdx+1) = x(1,mm) - x(1,nn)
               dx(mdx+2) = x(2,mm) - x(2,nn)
               dx(mdx+3) = x(3,mm) - x(3,nn)

            else

               write(6,*) "QI error: definition of div_ndx is incorrect"
               write(6,*) "div_ndx(", n, ") = ", div_ndx(n)
               stop 6

            endif

!  +---------------------------------------------------------------------------+
!  |  Compute mass-weighted norm of gradient of RC                             |
!  +---------------------------------------------------------------------------+

            gnorm = gnorm + ( gradRC(mdx+1,n)**2 + gradRC(mdx+2,n)**2 &
                            + gradRC(mdx+3,n)**2 ) / mass(mm)

         enddo

         ! for "imaginary-time" velocity correlation function.
! KFW 02142009 gradRC computed in evb_umb is the negative of the derivative (for umbrella force),
!              so we negate it here to be consistent with the vector of dx
!
         gradRC_norm(n) = sqrt( gnorm )

!        ftmp = ftmp * dot_product( gradRC(:,n), dx(:) )
         ftmp = ftmp * dot_product( -gradRC(:,n), dx(:) ) / gradRC_norm(n)

!        write(6,*) '>>> dx(:) = ', dx(:)
!        write(6,*) '>>> gradRC = ', gradRC(:,n), n

      enddo


!  +---------------------------------------------------------------------------+
!  |  F and G factors [ Eq. 2.29, JCP 120, 3086 (2004) ]                       |
!  |                                                                           |
!  |  F(x_1,x_2,...,x_P) = -mP/hbar^2/B^2 {sum(1,P/2) - sum(P/2+1,P)}          |
!  |       (x_k - x_k-1)^2 + 2/P {sum(1,P/2) - sum(P/2+1,P)} V(x_k)            |
!  |                                                                           |
!  |  G(x_1,x_2,...,x_P) = 2dP/B^2 - 4mP/hbar^2/B^3 sum(1,P) [x_k - x_(k-1)]^2 |
!  +---------------------------------------------------------------------------+


! +++++++++++++++++++++++++++++++
! write(6,*) 
! write(6,*) ' natom = ', natom
! write(6,*) ' nbead decrypt' 
! write(6,*) ' bead index, atom index'

!     do n = 1, ncopy
!        do m = 1, natomCL
!           mm = bead_dcrypt( m, n )
!           write(6,*) n, mm, mass(mm)
!        enddo
!     enddo

      F_sum1 = 0.0d0

      do m = 1, natomCL
         mm = bead_dcrypt( m, 1 )
         nn = bead_dcrypt( m, ncopy )
         F_sum1 = F_sum1 + mass(mm) * ( ( x(1,mm) - x(1,nn) )**2 &
                                      + ( x(2,mm) - x(2,nn) )**2 &
                                      + ( x(3,mm) - x(3,nn) )**2 )
      enddo

      do nbead = 2, ncopy/2
         do m = 1, natomCL
            mm = bead_dcrypt( m, nbead )
            nn = bead_dcrypt( m, nbead-1 )
            F_sum1 = F_sum1 + mass(mm) * ( ( x(1,mm) - x(1,nn) )**2 &
                                         + ( x(2,mm) - x(2,nn) )**2 &
                                         + ( x(3,mm) - x(3,nn) )**2 )
         enddo
      enddo


      F_sum2 = 0.0d0

      do nbead = ncopy/2+1, ncopy
         do m = 1, natomCL
            mm = bead_dcrypt( m, nbead )
            nn = bead_dcrypt( m, nbead-1 )
            F_sum2 = F_sum2 + mass(mm) * ( ( x(1,mm) - x(1,nn) )**2 &
                                         + ( x(2,mm) - x(2,nn) )**2 &
                                         + ( x(3,mm) - x(3,nn) )**2 )
         enddo
      enddo

!  +---------------------------------------------------------------------------+
!  |  2/P { sum(s,P/2-1) - sum(P/2+1,P-1) } V(r^s)                             |
!  +---------------------------------------------------------------------------+

      F_sum3 = 0.0d0

      do n = 1, ncopy/2-1
        F_sum3 = F_sum3 + vel0_bead(n)
      enddo

      F_sum4 = 0.0d0
      do n = ncopy/2+1, ncopy-1
        F_sum4 = F_sum4 + vel0_bead(n)
      enddo


! ++++++++++++++++++++++++++++++++

      G_sum = 0.0d0

      do m = 1, natomCL
         mm = bead_dcrypt( m, 1 )
         nn = bead_dcrypt( m, ncopy )
         G_sum = G_sum + mass(mm) * ( ( x(1,mm) - x(1,nn) )**2 &
                                    + ( x(2,mm) - x(2,nn) )**2 &
                                    + ( x(3,mm) - x(3,nn) )**2 )
      enddo

      do nbead = 2, ncopy
         do m = 1, natomCL
            mm = bead_dcrypt( m, nbead )
            nn = bead_dcrypt( m, nbead-1 )
            G_sum = G_sum + mass(mm) * ( ( x(1,mm) - x(1,nn) )**2 &
                                       + ( x(2,mm) - x(2,nn) )**2 &
                                       + ( x(3,mm) - x(3,nn) )**2 )
         enddo
      enddo

      ! for "imaginary-time" velocity correlation function.
      fv = - 0.25d0 * dble(ncopy) * coeff * ftmp

      Ftot = - coeff * ( F_sum1 - F_sum2 ) + coeff_F * ( F_sum3 - F_sum4 )

! The factor of 2 is to convert to 1/2 m V^2, so as to get the correct units
!kfw  Ftot = - coeff * ( F_sum1 - F_sum2 ) * 2.0d0 + coeff_F * ( F_sum3 - F_sum4 )

! The factor of 1/2 is to convert to 1/2 m V^2, so as to get the correct units
!kfw  Gtot = coeff_G1 - coeff_G2 * G_sum * 2.0d0
      Gtot = coeff_G1 - coeff_G2 * G_sum
!     Gtot = coeff_G1 - coeff_G2 * G_sum * 0.50d0

      f_v  = fv
      F_QI = Ftot
      G_QI = Gtot

!     write(6,*) '<<< coeff, F_sum1, F_sum2 = ', coeff, F_sum1, F_sum2
!     write(6,*) '<<< coeff_F, F_sum3, F_sum4 = ', coeff_F, F_sum3, F_sum4
!     write(6,*) '<<< vel0_bead = ', vel0_bead(:)

!   write(6,*) '>>> ncopy = ', ncopy

!   write(6,*) '>>> f_v = ', fv, coeff, ftmp
!   write(6,*) '>>>   F = ', F_QI
!   write(6,*) '>>>   G = ', Gtot


   endif

#  endif
#endif

   end subroutine qi_corrf_les



