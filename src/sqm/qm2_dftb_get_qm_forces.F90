! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-

#include "copyright.h"
#include "dprec.fh"
#include "def_time.h"
subroutine qm2_dftb_get_qm_forces(dxyzqm)

!     Gets the forces from the DFTB calculation
!
!     The gradient has already been calculated in the first call to 
! qm2_dftb_energy, and now we just have to convert it into forces, and
! store in the proper location.
!
!     Note that the DFTB routines are still in F77. So, this routine here
! actually just sets some values to call the DFTB routines. With time, I'll
! translate the DFBT code to F90.
!
!     Variables for qm-mm:
!
!     coord(3,qmmm_struct%nquant_nlink) - Cartesian coordinates of qm atoms.
!     qmmm_struct%nquant_nlink    - Total number of qm atoms. (Real + link)

      use qm2_dftb_module, only : izp_str, mcharge, mol, lmax
      use qmmm_module, only : qmmm_struct, qmmm_nml, qmmm_mpi
      use constants, only: AU_TO_KCAL, A_TO_BOHRS

      implicit none

      _REAL_, intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)

!Local
      integer i,j
      if (qmmm_mpi%commqmmm_master) then
        ! DISPERSION
        if (qmmm_nml%dftb_disper == 1) then
           call timer_start(TIME_QMMMDFTBDISPF)
           call  dispersion_grad(qmmm_struct%nquant_nlink,dxyzqm)
           call timer_stop(TIME_QMMMDFTBDISPF)
        endif

        call timer_start(TIME_QMMMDFTBREPULF)
        ! Calculate gradient of the repulsive energy
        call dftb_repulsivegrd(qmmm_struct%nquant_nlink,izp_str%izp,qmmm_struct%qm_coords,dxyzqm)
        call timer_stop(TIME_QMMMDFTBREPULF)

        call timer_start(TIME_QMMMDFTBHZEROF)
        ! Gradient due to H zero. (Non charge dependent part)
        call dftb_hzero_grad(qmmm_struct%nquant_nlink, izp_str%izp,lmax,qmmm_struct%qm_coords,dxyzqm)
        call timer_stop(TIME_QMMMDFTBHZEROF)

        call timer_start(TIME_QMMMDFTBGAMMAF)
        ! here are the contributions due to gamma if in scf mode - Charge Dependent Part
        call dftb_gammagrad(qmmm_struct%nquant_nlink,qmmm_struct%qm_coords, izp_str%izp,mcharge%uhubb,&
                            dxyzqm)
        call timer_stop(TIME_QMMMDFTBGAMMAF)

!===========================================
!     Puts forces into the proper matix
!===========================================
        do i=1,qmmm_struct%nquant_nlink
           do j = 1,3
              ! Convert to (kcal/mol)/Angstroms before 
              ! storing in Amber force array. 
              dxyzqm(j,i) = dxyzqm(j,i) * AU_TO_KCAL * A_TO_BOHRS
           end do
        end do

      end if  !(qmmm_mpi%commqmmm_master)

     ! part of qm_ewald.f. So, there's no need to calculate this here.
                     
      RETURN
end subroutine qm2_dftb_get_qm_forces

!============================================================================
subroutine dftb_gammagrad(nquant_nlink,qm_coords,atomtype,uhubb,gmgrd)

   use qm2_dftb_module, only: NNDIM, ks_struct
   use qmmm_module, only : qm2_struct

   implicit none

!! Passed in:
   integer, intent(in ) :: nquant_nlink          ! number of atoms in cell
   integer, intent(in ) :: atomtype(*)  ! list of atomic types
   _REAL_ , intent(in ) :: qm_coords(3,nquant_nlink)     ! atomic coordinates
   _REAL_ , intent(out) :: gmgrd(3,nquant_nlink)   ! gamma contribution to the gradient
   _REAL_ , intent(in ) :: uhubb(*)     ! hubbard parameters

!! Locals
   integer :: i,k,l
   _REAL_  :: deriv(3)
   _REAL_  :: tmpderiv(3)

   ! construct gamma derivatives
   call dftb_gammamatrix_deriv(nquant_nlink,NNDIM,qm_coords,atomtype,uhubb, &
                               ks_struct%derivx,ks_struct%derivy,ks_struct%derivz)

   do k=1,nquant_nlink

      deriv(1:3) = 0.0d0

      do i=1,nquant_nlink
         if (i > k) then
            tmpderiv(1) = ks_struct%derivx(i,k)
            tmpderiv(2) = ks_struct%derivy(i,k)
            tmpderiv(3) = ks_struct%derivz(i,k)
         else if (i == k) then
            tmpderiv(1) = 0.0d0
            tmpderiv(2) = 0.0d0
            tmpderiv(3) = 0.0d0
         else
            ! shortrange1(rrmuind - ri)  = -shortrange1(ri - rrmuind)
            tmpderiv(1) = -ks_struct%derivx(k,i)
            tmpderiv(2) = -ks_struct%derivy(k,i)
            tmpderiv(3) = -ks_struct%derivz(k,i)
         endif

         deriv(1) = deriv(1) + qm2_struct%scf_mchg(i)*tmpderiv(1)
         deriv(2) = deriv(2) + qm2_struct%scf_mchg(i)*tmpderiv(2)
         deriv(3) = deriv(3) + qm2_struct%scf_mchg(i)*tmpderiv(3)
      end do

      do l = 1,3
         gmgrd(l,k) = gmgrd(l,k) + qm2_struct%scf_mchg(k)*deriv(l)
      end do

   end do

   return

end subroutine dftb_gammagrad

!=============================================================================
! Build lower triang. matrices containg derivatives of long+shortrange expr.
!
! !!! NOTE THAT shortrange1(ri - rj) = -shortrange1(rj - ri) !!!
! !!! NOTE THAT phi1(ri - rj) = -phi1(rj - ri) !!!
!=============================================================================
subroutine dftb_gammamatrix_deriv(nquant_nlink,DIM,qm_coords,atomtype,u,gammamat1x,gammamat1y,gammamat1z)

   use constants, only : A_TO_BOHRS
   implicit none

!! Passed in:
   integer, intent(in ) :: nquant_nlink          ! Number of atoms
   integer, intent(in ) :: DIM          ! Dimension of the (square) matrices gammamat1[x,y,z]
   integer, intent(in ) :: atomtype(*)  ! 
   _REAL_ , intent(in ) :: qm_coords(3,nquant_nlink)     ! position of atoms
   _REAL_ , intent(in ) :: u(*)         ! hubbard parameters
   _REAL_ , intent(out) :: gammamat1x(DIM,DIM)
   _REAL_ , intent(out) :: gammamat1y(DIM,DIM)
   _REAL_ , intent(out) :: gammamat1z(DIM,DIM)

!!Locals
   integer :: i,j
   _REAL_ :: r(3),gdrv
   _REAL_ :: norm

   _REAL_  :: short_deriv(3), tol, basis(3,3)

   do i=1,nquant_nlink
      do j=1,(i-1)
         r(1:3)=(qm_coords(1:3,j)-qm_coords(1:3,i))*A_TO_BOHRS

         ! Lower matrix gammamat1x  contains (phi1+shortrange1)(ri-rj)x
         ! Lower matrix gammamat1y  contains (phi1+shortrange1)(ri-rj)y
         ! Lower matrix gammamat1z  contains (phi1+shortrange1)(ri-rj)z
         ! !!! NOTE THAT gamma1(ri - rj) = -gamma1(rj - ri) !!!i

         norm = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
         call GAM121(norm,u(atomtype(j)),u(atomtype(i)),gdrv)
         gammamat1x(i,j)= gdrv * r(1)
         gammamat1y(i,j)= gdrv * r(2)
         gammamat1z(i,j)= gdrv * r(3)

      end do
   end do

end subroutine dftb_gammamatrix_deriv

!
! USUALGRD
! ========
!
! Does the gradient of the Hamiltonian part of the SCC-DFTB energy.
!
subroutine dftb_hzero_grad(nquant_nlink,izp,lmax,qm_coords,grad)

   use constants, only : BOHRS_TO_A
   use qm2_dftb_module, only: NDIM,LDIM, ks_struct
   use qmmm_module, only : qmmm_nml
   implicit none

! Passed in:
   integer, intent(in ) :: nquant_nlink
   integer, intent(in ) :: izp(*)
   integer, intent(in ) :: lmax(*)
   _REAL_ , intent(inout) :: qm_coords(3,nquant_nlink)
   _REAL_ , intent(out) :: grad(3,nquant_nlink)

! Pointers
   ! to ks_struct
   integer, pointer :: ind(:)   ! (*)
   _REAL_ , pointer :: ev(:)    ! (*)
   _REAL_ , pointer :: occ(:)   ! (*)
   _REAL_ , pointer :: shift(:) ! (*)
   _REAL_ , pointer :: a(:,:)   ! (mdim,mdim)
   _REAL_ , pointer :: b(:,:)   ! (mdim,mdim)
   _REAL_ , pointer :: au(:,:)  ! (ldim,ldim)
   _REAL_ , pointer :: bu(:,:)  ! (ldim,ldim)
   _REAL_ , pointer :: auh(:,:) ! (ldim,ldim)
   _REAL_ , pointer :: buh(:,:) ! (ldim,ldim)

!!Locals:
   integer :: m,n,i,j,k,lj,lk,mj,mk
   integer :: mu,nu,izpj,izpk,indj,indk
   _REAL_  :: ocmcc,dgrh,xhelp,dgrs,dgr,dtmp
!!new locals:
   _REAL_  :: p(ndim,ndim),ep(ndim,ndim)
   integer :: iend(nquant_nlink)

!! Step size:
!! ----------
   ! Step size for "semi-numerical" derivative.
   ! The original DFTB code uses 1.0e-2, which is too big
   ! for accurate derivatives. On the other hand, the 
   ! distances in the integral files are in steps of 0.02, so 
   ! using deltax much smaller than that may not be completely
   ! meaningful. But since it gives MUCH, MUCH better results, I'll use
   ! it anyways.

   _REAL_, parameter :: deltax = 1.0d-5
   _REAL_, parameter :: rcdx = 1.0d0/deltax

   ! Since the H and S derivatives are done by:
   !
   !           f(x+delta) - f(x-delta)
   !   f'(x) = -----------------------
   !                  2 * delta
   !
   ! I believe that this should be 1/(2*delta).
   ! BUT, that doesn't work. I must be missing something here.
   !
!! Pointers
   ! to ks_struct
   ind   => ks_struct%ind  
   ev    => ks_struct%ev   
   occ   => ks_struct%occ  
   shift => ks_struct%shift
   a     => ks_struct%a    
   b     => ks_struct%b    
   au    => ks_struct%au   
   bu    => ks_struct%bu   
   auh   => ks_struct%auh  
   buh   => ks_struct%buh  

   p = 0.0d0
   ep= 0.0d0
   do m = 1,ndim
      do n = 1,m-1
         do i = 1,ndim
            ocmcc  = occ(i)*a(m,i)*a(n,i)
            p(m,n) = p(m,n) + ocmcc
            ep(m,n) = ep(m,n) + ocmcc*ev(i)
         end do
         p(n,m)=p(m,n)
         ep(n,m)=ep(m,n)
      end do
   end do

   do j = 1,nquant_nlink-1 
      iend(j)=ind(j+1)
   End Do
   iend(nquant_nlink)=ind(nquant_nlink)+ (lmax(izp(nquant_nlink))-1)**2 + 2*lmax(izp(nquant_nlink))-1
      
   do j = 1,nquant_nlink     ! Loop through all atoms 
      do k = 1,nquant_nlink  ! Loop through all atoms
         if(k /= j)then
            dtmp=0.5d0*(shift(k)+shift(j))
            ! loop X,Y,Z to get the derivative
            do i = 1,3
               ! saves the position of the atom
               xhelp = qm_coords(i,j)

               ! Plus step
               qm_coords(i,j) = xhelp + deltax*BOHRS_TO_A
               call slkmatrices_a_to_bohrs(k,j,qm_coords,au,bu,LDIM)

               ! Less step
               qm_coords(i,j) = xhelp - deltax*BOHRS_TO_A
               call slkmatrices_a_to_bohrs(k,j,qm_coords,auh,buh,LDIM)

               ! Restore coordinate
               qm_coords(i,j) = xhelp

                do nu=ind(j)+1,iend(j)
                    n=nu-ind(j)
                  do mu=ind(k)+1,iend(k)
                    m=mu-ind(k)
                           dgrs =  (bu(m,n)-buh(m,n))*rcdx
                           dgrh =  (au(m,n)-auh(m,n))*rcdx + dgrs*dtmp
                           grad(i,j) = grad(i,j) + dgrh*p(mu,nu) - dgrs*ep(mu,nu)
                     end do
                  end do
               end do
            endif
         end do
      end do


   !  End If
end subroutine dftb_hzero_grad

