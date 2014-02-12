! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"
#include "copyright.h"

! DFTB GB implementation by Gustavo Seabra (UFL) and Ross Walker (TSRI), 2005

subroutine hamilshift(qm_coords,atomtype,uhubb,niter,gammamat,shift,scf_mchg)

   !=======================================================
   ! get the hubbard contribution to the H matrix elements
   !=======================================================

   use qmmm_module, only : qmmm_struct, qm2_struct
   use qm2_dftb_module, only : mol, mcharge

   implicit none

   ! Passed in:
   _REAL_ , intent(in)  :: qm_coords(3,qmmm_struct%nquant_nlink)          ! atomic coordinates
   integer, intent(in)  :: atomtype(qmmm_struct%nquant_nlink)     ! 
   _REAL_ , intent(in)  :: uhubb(*)          ! list of hubbard parameters
   integer, intent(in)  :: niter             ! step in scf cycle if == 1: 
                                             !              build up gammamat
   _REAL_ , intent(out) :: gammamat(qmmm_struct%nquant_nlink,qmmm_struct%nquant_nlink)
   _REAL_ , intent(out) :: shift(qmmm_struct%nquant_nlink)        ! array contains shifts for hamilton matrix elements
   _REAL_ , intent(out) :: scf_mchg(qmmm_struct%nquant_nlink)

   ! Locals
   integer :: i,j
   _REAL_ :: tmpvalue

  
   scf_mchg(1:qmmm_struct%nquant_nlink) = mcharge%qzero(atomtype(1:qmmm_struct%nquant_nlink)) &
                                        - mol%qmat(1:qmmm_struct%nquant_nlink)
   ! build up gammamatrix before the first iteration. This matrix is
   ! then reused at the following iterations
   if (niter == 1) then
      call gammamatrix(qmmm_struct%nquant_nlink,qm_coords,atomtype,uhubb,gammamat)
   endif

   ! Calculate atomic hamilton shift (=sum over gamma * charges)
   do i=1,qmmm_struct%nquant_nlink
      do j=1,qmmm_struct%nquant_nlink
         ! gammamat is lower diagonal. All elements where j > i are zero.
         if (j > i) then
            tmpvalue = gammamat(j,i)
         else
            tmpvalue = gammamat(i,j)
         endif
         shift(i) = shift(i) - scf_mchg(j)*tmpvalue
      end do
   end do

   return

end subroutine HAMILSHIFT

