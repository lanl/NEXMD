! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"

subroutine gammamatrix(nquant_nlink, qm_coords, atomtype, u, gammamat)
!===========================================================================
! Build lower triangular Gamma matrix containing short range terms
!
!  INPUT Parameter:
!   INTEGER nquant_nlink          number of atoms
!   REAL*8 qm_coords(3,*)      position of atoms
!   REAL*8 u(*)          hubbard parameters
!
!  OUTPUT:
!   REAL*8 gammamat(*,*) matrix containing the values of the ewlad potential
!                       in the upper triangular part
!   !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
!   !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!
!
!============================================================================

    use constants, only : A_TO_BOHRS

    implicit none

!Passed in:
    integer, intent(in)  :: nquant_nlink                ! number of atoms
    _REAL_, intent(in)  :: qm_coords(3, nquant_nlink)           ! position of atoms
    integer, intent(in)  :: atomtype(*)        ! atom types
    _REAL_, intent(in)  :: u(*)               ! hubbard parameters
    _REAL_, intent(out) :: gammamat(nquant_nlink, nquant_nlink)  ! matrix containing the values of the gamma potential

!Locals
    integer :: i, j
    _REAL_  :: r(3)
    _REAL_  :: gval, norm
    external GAM12

    do i = 1, nquant_nlink
        do j = 1, i
            r(1:3) = (qm_coords(1:3, i) - qm_coords(1:3, j))*A_TO_BOHRS

            gval = 0.0d0

            ! Distance between the 2 atoms
            norm = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
            ! get value for Gamma
            call GAM12(norm, u(atomtype(i)), u(atomtype(j)), gval)

            gammamat(i, j) = gval

        end do
    end do
end subroutine gammamatrix
