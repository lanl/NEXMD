! <compile=optimized>
#include "dprec.fh"

!==================================================
! Calculates the energy shift due to the
! presence of the external charges.
!
! Returns: shiftE(j) = SUM_k {Q_k(classical)/r_jk}
!
! Where:
!       j = QM atom
!       k = MM atom
!
! include only REAL QM-MM interactions - skip link atoms.
!==================================================
subroutine externalshift(qmmm_mpi, qmmm_struct, qm_coords, izp, shiftE)

    use qm2_dftb_module, only : mcharge
    use qmmm_module, only : qmmm_mpi_structure
    use constants, only : A_TO_BOHRS
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

    ! Passed in:
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    integer, intent(in)  :: izp(qmmm_struct%nquant_nlink)     ! IZP for each atom
    _REAL_, intent(in)  :: qm_coords(3, qmmm_struct%nquant_nlink)     ! Coordinates
    _REAL_, intent(out) :: shiftE(qmmm_struct%nquant_nlink)  ! Energy shift for each atom

    !Locals
    integer :: i, j, k
    _REAL_  :: dif(3)
    _REAL_  :: r
    _REAL_  :: r2
    _REAL_  :: gamma

    do j = 1, qmmm_struct%nquant_nlink ! loop through real qm atoms
        shiftE(j) = 0.0d0
        do k = 1, qmmm_struct%qm_mm_pairs ! loop through external charges
            dif(1:3) = (qm_coords(1:3, j) - qmmm_struct%qm_xcrd(1:3, k))*A_TO_BOHRS
            r2 = dif(1)*dif(1) + dif(2)*dif(2) + dif(3)*dif(3)

            ! Due to the different charge magnitude in quantum and classical parts,
            ! a scaling may be necessary in this interaction. See eq. (4) in
            ! Elstner et al, J. Mol. Struc. (Theochem), 632, 24--41 (2003).

            ! WHEN CHANGING HERE, REMEMBER TO MAKE THE SAME CHANGES IN
            ! qm2_dftb_externalchgrad.

            ! Hubbard Parameter (for the scaling factor):
            ! uhub=mcharge%uhubb(izp(j))

            ! \zeta = 4.0 in Eq. (4):
            ! gamma = 1.0/sqrt(r2 + (0.5/uhub + 0.5/uhub)**2)

            ! \zeta = 0.4:
            ! gamma = 1.0/sqrt(r2 + 0.1 * (1.0/uhub)**2)

            ! In DFTB, \zeta in Eq. (4) is set to zero.
            ! This is similar to no scaling at all. (p.33)
            gamma = 1.0d0/sqrt(r2)
            shiftE(j) = shiftE(j) + gamma*qmmm_struct%qm_xcrd(4, k)
        end do
    end do
    return
end subroutine externalshift
