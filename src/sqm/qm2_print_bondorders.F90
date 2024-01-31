! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

! Calculates and prints the Bond Orders.
! Written by Johan Strumpfer (2010)

subroutine qm2_print_bondorders(qm2_params, qm2_struct, qmmm_struct)

    ! Requires a converged density matrix stored in qm2_struct%den_matrix
    use qmmm_module, only : qm2_structure
    use ElementOrbitalIndex, only : elementSymbol
    use qmmm_struct_module, only : qmmm_struct_type
    use qm2_params_module, only : qm2_params_type

    implicit none
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qm2_structure), intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    integer :: loop_count, orb_beg_i, orb_end_i, tri_i
    integer :: orb_beg_j, orb_end_j, tri_j, tri
    integer :: i, j, k, tri_k1, tri_k2, iqm, jqm
    _REAL_, dimension(:), pointer :: den_matrix2
    _REAL_ :: BO

    ! CALCULATE AND PRINT THE BOND ORDERS
    write (6, *) ''
    write (6, '("  QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER")')
    do iqm = 1, qmmm_struct%nquant_nlink
        orb_beg_i = qm2_params%orb_loc(1, iqm)
        orb_end_i = qm2_params%orb_loc(2, iqm)
        do jqm = 1, iqm - 1
            orb_beg_j = qm2_params%orb_loc(1, jqm)
            orb_end_j = qm2_params%orb_loc(2, jqm)
            BO = 0
            do i = orb_beg_i, orb_end_i
                do j = orb_beg_j, orb_end_j
                    tri = qm2_params%pascal_tri1(i) + j
                    BO = BO + qm2_struct%den_matrix(tri)*qm2_struct%den_matrix(tri)
                end do
            end do
            write (6, '("  QMMM:    ",I4," ",A4,"  ",I4," ",A4," ",f16.8)') &
                iqm, elementSymbol(qmmm_struct%iqm_atomic_numbers(iqm)), &
                jqm, elementSymbol(qmmm_struct%iqm_atomic_numbers(jqm)), BO
        end do
    end do

end subroutine qm2_print_bondorders
