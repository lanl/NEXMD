! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_fock_predict(qmmm_nml, num_qmmm_calls, hmatrix, matsize, fock_matrix, &
    fock_mat_final1, fock_mat_final2, &
    fock_mat_final3, fock_mat_final4)

!Subroutine by Ross Walker and Gustavo Seabra - This routine
!attempts to predict the Fock matrix based on Pulay et al CPL, 2004, 386, 272-278.

    use qmmm_nml_module, only : qmmm_nml_type

    implicit none

    !Passed in
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    integer, intent(in) :: num_qmmm_calls
    integer, intent(in) :: matsize
    _REAL_, intent(in) :: hmatrix(matsize)
    _REAL_, intent(out) :: fock_matrix(matsize)
    _REAL_, intent(in) :: fock_mat_final1(matsize)
    _REAL_, intent(in) :: fock_mat_final2(matsize)
    _REAL_, intent(in) :: fock_mat_final3(matsize)
    _REAL_, intent(in) :: fock_mat_final4(matsize)

    !This algorithm works by taking the previous final fock matrices from the last
    !4 MD steps and using these to predict the initial fock matrix for the next
    !MD step.

    !in parallel these arrays are all partial matrices but this is fine since
    !the fock matrix immediately gets reduced before actually being used.

    !If this is MD step 1 to 4 then there is nothing to do, we are using the standard
    !procedure for the SCF.
    if (num_qmmm_calls < 5) return

    !if we are step 5 or over we can attempt to build a Fock matrix.
    fock_matrix(1:matsize) = hmatrix(1:matsize)

    fock_matrix(1:matsize) = fock_matrix(1:matsize) + qmmm_nml%fockp_d4*fock_mat_final4(1:matsize) &
        + qmmm_nml%fockp_d3*fock_mat_final3(1:matsize) + qmmm_nml%fockp_d2*fock_mat_final2(1:matsize) &
        + qmmm_nml%fockp_d1*fock_mat_final1(1:matsize)
    return

end subroutine qm2_fock_predict

subroutine qm2_fock_store(qm2_struct, matsize, fock_matrix, hmatrix)

    use qmmm_module, only : qm2_structure

    implicit none

    !This routine takes a matrix from the end of the SCF and stores it in the latest
    !storage array for it and moves the other stored arrays down the chain.

    !Passed in
    type(qm2_structure), intent(inout) :: qm2_struct
    integer, intent(in) :: matsize
    _REAL_, intent(in) :: fock_matrix(matsize)
    _REAL_, intent(in) :: hmatrix(matsize)

    !Local
    _REAL_, dimension(:), pointer :: pointer_temp

!In parallel here each thread stores a it's part of the fock matrix - this ultimately
!gets reduced before diagonalization so does not cause problems. Although it obviously
!represents a serial bottleneck.

    !Step 1 move all current stored matrices down the chain.

!RCW: Just rotate pointers for speed.
    pointer_temp => qm2_struct%fock_mat_final4
    qm2_struct%fock_mat_final4 => qm2_struct%fock_mat_final3
    qm2_struct%fock_mat_final3 => qm2_struct%fock_mat_final2
    qm2_struct%fock_mat_final2 => qm2_struct%fock_mat_final1

    qm2_struct%fock_mat_final1 => pointer_temp

    !Step 2 store the current fock matrix (subtracting off the hmatrix) in fock_mat_final1
    qm2_struct%fock_mat_final1(1:matsize) = fock_matrix(1:matsize) - hmatrix(1:matsize)

    return

end subroutine qm2_fock_store
