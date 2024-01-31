! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"

! This file should contain the routines related to DFTB/GB.
!
!------------------------------------------------
! Gustavo Seabra (2006)
! Quantum Theory Project, University of Florida
! E-Mail: seabra@qtp.ufl.edu
! Url:    http://www.qtp.ufl.edu/~seabra
!------------------------------------------------
subroutine qm2_dftb_gb_shift(qm_gb, qmmm_mpi, qmmm_scratch, qmmm_struct, scf_mchg)
    ! Calculates the contribution from GB potential to the Shift vector. and
    ! updates the shift vector.
    !
    ! NOTE: The potentials are calculated in Hartrees/Angstrom. Multiplying by
    !       BOHRS_TO_A puts them in Hartrees/Bohrs.
    !
    ! NOTE: The '-' sign is here because we used the mulliken charges in the calculation of
    !       the potential, inplace of the electron density difference. But, this
    !       electron density difference is exactly the opposite as the mulliken charges.

!! Modules
    use qmmm_module, only : qm_gb_structure, qmmm_mpi_structure, qmmm_scratch_structure
    use qm2_dftb_module, only : ks_struct
    use constants, only : BOHRS_TO_A, zero
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

!! Passed in
    type(qm_gb_structure), intent(inout) :: qm_gb
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qmmm_scratch_structure), intent(inout) :: qmmm_scratch
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink) ! Mulliken charges per atom

#ifdef MPI
    include 'mpif.h'
#endif

!! Locals
    integer :: i, ier

    !Step 1 - calculate the potential at QM atoms due to QM atoms with current Mulliken charges
    !Parallel
    qm_gb%gb_qmpot(1:qmmm_struct%nquant_nlink) = zero
!Parallel
    call qmgb_calc_qm_pot(qm_gb, qmmm_mpi, qm_gb%gb_qmpot, qm_gb%qmqm_onefij, scf_mchg)
#ifdef MPI
    !Since only the master thread does most of DFTB at present we need to reduce
    !the gb_qmpot array from it's current distributed form. Hopefully
    !this can go away once the DFTB is properly parallelized.
# ifdef USE_MPI_IN_PLACE
    if (qmmm_mpi%commqmmm_master) then
        call mpi_reduce(MPI_IN_PLACE, qm_gb%gb_qmpot, qmmm_struct%nquant_nlink, &
            MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
    else
        call mpi_reduce(qm_gb%gb_qmpot, 0, qmmm_struct%nquant_nlink, &
            MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
    end if
# else
    call mpi_reduce(qm_gb%gb_qmpot, qmmm_scratch%matsize_red_scratch, qmmm_struct%nquant_nlink, &
        MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
    if (qmmm_mpi%commqmmm_master) &
        qm_gb%gb_qmpot(1:qmmm_struct%nquant_nlink) = qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
# endif
#endif

    !Step 2 - Add the mm potential and the qm potential to the shift.
!Only master does this at the moment
    if (qmmm_mpi%commqmmm_master) then
        ! The potential comes in electrons/A, and needs to be converted to electrons/Bohrs.
        do i = 1, qmmm_struct%nquant_nlink
            ks_struct%shift(i) = ks_struct%shift(i) + (qm_gb%gb_mmpot(i) + qm_gb%gb_qmpot(i))*BOHRS_TO_A
        end do
    end if

    return
end subroutine qm2_dftb_gb_shift
