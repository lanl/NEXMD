! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"

! This file should contain the routines related to DFTB/Ewald.
!
!------------------------------------------------
! Gustavo Seabra (2006)
! Quantum Theory Project, University of Florida
! E-Mail: seabra@qtp.ufl.edu
! Url:    http://www.qtp.ufl.edu/~seabra
!------------------------------------------------

subroutine qm2_dftb_ewald_shift(scf_mchg)
    ! Calculates the contribution from Ewald potential to the Shift vector. and
    ! updates the shift vector.
    !
    ! NOTE: we skip this if we are keeping the image charges fixed during the SCF (qm_ewald==2)
    !       but only after the first MD step has been done.
    !
    ! NOTE: The potentials are calculated in Hartrees/Angstrom. Multiplying by
    !       BOHRS_TO_A puts them in Hartrees/Bohrs.
    !
    ! NOTE: The '-' sign is here because we used the mulliken charges in the calculation of
    !       the potential, inplace of the electron density difference. But, this
    !       electron density difference is exactly the opposite as the mulliken charges.

!! Modules
    use qmmm_module, only : qmmm_struct, qmmm_nml, qmewald, qmmm_mpi, qmmm_scratch
    use qm2_dftb_module, only : ks_struct
    use constants, only : BOHRS_TO_A

    implicit none

!! Passed in
    _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink) ! Mulliken charges per atom

#ifdef MPI
    include 'mpif.h'
#endif

!! Locals
    integer :: i, ier

!!--
    ! Calculates the Ewald potential at the position of the atoms
    if (qmmm_nml%qm_ewald == 1 .OR. qmewald%ewald_startup) then
!Parallel
        call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg, & !qmewald%qmpot, &
            qmmm_struct%qm_coords, qmewald%kvec)

#ifdef MPI
        !Since only the master thread does most of DFTB at present we need to reduce
        !the qmpot array from it's current distributed form. Hopefully
        !this can go away once the DFTB is properly parallelized.
#ifdef USE_MPI_IN_PLACE
        if (qmmm_mpi%commqmmm_master) then
            call mpi_reduce(MPI_IN_PLACE, qmewald%qmpot, qmmm_struct%nquant_nlink, &
                MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
        else
            call mpi_reduce(qmewald%qmpot, 0, qmmm_struct%nquant_nlink, &
                MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
        end if
# else
        call mpi_reduce(qmewald%qmpot, qmmm_scratch%matsize_red_scratch, qmmm_struct%nquant_nlink, &
            MPI_DOUBLE_PRECISION, mpi_sum, 0, qmmm_mpi%commqmmm, ier)
        if (qmmm_mpi%commqmmm_master) &
            qmewald%qmpot(1:qmmm_struct%nquant_nlink) = qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
# endif
#endif

    end if

!Only master does this at the moment
    if (qmmm_mpi%commqmmm_master) then
        ! Puts the Ewald potential into Shift.
        ! The potential comes in au/A, and needs to be converted to au/Bohrs.
        do i = 1, qmmm_struct%nquant_nlink
            ks_struct%shift(i) = ks_struct%shift(i) - (qmewald%mmpot(i) + qmewald%qmpot(i))*BOHRS_TO_A
        end do
    end if

    !call timer_stop(TIME_QMMMENERGYSCFFOCKEWALD)
    return
end subroutine qm2_dftb_ewald_shift

!====================================================

subroutine qm2_dftb_ewald_corr(nquant_nlink, corr, mmpot, scf_mchg) !qmat)

    !This routine should be called after a call to qm2_helect. Up to this point
    !the energy for the Ewald sum only included half of the term from the MM atoms
    !and half from the QM atoms. The QM atoms should contribute only half but the
    !MM atoms should contribute full. This routine corrects EE for this.

    use constants, only : AU_TO_EV, BOHRS_TO_A, zero, half
    use qmmm_module, only : qm2_params

    implicit none

    !Passed in
    integer, intent(in)  :: nquant_nlink
    _REAL_, intent(out) :: corr   ! Correction in Hartree / Bohrs (a.u.)
    _REAL_, intent(in)  :: mmpot(nquant_nlink)
    _REAL_, intent(in)  :: scf_mchg(nquant_nlink)

    !Local
    _REAL_ :: etemp
    integer :: i, i1, ia, ib, i2

    corr = zero

    etemp = zero
    do i = 1, nquant_nlink
        etemp = etemp + mmpot(i)*scf_mchg(i) !qmat(i)
    end do

    corr = half*etemp*BOHRS_TO_A

    return

end subroutine qm2_dftb_ewald_corr
