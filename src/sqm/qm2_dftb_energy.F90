! <compile=optimized>

#include "copyright.h"
#include "dprec.fh"
subroutine qm2_dftb_energy(qm2_params, qmmm_nml, qmmm_scratch, qmmm_mpi, qm_gb, qmewald, qm2_struct, qmmm_struct, escf, scf_mchg)

! Calculates the self-consistent-charge (SCC) DFTB energy.
! The energy is returned in 'escf'.
!
!     Note that the DFTB routines are still in F77. So, this routine here
! actually just sets some values to call the DFTB routines. With time, I'll
! translate the DFTB code to F90.
!
!     Variables for qm-mm:
!
!     qmmm_struct%nquant_nlink    - Total number of qm atoms. (Real + link)
!
! This routine is called from qm2_energies.f

!In parallel all threads enter here.

    use qmmm_module, only : qm2_structure, qmmm_mpi_structure, qm_gb_structure, qm_ewald_structure, qmmm_scratch_structure
    use ElementOrbitalIndex, only : elementSymbol
    use qm2_dftb_module, only : disper, mol, lmax, izp_str, mcharge
    use constants, only : AU_TO_EV, AU_TO_KCAL, A_TO_BOHRS, BOHRS_TO_A
    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_nml_module, only : qmmm_nml_type
    use qm2_params_module, only : qm2_params_type

    implicit none
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qm2_structure), intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qm_gb_structure), intent(inout) :: qm_gb
    type(qm_ewald_structure), intent(inout) :: qmewald
    type(qmmm_scratch_structure), intent(inout) :: qmmm_scratch

    !Passed in
    _REAL_, intent(out)   :: escf
    _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

    !Locals
    !     Calculation results
    !
    _REAL_ total_e, geseatom

    integer k, i, j
    integer mm_link_atom

    !=======================
    !     Initializaton
    !=======================

    if (qmmm_mpi%commqmmm_master) then
        ! Calculation results in a.u.

        !==============================
        !     Set External charges
        !==============================

        if (qmmm_nml%verbosity > 0) then
            if (qmmm_struct%qm_mm_pairs > 0) then
                write (6, *) "QMMM SCC-DFTB: Classical atoms from nonbond list taken as point charges."
            else
                write (6, *) "QMMM SCC-DFTB: No external charges defined."
            end if
        end if

        !=======================
        !     Print stuff
        !=======================

        if (qmmm_nml%verbosity > 3) then
            write (6, '(/," QMMM SCC-DFTB: QM Region Input Cartesian Coordinates ")')
            write (6, '(" QMMM SCC-DFTB: ",4X,"NO.",2X,"TYP",2X,"L",2x"AT#",2X,"SYM",10X,"X",16X,"Y",16X,"Z",13X,"Charge")')
            do i = 1, qmmm_struct%nquant_nlink
                write (6, '(" QMMM SCC-DFTB: ",I6,3X,I2,2X,i2,1X,I3,4X,A2,2X,4F16.10)') i, &
                    izp_str%izp(i), &
                    lmax(izp_str%izp(i)), &
                    qmmm_struct%iqm_atomic_numbers(i), &
                    elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
                    (qmmm_struct%qm_coords(j, i), j=1, 3), scf_mchg(i) !( mcharge%qzero( izp_str%izp(i) ) - mol%qmat(i) )
            end do

            if (qmmm_struct%qm_mm_pairs > 0 .and. qmmm_nml%verbosity > 4) then
                write (6, *) 'QMMM SCC-DFTB: number of external charges', qmmm_struct%qm_mm_pairs
                write (6, *) 'QMMM SCC-DFTB: Coordinates of external charges (XYZ)'
                write (6, *)
                write (6, '(" QMMM SCC-DFTB: ",i3)') qmmm_struct%qm_mm_pairs
                do i = 1, qmmm_struct%qm_mm_pairs
                    write (6, '(" QMMM SCC-DFTB: ",4(2x,f10.6))') (qmmm_struct%qm_xcrd(j, i), j=1, 4)
                end do
            end if
        end if

        ! If (dispers), read dispersion parameters.
        !   this HAS to be done at every structure, because the dispersion
        !   parameters do depend on the arrangement of the neighbours
        if (qmmm_nml%dftb_disper == 1) then
            call dispersion_params(qmmm_struct, qmmm_struct%nquant_nlink, izp_str%izp)
        end if
    end if !qmmm_mpi%commqmmm_master

!In parallel all threads enter here.

    qmmm_struct%elec_eng = 0.0d0
    qmmm_struct%enuclr_qmqm = 0.0d0
    disper%edis = 0.0d0 ! Dispersion Energy, in a.u.
    escf = 0.0d0

    !do the SCF
    call qm2_dftb_scf(qm2_params, qmmm_nml, qmmm_scratch, qmmm_mpi, qm_gb, qmewald, qm2_struct, qmmm_struct, escf, &
        qmmm_struct%elec_eng, qmmm_struct%enuclr_qmqm, scf_mchg)

    return
end subroutine qm2_dftb_energy
