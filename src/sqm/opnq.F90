#include "copyright.h"
#include "dprec.fh"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the implementation of
! charge-dependent NB interactions (OPNQ) see
! T. Giese and D. York, J. Chem. Phys. 2007, v127, p194101
!
! Implemented by Taisung Lee (Rutgers, 2011)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opnq

    use EVDWMOD

    implicit none

! visuability declaration
    public:: Opnq_fock, Opnq_fock_atom_pair, Opnq_LJ_atom_pair

    private:: LJ2OPNQ, Initialize
    private:: MaxAtomicNumber_MM_opnq

    integer, parameter::MaxAtomicNumber_MM_opnq = 20

contains

    subroutine Opnq_fock(qmmm_opnq, qm2_params, qm2_struct, qmmm_struct, fock, density)
!***********************************************************************
!
!  This subroutine calculates the OPNQ contributions to the Fock matrix
!
!  Taisung Lee, Rutgers, 2011
!
!***********************************************************************
        use constants, only : A_TO_BOHRS, zero
        use ElementOrbitalIndex, only : MaxValenceOrbitals
        use qmmm_module, only : qm2_structure, qmmm_opnq_structure
        use qm2_params_module, only : qm2_params_type
        use qmmm_struct_module, only : qmmm_struct_type
        implicit none
        type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq
        type(qm2_params_type), intent(inout) :: qm2_params
        type(qm2_structure), intent(inout) :: qm2_struct
        type(qmmm_struct_type), intent(inout) :: qmmm_struct

        _REAL_, intent(inout) :: fock(:)
        _REAL_, intent(in) :: density(:)

! local
        integer::i, j, iqm, jmm, norbital
        _REAL_::fock_opnq, eOPNQ, LJ
        _REAL_::fock_opnq_pair, eOPNQ_pair, LJ_pair
        _REAL_, dimension(MaxValenceOrbitals)::fock_atom_diag, density_atom_diag

!  Check if initilization is necessary
        if (qmmm_opnq%opnq_initialized) then
            if (.not. allocated(qmmm_opnq%type_list_saved)) then
                qmmm_opnq%opnq_initialized = .false.
            else
                if (size(qmmm_opnq%type_list_saved) /= size(qmmm_opnq%MM_atomType)) then
                    qmmm_opnq%opnq_initialized = .false.
                else
                    do i = 1, size(qmmm_opnq%type_list_saved)
                        if (qmmm_opnq%type_list_saved(i) /= qmmm_opnq%MM_atomType(i)) then
                            qmmm_opnq%opnq_initialized = .false.
                            exit
                        end if ! (qmmm_opnq%type_list_saved /= qmmm_opnq%MM_atomType(i)
                    end do ! i=1, size(qmmm_opnq%type_list_saved)
                end if  ! (size(qmmm_opnq%type_list_saved) /= size(qmmm_opnq%MM_atomType)
            end if  !   (.not.allocated(qmmm_opnq%type_list_saved))
        end if !  (qmmm_opnq%opnq_initialized)

        if (.not. qmmm_opnq%opnq_initialized) call Initialize(qmmm_opnq)

        eOPNQ = zero
        LJ = zero

        do iqm = 1, qmmm_struct%nquant
            norbital = qm2_params%natomic_orbs(iqm)

            fock_opnq = zero
            do i = qm2_params%orb_loc(1, iqm), qm2_params%orb_loc(2, iqm)
                j = qm2_params%pascal_tri2(i)
                density_atom_diag(i) = density(j)
            end do ! i

            do jmm = 1, qmmm_struct%qm_mm_pairs
                call Opnq_fock_atom_pair(qmmm_opnq, qm2_params, qm2_struct, qmmm_struct, iqm, jmm, eOPNQ_pair, fock_opnq_pair)
                call Opnq_LJ_atom_pair(qmmm_opnq, qm2_params, qm2_struct, qmmm_struct, iqm, jmm, LJ_pair)
                eOPNQ = eOPNQ + eOPNQ_pair
                fock_opnq = fock_opnq + fock_opnq_pair
                LJ = LJ + LJ_pair
            end do ! jmm

            do i = qm2_params%orb_loc(1, iqm), qm2_params%orb_loc(2, iqm)
                j = qm2_params%pascal_tri2(i)
                fock(j) = fock(j) + fock_opnq
            end do ! i

        end do ! iqm

        qmmm_opnq%OPNQCorrection = eOPNQ
        qmmm_opnq%vdWCorrection = -LJ

        return

    end subroutine Opnq_fock

    subroutine Opnq_fock_atom_pair(qmmm_opnq, qm2_params, qm2_struct, qmmm_struct, iqm, jmm, eOPNQ_pair, fock_opnq_pair, dx, dy, dz)
!***********************************************************************
!
!  This subroutine calculates the OPNQ contributions to the Fock matrix
! for a single quantum atom due to an MM atom.
! the output unit is eV.
!
!  Taisung Lee, Rutgers, 2011
!
!***********************************************************************
        use constants, only : A_TO_BOHRS, AU_TO_EV, zero
        use qmmm_module, only : qm2_structure, qmmm_opnq_structure, MM_opnq
        use qm2_params_module, only : qm2_params_type
        use QM2_parameters, only : core_chg
        use opnq_switching, only : switchoff

        use qmmm_struct_module, only : qmmm_struct_type
        implicit none
        type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq
        type(qm2_params_type), intent(inout) :: qm2_params
        type(qm2_structure), intent(inout) :: qm2_struct
        type(qmmm_struct_type), intent(inout) :: qmmm_struct

        integer, intent(in)::iqm, jmm
        _REAL_, intent(out) :: fock_opnq_pair, eOPNQ_pair
        _REAL_, intent(out), optional::dx, dy, dz

!local
        integer::i, jmm_index, qmtype, mmtype, atomic_number
        _REAL_:: qm_charge, r2, rij, rijInAu, core_charge
        _REAL_::  ee, dEdQi, dEdQj, switching
        type(MM_opnq)::myOpnq

        fock_opnq_pair = zero
        eOPNQ_pair = zero
        if (present(dx) .and. present(dy) .and. present(dz)) then
            dx = zero; dy = zero; dz = zero
        end if

        qmType = qmmm_struct%qm_atom_type(iqm)
        if (qm2_params%qxd_supported(qmtype)) then

            ! calculate the effective charge for the qmatom
            call qm2_calc_mulliken(qm2_params, qm2_struct, iqm, qm_charge)

            jmm_index = qmmm_struct%qm_mm_pair_list(jmm)
            mmtype = qmmm_opnq%MM_atomType(jmm_index)
            if (qmmm_opnq%supported(mmtype)) then
                myOpnq = qmmm_opnq%MM_opnq_list_saved(mmtype)
                atomic_number = qmmm_opnq%atomic_number(mmtype)
                core_charge = core_chg(atomic_number)*1.d0

                r2 = sum((qmmm_struct%qm_xcrd(1:3, jmm) - qmmm_struct%qm_coords(1:3, iqm))**2)
                rij = sqrt(r2)
                rijInAu = rij*A_TO_BOHRS

                if (qmmm_opnq%switching) then
                    call switchoff(rij, qmmm_opnq%switch_cutoff1, qmmm_opnq%switch_cutoff2, switching)
                else
                    switching = 1.0d0
                end if

                call vdw_ij(qm_charge, qm2_params%qxd_s(qmtype), &  ! qm atom
                    qm2_params%qxd_z0(qmtype), qm2_params%qxd_zq(qmtype), &
                    qm2_params%qxd_d0(qmtype), qm2_params%qxd_dq(qmtype), &
                    qm2_params%qxd_q0(qmtype), qm2_params%qxd_qq(qmtype), &
                    qm2_params%qxd_neff(qmtype), qm2_params%core_chg(iqm), &
                    qmmm_struct%qm_xcrd(4, jmm), myopnq%s, myopnq%zeta, zero, &  ! mm atom
                    myopnq%alpha, zero, zero, zero, myopnq%neff, core_charge, &
                    rijInAu, ee, dedqi, dedqj)

                eOPNQ_pair = ee*switching*AU_TO_EV
                fock_opnq_pair = -dEdQi*switching*AU_TO_EV

                if (present(dx) .and. present(dy) .and. present(dz)) then
                    call vdw_ij_dri(qm_charge, qm2_params%qxd_s(qmtype), &  ! qm atom
                        qm2_params%qxd_z0(qmtype), qm2_params%qxd_zq(qmtype), &
                        qm2_params%qxd_d0(qmtype), qm2_params%qxd_dq(qmtype), &
                        qm2_params%qxd_q0(qmtype), qm2_params%qxd_qq(qmtype), &
                        qm2_params%qxd_neff(qmtype), qm2_params%core_chg(iqm), &
                        qmmm_struct%qm_xcrd(4, jmm), myopnq%s, myopnq%zeta, zero, &  ! mm atom
                        myopnq%alpha, zero, zero, zero, myopnq%neff, core_charge, &
                        qmmm_struct%qm_coords(1, iqm)*A_TO_BOHRS, &
                        qmmm_struct%qm_coords(2, iqm)*A_TO_BOHRS, &
                        qmmm_struct%qm_coords(3, iqm)*A_TO_BOHRS, &
                        qmmm_struct%qm_xcrd(1, jmm)*A_TO_BOHRS, &
                        qmmm_struct%qm_xcrd(2, jmm)*A_TO_BOHRS, &
                        qmmm_struct%qm_xcrd(3, jmm)*A_TO_BOHRS, &
                        dx, dy, dz)
                    dx = dx*AU_TO_EV*switching*A_TO_BOHRS
                    dy = dy*AU_TO_EV*switching*A_TO_BOHRS
                    dz = dz*AU_TO_EV*switching*A_TO_BOHRS
                end if

                continue

            end if

        end if ! (qm2_params%qxd_supported(qmtype))

        return

    end subroutine Opnq_fock_atom_pair

    subroutine Opnq_LJ_atom_pair(qmmm_opnq, qm2_params, qm2_struct, qmmm_struct, iqm, jmm, LJ_pair, dx, dy, dz)
!***********************************************************************
!
!  This subroutine calculates the classic LJ interactions on a single
! QM atom due to an MM atom.
! the output unit is eV.
!
!  Taisung Lee, Rutgers, 2011
!
!***********************************************************************
        use constants, only : KCAL_TO_EV, zero
        use qmmm_module, only : qm2_structure, qmmm_opnq_structure
        use opnq_switching, only : switchoff
        use qm2_params_module, only : qm2_params_type

        use qmmm_struct_module, only : qmmm_struct_type
        implicit none
        type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq
        type(qm2_params_type), intent(inout) :: qm2_params
        type(qm2_structure), intent(inout) :: qm2_struct
        type(qmmm_struct_type), intent(in) :: qmmm_struct

        integer, intent(in)::iqm, jmm
        _REAL_, intent(out) :: LJ_pair
        _REAL_, intent(out), optional::dx, dy, dz

!local
        integer::j, jmm_index, qmType, mmtype, mmtype_for_iqm
        _REAL_:: r2, rij
        _REAL_::temp1, temp2, temp3, switching

        LJ_pair = zero
        if (present(dx) .and. present(dy) .and. present(dz)) then
            dx = zero; dy = zero; dz = zero
        end if

        qmType = qmmm_struct%qm_atom_type(iqm)
        mmtype_for_iqm = qmmm_opnq%MM_atomType(qmmm_struct%iqmatoms(iqm))
        if (qm2_params%qxd_supported(qmtype)) then

            jmm_index = qmmm_struct%qm_mm_pair_list(jmm)
            mmtype = qmmm_opnq%MM_atomType(jmm_index)
            if (qmmm_opnq%supported(mmtype)) then

                r2 = sum((qmmm_struct%qm_xcrd(1:3, jmm) - qmmm_struct%qm_coords(1:3, iqm))**2)
                rij = sqrt(r2)

                ! classic LJ interaction
                temp1 = 0.5d0*(qmmm_opnq%LJ_r(mmtype_for_iqm) + qmmm_opnq%LJ_r(mmtype))
                temp2 = sqrt(qmmm_opnq%LJ_epsilon(mmtype_for_iqm)*qmmm_opnq%LJ_epsilon(mmtype))

                if (qmmm_opnq%switching) then
                    call switchoff(rij, qmmm_opnq%switch_cutoff1, qmmm_opnq%switch_cutoff2, switching)
                else
                    switching = 1.0d0
                end if

                LJ_pair = 4096.D0*(temp1**12)*temp2/(rij**12) - 128.D0*(temp1**6)*temp2/(rij**6)
                LJ_pair = LJ_pair*switching*KCAL_TO_EV

                if (present(dx) .and. present(dy) .and. present(dz)) then
                    temp3 = -12.D0*4096.D0*(temp1**12)*temp2/(rij**13) + 6.D0*128.D0*(temp1**6)*temp2/(rij**7)
                    dx = temp3*(qmmm_struct%qm_xcrd(1, jmm) - qmmm_struct%qm_coords(1, iqm))/rij
                    dy = temp3*(qmmm_struct%qm_xcrd(2, jmm) - qmmm_struct%qm_coords(2, iqm))/rij
                    dz = temp3*(qmmm_struct%qm_xcrd(3, jmm) - qmmm_struct%qm_coords(3, iqm))/rij
                    dx = -dx*switching*KCAL_TO_EV
                    dy = -dy*switching*KCAL_TO_EV
                    dz = -dz*switching*KCAL_TO_EV
                end if

            end if

        end if ! (qm2_params%qxd_supported(qmtype))

        return

    end subroutine Opnq_LJ_atom_pair

    subroutine LJ2OPNQ(atomic_number, sigma, epsilon, MM_entry)

        use constants, only : AU_TO_KCAL, A_TO_BOHRS
        use qmmm_module, only : MM_opnq
        implicit none

        integer, intent(in)::atomic_number
        _REAL_, intent(in)::sigma, epsilon
        type(MM_opnq), intent(out)::mm_entry

        _REAL_, parameter::A = 9.442292940115D0           !##############################
        _REAL_, parameter::B = 0.411072836373204D0        !#   These Mapping Coeff might
        _REAL_, parameter::C = 2.82083644658535D0         !# change in the future b/c
        _REAL_, parameter::D = 3.78925426936423D0         !# initial parameterization
        _REAL_, parameter::E = -0.0192103969646589D0      !# space could possibly change
        _REAL_, parameter::F = -0.724935124427059D0       !##############################

        _REAL_, parameter, dimension(MaxAtomicNumber_MM_opnq)::G_data = (/ &
            0.00395D0, 0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, & ! 1-5
            0.20636D0, 0.18738D0, 0.17208D0, 0.00000D0, 0.00000D0, & ! 6-10
            0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, & ! 11-15
            0.00000D0, 0.18944D0, 0.00000D0, 0.00000D0, 0.00000D0/) ! 16-20

        _REAL_, parameter, dimension(MaxAtomicNumber_MM_opnq)::neff_data = (/ &
            0.82400D0, 0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, & ! 1-5
            2.65700D0, 3.18700D0, 3.66300D0, 0.00000D0, 0.00000D0, & ! 6-10
            0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, 0.00000D0, & ! 11-15
            0.00000D0, 5.55100D0, 0.00000D0, 0.00000D0, 0.00000D0/) ! 16-20

        _REAL_::temp1, temp2, G, neff

        neff = neff_data(atomic_number)
        G = G_data(atomic_number)

        MM_entry%s = A*(epsilon**B)*(sigma**C)
        MM_entry%zeta = D*(epsilon**E)*(sigma**F)
        temp1 = 512.d0*(1.d0 - G)*(sigma*A_TO_BOHRS)**6.d0
        temp1 = temp1*epsilon/AU_TO_KCAL
        temp2 = (3.d0*sqrt(neff))

        MM_entry%alpha = (temp1/temp2)**(2.d0/3.d0)
        MM_entry%neff = neff_data(atomic_number)

        return

    end subroutine LJ2OPNQ

    subroutine Initialize(qmmm_opnq)

        use qmmm_module, only : qmmm_opnq_structure, MM_opnq
        implicit none
        type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq
        integer::i, natom, ntype
        type(MM_opnq)::temp

        natom = size(qmmm_opnq%MM_atomType)
        ntype = size(qmmm_opnq%LJ_r)

        if (allocated(qmmm_opnq%type_list_saved)) deallocate (qmmm_opnq%type_list_saved)
        allocate (qmmm_opnq%type_list_saved(natom))
        qmmm_opnq%type_list_saved = qmmm_opnq%MM_atomType

        if (allocated(qmmm_opnq%MM_opnq_list_saved)) deallocate (qmmm_opnq%MM_opnq_list_saved)
        allocate (qmmm_opnq%MM_opnq_list_saved(ntype))

        do i = 1, ntype
            temp%S = 0.d0
            temp%zeta = 0.d0
            temp%alpha = 0.d0
            temp%neff = 0.d0
            if (qmmm_opnq%atomic_number(i) >= 1 .and. &
                qmmm_opnq%atomic_number(i) <= MaxAtomicNumber_MM_opnq) then
                call LJ2OPNQ(qmmm_opnq%atomic_number(i), &
                    qmmm_opnq%LJ_r(i), qmmm_opnq%LJ_epsilon(i), temp)
            end if
            qmmm_opnq%MM_opnq_list_saved(i) = temp
        end do ! i=1, ntype

        qmmm_opnq%opnq_initialized = .true.

        return

    end subroutine Initialize

end module opnq
