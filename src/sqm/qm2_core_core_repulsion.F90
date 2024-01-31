! <compile=optimized>
#include "dprec.fh"
subroutine qm2_core_core_repulsion(qmmm_nml, qm2_params, qmmm_struct, iat, jat, rij, onerij, RI1, enuc)
! ----------------------------------------------------------------------
! PURPOSE: Calculate the core core repulsion energy between
!          atoms iat, jat
!
! E^{\rm core}_{i,j} = Z_i Z_j (s_i s_i|s_j s_j) [1 + f_{i,j}] + g_{i,j}
!
! MNDO   : f_{i,j} = e^{-\alpha_i R_{i,j}} + e^{-alpha_j R_{i,j}}
!          g_{i,j} = 0
!
! AM1/PM3: f_{i,j} = f^{\rm MNDO}_{i,j}
!          g_{i,j} = (Z_i Z_j)/R_{i,j}
!                              \sum_k a_k^i e^{-b_k^i (R_{i,j}-c_k^i)^2}
!
! PM6    : f_{i,j} = x_{ij}e^{-alpha_{ij}(R_ij+0.0003R_{ij}^6)}
!          with special cases for atom pairs of the type
!          OH, NH, CH, CC, SiO
!          g_{i,j} = g^{\rm AM1}_{i,j}
!                      + 10^{-8}[(Z_i^{1/3}+Z_j^{1/3})/R_{ij}]^{12}
!
! PM3-MAIS: f_{i,j} = f^{\rm MNDO}_{i,j}
!           g_{i,j} = sum_k a_k(ij) e^{-b_K(ij) [(R_ij-c_k(ij)^2]}
!
! The internuclear distance and its inverse are passed down for
! efficiency reasons. They have been calculated before for the
! evaluation of the ERI (s_i s_i|s_j s_j)
!
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : February 2010
! ----------------------------------------------------------------------

    use qm2_params_module, only : qm2_params_type
    use qmmm_nml_module, only : qmmm_nml_type
    use qmmm_module, only : EXPONENTIAL_CUTOFF
    use constants, only : zero, one, two, ten, ten_to_minus8, third, AU_TO_EV, BOHRS_TO_A
    use qmmm_struct_module, only : qmmm_struct_type

!DEBUG use qmmm_qmtheorymodule, only : String
    implicit none

    type(qmmm_struct_type), intent(in) :: qmmm_struct
    type(qmmm_nml_type), intent(in) :: qmmm_nml
    type(qm2_params_type), intent(in) :: qm2_params

    integer, intent(in)  :: iat, jat
    _REAL_, intent(in)  :: rij, onerij, RI1
    _REAL_, intent(out) :: enuc

    integer :: i, j
    integer :: iatyp, jatyp
    integer :: iatnum, jatnum, ijmin, ijmax
    _REAL_  :: alpi, alpj, scale, ccij, giijj
    _REAL_  :: anam1, tmp
    _REAL_  :: pddg_exp1, pddg_exp2, pddg_exp3, pddg_exp4, pddg_correction
    _REAL_  :: upc
    _REAL_, parameter :: threemin4 = 3.0d-04
    _REAL_, parameter :: ccfac = 9.28d0
    _REAL_, parameter :: ccexp = 5.98d0
    _REAL_, parameter :: siofac = 0.7d-03
    _REAL_, parameter :: sioexp = 2.9d0
    _REAL_ :: alpab, xab

    iatyp = qmmm_struct%qm_atom_type(iat)
    jatyp = qmmm_struct%qm_atom_type(jat)
    iatnum = qmmm_struct%iqm_atomic_numbers(iat)
    jatnum = qmmm_struct%iqm_atomic_numbers(jat)
    ijmin = min(iatnum, jatnum)
    ijmax = max(iatnum, jatnum)

    !   use rho_core for core-core interaction if necessary
    if (qmmm_nml%qmtheory%AM1D .or. qmmm_nml%qmtheory%MNDOD .or. qmmm_nml%qmtheory%PM6) then
        giijj = AU_TO_EV/sqrt(rij*rij/(BOHRS_TO_A*BOHRS_TO_A) &
            + (qm2_params%po(9, iatyp) + qm2_params%po(9, jatyp))**2)
    else
        giijj = RI1
    end if

    enuc = zero

    ! -----------------------------------------------------
    ! Calculate scaling factor for monopole interaction
    ! PM6: calculate also correction for unpolarizable core
    ! -----------------------------------------------------
    SWITCHMETHOD: if (qmmm_nml%qmtheory%PM6) then
        !DEBUG write (6,'(a)') 'PM6 scaling factor and correction for unpolarizable core'
        ! ---------------------------------------------------------
        ! Correction for unpolarizable core
        ! Equation (7) of Stewart, J Mol Mod 13 (2007) 1173
        ! Note: instead of the core charge as written in the paper,
        !       Stewart is using the atomic numbers here!
        ! ---------------------------------------------------------
        ! Note: Jimmy Stewart uses 0.3333d0 instead of one third
        !       For dinitrogen this leads to an error of
        !       0.0000765 eV = 0.00176 kcal/mol
        ! upc = onerij * (dble(iatnum)**0.3333d0 + dble(jatnum)**0.3333d0)
        upc = onerij*(dble(iatnum)**third + dble(jatnum)**third)
        upc = ten_to_minus8*upc**12

        enuc = enuc + upc

        ! -------------------------------------------------------------------
        ! PM6 scaling factor
        ! Equations (6), (9), (10), (11) of Stewart, J Mol Mod 13 (2007) 1173
        ! -------------------------------------------------------------------
        alpab = qm2_params%pm6_alpab(iatyp, jatyp)
        xab = qm2_params%pm6_xab(iatyp, jatyp)
        if ((ijmin == 1) .and. ((ijmax == 6) .or. (ijmax == 7) .or. (ijmax == 8))) then
            ! ---------------------------------------------------------
            ! CH, NH and OH special case, equation (9)
            ! Note: Stewart does not mention CH as special case
            !       in the paper, but his code treats CH like NH and OH
            ! ---------------------------------------------------------
            tmp = alpab*rij*rij
        else
            ! -------------------------------------------------
            ! General case, equation (6)
            ! -------------------------------------------------
            tmp = alpab*(rij + threemin4*(rij**6))
        end if
        ! ATTENTION: MOPAC2007 has this multiplied by a factor of two (see mndod.f90)
        ! This factor of two is necessary to match the MOPAC2009 PM6 results
        ! However, this factor of two is not documented in the PM6 paper!
        scale = one + two*xab*exp(-tmp)

        if ((iatnum == 6) .and. (jatnum == 6)) then
            ! -------------------------------------------
            ! CC scaling factor correction, equation (10)
            ! -------------------------------------------
            tmp = ccexp*rij
            tmp = exp(-tmp)
            scale = scale + ccfac*tmp
        else if ((ijmax == 14) .and. (ijmin == 8)) then
            ! --------------------------------------------
            ! SiO scaling factor correction, equation (11)
            ! --------------------------------------------
            tmp = (rij - sioexp)**2
            tmp = exp(-tmp)
            scale = scale - siofac*tmp
        end if

    else SWITCHMETHOD

        ! -------------------
        ! MNDO scaling factor
        ! -------------------
        alpi = exp(-qm2_params%cc_exp_params(iat)*rij)
        if (iatyp == jatyp) then
            alpj = alpi
        else
            alpj = exp(-qm2_params%cc_exp_params(jat)*rij)
        end if
        scale = one + alpi + alpj
        ! account for NH and OH atom pairs
        if ((ijmin == 1) .and. ((ijmax == 7) .or. (ijmax == 8))) then
            if (iatnum == 1) then
                scale = scale + (rij - 1.0D0)*alpj
            else
                scale = scale + (rij - 1.0D0)*alpi
            end if
        end if
    end if SWITCHMETHOD

    ! ---------------------------------------------
    ! Monopole term, multiplied with scaling factor
    ! ---------------------------------------------
    ccij = qm2_params%core_chg(iat)*qm2_params%core_chg(jat)
    tmp = ccij*giijj
    enuc = enuc + tmp*scale

    ! -----------------------------------------------
    ! AM1, PM3, PM6 and PM3-MAIS Gaussian corrections
    ! -----------------------------------------------
    if (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) then
        anam1 = zero
        do i = 1, qm2_params%num_fn(iatyp)
            tmp = rij - qm2_params%FN3(i, iatyp)
            tmp = qm2_params%FN2(i, iatyp)*tmp*tmp
            if (tmp < EXPONENTIAL_CUTOFF) then ! Skip doing the exponential if it is essentially zero
                anam1 = anam1 + qm2_params%FN1(i, iatyp)*exp(-tmp)
            end if
        end do
        if (iatyp == jatyp) then
            anam1 = anam1 + anam1
        else
            do i = 1, qm2_params%num_fn(jatyp)
                tmp = rij - qm2_params%FN3(i, jatyp)
                tmp = qm2_params%FN2(i, jatyp)*tmp*tmp
                if (tmp < EXPONENTIAL_CUTOFF) then ! Skip doing the exponential if it is essentially zero
                    anam1 = anam1 + qm2_params%FN1(i, jatyp)*exp(-tmp)
                end if
            end do
        end if
        anam1 = ccij*onerij*anam1

        ! apply the AM1-d/PhoT correction if existing (i.e., GNN/=1.0)
        if (qmmm_nml%qmtheory%AM1D) then
            anam1 = anam1*qm2_params%GNN(iatyp)*qm2_params%GNN(jatyp)
        end if

        enuc = enuc + anam1
    else if (qmmm_nml%qmtheory%PM3MAIS) then
        anam1 = zero
        do i = 1, 3
            tmp = rij - qm2_params%pm3mais_gamab(iatyp, jatyp, i)
            tmp = qm2_params%pm3mais_betab(iatyp, jatyp, i)*tmp*tmp
            if (tmp < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
                anam1 = anam1 + qm2_params%pm3mais_alpab(iatyp, jatyp, i)*exp(-tmp)
            end if
        end do
        enuc = enuc + anam1
    end if

    ! ---------------
    ! PDDG correction
    ! ---------------
    if (qmmm_struct%PDDG_IN_USE) then
        pddg_exp1 = exp(-ten*(rij - qm2_params%pddge1(iat) - qm2_params%pddge1(jat))**2)
        pddg_exp2 = exp(-ten*(rij - qm2_params%pddge1(iat) - qm2_params%pddge2(jat))**2)
        pddg_exp3 = exp(-ten*(rij - qm2_params%pddge2(iat) - qm2_params%pddge1(jat))**2)
        pddg_exp4 = exp(-ten*(rij - qm2_params%pddge2(iat) - qm2_params%pddge2(jat))**2)
        pddg_correction = qm2_params%pddg_term1(iatyp, jatyp)*pddg_exp1 + &
            qm2_params%pddg_term2(iatyp, jatyp)*pddg_exp2 + &
            qm2_params%pddg_term3(iatyp, jatyp)*pddg_exp3 + &
            qm2_params%pddg_term4(iatyp, jatyp)*pddg_exp4
        enuc = enuc + pddg_correction
    end if

end subroutine qm2_core_core_repulsion
