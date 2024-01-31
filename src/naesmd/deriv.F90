#include "dprec.fh"
#include "assert.fh"

! Subroutine for analytic and numerical derivatives of ground or excited states
subroutine deriv(sim, state) ! , xyz_in)
    use communism !for numerical derivatives
    use qm2_davidson_module
    use qmmm_module
    use constants, only : KCAL_TO_EV, EV_TO_KCAL

    implicit none

    type(simulation_t), target ::sim !communism module
    type(simulation_t), pointer::simpoint
    integer, intent(in), optional::state ! excited state where derivatives are calculated

    integer i, j, k, ihop, oldverbosity
    _REAL_::dxyz(sim%qmmm%nquant_nlink*3)
    _REAL_::dxyz_gs(sim%qmmm%nquant_nlink*3)
    _REAL_::dxyz1(3, sim%qmmm%nquant_nlink), dxyz1_test(3, sim%qmmm%nquant_nlink)
    _REAL_::xyz(3*sim%qmmm%nquant_nlink)
    _REAL_::charges2(sim%cosmo%nps), acharges2(sim%cosmo%numat), density2(sim%cosmo%lm61)
    _REAL_ :: Escf_right, Escf_left, E_ES_right, E_ES_left, h
    simpoint => sim
    !Collect coordinates in vector
    do i = 1, simpoint%qmmm%nquant_nlink
        do j = 1, 3
            xyz((i - 1)*3 + j) = simpoint%qmmm%qm_coords(j, i)
        end do
    end do

    !Determine state to for which to calculate derivatives
    if (present(state)) then ! if different state is specified on input
        ihop = state
    else
        ihop = simpoint%qmmm%state_of_interest !otherwise use the current state of interest
    end if

    if (simpoint%qmmm%ideriv .eq. 1) then !analytical derivatives

        !CALCULATE GROUND STATE DERIVATIVES !At some point these should be saved and
        !recalled for speed and convenience since deriv may be called multiple times for the
        !same rho_0 and thus repeated
        dxyz1 = 0.d0; dxyz1_test = 0.d0
        dxyz = 0.d0
        dxyz_gs = 0.d0

        ! Calculate ground state derivatives E_gr^x=E_nucl^x+E_el^x
        !   E_el^x=1/2 Tr((t^x+F^x) rho)
        ! qm2_get_exc_forces() is the SQM equivalent of DCART() in CEO

        call qm2_get_exc_forces(simpoint%qparams, simpoint%qnml, simpoint%rij, simpoint%qmpi, &
            simpoint%qm2, simpoint%qmmm, dxyz1, simpoint%qmmm%qm_coords)

        !add solvent part (nuclear and electronic for ground state with COSMO
        !surface derivatives)
        if ((simpoint%cosmo%solvent_model .gt. 0) .and. (simpoint%cosmo%solvent_model .ne. 10)) then
            if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then !have to specify simpoint%cosmo%ceps here because of division by 0 when eq 1
                call cosmo_1_tri(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2%den_matrix) !put simpoint%dav%den_matrix in the right place
                call diegrd(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test) !derivative
            elseif (simpoint%cosmo%potential_type .eq. 2) then
                call rcnfldgrad(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                    simpoint%qm2%den_matrix, simpoint%qm2%norbs)
            end if
            dxyz1 = dxyz1 + dxyz1_test
        end if

        do i = 1, simpoint%qmmm%nquant_nlink
            do j = 1, 3
                dxyz_gs((i - 1)*3 + j) = -dxyz1(j, i)*KCAL_TO_EV
            end do
        end do

        if (ihop > 0) then
            !CALCULATE EXCITED STATE DERIVATIVES

            !TERM 1: Tr(F^x rhoTZ)

            ! Need to call dealloc_rhotz() when finished or in polishing
            ! add allocation to the big allocation/dealloc regime

            ! Get excited state density
            call calc_rhotz(simpoint%qparams, simpoint%qnml, simpoint%qmpi, simpoint%cosmo, simpoint%qm2, &
                simpoint%dav, simpoint%qmmm, ihop, simpoint%dav%rhoTZ, .true.)
            call mo2sitef(simpoint%dav%Nb, simpoint%dav%vhf, simpoint%dav%rhoTZ, simpoint%dav%tz_scratch(1), &
                simpoint%dav%tz_scratch(simpoint%dav%Nb**2 + 1))
            call packing(simpoint%dav%Nb, simpoint%dav%tz_scratch(1), simpoint%dav%rhoTZ, 's')

            !Get transition density
            call getmodef(simpoint%dav%Np*simpoint%dav%Nh*2, simpoint%dav%Mx, simpoint%dav%Np, simpoint%dav%Nh, &
                ihop, simpoint%dav%v0, simpoint%dav%tz_scratch(1))
            call mo2sitef(simpoint%dav%Nb, simpoint%dav%vhf, simpoint%dav%tz_scratch(1), simpoint%dav%rhoLZ, &
                simpoint%dav%tz_scratch(simpoint%dav%Nb**2 + 1))
            dxyz1 = 0.d0; dxyz1_test = 0.d0;
            !calculate vacuum derivative for term 1
            call dcart1(simpoint%qnml, simpoint%qparams, simpoint%rij, simpoint%qmpi, simpoint%qm2, &
                simpoint%dav, simpoint%qmmm, dxyz1, &
                simpoint%qm2%den_matrix, simpoint%dav%rhoTZ, simpoint%qmmm%qm_coords)
            !add solvent part (symmetric only b/c symmetric matrix)
            if (simpoint%cosmo%solvent_model .gt. 0) then
                if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then !simpoint%cosmo%ceps.gt.1.0 because of singularity in cosmo subroutines
                    call cosmo_1_tri_2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%rhoTZ, density2, charges2, acharges2) !solvent and solute charges
                    call diegrd2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test, density2, charges2, acharges2) !derivative
                elseif (simpoint%cosmo%potential_type .eq. 2) then
                    call rcnfldgrad2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                        simpoint%qm2%den_matrix, simpoint%dav%rhoTZ, simpoint%dav%nb, .true.)
                end if
                dxyz1 = dxyz1 + dxyz1_test
            end if

            do i = 1, simpoint%qmmm%nquant_nlink
                do j = 1, 3
                    dxyz((i - 1)*3 + j) = dxyz((i - 1)*3 + j) - dxyz1(j, i)*KCAL_TO_EV
                end do
            end do

            !TERM 2: Tr(V^x(xi) xi^+)
            !Symmetric part
            dxyz1 = 0.d0; dxyz1_test = 0.d0; charges2 = 0.d0; acharges2 = 0.d0;
            !pack triangular of symmetric part and antisymmetric parts
            call packing(simpoint%dav%Nb, simpoint%dav%rhoLZ, simpoint%dav%tz_scratch(1), 's')
            call packing(simpoint%dav%Nb, simpoint%dav%rhoLZ, simpoint%dav%tz_scratch(simpoint%dav%nb**2 + 1), 'u')
            !vacuum derivative
            call dcart2(simpoint%qparams, simpoint%qnml, simpoint%rij, simpoint%qm2, simpoint%dav, &
                simpoint%qmmm, dxyz1, simpoint%dav%tz_scratch(1), simpoint%qmmm%qm_coords)
            call dcart2(simpoint%qparams, simpoint%qnml, simpoint%rij, simpoint%qm2, simpoint%dav, &
                simpoint%qmmm, dxyz1, simpoint%dav%tz_scratch(simpoint%dav%nb**2 + 1), &
                simpoint%qmmm%qm_coords)
            if (simpoint%cosmo%solvent_model .eq. 1) then
                if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then
                    simpoint%cosmo%qscnet(:, 1) = 0.d0; simpoint%cosmo%qdenet(:, 1) = 0.d0; !Clear Nuclear Charges
                    call cosmo_1_tri(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%tz_scratch(1)) !Fill Electronic Chrages
                    call diegrd(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test); !derivative
                elseif (simpoint%cosmo%potential_type .eq. 2) then
                    call rcnfldgrad_full(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                        simpoint%dav%rhoLZ, simpoint%dav%nb);
                end if
            end if
            dxyz1 = dxyz1 + 0.5*dxyz1_test

            do i = 1, simpoint%qmmm%nquant_nlink
                do j = 1, 3
                    dxyz((i - 1)*3 + j) = dxyz((i - 1)*3 + j) - dxyz1(j, i)*KCAL_TO_EV
                end do
            end do

            !STATE SPECIFIC SOLVENT TERMS
            dxyz1 = 0.d0; dxyz1_test = 0.d0; charges2 = 0.d0; acharges2 = 0.d0; density2 = 0.d0
            !Vertical Excitation Model
            if (simpoint%cosmo%solvent_model .eq. 2) then
                !Get unrelaxed difference density matrix for the state to calculate derivatives for
                call calc_rhotz(simpoint%qparams, simpoint%qnml, simpoint%qmpi, simpoint%cosmo, simpoint%qm2, simpoint%dav, &
                    simpoint%qmmm, ihop, simpoint%dav%rhoT, .false.)
                call mo2sitef(simpoint%dav%Nb, simpoint%dav%vhf, simpoint%dav%rhoT, simpoint%dav%tz_scratch(1), &
                    simpoint%dav%tz_scratch(simpoint%dav%Nb**2 + 1))
                call packing(simpoint%dav%Nb, simpoint%dav%tz_scratch(1), simpoint%dav%rhoT, 's')
                !Get relaxed or unrelaxed difference density matrix for the state specific state
                call calc_rhotz(simpoint%qparams, simpoint%qnml, simpoint%qmpi, simpoint%cosmo, simpoint%qm2, simpoint%dav, &
                    simpoint%qmmm, simpoint%qmmm%state_of_interest, simpoint%dav%rhoTZ, &
                    simpoint%cosmo%doZ);
                call mo2sitef(simpoint%dav%Nb, simpoint%dav%vhf, simpoint%dav%rhoTZ, simpoint%dav%tz_scratch(1), &
                    simpoint%dav%tz_scratch(simpoint%dav%Nb**2 + 1))
                call packing(simpoint%dav%Nb, simpoint%dav%tz_scratch(1), simpoint%dav%rhoTZ, 's')

                !calculate derivatives
                if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then !simpoint%cosmo%ceps.gt.1.0 because of glitch in cosmo subroutines
                    simpoint%cosmo%qscnet(:, 1) = 0.d0; simpoint%cosmo%qdenet(:, 1) = 0.d0; !Clear Nuclear Charges
                    call cosmo_1_tri(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%rhoTZ) !fill solvent charges
                    call cosmo_1_tri_2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%rhoT, density2, charges2, acharges2) !fill solute charges
                    call diegrd2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test, density2, charges2, acharges2) !derivative
                elseif (simpoint%cosmo%potential_type .eq. 2) then
                    call rcnfldgrad2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                        simpoint%dav%rhoTZ, simpoint%dav%rhoT, simpoint%dav%nb, .false.)
                end if
                dxyz1 = dxyz1 + 0.5*dxyz1_test

                do i = 1, simpoint%qmmm%nquant_nlink
                    do j = 1, 3
                        dxyz((i - 1)*3 + j) = dxyz((i - 1)*3 + j) - dxyz1(j, i)*KCAL_TO_EV
                    end do
                end do
            end if

            !State Specific Model
            if (simpoint%cosmo%solvent_model .eq. 4) then
                write (6, *) 'WARNING:DERIVATIVES FOR STATE SPECIFIC SOLVENT ARE NONVARIATIONAL'
                !Get relaxed density matrix for the state specific state
                simpoint%dav%rhoT = 0; simpoint%dav%rhoTZ = 0;
                call calc_rhotz(simpoint%qparams, simpoint%qnml, simpoint%qmpi, &
                    simpoint%cosmo, simpoint%qm2, simpoint%dav, simpoint%qmmm, &
                    simpoint%qmmm%state_of_interest, simpoint%dav%rhoT, simpoint%cosmo%doZ); !rhoT will be rhoTZ_k
                call mo2sitef(simpoint%dav%Nb, simpoint%dav%vhf, simpoint%dav%rhoT, simpoint%dav%tz_scratch(1), &
                    simpoint%dav%tz_scratch(simpoint%dav%Nb**2 + 1))
                simpoint%dav%rhoT = 0;
                call packing(simpoint%dav%Nb, simpoint%dav%tz_scratch(1), simpoint%dav%rhoT, 's')

                !Excited State Part
                if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then !simpoint%cosmo%ceps.gt.1.0 because of singularity in cosmo subroutines
                    call cosmo_1_tri(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%rhoT) !fill solvent charges
                    call cosmo_1_tri_2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%dav%rhoTZ, density2, charges2, acharges2) !fill solute charges
                    call diegrd2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test, density2, charges2, acharges2) !derivative
                elseif (simpoint%cosmo%potential_type .eq. 2) then
                    call rcnfldgrad2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                        simpoint%dav%rhoT, simpoint%dav%rhoTZ, simpoint%dav%nb, .false.)
                end if
                dxyz1 = dxyz1 + dxyz1_test

                do i = 1, simpoint%qmmm%nquant_nlink
                    do j = 1, 3
                        dxyz((i - 1)*3 + j) = dxyz((i - 1)*3 + j) - dxyz1(j, i)*KCAL_TO_EV
                    end do
                end do

                dxyz1 = 0.d0; dxyz1_test = 0.d0; charges2 = 0.d0; acharges2 = 0.d0; density2 = 0.d0
                !Ground State part
                if ((simpoint%cosmo%potential_type .eq. 3) .and. (simpoint%cosmo%ceps .gt. 1.0)) then !simpoint%cosmo%ceps.gt.1.0 because of singularity in cosmo subroutines
                    call cosmo_1_tri_2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2%den_matrix, &
                        density2, charges2, acharges2) !fill solute charges
                    call diegrd2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qmmm, dxyz1_test, &
                        density2, charges2, acharges2) !derivative
                elseif (simpoint%cosmo%potential_type .eq. 2) then
                    call rcnfldgrad2(simpoint%qparams, simpoint%qnml, simpoint%cosmo, simpoint%qm2, simpoint%qmmm, dxyz1_test, &
                        simpoint%dav%rhoT, simpoint%qm2%den_matrix, simpoint%dav%nb, .true.)
                end if
                dxyz1 = dxyz1 + dxyz1_test

                do i = 1, simpoint%qmmm%nquant_nlink
                    do j = 1, 3
                        dxyz_gs((i - 1)*3 + j) = dxyz_gs((i - 1)*3 + j) - dxyz1(j, i)*KCAL_TO_EV
                    end do
                end do

            end if

        end if !ihop>0

        !NUMERICAL DERIVATIVES
        !Currently wastes many resources by running full calculations for each
        !state. Need to store derivatives for all states since they are
        !calculated each time instead of calling deriv() for each state.

    elseif (simpoint%qmmm%ideriv .eq. 2) then
        h = simpoint%qmmm%numder_step !step size
        oldverbosity = simpoint%qnml%verbosity
        simpoint%qnml%verbosity = 0
        E_ES_left = 0; E_ES_right = 0;
        !Loop over the number of coordinates
        do i = 1, simpoint%qmmm%nquant_nlink*3
            !Modify the coordinates and calculate energies

            xyz(i) = xyz(i) + h !left
            do k = 1, simpoint%qmmm%nquant_nlink ! this can be eliminated by reorganizing xyz and dxyz
                do j = 1, 3
                    simpoint%qmmm%qm_coords(j, k) = xyz((k - 1)*3 + j)/0.529177249d0 !
                end do
            end do

            call do_sqm_davidson_update(simpoint, 0, vgs=Escf_left, rx=simpoint%qmmm%qm_coords(1, :) &
                , ry=simpoint%qmmm%qm_coords(2, :), rz=simpoint%qmmm%qm_coords(3, :))
            !Note: rx=___ etc. might already be uncessary because qm_coords might be used in the subroutine
            E_ES_left = sim%naesmd%Omega(ihop)

            xyz(i) = xyz(i) - 2*h !right
            do k = 1, simpoint%qmmm%nquant_nlink
                do j = 1, 3
                    simpoint%qmmm%qm_coords(j, k) = xyz((k - 1)*3 + j)/0.529177249d0 !
                end do
            end do

            call do_sqm_davidson_update(simpoint, 0, vgs=Escf_right, rx=simpoint%qmmm%qm_coords(1, :) &
                , ry=simpoint%qmmm%qm_coords(2, :), rz=simpoint%qmmm%qm_coords(3, :))
            E_ES_right = sim%naesmd%Omega(ihop)

            xyz(i) = xyz(i) + h !back to center
            do k = 1, simpoint%qmmm%nquant_nlink
                do j = 1, 3
                    simpoint%qmmm%qm_coords(j, k) = xyz((k - 1)*3 + j)/0.529177249d0 !
                end do
            end do

            call do_sqm_davidson_update(simpoint, 0, rx=simpoint%qmmm%qm_coords(1, :) &
                , ry=simpoint%qmmm%qm_coords(2, :), rz=simpoint%qmmm%qm_coords(3, :))

            !Calculate derivative
            dxyz_gs(i) = -(Escf_left - Escf_right)/(2*h)*27.2116 !GS
            if (ihop > 0) dxyz(i) = -(E_ES_left - E_ES_right)/(2*h)*27.2116 !ES
            write (6, *) 'ES:', E_ES_left, E_ES_right, 'dxyz(i)=', dxyz(i)
        end do

        simpoint%dav%verbosity = oldverbosity
        simpoint%qnml%verbosity = oldverbosity

    end if !deriv flag

    !WRITE RESULTS
    if (simpoint%qmpi%commqmmm_master .AND. simpoint%qnml%verbosity > 3) then
        !If verbosity level is greater than 3 we also print the force array on the
        !QM atoms
        write (6, '("QMMM: Forces on QM atoms ground state calculation (eV/A)")')
        write (6, '("QMMM: state=",2i3)') ihop, simpoint%dav%struct_opt_state
        write (6, '("QMMM: Atm ",i6,": ",3f20.14)') (j, -dxyz_gs(1 + 3*(j - 1)), &
            -dxyz_gs(2 + 3*(j - 1)), &
            -dxyz_gs(3 + 3*(j - 1)), j=1, simpoint%qmmm%nquant_nlink)
        if (ihop > 0) then
            write (6, '("QMMM: Forces on QM atoms excited state calculation (eV/A)")')
            write (6, '("QMMM: state=",2i3)') ihop, simpoint%dav%struct_opt_state
            write (6, '("QMMM: Atm ",i6,": ",3f20.14)') (j, -dxyz(1 + 3*(j - 1)), &
                -dxyz(2 + 3*(j - 1)), &
                -dxyz(3 + 3*(j - 1)), j=1, simpoint%qmmm%nquant_nlink)

            if (simpoint%qnml%verbosity > 4) then
                !Also print info in KJ/mol
                !write (6,'("QMMM:")')
                write (6, '("QMMM: Excited State Forces on QM atoms from SCF calculation (KJ/mol)")')
                write (6, '("QMMM: Atm ",i6,": ",3f20.14)') (j, &
                    -dxyz(1 + 3*(j - 1))*EV_TO_KCAL*4.184d0, &
                    -dxyz(2 + 3*(j - 1))*EV_TO_KCAL*4.184d0, &
                    -dxyz(3 + 3*(j - 1))*EV_TO_KCAL*4.184d0, j=1, simpoint%qmmm%nquant_nlink)
                write (6, '(/,''  QMMM: QM Region Cartesian Coordinates (*=link atom) '')')
                write (6, '(''  QMMM: QM_NO. MM_NO.'',2X,''ATOM'',9X,''X'',9X,''Y'',9X,''Z'')')
                do I = 1, simpoint%qmmm%nquant
                    write (6, '("  QMMM:",I6,2X,I7,6X,I2,3X,3F10.4)') I, &
                        simpoint%qmmm%iqmatoms(i), &
                        I, &
                        (simpoint%qmmm%qm_coords(J, I), J=1, 3)
                end do
            end if
        end if
    end if
    simpoint%deriv_forces = dxyz_gs + dxyz !incorporating into sim type
    return
end subroutine deriv
!
! Subroutine for analytic and numerical derivatives of ALL ground and excited states
! Used for Mean-Field Propagation
subroutine deriv_MF(sim, restart_flag)
    use communism
    use AIMC_type_module, only : AIMC_type

    implicit none
    interface
        subroutine deriv(sim, state)
            use communism !for numerical derivatives
            use qm2_davidson_module
            use qmmm_module
            use constants, only : KCAL_TO_EV, EV_TO_KCAL
            use naesmd_constants
            type(simulation_t), target ::sim
            integer, intent(in), optional::state
        end subroutine deriv
    end interface

    type(simulation_t), target ::sim !communism module
    type(simulation_t), pointer::simpoint

    type(naesmd_structure), pointer::namd
    type(AIMC_type), pointer::aimc
    _REAL_ vtemp(sim%excN), vgstemp
    _REAL_ NACR(sim%naesmd%natom*3)
    integer, intent(in) :: restart_flag !Flag for restarting (1) or not (0)
    integer ::state, state1, state2, n, i, j ! excited state where derivatives are calculated
    simpoint => sim
    namd => sim%naesmd
    aimc => sim%aimc
    n = simpoint%excN

    do state = 1, n
        call deriv(sim, state)
        simpoint%deriv_forces_state(state, :) = simpoint%deriv_forces(:)
    end do

    simpoint%deriv_forces = 0.0
    do state = 1, n
        simpoint%deriv_forces(:) = simpoint%deriv_forces(:) + namd%yg(state)**2.0d0*simpoint%deriv_forces_state(state, :)
    end do

    if (namd%dynam_type .eq. 'aimc') then
        aimc%FM(:) = simpoint%deriv_forces(:)
        aimc%imax = maxloc(abs(namd%yg(1:simpoint%excN)), 1)
        aimc%Fmax(:) = simpoint%deriv_forces_state(aimc%imax, :)
    end if

    call do_sqm_davidson_update(simpoint, 0, vmdqt=vtemp, vgs=vgstemp)
    do state1 = 1, n
        do state2 = state1 + 1, n
            call nacR_analytic_wrap(simpoint, state1, state2, NACR)
            NACR = sim%naesmd%sgn(state1, state2)*NACR
            if ((sim%naesmd%dynam_type .eq. 'aimc' .or. sim%naesmd%dynam_type .eq. 'mf' .or. sim%lprint .ge. 1) .and. restart_flag == 0) then
                if (sim%naesmd%icontw .eq. sim%naesmd%nstepw .or. sim%naesmd%tfemto .eq. 0.0d0) then
                    write (sim%outfile_6, 451) sim%naesmd%tfemto, state1, state2, (NACR(3*j - 2), NACR(3*j - 1), NACR(3*j), j=1, sim%naesmd%natom)
                    call flush (sim%outfile_6)
                end if
            end if
            simpoint%deriv_forces(:) = simpoint%deriv_forces(:) + 2.0*namd%yg(state1)*namd%yg(state2) &
                *cos(namd%yg(state2 + n) - namd%yg(state1 + n))*NACR(:)*(namd%vmdqt(state1) - namd%vmdqt(state2))*feVmdqt
        end do
    end do

    if (namd%dynam_type .eq. 'aimc') then
        aimc%FE(:) = simpoint%deriv_forces(:) - aimc%FM(:)
    end if

    return
451 format(F18.10, I5, I5, 10000(1X, F18.10))
end subroutine deriv_MF
!
