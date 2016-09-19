#include "dprec.fh"
!
!********************************************************************
!
!New Liouville functions
!
!by Kirill A Velizhanin (kirill@lanl.gov)
!
!********************************************************************
!
!
!********************************************************************
!
subroutine dav_wrap(sim_target)
    use communism !for sim hack
    use qmmm_module,only:qmmm_struct
    use qm2_davidson_module
    use cosmo_C, only: solvent_model

    implicit none

    _REAL_ f,ddot
    integer i
    type(simulation_t),target::sim_target
    type(simulation_t),pointer::sim

    sim=>sim_target
    ! Finding excited states

    !Vacuum, Linear Response Solvent have the single Davidson routine, Nonequilibrium State Specific
    !has iterative Davidson Wrapper, Equilibrium State Specific routine has scf and Davidson wrapper above this subroutine.
    if ((solvent_model.lt.2).or.(solvent_model.gt.3)) then
        call davidson();
    elseif ((solvent_model.eq.2).or.(solvent_model.eq.3)) then
        call solvent_scf_and_davidson_test(sim);
    end if

    ! Total energy of the ground state
    qm2ds%Eground=qmmm_struct%elec_eng+qmmm_struct%enuclr_qmmm &
        +qmmm_struct%enuclr_qmqm
    ! Total energy of required excited state
    qm2ds%Etot(:)=qm2ds%Eground+qm2ds%e0(:)
    ! cml-test Ereq may be able to be removed after implementation of struct_opt_state
    ! Total energy of the required state
    if(qm2ds%Mx>0) then
        qm2ds%Ereq=qm2ds%Etot(qm2ds%Mx)
    else
        qm2ds%Ereq=qm2ds%Eground
    end if

    ! Checking if some transition densities are suddenly changing signs
    ! correcting for this
    if(qm2ds%mdflag==2) then ! relevant only if molecular dynamics
        do i=1,qm2ds%Mx
            f=ddot(qm2ds%Ncis,qm2ds%v0_old(1,i),1,qm2ds%v0(1,qm2ds%kx(i)),1)

            if(f<0.d0) then
                qm2ds%v0(:,i)=-qm2ds%v0(:,i)
            end if
        end do
    end if

    qm2ds%v0_old(:,:)=qm2ds%v0(:,:)

    ! Output
    if (qm2ds%verbosity>0) then
        call outDavidson()
    endif
        
    if (qm2ds%calcxdens) then
        write(6,*)'Calculating cross densities and printing to CEO.out'
        call polarizab()
    endif

    qm2ds%has_been_run = .TRUE.

    return
endsubroutine dav_wrap
!
!********************************************************************
!
!Output of various results of Davidson
!
!********************************************************************
!
subroutine outDavidson()
    use qm2_davidson_module
    use qmmm_module, only : qmmm_nml
    use constants, only : ONE_AU
    implicit none

    integer i
    _REAL_ mu(3, qm2ds%Mx), alpha(3), ft

    ! Evaluation of transition dipole moments

    write(6,*) 'Frequencies (eV) and Oscillator strengths (unitless)'
    write(6,"(8x,'Omega',12x,'fx',14x,'fy',14x,'fz',10x,'ftotal')")

    call trans_dipole(mu, alpha) ! cml-test
    do i=1,qm2ds%Mx
        ft = (2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*(mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
        write(6,"(i4,5g25.15)") i,qm2ds%e0(i), &
            (2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(1,i)**2, &
            (2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(2,i)**2, &
            (2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(3,i)**2, ft
    end do

    write(6,*)
    write(6,*) 'Frequencies (eV) and Transition Dipole Moments (AU)'
    write(6,"(8x,'Omega',12x,'fx',14x,'fy',14x,'fz',10x,'ftotal')")
    do i=1,qm2ds%Mx
        ft = (mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
        write(6,"(i4,5g25.15)") i,qm2ds%e0(i), &
            mu(1,i),mu(2,i),mu(3,i), ft
    end do


    write(6,*)
    write(6,*)'Total energy of the ground state (eV,AU)'
    write(6,*) 0,qm2ds%Eground, qm2ds%Eground/ONE_AU

    write(6,*) 'Total energies of excited states (eV,AU)'
    do i=1,qm2ds%Mx
        write(6,*) i,qm2ds%Etot(i),qm2ds%Etot(i)/ONE_AU
    end do

    !Printing of transition densities to files
    !Print normal modes, etc. to file ! JAB
    if(qmmm_nml%printtd.gt.0) then
        open(76,file='normalmodes-ao.out')
        write(6,*) 'Printing Normal Modes in ao rep to file'
        call printNM_ao(76)
        close(76)
        if(qmmm_nml%printtd.gt.1) then
            open(76,file='normalmodes-mo.out')
            write(6,*) 'Printing Normal Modes in mo rep to file'
            call printNM_mo(76)
            close(76)
            if (qmmm_nml%printtd.gt.2) then
                !open(76,file='2eri.out')
                !write(6,*) 'Printing 2e Repulsion integrals to file'
                !call print2eri()
                !close(76)
                open(76,file='normalmodes-cfit.out')
                write(6,*) 'Printing Charges for Normal Modes to file'
                call printCfitNM(76)
                close(76)
            endif
        endif
    elseif(qmmm_nml%printtd.eq.-1) then !print all TD's to binary
        write(6,*) 'Printing Normal Modes to modes.b file'
        call printNM_binary(76,77)
    endif
    return
endsubroutine outDavidson
!
