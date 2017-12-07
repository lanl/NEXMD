#include "dprec.fh"
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!THIS SUBROUTINE IS THE NONEQUILIBRIUM STATE SPECIFIC SOLVENT SCF SELECTOR
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine solvent_scf_and_davidson_test(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);
    use qm2_davidson_module
    use cosmo_C, only : cosmo_C_structure
    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_module,only:qm2_structure

    implicit none
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct

    if (cosmo_c_struct%solvent_model.eq.3) then
        call calc_cosmo_3(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);
    else if (cosmo_c_struct%solvent_model.eq.2) then
        call calc_cosmo_2(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);
    end if
    return
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!This subroutine calculates the solvent energy <\xi|[V_S(T(\xi)),\xi]>
!It has to be called after all variables have been determined in other
!subroutines
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_excsolven(cosmo_c_struct, qm2ds, qmmm_struct,energy)
    use qm2_davidson_module
    use cosmo_C, only : cosmo_C_structure
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
     type(qm2_davidson_structure_type), intent(inout) :: qm2ds
     type(qmmm_struct_type), intent(inout) :: qmmm_struct

    _REAL_ :: energy,ddot
    _REAL_ :: tmp(qm2ds%Nb,qm2ds%Nb)
        
    tmp=0.d0; energy=0.d0
    !Calculate commutator [V_S(T(\xi)),\xi]
    call commutator(qm2ds%v2(1,qmmm_struct%state_of_interest),cosmo_c_struct%v_solvent_difdens,qm2ds%Nb,tmp,.false.)

    !Calculate dot product, scale to eV
    energy=ddot(qm2ds%Nb**2,qm2ds%v2(1,qmmm_struct%state_of_interest),1,tmp,1)
        !write(6,*)energy,qm2ds%xi(1),qm2ds%v2(1,qmmm_struct%state_of_interest)
end subroutine

!This subroutine calculates the solvent energy Tr(PV_S(T_k))
!It has to be called after all variables have been determined in other
!subroutines
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_gssolven(cosmo_c_struct,qm2_struct, qm2ds,qmmm_struct,energy)
    use qm2_davidson_module
    use qmmm_module,only:qm2_structure
    use cosmo_C, only : cosmo_C_structure ! rhotzpacked_k
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    _REAL_ :: energy,ddot

    energy=0.d0
    qm2ds%tz_scratch=0.d0
    call addfck(cosmo_c_struct,qmmm_struct, qm2ds%tz_scratch(1),qm2_struct%den_matrix);
    call unpacking(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%tz_scratch(qm2ds%Nb**2+1),'s')
    call unpacking(qm2ds%Nb,cosmo_c_struct%rhotzpacked_k,qm2ds%tz_scratch(1),'s')

    !Calculate dot product, scale to eV
    energy=ddot(qm2ds%Nb**2,qm2ds%tz_scratch(qm2ds%Nb**2+1),1,qm2ds%tz_scratch(1),1)
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Calculate the potential operator [[xi^T,V(X)],xi] for Z-vector and/or
! GS calculations
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_xicommutator(cosmo_c_struct,qm2ds,V_potential)
    use qm2_davidson_module
    use cosmo_C, only : cosmo_C_structure !xi_k
    implicit none
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    _REAL_ ::V_potential(qm2ds%Nb,qm2ds%Nb),tmp(qm2ds%Nb,qm2ds%Nb)
    tmp=0.d0;
    call commutator(cosmo_c_struct%xi_k,V_potential,qm2ds%Nb,tmp,.true.)!inner commutator
    call commutator(tmp,cosmo_c_struct%xi_k,qm2ds%Nb,V_potential,.false.) !second commutator with transpose
end subroutine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR [V_s(T+Z),xi] WHICH IS ADDED TO [F(xi),rho_0] in
!DAVIDSON
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_2(cosmo_c_struct,qm2_struct, qm2ds,qmmm_struct)
    use qm2_davidson_module
    use communism
    use cosmo_C, only: cosmo_C_structure !cosmo_c_struct%v_solvent_difdens, cosmo_c_struct%cosmo_scf_ftol,cosmo_c_struct%cosmo_scf_maxcyc,cosmo_c_struct%doZ,cosmo_c_struct%potential_type, &
        !cosmo_c_struct%linmixparam
    use constants, only : ZERO

    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_module,only:qm2_structure

    implicit none
     type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
     type(qm2_structure),intent(inout) :: qm2_struct
     type(qmmm_struct_type), intent(inout) :: qmmm_struct
     type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    _REAL_ :: lastxi(qm2ds%Nrpa)
    _REAL_ :: vsol_temp(qm2ds%Nb,qm2ds%Nb)
    integer verbosity_save
    integer i,k,soi_temp
    _REAL_ e0_0,e0_k,e0_k_1,f0,f1,ddot,energy
    logical calc_Z;

    !Initial Davidson Call (in vacuum) and T+Z calcualtion
    if (qm2ds%verbosity.lt.5) then
        verbosity_save=qm2ds%verbosity;
        qm2ds%verbosity=0;    !turn off davidson output
    endif

    !Write header for SCF iterations
    write(6,*)
    write(6,*)'Start Equilibrium Vertical Excitation Solvent Calculation'
    write(6,*)
    write(6,*)'SCF Step,  Excitation Energy,    DeltaE_sol, abs(error), error,  COSMO SCF Tolerance '
    write(6,*)'--------------------------------------------------'
        
    if(.not.qmmm_struct%qm_mm_first_call) then !From previous time step in dynamics
        lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Get old transition density
    endif
    soi_temp=qmmm_struct%state_of_interest

    call davidson(cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct);    !initial call in gas phase

    ! Tracking the transition density (checking for crossing) from last time
    ! step during dynamics
    if(.not.qmmm_struct%qm_mm_first_call) then
        f0=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,soi_temp),1))
        do i=1,qm2ds%Mx
            f1=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,i),1))
            if(qm2ds%verbosity.gt.4) write(6,*)'Overlaps=',f0,f1
            if(f0<f1) then
                write(6,*)'State crossing',qmmm_struct%state_of_interest,' to ',i
                write(6,*)'New state of interest is',i
                soi_temp=i
                f0=f1
            endif
        enddo
        qmmm_struct%state_of_interest=soi_temp
    endif
    lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition density

    if(qm2ds%verbosity.eq.5) call outDavidson(qm2_struct,qm2ds,qmmm_struct)

    !Initialize some variables
    e0_0 = qm2ds%e0(qmmm_struct%state_of_interest)
    e0_k = e0_0
    Calc_Z = cosmo_c_struct%doZ !Input to calculate Z or not
    qmmm_struct%qm_mm_first_call = .false. !After first iteration it's false
    qm2ds%eta(:) = 0.d0 !Clearing

    call calc_rhotz(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct,qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
    call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);

    cosmo_c_struct%v_solvent_difdens=0.d0;

    !Calculate Solvent Potential
    if(cosmo_c_struct%potential_type.eq.3) then !COSMO
        call VxiM(cosmo_c_struct,qm2ds,qm2ds%eta,cosmo_c_struct%v_solvent_difdens);
    elseif(cosmo_c_struct%potential_type.eq.2) then!ONSAGER
        call rcnfld(cosmo_c_struct,qm2_struct,qmmm_struct,cosmo_c_struct%v_solvent_difdens,qm2ds%eta,qm2ds%Nb)
    elseif(cosmo_c_struct%potential_type.eq.1) then!Straight Correlation
        call Vxi(qm2_struct,qm2ds,qmmm_struct,qm2ds%eta,cosmo_c_struct%v_solvent_difdens)
    endif
    cosmo_c_struct%v_solvent_difdens=cosmo_c_struct%linmixparam*cosmo_c_struct%v_solvent_difdens
    !Begin SCF loop
    do k=1,cosmo_c_struct%cosmo_scf_maxcyc
        call davidson(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);        !Calculate new excited states
        if(qm2ds%verbosity.eq.5) call outDavidson(qm2_struct,qm2ds,qmmm_struct)

        ! Tracking the transition density (checking for crossing)
        f0=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,soi_temp),1))
        do i=1,qm2ds%Mx
            f1=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,i),1))
            if(qm2ds%verbosity.gt.4) write(6,*)'Overlaps=',f0,f1
            if(f0<f1) then
                write(6,*)'State crossing due to solvent effect ',qmmm_struct%state_of_interest,' to ',i
                write(6,*)'New state of interest is',i
                soi_temp=i
                f0=f1
            end if
        end do
        !i=qmmm_struct%state_of_interest !for switching back at the end
        qmmm_struct%state_of_interest=soi_temp
        lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition densities

        !Initialize Variables for Solvent Potential
        vsol_temp=cosmo_c_struct%v_solvent_difdens
        call calc_rhotz(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct,qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
        call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);
        !Calculate Solvent Potential
        cosmo_c_struct%v_solvent_difdens=0.d0
        if(cosmo_c_struct%potential_type.eq.3) then !COSMO
            call VxiM(cosmo_c_struct,qm2ds,qm2ds%eta,cosmo_c_struct%v_solvent_difdens);
        elseif(cosmo_c_struct%potential_type.eq.2) then!ONSAGER
            call rcnfld(cosmo_c_struct,qm2_struct,qmmm_struct,cosmo_c_struct%v_solvent_difdens,qm2ds%eta,qm2ds%Nb)
        elseif(cosmo_c_struct%potential_type.eq.0) then!Straight Correlation
            call Vxi(qm2_struct,qm2ds,qmmm_struct,qm2ds%eta,cosmo_c_struct%v_solvent_difdens)
        endif
        cosmo_c_struct%v_solvent_difdens=cosmo_c_struct%linmixparam*&
		cosmo_c_struct%v_solvent_difdens+(1.0-cosmo_c_struct%linmixparam)*vsol_temp !mixing this and previous solution

        e0_k_1 = e0_k !Save last transition energy
        e0_k = qm2ds%e0(qmmm_struct%state_of_interest)
        !xi_abs_dif_sum=sum(abs(qm2ds%xi)-abs(xi_1))

        write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k ,cosmo_c_struct%cosmo_scf_ftol
        if  (abs( e0_k - e0_k_1 )< cosmo_c_struct%cosmo_scf_ftol) exit;        !Check for convergence
        !if  (abs( xi_abs_dif_sum )< cosmo_c_struct%cosmo_scf_ftol) exit; !Check for convergence
        if  (k==cosmo_c_struct%cosmo_scf_maxcyc) then
            write(6,*) '***SOLVENT SCF REACHED MAXIMUM ITERATIONS WITHOUT CONVERGENCE***'
            stop
        endif
    end do
        
    !qmmm_struct%qm_mm_first_call = .true.
    if(qm2ds%verbosity.eq.0) qm2ds%verbosity=verbosity_save
    call calc_excsolven(cosmo_c_struct,qm2ds,qmmm_struct,energy)
    if(qm2ds%verbosity>0) then
        write(6,*)
        write(6,*)'Final Results of Equilibrium Vertical Excitation Solvent Calculation '
        write(6,"('  ES nonlinear term energy=',e20.10,' eV')") energy
        !write(6,*)' i,   e0(i),    ferr(i),      ftol0'
        write(6,*)'-------------------------------------------------'
            !do k=1,qm2ds%Mx
            !        write(6,112) k,' +++',qm2ds%e0(k),qm2ds%ferr(k),qm2ds%ftol0
            !end do
            !write(6,*)'-------------------------------------------------'
            !write(6,*)
    end if
111 format (i3,' ',g24.16,5('        ',e10.3))
112 format (i3,a,g24.16,2(' ',e8.2))

end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR EQUILIBRIUM STATE SPECIFIC SOLVENT
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_4(sim_target)
    use qm2_davidson_module
    use communism
    use qmmm_module,only:qmmm_nml;
    !use cosmo_C, only : sim%cosmo%v_solvent_difdens, sim%cosmo%rhotzpacked_k,sim%cosmo%cosmo_scf_ftol,sim%cosmo%cosmo_scf_maxcyc,sim%cosmo%doZ,sim%cosmo%potential_type, &
    !    sim%cosmo%linmixparam,sim%cosmo%xi_k,sim%cosmo%v_solvent_difdens,sim%cosmo%solvent_model
    use constants, only : ZERO,AU_TO_EV

    implicit none
    type(simulation_t),target::sim_target
    type(simulation_t),pointer::sim
    integer verbosity_save,verbosity_save2,verbosity_save3
    logical verbosity_save4
    integer i,k,soi_temp
    _REAL_ e0_0,e0_k,e0_k_1,f0,f1,ddot,energy,gsenergy;
    logical calc_Z;
    _REAL_,allocatable:: rhotzpacked_k_temp(:),lastxi(:), &
        xi_k_temp(:),vsol_temp(:,:)

    sim=>sim_target

    call do_sqm_davidson_update(sim)
    e0_0 = (sim%naesmd%Omega(sim%qmmm%state_of_interest)+sim%naesmd%E0)*AU_TO_EV;
    e0_k = e0_0
    calc_Z=sim%cosmo%doZ
    !Initial Davidson Call (in vacuum) and T+Z calcualtion

    if(sim%dav%verbosity.lt.5) then
        verbosity_save=sim%dav%verbosity;
        verbosity_save2=qmmm_nml%verbosity;
        verbosity_save3=qmmm_nml%printdipole;
        verbosity_save4=qmmm_nml%printcharges;
        qmmm_nml%verbosity=0
        sim%dav%verbosity=0;        !turn off davidson output
        qmmm_nml%printdipole=0;
        qmmm_nml%printcharges=.false.;
    endif

    allocate(rhotzpacked_k_temp(sim%dav%nb*(sim%dav%nb+1)/2),lastxi(sim%dav%Nrpa),&
        xi_k_temp(sim%dav%nb**2))

    if(sim%cosmo%solvent_model.eq.5) allocate(vsol_temp(sim%dav%nb,sim%dav%nb))

    sim%qmmm%qm_mm_first_call = .false.
    sim%dav%eta(:)=0.d0 !Clearing
    sim%cosmo%rhotzpacked_k=0.d0

    lastxi=sim%dav%v0(:,sim%qmmm%state_of_interest) !Store old transition densities
    soi_temp=sim%qmmm%state_of_interest

    !Get relaxed or unrelaxed difference density
    call calc_rhotz(sim%cosmo,sim%qm2,sim%dav,sim%qmmm,sim%qmmm%state_of_interest,sim%dav%rhoTZ,calc_Z);
    call mo2sitef(sim%dav%Nb,sim%dav%vhf,sim%dav%rhoTZ,sim%dav%eta,sim%dav%tz_scratch);
    call packing(sim%dav%nb,sim%dav%eta,sim%cosmo%rhotzpacked_k, 's')

    sim%dav%eta(:) = 0.d0 !Clearing
        
    !Term for variational formulation
    if(sim%cosmo%solvent_model.eq.5) then !need solvent potential for excitation
        sim%cosmo%v_solvent_difdens=0.d0;
        !Calculate Solvent Potential
        if(sim%cosmo%potential_type.eq.3) then !COSMO
            call VxiM(sim%cosmo,sim%dav,sim%dav%eta,sim%cosmo%v_solvent_difdens);
        elseif(sim%cosmo%potential_type.eq.2) then!ONSAGER
            call rcnfld(sim%cosmo,sim%qm2,sim%qmmm,sim%cosmo%v_solvent_difdens,sim%dav%eta,sim%dav%Nb)
        elseif(sim%cosmo%potential_type.eq.1) then!Straight Correlation
            call Vxi(sim%qm2,sim%dav,sim%qmmm,sim%dav%eta,sim%cosmo%v_solvent_difdens)
        endif
        sim%cosmo%v_solvent_difdens=sim%cosmo%linmixparam*sim%cosmo%v_solvent_difdens
    endif

    !Get transition density k in AO
    call getmodef(2*sim%dav%Np*sim%dav%Nh,sim%dav%Mx,sim%dav%Np,sim%dav%Nh, &
        sim%qmmm%state_of_interest,sim%dav%v0,sim%dav%eta_tz)
    call mo2sitef(sim%dav%Nb,sim%dav%vhf,sim%dav%eta_tz,sim%cosmo%xi_k, &
        sim%dav%tz_scratch(sim%dav%Nb**2+1))
         
    !Linear mixing
    sim%cosmo%rhotzpacked_k=sim%cosmo%linmixparam*sim%cosmo%rhotzpacked_k
    sim%cosmo%xi_k=sim%cosmo%linmixparam*sim%cosmo%xi_k

    !Write header for SCF iterations
    write(6,*)
    write(6,*)'Start Equilibrium State Specific Solvent Calculation'
    if(sim%cosmo%solvent_model.eq.5) write(6,*)'******Using Varational Formulation******'
    write(6,*)
    write(6,*)'SCF Step,  Excited State Energy,    DeltaE_sol,      abs(error),error,  COSMO SCF Tolerance '
    write(6,*)'--------------------------------------------------'

    !Begin SCF loop
    do k=1,sim%cosmo%cosmo_scf_maxcyc
        call do_sqm_davidson_update(sim) !to include in the groundstate
        ! Tracking the transition density (checking for crossing)
        f0=abs(ddot(sim%dav%Ncis,lastxi,1,sim%dav%v0(1,soi_temp),1))
        do i=1,sim%dav%Mx
            f1=abs(ddot(sim%dav%Ncis,lastxi,1,sim%dav%v0(1,i),1))
            if(sim%dav%verbosity.gt.4) write(6,*)'Overlaps=',f0,f1
            if(f0<f1) then
                if(abs(f0-f1)>0.95) then
                    write(6,*)'State crossing due to solvent effect',sim%qmmm%state_of_interest,' to ',i
                    write(6,*)'New state of interest is',i
                    soi_temp=i
                    f0=f1
                else
                    write(6,*)'WARNING: STATES CANNOT BE FOLLOWED,&
                                               & TRY CHANGING THE LINEAR MIXING PARAMETER'
                endif
            end if
        end do
        sim%qmmm%state_of_interest=soi_temp
        lastxi=sim%dav%v0(:,sim%qmmm%state_of_interest) !Store old transition densities

        !Calculate new density for solvent potential
        rhotzpacked_k_temp=sim%cosmo%rhotzpacked_k
        call calc_rhotz(sim%cosmo, sim%qm2,sim%dav,sim%qmmm,sim%qmmm%state_of_interest,sim%dav%rhoTZ,calc_Z);
        call mo2sitef(sim%dav%Nb,sim%dav%vhf,sim%dav%rhoTZ,sim%dav%eta,sim%dav%tz_scratch);
        call packing(sim%dav%nb,sim%dav%eta,sim%cosmo%rhotzpacked_k, 's')
        sim%cosmo%rhotzpacked_k=sim%cosmo%linmixparam*sim%cosmo%rhotzpacked_k+(1.0-sim%cosmo%linmixparam)*rhotzpacked_k_temp
                
        !Potential excitation calculation in  variational formulation
        if(sim%cosmo%solvent_model.eq.5) then
            vsol_temp=sim%cosmo%v_solvent_difdens
            !Calculate Solvent Potential
            sim%cosmo%v_solvent_difdens=0.d0
            if(sim%cosmo%potential_type.eq.3) then !COSMO
                call VxiM(sim%cosmo,sim%dav,sim%dav%eta,sim%cosmo%v_solvent_difdens);
            elseif(sim%cosmo%potential_type.eq.2) then!ONSAGER
                call rcnfld(sim%cosmo,sim%qm2,sim%qmmm,sim%cosmo%v_solvent_difdens,sim%dav%eta,sim%dav%Nb)
            elseif(sim%cosmo%potential_type.eq.0) then!Straight Correlation
                call Vxi(sim%qm2,sim%dav,sim%qmmm,sim%dav%eta,sim%cosmo%v_solvent_difdens)
            endif
            sim%cosmo%v_solvent_difdens=sim%cosmo%linmixparam*sim%cosmo%v_solvent_difdens &
                +(1.0-sim%cosmo%linmixparam)*vsol_temp
        endif

        !get transition density k in AO
        xi_k_temp=sim%cosmo%xi_k
        call getmodef(2*sim%dav%Np*sim%dav%Nh,sim%dav%Mx,sim%dav%Np,sim%dav%Nh, &
            sim%qmmm%state_of_interest,sim%dav%v0,sim%dav%eta_tz)
        call mo2sitef(sim%dav%Nb,sim%dav%vhf,sim%dav%eta_tz,sim%cosmo%xi_k, &
            sim%dav%tz_scratch(sim%dav%Nb**2+1))
        sim%cosmo%xi_k=sim%cosmo%linmixparam*sim%cosmo%xi_k+(1.0-sim%cosmo%linmixparam)*xi_k_temp
                 
        e0_k_1 = e0_k !Save last transition energy
        e0_k = (sim%naesmd%Omega(sim%qmmm%state_of_interest)+sim%naesmd%E0)*AU_TO_EV;

        write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k ,sim%cosmo%cosmo_scf_ftol

        if  (abs( e0_k - e0_k_1 )< sim%cosmo%cosmo_scf_ftol) exit;        !Check for convergence
        if  (k==sim%cosmo%cosmo_scf_maxcyc) then
            write(6,*) 'WARNING! ***SOLVENT SCF REACHED MAXIMUM ITERATIONS WITHOUT CONVERGENCE***'
        endif
    end do


    !Printing out found eigenvalues, error and tolerance with solvent
    !qmmm_struct%qm_mm_first_call = .true.
    if(sim%dav%verbosity.eq.0) then
        sim%dav%verbosity=verbosity_save2 !hack
        qmmm_nml%verbosity=verbosity_save2
        qmmm_nml%printdipole=verbosity_save3
        qmmm_nml%printcharges=verbosity_save4
    endif

    !Calculate nonlinear term solvent energy
    sim%cosmo%v_solvent_difdens(1:sim%qm2%norbs,1:sim%qm2%norbs)=0.d0;    !Clearing
    !Save last transition density in AO
    call calc_rhotz(sim%cosmo, sim%qm2,sim%dav,sim%qmmm,sim%qmmm%state_of_interest,sim%dav%rhoTZ,calc_Z);
    call mo2sitef(sim%dav%Nb,sim%dav%vhf,sim%dav%rhoTZ,sim%dav%eta,sim%dav%tz_scratch);
    !Calculate Solvent Potential
    if(sim%cosmo%potential_type.eq.3) then !COSMO
        call VxiM(sim%cosmo,sim%dav,sim%dav%eta,sim%cosmo%v_solvent_difdens);
    elseif(sim%cosmo%potential_type.eq.2) then!ONSAGER
        call rcnfld(sim%cosmo,sim%qm2,sim%qmmm,sim%cosmo%v_solvent_difdens,sim%dav%eta,sim%dav%Nb)
    elseif(sim%cosmo%potential_type.eq.0) then!Straight Correlation
        call Vxi(sim%qm2,sim%dav,sim%qmmm,sim%dav%eta,sim%cosmo%v_solvent_difdens)
    endif

    call calc_excsolven(sim%cosmo,sim%dav,sim%qmmm,energy)
    call calc_gssolven(sim%cosmo,sim%qm2,sim%dav,sim%qmmm,gsenergy)
    if(sim%dav%verbosity>0) then
        write(6,*)
        write(6,*)'Final Results of Equilibrium State Specific Solvent Calculation'
        write(6,"('    GS nonlinear term energy=',e20.10,' eV')") gsenergy
        write(6,"('    ES nonlinear term energy=',e20.10,' eV')") energy
        write(6,*)'-------------------------------------------------'
        call do_sqm_davidson_update(sim) !to include in the groundstate
    end if

111 format (i3,' ',g24.16,4('        ',e10.3))

end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR [V_s(xi),xi] WHICH IS ADDED TO [F(xi),rho_0] in 
!Liouville operator routine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_3(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct)
    use qm2_davidson_module
    use qmmm_module,only:qm2_structure;
    use cosmo_C, only: cosmo_C_structure !cosmo_c_struct%v_solvent_xi, cosmo_c_struct%cosmo_scf_ftol,cosmo_c_struct%cosmo_scf_maxcyc,cosmo_c_struct%potential_type;
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
   _REAL_, DIMENSION(:), allocatable:: xi_k_1
    integer verbosity_save;
    integer k;
    _REAL_ e0_0,e0_k,e0_k_1,f,ddot;
 

    !Initial Davidson Call (in vacuum)
    call davidson(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);
    call mo2site(qm2ds,qm2ds%v0(1,qmmm_struct%state_of_interest),qm2ds%xi,qm2ds%xi_scratch);    !State of Interest to AO Basis
        
    !Initialize Variables for Solvent Potential
    allocate(xi_k_1(qm2_struct%norbs**2))
    cosmo_c_struct%v_solvent_xi(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.d0;    !Clearing
    xi_k_1=qm2ds%xi
        

    !Calculate Solvent Potential
    if(cosmo_c_struct%potential_type.eq.3) then!COSMO
        call VxiM(cosmo_c_struct,qm2ds,qm2ds%xi,cosmo_c_struct%v_solvent_xi);
    elseif(cosmo_c_struct%potential_type.eq.2) then!ONSAGER
        call rcnfld(cosmo_c_struct,qm2_struct,qmmm_struct,cosmo_c_struct%v_solvent_xi,qm2ds%xi,qm2ds%Nb)
    elseif(cosmo_c_struct%potential_type.eq.0) then!Straight Correlation
        call Vxi(qm2_struct,qm2ds,qmmm_struct,qm2ds%eta,cosmo_c_struct%v_solvent_xi)
    endif

    !First SCF step
    verbosity_save=qm2ds%verbosity;
    qm2ds%verbosity=0;    !turn off davidson output
    e0_0 = qm2ds%e0(qmmm_struct%state_of_interest);    !save vacuum energy
    e0_k_1 = e0_0 !initial energy
    call davidson(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);    !first davidson call with solvent potential
    e0_k = qm2ds%e0(qmmm_struct%state_of_interest);    !save first solventenergy

    !Write header for SCF iterations
    write(6,*)'Start state-specific COSMO SCF'
    write(6,*)
    write(6,*)'SCF Step,  Excitation Energy,    DeltaE_sol,      abs(error),     error,  COSMO SCF Tolerance '
    write(6,*)'--------------------------------------------------'
        
    !Write first SCF iteration results
    write(6,111)1, e0_k ,e0_k-e0_0,abs( e0_k - e0_k_1 ), e0_k_1-e0_k , cosmo_c_struct%cosmo_scf_ftol

    !Begin SCF loop
    do k=2,cosmo_c_struct%cosmo_scf_maxcyc
        if  (abs( e0_k - e0_k_1 )< cosmo_c_struct%cosmo_scf_ftol) exit;        !Check for convergence
 
        !Initialize Variables for Solvent Potential
        cosmo_c_struct%v_solvent_xi(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.d0;        !Clearing
                 
        !Save last transition density in AO Basis
        call mo2site(qm2ds,qm2ds%v0(1,qmmm_struct%state_of_interest),qm2ds%xi,qm2ds%xi_scratch);
                
        ! Checking if transition density changes sign and correcting for this
        f=ddot(qm2ds%Ncis,xi_k_1,1,qm2ds%xi,1)
        if(f<0.d0) then
            qm2ds%xi=(-1.d0)*qm2ds%xi;
        endif
                
        xi_k_1=qm2ds%xi!Save last transition density

        !Calculate Solvent Potential
        if(cosmo_c_struct%potential_type.eq.3) then!COSMO
            call VxiM(cosmo_c_struct,qm2ds,qm2ds%xi,cosmo_c_struct%v_solvent_xi);
        elseif(cosmo_c_struct%potential_type.eq.2) then!ONSAGER
            call rcnfld(cosmo_c_struct,qm2_struct,qmmm_struct,cosmo_c_struct%v_solvent_xi,qm2ds%xi,qm2ds%Nb)
        elseif(cosmo_c_struct%potential_type.eq.0) then!Straight Correlation
            call Vxi(qm2_struct,qm2ds,qmmm_struct,qm2ds%xi,cosmo_c_struct%v_solvent_xi)
        endif
                
        e0_k_1 = e0_k !Save last transition energy
        call davidson(cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct);        !Calculate new excited states
        e0_k = qm2ds%e0(qmmm_struct%state_of_interest);
        write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k , cosmo_c_struct%cosmo_scf_ftol

    end do
    qmmm_struct%qm_mm_first_call = .true.
    qm2ds%verbosity=verbosity_save
    deallocate(xi_k_1)

    !Printing out found eigenvalues, error and tolerance with solvent

    if(qm2ds%verbosity>0) then
        write(6,*)
        write(6,*)'    Final Results of last Davidson procedure  '
        write(6,*)
        write(6,*)' i,   e0(i),    ferr(i),      ftol0'
        write(6,*)'-------------------------------------------------'
        do k=1,qm2ds%Mx
            write(6,112) k,' +++ ',qm2ds%e0(k),qm2ds%ferr(k),qm2ds%ftol0
        end do
        write(6,*)'-------------------------------------------------'
        write(6,*)
    end if

111 format (i3,' ',g24.16,4('        ',e10.3))
112 format (i3,a,g24.16,2(' ',e8.2))

end subroutine

