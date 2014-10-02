#include "dprec.fh"
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!THIS SUBROUTINE IS THE NONEQUILIBRIUM STATE SPECIFIC SOLVENT SCF SELECTOR
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine solvent_scf_and_davidson_test();
  use cosmo_C
  implicit none

  if ((solvent_model.eq.3).or.(solvent_model.eq.5)) then
    call calc_cosmo_3();
  else if ((solvent_model.eq.2)) then
    call calc_cosmo_2();
  end if
  return
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!This subroutine calculates the solvent energy <\xi|[V_S(T(\xi)),\xi]>
!It has to be called after all variables have been determined in other
!subroutines
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_excsolven(energy)
        use qm2_davidson_module
        use qmmm_module,only:qm2_struct,qmmm_struct
        use cosmo_C, only : v_solvent_difdens
        implicit none
        _REAL_ :: energy,ddot
        _REAL_ :: tmp(qm2ds%Nb,qm2ds%Nb)
        
        tmp=0.d0; energy=0.d0        
        !Calculate commutator [V_S(T(\xi)),\xi]
        call commutator(qm2ds%v2(1,qmmm_struct%state_of_interest),v_solvent_difdens,qm2ds%Nb,tmp,.false.)

        !Calculate dot product, scale to eV
        energy=ddot(qm2ds%Nb**2,qm2ds%v2(1,qmmm_struct%state_of_interest),1,tmp,1)
        !write(6,*)energy,qm2ds%xi(1),qm2ds%v2(1,qmmm_struct%state_of_interest)
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR [V_s(T+Z),xi] WHICH IS ADDED TO [F(xi),rho_0] in
!DAVIDSON
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_2()
        use qm2_davidson_module
	use communism
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: v_solvent_difdens, cosmo_scf_ftol,cosmo_scf_maxcyc,doZ,potential_type,EF, &
                                linmixparam
        use constants, only : ZERO

        implicit none
        integer INFO;
        integer LWORK;
        _REAL_, DIMENSION(:), allocatable:: WORK,lastxi
        _REAL_ :: OPTIMALSIZE;
	_REAL_ :: xi_1(qm2ds%Nb**2),vsol_temp(qm2ds%Nb,qm2ds%Nb)
        integer verbosity_save,EFsave
        integer i,k,p,h,soi_temp
        _REAL_ e0_0,e0_k,e0_k_1,xi_abs_dif_sum,f0,f1,ddot,energy
        logical calc_Z;

        !Initial Davidson Call (in vacuum) and T+Z calcualtion
	if (qm2ds%verbosity.lt.5) then
                verbosity_save=qm2ds%verbosity;
                qm2ds%verbosity=0; !turn off davidson output
	endif

        call davidson(); !initial call in gas phase
        e0_0 = qm2ds%e0(qmmm_struct%state_of_interest)
        e0_k = e0_0
        allocate(lastxi(qm2ds%Nrpa))
        lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition densities
        soi_temp=qmmm_struct%state_of_interest !Initialize tracking variable

        Calc_Z = doZ

        qmmm_struct%qm_mm_first_call = .false.

        qm2ds%eta(:) = 0.d0 !Clearing

        call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
        call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);

        v_solvent_difdens=0.d0;

        !Calculate Solvent Potential
        if(potential_type.eq.3) then !COSMO
        call VxiM(qm2ds%eta,v_solvent_difdens);
        elseif(potential_type.eq.2) then!ONSAGER
        call rcnfld(v_solvent_difdens,qm2ds%eta,qm2ds%Nb)
        elseif(potential_type.eq.1) then!Straight Correlation
        call Vxi(qm2ds%eta,v_solvent_difdens)
        endif
        v_solvent_difdens=linmixparam*v_solvent_difdens

        !Write header for SCF iterations
        write(6,*)'Start Equilibrium Vertical Excitation Solvent Calculation'
        write(6,*)
        write(6,*)'SCF Step,  Excitation Energy,    DeltaE_sol,      abs(error),error,  COSMO SCF Tolerance '
        write(6,*)'--------------------------------------------------'

        !Begin SCF loop

        do k=1,cosmo_scf_maxcyc
                call davidson(); !Calculate new excited states

                ! Tracking the transition density (checking for crossing)
                f0=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,soi_temp),1))
                do i=1,qm2ds%Mx
                        f1=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,i),1))
                        if(qm2ds%verbosity.gt.4) write(6,*)'Overlaps=',f0,f1
                        if(f0<f1) then
                                write(6,*)'State crossing',qmmm_struct%state_of_interest,' to ',i
                                write(6,*)'New state of interest is',i
                                soi_temp=i
                                f0=f1
                        end if
                end do
                qmmm_struct%state_of_interest=soi_temp
                lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition densities

                !Initialize Variables for Solvent Potential
                vsol_temp=v_solvent_difdens
                call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
                call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);

                !Calculate Solvent Potential
                if(potential_type.eq.3) then !COSMO
                	call VxiM(qm2ds%eta,v_solvent_difdens);
                elseif(potential_type.eq.2) then!ONSAGER
                	call rcnfld(v_solvent_difdens,qm2ds%eta,qm2ds%Nb)
                elseif(potential_type.eq.0) then!Straight Correlation
                	call Vxi(qm2ds%eta,v_solvent_difdens)
                endif
                v_solvent_difdens=linmixparam*v_solvent_difdens+(1.0-linmixparam)*vsol_temp !mixing this and previous solution

                e0_k_1 = e0_k !Save last transition energy
                e0_k = qm2ds%e0(qmmm_struct%state_of_interest)
               !xi_abs_dif_sum=sum(abs(qm2ds%xi)-abs(xi_1))

                write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k ,cosmo_scf_ftol
                if  (abs( e0_k - e0_k_1 )< cosmo_scf_ftol) exit; !Check for convergence
               !if  (abs( xi_abs_dif_sum )< cosmo_scf_ftol) exit; !Check for convergence
                if  (k==cosmo_scf_maxcyc) then
                        write(6,*) 'WARNING! ***SOLVENT SCF REACHED MAXIMUM ITERATIONS WITHOUT CONVERGENCE***'
                endif
        end do
        
        qmmm_struct%qm_mm_first_call = .true.
        if(qm2ds%verbosity.eq.0) qm2ds%verbosity=verbosity_save
        call calc_excsolven(energy)
        if(qm2ds%verbosity>0) then
                write(6,*)
                write(6,*)'Final Results of Equilibrium Vertical Excitation Solvent Calculation '
                write(6,"('    Nonlinear term energy=',e20.10,' eV')") energy
                !write(6,*)' i,   e0(i),    ferr(i),      ftol0'
                write(6,*)'-------------------------------------------------'
                !do k=1,qm2ds%Mx
                !        write(6,112) k,' +++',qm2ds%e0(k),qm2ds%ferr(k),qm2ds%ftol0
                !end do
                !write(6,*)'-------------------------------------------------'
                !write(6,*)
        end if
111     format (i3,' ',g24.16,5('        ',e10.3))
112     format (i3,a,g24.16,2(' ',e8.2))

end subroutine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR EQUILIBRIUM STATE SPECIFIC SOLVENT
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_4(sim_target)
        use qm2_davidson_module
        use communism
        use qmmm_module,only:qm2_struct,qmmm_struct,qmmm_nml;
        use cosmo_C, only: EF,v_solvent_difdens, rhotzpacked_k,cosmo_scf_ftol,cosmo_scf_maxcyc,doZ,potential_type, &
                                linmixparam
        use constants, only : ZERO,AU_TO_EV

        implicit none
        type(simulation_t),target::sim_target
        type(simulation_t),pointer::sim
        integer verbosity_save,verbosity_save2,verbosity_save3,EFsave
	logical verbosity_save4
        integer i,k,p,h,soi_temp
        _REAL_ e0_0,e0_k,e0_k_1,f0,f1,ddot,energy;
        logical calc_Z;
        _REAL_,allocatable:: rhotzpacked_k_temp(:),lastxi(:)

	sim=>sim_target

	call do_sqm_davidson_update(sim)
        e0_0 = (sim%naesmd%Omega(qmmm_struct%state_of_interest)+sim%naesmd%E0)*AU_TO_EV;
        e0_k = e0_0
        calc_Z=doZ
        !Initial Davidson Call (in vacuum) and T+Z calcualtion

        if(qm2ds%verbosity.lt.5) then
        verbosity_save=qm2ds%verbosity;
        verbosity_save2=qmmm_nml%verbosity;
        verbosity_save3=qmmm_nml%printdipole;
        verbosity_save4=qmmm_nml%printcharges;
        qmmm_nml%verbosity=0
        qm2ds%verbosity=0; !turn off davidson output
        qmmm_nml%printdipole=0;
        qmmm_nml%printcharges=.false.;
        endif

        allocate(rhotzpacked_k_temp(qm2ds%nb*(qm2ds%nb+1)/2),lastxi(qm2ds%Nrpa))        

        qmmm_struct%qm_mm_first_call = .false.
        qm2ds%eta(:)=0.d0 !Clearing
	rhotzpacked_k=0.d0

        lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition densities
        soi_temp=qmmm_struct%state_of_interest

        !Get relaxed or unrelaxed difference density
        call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
        call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);
	call packing(qm2ds%nb,qm2ds%eta,rhotzpacked_k, 's')
        
        !Linear mixing
        rhotzpacked_k=linmixparam*rhotzpacked_k 

        !Write header for SCF iterations
        write(6,*)'Start Equilibrium State Specific Solvent Calculation'
        write(6,*)
        write(6,*)'SCF Step,  Excited State Energy,    DeltaE_sol,      abs(error),error,  COSMO SCF Tolerance '
        write(6,*)'--------------------------------------------------'

        !Begin SCF loop
        do k=1,cosmo_scf_maxcyc
                call do_sqm_davidson_update(sim) !to include in the groundstate
                ! Tracking the transition density (checking for crossing)
                f0=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,soi_temp),1))
                do i=1,qm2ds%Mx
                        f1=abs(ddot(qm2ds%Ncis,lastxi,1,qm2ds%v0(1,i),1))
                        if(qm2ds%verbosity.gt.4) write(6,*)'Overlaps=',f0,f1
                        if(f0<f1) then
                                write(6,*)'State crossing',qmmm_struct%state_of_interest,' to ',i
                                write(6,*)'New state of interest is',i
                                soi_temp=i
                                f0=f1
                        end if
                end do
                qmmm_struct%state_of_interest=soi_temp
                lastxi=qm2ds%v0(:,qmmm_struct%state_of_interest) !Store old transition densities

                !Calculate new density for solvent potential
                rhotzpacked_k_temp=rhotzpacked_k
                call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
                call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);
                call packing(qm2ds%nb,qm2ds%eta,rhotzpacked_k, 's')
                rhotzpacked_k=linmixparam*rhotzpacked_k+(1.0-linmixparam)*rhotzpacked_k_temp

                e0_k_1 = e0_k !Save last transition energy
                e0_k = (sim%naesmd%Omega(qmmm_struct%state_of_interest)+sim%naesmd%E0)*AU_TO_EV;

                write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k ,cosmo_scf_ftol

                if  (abs( e0_k - e0_k_1 )< cosmo_scf_ftol) exit; !Check for convergence
                if  (k==cosmo_scf_maxcyc) then
                        write(6,*) 'WARNING! ***SOLVENT SCF REACHED MAXIMUM ITERATIONS WITHOUT CONVERGENCE***'
                endif
        end do


!Printing out found eigenvalues, error and tolerance with solvent
        !qmmm_struct%qm_mm_first_call = .true.
	if(qm2ds%verbosity.eq.0) then
        qm2ds%verbosity=verbosity_save2 !hack
        qmmm_nml%verbosity=verbosity_save2
        qmmm_nml%printdipole=verbosity_save3
        qmmm_nml%printcharges=verbosity_save4
	endif

        !Calculate nonlinear term solvent energy
        v_solvent_difdens(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.d0;!Clearing
        !Save last transition density in AO
        call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,calc_Z);
        call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%eta,qm2ds%tz_scratch);
                !Calculate Solvent Potential
                if(potential_type.eq.3) then !COSMO
                        call VxiM(qm2ds%eta,v_solvent_difdens);
                elseif(potential_type.eq.2) then!ONSAGER
                        call rcnfld(v_solvent_difdens,qm2ds%eta,qm2ds%Nb)
                elseif(potential_type.eq.0) then!Straight Correlation
                        call Vxi(qm2ds%eta,v_solvent_difdens)
                endif

        call calc_excsolven(energy)
        if(qm2ds%verbosity>0) then
                write(6,*)
                write(6,*)'Final Results of Equilibrium State Specific Solvent Calculation'
                write(6,"('    Nonlinear term energy=',e20.10,' eV')") energy
                write(6,*)'-------------------------------------------------'
	        call do_sqm_davidson_update(sim) !to include in the groundstate
        end if

111     format (i3,' ',g24.16,4('        ',e10.3))
112     format (i3,a,g24.16,2(' ',e8.2))

end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!STATE SPECIFIC ROUTINE FOR [V_s(xi),xi] WHICH IS ADDED TO [F(xi),rho_0] in 
!Liouville operator routine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_cosmo_3()
        use qm2_davidson_module
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: v_solvent_xi, cosmo_scf_ftol,cosmo_scf_maxcyc,potential_type;

        implicit none
        integer INFO;
        integer LWORK;
        _REAL_, DIMENSION(:), allocatable:: WORK,xi_k_1
        _REAL_ :: OPTIMALSIZE;
        integer verbosity_save;
        integer k;
        _REAL_ e0_0,e0_k,e0_k_1,f,ddot;
 

        !Initial Davidson Call (in vacuum)
        call davidson();
        call mo2site(qm2ds%v0(1,qmmm_struct%state_of_interest),qm2ds%xi,qm2ds%xi_scratch); !State of Interest to AO Basis
        
        !Initialize Variables for Solvent Potential
        allocate(xi_k_1(qm2_struct%norbs**2))
        v_solvent_xi(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.d0; !Clearing
        xi_k_1=qm2ds%xi
        

        !Calculate Solvent Potential
        if(potential_type.eq.3) then!COSMO
        call VxiM(qm2ds%xi,v_solvent_xi);
        elseif(potential_type.eq.2) then!ONSAGER
        call rcnfld(v_solvent_xi,qm2ds%xi,qm2ds%Nb)
        elseif(potential_type.eq.0) then!Straight Correlation
        call Vxi(qm2ds%eta,v_solvent_xi)
        endif

        !First SCF step
        verbosity_save=qm2ds%verbosity; 
        qm2ds%verbosity=0; !turn off davidson output
        e0_0 = qm2ds%e0(qmmm_struct%state_of_interest); !save vacuum energy
        e0_k_1 = e0_0 !initial energy
        call davidson(); !first davidson call with solvent potential
        e0_k = qm2ds%e0(qmmm_struct%state_of_interest); !save first solventenergy

        !Write header for SCF iterations
        write(6,*)'Start state-specific COSMO SCF'
        write(6,*)
        write(6,*)'SCF Step,  Excitation Energy,    DeltaE_sol,      abs(error),     error,  COSMO SCF Tolerance '
        write(6,*)'--------------------------------------------------'
        
        !Write first SCF iteration results
        write(6,111)1, e0_k ,e0_k-e0_0,abs( e0_k - e0_k_1 ), e0_k_1-e0_k , cosmo_scf_ftol

        !Begin SCF loop
        do k=2,cosmo_scf_maxcyc
                if  (abs( e0_k - e0_k_1 )< cosmo_scf_ftol) exit; !Check for convergence
 
                !Initialize Variables for Solvent Potential
                v_solvent_xi(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.d0; !Clearing
                 
                !Save last transition density in AO Basis
                call mo2site(qm2ds%v0(1,qmmm_struct%state_of_interest),qm2ds%xi,qm2ds%xi_scratch); 
                
		!XL-BOXMD
                ! Checking if transition density changes sign and correcting for this
                f=ddot(qm2ds%Ncis,xi_k_1,1,qm2ds%xi,1)
                if(f<0.d0) then
                qm2ds%xi=(-1.d0)*qm2ds%xi;
                endif
                
                xi_k_1=qm2ds%xi!Save last transition density

                !Calculate Solvent Potential
                if(potential_type.eq.3) then!COSMO
                call VxiM(qm2ds%xi,v_solvent_xi);
                elseif(potential_type.eq.2) then!ONSAGER
                call rcnfld(v_solvent_xi,qm2ds%xi,qm2ds%Nb)
                elseif(potential_type.eq.0) then!Straight Correlation
                call Vxi(qm2ds%xi,v_solvent_xi)
                endif
                
                e0_k_1 = e0_k !Save last transition energy
                call davidson(); !Calculate new excited states
                e0_k = qm2ds%e0(qmmm_struct%state_of_interest);
                write(6,111)k, e0_k ,e0_k-e0_0,abs(e0_k-e0_k_1), e0_k_1-e0_k , cosmo_scf_ftol

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

111     format (i3,' ',g24.16,4('        ',e10.3))
112     format (i3,a,g24.16,2(' ',e8.2))

end subroutine
