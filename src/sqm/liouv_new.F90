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
subroutine dav_wrap()
	use qmmm_module,only:qmmm_struct
	use qm2_davidson_module
	use cosmo_C, only: solvent_model,ceps

	implicit none

	_REAL_ f,ddot
	integer i

	! Finding excited states

	!Vacuum, Linear Response Solvent have the single Davidson routine, Nonequilibrium State Specific
	!has iterative Davidson Wrapper, Equilibrium State Specific routine has scf and Davidson wrapper above this subroutine.
	if ((solvent_model.lt.2).or.(solvent_model.gt.3)) then
		call davidson();
	elseif ((solvent_model.eq.2).or.(solvent_model.eq.3)) then
		call solvent_scf_and_davidson_test();
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
	use constants, only : ONE_AU
	implicit none

	integer i
	_REAL_ mu(3, qm2ds%Mx), alpha(3), ft,mu_ground_state(3)

	! Evaluation of transition dipole moments

	write(6,*) 'Frequencies (eV) and Oscillator strengths (unitless)'
	write(6,"(8x,'Omega',12x,'fx',14x,'fy',14x,'fz',10x,'ftotal')")

	call trans_dipole(mu, alpha) ! cml-test
	do i=1,qm2ds%Mx
		ft = (2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*(mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
		write(6,"(i4,5g15.7)") i,qm2ds%e0(i), &
			(2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(1,i)**2, &
			(2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(2,i)**2, &
			(2.d0/3.d0)*(qm2ds%e0(i)/ONE_AU)*mu(3,i)**2, ft
	end do

	write(6,*) 
	write(6,*) 'Frequencies (eV) and Transitional Dipole Moments (AU)'
	write(6,"(8x,'Omega',12x,'fx',14x,'fy',14x,'fz',10x,'ftotal')")
	do i=1,qm2ds%Mx
		ft = (mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
		write(6,"(i4,5g30.20)") i,qm2ds%e0(i), &
			mu(1,i),mu(2,i),mu(3,i), ft
	end do



	write(6,*) ' Total energy of ground state (eV)'
	write(6,*) qm2ds%Eground

	write(6,*) ' Total energies of excited states (eV,AU)'
	do i=1,qm2ds%Mx
		write(6,*) i,qm2ds%Etot(i),qm2ds%Etot(i)/ONE_AU
	end do

	write(6,*) ' Required state (eV)',qm2ds%Ereq

	return
endsubroutine outDavidson
!
