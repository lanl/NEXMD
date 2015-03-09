#include "dprec.fh"
#include "assert.fh"

! Subroutine for analytic and numerical derivatives of ground or excited states
   subroutine deriv(sim,state) ! , xyz_in)
   use cosmo_C, only: doZ,lm61,nps,numat,qscnet,qscat,qdenet,phinet,gden,ipiden,ceps,solvent_model,potential_type,ediel,onsagE;
   use communism !for numerical derivatives
   use qm2_davidson_module
   use qmmm_module
   use constants, only : KCAL_TO_EV, EV_TO_KCAL 

   implicit none
   
   type(simulation_t),target ::sim !communism module
   type(simulation_t),pointer::simpoint
   integer,intent(in),optional::state ! excited state where derivatives are calculated
   
   integer i,j,k,ihop,oldverbosity
   _REAL_::dxyz(qmmm_struct%nquant_nlink*3),scratch(qmmm_struct%nquant_nlink**2)
   _REAL_::dxyz_gs(qmmm_struct%nquant_nlink*3)
   _REAL_::dxyz1(3,qmmm_struct%nquant_nlink),dxyz1_test(3,qmmm_struct%nquant_nlink)
   _REAL_::xyz(3*qmmm_struct%nquant_nlink)
   _REAL_::charges2(nps),acharges2(numat),density2(lm61)
   _REAL_ :: Escf_right,Escf_left,E_ES_right,E_ES_left,h
   _REAL_ :: tmp(qm2_struct%norbs,qm2_struct%norbs)
   simpoint=>sim
!write(6,*)'Calculating GS Derivatives'
!Collect coordinates in vector
   do i=1,qmmm_struct%nquant_nlink
      do j = 1,3
         xyz((i-1)*3+j)=qmmm_struct%qm_coords(j,i)
      end do
   end do
   
!Determine state to for which to calculate derivatives
   if(present(state)) then ! if different state is specified on input
      ihop=state
   else
      ihop=qmmm_struct%state_of_interest !otherwise use the current state of interest
   endif  

if(qmmm_struct%ideriv.eq.1) then !analytical derivatives

!CALCULATE GROUND STATE DERIVATIVES !At some point these should be saved and
!recalled for speed and convenience since deriv may be called multiple times for the
!same rho_0 and thus repeated
      dxyz1=0.d0; dxyz1_test=0.d0
      dxyz=0.d0
      dxyz_gs=0.d0

         ! Calculate ground state derivatives E_gr^x=E_nucl^x+E_el^x
         !   E_el^x=1/2 Tr((t^x+F^x) rho)
         ! qm2_get_exc_forces() is the SQM equivalent of DCART() in CEO 

      call qm2_get_exc_forces(dxyz1,qmmm_struct%qm_coords)

      !add solvent part (nuclear and electronic for ground state with COSMO
      !surface derivatives)
      if((solvent_model.gt.0).and.(solvent_model.ne.10)) then
        if((potential_type.eq.3).and.(ceps.gt.1.0)) then !have to specify ceps here because of division by 0 when eq 1
          call cosmo_1_tri(qm2_struct%den_matrix) !put qm2ds%den_matrix in the right place
          call diegrd(dxyz1_test) !derivative
        elseif(potential_type.eq.2) then
          call rcnfldgrad(dxyz1_test,qm2_struct%den_matrix,qm2_struct%norbs)   
          !write(6,*)'dxyz1=',dxyz1_test; dxyz1_test=0.d0  
          !call rcnfldgrad2(dxyz1_test,qm2_struct%den_matrix,qm2_struct%den_matrix,qm2ds%nb,.true.)
        endif     
      dxyz1=dxyz1+dxyz1_test
      endif

      do i=1,qmmm_struct%nquant_nlink
         do j=1,3
            dxyz_gs((i-1)*3+j)=-dxyz1(j,i)*KCAL_TO_EV
         end do
      end do

if (ihop>0) then
!CALCULATE EXCITED STATE DERIVATIVES
!Omega^x=Tr(F^x rhoTZ)+Tr(V^x(xi) xi^+)

     !TERM 1: Tr(F^x rhoTZ)

     ! Need to call dealloc_rhotz() when finished or in polishing
      ! add allocation to the big allocation/dealloc regime
      
      	! Get excited state density
      call calc_rhotz(ihop, qm2ds%rhoTZ,.true.)
      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%tz_scratch(1), &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))
      call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoTZ,'s')

	!Get transition density
      call getmodef(qm2ds%Np*qm2ds%Nh*2,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
         ihop,qm2ds%v0,qm2ds%tz_scratch(1))
      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%tz_scratch(1),qm2ds%rhoLZ, &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))
         dxyz1=0.d0; dxyz1_test=0.d0;

         !calculate vacuum derivative for term 1
         call dcart1(dxyz1,qm2_struct%den_matrix,qm2ds%rhoTZ,qmmm_struct%qm_coords)
         !add solvent part (symmetric only b/c symmetric matrix)  
         if(solvent_model.gt.0) then
                if((potential_type.eq.3).and.(ceps.gt.1.0)) then !ceps.gt.1.0 because of singularity in cosmo subroutines
                  call cosmo_1_tri_2(qm2ds%rhoTZ,density2,charges2,acharges2) !solvent and solute charges 
                  call diegrd2(dxyz1_test,density2,charges2,acharges2) !derivative
                elseif(potential_type.eq.2) then
                  call rcnfldgrad2(dxyz1_test,qm2_struct%den_matrix,qm2ds%rhoTZ,qm2ds%nb,.true.)
                endif
         dxyz1=dxyz1+0.5*dxyz1_test
         endif

         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do

!TERM 2: Tr(V^x(xi) xi^+) 
         !Symmetric part
         dxyz1=0.d0; dxyz1_test=0.d0; charges2=0.d0; acharges2=0.d0;
         !pack triangular of symmetric part and antisymmetric parts
         call packing(qm2ds%Nb, qm2ds%rhoLZ,qm2ds%tz_scratch(1), 's')
         call packing(qm2ds%Nb, qm2ds%rhoLZ,qm2ds%tz_scratch(qm2ds%nb**2+1),'u')
         !vacuum derivative
         call dcart2(dxyz1,qm2ds%tz_scratch(1),qmmm_struct%qm_coords)
         call dcart2(dxyz1,qm2ds%tz_scratch(qm2ds%nb**2+1),qmmm_struct%qm_coords)
         if(solvent_model.eq.1) then
         if((potential_type.eq.3).and.(ceps.gt.1.0)) then
            qscnet(:,1)=0.d0; qdenet(:,1)=0.d0; !Clear Nuclear Charges
            call cosmo_1_tri(qm2ds%tz_scratch(1)) !Fill Electronic Chrages
            call diegrd(dxyz1_test); !derivative
         elseif(potential_type.eq.2) then
            call rcnfldgrad_full(dxyz1_test,qm2ds%rhoLZ,qm2ds%nb); 
         end if
         end if
         dxyz1=dxyz1+0.5*dxyz1_test
         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do

!STATE SPECIFIC SOLVENT TERMS
         dxyz1=0.d0; dxyz1_test=0.d0; charges2=0.d0; acharges2=0.d0; density2=0.d0
        !Vertical Excitation Model
	if(solvent_model.eq.2) then
		!Get unrelaxed difference density matrix for the state to calculate derivatives for
		call calc_rhotz(ihop, qm2ds%rhoT,.false.)
      		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoT,qm2ds%tz_scratch(1), &
         		qm2ds%tz_scratch(qm2ds%Nb**2+1))
      		call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoT,'s')
                !Get relaxed or unrelaxed difference density matrix for the state specific state
                call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoTZ,doZ);
                call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
                call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoTZ,'s')
                
		!calculate derivatives
                if((potential_type.eq.3).and.(ceps.gt.1.0)) then !ceps.gt.1.0 because of glitch in cosmo subroutines
                  qscnet(:,1)=0.d0; qdenet(:,1)=0.d0; !Clear Nuclear Charges
                  call cosmo_1_tri(qm2ds%rhoTZ) !fill solvent charges
                  call cosmo_1_tri_2(qm2ds%rhoT,density2,charges2,acharges2) !fill solute charges 
                  call diegrd2(dxyz1_test,density2,charges2,acharges2) !derivative
                elseif(potential_type.eq.2) then
                  call rcnfldgrad2(dxyz1_test,qm2ds%rhoTZ,qm2ds%rhoT,qm2ds%nb,.false.)
                endif
	        dxyz1=dxyz1+0.5*dxyz1_test
                do i=1,qmmm_struct%nquant_nlink
                        do j=1,3
                                dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
                        end do
                end do
	endif

        !State Specific Model
	if(solvent_model.eq.4) then
                write(6,*)'WARNING:DERIVATIVES FOR STATE SPECIFIC SOLVENT ARE NONVARIATIONAL'
		!Get relaxed density matrix for the state specific state
                qm2ds%rhoT=0; qm2ds%rhoTZ=0;
	        call calc_rhotz(qmmm_struct%state_of_interest,qm2ds%rhoT,doZ); !rhoT will be rhoTZ_k
                call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoT,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
                qm2ds%rhoT=0;
                call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoT,'s')

		!Excited State Part
                if((potential_type.eq.3).and.(ceps.gt.1.0)) then !ceps.gt.1.0 because of singularity in cosmo subroutines
		  call cosmo_1_tri(qm2ds%rhoT) !fill solvent charges
                  call cosmo_1_tri_2(qm2ds%rhoTZ,density2,charges2,acharges2) !fill solute charges 
                  call diegrd2(dxyz1_test,density2,charges2,acharges2) !derivative
                elseif(potential_type.eq.2) then
                  call rcnfldgrad2(dxyz1_test,qm2ds%rhoT,qm2ds%rhoTZ,qm2ds%nb,.false.)
                endif
         	dxyz1=dxyz1+dxyz1_test
   
         	do i=1,qmmm_struct%nquant_nlink
            		do j=1,3
               			dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            		end do
         	end do

	        dxyz1=0.d0; dxyz1_test=0.d0; charges2=0.d0; acharges2=0.d0; density2=0.d0
		!Ground State part
		if((potential_type.eq.3).and.(ceps.gt.1.0)) then !ceps.gt.1.0 because of singularity in cosmo subroutines
                  call cosmo_1_tri_2(qm2_struct%den_matrix,density2,charges2,acharges2) !fill solute charges 
                  call diegrd2(dxyz1_test,density2,charges2,acharges2) !derivative
                elseif(potential_type.eq.2) then
                  call rcnfldgrad2(dxyz1_test,qm2ds%rhoT,qm2_struct%den_matrix,qm2ds%nb,.true.)
                endif
                dxyz1=dxyz1+dxyz1_test
      		do i=1,qmmm_struct%nquant_nlink
         		do j=1,3
            			dxyz_gs((i-1)*3+j)=dxyz_gs((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
        		end do
	     	end do

	endif

end if !ihop>0

!NUMERICAL DERIVATIVES
        !Currently wastes many resources by running full calculations for each
        !state. Need to store derivatives for all states since they are
        !calculated each time instead of calling deriv() for each state.

   elseif(qmmm_struct%ideriv.eq.2) then
   h=qmmm_struct%numder_step !step size
   oldverbosity=qmmm_nml%verbosity
   qmmm_nml%verbosity=0
   E_ES_left=0; E_ES_right=0;
   !Loop over the number of coordinates
   do i=1,qmmm_struct%nquant_nlink*3
   !write(6,*)'den_c',qm2_struct%den_matrix
   !write(6,*)'fock_c',qm2_struct%fock_matrix
   !Modify the coordinates and calculate energies

      xyz(i)=xyz(i)+h !left
      do k=1,qmmm_struct%nquant_nlink ! this can be eliminated by reorganizing xyz and dxyz
         do j = 1,3
            qmmm_struct%qm_coords(j,k)=xyz((k-1)*3+j)/0.529177249d0 !
         end do
      end do

      call do_sqm_davidson_update(simpoint,vgs=Escf_left,rx=qmmm_struct%qm_coords(1,:)&
           ,ry=qmmm_struct%qm_coords(2,:),rz=qmmm_struct%qm_coords(3,:))
        !Note: rx=___ etc. might already be uncessary because qm_coords might be used in the subroutine
        E_ES_left=sim%naesmd%Omega(ihop)

      xyz(i)=xyz(i)-2*h !right
      do k=1,qmmm_struct%nquant_nlink
         do j = 1,3
            qmmm_struct%qm_coords(j,k)=xyz((k-1)*3+j)/0.529177249d0 ! 
         end do
      end do

      call do_sqm_davidson_update(simpoint,vgs=Escf_right,rx=qmmm_struct%qm_coords(1,:)&
           ,ry=qmmm_struct%qm_coords(2,:),rz=qmmm_struct%qm_coords(3,:))
        E_ES_right=sim%naesmd%Omega(ihop)

      xyz(i)=xyz(i)+h !back to center
      do k=1,qmmm_struct%nquant_nlink
         do j = 1,3
            qmmm_struct%qm_coords(j,k)=xyz((k-1)*3+j)/0.529177249d0 ! 
         end do
      end do

      call do_sqm_davidson_update(simpoint,rx=qmmm_struct%qm_coords(1,:)&
           ,ry=qmmm_struct%qm_coords(2,:),rz=qmmm_struct%qm_coords(3,:)) 

      !Calculate derivative
      dxyz_gs(i)=-(Escf_left-Escf_right)/(2*h)*27.2116 !GS
      if(ihop>0) dxyz(i)= -(E_ES_left- E_ES_right)/(2*h)*27.2116 !ES
      write(6,*)'ES:',E_ES_left,E_ES_right,'dxyz(i)=',dxyz(i)
   enddo

   qm2ds%verbosity=oldverbosity
   qmmm_nml%verbosity=oldverbosity

   end if !deriv flag


!WRITE RESULTS 
   if(qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 3) then
      !If verbosity level is greater than 3 we also print the force array on the
      !QM atoms
      !write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms ground state calculation (eV)")')
      write (6,'("QMMM: state=",2i3)') ihop, qm2ds%struct_opt_state
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,-dxyz_gs(1+3*(j-1)), &
         -dxyz_gs(2+3*(j-1)), &
         -dxyz_gs(3+3*(j-1)), j=1,qmmm_struct%nquant_nlink)
	if(ihop>0) then
      !write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms excited state calculation (eV)")')
      write (6,'("QMMM: state=",2i3)') ihop, qm2ds%struct_opt_state
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,-dxyz(1+3*(j-1)),&
         -dxyz(2+3*(j-1)), &
         -dxyz(3+3*(j-1)), j=1,qmmm_struct%nquant_nlink)


      if (qmmm_nml%verbosity > 4) then
         !Also print info in KJ/mol
         !write (6,'("QMMM:")')
         write (6,'("QMMM: Excited State Forces on QM atoms from SCF calculation (KJ/mol)")')
         write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j, &
            -dxyz(1+3*(j-1))*EV_TO_KCAL*4.184d0, &
            -dxyz(2+3*(j-1))*EV_TO_KCAL*4.184d0, &
            -dxyz(3+3*(j-1))*EV_TO_KCAL*4.184d0,j=1,qmmm_struct%nquant_nlink)
         write(6,'(/,''  QMMM: QM Region Cartesian Coordinates (*=link atom) '')')
         write(6,'(''  QMMM: QM_NO. MM_NO.'',2X,''ATOM'',9X,''X'',9X,''Y'',9X,''Z'')')
           do I=1,qmmm_struct%nquant
              WRITE(6,'("  QMMM:",I6,2X,I7,6X,I2,3X,3F10.4)') I, &
                 qmmm_struct%iqmatoms(i), &
                 I, &
                 (qmmm_struct%qm_coords(J,I),J=1,3)
           end do
      end if
	endif
   end if
simpoint%deriv_forces=dxyz_gs+dxyz !incorporating into sim type
  return
   end subroutine deriv
!
