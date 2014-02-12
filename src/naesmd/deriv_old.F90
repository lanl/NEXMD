#include "dprec.fh"
#include "assert.fh"

! Subroutine for analytic and numerical derivatives of ground or excited states
   subroutine deriv(sim,state) ! , xyz_in)
   use cosmo_C, only: qscnet,lm61,qdenet,gden,ipiden,ceps,solvent_model,potential_type;
   use communism !for numerical derivatives
   use qm2_davidson_module
   use qmmm_module
   use constants, only : KCAL_TO_EV, EV_TO_KCAL 

   implicit none
   
   type(simulation_t),target ::sim !communism module
   type(simulation_t),pointer::simpoint
   integer,intent(in),optional::state ! excited state where derivatives are calculated
   
   integer i,j,k,ihop,oldverbosity
   _REAL_::dxyz(qmmm_struct%nquant_nlink*3)
   _REAL_::dxyz_gs(qmmm_struct%nquant_nlink*3)
   _REAL_::dxyz1(3,qmmm_struct%nquant_nlink)
   _REAL_::xyz(3*qmmm_struct%nquant_nlink)
   _REAL_ :: Escf_right,Escf_left,E_ES_right,E_ES_left,h

   simpoint=>sim

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
      dxyz1=0.d0
      dxyz=0.d0
      dxyz_gs=0.d0

         ! Calculate ground state derivatives E_gr^x=E_nucl^x+E_el^x
         !   E_el^x=1/2 Tr((t^x+F^x) rho)
         ! qm2_get_exc_forces() is the SQM equivalent of DCART() in CEO 
      call qm2_get_exc_forces(dxyz1,qmmm_struct%qm_coords)

      !add solvent part (nuclear and electronic for ground state with COSMO
      !surface derivatives)
      if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
         if (ihop>0) call addfck(qm2ds%tz_scratch(qm2ds%Nb**2+1),qm2_struct%den_matrix)!this makes sure we're using rho_0
         call diegrd_testing(dxyz1,.true.) !derivative with surface on
      endif
      write(6,*)'qscnetaftergsder:',qscnet(:,3)
         ! Convert from kcal/A to eV/A.
      do i=1,qmmm_struct%nquant_nlink
         do j=1,3
            dxyz_gs((i-1)*3+j)=-dxyz1(j,i)*KCAL_TO_EV
         end do
      end do

      if (ihop>0) then
!CALCULATE EXCITED STATE DERIVATIVES
!Omega^x=Tr(F^x rhoTZ)+Tr(V^x(xi) xi^+)

!TERM 1: Tr(F^x rhoTZ)

         !CAN COSMO GRAIDENT tr(V_s(rho_0)^(x)(rho_0+T+Z)) BE DONE WITH ABOVE GS FORCES FOR
         !SPEED? diegrd_testing does triangular matrix, but this is symmetric so
         !that's ok. Now doing only electronic part of F^x, but
         !is this correct? JAB Comment

      ! Need to call dealloc_rhotz() when finished or in polishing
      ! add allocation to the big allocation/dealloc regime
      
      ! Get excited state density
      call calc_rhotz(ihop, qm2ds%rhoTZ,.true.)
      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%tz_scratch(1), &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))
      call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoTZ,'s')
      call getmodef(qm2ds%Np*qm2ds%Nh*2,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
         ihop,qm2ds%v0,qm2ds%tz_scratch(1))
      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%tz_scratch(1),qm2ds%rhoLZ, &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))
         write(6,*)'qscnetafterrhotz:',qscnet(:,3)
         dxyz1=0.d0

         !calculate vacuum derivative for term 1
         call dcart1(dxyz1,qm2_struct%den_matrix,qm2ds%rhoTZ,qmmm_struct%qm_coords)
        
         !add solvent part (symmetric only b/c symmetric matrix)  
         if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
                   !do i=1,lm61
                       qdenet(i,1)=0.0 !clear nuclear part of densmatrix
                       qscnet(i,1)=0.0
                       !qscnet(i,3)=qscnet(i,2) !clear nuclear part of cavity charges
                   !enddo
                 !call cosmo_1(qm2_struct%den_matrix) !recalculate cavity charges
                 !replace once center charges from rho0 with those from rhoTZ
                    do i=1,lm61
                       qdenet(i,3)=gden(i)*qm2ds%rhoTZ(ipiden(i)) !use rhoTZ as densmatrix
                    enddo
                 call diegrd_testing(dxyz1,.true.) !derivative
         endif
                 !write(6,*)'Den_matrix size=',size(qm2_struct%den_matrix),'RhoTZ size=',size(qm2ds%rhoTZ),' Nb=',qm2ds%Nb,&
                 !       ' gden=',size(gden),' lm61=',lm61
                 !j=1; k=1;
                 !do i=1,qm2ds%Nb
                 !write(6,*) j,k,'Den_matrix:',qm2ds%rhoTZ(j:k)
                 !write(6,*) j,k,'gden:',gden(j:k)
                 !write(6,*) j,k,'ipiden:',ipiden(j:k)
                 !j=j+i; k=k+i+1;
                 !enddo
      
         !add derivatives from term 1 to total derivative and convert to eV/A
         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do


!TERM 2: Tr(V^x(xi) xi^+) !!Why are there no Tr(V^x(xi_s)xi^+_u) type terms??
         !Symmetric part
         dxyz1=0.d0
         !pack triangular of symmetric part
         call packing(qm2ds%Nb, qm2ds%rhoLZ,qm2ds%tz_scratch(1), 's')
         !vacuum derivative
         call dcart2(dxyz1, qm2ds%tz_scratch(1),qmmm_struct%qm_coords)
   
         if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
            !call addfck(f,p) with correct density matrix to fill once
            call addfck(qm2ds%rhoTZ,qm2ds%tz_scratch(1))
            call diegrd_testing(dxyz1,.true.) !derivative
         end if

         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do
         !call dcart2(dxyz1, qm2ds%tz_scratch(1),qmmm_struct%qm_coords)

         dxyz1=0.d0

         !pack triangular matrix of antisymmetric part
         call packing(qm2ds%Nb, qm2ds%rhoLZ, qm2ds%tz_scratch(qm2ds%Nb**2+1),'u')
         !vacuum derivative
         call dcart2(dxyz1,qm2ds%tz_scratch(qm2ds%Nb**2+1),qmmm_struct%qm_coords)

         dxyz1=-dxyz1 !for changing sign to give xi transposed
         if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.30)) then
                 !call addfck(f,p) with correct density matrix to fill once
                 !center charge matrix qdenet for diegrd
                 call addfck(qm2ds%rhoTZ,qm2ds%tz_scratch(qm2ds%Nb**2+1))
                 call diegrd_testing(dxyz1,.true.)
         end if
      
         
         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)+dxyz1(j,i)*KCAL_TO_EV 
            end do
         end do  

         end if !ihop>0


!NUMERICAL DERIVATIVES
        !Currently wastes many resources by running full calculations for each
        !state. Need to store derivatives for all states since they are
        !calculated each time instead of calling deriv() for each state.

   elseif(qmmm_struct%ideriv.eq.0) then
   h=qmmm_struct%numder_step !step size
write(6,*) 'DOING NUMERICAL DERIVATIVES!!!!!!!!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
   oldverbosity=qm2ds%verbosity
   qm2ds%verbosity=0
   qmmm_nml%verbosity=0
   E_ES_left=0; E_ES_right=0;
   !Loop over the number of coordinates
   do i=1,qmmm_struct%nquant_nlink*3
      !Modify the coordinates and calculate energies

      xyz(i)=xyz(i)+h !left

      do k=1,qmmm_struct%nquant_nlink ! this can be eliminated by reorganizing xyz and dxyz
         do j = 1,3
            qmmm_struct%qm_coords(j,k)=xyz((k-1)*3+j)/0.529177249d0 !
         end do
      end do

      call do_sqm_davidson_update(simpoint,vgs=Escf_left,rx=qmmm_struct%qm_coords(1,:)&
           ,ry=qmmm_struct%qm_coords(2,:),rz=qmmm_struct%qm_coords(3,:))
        !rx=___ etc. might already be uncessary because qm_coords might be used in the subroutine
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
      !write(6,*)'E_ES_left=',E_ES_left*27.2116,' E_ES_right=',E_ES_right*27.2116,' dxyz=',dxyz(i)*27.2116,' ISTATE=',ihop
   enddo
   !Call gs and es calculations one more time to make sure later calculations
      !use the correct densities, etc.

   do k=1,qmmm_struct%nquant_nlink
      do j = 1,3
         qmmm_struct%qm_coords(j,k)=xyz((k-1)*3+j)/0.529177249d0
      end do
   end do

      call do_sqm_davidson_update(simpoint,rx=qmmm_struct%qm_coords(1,:)&
           ,ry=qmmm_struct%qm_coords(2,:),rz=qmmm_struct%qm_coords(3,:))
   qm2ds%verbosity=oldverbosity
   qmmm_nml%verbosity=oldverbosity
   end if !deriv flag


!WRITE RESULTS 
   if(qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 3) then
      !If verbosity level is greater than 3 we also print the force array on the
      !QM atoms
      write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms ground state calculation (eV)")')
      write (6,'("QMMM: state=",2i3)') ihop, qm2ds%struct_opt_state
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,-dxyz_gs(1+3*(j-1)), &
         -dxyz_gs(2+3*(j-1)), &
         -dxyz_gs(3+3*(j-1)), j=1,qmmm_struct%nquant_nlink)

      write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms excited state calculation (eV)")')
      write (6,'("QMMM: state=",2i3)') ihop, qm2ds%struct_opt_state
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,-dxyz(1+3*(j-1)),&
         -dxyz(2+3*(j-1)), &
         -dxyz(3+3*(j-1)), j=1,qmmm_struct%nquant_nlink)


      if (qmmm_nml%verbosity > 4) then
         !Also print info in KJ/mol
         write (6,'("QMMM:")')
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
   end if
simpoint%deriv_forces=dxyz_gs+dxyz !incorporating into sim type
  return
   end subroutine deriv
!
