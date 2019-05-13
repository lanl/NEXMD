#include "dprec.fh"
#include "assert.fh"

! Subroutine for analytic derivativ
   subroutine deriv(dxyz,state) ! , xyz_in)
   use cosmo_C, only: ceps,solvent_model,potential_type;
   use qm2_davidson_module
   use qmmm_module
   use constants, only : KCAL_TO_EV, EV_TO_KCAL 

   implicit none

!-----------OLD VARS FROM CEO-------------------


   _REAL_ fPi,fbar,ff0,ff1,ff11,ddot
   parameter(fPi=3.1415926535898d0)
   parameter(fbar=0.05d0)  ! maximum dE in numerical derivative, eV.

   integer Na,Nm,N3,i,j,k,im,one,istate,ip,ih
   _REAL_::xx1(qmmm_struct%nquant_nlink)
   _REAL_::yy1(qmmm_struct%nquant_nlink)
   _REAL_::zz1(qmmm_struct%nquant_nlink)
   _REAL_,intent(out) :: dxyz(qmmm_struct%nquant_nlink*3)
   _REAL_::dxyz1(3,qmmm_struct%nquant_nlink)
   _REAL_::xyz(3*qmmm_struct%nquant_nlink)
   _REAL_ Omega,Omega_,EE0,EE0_,dEE0,f,dOmega,t,d,d1,ff
   
   integer mdflag
   character*20 datetime, keywr*40, machname*36
   _REAL_ Egr, E0_f
   integer, allocatable :: imodimp(:), jmodimp(:)
   integer imimp
   integer get_time
   _REAL_ tmark
   integer jm
   real time11,time12,time13,time23
   save time11,time12,time13,time23

!-----------END OLD VARS FROM CEO-------------------
!-----------NEW VARS-------------------
   
   integer,intent(in),optional::state ! excited state where derivatives are calculated

   _REAL_ :: Escf
   _REAL_, pointer :: rhogr(:)
   integer :: ier
   integer :: ihop
   integer :: ideriv ! Temporary until can be passed from input file
   integer :: iderivfl ! Temporary until we can find a use for it elsewhere
   logical :: calc_Z	

   ! These are dummy variables because sqm_energy needs them to run
   ! They don't appear to get assigned by sqm when calculating SPE
   ! These will probably need to be replaced by real values inside
   ! a struct when the time comes to use solvents and the like
   
   _REAL_ :: born_radii(1000), one_born_radii(1000)
   _REAL_ :: intdiel, extdiel, Arad

!-----------END NEW VARS-------------------

   N3=3*Na
   one=1
   istate=1
   ff0=0.0
   ff1=1.0
   ff11=-1.0
   Omega=0.0
   Omega_=0.0
!

! Cartesian coordinates are now located in qmmm_module
! But they are only there AFTER sqm is run outside this subroutine b/c it's passed in
! as an argument. Might want to take the sqm call outside or copy coords earlier
   xx1=qmmm_struct%qm_coords(1,:) ! xyz_in(1,:)
   yy1=qmmm_struct%qm_coords(2,:) ! xyz_in(2,:)
   zz1=qmmm_struct%qm_coords(3,:) ! xyz_in(3,:)

   do i=1,qmmm_struct%nquant_nlink
      do j = 1,3
         xyz((i-1)*3+j)=qmmm_struct%qm_coords(j,i) ! xyz_in(j,i)
      end do
   end do


! NEW CODE BLOCK 1
   ! The above block calculates the desired Hamiltonian and calcs the GS.
   ! SQM does this all at once
   ! No solvent effects have been added yet
   E0_f=qmmm_struct%elec_eng+qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm


! END NEW CODE BLOCK 1


! NEW CODE BLOCK 2
   rhogr=>qm2_struct%den_matrix
! END NEW CODE BLOCK 2

! REFORMED CODE BLOCK 3
   ! ihop is part of a common block in CEO and tracks the path
   ! of the non-adiabatic hops between excited states
   ! This is a trivial assignment until NAD is incorporated

   ihop=qm2ds%struct_opt_state ! this for standalone sqm
   write(6,*)"state=",state
   if(present(state)) then ! needed to explicitly specify the state
      ihop=state
   end if
   Write(6,*)"ihop=",ihop
   if (ihop.gt.0) then ! prepare for excited state gradients
      istate=ihop
! END REFORMED CODE BLOCK 3

! Calculate excited states and excited state density matrix T+Z in rhoTZ

! NEW CODE BLOCK 4
   ! Call to LZ() replaced by ex. st. den. matrix wrapper
   ! Need to call dealloc_rhotz() when finished or in polishing
   ! add allocation to the big allocation/dealloc regime
   ! Omega_f assignment was removed because this is taken care of
   ! in the Davidson procedure within
      calc_Z=.true.	
      call calc_rhotz(istate, qm2ds%rhoTZ,solvent_model,calc_Z)
! END NEW CODE BLOCK 4


! Get appropriate transition density and rhoTZ to AO
! NEW CODE BLOCK 5
      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,qm2ds%tz_scratch(1), &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))

      call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%rhoTZ,'s')
      call getmodef(qm2ds%Np*qm2ds%Nh*2,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
         istate,qm2ds%v0,qm2ds%tz_scratch(1))

      call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%tz_scratch(1),qm2ds%rhoLZ, &
         qm2ds%tz_scratch(qm2ds%Nb**2+1))

! END NEW CODE BLOCK 5

! REFORMED CODE BLOCK 6
   end if
   
   iderivfl=1    ! This is a flag for analytical derivatives. Not sure if needed in SQM
   ideriv=2

   if(ideriv>=2) then   ! Fast GS MOPAC derivatives
! END REFORMED CODE BLOCK 6

! NEW CODE BLOCK 7
      dxyz1=0.d0
      ! Calculate ground state derivatives E_gr^x=E_nucl^x+E_el^x
      !   E_el^x=1/2 Tr((t^x+F^x) rho)

      ! qm2_get_exc_forces() is the SQM equivalent of DCART() in CEO 
      call qm2_get_exc_forces(qm2_params,qmmm_nml,qm2_rij_eqns,qmmm_mpi,qm2_struct, dxyz1,qmmm_struct%qm_coords)
      !!JAB The function name should be changed? This is ground state derivative
      ! Convert from kcal/A to eV/A.
      ! Maps matrix onto vector. Only works when -Mbounds is not specified as compiler flag
      ! dxyz(1:qmmm_struct%nquant_nlink*3) = -dxyz1(1:qmmm_struct%nquant_nlink*3)*KCAL_TO_EV
      ! This will work with or without Mbounds. Less efficient.
      do i=1,qmmm_struct%nquant_nlink
         do j=1,3
            dxyz((i-1)*3+j)=-dxyz1(j,i)*KCAL_TO_EV
         end do
      end do
! END NEW CODE BLOCK 7

! NEW CODE BLOCK 8

!COSMO GRADIENT HERE? NEED TO IMPLEMENT COSMO GRADIENT IN THIS PART? OR CAN
!COSMO GRAIDENT tr(V_s(rho_0)^(x)(rho_0+T+Z)) BE DONE WITH ABOVE GS FORCES FOR
!SPEED? diegrd_testing does triangular matrix. use only symmetric part and only
!surface derivatives once! JAB Comment


      if (ihop>0) then
         dxyz1=0.d0
         ! Within dcart1(), dhc1() is called. However, between CEO and SQM, dhc1() is called
         ! a different number of times. The numbers returned are the same, however.
         ! This could be a discrepancy in the indexing or the looping which may be moot,
         ! but it may be worth investigating. CML 6/26/12
         call dcart1(qmmm_nml, qm2_params, qmmm_mpi, qm2_rij_eqns, qm2_struct, qm2ds, qmmm_struct, dxyz1,qm2_struct%den_matrix,qm2ds%rhoTZ, &
            qmmm_struct%qm_coords)
         ! Convert from kcal/A to eV/A.
         ! Maps matrix onto vector. Only works when -Mbounds is not specified as compiler flag
         ! dxyz(1:qmmm_struct%nquant_nlink*3) = dxyz(1:qmmm_struct%nquant_nlink*3)-dxyz1(1:qmmm_struct%nquant_nlink*3)*KCAL_TO_EV
        
         !Add solvent part (symmetric only b/c symmetric matrix)  
         if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
                 !call addfck(f,p) with correct density matrix to fill once
                 !center charge matrix qdenet for diegrd
                 call addfck(qm2_params, qmmm_nml,cosmo_c_struct, qm2ds%tz_scratch(qm2ds%Nb**2+1),qm2ds%tz_scratch(1))
                 call diegrd_testing(dxyz1,.true.)
         end if
         !same as above loop now commented out
         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do

! END NEW CODE BLOCK 8

! Term 2: Tr(V^x(xi) xi^+)

! NEW CODE BLOCK 9
      ! Symmetric part
         dxyz1=0.d0

         call packing(qm2ds%Nb, qm2ds%rhoLZ, qm2ds%tz_scratch(1), 's')

         call dcart2(qm2_params, qmmm_nml,qm2_rij_eqns, qm2_struct,qm2ds,qmmm_struct,dxyz1, qm2ds%tz_scratch(1), &
            qmmm_struct%qm_coords)

      ! Add solvent for symmetric part !!JAB
      ! Likely adding the SAS derivative twice for symmetric and asymmetric is
      ! incorrect. Using the wrong density matrix for this? 
      ! Fixed. Adding only once (flag in diegrd_testing)
      ! DOES THIS GIVE THE CORRECT FORM OF TERM 2 SCREENING FACTOR??
      ! Is this transposing xi for the trace?
      if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
                 !call addfck(f,p) with correct density matrix to fill once
                 !center charge matrix qdenet for diegrd
                 call addfck(qm2_params, qmmm_nml,cosmo_c_struct, qm2ds%tz_scratch(qm2ds%Nb**2+1),qm2ds%tz_scratch(1))
                 call diegrd_testing(dxyz1,.false.)
      end if

         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do

         ! Antisymmetric part
         dxyz1=0.d0

         call packing(qm2ds%Nb, qm2ds%rhoLZ, qm2ds%tz_scratch(1), 'u')

         call dcart2(qm2_params, qmmm_nml,qm2_rij_eqns,qm2_struct,qm2ds,qmmm_struct,dxyz1, qm2ds%tz_scratch(1), &
            qmmm_struct%qm_coords)

      ! Add solvent for antisymmetric part !!JAB
      if((ceps.gt.1.d0).and.(solvent_model.gt.0).and.(potential_type.eq.3)) then
                 !call addfck(f,p) with correct density matrix to fill once
                 !center charge matrix qdenet for diegrd
                 call addfck(qm2_params, qmmm_nml,cosmo_c_struct, qm2ds%tz_scratch(qm2ds%Nb**2+1),qm2ds%tz_scratch(1))
                 call diegrd_testing(dxyz1,.false.)
      end if

         do i=1,qmmm_struct%nquant_nlink
            do j=1,3
               dxyz((i-1)*3+j)=dxyz((i-1)*3+j)-dxyz1(j,i)*KCAL_TO_EV
            end do
         end do  
! END NEW CODE BLOCK 9

!		else	! "Standard CEO derivatives"

      end if
   end if

   if(qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 3) then	
      !If verbosity level is greater than 3 we also print the force array on the
      !QM atoms
      write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms excited state calculation (KCal/mol)")')
      write (6,'("QMMM: state=",2i3)') ihop, qm2ds%struct_opt_state
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,-dxyz(1+3*(j-1))*EV_TO_KCAL, &
         -dxyz(2+3*(j-1))*EV_TO_KCAL, &
         -dxyz(3+3*(j-1))*EV_TO_KCAL, j=1,qmmm_struct%nquant_nlink)

      if (qmmm_nml%verbosity > 4) then
         !Also print info in KJ/mol
         write (6,'("QMMM:")')
         write (6,'("QMMM: Forces on QM atoms from SCF calculation (KJ/mol)")')
         write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j, &
            -dxyz(1+3*(j-1))*EV_TO_KCAL*4.184d0, &
            -dxyz(2+3*(j-1))*EV_TO_KCAL*4.184d0, &
            -dxyz(3+3*(j-1))*EV_TO_KCAL*4.184d0,j=1,qmmm_struct%nquant_nlink)
      end if
   end if
   

   return
   end subroutine deriv
!
