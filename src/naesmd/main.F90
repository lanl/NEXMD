#include "dprec.fh"
#include "assert.fh"

#ifndef STDOUT
#define STDOUT 6
#endif


program MD_Geometry

    use qmmm_module,only:deallocate_qmmm
    use qm2_davidson_module
    use naesmd_constants
    use fewest_switches
    use writeoutput_module
    use langevin_temperature
    use cadiab_module
    use additional_subroutines
    use random
    use quantum_prop_add
    use verlet_module
    use communism
    use naesmd_module!,only: allocate_naesmd_module2
    use md_module
    use rksuite_90, only:setup, range_integrate, rk_comm_real_1d
    use cosmo_C, only : cosmo_C_structure
    implicit none
    !
    !--------------------------------------------------------------------
    !
    !  Molecular dynamics uses cartesian coordinates, analytical
    !  or numerical derivatives, and
    !  Newtonian (usual velocity-Verlet), Langevin thermostat velocity-Verlet,
    !  or Berendsen thermostat (all molecules collide with heat bath at
    !  once) algorithms.

    !  input.ceon file with CEO and MDQT parameters,
    !  as well as initial coordinates, velocities, and quantum coefficients;
    !
    !--------------------------------------------------------------------

    !The instances of all sim structures
    type(qm2_davidson_structure_type),target :: qm2ds_notp
    type(qmmm_struct_type), target :: qmmm_struct_notp
    type(qm2_structure), target :: qm2_struct_notp
    type(cosmo_C_structure), target :: cosmo_C_struct_notp
    type(naesmd_structure), target :: naesmd_struct_notp
    type(xlbomd_structure), target :: xlbomd_struct_notp
    type(md_structure), target :: md_struct_notp
    type(rk_comm_real_1d), target :: rk_struct_notp
    type(qm2_params_type),target :: qparams_struct_notp
    type(qmmm_nml_type),target :: qnml_struct_notp
    type(qm2_rij_eqns_structure),target:: rij_struct_notp
    type(qm_gb_structure),target :: gb_struct_notp
    type(qmmm_mpi_structure),target :: qmpi_struct_notp
    type(qmmm_opnq_structure),target :: opnq_struct_notp
    type(qmmm_div_structure), target :: div_struct_notp
    type(qmmm_vsolv_type), target :: vsolv_struct_notp
    type(qmmm_scratch_structure), target:: scratch_struct_notp
    type(qm_ewald_structure), target :: ewald_struct_notp

    !The sim structure, a set of pointers and some values
    type(simulation_t),target::sim_notp
    !The pointer to sim structure, more rapidly passed into subroutines
    type(simulation_t),pointer::sim
    integer :: Nsim = 1
   
    sim=>sim_notp
    call init0_simulation(sim,Nsim)
    call setp_simulation(sim,qmmm_struct_notp,qm2ds_notp, &
                qm2_struct_notp,naesmd_struct_notp,md_struct_notp,rk_struct_notp, &
                cosmo_c_struct_notp,xlbomd_struct_notp,qparams_struct_notp, &
                qnml_struct_notp,rij_struct_notp,gb_struct_notp,qmpi_struct_notp, &
                opnq_struct_notp,div_struct_notp,vsolv_struct_notp,scratch_struct_notp, &
                ewald_struct_notp)
         sim%id=0
    call nexmd_sim(sim)

contains


subroutine nexmd_sim(sim)
    type(simulation_t),pointer::sim

    integer i,j,k
    integer imdqt
    integer itime1,itime11,itime2,itime3,get_time
    _REAL_ time11,time12
    _REAL_ t_start,t_finish
    character*30 datetime, machname*36
    integer inputfdes
    character(100) ::  filename
    
   
    call get_date(datetime)
    call get_machine(machname)
    write (6,*)
    write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
    write (6,8) '| MD execution started at  ',datetime,' |'
    write (6,7) '| Computer: ',   machname,'|'
    write (6,*) '|________________________________________________|'
    write (6,*)
    itime1=get_time()
    
    call init_main(sim, inputfdes)

    !Put derivative variables into module

    !***********************************************************
    !Initialization and single point calculation or zeroeth
    !step of adiabatic or nonadiabatic dynamics
    !
    !Open input file, initialize variables, run first calculation
    !***********************************************************
    if(sim%Nsim.eq.1) then
        open (inputfdes,file='input.ceon',status='old')
    else
        write (filename, "(a6,i4.4,a5)") "input_", sim%id, ".ceon"
        open (inputfdes,file=trim(filename),status='old')
    endif 
    write(6,*)"init_sqm" 
    call init_sqm(sim,inputfdes,STDOUT) ! involves call to Davidson
    write(6,*)"init_sqm_done" 
    close(inputfdes)

    sim%naesmd%nbasis=sim%dav%Nb  ! this is number of atomic orbitals
    sim%nbasis=sim%naesmd%nbasis  ! not to confuse with Ncis or Nrpa !!!

    if(sim%excN.ne.0) call allocate_naesmd_module2(sim%naesmd, sim%dav%Ncis,sim%nbasis,sim%excN)

    ! calling the first time to check the quirality in what follows
    !  sim%naesmd%uumdqtflag =0 means that ceo is called by the first time,
    !     so m.o. matrix is stored
    !  sim%naesmd%uumdqtflag =1 means that the quirality of m.o. must be checked
    sim%naesmd%uumdqtflag=0
    !
    do i=1,sim%excN
        sim%naesmd%iorden(i)=i
        sim%naesmd%iordenhop(i)=0
    end do
     
    sim%naesmd%vgs=sim%naesmd%E0;
    sim%naesmd%vmdqt(1:sim%excN)=sim%naesmd%Omega(1:sim%excN)+sim%naesmd%vgs

    itime2=get_time()
    time12=real(itime2-itime1)/100

    !Calculate derivatives
    call cpu_time(t_start)
    if (sim%qmmm%ideriv.gt.0) then
        call deriv(sim,sim%naesmd%ihop)
    elseif (sim%naesmd%nstep.gt.0) then
        write(6,*)'Must choose derivative type greater than zero for dynamics'
        stop
    endif
    call cpu_time(t_finish)
    sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
    

    !Derivatives to forces
    call deriv2naesmd_accel(sim)

    !##########################################################
    !! Molecular dynamics with Quantum transition main loop
    !##########################################################
    sim%naesmd%icontw=1
    !Open output files
    if(Nsim.eq.1) then
	    call open_output(sim,sim%ibo,sim%naesmd%tfemto,sim%md%imdtype,sim%lprint)
	    sim%naesmd%ihop=sim%qmmm%state_of_interest 
    else
	    call open_output_multi(sim,sim%ibo,sim%naesmd%tfemto,sim%md%imdtype,sim%lprint)
	    sim%naesmd%ihop=sim%qmmm%state_of_interest 
    endif
    call writeoutputini(sim,sim%ibo,sim%naesmd%yg,sim%lprint)

    sim%rk_comm%tmax=sim%naesmd%nstep*sim%naesmd%dtmdqt
    do i =1,sim%excN
	sim%rk_comm%thresholds(i)=1.0d0
	sim%rk_comm%thresholds(i+sim%excN)=6.29d0
    enddo
    call setup(sim%rk_comm,sim%naesmd%tfemto,sim%naesmd%yg,sim%rk_comm%tmax,sim%rk_comm%rk_tol,sim%rk_comm%thresholds, &
	'M','R')
    do imdqt=1,sim%naesmd%nstep !Main loop
        !Classical propagation step - BOMD or NAESMD
        write(6,*)"Begin classical propagation step #",imdqt
        sim%naesmd%tfemto=sim%naesmd%tfemto+sim%naesmd%dtmdqt*convtf
        sim%qmmm%num_qmmm_calls=imdqt !necessary for BO dynamics

        if(sim%naesmd%state.eq.'exct'.and.sim%ibo.ne.1) then
            !NAESMD-propagate=====================================================================================
            sim%qmmm%state_of_interest=sim%naesmd%ihop;
            call verlet1(sim)
            write(6,*)'Begin nonadiabatic couplings and crossings calculations'
            call initialize(sim,sim%naesmd%yg)
            !****************************************************************
            ! calculation of the energies(=sim%naesmd%vmdqtold) and nacT(=sim%naesmd%cadiabold) values at the beginning
            ! of the classical step
            !*******************************************************
            call cadiaboldcalc(sim,imdqt)
            !**************************************************************
            ! calculation of the energies(=sim%naesmd%vmdqtnew) and nacT(=sim%naesmd%cadiabnew) values at the end
            ! of the classical step
            !*******************************************************
            call cadiabnewcalc(sim)
            !**************************************************************
            ! check crossing, if crossing takes place, sim%naesmd%nquantumstep=sim%naesmd%nstepcross*sim%naesmd%nquantumreal
            !*************************************************************** 
            call do_trivial_wrap(sim) 
            write(6,*)'End nonadiabatic couplings calculation'
            write(6,*)"End classical propagation step #",imdqt 
            !--------------------------------------------------------------------
            ! Loop for quantum propagation steps
            ! that implies CEO energy calculations
            !--------------------------------------------------------------------
            call quantum_propagation(sim, imdqt)         
            !--------------------------------------------------------------------
            ! last part of velocity verlet algorithm
            ! for ehrenfest should go after evalhop
            call verlet2(sim)
            !--------------------------------------------------------------------
            ! analyze the hopping
            !--------------------------------------------------------------------
            call evalhop(sim, sim%rk_comm, sim%lprint, sim%rk_comm%tend, sim%rk_comm%tmax, &
                         sim%rk_comm%rk_tol, sim%rk_comm%thresholds, &
                sim%Na, sim%naesmd%yg, sim%naesmd%cross)
            !--------------------------------------------------------------------

            !call check_ntotcoher(sim) 
             call decoherence_E0_and_C(sim)!currently this is not available 
        !--------------------------------------------------------------------
        else
            !BOMD-propagate======================================================================================
            call verlet(sim)
            sim%naesmd%ihop=sim%qmmm%state_of_interest
        end if

        ! evaluation of current kinetic energy
        call temperature(imdqt,sim%naesmd)
        if(imdqt.eq.1) then
            sim%naesmd%icont=1
            sim%naesmd%icontpdb=sim%naesmd%icontini
        end if

        if(sim%naesmd%icontw.ne.sim%naesmd%nstepw) then
            sim%naesmd%icontw=sim%naesmd%icontw+1
        else
            sim%naesmd%icontw=1
            call writeoutput(sim,sim%ibo,sim%naesmd%yg,sim%lprint,sim%naesmd%cross)
        end if
    end do

    !***********************************
    !      ttime=dtime(tarray)
    !      write(19,*) tarray(1)

779 format(F18.10,10000(1X,I3,3(1X,F18.10)))
886 format(10000(1X,F7.4))

    !##########################################################
    !! END of the Molecular dynamics with Quantum transition main loop
    !##########################################################

    itime11=get_time()
    time11=real(itime11-itime1)/100
    itime11=time11
    itime1=MOD(itime11,60)
    itime11=itime11/60
    itime2=MOD(itime11,60)
    itime11=itime11/60
    itime3=MOD(itime11,24)
    itime11=itime11/24

    call get_date(datetime)

    write (6,*)
    write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
    write (6,8) '| MD normal termination at ', datetime,' |'
    write (6,9) '| MD total CPU time	    ',  time11, ' seconds |'
    write(6,301) itime11,itime3,itime2,itime1
    write(6,*)

    ! Cumulative execution time of various code blocks
    write(6,*) ' SQM (ground sim%naesmd%state) took overall [s]:'
    write(6,'(g16.3)') sim%time_sqm_took
    write(6,*) ' Davidson (excited states) took overall [s]:'
    write(6,'(g16.3)') sim%time_davidson_took
    write(6,*) ' deriv (adiabatic forces) took overall [s]:'
    write(6,'(g16.3)') sim%time_deriv_took
    write(6,*) ' nacT (NA derivatives) took overall [s]:'
    write(6,'(g16.3)') sim%time_nact_took
    write (6,*) '|________________________________________________|'
    write (6,*)

7   format (' ',A12,' ',A35,A2)
8   format (' ',A27,'  ',A30,A2)
9   format (' ',A20,'    ',g16.5,A10)
301 format(' |     ',i2,' days ',i2,' hours ',i2,' minutes ',i2, &
        ' seconds     |')
302 format(A3,3f12.6)
303 format(i6,f9.3,g11.3,5g16.8)
304 format('STEP:',i5,f9.3,6g13.5,g10.2,f16.10)
305 format(a,i5,6g16.8,g12.4)
306 format(i7,10000g12.4)
307 format(a,g8.4,6g16.8,g12.4)
999 format(I3,1X,1000(1X,F18.10))

    call flush(6)

    stop

98  stop 'bad input'
end subroutine 

    !
    !********************************************************************
    !
    subroutine init_main(sim, inputfdes)
        use naesmd_constants
	use communism, only : simulation_t 
        implicit none
        type(simulation_t),pointer::sim
        _REAL_,allocatable::xx(:),yy(:),zz(:)
        integer, intent (inout) :: inputfdes
        integer ifind
        integer :: i,j,k,ii,jjj
        integer :: Na, Nm, N1, N2, N3
        integer slen
	integer :: itime1
        character*(150) txt
        character(100) ::  filename
        _REAL_ :: rk_tolerance
        logical moldyn_found
        ! variables of the moldyn namelist
	    integer natoms
	    integer bo_dynamics_flag,exc_state_init,n_exc_states_propagate
	    integer out_count_init
	    _REAL_ time_init,time_step
	    integer n_class_steps,n_quant_steps,quant_coeffs_reinit
	    _REAL_ num_deriv_step,therm_temperature
	    integer therm_type
	    _REAL_ berendsen_relax_const
	    integer heating,heating_steps_per_degree,out_data_steps
	    integer out_coords_steps
	    _REAL_ therm_friction
	    integer rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag
	    _REAL_ quant_step_reduction_factor
	    _REAL_ decoher_e0,decoher_c
	    integer decoher_type,dotrivial
	    _REAL_ deltared
	    integer iredpot,nstates
	    namelist /moldyn/ natoms,bo_dynamics_flag,exc_state_init, &
		n_exc_states_propagate,out_count_init,time_init, &
		rk_tolerance,time_step,n_class_steps,n_quant_steps, &
		num_deriv_step, &
		therm_temperature,therm_type, &
		berendsen_relax_const,heating, &
		heating_steps_per_degree,out_data_steps,out_coords_steps, &
		therm_friction,rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag, &
		quant_step_reduction_factor,decoher_e0,decoher_c,decoher_type,dotrivial,&
		iredpot,nstates,deltared
        !dtnact is the incremental time to be used at nact calculation
        sim%naesmd%dtnact=0.002d0



        ! definition of variables
        i=0
        sim%naesmd%ktbig(i)=char(48)//char(48)//char(48)//char(48)
        do i=1,9
            sim%naesmd%ktbig(i)=char(48)//char(48)//char(48)//char(48+i)
        end do

        ii=9
        do j=1,9
            do i=0,9
                ii=ii+1
                sim%naesmd%ktbig(ii)=char(48)//char(48)//char(48+j)//char(48+i)
            end do
        end do

        ii=99
        do k=1,9
            do j=0,9
                do i=0,9
                    ii=ii+1
                    sim%naesmd%ktbig(ii)=char(48)//char(48+k)//char(48+j)//char(48+i)
                end do
            end do
        end do

        ii=999
        do jjj=1,9
            do k=0,9
                do j=0,9
                    do i=0,9
                        ii=ii+1
                        sim%naesmd%ktbig(ii)=char(48+jjj)//char(48+k)//char(48+j)//char(48+i)
                    end do
                end do
            end do
        end do

        ! ido= Flag indicating the state of the computation
        ! use by divprk propagator
        !ido=1
        !idocontrol=1
        sim%naesmd%conthop=0
        sim%naesmd%conthop2=0
        !
        !--------------------------------------------------------------------
        !

        inputfdes=12*10000+sim%id
        ! copy input file
        if(sim%Nsim.eq.1) then
              open (inputfdes,file='input.ceon',status='old')
        else
              write (filename, "(a6,i4.4,a5)") "input_", sim%id, ".ceon"
              open (inputfdes,file=trim(filename),status='old')
        endif 
        j=0
        do while(readstring(inputfdes,sim%naesmd%cardini,slen).ge.0)
            j=j+1

            do i=1,slen
                sim%naesmd%txtinput(j)(i:i)=sim%naesmd%cardini(i:i)
            end do

            if(sim%naesmd%cardini(1:6).eq.'$COORD') then
                sim%naesmd%txtinicoord=j
            end if

            if(sim%naesmd%cardini(1:9).eq.'$ENDCOORD') then
                sim%naesmd%txtendcoord=j
            end if
        end do
        close(inputfdes)

        sim%naesmd%jend=j
        !
        !--------------------------------------------------------------------
        !
        !  Reading moldyn parameters
        !
        !--------------------------------------------------------------------
        !
        ! Default parameters for moldyn
        natoms=1000 !max number of atoms
        bo_dynamics_flag=1 ! 0-non-BO, 1-BO
        exc_state_init=0 ! initial excited state
        n_exc_states_propagate=0 ! total number of excited state to propagate
        out_count_init=0 ! iniit count for output files
        time_init=0.d0 ! initial time, fs
        rk_tolerance=1.d-7 ! tolerance for Runge-Kutta propagator
        time_step=0.1d0 ! classical time step, fs
        n_class_steps=1 ! number of classical steps
        n_quant_steps=4 ! number of quantum steps for each classical step
        num_deriv_step=1.d-3 ! finite step for numerical derivatives, A
        therm_temperature=300.d0 ! Thermostat temperature, K
        therm_type=1 ! Thermostate type (0-no thermostat, 1-Langevin,2-Berendsen)
        berendsen_relax_const=0.4d0 ! Berendsen bath relaxation constant
        heating=0 ! 0-equilibrium dynamics, 1-heating
        heating_steps_per_degree=100 ! number of steps per degree during heating
        out_data_steps=1 ! number of steps to write data
        out_coords_steps=10 ! number of steps to write coordinates (to restart)
        therm_friction=2.d1 ! friction coefficient, ps^{-1}
        rnd_seed=1 ! seed for random number generator
        out_data_cube=0 ! write view files to generate cubes, 0-no, 1-yes
        verbosity=-1 ! output verbosity (0-minimal,3-largest)
        moldyn_deriv_flag=1 ! derivatives (0-numer,1-analyt,2-fast GS analytical)
        quant_step_reduction_factor=0.025d0 ! ???
        decoher_e0=0.1d0 ! decoherence parameter E0
        decoher_c=0.1d0 ! decoherecne parameter C
        decoher_type=2 ! Type of decoherence: Reinitialize (0) Never, ! (1) At successful hops, (2) At successful plus frustrated hops... ! (3) Persico/Granucci, or (4) Truhlar [2]
        dotrivial=1 !do trivial unavoided crossing routine (1) or not (0)  
        iredpot=0 !don't reduce the number of potentials during dynamics
        nstates=2 !do 2 states higher by default for iredpot
        deltared=1 !do 1 eV higher in energy for iredpot
 
        if(sim%Nsim.eq.1) then
              open (inputfdes,file='input.ceon',status='old')
        else
              write (filename, "(a6,i4.4,a5)") "input_", sim%id, ".ceon"
              open (inputfdes,file=trim(filename),status='old')
        endif 
               
        ! looking for namelist
        call nmlsrc('moldyn',inputfdes,ifind)
        if(ifind.ne.0) moldyn_found=.true.

        !Read qmmm namelist
        rewind inputfdes
        if(moldyn_found) then
            read(inputfdes,nml=moldyn)
        else
            write(6,*)'Could not find moldyn namelist'
            stop
        endif
        
        sim%naesmd%iredpot=iredpot
        sim%naesmd%nstates=nstates 
        sim%naesmd%deltared=deltared 
        sim%rk_comm%rk_tol=rk_tolerance 
        sim%dotrivial=dotrivial
        !Input checks for features under development
        if(decoher_type>2) then
            write(6,*) "decoher_type under development. Do not use."
            stop
        endif
        
        if(therm_type==2) then
            write(6,*) "Berendsen thermostat under development. Do not use."
            stop
        endif
        
        if(heating>0) then
            write(6,*) "heating under development. Do not use."
        endif
        !End input checks
        
        ! Set quant_coeffs_reinit
        if(decoher_type==1) then
            quant_coeffs_reinit=1 ! reinit of quantum coeffs after a hop (1-yes,0-no)
        else if(decoher_type==2) then
            quant_coeffs_reinit=2 ! reinit after frustrated and actual hops
        else
            quant_coeffs_reinit=0
        endif
        
        !This happens if verbosity is not set durring read namelist
        if(verbosity==-1) then
            if(n_class_steps>1) then
                verbosity=2
            else
                verbosity=3
            endif
        endif

        ! Tentatively, loading new parameter from the moldyn namelist to
        ! old parameters
        sim%ibo=bo_dynamics_flag
   
        sim%md%imdtype=exc_state_init
        sim%naesmd%ihop=sim%md%imdtype
        if(sim%md%imdtype==0) then
            sim%naesmd%state='fund'
        else
            sim%naesmd%state='exct'
        end if

        sim%naesmd%npot=n_exc_states_propagate
   
        sim%naesmd%icontini=out_count_init
        sim%naesmd%tfemto=time_init

        sim%naesmd%dtmdqt=time_step
       ! h=sim%naesmd%dtmdqt
        sim%naesmd%dtmdqt=sim%naesmd%dtmdqt/convtf

        sim%naesmd%nstep=n_class_steps

        sim%naesmd%nquantumreal=n_quant_steps
        sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumreal)

        sim%naesmd%decorhop=quant_coeffs_reinit
        sim%qmmm%numder_step=num_deriv_step
        sim%qmmm%ideriv=moldyn_deriv_flag

        sim%naesmd%temp0=therm_temperature
        sim%md%ttt=sim%naesmd%temp0

        therm_type=therm_type
        if(therm_type==0) sim%naesmd%ensemble ='energy'
        if(therm_type==1) sim%naesmd%ensemble ='langev'
        if(therm_type==2) sim%naesmd%ensemble ='temper'

        sim%naesmd%tao=berendsen_relax_const

        if(heating==1) then
            sim%naesmd%prep='heat'
        else
            sim%naesmd%prep='equi'
        end if

        sim%naesmd%istepheat=heating_steps_per_degree
        sim%naesmd%nstepw=out_data_steps
        sim%naesmd%nstepcoord=out_coords_steps
   
        sim%naesmd%friction=therm_friction/1.d3*convtf !convter therm_friction (in 1/ps) to sim%naesmd%friction (1/AU)

        sim%naesmd%iseedmdqt=rnd_seed
        write(6,*)"rnd_seed=",rnd_seed
        sim%naesmd%iview=out_data_cube
        sim%lprint=verbosity
        sim%md%ideriv=moldyn_deriv_flag

        sim%naesmd%nstepcross=int(1.d0/quant_step_reduction_factor)
        sim%constcoherE0=decoher_e0
        sim%constcoherC=decoher_c

        sim%cohertype=decoher_type

        sim%naesmd%deltared=sim%naesmd%deltared/feVmdqt !for reducing number of potentials   

        write(6,*) '!!!!!!-----MD INPUT-----!!!!!!'
        write(6,*)

        if (sim%md%imdtype.eq.0) then
            sim%qmmm%state_of_interest=sim%md%imdtype
            write(6,*)'Ground state MD,      imdtype=',sim%md%imdtype
    
        else
            sim%qmmm%state_of_interest=sim%md%imdtype
            write(6,*)'Excited state MD,      state=',sim%md%imdtype
            write(6,*)'Number of states to propagate',sim%naesmd%npot
        end if

        write(6,*)'Initial count for out files   ',sim%naesmd%icontini

        if(sim%naesmd%ensemble.eq.'langev') then
            write(6,*)'MD using Langevin '
            !write(6,*)'with friction coefficient [1/ps]=',friction*1.0D3

        else
            if(sim%naesmd%ensemble.eq.'temper') then
                write(6,*)'MD using Berendsen thermostat'
                write(6,*)'with bath relaxation constant [ps]=',sim%naesmd%tao
            else
                write(6,*)'MD at constant energy'
            end if
        end if

        write(6,*)'Starting time, fs        ',sim%naesmd%tfemto
        write(6,*)'Time step, fs		      ',sim%naesmd%dtmdqt*convtf
        write(6,*)'Number of classical steps (CS) ',sim%naesmd%nstep
        write(6,*)'Quantum steps/per CS         ',sim%naesmd%nquantumreal
        !	write(6,*)	'Vibr. modes to consider      ',   symm
        write(6,*)'Displacement for deriv. [A]   ',sim%qmmm%numder_step

        if(therm_type.eq.0) then
            write(6,*)'Newtonian dynamics,           therm_type=',therm_type
        else if(therm_type.eq.1) then
            write(6,*)'Langevin thermostat dynamics, therm_type=',therm_type
        else if(therm_type.eq.2) then
            write(6,*)'Temperature thermostat,       therm_type=',therm_type
        end if

        write(6,*)'Temperature, K               ',sim%naesmd%temp0
        write(6,*)'Bath relaxation constant     ',sim%naesmd%tao

        if(sim%naesmd%prep.eq.'equi') write(6,*)'Equilibrated dynamics'
        if(sim%naesmd%prep.eq.'heat') write(6,*)'Heating dynamics from T=0'

        write(6,*)'Number of steps per degree K  ',sim%naesmd%istepheat
        write(6,*)'Number of steps to write data ',sim%naesmd%nstepw
        write(6,*)'Number of steps to write coords',sim%naesmd%nstepcoord
        write(6,*)'Friction coefficient [1/ps]   ',sim%naesmd%friction/convtf*1d3 !convert back to ps
        write(6,*)'Seed for random generator    ',sim%naesmd%iseedmdqt

        if(sim%naesmd%iview.eq.1) write(6,*)'Will write files to generate cubes'

        write(6,*)'NAESMD verbosity, ',sim%lprint

        if(sim%md%ideriv.eq.0) then
            write(6,*)'MD will use numerical deriv, ideriv=',sim%md%ideriv
        else if(sim%md%ideriv.ge.2) then
            write(6,*)'MD will use fast GS analytic deriv,  sim%md%ideriv=',sim%md%ideriv
        else
            write(6,*)'MD will use analytic deriv,  ideriv=',sim%md%ideriv
        end if

        write(6,*) 'Cartesian coordinates are used'


        !***************************************************
        ! Initial allocations
        ! **************************************************
        call allocate_naesmd_module_init(sim%naesmd,natoms)
        call allocate_md_module_init(sim%md,natoms)
        allocate(xx(natoms),yy(natoms),zz(natoms))

        ! reading atoms and modes:
        rewind (inputfdes)
       Na=0
        ! read the cartesian coordinates
        !****************************************************************************
        do
            read(inputfdes,'(a)',err=29) txt
            if ( &
                txt(1:6).eq.'$COORD' &
                .or.txt(1:6).eq.'$coord' &
                .or.txt(1:6).eq.'&COORD' &
                .or.txt(1:6).eq.'&coord') exit ! breaking infinite loop
        end do

        read(inputfdes,'(a)',err=29) txt

        do while( &
            txt(1:9).ne.'$ENDCOORD' &
            .and.txt(1:9).ne.'$endcoord' &
            .and.txt(1:9).ne.'&ENDCOORD' &
            .and.txt(1:9).ne.'&endcoord')

            Na=Na+1 ! counting atoms

            read(txt,*,err=29) sim%md%atoms(Na),sim%md%r0(3*Na-2),sim%md%r0(3*Na-1),sim%md%r0(3*Na)
            sim%md%atmass(Na)=2*sim%md%atoms(Na)
            sim%naesmd%atomtype(Na)=sim%md%atoms(Na)
            !         if(sim%naesmd%atomtype(Na).eq.1) sim%naesmd%massmdqt(Na)=1.0079d0*convm
            if(sim%naesmd%atomtype(Na).eq.1) sim%naesmd%massmdqt(Na)=1.00d0*convm
            if(sim%naesmd%atomtype(Na).eq.1) sim%naesmd%atomtype2(Na)='H '
            if(sim%naesmd%atomtype(Na).eq.2) sim%naesmd%massmdqt(Na)=4.0026d0*convm
            if(sim%naesmd%atomtype(Na).eq.2) sim%naesmd%atomtype2(Na)='He'
            if(sim%naesmd%atomtype(Na).eq.3) sim%naesmd%massmdqt(Na)=6.941d0*convm
            if(sim%naesmd%atomtype(Na).eq.3) sim%naesmd%atomtype2(Na)='Li'
            if(sim%naesmd%atomtype(Na).eq.4) sim%naesmd%massmdqt(Na)=9.0122d0*convm
            if(sim%naesmd%atomtype(Na).eq.4) sim%naesmd%atomtype2(Na)='Be'
            if(sim%naesmd%atomtype(Na).eq.5) sim%naesmd%massmdqt(Na)=10.811d0*convm
            if(sim%naesmd%atomtype(Na).eq.5) sim%naesmd%atomtype2(Na)='B '
            if(sim%naesmd%atomtype(Na).eq.6) sim%naesmd%massmdqt(Na)=12.011d0*convm
            if(sim%naesmd%atomtype(Na).eq.6) sim%naesmd%atomtype2(Na)='C '
            if(sim%naesmd%atomtype(Na).eq.7) sim%naesmd%massmdqt(Na)=14.007d0*convm
            if(sim%naesmd%atomtype(Na).eq.7) sim%naesmd%atomtype2(Na)='N '
            if(sim%naesmd%atomtype(Na).eq.8) sim%naesmd%massmdqt(Na)=15.9994d0*convm
            if(sim%naesmd%atomtype(Na).eq.8) sim%naesmd%atomtype2(Na)='O '
            if(sim%naesmd%atomtype(Na).eq.9) sim%naesmd%massmdqt(Na)=18.998d0*convm
            if(sim%naesmd%atomtype(Na).eq.9) sim%naesmd%atomtype2(Na)='F '
            if(sim%naesmd%atomtype(Na).eq.10) sim%naesmd%massmdqt(Na)=20.180d0*convm
            if(sim%naesmd%atomtype(Na).eq.10) sim%naesmd%atomtype2(Na)='Ne'
            if(sim%naesmd%atomtype(Na).eq.14) sim%naesmd%massmdqt(Na)=28.086d0*convm
            if(sim%naesmd%atomtype(Na).eq.14) sim%naesmd%atomtype2(Na)='Si'
            if(sim%naesmd%atomtype(Na).eq.15) sim%naesmd%massmdqt(Na)=30.974d0*convm
            if(sim%naesmd%atomtype(Na).eq.15) sim%naesmd%atomtype2(Na)='P '
            if(sim%naesmd%atomtype(Na).eq.16) sim%naesmd%massmdqt(Na)=32.065d0*convm
            if(sim%naesmd%atomtype(Na).eq.16) sim%naesmd%atomtype2(Na)='S '
            if(sim%naesmd%atomtype(Na).eq.17) sim%naesmd%massmdqt(Na)=35.453d0*convm
            if(sim%naesmd%atomtype(Na).eq.17) sim%naesmd%atomtype2(Na)='Cl'

            read(inputfdes,'(a)',err=29) txt
        end do
        !--------------------------------------------------------------------
        !Allocate variables
        write(6,*)'Allocating variables'
        sim%naesmd%natom=Na
        sim%Na=Na
        allocate(sim%coords(Na*3))
        sim%excN=sim%naesmd%npot
        call allocate_naesmd_module(sim%naesmd,Na,sim%naesmd%npot)
        call allocate_md_module(sim%md,Na)
        if(sim%naesmd%npot.gt.0) then
                allocate(sim%rk_comm%thresholds(2*sim%excN))
                allocate(sim%naesmd%yg(2*sim%excN))
                allocate(sim%naesmd%ygprime(2*sim%excN))
                allocate(sim%naesmd%cross(sim%excN))
                allocate(sim%naesmd%Omega(sim%excN))
        endif

        !--------------------------------------------------------------------
        !
        !  Reading velocities
        !
        !--------------------------------------------------------------------

        rewind (inputfdes)
        do
            read(inputfdes,'(a)',err=29) txt
            if( &
                txt(1:6).eq.'$VELOC' &
                .or.txt(1:6).eq.'$veloc' &
                .or.txt(1:6).eq.'&VELOC' &
                .or.txt(1:6).eq.'&veloc') exit ! exiting infinite loop
        end do

        i=1
        do
            read(inputfdes,'(a)',err=29) txt
            if( &
                txt(1:9).eq.'$ENDVELOC' &
                .or.txt(1:9).eq.'$endveloc' &
                .or.txt(1:9).eq.'&ENDVELOC' &
                .or.txt(1:9).eq.'&endveloc') exit

            read(txt,*,err=29) sim%naesmd%vx(i),sim%naesmd%vy(i),sim%naesmd%vz(i)
            i=i+1
        end do

        sim%naesmd%vx(1:sim%naesmd%natom)=sim%naesmd%vx(1:sim%naesmd%natom)/convl*convt
        sim%naesmd%vy(1:sim%naesmd%natom)=sim%naesmd%vy(1:sim%naesmd%natom)/convl*convt
        sim%naesmd%vz(1:sim%naesmd%natom)=sim%naesmd%vz(1:sim%naesmd%natom)/convl*convt
        !
        !--------------------------------------------------------------------
        ! Reading quantum coeffients
        !
        !--------------------------------------------------------------------
        !
        rewind (inputfdes)
        do
            read(inputfdes,'(a)',err=29) txt
            if( &
                txt(1:6).eq.'$COEFF' &
                .or.txt(1:6).eq.'$coeff' &
                .or.txt(1:6).eq.'&COEFF' &
                .or.txt(1:6).eq.'&coeff') exit ! leaving the infinite loop
        end do

        i=1

        if(sim%excN.gt.0) then
          do
            read(inputfdes,'(a)',err=29) txt
            if( &
                txt(1:9).eq.'$ENDCOEFF' &
                .or.txt(1:9).eq.'$endcoeff' &
                .or.txt(1:9).eq.'&ENDCOEFF' &
                .or.txt(1:9).eq.'&endcoeff') exit

            read(txt,*,err=29) sim%naesmd%yg(i),sim%naesmd%yg(i+sim%excN)

            sim%naesmd%yg(i)=dsqrt(sim%naesmd%yg(i))
            sim%naesmd%yg(i+sim%excN)=dasin(sim%naesmd%yg(i+sim%excN))

            i=i+1
            if(i.gt.sim%excN) exit ! too many coefficients - leaving the loop
          end do
        endif

        close (inputfdes)
        !
        !--------------------------------------------------------------------
        !
        write(6,*)' Finish reading coordinates '
        !
        !--------------------------------------------------------------------
        !
        Nm=3*Na
        sim%md%v(1:Nm,1:Nm)=0.d0
        do i=1,Na
            if (sim%md%atoms(i).eq.1) then
                sim%md%fo(3*i-2)=3000
                sim%md%fo(3*i-1)=3000
                sim%md%fo(3*i)=3000
            else
                sim%md%fo(3*i-2)=1000
                sim%md%fo(3*i-1)=1000
                sim%md%fo(3*i)=1000
            end if
   
            sim%md%fm(3*i-2)=sim%md%atmass(i)
            sim%md%fm(3*i-1)=sim%md%atmass(i)
            sim%md%fm(3*i)=sim%md%atmass(i)
            sim%md%v(3*i-2,3*i-2)=1.0
            sim%md%v(3*i-1,3*i-1)=1.0
            sim%md%v(3*i,3*i)=1.0
        end do

        ! initialize variables:
        N3=3*Na
        N2=2*Nm
        N1=Nm+1

        ! compute initial cartesian coordinates:
        do j=3,N3,3
            xx(j/3)=sim%md%r0(j-2)
            yy(j/3)=sim%md%r0(j-1)
            zz(j/3)=sim%md%r0(j)
        end do

        ! define cartesian coordinates for mdqt
        ! transform coordiantes from angstrom to atomic units
        sim%naesmd%rx(1:sim%naesmd%natom)=xx(1:sim%naesmd%natom)/convl
        sim%naesmd%ry(1:sim%naesmd%natom)=yy(1:sim%naesmd%natom)/convl
        sim%naesmd%rz(1:sim%naesmd%natom)=zz(1:sim%naesmd%natom)/convl

        ! position of the center of mass
        sim%naesmd%xcmini=0.d0
        sim%naesmd%ycmini=0.d0
        sim%naesmd%zcmini=0.d0

        sim%naesmd%masstot=0.d0 ! kav: - zeroing was not here, needed?
        do j=1,sim%naesmd%natom
            sim%naesmd%masstot=sim%naesmd%masstot+sim%naesmd%massmdqt(j)
        end do
        do j=1,sim%naesmd%natom
            sim%naesmd%xcmini=sim%naesmd%xcmini+sim%naesmd%rx(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            sim%naesmd%ycmini=sim%naesmd%ycmini+sim%naesmd%ry(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            sim%naesmd%zcmini=sim%naesmd%zcmini+sim%naesmd%rz(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
        end do
   
        !Remove rotation and translation from initial velocity
        write(6,*)'Rescaling velocity'
        call rescaleveloc(sim%naesmd%rx,sim%naesmd%ry,sim%naesmd%rz, &
		sim%naesmd%vx,sim%naesmd%vy,sim%naesmd%vz,sim%naesmd%massmdqt,sim%naesmd%natom)

        ! compute kinetic energy, for cartesian option
        sim%naesmd%kin=0.d0
        do i=1,sim%naesmd%natom
            sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(i)* &
		(sim%naesmd%vx(i)**2+sim%naesmd%vy(i)**2+sim%naesmd%vz(i)**2)/2
        end do

        write(6,*)
        write(6,*)'Initial Geometry'

        do i=1,Na
            write(6,100) sim%md%atoms(i),xx(i),yy(i),zz(i)
        end do
100     format(I5,'     ',3F12.6)


!        if(sim%naesmd%tfemto.eq.0.d0) then
!            open (9,file='coords.xyz')
!            write (9,*) Na
!            write (9,*) 'Input Geometry'
!
!            do j = 1,Na
!                write(9,302) ELEMNT(sim%md%atoms(j)),xx(j),yy(j),zz(j)
!            end do
!            close (9)
!        end if
        sim%coords(1:3*Na)=sim%md%r0(1:3*Na)
        allocate(sim%deriv_forces(3*Na))

        return

7       format (' ',A12,' ',A35,A2)
8       format (' ',A27,'  ',A30,A2)
9       format (' ',A20,'    ',g16.5,A10)
301     format(' |     ',i2,' days ',i2,' hours ',i2,' minutes ',i2, &
            ' seconds     |')
302     format(A3,3f12.6)
303     format(i6,f9.3,g11.3,5g16.8)
304     format('STEP:',i5,f9.3,6g13.5,g10.2,f16.10)
305     format(a,i5,6g16.8,g12.4)
306     format(i7,10000g12.4)
307     format(a,g8.4,6g16.8,g12.4)
999     FORMAT(I3,1X,1000(1X,F18.10))
29      write(6,*) txt
98      stop 'bad input'
    end subroutine





  
end
