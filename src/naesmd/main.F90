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
    use cosmo_C, only : cosmo_C_structure
    use clone_module
    use AIMC_type_module
    implicit none
    !
    !--------------------------------------------------------------------
    !
    !  Molecular dynamics uses cartesian coordinates, analytical
    !  or numerical derivatives, and
    !  Newtonian (usual velocity-Verlet), Langevin thermostat velocity-Verlet,
    !  or Berendsen thermostat (all molecules collide with heat bath at
    !  once) algorithms.
    !
    !  input.ceon file with CEO and MDQT parameters,
    !  as well as initial coordinates, velocities, and quantum coefficients;
    !
    !  reastart.out instead for restarting
    !
    !--------------------------------------------------------------------

    integer :: Nsim = 1, Nsim_max=100
    integer :: restart_flag

    call calculate_Nsim(Nsim,restart_flag)
    call simulations(Nsim,Nsim_max,restart_flag)

contains

subroutine simulations(Nsim,Nsim_max,restart_flag)
    use AIMC, only: AIMC_clone_check, AIMC_clone
    use dropout_module
    use naesmd_constants
    integer, intent(inout) :: Nsim
    integer, intent(in) :: Nsim_max
    integer, intent(in) :: restart_flag
    integer :: i,imdqt,j,k,l
    _REAL_ :: tini, tfin !For integrating MCE equations
    _REAL_ :: time
    integer,dimension(:) :: drops(Nsim_max) !For stroing dropout indexes

    !The instances of all sim structures
    type(qm2_davidson_structure_type),target :: qm2ds_notp(Nsim_max)
    type(qmmm_struct_type), target :: qmmm_struct_notp(Nsim_max)
    type(qm2_structure), target :: qm2_struct_notp(Nsim_max)
    type(cosmo_C_structure), target :: cosmo_C_struct_notp(Nsim_max)
    type(naesmd_structure), target :: naesmd_struct_notp(Nsim_max)
    type(xlbomd_structure), target :: xlbomd_struct_notp(Nsim_max)
    type(md_structure), target :: md_struct_notp(Nsim_max)
    type(qm2_params_type),target :: qparams_struct_notp(Nsim_max)
    type(qmmm_nml_type),target :: qnml_struct_notp(Nsim_max)
    type(qm2_rij_eqns_structure),target:: rij_struct_notp(Nsim_max)
    type(qm_gb_structure),target :: gb_struct_notp(Nsim_max)
    type(qmmm_mpi_structure),target :: qmpi_struct_notp(Nsim_max)
    type(qmmm_opnq_structure),target :: opnq_struct_notp(Nsim_max)
    type(qmmm_div_structure), target :: div_struct_notp(Nsim_max)
    type(qmmm_vsolv_type), target :: vsolv_struct_notp(Nsim_max)
    type(qmmm_scratch_structure), target:: scratch_struct_notp(Nsim_max)
    type(qm_ewald_structure), target :: ewald_struct_notp(Nsim_max)
    type(AIMC_type), target :: aimc_struct_notp(Nsim_max)

    !The sim structure, a set of pointers and some values
    type(simulation_t),target::sim_notp(Nsim_max)
    !The pointer to sim structure, more rapidly passed into subroutines
    type(sim_pointer_array) ::sims(Nsim_max)

    !The MCE structure for solving the nuclear TDS
    type(MCE) :: nuclear
    logical file_exists2
    integer dps
    _REAL_ junk_time
    integer, allocatable :: indexes(:)

    dps = 0
    inquire(FILE="dropped.out", EXIST=file_exists2)
    if(file_exists2) then
        call system('wc dropped.out > dropped.tmp')
        open(1,file='dropped.tmp')
            read(1,*) dps
        close(1)
        call system('rm dropped.tmp')
        if(dps.gt.0) then
            allocate(indexes(dps))
            open(1,file='dropped.out')
                do i = 1, dps
                    read(1,*) junk_time, indexes(i)
                enddo
            close(1)
        else
            allocate(indexes(1))
            indexes(1) = -1
        endif
    else
        allocate(indexes(1))
        indexes(1) = -1
    endif
    Nsim = Nsim - dps

    drops = 0
    do i = 0, dps - 1
        if (any(indexes == i)) drops(i + 1) = 1
    enddo

    do i=1, Nsim_max
            sims(i)%sim=>sim_notp(i)
            sims(i)%sim%id=i-1
            call init0_simulation(sims(i)%sim,Nsim)
            call setp_simulation(sims(i)%sim,qmmm_struct_notp(i),qm2ds_notp(i), &
                        qm2_struct_notp(i),naesmd_struct_notp(i),md_struct_notp(i), &
                        cosmo_c_struct_notp(i),xlbomd_struct_notp(i),qparams_struct_notp(i), &
                        qnml_struct_notp(i),rij_struct_notp(i),gb_struct_notp(i),qmpi_struct_notp(i), &
                        opnq_struct_notp(i),div_struct_notp(i),vsolv_struct_notp(i),scratch_struct_notp(i), &
                        ewald_struct_notp(i),aimc_struct_notp(i))
    enddo
    
    j=0
    do i=1, Nsim + dps
        if(.not.any(sims(i)%sim%id == indexes)) then
            j = j + 1
            call initiate_sim(sims(j)%sim,Nsim,Nsim_max,restart_flag,nuclear,i-1)
            if(sims(j)%sim%naesmd%tfemto.gt.0.0d0.and.(sims(j)%sim%naesmd%dynam_type.eq.'aimc'.or.sims(j)%sim%naesmd%dynam_type.eq.'mf')) then
                call check_sign_for_restart(sims(j)%sim)
            endif
        endif
    enddo 
    if(sims(1)%sim%naesmd%dynam_type.eq.'aimc'.and.restart_flag.eq.1) then
        open(1,file='Pha_last.out')
            read(1,999) j, time, (sims(i)%sim%naesmd%pha,i=1,Nsim)
        close(1)
    endif
    
!Parallel Trajectories (need to make this makefile option) 
     do imdqt=1,sims(Nsim)%sim%naesmd%nstep !Main loop (ONLY LAST INPUT FILE NSTEP is used)
         if (sims(1)%sim%naesmd%dynam_type.eq.'aimc'.and.sims(1)%sim%naesmd%tfemto.eq.0.0d0) then
            call write_mce_data_0(sims, nuclear, Nsim)
         endif
         do i=1,Nsim
            call nexmd_sim_step(sims(i)%sim,imdqt,nuclear,Nsim)
         enddo
!Dropping trajectories at S0/S1 conincal intersection
         do i = 1, Nsim
            if (sims(i)%sim%excN.gt.0.and.(sims(i)%sim%naesmd%vmdqt(1)-sims(i)%sim%naesmd%vgs)*feVmdqt.lt.sims(i)%sim%naesmd%S0_S1_threshold) then
                call dropout(sims,nuclear,i,Nsim,drops,Nsim_max)
            endif
         enddo
!Dropping trajectories at S0/S1 conincal intersection
         if(sims(1)%sim%naesmd%dynam_type.eq.'aimc') then
            call propagate_sE(sims,nuclear,Nsim)
            if(sims(1)%sim%naesmd%icontw.ne.sims(1)%sim%naesmd%nstepw) then
                call update_nuclear(nuclear,sims,Nsim,0,0)
            else
                if(sims(i)%sim%naesmd%icont.ne.sims(i)%sim%naesmd%nstepcoord) then
                    call update_nuclear(nuclear,sims,Nsim,1,0)
                else
                    call update_nuclear(nuclear,sims,Nsim,1,1)
                endif
            endif
         endif
         do i=1, Nsim
            if(sims(i)%sim%naesmd%dynam_type.eq.'aimc') then
                if( AIMC_clone_check(sims(i)%sim).and.(Nsim.lt.Nsim_max) ) then
                        Nsim=Nsim+1
                        sims(Nsim)%sim%id=Nsim-1
                        sims(i)%sim%aimc%nclones=sims(i)%sim%aimc%nclones+1
                        call clone_sim(sims(i)%sim,sims(Nsim)%sim) !exact copy of sim
                        call AIMC_clone(sims(i)%sim, sims(Nsim)%sim, nuclear, Nsim, i) !AIMC cloning procedure
                        sims(Nsim)%sim%aimc%new_clone=.true.
                        call open_output_multi(sims(Nsim)%sim,sims(Nsim)%sim%ibo,&
                                sims(Nsim)%sim%naesmd%tfemto,sims(Nsim)%sim%md%imdtype,sims(Nsim)%sim%lprint,sum(drops))
                endif
            endif
            if(sims(i)%sim%naesmd%icontw.ne.sims(i)%sim%naesmd%nstepw) then 
                continue
            else 
                call writeoutput(sims(i)%sim,sims(i)%sim%ibo,sims(i)%sim%naesmd%yg,sims(i)%sim%naesmd%yg_new,sims(i)%sim%lprint,sims(i)%sim%naesmd%cross,Nsim)
            end if
         enddo
         do i=1,Nsim
             if(sims(i)%sim%naesmd%icontw.ne.sims(i)%sim%naesmd%nstepw) then
                 sims(i)%sim%naesmd%icontw=sims(i)%sim%naesmd%icontw+1
             else
                 sims(i)%sim%naesmd%icontw=1
                 if(sims(i)%sim%naesmd%icont.ne.sims(i)%sim%naesmd%nstepcoord) then
                     sims(i)%sim%naesmd%icont=sims(i)%sim%naesmd%icont+1
                 else
                     call write_restart(sims(i)%sim, Nsim, drops, Nsim_max) !Restart written here to achieve syncronization
                     print *, 'Restart written for tr ', sims(i)%sim%id, ' at ', sims(i)%sim%naesmd%tfemto, 'fs'
                     sims(i)%sim%naesmd%icont=1
                     sims(i)%sim%naesmd%icontpdb=sims(i)%sim%naesmd%icontpdb+1
                 endif
             endif
        enddo
    enddo
    do i=1, Nsim
        call finish_sim(sims(i)%sim,Nsim)
    enddo 
    stop

999     FORMAT(I3,1000(1X,F18.10))

end subroutine simulations


subroutine initiate_sim(sim,Nsim,Nsim_max,restart_flag,nuclear,tr_number)
    type(simulation_t), pointer::sim
    type(MCE) :: nuclear
    character*30 datetime, machname*36
    integer inputfdes
    integer i,j,k,l
    integer, intent(inout) :: Nsim
    integer, intent(in) :: Nsim_max
    integer, intent(in) :: restart_flag
    integer nstep0
    _REAL_ t_start,t_finish 
    character(100) ::  filename
    integer drops
    integer, intent(in) :: tr_number

    call get_date(datetime)
    call get_machine(machname)
    write (6,*)
    write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
    write (6,8) '| MD execution started at  ',datetime,' |'
    write (6,7) '| Computer: ',   machname,'|'
    write (6,*) '|________________________________________________|'
    write (6,*)
    call cpu_time(sim%time1) 

    filename = ''
    
    if(restart_flag.eq.1) then
        call check_files(Nsim,nstep0,tr_number,sim%id)
        write(6,*) 'All files checked and ready for restart for trajectory: ', tr_number
    endif
    call init_main(sim, inputfdes, Nsim, Nsim_max, restart_flag, nuclear, tr_number)
    if(restart_flag.eq.1) sim%naesmd%nstep0=nstep0

    !Put derivative variables into module

    !***********************************************************
    !Initialization and single point calculation or zeroeth
    !step of adiabatic or nonadiabatic dynamics
    !
    !Open input file, initialize variables, run first calculation
    !***********************************************************
    if(Nsim.eq.1.and.restart_flag.eq.0) then
        open (inputfdes,file='input.ceon',status='old')
    elseif(Nsim.eq.1.and.restart_flag.eq.1) then
        open(inputfdes,file='restart.out',status='old')
    else
        write (filename, "(a8,i4.4,a4)") "restart_", sim%id, ".out"
        open (inputfdes,file=trim(adjustl(filename)),status='old')
    endif 

    !Open output files
    if((Nsim.eq.1).and.sim%naesmd%dynam_type.ne.'aimc') then
            call open_output(sim,sim%ibo,sim%naesmd%tfemto,sim%md%imdtype,sim%lprint)
            sim%naesmd%ihop=sim%qmmm%state_of_interest 
    else
            drops = 0
            call open_output_multi(sim,sim%ibo,sim%naesmd%tfemto,sim%md%imdtype,sim%lprint,tr_number-sim%id)
            sim%naesmd%ihop=sim%qmmm%state_of_interest 
    endif

    call init_sqm(sim,inputfdes,STDOUT) ! involves call to Davidson
    close(inputfdes)

    if((((sim%naesmd%dynam_type.eq.'mf').or.(sim%naesmd%dynam_type.eq.'aimc'))).and.(sim%cosmo%solvent_model>0)) then
        write(6,*) 'Mean-field not yet compatible with solvent models'
        stop
    endif    

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
        sim%naesmd%torden(i)=i
        sim%naesmd%iordenhop(i)=0
    end do
     
    sim%naesmd%vgs=sim%naesmd%E0;
    sim%naesmd%vmdqt(1:sim%excN)=sim%naesmd%Omega(1:sim%excN)+sim%naesmd%vgs


    !Calculate derivatives
    call cpu_time(t_start)
    if (sim%qmmm%ideriv.gt.0) then
        if(((sim%naesmd%dynam_type.eq.'mf').or.(sim%naesmd%dynam_type.eq.'aimc'))) then
                call deriv_MF(sim,restart_flag)
        else
                call deriv(sim,sim%naesmd%ihop)
        endif
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

    if (sim%naesmd%tfemto==0.0d0) call writeoutputini(sim,sim%ibo,sim%naesmd%yg,sim%naesmd%yg_new,sim%lprint)

7   format (' ',A12,' ',A35,A2)
8   format (' ',A27,'  ',A30,A2)

    call flush(6)
    return

end subroutine initiate_sim

subroutine nexmd_sim_step(sim, imdqt, nuclear, Nsim)
    type(simulation_t),pointer::sim
    integer, intent(in) :: imdqt
    integer ii
    type(MCE) :: nuclear
    integer, intent(in) :: Nsim
    
        !Classical propagation step - BOMD or NAESMD
        write(6,*)"Begin classical propagation step #",imdqt, " sim #: ", sim%id
        sim%naesmd%tfemto=sim%naesmd%tfemto+sim%naesmd%dtmdqt*convtf
        sim%qmmm%num_qmmm_calls=imdqt !necessary for BO dynamics

        if(sim%naesmd%state.eq.'exct'.and.sim%ibo.ne.1) then
            !NAESMD-propagate=====================================================================================
            sim%qmmm%state_of_interest=sim%naesmd%ihop;
            call verlet1(sim)
            write(6,*)'Begin nonadiabatic couplings and crossings calculations'
            call initialize(sim,sim%naesmd%yg_new)
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
            write(6,*)"End classical propagation step #",imdqt, " sim #: ", sim%id
            !--------------------------------------------------------------------
            ! Loop for quantum propagation steps
            ! that implies CEO energy calculations
            !--------------------------------------------------------------------
            call quantum_propagation(sim, imdqt, nuclear, Nsim)         
            !--------------------------------------------------------------------
            ! last part of velocity verlet algorithm
            do ii=1,sim%excN
                sim%naesmd%yg(ii)=sqrt(sim%naesmd%yg_new(ii)**2+sim%naesmd%yg_new(ii+sim%excN)**2)
                sim%naesmd%yg(ii+sim%excN)=atan2(sim%naesmd%yg_new(ii+sim%excN),sim%naesmd%yg_new(ii))+sim%naesmd%yg_new(ii+2*sim%excN)
            enddo
            call verlet2(sim)
            if(sim%naesmd%dynam_type.eq.'tsh') then
                    !--------------------------------------------------------------------
                    ! analyze the hopping
                    !--------------------------------------------------------------------
                    call evalhop(sim,sim%lprint,sim%Na,sim%naesmd%yg,sim%naesmd%yg_new,sim%naesmd%cross)
                    !--------------------------------------------------------------------

                     call decoherence_E0_and_C(sim)!currently this is not available
            else if (sim%dotrivial.eq.1) then
                     call just_trivial(sim,sim%naesmd%yg,sim%naesmd%yg_new,sim%naesmd%cross,nuclear,Nsim)
            endif 
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


    !***********************************
    !      ttime=dtime(tarray)
    !      write(19,*) tarray(1)

    return
end subroutine nexmd_sim_step

subroutine finish_sim(sim,Nsim)
    use writeoutput_module

    type(simulation_t),pointer::sim

    integer itime11,itime2,itime3
    integer, intent(in) :: Nsim
    _REAL_ time11,time12
    character*30 datetime

    !##########################################################
    !! END of the Molecular dynamics with Quantum transition main loop
    !##########################################################

    call cpu_time(time11)
    time11=time11-sim%time1
    itime11=time11
    sim%itime1=MOD(itime11,60)
    itime11=itime11/60
    itime2=MOD(itime11,60)
    itime11=itime11/60
    itime3=MOD(itime11,24)
    itime11=itime11/24

    call get_date(datetime)

    write (6,*)
    write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
    write (6,8) '| MD normal termination at ', datetime,' |'
    write (6,9) '| MD total CPU time            ',  time11, ' seconds |'
    write(6,301) itime11,itime3,itime2,sim%itime1
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


8   format (' ',A27,'  ',A30,A2)
9   format (' ',A20,'    ',g16.5,A10)
301 format(' |     ',i2,' days ',i2,' hours ',i2,' minutes ',i2, &
        ' seconds     |')

    return
end subroutine finish_sim

    !
    !********************************************************************
    !
    subroutine init_main(sim, inputfdes, Nsim, Nsim_max, restart_flag, nuclear, tr_number)
        use naesmd_constants
        use communism
        implicit none
        type(simulation_t),pointer::sim
        type(MCE) :: nuclear
        integer, intent(in) :: Nsim
        integer, intent(in) :: Nsim_max
        integer, intent(in) :: restart_flag
        _REAL_,allocatable::xx(:),yy(:),zz(:)
        integer, intent (inout) :: inputfdes
        integer, intent(in) :: tr_number
        integer ifind
        integer :: i,j,k,l,m,ii,jjj,ios
        integer :: Na, Nm, N1, N2, N3
        integer slen
        integer :: itime1
        character*(150) txt
        character(100) ::  filename
        character(100) :: NAMD_type
        _REAL_ :: AIMC_dclone_1,AIMC_dclone_2,AIMC_dclone_3,AIMC_max_clone,AIMC_force_pop_min
        integer nclones0
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
            integer ifixed
            integer :: nmc, npc
            integer :: printTdipole
            integer :: printTDM
            _REAL_ :: S0_S1_threshold
            namelist /moldyn/ natoms,bo_dynamics_flag,NAMD_type,exc_state_init, &
                n_exc_states_propagate,out_count_init,time_init, &
                time_step,n_class_steps,n_quant_steps, &
                num_deriv_step, &
                therm_temperature,therm_type, &
                berendsen_relax_const,heating, &
                heating_steps_per_degree,out_data_steps,out_coords_steps, &
                therm_friction,rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag, &
                quant_step_reduction_factor,decoher_e0,decoher_c,decoher_type,dotrivial,&
                iredpot,nstates,deltared,ifixed,AIMC_dclone_1,AIMC_dclone_2,AIMC_dclone_3,&
                AIMC_max_clone,AIMC_force_pop_min,nmc,npc,printTdipole,printTDM,nclones0,S0_S1_threshold
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

        sim%naesmd%conthop=0
        sim%naesmd%conthop2=0
        !
        !--------------------------------------------------------------------
        !

        inputfdes=12*10000+sim%id
        filename = ''
        ! copy input file
        if(Nsim.eq.1.and.restart_flag.eq.0) then
            open(inputfdes,file='input.ceon',status='old')
        elseif(Nsim.eq.1.and.restart_flag.eq.1) then
            open(inputfdes,file='restart.out',status='old')
        else
            write (filename, "(a8,i4.4,a4)") "restart_", tr_number, ".out"
            open (inputfdes,file=filename(1:16),status='old')
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
        ifixed=0 ! do not use fixed atoms
        NAMD_type='tsh' ! 'tsh' for SH, 'mf' for EHR, 'aimc' for cloning
        AIMC_dclone_1=1.5
        AIMC_dclone_2=0.2617
        AIMC_dclone_3=0.005
        AIMC_max_clone=4
        AIMC_force_pop_min=0.0
        nclones0=0 !initial count for the number of clones
        nmc=0 ! number of nuclear directions to freeze
        npc=0 ! number of pairs of atoms to freeze the distance between them
        printTdipole=0 !to print transition dipole moments in separate files
        printTDM=0 !to print complete TDM
        S0_S1_threshold=0.2d0 !for dropping trajectory when S0/S1 gap is lower than "S0_S1_threshold" (eV)

        filename = ''
        if(Nsim.eq.1.and.restart_flag.eq.0) then
            open (inputfdes,file='input.ceon',status='old')
        elseif(Nsim.eq.1.and.restart_flag.eq.1) then
            open(inputfdes,file='restart.out',status='old')
        else
            write (filename, "(a8,i4.4,a4)") "restart_", tr_number, ".out"
            open (inputfdes,file=filename(1:16),status='old')
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

        if(ifixed.ne.0) then
          open(617,file='fixed_atoms')
          do i=1,ifixed
             read(617,*) sim%naesmd%ifxd(i) 
          enddo
          close(617)
        endif
        
        sim%naesmd%iredpot=iredpot
        sim%naesmd%nstates=nstates 
        sim%naesmd%deltared=deltared 
        sim%dotrivial=dotrivial
        !Input checks for features under development
        if(decoher_type>2) then
            write(6,*) "decoher_type under development. Do not use."
            stop
        endif
        call int_legal_range('moldyn: (decoher_type ) ', decoher_type,0,2) 
        call int_legal_range('moldyn: (trivial crossing ) ', dotrivial,0,1)
        call float_legal_range('moldyn: (quant_step_reduction_factor ) ',quant_step_reduction_factor,0.0d-16,1.0d1)
        call float_legal_range('moldyn: (therm_friction ) ',therm_friction,0.0d0,1.0d6)
        call float_legal_range('moldyn: (therm_temperature ) ',therm_temperature,0.0d0,1.0d7)
        call int_legal_range('moldyn: (therm_type ) ', therm_type,0,2)
        call int_legal_range('moldyn: (number of atoms ) ', natoms,0,999999999)
        call int_legal_range('moldyn: (bo_dynamics_flag ) ', bo_dynamics_flag,0,1)
        if (NAMD_type.eq.'tsh') call int_legal_range('moldyn: (exc_state_init ) ', exc_state_init,0,999999999)
        call int_legal_range('moldyn: (n_exc_states_propagate ) ', n_exc_states_propagate,0,999999999)

        call float_legal_range('moldyn: (initial time ) ',time_init,0.0d0,1.0d21)
        call float_legal_range('moldyn: (time step ) ',time_step,0.0d-6,1.0d21)

        call int_legal_range('moldyn: (classical steps ) ', n_class_steps,0,999999999)
        call int_legal_range('moldyn: (quantum steps ) ', n_quant_steps,1,999999999)

        call int_legal_range('moldyn: (moldyn_deriv_flag ) ', moldyn_deriv_flag,0,2)

        call float_legal_range('moldyn: (Displacement for numerical derivatives) ',num_deriv_step,0.0d-16,1.0d0)

        call int_legal_range('moldyn: (verbosity ) ', verbosity,0,3)

        call int_legal_range('moldyn: (Number of steps to write data  ) ', out_data_steps,0,999999999)
        call int_legal_range('moldyn: (Number of steps to write the restart file ) ', out_coords_steps,0,999999999)
        call int_legal_range('moldyn: (view files to generate cubes ) ', out_data_cube,0,1)
        
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
        sim%naesmd%dynam_type=NAMD_type
        if(n_exc_states_propagate.eq.0.and.NAMD_type.eq.'aimc') then
            write(6,*) "Can't run ground state dynamics (n_exc_states_propagate = 0) and AIMC (NAMD_type = 'aimc'). Fix and relaunch!"
            STOP;
        endif
        if((NAMD_type.eq.'mf').or.(NAMD_type.eq.'aimc')) sim%md%imdtype=-1
        if((sim%ibo.eq.0).AND.NAMD_type.ne.'tsh') then
            if((NAMD_type.ne.'mf').AND.(NAMD_type.ne.'aimc')) then
               write(6,*) 'only Surface hopping-"NAMD_type=tsh",Ehrenfest-"NAMD_type=mf"'
               write(6,*) 'and Ab-Initio Multiple Cloning- "NAMD_type=aimc" are allowed'
               stop
            endif
        endif
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
        sim%naesmd%dtmdqt=sim%naesmd%dtmdqt/convtf

        sim%naesmd%nstep=n_class_steps
        if(restart_flag.eq.0) sim%naesmd%nstep0=sim%naesmd%nstep
        if((NAMD_type.eq.'aimc').and.(1.ne.n_quant_steps)) then
             n_quant_steps=1
             write(6,*) 'AIMC requires n_quant_steps Equal 1, for now'
        endif

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
        sim%naesmd%fix = ifixed 

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

        sim%aimc%delta_clone_1=AIMC_dclone_1
        sim%aimc%delta_clone_2=AIMC_dclone_2
        sim%aimc%delta_clone_3=AIMC_dclone_3
        sim%aimc%force_pop_min=AIMC_force_pop_min
        sim%aimc%nclones_max=AIMC_max_clone
        sim%cohertype=decoher_type

        sim%naesmd%nmc=nmc
        sim%naesmd%npc=npc
        sim%naesmd%printTdipole=printTdipole
        sim%naesmd%printTDM=printTDM
        sim%naesmd%S0_S1_threshold=S0_S1_threshold
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

        else
            if(sim%naesmd%ensemble.eq.'temper') then
                write(6,*)'MD using Berendsen thermostat'
                write(6,*)'with bath relaxation constant [ps]=',sim%naesmd%tao
            else
                write(6,*)'MD at constant energy'
            end if
        end if

        if (sim%naesmd%printTdipole.gt.0) write(6,*) 'Printing transition dipole moments: ', sim%naesmd%printTdipole
        if (sim%naesmd%printTDM.gt.0) write(6,*) 'Printing the complete TDM', sim%naesmd%printTDM

        write(6,*)'Starting time, fs        ',sim%naesmd%tfemto
        write(6,*)'Time step, fs                      ',sim%naesmd%dtmdqt*convtf
        write(6,*)'Number of classical steps (CS) ',sim%naesmd%nstep
        write(6,*)'Quantum steps/per CS         ',sim%naesmd%nquantumreal
        write(6,*)'Displacement for deriv. [A]   ',sim%qmmm%numder_step

        if(therm_type.eq.0) then
            write(6,*)'Newtonian dynamics,           therm_type=',therm_type
        else if(therm_type.eq.1) then
            write(6,*)'Langevin thermostat dynamics, therm_type=',therm_type
        else if(therm_type.eq.2) then
            write(6,*)'Temperature thermostat,       therm_type=',therm_type
        end if

        if(nmc.gt.0.and.npc.gt.0) then
            write(6,*)'Fatal error: NEXMD can not freeze normal modes (nmc > 0) and distances between atoms (npc > 0) at the same time!'
            stop
        endif

        if(nmc.gt.0) then
            write(6,*)'Freezing a total of ', nmc, ' normal modes'
        else if(npc.gt.0) then
            write(6,*)'Freezing a total of ', npc, 'distances between atoms'
        else
            write(6,*)'MD without freezing'
        endif
        
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
            if(sim%naesmd%atomtype(Na).eq.18) sim%naesmd%massmdqt(Na)=39.948d0*convm
            if(sim%naesmd%atomtype(Na).eq.18) sim%naesmd%atomtype2(Na)='Ar'
            if(sim%naesmd%atomtype(Na).eq.19) sim%naesmd%massmdqt(Na)=39.098d0*convm
            if(sim%naesmd%atomtype(Na).eq.19) sim%naesmd%atomtype2(Na)='K'
            if(sim%naesmd%atomtype(Na).eq.20) sim%naesmd%massmdqt(Na)=40.078d0*convm
            if(sim%naesmd%atomtype(Na).eq.20) sim%naesmd%atomtype2(Na)='Ca'
            if(sim%naesmd%atomtype(Na).eq.21) sim%naesmd%massmdqt(Na)=44.956d0*convm
            if(sim%naesmd%atomtype(Na).eq.21) sim%naesmd%atomtype2(Na)='Sc'
            if(sim%naesmd%atomtype(Na).eq.22) sim%naesmd%massmdqt(Na)=47.867d0*convm
            if(sim%naesmd%atomtype(Na).eq.22) sim%naesmd%atomtype2(Na)='Ti'
            if(sim%naesmd%atomtype(Na).eq.23) sim%naesmd%massmdqt(Na)=50.942d0*convm
            if(sim%naesmd%atomtype(Na).eq.23) sim%naesmd%atomtype2(Na)='V'
            if(sim%naesmd%atomtype(Na).eq.24) sim%naesmd%massmdqt(Na)=51.996d0*convm
            if(sim%naesmd%atomtype(Na).eq.24) sim%naesmd%atomtype2(Na)='Cr'
            if(sim%naesmd%atomtype(Na).eq.25) sim%naesmd%massmdqt(Na)=54.938d0*convm
            if(sim%naesmd%atomtype(Na).eq.25) sim%naesmd%atomtype2(Na)='Mn'
            if(sim%naesmd%atomtype(Na).eq.26) sim%naesmd%massmdqt(Na)=55.845d0*convm
            if(sim%naesmd%atomtype(Na).eq.26) sim%naesmd%atomtype2(Na)='Fe'
            if(sim%naesmd%atomtype(Na).eq.27) sim%naesmd%massmdqt(Na)=58.933d0*convm
            if(sim%naesmd%atomtype(Na).eq.27) sim%naesmd%atomtype2(Na)='Co'
            if(sim%naesmd%atomtype(Na).eq.28) sim%naesmd%massmdqt(Na)=58.693d0*convm
            if(sim%naesmd%atomtype(Na).eq.28) sim%naesmd%atomtype2(Na)='Ni'
            if(sim%naesmd%atomtype(Na).eq.29) sim%naesmd%massmdqt(Na)=63.546d0*convm
            if(sim%naesmd%atomtype(Na).eq.29) sim%naesmd%atomtype2(Na)='Cu'
            if(sim%naesmd%atomtype(Na).eq.30) sim%naesmd%massmdqt(Na)=65.38d0*convm
            if(sim%naesmd%atomtype(Na).eq.30) sim%naesmd%atomtype2(Na)='Zn'
            if(sim%naesmd%atomtype(Na).eq.31) sim%naesmd%massmdqt(Na)=69.723d0*convm
            if(sim%naesmd%atomtype(Na).eq.31) sim%naesmd%atomtype2(Na)='Ga'
            if(sim%naesmd%atomtype(Na).eq.32) sim%naesmd%massmdqt(Na)=72.640d0*convm
            if(sim%naesmd%atomtype(Na).eq.32) sim%naesmd%atomtype2(Na)='Ge'
            if(sim%naesmd%atomtype(Na).eq.33) sim%naesmd%massmdqt(Na)=74.922d0*convm
            if(sim%naesmd%atomtype(Na).eq.33) sim%naesmd%atomtype2(Na)='As'
            if(sim%naesmd%atomtype(Na).eq.34) sim%naesmd%massmdqt(Na)=78.960d0*convm
            if(sim%naesmd%atomtype(Na).eq.34) sim%naesmd%atomtype2(Na)='Se'
            if(sim%naesmd%atomtype(Na).eq.35) sim%naesmd%massmdqt(Na)=79.904d0*convm
            if(sim%naesmd%atomtype(Na).eq.35) sim%naesmd%atomtype2(Na)='Br'
            if(sim%naesmd%atomtype(Na).eq.36) sim%naesmd%massmdqt(Na)=83.798d0*convm
            if(sim%naesmd%atomtype(Na).eq.36) sim%naesmd%atomtype2(Na)='Kr'
            if(sim%naesmd%atomtype(Na).eq.37) sim%naesmd%massmdqt(Na)=85.468d0*convm
            if(sim%naesmd%atomtype(Na).eq.37) sim%naesmd%atomtype2(Na)='Rb'
            if(sim%naesmd%atomtype(Na).eq.38) sim%naesmd%massmdqt(Na)=87.620d0*convm
            if(sim%naesmd%atomtype(Na).eq.38) sim%naesmd%atomtype2(Na)='Sr'
            if(sim%naesmd%atomtype(Na).eq.39) sim%naesmd%massmdqt(Na)=88.906d0*convm
            if(sim%naesmd%atomtype(Na).eq.39) sim%naesmd%atomtype2(Na)='Y'
            if(sim%naesmd%atomtype(Na).eq.40) sim%naesmd%massmdqt(Na)=91.224d0*convm
            if(sim%naesmd%atomtype(Na).eq.40) sim%naesmd%atomtype2(Na)='Zr'
            if(sim%naesmd%atomtype(Na).eq.41) sim%naesmd%massmdqt(Na)=92.906d0*convm
            if(sim%naesmd%atomtype(Na).eq.41) sim%naesmd%atomtype2(Na)='Nb'
            if(sim%naesmd%atomtype(Na).eq.42) sim%naesmd%massmdqt(Na)=95.950d0*convm
            if(sim%naesmd%atomtype(Na).eq.42) sim%naesmd%atomtype2(Na)='Mo'
            if(sim%naesmd%atomtype(Na).eq.43) sim%naesmd%massmdqt(Na)=98.000d0*convm
            if(sim%naesmd%atomtype(Na).eq.43) sim%naesmd%atomtype2(Na)='Tc'
            if(sim%naesmd%atomtype(Na).eq.44) sim%naesmd%massmdqt(Na)=101.07d0*convm
            if(sim%naesmd%atomtype(Na).eq.44) sim%naesmd%atomtype2(Na)='Ru'
            if(sim%naesmd%atomtype(Na).eq.45) sim%naesmd%massmdqt(Na)=102.91d0*convm
            if(sim%naesmd%atomtype(Na).eq.45) sim%naesmd%atomtype2(Na)='Rh'
            if(sim%naesmd%atomtype(Na).eq.46) sim%naesmd%massmdqt(Na)=106.42d0*convm
            if(sim%naesmd%atomtype(Na).eq.46) sim%naesmd%atomtype2(Na)='Pd'
            if(sim%naesmd%atomtype(Na).eq.47) sim%naesmd%massmdqt(Na)=107.87d0*convm
            if(sim%naesmd%atomtype(Na).eq.47) sim%naesmd%atomtype2(Na)='Ag'
            if(sim%naesmd%atomtype(Na).eq.48) sim%naesmd%massmdqt(Na)=112.41d0*convm
            if(sim%naesmd%atomtype(Na).eq.48) sim%naesmd%atomtype2(Na)='Cd'
            if(sim%naesmd%atomtype(Na).eq.49) sim%naesmd%massmdqt(Na)=114.82d0*convm
            if(sim%naesmd%atomtype(Na).eq.49) sim%naesmd%atomtype2(Na)='In'
            if(sim%naesmd%atomtype(Na).eq.50) sim%naesmd%massmdqt(Na)=118.71d0*convm
            if(sim%naesmd%atomtype(Na).eq.50) sim%naesmd%atomtype2(Na)='Sn'
            if(sim%naesmd%atomtype(Na).eq.51) sim%naesmd%massmdqt(Na)=121.76d0*convm
            if(sim%naesmd%atomtype(Na).eq.51) sim%naesmd%atomtype2(Na)='Sb'
            if(sim%naesmd%atomtype(Na).eq.52) sim%naesmd%massmdqt(Na)=127.60d0*convm
            if(sim%naesmd%atomtype(Na).eq.52) sim%naesmd%atomtype2(Na)='Te'
            if(sim%naesmd%atomtype(Na).eq.53) sim%naesmd%massmdqt(Na)=126.90d0*convm
            if(sim%naesmd%atomtype(Na).eq.53) sim%naesmd%atomtype2(Na)='I'
            if(sim%naesmd%atomtype(Na).eq.54) sim%naesmd%massmdqt(Na)=131.29d0*convm
            if(sim%naesmd%atomtype(Na).eq.54) sim%naesmd%atomtype2(Na)='Xe'
            if(sim%naesmd%atomtype(Na).eq.55) sim%naesmd%massmdqt(Na)=132.91d0*convm
            if(sim%naesmd%atomtype(Na).eq.55) sim%naesmd%atomtype2(Na)='Cs'
            if(sim%naesmd%atomtype(Na).eq.56) sim%naesmd%massmdqt(Na)=137.33d0*convm
            if(sim%naesmd%atomtype(Na).eq.56) sim%naesmd%atomtype2(Na)='Ba'
            if(sim%naesmd%atomtype(Na).eq.78) sim%naesmd%massmdqt(Na)=195.08d0*convm
            if(sim%naesmd%atomtype(Na).eq.78) sim%naesmd%atomtype2(Na)='Pt'
            if(sim%naesmd%atomtype(Na).eq.79) sim%naesmd%massmdqt(Na)=196.97d0*convm
            if(sim%naesmd%atomtype(Na).eq.79) sim%naesmd%atomtype2(Na)='Au'
            if(sim%naesmd%atomtype(Na).eq.80) sim%naesmd%massmdqt(Na)=200.59d0*convm
            if(sim%naesmd%atomtype(Na).eq.80) sim%naesmd%atomtype2(Na)='Hg'
            if(sim%naesmd%atomtype(Na).eq.82) sim%naesmd%massmdqt(Na)=207.20d0*convm
            if(sim%naesmd%atomtype(Na).eq.82) sim%naesmd%atomtype2(Na)='Pb'
            if(sim%naesmd%atomtype(Na).eq.83) sim%naesmd%massmdqt(Na)=208.98d0*convm
            if(sim%naesmd%atomtype(Na).eq.83) sim%naesmd%atomtype2(Na)='Bi'
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

        !Reading normal modes to freeze
        if (sim%naesmd%nmc.gt.0) then
            rewind (inputfdes)
            do  
                read(inputfdes,'(a)',err=29) txt 
                if( &
                    txt(1:6).eq.'$MODES' &
                    .or.txt(1:6).eq.'$modes' &
                    .or.txt(1:6).eq.'&MODES' &
                    .or.txt(1:6).eq.'&modes') exit ! exiting infinite loop
            end do
            i=1 
            do  
                read(inputfdes,'(a)',err=29) txt 
                if( &
                    txt(1:9).eq.'$ENDMODES' &
                    .or.txt(1:9).eq.'$endmodes' &
                    .or.txt(1:9).eq.'&ENDMODES' &
                    .or.txt(1:9).eq.'&endmodes') exit
                read(txt,*,err=29) sim%naesmd%mc(i)
                i=i+1
            end do
            write(6,*) 'Normal modes to freeze:'
            do i=1,sim%naesmd%nmc
                write(6,*) sim%naesmd%mc(i)
            enddo
            !Reading all normal modes for freezing
            open(34,file='nma_modes.out')
            do i=1,3*sim%naesmd%natom
               print *, 'Reading degree # ', i
               read(34,*) (sim%naesmd%enm(i,j),j=1,sim%naesmd%natom*3)
!               read(34,198) (sim%naesmd%enm(i,j),j=1,sim%naesmd%natom*3)
            enddo
            close(34)
            !Reading equilibrium positions for normal mode freezing
            open(34,file='reference.xyz')
            read(34,*)
            read(34,*)
            do i=1,sim%naesmd%natom
               read(34,*) txt, sim%naesmd%xbf0(i), sim%naesmd%ybf0(i),&
                          sim%naesmd%zbf0(i) 
               sim%naesmd%xbf0(i) = sim%naesmd%xbf0(i)/convl
               sim%naesmd%ybf0(i) = sim%naesmd%ybf0(i)/convl
               sim%naesmd%zbf0(i) = sim%naesmd%zbf0(i)/convl
            enddo
            close(34)
            !Calculating total mass
            sim%naesmd%masatotal = sum(sim%naesmd%massmdqt)
            do i=1, sim%naesmd%natom
               sim%naesmd%massqrt(i) = sqrt(sim%naesmd%massmdqt(i))
            enddo
        else if (sim%naesmd%npc.gt.0) then
            !Reading pairs of distances to freeze
            rewind (inputfdes)
            do
                read(inputfdes,'(a)',err=29) txt
                if( &
                    txt(1:6).eq.'$PAIRS' &
                    .or.txt(1:6).eq.'$pairs' &
                    .or.txt(1:6).eq.'&PAIRS' &
                    .or.txt(1:6).eq.'&pairs') exit ! exiting infinite loop
            end do
            i=1
            do
                read(inputfdes,'(a)',err=29) txt
                if( &
                    txt(1:9).eq.'$ENDPAIRS' &
                    .or.txt(1:9).eq.'$endpairs' &
                    .or.txt(1:9).eq.'&ENDPAIRS' &
                    .or.txt(1:9).eq.'&endpairs') exit
                read(txt,*,err=29) sim%naesmd%dtc(i,1),sim%naesmd%dtc(i,2)
                i=i+1
            end do
            write(6,*) 'Pairs of atoms to freeze:'
            do i=1,sim%naesmd%npc
                write(6,*) sim%naesmd%dtc(i,1),sim%naesmd%dtc(i,2)
            enddo
        endif

        if(sim%naesmd%npot.gt.0) then
            allocate(sim%naesmd%yg(2*sim%excN))
            allocate(sim%naesmd%ygprime(2*sim%excN))
            allocate(sim%naesmd%yg_new(3*sim%excN))
            allocate(sim%naesmd%ygprime_new(3*sim%excN))
            allocate(sim%naesmd%cross(sim%excN))
            allocate(sim%naesmd%Omega(sim%excN))
            sim%naesmd%Omega = 0.0d0
            allocate(sim%naesmd%sgn(sim%excN,sim%excN))
        endif
        if(sim%naesmd%dynam_type.eq.'aimc'.and.sim%id.eq.0) then
            !Output files
!            nuclear%outfile_0=400 !for debugging
            nuclear%outfile_1=401 !for populations 
            nuclear%outfile_2=402 !for nuclear coefficients
            nuclear%outfile_3=403 !for electronic overlaps
            nuclear%outfile_4=404 !for last Heff (for restart)
            nuclear%outfile_5=405 !for last Pha (for restart)
            nuclear%outfile_6=406 !for dropped trajectories (for restart)
            call open_output_mce(nuclear,sim%naesmd%tfemto)
            !Effective Hamiltonians
            allocate(nuclear%Heff(Nsim_max,Nsim_max),nuclear%Heff_old(Nsim_max,Nsim_max))
            nuclear%Heff=0.0d0
            nuclear%Heff_old=0.0d0
            allocate(nuclear%cloned(Nsim_max))
            nuclear%cloned=.false.
            !Nuclear cofficients D
            allocate(nuclear%D(Nsim))
            nuclear%D = (0.0d0, 0.0d0)
            nuclear%D(1) = (1.0d0, 0.0d0)
            !Electronic populations
            allocate(nuclear%pop(sim%excN))
            !Defoult initial values for the electronic overaps
            allocate(nuclear%sE(Nsim_max,Nsim_max,sim%excN,sim%excN))
            nuclear%sE = 0.0d0
            do i = 1, sim%excN
                do j = 1, Nsim
                    nuclear%sE(j,j,i,i) = 1.0d0
                enddo
            enddo
        endif
        if(sim%naesmd%dynam_type.eq.'aimc') then
            !Initial defoult value for the nuclear phase (gamma)
            sim%naesmd%pha = 0.0d0
            !Gaussian widths for MCE (1/Bohr^2)
            allocate(sim%naesmd%w(3*sim%naesmd%natom))
            do i = 1, sim%naesmd%natom
                if(sim%naesmd%atomtype(i).eq.1) then !H
                    sim%naesmd%w(3*i-2) = 4.7d0
                    sim%naesmd%w(3*i-1) = 4.7d0
                    sim%naesmd%w(3*i) = 4.7d0
                elseif(sim%naesmd%atomtype(i).eq.5) then !B
                    sim%naesmd%w(3*i-2) = 15.2d0
                    sim%naesmd%w(3*i-1) = 15.2d0
                    sim%naesmd%w(3*i) = 15.2d0
                elseif(sim%naesmd%atomtype(i).eq.6) then !C
                    sim%naesmd%w(3*i-2) = 22.7d0
                    sim%naesmd%w(3*i-1) = 22.7d0
                    sim%naesmd%w(3*i) = 22.7d0
                elseif(sim%naesmd%atomtype(i).eq.7) then !N
                    sim%naesmd%w(3*i-2) = 19.0d0
                    sim%naesmd%w(3*i-1) = 19.0d0
                    sim%naesmd%w(3*i) = 19.0d0
                elseif(sim%naesmd%atomtype(i).eq.8) then !O
                    sim%naesmd%w(3*i-2) = 12.2d0
                    sim%naesmd%w(3*i-1) = 12.2d0
                    sim%naesmd%w(3*i) = 12.2d0
                elseif(sim%naesmd%atomtype(i).eq.9) then !F
                    sim%naesmd%w(3*i-2) = 8.5d0
                    sim%naesmd%w(3*i-1) = 8.5d0
                    sim%naesmd%w(3*i) = 8.5d0
                elseif(sim%naesmd%atomtype(i).eq.12) then !Mg (Calculated for the Chlorophyll monomer)
                    sim%naesmd%w(3*i-2) = 10.3
                    sim%naesmd%w(3*i-1) = 10.3
                    sim%naesmd%w(3*i) = 10.3
                elseif(sim%naesmd%atomtype(i).eq.16) then !S
                    sim%naesmd%w(3*i-2) = 16.7d0
                    sim%naesmd%w(3*i-1) = 16.7d0
                    sim%naesmd%w(3*i) = 16.7d0
                elseif(sim%naesmd%atomtype(i).eq.17) then !Cl
                    sim%naesmd%w(3*i-2) = 7.4d0
                    sim%naesmd%w(3*i-1) = 7.4d0
                    sim%naesmd%w(3*i) = 7.4d0
                else
                    print *, 'Fatal error, Gaussian width mising for atomic number: ', sim%naesmd%atomtype(i)
                    stop
                endif
            enddo
        if(restart_flag.eq.1) call read_aimc_for_restart(nuclear,Nsim,sim%excN)
        endif

        do i=1,sim%excN
          do j=1, sim%excN
            sim%naesmd%sgn(i,j)=1.0d0
          enddo
        enddo

        !--------------------------------------------------------------------
        !
        !  Reading velocities
        !
        !--------------------------------------------------------------------

        sim%naesmd%vx = 0.0d0
        sim%naesmd%vy = 0.0d0
        sim%naesmd%vz = 0.0d0
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

            read(txt,*,err=29) sim%naesmd%yg(i), sim%naesmd%yg(i+sim%excN)!, sim%naesmd%yg_new(i), sim%naesmd%yg_new(i+sim%excN), sim%naesmd%yg_new(i+2*sim%excN)
            sim%naesmd%yg_new(i)=sim%naesmd%yg(i)*dcos(sim%naesmd%yg(i+sim%excN))
            sim%naesmd%yg_new(i+sim%excN)=sim%naesmd%yg(i)*dsin(sim%naesmd%yg(i+sim%excN))
            sim%naesmd%yg_new(i+2*sim%excN)=0.0d0
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
198     format(2000(e18.10,1x))


        sim%coords(1:3*Na)=sim%md%r0(1:3*Na)
        allocate(sim%deriv_forces(3*Na))
        allocate(sim%deriv_forces_state(sim%excN,3*Na))

        !Reading input.ceon lines for restart 
        open(10000-1,file='input.ceon',status='old')
        j=1
        do
          read(10000-1, '(A)', iostat=ios) txt
          if (ios /= 0 .or. txt(2:6).eq.'coord') exit
          j = j + 1
        end do
        allocate(sim%input_line(j)) 
        rewind(10000-1)
        do i = 1, j
          read(10000-1, '(A)', iostat=ios) sim%input_line(i)
        end do
        close(10000-1)

        if(sim%naesmd%dynam_type.eq.'aimc') then
                allocate(sim%aimc%FM(3*Na))
                allocate(sim%aimc%FE(3*Na))
                allocate(sim%aimc%Fmax(3*Na))
                sim%aimc%nclones=nclones0
        endif

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
