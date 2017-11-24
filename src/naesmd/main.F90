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

    !integer ido,neq,idocontrol,lprint
    integer neq,lprint
    integer Na,Nm,N1,N2,N3,natoms
    integer i,j,k
    integer ii,jjj,ibo
    integer imdqt,iimdqt
    _REAL_ h,d,t,Ek,Etot
    _REAL_,allocatable::fosc(:)
    integer ither,win
    _REAL_ tini,tend,tmax,tgot!,toldivprk,param(50)
    integer slen
    integer ifind
    logical moldyn_found
    type(qm2_davidson_structure_type),target :: qm2ds_notp
    type(qmmm_struct_type), target :: qmmm_struct_notp
    type(qm2_structure), target :: qm2_struct_notp
    type(naesmd_structure), target :: naesmd_struct_notp
    type(md_structure), target :: md_struct_notp
    type(rk_comm_real_1d), target :: rk_struct_notp
    type(simulation_t),target::sim_notp
    type(simulation_t),pointer::sim
    ! variables of the moldyn namelist
    integer bo_dynamics_flag,exc_state_init,n_exc_states_propagate
    integer out_count_init
    _REAL_ time_init,rk_tolerance,time_step
    _REAL_,allocatable::thresholds(:)
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
    integer cohertype
    _REAL_, dimension(:), allocatable :: yg, ygprime
    integer rk_flag
    integer,allocatable :: cross(:)
    integer crosstot
    _REAL_ ntotcoher
    _REAL_ constcoherE0,constcoherC
    character*(150) txt
    external diff
    integer itime1,itime11,itime2,itime3,get_time
    _REAL_ time11,time12
    character*30 datetime, machname*36
     external :: fcn ! function call for interpolation
    _REAL_ t_start,t_finish
    integer inputfdes

    sim=>sim_notp
   
    call init0_simulation(sim)
    call setp_simulation(sim,qmmm_struct_notp,qm2ds_notp,&
		qm2_struct_notp,naesmd_struct_notp, md_struct_notp, rk_struct_notp)
    call init_main(sim, sim%naesmd, sim%md)

    !Put derivative variables into module
    sim%qmmm%ideriv=moldyn_deriv_flag
    sim%qmmm%numder_step=num_deriv_step

    open (37,file='debug.out')

    !***********************************************************
    !Initialization and single point calculation or zeroeth
    !step of adiabatic or nonadiabatic dynamics
    !
    !Open input file, initialize variables, run first calculation
    !***********************************************************
    open (inputfdes,file='input.ceon',status='old')
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
     
    !sim%naesmd%uumdqtflag=1
    if(sim%md%imdtype.eq.0) then
        Etot=sim%naesmd%E0+Ek
    else
        Etot=sim%naesmd%E0+sim%naesmd%Omega(sim%md%imdtype)+Ek
    endif
    sim%naesmd%vgs=sim%naesmd%E0;
    sim%naesmd%vmdqt(1:sim%excN)=sim%naesmd%Omega(1:sim%excN)+sim%naesmd%vgs
    t=0.0

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
    time11=0.0

    !##########################################################
    !! Molecular dynamics with Quantum transition main loop
    !##########################################################
    i=1
    sim%naesmd%icontw=1
    !Open output files
    call open_output(ibo,sim%naesmd%tfemto,sim%md%imdtype,lprint)
    sim%naesmd%ihop=sim%qmmm%state_of_interest !FIXME change sim%naesmd%ihop to module variable
    call writeoutputini(sim,ibo,yg,lprint)

    tmax=sim%naesmd%nstep*sim%naesmd%dtmdqt
    do i =1,sim%excN
	thresholds(i)=1.0d0	
	thresholds(i+sim%excN)=6.29d0	
    enddo
    call setup(sim%rk_comm,sim%naesmd%tfemto,yg,tmax,rk_tolerance,thresholds, &
	'M','R')
    do imdqt=1,sim%naesmd%nstep !Main loop
        !Classical propagation step - BOMD or NAESMD
        write(6,*)"Begin classical propagation step #",imdqt
        sim%naesmd%tfemto=sim%naesmd%tfemto+sim%naesmd%dtmdqt*convtf
        sim%qmmm%num_qmmm_calls=imdqt !necessary for BO dynamics

        if(sim%naesmd%state.eq.'exct'.and.ibo.ne.1) then
            !NAESMD-propagate
            sim%qmmm%state_of_interest=sim%naesmd%ihop;
            call verlet1(sim)
        else
            !BOMD-propagate
            call verlet(sim)
            sim%naesmd%ihop=sim%qmmm%state_of_interest
        end if

        if(sim%naesmd%state.eq.'exct'.and.ibo.ne.1) then
            write(6,*)'Begin nonadiabatic couplings and crossings calculations'
            call initialize(sim,yg)
	    
            !*******************************************************
            ! The analytic NAC for t.
            ! are calculated inside of cadiaboldcalc,cadiabmiddlecalc, and cadiabnewcalc
            ! Input should be xxp,yyp,zzp 'plus' xyz at t+dt
            ! xxm,yym,zzm 'minus' - xyz at t-dt, dt should correspond ~10^-4 10^-5 A shift
            ! Vectors and frequencies from the above ceo, at xx,yy,zz geometry
            ! The routine will output sim%naesmd%cadiab_analt - an array of NA couplings for T
            ! from sim%naesmd%state sim%naesmd%ihop to all other states sim%naesmd%npot
            ! xxp,yyp, and zzp are xyz at t + sim%naesmd%dtnact
            ! xxm,yym, and zzm are xyz at t - sim%naesmd%dtnact
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
            
            if(dotrivial.eq.1) then
            call checkcrossing(sim,cross,lprint)
            else
               write(6,*)'WARNING:TRIVIAL UNAVOIDED CROSSING DETECTION IS OFF'
               do i=1,sim%naesmd%npot
                  cross(i)=0
               enddo
            endif

            sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
            sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)
            !*******************************************
            ! loop to detect the crossing point
            !******************************************
            crosstot=0
            do i=1,sim%excN
                if(cross(i).eq.1) crosstot=1
            end do
            if(crosstot.eq.1) then
                if (lprint.gt.1) then
                    write(6,*)'there is crossing'
                    write(6,*)'cross(:)=',cross(1:sim%excN)
                endif
                sim%naesmd%cadiabhop=sim%naesmd%cadiabnew(sim%qmmm%state_of_interest,sim%naesmd%iorden(sim%qmmm%state_of_interest))
                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal*sim%naesmd%nstepcross
                sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)
                do j=1,sim%excN
                    sim%naesmd%lowvalue(j)=1000.0d0
                end do
                do iimdqt=1,sim%naesmd%nquantumstep
                    sim%naesmd%tfemtoquantum=sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                        +iimdqt*sim%naesmd%dtquantum*convtf
                    call vmdqtmiddlecalc(sim,iimdqt,Na)
                    call checkcrossingmiddle(sim,cross)
                    if(lprint.ge.2) then
                        do j=1,sim%excN
                            if(cross(j).eq.1) then
                                write(101,888) sim%naesmd%tfemto,sim%naesmd%tfemtoquantum, &
                                    cross(j),j,sim%naesmd%iordenhop(j),sim%naesmd%scpr(j,sim%naesmd%iordenhop(j)), &
                                    sim%naesmd%scprreal(j,sim%naesmd%iordenhop(j)), &
                                    sim%naesmd%vmdqtmiddle(j)-sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j))
                                call flush(101)
                            end if
                        end do
                    end if
                end do
                do win=1,3
                    do j=1,sim%excN
                        if(cross(j).eq.1) then
                            if(j.lt.sim%naesmd%iordenhop(j).or.j.eq.sim%naesmd%ihop) then
                                if(j.ne.sim%naesmd%iordenhop(j)) then
                                    call vmdqtlowvalue(sim,win,sim%naesmd%lowvaluestep(j),Na)
                                    call checkcrossingmiddle(sim,cross)
                                    if(lprint.ge.2) then
                                        write(101,888) sim%naesmd%tfemto,sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                                            +sim%naesmd%lowvaluestep(j)*sim%naesmd%dtquantum*convtf, &
                                            cross(j),j,sim%naesmd%iordenhop(j),sim%naesmd%scpr(j,sim%naesmd%iordenhop(j)), &
                                            sim%naesmd%scprreal(j,sim%naesmd%iordenhop(j)),sim%naesmd%lowvalue(j)
                                        call flush(101)
                                    end if
                                end if
                            end if
                        end if
                    end do
                end do
                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)
            end if

            !
            !  remove the couplings if cross=2
            !
            do i=1,sim%excN
                if(i.eq.sim%qmmm%state_of_interest) then

                    if(cross(i).eq.2) then
                        sim%naesmd%nquantumstep=sim%naesmd%nquantumreal

                        if(sim%naesmd%conthop.gt.0) then
                            if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop=0
                        end if

                        if(sim%naesmd%conthop2.gt.0) then
                            if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop2=0
                        end if

                        if(sim%naesmd%conthop.eq.0) then
                            sim%naesmd%cadiabold(i,sim%naesmd%iorden(i))=0.0d0
                            sim%naesmd%cadiabold(sim%naesmd%iorden(i),i)=0.0d0
                            sim%naesmd%cadiabnew(i,sim%naesmd%iorden(i))=0.0d0
                            sim%naesmd%cadiabnew(sim%naesmd%iorden(i),i)=0.0d0
                        end if
                    end if
                else
                    if(i.ne.sim%naesmd%iorden(sim%qmmm%state_of_interest)) then
                        if(i.lt.sim%naesmd%iorden(i)) then
                            if(cross(i).eq.2) then
                                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                                sim%naesmd%cadiabold(i,sim%naesmd%iorden(i))=0.0d0
                                sim%naesmd%cadiabold(sim%naesmd%iorden(i),i)=0.0d0
                                sim%naesmd%cadiabnew(i,sim%naesmd%iorden(i))=0.0d0
                                sim%naesmd%cadiabnew(sim%naesmd%iorden(i),i)=0.0d0
                            end if
                        end if
                    end if
                end if
            end do
            write(6,*)'End nonadiabatic couplings calculation'
            !--------------------------------------------------------------------
            ! Loop for quantum propagation steps
            ! that implies CEO energy calculations
            !--------------------------------------------------------------------
            write(6,*)"End classical propagation step #",imdqt
         
            do iimdqt=1,sim%naesmd%nquantumstep
                write(6,*)'Begin quantum step #',iimdqt
                sim%naesmd%tfemtoquantum=sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                    +iimdqt*sim%naesmd%dtquantum*convtf
                ! Definition of initial and final time for the quantum propagator
                if(imdqt.eq.1.and.iimdqt.eq.1) then
                    tini=0.0d0
                else
                    tini=tend
                end if

                tend=sim%naesmd%tfemtoquantum/convtf
                sim%naesmd%tini0=tini
                !--------------------------------------------------------------------
                ! Calculation of the coordinates at the middle of the step.
                ! This intermediate structure is obtained using Newton equation
                ! with constant velocities and accelerations.
                ! and calculation of the energies(=sim%naesmd%vmdqtmiddle) and nacT(=sim%naesmd%cadiabmiddle) values
                ! at intermediate times of the classical step
                !--------------------------------------------------------------------
                call cadiabmiddlecalc(sim,iimdqt,Na,cross)

                if(lprint.ge.3) then
                    if(iimdqt.ne.sim%naesmd%nquantumstep) then
                        write(96,889) sim%naesmd%tfemtoquantum, &
                            sim%naesmd%vgs*feVmdqt,(sim%naesmd%vmdqtmiddle(j)*feVmdqt,j=1,sim%excN)

                        write(93,889) sim%naesmd%tfemtoquantum, &
                            ((sim%naesmd%cadiabmiddle(j,k),k=1,sim%excN),j=1,sim%excN)
                        call flush(96)
                        call flush(93)
                    end if
                end if
                !--------------------------------------------------------------------
                !
                !  Calculation of parameters to fit the values of sim%naesmd%cadiab and
                !  sim%naesmd%vmdqt during propagation
                !
                !--------------------------------------------------------------------
                call fitcoef(sim)
                !--------------------------------------------------------------------
                ! Runge-Kutta-Verner propagator
                !--------------------------------------------------------------------
                !write(6,*)'Propagator about to be called:',ido,neq,tini,tend,toldivprk
                write(6,*)'Propagator about to be called:', tini,tend,rk_tolerance
                !write(6,*)'param:',param
                write(6,*)'yg:',yg
                write(6,*)'T_start / End :', imdqt, iimdqt, sim%rk_comm%t_start, tend, tmax
	 	tend=min(tend,tmax)	
		call range_integrate(sim%rk_comm,fcn,tend,tgot,yg,ygprime,sim%naesmd,rk_flag)
                write(6,*)'T_got :',tgot, sim%rk_comm%t_start
		tend=tgot
		!call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
                ! Check the norm
                call checknorm(sim%rk_comm,sim%excN,tend,tmax,rk_tolerance,thresholds,yg)
                !call checknorm(sim,ido,neq,tini,tend,toldivprk,param,yg,idocontrol)
                !******************************************************
                ! values for hop probability
                do k=1,sim%excN
                    do j=1,sim%excN
                        sim%naesmd%vnqcorrhop(k,j)=-1.0d0*yg(j) &
                            *dcos(yg(j+sim%excN)-yg(k+sim%excN))*sim%naesmd%cadiabmiddle(k,j)
                        sim%naesmd%vnqcorrhop(k,j)=sim%naesmd%vnqcorrhop(k,j)*2.0d0*yg(k)
                    end do
                end do

                do k=1,sim%excN
                    do j=1,sim%excN
                        sim%naesmd%vnqcorrhoptot(k,j)=sim%naesmd%vnqcorrhoptot(k,j) &
                            +sim%naesmd%vnqcorrhop(k,j)*sim%naesmd%dtquantum
                    end do
                end do
                write(6,*)'End quantum step #',iimdqt
            end do

            write(6,*)'Now doing some other things'
            !--------------------------------------------------------------------
            ! last part of velocity verlet algorithm
            ! for ehrenfest should go after evalhop
            call verlet2(sim)
            !--------------------------------------------------------------------
            ! analyze the hopping
            !--------------------------------------------------------------------
	    write(6,*) "Here?"
            call evalhop(sim, sim%rk_comm, lprint, tend, tmax, rk_tolerance, thresholds, &
                Na, yg, cross)
	    write(6,*) "no"
            write(6,*)'Now we done some other things'
            !call evalhop(sim,lprint,ido,neq,tini,tend,toldivprk, &
            !    param,Na,yg,cross,idocontrol)
            !--------------------------------------------------------------------
            if(sim%naesmd%icontw.eq.sim%naesmd%nstepw) then
                if(sim%naesmd%state.eq.'exct') then
                    if(lprint.ge.1) then
                        ntotcoher=0.0d0

                        do j=1,sim%excN
                            ntotcoher=ntotcoher+yg(j)**2
                        end do
!BTN: removed file coeff-n-before.out. grep this line to undo
!                        write(105,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN), &
!                            ntotcoher
                        call flush(105)
                    end if
                end if
            end if

            if(constcoherE0.ne.0.d0.and.constcoherC.ne.0.d0) then
                if(sim%naesmd%conthop.ne.1.and.sim%naesmd%conthop2.ne.1) then
                    !BTN: I have not had time to investigate this fully, but the simulation object gets garbled when passed to coherence. 
                    !My guess is that one of the variables passed (though it does not appear to be sim) is not defined identically correctly inside coherence.
                    !For now, just exit.
                    write(6,*) 'Your choices for decoher_e0 and decoher_c &
                        are not currently available'
                    stop
                    call coherence(sim,sim%rk_comm, Na, tend,tmax, rk_tolerance, &
                        thresholds,yg,constcoherE0,constcoherC,cohertype)
                end if
            end if
           !write(6,*)'Now finished with the other things'
        !--------------------------------------------------------------------
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
            call writeoutput(sim,ibo,yg,lprint,cross)
        end if
    end do

    !  final call to divprk to release workspace
    !  which was automatically allocated by the initial call with IDO=1.
    if(sim%naesmd%state.eq.'exct'.and.ibo.ne.1) then
        if(1==0) then ! kav: to skip it whatsoever
                       ! jakb: not sure if this should be skipped or not
            !if(idocontrol.eq.0) then
            !        ido=3
            !        call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
            !        ido=1
            !        idocontrol=1
            !endif
        end if
    end if

    !***********************************
    !      ttime=dtime(tarray)
    !      write(19,*) tarray(1)

779 format(F18.10,10000(1X,I3,3(1X,F18.10)))
886 format(10000(1X,F7.4))
888 format(F7.4,1X,F7.4,1X,I3,1X,I3,1X,I3,4(1X,F18.10))
889 format(10000(1X,F18.10))

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

29  write(6,*) txt
98  stop 'bad input'




contains
    !
    !********************************************************************
    !
    subroutine init_main(sim, naesmd_struct,md_struct)
        use naesmd_module
        use md_module
        use naesmd_constants
        use md_module

        implicit none
        type(simulation_t),pointer::sim
	type(naesmd_structure), intent(inout) :: naesmd_struct	
	type(md_structure), intent(inout) :: md_struct	
        _REAL_,allocatable::xx(:),yy(:),zz(:)

        !dtnact is the incremental time to be used at nact calculation
        naesmd_struct%dtnact=0.002d0

        ! Setting divprk propagator parameters
        !call sset(50,0.0d0,param,1)
        !param(10)=1.d0

        ! Maximum number of steps allowed
        !param(4)=5.d6

        call get_date(datetime)
        call get_machine(machname)
        itime1=get_time()

        ! definition of variables
        i=0
        naesmd_struct%ktbig(i)=char(48)//char(48)//char(48)//char(48)
        do i=1,9
            naesmd_struct%ktbig(i)=char(48)//char(48)//char(48)//char(48+i)
        end do

        ii=9
        do j=1,9
            do i=0,9
                ii=ii+1
                naesmd_struct%ktbig(ii)=char(48)//char(48)//char(48+j)//char(48+i)
            end do
        end do

        ii=99
        do k=1,9
            do j=0,9
                do i=0,9
                    ii=ii+1
                    naesmd_struct%ktbig(ii)=char(48)//char(48+k)//char(48+j)//char(48+i)
                end do
            end do
        end do

        ii=999
        do jjj=1,9
            do k=0,9
                do j=0,9
                    do i=0,9
                        ii=ii+1
                        naesmd_struct%ktbig(ii)=char(48+jjj)//char(48+k)//char(48+j)//char(48+i)
                    end do
                end do
            end do
        end do

        ! ido= Flag indicating the state of the computation
        ! use by divprk propagator
        !ido=1
        !idocontrol=1
        naesmd_struct%conthop=0
        naesmd_struct%conthop2=0
        !
        !--------------------------------------------------------------------
        !
        write (6,*)
        write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
        write (6,8) '| MD execution started at  ',datetime,' |'
        write (6,7) '| Computer: ',   machname,'|'
        write (6,*) '|________________________________________________|'
        write (6,*)

        ! copy input file
        open (12,file='input.ceon',status='old')
        j=0
        do while(readstring(12,naesmd_struct%cardini,slen).ge.0)
            j=j+1

            do i=1,slen
                naesmd_struct%txtinput(j)(i:i)=naesmd_struct%cardini(i:i)
            end do

            if(naesmd_struct%cardini(1:6).eq.'$COORD') then
                naesmd_struct%txtinicoord=j
            end if

            if(naesmd_struct%cardini(1:9).eq.'$ENDCOORD') then
                naesmd_struct%txtendcoord=j
            end if
        end do
        close(12)

        naesmd_struct%jend=j
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
 
        inputfdes=12
        open (inputfdes,file='input.ceon',status='old')

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
        
        naesmd_struct%iredpot=iredpot
        naesmd_struct%nstates=nstates 
        naesmd_struct%deltared=deltared 
        
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
        ibo=bo_dynamics_flag
   
        md_struct%imdtype=exc_state_init
        naesmd_struct%ihop=md_struct%imdtype
        if(md_struct%imdtype==0) then
            naesmd_struct%state='fund'
        else
            naesmd_struct%state='exct'
        end if

        naesmd_struct%npot=n_exc_states_propagate
        neq=2*naesmd_struct%npot !Number of equations for divprk -- could use excN instead
   
        naesmd_struct%icontini=out_count_init
        naesmd_struct%tfemto=time_init
        !toldivprk=rk_tolerance

        naesmd_struct%dtmdqt=time_step
        h=naesmd_struct%dtmdqt
        naesmd_struct%dtmdqt=naesmd_struct%dtmdqt/convtf

        naesmd_struct%nstep=n_class_steps

        naesmd_struct%nquantumreal=n_quant_steps
        naesmd_struct%dtquantum=naesmd_struct%dtmdqt/dfloat(naesmd_struct%nquantumreal)

        naesmd_struct%decorhop=quant_coeffs_reinit
        d=num_deriv_step

        naesmd_struct%temp0=therm_temperature
        md_struct%ttt=naesmd_struct%temp0

        ither=therm_type
        if(ither==0) naesmd_struct%ensemble ='energy'
        if(ither==1) naesmd_struct%ensemble ='langev'
        if(ither==2) naesmd_struct%ensemble ='temper'

        naesmd_struct%tao=berendsen_relax_const

        if(heating==1) then
            naesmd_struct%prep='heat'
        else
            naesmd_struct%prep='equi'
        end if

        naesmd_struct%istepheat=heating_steps_per_degree
        naesmd_struct%nstepw=out_data_steps
        naesmd_struct%nstepcoord=out_coords_steps
   
        naesmd_struct%friction=therm_friction/1.d3*convtf !convter therm_friction (in 1/ps) to naesmd_struct%friction (1/AU)

        naesmd_struct%iseedmdqt=rnd_seed
        write(6,*)"rnd_seed=",rnd_seed
        naesmd_struct%iview=out_data_cube
        lprint=verbosity
        md_struct%ideriv=moldyn_deriv_flag

        naesmd_struct%nstepcross=int(1.d0/quant_step_reduction_factor)
        constcoherE0=decoher_e0
        constcoherC=decoher_c

        cohertype=decoher_type

        naesmd_struct%deltared=naesmd_struct%deltared/feVmdqt !for reducing number of potentials   

        write(6,*) '!!!!!!-----MD INPUT-----!!!!!!'
        write(6,*)

        if (md_struct%imdtype.eq.0) then
            sim%qmmm%state_of_interest=md_struct%imdtype
            write(6,*)'Ground state MD,      imdtype=',md_struct%imdtype
    
        else
            sim%qmmm%state_of_interest=md_struct%imdtype
            write(6,*)'Excited state MD,      state=',md_struct%imdtype
            write(6,*)'Number of states to propagate',naesmd_struct%npot
        end if

        write(6,*)'Initial count for out files   ',naesmd_struct%icontini

        if(naesmd_struct%ensemble.eq.'langev') then
            write(6,*)'MD using Langevin '
            !write(6,*)'with friction coefficient [1/ps]=',friction*1.0D3

        else
            if(naesmd_struct%ensemble.eq.'temper') then
                write(6,*)'MD using Berendsen thermostat'
                write(6,*)'with bath relaxation constant [ps]=',naesmd_struct%tao
            else
                write(6,*)'MD at constant energy'
            end if
        end if

        write(6,*)'Starting time, fs        ',naesmd_struct%tfemto
        write(6,*)'Time step, fs		      ',naesmd_struct%dtmdqt*convtf
        write(6,*)'Number of classical steps (CS) ',naesmd_struct%nstep
        write(6,*)'Quantum steps/per CS         ',naesmd_struct%nquantumreal
        !	write(6,*)	'Vibr. modes to consider      ',   symm
        write(6,*)'Displacement for deriv. [A]   ',d

        if(ither.eq.0) then
            write(6,*)'Newtonian dynamics,           ither=',ither
        else if(ither.eq.1) then
            write(6,*)'Langevin thermostat dynamics, ither=',ither
        else if(ither.eq.2) then
            write(6,*)'Temperature thermostat,       ither=',ither
        end if

        write(6,*)'Temperature, K               ',naesmd_struct%temp0
        write(6,*)'Bath relaxation constant     ',naesmd_struct%tao

        if(naesmd_struct%prep.eq.'equi') write(6,*)'Equilibrated dynamics'
        if(naesmd_struct%prep.eq.'heat') write(6,*)'Heating dynamics from T=0'

        write(6,*)'Number of steps per degree K  ',naesmd_struct%istepheat
        write(6,*)'Number of steps to write data ',naesmd_struct%nstepw
        write(6,*)'Number of steps to write coords',naesmd_struct%nstepcoord
        write(6,*)'Friction coefficient [1/ps]   ',naesmd_struct%friction/convtf*1d3 !convert back to ps
        write(6,*)'Seed for random generator    ',naesmd_struct%iseedmdqt

        if(naesmd_struct%iview.eq.1) write(6,*)'Will write files to generate cubes'

        write(6,*)'NAESMD verbosity, ',lprint

        if(md_struct%ideriv.eq.0) then
            write(6,*)'MD will use numerical deriv, ideriv=',md_struct%ideriv
        else if(md_struct%ideriv.ge.2) then
            write(6,*)'MD will use fast GS analytic deriv,  md_struct%ideriv=',md_struct%ideriv
        else
            write(6,*)'MD will use analytic deriv,  ideriv=',md_struct%ideriv
        end if

        write(6,*) 'Cartesian coordinates are used'


        !***************************************************
        ! Initial allocations
        ! **************************************************
        call allocate_naesmd_module_init(naesmd_struct,natoms)
        call allocate_md_module_init(md_struct,natoms)
        allocate(xx(natoms),yy(natoms),zz(natoms))

        ! reading atoms and modes:
        rewind (12)
        Na=0
        ! read the cartesian coordinates
        !****************************************************************************
        do
            read(12,'(a)',err=29) txt
            if ( &
                txt(1:6).eq.'$COORD' &
                .or.txt(1:6).eq.'$coord' &
                .or.txt(1:6).eq.'&COORD' &
                .or.txt(1:6).eq.'&coord') exit ! breaking infinite loop
        end do

        read(12,'(a)',err=29) txt

        do while( &
            txt(1:9).ne.'$ENDCOORD' &
            .and.txt(1:9).ne.'$endcoord' &
            .and.txt(1:9).ne.'&ENDCOORD' &
            .and.txt(1:9).ne.'&endcoord')

            Na=Na+1 ! counting atoms

            read(txt,*,err=29) md_struct%atoms(Na),md_struct%r0(3*Na-2),md_struct%r0(3*Na-1),md_struct%r0(3*Na)
            md_struct%atmass(Na)=2*md_struct%atoms(Na)
            naesmd_struct%atomtype(Na)=md_struct%atoms(Na)
            !         if(naesmd_struct%atomtype(Na).eq.1) naesmd_struct%massmdqt(Na)=1.0079d0*convm
            if(naesmd_struct%atomtype(Na).eq.1) naesmd_struct%massmdqt(Na)=1.00d0*convm
            if(naesmd_struct%atomtype(Na).eq.1) naesmd_struct%atomtype2(Na)='H '
            if(naesmd_struct%atomtype(Na).eq.2) naesmd_struct%massmdqt(Na)=4.0026d0*convm
            if(naesmd_struct%atomtype(Na).eq.2) naesmd_struct%atomtype2(Na)='He'
            if(naesmd_struct%atomtype(Na).eq.3) naesmd_struct%massmdqt(Na)=6.941d0*convm
            if(naesmd_struct%atomtype(Na).eq.3) naesmd_struct%atomtype2(Na)='Li'
            if(naesmd_struct%atomtype(Na).eq.4) naesmd_struct%massmdqt(Na)=9.0122d0*convm
            if(naesmd_struct%atomtype(Na).eq.4) naesmd_struct%atomtype2(Na)='Be'
            if(naesmd_struct%atomtype(Na).eq.5) naesmd_struct%massmdqt(Na)=10.811d0*convm
            if(naesmd_struct%atomtype(Na).eq.5) naesmd_struct%atomtype2(Na)='B '
            if(naesmd_struct%atomtype(Na).eq.6) naesmd_struct%massmdqt(Na)=12.011d0*convm
            if(naesmd_struct%atomtype(Na).eq.6) naesmd_struct%atomtype2(Na)='C '
            if(naesmd_struct%atomtype(Na).eq.7) naesmd_struct%massmdqt(Na)=14.007d0*convm
            if(naesmd_struct%atomtype(Na).eq.7) naesmd_struct%atomtype2(Na)='N '
            if(naesmd_struct%atomtype(Na).eq.8) naesmd_struct%massmdqt(Na)=15.9994d0*convm
            if(naesmd_struct%atomtype(Na).eq.8) naesmd_struct%atomtype2(Na)='O '
            if(naesmd_struct%atomtype(Na).eq.9) naesmd_struct%massmdqt(Na)=18.998d0*convm
            if(naesmd_struct%atomtype(Na).eq.9) naesmd_struct%atomtype2(Na)='F '
            if(naesmd_struct%atomtype(Na).eq.10) naesmd_struct%massmdqt(Na)=20.180d0*convm
            if(naesmd_struct%atomtype(Na).eq.10) naesmd_struct%atomtype2(Na)='Ne'
            if(naesmd_struct%atomtype(Na).eq.14) naesmd_struct%massmdqt(Na)=28.086d0*convm
            if(naesmd_struct%atomtype(Na).eq.14) naesmd_struct%atomtype2(Na)='Si'
            if(naesmd_struct%atomtype(Na).eq.15) naesmd_struct%massmdqt(Na)=30.974d0*convm
            if(naesmd_struct%atomtype(Na).eq.15) naesmd_struct%atomtype2(Na)='P '
            if(naesmd_struct%atomtype(Na).eq.16) naesmd_struct%massmdqt(Na)=32.065d0*convm
            if(naesmd_struct%atomtype(Na).eq.16) naesmd_struct%atomtype2(Na)='S '
            if(naesmd_struct%atomtype(Na).eq.17) naesmd_struct%massmdqt(Na)=35.453d0*convm
            if(naesmd_struct%atomtype(Na).eq.17) naesmd_struct%atomtype2(Na)='Cl'

            read(12,'(a)',err=29) txt
        end do
        !--------------------------------------------------------------------
        !Allocate variables
        write(6,*)'Allocating variables'
        naesmd_struct%natom=Na
        sim%Na=Na
        allocate(sim%coords(Na*3))
        sim%excN=naesmd_struct%npot
        call allocate_naesmd_module(naesmd_struct,Na,naesmd_struct%npot)
        call allocate_md_module(md_struct,Na)
        if(naesmd_struct%npot.gt.0) then
                allocate(thresholds(2*sim%excN))
                allocate(yg(2*sim%excN))
                allocate(ygprime(2*sim%excN))
                allocate(cross(sim%excN))
                allocate(fosc(sim%excN))
                allocate(naesmd_struct%Omega(sim%excN))
        endif

        !--------------------------------------------------------------------
        !
        !  Reading velocities
        !
        !--------------------------------------------------------------------

        rewind (12)
        do
            read(12,'(a)',err=29) txt
            if( &
                txt(1:6).eq.'$VELOC' &
                .or.txt(1:6).eq.'$veloc' &
                .or.txt(1:6).eq.'&VELOC' &
                .or.txt(1:6).eq.'&veloc') exit ! exiting infinite loop
        end do

        i=1
        do
            read(12,'(a)',err=29) txt
            if( &
                txt(1:9).eq.'$ENDVELOC' &
                .or.txt(1:9).eq.'$endveloc' &
                .or.txt(1:9).eq.'&ENDVELOC' &
                .or.txt(1:9).eq.'&endveloc') exit

            read(txt,*,err=29) naesmd_struct%vx(i),naesmd_struct%vy(i),naesmd_struct%vz(i)
            i=i+1
        end do

        naesmd_struct%vx(1:naesmd_struct%natom)=naesmd_struct%vx(1:naesmd_struct%natom)/convl*convt
        naesmd_struct%vy(1:naesmd_struct%natom)=naesmd_struct%vy(1:naesmd_struct%natom)/convl*convt
        naesmd_struct%vz(1:naesmd_struct%natom)=naesmd_struct%vz(1:naesmd_struct%natom)/convl*convt
        !
        !--------------------------------------------------------------------
        ! Reading quantum coeffients
        !
        !--------------------------------------------------------------------
        !
        rewind (12)
        do
            read(12,'(a)',err=29) txt
            if( &
                txt(1:6).eq.'$COEFF' &
                .or.txt(1:6).eq.'$coeff' &
                .or.txt(1:6).eq.'&COEFF' &
                .or.txt(1:6).eq.'&coeff') exit ! leaving the infinite loop
        end do

        i=1

        if(sim%excN.gt.0) then
          do
            read(12,'(a)',err=29) txt
            if( &
                txt(1:9).eq.'$ENDCOEFF' &
                .or.txt(1:9).eq.'$endcoeff' &
                .or.txt(1:9).eq.'&ENDCOEFF' &
                .or.txt(1:9).eq.'&endcoeff') exit

            read(txt,*,err=29) yg(i),yg(i+sim%excN)

            yg(i)=dsqrt(yg(i))
            yg(i+sim%excN)=dasin(yg(i+sim%excN))

            i=i+1
            if(i.gt.sim%excN) exit ! too many coefficients - leaving the loop
          end do
        endif

        close (12)
        !
        !--------------------------------------------------------------------
        !
        write(6,*)' Finish reading coordinates '
        !
        !--------------------------------------------------------------------
        !
        Nm=3*Na
        md_struct%v(1:Nm,1:Nm)=0.d0
        do i=1,Na
            if (md_struct%atoms(i).eq.1) then
                md_struct%fo(3*i-2)=3000
                md_struct%fo(3*i-1)=3000
                md_struct%fo(3*i)=3000
            else
                md_struct%fo(3*i-2)=1000
                md_struct%fo(3*i-1)=1000
                md_struct%fo(3*i)=1000
            end if
   
            md_struct%fm(3*i-2)=md_struct%atmass(i)
            md_struct%fm(3*i-1)=md_struct%atmass(i)
            md_struct%fm(3*i)=md_struct%atmass(i)
            md_struct%v(3*i-2,3*i-2)=1.0
            md_struct%v(3*i-1,3*i-1)=1.0
            md_struct%v(3*i,3*i)=1.0
        end do

        ! initialize variables:
        N3=3*Na
        N2=2*Nm
        N1=Nm+1

        ! compute initial cartesian coordinates:
        do j=3,N3,3
            xx(j/3)=md_struct%r0(j-2)
            yy(j/3)=md_struct%r0(j-1)
            zz(j/3)=md_struct%r0(j)
        end do

        ! define cartesian coordinates for mdqt
        ! transform coordiantes from angstrom to atomic units
        naesmd_struct%rx(1:naesmd_struct%natom)=xx(1:naesmd_struct%natom)/convl
        naesmd_struct%ry(1:naesmd_struct%natom)=yy(1:naesmd_struct%natom)/convl
        naesmd_struct%rz(1:naesmd_struct%natom)=zz(1:naesmd_struct%natom)/convl

        ! position of the center of mass
        naesmd_struct%xcmini=0.d0
        naesmd_struct%ycmini=0.d0
        naesmd_struct%zcmini=0.d0

        naesmd_struct%masstot=0.d0 ! kav: - zeroing was not here, needed?
        do j=1,naesmd_struct%natom
            naesmd_struct%masstot=naesmd_struct%masstot+naesmd_struct%massmdqt(j)
        end do
        do j=1,naesmd_struct%natom
            naesmd_struct%xcmini=naesmd_struct%xcmini+naesmd_struct%rx(j)*naesmd_struct%massmdqt(j)/naesmd_struct%masstot
            naesmd_struct%ycmini=naesmd_struct%ycmini+naesmd_struct%ry(j)*naesmd_struct%massmdqt(j)/naesmd_struct%masstot
            naesmd_struct%zcmini=naesmd_struct%zcmini+naesmd_struct%rz(j)*naesmd_struct%massmdqt(j)/naesmd_struct%masstot
        end do
   
        !Remove rotation and translation from initial velocity
        write(6,*)'Rescaling velocity'
        call rescaleveloc(naesmd_struct%rx,naesmd_struct%ry,naesmd_struct%rz, &
		naesmd_struct%vx,naesmd_struct%vy,naesmd_struct%vz,naesmd_struct%massmdqt,naesmd_struct%natom)

        ! compute kinetic energy, for cartesian option
        naesmd_struct%kin=0.d0
        do i=1,naesmd_struct%natom
            naesmd_struct%kin=naesmd_struct%kin+naesmd_struct%massmdqt(i)* &
		(naesmd_struct%vx(i)**2+naesmd_struct%vy(i)**2+naesmd_struct%vz(i)**2)/2
        end do

        write(6,*)
        write(6,*)'Initial Geometry'

        do i=1,Na
            write(6,100) md_struct%atoms(i),xx(i),yy(i),zz(i)
        end do
100     format(I5,'     ',3F12.6)


!        if(naesmd_struct%tfemto.eq.0.d0) then
!            open (9,file='coords.xyz')
!            write (9,*) Na
!            write (9,*) 'Input Geometry'
!
!            do j = 1,Na
!                write(9,302) ELEMNT(md_struct%atoms(j)),xx(j),yy(j),zz(j)
!            end do
!            close (9)
!        end if
        sim%coords(1:3*Na)=md_struct%r0(1:3*Na)
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
