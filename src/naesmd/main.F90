#include "dprec.fh"
#include "assert.fh"

#ifndef STDOUT
#define STDOUT 6
#endif


program MD_Geometry

    use qmmm_module,only:qmmm_struct,deallocate_qmmm
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

    integer ido,neq,idocontrol,lprint
    integer Na,Nm,N1,N2,N3,natoms
    integer i,j,k
    integer ii,jjj,ibo
    integer imdqt,iimdqt
    _REAL_ h,d,t,Ek,Etot
    _REAL_,allocatable::fosc(:)
    _REAL_,target,allocatable::Omega(:)
    _REAL_,target::E0
    _REAL_,pointer::pOmega(:),pE0
    integer ither,win
    _REAL_ tini,tend,toldivprk,param(50)
    integer slen
    integer ifind
    logical moldyn_found
    type(simulation_t),target::sim_notp
    type(simulation_t),pointer::sim
    ! variables of the moldyn namelist
    integer bo_dynamics_flag,exc_state_init,n_exc_states_propagate
    integer out_count_init
    _REAL_ time_init,rk_tolerance,time_step
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
    namelist /moldyn/ natoms,bo_dynamics_flag,exc_state_init, &
        n_exc_states_propagate,out_count_init,time_init, &
        rk_tolerance,time_step,n_class_steps,n_quant_steps, &
        quant_coeffs_reinit,num_deriv_step, &
        therm_temperature,therm_type, &
        berendsen_relax_const,heating, &
        heating_steps_per_degree,out_data_steps,out_coords_steps, &
        therm_friction,rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag, &
        quant_step_reduction_factor,decoher_e0,decoher_c,decoher_type,dotrivial
    integer cohertype
    _REAL_,allocatable :: yg(:)
    integer,allocatable :: cross(:)
    integer crosstot
    _REAL_ ntotcoher
    _REAL_ constcoherE0,constcoherC
    character*(150) txt
    external diff
    integer itime1,itime11,itime2,itime3,get_time
    _REAL_ time11,time12
    character*30 datetime, machname*36
    external fcn ! function call for interpolation
    _REAL_ t_start,t_finish
    integer inputfdes

    sim=>sim_notp
   
    call init0_simulation(sim)
    call init_main()

    !Put derivative variables into module
    qmmm_struct%ideriv=moldyn_deriv_flag
    qmmm_struct%numder_step=num_deriv_step


    !***********************************************************
    !Initialization and single point calculation or zeroeth
    !step of adiabatic or nonadiabatic dynamics
    !
    !Open input file, initialize variables, run first calculation
    !***********************************************************
    open (inputfdes,file='input.ceon',status='old')
    call init_sqm(sim,inputfdes,STDOUT) ! involves call to Davidson
    close(inputfdes)

    nbasis=sim%dav%Nb  ! this is number of atomic orbitals
    sim%nbasis=nbasis  ! not to confuse with Ncis or Nrpa !!!

    if(sim%excN.ne.0) call allocate_naesmd_module2(sim%dav%Ncis,sim%nbasis,sim%excN)

    ! calling the first time to check the quirality in what follows
    !  uumdqtflag =0 means that ceo is called by the first time,
    !     so m.o. matrix is stored
    !  uumdqtflag =1 means that the quirality of m.o. must be checked
    uumdqtflag=0
    !
    do i=1,sim%excN
        iorden(i)=i
        iordenhop(i)=0
    end do
     
    !uumdqtflag=1
    if(imdtype.eq.0) then
        Etot=E0+Ek
    else
        Etot=E0+Omega(imdtype)+Ek
    endif
    vgs=E0;
    vmdqt(1:sim%excN)=Omega(1:sim%excN)+vgs
    t=0.0

    itime2=get_time()
    time12=real(itime2-itime1)/100

    !Calculate derivatives
    call cpu_time(t_start)
    if (qmmm_struct%ideriv.gt.0) then
        call deriv(sim,ihop)
    elseif (nstep.gt.0) then
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
    icontw=1
    !Open output files
    call open_output(ibo,tfemto,imdtype,lprint)
    ihop=qmmm_struct%state_of_interest !FIXME change ihop to module variable
    if(tfemto.eq.0.d0) call writeoutputini(sim,ibo,yg,lprint)

    do imdqt=1,nstep !Main loop
        !Classical propagation step - BOMD or NAESMD
        write(6,*)"Begin classical propagation step #",imdqt
        tfemto=tfemto+dtmdqt*convtf
        qmmm_struct%num_qmmm_calls=imdqt !necessary for BO dynamics

        if(state.eq.'exct'.and.ibo.ne.1) then
            !NAESMD-propagate
            qmmm_struct%state_of_interest=ihop;
            call verlet1(sim)
        else
            !BOMD-propagate
            call verlet(sim)
            ihop=qmmm_struct%state_of_interest
        end if

        if(state.eq.'exct'.and.ibo.ne.1) then
            write(6,*)'Begin nonadiabatic couplings and crossings calculations'
            call initialize(sim,yg)
            !*******************************************************
            ! The analytic NAC for t.
            ! are calculated inside of cadiaboldcalc,cadiabmiddlecalc, and cadiabnewcalc
            ! Input should be xxp,yyp,zzp 'plus' xyz at t+dt
            ! xxm,yym,zzm 'minus' - xyz at t-dt, dt should correspond ~10^-4 10^-5 A shift
            ! Vectors and frequencies from the above ceo, at xx,yy,zz geometry
            ! The routine will output cadiab_analt - an array of NA couplings for T
            ! from state ihop to all other states npot
            ! xxp,yyp, and zzp are xyz at t + dtnact
            ! xxm,yym, and zzm are xyz at t - dtnact
            !****************************************************************
            ! calculation of the energies(=vmdqtold) and nacT(=cadiabold) values at the beginning
            ! of the classical step
            !*******************************************************
            call cadiaboldcalc(sim,imdqt)
            !**************************************************************
            ! calculation of the energies(=vmdqtnew) and nacT(=cadiabnew) values at the end
            ! of the classical step
            !*******************************************************
            call cadiabnewcalc(sim)
            !**************************************************************
            ! check crossing, if crossing takes place, nquantumstep=nstepcross*nquantumreal
            !***************************************************************
            
            if(dotrivial.eq.1) then
            call checkcrossing(sim,cross,lprint)
            else
               write(6,*)'WARNING:TRIVIAL UNAVOIDED CROSSING DETECTION IS OFF'
               do i=1,npot
                  cross(i)=0
               enddo
            endif

            nquantumstep=nquantumreal
            dtquantum=dtmdqt/dfloat(nquantumstep)
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
                cadiabhop=cadiabnew(qmmm_struct%state_of_interest,iorden(qmmm_struct%state_of_interest))
                nquantumstep=nquantumreal*nstepcross
                dtquantum=dtmdqt/dfloat(nquantumstep)
                do j=1,sim%excN
                    lowvalue(j)=1000.0d0
                end do
                do iimdqt=1,nquantumstep
                    tfemtoquantum=tfemto-dtmdqt*convtf &
                        +iimdqt*dtquantum*convtf
                    call vmdqtmiddlecalc(sim,iimdqt,Na)
                    call checkcrossingmiddle(sim,cross)
                    if(lprint.ge.2) then
                        do j=1,sim%excN
                            if(cross(j).eq.1) then
                                write(101,888) tfemto,tfemtoquantum, &
                                    cross(j),j,iordenhop(j),scpr(j,iordenhop(j)), &
                                    scprreal(j,iordenhop(j)), &
                                    vmdqtmiddle(j)-vmdqtmiddle(iordenhop(j))
                                call flush(101)
                            end if
                        end do
                    end if
                end do
                do win=1,3
                    do j=1,sim%excN
                        if(cross(j).eq.1) then
                            if(j.lt.iordenhop(j).or.j.eq.ihop) then
                                if(j.ne.iordenhop(j)) then
                                    call vmdqtlowvalue(sim,win,lowvaluestep(j),Na)
                                    call checkcrossingmiddle(sim,cross)
                                    if(lprint.ge.2) then
                                        write(101,888) tfemto,tfemto-dtmdqt*convtf &
                                            +lowvaluestep(j)*dtquantum*convtf, &
                                            cross(j),j,iordenhop(j),scpr(j,iordenhop(j)), &
                                            scprreal(j,iordenhop(j)),lowvalue(j)
                                        call flush(101)
                                    end if
                                end if
                            end if
                        end if
                    end do
                end do
                nquantumstep=nquantumreal
                dtquantum=dtmdqt/dfloat(nquantumstep)
            end if

            !
            !  remove the couplings if cross=2
            !
            do i=1,sim%excN
                if(i.eq.qmmm_struct%state_of_interest) then

                    if(cross(i).eq.2) then
                        nquantumstep=nquantumreal

                        if(conthop.gt.0) then
                            if(iordenhop(i).ne.ihopprev) conthop=0
                        end if

                        if(conthop2.gt.0) then
                            if(iordenhop(i).ne.ihopprev) conthop2=0
                        end if

                        if(conthop.eq.0) then
                            cadiabold(i,iorden(i))=0.0d0
                            cadiabold(iorden(i),i)=0.0d0
                            cadiabnew(i,iorden(i))=0.0d0
                            cadiabnew(iorden(i),i)=0.0d0
                        end if
                    end if
                else
                    if(i.ne.iorden(qmmm_struct%state_of_interest)) then
                        if(i.lt.iorden(i)) then
                            if(cross(i).eq.2) then
                                nquantumstep=nquantumreal
                                cadiabold(i,iorden(i))=0.0d0
                                cadiabold(iorden(i),i)=0.0d0
                                cadiabnew(i,iorden(i))=0.0d0
                                cadiabnew(iorden(i),i)=0.0d0
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
         
            do iimdqt=1,nquantumstep
                write(6,*)'Begin quantum step #',iimdqt
                tfemtoquantum=tfemto-dtmdqt*convtf &
                    +iimdqt*dtquantum*convtf
                ! Definition of initial and final time for the quantum propagator
                if(imdqt.eq.1.and.iimdqt.eq.1) then
                    tini=0.0d0
                else
                    tini=tend
                end if

                tend=tfemtoquantum/convtf
                tini0=tini
                !--------------------------------------------------------------------
                ! Calculation of the coordinates at the middle of the step.
                ! This intermediate structure is obtained using Newton equation
                ! with constant velocities and accelerations.
                ! and calculation of the energies(=vmdqtmiddle) and nacT(=cadiabmiddle) values
                ! at intermediate times of the classical step
                !--------------------------------------------------------------------
                call cadiabmiddlecalc(sim,iimdqt,Na,cross)

                if(lprint.ge.3) then
                    if(iimdqt.ne.nquantumstep) then
                        write(96,889) tfemtoquantum, &
                            vgs*feVmdqt,(vmdqtmiddle(j)*feVmdqt,j=1,sim%excN)

                        write(93,889) tfemtoquantum, &
                            ((cadiabmiddle(j,k),k=1,sim%excN),j=1,sim%excN)
                        call flush(96)
                        call flush(93)
                    end if
                end if
                !--------------------------------------------------------------------
                !
                !  Calculation of parameters to fit the values of cadiab and
                !  vmdqt during propagation
                !
                !--------------------------------------------------------------------
                call fitcoef(sim)
                !--------------------------------------------------------------------
                ! Runge-Kutta-Verner propagator
                !--------------------------------------------------------------------
                call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
                ! Check the norm
                call checknorm(sim,ido,neq,tini,tend,toldivprk,param,yg)
                !******************************************************
                ! values for hop probability
                do k=1,sim%excN
                    do j=1,sim%excN
                        vnqcorrhop(k,j)=-1.0d0*yg(j) &
                            *dcos(yg(j+sim%excN)-yg(k+sim%excN))*cadiabmiddle(k,j)
                        vnqcorrhop(k,j)=vnqcorrhop(k,j)*2.0d0*yg(k)
                    end do
                end do

                do k=1,sim%excN
                    do j=1,sim%excN
                        vnqcorrhoptot(k,j)=vnqcorrhoptot(k,j) &
                            +vnqcorrhop(k,j)*dtquantum
                    end do
                end do
                write(6,*)'End quantum step #',iimdqt
            end do

            !write(6,*)'Now doing some other things'
            !--------------------------------------------------------------------
            ! last part of velocity verlet algorithm
            ! for ehrenfest should go after evalhop
            call verlet2(sim)
            !--------------------------------------------------------------------
            ! analyze the hopping
            !--------------------------------------------------------------------
            call evalhop(sim,lprint,ido,neq,tini,tend,toldivprk, &
                param,Na,yg,cross,idocontrol)
            !--------------------------------------------------------------------
            if(icontw.eq.nstepw) then
                if(state.eq.'exct') then
                    if(lprint.ge.1) then
                        ntotcoher=0.0d0

                        do j=1,sim%excN
                            ntotcoher=ntotcoher+yg(j)**2
                        end do

                        write(105,999) ihop,tfemto,(yg(j)**2,j=1,sim%excN), &
                            ntotcoher
                        call flush(105)
                    end if
                end if
            end if

            if(constcoherE0.ne.0.d0.and.constcoherC.ne.0.d0) then
                if(conthop.ne.1.and.conthop2.ne.1) then
                    call coherence(sim,Na,ido,neq,tini,tend,toldivprk, &
                        param,yg,constcoherE0,constcoherC,cohertype,idocontrol)
                end if
            end if
           !write(6,*)'Now finished with the other things'
        !--------------------------------------------------------------------
        end if


        ! evaluation of current kinetic energy
        call temperature(imdqt)
        if(imdqt.eq.1) then
            icont=1
            icontpdb=icontini
        end if

        if(icontw.ne.nstepw) then
            icontw=icontw+1
        else
            icontw=1
            call writeoutput(sim,ibo,yg,lprint,cross)
        end if
    end do

    !  final call to divprk to release workspace
    !  which was automatically allocated by the initial call with IDO=1.
    if(state.eq.'exct'.and.ibo.ne.1) then
        if(1==0) then ! kav: to skip it whatsoever
                       ! jakb: not sure if this should be skipped or not
            ido=3
            call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
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
    write(6,*) ' SQM (ground state) took overall [s]:'
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
    subroutine init_main()
        use naesmd_module
        use md_module

        implicit none
        _REAL_,allocatable::xx(:),yy(:),zz(:)

        ! dtnact is the incremental time to be used at nact calculation
        dtnact=0.002d0

        ! Setting divprk propagator parameters
        call sset(50,0.0d0,param,1)
        param(10)=1.d0

        ! Maximum number of steps allowed
        param(4)=5.d6

        call get_date(datetime)
        call get_machine(machname)
        itime1=get_time()

        ! definition of variables
        i=0
        ktbig(i)=char(48)//char(48)//char(48)//char(48)
        do i=1,9
            ktbig(i)=char(48)//char(48)//char(48)//char(48+i)
        end do

        ii=9
        do j=1,9
            do i=0,9
                ii=ii+1
                ktbig(ii)=char(48)//char(48)//char(48+j)//char(48+i)
            end do
        end do

        ii=99
        do k=1,9
            do j=0,9
                do i=0,9
                    ii=ii+1
                    ktbig(ii)=char(48)//char(48+k)//char(48+j)//char(48+i)
                end do
            end do
        end do

        ii=999
        do jjj=1,9
            do k=0,9
                do j=0,9
                    do i=0,9
                        ii=ii+1
                        ktbig(ii)=char(48+jjj)//char(48+k)//char(48+j)//char(48+i)
                    end do
                end do
            end do
        end do

        ! ido= Flag indicating the state of the computation
        ! use by divprk propagator
        ido=1

        conthop=0
        conthop2=0
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
        do while(readstring(12,cardini,slen).ge.0)
            j=j+1

            do i=1,slen
                txtinput(j)(i:i)=cardini(i:i)
            end do

            if(cardini(1:6).eq.'$COORD') then
                txtinicoord=j
            end if

            if(cardini(1:9).eq.'$ENDCOORD') then
                txtendcoord=j
            end if
        end do
        close(12)

        jend=j
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
        quant_coeffs_reinit=0 ! reinit of quantum coeffs after a hop (1-yes,0-no)
        num_deriv_step=1.d-3 ! finite step for numerical derivatives, A
        therm_temperature=300.d0 ! Thermostat temperature, K
        therm_type=0 ! Thermostate type (0-no thermostat, 1-Langevin,2-Berendsen)
        berendsen_relax_const=0.4d0 ! Berendsen bath relaxation constant
        heating=0 ! 0-equilibrium dynamics, 1-heating
        heating_steps_per_degree=100 ! number of steps per degree during heating
        out_data_steps=1 ! number of steps to write data
        out_coords_steps=10 ! number of steps to write coordinates (to restart)
        therm_friction=2.d0 ! friction coefficient, ps^{-1}
        rnd_seed=1 ! seed for random number generator
        out_data_cube=0 ! write view files to generate cubes, 0-no, 1-yes
        verbosity=2 ! output verbosity (0-minimal,3-largest)
        moldyn_deriv_flag=1 ! derivatives (0-numer,1-analyt,2-fast GS analytical)
        quant_step_reduction_factor=0.1d0 ! ???
        decoher_e0=0.1d0 ! decoherence parameter E0
        decoher_c=0.1d0 ! decoherecne parameter C
        decoher_type=0 ! decoherence type, Persico/Granucci(0), or Truhlar (1)
        dotrivial=1 !do trivial unavoided crossing routine (1) or not (0)  
 
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

        ! Tentatively, loading new parameter from the moldyn namelist to
        ! old parameters
        ibo=bo_dynamics_flag
   
        imdtype=exc_state_init
        ihop=imdtype
        if(imdtype==0) then
            state='fund'
        else
            state='exct'
        end if

        npot=n_exc_states_propagate
        neq=2*npot !Number of equations for divprk -- could use 2*sim%excN instead
   
        icontini=out_count_init
        tfemto=time_init
        toldivprk=rk_tolerance

        dtmdqt=time_step
        h=dtmdqt
        dtmdqt=dtmdqt/convtf

        nstep=n_class_steps

        nquantumreal=n_quant_steps
        dtquantum=dtmdqt/dfloat(nquantumreal)

        decorhop=quant_coeffs_reinit
        d=num_deriv_step

        temp0=therm_temperature
        ttt=temp0

        ither=therm_type
        if(ither==0) ensemble ='energy'
        if(ither==1) ensemble ='langev'
        if(ither==2) ensemble ='temper'

        tao=berendsen_relax_const

        if(heating==1) then
            prep='heat'
        else
            prep='equi'
        end if

        istepheat=heating_steps_per_degree
        nstepw=out_data_steps
        nstepcoord=out_coords_steps
   
        friction=therm_friction
        friction=friction/1.d3
        friction=friction*convtf

        iseedmdqt=rnd_seed
        write(6,*)"rnd_seed=",rnd_seed
        iview=out_data_cube
        lprint=verbosity
        ideriv=moldyn_deriv_flag

        nstepcross=int(1.d0/quant_step_reduction_factor)
        constcoherE0=decoher_e0
        constcoherC=decoher_c

        cohertype=decoher_type

        write(6,*) '!!!!!!-----MD INPUT-----!!!!!!'
        write(6,*)

        if (imdtype.eq.0) then
            qmmm_struct%state_of_interest=imdtype
            write(6,*)'Ground state MD,      imdtype=',imdtype
    
        else
            qmmm_struct%state_of_interest=imdtype
            write(6,*)'Excited state MD,      state=',imdtype
            write(6,*)'Number of states to propagate',npot
        end if

        write(6,*)'Initial count for out files   ',icontini

        if(ensemble.eq.'langev') then
            write(6,*)'MD using Langevin '
            write(6,*)'with friction coefficient [1/ps]=',friction*1.0D3

        else
            if(ensemble.eq.'temper') then
                write(6,*)'MD using Berendsen thermostat'
                write(6,*)'with bath relaxation constant [ps]=',tao
            else
                write(6,*)'MD at constant energy'
            end if
        end if

        write(6,*)'Starting time, fs        ',tfemto
        write(6,*)'Time step, fs		      ',dtmdqt*convtf
        write(6,*)'Number of classical steps (CS) ',nstep
        write(6,*)'Quantum steps/per CS         ',nquantumreal
        !	write(6,*)	'Vibr. modes to consider      ',   symm
        write(6,*)'Displacement for deriv. [A]   ',d

        if(ither.eq.0) then
            write(6,*)'Newtonian dynamics,           ither=',ither
        else if(ither.eq.1) then
            write(6,*)'Langevin thermostat dynamics, ither=',ither
        else if(ither.eq.2) then
            write(6,*)'Temperature thermostat,       ither=',ither
        end if

        write(6,*)'Temperature, K               ',temp0
        write(6,*)'Bath relaxation constant     ',tao

        if(prep.eq.'equi') write(6,*)'Equilibrated dynamics'
        if(prep.eq.'heat') write(6,*)'Heating dynamics from T=0'

        write(6,*)'Number of steps per degree K  ',istepheat
        write(6,*)'Number of steps to write data ',nstepw
        write(6,*)'Number of steps to write coords',nstepcoord
        write(6,*)'Friction coefficient [1/ps]   ',friction
        write(6,*)'Seed for random generator    ',iseedmdqt

        if(iview.eq.1) write(6,*)'Will write files to generate cubes'

        write(6,*)'NAESMD verbosity, ',lprint

        if(ideriv.eq.0) then
            write(6,*)'MD will use numerical deriv, ideriv=',ideriv
        else if(ideriv.ge.2) then
            write(6,*)'MD will use fast GS analytic deriv,  ideriv=',ideriv
        else
            write(6,*)'MD will use analytic deriv,  ideriv=',ideriv
        end if

        write(6,*) 'Cartesian coordinates are used'


        !***************************************************
        ! Initial allocations
        ! **************************************************
        call allocate_naesmd_module_init(natoms)
        call allocate_md_module_init(natoms)
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

            read(txt,*,err=29) atoms(Na),r0(3*Na-2),r0(3*Na-1),r0(3*Na)
            atmass(Na)=2*atoms(Na)
            atomtype(Na)=atoms(Na)
            !         if(atomtype(Na).eq.1) massmdqt(Na)=1.0079d0*convm
            if(atomtype(Na).eq.1) massmdqt(Na)=1.00d0*convm
            if(atomtype(Na).eq.1) atomtype2(Na)='H '
            if(atomtype(Na).eq.2) massmdqt(Na)=4.0026d0*convm
            if(atomtype(Na).eq.2) atomtype2(Na)='He'
            if(atomtype(Na).eq.3) massmdqt(Na)=6.941d0*convm
            if(atomtype(Na).eq.3) atomtype2(Na)='Li'
            if(atomtype(Na).eq.4) massmdqt(Na)=9.0122d0*convm
            if(atomtype(Na).eq.4) atomtype2(Na)='Be'
            if(atomtype(Na).eq.5) massmdqt(Na)=10.811d0*convm
            if(atomtype(Na).eq.5) atomtype2(Na)='B '
            if(atomtype(Na).eq.6) massmdqt(Na)=12.011d0*convm
            if(atomtype(Na).eq.6) atomtype2(Na)='C '
            if(atomtype(Na).eq.7) massmdqt(Na)=14.007d0*convm
            if(atomtype(Na).eq.7) atomtype2(Na)='N '
            if(atomtype(Na).eq.8) massmdqt(Na)=15.9994d0*convm
            if(atomtype(Na).eq.8) atomtype2(Na)='O '
            if(atomtype(Na).eq.9) massmdqt(Na)=18.998d0*convm
            if(atomtype(Na).eq.9) atomtype2(Na)='F '
            if(atomtype(Na).eq.10) massmdqt(Na)=20.180d0*convm
            if(atomtype(Na).eq.10) atomtype2(Na)='Ne'
            if(atomtype(Na).eq.14) massmdqt(Na)=28.086d0*convm
            if(atomtype(Na).eq.14) atomtype2(Na)='Si'
            if(atomtype(Na).eq.15) massmdqt(Na)=30.974d0*convm
            if(atomtype(Na).eq.15) atomtype2(Na)='P '
            if(atomtype(Na).eq.16) massmdqt(Na)=32.065d0*convm
            if(atomtype(Na).eq.16) atomtype2(Na)='S '
            if(atomtype(Na).eq.17) massmdqt(Na)=35.453d0*convm
            if(atomtype(Na).eq.17) atomtype2(Na)='Cl'

            read(12,'(a)',err=29) txt
        end do
        !--------------------------------------------------------------------
        !Allocate variables
        write(6,*)'Allocating variables'
        natom=Na
        sim%Na=Na
        allocate(sim%coords(Na*3))
        sim%excN=npot
        call allocate_naesmd_module(Na,npot)
        call allocate_md_module(Na)
        if(npot.gt.0) then
                allocate(yg(2*sim%excN))
                allocate(cross(sim%excN))
                allocate(fosc(sim%excN))
                allocate(Omega(sim%excN))
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

            read(txt,*,err=29) vx(i),vy(i),vz(i)
            i=i+1
        end do

        vx(1:natom)=vx(1:natom)/convl*convt
        vy(1:natom)=vy(1:natom)/convl*convt
        vz(1:natom)=vz(1:natom)/convl*convt
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
        v(1:Nm,1:Nm)=0.d0
        do i=1,Na
            if (atoms(i).eq.1) then
                fo(3*i-2)=3000
                fo(3*i-1)=3000
                fo(3*i)=3000
            else
                fo(3*i-2)=1000
                fo(3*i-1)=1000
                fo(3*i)=1000
            end if
   
            fm(3*i-2)=atmass(i)
            fm(3*i-1)=atmass(i)
            fm(3*i)=atmass(i)
            v(3*i-2,3*i-2)=1.0
            v(3*i-1,3*i-1)=1.0
            v(3*i,3*i)=1.0
        end do

        ! initialize variables:
        N3=3*Na
        N2=2*Nm
        N1=Nm+1

        ! compute initial cartesian coordinates:
        do j=3,N3,3
            xx(j/3)=r0(j-2)
            yy(j/3)=r0(j-1)
            zz(j/3)=r0(j)
        end do

        ! define cartesian coordinates for mdqt
        ! transform coordiantes from angstrom to atomic units
        rx(1:natom)=xx(1:natom)/convl
        ry(1:natom)=yy(1:natom)/convl
        rz(1:natom)=zz(1:natom)/convl

        ! position of the center of mass
        xcmini=0.d0
        ycmini=0.d0
        zcmini=0.d0

        masstot=0.d0 ! kav: - zeroing was not here, needed?
        do j=1,natom
            masstot=masstot+massmdqt(j)
        end do
        do j=1,natom
            xcmini=xcmini+rx(j)*massmdqt(j)/masstot
            ycmini=ycmini+ry(j)*massmdqt(j)/masstot
            zcmini=zcmini+rz(j)*massmdqt(j)/masstot
        end do
   
        !Remove rotation and translation from initial velocity
        write(6,*)'Rescaling velocity'
        call rescaleveloc(rx,ry,rz,vx,vy,vz,massmdqt,natom)

        ! compute kinetic energy, for cartesian option
        kin=0.d0
        do i=1,natom
            kin=kin+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)/2
        end do

        write(6,*)
        write(6,*)'Initial Geometry'

        do i=1,Na
            write(6,100) atoms(i),xx(i),yy(i),zz(i)
        end do
100     format(I5,'     ',3F12.6)


        if(tfemto.eq.0.d0) then
            open (9,file='coords.xyz')
            write (9,*) Na
            write (9,*) 'Input Geometry'

            do j = 1,Na
                write(9,302) ELEMNT(atoms(j)),xx(j),yy(j),zz(j)
            end do
            close (9)
        end if
        sim%coords(1:3*Na)=r0(1:3*Na)
        allocate(sim%deriv_forces(3*Na))
        pOmega=>Omega
        pE0=>E0

        call init_naesmd_space_globs(sim%naesmd, pOmega, pE0)

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
