#include "dprec.fh"
#include "assert.fh"

module additional_subroutines
    use naesmd_constants
    use communism
contains
    !----- SUBROUTINES

    integer function readstring(file,card,flen)
        integer       file, flen
        character*200 card
        !if(file.gt.200)STOP 'ERROR: file number too large'
        read(file,'(a)',err=100,end=100)card
        flen=199
        do while(card(flen:flen).eq.' ')
            flen=flen-1
            if (flen.eq.0) goto 100
        enddo
        readstring=flen
        return
100     readstring=-1
        return
    end function

    !
    !********************************************************************
    !
    !********************************************************************
    !
    subroutine open_output_multi(sim,ibo,tfemto,imdtype,lprint)
        implicit none
        type(simulation_t),pointer::sim

        _REAL_ tfemto
        integer imdtype,ibo,lprint
        character*10 file_status
        character*10 file_access
        character(100) ::  filename

        if(tfemto.eq.0.d0) then
            file_access='sequential'
            file_status='unknown'
        else
            file_access='append'
            file_status='old'
        end if

        ! xyz.out: time(fs), (atomic_number(i), xyz-coordinates(i)(Amstrong), i=1,number_of_atoms)
        ! veloc.out: time(fs), (atomic_number(i), xyz-veloc(i)(Amstrong/ps), i=1,number_of_atoms)
        ! energy-ev.out: time(fs),kinetic energy(t),kinetic energy(t)-kinetic energy(0),
        !                potential energy(t),potential energy(t)-potential energy(0),
        !                total energy(t),total energy(t)-total energy(0). Units=eV
        ! energy-au.out: time(fs),kinetic energy, potential energy, total energy in atomic units
        ! pes.out: time(fs), energy(ground state, excited state(n), n=1,# of states), units=eV
        ! coeff-n.out: current state where the nuclei are propagated(it changes at each hop), time(fs),
        !              square of the module of each quantum coefficients(yg(i)**2,i=1,# of excited states),
        !              sum of them (it should be equal 1)
        ! coeff-n-before.out: current state where the nuclei are propagated(it changes at each hop), time(fs),
        !              square of the module of each quantum coefficients(yg(i)**2,i=1,# of excited states)
        !              before the decoherence correction
        !              sum of them (it should be equal 1)
        ! coeff-q.out: time(fs), fase of each quantum coefficient(yg(i+npot)),
        !              expressed as sin(yg(i+npot),i=1,# of excited states)
        ! nact.out: time(fs), nonadiabatic couplings cadiab(j,k)=(d_jk(nonadiabatic coupling vector between
        !           states j and k)*(nuclear velocities vector). It is written as ((cadiab(j,k),k=1,# of
        !           excited states),j=1,# of excited states)
        ! temperature.out: time(fs), temperature, ideal temperature
        ! cm.out: time(fs), ix,y,z coordinates of the center of mass at t-coordinates of the center of mass at 0
        ! fx.out, fy.out,fz.out: forces in x,y,z for each atom, in ev/amstrong
        ! nacr.out: atom_#,xyz-coordinates of the nonadiabatic coupling vector on this atom,
        !           it is written each time the probability to hop is greater that the random number, despite
        !           that the hop can be forbidden by energy requirements
        ! hops-trial.out: time(fs), #_of_state_at_which_to_system_will_try_to_hop,
        !                 parameter_that_indicates_that_the_hop_has_been_done
        !                 (0 for successful hops and 1 for forbidden hops)
        !                 it is written each time the probability to hop is greater that the random number, despite
        !                 that the hop can be forbidden by energy requirements
        ! hops.out: time(fs), #_of_state_before_hop, #_of_state_after_hop,total_energy_before,
        !           total_energy_after (should be equal), units=a.u. It is written each time the hop takes place.
        ! transition-densities.out: time(fs),(transition_density(i), i=1, number_of_atomic_orbitals)
        ! transition-densities-tot.out: time(fs),(transition_density(i,k), i=1, number_of_atomic_orbitals,k=1,4)
        ! order.out: time(fs),diabatic order of states respect to the initial order at t=0
        ! cross-steps.out: time(fs), overlap between the hypothetical crossing states. It is written when the overlap between
        ! the states is < 0.9 and therefore a reduction of the quantum step is required.

        sim%outfile_23=34*10000+sim%id
        write (filename, "(a3,i4.4,a4)") "34_", sim%id, ".out"
        OPEN(sim%outfile_23,FILE=trim(filename) ,action='write',STATUS=file_status,ACCESS=file_access)
        sim%outfile_24=110*10000+sim%id
        write (filename, "(a4,i4.4,a4)") "110_", sim%id, ".out"
        OPEN(sim%outfile_24,FILE=trim(filename) ,action='write',STATUS=file_status,ACCESS=file_access) 
        sim%outfile_20=202*10000+sim%id
        write (filename, "(a9,i4.4,a4)") "velocity_", sim%id, ".out"
        OPEN(sim%outfile_20,FILE=trim(filename) ,action='write',STATUS=file_status,ACCESS=file_access)
        sim%outfile_21=203*10000+sim%id
        write (filename, "(a12,i4.4,a4)") "coefficient_", sim%id, ".out"
        if(ibo==0) OPEN(sim%outfile_21,FILE=trim(filename) ,action='write',STATUS=file_status,ACCESS=file_access)
        sim%outfile_22=9*10000+sim%id
        write (filename, "(a7,i4.4,a4)") "coords_", sim%id, ".xyz"
        OPEN(sim%outfile_22,file=trim(filename) ,action='write',STATUS=file_status,ACCESS=file_access)
        if(lprint.ge.0) then
            !
            ! use for vibration***************************
            write (filename, "(A4,I4.4,A4)") "xyz_", sim%id, ".out"
            sim%outfile_18=99*10000+sim%id 
            !open(sim%outfile_18,file=trim(filename),status=file_status,access=file_access)
            write (filename, "(A6,I4.4,A4)") "veloc_", sim%id, ".out"
            sim%outfile_19=80*10000+sim%id 
            !open(sim%outfile_19,file=trim(filename),status=file_status,access=file_access)
            !*****************************************
            !

            write (filename, "(A10,I4.4,A4)") "energy-ev_", sim%id, ".out"
            sim%outfile_1=98*10000+sim%id 
            open(sim%outfile_1,file=trim(filename),status=file_status,access=file_access)

            ! creating header
            ! BTN: only print header when not restarting
            if(tfemto.eq.0.d0) then
              write(sim%outfile_1,'(a)') '##     time(fs)           T                  &
	      &T-T0               U                  U-U0               E                  E-E0'
	          end if
            write (filename, "(A12,I4.4,A4)") "temperature_", sim%id, ".out"
            sim%outfile_2=92*10000+sim%id  
            open(sim%outfile_2,file=trim(filename),status=file_status,access=file_access)

            ! creating header
            ! BTN: only print header when not restarting
            if(tfemto.eq.0.d0) then
              write(sim%outfile_2,'(a)') '## time(fs)    Temp(current)       Temp(thermostat)'
            end if 

            if(imdtype.gt.0.and.ibo.ne.1) then
                write (filename, "(A5,I4.4,A4)") "hops_", sim%id, ".out"
                sim%outfile_3=30*10000+sim%id 
                open(sim%outfile_3,file=trim(filename),status=file_status,access=file_access)
      
                ! headaer
                write(sim%outfile_3,*) '## time [fs] type_of_hop state_before state_after energy_before [a.u.] energy_after [a.u]'
            end if
        end if

        if(lprint.ge.1) then
            write(6,*)'LPRINT:',lprint,imdtype
            if (imdtype.gt.0) then 
                write (filename, "(A4,I4.4,A4)") "pes_", sim%id, ".out"
                sim%outfile_4=96*10000+sim%id 
                open(sim%outfile_4,file=trim(filename),status=file_status,access=file_access)
            end if
            if (imdtype.gt.0.)  then
                write (filename, "(A21,I4.4,A4)") "transition-densities_", sim%id, ".out"
                sim%outfile_5=89*10000+sim%id 
                open(sim%outfile_5,file=trim(filename),status=file_status,access=file_access)
            endif
            if (imdtype.gt.0.and.ibo.ne.1) then
                write (filename, "(A5,I4.4,A4)") "nacr_", sim%id, ".out"
                sim%outfile_6=29*10000+sim%id 
                open(sim%outfile_6,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A5,I4.4,A4)") "nact_", sim%id, ".out"
                sim%outfile_7=93*10000+sim%id 
                open(sim%outfile_7,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A8,I4.4,A4)") "coeff-n_", sim%id, ".out"
                sim%outfile_8=95*10000+sim%id 
                open(sim%outfile_8,file=trim(filename),status=file_status,access=file_access)
!BTN: removed file coeff-n-before.out. grep this line to undo
                write (filename, "(A15,I4.4,A4)") "coeff-n-before_", sim%id, ".out"
                sim%outfile_9=105*10000+sim%id 
!                open(sim%outfile_9,file=trim(filename),status=file_status,access=file_access)
            endif
        endif

        if(lprint.ge.2) then
            if (imdtype.gt.0..and.ibo.ne.1) then
                write (filename, "(A11,I4.4,A4)") "hops-trial_", sim%id, ".out"
                sim%outfile_10=33*10000+sim%id 
                open(sim%outfile_10,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A6,I4.4,A4)") "order_", sim%id, ".out"
                sim%outfile_11=100*10000+sim%id 
                open(sim%outfile_11,file=trim(filename),status=file_status,access=file_access)
            endif
                write (filename, "(A3,I4.4,A4)") "fy_", sim%id, ".out"
                sim%outfile_13=84*10000+sim%id 
            !    open(sim%outfile_13,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A3,I4.4,A4)") "fz_", sim%id, ".out"
                sim%outfile_14=83*10000+sim%id 
             !   open(sim%outfile_14,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A12,I4.4,A4)") "cross-steps_", sim%id, ".out"
                sim%outfile_15=101*10000+sim%id 
                open(sim%outfile_15,file=trim(filename),status=file_status,access=file_access)
        endif

        if(lprint.ge.3) then
            if (imdtype.gt.0.and.ibo.ne.1) then
                write (filename, "(A8,I4.4,A4)") "coeff-q_", sim%id, ".out"
                sim%outfile_16=94*10000+sim%id 
                open(sim%outfile_16,file=trim(filename),status=file_status,access=file_access)
            endif
                write (filename, "(A3,I4.4,A4)") "cm_", sim%id, ".out"
                sim%outfile_17=91*10000+sim%id 
                open(sim%outfile_17,file=trim(filename),status=file_status,access=file_access)
                write (filename, "(A7,I4.4,A4)") "forces_", sim%id, ".out"
                sim%outfile_12=85*10000+sim%id 
                open(sim%outfile_12,file=trim(filename),status=file_status,access=file_access)
        endif

        return
    end SUBROUTINE

    subroutine open_output(sim,ibo,tfemto,imdtype,lprint)
        implicit none
        type(simulation_t),pointer::sim

        _REAL_ tfemto
        integer imdtype,ibo,lprint
        character*10 file_status
        character*10 file_access

        if(tfemto.eq.0.d0) then
            file_access='sequential'
            file_status='unknown'
        else
            file_access='append'
            file_status='old'
        end if
        sim%outfile_1=98
        sim%outfile_2=92
        sim%outfile_3=30
        sim%outfile_4=96
        sim%outfile_5=89
        sim%outfile_6=29
        sim%outfile_7=93
        sim%outfile_8=95
        sim%outfile_9=105
        sim%outfile_10=33
        sim%outfile_11=100
        sim%outfile_12=85
        sim%outfile_13=84
        sim%outfile_14=83
        sim%outfile_15=101
        sim%outfile_16=94
        sim%outfile_17=91
        sim%outfile_18=99
        sim%outfile_19=80
        sim%outfile_20=202
        sim%outfile_21=203
        sim%outfile_22=9
        sim%outfile_23=34
        sim%outfile_24=110
        OPEN(sim%outfile_23,FILE='34.out' ,action='write',STATUS=file_status,ACCESS=file_access)
        OPEN(sim%outfile_24,FILE='110.out' ,action='write',STATUS=file_status,ACCESS=file_access)
       
        OPEN(sim%outfile_20,FILE='velocity.out' ,action='write',STATUS=file_status,ACCESS=file_access)
        if(ibo==0) OPEN(sim%outfile_21,FILE='coefficient.out' ,action='write',STATUS=file_status,ACCESS=file_access)
        OPEN(sim%outfile_22,file='coords.xyz' ,action='write',STATUS=file_status,ACCESS=file_access)

        ! xyz.out: time(fs), (atomic_number(i), xyz-coordinates(i)(Amstrong), i=1,number_of_atoms)
        ! veloc.out: time(fs), (atomic_number(i), xyz-veloc(i)(Amstrong/ps), i=1,number_of_atoms)
        ! energy-ev.out: time(fs),kinetic energy(t),kinetic energy(t)-kinetic energy(0),
        !                potential energy(t),potential energy(t)-potential energy(0),
        !                total energy(t),total energy(t)-total energy(0). Units=eV
        ! energy-au.out: time(fs),kinetic energy, potential energy, total energy in atomic units
        ! pes.out: time(fs), energy(ground state, excited state(n), n=1,# of states), units=eV
        ! coeff-n.out: current state where the nuclei are propagated(it changes at each hop), time(fs),
        !              square of the module of each quantum coefficients(yg(i)**2,i=1,# of excited states),
        !              sum of them (it should be equal 1)
        ! coeff-n-before.out: current state where the nuclei are propagated(it changes at each hop), time(fs),
        !              square of the module of each quantum coefficients(yg(i)**2,i=1,# of excited states)
        !              before the decoherence correction
        !              sum of them (it should be equal 1)
        ! coeff-q.out: time(fs), fase of each quantum coefficient(yg(i+npot)),
        !              expressed as sin(yg(i+npot),i=1,# of excited states)
        ! nact.out: time(fs), nonadiabatic couplings cadiab(j,k)=(d_jk(nonadiabatic coupling vector between
        !           states j and k)*(nuclear velocities vector). It is written as ((cadiab(j,k),k=1,# of
        !           excited states),j=1,# of excited states)
        ! temperature.out: time(fs), temperature, ideal temperature
        ! cm.out: time(fs), ix,y,z coordinates of the center of mass at t-coordinates of the center of mass at 0
        ! fx.out, fy.out,fz.out: forces in x,y,z for each atom, in ev/amstrong
        ! nacr.out: atom_#,xyz-coordinates of the nonadiabatic coupling vector on this atom,
        !           it is written each time the probability to hop is greater that the random number, despite
        !           that the hop can be forbidden by energy requirements
        ! hops-trial.out: time(fs), #_of_state_at_which_to_system_will_try_to_hop,
        !                 parameter_that_indicates_that_the_hop_has_been_done
        !                 (0 for successful hops and 1 for forbidden hops)
        !                 it is written each time the probability to hop is greater that the random number, despite
        !                 that the hop can be forbidden by energy requirements
        ! hops.out: time(fs), #_of_state_before_hop, #_of_state_after_hop,total_energy_before,
        !           total_energy_after (should be equal), units=a.u. It is written each time the hop takes place.
        ! transition-densities.out: time(fs),(transition_density(i), i=1, number_of_atomic_orbitals)
        ! transition-densities-tot.out: time(fs),(transition_density(i,k), i=1, number_of_atomic_orbitals,k=1,4)
        ! order.out: time(fs),diabatic order of states respect to the initial order at t=0
        ! cross-steps.out: time(fs), overlap between the hypothetical crossing states. It is written when the overlap between
        ! the states is < 0.9 and therefore a reduction of the quantum step is required.

        if(lprint.ge.0) then
            !
            ! use for vibration***************************
            !          OPEN(99, FILE= 'xyz.out', status=file_status,
            !     $access=file_access)
            !        OPEN(80, FILE= 'veloc.out', status=file_status,
            !     $access=file_access)
            !*****************************************
            !

            open(98,file='energy-ev.out',status=file_status,access=file_access)

            ! creating header
            ! BTN: only print header when not restarting
            if(tfemto.eq.0.d0) then
              write(98,'(a)') '##     time(fs)           T                  &
	      &T-T0               U                  U-U0               E                  E-E0'
	          end if
 
            open(92,file='temperature.out',status=file_status,access=file_access)

            ! creating header
            ! BTN: only print header when not restarting
            if(tfemto.eq.0.d0) then
              write(92,'(a)') '## time(fs)    Temp(current)       Temp(thermostat)'
            end if 

            if(imdtype.gt.0.and.ibo.ne.1) then
                open(30,file='hops.out',status=file_status, access=file_access)
      
                ! headaer
                write(30,*) '## time [fs] type_of_hop state_before state_after energy_before [a.u.] energy_after [a.u]'
            end if
        end if

        if(lprint.ge.1) then
            write(6,*)'LPRINT:',lprint,imdtype
            if (imdtype.gt.0) OPEN(96, FILE= 'pes.out',  &
                status=file_status, access=file_access)
            if (imdtype.gt.0.)  &
                OPEN(89, FILE= 'transition-densities.out', &
                status=file_status, access=file_access)
            !          if (imdtype.gt.0.)
            !     $ OPEN(77, FILE= 'transition-densities-tot.out',
            !     $ status=file_status, access=file_access)
            if (imdtype.gt.0.and.ibo.ne.1) OPEN(29, FILE= 'nacr.out',   &
                status=file_status, access=file_access)
            if (imdtype.gt.0.and.ibo.ne.1) OPEN(93, FILE= 'nact.out',   &
                status=file_status, access=file_access)
            if (imdtype.gt.0.and.ibo.ne.1) OPEN(95, FILE= 'coeff-n.out', &
                status=file_status, access=file_access)
!BTN: removed file coeff-n-before.out. grep this line to undo
!            if (imdtype.gt.0.and.ibo.ne.1)  &
!                OPEN(105, FILE= 'coeff-n-before.out', &
!                status=file_status, access=file_access)
        endif

        if(lprint.ge.2) then
            if (imdtype.gt.0..and.ibo.ne.1) OPEN(33,FILE='hops-trial.out', &
                status=file_status, access=file_access)
            if (imdtype.gt.0..and.ibo.ne.1) OPEN(100,FILE='order.out', &
                status=file_status, access=file_access)
!            OPEN(84, FILE= 'fy.out', status=file_status, &
!                access=file_access)
!            OPEN(83, FILE= 'fz.out', status=file_status, &
!                access=file_access)
            OPEN(101, FILE= 'cross-steps.out', status=file_status, &
                access=file_access)
        endif

        if(lprint.ge.3) then
            if (imdtype.gt.0.and.ibo.ne.1) OPEN(94, FILE= 'coeff-q.out',   &
                status=file_status, access=file_access)
            OPEN(91, FILE= 'cm.out', status=file_status, &
                access=file_access)
            OPEN(85, FILE= 'forces.out', status=file_status, &
                access=file_access)
        endif

        return
    end SUBROUTINE
    SUBROUTINE checknorm(rk_comm,nco,t,tmax,rk_tolerance,thresholds,yg)
       use rksuite_90, only: rk_comm_real_1d, setup 
	 IMPLICIT NONE
        integer nco,  k
        _REAL_ t, tmax,rk_tolerance,norm,normdiff
        _REAL_, dimension(:) :: yg, thresholds
        type(rk_comm_real_1d) :: rk_comm 

        norm=0.0d0
        do k = 1,nco
            norm=norm+yg(k)*yg(k)
        enddo
        normdiff=dabs(norm-1.0d0)
        if(normdiff.ge.1.0d-5) then
            do k = 1,nco
                yg(k)=yg(k)/dsqrt(norm)
            enddo
	    write(45,*) t*convtf,normdiff
	    call    setup(rk_comm, t, yg, tmax, rk_tolerance, thresholds, 'M','R')
        endif

        RETURN
    END SUBROUTINE

    ! Subroutine to calculate the coefficients to fit to a linear eq. the values of
    ! sim%naesmd%cadiab and sim%naesmd%vmdqt during propagation

    SUBROUTINE fitcoef(sim)
        use communism
        IMPLICIT NONE
        type(simulation_t), pointer :: sim
        integer k,j

        do k=1,sim%excN
            do j=1,sim%excN
                sim%naesmd%bcoeffcadiab(k,j)=(sim%naesmd%cadiabmiddle(k,j)- &
                    sim%naesmd%cadiabmiddleold(k,j))/sim%naesmd%dtquantum
            enddo
        enddo
        
        do k=1,sim%excN
            sim%naesmd%bcoeffvmdqt(k)=(sim%naesmd%vmdqtmiddle(k)- &
                sim%naesmd%vmdqtmiddleold(k))/sim%naesmd%dtquantum
        enddo

        RETURN
    END SUBROUTINE

    ! Subroutine that initialize the values of
    ! sim%naesmd%vnqcorrhoptot that are integrated in time during each quantum propagation
    ! it also store the value of yg(sim%naesmd%ihop) at t to be used in the hopping evaluation

    SUBROUTINE initialize(sim,yg)
        use communism
        IMPLICIT NONE

        type(simulation_t), pointer :: sim
        integer k,j
        double precision yg(sim%excN)!yg(nmaxpot) 

        sim%naesmd%nqold=yg(sim%naesmd%ihop)

        do k=1,sim%excN
            do j=1,sim%excN
                sim%naesmd%vnqcorrhoptot(k,j)=0.0d0
            enddo
        enddo

        sim%naesmd%tfemtoquantum=0.0d0

        RETURN
    END SUBROUTINE

end module
