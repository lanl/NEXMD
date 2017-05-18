#include "dprec.fh"
#include "assert.fh"

module additional_subroutines
    use naesmd_constants

contains
    !----- SUBROUTINES

    integer function readstring(file,card,flen)
        integer       file, flen
        character*200 card
        if(file.gt.200)STOP'ERROR: file number too large'
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
    subroutine open_output(ibo,tfemto,imdtype,lprint)
        implicit none

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
            write(98,'(a)') '##     time(fs)           T                  &
	      &T-T0               U                  U-U0               E                  E-E0'
 
            open(92,file='temperature.out',status=file_status,access=file_access)

            ! creating header
            write(92,'(a)') '## time(fs)    Temp(current)       Temp(thermostat)'

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
            if (imdtype.gt.0.and.ibo.ne.1)  &
                OPEN(105, FILE= 'coeff-n-before.out', &
                status=file_status, access=file_access)
        endif

        if(lprint.ge.2) then
            if (imdtype.gt.0..and.ibo.ne.1) OPEN(33,FILE='hops-trial.out', &
                status=file_status, access=file_access)
            if (imdtype.gt.0..and.ibo.ne.1) OPEN(100,FILE='order.out', &
                status=file_status, access=file_access)
            OPEN(85, FILE= 'fx.out', status=file_status, &
                access=file_access)
            OPEN(84, FILE= 'fy.out', status=file_status, &
                access=file_access)
            OPEN(83, FILE= 'fz.out', status=file_status, &
                access=file_access)
            OPEN(101, FILE= 'cross-steps.out', status=file_status, &
                access=file_access)
        endif

        if(lprint.ge.3) then
            if (imdtype.gt.0.and.ibo.ne.1) OPEN(94, FILE= 'coeff-q.out',   &
                status=file_status, access=file_access)
            OPEN(91, FILE= 'cm.out', status=file_status, &
                access=file_access)
        endif

        return
    end SUBROUTINE

    SUBROUTINE checknorm(sim,ido,neq,tini,tend,toldivprk,param,yg,idocontrol)
        use naesmd_module
        use md_module
        use communism
        IMPLICIT NONE
        type(simulation_t), pointer :: sim
        integer k,ido,neq, idocontrol
        double precision tini,tend,toldivprk,norm,normdiff,param(50)
        double precision yg(sim%excN)
        external fcn

        norm=0.0d0
        do k = 1,sim%excN
            norm=norm+yg(k)*yg(k)
        enddo
        normdiff=dabs(norm-1.0d0)
        if(normdiff.ge.1.0d-5) then
            do k = 1,sim%excN
                yg(k)=yg(k)/dsqrt(norm)
            enddo
            if(idocontrol.eq.0) then
                    ido=3
                    write(45,*) tini*convtf,normdiff
                    call divprk(ido,neq,fcn,tini,tend, &
                        toldivprk,param,yg)
                    ido=1
                    idocontrol=1
            endif
        endif

        RETURN
    END SUBROUTINE

    ! Subroutine to calculate the coefficients to fit to a linear eq. the values of
    ! cadiab and vmdqt during propagation

    SUBROUTINE fitcoef(sim)
        use naesmd_module
        use md_module
        use communism
        IMPLICIT NONE
        type(simulation_t), pointer :: sim
        integer k,j

        do k=1,sim%excN
            do j=1,sim%excN
                bcoeffcadiab(k,j)=(cadiabmiddle(k,j)- &
                    cadiabmiddleold(k,j))/dtquantum
            enddo
        enddo
        
        do k=1,sim%excN
            bcoeffvmdqt(k)=(vmdqtmiddle(k)- &
                vmdqtmiddleold(k))/dtquantum
        enddo

        RETURN
    END SUBROUTINE

    ! Subroutine that initialize the values of
    ! vnqcorrhoptot that are integrated in time during each quantum propagation
    ! it also store the value of yg(ihop) at t to be used in the hopping evaluation

    SUBROUTINE initialize(sim,yg)
        use naesmd_module
        use md_module
        use communism
        IMPLICIT NONE

        type(simulation_t), pointer :: sim
        integer k,j
        double precision yg(sim%excN)!yg(nmaxpot) 

        nqold=yg(ihop)

        do k=1,sim%excN
            do j=1,sim%excN
                vnqcorrhoptot(k,j)=0.0d0
            enddo
        enddo

        tfemtoquantum=0.0d0

        RETURN
    END SUBROUTINE

end module
