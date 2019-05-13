#include "dprec.fh"
#include "assert.fh"

!***************************************************
!** Verlet velocity algorithm                     **
!***************************************************
module verlet_module
    use naesmd_constants
    use langevin_temperature
    use communism
    implicit none

contains
    !
    !--------------------------------------------------------------------
    !
    subroutine verlet(sim)
        implicit none

        type(simulation_t),pointer::sim

        integer i
        _REAL_ dt_2,dt2_2
        _REAL_ t_start,t_finish

        dt_2=sim%naesmd%dtmdqt/2.d0 ! dt/2
        dt2_2=sim%naesmd%dtmdqt**2/2 ! dt^2/2
        !
        !--------------------------------------------------------------------
        !
        !  Get frictional and random terms for position and velocity
        !
        !--------------------------------------------------------------------
        !
        call sdterm(sim%naesmd)
        !
        !--------------------------------------------------------------------
        !
        !  Store the current atom positions, then find new atom
        !  positions and half-step velocities via Verlet recursion
        !
        !--------------------------------------------------------------------
        !
        sim%naesmd%rxold=sim%naesmd%rx
        sim%naesmd%ryold=sim%naesmd%ry
        sim%naesmd%rzold=sim%naesmd%rz
        sim%naesmd%vxold=sim%naesmd%vx
        sim%naesmd%vyold=sim%naesmd%vy
        sim%naesmd%vzold=sim%naesmd%vz
        sim%naesmd%axold=sim%naesmd%ax
        sim%naesmd%ayold=sim%naesmd%ay
        sim%naesmd%azold=sim%naesmd%az
        do i=1,sim%naesmd%natom
            if(sim%naesmd%ensemble.eq.'langev') then
                sim%naesmd%rx(i)=sim%naesmd%rx(i)+sim%naesmd%vx(i)*sim%naesmd%vfric(i)+sim%naesmd%ax(i)*sim%naesmd%afric(i)  &
                    +sim%naesmd%prand(1,i)
                sim%naesmd%ry(i)=sim%naesmd%ry(i)+sim%naesmd%vy(i)*sim%naesmd%vfric(i)+sim%naesmd%ay(i)*sim%naesmd%afric(i) &
                    +sim%naesmd%prand(2,i)
                sim%naesmd%rz(i)=sim%naesmd%rz(i)+sim%naesmd%vz(i)*sim%naesmd%vfric(i)+sim%naesmd%az(i)*sim%naesmd%afric(i) &
                    +sim%naesmd%prand(3,i)

                sim%naesmd%vx(i)=sim%naesmd%vx(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%ax(i)*sim%naesmd%vfric(i)
                sim%naesmd%vy(i)=sim%naesmd%vy(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%ay(i)*sim%naesmd%vfric(i)
                sim%naesmd%vz(i)=sim%naesmd%vz(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%az(i)*sim%naesmd%vfric(i)

            end if

            if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
                ! no sim%naesmd%friction or random force, ballistic propagation
                sim%naesmd%rx(i)=sim%naesmd%rx(i)+sim%naesmd%vx(i)*sim%naesmd%dtmdqt+sim%naesmd%ax(i)*dt2_2
                sim%naesmd%ry(i)=sim%naesmd%ry(i)+sim%naesmd%vy(i)*sim%naesmd%dtmdqt+sim%naesmd%ay(i)*dt2_2
                sim%naesmd%rz(i)=sim%naesmd%rz(i)+sim%naesmd%vz(i)*sim%naesmd%dtmdqt+sim%naesmd%az(i)*dt2_2

                sim%naesmd%vx(i)=sim%naesmd%vx(i)+sim%naesmd%ax(i)*dt_2
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+sim%naesmd%ay(i)*dt_2
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+sim%naesmd%az(i)*dt_2
            end if
        end do
        !
        !--------------------------------------------------------------------
        !
        !     Get the potential energy, atomic forces and accelerations
        !     Newton second law to get the next accelerations;
        !     the accelerations and forces are calculated at deriv.f
        if((sim%cosmo%solvent_model.eq.4).or.(sim%cosmo%solvent_model.eq.5)) then
            call calc_cosmo_4(sim)
        else
            if(sim%excn>0) then
                call do_sqm_davidson_update(sim,sim%naesmd%cmdqt,sim%naesmd%vmdqt,sim%naesmd%vgs, &
                    statelimit=sim%excn+1)
            else
                call do_sqm_davidson_update(sim,sim%naesmd%cmdqt,sim%naesmd%vmdqt,sim%naesmd%vgs, &
                    statelimit=sim%qmmm%state_of_interest)
            end if
        endif
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call cpu_time(t_start)
        call deriv(sim,sim%naesmd%ihop)
        call cpu_time(t_finish)
        sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
        call deriv2naesmd_accel(sim)
        !
        !
        !     find the full-step velocities using the Verlet recursion
        do i=1,sim%naesmd%natom
            if(sim%naesmd%ensemble.eq.'langev') then
                sim%naesmd%vx(i)=sim%naesmd%vx(i)+0.5d0*sim%naesmd%ax(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(1,i)
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+0.5d0*sim%naesmd%ay(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(2,i)
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+0.5d0*sim%naesmd%az(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(3,i)
            end if

            if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
                ! ballistic
                sim%naesmd%vx(i)=sim%naesmd%vx(i)+sim%naesmd%ax(i)*dt_2
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+sim%naesmd%ay(i)*dt_2
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+sim%naesmd%az(i)*dt_2
            end if
        end do

        !remove rotation/translation by rescaling velocities
        call rescaleveloc(sim%naesmd%rx,sim%naesmd%ry,sim%naesmd%rz,&
		sim%naesmd%vx,sim%naesmd%vy,sim%naesmd%vz,sim%naesmd%massmdqt,sim%naesmd%natom)

        !
        !  kinetic energy calculation
        !
        sim%naesmd%kin=0.d0

        do i=1,sim%naesmd%natom
            sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(i)*(sim%naesmd%vx(i)**2+sim%naesmd%vy(i)**2+sim%naesmd%vz(i)**2)/2
        end do


        return
    end subroutine
    !
    !--------------------------------------------------------------------
    !
    subroutine verlet1(sim)

        implicit none

        type(simulation_t),pointer :: sim

        integer Na,Nm
        integer i
        _REAL_ dt_2,dt2_2


        dt_2=sim%naesmd%dtmdqt/2.d0 ! dt/2
        dt2_2=sim%naesmd%dtmdqt**2/2 ! dt^2/2
        !
        !--------------------------------------------------------------------
        !
        !  Get frictional and random terms for position and velocity
        !
        !--------------------------------------------------------------------
        !
        call sdterm(sim%naesmd)
        !
        !--------------------------------------------------------------------
        !
        !  Store the current atom positions, then find new atom
        !  positions and half-step velocities via Verlet recursion
        !
        !--------------------------------------------------------------------
        !
        sim%naesmd%rxold=sim%naesmd%rx
        sim%naesmd%ryold=sim%naesmd%ry
        sim%naesmd%rzold=sim%naesmd%rz
        sim%naesmd%vxold=sim%naesmd%vx
        sim%naesmd%vyold=sim%naesmd%vy
        sim%naesmd%vzold=sim%naesmd%vz
        sim%naesmd%axold=sim%naesmd%ax
        sim%naesmd%ayold=sim%naesmd%ay
        sim%naesmd%azold=sim%naesmd%az

        do i=1,sim%naesmd%natom
            if(sim%naesmd%ensemble.eq.'langev') then
                sim%naesmd%rx(i)=sim%naesmd%rx(i)+sim%naesmd%vx(i)*sim%naesmd%vfric(i)+sim%naesmd%ax(i)*sim%naesmd%afric(i)  &
                    +sim%naesmd%prand(1,i)
                sim%naesmd%ry(i)=sim%naesmd%ry(i)+sim%naesmd%vy(i)*sim%naesmd%vfric(i)+sim%naesmd%ay(i)*sim%naesmd%afric(i) &
                    +sim%naesmd%prand(2,i)
                sim%naesmd%rz(i)=sim%naesmd%rz(i)+sim%naesmd%vz(i)*sim%naesmd%vfric(i)+sim%naesmd%az(i)*sim%naesmd%afric(i) &
                    +sim%naesmd%prand(3,i)

                sim%naesmd%vx(i)=sim%naesmd%vx(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%ax(i)*sim%naesmd%vfric(i)
                sim%naesmd%vy(i)=sim%naesmd%vy(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%ay(i)*sim%naesmd%vfric(i)
                sim%naesmd%vz(i)=sim%naesmd%vz(i)*sim%naesmd%pfric(i)+0.5d0*sim%naesmd%az(i)*sim%naesmd%vfric(i)

            end if

            if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
                ! no sim%naesmd%friction or random force, ballistic propagation
                sim%naesmd%rx(i)=sim%naesmd%rx(i)+sim%naesmd%vx(i)*sim%naesmd%dtmdqt+sim%naesmd%ax(i)*dt2_2
                sim%naesmd%ry(i)=sim%naesmd%ry(i)+sim%naesmd%vy(i)*sim%naesmd%dtmdqt+sim%naesmd%ay(i)*dt2_2
                sim%naesmd%rz(i)=sim%naesmd%rz(i)+sim%naesmd%vz(i)*sim%naesmd%dtmdqt+sim%naesmd%az(i)*dt2_2

                sim%naesmd%vx(i)=sim%naesmd%vx(i)+sim%naesmd%ax(i)*dt_2
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+sim%naesmd%ay(i)*dt_2
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+sim%naesmd%az(i)*dt_2
            end if
        end do

        call flush(130)
        call flush(131)
        call flush(132)

        return
    end subroutine
    !
    !--------------------------------------------------------------------
    !
    subroutine verlet2(sim)
        implicit none

        type(simulation_t),pointer :: sim

        integer Na,Nm
        integer i
        _REAL_ dt_2,dt2_2
        _REAL_ t_start,t_finish


        dt_2=sim%naesmd%dtmdqt/2.d0 ! dt/2
        dt2_2=sim%naesmd%dtmdqt**2/2 ! dt^2/2
        !
        !--------------------------------------------------------------------
        !
        !     Get the potential energy, atomic forces and accelerations
        !     Newton second law to get the next accelerations;
        !     the forces are calculated by deriv

        if(sim%excn>0) then
            call do_sqm_davidson_update(sim,sim%naesmd%cmdqt,sim%naesmd%vmdqt,sim%naesmd%vgs, &
                statelimit=sim%excn+1)
        else
            call do_sqm_davidson_update(sim,sim%naesmd%cmdqt,sim%naesmd%vmdqt,sim%naesmd%vgs, &
                statelimit=sim%qmmm%state_of_interest)
        end if

        call cpu_time(t_start)
        call deriv(sim,sim%naesmd%ihop)
        call cpu_time(t_finish)
        sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
 
        call deriv2naesmd_accel(sim)

        !  find the full-step velocities using the Verlet recursion
        do i=1,sim%naesmd%natom
            if(sim%naesmd%ensemble.eq.'langev') then
                sim%naesmd%vx(i)=sim%naesmd%vx(i)+0.5d0*sim%naesmd%ax(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(1,i)
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+0.5d0*sim%naesmd%ay(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(2,i)
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+0.5d0*sim%naesmd%az(i)*sim%naesmd%vfric(i)+sim%naesmd%vrand(3,i)

            end if

            if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
                ! ballistic
                sim%naesmd%vx(i)=sim%naesmd%vx(i)+sim%naesmd%ax(i)*dt_2
                sim%naesmd%vy(i)=sim%naesmd%vy(i)+sim%naesmd%ay(i)*dt_2
                sim%naesmd%vz(i)=sim%naesmd%vz(i)+sim%naesmd%az(i)*dt_2
            end if
        end do

        !remove translation/rotation
        call rescaleveloc(sim%naesmd%rx,sim%naesmd%ry,sim%naesmd%rz, &
		sim%naesmd%vx,sim%naesmd%vy,sim%naesmd%vz,sim%naesmd%massmdqt,sim%naesmd%natom)


        call flush(133)
        !
        !     kinetic energy calculation
        !
        sim%naesmd%kin=0.d0

        do i=1,sim%naesmd%natom
            sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(i)*(sim%naesmd%vx(i)**2+sim%naesmd%vy(i)**2+sim%naesmd%vz(i)**2)/2
        end do


        return
    end subroutine
!
end module
!
