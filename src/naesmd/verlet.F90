#include "dprec.fh"
#include "assert.fh"

!***************************************************
!** Verlet velocity algorithm                     **
!***************************************************
module verlet_module
   use naesmd_constants
   use langevin_temperature
   use communism
   use cosmo_C,only:solvent_model
   use qmmm_module, only:qmmm_struct
   implicit none

   contains
!
!--------------------------------------------------------------------
!
   subroutine verlet(sim)
   implicit none

   type(simulation_t),pointer::sim

   integer Na,Nm,temp
   integer k,i,j

   include 'sizes'
   _REAL_ E0,d
   include 'md.par'
   include 'parH.par'
   !real*8 Omega(Mx_M),fosc(Mx_M)
   include 'md.cmn'
   include 'common'

   _REAL_ dt_2,dt2_2
   _REAL_ t_start,t_finish 

   Na=sim%naesmd%Na
   Nm=sim%naesmd%Nm

   dt_2=dtmdqt/2.d0 ! dt/2
   dt2_2=dtmdqt**2/2 ! dt^2/2
!
!--------------------------------------------------------------------
!
!  Get frictional and random terms for position and velocity
!
!--------------------------------------------------------------------
!
   call sdterm
!
!--------------------------------------------------------------------
!         
!  Store the current atom positions, then find new atom
!  positions and half-step velocities via Verlet recursion
!
!--------------------------------------------------------------------
!
   rxold=rx
   ryold=ry
   rzold=rz
   vxold=vx
   vyold=vy
   vzold=vz
   axold=ax
   ayold=ay
   azold=az
   do i=1,natom
      if(ensemble.eq.'langev') then
         rx(i)=rx(i)+vx(i)*vfric(i)+ax(i)*afric(i)  &
            +prand(1,i)
         ry(i)=ry(i)+vy(i)*vfric(i)+ay(i)*afric(i) &
            +prand(2,i)
         rz(i)=rz(i)+vz(i)*vfric(i)+az(i)*afric(i) &
            +prand(3,i)

         vx(i)=vx(i)*pfric(i)+0.5d0*ax(i)*vfric(i)
         vy(i)=vy(i)*pfric(i)+0.5d0*ay(i)*vfric(i)
         vz(i)=vz(i)*pfric(i)+0.5d0*az(i)*vfric(i)

      end if

      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         ! no friction or random force, ballistic propagation
         rx(i)=rx(i)+vx(i)*dtmdqt+ax(i)*dt2_2
         ry(i)=ry(i)+vy(i)*dtmdqt+ay(i)*dt2_2
         rz(i)=rz(i)+vz(i)*dtmdqt+az(i)*dt2_2

         vx(i)=vx(i)+ax(i)*dt_2
         vy(i)=vy(i)+ay(i)*dt_2
         vz(i)=vz(i)+az(i)*dt_2
      end if
   end do
!
!--------------------------------------------------------------------
!
!     Get the potential energy, atomic forces and accelerations
!     Newton second law to get the next accelerations;
!     the accelerations and forces are calculated at deriv.f

   if((solvent_model.eq.4).or.(solvent_model.eq.5)) then
        call calc_cosmo_4(sim)
   else
        call do_sqm_davidson_update(sim,cmdqt,vmdqt,vgs,statelimit=qmmm_struct%state_of_interest)
   endif
   ihop=qmmm_struct%state_of_interest

   call cpu_time(t_start)
   call deriv(sim,ihop)
   call cpu_time(t_finish)
   sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start

   call deriv2naesmd_accel(sim)
! 
!
!     find the full-step velocities using the Verlet recursion
   do i=1,natom
      if(ensemble.eq.'langev') then
         vx(i)=vx(i)+0.5d0*ax(i)*vfric(i)+vrand(1,i)
         vy(i)=vy(i)+0.5d0*ay(i)*vfric(i)+vrand(2,i)
         vz(i)=vz(i)+0.5d0*az(i)*vfric(i)+vrand(3,i)
      end if

      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         ! ballistic
         vx(i)=vx(i)+ax(i)*dt_2
         vy(i)=vy(i)+ay(i)*dt_2
         vz(i)=vz(i)+az(i)*dt_2
      end if
   end do

   !remove rotation/translation by rescaling velocities
   call rescaleveloc(rx,ry,rz,vx,vy,vz,massmdqt,natom)

!
!  kinetic energy calculation
!
   kin=0.d0

   do i=1,natom
      kin=kin+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)/2
   end do

889   format(10000(1x,f18.10))

   return
   end subroutine
!
!--------------------------------------------------------------------
!
   subroutine verlet1(sim)
   implicit none

   type(simulation_t),pointer :: sim

   integer Na,Nm,temp
   integer k,i,j

   include 'sizes'
   real(8) E0,d
   include 'md.par'
   include 'parH.par'
   !real*8 Omega(Mx_M),fosc(Mx_M)
   include 'md.cmn'
   include 'common'

   real(8) dt_2,dt2_2 

   Na=sim%naesmd%Na
   Nm=sim%naesmd%Nm

   dt_2=dtmdqt/2.d0 ! dt/2
   dt2_2=dtmdqt**2/2 ! dt^2/2
!
!--------------------------------------------------------------------
!
!  Get frictional and random terms for position and velocity
!
!--------------------------------------------------------------------
!
   call sdterm
!
!--------------------------------------------------------------------
!         
!  Store the current atom positions, then find new atom
!  positions and half-step velocities via Verlet recursion
!
!--------------------------------------------------------------------
!
   rxold=rx
   ryold=ry
   rzold=rz
   vxold=vx
   vyold=vy
   vzold=vz
   axold=ax
   ayold=ay
   azold=az

   !write(130,889) (ry(i),rz(i),i=1,natom)
   !write(132,889) (vx(i),vy(i),vz(i),i=1,natom)

   do i=1,natom
      if(ensemble.eq.'langev') then
         rx(i)=rx(i)+vx(i)*vfric(i)+ax(i)*afric(i)  &
            +prand(1,i)
         ry(i)=ry(i)+vy(i)*vfric(i)+ay(i)*afric(i) &
            +prand(2,i)
         rz(i)=rz(i)+vz(i)*vfric(i)+az(i)*afric(i) &
            +prand(3,i)

         vx(i)=vx(i)*pfric(i)+0.5d0*ax(i)*vfric(i)
         vy(i)=vy(i)*pfric(i)+0.5d0*ay(i)*vfric(i)
         vz(i)=vz(i)*pfric(i)+0.5d0*az(i)*vfric(i)

      end if

      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         ! no friction or random force, ballistic propagation
         rx(i)=rx(i)+vx(i)*dtmdqt+ax(i)*dt2_2
         ry(i)=ry(i)+vy(i)*dtmdqt+ay(i)*dt2_2
         rz(i)=rz(i)+vz(i)*dtmdqt+az(i)*dt2_2

         vx(i)=vx(i)+ax(i)*dt_2
         vy(i)=vy(i)+ay(i)*dt_2
         vz(i)=vz(i)+az(i)*dt_2
      end if
   end do

   !write(131,889) (ry(i),rz(i),i=1,natom)
   call flush(130)
   call flush(131)
   call flush(132)

889   format(10000(1x,f18.10))

   return
   end subroutine
!
!--------------------------------------------------------------------
!
   subroutine verlet2(sim)
   implicit none

   type(simulation_t),pointer :: sim

   integer Na,Nm,temp
   integer k,i,j

   include 'sizes'
   _REAL_ E0,d
   include 'md.par'
   include 'parH.par'
   include 'md.cmn'
   include 'common'

   _REAL_ dt_2,dt2_2
   _REAL_ t_start,t_finish 

   Na=sim%naesmd%Na
   Nm=sim%naesmd%Nm

   dt_2=dtmdqt/2.d0 ! dt/2
   dt2_2=dtmdqt**2/2 ! dt^2/2
!
!--------------------------------------------------------------------
!
!     Get the potential energy, atomic forces and accelerations
!     Newton second law to get the next accelerations;
!     the forces are calculated by deriv

   call do_sqm_davidson_update(sim,cmdqt,vmdqt,vgs,statelimit=qmmm_struct%state_of_interest)

   call cpu_time(t_start)
   call deriv(sim,ihop)
   call cpu_time(t_finish)
   sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
 
   call deriv2naesmd_accel(sim)

!  find the full-step velocities using the Verlet recursion
   do i=1,natom
      if(ensemble.eq.'langev') then
         vx(i)=vx(i)+0.5d0*ax(i)*vfric(i)+vrand(1,i)
         vy(i)=vy(i)+0.5d0*ay(i)*vfric(i)+vrand(2,i)
         vz(i)=vz(i)+0.5d0*az(i)*vfric(i)+vrand(3,i)

!            if(i.eq.47.or.i.eq.73) then
!              vx(i)=0.0d0
!              vy(i)=0.0d0
!              vz(i)=0.0d0
!            endif
      end if

      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         ! ballistic
         vx(i)=vx(i)+ax(i)*dt_2
         vy(i)=vy(i)+ay(i)*dt_2
         vz(i)=vz(i)+az(i)*dt_2
      end if
   end do

   !remove translation/rotation
   call rescaleveloc(rx,ry,rz,vx,vy,vz,massmdqt,natom)

   !write(133,889) (vx(i),vy(i),vz(i),i=1,natom)
   call flush(133)
!
!     kinetic energy calculation
!
   kin=0.d0

   do i=1,natom
      kin=kin+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)/2
   end do

889   format(10000(1x,f18.10))

   return
   end subroutine
!
end module
!
