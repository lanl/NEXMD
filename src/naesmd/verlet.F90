#include "dprec.fh"
#include "assert.fh"

!***************************************************
!** Verlet velocity algorithm                     **
!***************************************************
!Na,Nm,d,E0,Omega)
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
write(6,*)'vy=',vy   
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

!         do j=1,natom
!             rx(j)=rx(j)*convl
!             ry(j)=ry(j)*convl
!             rz(j)=rz(j)*convl
!         enddo
! for classical path
!        temp=ihop
!        ihop=0
!      call deriv(Na,Nm,d,tfemto,E0,Omega)
!        ihop=temp
!         do j=1,natom
!             rx(j)=rx(j)/convl
!             ry(j)=ry(j)/convl
!             rz(j)=rz(j)/convl
!         enddo

   ! kav: the following block is subsituted by
   !call do_sqm_and_davidson(sim)

   !call deriv(sim%deriv_forces) 

   !call deriv2naesmd_forces(sim)
   !call dav2naesmd_Omega(sim)

   !vgs=sim%naesmd%E0

   !do i=1,npot
      !vmdqt(i)=sim%naesmd%Omega(i)+vgs
   !end do

   if((solvent_model.eq.4).or.(solvent_model.eq.5)) then
        call calc_cosmo_4(sim)
   else
   	call do_sqm_davidson_update(sim,cmdqt,vmdqt,vgs)
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
!
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

   !remove rotation/translation by rescaling velocities
   call rescaleveloc(rx,ry,rz,vx,vy,vz,massmdqt,natom)

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
write(6,*)'vy1:',vy
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
!     Get the potential energy, atomic forces and accelerations
!     Newton second law to get the next accelerations;
!     the accelerations and forces are calculated at deriv.f

!         do j=1,natom
!             rx(j)=rx(j)*convl
!             ry(j)=ry(j)*convl
!             rz(j)=rz(j)*convl
!         enddo
! for classical path
!        temp=ihop
!        ihop=0
!      call deriv(Na,Nm,d,tfemto,E0,Omega)
!        ihop=temp
!         do j=1,natom
!             rx(j)=rx(j)/convl
!             ry(j)=ry(j)/convl
!             rz(j)=rz(j)/convl
!         enddo
write(6,*)'vy2:',vy
   ! kav: this was commented out for some reason
   ! by Grisha, but without it does not work
   call do_sqm_davidson_update(sim,cmdqt,vmdqt,vgs)

   call cpu_time(t_start)
   call deriv(sim,ihop)
   call cpu_time(t_finish)
   sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
 
   call deriv2naesmd_accel(sim)

   ! kav: I do not understand how the following works,
   ! this subroutine is defined in deriv.F90 with multiple
   ! arguments, not just one used here   

   !call dav2naesmd_Omega(sim)

   !vgs=sim%naesmd%E0

   !do i=1,npot
   !   vmdqt(i)=sim%naesmd%Omega(i)+vgs
   !end do
! 
!
!     find the full-step velocities using the Verlet recursion
!
   !write(125,889) tfemto,(ax(i),i=1,natom)
   !stop ! kav: checking

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
