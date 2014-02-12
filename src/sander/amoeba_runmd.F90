#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_runmd
implicit none
private
integer, save :: do_flag

_REAL_,save, allocatable :: acceleration(:,:),old_acceleration(:,:), &
                            old_crd(:,:)
_REAL_,parameter :: kcal_to_gm_x_Ang2_per_ps2 = 4.184d+2
_REAL_,parameter :: ideal_gas_const_kcal_per_mol_K = 1.9872065d-3
_REAL_, parameter :: sander_to_tinker_time_convert = 20.455d0
_REAL_, parameter :: tinker_to_sander_time_convert = 1.d0 / 20.455d0

public AM_RUNMD_get_coords,AM_RUNMD_init,AM_RUNMD,AM_RUNMD_scale_cell, &
       AM_RUNMD_get_ucell_info

contains
!------------------------------------------------------------
subroutine AM_RUNMD_init(natom)
   use amoeba_mdin, only : beeman_integrator
   integer,intent(in) :: natom

   integer :: ier
   if ( beeman_integrator == 1 )then
      allocate(acceleration(3,natom),stat=ier)
      REQUIRE( ier==0 )
      allocate(old_acceleration(3,natom),stat=ier)
      REQUIRE( ier==0 )
      allocate(old_crd(3,natom),stat=ier)
      REQUIRE( ier==0 )
   endif
end subroutine AM_RUNMD_init
!------------------------------------------------------------
subroutine AM_RUNMD_get_coords(natom,tt,irest,ntb,crd,vel)
   use amoeba_mdin, only : beeman_integrator,am_nbead
   use nblist, only: a,b,c,alpha,beta,gamma
   use file_io_dat
   integer,intent(in) :: natom,irest,ntb
   _REAL_,intent(out) :: tt,crd(3,*),vel(3,*)

!  (DAC note: this routine is only called once we have determined that we
!  are reading a "new-style" coordinate file.  Therefore, earlier checks
!  on the value of "iok" can be (and have been) removed.)

   integer dim1,num_list,j,k,nf,ibead
   integer :: iok,ionerr
   character(len=80) :: fmt,fmtin,dtype

   ! open inpcrd file
   nf = 9  !auto assigned at some point??
   call amopen(nf,inpcrd,'O','F','R')
   call nxtsec_crd_reset()
   fmtin = '(E16.8)'
   dtype = 'ATOMIC_COORDS_SIMULATION_TIME'
   ionerr = 1 ! not fatal if missing ---should be there
   call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)tt
   fmtin = '(I8)'
   dtype = 'ATOMIC_COORDS_NUM_LIST'
   ionerr = 1 ! not fatal if missing ---should be there
   call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)num_list

   
   if ( num_list == natom*am_nbead )then
      dim1 = 3
      ionerr = 0 !fatal if missing
      fmtin = '(5E16.8)'
      dtype = 'ATOMIC_COORDS_LIST'
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)((crd(j,k),j=1,dim1),k=1,num_list)
   else if( num_list == natom ) then
      dim1 = 3
      ionerr = 0 !fatal if missing
      fmtin = '(5E16.8)'
      dtype = 'ATOMIC_COORDS_LIST'
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)((crd(j,k),j=1,dim1),k=1,num_list)
      do ibead=1,am_nbead-1
         crd(1:3,ibead*natom+1:ibead*natom+natom)=crd(1:3,1:natom)
      end do
   else
      write(6,*)'AMOEBA_get_coords: wrong number of atomic coords'
      call mexit(6,1)
   endif
 
   if ( irest /= 0 )then !look for velocities
      fmtin = '(I8)'
      dtype = 'ATOMIC_VELOCITIES_NUM_LIST'
      ionerr = 1 ! not fatal if missing ---should be there
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)num_list
      if ( num_list /= natom )then
         write(6,*)'AMOEBA_get_coords: wrong number of atomic velocities'
         call mexit(6,1)
      endif
      dim1 = 3
      ionerr = 0 !fatal if missing
      fmtin = '(5E16.8)'
      dtype = 'ATOMIC_VELOCITIES_LIST'
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)((vel(j,k),j=1,dim1),k=1,num_list)
      do k = 1,num_list
         do j = 1,3
            vel(j,k) = sander_to_tinker_time_convert*vel(j,k)
         enddo
      enddo
   endif ! irest /= 0
   
   if ( irest /= 0 .and. beeman_integrator == 1 )then !find accel, old_accel
      ! first the accelerations
      fmtin = '(I8)'
      dtype = 'ATOMIC_ACCELERATIONS_NUM_LIST'
      ionerr = 1 ! not fatal if missing ---should be there
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)num_list
      if ( num_list /= natom )then
         write(6,*)'AMOEBA_get_coords: wrong number of atomic accelerations'
         call mexit(6,1)
      endif
      dim1 = 3
      ionerr = 0 !fatal if missing
      fmtin = '(5E16.8)'
      dtype = 'ATOMIC_ACCELERATIONS_LIST'
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)((acceleration(j,k),j=1,dim1),k=1,num_list)
      do k = 1,num_list
         do j = 1,3
            acceleration(j,k) = sander_to_tinker_time_convert**2* &
                                acceleration(j,k)
         enddo
      enddo

      ! next the old accelerations
      fmtin = '(I8)'
      dtype = 'OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
      ionerr = 1 ! not fatal if missing ---should be there
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)num_list
      if ( num_list /= natom )then
         write(6,*)'AMOEBA_get_coords: ', &
                   'wrong number of old atomic accelerations'
         call mexit(6,1)
      endif
      dim1 = 3
      ionerr = 0 !fatal if missing
      fmtin = '(5E16.8)'
      dtype = 'OLD_ATOMIC_ACCELERATIONS_LIST'
      call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
      read(nf,fmt)((old_acceleration(j,k),j=1,dim1),k=1,num_list)
      do k = 1,num_list
         do j = 1,3
            old_acceleration(j,k) = sander_to_tinker_time_convert**2* &
                                    old_acceleration(j,k)
         enddo
      enddo
   endif ! irest /= 0 .and. beeman_integrator == 1
   fmtin = '(3e16.8)'
   dtype = 'UNIT_CELL_PARAMETERS'
   ionerr = 1 ! not fatal if missing
   call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( ntb > 0 )then
      read(nf,fmt)a,b,c,alpha,beta,gamma
      call fill_ucell(a,b,c,alpha,beta,gamma)
      write(6,9129)a,b,c,alpha,beta,gamma
   endif
   call nxtsec_crd_reset()
   close(nf)
   9129 format(t2,'NEW EWALD BOX PARAMETERS from inpcrd file:', &
         /5x,'A     =',f10.5,'  B    =',f10.5,'  C     =',f10.5,/, &
         /5x,'ALPHA =',f10.5,'  BETA =',f10.5,'  GAMMA =',f10.5,/)

end subroutine AM_RUNMD_get_coords
!------------------------------------------------------------
subroutine AM_RUNMD_write_restrt(natom,t,title,restrt,crd,vel,owrite)
   use nblist, only: a,b,c,alpha,beta,gamma
   integer, intent(in) :: natom
   _REAL_,intent(in) :: t
   character(len=*),intent(in) :: title,restrt
   character, intent(in)       :: owrite
   _REAL_,intent(in) :: crd(3,*),vel(3,*)

   logical first
   save first
   data first/.true./
   integer :: values(8),nf,j,n
   character(len=12) :: date,time,zone
   character(len=8)word,word1
   character(len=16)word2

   nf = 16
   if( first ) then
      call amopen(nf,restrt,owrite,'F','W')
      first = .false.
   else
      call amopen(nf,restrt,'O','F','W')
   end if

   call date_and_time(date,time,zone,values)
   write(nf,'(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
           '%VERSION  VERSION_STAMP = V0001.000  DATE = ', &
           values(2),'/',values(3),'/',values(1)-2000, &
           '  ',values(5),':',values(6),':',values(7)
   write(nf,'(a)')'%FLAG TITLE'
   write(nf,'(a)')'%FORMAT(a)'
   write(nf,'(a)')title(1:len_trim(title))
   write(nf,'(a)')'%FLAG ATOMIC_COORDS_SIMULATION_TIME'
   write(nf,'(a)')'%FORMAT(E16.8)'
   write(nf,'(E16.8)')t
!-----------------------------------------------------------
   write(nf,'(a)')'%FLAG ATOMIC_COORDS_NUM_LIST'
   write(nf,'(a)')'%FORMAT(i8)'
   write(nf,'(I8)')natom
   write(word,'(I8)')natom
   word1 = adjustl(word)
   word2 = '(3,'//word1(1:len_trim(word1))//')'
   write(nf,'(a)')'%FLAG ATOMIC_COORDS_LIST'
   write(nf,'(a)')'%COMMENT   dimension = '//word2
   write(nf,'(a)')'%FORMAT(3e20.12)'
   write(nf,'(3e20.12)')((crd(j,n),j=1,3),n=1,natom)
   write(nf,'(a)')'%FLAG ATOMIC_VELOCITIES_NUM_LIST'
   write(nf,'(a)')'%FORMAT(i8)'
   write(nf,'(I8)')natom
   write(nf,'(a)')'%FLAG ATOMIC_VELOCITIES_LIST'
   write(nf,'(a)')'%COMMENT   dimension = '//word2
   write(nf,'(a)')'%FORMAT(3e20.12)'
   write(nf,'(3e20.12)') &
          ((tinker_to_sander_time_convert*vel(j,n),j=1,3),n=1,natom)
   write(nf,'(a)')'%FLAG ATOMIC_ACCELERATIONS_NUM_LIST'
   write(nf,'(a)')'%FORMAT(i8)'
   write(nf,'(I8)')natom
   write(nf,'(a)')'%FLAG ATOMIC_ACCELERATIONS_LIST'
   write(nf,'(a)')'%COMMENT   dimension = '//word2
   write(nf,'(a)')'%FORMAT(3e20.12)'
   write(nf,'(3e20.12)') &
           ((tinker_to_sander_time_convert**2*acceleration(j,n),j=1,3), &
                                                                n=1,natom)
   write(nf,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
   write(nf,'(a)')'%FORMAT(i8)'
   write(nf,'(I8)')natom
   write(nf,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_LIST'
   write(nf,'(a)')'%COMMENT   dimension = '//word2
   write(nf,'(a)')'%FORMAT(3e20.12)'
   write(nf,'(3e20.12)') &
           ((tinker_to_sander_time_convert**2*old_acceleration(j,n), &
                                                          j=1,3),n=1,natom)
!-----------------------------------------------------------
  write(nf,'(a)')'%FLAG UNIT_CELL_PARAMETERS'
  write(nf,'(a,a)')'%COMMENT lengths a,b,c; ', &
                          'then angles alpha,beta,gamma'
  write(nf,'(a)')'%FORMAT(3e20.12)'
  write(nf,'(3e20.12)')a,b,c,alpha,beta,gamma
  close(nf)
end subroutine AM_RUNMD_write_restrt
!------------------------------------------------------------
subroutine AM_RUNMD( ix,ih,ipairs, &
                     winv,mass,xx, &
                     crd,vel,frc,qsetup)
   use nblist, only: volume
   use file_io_dat
   use state
   integer,intent(in) :: ix(*),ipairs(*)
    character(len=4),intent(in) :: ih(*)
   _REAL_,intent(in) :: winv(*),mass(*),xx(*)
   _REAL_,intent(inout) :: crd(3,*),vel(3,*),frc(3,*)
   logical,intent(in) :: qsetup

#  include "md.h"
#  include "memory.h"
   integer :: nstep,nitp,nits,nren
   _REAL_          :: vir(4),onefac(3),degrees_of_freedom
   type(state_rec) :: ener,enert,enert2,enert_old,enert2_old
   logical do_list_update

   degrees_of_freedom = dble(3*natom - 3)  !num crds minus overall translation
                                           ! FIXME change if constraints
   onefac(1) = 2.d0 / (degrees_of_freedom*ideal_gas_const_kcal_per_mol_K)
   onefac(2) = 1.d-6
   onefac(3) = 1.d-6
   ! clear the energy array
   ener         = null_state_rec
   ener%volume  = volume
   ener%density = tmass / (0.602204d0*volume)
   nren = 51
   do_list_update = .FALSE.  !only important for psander--> fix later

   nitp = 0
   nits = 0
   enert      = null_state_rec
   enert2     = null_state_rec
   enert_old  = null_state_rec
   enert2_old = null_state_rec
   if ( ntpr <= 0 )ntpr = nstlim
   call amopen(7,mdinfo,'U','F','W')
   do nstep = 1,nstlim
      call AM_RUNMD_beeman_md_step(&
                                   ix,ih,ipairs,ntt,ntp, &
                                   dt,winv,mass,xx,temp0,tautp, &
                                   pres0,taup, &
                                   crd,vel, &
                                   frc,vir,ener,volume,tmass, &
                                   qsetup,do_list_update,nstep)
      call AM_RUNMD_update_enert(nren,ener,enert,enert2)
      t = t + dt
      if ( nstep == 1 .or. nstep == nstlim .or. mod(nstep,ntpr) == 0 )then
         rewind(7)
         call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)
         call amflsh(7)
      endif
      if ( ntwr > 0 )then
         if ( mod(nstep,ntwr)==0 )then
            call AM_RUNMD_write_restrt(natom,t,title,restrt,crd,vel,owrite)
         endif
      endif
      if ( nscm > 0 )then
        if ( mod(nstep,nscm) == 0 )then
           call AM_RUNMD_remove_com(natom,mass,tmass,onefac,vel)
        endif
      endif
      if ( ntwx > 0 )then
         if ( mod(nstep,ntwx)==0 )then
            call AM_RUNMD_crd_archive(natom,crd,MDCRD_UNIT)
         endif
      endif
      if ( ntwv > 0 )then
         if ( mod(nstep,ntwv)==0 )then
            call AM_RUNMD_vel_archive(natom,vel,MDVEL_UNIT)
         endif
      endif
      if ( ntave > 0 )then
        if ( mod(nstep,ntave) == 0 )then
           call AM_RUNMD_running_ave(nren,nstep,ntave,onefac,t, &
                  enert,enert2,enert_old,enert2_old)
        endif
      endif
   enddo
   call AM_RUNMD_overall_ave(nren,nstlim,onefac,t,enert,enert2)
   call AM_RUNMD_write_restrt(natom,t,title,restrt,crd,vel,owrite)
   
end subroutine AM_RUNMD
!------------------------------------------------------------
subroutine AM_RUNMD_beeman_md_step( &
                                   ix,ih,ipairs,ntt,ntp, &
                                   delt,winv,mass,xx,temp0,tautp, &
                                   pres0,taup, &
                                   crd,vel, &
                                   frc,vir,ener,volume,tmass, &
                                   qsetup,do_list_update,nstep)
   use state
   use nblist, only: fill_tranvec
   integer,intent(in) :: ix(*),ipairs(*),ntt,ntp,nstep
   character(len=4),intent(in) :: ih(*)
   _REAL_,intent(in) :: delt,winv(*),mass(*),xx(*), &
                        temp0,tautp,pres0,taup,tmass
   _REAL_,intent(inout) :: crd(3,*),vel(3,*)
   _REAL_,intent(out)   :: frc(3,*),vir(4)
   type(state_rec),intent(out) :: ener
   _REAL_ :: volume
   logical,intent(in) :: qsetup,do_list_update

   _REAL_ :: kinetic_energy(3,3),eksum,degrees_of_freedom, &
             convert,temp,pressure
#  include "memory.h"

   call AM_RUNMD_beeman_1st_crd_update(natom,delt, &
                                 acceleration,old_acceleration, &
                                 old_crd,crd,vel)

! call rattle here if used---future
! get the forces
   call force(xx,ix,ih,ipairs,crd,frc,ener,vir, & 
            xx(l96), xx(l97), xx(l98), xx(l99), &
            qsetup,do_list_update,nstep)
! get the next accelerations from forces
   call AM_RUNMD_beeman_full_step_vels(natom,delt,  &
                                 kcal_to_gm_x_Ang2_per_ps2,winv,  &
                                 acceleration,old_acceleration,frc,vel)
! call rattle here again if used---future
! get kin ene here
   call AM_RUNMD_kinetic_energy_calc(natom,kcal_to_gm_x_Ang2_per_ps2, &
                                     mass,vel,kinetic_energy,eksum)
   ener%kin%tot = eksum
   ener%tot     = ener%kin%tot + ener%pot%tot
   degrees_of_freedom = dble(3*natom - 3)!num crds minus overall translation
                                           ! FIXME change if constraints
   convert = 2.d0 / (degrees_of_freedom*ideal_gas_const_kcal_per_mol_K)
   temp = eksum*convert
   if ( ntt > 0 )then
      call AM_RUNMD_berendson_scaling(natom,vel,delt,temp,temp0,tautp)
   endif
   call AM_RUNMD_pressure(eksum,vir,volume,pressure)
   ener%pres(4) = pressure
   ener%volume  = volume
   ener%density = tmass / (0.602204d0*volume)

   if ( ntp > 0 )then !atom based scaling
      call AM_RUNMD_scale_cell(natom,pressure,delt,pres0,taup,crd)
      call fill_tranvec()
   endif
end subroutine AM_RUNMD_beeman_md_step
!------------------------------------------------------------
subroutine AM_RUNMD_beeman_1st_crd_update(numatoms,delt,accel,old_accel, &
                                 oldcrd,crd,vel)
  implicit none
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: delt,accel(3,*),old_accel(3,*)
  _REAL_,intent(out) :: oldcrd(3,*)
  _REAL_,intent(inout) :: crd(3,*),vel(3,*)
!
!     store the current atom positions, then find new atom
!     positions and half-step velocities via Beeman recursion
!

  integer :: n,j
  _REAL_ :: term(3),delt_8,delt2_8

  delt_8 = delt / 8.0d0
  delt2_8 = delt * delt_8

  do n = 1,numatoms
    do j = 1,3
      oldcrd(j,n) = crd(j,n)
      term(j) = 5.d0*accel(j,n) - old_accel(j,n)
      crd(j,n) = crd(j,n) + vel(j,n)*delt + term(j)*delt2_8
      vel(j,n) = vel(j,n) + term(j)*delt_8
    enddo
  enddo
end subroutine AM_RUNMD_beeman_1st_crd_update
!------------------------------------------------------------
subroutine AM_RUNMD_beeman_full_step_vels(numatoms,delt,factor,winv,  &
                                 accel,old_accel,frc,vel)
  implicit none
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: delt,factor,winv(*),frc(3,*)
  _REAL_,intent(out) :: old_accel(3,*)
  _REAL_,intent(inout) :: vel(3,*),accel(3,*)
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the Beeman recursion
!
  integer :: n,j
  _REAL_ :: delt_8

  delt_8 = delt / 8.0d0
  do n = 1,numatoms
    do j = 1,3
      old_accel(j,n) = accel(j,n)
      accel(j,n) = factor*frc(j,n)*winv(n)
      vel(j,n) = vel(j,n) + (3.d0*accel(j,n) + old_accel(j,n))*delt_8
    enddo
  enddo
  return
end subroutine AM_RUNMD_beeman_full_step_vels
!------------------------------------------------------------
subroutine AM_RUNMD_kinetic_energy_calc( &
                      numatoms,factor,mass,vel,kinetic_energy,eksum)
  implicit none
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: factor,mass(*),vel(3,*)
  _REAL_,intent(out) :: kinetic_energy(3,3),eksum

  integer j,k,n
  _REAL_ term
  do k = 1,3
    do j = 1,3
      kinetic_energy(j,k) = 0.d0
    enddo
  enddo
  do n = 1,numatoms
    do k = 1,3
      term = mass(n)*vel(k,n)
      do j = 1,3
        kinetic_energy(j,k) = kinetic_energy(j,k) + term*vel(j,n)
      enddo
    enddo
  enddo
! convert units & multiply by 1/2
  term = 0.5d0 / factor
  do k = 1,3
    do j = 1,3
      kinetic_energy(j,k) = term*kinetic_energy(j,k)
    enddo
  enddo
  eksum = kinetic_energy(1,1)+kinetic_energy(2,2)+kinetic_energy(3,3)
  return
end subroutine AM_RUNMD_kinetic_energy_calc
!------------------------------------------------------------
subroutine AM_RUNMD_berendson_scaling(numatoms,vel,delt,current_temp, &
                         target_temp,tautp)
   implicit none
   integer,intent(in) :: numatoms
   _REAL_,intent(in) :: delt,current_temp,target_temp,tautp
   _REAL_,intent(inout) :: vel(3,*)

   integer :: n
   _REAL_ :: scale_factor

   scale_factor = sqrt(1.0d0 + (delt/tautp)*(target_temp/current_temp-1.0d0))
   do n = 1,numatoms
     vel(1,n) = scale_factor*vel(1,n)
     vel(2,n) = scale_factor*vel(2,n)
     vel(3,n) = scale_factor*vel(3,n)
   enddo
end subroutine AM_RUNMD_berendson_scaling
!------------------------------------------------------------
subroutine AM_RUNMD_remove_com(numatoms,mass,tmass,onefac,vel)
   implicit none
   integer,intent(in) :: numatoms
   _REAL_,intent(in) :: mass(*),tmass,onefac(3)
   _REAL_,intent(inout) :: vel(3,*)

   integer :: n
   _REAL_ :: vcmx,vcmy,vcmz,amass,tmassinv,vel2,atempdrop

   vcmx = 0.d0
   vcmy = 0.d0
   vcmz = 0.d0
   do n = 1,numatoms
      amass = mass(n)
      vcmx = vcmx + amass*vel(1,n)
      vcmy = vcmy + amass*vel(2,n)
      vcmz = vcmz + amass*vel(3,n)
   enddo
   tmassinv = 1.d0 / tmass
   vcmx = vcmx * tmassinv
   vcmy = vcmy * tmassinv
   vcmz = vcmz * tmassinv
   vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz
   atempdrop = 0.5d0 * tmass * vel2 * onefac(1)
   vel2 = sqrt(vel2)
   write (6,'(a,f15.6,f9.2,a)')'check COM velocity, temp: ',vel2,atempdrop, &
                               '(Removed)'
   do n = 1,numatoms
      vel(1,n) = vel(1,n) - vcmx
      vel(2,n) = vel(2,n) - vcmy
      vel(3,n) = vel(3,n) - vcmz
   enddo
end subroutine AM_RUNMD_remove_com
!------------------------------------------------------------
subroutine AM_RUNMD_crd_archive(numatoms,crd,MDCRD_UNIT)
   implicit none
   integer,intent(in) :: numatoms,MDCRD_UNIT
   _REAL_,intent(in) :: crd(3,*)

#  include "box.h"
   logical :: loutfm

   loutfm = .TRUE.
   call corpac(crd,1,3*numatoms,MDCRD_UNIT,loutfm)
   if ( ntb > 0 )then
     call corpac(box,1,3,MDCRD_UNIT,loutfm)
   endif
end subroutine AM_RUNMD_crd_archive
!------------------------------------------------------------
subroutine AM_RUNMD_vel_archive(numatoms,vel,MDVEL_UNIT)
   implicit none
   integer,intent(in) :: numatoms,MDVEL_UNIT
   _REAL_,intent(in) :: vel(3,*)

   logical :: loutfm
   integer :: ier,j,n
   _REAL_, allocatable :: svel(:,:)

   allocate(svel(3,numatoms),stat=ier)
   REQUIRE( ier==0 )
   loutfm = .TRUE.
   do n = 1,numatoms
      do j = 1,3
         svel(j,n) = tinker_to_sander_time_convert*vel(j,n)
      enddo
   enddo
   call corpac(svel,1,3*numatoms,MDVEL_UNIT,loutfm)
   deallocate(svel)
end subroutine AM_RUNMD_vel_archive
!------------------------------------------------------------
subroutine AM_RUNMD_get_ucell_info(inpcrd,a,b,c,alpha,beta,gam)
   implicit none
   character(len=*),intent(in) :: inpcrd
   _REAL_,intent(out) :: a,b,c,alpha,beta,gam

   integer nf
   integer :: iok,ionerr
   character(len=80) :: fmt,fmtin,dtype

   nf = 9  !auto assigned at some point??
   call amopen(nf,inpcrd,'O','F','R')
   call nxtsec_crd_reset()

   fmtin = '(3e16.8)'
   dtype = 'UNIT_CELL_PARAMETERS'
   ionerr = 1 ! not fatal if missing
   call nxtsec_crd(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found in prmtop
      read(nf,fmt)a,b,c,alpha,beta,gam
   else
      write(6,*)'AMOEBA_get_ucell_info: UNIT_CELL_PARAMETERS missing'
      call mexit(6,1)
   endif
   call nxtsec_crd_reset()
   close(nf)
end subroutine AM_RUNMD_get_ucell_info
!------------------------------------------------------------
subroutine AM_RUNMD_update_enert(nren,ener,enert,enert2)
   use state
   implicit none
   integer,intent(in) :: nren
   type(state_rec), intent(in)    :: ener
   type(state_rec), intent(inout) :: enert,enert2

    enert  = enert  + ener
    enert2 = enert2 + (ener * ener)
end subroutine AM_RUNMD_update_enert
!------------------------------------------------------------
subroutine AM_RUNMD_running_ave(nren,nstep,ntave,onefac,t, &
                  enert,enert2,enert_old,enert2_old)
   use state
   implicit none
   integer,intent(in) :: nren,nstep,ntave

   _REAL_  :: onefac(3),t
   type(state_rec),intent(in)    :: enert, enert2
   type(state_rec),intent(inout) :: enert_old, enert2_old

   type(state_rec)               :: enert_tmp,enert2_tmp
   integer :: m

   enert_tmp     = enert - enert_old
   enert2_tmp    = enert2 - enert2_old
   enert_old     = enert
   enert2_old    = enert2
   enert_tmp     = enert_tmp/ntave
   enert2_tmp    = enert2_tmp/ntave - enert_tmp*enert_tmp
   call zero_neg_values_state(enert2_tmp)
   enert2_tmp    = sqrt(enert2_tmp)

   write(6,'(/5x,a,i7,a,/)')' A V E R A G E S   O V E R ',ntave,' S T E P S'
   call prntmd(nstep,0,0,t,enert_tmp,onefac,0,.false.)
   write(6,'(/5x,a,/)')' R M S  F L U C T U A T I O N S'
   call prntmd(nstep,0,0,t,enert2_tmp,onefac,0,.false.)
end subroutine AM_RUNMD_running_ave
!------------------------------------------------------------
subroutine AM_RUNMD_overall_ave(nren,nstlim,onefac,t,enert,enert2)
   use state
   implicit none
   integer,intent(in) :: nren,nstlim
   _REAL_,intent(in) :: onefac(3),t
   type(state_rec) ,intent(inout) :: enert,enert2

   integer :: m

   enert  = enert/nstlim
   enert2 = enert2/nstlim - enert*enert
   call zero_neg_values_state(enert2)
   enert2 =  sqrt(enert2)

   write(6,'(/5x,a,i7,a,/)')' A V E R A G E S   O V E R ',nstlim,' S T E P S'
   call prntmd(nstlim,0,0,t,enert,onefac,0,.false.)
   write(6,'(/5x,a,/)')' R M S  F L U C T U A T I O N S'
   call prntmd(nstlim,0,0,t,enert2,onefac,0,.false.)
end subroutine AM_RUNMD_overall_ave
!------------------------------------------------------------
subroutine AM_RUNMD_pressure(eksum,vir,volume,pressure)
   implicit none
   _REAL_,intent(in) :: eksum,vir(4),volume
   _REAL_,intent(out) :: pressure

   _REAL_,parameter :: pressure_constant = 6.85695d+4

   pressure = (pressure_constant/volume)*(2.d0*eksum - vir(4)) / 3.d0
  
end subroutine AM_RUNMD_pressure
!------------------------------------------------------------
subroutine AM_RUNMD_scale_cell(numatoms,pressure,dt,pres0,taup,crd)
   use nblist, only: oldrecip,ucell,sphere,cutlist
   use constants, only : one,third
   use amoeba_mdin, only : compress
   implicit none
   integer,intent(in) :: numatoms
   _REAL_,intent(in) :: pressure,dt,pres0,taup
   _REAL_,intent(inout) :: crd(3,*)

   _REAL_ :: scale_factor(3),f1,f2,f3
   integer :: n

   scale_factor(1) = (one + (dt*compress/taup)*(pressure-pres0))**third
   scale_factor(2) = scale_factor(1)
   scale_factor(3) = scale_factor(1)
   ! rescale the unit cell
   call redo_ucell(scale_factor)
   do n = 1,numatoms
      !     ...get fractional coordinates with old cell parameters
         
      f1 = crd(1,n)*oldrecip(1,1)+crd(2,n)*oldrecip(2,1)+ &
           crd(3,n)*oldrecip(3,1)
      f2 = crd(1,n)*oldrecip(1,2)+crd(2,n)*oldrecip(2,2)+ &
           crd(3,n)*oldrecip(3,2)
      f3 = crd(1,n)*oldrecip(1,3)+crd(2,n)*oldrecip(2,3)+ &
           crd(3,n)*oldrecip(3,3)
         
      !     ...use these with the new cell parameters to get new
      !     cartesian coordinates
         
      crd(1,n) = f1*ucell(1,1)+f2*ucell(1,2)+f3*ucell(1,3)
      crd(2,n) = f1*ucell(2,1)+f2*ucell(2,2)+f3*ucell(2,3)
      crd(3,n) = f1*ucell(3,1)+f2*ucell(3,2)+f3*ucell(3,3)
   enddo
end subroutine AM_RUNMD_scale_cell
!------------------------------------------------------------
end module amoeba_runmd
