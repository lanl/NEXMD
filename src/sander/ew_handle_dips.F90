! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!-------------------------------------------------------------------
! HANDLE INDUCED
! main driver for induced dipoles
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine handle_induced here]
subroutine handle_induced(crd,numatoms,iac,ico,charge, &
      cn1,cn2,asol,bsol,eelt,epol,frc,x,ix,ipairs, &
!      pol, &
! Modified by WJM, YD
      pol, dampfactor,pol2, &
!
      xr,virvsene,ipres,ibgwat,nres,aveper,aveind,avetot,emtot, &
      dipole_rms, dipiter,dipole_temp,dt,ntb &
      ,cn3,cn4,cn5 &
      )

   implicit none
   integer numatoms,iac(*),ico(*),ibgwat,nres,ipres(*)
   _REAL_ crd(3,*),charge(*),cn1(*),cn2(*),asol(*), &
         bsol(*),eelt,epol,frc(3,*),xr(3,*),pol(*),virvsene, &
         aveper,aveind,avetot,emtot, &
         dipole_rms, dipiter,dipole_temp,dt
   _REAL_ x(*)
! Modified by WJM, YD
   _REAL_ dampfactor(*), pol2(*)
!
   _REAL_ cn3(*),cn4(*),cn5(*)
   integer ix(*),ipairs(*), ntb
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "ew_frc.h"

#  include "def_time.h"

   integer issetup,iter,initdip
   save issetup,initdip
   data issetup/1/
   data initdip/1/

   iter = 1
   if ( issetup == 1 )then
      if ( indmeth >= 0 .and. indmeth <= 2 )then
         call zero_array(x(leold1),3*numatoms)
      end if
      if ( indmeth >= 1 .and. indmeth <= 2 )then
         call zero_array(x(leold2),3*numatoms)
      end if
      if ( indmeth == 2 )then
         call zero_array(x(leold3),3*numatoms)
      end if
      if ( indmeth == 3 .and. irstdip == 0 )then
         call zero_array(x(linddip),3*numatoms)
         call zero_array(x(ldipvel),3*numatoms)
      end if
      iquench = 0
      issetup = 0
   end if
   if ( indmeth == 3 .and. irstdip == 1 )initdip = 0
   if ( indmeth == 3 .and. nquench > 0 )then
      if ( mod(iquench,nquench) == 0 )initdip = 1
   end if
 
   if ( indmeth < 3 .or. initdip == 1 )then
      
      if ( indmeth < 3 )then
         
         !         ---predict induced dipoles from previous time steps:
         
         call dip_init(numatoms,x(linddip),pol, &
               x(leold1),x(leold2),x(leold3),indmeth)
      end if
      
      do iter = 1,maxiter
         call ewald_force(crd,numatoms,iac,ico,charge, &
               cn1,cn2,asol,bsol,eelt,epol,frc,x,ix,ipairs, &
!               xr,virvsene,pol,.false.)
! Modified by WJM, YD
                xr,virvsene,pol, dampfactor, pol2, .false. &
!
                ,cn3,cn4,cn5 &
                )
         if ( verbose >= 1 ) write(6,66)iter,diprms,diptol
         66 format(1x,'Induced dipoles: iter = ',i3,' RMS = ', &
               e12.3, 'TOL = ',e12.3)
         if ( diprms < diptol )goto 100
         call dip_iter(numatoms,x(linddip),pol,x(lfield))
      end do
      100 continue
      initdip = 0
      if ( indmeth < 3 )then
         call dip_fin(numatoms,x(lfield), &
               x(leold1),x(leold2),x(leold3),indmeth)
      else
         call zero_array(x(ldipvel),3*numatoms)
      end if
      
   else
      
      !         ---Car-Parrinello method. No iterations
      
      call ewald_force(crd,numatoms,iac,ico,charge, &
            cn1,cn2,asol,bsol,eelt,epol,frc,x,ix,ipairs, &
!            xr,virvsene,pol,.false.)
! Modified by WJM, YD
             xr,virvsene,pol,dampfactor,pol2,.false.&
!
             ,cn3,cn4,cn5 &
            )
   end if  ! ( indmeth < 3 .or. initdip == 1 )
   if ( indmeth == 3 .and. nquench > 0 ) iquench = iquench + 1
   dipole_rms = diprms
   dipiter = iter
   dipole_temp = diptemp
   call dipole_stat(crd,charge,x(linddip),ibgwat,nres, &
         ipres,aveper,aveind,avetot,emtot)
   
   return
end subroutine handle_induced 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dip_init here]
subroutine dip_init(numatoms,ind_dip,pol, &
      eold1,eold2,eold3,indmeth)
   implicit none
   integer numatoms,indmeth
   _REAL_ eold1(3,*),eold2(3,*),eold3(3,*)
   _REAL_ ind_dip(3,*),pol(*)
   integer n,istart,iend
#ifdef MPI
#  include "parallel.h"
#endif
#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = numatoms
#endif
   
   ! first order prediction
   ! eold1 is field from last MD step, eold2 from 2 MD steps ago
   ! eold3 is from 3 MD steps ago
   
   call timer_start(TIME_DIPUP)
   if ( indmeth == 0 )then
      do n = istart,iend
         ind_dip(1,n) = pol(n)*eold1(1,n)
         ind_dip(2,n) = pol(n)*eold1(2,n)
         ind_dip(3,n) = pol(n)*eold1(3,n)
      end do
   else if ( indmeth == 1 )then
      do n = istart,iend
         ind_dip(1,n) = pol(n)*(2.d0*eold1(1,n) - eold2(1,n))
         ind_dip(2,n) = pol(n)*(2.d0*eold1(2,n) - eold2(2,n))
         ind_dip(3,n) = pol(n)*(2.d0*eold1(3,n) - eold2(3,n))
      end do
   else if ( indmeth == 2 )then
      do n = istart,iend
         ind_dip(1,n) = pol(n)*(3.d0*eold1(1,n) - &
               3.d0*eold2(1,n) + eold3(1,n))
         ind_dip(2,n) = pol(n)*(3.d0*eold1(2,n) - &
               3.d0*eold2(2,n) + eold3(2,n))
         ind_dip(3,n) = pol(n)*(3.d0*eold1(3,n) - &
               3.d0*eold2(3,n) + eold3(3,n))
      end do
   end if
   call timer_stop(TIME_DIPUP)

   return
end subroutine dip_init 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dip_iter here]
subroutine dip_iter(numatoms,ind_dip,pol,field)
   implicit none
   integer numatoms
   _REAL_ field(3,*)
   _REAL_ ind_dip(3,*),pol(*)
   integer n,istart,iend
#ifdef MPI
#  include "parallel.h"
#endif

   _REAL_ maxf,es,maxd,ds
   integer nf,nd

   nf = 1
   nd = 1
   maxf = 0.d0
   maxd = 0.d0
#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = numatoms
#endif
   call timer_start(TIME_DIPUP)
   do n = istart,iend
      ind_dip(1,n) = pol(n)*field(1,n)
      ind_dip(2,n) = pol(n)*field(2,n)
      ind_dip(3,n) = pol(n)*field(3,n)
   end do
   call timer_stop(TIME_DIPUP)
   return
end subroutine dip_iter 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dip_fin here]
subroutine dip_fin(numatoms,field,eold1,eold2,eold3,indmeth)
   implicit none
   integer numatoms,indmeth
   _REAL_ eold1(3,*),eold2(3,*),eold3(3,*)
   _REAL_ field(3,*)
   integer n,istart,iend
#ifdef MPI
#  include "parallel.h"
#endif
#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = numatoms
#endif
   call timer_start(TIME_DIPUP)
   do n = istart,iend
      if ( indmeth == 2 )then
         eold3(1,n) = eold2(1,n)
         eold3(2,n) = eold2(2,n)
         eold3(3,n) = eold2(3,n)
      end if
      if ( indmeth >= 1 .and. indmeth <= 2)then
         eold2(1,n) = eold1(1,n)
         eold2(2,n) = eold1(2,n)
         eold2(3,n) = eold1(3,n)
      end if
      if ( indmeth >= 0 .and. indmeth <= 2)then
         eold1(1,n) = field(1,n)
         eold1(2,n) = field(2,n)
         eold1(3,n) = field(3,n)
      end if
   end do
   call timer_stop(TIME_DIPUP)
   return
end subroutine dip_fin 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_dipinfo here]
subroutine get_dipinfo(numatoms,pol,field,ind_dip,dip_vel, &
      dipself,dipkine,diprms,dipndf,indmeth)
   implicit none
   integer numatoms
   _REAL_ pol(*),field(3,*),ind_dip(3,*),dip_vel(3,*)
   _REAL_ dipself,dipkine,diprms,dipndf
   _REAL_ small
   integer ic,n,indmeth,istart,iend
#ifdef MPI
#  include "parallel.h"
#endif
   small = 1.d-6
#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = numatoms
#endif

   call timer_start(TIME_DIPUP)
   ic = 0
   dipkine = 0.d0
   diprms = 0.d0
   dipself = 0.d0
   do n = istart,iend
      if ( pol(n) > small )then
         ic = ic + 1
         if ( indmeth == 3 )then
            dipkine = dipkine + dip_vel(1,n)**2 + &
                  dip_vel(2,n)**2 + dip_vel(3,n)**2
         end if
         diprms = diprms + (pol(n)*field(1,n)-ind_dip(1,n))**2 + &
               (pol(n)*field(2,n)-ind_dip(2,n))**2 + &
               (pol(n)*field(3,n)-ind_dip(3,n))**2
         dipself = dipself + 0.5d0*(ind_dip(1,n)*ind_dip(1,n)+ &
               ind_dip(2,n)*ind_dip(2,n)+ind_dip(3,n)*ind_dip(3,n)) / &
               pol(n)
      end if
   end do
   dipndf = ic
   call timer_stop(TIME_DIPUP)
   return
end subroutine get_dipinfo 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dipole_stat here]
subroutine dipole_stat(crd,charge,dipole,ibgwat,nres,ipres, &
      aveper,aveind,avetot,emtot)
   use constants, only : INV_AMBER_ELECTROSTATIC
   implicit none
   _REAL_ crd(3,*),charge(*),dipole(3,*), &
         aveper,aveind,avetot
   integer ibgwat,nres,ipres(*)

   _REAL_ pper,pind,ptot,px,py,pz,pindx,pindy,pindz, &
         pi_mag,pp_mag,emx,emy,emz,emtot,term1,emsq
   integer j,i,llim,iulim
   integer kstep
   save kstep
   data kstep/0/

   call timer_start(TIME_DIPUP)
   pper = 0.0d0
   pind = 0.0d0
   ptot = 0.0d0
   emx = 0.0d0
   emy = 0.0d0
   emz = 0.0d0
   emtot = 0.d0

   if (ibgwat > 0) then
      do j = ibgwat,nres
         px = 0.0d0
         py = 0.0d0
         pz = 0.0d0
         pindx = 0.0d0
         pindy = 0.0d0
         pindz = 0.0d0
         llim = ipres(j)
         iulim = ipres(j+1)-1
         do i = llim,iulim
            px = px + charge(i)*crd(1,i)
            py = py + charge(i)*crd(2,i)
            pz = pz + charge(i)*crd(3,i)
            pindx = pindx + dipole(1,i)
            pindy = pindy + dipole(2,i)
            pindz = pindz + dipole(3,i)
         end do
         emx = emx + px+pindx
         emy = emy + py+pindy
         emz = emz + pz+pindz
         pi_mag = sqrt(pindx**2 + pindy**2 + pindz**2)
         pp_mag = sqrt(px**2 + py**2 + pz**2)
         pper = pper + pp_mag
         pind = pind +  pi_mag
         ptot = ptot + sqrt( (px+pindx)**2 + (py+pindy)**2 &
               + (pz+pindz)**2)
      end do
   end if
   kstep = kstep + 1
   term1 = 4.8d0*INV_AMBER_ELECTROSTATIC
   emsq = ((emx*term1)**2 + (emy*term1)**2 + (emz*term1)**2)
   emtot = emtot + emsq/nres
   aveper = (pper / nres) * term1
   aveind = (pind / nres) * term1
   avetot = (ptot / nres) * term1
   call timer_stop(TIME_DIPUP)

   return
end subroutine dipole_stat 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_torque here]
subroutine get_torque(torque,dipole,field,nstart,ntop)
   implicit none
   integer nstart,ntop
   _REAL_ torque(3,*),dipole(3,*),field(3,*)
   integer n
   
   !     ---torque = dipole x field
   
   do n = nstart,ntop
      torque(1,n) = dipole(2,n)*field(3,n)-dipole(3,n)*field(2,n)
      torque(2,n) = dipole(3,n)*field(1,n)-dipole(1,n)*field(3,n)
      torque(3,n) = dipole(1,n)*field(2,n)-dipole(2,n)*field(1,n)
   end do
   return
end subroutine get_torque 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine add_induced here]
subroutine add_induced(numatoms,ind_dip,dipole)
   implicit none
   integer numatoms
   _REAL_ ind_dip(3,*),dipole(3,*)
   integer n,i
   call timer_start(TIME_DIPUP)
   do n = 1,numatoms
      dipole(1,n) = dipole(1,n) + ind_dip(1,n)
      dipole(2,n) = dipole(2,n) + ind_dip(2,n)
      dipole(3,n) = dipole(3,n) + ind_dip(3,n)
   end do
   call timer_stop(TIME_DIPUP)
   return
end subroutine add_induced 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cp_dips here]
subroutine cp_dips(numatoms,pol,x,dt)
   implicit none
   _REAL_ pol(*),x(*),dt
   integer istart,iend,numatoms
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "debug.h"
#ifdef MPI
#  include "parallel.h"
#endif
   if ( do_debugf /= 0 )return

#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = numatoms
#endif
   call timer_start(TIME_DIPUP)
   call do_cp_dip(istart,iend,pol,x(lfield),x(linddip), &
         x(ldipvel),diptemp,dipmass,diptau,nttdip,dt)
   call timer_stop(TIME_DIPUP)
end subroutine cp_dips 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_cp_dip here]
subroutine do_cp_dip(istart,iend,pol,field,ind_dip,dip_vel, &
      diptemp,dipmass,diptau,nttdip,dt)
   use constants, only : AMBER_ELECTROSTATIC2
   implicit none
   integer istart,iend,nttdip
   _REAL_ pol(*),field(3,*),ind_dip(3,*), &
         dip_vel(3,*),dipself
   _REAL_ diptemp,dipmass,diptau,dt
   _REAL_ small
   _REAL_ dtx,boltz2,dttp,fac,rndf, &
         massinv,fac2,ffac,wfac,scale_t
   integer n
   small = 1.d-6
   
   dttp = dt / diptau
   dtx = dt*20.455d0
   massinv = AMBER_ELECTROSTATIC2/dipmass
   if ( nttdip == 1 )then
      scale_t = sqrt(1.d0 - dttp)
   else
      scale_t = 1.d0
   end if
   
   !     ---update the velocities
   
   do n = istart,iend
      if ( pol(n) > small )then
         wfac = dtx*massinv
         ffac = 1.d0/pol(n)
         
         !         ---force is field - dip/pol
         
         dip_vel(1,n) = (dip_vel(1,n) + &
               wfac*(field(1,n)-ind_dip(1,n)*ffac))*scale_t
         dip_vel(2,n) = (dip_vel(2,n) + &
               wfac*(field(2,n)-ind_dip(2,n)*ffac))*scale_t
         dip_vel(3,n) = (dip_vel(3,n) + &
               wfac*(field(3,n)-ind_dip(3,n)*ffac))*scale_t
      end if
   end do
   
   !     ---now update the dipoles
   
   do n = istart,iend
      if ( pol(n) > small )then
         ind_dip(1,n) = ind_dip(1,n) + dtx*dip_vel(1,n)
         ind_dip(2,n) = ind_dip(2,n) + dtx*dip_vel(2,n)
         ind_dip(3,n) = ind_dip(3,n) + dtx*dip_vel(3,n)
      end if
   end do
   
   return
end subroutine do_cp_dip 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_dips here]
subroutine get_dips(x,nr)
   use file_io_dat
   implicit none
   _REAL_ x(*)
   integer nr,irest
   _REAL_ tt
   integer i,nr3,natom
#  include "ew_mpole.h"
   
   if ( indmeth /= 3 )return
   
   if ( irstdip /= 0 ) then
      write(6,9108)
      call amopen(19,inpdip,'O','F','R')
      nr3 = 3*nr
      read(19,9008) title1
      read(19,9018) natom,tt
      if( natom /= nr ) then
         write(6,9118)
         call mexit(6, 1)
      end if
      read(19,9028) (x(linddip+i-1),i=1,nr3)
      read(19,9028) (x(ldipvel+i-1),i=1,nr3)
      close(19)
      irstdip = 1
      write(6,9008) title1
      write(6,9009) tt
   end if
   
   !       ---open the restart file for dipoles
   
   call amopen(20,rstdip,owrite,'F','W')
   return
   
   9108 format(/,'   3.  induced dipoles and velocities',/)
   9118 format(/2x,'FATAL: NATOM mismatch in dipole and ', &
         'topology files')
   9008 format(a80)
   9009 format(t2,'begin time read from input dipoles =',f10.3,' ps'/)
   9018 format(i5,5e15.7)
   9028 format(6f12.7)
end subroutine get_dips 
!-------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wrt_dips here]
subroutine wrt_dips(dip,dpv,nr,tt,title)
   implicit none
   _REAL_ dip(*),dpv(*)
   _REAL_ tt
   integer nr
   character(len=80) title
   integer i,nr3
#  include "ew_mpole.h"
   nr3 = 3*nr
   write(20,9008) title
   write(20,9018) nr,tt
   write(20,9028) (dip(i),i=1,nr3)
   write(20,9028) (dpv(i),i=1,nr3)
   9008 format(a80)
   9018 format(i5,5e15.7)
   9028 format(6f12.7)
   rewind(20)
   return
end subroutine wrt_dips 
