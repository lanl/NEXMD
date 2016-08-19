#include "dprec.fh"
#include "assert.fh"

! Subroutine to calculate the probability to hop and the velocity adjustment after an effective hop
module fewest_switches
use naesmd_constants
use random
use langevin_temperature
use communism
use naesmd_space_module
implicit none
contains
   subroutine evalhop(sim, lprint,ido,neq,tini,tend,toldivprk, &
      param,Na,Nm,atm2,mdflag, &
      d,E0,Omega,fosc,yg,cross,idocontrol,ibo)
   implicit none
   type(simulation_t), pointer :: sim
   integer Na,Nm,lprint,ibo
   integer mdflag
   integer k,i,j,icheck,itest,ini,ihopavant
   integer atm2(Na)
   integer ido,neq,idocontrol
   include 'sizes'
   _REAL_ tini,tend,toldivprk,param(50)
   _REAL_ g(nmaxpot),gacum(nmaxpot)
   _REAL_ iseedhop,eavant, eapres 
   _REAL_ E0,d
   include 'md.par'
   include 'parH.par'
   _REAL_ Omega(Mx_M),fosc(Mx_M)
   _REAL_ xx(Na),yy(Na),zz(Na)
   include 'common'
   _REAL_ yg(nmaxpot),ytemp,ytemp2
   _REAL_ t_start,t_finish 
   integer cross(nmaxpot),crosstemp,ininonhop
   external fcn
   idocontrol=0
! conthop is used to not allow crossing inmediately after a hop
! conthop2 is used to not allow hoppings inmediately after a crossing
   if(conthop.eq.3) conthop=0
   if(conthop2.eq.3) conthop2=0
   if(conthop.gt.0) conthop=conthop+1
   if(conthop2.gt.0) conthop2=conthop2+1
   if(conthop.gt.0) then
      if(cross(ihop).eq.2) then
            if(iordenhop(ihop).ne.ihopprev) conthop=0
      end if
   end if
   iseedhop=rranf1(iseedmdqt)
   ihopavant=ihop
   eavant=vmdqtnew(ihop)
   do j=1,natom
      eavant=eavant+0.5d0*massmdqt(j)*(vx(j)**2+vy(j)**2+vz(j)**2)
   end do
   do j=1,npot
      g(j)=0.d0
      gacum(j)=0.d0
   end do
! g(j) is the probability to hop from the 
! current state to the state j
   do j=1,npot
      if(j.ne.ihop) then
         g(j)=vnqcorrhoptot(j,ihop)/(nqold**2)
         if(g(j).lt.0.0d0) g(j)=0.0d0
      end if
   end do
   icheck=0
   do j=1,npot
      gacum(j)=0.0d0
      if(j.ne.ihop) then
         do k=1,j
            if(k.ne.ihop) gacum(j)=gacum(j)+g(k)
         end do
      end if
   end do
   itest=0
   do j=1,npot
      if(j.ne.ihop) then
         if(iseedhop.le.gacum(j).and.itest.eq.0) then
            icheck=j
            itest=1
         end if
      end if
   end do
   crosstemp=cross(ihop)
   if(icheck.ne.0.and.cross(ihop).ne.2) then
      ini=0
      ! adjustment of velocities
      if(conthop2.gt.0) then
         if(ihopprev.ne.icheck) conthop2=0
      end if
      if(conthop2.eq.0) then
         call veladjustment(sim, lprint,Na,Nm,atm2,mdflag, &
            icheck,ini,d,E0,Omega,fosc)
      if(lprint.ge.2) then
         write(33,*) tfemto,icheck,ini
         call flush(33)
      end if
   else
      ini=1
   end if
   if(ini.eq.0) then
      ihop=icheck
      call naesmd2qmmm_r(sim)
      call cpu_time(t_start)
      call deriv(sim,ihop)
      call cpu_time(t_finish)
      sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
      call do_sqm_davidson_update(sim,vmdqt=vmdqtnew,vgs=vgs)  
      if(decorhop.eq.1) then
         do j=1,npot
            yg(j)=0.0d0
            yg(j+npot)=rranf1(iseedmdqt)
         end do
         yg(ihop)=1.d0
         ido=3
         call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
         ido=1
         idocontrol=1
      end if
      conthop=1
   end if
        endif
   if(crosstemp.eq.2.and.conthop.eq.0) then
      conthop2=1
! ihopprev keep the value of the previous state from where we hop
! in order to allow new hops for the next 2 steps only in case that
! we do not want to hop again to the same state
      ihopprev=ihop
      ini=0
      ihop=iordenhop(ihop)
      icheck=ihop
! after the hop, we reinitialize the variables
      call naesmd2qmmm_r(sim)
      call cpu_time(t_start)
      call deriv(sim,ihop)
      call cpu_time(t_finish)
      sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
      do j=1,natom
         xx(j)=rx(j)
         yy(j)=ry(j)
         zz(j)=rz(j)
      end do
      call do_sqm_davidson_update(sim,vmdqt=vmdqtnew,vgs=vgs)        
      ytemp=yg(ihopavant)
      ytemp2=yg(ihopavant+npot)
      yg(ihopavant)=yg(ihop)
      yg(ihopavant+npot)=yg(ihop+npot)
      yg(ihop)=ytemp
      yg(ihop+npot)=ytemp2
      if(idocontrol.eq.0) then
         ido=3
         call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
         ido=1
         idocontrol=1
      end if
   end if
! Evaluation of other crossings that do not involve the ihop state
!*************************************************************
   ininonhop=1
   do i=1,npot
   if(i.lt.iorden(i).and.i.ne.ihopavant.and.i.ne.iorden(ihopavant)) then
      if(cross(i).eq.2) then
         ininonhop=0
         icheck=i
! after the hop, we reinicialize the variables
         ytemp=yg(i)
         ytemp2=yg(i+npot)
         yg(i)=yg(iorden(i))
         yg(i+npot)=yg(iorden(i)+npot)
         yg(iorden(i))=ytemp
         yg(iorden(i)+npot)=ytemp2
         if(idocontrol.eq.0) then
            ido=3
            call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
            ido=1
            idocontrol=1
         end if
      end if
   end if
   enddo
! Check the conservation of the total energy in the hop
! and recalculate the kin to be printed in writeoutput.f
   if(icheck.ne.0) then
      if(ini.eq.0.or.ininonhop.eq.0) then
         eapres=vmdqtnew(ihop)
         kin=0.0d0
         do j=1,natom
            kin=kin+massmdqt(j)*(vx(j)**2+vy(j)**2+vz(j)**2)/2
         end do
         eapres=eapres+kin
         do j=1,npot
            if(j.ne.ihop) then
               if(j.lt.iorden(j)) then
                  if(cross(j).ne.0) then
                     write(30,887) tfemto,cross(j),j,iorden(j),eavant,eapres
                     call flush(30)
                  end if
               end if
            else
               if(ihop.ne.ihopavant) then
                  write(30,887) tfemto,cross(ihop),ihopavant,ihop,eavant,eapres
                  call flush(30)
               end if
            end if
         end do
      end if
   end if
!***********************************
! end analyze the hopping 
!**********************************
889   FORMAT(I3,10000(1X,F18.10))
888   FORMAT(10000(1X,F18.10))
887   FORMAT(F18.10,1X,I2,1X,I3,1X,I3,10000(1X,F18.10))
   return
   end subroutine
! At the point of hop, in general, the value of the potential energy in the new
! surface is different to the one in the older.
! In order to conserve the energy, we adjust the velocities
   subroutine veladjustment(sim, lprint,Na,Nm,atm2,mdflag,icheck,ini,d, &
      E0,Omega,fosc) 
   implicit none
   type(simulation_t),pointer::sim
   integer Na,Nm,lprint
   integer mdflag
   integer atm2(Na)
   integer i,j,icheck,ini,ihoptemp 
   include 'sizes'
   include 'md.par'
   include 'md.cmn'
   double precision dij(nmax*3),vicheck
   double precision alpha,racine,ctehop1,dctehop1 
   double precision vtemp(nmaxpot),vgstemp
   real*8 xx(Na),yy(Na),zz(Na)
   real*8 E0,d
   include 'parH.par'
   real*8 Omega(Mx_M),fosc(Mx_M)
   include 'common'
!********************************************************
! adjustment of velocities
!********************************************************
! Added by ST: calculate here NAC <psi| d psi/dR> in one step:
!   Current energy calculation
   call do_sqm_davidson_update(sim,vmdqt=vtemp,vgs=vgstemp)        
!    Feed here energies and wavefunctions and geometry, get back dij  
! if necessary here the signs of the CI coefficient matrix can be checked right here
! analytical calculation of nacR
   call nacR_analytic_wrap(sim, ihop, icheck, dij)
! end of the calculation of the non-adiabatic coupling vector dij
!*********************************
   if(lprint.ge.1) then
      j=1
      do i=1,natom
         write(29,*) i,dij(j),dij(j+1),dij(j+2)
         j=j+3
      end do
      call flush(29)
   end if
! calculation of the current energy
! and the velocities adjustment
   ihoptemp=icheck
   vicheck=vtemp(ihoptemp)
   alpha=vmdqtnew(ihop)-vicheck
   racine = 0.0d0

   j=1
   do i=1,natom
      racine=racine+vx(i)*dij(j)+vy(i)*dij(j+1) &
         +vz(i)*dij(j+2)
      j=j+3
   end do
   racine=racine**2
   j=1
   do i=1,natom
      racine=racine+2.0d0*alpha/massmdqt(i) &
         *(dij(j)**2+dij(j+1)**2+dij(j+2)**2)
      j=j+3
   end do
   if(racine.le.0.0d0) then
      ini=1
      goto 4321
   end if
   ctehop1=0.0d0
   j=1
   do i=1,natom
      ctehop1=ctehop1+ &
         +vx(i)*dij(j)+vy(i)*dij(j+1)+vz(i)*dij(j+2)
      j=j+3
   end do
   ctehop1=ctehop1+dsqrt(racine)
   dctehop1=0.d0
   j=1
   do i=1,natom
      dctehop1=dctehop1+1.0d0/massmdqt(i) &
        *(dij(j)**2+dij(j+1)**2+dij(j+2)**2)
      j=j+3
   end do
   ctehop1=ctehop1/dctehop1

! option to adjust the velocities in the direction of
! the nonadiabatic coupling vector
   j=1
   do i=1,natom
      vx(i)=vx(i)-ctehop1*dij(j)/massmdqt(i)
      vy(i)=vy(i)-ctehop1*dij(j+1)/massmdqt(i)
      vz(i)=vz(i)-ctehop1*dij(j+2)/massmdqt(i)
      j=j+3
   end do
4321     continue
!********************************************************
! end of adjustment of velocities
!********************************************************
889   format(10000(1x,f18.10))
   return
   end subroutine
end module
