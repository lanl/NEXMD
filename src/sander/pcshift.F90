! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
!============================================================================

!                          PROGRAM PCSHIFT
!   Pseudocontact shift constraints for energy minimization and MD

!                             written by
!              Giovanni Gori-Savellini and Andrea Romagnoli
!             Department of Chemistry, University of Florence
!                   e-mail: nanni@risc1.lrm.fi.cnr.it
!                           roma@risc1.lrm.fi.cnr.it

!============================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pcshift here]
subroutine pcshift(natom,x,f)
   
   use file_io_dat
   implicit none
   integer:: i, ihelp, ihep, iin, itmp, m, n, natom, nk, nres
   _REAL_ :: f, prob, rf, rms, vett, x, xrfac6, xrfact, z
#  include "nmr.h"
#  include "md.h"
#  include "pcshift.h"
#  include "def_time.h"
   integer first
   dimension vett(np)
   dimension x(3,*),f(3,*)
   save first
   data first /1/
   character(len=3) nmpmc
   namelist /pcshf/ nprot,nfe,obs,wt,iprot,tolpro,mltpro, &
         optphi,opttet,optomg, &
         opta1,opta2,optkon,nmpmc
   
   !=====================================================================
   !   First time through, set up a number of arrays, and read input data:
   !=====================================================================
   
   if (first > 0) then
      nstampa=0
      nstampaten=0
      ! If restraint input has been redirected, open the appropriate file
      iin = 7
      if (iredir(7) /= 0) then
         call amopen(37, redir(7)(1:iredir(7)),'O','F','R')
         iin = 37
         write(6,10) redir(7)(1:iredir(7))
         10 format(' Pseudocontact shifts will be read from file: ',a)
      end if
      
      ! read the namelist pcshf
      
      read(iin,nml=pcshf,err=30)
      write(6,*) 'namelist no reports error '
      goto 40
      30 write(6,*) 'namelist reports error reading &pcshf'
      call mexit(6,1)
      40 continue
      if(nprot > mshfd.or. nprot < 1) then
         write(6,70) nprot, mshf
         70 format('nprot out of range:',2i5)
         call mexit(6,1)
      end if
      itmp=0
      write(6,*) 'Number of Total protons ', nprot
      write(6,*) 'Number of paramagnetic center ', nfe
      do i=1,natom
         if (resat(i)(1:3) == nmpmc)then
            read(resat(i)(10:13),'(i4)')nres
            itmp=itmp+1
            ippmc(itmp)=i
            
            !  assign the coordinate of paramagnetic center(s)
            
            cmx(itmp)=x(1,i)
            cmy(itmp)=x(2,i)
            cmz(itmp)=x(3,i)
            write(6,*) 'Coord. of paramagnetic center(s) : ',nmpmc
            write(6,*) cmx(itmp),cmy(itmp),cmz(itmp)
         end if
      end do
      
      ! minimizz. iniziale del tensore
      
      ! assign the coordinate of protons
      
      do i=1,nprot
         ihelp = iprot(i)
         coox(i) = x(1,ihelp)
         cooy(i) = x(2,ihelp)
         cooz(i) = x(3,ihelp)
      end do
      
      ! initialize variables
      
      do n=1,3
         do m=1,natom
            d(n,m) = 0.0d0
         end do
      end do
      
      ! inizialize tensor(s) parameters
      
      do nk=1,nfe
         vett((nk-1)*mpar+1)=optphi(nk)
         vett((nk-1)*mpar+2)=opttet(nk)
         vett((nk-1)*mpar+3)=optomg(nk)
         vett((nk-1)*mpar+4)=opta1(nk)
         vett((nk-1)*mpar+5)=opta2(nk)
      end do
      
      ! evaluate initial shifts and derivates
      
      call calcresid(vett)
      call dispshift
      epcshf = 0.0d0
      call derivshift(vett)

      first = 0
   end if
   
   !  End of "first" loop
   
   call timer_start(TIME_SHFDER)
   
   
   !=====================================================================
   !=====================================================================
   !     aggiornamento del tensore
   
   do nk=1,nfe
      vett((nk-1)*mpar+1)=optphi(nk)
      vett((nk-1)*mpar+2)=opttet(nk)
      vett((nk-1)*mpar+3)=optomg(nk)
      vett((nk-1)*mpar+4)=opta1(nk)
      vett((nk-1)*mpar+5)=opta2(nk)
   end do
   
   ! read the new coordinates of protons & paramagnetic center(s)
   
   do i=1,nfe
      cmx(i) = x(1,ippmc(i))
      cmy(i) = x(2,ippmc(i))
      cmz(i) = x(3,ippmc(i))
   end do
   do i=1,nprot
      ihep = iprot(i)
      coox(i) = x(1,ihep)
      cooy(i) = x(2,ihep)
      cooz(i) = x(3,ihep)
   end do
   
   ! final result
   
   if (natom < 0) then
      rewind(53)
      write(53,*) ' Final pseudocontact shifts'
      call dispshift
      return
   end if
   
   !  estimate pseudocontact shifts & derivates
   
   epcshf = 0.0d0
   call derivshift(vett)
   !                       print intermediate results
   nstampa=nstampa+1
   if (nstampa > 9) then
      call ptrir
      nstampa=0
   end if
   nstampaten=nstampaten+1
   if (nstampaten > 999) then
      rewind(53)
      call dispshift
      nstampaten=0
   end if

   
   ! --  update the force array
   
   do n=1,3
      do m=1,natom
         f(n,m) = f(n,m) + d(n,m)
         d(n,m) = 0.d0
      end do
   end do

   call timer_stop_start(TIME_SHFDER,TIME_KMAT)
   if (iprint /= 0) then
      write(6,44) epcshf
      44 format(40x,'Total pcshift constraint:',f8.2)
      write(6,390)
      390 format(/21x,'#  Pearson r           rms error'/18x,35('-'))
      call pearsn(obsp,shavp,ipear,rf,prob,z,rms,xrfact,xrfac6,iuse)
      write(6,'(a18,i5,3f10.5)') 'pcshift correlation:',iuse,rf,z,rms
   end if

   call timer_stop(TIME_KMAT)
   
   return
end subroutine pcshift 
!-------------------------------------------------------------
!                     CALCRESID
! evaluate pseudocontact shift and residuals

! IN:  vett = tensor(s) parameters

! OUT: resid = total residue
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine calcresid here]
subroutine calcresid(vett)
   use constants, only : zero,RAD_TO_DEG 
   implicit none
   integer:: i, j, k, l, m, nprotml
   _REAL_ :: a1d, a2d, app, axx, axy, axz, ayx, ayy, &
        ayz, azx, azy, azz, calcolo, &
        costmult, g1, g2, o, p, r, scalx, scaly, &
        scalz, t, tmp2, vett, xapp, yapp, zapp

#  include "pcshift.h"
   dimension vett(np)
   
   !                 calcolo dell errore e delle violazioni
   
   shift(1:nprot)=zero
   resid=zero
   iviolation=0
   do k=1,nfe
      p=vett((k-1)*mpar+1)
      t=vett((k-1)*mpar+2)
      o=vett((k-1)*mpar+3)
      a1d=vett((k-1)*mpar+4)
      a2d=vett((k-1)*mpar+5)
      axx=cos(p)*cos(o)
      axy=sin(p)*cos(o)
      axz=sin(o)
      ayx=-cos(t)*sin(p)-sin(o)*cos(p)*sin(t)
      ayy=cos(t)*cos(p)-sin(o)*sin(p)*sin(t)
      ayz=sin(t)*cos(o)
      azx=sin(t)*sin(p)-sin(o)*cos(p)*cos(t)
      azy=-sin(t)*cos(p)-sin(o)*sin(p)*cos(t)
      azz=cos(t)*cos(o)
      i=1
      m=1
      do while (i <= nprot)
         calcolo=zero
         nprotml=mltpro(i)
         if (nprotml > 1) then
            costmult=1/dble(nprotml)
            app=zero
            do l=1,nprotml
               xapp=coox(i)-cmx(k)
               yapp=cooy(i)-cmy(k)
               zapp=cooz(i)-cmz(k)
               r=sqrt(xapp**2+yapp**2+zapp**2)
               scalz=(xapp*azx+yapp*azy+zapp*azz)
               scalx=(xapp*axx+yapp*axy+zapp*axz)
               scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
               g1=(sqrt(3.d0)*scalz-r)*(sqrt(3.d0)*scalz+r)
               g2=(scalx-scaly)*(scalx+scaly)
               app=app+(a1d*g1+1.5d0*a2d*g2)*costmult/r**5
               i=i+1
            end do
            do l=1,nprotml
               shift(i-l)=shift(i-l)+app
            end do
         else
            xapp=coox(i)-cmx(k)
            yapp=cooy(i)-cmy(k)
            zapp=cooz(i)-cmz(k)
            r=sqrt(xapp**2+yapp**2+zapp**2)
            scalz=(xapp*azx+yapp*azy+zapp*azz)
            scalx=(xapp*axx+yapp*axy+zapp*axz)
            scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
            g1=(sqrt(3.d0)*scalz-r)*(sqrt(3.d0)*scalz+r)
            g2=(scalx-scaly)*(scalx+scaly)
            shift(i)=shift(i)+(a1d*g1+1.5d0*a2d*g2)/r**5
            i=i+1
         end if
      end do
   end do

   do j=1,nprot
      tmp2=abs(shift(j)-obs(j))-tolpro(j)
      if (tmp2 > zero) then
         iviolation=iviolation+1
         resid=resid+tmp2**2*wt(j)/dble(mltpro(j))
      end if
   end do
   if (nstampa == 0) then
      write(6,'(1x,A)') 'TENSOR PARAMETERS'
      write (6,31) (optphi(i)*RAD_TO_DEG,opttet(i)*RAD_TO_DEG &
            ,optomg(i)*RAD_TO_DEG,i=1,nfe)
      31 format(1x,'PHI= ',f7.3,1x,'TETA= ',f7.3,1x,'OMEGA= ',f7.3)
      write (6,32) (opta1(i),opta2(i),i=1,nfe)
      32 format (1x,'A1= ',f9.3,1x,'A2= ',f9.3)
      write (6,30) iviolation,resid
      30 format (1x,' VIOLATIONS ',i6,1x,' RESIDUALS ',f10.5)
   else
      write (6,40) iviolation,resid
      40 format (15x,' VIOLATIONS ',i6,1x,' RESIDUALS ',f10.5)
   end if
   return
end subroutine calcresid 
!--------------------------------------------------------------
!                             PTRIR
! print intermediate results
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ptrir here]
subroutine ptrir
   use constants, only: zero, RAD_TO_DEG 
   implicit none
   integer:: i, j
   _REAL_ :: tmp2

#  include "nmr.h"
#  include "pcshift.h"
   resid=zero
   iviolation=0
   do j=1,nprot
      tmp2=abs(shift(j)-obs(j))-tolpro(j)
      if (tmp2 > zero) then
         iviolation=iviolation+1
         resid=resid+tmp2**2*wt(j)/dble(mltpro(j))
      end if
   end do

   if (nstampa == 0) then
      write(6,'(1x,A)') 'TENSOR PARAMETERS'
      write (6,31) (optphi(i)*RAD_TO_DEG,opttet(i)*RAD_TO_DEG &
            ,optomg(i)*RAD_TO_DEG,i=1,nfe)
      31 format(1x,'PHI= ',f7.3,1x,'TETA= ',f7.3,1x,'OMEGA= ',f7.3)
      write (6,32) (opta1(i),opta2(i),i=1,nfe)
      32 format (1x,'A1= ',f9.3,1x,'A2= ',f9.3)
      write (6,30) iviolation,resid
      30 format (1x,' VIOLATIONS ',i6,1x,' RESIDUALS ',f10.5)
   else
      write (6,40) iviolation,resid
      40 format (15x,' VIOLATIONS ',i6,1x,' RESIDUALS ',f10.5)
   end if
   return
end subroutine ptrir 
!--------------------------------------------------------------
!                       DERIVSHIFT
! evaluate the derivates

! IN: vett = tensor(s) parameters

! OUT: d   = derivates
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine derivshift here]
subroutine derivshift(vett)
   
   implicit none
   integer:: i, k, l, m, nprotml, num
   _REAL_ :: a1d, a2d, app, axx, axy, axz, ayx, ayy, &
        ayz, azx, azy, azz, bx, bxy, bxz, by, byz, bz, cader, cadx, cady, &
        cadz, calcolo, capdx, capdy, capdz, &
        costmult, cx, cxy, cxz, cy, cyz, cz, den, eshfi, fac, &
        g1, g2, o, p, r, scalx, scaly, scalz, &
        sgnerr, shav, shi, stdev, t, tmp, tmp2, &
        vett, xapp, yapp, zapp

#  include "nmr.h"
#  include "pcshift.h"
   dimension vett(np)
   dimension cadx(mshfd),cady(mshfd),cadz(mshfd), &
         capdx(mshfd),capdy(mshfd),capdz(mshfd)

   do i=1,nprot
      shift(i)=0.0d0
   end do
   ! evaluate the pseudocontact shifts
   do k=1,nfe
      p=vett((k-1)*mpar+1)
      t=vett((k-1)*mpar+2)
      o=vett((k-1)*mpar+3)
      a1d=vett((k-1)*mpar+4)
      a2d=vett((k-1)*mpar+5)
      axx=cos(p)*cos(o)
      axy=sin(p)*cos(o)
      axz=sin(o)
      ayx=-cos(t)*sin(p)-sin(o)*cos(p)*sin(t)
      ayy=cos(t)*cos(p)-sin(o)*sin(p)*sin(t)
      ayz=sin(t)*cos(o)
      azx=sin(t)*sin(p)-sin(o)*cos(p)*cos(t)
      azy=-sin(t)*cos(p)-sin(o)*sin(p)*cos(t)
      azz=cos(t)*cos(o)

      tmp=0.0d0
      tmp2=0.0d0
      i=1
      m=1
      stdev=0.0d0
      shi=0.0d0
      ipear=0
      do while (i <= nprot)
         calcolo=0.0d0
         cader=0.0d0
         nprotml=mltpro(i)
         costmult=1/dble(nprotml)
         if (nprotml > 1) then
            costmult=1/dble(nprotml)
            app=0.0d0
            do l=1,nprotml
               xapp=coox(i)-cmx(k)
               yapp=cooy(i)-cmy(k)
               zapp=cooz(i)-cmz(k)
               r=sqrt(xapp**2+yapp**2+zapp**2)
               scalz=(xapp*azx+yapp*azy+zapp*azz)
               scalx=(xapp*axx+yapp*axy+zapp*axz)
               scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
               g1=(sqrt(3.d0)*scalz-r)*(sqrt(3.d0)*scalz+r)
               g2=(scalx-scaly)*(scalx+scaly)
               app=app+(a1d*g1+1.5d0*a2d*g2)*costmult/r**5
               i=i+1
            end do
            do l=1,nprotml
               shift(i-l)=shift(i-l)+app
            end do
         else
            xapp=coox(i)-cmx(k)
            yapp=cooy(i)-cmy(k)
            zapp=cooz(i)-cmz(k)
            r=sqrt(xapp**2+yapp**2+zapp**2)
            scalz=(xapp*azx+yapp*azy+zapp*azz)
            scalx=(xapp*axx+yapp*axy+zapp*axz)
            scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
            g1=(sqrt(3.d0)*scalz-r)*(sqrt(3.d0)*scalz+r)
            g2=(scalx-scaly)*(scalx+scaly)
            shift(i)=shift(i)+(a1d*g1+1.5d0*a2d*g2)/r**5
            i=i+1
         end if
      end do
   end do
   
   ! evaluate the derivates
   
   do i=1,nprot
      cadx(i)=0.0d0
      cady(i)=0.0d0
      cadz(i)=0.0d0
   end do
   i=1
   do k=1,nfe
      do while (i <= nprot)
         shi=shift(i)
         nprotml=mltpro(i)
         costmult=1.0d0/nprotml
         xapp=coox(i)-cmx(k)
         yapp=cooy(i)-cmy(k)
         zapp=cooz(i)-cmz(k)
         r=sqrt(xapp**2+yapp**2+zapp**2)
         bx=a1d*(3*azx**2-1)+1.5d0*a2d*(axx**2-ayx**2)
         by=a1d*(3*azy**2-1)+1.5d0*a2d*(axy**2-ayy**2)
         bz=a1d*(3*azz**2-1)+1.5d0*a2d*(axz**2-ayz**2)
         bxy=6*a1d*azx*azy+2*(axx*axy-ayx*ayy)
         bxz=6*a1d*azx*azz+2*(axx*axz-ayx*ayz)
         byz=6*a1d*azy*azz+2*(axy*axz-ayy*ayz)

         num=xapp**2*bx+yapp**2*by+zapp**2*bz+ &
               xapp*yapp*bxy+xapp*zapp*bxz+yapp*zapp*byz

         cx=axx**2+ayx**2+azx**2
         cy=axy**2+ayy**2+azy**2
         cz=axz**2+ayz**2+azz**2
         cxy=2*(axx*axy+ayx*ayy+azx*azy)
         cxz=2*(axx*axz+ayx*ayz+azx*azz)
         cyz=2*(axy*axz+ayy*ayz+azy*azz)

         den=r**5

         capdx(i)=num*(-2.5d0)*(2*xapp)/r**7+ &
               (2*xapp*bx+yapp*bxy+zapp*bxz)/den

         capdy(i)=num*(-2.5d0)*(2*yapp)/r**7+ &
               (2*yapp*by+xapp*bxy+zapp*byz)/den

         capdz(i)=num*(-2.5d0)*(2*zapp)/r**7+ &
               (2*zapp*bz+xapp*bxz+yapp*byz)/den

         if ((shi-obs(i)) >= 0.0) then
            sgnerr=-(abs(shi-obs(i)))
         else
            sgnerr=abs(shi-obs(i))
         end if
         cadx(i)=cadx(i)+2*capdx(i)*sgnerr*costmult
         cady(i)=cady(i)+2*capdy(i)*sgnerr*costmult
         cadz(i)=cadz(i)+2*capdz(i)*sgnerr*costmult
         i=i+1
      end do
   end do
   
   ! evaluate energy contributions
   
   i=1
   do while (i <= nprot)
      shi=shift(i)
      nprotml=mltpro(i)
      costmult=1.0d0/nprotml
      if (wt(i) > 0.0) then
         shav = abs(shi-obs(i))-tolpro(i)
         if (ipnlty == 1 .or. ipnlty == 3) then
            if (shav > 0.0) then
               eshfi = wt(i)*shav*costmult
               fac= wt(i)
            else
               eshfi = 0.0d0
               fac = 0.0d0
            end if
         else if (ipnlty == 2) then
            if (shav > 0.0) then
               eshfi = wt(i)*shav**2*costmult
               fac= wt(i)
            else
               eshfi = 0.0d0
               fac = 0.0d0
            end if
         else
            write(6,*) 'bad values for ipnlty:', ipnlty
            call mexit(6,1)
         end if
         ipear = ipear + 1
         shavp(ipear) = shi
         obsp(ipear) = obs(i)
         epcshf=epcshf+eshfi*optkon
         
         ! evaluate the pseudocontact force term
         
         d(1,iprot(i))=optkon*fac*cadx(i)
         d(2,iprot(i))=optkon*fac*cady(i)
         d(3,iprot(i))=optkon*fac*cadz(i)
      end if  ! (wt(i) > 0.0)
      i=i+1
   end do
   return
end subroutine derivshift 
!-------------------------------------------------------------
!                      DISPSHIFT

! OUT: write the violated pcshifts in the scratch file fort.53

!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dispshift here]
subroutine dispshift
   
   implicit none
   integer:: i
   _REAL_ :: err, tmp2
#  include "nmr.h"
#  include "pcshift.h"
   write(53,*) '           Visualizzazione degli shifts'
   write(53,*) '    Atom                Oss.       Calc.        Err.'
   do i=1,nprot
      tmp2=abs(shift(i)-obs(i))-tolpro(i)
      if( tmp2 > 0.0 ) then
         err = tmp2**2*wt(i)/dble(mltpro(i))
         if (err > pencut) then
            write(53,12) iprot(i),resat(iprot(i)), &
                  obs(i),shift(i),err
            12 format(1x,i4,1x,a13,1x,f9.5,1x,f9.5,1x,f9.5)
         end if
      end if
   end do
   return
end subroutine dispshift 
