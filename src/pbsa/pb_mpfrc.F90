! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of multipole expansion reaction field interactions
subroutine pb_mpfrc(natom,atmfirst,atmlast,lmax,rdiel,xctr,yctr,zctr,epsin,epsout,cg,x,f,eel)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
    
#  include "pb_constants.h"

   ! passed variables
    
   integer natom,atmfirst,atmlast,lmax
   _REAL_ rdiel, xctr, yctr, zctr
   _REAL_ epsin, epsout, cg(natom), x(3,natom)
   _REAL_ f(3,natom), eel
    
   ! local variables
    
   integer k,l,m,i,j,lmpol,lstl(lmax*(lmax+3)/2+1),lstm(lmax*(lmax+3)/2+1)
   _REAL_ acofftmp,r2,factortmp,ddx,ddy,ddz,rst,ar,ebulk,factconv
   _REAL_ aqlmoutrtmp,aqlmoutitmp,aqlmoutr(lmax*(lmax+3)/2+1),aqlmouti(lmax*(lmax+3)/2+1)
   _REAL_ cmpk(natom,lmax),smpk(natom,lmax)
   _REAL_ acoff(0:lmax),factor(0:lmax,0:lmax)
   _REAL_ plmk(natom,0:lmax+1,0:lmax),dplmk(natom,0:lmax+1,0:lmax),acr(natom,0:lmax)
   _REAL_ dact(natom),dacp(natom),dasp(natom),da1ost(natom),dast(natom)
   _REAL_ act(natom),ast(natom),acp(natom),asp(natom)
   _REAL_ plmktmp(natom),dplmktmp(natom),fx(natom),fy(natom),fz(natom),cmpktmp(natom),smpktmp(natom)
   _REAL_ dx,dy,dz,dxr,dxt,dxp,dyr,dyt,dyp,dzr,dzt,drq,dtq,dpq,dl,dm
    
   eel = ZERO; fx = ZERO; fy = ZERO; fz = ZERO
    
   ! calculate r, theta, phi of every particle
    
   do i = atmfirst, atmlast
     ddx = x(1,i)-xctr; ddy = x(2,i)-yctr; ddz = x(3,i)-zctr
     r2 = ddx**2+ddy**2+ddz**2; ar = sqrt(r2)
     acr(i,0) = cg(i)
     do l = 1, lmax
        acr(i,l) = acr(i,l-1)*ar
     end do
     if ( ddx == ZERO .and. ddy == ZERO .and. ddz == ZERO ) then
        dact(i) = ONE
        dast(i) = ZERO
        dacp    = ONE
        dasp    = ZERO
     else if ( ddx == ZERO .and. ddy == ZERO .and. ddz /= ZERO) then
        dact(i) = ddz/ar
        dast(i) = sqrt(ONE-ddz**2/r2)
        dacp(i) = ONE
        dasp(i) = ZERO
     else
       dact(i)  = ddz/ar
       dast(i)  = sqrt(ONE-ddz**2/r2)
       rst      = sqrt(r2-ddz**2)
       dacp(i)  = ddx/rst
       dasp(i)  = ddy/rst
     end if
     if ( dast(i) == ZERO ) then
       da1ost(i) = ZERO
     else
       da1ost(i) = ONE/dast(i)
     end if
     act(i) = dact(i)
     ast(i) = dast(i)
     acp(i) = dacp(i)
     asp(i) = dasp(i)
   end do

   ! calculate factorial

   do l = 0, lmax
      factor(l,0) = ONE
      factor(l,0) = sqrt( dble(2*l+1)/( FOURPI*factor(l,0) ) )
      do m = 1, l
         factor(l,m) = ONE
         do j = (l-m+1), (l+m)
            factor(l,m) = factor(l,m)*dble(j)
         end do
         factor(l,m) = sqrt( dble(2*l+1)/( FOURPI*factor(l,m) ) )
      end do
   end do
     
   ! calculate legendre polynomial
    
   call plgndr(natom,atmfirst,atmlast,lmax,dact,plmk) 
    
   ! calculate derivative of legendre polynomial
    
   call dplgndr(natom,atmfirst,atmlast,lmax,dact,plmk,dplmk) 
    
   ! calculate cosmphi and sinmphi
    
   call cosmphi(natom,atmfirst,atmlast,lmax,dacp,cmpk)
   call sinmphi(natom,atmfirst,atmlast,lmax,dacp,dasp,smpk)
    
   ! calculate qlm for input charges
    
   lmpol=lmax*(lmax+3)/2+1
   call genlist(lmax,lstl,lstm)

   aqlmoutr = ZERO; aqlmouti = ZERO
   do i = 1, lmax+1
      l = lstl(i)
      do k = atmfirst, atmlast
         aqlmoutr(i) = aqlmoutr(i) + acr(k,l)*factor(l,0)*plmk(k,l,0)
      end do
   end do
   do i = lmax+2, lmpol
      l = lstl(i)
      m = lstm(i)
      do k = atmfirst, atmlast
         aqlmoutr(i) = aqlmoutr(i) + acr(k,l)*factor(l,m)*plmk(k,l,m)*cmpk(k,m)
         aqlmouti(i) = aqlmouti(i) - acr(k,l)*factor(l,m)*plmk(k,l,m)*smpk(k,m)
      end do  
   end do

   ! calculate the coffient

   ebulk = epsout/epsin
   factconv = (18.2223d0)**2

   do l = 0, lmax
      acoff(l) = FOURPI*(ONE-ebulk)/(dble(2*l+1)*(rdiel**(2*l+1))*(ebulk+dble(l)/dble(l+1)))
   end do

   ! calculate the reaction field energy for input charges

   eel = ZERO
   do i = 1, lmax+1
      l = lstl(i)
      eel = eel + acoff(l)*aqlmoutr(i)*aqlmoutr(i)
   end do
   do i = lmax+2, lmpol
      l = lstl(i)
      eel = eel + acoff(l)*TWO*( aqlmoutr(i)*aqlmoutr(i)+aqlmouti(i)*aqlmouti(i) )
   end do
   eel = HALF*factconv*eel
    
   ! calculate reaction field force
   ! first computer all contributes of term with m=0 
   ! when l=0, dr,dt,dp=0, so l is from 1 to lmax, i from 2 to lmpol
    
   do i = 2, lmax+1
      l = lstl(i)
      dl = dble(l)
      factortmp = factor(l,0)
      acofftmp = acoff(l)
      aqlmoutrtmp = aqlmoutr(i)
      do k = atmfirst, atmlast
         plmktmp(k) = plmk(k,l,0)*factortmp
         dplmktmp(k) = dplmk(k,l,0)*factortmp
         cmpktmp(k) = acofftmp*acr(k,l-1)*aqlmoutrtmp
      end do
       
      do k = atmfirst, atmlast
         drq = dl*plmktmp(k)*cmpktmp(k)
         dtq = -ast(k)*dplmktmp(k)*cmpktmp(k)
            
         dxr = (ast(k)*acp(k))*drq
         dxt = (act(k)*acp(k))*dtq
         dx = dxr+dxt
            
         dyr = (ast(k)*asp(k))*drq
         dyt = (act(k)*asp(k))*dtq
         dy = dyr+dyt
            
         dzr = (act(k))*drq
         dzt = (ast(k))*dtq
         dz = dzr-dzt
            
         fx(k) = fx(k) - dx
         fy(k) = fy(k) - dy
         fz(k) = fz(k) - dz
      end do
    
   end do  ! i = 2, lmax+1
    
   ! terms with m /= 0
    
   do i = lmax+2, lmpol
      l = lstl(i)
      m = lstm(i)
      dl = dble(l)
      dm = dble(m)
      factortmp = factor(l,m)
      acofftmp = acoff(l)
      aqlmoutrtmp = aqlmoutr(i)
      aqlmoutitmp = aqlmouti(i)

      do k = atmfirst, atmlast
         plmktmp(k) = plmk(k,l,m)*factortmp
         dplmktmp(k) = dplmk(k,l,m)*factortmp
         cmpktmp(k) = ( cmpk(k,m)*aqlmoutrtmp-smpk(k,m)*aqlmoutitmp )*acofftmp*acr(k,l-1)
         smpktmp(k) = ( smpk(k,m)*aqlmoutrtmp+cmpk(k,m)*aqlmoutitmp )*da1ost(k)*acofftmp*acr(k,l-1)
      end do
       
      do k = atmfirst, atmlast
         drq = 2.0*dl*cmpktmp(k)*plmktmp(k)
         dpq = -2.0*dm*smpktmp(k)*plmktmp(k)
         dtq = -2.0*cmpktmp(k)*dplmktmp(k)*ast(k)
          
         dxr = (ast(k)*acp(k))*drq
         dxt = (act(k)*acp(k))*dtq
         dxp = (asp(k))*dpq
         dx = dxr+dxt-dxp
          
         dyr = (ast(k)*asp(k))*drq
         dyt = (act(k)*asp(k))*dtq
         dyp = (acp(k))*dpq
         dy = dyr+dyt+dyp
          
         dzr = (act(k))*drq
         dzt = (ast(k))*dtq
         dz = dzr-dzt
          
         fx(k) = fx(k) - dx
         fy(k) = fy(k) - dy
         fz(k) = fz(k) - dz
      end do
       
   end do  ! i = lmax+2, lmpol
    
   do k = atmfirst, atmlast
      f(1,k) = f(1,k) + fx(k)*factconv
      f(2,k) = f(2,k) + fy(k)*factconv
      f(3,k) = f(3,k) + fz(k)*factconv
   end do
    
contains    
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute p(l,m)
subroutine plgndr(natm,atmfirst,atmlast,lmax,x,ap)
    
   implicit none

   ! passed variables
    
   integer natm,lmax,atmfirst,atmlast
   _REAL_ x(natm),ap(natm,0:lmax+1,0:lmax)
    
   ! local variables
    
   integer l,m,k
   _REAL_ fact,somx2
    
   ! compute p(m,m)
    
   fact = ONE
   do m = 1, lmax
      do k = atmfirst, atmlast
         ap(k,0,0) = ONE
         somx2 = sqrt((ONE-x(k))*(ONE+x(k)))
         ap(k,m,m) = ap(k,m-1,m-1)*fact*somx2
      end do
      fact = fact+TWO
   end do
    
   ! compute p(m+1,m)
    
   do m = 0, lmax
      do k = atmfirst, atmlast
         ap(k,m+1,m) = x(k)*dble(2*m+1)*ap(k,m,m)
      end do
   end do
    
   ! compute p(l,m)
    
   do l = 2, lmax
      do m = 0, l-2
         do k = atmfirst, atmlast
            ap(k,l,m) = (x(k)*dble(2*l-1)*ap(k,l-1,m)-dble(l+m-1)*ap(k,l-2,m))/dble(l-m)
         end do
      end do
   end do
    
end subroutine plgndr
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute derivative of p(l,m)
subroutine dplgndr(natm,atmfirst,atmlast,lmax,x,ap,adp)
    
   implicit none

   ! passed variables
    
   integer natm,lmax,atmfirst,atmlast
   _REAL_ x(natm),ap(natm,0:lmax+1,0:lmax),adp(natm,0:lmax+1,0:lmax)
    
   ! local variables
    
   integer l,m
   _REAL_ fact(natm)
    
   do k = atmfirst, atmlast
      fact(k) = ONE/sqrt(ONE-x(k)*x(k))
   end do
    
   do l = 0, lmax
      do m =0, l-1
         do k = atmfirst, atmlast
            adp(k,l,m) = fact(k)*(ap(k,l,m+1)-fact(k)*m*x(k)*ap(k,l,m))
         end do 
      end do
      do k = atmfirst, atmlast
         adp(k,l,l) = -fact(k)*fact(k)*l*x(k)*ap(k,l,l)
      end do
   end do
    
end subroutine dplgndr
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generate list of qlm, l, m
subroutine genlist(lmax,lstl,lstm)
    
   implicit none

   ! passed variables
    
   integer lmax,lstl(*),lstm(*)
    
   ! local variables
    
   integer l,m,index1

   ! always the same order in spherical harmonics

   index1 = 0
   do l = 0, lmax
      index1 = index1+1
      lstl(index1) = l
      lstm(index1) = 0
   enddo
   do l = 1, lmax
      do m = 1, l
         index1 = index1+1
         lstl(index1) = l
         lstm(index1) = m
      end do
   end do
    
end subroutine genlist
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ cos(m*phi) calculation (m > 0)
subroutine cosmphi(natm,atmfirst,atmlast,lmax,acp,cmpk)

   implicit none

   ! passed variables

   integer natm,lmax,atmfirst,atmlast
   _REAL_ acp(natm),cmpk(natm,lmax)

   ! local variables

   integer k,m

   do k = atmfirst, atmlast
      cmpk(k,1) = acp(k)
      cmpk(k,2) = TWO*acp(k)*acp(k)-ONE
      cmpk(k,3) = (TWO*cmpk(k,2)-ONE)*acp(k)
   end do
   do m = 4,lmax
      do k = atmfirst, atmlast
         cmpk(k,m) = TWO*acp(k)*cmpk(k,m-1)-cmpk(k,m-2)
      end do
   end do
    
end subroutine cosmphi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ sin(m*phi) calculation (m > 0)
subroutine sinmphi(natm,atmfirst,atmlast,lmax,acp,asp,smpk)
    
   implicit none

   ! passed variables
    
   integer natm,lmax,atmfirst,atmlast
   _REAL_ acp(natm),asp(natm),smpk(natm,lmax)
    
   ! local variables
    
   integer k,m
    
   do k = atmfirst, atmlast
      smpk(k,1) = asp(k)
      smpk(k,2) = TWO*acp(k)*asp(k)
      smpk(k,3) = TWO*acp(k)*smpk(k,2)-asp(k)
   end do
   do m = 4, lmax
      do k = atmfirst, atmlast
         smpk(k,m) = TWO*acp(k)*smpk(k,m-1)-smpk(k,m-2)
      end do
   end do
    
end subroutine sinmphi

end subroutine pb_mpfrc
