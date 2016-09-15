#include "dprec.fh"
#include "assert.fh"
!
!********************************************************************      
!********************************************************************      
!   
   subroutine fcn(n,x,yg,yprime)
   implicit none
   integer n,k,j
   _REAL_ x,yg(n),yprime(n)
 
   call interpolate(n,x)
   call vqcalc(n,yg,yprime)

   return
   end subroutine
!
!********************************************************************
!
!  Subroutine to interpolate the values of the terms <phi_k | d phi_j/ d t >
!  and the adiabatic energies
!
!********************************************************************
!
   subroutine interpolate(n,x)
   use naesmd_module
   use md_module
   implicit none
   integer n,k,j
   _REAL_ x
 
   do k=1,n/2
      do j=1,n/2
         cadiab(k,j)=cadiabmiddleold(k,j) &
            +bcoeffcadiab(k,j)*(x-tini0) 
      end do
   end do

   do k=1,n/2
      vmdqt(k)=vmdqtmiddleold(k)+bcoeffvmdqt(k)*(x-tini0) 
   end do
        
889   FORMAT(10000(1X,F18.10))

   return
   end subroutine
!
!********************************************************************
!
!  Subroutine to calculate de d(yg)/dt and dqq/dt
!  That is, the derivatives of the modulus and phase of the
!  electronic coefficients
!
!********************************************************************
!

   subroutine vqcalc(n,yg,yprime)
   use naesmd_module
   use md_module
   implicit none
   integer n,k,j
   _REAL_ x,yg(n),yprime(n)
   _REAL_ norm 

   do k=1,n/2
      yprime(k)=0.d0
      yprime(k+n/2)=-1.0d0*vmdqt(k)

      do j=1,n/2
!                vnqcorrhop(k,j) = -2.0d0*dsqrt(yg(j)*yg(k))*
!     $dcos(yg(j+n)-yg(k+n))*cadiab(k,j)
         vnqcorrhop(k,j)=-1.0d0*yg(j)* &
            dcos(yg(j+n/2)-yg(k+n/2))*cadiab(k,j)

         yprime(k)=yprime(k)+vnqcorrhop(k,j) 
         vnqcorrhop(k,j)=vnqcorrhop(k,j)*2.0d0*yg(k)
      end do

      if(dabs(yg(k)).lt.1.0d-7) then
         yprime(k+n/2)=0.0d0
      else
         do j=1,n/2
!                vqqcorr(k) = vqqcorr(k) - dsqrt(yg(j)/yg(k))*
!     $dsin(yg(j+n)-yg(k+n))*cadiab(k,j) 
            yprime(k+n/2)=yprime(k+n/2)-yg(j)/yg(k)* &
               dsin(yg(j+n/2)-yg(k+n/2))*cadiab(k,j) 
         end do
      end if
   end do

!        norm=0.0d0
!        do k = 1,n
!            norm=norm+yg(k)*yprime(k)
!        enddo
!        write(10,*) norm
!        call flush(10)

889   FORMAT(10000(1X,F18.10))

   return
   end subroutine
!
