! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.fh"

! geometrical transformation of ss-overlapp and hamilton matrix elements

subroutine skss(x,x2,i,j,r2,iovpar,em,ne)
   implicit none

!! Passed in
   _REAL_ , intent(in ) :: x(6)
   _REAL_ , intent(in ) :: x2(6)
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   _REAL_ , intent(out) :: em(1,1)
   integer, intent(in ) :: ne


!! Locals
   integer :: id
   _REAL_ :: parm(13)

   id=iovpar(i,j,r2,parm)
   em(1,1)=parm(10)

   return
end subroutine skss


! geometrical transformation of sp-overlapp and hamilton matrix elements
subroutine sksp(x,x2,i,j,r2,iovpar,em,emt,ne)
   implicit none 

!! Passed in
   _REAL_ , intent(in ) :: x(6)
   _REAL_ , intent(in ) :: x2(6)
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,3)
   _REAL_ , intent(out) :: emt(ne,3)

!! Locals
   integer :: id,l
   _REAL_  :: parm(13)

   id=iovpar(i,j,r2,parm)
   do l=1,3
      em(1,l)  = x(l) * parm(9)
      emt(l,1) = -em(1,l)
   end do

   return
end subroutine sksp


! geometrical transformation of sd-overlapp and hamilton matrix elements

subroutine sksd(x,x2,i,j,r2,iovpar,em,emt,ne)
   implicit none

!! Passed in
   _REAL_ , intent(in ) :: x(6)
   _REAL_ , intent(in ) :: x2(6)
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,5)
   _REAL_ , intent(out) :: emt(ne,5)

!! Locals
   integer :: id, l
   _REAL_  :: es(5), r3, d4, d5
   _REAL_  :: parm(13)

!--

   r3=sqrt(3.0)
   d4=x2(3)-0.5*(x2(1)+x2(2))
   d5=x2(1)-x2(2)
   id=iovpar(i,j,r2,parm)

   do l=1,3
      es(l)=r3*x(l)*x(l+1)
   end do
   es(4)=0.5*r3*d5
   es(5)=d4
   do l=1,5
      em(1,l)=es(l)*parm(8)
      emt(l,1)=em(1,l)
   end do
   return
end subroutine sksd


! geometrical transformation of pp-overlapp and hamilton matrix elements
subroutine skpp(x,x2,i,j,r2,iovpar,em,ne)
   implicit none

!! Passed in:
   _REAL_ , intent(in ) :: x(6)           
   _REAL_ , intent(in ) :: x2(6)          
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,3)        

!! Locals
   integer :: id, l, ir, is, ii, k
   _REAL_ :: dm(6),parm(13),epp(6), hp

!--

   id=iovpar(i,j,r2,parm)
   do l=1,3
      epp(l)=x2(l)
      epp(l+3)=x(l)*x(l+1)
30 enddo
   do l=1,3
      hp=epp(l)
      dm(l)=hp*parm(6)+(1.0-hp)*parm(7)
   enddo
   do l=4,6
      dm(l)=epp(l)*(parm(6)-parm(7))
   enddo
   do ir=1,3
      do is=1,ir
         ii=ir-is
         k=3*ii-(ii*(ii-1))/2+is
         em(is,ir)=dm(k)
         em(ir,is)=dm(k)
      enddo
   enddo
   return
end subroutine skpp


! geometrical transformation of pd-overlapp and hamilton matrix elements

subroutine skpd(x,x2,i,j,r2,iovpar,em,emt,ne)
   implicit none

!! Passed in
   _REAL_ , intent(in ) :: x(6)
   _REAL_ , intent(in ) :: x2(6)
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,5)
   _REAL_ , intent(out) :: emt(ne,5)


!! Locals
   integer :: id, l, m, ir, is, k
   _REAL_  :: r3, d3, d4, d5, d6
   _REAL_  :: parm(13),epd(13,2),dm(15)


   r3=sqrt(3.0)
   d3=x2(1)+x2(2)
   d4=x2(3)-0.5*d3
   d5=x2(1)-x2(2)
   d6=x(1)*x(2)*x(3)
   id=iovpar(i,j,r2,parm)

   do l=1,3
      epd(l,1)=r3*x2(l)*x(l+1)
      epd(l,2)=x(l+1)*(1.0-2.0*x2(l))
      epd(l+4,1)=r3*x2(l)*x(l+2)
      epd(l+4,2)=x(l+2)*(1.0-2.0*x2(l))
      epd(l+7,1)=0.5*r3*x(l)*d5
      epd(l+10,1)=x(l)*d4
   end do

   epd(4,1)=r3*d6
   epd(4,2)=-2.0*d6
   epd(8,2)=x(1)*(1.0-d5)
   epd(9,2)=-x(2)*(1.0+d5)
   epd(10,2)=-x(3)*d5
   epd(11,2)=-r3*x(1)*x2(3)
   epd(12,2)=-r3*x(2)*x2(3)
   epd(13,2)=r3*x(3)*d3
   do l=1,15
      dm(l)=0.0
   end do
   do m=1,2
      dm(1)=dm(1)+epd(1,m)*parm(m+3)
      dm(2)=dm(2)+epd(6,m)*parm(m+3)
      dm(3)=dm(3)+epd(4,m)*parm(m+3)
      dm(5)=dm(5)+epd(2,m)*parm(m+3)
      dm(6)=dm(6)+epd(7,m)*parm(m+3)
      dm(7)=dm(7)+epd(5,m)*parm(m+3)
      dm(9)=dm(9)+epd(3,m)*parm(m+3)

      do l=8,13
         dm(l+2)=dm(l+2)+epd(l,m)*parm(m+3)
      end do

   end do

   dm(4)=dm(3)
   dm(8)=dm(3)

   do ir=1,5
      do is=1,3
         k=3*(ir-1)+is
         emt(ir,is)=-dm(k)
         em(is,ir)=dm(k)
      enddo
   enddo
   return
end subroutine skpd


! geometrical transformation of dd-overlapp and hamilton matrix elements
subroutine skdd(x,x2,i,j,r2,iovpar,em,ne)
   implicit none

!! Passed in
   _REAL_ , intent(in ) :: x(6)
   _REAL_ , intent(in ) :: x2(6)
   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, external    :: iovpar
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,5)


!! Locals
   integer :: id, l, m, ir, is, ii, k
   _REAL_  :: r3, d3, d4, d5
   _REAL_  :: parm(13),e(15,3),dm(15), dd(3)

!--

   r3=sqrt(3.0)
   d3=x2(1)+x2(2)
   d4=x2(3)-0.5*d3
   d5=x2(1)-x2(2)
   id=iovpar(i,j,r2,parm)

   do l=1,3
      e(l,1)=x2(l)*x2(l+1)
      e(l,2)=x2(l)+x2(l+1)-4.0*e(l,1)
      e(l,3)=x2(l+2)+e(l,1)
      e(l,1)=3.0*e(l,1)
   end do

   e(4,1)=d5*d5
   e(4,2)=d3-e(4,1)
   e(4,3)=x2(3)+0.25*e(4,1)
   e(4,1)=0.75*e(4,1)
   e(5,1)=d4*d4
   e(5,2)=3.0*x2(3)*d3
   e(5,3)=d3*d3*0.75
   dd(1)=x(1)*x(3)
   dd(2)=x(2)*x(1)
   dd(3)=x(3)*x(2)

   do l=1,2
      e(l+5,1)=3.0*x2(l+1)*dd(l)
      e(l+5,2)=dd(l)*(1.0-4.0*x2(l+1))
      e(l+5,3)=dd(l)*(x2(l+1)-1.0)
   end do

   e(8,1)=dd(1)*d5*1.5
   e(8,2)=dd(1)*(1.0-2.0*d5)
   e(8,3)=dd(1)*(0.5*d5-1.0)
   e(9,1)=d5*0.5*d4*r3
   e(9,2)=-d5*x2(3)*r3
   e(9,3)=d5*0.25*(1.0+x2(3))*r3
   e(10,1)=x2(1)*dd(3)*3.0
   e(10,2)=(0.25-x2(1))*dd(3)*4.0
   e(10,3)=dd(3)*(x2(1)-1.0)
   e(11,1)=1.5*dd(3)*d5
   e(11,2)=-dd(3)*(1.0+2.0*d5)
   e(11,3)=dd(3)*(1.0+0.5*d5)
   e(13,3)=0.5*d5*dd(2)
   e(13,2)=-2.0*dd(2)*d5
   e(13,1)=e(13,3)*3.0
   e(12,1)=d4*dd(1)*r3
   e(14,1)=d4*dd(3)*r3
   e(15,1)=d4*dd(2)*r3
   e(15,2)=-2.0*r3*dd(2)*x2(3)
   e(15,3)=0.5*r3*(1.0+x2(3))*dd(2)
   e(14,2)=r3*dd(3)*(d3-x2(3))
   e(14,3)=-r3*0.5*dd(3)*d3
   e(12,2)=r3*dd(1)*(d3-x2(3))
   e(12,3)=-r3*0.5*dd(1)*d3

   do l=1,15
      dm(l)=0.0
      do m=1,3
         dm(l)=dm(l)+e(l,m)*parm(m)
      enddo
   enddo

   do ir=1,5
      do is=1,ir
         ii=ir-is
         k=5*ii-(ii*(ii-1))/2+is
         em(ir,is)=dm(k)
         em(is,ir)=dm(k)
      enddo
   enddo

   return
end subroutine skdd
