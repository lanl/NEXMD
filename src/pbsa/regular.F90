#define _REAL_ double precision

! ----- for a regular point (i0,j0,k0), find corresponding coefficients


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine regula here]
subroutine regula(l,m,n,hx,hy,hz,b_in,b_out,i0,j0,k0,info,x,y,z, coe1,rhs)
   implicit none

  !Passed variables
   integer   l,m,n,i0,j0,k0,info
   _REAL_    hx,hy,hz,b_in,b_out,rhs
   _REAL_    x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_    coe1(7)
  
  !Local variables
   
   _REAL_    bx1,bx2,by1,by2,bz1,bz2,fkk
   _REAL_    fb_in,fb_out,ff_in,ff_out,fk_in,fk_out

   if (info <= 3) then
      bx1 = fb_in(b_in,b_out,0.5*(x(i0-1)+x(i0)),y(j0),z(k0))
      bx2 = fb_in(b_in,b_out,0.5*(x(i0+1)+x(i0)),y(j0),z(k0))
      by1 = fb_in(b_in,b_out,x(i0),0.5*(y(j0-1)+y(j0)),z(k0))
      by2 = fb_in(b_in,b_out,x(i0),0.5*(y(j0+1)+y(j0)),z(k0))
      bz1 = fb_in(b_in,b_out,x(i0),y(j0),0.5*(z(k0-1)+z(k0)))
      bz2 = fb_in(b_in,b_out,x(i0),y(j0),0.5*(z(k0+1)+z(k0)))
      fkk = fk_in(x(i0),y(j0),z(k0))
      rhs = ff_in(x(i0),y(j0),z(k0))
   else
      bx1 = fb_out(b_in,b_out,0.5*(x(i0-1)+x(i0)),y(j0),z(k0))
      bx2 = fb_out(b_in,b_out,0.5*(x(i0+1)+x(i0)),y(j0),z(k0))
      by1 = fb_out(b_in,b_out,x(i0),0.5*(y(j0-1)+y(j0)),z(k0))
      by2 = fb_out(b_in,b_out,x(i0),0.5*(y(j0+1)+y(j0)),z(k0))
      bz1 = fb_out(b_in,b_out,x(i0),y(j0),0.5*(z(k0-1)+z(k0)))
      bz2 = fb_out(b_in,b_out,x(i0),y(j0),0.5*(z(k0+1)+z(k0)))
      fkk = fk_out(x(i0),y(j0),z(k0))
      rhs = ff_out(x(i0),y(j0),z(k0))
   end if

   coe1(2) = bx1/hx/hx
   coe1(3) = bx2/hx/hx
   coe1(4) = by1/hy/hy
   coe1(5) = by2/hy/hy
   coe1(6) = bz1/hz/hz
   coe1(7) = bz2/hz/hz
   coe1(1) = - (coe1(2)+coe1(3)+coe1(4)+coe1(5)+coe1(6) &
         +coe1(7))+fkk

   return
end subroutine regula 



