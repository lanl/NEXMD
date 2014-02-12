#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine curv here]
subroutine curv(l,m,n,h,hx,hy,hz,x1,y1,z1,i0,j0,k0,x,y,z,phi,t, xyy,xzz,xyz)
   implicit none

   !       ----------------------------------------------
   !       Compute curvature at the projection (x1,y1,z1)
   !       ----------------------------------------------
   !Passed variables:
   integer l,m,n,i0,j0,k0
   _REAL_  h,hx,hy,hz,x1,y1,z1,x(0:l+1), y(0:m+1), z(0:n+1),&
           phi(0:l+1,0:m+1,0:n+1),t(3,3),xyy,xzz,xyz
   
   !Local variables:
   _REAL_ ph0,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz,phieps,&
          tmp1,tmp2,tmp3,tmp

   call grtopr(l,m,n,x,y,z,h,hx,hy,hz,x1,y1,z1,i0,j0,k0,phi,1,1,1, &
         ph0,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz)

   phieps = t(1,1)*phx+t(1,2)*phy+t(1,3)*phz

   if ( abs(phieps) < 1.0e-20) then
      write(*,*) "   Error: phieps = 0.0 in curv()!"
      stop
   end if

   tmp1 = t(2,1)*phxx+t(2,2)*phxy+t(2,3)*phxz
   tmp2 = t(2,1)*phxy+t(2,2)*phyy+t(2,3)*phyz
   tmp3 = t(2,1)*phxz+t(2,2)*phyz+t(2,3)*phzz

   tmp = t(2,1)*tmp1+t(2,2)*tmp2+t(2,3)*tmp3
   xyy = -tmp/phieps

   tmp1 = t(3,1)*phxx+t(3,2)*phxy+t(3,3)*phxz
   tmp2 = t(3,1)*phxy+t(3,2)*phyy+t(3,3)*phyz
   tmp3 = t(3,1)*phxz+t(3,2)*phyz+t(3,3)*phzz

   tmp = t(3,1)*tmp1+t(3,2)*tmp2+t(3,3)*tmp3
   xzz = -tmp/phieps

   tmp = t(2,1)*tmp1+t(2,2)*tmp2+t(2,3)*tmp3
   xyz = -tmp/phieps

   return
end subroutine curv 

