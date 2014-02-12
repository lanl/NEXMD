#define _REAL_ double precision


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine transf here]
subroutine transf(l,m,n,h,hx,hy,hz,xi,yi,zi,i0,j0,k0,x,y,z,phi, t)
   implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  transf defines a local coordinate transfomation matrix    *                                 *
   !       *                                                            *
   !       **************************************************************
   
   !Passed variables:
   integer l,m,n,i0,j0,k0
   _REAL_ h,hx,hy,hz,xi,yi,zi,x(0:l+1),y(0:m+1), z(0:n+1),&
          phi(0:l+1,0:m+1,0:n+1),t(3,3)
   
   !Local variables:
  
   _REAL_ ph0,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz,t11,t12,&
          t13,tp1,t21,t22,t23,tp2,t31,t32,t33,tp3

   call grtopr(l,m,n,x,y,z,h,hx,hy,hz,xi,yi,zi,i0,j0,k0,phi,1,1,0, &
         ph0,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz)
   !       print *,xi,yi,zi
   !       print *,'transf',ph0
   !       print *,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz

   if (phx == 0.0 .and. phz == 0) then

      t11 = phx
      t12 = phy
      t13 = phz
      tp1 = sqrt(t11*t11 + t12*t12 + t13*t13)

      if (abs(tp1) < 1.0e-20) then
         write(*,*) "   Error: tp1 = 0.0 in transf()!"
         stop
      end if

      t(1,1) = t11/tp1
      t(1,2) = t12/tp1
      t(1,3) = t13/tp1

      t21 = t12
      t22 = -t11
      t23 = 0.0
      tp2 = sqrt(t21*t21 + t22*t22 + t23*t23)

      if (abs(tp2) < 1.0e-20) then
         write(*,*) "   Error: tp2 = 0.0 in transf()!"
         stop
      end if

      t(2,1) = t21/tp2
      t(2,2) = t22/tp2
      t(2,3) = t23/tp2

      t31 = t11*t13
      t32 = t12*t13
      t33 = -t11*t11-t12*t12
      tp3 = sqrt(t31*t31 + t32*t32 + t33*t33)

      if (abs(tp3) < 1.0e-20) then
         write(*,*) "   Error: tp3 = 0.0 in transf()!"
         stop
      end if

      t(3,1) = t31/tp3
      t(3,2) = t32/tp3
      t(3,3) = t33/tp3
   else
      t11 = phx
      t12 = phy
      t13 = phz
      tp1 = sqrt(t11*t11 + t12*t12 + t13*t13)

      if (abs(tp1) < 1.0e-20) then
         write(*,*) "   Error: tp1 = 0.0 in transf()!"
         stop
      end if

      t(1,1) = t11/tp1
      t(1,2) = t12/tp1
      t(1,3) = t13/tp1

      t21 = t13
      t22 = 0.0
      t23 = -t11
      tp2 = sqrt(t21*t21 + t22*t22 + t23*t23)

      if (abs(tp2) < 1.0e-20) then
         write(*,*) "   Error: tp2 = 0.0 in transf()!"
         stop
      end if

      t(2,1) = t21/tp2
      t(2,2) = t22/tp2
      t(2,3) = t23/tp2

      t31 = -t11*t12
      t32 = t11*t11 + t13*t13
      t33 = -t12*t13
      tp3 = sqrt(t31*t31 + t32*t32 + t33*t33)

      if (abs(tp3) < 1.0e-20) then
         write(*,*) "   Error: tp3 = 0.0 in transf()!"
         stop
      end if

      t(3,1) = t31/tp3
      t(3,2) = t32/tp3
      t(3,3) = t33/tp3
   end if  ! (phx == 0.0 .and. phz == 0)

   return
end subroutine transf 

! ----- End of transf()

