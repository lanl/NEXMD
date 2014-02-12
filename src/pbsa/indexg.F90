#define _REAL_ double precision

! ----- decide the total number of irregular points
!       for each irregular point, find its corresponding projection
!       point, and local transfomation matrix, etc.


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexg here]
subroutine indexg(l,m,n,h,hx,hy,hz,maxirr,x,y,z,index,phi,index2,cirreg)
  implicit none

   !common /lmn/l, m, n, nirreg
   integer l,m,n,nirreg,i,j,k,i0,j0,k0,maxirr
  
   _REAL_ x1,y1,z1,xyy,xzz,xyz
   _REAL_ x(0:l+1), y(0:m+1), z(0:n+1)
   integer  index(l,m,n),index2(l,m,n)
   _REAL_ phi(0:l+1,0:m+1,0:n+1)
   _REAL_ cirreg(maxirr,15)
   _REAL_ hx,hy,hz,h
   _REAL_ t(3,3)

   !                 print *,'WJ maxirr',maxirr
   nirreg = 0      !# number of irregular points

   !                 print *,l,m,n
   !                 print *,'WJ nirreg',nirreg
   !       print *,index(1,1,1:18)
  !write(8101,*) l,m,n,maxirr,x,y,z
  !write(8102,*) index
  !write(8103,*) phi
   h = x(1)-x(0)
   index2 = 0
   do k=1,n
      do j=1,m
         do i=1,l
            if (index(i,j,k) < 5 .and. index(i,j,k) > 1) then
               nirreg = nirreg + 1

               if (nirreg > maxirr) then
                  write(*,*) "   maxirr is defined too small!"
                  stop
               end if

               ! --------------- find the projection of the irregular point
               !                 print *,'nirreg=',nirreg
               !                print *,i,j,k,index(i,j,k)
               !                print *,'WJ project ...'
               call project(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi,x1,y1,z1)
               !                print *,x1,y1,z1
               !                 write (81,*),x1,y1,z1

               ! --------------- find the local coordinate transformation matrix

!                                print *,'WJ transf ...',i,j,k
               i0 = nint((x1 - x(0))/h)
               j0 = nint((y1 - y(0))/h)
               k0 = nint((z1 - z(0))/h)
               call transf(l,m,n,h,hx,hy,hz,x1,y1,z1,i0,j0,k0,x,y,z,phi, t)

               ! --------------- find the curvertures at the projection

!                                print *,'WJ curv ...'
               call curv(l,m,n,h,hx,hy,hz,x1,y1,z1,i0,j0,k0,x,y,z,phi,t, xyy,xzz,xyz)
               !WJstop
               !                 stop

               cirreg(nirreg, 1)  = x1
               cirreg(nirreg, 2)  = y1
               cirreg(nirreg, 3)  = z1
               cirreg(nirreg, 4)  = xyy
               cirreg(nirreg, 5)  = xzz
               cirreg(nirreg, 6)  = xyz
               cirreg(nirreg, 7)  = t(1,1)
               cirreg(nirreg, 8)  = t(1,2)
               cirreg(nirreg, 9)  = t(1,3)
               cirreg(nirreg, 10) = t(2,1)
               cirreg(nirreg, 11) = t(2,2)
               cirreg(nirreg, 12) = t(2,3)
               cirreg(nirreg, 13) = t(3,1)
               cirreg(nirreg, 14) = t(3,2)
               cirreg(nirreg, 15) = t(3,3)

               index2(i,j,k) = nirreg
            end if  ! (index(i,j,k) < 5 .and. index(i,j,k) > 1)
         end do  !  k=1,n
      end do  !  j=1,m
   end do  !  i=1,l

   ! write(*,*) "   The number of irregular points is ", nirreg
   ! write(*,*)

   return
end subroutine indexg 

! ----- End of indexg()
