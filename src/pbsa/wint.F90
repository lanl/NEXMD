#define _REAL_ double precision


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wint here]
subroutine wint(l,m,n,h,hx,hy,hz,xs,ys,zs,maxirr,nq,index,index2,cirreg, wp, &
      wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe, &
      w0p,wyp,wzp,wyyp,wzzp,wyzp)
   implicit none

   !Passed variables
   integer  l,m,n,maxirr,nq
   _REAL_   xs,ys,zs,hx,hy,hz,h

   integer  index(l,m,n), index2(l,m,n)
   _REAL_   cirreg(maxirr, 15)
   _REAL_   wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_   wyycoe(maxirr, nq)
   _REAL_   wzzcoe(maxirr, nq), wyzcoe(maxirr, nq)
   _REAL_   wp(maxirr)
   _REAL_   w0p(maxirr),wyp(maxirr),wzp(maxirr)
   _REAL_   wyyp(maxirr),wzzp(maxirr),wyzp(maxirr)

  !Local variables
  integer   i,j,k,ir
  _REAL_    w0,wy,wz,wyy,wzz,wyz

   do i=1,l
      do j=1,m
         do k=1,n
            if (index(i,j,k) > 1 .and. index(i,j,k) < 5) then
               ir = index2(i, j, k)
               call wint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg, wp, &
                     wcoe,wycoe,wzcoe, wyycoe,wzzcoe,wyzcoe, &
                     w0,wy,wz,wyy,wzz,wyz)
               w0p(ir)  = w0
               wyp(ir)  = wy
               wzp(ir)  = wz
               wyyp(ir) = wyy
               wzzp(ir) = wzz
               wyzp(ir) = wyz
            end if
         end do
      end do
   end do

   return
end subroutine wint 




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wint2 here]
subroutine wint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg, wp, &
      wcoe,wycoe,wzcoe, wyycoe,wzzcoe,wyzcoe, &
      w0,wy,wz,wyy,wzz,wyz)

   implicit none

   !Passed variables

   integer,parameter :: nqq = 27 
   integer   l,m,n,ir,maxirr,nq
   _REAL_    xs,ys,zs,h,hx,hy,hz

   integer   index2(l,m,n)
   _REAL_    cirreg(maxirr, 15)
   _REAL_    wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_    wyycoe(maxirr, nq)
   _REAL_    wzzcoe(maxirr, nq), wyzcoe(maxirr, nq)
   _REAL_    wp(maxirr)
   _REAL_    w0,wy,wz,wyy,wzz,wyz

  !Local variables

   integer   if0,i0,j0,k0,nsub,i1,j1,k1,i2,j2,k2,ip,ii,idis,nn2
   _REAL_    hmax,alf,x1,y1,z1,x2,y2,z2,x3,y3,z3,dis,xyy,xzz,xyz
   _REAL_    t(3,3)
   integer   itemp(nqq)

   hmax=h
   if (nq /= nqq) then
      write(*,*) "   Please change nqq in wint2() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 10.1*hmax
   if0 = int(alf/hmax) + 1

   w0  = 0.0
   wy  = 0.0
   wz  = 0.0
   wyy = 0.0
   wzz = 0.0
   wyz = 0.0

   x1 = cirreg(ir, 1)
   y1 = cirreg(ir, 2)
   z1 = cirreg(ir, 3)

   xyy = cirreg(ir, 4)
   xzz = cirreg(ir, 5)
   xyz = cirreg(ir, 6)

   t(1,1) = cirreg(ir, 7)
   t(1,2) = cirreg(ir, 8)
   t(1,3) = cirreg(ir, 9)
   t(2,1) = cirreg(ir, 10)
   t(2,2) = cirreg(ir, 11)
   t(2,3) = cirreg(ir, 12)
   t(3,1) = cirreg(ir, 13)
   t(3,2) = cirreg(ir, 14)
   t(3,3) = cirreg(ir, 15)

   i0 = nint((x1-xs)/hx)
   j0 = nint((y1-ys)/hy)
   k0 = nint((z1-zs)/hz)


   nsub = 0
   do i2=0, if0          !# Starting from the center
      do i1 = i0-i2, i0+i2
         do j1 = j0-i2, j0+i2
            do k1 = k0-i2, k0+i2
               if (i1 < 1 .or. j1 < 1 .or. k1 < 1 &
                     .or. i1 > l .or. j1 > m .or. k1 > n) goto 555

               idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
               nn2 = index2(i1,j1,k1)

               if (idis == i2 .and. nn2 > 0) then

                  x2 = cirreg(nn2, 1)
                  y2 = cirreg(nn2, 2)
                  z2 = cirreg(nn2, 3)

                  do ii = 1, nsub
                     ip = itemp(ii)
                     x3 = cirreg(ip, 1)
                     y3 = cirreg(ip, 2)
                     z3 = cirreg(ip, 3)
                     call distan(x2,y2,z2,x3,y3,z3, dis)
                     if (dis < 0.1*hmax) then
                        goto 555
                     end if
                  end do

                  call distan(x1,y1,z1,x2,y2,z2, dis)
                  if ( dis < alf ) then
                     nsub = nsub + 1
                     itemp(nsub) = nn2
                     ip = index2(i1,j1,k1)
                     w0  = w0 + wcoe(ir, nsub)*wp(ip)
                     wy = wy + wycoe(ir, nsub)*wp(ip)
                     wz = wz + wzcoe(ir, nsub)*wp(ip)
                     wyy = wyy + wyycoe(ir, nsub)*wp(ip)
                     wzz = wzz + wzzcoe(ir, nsub)*wp(ip)
                     wyz = wyz + wyzcoe(ir, nsub)*wp(ip)
                     if ( nsub .eq. nq) goto 577
                  end if
               end if
               555 continue
            end do  ! k1 = k0-i2, k0+i2
         end do  ! j1 = j0-i2, j0+i2
      end do  ! i1 = i0-i2, i0+i2
   end do  ! i2=0, if0
   577 continue

   return
end subroutine wint2 




