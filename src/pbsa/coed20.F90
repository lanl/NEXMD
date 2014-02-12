#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine coed20 here]
subroutine coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,maxirr,nq,index2,cirreg, &
      wcoe, wxcoe, wycoe,wzcoe,wxxcoe,wyycoe, &
      wzzcoe,wxycoe,wxzcoe,wyzcoe)
   use iim_use
   implicit none
   integer,parameter :: nqq = 27 

   !Passed variables

   integer l,m,n,nirreg,maxirr,nq
   integer index2(l,m,n)
   _REAL_  xs,ys,zs,h,hx,hy,hz
   _REAL_  cirreg(maxirr, 15)
   _REAL_  wcoe(maxirr,nq),wxcoe(maxirr, nq),wycoe(maxirr, nq)
   _REAL_  wzcoe(maxirr,nq),wxxcoe(maxirr, nq),wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq),wxycoe(maxirr, nq)
   _REAL_  wxzcoe(maxirr, nq),wyzcoe(maxirr, nq)

   !Local variables
   integer i0,j0,k0,ir,isvd,i,j,k,i1,j1,k1,ip,job,ii,&
           ij,jj,nsub,nn2,i2,if0,idis,inf
   _REAL_  alf,x1,y1,z1,x2,y2,z2,x3,y3,z3,xyy,xzz,xyz,dis,hmax
   _REAL_  w1(20, nqq), w4(20)
   _REAL_  sd(21), uw(20,20), v(nqq, nqq)
   _REAL_  ew(nqq)
   _REAL_  ew1(nqq), ew2(nqq), ew3(nqq)
   _REAL_  ew4(nqq), ew5(nqq), ew6(nqq)
   _REAL_  ew7(nqq), ew8(nqq), ew9(nqq), ew10(nqq)
   _REAL_  w31(nqq), w32(nqq), w33(nqq)
   _REAL_  w34(nqq), w35(nqq), w36(nqq)
   _REAL_  w37(nqq), w38(nqq), w39(nqq), w310(nqq)
   _REAL_  sd2(20), work(1000)
   _REAL_  t(3,3)
   _REAL_  tempx(3),tempy(3)
   integer itemp(nqq)


   ! ----- set SVD subroutine: isvd=1  call ssvdc()
   !                           isvd=2  call sgecvd()
   !do k = 1,n
   !  do j = 1,m
   !    do i = 1,l
   !  write(17,*) i,j,k,index2(i,j,k)
   !    end do
   !    end do
   !    end do

   
   hmax=h
   isvd = 1

   if (nq /= nqq) then
      write(*,*) "   Please change nqq in coedef() s.t. nqq=nq in main()!"
      stop
   end if
   !alf = 10.1*hmax
   !alf = 14.1*hmax !for 0.25 space
    alf = 18.1*hmax
   if0 = int(alf/hmax) + 1

   do i=1, nirreg
      do j=1, nq
         wcoe(i,j)   = 0.0
         wxcoe(i,j)  = 0.0
         wycoe(i,j)  = 0.0
         wzcoe(i,j)  = 0.0
         wxxcoe(i,j) = 0.0
         wyycoe(i,j) = 0.0
         wzzcoe(i,j) = 0.0
         wxycoe(i,j) = 0.0
         wxzcoe(i,j) = 0.0
         wyzcoe(i,j) = 0.0
      end do
   end do

   ! ----- for each projection point

   do ir = 1, nirreg

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

      ! -------- initialize

      do jj=1, nq
         do ii=1, nq
            v(ii, jj) = 0.0      !# the V matrix in SVD: U*S*V
         end do
         do ii=1, 20
            w1(ii, jj) = 0.0     !# the coefficient matrix to be SVDed
         end do
      end do

      ! -------- Generate the matrix
      

      nsub = 0
      do i2=0, if0          !# Starting from the center
         do i1 = i0-i2, i0+i2
            do j1 = j0-i2, j0+i2
               do k1 = k0-i2, k0+i2
                  if (i1 < 1 .or. j1 < 1 .or. k1 < 1 &
                        .or. i1 > l .or. j1 > m .or. k1 > n) goto 55

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
                           goto 55
                        end if
                     end do

                     call distan(x1,y1,z1,x2,y2,z2, dis)
                     if ( dis < alf ) then
                        nsub = nsub + 1
                        itemp(nsub) = nn2
                        tempx(1) = x2 - x1
                        tempx(2) = y2 - y1
                        tempx(3) = z2 - z1
                        call matvec(3,3,t, tempx, tempy)
                        w1(1, nsub)  = 1.0
                        w1(2, nsub)  = tempy(1)
                        w1(3, nsub)  = tempy(2)
                        w1(4, nsub)  = tempy(3)
                        w1(5, nsub)  = 0.5*tempy(1) * tempy(1)
                        w1(6, nsub)  = 0.5*tempy(2) * tempy(2)
                        w1(7, nsub)  = 0.5*tempy(3) * tempy(3)
                        w1(8, nsub)  = tempy(1) * tempy(2)
                        w1(9, nsub)  = tempy(1) * tempy(3)
                        w1(10, nsub) = tempy(2) * tempy(3)
                        w1(11, nsub) = tempy(1)*tempy(1)*tempy(1)
                        w1(12, nsub) = tempy(1)*tempy(1)*tempy(2)
                        w1(13, nsub) = tempy(1)*tempy(1)*tempy(3)
                        w1(14, nsub) = tempy(1)*tempy(2)*tempy(2)
                        w1(15, nsub) = tempy(1)*tempy(2)*tempy(3)
                        w1(16, nsub) = tempy(1)*tempy(3)*tempy(3)
                        w1(17, nsub) = tempy(2)*tempy(2)*tempy(2)
                        w1(18, nsub) = tempy(2)*tempy(2)*tempy(3)
                        w1(19, nsub) = tempy(2)*tempy(3)*tempy(3)
                        w1(20, nsub) = tempy(3)*tempy(3)*tempy(3)
                        if ( nsub .eq. nq) goto 77
                     end if
                  end if  ! (idis == i2 .and. nn2 > 0)
                  55 continue
               end do  ! k1 = k0-i2, k0+i2
            end do  ! j1 = j0-i2, j0+i2
         end do ! i1 = i0-i2, i0+i2
      end do ! i2=0, if0
      77 continue

      if (nsub < 22) then
         write(*,*) "   nsub = ",nsub,"alf is too small in coede20f()!"
         stop
      end if

      ! -------- Call least square routine

      if (isvd == 1) then
         job = 11
         call dsvdc(w1,20,20,nsub,sd,ew,uw,20,v,nq,w4,job,inf)
      else
         call dgesvd('A','A',20,nsub,w1,20,sd2,uw,20,v,nq, &
               work, 1000, inf)
         do ij = 1, 20
            sd(ij) = sd2(ij)
         end do
      end if

      if (inf /= 0) then
         write(*,*) inf, " - ssvdc() or sgesvd() failed in coedef!"
         stop
      end if

      do i1=1, 20
         if (abs(sd(i1)) > 1.0e-14 ) then
            ew1(i1)=uw(1, i1)/sd(i1)
            ew2(i1)=uw(2, i1)/sd(i1)
            ew3(i1)=uw(3, i1)/sd(i1)
            ew4(i1)=uw(4, i1)/sd(i1)
            ew5(i1)=uw(5, i1)/sd(i1)
            ew6(i1)=uw(6, i1)/sd(i1)
            ew7(i1)=uw(7, i1)/sd(i1)
            ew8(i1)=uw(8, i1)/sd(i1)
            ew9(i1)=uw(9, i1)/sd(i1)
            ew10(i1)=uw(10, i1)/sd(i1)
         else
            ew1(i1)= 0.0
            ew2(i1)= 0.0
            ew3(i1)= 0.0
            ew4(i1)= 0.0
            ew5(i1)= 0.0
            ew6(i1)= 0.0
            ew7(i1)= 0.0
            ew8(i1)= 0.0
            ew9(i1)= 0.0
            ew10(i1)= 0.0
         end if
      end do


      ! -------- w31(i),w32(i),......,w36(i) are the solutions

      do i1=1, nsub
         w31(i1) = 0.0
         w32(i1) = 0.0
         w33(i1) = 0.0
         w34(i1) = 0.0
         w35(i1) = 0.0
         w36(i1) = 0.0
         w37(i1) = 0.0
         w38(i1) = 0.0
         w39(i1) = 0.0
         w310(i1) = 0.0

         if (isvd == 1) then
            do j1=1,20
               w31(i1) = w31(i1) + v(i1,j1)*ew1(j1)
               w32(i1) = w32(i1) + v(i1,j1)*ew2(j1)
               w33(i1) = w33(i1) + v(i1,j1)*ew3(j1)
               w34(i1) = w34(i1) + v(i1,j1)*ew4(j1)
               w35(i1) = w35(i1) + v(i1,j1)*ew5(j1)
               w36(i1) = w36(i1) + v(i1,j1)*ew6(j1)
               w37(i1) = w37(i1) + v(i1,j1)*ew7(j1)
               w38(i1) = w38(i1) + v(i1,j1)*ew8(j1)
               w39(i1) = w39(i1) + v(i1,j1)*ew9(j1)
               w310(i1) = w310(i1) + v(i1,j1)*ew10(j1)
            end do
         else
            do j1=1,20
               w31(i1) = w31(i1) + v(j1,i1)*ew1(j1)
               w32(i1) = w32(i1) + v(j1,i1)*ew2(j1)
               w33(i1) = w33(i1) + v(j1,i1)*ew3(j1)
               w34(i1) = w34(i1) + v(j1,i1)*ew4(j1)
               w35(i1) = w35(i1) + v(j1,i1)*ew5(j1)
               w36(i1) = w36(i1) + v(j1,i1)*ew6(j1)
               w37(i1) = w37(i1) + v(j1,i1)*ew7(j1)
               w38(i1) = w38(i1) + v(j1,i1)*ew8(j1)
               w39(i1) = w39(i1) + v(j1,i1)*ew9(j1)
               w310(i1) = w310(i1) + v(j1,i1)*ew10(j1)
            end do
         end if

         wcoe(ir,i1)   = w31(i1)
         wxcoe(ir,i1)  = w32(i1)
         wycoe(ir,i1)  = w33(i1)
         wzcoe(ir,i1)  = w34(i1)
         wxxcoe(ir,i1) = w35(i1)
         wyycoe(ir,i1) = w36(i1)
         wzzcoe(ir,i1) = w37(i1)
         wxycoe(ir,i1) = w38(i1)
         wxzcoe(ir,i1) = w39(i1)
         wyzcoe(ir,i1) = w310(i1)

      end do  !  i1=1, nsub


      if(ir < 5) then

         !        call coec20(ir, maxirr,nq,index2,cirreg, wcoe)
         !        call coec20(ir, maxirr,nq,index2,cirreg, wycoe)
         !        call coec20(ir, maxirr,nq,index2,cirreg, wzcoe)
         !        call coec20(ir, maxirr,nq,index2,cirreg, wyycoe)
         !        call coec20(ir, maxirr,nq,index2,cirreg, wzzcoe)
         !        call coec20(ir, maxirr,nq,index2,cirreg, wyzcoe)

      end if

   end do  ! ir = 1, nirreg

   return
end subroutine coed20 





!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine coec20 here]
subroutine coec20(h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg, coe)
   use iim_use
   implicit _REAL_ (a-h,o-z) !Not called.
   parameter ( nqq = 27 )

   common /lmn/l, m, n, nirreg
   !common /hxyz/hx,hy,hz,hmax
   !common /xyz/xs, xf, ys, yf, zs, zf
   _REAL_ xs,ys,zs

   dimension index2(l,m,n)
   dimension cirreg(maxirr, 15)
   dimension coe(maxirr,nq)

   dimension t(3,3)
   dimension tempx(3),tempy(3)
   dimension itemp(nqq)

   if (nq /= nqq) then
      write(*,*) "   Please change nqq in coec20() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 10.1*hmax
   !XP if0 = int(alf/hmax) + 1
   if0 = int(alf/hmax) + 4


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

   s1 = 0.0
   s2 = 0.0
   s3 = 0.0
   s4 = 0.0
   s5 = 0.0
   s6 = 0.0
   s7 = 0.0 
   s8 = 0.0 
   s9 = 0.0
   s10 = 0.0


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
                     tempx(1) = x2 - x1
                     tempx(2) = y2 - y1
                     tempx(3) = z2 - z1
                     call matvec(3,3,t, tempx, tempy)

                     s1 = s1 + coe(ir, nsub)
                     s2 = s2 + coe(ir, nsub)*tempy(1)
                     s3 = s3 + coe(ir, nsub)*tempy(2)
                     s4 = s4 + coe(ir, nsub)*tempy(3)
                     s5 = s5 + coe(ir, nsub)*0.5*tempy(1)*tempy(1)
                     s6 = s6 + coe(ir, nsub)*0.5*tempy(2) * tempy(2)
                     s7 = s7 + coe(ir, nsub)*0.5*tempy(3) * tempy(3)
                     s8 = s8 + coe(ir, nsub)*tempy(1)*tempy(2)
                     s9 = s9 + coe(ir, nsub)*tempy(1)*tempy(3)
                     s10 = s10 + coe(ir, nsub)*tempy(2) * tempy(3)
                     if ( nsub .eq. nq) goto 577
                  end if
               end if  ! (idis == i2 .and. nn2 > 0)
               555 continue
            end do  ! k1 = k0-i2, k0+i2
         end do  ! j1 = j0-i2, j0+i2
      end do  ! i1 = i0-i2, i0+i2
   end do  ! i2=0, if0
   577 continue

   write(*,*) s1
   write(*,*) s2,s3,s4
   write(*,*) s5,s6,s7
   write(*,*) s8,s9,s10
   write(*,*) "-------------------------------"

   return
end subroutine coec20 
