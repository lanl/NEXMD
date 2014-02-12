#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grtopr here]
subroutine grtopr(l,m,n,x,y,z,h,hx,hy,hz,xx,yy,zz,i0,j0,k0,qg,inf1,inf2,inf3, &
      q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)
   implicit none
   integer,parameter :: nq = 27
   
   !       **************************************************************
   !       *                                                            *
   !       *  GrToPr find:                                              *
   !       *     if inf1=1: the value q0 of a quantity q at (xx,yy,zz)  *
   !       *     if inf2=1: the value qx,qy,qz of  q at (xx,yy,zz)      *
   !       *     if inf3=1: the value qxx,qyy,qzz,qxy,qxz,qyz of q at   *
   !       *                (xx,yy,zz)                                  *
   !       *  based upon its values at grid points qg(*,*,*).           *
   !       *                                                            *
   !       *   im=1: 8-pint bilinear interpolation scheme used.         *
   !       *         (only used in case: inf1=1 only)                   *
   !       *   im=2: General interpolation scheme used.                 *
   !       *                                                            *
   !       **************************************************************

   !Passed variables:
   integer l,m,n,i0,j0,k0,inf1,inf2,inf3
   _REAL_  x(0:l+1),y(0:m+1),z(0:n+1),qg(0:l+1,0:m+1,0:n+1),h,hx,hy,hz,&
           q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz,xx,yy,zz

   !Local variables:
   integer im,isvd,jj,i1,j1,k1,i,ii,j,ij,nsub,job,inf
   _REAL_ x0b,x1b,y0b,y1b,x1,y1,z1,x0,y0,z0,z0b,z1b
   _REAL_ w1(10, nq), w4(10)
   _REAL_ sd(11), uw(10,10), v(nq, nq)
   _REAL_ ew(nq)
   _REAL_ ew1(nq), ew2(nq), ew3(nq)
   _REAL_ ew4(nq), ew5(nq), ew6(nq)
   _REAL_ ew7(nq), ew8(nq), ew9(nq), ew10(nq)
   _REAL_ w31(nq), w32(nq), w33(nq)
   _REAL_ w34(nq), w35(nq), w36(nq)
   _REAL_ w37(nq), w38(nq), w39(nq), w310(nq)
   _REAL_ qcoe(10,nq)
   _REAL_ sd2(10), work(1000)

   ! ----- set SVD subroutine: isvd=1  call ssvdc()
   !                           isvd=2  call sgecvd()

   isvd = 2
   im = 2

   !       i0 = nint((xx - x(0))/hx)
   !       j0 = nint((yy - y(0))/hy)
   !       k0 = nint((zz - z(0))/hz)


   if (im == 1 .and. inf2 == 0 .and. inf3 == 0) then
      i1 = i0 + 1
      j1 = j0 + 1
      k1 = k0 + 1
      x0 = x(i0)
      y0 = y(j0)
      z0 = z(k0)
      x1 = x0 + hx
      y1 = y0 + hy
      z1 = z0 + hz

      x0b = 2.0*(x1 - xx)/hx
      x1b = 2.0*(xx - x0)/hx
      y0b = 2.0*(y1 - yy)/hy
      y1b = 2.0*(yy - y0)/hy
      z0b = 2.0*(z1 - zz)/hz
      z1b = 2.0*(zz - z0)/hz

      q0 =   x0b * y0b * z0b * qg(i0,j0,k0) &
            + x0b * y0b * z1b * qg(i0,j0,k1) &
            + x0b * y1b * z0b * qg(i0,j1,k0) &
            + x0b * y1b * z1b * qg(i0,j1,k1) &
            + x1b * y0b * z0b * qg(i1,j0,k0) &
            + x1b * y0b * z1b * qg(i1,j0,k1) &
            + x1b * y1b * z0b * qg(i1,j1,k0) &
            + x1b * y1b * z1b * qg(i1,j1,k1)

      q0 = q0/8.0

   else

      do i=1, 10
         do j=1, nq
            qcoe(i,j) = 0.0
         end do
      end do

      ! -------- initialize

      do jj=1, nq
         do ii=1, nq
            v(ii, jj) = 0.0      !# the V matrix in SVD: U*S*V
         end do
         do ii=1, 10
            w1(ii, jj) = 0.0     !# the coefficient matrix to be SVDed
         end do
      end do

      ! -------- Generate the matrix

      nsub = 0
      do i1 = i0-1, i0+1
         do j1 = j0-1, j0+1
            do k1 = k0-1, k0+1
               if (i1 < 1 .or. j1 < 1 .or. k1 < 1 &
                     .or. i1 > l .or. j1 > m .or. k1 > n) goto 55

               nsub = nsub + 1
               w1(1, nsub)  = 1.0
               w1(2, nsub)  = x(i1) - xx
               w1(3, nsub)  = y(j1) - yy
               w1(4, nsub)  = z(k1) - zz
               w1(5, nsub)  = 0.5*(x(i1) - xx)*(x(i1) - xx)
               w1(6, nsub)  = 0.5*(y(j1) - yy)*(y(j1) - yy)
               w1(7, nsub)  = 0.5*(z(k1) - zz)*(z(k1) - zz)
               w1(8, nsub)  = (x(i1) - xx)*(y(j1) - yy)
               w1(9, nsub)  = (x(i1) - xx)*(z(k1) - zz)
               w1(10, nsub) = (y(j1) - yy)*(z(k1) - zz)
               55 continue
            end do
         end do
      end do
      if (nsub < 12) then
         write(*,*) "   nsub is too small in GrToPr()!"
         !            print *,xx,yy,zz
         stop
      end if

      ! -------- Call least square routine

      if (isvd == 1) then
         job = 11
         call dsvdc(w1,10,10,nsub,sd,ew,uw,10,v,nq,w4,job,inf)
      else
         call dgesvd('A','A',10,nsub,w1,10,sd2,uw,10,v,nq, &
               work, 1000, inf)
         do ij = 1, 10
            sd(ij) = sd2(ij)
         end do
      end if

      if (inf /= 0) then
         write(*,*) inf, " - ssvdc() or sgesvd() failed in GrToPr!"
         stop
      end if


      do i1=1, 10
         if (abs(sd(i1)) > 1.0e-14) then
            if (inf1 == 1) then
               ew1(i1)=uw(1, i1)/sd(i1)
            end if
            if (inf2 == 1) then
               ew2(i1)=uw(2, i1)/sd(i1)
               ew3(i1)=uw(3, i1)/sd(i1)
               ew4(i1)=uw(4, i1)/sd(i1)
            end if
            if (inf3 == 1) then
               ew5(i1)=uw(5, i1)/sd(i1)
               ew6(i1)=uw(6, i1)/sd(i1)
               ew7(i1)=uw(7, i1)/sd(i1)
               ew8(i1)=uw(8, i1)/sd(i1)
               ew9(i1)=uw(9, i1)/sd(i1)
               ew10(i1)=uw(10, i1)/sd(i1)
            end if
         else
            if (inf1 == 1) then
               ew1(i1)= 0.0
            end if
            if (inf2 == 1) then
               ew2(i1)= 0.0
               ew3(i1)= 0.0
               ew4(i1)= 0.0
            end if
            if (inf3 == 1) then
               ew5(i1)= 0.0
               ew6(i1)= 0.0
               ew7(i1)= 0.0
               ew8(i1)= 0.0
               ew9(i1)= 0.0
               ew10(i1)= 0.0
            end if
         end if  ! (abs(sd(i1)) > 1.0e-14)
      end do  !  i1=1, 10


      ! -------- w3(i) is the solution

      do i1=1, nsub
         if (inf1 == 1) then
            w31(i1) = 0.0
         end if
         if (inf2 == 1)  then
            w32(i1) = 0.0
            w33(i1) = 0.0
            w34(i1) = 0.0
         end if
         if (inf3 == 1) then
            w35(i1) = 0.0
            w36(i1) = 0.0
            w37(i1) = 0.0
            w38(i1) = 0.0
            w39(i1) = 0.0
            w310(i1) = 0.0
         end if

         do j1=1,10
            if (isvd == 1) then
               if (inf1 == 1) then
                  w31(i1) = w31(i1) + v(i1,j1)*ew1(j1)
               end if
               if (inf2 == 1) then
                  w32(i1) = w32(i1) + v(i1,j1)*ew2(j1)
                  w33(i1) = w33(i1) + v(i1,j1)*ew3(j1)
                  w34(i1) = w34(i1) + v(i1,j1)*ew4(j1)
               end if
               if (inf3 == 1) then
                  w35(i1) = w35(i1) + v(i1,j1)*ew5(j1)
                  w36(i1) = w36(i1) + v(i1,j1)*ew6(j1)
                  w37(i1) = w37(i1) + v(i1,j1)*ew7(j1)
                  w38(i1) = w38(i1) + v(i1,j1)*ew8(j1)
                  w39(i1) = w39(i1) + v(i1,j1)*ew9(j1)
                  w310(i1) = w310(i1) + v(i1,j1)*ew10(j1)
               end if
            else
               if (inf1 == 1) then
                  w31(i1) = w31(i1) + v(j1,i1)*ew1(j1)
               end if
               if (inf2 == 1) then
                  w32(i1) = w32(i1) + v(j1,i1)*ew2(j1)
                  w33(i1) = w33(i1) + v(j1,i1)*ew3(j1)
                  w34(i1) = w34(i1) + v(j1,i1)*ew4(j1)
               end if
               if (inf3 == 1) then
                  w35(i1) = w35(i1) + v(j1,i1)*ew5(j1)
                  w36(i1) = w36(i1) + v(j1,i1)*ew6(j1)
                  w37(i1) = w37(i1) + v(j1,i1)*ew7(j1)
                  w38(i1) = w38(i1) + v(j1,i1)*ew8(j1)
                  w39(i1) = w39(i1) + v(j1,i1)*ew9(j1)
                  w310(i1) = w310(i1) + v(j1,i1)*ew10(j1)
               end if
            end if  ! (isvd == 1)
         end do  !  j1=1,10

         if (inf1 == 1) then
            qcoe(1,i1) = w31(i1)
         end if
         if (inf2 == 1) then
            qcoe(2,i1) = w32(i1)
            qcoe(3,i1) = w33(i1)
            qcoe(4,i1) = w34(i1)
         end if
         if (inf3 == 1) then
            qcoe(5,i1) = w35(i1)
            qcoe(6,i1) = w36(i1)
            qcoe(7,i1) = w37(i1)
            qcoe(8,i1) = w38(i1)
            qcoe(9,i1) = w39(i1)
            qcoe(10,i1) = w310(i1)
         end if

      end do  !  i1=1, nsub

      q0  = 0
      qx  = 0
      qy  = 0
      qz  = 0
      qxx = 0
      qyy = 0
      qzz = 0
      qxy = 0
      qxz = 0
      qyz = 0

      nsub = 0
      do i1 = i0-1, i0+1
         do j1 = j0-1, j0+1
            do k1 = k0-1, k0+1
               if (i1 < 1 .or. j1 < 1 .or. k1 < 1 &
                     .or. i1 > l .or. j1 > m .or. k1 > n) goto 155
               if (qg(i1,j1,k1) > 999.0 ) then
                  print *,'unexcepted grid',i1,j1,k1,qg(i1,j1,k1)
                  print *,i0,j0,k0
                  print *,xx,yy,zz
                  stop
               end if

               nsub = nsub + 1
               if (inf1 == 1) then
                  q0  =  q0  + qcoe(1,nsub)*qg(i1,j1,k1)
               end if
               if (inf2 == 1) then
                  qx  =  qx  + qcoe(2,nsub)*qg(i1,j1,k1)
                  qy  =  qy  + qcoe(3,nsub)*qg(i1,j1,k1)
                  qz  =  qz  + qcoe(4,nsub)*qg(i1,j1,k1)
               end if
               if (inf3 == 1) then
                  qxx =  qxx + qcoe(5,nsub)*qg(i1,j1,k1)
                  qyy =  qyy + qcoe(6,nsub)*qg(i1,j1,k1)
                  qzz =  qzz + qcoe(7,nsub)*qg(i1,j1,k1)
                  qxy =  qxy + qcoe(8,nsub)*qg(i1,j1,k1)
                  qxz =  qxz + qcoe(9,nsub)*qg(i1,j1,k1)
                  qyz =  qyz + qcoe(10,nsub)*qg(i1,j1,k1)
               end if
               155 continue
            end do
         end do
      end do

   end if

   return
end subroutine grtopr 

