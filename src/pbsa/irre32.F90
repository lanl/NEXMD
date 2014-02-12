#define _REAL_ double precision
! ----- for an irregular point (i0,j0,k0), find corresponding coefficients


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine irre32 here]
subroutine irre32(l,m,n,h,hx,hy,hz,ifail, &
      i0,j0,k0,info,bmax,b_in,b_out,x,y,z,phi,index, &
      q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp, &
      nq,maxirr,index2,cirreg,wp,qp,xx,rhs)
use iim_use

   implicit none
   integer,parameter :: im=4, ime=4, immax=im, ns=27, in=ns, &
         inmax=in, imnn=im+2*in, &
         lwar=3*inmax*inmax/2+10*inmax+2*immax+1, &
         liwar=in

   !Passed variables
   _REAL_   h,hx,hy,hz,bmax,b_in,b_out,rhs
   integer  l,m,n,ifail,i0,j0,k0,info,nq,maxirr
   integer  index(l,m,n),index2(l,m,n)
   _REAL_   x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_   phi(0:l+1,0:m+1, 0:n+1)
   _REAL_   cirreg(maxirr, 15)
   _REAL_   wp(maxirr), qp(maxirr)
   _REAL_   q0p(maxirr),qyp(maxirr),qzp(maxirr)
   _REAL_   w0p(maxirr),wyp(maxirr),wzp(maxirr)
   _REAL_   wyyp(maxirr),wzzp(maxirr),wyzp(maxirr)

   !Local variables
   integer  nc,kr,kctr,ndis,ji,i,j,k,ir,iprint,ii,iout,ijk
   _REAL_   xyy,xzz,xyz,x1,y1,z1,ujmp,wy,wz,wyy,wzz,wyz,q,qy,qz,&
            fk_out,ff_out,fk_in,ff_in,tmp,temp,tmp1,tmp2,tmp3,rho,&
            rho1,hxx,hyy,hzz,hxyz,fb_out,fb_in,fmin,h1,bl2jmp,bl3jmp,&
            fjmp,corr,eps,fkk_in,fkk_out,fkk,fkjmp
            
            
   _REAL_   t(3,3),a(20)
   _REAL_   xloc(3),tmpx(3)
   _REAL_   bl_out(3),bl_in(3)
   _REAL_   cc(inmax,inmax), ca(immax,inmax), caa(immax,inmax)
   _REAL_   cd(inmax),cb(immax),xx(in),xl(in),xu(in),uu(imnn)
   _REAL_   war(lwar),iwar(liwar)
   _REAL_   tmpcoe(7)

   ir = index2(i0, j0, k0)

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

   ! ----- find jump conditions

   call fjmps(x1,y1,z1, fjmp)

   call fkjmps(x1,y1,z1, fkjmp)

   ujmp = wp(ir)
   wy   = wyp(ir)
   wz   = wzp(ir)
   wyy  = wyyp(ir)
   wzz  = wzzp(ir)
   wyz  = wyzp(ir)

   q =  qp(ir)
   qy = qyp(ir)
   qz = qzp(ir)


   call betas(hx,hy,hz,x1,y1,z1,t,bl_in,bl_out,b_in,b_out)
   fkk_in = fk_in(x1,y1,z1)
   fkk_out = fk_out(x1,y1,z1)

   hxx=hx*hx
   hyy=hy*hy
   hzz=hz*hz
   h1=fmin(hxx,fmin(hyy,hzz))

   eps=0.1e-11

   do i=1,immax
      do j=1,inmax
         ca(i,j) = 0.0
      end do
      cb(i) = 0.0
   end do

   nc = 0
   kr = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            kr = kr +1
            if (i == i0 .and. j == j0 .and. k == k0) then
               !               kr = kr -1
               kctr = nc
               !               ca(kr+ime,nc) = -1.0d0
               !               cb(kr+ime) = -1.0d-20
            else
               !               ca(kr+ime,nc) = 1.0d0
               !               cb(kr+ime) = -100.0d0
            end if
            xl(nc) = 0.0
            xu(nc) = 2.0*bmax/h1
         end do
      end do
   end do

   xl(kctr)=-12.0*bmax/h1
   xu(kctr)=0.0

   rho = b_in/b_out
   rho1 = 1 - rho
   bl2jmp = bl_out(2)-bl_in(2)
   bl3jmp = bl_out(3)-bl_in(3)
   temp = fkjmp/b_out

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1

            ! ----------- local coordinates transformation

            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            ! ----------- form the coefficients of equality constraints

            if (index(i,j,k) <= 3) then
               ca(1,nc)  = 1.0
               ca(2,nc)  = xloc(1)
               ca(3,nc)  = xloc(2)
               ca(4,nc)  = xloc(3)
            else
               ca(1,nc)  = 1.0-0.5*temp*xloc(1)*xloc(1)
               ca(2,nc)  = rho*xloc(1) &
                     - 0.5*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                     - bl_in(1)/b_out + rho*bl_out(1)/b_out) &
                     + 0.5*xloc(2)*xloc(2)*xyy*rho1 &
                     + 0.5*xloc(3)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*(bl_in(2)/b_out &
                     - rho*bl_out(2)/b_out) &
                     + xloc(1)*xloc(3)*(bl_in(3)/b_out &
                     - rho*bl_out(3)/b_out) &
                     + xloc(2)*xloc(3)*xyz*rho1
               ca(3,nc)  = - 0.5*xloc(1)*xloc(1)*bl2jmp/b_out &
                     + xloc(1)*xloc(2)*xyy*rho1 &
                     + xloc(1)*xloc(3)*xyz*rho1 &
                     + xloc(2)
               ca(4,nc)  = - 0.5*xloc(1)*xloc(1)*bl3jmp/b_out &
                     + xloc(1)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*xyz*rho1 &
                     + xloc(3)
            end if
         end do  !  k=k0-1,k0+1
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1

   if (info > 3) then
      cb(1) = fkjmp
   end if

   cb(2)  = -bl_in(1)
   cb(3)  = -bl_in(2)
   cb(4)  = -bl_in(3)

   do i=1,inmax
      do j=1,inmax
         cc(i,j) = 0.0
      end do
      cc(i,i) = 1.0
   end do

   hxyz=hx*hx+hy*hy+hz*hz
   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            ndis=abs(i-i0)+abs(j-j0)+abs(k-k0)
            if (ndis <= 1) then
               if (phi(i,j,k) <= 0.0) then
                  cd(nc)=-3.0*fb_in(b_in,b_out,x(i),y(j),z(k))/hxyz
               else
                  cd(nc)=-3.0*fb_out(b_in,b_out,x(i),y(j),z(k))/hxyz
               end if
            else
               cd(nc)=0.0
            end if

            if (ndis == 0) cd(nc)=-6.0*cd(nc)
         end do
      end do
   end do

   do ijk=1, in
      !           ca(11, ijk) = 1.0
   end do
   !        cb(11) = 0.0

   if (info > 3) then
      call cpymat(immax,inmax,ca,caa)
   end if

   iprint = 1
   iwar(1) = 1


   ! ----- solve quadratic programming

   call ql0001(im,ime,immax,in,inmax,imnn,cc,cd,ca,cb,xl,xu, &
         xx,uu,iout,ifail,iprint,war,lwar,iwar,liwar,eps)

   tmp=0.0
   do ii=1,4
      tmp1=0.0
      do ji=1,inmax
         tmp1=tmp1+ca(ii,ji)*xx(ji)
      end do
      if(abs(tmp1+cb(ii)) > tmp) tmp = abs(tmp1+cb(ii))
   end do

   if(ifail > 10) then
      write(*,*) "irre32-Warning!",i0,j0,k0,"Inconsistent Constraints"
      do i=1, in
         xx(i) = 0.0
      end do
      call regula(l,m,n,hx,hy,hz,b_in,b_out,i0,j0,k0,info,x,y,z, tmpcoe,rhs)
      xx(5)  = tmpcoe(2)
      xx(23) = tmpcoe(3)
      xx(11) = tmpcoe(4)
      xx(17) = tmpcoe(5)
      xx(13) = tmpcoe(6)
      xx(15) = tmpcoe(7)
      xx(14) = tmpcoe(1)
      return
   end if

   do i=1,20
      a(i) = 0.0
   end do

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1

            nc = nc + 1
            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            if (index(i,j,k) <= 3) then
               a(1) =a(1) +xx(nc)
               a(3) =a(3) +xx(nc)*xloc(1)
               a(5) =a(5) +xx(nc)*xloc(2)
               a(7) =a(7) +xx(nc)*xloc(3)
               a(9) =a(9) +xx(nc)*xloc(1)*xloc(1)*0.5
               a(11)=a(11)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(13)=a(13)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(15)=a(15)+xx(nc)*xloc(1)*xloc(2)
               a(17)=a(17)+xx(nc)*xloc(1)*xloc(3)
               a(19)=a(19)+xx(nc)*xloc(2)*xloc(3)
            else
               a(2) =a(2) +xx(nc)
               a(4) =a(4) +xx(nc)*xloc(1)
               a(6) =a(6) +xx(nc)*xloc(2)
               a(8) =a(8) +xx(nc)*xloc(3)
               a(10)=a(10)+xx(nc)*xloc(1)*xloc(1)*0.5
               a(12)=a(12)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(14)=a(14)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(16)=a(16)+xx(nc)*xloc(1)*xloc(2)
               a(18)=a(18)+xx(nc)*xloc(1)*xloc(3)
               a(20)=a(20)+xx(nc)*xloc(2)*xloc(3)
            end if
         end do
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1

   tmp  = 1.0/b_out
   tmp1 =   a(4) + a(10)*(xyy+xzz-bl_out(1)*tmp) &
         - a(12)*xyy - a(14)*xzz - a(16)*bl_out(2)*tmp &
         - a(18)*bl_out(3)*tmp - a(20)*xyz
   tmp2 =   a(6) - a(10)*bl_out(2)*tmp &
         + a(16)*xyy + a(18)*xyz
   tmp3 =   a(8) - a(10)*bl_out(3)*tmp &
         + a(16)*xyz + a(18)*xzz

   corr =   a(10)*((fjmp-fkk_out*ujmp)*tmp-wyy-wzz) &
         + a(12)*wyy + a(14)*wzz + a(16)*qy*tmp &
         + a(18)*qz*tmp + a(20)*wyz + a(2)*ujmp &
         + tmp*tmp1*q + tmp2*wy + tmp3*wz

   !.......if the grid point is on (+) side

   if (info > 3) then
      corr = corr + fkk_out*ujmp - fjmp
   end if

   if (info <= 3) then
      fkk = fk_in(x(i0),y(j0),z(k0))
      rhs = ff_in(x(i0),y(j0),z(k0))
   else
      fkk = fk_out(x(i0),y(j0),z(k0))
      rhs = ff_out(x(i0),y(j0),z(k0))
   end if

   xx(kctr) = xx(kctr)+fkk
   rhs = rhs + corr

   return
end subroutine irre32 
