#define _REAL_ double precision
  ![ Calculate the unin for the target pnt on the interface]
      subroutine interpx(l,m,n,bi,bo,n3,i0,j0,k0,sgni,sgnj,sgnk,&
          nbnd,xp,yp,zp,uu,dudx,dudy,dudz, &
                 wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp, & 
                 xyy,xzz,xyz,t,u,phi,xs,ys,zs,h,&
                 index2,nn1)
     use iim_use
      implicit none
   
     integer, parameter :: nq=27

      !common /lmn/ l1,m1,n1,ndum
      !common /hxyz/hx,hy,hz,hmax
      !common /para/bi,bo

   ! passed variables

   ! l,m,n                   box dimension
   ! n3                      when n3=10 2nd order
   ! i0,j0,k0                projection pnts in grid crds
   ! sgni,sgnj,sgnk          dummy variables here
   ! nbnd                    initial dimension without clustering
   ! xp,yp,zp                projection crds
   ! uu,dudx,dudy,dudz       potential and derivatives of the pnt on interface
   ! phi,u                   lvlset function; potential of whole grids
   ! wp,qp                   [u]and[betaUn]for the interface
   ! wp,wyp,wzp,wyyp,        jump conditions and derivatives of [u]
   ! wyzp,wzzp
   ! qp,qyp,qzp              jump conditions and derivatives of [beta*u]
   ! xyy,xzz,xyz             derivative of the geometry on the pnt on interface
   ! t                       crds trans matrix
   ! index2                  the index for the irregular point
   ! nn1                     the number for the target pnt on interface 
   
   ! local variables

   ! xloc                    local crds distance 
   ! rho,rho1                ratio of the dielectric coefficients 

   !passed variables
   _REAL_  q0p(nbnd),qyp(nbnd),qzp(nbnd),w0p(nbnd),wyp(nbnd),wzp(nbnd),&
            wyyp(nbnd),wzzp(nbnd),wyzp(nbnd),&
            wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
            wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
            wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq),&
            sss1(nbnd),sss2(nbnd),c(l,m,n,7),c2(nbnd, 27),&
            wp(1:nbnd),qp(1:nbnd)
    integer index2(l,m,n),nn1




      integer l,m,n,ndum,l1,m1,n1,nbnd,ir
      _REAL_ hx,hy,hz,hmax
 
      _REAL_ x(0:l+1), y(0:m+1), z(0:n+1)
      _REAL_ phi(0:l+1,0:m+1, 0:n+1)
      _REAL_ u(1:l,1:m, 1:n)
      _REAL_ uu,dudx,dudy,dudz
      _REAL_ xp,yp,zp
      _REAL_ x1,y1,z1,xs,ys,zs,h
      integer i0,j0,k0,sgni,sgnj,sgnk
      !cq
      _REAL_ ds,dr,ndir(3),alf
      _REAL_ t(3,3),a(20)
      _REAL_ xloc(3),tmpx(3)
      _REAL_ tmp1(3),tmp2(3),tmpa(3,3),tmpb(3,3),tmpc(3,3)
      _REAL_ bi,bo
      _REAL_ bl_out(3),bl_in(3)
      _REAL_ xyy,xzz,xyz
      _REAL_ fjmp,fkjmp,ujmp
      _REAL_ rho,rho1
      _REAL_ bl2jmp,bl3jmp,temp
      integer nc,i,j,k,job,n3,nsub,info
      _REAL_ wy,wz,wyy,wyz,wzz,q,qy,qz
      _REAL_ fw,fq,fbi,fbo

      integer, parameter :: nn = 100
      _REAL_,dimension(1:n3,1:nn) :: ca
      _REAL_,dimension(1:nn) :: b
      _REAL_,dimension(1:nn) :: ew,w4
      _REAL_,dimension(1:nn) :: w3ux,w3uy,w3uz,w3u
      _REAL_,dimension(1:nn) :: ewux,ewuy,ewuz,ewu
      _REAL_,dimension(1:n3,1:n3) :: uw
      _REAL_,dimension(1:nn,1:nn) :: vl
      _REAL_,dimension(1:n3+1) :: sd
        
   
      l1 = l; m1 = m; n1 = n
      hx = h; hy = h; hz = h; hmax = h
      ! XP: build a larger grid, can solve problem concerning the membrance case
      ! of Wes using PBC
      do i = 0, l+1
         x(i) = xs + i*hx
      end do
      do j = 0, m+1
         y(j) = ys + j*hy
      end do
      do k = 0, n+1
         z(k) = zs + k*hz
      end do
      x1=xp ! (xp,yp,zp) is already in the lab frame
      y1=yp
      z1=zp
      !x1=xs+xp*h ! (xp,yp,zp) is already in the lab frame
      !y1=ys+yp*h
      !z1=zs+zp*h

   !  print *,l,m,n,n3
   !  print *,i0,j0,k0
   !  print *,x1,y1,z1
   !  print *,c
   !  print *,c1
   !  print *,c2
   !  print *,en
   !  print *,en1

      !call fjmps(x1,y1,z1, fjmp)
      fjmp = 0.d0
           
      !call fkjmps(x1,y1,z1, fkjmp)
      fkjmp = 0.d0
   ! get ir
      ir=nn1
       
      ! ir=index2(i0,j0,k0)
      ujmp = wp(ir)
      wy   = wyp(ir)
      wz   = wzp(ir)
      wyy  = wyyp(ir)
      wzz  = wzzp(ir)
      wyz  = wyzp(ir)
   
      q =  qp(ir)
      qy = qyp(ir)
      qz = qzp(ir)
      
      bl_out=0.0d0;bl_in=0.0d0
 
      rho = bi/bo
      rho1 = 1 - rho
      bl2jmp = bl_out(2)-bl_in(2)
      bl3jmp = bl_out(3)-bl_in(3)
      temp = fkjmp/bo
   
   !  ds = 1.d0/sqrt((atmctr(1)-x1)**2+(atmctr(2)-y1)**2 &
   !       +(atmctr(3)-z1)**2)
   !  ndir(1) = (atmctr(1)-x1)*ds
   !  ndir(2) = (atmctr(2)-y1)*ds
   !  ndir(3) = (atmctr(3)-z1)*ds

      alf = 1.6d0*h

      nc = 0
   !  if ( i0 == xp .or. j0 == yp .or. k0 == zp ) write(10,*) 'xyzp',xp,yp,zp
   !  if ( i0 == xp .or. j0 == yp .or. k0 == zp ) write(10,*) 'ijk0',i0,j0,k0
      do i=i0-1,i0+1
         do j=j0-1,j0+1
            do k=k0-1,k0+1
 !              if ( (i-i0)**2 + (j-j0)**2 + (k-k0)**2 > 1 )  cycle
   !           if ( i < floor(xp+sgni*1.d-6) .or. i > ceiling(xp+sgni*1.d-6) ) cycle
   !           if ( j < floor(yp+sgnj*1.d-6) .or. j > ceiling(yp+sgnj*1.d-6) ) cycle
   !           if ( k < floor(zp+sgnk*1.d-6) .or. k > ceiling(zp+sgnk*1.d-6) ) cycle

   !           if ( i0 == xp .or. j0 == yp .or. k0 == zp ) write(10,*) i,j,k,phi(i,j,k)
   !        print "(3f10.5,2f20.10)",x(i),y(j),z(k),phi(i,j,k),u(i,j,k)
               tmpx(1) = x(i) - x1
               tmpx(2) = y(j) - y1
               tmpx(3) = z(k) - z1

               ds=sqrt(tmpx(1)**2+tmpx(2)**2+tmpx(3)**2)
      
               !XP: to make sure more points are used to interpolate.
               if ( ds > alf) cycle
               nc = nc + 1

   !           if ( ds == 0.d0 ) then
   !              ds = 1.d0+dcos(3.1415926535897932d0*ds/alf)
   !              ds = exp(-ds**2/0.1d0)
   !           else
   !              not tested for outside interpolation
   !              dr = abs(dot_product(tmpx,ndir))/ds
   !              ds = (1.d0+dr)* &
   !                   (1.d0+dcos(3.1415926535897932d0*ds/alf))
   !              ds = exp(-(dr-1)**2/0.1d0)*exp(-ds**2/0.1d0)
   !           end if

               ds = 1.d0+dcos(3.1415926535897932d0*ds/alf)
               ds = 1.d0
               call matvec(3,3,t,tmpx,xloc)
   !           print *,'coord',tmpx
  
            
               if (phi(i,j,k) .le. 0.d0) then       
                 ca(1,nc)  = 1.0d0*ds
                 ca(2,nc)  = xloc(1)*ds
                 ca(3,nc)  = xloc(2)*ds
                 ca(4,nc)  = xloc(3)*ds
                 if ( n3 == 10 ) then
                    ca(5,nc)  = 0.5d0*xloc(1)*xloc(1)*ds
                    ca(6,nc)  = 0.5d0*xloc(2)*xloc(2)*ds
                    ca(7,nc)  = 0.5d0*xloc(3)*xloc(3)*ds
                    ca(8,nc)  = xloc(1)*xloc(2)*ds
                    ca(9,nc)  = xloc(1)*xloc(3)*ds
                    ca(10,nc) = xloc(2)*xloc(3)*ds
                 end if
                 b(nc)     = u(i,j,k)*ds
               else
                 ca(1,nc)  = 1.0d0*ds!-0.5*temp*xloc(1)*xloc(1)
           !     ca(2,nc)  = rho*xloc(1)*ds
           !     ca(3,nc)  = xloc(2)*ds
           !     ca(4,nc)  = xloc(3)*ds
                     ca(2,nc)  = rho*xloc(1) &
                               - 0.5*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                               - bl_in(1)/bo + rho*bl_out(1)/bo) &     !XP: bl_in()
                               + 0.5*xloc(2)*xloc(2)*xyy*rho1 &        !1st der
                               + 0.5*xloc(3)*xloc(3)*xzz*rho1 &
                               + xloc(1)*xloc(2)*(bl_in(2)/bo &
                               - rho*bl_out(2)/bo) &
                               + xloc(1)*xloc(3)*(bl_in(3)/bo &     
                               - rho*bl_out(3)/bo) &
                               + xloc(2)*xloc(3)*xyz*rho1
                     ca(3,nc)  = - 0.5*xloc(1)*xloc(1)*bl2jmp/bo &
                               + xloc(1)*xloc(2)*xyy*rho1 &
                               + xloc(1)*xloc(3)*xyz*rho1 &
                               + xloc(2)     
                     ca(4,nc)  = - 0.5*xloc(1)*xloc(1)*bl3jmp/bo &
                               + xloc(1)*xloc(3)*xzz*rho1 &
                               + xloc(1)*xloc(2)*xyz*rho1 &
                               + xloc(3)     
                 if ( n3 == 10 ) then
                    ca(5,nc)  = 0.5d0*xloc(1)*xloc(1)*rho*ds
                    ca(6,nc)  = 0.5d0*xloc(2)*xloc(2)*ds &
                              + 0.5d0*(rho-1)*xloc(1)*xloc(1)*ds
                    ca(7,nc)  = 0.5d0*xloc(3)*xloc(3)*ds &
                              + 0.5d0*(rho-1)*xloc(1)*xloc(1)*ds
                    ca(8,nc)  = xloc(1)*xloc(2)*rho*ds
                    ca(9,nc)  = xloc(1)*xloc(3)*rho*ds
                    ca(10,nc) = xloc(2)*xloc(3)*ds
                 end if
                   b(nc)    = u(i,j,k)*ds - ujmp*ds &
                           - xloc(1)*q/bo*ds &
                           - xloc(2)*wy*ds &
                           - xloc(3)*wz*ds &
                              - xloc(2)*xloc(3)*(-q/bo*xyz+wyz) &
                              - 0.5*xloc(2)*xloc(2)*(-q/bo*xyy+wyy) &
                              - 0.5*xloc(3)*xloc(3)*(-q/bo*xzz+wzz) &
                              - xloc(1)*xloc(2)*(wy*xyy+wz*xyz+qy/bo) &
                              - xloc(1)*xloc(3)*(wy*xyz+wz*xzz+qz/bo) &
                              - 0.5*xloc(1)*xloc(1)* &
                                             (q/bo*(xyy+xzz)-wyy-wzz)
!                print *,'xloc',xloc(1:3),bo,ds
!                print *,'out',u(i,j,k),ujmp,q,wy,wz
               end if
!              print *,'ca',phi(i,j,k),ca(1:4,nc),b(nc)
!  !           write(100,"(5f15.8)") ca(1:4,nc),b(nc)
            end do
         end do
      end do
 
   !print *,'nc=',nc
   if ( nc < 4 ) then
   write(6,*) 'pb bomb: nc < 4'
   stop
   end if

       job = 11
   !   print *,'n2=',n2
   !   print *,'n3=',n3
   !   n3 = 10
   !   nsub = 27
   !   nn = 27
       nsub = nc 
   !   print *,'nsub=',nsub
       call dsvdc(ca,n3,n3,nsub,sd,ew,uw,n3,vl,nn,w4,job,info)

       do i=1,n3
          if ( sd(i) .gt. 1.0d-12) then
          ewu(i) = uw(1,i)/sd(i)
          ewux(i)= uw(2,i)/sd(i)
          ewuy(i)= uw(3,i)/sd(i)
          ewuz(i)= uw(4,i)/sd(i)
          else
          ewu(i) = 0.0d0
          ewux(i)= 0.0d0
          ewuy(i)= 0.0d0
          ewuz(i)= 0.0d0
          end if
       enddo
   !   do i=1,n3-1
   !      ewu(i) = uw(1,i)/sd(i)
   !      ewux(i)= uw(2,i)/sd(i)
   !      ewuy(i)= uw(3,i)/sd(i)
   !      ewuz(i)= uw(4,i)/sd(i)
   !   enddo
   !   if(sd(n3) .gt. 1.0d-12) then
   !      ewu(n3) = uw(1,n3)/sd(n3)
   !      ewux(n3) = uw(2,n3)/sd(n3)
   !      ewuy(n3) = uw(3,n3)/sd(n3)
   !      ewuz(n3) = uw(4,n3)/sd(n3)
   !   else
   !      ewu(n3) = 0.0d0
   !      ewux(n3) = 0.0d0
   !      ewuy(n3) = 0.0d0
   !      ewuz(n3) = 0.0d0
   !   endif
   !   print *,v

   !   do i=1,n2
       do i=1,nsub
          w3u(i)  = 0.0d0
          w3ux(i) = 0.0d0
          w3uy(i) = 0.0d0
          w3uz(i) = 0.0d0
          do j=1,n3
             w3u(i)  = w3u(i)  + vl(i,j)*ewu(j)
             w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
             w3uy(i) = w3uy(i) + vl(i,j)*ewuy(j)
             w3uz(i) = w3uz(i) + vl(i,j)*ewuz(j)
          enddo
   !      print *,i,v
       enddo
   !   print *,v

       uu   = 0.0d0
       dudx = 0.0d0
       dudy = 0.0d0
       dudz = 0.0d0
   !   do i = 1, n2
       do i = 1, nsub
          uu   = uu   +  w3u(i) *b(i)
          dudx = dudx +  w3ux(i)*b(i)
          dudy = dudy +  w3uy(i)*b(i)
          dudz = dudz +  w3uz(i)*b(i)
       end do
!      print *, 'local',uu,dudx,dudy,dudz

   !   tmp1(1) = dudx
   !   tmp1(2) = dudy
   !   tmp1(3) = dudz
   !   call  transp(3,3,t,tmpa)
   !   call  matvec(3,3,tmpa,tmp1,tmp2)
   !   dudx = tmp2(1)
   !   dudy = tmp2(2)
   !   dudz = tmp2(3)

       end subroutine
