#define _REAL_ double precision

! This will be merged into pb_iimdrv()
!
subroutine IIM(gox,goy,goz,l,m,n,lvlst,insas,nbnd,iepsav,epsin,epsout,bv,u,u0,accept,h,&
               atmfirst,atmlast,bcopt,solvopt)
    
!    xs <- gox, ys <- goy, zs <- goz
!    l  <- xm , m  <- ym , n  <- zm
!    bi <- epsin, bo <- epsout
!    u  <- phi, u0 <- xs
    
   implicit none
    
   ! passed variables
      
   _REAL_  gox,goy,goz,h,accept,epsout,epsin
   integer l,m,n,nbnd,atmfirst,atmlast,bcopt,solvopt
   integer insas(0:l+1,0:m+1,0:n+1),iepsav(1:4,1:l*m*n)
   _REAL_  bv(l,m,n),u0(l,m,n),u(l,m,n)
   _REAL_  lvlst(0:l+1,0:m+1,0:n+1)
    
   ! common variables
    
   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo
   integer  nirreg
  !common /lmn/ l, m, n, nirreg
  !common /xyz/ xs, xf, ys, yf, zs, zf
  !common /hxyz/ hx,hy,hz,hmax
  !common /para/ bi,bo

#  include "pb_constants.h"

   ! local variables
      
   _REAL_, parameter :: eps0 = 8.8542D-12 / (1.6022D-19)**2 /  &
                               (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer, parameter :: nq=27
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ beta_max,rhs
   _REAL_ coe1(7), coe2(27)
   _REAL_ fw
   
    integer it,jt,kt,iit,it1,jt1,kt1,iit1
   _REAL_,allocatable :: x(:),y(:),z(:),cirreg(:,:)
   integer,allocatable :: index(:,:,:),index2(:,:,:)
   _REAL_,allocatable :: phi(:,:,:)
   _REAL_,allocatable :: wp(:),qp(:)
   _REAL_,allocatable :: q0p(:),qyp(:),qzp(:)
   _REAL_,allocatable :: w0p(:),wyp(:),wzp(:)
   _REAL_,allocatable :: wyyp(:),wzzp(:),wyzp(:)
   _REAL_,allocatable :: wcoe(:,:),wxcoe(:, :),wycoe(:, :)
   _REAL_,allocatable :: wzcoe(:,:),wxxcoe(:, :),wyycoe(:, :)
   _REAL_,allocatable :: wzzcoe(:, :),wxycoe(:, :)
   _REAL_,allocatable :: wxzcoe(:, :),wyzcoe(:, :)
   _REAL_,allocatable :: sss1(:),sss2(:)
   _REAL_,allocatable :: c(:,:,:,:), c2(:, :)
   _REAL_,allocatable :: f(:,:,:)

   ! interfacing variables 

   nirreg = nbnd
   xs = gox; ys = goy; zs = goz
   hx = h; hy = h; hz = h; hmax = h
   bi = epsin; bo = epsout

   ! allocating working arrays

   allocate (cirreg(1:nbnd,1:15)) ! geometrical parameters at irregular points
   allocate (index(1:l,1:m,1:n), index2(1:l,1:m,1:n))
   allocate (x(0:l+1), y(0:m+1), z(0:n+1))
   allocate (phi(0:l+1,0:m+1,0:n+1))
   allocate (wp(nbnd), qp(nbnd))  ! jump conditions at irregular points
   allocate (q0p(nbnd),qyp(nbnd),qzp(nbnd)) ! tangential derivatives of field jump conditions
   allocate (w0p(nbnd),wyp(nbnd),wzp(nbnd)) ! tangential derivatives of potential jump conditions
   allocate (wyyp(nbnd),wzzp(nbnd),wyzp(nbnd)) ! second derivatives of potential jump conditions
   allocate (wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd, nq)) ! coefficients obtained by SVD
   allocate (wzcoe(nbnd,nq),wxxcoe(nbnd, nq),wyycoe(nbnd, nq))
   allocate (wzzcoe(nbnd, nq),wxycoe(nbnd, nq))
   allocate (wxzcoe(nbnd, nq),wyzcoe(nbnd, nq))
   allocate (sss1(nbnd),sss2(nbnd)) !dummy arrays
   allocate (c(l,m,n,7), c2(nbnd, 27))
   allocate (f(l,m,n))
    
   ! setting coordinates of grid points
    
   do i = 0, l+1
      x(i) = xs + i*hx
     if(x(i)==1.0d0) then
       it=i
    end if
    if(x(i)==1.5d0) it1=i
   end do
   do j = 0, m+1
      y(j) = ys + j*hy
    if(y(j)==0.0d0) then
      jt=j
   end if
   end do
   do k = 0, n+1
      z(k) = zs + k*hz
   if(z(k)==0.0d0) then
    kt=k
   end if
   end do
   !write(99,*) x(1:l),y(1:m),z(1:n)
   
   ! set up the level set function

   phi(1:l,1:m,1:n) = lvlst(1:l,1:m,1:n)
   !do k =1, n; do j =1, m; do i = 1, l
   !   write(99,*) i,j,k,phi(i,j,k)
   !end do;end do;end do
   
   ! setting up index
   ! index = 1 for interior regular grid points;
   ! index = 2 for interior irregular grid points;
   ! index = 3 for boundary grid points right on the interface;
   ! index = 4 for exterior irregular grid points;
   ! index = 5 for exterior regular grid points;
   do k = 1, n
      do j = 1, m
         do i = 1, l
            !write(93,*),i,j,k,insas(i,j,k),phi(i,j,k)
            if ( insas(i,j,k) > 0 ) then
               index(i,j,k) = 1
            else
               index(i,j,k) = 5
            end if
         end do
      end do
   end do
   
   do ii = 1, nbnd
      i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)
      !write(86,*) i,j,k,insas(i,j,k)
      if ( phi(i,j,k) > 0 ) then
         index(i,j,k) = 4
      else if ( phi(i,j,k) < 0 ) then
         index(i,j,k) = 2
      else
         index(i,j,k) = 3
      end if
   end do
   !do k =1, n; do j =1, m; do i = 1, l
   !   write(99,*) i,j,k,index(i,j,k)
   !end do;end do;end do
   
   ! setting up index2, the address index of irregular grid points
   !            cirreg, projection points, its transformation matrix & curvature

   call indexg(l,m,n,h,hx,hy,hz,nbnd,x,y,z,index,phi, index2,cirreg)
   !do k =1, n; do j =1, m; do i = 1, l
   !   write(99,*) i,j,k,index2(i,j,k)
   !end do;end do;end do
   !do i = 1, nbnd
   !   write(99,*) i
   !   write(99,*) cirreg(i,1:15)
   !end do
   
   ! setting up the jump conditions at the projection points of the irreular (boundary) grid points

   call prodis(l,m,n,nirreg,bcopt, bi,bo,atmfirst, atmlast, nbnd, cirreg, wp, qp)
   !do i = 1, nbnd
   !   write(98,*) i,wp(i),qp(i)
   !end do
    
   ! calculate the first and second derivatives of the jump conditions in the
   ! surface tangential directions in the local coordinate system
   ! step 1: this is the first derivatives
   call coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
               wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)
   call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,qp, &
             wcoe,wycoe,wzcoe,q0p,qyp,qzp)
   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)
    
   ! step 2: this is the second derivatives
    
   call coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
              wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)
   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,sss1,sss2,wyyp,wzzp,wyzp)
   !write(200,*) wcoe,wycoe,wzcoe,q0p,qyp,qzp,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp
   !write(200,*) q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp
    
   ! setting up linear system coefficient matrix
    
   beta_max = max(bi,bo)
   nz_num = 0
   do k=1,n
      do j=1,m
         do i=1,l
            !write(200,*),i,j,k,index(i,j,k),nz_num
            if ( index(i,j,k).eq.1 ) then
               do ii = 2, 7
                  c(i,j,k,ii) = epsin/h/h
               end do
               if ( i == 1 ) c(i,j,k,2) = 0.d0
               if ( i == l ) c(i,j,k,3) = 0.d0
               if ( j == 1 ) c(i,j,k,4) = 0.d0
               if ( j == m ) c(i,j,k,5) = 0.d0
               if ( k == 1 ) c(i,j,k,6) = 0.d0
               if ( k == n ) c(i,j,k,7) = 0.d0
               c(i,j,k,1) = 6.d0*epsin/h/h
               ii = i+(j-1)*l+(k-1)*l*m
               f(i,j,k)   = bv(i,j,k)/h/h/h*FOURPI
               do ii = 1 , 7
                  if ( abs(c(i,j,k,ii)) > 1.d-10 ) nz_num=nz_num+1
               end do
            else if (index(i,j,k).eq.5 ) then
               do ii = 2, 7
                  c(i,j,k,ii) = epsout/h/h
               end do
               if ( i == 1 ) c(i,j,k,2) = 0.d0
               if ( i == l ) c(i,j,k,3) = 0.d0
               if ( j == 1 ) c(i,j,k,4) = 0.d0
               if ( j == m ) c(i,j,k,5) = 0.d0
               if ( k == 1 ) c(i,j,k,6) = 0.d0
               if ( k == n ) c(i,j,k,7) = 0.d0
               c(i,j,k,1) = 6.d0*epsout/h/h
               ii = i+(j-1)*l+(k-1)*l*m
               f(i,j,k)   = bv(i,j,k)/h/h/h*FOURPI
               !if (i /= 1 .and. i /= l .and. &
               !    j /= 1 .and. j /= m .and. &
               !    k /= 1 .and. k /= n)
               do ii = 1, 7
                  if ( c(i,j,k,ii ) > 1.d-10 ) nz_num=nz_num+1
               end do
!print *,c(i,j,k,1:7)
!print *,"==========",i,j,k,"========="
!print *,c(i,j,k,1:7),f(i,j,k)
!              do nc=1,7
!                 c(i,j,k,nc) = coe1(nc)
!              end do          
!              f(i,j,k) = rhs
!print *,i,j,k, f(i,j,k)
            else
               ii = i+(j-1)*l+(k-1)*l*m
               !write(100,*), i,j,k, sngl(bv(ii)), index(i,j,k)
               !write(200,*) i,j,k,index(i,j,k),'using irre31'
               call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
                           i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                           qyp,qzp,wyp,wzp,wyyp,wzzp,wyzp,  &  
                           nq,nbnd,index2,cirreg,wp,qp,coe2,rhs)
               if (IFAIL.gt.10) then
                  call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
                              i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                              q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp, &
                              nq,nbnd,index2,cirreg,wp,qp,coe2,rhs)
               endif
               ir = index2(i,j,k)
!print *,"=== i ====",i,j,k,"========="
!print *,insas(i,j,k),index(i,j,k)
!print *,x(i),y(j),z(k)
!print *,cirreg(ir,1:15)
!print *,'----------------------------'
!print *,w0p(ir),wyp(ir),wzp(ir)
!print *,wyyp(ir),wzzp(ir),wyzp(ir)
!print *,q0p(ir),qyp(ir),qzp(ir)
!print *,'----------------------------'

               do nc = 1, 27
                  c2(ir,nc) = coe2(nc)
               end do
               f(i,j,k) = rhs
               do ii = 1, 27
                  if ( abs(c2(ir,ii)) > 1.d-10 ) nz_num=nz_num+1
               end do
!print *,c2(ir,1:27),f(i,j,k)
            endif
         end do
      end do
   end do

   ! cleaning up working arrays
   
   deallocate (x, y, z)
   deallocate (phi)
   deallocate (cirreg) ! geometrical parameters at irregular points
   deallocate (wp, qp) ! jump condition at irregular points
   deallocate (q0p,qyp,qzp) ! tangential derivatives of field jump conditions
   deallocate (w0p,wyp,wzp) ! tangential derivatives of potential jump conditions
   deallocate (wyyp,wzzp,wyzp) ! second derivatives of potential jump conditions
   deallocate (wcoe,wxcoe,wycoe) ! coefficients obtained by SVD
   deallocate (wzcoe,wxxcoe,wyycoe)
   deallocate (wzzcoe,wxycoe)
   deallocate (wxzcoe,wyzcoe)
   deallocate (sss1,sss2) !dummy arrays

   ! entering the linear system solver

   if ( solvopt == 1 ) then
      ! algebraic multigrid
      !call amg(l,m,n,nbnd,c,c2,index,index2,f,u,u0,accept) ! should be modified
    write(6,*) "AMG Driver is retired."; call mexit(6,1) 
   else if ( solvopt == 2 ) then
      ! ILU preconditioned GMRES
      call gmres(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,accept)
   else if ( solvopt == 3 ) then
      ! ILU preconditioned BiCG
      call bicg(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,accept)
   else
      write(6,*) "Unsupported Solver option"; call mexit(6,1)
   end if

   ! converting back to the PBSA unit for potential
 

      !write(10,*) "test point potential",u(it,jt,kt),"analytical result",79.0d0/80.0d0/1.2d0
      !write(10,*) "test point1 potential",u(it1,jt,kt),"analytical result",-1.0d0/80.0d0/1.5d0
  
   do k=1, n
      do j=1, m
         do i=1, l
!print *,i,j,k,index(i,j,k)
!          if (index(i,j,k) < 3) u(i,j,k) = u(i,j,k) + fw(x(i),y(j),z(k))
!         if(index(i,j,k)==4   ) then
!           write(8818,*) abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2),&
!            l_green(abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2)),&
!                 -l_green(abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2))/80.0d0/0.5d0,u(i,j,k),&
!               abs(-l_green(abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2))/80.0d0/0.5d0-u(i,j,k))/abs(u(i,j,k))&
!           *100,"%"
!           end if
       u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
!          write (85,*),i,j,k,u(i,j,k) / 4.0 / pi / eps0
!          write (85,*),i,j,k,u(i,j,k)
         end do
      end do
   end do  
   
   deallocate (index, index2)
   deallocate (c, c2)
   deallocate (f)
contains

function l_green(i,j,k)

   implicit none
   integer i,j,k
   _REAL_ l_green
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green

   if ( i <= 20  .and. j <= 20 .and. k <= 20 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green


end subroutine IIM
