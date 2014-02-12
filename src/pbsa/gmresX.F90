#define _REAL_ double precision

! The following process is the GMRES(mm) iterative method for  the algebraic 
! system  A*x = b  with A being any non-symmetric matrix. The original code 
! is under FIIM LEVEL from Zhilin and is transformed form 2D to 3D now. nbnd 
! is number of all irregular points,while ctn is the reduced one after clustering now.

! The algorithm for GMRES is described as below:
! To solve A*X=f; x0 is the guess vector for solution r0=f-Ax0 and v1=r0/||r0||
! The key aim is to use Schmidt orthogonalization method to orthogonalize the
! Krylov subspace {v1,Av1,A^2V1,...,A^(k-1)v1}
! Iterate: for j=1,2,...,m  ;
! Vj+1=Avj-sum((Avj*vi)*vi,i=1,2,..,j, vi is normalized.
! vj+1=Vj+1/||Vj+1||  Got the nomalized and orthogonalized basis.
! solution x=x0+Vm*ym, please refer to related material to know the math of
! implicit residual.
! And A*x is done my subroutine matvec.

! Implicit residual, methods to calculat the residual without needing to
! calculate ym at the time:
! AVk=Vk+1 H'  res=||f-A(x0+vmym||=||r0-Vk+1 H'ym||
! => res=||Vk+1(beta*e1-H'ym)||=||beta*e1-H'ym|| 
! Please be advised that only the last element of the last row of H' is nonzero,
! which makes that when Q applys, res=||Q(beta*e1-H'ym)||=||gk-R*ym|| last row
! of R is zero, and the minimum norm would be the others are zero, thus the last
! element of gk is the residual, the implicit one.

subroutine gmresx(mm,l,m,n,nbnd,ctn,imax,h,phi,x,y,z,xs,ys,zs,&
                   index,index2,& 
                   wp,qp,& 
                   cirreg,&
                   unj,u,f,x0,bf,&
                   tol,iter,error,bout,bin,&
                   q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
                   wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
                   wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
                   accept,bv,solvopt,pbverbose,norm,inorm)          
   
   !  Passed variables

   ! Input:
   !  pbverbose              logical flag variable, .true. to print detail
   !  mm	             INTEGER number of iterations between restarts
   !  l,m,n                  dimension of the 3-D box
   !  [x-z]s,x-z             starting crds and 3-D crds array of the box
   !  nbnd,ctn               original dimension of irregular ptn; clustered one
   !  imax                   max times of restarts of the GMRES algorithm
   !  h                      grid spacing
   !  phi                    level set function   
   !  index                  the index labeling for all grid points
   !  index2                 index of number for irregular points
   !  cirreg                 geometrical parameters at irregular points
   !  unj                    jump of 1st derivative of potential on interface
   !  wp,qp                  [u]and[betaUn]for the interface
   !  x0	             REAL initial guess vector
   !  bf f                   REAL right hand side vector & variable to store
   !  tol                    determine resolution of the solution
   !  bout,bin               dielectric coefficients of outside and inside
   !  w0p,wyp,wzp,w0p,       tangential derivatives of jump conditions
   !  wyp,wzp                
   !  wcoe,wxcoe,wycoe,wzcoe coefficients obtained by SVD
   !  wxxcoe,wyycoe wzzcoe
   !  wxycoe,wxzcoe,wyzcoe
   !  accept                 convergence criteria of the linear solvers 
   !  bv                     initial rhs of the linear eqns of the whole grids
   !  solvopt                solvopt for linear solver of whole grid potential

   ! Output:
   !  u                      potential of whole grids
   !  x0                     The solution vector
   !  iter          	     The number of iterations
   !  error	             The 2-norm of the residue of the solution
   
   ! Local variables 
   
   !  hg                     stores the Hessenberg matrix.
   !  v                      stores the modified Gram-Schmidt orthonormal
   !		             which forms the Krylov subspace.
   !  ser1,ser2              implicit residual; explicit residual
   !  RSD,MAXD,ABSD          all kinds of criteria for convergence,debuging use 

   ! n2<-nind3  index<-index2    index2<-index3

   implicit none 
   integer, parameter :: nq=27 

    !passed varibles for IIM 

    logical pbverbose
    _REAL_  accept,tol,bin,bout,h, ssum,nsum,xs,ys,zs
    integer mm,m,n,l,n2,nbnd,imax,solvopt,ctn
    _REAL_  q0p(nbnd),qyp(nbnd),qzp(nbnd),w0p(nbnd),wyp(nbnd),wzp(nbnd),&
            wyyp(nbnd),wzzp(nbnd),wyzp(nbnd)
    _REAL_  wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
            wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
            wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq)
   _REAL_   sss1(nbnd),sss2(nbnd)
   _REAL_   c(l,m,n,7),c2(nbnd, 27)
   _REAL_   x0(nbnd),bf(nbnd),r(nbnd),s(nbnd), &
            x11(nbnd),vk(nbnd),w(nbnd),hj(mm,2),yy(nbnd),&
            x110(nbnd),w0(nbnd),yy0(nbnd)
   _REAL_   x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_   f(l,m,n),u(1:l,1:m,1:n),phi(0:l+1,0:m+1,0:n+1),bv(l,m,n)
   _REAL_   unj(nbnd)
   _REAL_   cirreg(1:nbnd,1:15)
   _REAL_   wp(1:nbnd),qp(1:nbnd)
   integer  index(1:l,1:m,1:n),index2(1:l,1:m,1:n)
   _REAL_   norm, inorm

  !local variables

   integer idimf,infon,i1,j2,i2,k,k1,k2,mi,n1,j,i,iter,j1,icall
   _REAL_  drl,ser0,a,b,d,dr1,dsqrt,&
           hit,dt,st,ser1,serr,error,ser2,RSD0,MAXD0,ABSD0,RSD,MAXD,ABSD
   _REAL_,allocatable :: hg(:,:),v(:,:)
   
  allocate (hg(nbnd,mm),v(nbnd,mm))
  write(6,*) 'The gmres dimension is ',ctn
  iter = 0

  ! Loop 2, this is the restart of the algorthm with a guess vector x0

   do j=1,imax

  ! Initializaiton for variables used each start of the algorithm

   ser0 = 0.0;hg=0.0d0;v=0.0d0;w0=0.0d0

   icall = 20  ! Flag for debuging

   
  ! first time to call matvec to perform r0=Ax0-b, get residual r0 and calculate v1

     call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
              wp,qp,&
              cirreg,&
              unj,x0,yy0,bf,bin,bout,index,index2,&
              q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
              wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
              wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
              accept,bv,icall,solvopt,norm,inorm )  
 
     do i=1,ctn       
      yy(i)=yy0(i)
     end do
      
   call residX(ctn,yy,bf,r)       !chng resid to residX  note we use 2 sides

 ! Initial values of RSD MAXD and ABSD values, used for comparative convergence
 ! criteria

 ! if(j==1) then 
 ! RSD0=RSD(ctn,r);MAXD0=MAXD(ctn,r);ABSD0=ABSD(ctn,r)    
 ! write(6,*) 'Initial ones:','RSD0=',RSD0,'MAXD0',MAXD0,'ABSD0',ABSD0
 ! end if

   dr1 = dsqrt(dot_product(r(1:ctn),r(1:ctn)))

   ! if (pbverbose) write(6,*) 'norm of v1',dr1

   ! if x0 satisfy the eqn..  already get the desired potential u, return
  !write(10,*) 'norm of v1',dr1
   if(dr1 <= 1e-12) return

   ! Get v1
   do i=1,ctn
       v(i,1) = r(i)/dr1
       s(i) = 0.0d0
   end do

   s(1) = dr1

      !--Loop 1-- For i=1,2, ..., mm -----calculate w v and matrix h-------------------------

    do i=1,mm-1
     iter = iter + 1

   do i1=1,ctn
      vk(i1) = v(i1,i)
   end do

   ! do Avi, calculate the next basis vector of the Krylov subspace

   call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
                wp,qp,&
                cirreg,&
                unj,vk,w0,bf,bin,bout,index,index2,&            
                q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
                wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
                wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
                accept,bv,icall,solvopt,norm,inorm)             !input vk,output w
   do k=1,ctn
      w(k)=w0(k)
   end do 
   
   ! Calculate the projs of the new  basis in direction of former basises

   do k=1,i

      do i1=1,ctn
         vk(i1) = v(i1,k)
      end do

      hg(k,i) = dot_product(w(1:ctn),vk(1:ctn)) ! basis in vk direction
           
   !orthogonalizing..

      do i1=1,ctn
         w(i1) = w(i1) - hg(k,i)*vk(i1)
      end do

   end do

   ! norm of the next orthogonalized basis, the i+1 th one

          hg(i+1,i) = dsqrt(dot_product(w(1:ctn),w(1:ctn)))

   ! If hg == 0, then we got the complete basis of the space.

   if( abs(hg(i+1,i)) <= 1e-14) then
         infon = 1
    write(6,*) 'info=1, Got complete basis of space.'
   else
    infon = 0
  
   ! the nomalized new vector of the i+1 th.

    do i1=1,ctn
       v(i1,i+1) = w(i1)/hg(i+1,i)
    end do

   end if
    !----------------------------------------------------------------
    ! Convert matrix h to get the residual implicitly
    !
    !----- Apply J_1, j_2, ..., J_{i-1} on (hg_{1,i}, ..., hg_{i+1,i} ---------
    !   Suppose J_i =
    !                  | I 0 |    p = | cos(\alf)  sin(\alf) |
    !                  | 0 P |            | -sin(\alf) cos(\alf) |

    !       cos(\alf) = hg(i,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)
    !       sin(\alf) = hg(i+1,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)

    !------- Form J_i so that the (i+1)th component of J_i*h(:,i) is zero.

       !write(*,*) '##########below to see whether AV1=V2H############'

      ! do i1=1,ctn
      !  w0(i1)=w0(i1)-hg(1,1)*v(i1,1)-hg(2,1)*v(i1,2)
      ! end do
      !write(*,*) 'This should be zero',dsqrt(dot_product(w0(1:ctn),w0(1:ctn)))
         
      !write(*,*) ' ################# test results for iter=1##########'
      !write(*,*) 'norm v1 is',dsqrt(dot_product(v(1:ctn,1),v(1:ctn,1)))
      !write(*,*) 'norm v2 is',dsqrt(dot_product(v(1:ctn,2),v(1:ctn,2)))
      ! h11=hg(1,1);h21=hg(2,1)
      !write(*,*) 'origial hg(1,1)',hg(1,1),'hg(2,1)',hg(2,1)

    do k=1,i-1
       hit = hj(k,1)* hg(k,i) + hj(k,2)*hg(k+1,i)
       hg(k+1,i) = -hj(k,2)*hg(k,i) + hj(k,1)*hg(k+1,i)
       hg(k,i) = hit
    end do
      !s11=s(1); s22=s(2)
      !write(*,*) 'original s(1)',s(1),'s(2)',s(2),'s(1) should be v1,s(2)0'
    if(infon == 0) then
        dt = dsqrt(hg(i,i)*hg(i,i)+hg(i+1,i)*hg(i+1,i))
        hj(i,1) = hg(i,i)/dt
        hj(i,2) = hg(i+1,i)/dt
        st = hj(i,1)*s(i) + hj(i,2)*s(i+1)
        s(i+1) = -hj(i,2)*s(i) + hj(i,1)*s(i+1)
        s(i) = st
        hit = hj(i,1)* hg(i,i) + hj(i,2)*hg(i+1,i)
        hg(i+1,i) = -hj(i,2)*hg(i,i) + hj(i,1)*hg(i+1,i) ! found by XP
        hg(i,i) = hit

    end if
      ! Debugging 
     
      !write(*,*) 'cos',hj(1,1),'sin',hj(1,2),'norm',hj(1,1)**2+hj(1,2)**2
      !write(*,*) 'new hg(1,1)',hg(1,1),'hg(2,1)',hg(2,1)
      !write(*,*) 'new s(1)',s(1),'s(2)',s(2)
      !write(*,*) 'so ser1=',s(2)

         ser1 = abs(s(i+1))

      !begine explicit test by XP

       mi=i
       yy(mi) = s(mi)/(hg(mi,mi))!+1.0d-14)

      do k=mi-1,1,-1
         yy(k) = s(k)
         do j1 = k+1,mi
            yy(k) = yy(k) - hg(k,j1)*yy(j1)
         end do
         yy(k) = yy(k)/hg(k,k)   ! The coefficients for each basises.
      end do

      do i2=1,ctn
         x11(i2) = x0(i2)

         do k=1,mi
            x11(i2) = x11(i2) + yy(k)*v(i2,k)
         end do
      end do

     do i2=1,ctn
         x110(i2)=x11(i2)
     end do 
 
    ! calculate the real residual, explicit residual

    call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
               wp,qp,&
               cirreg,&
               unj,x110,w,bf,bin,bout,index,index2,&    
               q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
               wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
               wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
               accept,bv,icall,solvopt,norm,inorm)       !intput x110 out put w
        do i2=1,ctn
         x11(i2)=x110(i2)    
        end do           
   
    call residX(ctn,w,bf,r)
       ser2 = dsqrt(dot_product(r(1:ctn),r(1:ctn)))

         
   !Judge residual if OK goto calculate it explicitly if not go on loop1

   !XP: Here we use the explicit residual as the criteria
        !if (ser1 > 1.0d-15) then
        !   serr = abs((ser1-ser0)/ser1)
        !else
        !   serr = 0.0
        !end if

  !write(10,*) "ser1=",ser1,"ser2=",ser2,"tol=",tol,"iter=",iter;flush(10)
     !write(6,*) "i,j",i,j
     !write(6,*) "ser1=",ser1,"ser2=",ser2,"tol=",tol,"iter=",iter
     !write(6,*) "relative:","RSD",RSD(ctn,r)/RSD0,"MAXD",MAXD(ctn,r)/MAXD0,"ABSD",ABSD(ctn,r)/ABSD0
     !write(6,*) "absolute:","RSD",RSD(ctn,r),"MAXD",MAXD(ctn,r),"ABSD",ABSD(ctn,r)
  
     !Choices for criteria for convergence.

      !if(ser1 <= tol .or. serr<=tol ) then  Notes: serr= ser1/(last ser1),
      if(ser2 <= tol  ) then
      !if(RSD(ctn,r)/RSD0 <= tol  ) then
      !if( ser1< tol ) then
      !write(10,*) 'Residual=',ser2,'ITER NUMBER=',iter,'j',j
      !write(6,*) 'Residual=',ser2,'ITER NUMBER=',iter,'j',j
            serr = tol - 1.0e-15
            mi = i
            goto 100  
        end if
       
      ! restart if the explicit and implicit residuals are far different

        if(ser1/ser2 > 10.0d0 .or. ser1/ser2 < 1.0e-1) then
          mi=i
        !write(10,*) ' ex&im '
      !write(6,*) ' ex&im ','j',j
         goto 100
       end if
        ser0 = ser1
      end do  
      mi = mm - 1

      !end loop 1 for i=1,2,...mm-----------------------------------

      ! Update(x_i,i) Calculate y x and call matvec to calculate the residual

      100 yy(mi) = s(mi)/(hg(mi,mi)+1.0d-14)
      !100    yy(mi) = s(mi)/hg(mi,mi)

      do k=mi-1,1,-1
         yy(k) = s(k)
         do j1 = k+1,mi
            yy(k) = yy(k) - hg(k,j1)*yy(j1)
         end do
         yy(k) = yy(k)/hg(k,k)
      end do

      !do i=1,nbnd
      do i=1,ctn
         x11(i) = x0(i)
         do k=1,mi
            x11(i) = x11(i) + yy(k)*v(i,k)
         end do
      end do
 
  !   write(*,*) 'XP: Following are used to test AVk=Vk+1Hk'
  !  
  !   write(*,*) ' XP: Following are tests for solution:'
 
  !stest=0.0d0     
  !   do i2=1,mi
  !    stest(i2)=s(i2)              
  !     do j2=1,mi
  !    stest(i2)=stest(i2)-hg(i2,j2)*yy(j2)              
  !     end do
  !  write(*,*) 'stest',i2,stest(i2)
  !   end do 
 
     do i=1,ctn
         x110(i)=x11(i)
     end do 
    icall = 100
    call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
               wp,qp,&
               cirreg,&
               unj,x110,w,bf,bin,bout,index,index2,&    
               q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
               wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
               wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
               accept,bv,icall,solvopt,norm,inorm)       !intput x110 out put w

        do i=1,ctn
         x11(i)=x110(i)    
       
        end do           
    !explicit residual

    call residX(ctn,w,bf,r)
    s(mi+1) = dsqrt(dot_product(r(1:ctn),r(1:ctn)))
 
    !update x0 for potential restart.
    do k=1,ctn
       x0(k) = x11(k)
    end do

   !-----Judge if residual OK reture, or back to loop 2 j=1 imax------------
 
    !  write(6,*) 'abs(s)=',abs(s(mi+1)),'tol=',tol,'j',j
    !  write(6,*) "relative:","RSD",RSD(ctn,r)/RSD0,"MAXD",MAXD(ctn,r)/MAXD0,"ABSD",ABSD(ctn,r)/ABSD0
    !  write(6,*) "absolute:","RSD",RSD(ctn,r),"MAXD",MAXD(ctn,r),"ABSD",ABSD(ctn,r)
    !  !!write(10,*) abs(s(mi+1)),iter

      if( abs(s(mi+1)) < tol ) then
      !if( RSD(ctn,r)/RSD0 < tol ) then
      !write(6,*) 'return'
       deallocate(hg,v)
     return
     end if
     
   !avoid dead lock, when convergence can't be minimized, just restart.

   if (j > 4) then
   ! write(6,*) 'Possible dead lock, possible failure'
   ! write(6,*) 'j=',j,'res=',s(mi+1)
     deallocate(hg,v)
      return
     end if

      error = s(mi+1)

   end do  ! j=1,imax
   !write(6,*) 'j=',j,'imax=',imax
   deallocate(hg,v)
   return
end subroutine gmresx

!*******************************************************************************


subroutine residX(n,x,b,r)


   !+ [The process is to perform y=Ax.]

   implicit none 

   integer i,j,n
   _REAL_  x(n), b(n), r(n)

   do i=1,n
      r(i) = b(i) - x(i)
   end do
   return
end subroutine residX

function RSD(n,x)

   implicit none
 
   integer n,i
   _REAL_ x(n),RSD
   
   RSD=0.0d0
   do i=1,n
   RSD=RSD+x(i)*x(i)
   end do
   RSD=sqrt(RSD)
end function RSD
  
function MAXD(n,x)
 
  implicit none
  integer n,i
  _REAL_  x(n),MAXD
  MAXD=0.0d0
  do i=1,n
  if( MAXD < abs(x(i))) MAXD=x(i)
  end do
end function MAXD

function ABSD(n,x)

  implicit none
  integer n,i
  _REAL_  x(n),ABSD

  ABSD=0.0d0
  do i=1,n
  ABSD=ABSD+abs(x(i))
  end do
end function ABSD 

