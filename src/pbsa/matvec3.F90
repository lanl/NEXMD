#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Caculate A*X with the Un inputted,which is [beta*un] + bf  ]
subroutine matvec3(xm,ym,zm,h,nbnd,nind3,x,y,z,xs,ys,zs,phi,u,f,&
             wp,qp,&
             cirreg,&
             unj,unjf,fvec,bf,epsin,epsout,index,index2,&
             q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp,&
             wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe,&
             wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2,&
             accept,bv,icall,solvopt,norm,inorm )      

  ! input unj, output fvec
  ! nind3 is the  reduced number of  all the irregular points
  ! note unjf and  bf are passed in u&u0 are not
  ! note this is two sides aug approaches.

 use aug_solver   ! change interface names, 'personalized' solver used.
   implicit none 
   integer, parameter :: nq=27
   logical alive
   integer :: status = 0

   ! passed variables

   ! nbnd,nind3              the initial dimension; reduced dimension of aug
   ! x(),y(),z()             coordinates of grid points
   ! xs,ys,zs                staring crds of the box
   ! phi,u                   lvlset function; potential of the whole grids
   ! cirreg                  geometrical parameters at irregular points
   ! index()                 flag of all grid points with 1-5
   ! index2()                number index for irregular grid points
   ! wp,qp                   [u]and[betaUn]for the interface
   ! unj(),unjf()            [Un] for irregular points, unjf is the input
   ! bf(), fvec()            rhs of Ag=b; output of matrix-vector multiply
   ! w0p,wyp,wzp,w0p,        tangential derivatives of jump conditions
   ! wyp,wzp                
   ! wcoe,wxcoe,wycoe,wzcoe  coefficients obtained by SVD
   ! wxxcoe,wyycoe wzzcoe
   ! wxycoe,wxzcoe,wyzcoe
   ! bv                      initial rhs of the linear eqns of the whole grids
   ! solvopt,icall           solvopt for linear solvers above, debuging variable
   ! c,c2                    linear coefficients
   ! accept                  convergence criteria of the linear solvers 
   ! sss1,sss2               dummy arrays
   
   ! local variables

   ! xp,yp,zp                projection crds on interface of the irregular pnts
   ! eps[x-z]                the coefs of grid pnts, which is 1 for AUG method
   ! f(i,j,k)                rhs of the eqns for the whole grid for AUG method
   ! xyy,xzz,xyz             derivatives of the geometry on projection pnts
   ! t                       trans matrix between local crds and grid crds
   ! unin                    derivative of u from inner side on proj pnts
                  
   integer xm,ym,zm,nbnd,nind3,n3,n1,icall
   integer index(1:xm,1:ym,1:zm),index2(1:xm,1:ym,1:zm)
   _REAL_ xyz,xyy,xzz,gox,goy,goz,h,accept,epsout,epsin
   _REAL_  q0p(nbnd),qyp(nbnd),qzp(nbnd),w0p(nbnd),wyp(nbnd),wzp(nbnd),&
           wyyp(nbnd),wzzp(nbnd),wyzp(nbnd)
   _REAL_  wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
           wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
           wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq)
   _REAL_  sss1(nbnd),sss2(nbnd)
   _REAL_  c(xm,ym,zm,7),c2(nbnd, 27)
   _REAL_  wp(nbnd),qp(nbnd)
   _REAL_  t(3,3)
   _REAL_  bv(xm,ym,zm)
   _REAL_  x(0:xm+1),y(0:ym+1),z(0:zm+1)
   _REAL_  phi(0:xm+1,0:ym+1,0:zm+1)
   _REAL_  u(1:xm,1:ym,1:zm),f(xm,ym,zm)
   _REAL_  unj(1:nbnd),cirreg(1:nbnd,1:15),&
           unjf(1:nbnd),fvec(1:nbnd),bf(1:nbnd)

   ! variables to use mg
   _REAL_  inorm,norm,dummy
   integer itn , solvopt
   _REAL_  epsx,epsy,epsz
   _REAL_  iv(1:xm*ym*zm)
   _REAL_  xso(xm*ym*zm+2*xm*ym)

   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo
   _REAL_ uu,dudx,dudy,dudz,unin,unout 
   integer l, m, n, nirreg

   !sgni,j,k is not used for now... XL
   integer :: sgni=0  
   integer :: sgnj=0
   integer :: sgnk=0
   integer i0,j0,k0,nn1

   !WMBS - for use in neutralizing charge on f for periodic solvers
   integer cntirreg !count of the number of irregular grid nodes
   _REAL_ fsum !holds the average charge per irregular grid node to
               !be neutralized
   !!!!!!

   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ beta_max,rhs
   _REAL_  coe2(27)
   _REAL_ fw

   _REAL_ xp,yp,zp
  
   l = xm; m = ym; n = zm; nirreg = nind3
   hx=h;hy=h;hz=h

   epsx=1.0d0;epsy=1.0d0;epsz=1.0d0   
   bi=epsin;bo=epsout
  
   cntirreg = 0

   do i = 1, nind3
      unj(i) = unjf(i)
   end do
 
  !start to check
  !write(1001,*) cirreg
  !write(1001,*) index2
  !write(1001,*) unjf
       
   ! calculate the first and second derivatives of the jump conditions in the
   ! surface tangential directions in the local coordinate system
   ! step 1: this is the first derivatives of w g for  irregular points only 

   call coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
             wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,unj, &
             wcoe,wycoe,wzcoe,q0p,qyp,qzp)
    
   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)

   ! write(200,*) wcoe,wycoe,wzcoe,q0p,qyp,qzp,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp
   ! write(199,*) q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp

   ! step 2: this is the second derivatives ( when secodary order, should be used )

   call coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
              wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,sss1,sss2,wyyp,wzzp,wyzp)

 !  array check
 !  write(71,*) 'qp',unj
 !  write(71,*) 'q0p',q0p
 !  write(71,*) 'qyp',qyp
 !  write(71,*) 'qzp',qzp
 ! do i=1,nind3
 ! write(721,*) 'wp',i,wp(i)
 ! write(722,*) 'w0p',i,w0p(i)
 ! write(723,*) 'wyp',i,wyp(i)
 ! write(724,*) 'wzp',i,wzp(i)
 ! write(725,*) 'wyyp',i,wyyp(i)
 ! write(726,*) 'wzzp',i,wzzp(i)
 ! write(727,*) 'wyzp',i,wyzp(i)
 ! end do 
 ! q0p=qp;qyp=0.0d0;qzp=0.0d0
 ! wp=-0.5d0;w0p=-0.5d0;wyp=0.0d0;
 ! wzp=0.0d0;wyyp=0.0d0;wyzp=0.0d0;wzzp=0.0d0
 ! write(200,*) wcoe,wycoe,wzcoe
 ! write(201,*) q0p,qyp,qzp
 ! write(202,*) wyycoe,wzzcoe,wyzcoe
 ! write(203,*) w0p,wyp,wzp,wyyp,wzzp,wyzp
 ! write(204,*) q0p,qyp,qzp
 ! write(205,*) w0p,wyp,wzp
 ! write(206,*) wyyp,wzzp,wyzp
   ! setting up linear system coefficient matrix
   ! the current version uses irregular points on both sides
 
   beta_max=1.0
   nz_num = 0

   do k = 1, n
      do j= 1, m
         do i= 1, l

            ! inside regular points

            if ( index(i,j,k) == 1 ) then

            ! set up the 7-band laplassian operator in uniform beta_max
    
            ! write(401,*) bv(i,j,k)  check bv

            do ii = 2, 7
               c(i,j,k,ii) =1.0d0/h/h             !epsin/h/h
            end do
            if ( i == 1 ) c(i,j,k,2) = 0.d0
            if ( i == l ) c(i,j,k,3) = 0.d0
            if ( j == 1 ) c(i,j,k,4) = 0.d0
            if ( j == m ) c(i,j,k,5) = 0.d0
            if ( k == 1 ) c(i,j,k,6) = 0.d0
            if ( k == n ) c(i,j,k,7) = 0.d0
            c(i,j,k,1) = 6.0d0/h/h                 !epsin/h/h  XP: this is still for general IIM, epsin /= epsout
            f(i,j,k) = bv(i,j,k)/h/h/h/epsin*FOURPI
            do ii = 1 , 7
               if ( abs(c(i,j,k,ii)) > 1.d-10 ) nz_num = nz_num + 1
            end do
            
            ! outside regular points
 
            else if (index(i,j,k) == 5 ) then

            !write(402,*) bv(i,j,k)   check bv
         
               do ii = 2, 7
                  c(i,j,k,ii) = 1.0d0/h/h          !epsout/h/h
               end do
               if ( i == 1 ) c(i,j,k,2) = 0.d0
               if ( i == l ) c(i,j,k,3) = 0.d0
               if ( j == 1 ) c(i,j,k,4) = 0.d0
               if ( j == m ) c(i,j,k,5) = 0.d0
               if ( k == 1 ) c(i,j,k,6) = 0.d0
               if ( k == n ) c(i,j,k,7) = 0.d0
               c(i,j,k,1) = 6.0d0/h/h              !epsout/h/h ! XP: this is still for general IIM, epsin /= epsout

               f(i,j,k) = bv(i,j,k)/h/h/h/epsout*FOURPI
               do ii = 1, 7
                  if ( c(i,j,k,ii ) > 1.d-10 ) nz_num=nz_num+1
               end do
             
              !write(61,*) i,j,k,f(i,j,k)  check rhs
   !print *,c(i,j,k,1:7)
   !print *,"==========",i,j,k,"========="
   !print *,c(i,j,k,1:7),f(i,j,k)
   !              do nc=1,7
   !                 c(i,j,k,nc) = coe1(nc)
   !              end do          
   !              f(i,j,k) = rhs
   !print *,i,j,k, f(i,j,k)
 
   !   write(200,"(7f20.6)") c(i,j,k,1:7)

            ! irregular points

            else 

            cntirreg = cntirreg+1  !WS
     ! interface variables check
     !write(1001,*) i,j,k,x,y,z,phi,index
     !write(1002,*) q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp
     !write(1003,*) nq,nbnd
     !write(1004,*) index2
     !write(1005,*) cirreg
     !write(1006,*) wp,unj   
     !XP: set IIM as 80.0 same result from updated: set 1 80
     !write(1007,"(7f20.6)") coe2(1:7)/beta_max
        
       coe2=0.0d0 

     ! XP: Note the interface has changed. b_in, b_out, and bi bo are all passed in.

     bo=bi 
             call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
                           i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                           qyp,qzp,wyp,wzp,wyyp,wzzp,wyzp,  &
                           nq,nbnd,index2,cirreg,wp,unj,coe2,rhs)

               if (IFAIL.gt.10) call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
                            i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                            q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp, &
                            nq,nbnd,index2,cirreg,wp,unj,coe2,rhs)
               ir = index2(i,j,k)
     bo=epsout
   !  interface variables check
   !write(30,*) i,j,k,rhs
   !print *,"=== i ====",i,j,k,"========="
   !print *,x(i),y(j),z(k)
   !print *,cirreg(ir,1:15)
   !print *,'----------------------------'
   !print *,w0p(ir),wyp(ir),wzp(ir)
   !print *,wyyp(ir),wzzp(ir),wyzp(ir)
   !print *,q0p(ir),qyp(ir),qzp(ir)
   !print *,'----------------------------'

                     
               f(i,j,k) = -rhs            !*(-6.0d0)/coe2(14)
        
   !     Out of the 27 neighbors, 7 is nonzero and contribute
   !     coe2(14)=-6.0d0*4.0d0
   !     coe2(5)=1.0d0*4.0d0
   !     coe2(11)=1.0d0*4.0d0
   !     coe2(13)=1.0d0*4.0d0
   !     coe2(15)=1.0d0*4.0d0
   !     coe2(17)=1.0d0*4.0d0
   !     coe2(23)=1.0d0*4.0d0

 
               do nc = 1, 27
                  c2(ir,nc) = coe2(nc)
               end do
               do ii = 1, 27
                  if ( abs(c2(ir,ii)) > 1.d-10 ) nz_num = nz_num + 1
               end do
      !  variable check
      !         write(201,"(7f20.6)") c2(ir,1:7)!/beta_max
      !         write(202,"(20f20.6)") c2(ir,8:27)!/beta_max

            end if
         end do
      end do
   end do
    !Outputs check for the setup of the linear eqns:
    !  entering the linear system solver
    !  if(icall == 2 ) then 
    !  write(7001,*) l,m,n,nbnd,nz_num 
    !  write(7002,"(7f20.6)") c 
    !  write(7007,"(7f20.6)") c2              
    !  write(7003,*) index,index2    
    !  write(7004,*) f               !error
    !  write(7006,"(7f20.6)") u               
    !  write(7005,*) u0              
    !  write(7005,*) accept,icall          
    !  write(7008,*) "h=",h,"epsin=",epsin,"epsout=",epsout
    !  end if

    !  if(icall == 1 ) then 
    !  write(6009,*) bv
      !write(6001,*) l,m,n,nbnd,nz_num 
      !write(6002,"(7f20.6)") c 
      !write(6007,"(7f20.6)") c2              
      !write(6003,*) index,index2    
      !write(6004,*) f               !error
      !write(6006,"(7f20.6)") u               
      !write(6005,*) u0              
      !write(6005,*) accept          
      !write(6008,*) "h=",h,"epsin=",epsin,"epsout=",epsout
    !  end if
 
   !bicg solver which can solve linear eqns when epsin/=epsout

!  write(10,*) "here we are before bicg"
!      call bicg(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,accept)! caculate u
!             write (88,"(f20.10)") u
 
     ! Use multigrid to solve it efficiently!!
          
             f=f*h*h         

             dummy=0.0d0;iv=0.0d0; !use iccg

             !WMBS - Need to ensure neutral system if using periodic cg as
             !solver. Any other periodic solvers will also need to do the same
             !since electrical potential is ill-defined for a periodic monopole
             !system.
           !if ( solvopt==5) then
             !   if (abs(sum(f)/(xm*ym*zm)) > 1.0d-8) then 
             !           write(6,*) "Non-neutral total system charge of: "
             !           write(6,*) CHAR(09),sum(f)
             !           write(6,*) "Adding charge of ",sum(f)/(xm*ym*zm)," to all"&
             !                   ,xm*ym*zm," Nodes"
             !           f = f - sum(f)/(xm*ym*zm)
             !           write(6,*) "New system net charge is: ",sum(f)
             !           flush(6) 
            !end if

           !  u = u*h
             !write(6,*) 'writting u for fft'; flush(6)
             !write(17,*) 'fft'
            !do i=1,xm; do j=1,ym; do k=1,zm
            !    write(17,*) i,j,k,u(i,j,k)
            ! end do; end do; end do
            ! flush(17)
          !end if
             !call init_param(l,m,n,l*m,l*m*n,10000,1.95d0,accept,0.0d0,1.0d0,h,dummy)
 
             call init_param(l,m,n,l*m,l*m*n,10000,dummy,accept,&
                  0.0d0,1.0d0,h,1.9d0)
             ! call init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)
 
             call allocate_array(solvopt) 
             
             xso(:) = 0.d0
             call init_array(solvopt,epsx,epsy,epsz,f,iv,xso)
              
          if(solvopt/=7 ) then !.and. solvopt/=5) then
           call pb_iccg(u,xso)
             !write(6,*) 'writting u for iccg'; flush(6)
             !write(16,*) 'iccg'
             !do i=1,xm; do j=1,ym; do k=1,zm
             !   write(16,*) i,j,k,u(i,j,k)
             !end do; end do; end do
             !flush(16)
            !call pb_sor(u,xso) 
             !write (88,"(f20.10)") u
          !else if ( solvopt==2) then
          !   call pb_mg(u,xso)
        ! else if ( solvopt==5) then
          !WMBS - Caveat emptor!
          !     This solver is slow to begin with.
          !     Convergence for gmres may become excruciatingly slow.
        !    call pb_pcg(u,xso)
             !u = u*h
             !write(6,*) 'writting u for fft'; flush(6)
             !write(17,*) 'fft'
            !do i=1,xm; do j=1,ym; do k=1,zm
            !    write(17,*) i,j,k,u(i,j,k)
            ! end do; end do; end do
            ! flush(17)
          else if ( solvopt==7) then
             !WMBS - Need to neutralize any excess charge before calling fft.
             !FFT will do this internally, but it won't show up in f since it
             !will switch to an internal charge grid.
             if (abs(sum(f)/(xm*ym*zm)) > 1.0d-8) then 
               ! write(6,*) "Non-neutral total system charge of: "
               ! write(6,*) CHAR(09),sum(f)
               ! write(6,*) "Adding charge of ",sum(f)/(xm*ym*zm)," to all"&
               !                 ,xm*ym*zm," Nodes"
                fsum = sum(f)/(xm*ym*zm)
                do i=1,l; do j=1,m; do k=1,n
                    !if (.not. (index(i,j,k) == 1 .or. index(i,j,k) == 5) ) then
                        f(i,j,k) = f(i,j,k) - fsum
                    !end if    
                end do; end do; end do;
               ! write(6,*) "New system net charge is: ",sum(f)
                flush(6) 
             end if
             call pb_fftsolv(f,2,xm,ym,zm,u,h,64)
             u = u*h
             !write(6,*) 'writting u for fft'; flush(6)
             !write(17,*) 'fft'
            !do i=1,xm; do j=1,ym; do k=1,zm
            !    write(17,*) i,j,k,u(i,j,k)
            ! end do; end do; end do
            ! flush(17)
          end if
               itn = l_itn
              inorm = l_inorm                                                             
              norm = l_norm  
           
              !itn = l_itn
              !inorm = l_inorm
              !norm = l_norm
           
              call deallocate_array(solvopt)

!       do k=1,n
!        do j=1,m
!          do i=1,l
!          if(index(i,j,k)<=3)  then
!            write(288,*) u(i,j,k)
!           end if
!           u(i,j,k)=0.49375d0
!             end do 
!           end do 
!        end do 
   


    !Below is the green setup of u
    
    !     write(187,*) "nbnd=",nbnd
    !
    !     do k=1,n
    !       do j=1,m
    !         do i=1,l
    !        if(index(i,j,k)<=3) then
    !          u(i,j,k)= l_green(abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2))/0.5d0&
    !          +0.5d0/80.d0-0.5d0
    !        else
    !          u(i,j,k)= l_green(abs(i-(l+1)/2),abs(j-(m+1)/2),abs(k-(n+1)/2))/0.5d0/bo
    !          
    !        end if
    !            end do 
    !          end do 
    !       end do 



     ! Below is the read in of u
     
     !   inquire(file=filename,exist=alive)
     !   if(alive) then
     !   write(6,*) filename,"exist"
     !   end if
     !
     !   if(alive) then
     !   open(unit=10,file=filename,access="sequential",status="old")
     !   ii=0
     !   do while(.true.)
     !      ii=ii+1
     !      read(unit=10,fmt=*,iostat=status) u
     !   if(status/=0) exit
     !      do k=1,n
     !       do j=1,m
     !         do i=1,l
     !   write(818,*) u(i,j,k)
     !            end do 
     !          end do 
     !       end do 
     !
     !  end do
     !  else
     !   write(*,*) TRIM(filename)," doesn't exit!"
     !  end if

     !Below is the analytically setup of u

     ! write(12,*) 'b interp'
     !do k=1,n
     !   do j=1,m
     !      do i=1,l
     !      if(index(i,j,k) <= 3) then
     !      write(88,*) i,j,k,u(i,j,k)
     !   end if
     !   if(index(i,j,k) <= 3) then
     !    u(i,j,k)=1.0d0/2-1.0d0/2/80
     !    else 
     !    u(i,j,k)=-1.0d0/sqrt((i-25.0d0)*(i-25.0d0)+&
     !    (j-25.0d0)*(j-25.0d0)+(k-25.0d0)*(k-25.0d0))/0.25/80.0d0
     !    end if
     !    write(84,*) i,j,k,u(i,j,k)
     !    end do 
     !  end do
     !end do
  
   n3=10! second order interplation
   do k=1,n
      do j=1,m
         do i=1,l
            nn1=index2(i,j,k)

            if (nn1 .gt. 0) then

            !projection crds and geometrical info 
   
               xp=cirreg(nn1,1)
               yp=cirreg(nn1,2)
               zp=cirreg(nn1,3)

               xyy = cirreg(nn1, 4)
               xzz = cirreg(nn1, 5)
               xyz = cirreg(nn1, 6)
    
               t(1,1) = cirreg(nn1, 7)
               t(1,2) = cirreg(nn1, 8)
               t(1,3) = cirreg(nn1, 9)
               t(2,1) = cirreg(nn1, 10)
               t(2,2) = cirreg(nn1, 11)
               t(2,3) = cirreg(nn1, 12)
               t(3,1) = cirreg(nn1, 13)
               t(3,2) = cirreg(nn1, 14)
               t(3,3) = cirreg(nn1, 15)

               !the nearest grid point

               i0 = nint((xp - x(0))/h)
               j0 = nint((yp - y(0))/h)
               k0 = nint((zp - z(0))/h)
    bo=bi 
               call interpx(l,m,n,bi,bo,n3,i0,j0,k0,sgni,sgnj,sgnk,&
                       nbnd,xp,yp,zp,uu,dudx,dudy,dudz,&
                       wp,wyp,wzp,wyyp,wzzp,wyzp,unj,qyp,qzp, &
                       xyy,xzz,xyz,t,u,phi,x(0),y(0),z(0),h,&
                       index2,nn1)
    bo=epsout

           ! call interp output Un^-, dudx in the local crds
                              
               unin = dudx 
               unout = unin + unj(nn1)

            !  Preconditioning, which results in smaller residual
            !  unout = (qp(nn1)-unj(nn1))/79.0d0 , as we are using 1:80

               fvec(nn1)=epsout*unout-epsin*unin-qp(nn1)
               fvec(nn1)=fvec(nn1)+bf(nn1)
            end if
         end do
      end do 
   end do
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
end subroutine matvec3
