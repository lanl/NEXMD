#include "dprec.fh"
!
!********************************************************************
!
!  Davidson diagonalization
!
!  Converted to fortran 90 (from fortran 77)
!  and interfaced with sqm(AmberTools) by
!  Kirill A. Velizhanin (kirill@lanl.gov)
!
!********************************************************************	
!
!********************************************************************
!   
   subroutine davidson(qmmm_struct)
   use qmmm_module,only:qm2_struct !cml-test
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type


   implicit none

     type(qmmm_struct_type), intent(inout) :: qmmm_struct

   _REAL_ random,rranset,rranf
   _REAL_ ferr1,ferr2,fn
   _REAL_ f,f1,f2,ddot,f0,f11

   integer i,j,lprint
   integer iseed,ip,ih,j0,Mx0,imin,one,iloops
   integer kflag,nd,nd1,nd1_old,j1

   integer,save::istore=0 ! zero initially
   integer,save::istore_M=0

   integer Lt,Mb
!
!--------------------------------------------------------------------
!
!  Initialize Davidson
!
!--------------------------------------------------------------------
!
   Lt=qm2ds%Nb*(qm2ds%Nb+1)/2
   Mb=qm2ds%Nb**2
   one=1 
   lprint=qm2ds%verbosity       

   nd=qm2ds%nd 

   iloops=0

   if(qm2ds%dav_guess==0) then
      istore=0 ! overwriting istore to not use guess
   end if

! Split Mx into batches of j1 size
   j1=nd/qm2ds%idavfrac
   j1=min(j1,qm2ds%Mx)
   nd1=min(j1+2,nd/2)
   if (qm2ds%mdflag.ge.0) j1=qm2ds%Mx ! No shift is allowed for MD points
      
!	if (irflag.gt.1.and.mdflag.lt.0) then
   !qm2ds%irflag=5
   !write(6,*)"qm2ds%irflag=",qm2ds%irflag

   if (qm2ds%irflag.gt.1) then
!     Calculations will start from irflag state for ceo
      j0=qm2ds%irflag-1
      !write(6,*)"begin at ", j0,qm2ds%Mx,qm2ds%irflag
      if (j0.gt.qm2ds%Mx) then
         write(6,*) 'Looks like irflag= ', qm2ds%irflag
         write(6,*) 'Show that states Mx= ',qm2ds%Mx
         write(6,*) 'are already calculated, exiting'

         stop 
      end if

!     Read irflag-1 states 
      open (10,file='modes.b',form='unformatted',status='old')
      open (12,file='ee.b',status='old')

      do j=1,j0
         read (10) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
         read (12,*) qm2ds%e0(j)
      enddo
   
      close(10)
      close(12)

      if (j0.eq.qm2ds%Mx) goto 70

   else    ! Assume start from the beginning 
      j0=0
   end if

   qm2ds%Mj=j0

   if (lprint.gt.1) then
      write(6,*) 'Davidson parameters'
      write(6,*) 'mdflag=',qm2ds%mdflag
      write(6,*) 'irflag=',qm2ds%irflag
      write(6,*) 'M4=',qm2ds%Ncis
      write(6,*) 'Mx=',qm2ds%Mx
      write(6,*) 'Mj=',qm2ds%Mj
      write(6,*) 'nd=',nd 
      write(6,*) 'nd1=',nd1
      write(6,*) 'j1=',j1
      write(6,*)'istore=',istore
      write(6,*)'istore_M=',istore_M
   end if

!	if (irflag.gt.1.and.Nb.gt.100.and.mdflag.gt.0) then
!	   istore=min(Mx,Mx_ss)
! --- read stored modes from the previous calculations in MD run
!         open (10,file='modes-ao.b',form='unformatted',status='old')       
!         do j=1,istore
!	    read (10) (v2(i,j),i=1,Nb*Nb) 
!	   enddo
!	   close(10)	  
!	endif	 

   if(istore.gt.0) then ! MD point only!!!!	 
!     recover excited state vectors from AO representation in v2
!     recovered state vectors are put to v0
      do j=1,istore
         call site2mo(qm2ds%rrwork,qm2ds%v2(1,j),qm2ds%v0(1,j))
      end do
   endif

! kav: the lines below are commented out since they are not present in original davidson
! 
! CML TEST 7/15/12
   !if(qm2ds%mdflag /= 0) istore = qm2ds%Mx ! If first run, don't try and recover vectors. Do it on future calls.
   !if (qm2ds%mdflag /= 0) istore = 0 ! If first run, don't try and recover vectors. Do it on future calls. ! KGB
!
!--------------------------------------------------------------------
!  Begin big loop
!--------------------------------------------------------------------
! find many vectors:
10 continue
   if (j0.lt.qm2ds%Mx) then 

   iloops=iloops+1
!	if (iloops.eq.2) j1=j1/2 ! After the first one convergence is slow
   j1=min(j1,(qm2ds%Mx-j0))
   nd1=min(j1+2,nd/2)

   if (lprint.gt.3) then
      write(6,*)
      write(6,*) 'Davidson batch ', iloops
      write(6,*) 'So far found ', j0, ' states'
      write(6,*) 'out of requested ',qm2ds%Mx, ' states'
      write(6,*) 'This batch will seek',j1,' vectors'

      if(j0.gt.0) write(6,*) 'Shift is',qm2ds%fs+qm2ds%e0(j0),' eV'
   end if 

! order quasidiagonal:
   i=0
   do ip=1,qm2ds%Np
    do ih=1,qm2ds%Nh
       i=i+1
       qm2ds%rrwork(i)=qm2ds%ehf(ih+qm2ds%Np)-qm2ds%ehf(ip)
    end do
   end do

   call rrdpsort(qm2ds%rrwork,qm2ds%Ncis,qm2ds%ix,1) 
!  Account for found vectors
   do j=1,j0
    do i=1,qm2ds%Ncis
       qm2ds%rrwork(i)=qm2ds%rrwork(i)+(qm2ds%fs+qm2ds%e0(j0)) &
          *abs(qm2ds%v0(i,j)**2-qm2ds%v0(qm2ds%Ncis+i,j)**2)
    end do
   end do

! try to find vector new vectors in the batch:
!write(6,*) 'Entering davidson0'
   call davidson0(qmmm_struct,qm2ds%Ncis,lprint,qm2ds%ftol0,qm2ds%ftol1,qm2ds%ferr, &
      qm2ds%Np,qm2ds%Nh,j0,j1, &
      qm2ds%e0,qm2ds%v0,kflag,qm2ds%ix, &
      qm2ds%rrwork(1),qm2ds%rrwork(2*qm2ds%Ncis+1), &
      qm2ds%rrwork(4*qm2ds%Ncis+1), &
      nd,nd1,qm2ds%vexp1,qm2ds%vexp,qm2ds%ray,qm2ds%rayv,qm2ds%rayvL, &
      qm2ds%rayvR,qm2ds%raye,qm2ds%raye1, &
      qm2ds%ray1,qm2ds%ray1a,qm2ds%ray2,qm2ds%idav,istore)

! Printing out found eigenvalues, error and tolerance
   if(lprint>0) then
      write(6,*)' i, e0(i), ferr(i), ftol0'
      do i=1,j0
         write(6,111) i,' +++ ',qm2ds%e0(i),qm2ds%ferr(i),qm2ds%ftol0
      end do
      write(6,*)'-------------------------------------------------'
      write(6,*) 
   end if

   if (qm2ds%mdflag.le.-3) qm2ds%Mj=j0
111   format (i3,a,g24.16,2(' ',e8.2))
   call flush(6)

   
! Write vectors only for BIG sizes in the case of crash/restart       	  
!   write(6,*) "Davidson Eigenvectors"
!   do j=1,j0
!     write(6, '(40e18.9)') (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
!   end do
   
   if(qm2ds%mdflag.lt.0.and.qm2ds%Nb.gt.100) then
      open (10,file='modes.b',form='unformatted')
      open (12,file='ee.b')

      do j=1,j0
         write (10) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
         write (12,*)  qm2ds%e0(j)
      end do

      close(10)
      close(12)
   end if

   goto 10
   end if
       
!
!--------------------------------------------------------------------
!
! End big loop
!
!--------------------------------------------------------------------
!
 70   continue

! end find many vectors
!--------------------------------------------------------------------
!
   if(qm2ds%mdflag.ne.0) then   ! initial MD point, store vector
!  normalize vectors
      do j=1,qm2ds%Mx
         fn=0
! CIS
         if(qm2ds%idav==1) then
            do i=1,qm2ds%Ncis
               qm2ds%v0(i+qm2ds%Ncis,j)=0.0
            end do
         end if
! end CIS

         fn=ddot(qm2ds%Ncis,qm2ds%v0(1,j),one,qm2ds%v0(1,j),one) &
            -ddot(qm2ds%Ncis,qm2ds%v0(qm2ds%Ncis+1,j),one, &
            qm2ds%v0(qm2ds%Ncis+1,j),one)

!         do i = 1,M4
!            fn = fn + v0(i,j)**2 - v0(i+M4,j)**2
!         enddo

         f=1/sqrt(abs(fn))
         call dscal(qm2ds%Nrpa,f,qm2ds%v0(1,j),one)

!         do i = 1,M2
!            v0(i,j) = f*v0(i,j)
!         enddo
!      write(6,*) j, f
      end do

!     write to the hard disk
      call rrdpsort(qm2ds%e0,qm2ds%Mx,qm2ds%kx,2)
!	  do i=1,M4
!	     xi_old(i)=v0(i,2)+v0(i+M4,2)
!	  enddo
!	  f2= ddot(M4,xi_old,one,xi_old,one)
!	  f2 = 1/sqrt(abs(f2))
!	  do i=1,M4
!	    xi_old(i)=f2*xi_old(i)
!	  enddo

      if (qm2ds%mdflag.lt.0) then
         open (10,file='modes.b',form='unformatted')
         open (12,file='ee.b')
!        open (11,file='modes.in',access='append')
!        write(11,210)
!        write(11,210) '$MODES'
!
         do j=1,qm2ds%Mx
            write (10) (qm2ds%v0(i,qm2ds%kx(j)),i=1,qm2ds%Nrpa)
            write (12,*)  qm2ds%e0(j)
!	         write(11,200) j,e0(j),1
         end do
 
         close(10)
         close(12)
!        write(11,210) '$ENDMODES'           
!        close(11)
!200     format(i7,g20.12,i4)
!210     format(A)
      end if
   end if

   if (qm2ds%mdflag.ge.0) then ! MD point only!!!!
!     store excited state vectors in AO representation in v2
!     istore=min(Mx,Mx_ss)
      if (istore.le.qm2ds%Mx) istore=qm2ds%Mx
      if (istore_M.le.qm2ds%Mx) istore_M=qm2ds%Mx
      if (istore.le.istore_M) istore=istore_M

      do j=1,qm2ds%Mx
         call mo2site(qm2ds%v0(1,j),qm2ds%v2(1,j),qm2ds%rrwork)
      end do
!	 if (Nb.gt.100) then
! --- write stored modes from the previous calculations in MD run
!         open (10,file='modes-ao.b',form='unformatted')       
!         do j=1,istore
!	    write (10) (v2(i,j),i=1,Nb*Nb) 
!	   enddo
!	   close(10)
!	 endif  	  	 
   end if

!      if (mdflag.ne.0) then
! Calculate the densities of the excited states.
! This is exact for CIS (davcis) and approximate for RPA (davrpa)
! In the matrix form it is given as [[xi^+,rho],eta]
!        f0=0.0
!          f1=1.0
!          f11=-1.0
!	  l=0
!         do j=1,Mx_s
!	  do k=1,Mx_s
!	   l=l+1
! Make [xi^+,rho]
!          do i = 1,M4
!            rrwork(Mb+i) = - v0(i+M4,j)
!                rrwork(Mb+M4+i) = v0(i,j)
!          enddo
! Here rrwork(Mb+1) contains [xi^+,rho] with 2*M4 elements

! Now expand matrices to square form
!         call expinter(Nb,Np,Nh,M4,rrwork(Mb+1),eta)
!         call expinter(Nb,Np,Nh,M4,v0(1,k),xi)

!        write(6,*) 'States', j, k
!	write (6,*) '[xi^+,rho]'
!        call prmat(Nb,eta)
!        write (6,*) 'xi'
!        call prmat(Nb,xi)

! Finally commute [eta,xi]
!         call dgemm('N','N',Nb,Nb,Nb,f1,eta,Nb,xi,Nb,f0,rrwork(1),Nb)
!         call dgemm('N','N',Nb,Nb,Nb,f11,xi,Nb,eta,Nb,f1,rrwork(1),Nb)
! Here rrwork(1) has excited state density (delta rho) of state j and k in MO
!         call copying(Mb,rrwork(1),v2(1,l))

!	write (6,*) '[[xi^+,rho],eta]'
!        call prmat(Nb,rrwork(1))       
!         enddo
!	enddo
!       endif
   return     
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine davidson0(qmmm_struct,M4,lprint,ftol0,ftol1,ferr, &
      Np,Nh,j0,j1,e0,v0,kflag,iee2,ee2,eta,xi, &
      nd,nd1,vexp1,vexp,ray,rayv,rayvL,rayvR,raye,raye1, &
      ray1,ray1a,ray2,idav,istore)
   use qm2_davidson_module   
  use cosmo_C,only:solvent_model 
  use qmmm_struct_module, only : qmmm_struct_type

   implicit none
     type(qmmm_struct_type), intent(inout) :: qmmm_struct

   !logical check_symmetry; !!JAB Testing

   integer Np,Nh,M4,lprint,kflag,nd,nd1,nd1_old,j0,j1
   integer one,i,j,k,n,m,icount,idav,istore,iloop,itarget
   integer info,iee2(qm2ds%Ncis)
   !integer l,u,c
   _REAL_ ftol0,ftol1,ferr(j1+j0),ddot,fn,f2m
   _REAL_ f0,f1,f2,f3,f4,tresh2
   _REAL_ e0(M4),v0(qm2ds%Nrpa,qm2ds%Mx),ee2(qm2ds%Ncis)
   _REAL_ eta(qm2ds%Nrpa),xi(qm2ds%Nrpa)
!  Davidson expansion vectors:
   _REAL_ vexp1(qm2ds%Ncis,nd),vexp(qm2ds%Nrpa,nd)

! Davidson Rayleigh matrices:
   _REAL_ ray1(nd,nd),ray2(nd,nd),ray1a(nd,nd)
   _REAL_ ray(nd,nd),raye(nd),raye1(nd)
   _REAL_ rayv(nd,nd),rayvL(nd,nd),rayvR(nd,nd)
   !_REAL_ f11;

   if (lprint.gt.4) write(6,*)' Entering davidson0'

   call clearing(nd*nd,ray1(1,1))
   call clearing(nd*nd,ray1a(1,1))
   call clearing(nd*nd,ray2(1,1))
   call clearing(nd*nd,ray(1,1))
   call clearing(nd*nd,rayv(1,1))
   call clearing(nd*nd,rayvL(1,1))
   call clearing(nd*nd,rayvR(1,1))
   call clearing(nd,raye)
   call clearing(nd,raye1)

   one=1

!	if (lprint.gt.1) write(6,*) 'Start Davidson'
   icount=0
   nd1_old=0
   m=0
   n=0
   iloop=0

   if (istore.gt.0) then    ! Use vectors from the previous step
      j1=min(j1,istore)
      goto 70 
   end if 

! **** Assign Davidson trial vectors

85 continue

   tresh2=qm2ds%tresh*0.9
   do j=1,nd   
      call clearing (2*M4,vexp(1,j))
      call clearing (M4,vexp1(1,j))
   end do 

! vexp and vexp1 are now zero !!JAB
! assign trial vectors based on the MO
   itarget =0
   do j=1,nd1
80    continue

      itarget=itarget+1
      if (itarget.ge.Np*Nh) goto 85 ! Restart vectors
      f1=0.0

      do i=1,j0
         f1=f1+v0(iee2(itarget),i)**2+v0(iee2(itarget)+M4,i)**2
      end do

      if (f1.ge.tresh2) goto 80  ! MO pair is not accepted
      if (lprint.gt.2) write(6,*)'++ START ',itarget,iee2(itarget)

      vexp1(iee2(itarget),j) = vexp1(iee2(itarget),j)+1.0
   end do

!  Orthogonolize trial vectors (vexp1) to found eigenvectors, some strange
!  function XI(V0) !!JAB
 
  !!ORTHONORMALIZE WITH XI=X+Y??
   
   do i=1,j0
      !!MAKE XI
      do k=1,M4
         xi(k)=v0(k,i)+v0(k+M4,i) !xi=X+Y??
      end do
      !!NORMALIZE
      f2=ddot(M4,xi,one,xi,one) !xi.xi
      f2=1.0D0/sqrt(abs(f2)) !normalization constant

      call dscal(M4,f2,xi,one) !normalize xi
     
      do j=1,nd1
         !!ORTHOGONALIZE
         f1=-ddot(M4,xi,one,vexp1(1,j),one) !xi.vexp1
         call daxpy(M4,f1,xi,one,vexp1(1,j),one) !vexp1=vexp1+xi.vexp1*xi
      end do

      f2=ddot(M4,vexp(1,j),one,vexp1(1,j),one) !vexp.vexp1 but vexp=0 if started from 85??
      f2=1.0D0/sqrt(abs(f2)) !normalization constant
      call dscal(M4,f2,vexp1(1,j),one) !normalize

      !!SAME AS ABOVE BUT WITH XI=X-Y
      do k=1,M4
         xi(k)=v0(k,i)-v0(k+M4,i) !xi=X-Y??
      end do

      f2=ddot(M4,xi,one,xi,one) !xi.xi
      f2=1.0D0/sqrt(abs(f2)) 
      call dscal(M4,f2,xi,one) !normalize xi

      do j=1,nd1   
         f1=-ddot(M4,xi,one,vexp1(1,j),one) 
         call daxpy(M4,f1,xi,one,vexp1(1,j),one)
      end do

      f2=ddot(M4,vexp(1,j),one,vexp1(1,j),one)
      f2=1.0D0/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do

!  ORTHOGONALIZE AND NORMALIZE TRIAL VECTORS
90     continue
   do j=1,nd1

      !!NORMALIZE
      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)

      do i=1,j-1
         
         !!ORTHOGONALIZE
         f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh1) then !!REMOVE VERY NONORTHOGONAL VECTORS
            call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
            nd1=nd1-1 !!NUMBER OF VECTORS
            if(lprint.gt.2) write(6,*) 'Trial removed #',j,i,f1,nd1
            goto 90
         end if
         
         call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one) !!ORTHOGONALIZE
      end do

      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one) !!NORMALIZE AGAIN
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do 

!  CHECK ORTHOGONALITY OF EXPANSIONS
   if(lprint.gt.3) write(6,*) 'Check expansion to expansion'
   do j=1,nd1
      do i=1,nd1
         f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
!	  if (i.eq.j) write(6,*) i,j,f1
!        if (f1.gt.1.0E-15.and.i.ne.j) write(6,*) i,j,f1
         if(f1.gt.1.0E-15.and.i.ne.j) goto 90
      end do
   end do

   if (nd1.le.j1) j1=nd1
! **** End assign Davidson trial vectors
10 continue
! **** Write some things, check the number of iterations, etc.
   icount=icount+1 !iteration counter
   if(lprint.gt.4) write(6,*) 'COUNT=',icount,'Exp=',nd1
   if(icount.gt.qm2ds%icount_M) then
         write(6,*) "Number of davidson iterations exceeded, exiting"
         stop
   endif
   if(lprint.gt.4) write(6,*) 'nd1,nd1_old,j0,j1',nd1,nd1_old,j0,j1

! **** Davidson Restart
   if(nd1.gt.nd) then !if the number of Krylov subspace expansions is  
      iloop=iloop+1
      if(lprint.gt.4) write(6,*)
      write(6,*) 'Davidson not converged, expansion=',nd
   end if

70     continue

   if(nd1.gt.nd.or.istore.gt.0) then  ! Use vectors from the previous step
      istore=0
      if (lprint.gt.4) write(6,*) 'Restart Davidson with guesses from the previous step'
!	   do j=1,j1
!	    write(6,*) 'Mode',j, e0(j)
!	    write(6,*) 'MO_p-h'
!	    call prarr(M4,v0(1,j)) 
!	    write(6,*) 'MO_h-p'
!	    call prarr(M4,v0(M4+1,j))
!         enddo
      if(lprint.gt.4) write(6,*) 'Restart loop = ',iloop
      icount=0
      nd1_old=0
      m=0
      n=0

      if(idav.eq.2) nd1=min((nd-2),2*j1)  ! RPA
      if(idav.eq.1) nd1=min((nd-2),j1)  ! CIS
      if(lprint.gt.1) write(6,*) 'Currently have ',j1,' states'
      if(lprint.gt.1) write(6,*) 'With ',nd1,' initial guesses'

      !NEW EXPANSION VECTORS
      if(idav.eq.2) then ! RPA 
         do j=1,j1
          do i=1,M4
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0)+v0(i+M4,j+j0)
             if((j+j1).le.nd1) vexp1(i,j+j1)=v0(i,j+j0)-v0(i+M4,j+j0)
          end do
         end do
      end if

      if(idav.eq.1) then ! CIS 
         do j=1,j1
          do i=1,M4
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0) 
          end do
         end do
      end if

!	   do j=1,nd1
!	    write(6,*) '#',j
!	    call prarr(M4,vexp1(1,j)) 
!	   enddo

!  Orthogonalize and normalize new expansion vectors
12    continue
      do j=1,nd1
         f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(M4,f2,vexp1(1,j),one)

         do i=1,j-1
            f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
            if (abs(f1).gt.qm2ds%tresh1) then
               call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
               nd1=nd1-1
               if(lprint.gt.2) write(6,*) 'New expansion removed #',j,i,f1,nd1
               goto 12
            end if

            call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
         end do

         f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(M4,f2,vexp1(1,j),one)
      end do

!	 write(6,*) 'Expansion vectors left', nd1

!  Check
!	 write(6,*) 'Check expansion to expansion'
      do j=1,nd1
       do i=1,nd1
          f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
          if (f1.gt.1.0E-15.and.i.ne.j.and.lprint.gt.0) write(6,*) i,j,f1
          if(f1.gt.1.0E-15.and.i.ne.j) goto 12
       end do
      end do

      goto 10
   end if

! **** End Davidson Restart	
! Next is calculated vexp=L(vexp1), so vexp becomes the Lxi of the guess xi !!JAB
! Vexp1 was the orthonormalized guess expansion vectors from the previous
! section. Throughout comments using A and B are related to  L(xi)=L([X;Y]) similarly 
! to [A,B;B,A][X;Y] from RPA eigenvalue problem, but definately not the band 
! ABBA as one might mistakenly assume !!JAB

   do i=nd1_old+1,nd1
      call clearing (2*M4,eta)
      call dcopy(M4,vexp1(1,i),one,eta,one)
      call Lxi_testing(qmmm_struct,eta,vexp(1,i),solvent_model)
      !call Lxi(eta,vexp(1,i),solvent_model)	
       
! CIS - set Y=0
      if(idav.eq.1) then
         do j=1,M4
            vexp(M4+j,i)=0.0
         enddo
      endif

! form vexp(1,i) having Ab and Bb to (A+B)b and (A-B)b	  
! thus vexp was L([X;Y], now it is [LX+LY;LX-LY] !!JAB
      do j=1,M4
         f1=vexp(j,i)
         vexp(j,i)=vexp(j,i)+vexp(M4+j,i)
         vexp(M4+j,i)=f1-vexp(M4+j,i)
      end do
   end do

!       stop
! **** Operations in Krylov space
! form ray1=b(A-B)b and ray2=b(A+B)b	
!That is, these are formed from vexp=L(xi) and the new guess expansions Vexp1
   do i=1,nd1
      do j=nd1_old+1,nd1
         ray1(i,j)=ddot(M4,vexp1(1,i),one,vexp(1,j),one)
         ray2(i,j)=ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
      end do
   end do
  
   if(nd1_old.ne.0) then
      do i=nd1_old+1,nd1
       do j=1,nd1
          ray1(i,j)=ddot(M4,vexp1(1,i),one,vexp(1,j),one)
          ray2(i,j)=ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
       end do
      end do
   end if

   nd1_old=nd1

     !!JAB Testing
     !do i=1,nd1 !!JAB Checking matrix elements
!	do j=1,nd1
!	write(6,*) 'Ray1',i,j,ray1(i,j)
!	enddo 
 !    enddo
     !!JAB Testing

! form ray1a=sqrt(b(A-B)b)
   f1=1.0
   f0=0.0

   call dcopy(nd*nd,ray1,one,rayvR,one)

   !call symmetr(nd,rayvR) !!JAB Test

   !if (qm2ds%verbosity > 4) 
   !print*,'1 qm2ds%Nrpa=',qm2ds%Nrpa

   !!JAB Testing
   !Check matrix symmetry as required by dsyev eigenproblem routine
        !if(.not.check_symmetry(rayvR,nd)) then
                        
        !write(6,*) 'M4',M4,'qm2ds%nb',qm2ds%nb,'size(vexp(:,1))',size(vexp(:,1))
        !write(6,*) 'TEST nd1    nd1     orb     orb        rayvR(i.j)      rayvR(j.i)      vexp    vexp1'
        !do i=1,nd1
        !do j=1,nd1
        !c=0
        !do k=1,qm2ds%nb
        !do l=1,qm2ds%nb/2
        !u=(k-1)*qm2ds%nb/2+l
        !c=c+1
        !write(6,"(5I4,8g15.5)")&
        ! i,j,k,l,c,rayvR(i,j),rayvR(j,i),vexp1(u,i),vexp(u,j),vexp1(u,j),vexp(u,i),&
        ! vexp1(u,i)*vexp(u,j),vexp1(u,j)*vexp(u,i)
        !end do
        !end do
        !end do
        !end do
        !write(6,*)"Error: rayvR non-symmetric. Something is wrong&
        !                 with the Liouville equation.(x_n|L(x_m)).ne.(x_m|L(x_n))"
        !stop
        !end if
   !!End JAB Testing
   write(6,*)'info00',info 
 
   call dsyev ('v','u',nd1,rayvR,nd,raye,xi,qm2ds%Nrpa,info) 
        !Eigenvalues of ray1 in raye and eigenvectors in rayvR
   write(6,*)'info000',info
   write(6,*)'info Nrpa,nd,nd1',qm2ds%Nrpa,nd,nd1
   write(6,*)'info shapes',shape(rayvR),shape(raye),shape(xi)
 
   do j=1,nd1
      raye1(j)=Sign(Sqrt(Abs(raye(j))),raye(j)) !make raye1 sqrt(eigenvalues) with same sign as eigenvalues
   end do

   do j=1,nd1
    do i=1,nd1
       rayv(i,j)=rayvR(j,i)*raye1(i)
    end do
   end do
 
   call dgemm('N','N',nd1,nd1,nd1,f1,rayvR,nd,rayv,nd,f0,ray1a,nd) !

!	write(6,*) 'ray1 eigenvalues'
!	call prarr(nd1,raye)
!	call prarr(nd1,raye1)
!	write(6,*) 'ray1a = sqrt(ray1)'
!	do i=1,nd1
!	 call prarr(nd1,ray1a(1,i))
!	enddo 

! form ray=Sqrt[b(A-B)b]*b(A+B)b*Sqrt[b(A-B)b] and Diagonalize

   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,ray1a,nd,f0,rayv,nd)
   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,ray,nd)
   call dcopy(nd*nd,ray,one,rayv,one)
   call symmetr(nd,rayv)

   !!JAB Testing
   !Check matrix symmetry as required by dsyev eigenproblem routine
   !     if(.not.check_symmetry(rayv,nd)) then
   !              write(6,*) "rayv symmetric"
   !     else
   !              write(6,*) "Error: rayv non-symmetric"
   !              !stop
   !     end if
   !!End JAB Testing

! find eigenvalues and eigenvectors of ray
   write(6,*)'info0',info
   call dsyev ('v','u',nd1,rayv,nd,raye,xi,qm2ds%Nrpa,info)
   do j=1,nd1
      raye(j) = Sign(Sqrt(Abs(raye(j))),raye(j))
   enddo

!      write(6,*) 'ray'
!	do i=1,nd1
!	 call prarr(nd1,ray(1,i))
!	enddo	
!	write(6,*) 'ray eigenvalues'
!	call prarr(nd1,raye)	
!	call prarr(nd1,raye1)
!      write(6,*) 'rayv'
!	do i=1,nd1
!	 call prarr(nd1,rayv(1,i))
!	enddo

! current EigenVectors rayv = Sqrt(1/b(A-B)b))|X+Y>
! Solve for Right EigenVector  rayvR = |X+Y>=ray1a*rayv
   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,rayvR,nd)
! Solve for Left EigenVector  
! rayvL = |X-Y> = 1/E b(A+B)b|X+Y> =1/raye * ray2 * rayvR
   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,rayvR,nd,f0,rayvL,nd)

   do j=1,nd1
    do i=1,nd1
       rayvL(i,j)=rayvL(i,j)/raye(j)
    end do
   end do
!      write(6,*) 'right rayvR'
!	do i=1,nd1
!	 call prarr(nd1,rayvR(1,i))
!	enddo
!      write(6,*) 'left rayvL'
!	do i=1,nd1
!	 call prarr(nd1,rayvL(1,i))
!	enddo

   if(lprint.gt.2) then 
      write(6,*)'info',info
      write(6,910) (raye(i),i=1,j1)
   end if 

   if(raye(1).le.0.1) goto 100           ! RPA bad behavior
! **** End operations in Krylov space

!  Form approximate eigenvectors
   do i=j0+1,j0+j1
      call clearing(2*M4,v0(1,i))
   end do

   do k=j0+1,j0+j1

 !new eigenvectors from converged onward (j0 is the number of converged vectors)
 !vexp1 is calculated at the beginning of the routine. The previous Krylov space
 !calculations result in f1 and f2 below which are used to mix vexp1 with
 !unconverged vectors in v0
      do i=1,nd1
         f1=rayvR(i,(k-j0))
         f2=rayvL(i,(k-j0))
         call daxpy(M4,f1,vexp1(1,i),one,v0(1,k),one) !V0=V0+f1*Vexp1
         call daxpy(M4,f2,vexp1(1,i),one,v0(1+M4,k),one)
      end do

      do i=1,M4 !Whatever is going on here, it looks like some X-Y, X+Y stuff
         f3=v0(i,k)
         v0(i,k)=-v0(i,k)-v0(M4+i,k)
         v0(M4+i,k)=v0(M4+i,k)-f3
      end do
       
! CIS : Y=0
      if(idav.eq.1) then
         do i=1,M4
            v0(i+M4,k)=0.0
         end do
      end if
   end do

!  Orthogonalize and normalize approximate eigenvectors
   do j=j0+1,j0+j1
      do i=j0+1,j-1
         f1=-ddot(M4,v0(1,i),one,v0(1,j),one) &
            +ddot(M4,v0(1+M4,i),one,v0(1+M4,j),one)
         call daxpy(2*M4,f1,v0(1,i),one,v0(1,j),one)
      end do

      f2=ddot(M4,v0(1,j),one,v0(1,j),one) &
         -ddot(M4,v0(1+M4,j),one,v0(1+M4,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(2*M4,f2,v0(1,j),one)

! 	  write(6,*) 'Mode',j, raye(j)
!  	  write(6,*) 'MO_p-h'
!  	  call prarr(M4,v0(1,j)) 
!  	  write(6,*) 'MO_h-p'
!  	  call prarr(M4,v0(M4+1,j))	  	   
   end do
 
! FIND eigenvalue, residual vectors and residual norm:
! print results of this iteration and check for converged vectors

   if (lprint.gt.1) write(6,*)'eigenvalues and residual norm'
   n=0
   m=0

   do j=j0+1,j0+j1
      
      call Lxi_testing(qmmm_struct,v0(1,j),eta,solvent_model) !L(xi) output in eta !!JAB
      !call Lxi(v0(1,j),eta,solvent_model)	
      !write(6,*)"loc(v0(1,j))=",j,loc(v0(1,j))	
      f1=ddot(M4,v0(1,j),one,eta(1),one) &
         -ddot(M4,v0(1+M4,j),one,eta(1+M4),one)
! CIS
      if (idav.eq.1) f1=ddot(M4,v0(1,j),one,eta(1),one) !E=<xi|L(xi)> now in f1 !!JAB
      call dcopy(2*M4,eta,one,xi,one)
      call daxpy(2*M4,-f1,v0(1,j),one,xi,one)
      f2=ddot(M4,xi(1),one,xi(1),one)
      f3=ddot(M4,xi(1+M4),one,xi(1+M4),one)
! CIS
      if(idav.eq.1) then !Clear Y after Liouville equation
         f3=0.0 
         call clearing(M4,xi(1+M4))
         call clearing(M4,eta(1+M4))
      endif

      f2=f2+f3 !f(x) + f(y)
      if(f2.le.ftol0) then ! Converged vector
         n=n+1
         call dcopy(2*M4,v0(1,j0+n),one,v0(1,j),one) !Move converged vector to beginning
         e0(j0+n)= f1 !move converged energy to beginning
         ferr(j0+n)=abs(f2)+abs(f3) !

         if(lprint.gt.3) then
            write(6,"(3i5,1x,2f14.9,2g10.3,2x,A)") &
               j,n,m,raye(j-j0),f1,f2,f3,' Converged!'
         end if
! CIS
      else if(idav.eq.1) then       

 ! FORM perturbed residual vectors for new initial guesses
         if((nd1+m).eq.nd) goto 45
         m=m+1
        
         !Not obvious how this is done, but f1/(f1-ee2(i) is some sort of energy
         !weighting and eta is L(xi) v0
         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)-f1*v0(i,j))/(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)

         if(lprint.gt.3) then 
                write(6,"(3i5,1x,2f14.9,2g10.3)") &
               j,n,m,raye(j-j0),f1,f2,f3
         end if
      else
         if((nd1+m).eq.nd) goto 45
         m=m+1

         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)-eta(i+M4)-f1*(v0(i,j)-v0(i+M4,j))) &
               /(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)

         if((nd1+m).eq.nd) goto 45
         m=m+1

         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)+eta(i+M4)-f1*(v0(i,j)+v0(i+M4,j))) &
               /(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)
         if(lprint.gt.1) write(6,"(3i5,1x,2f14.9,2g10.3)") &
            j,n,m,raye(j-j0),f1,f2,f3
      end if          
   end do
!          
45 continue
        
   if((j1-n).eq.0) then
      if(lprint.gt.4) write(6,*) 'All vectors found after loop' &
         ,iloop, ', Expansion ', nd1
      goto 100
   end if

   if(m.eq.0) then ! Restart Davidson
      nd1=nd+1
      goto 10
   end if

   if(lprint.gt.4) write(6,*) 'New perturbed m=',m
  
!  Orthogonalize and normalize residual vectors to found eigenvectors
!	 do j=nd1+1,nd1+m
!	  do i=1,j0
!	   f1= -0.5*(ddot(M4,v0(1,i),one,vexp1(1,j),one)
!     +       +ddot(M4,v0(1,i),one,vexp1(1,j),one))
!         f1= f1/ddot(M4,v0(1,i),one,v0(1,i),one)
!	   call daxpy(M4,f1,v0(1,i),one,vexp1(1,j),one) 
!	  enddo
!	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
!        f2 = 1/sqrt(abs(f2))
!	   call dscal(M4,f2,vexp1(1,j),one)
!	 enddo

!  Orthogonalize and normalize residual vectors to expansion vectors
15 continue

   do j=1+nd1,nd1+m
      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)

      do i=1,j-1
         f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh) then
            call dcopy(M4,vexp1(1,nd1+m),one,vexp1(1,j),one) 
            m=m-1     
            if(lprint.gt.2) write(6,*) 'exp removed #',j,i,f1,nd1+m
            goto 15
         end if

         call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
      end do

      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do
!  Check
!	 write(6,*) 'Check expansion to expansion'
   do j=1+nd1,nd1+m
    do i=1,j-1
       f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
!      if (f1.gt.1.0E-15) write(6,*) i,j,f1
       if(f1.gt.1.0E-15.and.i.ne.j) goto 15          
    end do
   end do

!           do j=1,nd1+m
!	     write(6,*) '#',j
!	     call prarr(M4,vexp1(1,j)) 
!           enddo

   if(lprint.gt.1) write(6,*) 'Perturbed left m=',m

   if(m.eq.0) then      
      write(6,*) 'Run out of expansion vectors'
!           do j=1,nd1
!	     write(6,*) '#',j
!	     write(6,*) 'MO_p-h'
!	     call prarr(M4,vexp(1,j)) 
!	     write(6,*) 'MO_h-p'
!	     call prarr(M4,vexp(M4+1,j))
!           enddo
      goto 100 
   end if

   if(iloop.gt.qm2ds%iloop_M) then
      write(6,*)'Davidson not converged after ',iloop,' loops'
      goto 100 
   end if

   if((nd1+m).gt.(nd-1).and.nd1.ne.nd) then
      nd1=nd
   else
      nd1=nd1+m
   endif

!  Check
!	 write(6,*) 'Check expansion to found'
!	 do j=1,nd1
!	  do i=1,j0
!	   f1= ddot(M4,v0(1,i),one,vexp(1,j),one)
!        if (f1.gt.1.0E-15) write(6,*) i,j,f1
!        enddo
!	 enddo

!  Check
!	 write(6,*) 'Check expansion to expansion'
!	 do j=1,nd1
!	  do i=1,nd1
!	   f1= ddot(M4,vexp(1,i),one,vexp(1,j),one)
!        if (f1.gt.1.0E-15) write(6,*) i,j,f1
!        enddo
!	 enddo

!      write(6,*) 'Expansion vectors' 
!	do j=1,nd1+m
!	 write(6,*) '#',j
!	 write(6,*) 'MO_p-h'
!	 call prarr(M4,vexp(1,j)) 
!	 write(6,*) 'MO_h-p'
!	 call prarr(M4,vexp(M4+1,j))
!	enddo
 
!	 stop
   goto 10

100   continue
   j0=j0+n
   if (n.eq.0) then
      write(6,*)'Could not go further ',j0,' vectors'
      kflag=1
   end if
    
   if(lprint.gt.0) then
      write(6,*) '@@@@ Davidson subroutine Found vectors',j0
!           do j=1,j0
!	      write(6,*) 'Mode',j, e0(j)
!	      write(6,*) 'MO_p-h'
!	      call prarr(M4,v0(1,j)) 
!	      write(6,*) 'MO_h-p'
!	      call prarr(M4,v0(M4+1,j))
!           enddo
   end if

910   format(' ', 100g11.3)

   return
   end
!
