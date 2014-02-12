	subroutine davidson(mdflag)
	
	implicit none
	
      include 'parH.par'   
      include 'const.par'
      include 'int.var'
      include 'hf.var'      
      include 'rho.var'
      include 'lanczos.var'
      include 'davidson.var'
      include 'etaxi.var'
      include 'modes.var'
      include 'parvar.var'
      include 'hf.cmn'
      include 'rho.cmn'
      include 'lanczos.cmn'
      include 'etaxi.cmn'
      include 'davidson.cmn'
      include 'modes.cmn'
      include 'parvar.cmn' 
      
      real*8 random,rranset,rranf
      real*8 ferr1,ferr2,fn,ferr(Mx_M)
      real*8 f,f1,f2,fo1,fo2,ddot,f0,f11
	
      integer iseed,ip,ih,j0,Mx0,imin,one,iloops
      integer kflag,nd,nd1,nd1_old,j1,mdflag,istore,istore_M
	save istore,istore_M
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize Davidson:

	one=1 
!	lprint = 1	       

! set Davidson parameters
!      if (Mx.ge.irflag) Mx=irflag-1
      nd=M4_M*Nd_M/(M4+5)
	if (nd.gt.Nd_MM) nd=nd_MM
      iloops=0
      Mx=min(M4,Mx)
! Split Mx into batches of j1 size
      j1=nd/idavfrac
	j1=min(j1,Mx)
	nd1=min(j1+2,nd/2)
	if (mdflag.ge.0) j1=Mx ! No shift is allowed for MD points
      
!	if (irflag.gt.1.and.mdflag.lt.0) then
 	if (irflag.gt.1) then
! Calculations will start from irflag state for ceo
	 j0=irflag-1
	 print *, j0,Mx,irflag
	 if (j0.gt.Mx) then
	   print *, 'Looks like irflag= ', irflag
	   print *, 'Show that states Mx= ', Mx
	   print *, 'are already calculated, exiting'
	   stop 
	  endif	
!  Read irflag-1 states 
         open (10,file='modes.b',form='unformatted',status='old')
	   open (12,file='ee.b',status='old')
         do j=1,j0
          read (10) (v0(i,j),i=1,M2)
	    read (12,*) e0(j)
         enddo   
         close(10)
	   close(12)
	 if (j0.eq.Mx) goto 70
	else    ! Assume start from the beginning 
	 j0=0
	endif
	 Mj=j0
			
	if (lprint.gt.1) then
	 print *, 'Davidson parameters'
	 print *, 'mdflag=',mdflag
	 print *, 'irflag=',irflag
	 print *, 'M4=',M4
	 print *, 'Mx=',Mx
	 print *, 'Mj=',Mj
	 print *, 'nd=',nd 
	 print *, 'nd1=',nd1
	 print *, 'j1=',j1
	endif	 	

!	if (irflag.gt.1.and.Nb.gt.100.and.mdflag.gt.0) then
!	   istore=min(Mx,Mx_ss)
! --- read stored modes from the previous calculations in MD run
!         open (10,file='modes-ao.b',form='unformatted',status='old')       
!         do j=1,istore
!	    read (10) (v2(i,j),i=1,Nb*Nb) 
!	   enddo
!	   close(10)	  
!	endif	 

      if (istore.gt.0) then ! MD point only!!!!	 
!      recover excited state vectors from AO representation in v2
       do j=1,istore
	   call site2mo (rrwork,v2(1,j),v0(1,j))
	 enddo
	endif

! ++++ begin big loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
! find many vectors:
10	continue
      if (j0.lt.Mx) then 

	    iloops=iloops+1
!	    if (iloops.eq.2) j1=j1/2 ! After the first one convergence is slow
	    j1=min(j1,(Mx-j0))
	    nd1=min(j1+2,nd/2)
	    if (lprint.gt.0) then
	     print *
	     print *, 'Davidson batch ', iloops
	     print *, 'So far found ', j0, ' states'
	     print *, 'out of requested ', Mx, ' states'
	     print *, 'This batch will seek',j1,' vectors'
	     if(j0.gt.0) print *, 'Shift is',fs+e0(j0),' eV'
	    endif 

! order quasidiagonal:
        i = 0
        do ip = 1,Np
         do ih = 1,Nh
           i = i + 1
           rrwork(i) = ee(ih+Np) - ee(ip)
         enddo
        enddo	
        call rrdpsort (rrwork,M4,ix,1) 
! Account for found vectors
	   do j=1,j0
          do i = 1,M4
           rrwork(i)=rrwork(i)+(fs+e0(j0))*abs(v0(i,j)**2-v0(M4+i,j)**2)
          enddo
	   enddo
! try to find vector new vectors in the batch:

	   call davidson0 (M4,lprint,ftol0,ftol1,ferr,Np,Nh,j0,j1, &
      	 e0,v0,kflag,ix,rrwork(1),rrwork(2*M4+1),rrwork(4*M4+1), &
             nd,nd1,vexp1,vexp,ray,rayv,rayvL,rayvR,raye,raye1, &
             ray1,ray1a,ray2,idav,istore)

       
! write vector:
	   do i=1,j0
          if (lprint.gt.0) print 111,i,' +++ ',e0(i),ftol0,ferr(i)
         enddo
	   if (mdflag.le.-3) Mj=j0
 111     format (i3,a,g24.16,2(' ',e8.2))
         call flush(6)

! Write vectors only for BIG sizes in the case of crash/restart       	  
	  if (mdflag.lt.0.and.Nb.gt.100) then
         open (10,file='modes.b',form='unformatted')
	   open (12,file='ee.b')
         do j=1,j0
          write (10) (v0(i,j),i=1,M2)
	    write (12,*)  e0(j)
         enddo
         close(10)
	   close(12)
	  endif
	  
       goto 10
       endif
! ++++ End big loop
 70     continue

! end find many vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (mdflag.ne.0) then   ! initial MD point, store vector
!     normalize vectors
      do j = 1,Mx
         fn = 0
! CIS
      if (idav.eq.1) then
         do i = 1,M4
            v0(i+M4,j) = 0.0
         enddo
      endif
! end CIS
        fn= ddot(M4,v0(1,j),one,v0(1,j),one)- &
            ddot(M4,v0(M4+1,j),one,v0(M4+1,j),one)
!         do i = 1,M4
!            fn = fn + v0(i,j)**2 - v0(i+M4,j)**2
!         enddo
         f = 1/sqrt(abs(fn))
	   call dscal (M2,f,v0(1,j),one)
!         do i = 1,M2
!            v0(i,j) = f*v0(i,j)
!         enddo
!      print *, j, f
      enddo

!    write to the hard disk
      call rrdpsort (e0,Mx,kx,2)
!	  do i=1,M4
!	     xi_old(i)=v0(i,2)+v0(i+M4,2)
!	  enddo
!	  f2= ddot(M4,xi_old,one,xi_old,one)
!	  f2 = 1/sqrt(abs(f2))
!	  do i=1,M4
!	    xi_old(i)=f2*xi_old(i)
!	  enddo
       	  
	  if (mdflag.lt.0) then
         open (10,file='modes.b',form='unformatted')
	   open (12,file='ee.b')
!        open (11,file='modes.in',access='append')
!        write(11,210)
!        write(11,210) '$MODES'
!
         do j=1,Mx
          write (10) (v0(i,kx(j)),i=1,M2)
	    write (12,*)  e0(j)
!	    write(11,200) j,e0(j),1
         enddo
	   
         close(10)
	   close(12)
!        write(11,210) '$ENDMODES'           
!        close(11)
!200     format(i7,g20.12,i4)
!210     format(A)
	  endif
       endif

      if (mdflag.ge.0) then ! MD point only!!!!
!      store excited state vectors in AO representation in v2
!       istore=min(Mx,Mx_ss)
       if (istore.le.Mx) istore=Mx
       if (istore_M.le.Mx) istore_M=Mx
       if (istore.le.istore_M) istore=istore_M
       do j=1,Mx
	   call mo2site (v0(1,j),v2(1,j),rrwork)
	 enddo
!	 if (Nb.gt.100) then
! --- write stored modes from the previous calculations in MD run
!         open (10,file='modes-ao.b',form='unformatted')       
!         do j=1,istore
!	    write (10) (v2(i,j),i=1,Nb*Nb) 
!	   enddo
!	   close(10)
!	 endif  	  	 
	endif

 	 
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

!        print *, 'States', j, k
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


	   subroutine davidson0 (M4,lprint,ftol0,ftol1,ferr, &
      	Np,Nh,j0,j1,e0,v0,kflag,iee2,ee2,eta,xi, &
            nd,nd1,vexp1,vexp,ray,rayv,rayvL,rayvR,raye,raye1, &
             ray1,ray1a,ray2,idav,istore)
     
	implicit none	
      include 'parH.par'
			
	integer Np,Nh,M4,lprint,kflag,nd,nd1,nd1_old,j0,j1
	integer one,i,j,k,n,m,icount,idav,istore,iloop,itarget
	integer info,ix(2*Nd_M),iee2(M4)
	real*8 ftol0,ftol1,ferr(j1+j0),ddot,fn,fm
	real*8 f0,f1,f2,f3,f4,tresh2
	real*8 e0(M4),v0(2*M4_M,Mx_M),ee2(M4),eta(2*M4),xi(2*M4)
! Davidson expansion vectors:
      real*8 vexp1(M4,nd),vexp(2*M4,nd)

! Davidson Rayleigh matrices:
      real*8 ray1(nd,nd),ray2(nd,nd),ray1a(nd,nd)
	real*8 ray(nd,nd),raye(nd),raye1(nd)
	real*8 rayv(nd,nd),rayvL(nd,nd),rayvR(nd,nd)
	
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
		
!	if (lprint.gt.1) print *, 'Start Davidson'
	icount=0
	nd1_old=0
	m=0
	n=0
	iloop=0
	
      if (istore.gt.0) then    ! Use vectors from the previous step
        j1=min(j1,istore)
	  goto 70 
	endif  
	
! **** Assign Davidson trial vectors

85      continue	
        tresh2=tresh*0.9 	
 	  do j=1,nd   
	   call clearing (2*M4,vexp(1,j))
	   call clearing (M4,vexp1(1,j))
        enddo 
	 
! assign trial vectors based on the MO
        itarget =0			
	  do j=1,nd1
80       continue
	   itarget=itarget+1
	   if (itarget.ge.Np*Nh) goto 85 ! Restart vectors
	   f1=0.0	   	 
         do i=1,j0
	    f1=f1+v0(iee2(itarget),i)**2+v0(iee2(itarget)+M4,i)**2
	   enddo
	   if (f1.ge.tresh2) goto 80  ! MO pair is not accepted
	   if (lprint.gt.2) print *,'++ START ',itarget,iee2(itarget)
	   vexp1(iee2(itarget),j) = vexp1(iee2(itarget),j)+1.0
	  enddo

! Orthogonolize trial vectors to found eigenvectors
	 do i=1,j0
	    do k=1,M4
	     xi(k)=v0(k,i)+v0(k+M4,i)
	    enddo
	    f2= ddot(M4,xi,one,xi,one)	  
	    f2 = 1/sqrt(abs(f2))
	    call dscal(M4,f2,xi,one) 
	   do j=1,nd1	   
	    f1= -ddot(M4,xi,one,vexp1(1,j),one)
	    call daxpy(M4,f1,xi,one,vexp1(1,j),one)
	   enddo
	   f2= ddot(M4,vexp(1,j),one,vexp1(1,j),one)
         f2 = 1/sqrt(abs(f2))
	   call dscal(M4,f2,vexp1(1,j),one)
	    do k=1,M4
	     xi(k)=v0(k,i)-v0(k+M4,i)
	    enddo
	    f2= ddot(M4,xi,one,xi,one)	  
	    f2 = 1/sqrt(abs(f2))
	    call dscal(M4,f2,xi,one) 
	   do j=1,nd1	   
	    f1= -ddot(M4,xi,one,vexp1(1,j),one)
	    call daxpy(M4,f1,xi,one,vexp1(1,j),one)
	   enddo
	   f2= ddot(M4,vexp(1,j),one,vexp1(1,j),one)
         f2 = 1/sqrt(abs(f2))
	   call dscal(M4,f2,vexp1(1,j),one)
	 enddo        

!  Orthogonalize and normalize trial vectors
90     continue
	 do j=1,nd1
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)	  
	  do i=1,j-1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if (abs(f1).gt.tresh1) then
	     call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
	     nd1=nd1-1	     
	     if (lprint.gt.2) print *, 'Trial removed #',j,i,f1,nd1
	     goto 90
	   endif
	   call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
	  enddo
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)
	 enddo	 

!  Check
!	 print *, 'Check expansion to expansion'
	 do j=1,nd1
	  do i=1,nd1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
!	  if (i.eq.j) print *, i,j,f1
!        if (f1.gt.1.0E-15.and.i.ne.j) print *, i,j,f1
	  if (f1.gt.1.0E-15.and.i.ne.j) goto 90
        enddo
	 enddo
	 if (nd1.le.j1) j1=nd1
! **** End assign Davidson trial vectors

10	continue
     
	icount=icount+1
	if (lprint.gt.1) print *, 'COUNT=',icount,'Exp=',nd1
	if (icount.gt.icount_M) stop
	
	if (lprint.gt.1) print *, 'nd1,nd1_old,j0,j1',nd1,nd1_old,j0,j1
! **** Davidson Restart
	if (nd1.gt.nd) then 	 
	 iloop=iloop+1
	 if (lprint.gt.0) print *
	  print*, 'Davidson not converged, expansion=',nd
      endif
70     continue
      if (nd1.gt.nd.or.istore.gt.0) then  ! Use vectors from the previous step
	 istore=0
       if (lprint.gt.1) print*, 'Restart Davidson with improved guesses'
!	   do j=1,j1
!	    print *, 'Mode',j, e0(j)
!	    print *, 'MO_p-h'
!	    call prarr(M4,v0(1,j)) 
!	    print *, 'MO_h-p'
!	    call prarr(M4,v0(M4+1,j))
!         enddo
	 if (lprint.gt.0) print*, 'Restart loop = ',iloop
	 icount=0
	 nd1_old=0
	 m=0
	 n=0
	 if (idav.eq.2) nd1=min((nd-2),2*j1)	 ! RPA
	 if (idav.eq.1) nd1=min((nd-2),j1)	 ! CIS
	 if (lprint.gt.1) print *, 'Looking after ',j1,' states'
	 if (lprint.gt.1) print *, 'With ',nd1,' initial guesses'	
	 if (idav.eq.2) then                       ! RPA 
        do j=1,j1
	   do i=1,M4
	    if (j.le.nd1) vexp1(i,j)=v0(i,j+j0)+v0(i+M4,j+j0)
	    if ((j+j1).le.nd1) vexp1(i,j+j1)=v0(i,j+j0)-v0(i+M4,j+j0)
	   enddo
	  enddo
	 endif
	 if (idav.eq.1) then                       ! CIS 
        do j=1,j1
	   do i=1,M4
	    if (j.le.nd1) vexp1(i,j)=v0(i,j+j0)
	   enddo
	  enddo
	 endif	
!	   do j=1,nd1
!	    print *, '#',j
!	    call prarr(M4,vexp1(1,j)) 
!	   enddo

!  Orthogonalize and normalize new expansion vectors
12     continue
	 do j=1,nd1
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)	  
	  do i=1,j-1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if (abs(f1).gt.tresh1) then
	     call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
	     nd1=nd1-1	     
	     if (lprint.gt.2) print *, 'New expansion removed #',j,i,f1,nd1
	     goto 12
	   endif
	   call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
	  enddo
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)
	 enddo	 
!	 print *, 'Expansion vectors left', nd1

!  Check
!	 print *, 'Check expansion to expansion'
	 do j=1,nd1
	  do i=1,nd1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
!	  if (i.eq.j) print *, i,j,f1
!        if (f1.gt.1.0E-15.and.i.ne.j) print *, i,j,f1
	  if (f1.gt.1.0E-15.and.i.ne.j) goto 12
        enddo
	 enddo

	  goto 10
	endif	
! **** End Davidson Restart	

	do i=nd1_old+1,nd1
	 call clearing (2*M4,eta)
	 call dcopy(M4,vexp1(1,i),one,eta,one)
	 call Lxi (eta,vexp(1,i))
! CIS
       if (idav.eq.1) then
        do j=1,M4
	   vexp(M4+j,i) = 0.0
	  enddo
       endif
! form vexp(1,i) having Ab and Bb to (A+B)b and (A-B)b	  
        do j=1,M4
	   f1=vexp(j,i)
	   vexp(j,i)=vexp(j,i)+vexp(M4+j,i)
	   vexp(M4+j,i)=f1-vexp(M4+j,i)
	  enddo
	enddo
!       stop
! **** Operations in Krylov space
! form ray1=b(A-B)b and ray2=b(A+B)b	
	do i=1,nd1
	 do j=nd1_old+1,nd1
	  ray1(i,j) = ddot(M4,vexp1(1,i),one,vexp(1,j),one)
	  ray2(i,j) = ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
	 enddo
	enddo
  
	if (nd1_old.ne.0) then
	do i=nd1_old+1,nd1
	 do j=1,nd1
	  ray1(i,j) = ddot(M4,vexp1(1,i),one,vexp(1,j),one)
	  ray2(i,j) = ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
	 enddo
	enddo
	endif
	
	nd1_old=nd1
	
!      print *, 'ray1'
!	do i=1,nd1
!	 call prarr(nd1,ray1(1,i))
!	enddo        
!      print *, 'ray2'
!	do i=1,nd1
!	 call prarr(nd1,ray1(1,i))
!	enddo 

! form ray1a=sqrt(b(A-B)b)
	f1=1.0
	f0=0.0
      call dcopy(nd*nd,ray1,one,rayvR,one)
      call dsyev ('v','u',nd1,rayvR,nd,raye,xi,2*M4_M,info)
	    	
	do j=1,nd1
         raye1(j) = Sign(Sqrt(Abs(raye(j))),raye(j))
	enddo	
	do j=1,nd1
        do i=1,nd1
         rayv(i,j) = rayvR(j,i)*raye1(i)
	  enddo
	enddo
	call dgemm('N','N',nd1,nd1,nd1,f1,rayvR,nd,rayv,nd,f0,ray1a,nd)

!	print *, 'ray1 eigenvalues'
!	call prarr(nd1,raye)
!	call prarr(nd1,raye1)
!	print *, 'ray1a = sqrt(ray1)'
!	do i=1,nd1
!	 call prarr(nd1,ray1a(1,i))
!	enddo 

! form ray=Sqrt[b(A-B)b]*b(A+B)b*Sqrt[b(A-B)b] and Diagonalize
      call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,ray1a,nd,f0,rayv,nd)
	call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,ray,nd)
	call dcopy(nd*nd,ray,one,rayv,one)
	call symmetr(nd,rayv)
	call dsyev ('v','u',nd1,rayv,nd,raye,xi,2*M4_M,info)
	do j=1,nd1
         raye(j) = Sign(Sqrt(Abs(raye(j))),raye(j))
	enddo	      

!      print *, 'ray'
!	do i=1,nd1
!	 call prarr(nd1,ray(1,i))
!	enddo	
!	print *, 'ray eigenvalues'
!	call prarr(nd1,raye)	
!	call prarr(nd1,raye1)
!      print *, 'rayv'
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
        rayvL(i,j) = rayvL(i,j)/raye(j)
	 enddo
	enddo	
!      print *, 'right rayvR'
!	do i=1,nd1
!	 call prarr(nd1,rayvR(1,i))
!	enddo
!      print *, 'left rayvL'
!	do i=1,nd1
!	 call prarr(nd1,rayvL(1,i))
!	enddo
				
      if (lprint.gt.2) print *, 'info', info
	if (lprint.gt.2) write (6,910) (raye(i),i=1,j1)
	
      if (raye(1).le.0.1) goto 100           ! RPA bad behavior
! **** End operations in Krylov space

!  Form approximate eigenvectors
      do i=j0+1,j0+j1
	 call clearing(2*M4,v0(1,i))
	enddo
	
	do k=j0+1,j0+j1
	 do i=1,nd1
	   f1 = rayvR(i,(k-j0))
	   f2 = rayvL(i,(k-j0))
	   call daxpy(M4,f1,vexp1(1,i),one,v0(1,k),one)
	   call daxpy(M4,f2,vexp1(1,i),one,v0(1+M4,k),one)
!!	  do j=1,M4
!!	   v0(j,k)=v0(j,k)+rayv(i,ix(k-j0))*vexp(j,i)
!!	   v0(j+M4,k)=v0(j+M4,k)-rayv(i+nd1,ix(k-j0))*vexp(j,i)
!	   v0(j+M4,k)=v0(j+M4,k)+rayv(i+nd1,ix(k-j0))*vexp(j+M4,i)
!	   v0(j,k)=v0(j,k)-rayv(i,ix(k-j0))*vexp(j+M4,i)
!!	  enddo
	 enddo
          do i=1,M4
           f3=v0(i,k)
	     v0(i,k)=-v0(i,k)-v0(M4+i,k)
	     v0(M4+i,k)=v0(M4+i,k)-f3
          enddo	 
       
! CIS : Y=0
       if (idav.eq.1) then
          do i=1,M4
           v0(i+M4,k)=0.0
          enddo
         endif
        enddo
		 
!  Orthogonalize and normalize approximate eigenvectors to found ones
	 do j=j0+1,j0+j1
	  do i=j0+1,j-1
	   f1= -ddot(M4,v0(1,i),one,v0(1,j),one) &
      	 +ddot(M4,v0(1+M4,i),one,v0(1+M4,j),one)
         call daxpy(2*M4,f1,v0(1,i),one,v0(1,j),one)
	  enddo
	 f2= ddot(M4,v0(1,j),one,v0(1,j),one) &
            -ddot(M4,v0(1+M4,j),one,v0(1+M4,j),one)
	 f2 = 1/sqrt(abs(f2))
	  call dscal(2*M4,f2,v0(1,j),one)

! 	  print *, 'Mode',j, raye(j)
!  	  print *, 'MO_p-h'
!  	  call prarr(M4,v0(1,j)) 
!  	  print *, 'MO_h-p'
!  	  call prarr(M4,v0(M4+1,j))	  	   
	enddo
 
! find eigenvalue, residual vectors and residual norm:
 	 if (lprint.gt.1) print *, 'eigenvalues and residual norm'
	 n=0
	 m=0
        do j = j0+1,j0+j1
         call Lxi (v0(1,j),eta)
	    f1 = ddot(M4,v0(1,j),one,eta(1),one)- &
      	   ddot(M4,v0(1+M4,j),one,eta(1+M4),one)
! CIS
         if (idav.eq.1) f1 = ddot(M4,v0(1,j),one,eta(1),one)
	    call dcopy(2*M4,eta,one,xi,one)
	    call daxpy(2*M4,-f1,v0(1,j),one,xi,one)
	    f2 = ddot(M4,xi(1),one,xi(1),one)
	    f3 = ddot(M4,xi(1+M4),one,xi(1+M4),one)
! CIS
         if (idav.eq.1) then
	       f3=0.0
             call clearing(M4,xi(1+M4))
		 call clearing(M4,eta(1+M4))
           endif
	    f2=f2+f3
    	   
	   if (f2.le.ftol0) then	       ! Converged vector
           n=n+1
	      call dcopy(2*M4,v0(1,j0+n),one,v0(1,j),one)
		 e0(j0+n)= f1
		 ferr(j0+n)=abs(f2)+abs(f3)
             if (lprint.gt.1) print "(3i5,1x,2f14.9,2g10.3,2x,A)",  &
              j,n,m,raye(j-j0),f1,f2,f3, ' Converged!'
! CIS
	   else if (idav.eq.1)  then        ! Form perturbed residual vectors
	    if ((nd1+m).eq.nd) goto 45
	     m=m+1
	    do i = 1,M4
	     vexp1(i,nd1+m)=(eta(i)-f1*v0(i,j))/(f1-ee2(i))
	    enddo
	    f4= ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
	     f4 = 1/sqrt(abs(f4))
	     call dscal(M4,f4,vexp1(1,nd1+m),one)
           if (lprint.gt.1) print "(3i5,1x,2f14.9,2g10.3)", &
              j,n,m,raye(j-j0),f1,f2,f3
	   else
	    if ((nd1+m).eq.nd) goto 45
	     m=m+1
	    do i = 1,M4
	     vexp1(i,nd1+m)=(eta(i)-eta(i+M4)-f1*(v0(i,j)-v0(i+M4,j))) &
                          /(f1-ee2(i))
	    enddo
	    f4= ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
	     f4 = 1/sqrt(abs(f4))
	     call dscal(M4,f4,vexp1(1,nd1+m),one)

	    if ((nd1+m).eq.nd) goto 45
	     m=m+1
	    do i = 1,M4
	     vexp1(i,nd1+m)=(eta(i)+eta(i+M4)-f1*(v0(i,j)+v0(i+M4,j))) &
                          /(f1-ee2(i))
	    enddo
	    f4= ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
	     f4 = 1/sqrt(abs(f4))
	     call dscal(M4,f4,vexp1(1,nd1+m),one)
	   if (lprint.gt.1) print "(3i5,1x,2f14.9,2g10.3)", &
              j,n,m,raye(j-j0),f1,f2,f3
	   endif   	        
	   enddo          
45      continue
        
	  if ((j1-n).eq.0) then
	   if (lprint.gt.0) print *, 'All vectors found after loop' &
                               ,iloop, ', Expansion ', nd1
	   goto 100
	  endif
        if (m.eq.0) then ! Restart Davidson
	   nd1=nd+1
	   goto 10
	  endif
	  if (lprint.gt.1) print *, 'New perturbed m=',m


	   
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

!  Orthogonalize and normilize residual vectors to expansion vectors
15     continue
	 do j=1+nd1,nd1+m
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)
	  do i=1,j-1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if (abs(f1).gt.tresh) then
	     call dcopy(M4,vexp1(1,nd1+m),one,vexp1(1,j),one) 
	     m=m-1	     
	     if (lprint.gt.2) print *, 'exp removed #',j,i,f1,nd1+m
	     goto 15
	   endif
	   call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
	  enddo
	  f2= ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
        f2 = 1/sqrt(abs(f2))
	  call dscal(M4,f2,vexp1(1,j),one)
	 enddo
!  Check
!	 print *, 'Check expansion to expansion'
	 do j=1+nd1,nd1+m
	  do i=1,j-1
	   f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
!        if (f1.gt.1.0E-15) print*, i,j,f1
	   if (f1.gt.1.0E-15.and.i.ne.j) goto 15          
        enddo
	 enddo

!           do j=1,nd1+m
!	     print *, '#',j
!	     call prarr(M4,vexp1(1,j)) 
!           enddo

      if (lprint.gt.1) print *, 'Perturbed left m=',m
	
	 if (m.eq.0) then      
	   print *, 'Run out of expansion vectors'
!           do j=1,nd1
!	     print *, '#',j
!	     print *, 'MO_p-h'
!	     call prarr(M4,vexp(1,j)) 
!	     print *, 'MO_h-p'
!	     call prarr(M4,vexp(M4+1,j))
!           enddo
	   goto 100 
	 endif
	  if (iloop.gt.iloop_M) then
	   print *, 'Davidson not converged after ', iloop, ' loops'
	   goto 100 
	  endif

	 if ((nd1+m).gt.(nd-1).and.nd1.ne.nd) then
	  nd1=nd
	 else
	  nd1=nd1+m
	 endif
!  Check
!	 print *, 'Check expansion to found'
!	 do j=1,nd1
!	  do i=1,j0
!	   f1= ddot(M4,v0(1,i),one,vexp(1,j),one)
!        if (f1.gt.1.0E-15) print*, i,j,f1
!        enddo
!	 enddo

!  Check
!	 print *, 'Check expansion to expansion'
!	 do j=1,nd1
!	  do i=1,nd1
!	   f1= ddot(M4,vexp(1,i),one,vexp(1,j),one)
!        if (f1.gt.1.0E-15) print*, i,j,f1
!        enddo
!	 enddo

!      print *, 'Expansion vectors' 
!	do j=1,nd1+m
!	 print *, '#',j
!	 print *, 'MO_p-h'
!	 call prarr(M4,vexp(1,j)) 
!	 print *, 'MO_h-p'
!	 call prarr(M4,vexp(M4+1,j))
!	enddo
	 
!	 stop
	 goto 10

 100  continue
          j0=j0+n
          if (n.eq.0) then
	     print *, 'Could not go further ',j0,' vectors'
	     kflag=1
	    endif
	    
	    if (lprint.gt.0) then
	     print *, '@@@@ Found vectors',j0
!           do j=1,j0
!	      print *, 'Mode',j, e0(j)
!	      print *, 'MO_p-h'
!	      call prarr(M4,v0(1,j)) 
!	      print *, 'MO_h-p'
!	      call prarr(M4,v0(M4+1,j))
!           enddo
	    endif
 910  format(' ', 100g11.3)	
	return
	end
	
