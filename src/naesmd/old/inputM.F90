!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         This program read input parameters
!         and Cartesian coordinates.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine input1(keywr,Nat1,xx1,yy1,zz1,atm1,mdflag)
      implicit none
      
	integer Nat1,atm1(Nat1),mdflag
	real*8 xx1(Nat1),yy1(Nat1),zz1(Nat1)
	
      include 'parH.par'    
      include 'const.par'
      include 'scf.par' 
      include 'int.var'
      include 'real.var'
      include 'cart.var'
      include 'coulombM.var'
      include 'envir.var'
      include 'etaxi.var'
      include 'hf.var'      
      include 'rho.var'
      include 'lanczos.var'
      include 'modes.var'
      include 'parvar.var'
	include 'dipole.var'
	include 'derivfl.cmn'
	include 'dipole.cmn'
      include 'cart.cmn'
      include 'coulombM.cmn'
      include 'etaxi.cmn'
      include 'hf.cmn'
      include 'rho.cmn'
      include 'lanczos.cmn'
      include 'modes.cmn'
      include 'envir.cmn'
      include 'parvar.cmn'   
      character*(150) txt, keywr*40
	logical first
      data first /.true./
      save first
   
      Nat=Nat1
	do j=1,Nat1
	 atm(j)=atm1(j)
	 xx(j)=xx1(j)
	 yy(j)=yy1(j)
	 zz(j)=zz1(j)
	enddo
      j=0
      k=0

   print*,'kav - inputM'
	
! --- if this is the first time in this routine, load initial parameters
      if (first) then
!      open (11,file='modes.in')
      open (12,file='input.ceon',status='old')
	
10    continue      
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)
      if (txt(1:6).ne.'$PARAM') goto 10
      read(12,'(a)',err=29) txt            ! Model Hamiltonian type
!      write(11,'(a)') txt(1:80)           
	keywr=txt(1:15)
	if (index(keywr,'INDO').NE.0.AND.index(keywr,'MINDO').EQ.0) then
	  print *, keywr,' hamiltonian requested'
	  print *, 'Use ceoDavZ program for ', keywr
	  stop
	endif
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)  
      read(txt,*,err=29) inttyp          ! Type of Computation
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29)(intf(i),i=1,6)  ! Integral factors (if inttyp=2)
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) irflag          ! flag for ground state density matrix	
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)
      read(txt,*,err=29) rtol            ! requested convergence for rho
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) Mx              ! Modes to be computed
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) ftol,ftol0,ftol1 ! Lanczos convergence parameters
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) Nrestart,lprint  ! Lanczos max restarts and print level
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) Charg           ! Electronic charge
      read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) iflagE,Ex,Ey,Ez    ! External field
	read(12,'(a)',err=29) txt
!      write(11,'(a)') txt(1:80)   
      read(txt,*,err=29) iflagS,eps,a0,b0,c0   ! Onsager parameters
	
!15    continue      
!      read(*,'(a)',err=29) txt
!       if (txt(1:7).ne.'$ENDPAR') then
!        write(11,'(a)') txt(1:80)
!        goto 15  
!       endif    
! 
      if (mdflag.le.-3) then    
20     continue      
       read(12,'(a)',err=29) txt
       if (txt(1:6).ne.'$COORD') goto 20
          
       read(12,'(a)',err=29) txt
        do while (txt(1:9).ne.'$ENDCOORD')
        Nat=Nat+1
        read(txt,*,err=29) atm(Nat),xx(Nat),yy(Nat),zz(Nat)
        read(12,'(a)',err=29) txt
        enddo
       
!      if (iflagE.gt.1) then
!25    continue      
!      read(12,'(a)',err=29) txt
!      if (txt(1:6).ne.'$FIELD') goto 25
!       
!       i=0   
!       read(*,'(a)',err=29) txt
!       do while (txt(1:9).ne.'$ENDFIELD')
!       i=i+1
!       read(txt,*,err=29) Nfd(i),EEx(i),EEy(i),EEz(i)
!       read(*,'(a)',err=29) txt
!       enddo
!       
!       if (i.ne.Nat) then
!       print "(A,5x,i6,5x,A)",'Field is given only for',i,'atoms'
!       endif
!      endif

!      close(11)
      endif
	close(12)	
	
	endif
210   format(A)
   
!      if (irflag.eq.0) then
! compute MOPAC hamiltonian parameters 
        Np=Charg
      ! FIXME
      ! call MOPAC(ITYPE,ID,keywr,atm,xx,yy,zz,Nat, &
            ! Nb,Nc,Np,Enucl,Atheat,tt,dip0x,dip0y,dip0z, &
            ! nfirst,nmidle,nlast,atmass) 

        Lt=Nb*(Nb+1)/2
        
	  if (iderivfl.eq.0) then ! We are not in analytic derivatives
        do i=1,Lt
         dip0x(i) = dip0x(i)/fDeb
         dip0y(i) = dip0y(i)/fDeb
         dip0z(i) = dip0z(i)/fDeb
        enddo
	  endif

      if (mdflag.lt.0) then      
!  Unformatted write to the disk of dipole matrices
	call wrb_arr(dip0x,Lt,'dip0x.b','w','s')
	call wrb_arr(dip0y,Lt,'dip0y.b','w','s')
	call wrb_arr(dip0z,Lt,'dip0z.b','w','s')
	endif
	     
!      print *, 'F=',Ex,Ey,Ez
!  Add external field contribution to the one-electron matrix
      do i=1,Nat
       if (NFIRST(i).eq.NLAST(i)) then
       Ntype(i)=1
       else
       Ntype(i)=2
       endif
      enddo

        if (iflagE.gt.0) then
         j=0
	 f=0.0
         do i=1,Nat
	    if (Ntype(i).eq.1) then
	     j=j+1
	     n=0
	     f=f-(atm(i)-n)*(dip0x(j*(j-1)/2+j)*Ex+ &
              dip0y(j*(j-1)/2+j)*Ey+dip0z(j*(j-1)/2+j)*Ez)
	    elseif (Ntype(i).eq.2) then 
	     n=2
	     if (atm(i).gt.10) n=10
	     j=j+4
	     do m=j-3,j
	     f=f-(atm(i)-n)/4.0*(dip0x(m*(m-1)/2+m)*Ex+ &
              dip0y(m*(m-1)/2+m)*Ey+dip0z(m*(m-1)/2+m)*Ez)
	     enddo
	    else
	     n=10
	     j=j+9
	     do m=j-8,j
	     f=f-(atm(i)-n)/9.0*(dip0x(m*(m-1)/2+m)*Ex+ &
              dip0y(m*(m-1)/2+m)*Ey+dip0z(m*(m-1)/2+m)*Ez)
	     enddo
	    endif
         enddo
	  do i=1,Lt
	   tt(i)=tt(i)+dip0x(i)*Ex+dip0y(i)*Ey+dip0z(i)*Ez
	  enddo
!  f is a nuclear piece to energy due to field
	Enucl=Enucl+f/feV
        endif

      if (mdflag.lt.0) then
!  Unformatted write to the disk of one and two-electron matrices
	call wrb_arr(tt,Lt,'oneel.b','w','s')
      call wrb_arr2(Nc,WJ,WK,GSS,GSP,GPP,GP2, &
              HSP,GSD,GPD,GDD,'vvv.b','w','s')
      endif
	    
220     format(i6,20x,A) 
230     format(i7,3g15.7,3i6)   
      
!       if (irflag.ne.0) then ! read previous scf results	
!        call wrb_arr (ee,Nb,'fockd.b','r','s')
! 	  call wrb_arr (uu,Nb*Nb,'molor.b','r','s')
!        call wrb_arr2(Nc,WJ,WK,GSS,GSP,GPP,GP2,
!     +        HSP,GSD,GPD,GDD,'vvv.b','r','s')
!       endif

!  form one-electron matrix to check if converged SCF is given:
!	do i=1,Nb*Nb
!	 eta(i)=0.0
!	enddo
!	do i=1,Nb
!	 eta(Nb*(i-1)+i)=ee(i)
!	enddo
!	call transp(Nb,uu,rrwork)
!        call multiple(Nb,eta,rrwork,xi)
!	call multiple(Nb,uu,xi,rrwork)
!
!	call wrb_arr(eta,Lt,'rhogr.b','r','s')
!        call Vxi_pack(eta,xi)
!
!         do i = 1,Nb
!            do j = 1,i
!            tt(i*(i-1)/2+j) = rrwork(Nb*(i-1)+j)- xi(i*(i-1)/2+j)
!            enddo
!         enddo
	        
!	call wrb_arr(xi,Lt,'dip0x.b','r','s')
!	print *, 'dipx'
!	call prmat1(10,xi)
!	call wrb_arr(xi,Lt,'dip0y.b','r','s')
!	print *, 'dipy'
!	call prmat1(10,xi)
!	call wrb_arr(xi,Lt,'dip0z.b','r','s')
!	print *, 'dipz'
!	call prmat1(10,xi)
!	print *, 'Energies'
!	print "(100g13.3)", (ee(i),i=1,Nb)
!	print *, 'fock1'
!	call prmat(Nb,rrwork)
!	call wrb_arr(rrwork,Lt,'focks.b','r','s')
!	print *, 'fock'
!	call prmat1(Nb,rrwork)
!	print *, 'tt1'
!	call prmat1(Nb,tt)
!	call wrb_arr(tt,Lt,'oneel.b','r','s')
!	print *, 'tt'
!	call prmat1(Nb,tt)
!        print *, 'Nfirst, Nmidle,Nlast'
!	do i=1, Nat
!	 print *, NFIRST(i),NMIDLE(i),NLAST(i)
!	enddo
!	print *, 'coulomb'
!	do i=1,Nc
!	 print "(3g12.2)", W(i),WJ(i),WK(i)
!        enddo
      
!     compute nuclear repulsion energy
!      Enucl=0
!      do i=1,Nat
!        if (Ntype(i).eq.1) then
!        k=0
!        elseif (Ntype(i).eq.2) then 
!        k=2
!        if (atm(i).gt.10) k=10
!        else
!        k=10
!        endif
!       do j=i+1,Nat
!        if (Ntype(j).eq.1) then
!        m=0
!        elseif (Ntype(j).eq.2) then 
!        m=2
!        if (atm(i).gt.10) m=10
!        else
!        m=10
!        endif
!      f=(xx(i)-xx(j))**2+(yy(i)-yy(j))**2+(zz(i)-zz(j))**2
!      Enucl=Enucl+(atm(i)-k)*(atm(j)-m)/sqrt(f)
!       enddo
!      enddo
!      Enucl=Enucl*bohrs

! assign the rest of the variables
      Np=Np/2     
      Nh=Nb-Np
      Lt = Nb*(Nb+1)/2
      Mb = Nb*Nb
      M4 = Np*Nh
      M2 = 2*M4
	
      if (first) then
      if (iflagS.eq.2) then
! Compute Onsager ellipsoid parameters if requested    
!
!     initialize variables for integratiion
! 
      delts=.0002                        
      aaa=0
      bbb=0
      ccc=0
      s=0
      i=0
      q=a0
      if (a0.gt.b0) q=b0
      if (b0.gt.c0) q=c0
!      
!
!     perform integral 
!
!        
90       aa=aaa
         bb=bbb
         cc=ccc
         aaa=aaa+1/(s+a0**2)/sqrt((s+a0**2)*(s+b0**2)*(s+c0**2))*delts
         bbb=bbb+1/(s+b0**2)/sqrt((s+a0**2)*(s+b0**2)*(s+c0**2))*delts
         ccc=ccc+1/(s+c0**2)/sqrt((s+a0**2)*(s+b0**2)*(s+c0**2))*delts
         if (s.ge.10*q) delts = .05
         s=s+delts
         i=i+1
         if (i.le.1.or.abs(1-aa/aaa).gt.1E-8.or.abs(1-bb/bbb).gt. &
         1E-8.or.abs(1-cc/ccc).gt.1E-8) go to 90
      Aa1=aaa*a0*b0*c0/2
      Ab=bbb*a0*b0*c0/2
      Ac=ccc*a0*b0*c0/2
!
!     normalize integrals so Aa1+Ab+Ac=1     
!
      tot=Aa1+Ab+Ac
      Aan=Aa1/tot
      Abn=Ab/tot
      Acn=Ac/tot
!
!     calculate tensor elements
! 
      gam1=-3/a0/b0/c0*Aan*(1-Aan)*(eps-1)/(eps+(1-eps)*Aan)*14.405
      gam2=-3/a0/b0/c0*Abn*(1-Abn)*(eps-1)/(eps+(1-eps)*Abn)*14.405
      gam3=-3/a0/b0/c0*Acn*(1-Acn)*(eps-1)/(eps+(1-eps)*Acn)*14.405
!
      endif
	endif
	
	first=.false.
	return
       
29    print *, txt
      stop 'Bad input file'
      
      end
      






