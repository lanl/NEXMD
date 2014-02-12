! Subroutine for analytic derivatives
      subroutine deriv(Na,Nm,d,t,E0_f,Omega_f)
      use naesmd_constants
      implicit none

      include 'md.par'

      include 'sizes'
	
	real*8 fPi,fbar,ff0,ff1,ff11,ddot
      parameter (fPi = 3.1415926535898d0)
      parameter (fbar = 0.05d0)  ! maximum dE in numerical derivative, eV.
	
	integer Na,Nm,N3,i,j,k,im,one,istate,ip,ih	
      real*8 xx1(Na_M),yy1(Na_M),zz1(Na_M)
	real*8 dxyz(Nm),dxyz1(3,Na_M),xyz(3,Na_M)
!	real*8 dgs1(Nm),desa1(Nm),desb1(Nm)
!	real*8 dgs2(Nm),desa2(Nm),desb2(Nm)
      real*8 r(N3_M),disp0(N3_M)            ! atomic cartesian coordinates
	real*8 Omega,Omega_,EE0,EE0_,dEE0,f,dOmega,t,d,d1,ff
	
	integer mdflag
      include 'parH.par'    
      include 'const.par'
      include 'scf.par' 
      include 'resp.var'	
      include 'etaxi.var'
      include 'envir.var'
      include 'hf.var'      
      include 'rho.var'
      include 'lanczos.var'
	include 'davidson.var'
      include 'modes.var'
      include 'parvar.var'
      include 'coulombM.var'
	include 'derivfl.cmn'
      include 'coulombM.cmn'
      include 'resp.cmn'
      include 'envir.cmn'
      include 'etaxi.cmn'
	include 'davidson.cmn'
      include 'hf.cmn'
      include 'rho.cmn'
      include 'lanczos.cmn'
      include 'modes.cmn'
      include 'parvar.cmn' 
      include 'md.cmn'
      character*20 datetime, keywr*40, machname*36
      common /keywr/ keywr
	real*8 E0_f,Omega_f(Mx_M)
	real*8 rhogr(Nb_M*(Nb_M+1)/2),rhoTZ(Mb_M),rhoLZ(Mb_M),Egr
	integer imodimp(N3_M), jmodimp(N3_M),imimp
      integer itime1,itime11,itime2,itime3,get_time
      integer jm
      double precision ftemp(nmax*3)
      real time11,time12,time13,time23
	save time11,time12,time13,time23
	
      include 'common'

      print*,'kav: in deriv'

      N3 = 3*Na
      one=1
      istate=1
      ff0=0.0
      ff1=1.0
      ff11=-1.0
      Omega=0.0
      Omega_=0.0

! find current cartesian coordinates:
      do j = 1,natom
        xx1(j) = rx(j)
        yy1(j) = ry(j)
        zz1(j) = rz(j) 
      enddo

      do j=1,Nm
         imodimp(j)=1
      enddo

! find analytic derivatives:

! First of all calculate curent Hamiltonian
      mdflag=0
      call input1(keywr,Na,xx1,yy1,zz1,atoms,mdflag)
!	print *, 'tt'
!	call prmat1(Nb,tt)	
	call dcopy(Lt,tt,one,rhoLZ,one) 		
! Now calculate calculate ground state:
      call scf(tt,ee,uu,eta,xi,Nit,rdamp,iflagS,mdflag)
      E0_f=(Enucl+Eelec+Esol)*feV
!	print *, 'E_elect=', EelecE
!	print *, 'Enucl=', Enucl*feV
!	print *, 'E0_f=', E0_f
! Store ground state density (AO) in in rhogr for future
	call dcopy(Lt,tt,one,rhogr,one)      
!	print *, 'Rhogr'
!	call prmat1(Nb,rhogr)
!	call Vxi_pack(rhogr,xi) 
!	print *, 'V(rho)'
!	call prmat1(Nb,xi)
!	call summing(Lt,xi,rhoLZ)
!	print *, 'F'
!	call prmat1(Nb,rhoLZ)
			
!      if (imdtype.gt.0) then ! prepare for excited state gradients
!      istate=imdtype
      if (ihop.gt.0) then ! prepare for excited state gradients
      istate=ihop
! Calculate excited states and excited state density matrix T+Z in rhoTZ
       call LZ(rhogr,rhoTZ,rhoLZ,istate)
       do j=1,istate
        Omega_f(j)=e0(kx(j))
       enddo
	 
! Get appropriate transition density and rhoTZ to AO
       call mo2sitef (Nb,uu,rhoTZ,rrwork(1),rrwork(Mb+1))
!	 print *, 'T+Z in MO'
!       call prmat(Nb,rhoTZ)
!	 print *, 'T+Z in AO'
!       call prmat(Nb,rrwork(1))
       call packing(Nb,rrwork(1),rhoTZ,'s')
!	 print *, 'T+Z in AO'
!       call prmat1(Nb,rhoTZ)

      call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,rrwork(1))
      call mo2sitef (Nb,uu,rrwork(1),rhoLZ,rrwork(Mb+1))

! Some sanity checks
!      call Iminus2rho(Nb,Np,rrwork(1),xi)
!	ff=ddot(Mb,rrwork(1),one,xi,one)
!	print *, 'ff= ', ff
!	call unpacking(Nb,rhogr,xi,'s')
!	call commut(Nb,rhoLZ,xi,eta)
!	call transp1(Nb,rhoLZ)
!	call multiple(Nb,rhoLZ,eta,xi)
!	call prmat(Nb,xi)
!	call trace(Nb,xi,ff)
!	print *, 'ff1= ', ff/2.0	
!	 print *, 'xi in MO', kx(istate)
!       call prmat(Nb,rrwork(1))
!	 print *, 'xi in AO', kx(istate)
!       call prmat(Nb,rhoLZ)
!       call site2mof(Nb,uu,rhoLZ,rrwork(1),rrwork(Mb+1))
!	 call continter(Nb,Np,Nh,M4,rrwork(Mb+1),rrwork(1))
!       call Lxi(rrwork(Mb+1),rrwork(1))
! 	 call expinter(Nb,Np,Nh,M4,rrwork(1),rrwork(Mb+1))	
!      do i=1,Mb
!	 rrwork(Mb+i)=rrwork(Mb+i)/e0(kx(istate))
!	enddo
!	 print *, 'xi in MO', kx(istate),e0(kx(istate))
!       call prmat(Nb,rrwork(Mb+1))
       
      endif

      iderivfl=1	 	
!****************************************
! start the bucle for derivatives calculations
!**********************************************

! Figure type of derivetives to make:

      if (ideriv.ge.2) then  ! Fast GS MOPAC derivatives
! Assign coordinates and zero derivatives
       do j = 1, Na
         xyz(1,j) = rx(j) 
         xyz(2,j) = ry(j) 
         xyz(3,j) = rz(j)
         dxyz1(1,j) = 0.0 
         dxyz1(2,j) = 0.0
         dxyz1(3,j) = 0.0  
        enddo 

! Calculate ground state derivatives E_gr^x=E_nucl^x+E_el^x
!   E_el^x=1/2 Tr((t^x+F^x) rho)  
	
	 call DCART(xyz,dxyz1,rhogr)
! Convert from kcal/A to eV/A
       do j = 3,N3,3
         dxyz(j-2) = -dxyz1(1,j/3)/23.061
         dxyz(j-1) = -dxyz1(2,j/3)/23.061
         dxyz(j) = -dxyz1(3,j/3)/23.061
       enddo 
	 

! Calculate excited state derivatives Omega^x=Tr(F^x rhoTZ)+Tr(V^x(xi) xi^+)
      if (ihop.gt.0) then 
! Term 1: Tr(F^x rhoTZ)
	 call DCART1(xyz,dxyz1,rhogr,rhoTZ)
! Convert from kcal/A to eV/A
       do j = 3,N3,3
         dxyz(j-2) = dxyz(j-2)-dxyz1(1,j/3)/23.061
         dxyz(j-1) = dxyz(j-1)-dxyz1(2,j/3)/23.061
         dxyz(j) = dxyz(j)    -dxyz1(3,j/3)/23.061
       enddo 	 	

! Term 2: Tr(V^x(xi) xi^+) 
! Symmetric part 
	  call packing(Nb,rhoLZ,rrwork(1),'s')
	  call DCART2(xyz,dxyz1,rrwork(1))
! Convert from kcal/A to eV/A
       do j = 3,N3,3
         dxyz(j-2) = dxyz(j-2)-dxyz1(1,j/3)/23.061
         dxyz(j-1) = dxyz(j-1)-dxyz1(2,j/3)/23.061
         dxyz(j) = dxyz(j)	-dxyz1(3,j/3)/23.061
       enddo 

! Antisymmetric part	  
	  call packing(Nb,rhoLZ,rrwork(1),'u')	
	  call DCART2(xyz,dxyz1,rrwork(1))
! Convert from kcal/A to eV/A
       do j = 3,N3,3
         dxyz(j-2) = dxyz(j-2)-dxyz1(1,j/3)/23.061
         dxyz(j-1) = dxyz(j-1)-dxyz1(2,j/3)/23.061
         dxyz(j) = dxyz(j)	-dxyz1(3,j/3)/23.061
       enddo 	 	
	 
      endif
      else ! Standard CEO analytic derivatives
      do im = 1,Nm
         itime1=get_time()
            ff=1.0
10      continue 
            d1=d*ff       
            mdflag=0
! Make geometry increments 
! Find differences Vc^+ - Vc^- and t^+ - t^-

! Increment -
            do j = 3,N3,3
               xx1(j/3) = rx(j/3) - d1*v(j-2,im)
               yy1(j/3) = ry(j/3) - d1*v(j-1,im)
               zz1(j/3) = rz(j/3) - d1*v(j,im)
            enddo
! Find Hamiltonian 	
            call input1(keywr,Na,xx1,yy1,zz1,atoms,mdflag)
! Write hamiltonian
            call hamilderiv(tt,W,WJ,WK,vexp,'w')
! Nuclear energy
            EE0=(Enucl+Esol)*feV

! Increment +
            do j = 3,N3,3
               xx1(j/3) = rx(j/3) + d1*v(j-2,im)
               yy1(j/3) = ry(j/3) + d1*v(j-1,im)
               zz1(j/3) = rz(j/3) + d1*v(j,im)
            enddo
! Find Hamiltonian 	
	    call input1(keywr,Na,xx1,yy1,zz1,atoms,mdflag)
! Read - hamiltonian and substract it from +
! i.e WJ, WK, and tt are dWJ, dWK, and dtt
            call hamilderiv(tt,W,WJ,WK,vexp,'r')
! Difference in nuclear energy
            EE0=(Enucl+Esol)*feV - EE0
            itime11=get_time()
	  	 
!   Ground state derivative is E_gr^x=E_nucl^x+E_el^x
!   E_el^x=1/2 Tr((t^x+F^x) rho) - only derivatives of hamiltonian are needed	
! Electronic piece
            call dcopy(Lt,rhogr,one,eta,one)
            call Vxi_pack(eta,xi) 
!       print *, 'j ', j      
!       print *,'V^x rho'
!       call prmat1(Nb,xi)	
!	     print *, 'xi'
!            call prmat(Nb,xi)
            Eelec = 0
!            do i = 1,Nb
!            do j = 1,i-1
!             Eelec = Eelec+2*eta(i*(i-1)/2+j)*
!     +           (xi(i*(i-1)/2+j)+2*tt(i*(i-1)/2+j))
!            enddo
!             Eelec = Eelec+eta(i*(i-1)/2+i)*
!     +           (xi(i*(i-1)/2+i)+2*tt(i*(i-1)/2+i))
!           enddo 
            call dcopy(Lt,tt,one,rrwork,one)
	    call daxpy(Lt,0.5d0,xi,one,rrwork,one)
            do i = 1,Nb
               Eelec = Eelec-eta(i*(i-1)/2+i)*rrwork(i*(i-1)/2+i)
            enddo 
            Eelec=2*Eelec+4*ddot(Lt,eta,one,rrwork,one)

            dEE0=EE0+Eelec
		
            itime2=get_time()
	    time13=time13+real(itime2-itime11)/100
	    time11=time11+real(itime2-itime1)/100
!   Excited state derivatives Omega^x=Tr(F^x rhoTZ)+Tr(V^x(xi) xi^+)
!   Only derivatives of hamiltonian are needed	
!            if (imdtype.gt.0) then 
            if (ihop.gt.0) then 
! First piece Tr(F^x rhoTZ)
! Note V(rho) is in xi and tt is given by input1 above
                Omega = 0
!               do i = 1,Nb
!                 do j = 1,i-1
!                    Omega = Omega+2*rhoTZ(i*(i-1)/2+j)*
!     +           (xi(i*(i-1)/2+j)+tt(i*(i-1)/2+j))
!                 enddo
!                 Omega = Omega+rhoTZ(i*(i-1)/2+i)*
!     +           (xi(i*(i-1)/2+i)+tt(i*(i-1)/2+i))
!               enddo            
                call daxpy(Lt,1.0d0,xi,one,tt,one)
                do i = 1,Nb
                   Omega = Omega-rhoTZ(i*(i-1)/2+i)*tt(i*(i-1)/2+i)
                enddo 
                Omega=Omega+2*ddot(Lt,rhoTZ,one,tt,one)
!	 print *, 'f= ', Omega/(2.0*d1)      
!        do j = 3,N3,3
!           print *, xx1(j/3)-rx(j/3),yy1(j/3)-ry(j/3),zz1(j/3)-rz(j/3)
!        enddo
!             print *,'F^x'
!             call prmat1(Nb,tt)
		    
                itime11=get_time()
	   	   
! Second piece Tr(V^x(xi) xi^+)
                call Vxi (rhoLZ,rrwork(1))
!		    print *, 'f= ', Omega  
!		    print *, 'f1= ', ddot(Mb,rhoLZ,one,rrwork(1),one)	
                Omega = Omega+ddot(Mb,rhoLZ,one,rrwork(1),one)
		    		
	        dOmega = Omega
                itime1=get_time()
                time12=time12+real(itime1-itime2)/100	
	        time23=time23+real(itime1-itime11)/100		   
	    else
	        dOmega = 0.0
	    endif      

!        print *, 'im= ', im,dEE0,dOmega

            if (abs(dEE0).gt.fbar.or.abs(dOmega).gt.fbar) then
                write(13,302) 'WARNING: dE is too big: ',im,dEE0,dOmega
	        ff=ff*0.5
	        goto 10
	    endif

            if (ihop.gt.0) then ! Excited state 
                dxyz(im)=-(dEE0+dOmega)/d1/2.0         !   Derrivative (eV/A)
             else	                             ! Ground state 	
               dxyz(im)=-dEE0/d1/2.0                  !   Derrivative (eV/A)  
            endif
!	    write(13,301) im,a(im),dEE0,dOmega
      enddo
	endif
	
       do im=1,Nm
	   if (imodimp(im).eq.1) then
! modified by Seba**********************
	    ftemp(im)= dxyz(im) 
!****************************************
         endif
	 enddo
!*******************************************************
! end of the bucle for derivatives calculation
!*******************************************************
! define the current cartesian forces and accelerations
      do jm = 3,N3,3
         fxmdqt(jm/3) = ftemp(jm-2)
         fymdqt(jm/3) = ftemp(jm-1)
         fzmdqt(jm/3) = ftemp(jm)
      enddo
      do jm=1,natom
         fxmdqt(jm)=fxmdqt(jm)*convl/feVmdqt
         fymdqt(jm)=fymdqt(jm)*convl/feVmdqt
         fzmdqt(jm)=fzmdqt(jm)*convl/feVmdqt
      enddo
      do jm=1,natom
         ax(jm) = fxmdqt(jm) / massmdqt(jm)
         ay(jm) = fymdqt(jm) / massmdqt(jm)
         az(jm) = fzmdqt(jm) / massmdqt(jm)
      enddo

	  close (13)
!****************************************
	 
!      print *, 'dxyz'	
!      do j = 3,N3,3
!	  print "(3e12.3)", dxyz(j-2),dxyz(j-1),dxyz(j) 
!       enddo	 
	   
	iderivfl=0	
	  
! find current cartesian coordinates again and calculate current Hamiltonian:
      do j = 3,N3,3
        xx1(j/3) = r(j-2)
        yy1(j/3) = r(j-1)
        zz1(j/3) = r(j) 
       enddo
      do j = 1,natom
        xx1(j) = rx(j)
        yy1(j) = ry(j)
        zz1(j) = rz(j) 
      enddo
      mdflag=0
      call input1(keywr,Na,xx1,yy1,zz1,atoms,mdflag)

301   format('diff:',i5,4g17.9)
302   format(a,i3,2g13.5)
      
	print *, 'Time GS ', time11,' Vxi ',time13
	print *, 'Time ES ', time12,' Vxi ',time23

      return
	
      end

! Subroutine for to solve for relaxed excited state density
      subroutine LZ(rhogr,rhoT,rhoLZ,istate)
      implicit none


	real*8 f,t,d,d1,ff,ff0,ff1,ff11	
	integer i,ii,j,k,im,one,istate,mdflag,ip,ih	
	
      include 'parH.par'    
      include 'const.par'
      include 'scf.par' 
      include 'resp.var'	
      include 'etaxi.var'
      include 'envir.var'
      include 'hf.var'      
      include 'rho.var'
      include 'lanczos.var'
      include 'modes.var'
      include 'parvar.var'
      include 'coulombM.var'
      include 'coulombM.cmn'
      include 'resp.cmn'
      include 'envir.cmn'
      include 'etaxi.cmn'
      include 'hf.cmn'
      include 'rho.cmn'
      include 'lanczos.cmn'
      include 'modes.cmn'
      include 'parvar.cmn' 

	real*8 rhogr(Nb_M*(Nb_M+1)/2),rhoT(Mb_M),rhoLZ(Mb_M)
! temporary, later delete
!	real*8 tmp1(Mb_M),tmp2(Mb_M),tmp3(Mb_M),tmp4(Mb_M),tmp5(Mb_M)
	
	one=1
	mdflag=0
	ff0=0.0
	ff1=1.0
	ff11=-1.0
      
! First calculate excited states energies
!	 Mx=istate+1   ! the safest way to ensure state in question
       Mx=istate
       call davidson(mdflag)
	 call rrdpsort (e0,Mx,kx,2)
!	 do i=1,Mx
!	  print *, i, kx(i), e0(kx(i))
!	  call getmodef(M2_M,Mx_M,Np,Nh,kx(i),v0,rrwork)
!	  call prmat(Nb,rrwork)
!	 enddo
! Calculate unrelaxed density of a given excited state
! T=[[xi^+ rho] xi]= (I-2rho)(xi^+ xi+xi xi^+)
       call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,rrwork)
	 call dgemm('N','T',Nb,Nb,Nb,ff1,rrwork,Nb,rrwork,Nb,ff0,eta,Nb)
       call dgemm('T','N',Nb,Nb,Nb,ff1,rrwork,Nb,rrwork,Nb,ff1,eta,Nb)
       call Iminus2rho(Nb,Np,eta,rhoT) 
!       print *, 'delta rho MO'
!	 call prmat(Nb,rhoT)
! Store unrelaxed excited state density (MO) in in rhogT for future
	 
! Start putting together left hand side of LZ-equation
! Form [V(rhoT), rho]=(I-2rho)V(rhoT)  
       call mo2sitef (Nb,uu,rhoT,eta,rrwork(Mb+1))
       
!	 print *, 'delta rho AO'
!	 call prmat(Nb,eta)
	 
	 call Vxi (eta,rrwork(1))
! Commutator version
!	 call unpacking(Nb,rhogr,tmp1,'s')
!       call site2mof (Nb,uu,tmp1,tmp2,rrwork(Mb+1))
!	 print *, 'rho(AO)'
!	 call prmat(Nb,tmp1)
!	 print *, 'rho(MO)'
!	 call prmat(Nb,tmp2)
!       call commut (Nb,rrwork(1),tmp1,tmp2)
!	 call site2mof (Nb,uu,tmp2,tmp1,rrwork(Mb+1))
!	 print *, '[V(rhoT), rho](MO)'
!	 call prmat(Nb,tmp1)

! Normal version	 
	 call site2mof (Nb,uu,rrwork(1),xi,rrwork(Mb+1))
	 call Iminus2rho(Nb,Np,xi,rhoLZ) 
	 call project(Nb,Np,Nh,rhoLZ)
!       print *, 'V(delta rho AO)'
!	 call prmat(Nb,rrwork(1))
!       print *, 'V(delta rho MO)'
!	 call prmat(Nb,xi)
!	 print *, '[V(rhoT), rho](MO)'
!	 call prmat(Nb,rhoLZ)	 

! Form [[xi^+, rho], V(xi)], rho] +cc= (I-2rho)[(I-2rho)xi^+, V(xi)]+cc
! Normal version
      call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,xi)
      call mo2sitef (Nb,uu,xi,rrwork(1),rrwork(Mb+1))
	call Vxi (rrwork(1),eta)
      call site2mof (Nb,uu,eta,rrwork(1),rrwork(Mb+1))
!       print *, 'V(xi), MO'
!	 call prmat(Nb,rrwork(1))
	call transp1(Nb,xi)
      call Iminus2rho(Nb,Np,xi,eta)
      call dgemm('N','N',Nb,Nb,Nb,ff1,eta,Nb,rrwork(1),Nb,ff0,xi,Nb)
      call dgemm('N','N',Nb,Nb,Nb,ff11,rrwork(1),Nb,eta,Nb,ff1,xi,Nb)   
      call symmetr(Nb,xi)	 
	call Iminus2rho(Nb,Np,xi,eta)
      call project(Nb,Np,Nh,eta)
!	 print *, '(I-2rho)V(rhoT) (MO)'
!	 call prmat(Nb,rhoLZ)	 
!	 print *, '(I-2rho)[(I-2rho)xi^+, V(xi)] (MO, symm)' 
!	 call prmat(Nb,eta)
	 	
! Commutator version
!	call unpacking(Nb,rhogr,tmp1,'s')
!      call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,xi)
!      call mo2sitef (Nb,uu,xi,rrwork(1),rrwork(Mb+1))
!	call Vxi (rrwork(1),eta)
!      call transp1(Nb,rrwork(1))
!	 print *, 'xi, AO' 
!	 call prmat(Nb,rrwork(1))
!       print *, 'V(xi), AO'
!	 call prmat(Nb,eta)	 
!	call commut(Nb,rrwork(1),tmp1,tmp2)
!	call commut(Nb,tmp2,eta,xi)
!      call commut(Nb,xi,tmp1,tmp2)
!	call site2mof (Nb,uu,tmp2,rrwork(1),rrwork(Mb+1))
!	 print *, '[[[xi^+, rho], V(xi)],rho]' 
!	 call prmat(Nb,rrwork(1))

! Finally put together left hand side of LZ-equation
! 2 is coming from the complex comjugate for eta
      do i=1,Mb
	 rhoLZ(i)= rhoLZ(i) - 2*eta(i)
	enddo
!	 print *, 'LZ (MO)' 
!	 call prmat(Nb,rhoLZ)

! Now solve for Z equation LZ=rhoLZ

! Start with 0-order
!      call continter(Nb,Np,Nh,M4,tmp1,rhoLZ)
!	f=1.0
! Start loop here     
!      do im=1,2
!	 print *
!	 print *, 'It ', im
!	 print *, 'Delta LZ',f
!	 call expinter(Nb,Np,Nh,M4,xi,tmp2)
!	 call prmat(Nb,tmp2)
!       i = 0
!       do ip = 1,Np
!         do ih = Np+1,Nb
!            i = i + 1
!            ff = 1/(ee(ih) - ee(ip))
!            tmp2(i) = ff*tmp1(i)
!		tmp2(i+M4) = -ff*tmp1(i+M4)
!         enddo
!       enddo
!	 call Lxi (tmp2,xi)   
!	 call substract(M2,xi,tmp1)
! Check for convergency
!	 f=0.0
!	 do j=1,M2
!	  f=f+abs(tmp1(j))
!	 enddo
! If converged, exit       
!	enddo

! Start with 0-order
      call continter(Nb,Np,Nh,M4,rrwork(3*M2+1),rhoLZ)
	call dcopy(M2,rrwork(3*M2+1),one,rrwork(1),one)	
	f=1.0
! Start loop here     
      do im=1,Ni_M
!	 print *, 'It ', im,' Delta LZ ',f
       i = 0
       do ip = 1,Np
         do ih = Np+1,Nb
            i = i + 1
            ff = 1/(ee(ih) - ee(ip))
            eta(i) = ff*rrwork(3*M2+i)
		eta(i+M4) = -ff*rrwork(3*M2+M4+i)
         enddo
       enddo
	 if (im.eq.1) then
	  call dcopy(M2,eta,one,rrwork(M2+1),one)
	 else
	  call summing(M2,eta,rrwork(M2+1))
	 endif 	 
	 call Lxi(rrwork(M2+1),rrwork(2*M2+1))   
! Check for convergency
	 f=0.0
	 do j=1,M2
	  rrwork(3*M2+j)=rrwork(j)-rrwork(2*M2+j)
	  f=f+abs(rrwork(3*M2+j)**2)
	 enddo
! If converged, exit
       if (f.lt.ftol0**2) then ! converged
!	  print *, 'rhoLZ'
!	  call expinter(Nb,Np,Nh,M4,rrwork(1),eta)
!	  call prmat(Nb,eta)
!	  print *, 'rhoLZ obtained'
! 	  call expinter(Nb,Np,Nh,M4,rrwork(2*M2+1),eta)
!	  call prmat(Nb,eta)
!	  print *, 'Z -matrix obtained'     
	  call expinter(Nb,Np,Nh,M4,rrwork(M2+1),eta)
!	  call prmat(Nb,eta)
!	  print *, 'T -matrix obtained'     
!	  call prmat(Nb,rhoT)	
!	  call summing(Mb,eta,rhoT)
         do j=1,Mb
	    rhoT(j)=rhoT(j)+eta(j)
	   enddo
!	  print *, 'T+Z -matrix obtained'     
!	  call prmat(Nb,rhoT)	
	  goto 20
	 endif 
! Quit if convergence is not achieved
	 if (im.eq.Ni_M) then 
	  print *, 'Eq. for Z: LZ=rhoLZ did not converged'
	  print *, 'Achieved convergence= ', f
	  stop
	 endif	  	  
	enddo

20    continue

! Now calculate Z as Z= sum_i V_-11i / omega_i xi_i^+
! Note, possible only for small molecules where Mx=Np*Nh 
! i.e. all staes are needed
!       call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,xi)
!       call getmodef(M2_M,Mx_M,Np,Nh,kx(istate),v0,eta)	
!	 call transp1(Nb,eta)
!	 call clearing(Mb,tmp5)
!	 do i=-Mx,Mx
!	 ff = 0.0
!	  if (i.ne.0) then  
!	   ii=abs(i)
!	    call getmodef(M2_M,Mx_M,Np,Nh,kx(ii),v0,tmp1)	      
!	    if (i.gt.0) call transp1(Nb,tmp1)
!	    call vv3(Nb,Np,uu,xi,eta,tmp1,tmp2,tmp3,tmp4,f)
!	    ff=ff+f/e0(ii)
!	    call vv3(Nb,Np,uu,tmp1,eta,xi,tmp2,tmp3,tmp4,f)
!	    ff=ff+f/e0(ii)
!	    call vv3(Nb,Np,uu,xi,tmp1,eta,tmp2,tmp3,tmp4,f)
!	    ff=ff+f/e0(ii)
!	   print *, 'ff= ',ff
!	   do j=1,Mb
!	    tmp5(j)=tmp5(j)+ff*tmp1(j)
!	   enddo
!	  endif
!	 enddo
!	  print *, 'Z -matrix another way'     
!	  call prmat(Nb,tmp5)		 

      return
      end

! Subroutine to calculate Coulomb coupling V_ijk
! See right above
      subroutine vv3(Nb,Np,uu,xi1,xi2,xi3,tmp1,tmp2,tmp3,f)
      implicit none

	real*8 f,ff,f0,f1,f11	
	integer i,j,k,Nb,Np
	real*8 xi1(Nb,Nb),xi2(Nb,Nb),xi3(Nb,Nb),uu(Nb,Nb) 
      real*8 tmp1(Nb,Nb),tmp2(Nb,Nb),tmp3(Nb,Nb) 

	f0=0.0
	f1=1.0
	f11=-1.0	
       call mo2sitef (Nb,uu,xi3,tmp2,tmp3)
	 call Vxi (tmp2,tmp1)
       call site2mof (Nb,uu,tmp1,tmp2,tmp3)	
	 
	 call dgemm('N','N',Nb,Nb,Nb,f1,xi1,Nb,xi2,Nb,f0,tmp3,Nb)
	 call dgemm('N','N',Nb,Nb,Nb,f1,xi2,Nb,xi1,Nb,f1,tmp3,Nb)
	 call Iminus2rho(Nb,Np,tmp3,tmp1)
	 
	 call dgemm('N','N',Nb,Nb,Nb,f1,tmp1,Nb,tmp2,Nb,f0,tmp3,Nb)  
	 call trace (Nb,tmp3,f)
	 f=f/2.0
	 	
      return
      end

! Subroutine for derivatives of Coulomb and hopping matrices
      subroutine hamilderiv(tt,W,WJ,WK,vexp,rw)
      implicit none
 
      include 'parH.par'   
      include 'const.par'
      include 'int.var'      
      include 'parvar.var'            
      include 'parvar.cmn'

	character rw*1      
	real*8 tt(Lt),W(Nc),WJ(Nc),WK(Nc),vexp(2*M2_M*Nd_M),f
	integer one

	one=1

      if (rw.eq.'w') then   ! write Coulomb and hopping
	 i=1
	 call dcopy(Lt,tt,one,vexp(i),one)
	 i=i+Lt
	 call dcopy(Nc,W,one,vexp(i),one)
	 i=i+Nc
	 call dcopy(Nc,WK,one,vexp(i),one)
	endif

      if (rw.eq.'r') then   ! read Coulomb and hopping
	 f=-1.0
	 i=1
	 call daxpy(Lt,f,vexp(i),one,tt,one)
	 i=i+Lt
	 call daxpy(Nc,f,vexp(i),one,W,one)
	 i=i+Nc
	 call daxpy(Nc,f,vexp(i),one,WK,one)
!	 print *, rw
!	 call prmat0(10,200,WJ)
!       do j=1,Nc
!	  print *, j,Nc,W(j),vexp(Lt+j)
!	 enddo
	 endif	

	return
	end
