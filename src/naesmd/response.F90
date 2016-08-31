      subroutine response1a (Nb,Nb_M,Mx,M2,e0, &
          dip,eta,xi,rrwork,mu,polar)

      implicit none
      
      integer Nb,Nb_M,Mx,M2,i,j,one
      real*8 xi(Nb_M*Nb_M),eta(Nb_M*Nb_M),dip(Nb_M*(Nb_M+1)/2)
      real*8 rrwork(*)
      real*8 polar,mu(Mx),e0(Mx),fn,f,ddot

      polar=0.0
	one=1
	   
      open (10,file='modes.b',form='unformatted')
      do j = 1,Mx
        read (10) (rrwork(i),i=1,M2)
         call mo2site(rrwork,xi,eta)
         call unpacking(Nb,dip,eta,'s')	   
         mu(j) = ddot(Nb*Nb,xi,one,eta,one)*sqrt(2.0)      
         polar=polar+2*mu(j)**2/e0(j)
       enddo
       close(10)      
           
      return
      end

      subroutine linear ()

      implicit none
      
      include 'parH.par'      
      include 'const.par'
      include 'int.var'
      include 'real.var'
      include 'cart.var'
      include 'hf.var'
      include 'resp.var'
      include 'etaxi.var'      
      include 'rho.var'
      include 'lanczos.var'
      include 'modes.var'
      include 'parvar.var'
      include 'cart.cmn'
      include 'hf.cmn'
      include 'resp.cmn'
      include 'etaxi.cmn'
      include 'rho.cmn'
      include 'lanczos.cmn'
      include 'modes.cmn'
      include 'parvar.cmn' 
      real*8 fx,fy,fz,ft
      

      	call wrb_arr(tt,Lt,'dip0x.b','r','s')
	call copying(Lt,tt,eta)
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mux,alphax)
      	call wrb_arr(xi,Lt,'rho1x.b','w','s')
      	      	
      	call wrb_arr(tt,Lt,'dip0y.b','r','s')
	call copying(Lt,tt,eta)
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,muy,alphay)
      	call wrb_arr(xi,Lt,'rho1y.b','w','s')
              
      	call wrb_arr(tt,Lt,'dip0z.b','r','s')
	call copying(Lt,tt,eta)
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,muz,alphaz)
      	call wrb_arr(xi,Lt,'rho1z.b','w','s')
      	
      	do j=1,Mx
      	ftot(j)=2*e0(j)*(mux(j)**2+muy(j)**2+muz(j)**2)
      	enddo
      	
	open (11,file='freq.dat')
        print *, 'Frequencies (eV) and Oscillator strength (Unitless)'
	print 110 
	do j=1,Mx
	fx=2*e0(j)*mux(j)**2
	fy=2*e0(j)*muy(j)**2
	fz=2*e0(j)*muz(j)**2
	ft=ftot(j)
	print "(i4,5g15.7)",j,e0(j), fx, fy, fz, ft
	write(11,"(i4,4g15.7)") j,e0(j), mux(j), muy(j), muz(j)
        enddo
        close (11)
	
        alpha=(alphaz+alphay+alphax)/3
        f=1.441d-23
        
        print *
        print *, 'Linear polarizabilities, (e*A**2/V) and esu'    	
	print "(5x,'alpha_xx',2g15.7)",alphax,alphax*f
	print "(5x,'alpha_yy',2g15.7)",alphay,alphay*f
	print "(5x,'alpha_zz',2g15.7)",alphaz,alphaz*f
	print *,'    alpha=(alpha_zz+alpha_yy+alpha_xx)/3'
	print "(5x,'alpha   ',2g15.7)",alpha,alpha*f

	return
	
        entry quadratic ()
        
        call wrb_arr(xi,Lt,'rho1x.b','r','s')
	call wrb_arr(eta,Lt,'rhogr.b','r','s')
	call wrb_arr(tt,Lt,'dip0x.b','r','s')
	call quadr2(Nb,tt,eta,xi,rrwork,v0,beta1x)
        call wrb_arr(xi,Lt,'rho2x.b','w','s')
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu2x,beta2x)
        call wrb_arr(tt,Lt,'rho2x.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho2x.b','w','s')
        
        call wrb_arr(xi,Lt,'rho1y.b','r','s')
        call wrb_arr(eta,Lt,'rhogr.b','r','s')
        call wrb_arr(tt,Lt,'dip0y.b','r','s')
        call quadr2(Nb,tt,eta,xi,rrwork,v0,beta1y)
        call wrb_arr(xi,Lt,'rho2y.b','w','s')
        call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu2y,beta2y)
        call wrb_arr(tt,Lt,'rho2y.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho2y.b','w','s')
        
        call wrb_arr(xi,Lt,'rho1z.b','r','s')
        call wrb_arr(eta,Lt,'rhogr.b','r','s')
        call wrb_arr(tt,Lt,'dip0z.b','r','s')
        call quadr2(Nb,tt,eta,xi,rrwork,v0,beta1z)
        call wrb_arr(xi,Lt,'rho2z.b','w','s')
        call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu2z,beta2z)
        call wrb_arr(tt,Lt,'rho2z.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho2z.b','w','s')

        do j=1,Mx
        ftot2(j)=2*e0(j)*(mu2x(j)**2+mu2y(j)**2+mu2z(j)**2)
        enddo

        print *, 'Frequencies (eV) and Oscillator strength (Unitless)'
        print 110 
        do j=1,Mx
        fx=2*e0(j)*mu2x(j)**2
        fy=2*e0(j)*mu2y(j)**2
        fz=2*e0(j)*mu2z(j)**2
        ft=ftot2(j)
        print "(i4,5g15.7)",j,e0(j), fx, fy, fz, ft
        enddo

        betax=beta1x+beta2x
        betay=beta1y+beta2y
        betaz=beta1z+beta2z
        beta=(betaz+betay+betax)/3
        f=4.323d-29

        print *
        print*, beta1x,beta2x
	print*, beta1y,beta2y
	print*, beta1z,beta2z
        print *, 'Quadratic polarizabilities, (e*A**3/V**2) and esu'    	
	print "(5x,'beta_xx',2g15.7)",betax,betax*f
	print "(5x,'beta_yy',2g15.7)",betay,betay*f
	print "(5x,'beta_zz',2g15.7)",betaz,betaz*f
	print *,'    beta=(beta_zz+beta_yy+beta_xx)/3'
	print "(5x,'beta   ',2g15.7)",beta,beta*f
	
	return
	
        entry cubic ()

	call wrb_arr(rrwork,Lt,'rho2x.b','r','s')        
        call wrb_arr(xi,Lt,'rho1x.b','r','s')
	call wrb_arr(eta,Lt,'rhogr.b','r','s')
	call wrb_arr(tt,Lt,'dip0x.b','r','s')
	call cubic3(Nb,tt,eta,xi,rrwork,v0,gamma1x)
        call wrb_arr(xi,Lt,'rho3x.b','w','s')
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu3x,gamma2x)
        call wrb_arr(tt,Lt,'rho3x.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho3x.b','w','s')

	call wrb_arr(rrwork,Lt,'rho2y.b','r','s')        
        call wrb_arr(xi,Lt,'rho1y.b','r','s')
	call wrb_arr(eta,Lt,'rhogr.b','r','s')
	call wrb_arr(tt,Lt,'dip0y.b','r','s')
	call cubic3(Nb,tt,eta,xi,rrwork,v0,gamma1y)
        call wrb_arr(xi,Lt,'rho3y.b','w','s')
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu3y,gamma2y)
        call wrb_arr(tt,Lt,'rho3y.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho3y.b','w','s')
        
	call wrb_arr(rrwork,Lt,'rho2z.b','r','s')        
        call wrb_arr(xi,Lt,'rho1z.b','r','s')
	call wrb_arr(eta,Lt,'rhogr.b','r','s')
	call wrb_arr(tt,Lt,'dip0z.b','r','s')
	call cubic3(Nb,tt,eta,xi,rrwork,v0,gamma1z)
        call wrb_arr(xi,Lt,'rho3z.b','w','s')
      	call response1(Nb,Nb_M,Mx,M2,e0, &
          tt,eta,xi,rrwork,mu3z,gamma2z)
        call wrb_arr(tt,Lt,'rho3z.b','r','d')
        call summing(Lt,xi,tt)
        call wrb_arr(tt,Lt,'rho3z.b','w','s')
	
      	do j=1,Mx
      	ftot3(j)=2*e0(j)*(mu3x(j)**2+mu3y(j)**2+mu3z(j)**2)
      	enddo
		
        print *, 'Frequencies (eV) and Oscillator strength (Unitless)'
	print 110 
	do j=1,Mx
	fx=2*e0(j)*mu3x(j)**2
	fy=2*e0(j)*mu3y(j)**2
	fz=2*e0(j)*mu3z(j)**2
	ft=ftot3(j)
	print "(i4,5g15.7)",j,e0(j), fx, fy, fz, ft
        enddo

	gammax=gamma1x+gamma2x
	gammay=gamma1y+gamma2y
	gammaz=gamma1z+gamma2z
        gamma=(gammaz+gammay+gammax)/3
        f=1.297d-34

        print *
	print *, gamma1x, gamma2x
	print *, gamma1y, gamma2y
	print *, gamma1z, gamma2z
        print *, 'Cubic polarizabilities, (e*A**4/V**3) and esu'    	
	print "(5x,'gamma_xx',2g15.7)",gammax,gammax*f
	print "(5x,'gamma_yy',2g15.7)",gammay,gammay*f
	print "(5x,'gamma_zz',2g15.7)",gammaz,gammaz*f
	print *,'    gamma=(gamma_zz+gamma_yy+gamma_xx)/3'
	print "(5x,'gamma   ',2g15.7)",gamma,gamma*f

110     format(8x,'Omega',12x,'fx',14x,'fy',14x,'fz',10x,'ftotal')
       return
       end
      
      subroutine response1 (Nb,Nb_M,Mx,M2,e0, &
          dip,eta,xi,rrwork,mu,polar)

      implicit none
      
      integer Nb,Nb_M,Mx,M2,i,Lt,Lt_M
      real*8 xi(Nb_M*Nb_M),eta(Nb_M*Nb_M),dip(Nb_M*(Nb_M+1)/2)
      real*8 rrwork(*)
      real*8 polar,mu(Mx),e0(Mx),fn,f

      Lt=Nb*(Nb+1)/2
      Lt_M=Nb_M*(Nb_M+1)/2
      
      do i=1,Lt
      rrwork(i)=eta(i)
      enddo
      
      call response11 (Nb,Nb_M,Mx,M2,e0,dip,rrwork(1), &
          eta,xi,rrwork(Lt_M+1),rrwork(2*Lt_M+1),mu,polar)
          
      do i=1,Lt
      xi(i)=rrwork(Lt_M+i)
      enddo
            
      return
      end


      subroutine response11 (Nb,Nb_M,Mx,M2,e0,dip,dip1, &
          eta1,xi1,drho,mode,mu,polar)

      implicit none
      
      integer Nb,Nb_M,Mx,M2,i,j,k,n,m
      real*8 xi1(Nb,Nb),eta1(Nb,Nb),dip(Nb*(Nb+1)/2)
      real*8 mode(*),drho(Nb*(Nb+1)/2),dip1(Nb*(Nb+1)/2)
      real*8 polar,mu(Mx),e0(Mx),fn,f


         do m = 1,Nb*(Nb+1)/2
            drho(m) = 0.0
         enddo
         
         do m=1,Mx
          mu(m)=0.0
         enddo 

	fn=1/sqrt(2.0)
      open (10,file='modes.b',form='unformatted')
      do j = 1,Mx
        read (10) (mode(i),i=1,M2)
         call mo2site(mode,xi1,eta1)

         do m = 1,Nb
            do n = 1,m-1
             mu(j) = mu(j)+2*(xi1(n,m)+xi1(m,n))*dip1(m*(m-1)/2+n)*fn
            enddo
             mu(j) = mu(j)+2*xi1(m,m)*dip1(m*(m-1)/2+m)*fn
         enddo
         
         do m = 1,Nb
          do n = 1,m
            k=m*(m-1)/2+n
           drho(k) = drho(k)+mu(j)*(xi1(n,m)+xi1(m,n))*fn/e0(j)
          enddo
         enddo         
!         polar=polar+2*mu(j)**2/e0(j)
       enddo
       close(10)

       polar=0.0
       do m = 1,Nb
        do n = 1,m-1
         polar = polar+4*drho(m*(m-1)/2+n)*dip(m*(m-1)/2+n)
        enddo
         polar = polar+2*drho(m*(m-1)/2+m)*dip(m*(m-1)/2+m)
       enddo
              
      return
      end
 
      subroutine quadr2(Nb,tt,eta,xi,rrwork,v0,beta1)

      implicit none
      integer Nb,Lt,Mb,i   
      real*8 tt(*),rrwork(*),v0(*)   
      real*8 xi(*),eta(*),beta1
      
      Lt=Nb*(Nb+1)/2
      Mb=Nb*Nb
      
      do i=1,Lt
      rrwork(i)=eta(i)
      rrwork(i+Lt)=xi(i)
      enddo
    
      call quadr22(Nb,Lt,tt,rrwork(1),rrwork(Lt+1),rrwork(2*Lt+1), &
       eta,xi,v0(1),v0(Mb+1),v0(2*Mb+1),v0(3*Mb+1),beta1)
     
      do i=1,Lt
      xi(i)=rrwork(2*Lt+i)
      eta(i)=rrwork(Lt+i)
      enddo     
  	
      return
      end

      subroutine quadr22(Nb,Lt,dip,rho,ksi1,w2, &
       rhoa,ksi1a,dipa,temp1,temp2,temp3,beta1)
     
      implicit none
      integer Nb,Lt,i,j
      real*8 beta1,dip(Lt),rho(Lt),ksi1(Lt),w2(Lt)
      real*8 rhoa(Nb,Nb),ksi1a(Nb,Nb),temp1(Nb,Nb)
      real*8 temp2(Nb,Nb),temp3(Nb,Nb)
      real*8 dipa(Nb,Nb)

       call unpacking (Nb,rho,rhoa,'s')
       call unpacking (Nb,ksi1,ksi1a,'s')
       call unpacking (Nb,dip,dipa,'s')
	
!     ****** assign W2=(I-2rho)*ksi1^2
      call multiple(Nb,ksi1a,ksi1a,temp1)
      call multiple(Nb,temp1,rhoa,temp2)
      
      do i=1,Nb
         do j=1,Nb
            temp2(i,j)=temp1(i,j)-2.0*temp2(i,j)
         enddo
      enddo	 

      call packing (Nb,temp2,w2,'s')
      call unpacking (Nb,w2,temp2,'s')      
      
!!     ****** calculate beta1=Tr(dip,W2)
      call multiple(Nb,dipa,temp2,temp1)
      call trace(Nb,temp1,beta1)
 
!     *******calculate A=[V(ksi),ksi]+[V(W2),rho]-[dip1,ksi1]  

      call Vxi_pack (ksi1,rho)
      call Vxi_pack (w2,ksi1)
      call unpacking (Nb,rho,temp1,'s') 
      call unpacking (Nb,ksi1,temp2,'s')                 
         
      call commut(Nb,temp1,ksi1a,temp3)
      call commut(Nb,temp2,rhoa,temp1)
      call commut(Nb,dipa,ksi1a,temp2)

      do i=1,Nb
       do j=1,Nb
        temp1(i,j)=temp2(i,j)-(temp3(i,j)+temp1(i,j))
       enddo
      enddo

!     ****** calculate new dipole dip2=[A,rho]
      call commut(Nb,temp1,rhoa,temp2)
      call commut(Nb,temp2,rhoa,temp1)
      call commut(Nb,temp1,rhoa,temp2)
      
      call packing (Nb,temp2,ksi1,'s')
      
      return
      end

      subroutine cubic3(Nb,tt,eta,xi,rrwork,v0,gamma1)

      implicit none
      integer Nb,Lt,Mb,i   
      real*8 tt(*),rrwork(*),v0(*)   
      real*8 xi(*),eta(*),gamma1
      
      Lt=Nb*(Nb+1)/2
      Mb=Nb*Nb
      
      do i=1,Lt
      rrwork(i+Lt)=eta(i)
      rrwork(i+2*Lt)=xi(i)
      enddo
    
      call cubic33(Nb,Lt,tt,rrwork(Lt+1),rrwork(2*Lt+1),rrwork(1), &
       eta,xi,v0(1),v0(Mb+1),v0(2*Mb+1),v0(3*Mb+1),v0(4*Mb+1),gamma1)
     
      do i=1,Lt
      eta(i)=rrwork(2*Lt+i)
      xi(i)=rrwork(i)
      enddo     
  	
      return
      end

      subroutine cubic33(Nb,Lt,dip,rho,ksi1,rho2, &
       rhoa,ksi1a,rho2a,dipa,temp1,temp2,temp3,gamma1)
     
      implicit none
      integer Nb,Lt,i,j
      real*8 gamma1,dip(Lt),rho(Lt),ksi1(Lt),rho2(Lt)
      real*8 rhoa(Nb,Nb),ksi1a(Nb,Nb),rho2a(Nb,Nb)
      real*8 temp1(Nb,Nb),temp2(Nb,Nb),temp3(Nb,Nb)
      real*8 dipa(Nb,Nb)

       call unpacking (Nb,rho,rhoa,'s')
       call unpacking (Nb,ksi1,ksi1a,'s')
       call unpacking (Nb,rho2,rho2a,'s')
       call unpacking (Nb,dip,dipa,'s')
 
!     ****** assign W3=(I-2rho)*(ksi2ksi1+ksi1ksi2)
      call commut(Nb,rho2a,rhoa,temp1)
      call commut(Nb,temp1,rhoa,temp3)
      call multiple(Nb,ksi1a,temp3,temp1)
      call multiple(Nb,temp3,ksi1a,temp2)
      do i=1,Nb
          do j=1,Nb
            temp1(i,j)=temp1(i,j)+temp2(i,j)
          enddo
      enddo
      call multiple(Nb,temp1,rhoa,temp3)
      do i=1,Nb
         do j=1,Nb
            temp3(i,j)=temp1(i,j)-2.0*temp3(i,j)
          enddo
      enddo

!     ****** calculate gamma1=Tr(dip,W3)
      call multiple(Nb,dipa,temp3,temp1)
      call trace(Nb,temp1,gamma1)

!     *******calculate A=[dip,rho2]-
!     [V(rho2),ksi1]+[V(ksi1),rho2]+[V(W3),rho]

      call Vxi_pack (ksi1,rho)
      call Vxi_pack (rho2,ksi1)                  
      call packing (Nb,temp3,rho2,'s')
      call unpacking (Nb,rho,temp2,'s')
      call unpacking (Nb,ksi1,temp1,'s')
      call commut(Nb,dipa,rho2a,temp3)
      call commut(Nb,temp1,ksi1a,dipa)
      call commut(Nb,temp2,rho2a,temp1)
      call Vxi_pack (rho2,rho)
      call unpacking (Nb,rho,temp2,'s')    
      call commut(Nb,temp2,rhoa,rho2a)
      
      do i=1,Nb
         do j=1,Nb
      temp1(i,j)=temp3(i,j)-(temp1(i,j)+dipa(i,j)+rho2a(i,j))
         enddo
      enddo
      
!     ****** calculate new dipole dip3=[A,rho]
      call commut(Nb,temp1,rhoa,temp2)
      call commut(Nb,temp2,rhoa,temp1)
      call commut(Nb,temp1,rhoa,temp2)
      call packing (Nb,temp2,ksi1,'s')

      return
      end
