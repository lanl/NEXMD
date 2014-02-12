      subroutine Lxi (u1,v1)

      include 'parH.par'
      include 'int.var'
      include 'hf.var'
      include 'etaxi.var'
      include 'parvar.var'
      include 'modes.var'
      include 'modes.cmn'
      include 'parvar.cmn'
      include 'hf.cmn'
      include 'etaxi.cmn'

      real*8 u1(M2_M),v1(M2_M)
	real*8 temp1(M2_M),temp2(M2_M)
      real*8 f,f1,f2,fs1,ddot
	integer one
      parameter (one = 1)
	
	fs1=0

! multiply:
      call mo2site (u1,xi,eta)
      call Vxi (xi,eta)
      call site2mo (xi,eta,v1)
      i = 0
      do p = 1,Np
         do h = Np+1,Nb
            i = i + 1
            f = ee(h) - ee(p)
            v1(i) = v1(i) + f*u1(i)
            v1(i+M4) = - (v1(i+M4) + f*u1(i+M4))
         enddo
      enddo

! shift:
      do i=1,2*M4
	 temp1(i)= v0(i,1)
	 temp2(i)= u1(i)
	enddo
	if (Mj.gt.0) then
	 fs1=fs+e0(Mj)
       do j = 1,Mj
	   f1 = fs1*(ddot(M4,v0(1,j),one,u1(1),one)-ddot(M4,v0(M4+1,j),one,u1(M4+1),one))
	   f2 = fs1*(ddot(M4,v0(M4+1,j),one,u1(1),one)-ddot(M4,v0(1,j),one,u1(M4+1),one))
         call daxpy(M4,f1,v0(1,j),one,v1(1),one)
	   call daxpy(M4,f2,v0(1+M4,j),one,v1(1),one)
         call daxpy(M4,f1,v0(1+M4,j),one,v1(1+M4),one)
	   call daxpy(M4,f2,v0(1,j),one,v1(1+M4),one)
!! 	    f1 = 0
!! 	    f2 = 0
!! 	    do i = 1,M4
!! 		 f1 = f1 + v0(i,j)*u1(i) - v0(i+M4,j)*u1(i+M4)
!! 		 f2 = f2 + v0(i+M4,j)*u1(i) - v0(i,j)*u1(i+M4)
!! 	    enddo
!! 	    f1 = fs*f1
!! 	    f2 = fs*f2
!! 	    do i = 1,M4
!! 		 v1(i) = v1(i) + f1*v0(i,j) + f2*v0(i+M4,j)
!! 		 v1(i+M4) = v1(i+M4) + f1*v0(i+M4,j) + f2*v0(i,j)
!! 	    enddo
       enddo
	endif
	return
      end

      subroutine site2mo (zz,xi,v1)
      
      implicit none

      include 'parH.par'
      include 'int.var'
      include 'rho.var'
      include 'parvar.var'
      include 'parvar.cmn'
      include 'rho.cmn'
      
      real*8 f,fu,fv,f0,f1
      parameter (f0 = 0.0)
      parameter (f1 = 1.0)

      real*8 xi(Nb,Nb),v1(M2)

      real*8 zz(Nb,Nb)
      
      call dgemm ('T','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
      call dgemm ('T','N',Nh,Np,Nb,f1,uu(Nb*Np+1),Nb,zz,Nb,f0,v1,Nh)
      call dgemm ('T','N',Nh,Np,Nb,f1,zz(1,Np+1),Nb,uu,Nb,f0, &
       v1(M4+1),Nh)      

!      do l = 1,Nb
!         do j = 1,Nb
!            f = 0
!            do i = 1,Nb
!               f = f + uu(i,l)*xi(i,j)
!            enddo
!            zz(j,l) = f
!         enddo
!      enddo
!      i = 0
!      do p = 1,Np
!         do h = Np+1,Nb
!            i = i + 1
!            fu = 0
!            fv = 0
!            do j = 1,Nb
!               fu = fu + uu(j,h)*zz(j,p)
!            enddo
!            do j = 1,Nb
!               fv = fv + uu(j,p)*zz(j,h)
!            enddo
!            v1(i) = fu
!            v1(i+M4) = fv
!         enddo
!      enddo

      return


      entry mo2site (v1,xi,zz)

      call dgemm ('T','T',Np,Nb,Nh,f1,v1(1),Nh,uu(Nb*Np+1),Nb, &
       f0,zz,Nb)
      call dgemm ('N','T',Nh,Nb,Np,f1,v1(M4+1),Nh,uu,Nb, &
       f0,zz(Np+1,1),Nb)
      call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,xi,Nb)

!      do j = 1,Nb
!         i = 0
!         do p = 1,Np
!            fu = 0
!            fv = 0
!            do h = Np+1,Nb
!               i = i + 1
!               fu = fu + ww(h,j)*v1(i)
!               fv = fv + ww(h,j)*v1(i+M4)
!            enddo
!            zz(p,j) = fu
!            zz(p+Np,j) = fv
!         enddo
!      enddo
!      do j = 1,Nb
!         do i = 1,Nb
!            f = 0
!            do p = 1,Np
!               f = f + ww(p,i)*zz(p,j) + ww(p,j)*zz(p+Np,i)
!            enddo
!            xi(i,j) = f
!         enddo
!      enddo
      
      return
      end

      subroutine site2moph (Nb,Np,Nh,M4,uu,xi,v1,zz)
      
      implicit none
      integer Nb,Np,Nh,M4
      real*8 uu(Nb*Nb),xi(Nb,Nb),zz(Nb,Nb),v1(2*M4)
      real*8 f0,f1
      
      parameter (f0 = 0.0)
      parameter (f1 = 1.0)
        
      call dgemm ('T','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
      call dgemm ('T','N',Nh,Np,Nb,f1,uu(Nb*Np+1),Nb,zz,Nb,f0,v1,Nh)
      call dgemm ('T','N',Nh,Np,Nb,f1,zz(1,Np+1),Nb,uu,Nb,f0, &
       v1(M4+1),Nh)      

      return
        
        entry mo2siteph (Nb,Np,Nh,M4,uu,v1,xi,zz)

      call dgemm ('T','T',Np,Nb,Nh,f1,v1(1),Nh,uu(Nb*Np+1),Nb, &
       f0,zz,Nb)
      call dgemm ('N','T',Nh,Nb,Np,f1,v1(M4+1),Nh,uu,Nb, &
       f0,zz(Np+1,1),Nb)
      call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,xi,Nb)

      return
        end

      subroutine site2mof (Nb,uu,xi,eta,zz)
      
      implicit none
      integer Nb
      real*8 uu(Nb,Nb),xi(Nb,Nb),zz(Nb,Nb),eta(Nb,Nb)
      real*8 f0,f1
      
      parameter (f0 = 0.0)
      parameter (f1 = 1.0)

      call dgemm ('N','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
      call dgemm ('T','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)

      return
        
        entry mo2sitef (Nb,uu,xi,eta,zz)

      call dgemm ('N','T',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
      call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)

      return
        end




