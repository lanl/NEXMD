#include "dprec.fh"

!---------------------------------------------------
!       This subroutine computes (I-2rho)*xi
!       and xi*(I-2rho) for real and complex matrix
!---------------------------------------------------

subroutine Iminus2rho(Nb,Np,temp1,temp2) 
          
	implicit none

	integer Nb,Np,j,i
	_REAL_  temp1(Nb,Nb),temp2(Nb,Nb)
        
! Particle-particle and Hole-particle
	do i=1,Np
		do j=1,Nb
			temp2(i,j)=-temp1(i,j)
		enddo
	enddo
	do i=Np+1,Nb
		do j=1,Nb
			temp2(i,j)=temp1(i,j)
		enddo
	enddo
         
	return
        
entry Iminus2rho1(Nb,Np,temp1,temp2)

	do i=1,Np
		do j=1,Nb
			temp2(j,i)=-temp1(j,i)
		enddo
	enddo
	do i=Np+1,Nb
		do j=1,Nb
			temp2(j,i)=temp1(j,i)
		enddo
	enddo   
        
	return
end subroutine Iminus2rho

subroutine Iminus2rhoz(Nb,Np,ztmp1,ztmp2) 
        
	implicit none

	integer Nb,Np,j,i
	complex*16  ztmp1(Nb,Nb),ztmp2(Nb,Nb)
        
! Particle-particle and Hole-particle
	do i=1,Np
		do j=1,Nb
			ztmp2(i,j)=-ztmp1(i,j)
		enddo
	enddo
	do i=Np+1,Nb
		do j=1,Nb
			ztmp2(i,j)=ztmp1(i,j)
		enddo
	enddo
        
	return
        
entry Iminus2rho1z(Nb,Np,ztmp1,ztmp2)

	do i=1,Np
		do j=1,Nb
			ztmp2(j,i)=-ztmp1(j,i)
		enddo
	enddo
	do i=Np+1,Nb
		do j=1,Nb
			ztmp2(j,i)=ztmp1(j,i)
		enddo
	enddo   
	return
end subroutine Iminus2rhoz

subroutine getmode(M2_M,Mx_M,M2,i,v0,tmp)
	implicit none
	integer M2_M,Mx_M,M2,i,j
	_REAL_ v0(M2_M,Mx_M),tmp(M2)

	do j=1,M2
	  tmp(j)=v0(j,i)
	enddo
	
           return
end subroutine getmode
       
subroutine getmodez(M2_M,Mx_M,M2,i,v0,tmp)
	implicit none
	integer M2_M,Mx_M,M2,i,j
	_REAL_ v0(M2_M,Mx_M)
	complex*16 tmp(M2)

	do j=1,M2
		tmp(j)=dcmplx(v0(j,i),0.0d0) 
	enddo
	
	return
end subroutine getmodez

subroutine getmodef(M2_M,Mx_M,Np,Nh,i,v0,temp)
      	implicit none
      	integer M2_M,Mx_M,i,j,k,jj,Np,Nh
       	_REAL_ v0(M2_M,Mx_M),temp(Np+Nh,Np+Nh)

      	do j=1,Np
         	do k=1,Np
          		temp(j,k)=0.0
         	enddo
       	enddo
      	do j=Np+1,Nh+Np
         	do k=Np+1,Nh+Np
          		temp(j,k)=0.0
         	enddo
      	enddo
         
       	jj=0    
      	do j=1,Np
         	do k=1,Nh
          		jj=jj+1
          		temp(j,Np+k)=v0(jj,i)
        		temp(Np+k,j)=v0(jj+Np*Nh,i)
         	enddo
       	enddo   
      	return
end subroutine getmodef

!---------------------------------------------------
!       This subrotine prints a square matrix
!---------------------------------------------------
subroutine prmat(n,d)

	implicit none
	integer n,i,j
	_REAL_ d(n,n)

	do i=1,n
		write(6,910) (d(i,j),j=1,n)
	enddo

	910  format(' ', 40g11.3)
! 910  format(' ', 40ES14.5)
! 910  format(' ', 20f11.7)
	return
end subroutine prmat

!---------------------------------------------------
!       This subrotine prints an nxm matrix
!---------------------------------------------------
subroutine prmat0(n,m,d)

      implicit none
      integer n,m,i,j
      _REAL_ d(n,m)
       do i=1,m
         write(6,910) (d(i,j),j=1,n)
       enddo
 910  format(' ', 10g11.3)
      return
end subroutine prmat0

!---------------------------------------------------
!       This subrotine print triangular matrix
!---------------------------------------------------
subroutine prmat1(n,d)
             
      implicit none
      integer n,i,j
      _REAL_ d(*)
      do i=1,n
         write(6,910) (d(i*(i-1)/2+j),j=1,i)
      enddo 
!910  format(' ', 10g10.3)
 910  format(' ', 20f11.7)
      return
end subroutine prmat1

!---------------------------------------------------
!       This subroutine project matrix to interband
!       space
!---------------------------------------------------

      subroutine project(Nb,Np,Nh,tmp)

      implicit none
      integer Nb,Np,Nh,M4,j,i,k
      _REAL_  tmp(Nb,Nb)
! Particle-particle
      do i=1,Np
         do j=1,Np
          tmp(i,j)=0.0
         enddo
        enddo
! Hole-hole
      do i=Np+1,Nb
         do j=Np+1,Nb
          tmp(i,j)=0.0
         enddo
      enddo

        return
end subroutine project
!---------------------------------------------------
!       This subrotine transpose matrix
!---------------------------------------------------

subroutine transp(n,d1,d2)
	implicit none
	integer n,i,j
	_REAL_ d1(n,n),d2(n,n)

	do i=1,n
		do j=1,n
			d2(j,i)=d1(i,j)
		enddo
	enddo
        
	return
end subroutine transp

subroutine transp1(n,d1)
	implicit none
	integer n,i,j
	_REAL_ d1(n,n),f

	do i=1,n
		do j=1,i
			f=d1(i,j)
			d1(i,j)=d1(j,i)
			d1(j,i)=f
		enddo
	enddo
        
	return
end subroutine transp1

subroutine transp1z(n,d1)
	implicit none
	integer n,i,j
	complex*16 d1(n,n),f

	do i=1,n
		do j=1,i
			f=d1(i,j)
			d1(i,j)=d1(j,i)
			d1(j,i)=f
		enddo
	enddo
        
	return
end subroutine transp1z

!---------------------------------------------------
!       This subroutine expand interband matrix
!       from 2*M4 to full Mb and back
!---------------------------------------------------
subroutine expinter(Nb,Np,Nh,M4,tmp1,tmp2)

	implicit none
	integer Nb,Np,Nh,M4,j,i,k
	_REAL_  tmp1(2*M4),tmp2(Nb,Nb)

! Particle-particle and Hole-hole
	call clearing(Nb*Nb,tmp2)
!c Particle-hole
	do i=Np+1,Nb
		do j=1,Np
			tmp2(j,i)=tmp1(Nh*(j-1)+i-Np)
		enddo
	enddo
! Hole-particle
	do i=1,Np
		do j=Np+1,Nb
			tmp2(j,i)=tmp1(Nh*(i-1)+j-Np+M4)
		enddo
	enddo

	return

entry continter(Nb,Np,Nh,M4,tmp1,tmp2)

! Particle-hole
	do i=Np+1,Nb
		do j=1,Np
			tmp1(Nh*(j-1)+i-Np)=tmp2(j,i)
		enddo
	enddo
! Hole-particle
	do i=1,Np
		do j=Np+1,Nb
			tmp1(Nh*(i-1)+j-Np+M4)=tmp2(j,i)
		enddo
	enddo

	return
end subroutine expinter


!---------------------------------------------------
!       This subrotine copy matrix to matrix
!---------------------------------------------------
subroutine copying(n,ini,fin)
	implicit none
	integer n,i
	_REAL_ ini(n),fin(n)

	do i=1,n
		fin(i)=ini(i)
	enddo

	return

!---------------------------------------------------
!       This subrotine substract two matrices
!---------------------------------------------------
entry substract(n,ini,fin)

	do i=1,n
		fin(i)=fin(i)-ini(i)
	enddo

	return

!---------------------------------------------------
!      sum two matrices (should be deprecated in F90)
!---------------------------------------------------
entry summing(n,ini,fin)

	do i=1,n
		fin(i)=fin(i)+ini(i)
	enddo

	return

end subroutine copying

!---------------------------------------------------
!    copy triangular matrix to square matrix
!---------------------------------------------------
subroutine unpacking(n,ini,fin,su)
	implicit none
      	integer n,i,j
      	character*1 su
      	_REAL_ ini(n*(n+1)/2),fin(n,n)

      	if (su.eq.'s') then
       		do i=1,n
      			do j=1,i
   				fin(j,i)=ini(i*(i-1)/2+j)
     				 fin(i,j)=fin(j,i)
       			enddo
       		enddo
      	elseif (su.eq.'u') then
       		do i=1,n
     			do j=1,i
    				fin(j,i)=ini(i*(i-1)/2+j)
   				fin(i,j)=-fin(j,i)
    			enddo
      			fin(i,i)=0.0
     		enddo
     	else
       		write(6,*)  'Unrecognized flag to unpacking'
     	endif
     	return

!---------------------------------------------------
!       This subrotine copy square matrix to triangular matrix
!---------------------------------------------------
entry packing(n,fin,ini,su)

      if (su.eq.'s') then
       do i=1,n
         do j=1,i
            ini(i*(i-1)/2+j)=0.5*(fin(j,i)+fin(i,j))
        enddo
       enddo
      elseif (su.eq.'u') then
       do i=1,n
         do j=1,i
            ini(i*(i-1)/2+j)=0.5*(fin(j,i)-fin(i,j))
         enddo
            ini(i*(i-1)/2+i)=0.0
       enddo
      else
       write(6,*)  'Unrecognized flag to packing'
      endif

      return
end subroutine unpacking

!---------------------------------------------------
!       This subrotine write-read array to disk
!---------------------------------------------------
subroutine wrb_arr (arr,idim,fname,rw,sd)

      implicit none

      character fname*(*),rw*1,sd*1
      integer idim,i
      real*8 arr(idim)

      open (10,file=fname,form='unformatted')

      if (rw.eq.'w') then
      write (10) (arr(i),i=1,idim)
      elseif (rw.eq.'r') then
      read (10) (arr(i),i=1,idim)
      else
      write(6,*)  'Cannot read or write array'
      endif

      if (sd.eq.'d') then
      close (10,status='delete')
      else
      close (10)
      endif

    return
end subroutine wrb_arr

