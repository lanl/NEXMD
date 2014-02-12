!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine write-read array to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		
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
      print *, 'Cannot read or write array'
      endif
      
      if (sd.eq.'d') then
      close (10,status='delete')
      else
      close (10)
      endif

!      open (10,file=fname)
!      do i=1,idim
!      write (10,*) arr(i)
!      enddo
!      close (10)

	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine write-read coulomb matrix to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		
      subroutine wrb_arr1(arr,ar1,ar2,ar3,ar4,idim,fname,rw,sd)
      
      implicit none
      
      character fname*(*),rw*1,sd*1
      integer idim,i
      integer ar1(idim),ar2(idim),ar3(idim),ar4(idim)
      real*8 arr(idim)	
      
      open (10,file=fname,form='unformatted')
            
      if (rw.eq.'w') then
      write (10) (arr(i), i=1,idim) 
      write (10) (ar1(i), i=1,idim)
	write (10) (ar2(i), i=1,idim)
	write (10) (ar3(i), i=1,idim)
	write (10) (ar4(i), i=1,idim)
	elseif (rw.eq.'r') then
      read (10) (arr(i), i=1,idim) 
      read (10) (ar1(i), i=1,idim)
	read (10) (ar2(i), i=1,idim)
	read (10) (ar3(i), i=1,idim)
	read (10) (ar4(i), i=1,idim)	 
      else 
      print *, 'Cannot read or write arrays'
      endif
      
      if (sd.eq.'d') then
      close (10,status='delete')
      else
      close (10)
      endif
      
	return
	end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine write-read coulomb matrix to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		
      subroutine wrb_arr2(idim,WJ,WK,GSS,GSP,GPP,GP2,
     + HSP,GSD,GPD,GDD,fname,rw,sd)
      
      implicit none
      
      character fname*(*),rw*1,sd*1
      integer idim,i
      real*8 WJ(idim),WK(idim)                    ! Coulomb matrix elements
      real*8 GSS(107),GSP(107),GPP(107),GP2(107),
     1       HSP(107),GSD(107),GPD(107),GDD(107)
      
      open (11,file=fname,form='unformatted')
            
      if (rw.eq.'w') then
	   write (11) (WJ(i),i=1,idim)
	   write (11) (WK(i),i=1,idim)
	   write (11) (GSS(i),i=1,107)
	   write (11) (GSP(i),i=1,107)
	   write (11) (GPP(i),i=1,107)
	   write (11) (GP2(i),i=1,107)
	   write (11) (HSP(i),i=1,107)
	   write (11) (GSD(i),i=1,107)
	   write (11) (GPD(i),i=1,107)
	   write (11) (GDD(i),i=1,107)
      elseif (rw.eq.'r') then
	   read (11) (WJ(i),i=1,idim)
	   read (11) (WK(i),i=1,idim)
	   read (11) (GSS(i),i=1,107)
	   read (11) (GSP(i),i=1,107)
	   read (11) (GPP(i),i=1,107)
	   read (11) (GP2(i),i=1,107)
	   read (11) (HSP(i),i=1,107)
	   read (11) (GSD(i),i=1,107)
	   read (11) (GPD(i),i=1,107)
	   read (11) (GDD(i),i=1,107)
      else 
      print *, 'Cannot read or write arrays'
      endif
      
      if (sd.eq.'d') then
      close (11,status='delete')
      else
      close (11)
      endif
      
	return
	end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine write-read matrix to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		
      subroutine wrb_matr (matr,idim1,idim2,ldim,fname,rw,sd)
      
      implicit none
      
      character fname*(*),rw*1,sd*1
      integer idim1,idim2,ldim,i,j
      real*8 matr(ldim,idim2)
      
      if (idim1.gt.ldim) then 
      print *, 'Cannot read or write matrix'
      goto 100
      endif	
      
      open (10,file=fname,form='unformatted')
            
      if (rw.eq.'w') then
      do j=1,idim1
      write (10) (matr(j,i),i=1,idim2)
      enddo
      elseif (rw.eq.'r') then
      do j=1,idim1
      read (10) (matr(j,i),i=1,idim2)
      enddo 
      else 
      print *, 'Cannot read or write matrix'
      endif
      
      if (sd.eq.'d') then
      close (10,status='delete')
      else
      close (10)
      endif

100    continue

	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine print square matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      subroutine prmat(n,d)
      
      implicit none
      integer n,i,j
      real*8 d(n,n)
       do i=1,n
         write(6,910) (d(i,j),j=1,n)
       enddo
 910  format(' ', 40g11.3)
! 910  format(' ', 20f11.7)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine print nxm matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      subroutine prmat0(n,m,d)
      
      implicit none
      integer n,m,i,j
      real*8 d(n,m)
       do i=1,m
         write(6,910) (d(i,j),j=1,n)
       enddo
 910  format(' ', 10g11.3)
      return
      end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine print triangular matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      subroutine prmat1(n,d)
      
      implicit none
      integer n,i,j
      real*8 d(*)
      do i=1,n
         write(6,910) (d(i*(i-1)/2+j),j=1,i)
      enddo
! 910  format(' ', 10g10.3)
 910  format(' ', 20f11.7)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine print array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      subroutine prarr(n,d)
      
      implicit none
      integer n,i,j
      real*8 d(n)
         write(6,910) (d(j),j=1,n)
 910  format(' ', 10g11.3)

      return
      end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine copy matrix to matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
      subroutine copying(n,ini,fin)
      implicit none
      integer n,i
      real*8 ini(n),fin(n)
      
      do i=1,n
         fin(i)=ini(i)
      enddo
      
      return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine substract two matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      entry substract(n,ini,fin)

      do i=1,n
         fin(i)=fin(i)-ini(i)
      enddo
      
      return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine sum two matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	
      entry summing(n,ini,fin)

      do i=1,n
         fin(i)=fin(i)+ini(i)
      enddo
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine copy triangular matrix to square matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
      subroutine unpacking(n,ini,fin,su)
      implicit none
      integer n,i,j
      character*1 su
      real*8 ini(n*(n+1)/2),fin(n,n)

      if (su.eq.'s') then
       do i=1,n
         do j=1,i
            fin(j,i)=ini(i*(i-1)/2+j)
            fin(i,j)=fin(j,i)
        enddo
       enddo
      endif 
      if (su.eq.'u') then
       do i=1,n
         do j=1,i
            fin(j,i)=ini(i*(i-1)/2+j)
            fin(i,j)=-fin(j,i)
        enddo
	 fin(i,i)=0.0
       enddo 
      endif      
      goto 10
      
      print *, 'Unrecognized flag to unpacking'

          
10    continue
     
      return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine copy square matrix to triangular matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
      entry packing(n,fin,ini,su)

      if (su.eq.'s') then
       do i=1,n
         do j=1,i
            ini(i*(i-1)/2+j)=0.5*(fin(j,i)+fin(i,j))
        enddo
       enddo
      endif 
      if (su.eq.'u') then
       do i=1,n
         do j=1,i
            ini(i*(i-1)/2+j)=0.5*(fin(j,i)-fin(i,j))
	 enddo
	    ini(i*(i-1)/2+i)=0.0
       enddo 
      endif
      goto 20
      
      print *, 'Unrecognized flag to packing'
           
20    continue
           
      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This subroutine commutes two matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine commut(n,d1,d2,dd)
      implicit none
      integer n,i,j,k
      real*8 d1(n,n),d2(n,n),dd(n,n),tr
      
       do i=1,n
        do j=1,n
          tr=0.0
           do k=1,n
            tr=tr+d1(i,k)*d2(k,j)-d2(i,k)*d1(k,j)
           enddo 
         dd(i,j)=tr
        enddo
       enddo

       return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This subroutine anticommutes two matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      entry anticommut(n,d1,d2,dd)

       do i=1,n
        do j=1,n
          tr=0.0
           do k=1,n
            tr=tr+d1(i,k)*d2(k,j)+d2(i,k)*d1(k,j)
           enddo 
         dd(i,j)=tr
        enddo
       enddo

       return
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This subroutine multiples two matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      entry multiple(n,d1,d2,dd)

       do i=1,n
        do j=1,n
          tr=0.0
           do k=1,n
            tr=tr+d1(i,k)*d2(k,j)
           enddo 
         dd(i,j)=tr
        enddo
       enddo

       return
       end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This subroutine trace matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine trace(n,dd,tr)
      implicit none
      integer n,k
      real*8 dd(n,n),tr

          tr=0.0
        do k=1,n
           tr=tr+2*dd(k,k)
        enddo
	
        return
	end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine transpose matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine transp(n,d1,d2)
      implicit none
      integer n,i,j
      real*8 d1(n,n),d2(n,n)

        do i=1,n
          do j=1,n
            d2(j,i)=d1(i,j)
          enddo
	enddo
	
           return
            end

      subroutine transp1(n,d1)
      implicit none
      integer n,i,j
      real*8 d1(n,n),f

        do i=1,n
          do j=1,i
	      f=d1(i,j)
		d1(i,j)=d1(j,i)
            d1(j,i)=f
          enddo
	enddo
	
           return
            end

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
            end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine zeros array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine clearing(n,d)
      implicit none
      integer n,i
      real*8 d(n)
      
	do i=1,n
	 d(i)=0.0
      enddo
	
      return
       end	

      subroutine clearingz(n,d)
      implicit none
      integer n,i
      complex*16 d(n)
     
	do i=1,n
	 d(i)=(0.0,0.0)
      enddo
	
      return
       end
	 	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subrotine get i-th mode from the mode array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getmode(M2_M,Mx_M,M2,i,v0,tmp)
      implicit none
      integer M2_M,Mx_M,M2,i,j
	real*8 v0(M2_M,Mx_M),tmp(M2)

      do j=1,M2
	  tmp(j)=v0(j,i)
	enddo
	
           return
            end
       
      subroutine getmodez(M2_M,Mx_M,M2,i,v0,tmp)
      implicit none
      integer M2_M,Mx_M,M2,i,j
	real*8 v0(M2_M,Mx_M)
	complex*16 tmp(M2)

      do j=1,M2
	  tmp(j)=dcmplx(v0(j,i),0.0d0) 
	enddo
	
           return
            end

      subroutine getmodef(M2_M,Mx_M,Np,Nh,i,v0,temp)
      implicit none
      integer M2_M,Mx_M,i,j,k,jj,Np,Nh
        real*8 v0(M2_M,Mx_M),temp(Np+Nh,Np+Nh)

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
            end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subroutine expand interband matrix
!       from 2*M4 to full Mb and back
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine expinter(Nb,Np,Nh,M4,tmp1,tmp2)

      implicit none
      integer Nb,Np,Nh,M4,j,i,k
        real*8  tmp1(2*M4),tmp2(Nb,Nb)

! Particle-particle and Hole-hole
        call clearing(Nb*Nb,tmp2)
! Particle-hole
      do i=Np+1,Nb
         do j=1,Np
!         tmp2(j,i)=tmp1(Np*(i-Np-1)+j)
! tmp1 needs to be transposed
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
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subroutine project matrix to interband 
!       space 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      subroutine project(Nb,Np,Nh,tmp)

      implicit none
      integer Nb,Np,Nh,M4,j,i,k
      real*8  tmp(Nb,Nb)
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
        end        	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subroutine symmetrize square matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      subroutine symmetr(Nb,tmp)

      implicit none
      integer Nb,j,i,k
      real*8  tmp(Nb,Nb)

      do i=1,Nb
         do j=1,i
          tmp(i,j)=0.5*(tmp(i,j)+tmp(j,i))
	    tmp(j,i)=tmp(i,j)
         enddo
      enddo

        return
        end        
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This subroutine computes (I-2rho)*xi
!       and xi*(I-2rho) for real and complex matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Iminus2rho(Nb,Np,temp1,temp2) 
          
      implicit none
      integer Nb,Np,j,i
        real*8  temp1(Nb,Nb),temp2(Nb,Nb)
        
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
        end

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
        end
