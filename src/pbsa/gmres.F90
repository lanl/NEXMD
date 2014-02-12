#define _REAL_ double precision

      subroutine gmres(xm,ym,zm,nbnd,nz_num,c,c2,index,index2,bv,phi,xs,EPS1)

!     EPS1 <- accept
!     bv <- f
!     phi <- u
!     xs  <- u0

      implicit none
      integer xm,ym,zm,xmymzm,nbnd,nz_num
      _REAL_ :: c(xm,ym,zm,7),c2(nbnd,27)
      _REAL_ :: bv(1:xm*ym*zm)
      _REAL_ :: phi(1:xm*ym*zm),xs(1:xm*ym*zm)
      integer :: index(xm,ym,zm),index2(xm,ym,zm)
      _REAL_ EPS1

      _REAL_,allocatable :: a(:),u(:),f(:)
      integer,allocatable :: ia(:),ja(:),ig(:)
      integer NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,NNA
      integer MATRIX,IFIRST,ISWTCH,IOUT,IPRINT,IERR

      _REAL_ time1,time2,tol_abs,tol_rel
      integer ist,itrmax,mr

      integer n,nelt,nsave,isym,itol,itmax,iunit,iter,lierr,lenw,leniw
      _REAL_ tol,err
      integer,allocatable :: iwork(:)
      _REAL_,allocatable :: rwork(:)

!     print *,'entered amg'
      
      xmymzm = xm*ym*zm

!       read (90,*),NDA
!       allocate (A(NDA))
!       read (90,*),A

!       read (90,*),NDIA
!       allocate (IA(NDIA))
!       read (90,*),IA

!       read (90,*),NDJA
!       allocate (JA(NDJA))
!       read (90,*),JA

!       read (90,*),NDF
!       allocate (F(NDF))
!       read (90,*),F

!       read (90,*),NDU
!       allocate (U(NDU))
!       read (90,*),U

!       read (90,*),NDIG
!       allocate (IG(NDIG))
!       read (90,*),IG

!     print *,nda,ndia,ndja,ndu,ndf,ndig
      allocate ( a(1:nz_num), stat=ist)
      if ( ist /= 0 ) then 
         print *,'allocate a(1:nz_num) failed',nda,ist
         stop
      end if
      allocate ( ja(1:nz_num), stat=ist)
      if ( ist /= 0 ) then 
         print *,'allocate ja(1:nz_num) failed',ndja,ist
         stop
      end if
      allocate ( ia(1:nz_num), stat=ist)
      if ( ist /= 0 ) then 
         print *,'allocate ia(1:nz_num) failed',ndia,ist
         stop
      end if
      allocate ( u(1:xmymzm)  , stat=ist)
      if ( ist /= 0 ) then 
         print *,'allocate u(1:xmymzm)   failed',ndu,ist
         stop
      end if
      allocate ( f(1:xmymzm)  , stat=ist)
      if ( ist /= 0 ) then 
         print *,'allocate f(1:xmymzm)   failed',ndf,ist
         stop
      end if
!     print *,'allocation success!'

!     print *,'entering seta'
      call SETA(xm,ym,zm,nbnd,nz_num, &
           index,index2,c,c2,bv,a, f, ia, ja)

      u(1:xmymzm) = xs(1:xmymzm)

      call wallclock(time1)

      n = xmymzm
      nelt = nz_num
      isym = 0 
      nsave = 10
      itol = 0
      tol = EPS1
      itmax = 1000 ! itmax = nsave*nrmax
      iunit = 0
      !lenw = 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NL+NU
      lenw =  1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+ nelt
      !leniw = NL+NU+4*N+32
      leniw = nz_num+4*n+32
      allocate(rwork(lenw),iwork(leniw))

      call dslugm(n, f, u, nelt, ia, ja, a, isym, nsave, itol, &
           tol, itmax,iter,err,lierr,iunit,rwork,lenw, iwork,leniw)
      call wallclock(time2)

      xs(1:xmymzm) = u(1:xmymzm)
      phi(1:xmymzm) = u(1:xmymzm)

      deallocate(rwork,iwork)
      deallocate(a,ia,ja)
      deallocate(u)
      deallocate(f)

      contains

      subroutine SETA(l,m,n,maxirr,nz_num,index,index2,c,c2, &
                 ff,a,f,ia,ja)

!     dimension uu(0:l+1,0:m+1,0:n+1)	
      implicit none
      integer l,m,n,maxirr,nz_num
      integer :: index(l,m,n), index2(l,m,n)
      _REAL_ :: c(l,m,n,7), c2(maxirr,27)
      _REAL_ :: ff(l,m,n)
      _REAL_ :: a(nz_num),f(l*m*n)
      integer :: ia(nz_num), ja(nz_num)

      integer ncount,i,j,k,ne,ii,ir,nc1,i0,j0,k0,ndis,ne1
!     integer npos
      
!     print *,'entered seta'
!     print *,"l,m,n",l,m,n,nz_num
      ncount = 0
      
      do k=1, n
         do j=1, m
            do i=1, l
!              print *,i,j,k
               if (index(i,j,k).eq.1 .or. index(i,j,k).eq.5) then
!                 print *,'entered if'
                  ne1 = npos(l,m,n,i,j,k)
                  ncount = ncount + 1
                  a(ncount) = c(i,j,k,1)
                  ia(ncount) = ne1
                  ja(ncount) = ne1
                  f(ne1) = ff(i,j,k)
!                 print *,ncount,ia(ncount),ja(ncount)

                  if ( abs(c(i,j,k,2)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i-1, j, k)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,2)
                     ia(ncount) = ne1
                     ja(ncount) = ne
                  end if 

                  if ( abs(c(i,j,k,3)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i+1, j, k)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,3)
                     ia(ncount) = ne1
                     ja(ncount)=ne
                  end if 

                  if ( abs(c(i,j,k,4)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i, j-1, k)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,4)
                     ia(ncount) = ne1
                     ja(ncount)=ne
                  end if 

                  if ( abs(c(i,j,k,5)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i, j+1, k)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,5)
                     ia(ncount) = ne1
                     ja(ncount)=ne
                  end if 

                  if ( abs(c(i,j,k,6)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i, j, k-1)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,6)
                     ia(ncount) = ne1
                     ja(ncount)=ne
                  end if 

                  if ( abs(c(i,j,k,7)) > 1.d-10 ) then 
                     ne = npos(l,m,n,i, j, k+1)
                     ncount = ncount + 1
                     a(ncount) = -c(i,j,k,7)
                     ia(ncount) = ne1
                     ja(ncount)=ne
                  end if

! WJ
!                 ne = npos(l,m,n,i, j, k)
!                 write(190,*) '++++',i,j,k,ne,'+++++++++'
!                 write(190,*) c(i,j,k,1:7)
!                 do ii=ia(ne),ncount
!                    write(190,*) a(ii),ja(ii)
!                 end do
!                 write(190,*) f(ne)
 
                
               else
                   
                  ir = index2(i,j,k)
                  ne1 = npos(l,m,n,i,j,k)
                  ncount = ncount + 1
                  a(ncount) = -c2(ir,14)
                  ia(ncount) = ne1
                  ja(ncount) = ne1
                  f(ne1) = -ff(i,j,k)
                  
                  nc1 = 0
                  do i0 = i-1, i+1
                  do j0 = j-1, j+1
                  do k0 = k-1, k+1
                     ndis=abs(i-i0)+abs(j-j0)+abs(k-k0)
                     if (ndis .eq. 0) then 
                       nc1 = nc1 + 1
                     else
                       nc1 = nc1 + 1
                       if ( abs(c2(ir,nc1)) > 1.d-10 ) then
                       ne = npos(l,m,n,i0,j0,k0)
                       ncount = ncount + 1
                       a(ncount) = -c2(ir,nc1)
                       ia(ncount) = ne1
                       ja(ncount) = ne
                       end if
                     endif
                  end do
                  end do
                  end do
! WJ
!                 ne = npos(l,m,n,i, j, k)
!                 write (190,'(3i5,$)'),i,j,k
!                 write (190,*),"IRREGULAR"
!                 do ii=ia(ne),ncount
!                    write(190,'(f20.10,$)') a(ii)
!                 end do
!                 write(190,'(f20.10)') f(ne)
!
                
               endif
!     write(102,*) ncount
!     flush(102)
            end do
         end do
      end do
      
      return

      end subroutine
        
      integer function npos(l,m,n,i,j,k)
      implicit none
         integer l,m,n,i,j,k
         npos = i + (j-1)*l + (k-1)*l*m 
         return
      end function npos
 
      end subroutine gmres
