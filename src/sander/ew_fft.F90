! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"


!-------------------------------------------------------------------
!     --- GET_FFTDIMS ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_fftdims here]
subroutine get_fftdims(nfft1,nfft2,nfft3, &
      nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork, &
      sizfftab,sizffwrk)

   use trace
   implicit none
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
         nfftable,nffwork,sizfftab,sizffwrk
   integer n,nfftmax
#ifdef MPI
   integer numtasks0
#  include "parallel.h"
#  include "ew_parallel.h"
#endif
   call trace_enter( 'get_fftdims' )
   nfftdim1 = nfft1
   n = nfft1/2
   if ( nfft1 == 2*n )then
      nfftdim1 = nfft1/2+2
   else
      call sander_bomb("GET_FFTDIMS (ew_fft.f)", &
            "For RealComplex FFT","nfft1 must be even")
   end if
   nfftdim2 = nfft2
   n = nfft2/2
   if ( nfft2 == 2*n )nfftdim2 = nfft2+1
   nfftdim3 = nfft3
   n = nfft3/2
   if ( nfft3 == 2*n )nfftdim3 = nfft3+1
   nfftmax = max(nfft1,nfft2,nfft3)

#ifdef MPI
   nfftable = 4*nfftmax + 15
   sizfftab = 3*nfftable

   !       space work for the individual 1D ffts
   nffwork = 2*nfftdim2

   
   !       space for the x-z slabs
   !       The x-y slabs are tranposed into the xz slabs in
   !       the work array.
   
   numtasks0 = numtasks
   numtasks = num_recip
   nffwork = nffwork + nfftdim1*nfft3*2*(nfft2/numtasks + 1)
   
   !         set up division of fft grid among processors
   
   call par_fft_setup(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
         indz,ntxyslab,ntxzslab,mxyslabs,mxzslabs, &
         nxyslab, nxzslab, mxystart,mxzstart, num_recip)
   
   !       space for the tmp arrays used in transpose communications
   !          assume that node 0 has the largest nxyslab and nxzslab
   
   if ( numtasks > 1)then
      ind_tr_tmp=nffwork+1
      ind_tr_tmp1=ind_tr_tmp+2*nfftdim1*nxyslab(0)*nxzslab(0)
      nffwork=nffwork+4*nfftdim1*nxyslab(0)*nxzslab(0)
   end if
   
   sizffwrk  = nffwork
   numtasks = numtasks0

#else /* not MPI */
   !-----------------------------------------------------------

   nfftable = 4*nfftmax + 15
   sizfftab = 3*nfftable

   !       space work for the individual 1D ffts

   nffwork = 2*nfftdim2

   nffwork = nffwork + nfftdim1*nfft3*2*(nfft2+ 1)
   
   sizffwrk  = nffwork

#endif /* MPI */

   call trace_exit( 'get_fftdims' )
   return
end subroutine get_fftdims 

!-------------------------------------------------------------------
!     --- FFT_BACK RC ---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft_backrc here]
subroutine fft_backrc(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork, tmpy, alpha,beta)

   use trace
   implicit none

   _REAL_  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork
   _REAL_  tmpy(*), &
         alpha(*),beta(*)

   integer isign
   _REAL_  scale
   
   call trace_enter( 'fft_backrc' )
   isign = -1
   scale = 1.d0
   call fft3d0rc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta )

   call trace_exit( 'fft_backrc' )
   return
end subroutine fft_backrc 

!-------------------------------------------------------------------

!     --- FFT_FORWARD RC---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft_forwardrc here]
subroutine fft_forwardrc(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork,tmpy,alpha,beta)

   use trace
   implicit none
   
   _REAL_  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork
   _REAL_  tmpy(*), &
         alpha(*),beta(*)

   integer isign
   _REAL_  scale
   
   call trace_enter( 'fft_forwardrc' )
   isign = 1
   scale = 1.d0
   call fft3d_zxyrc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,array,nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta )
   
   call trace_exit( 'fft_forwardrc' )
   return
end subroutine fft_forwardrc 

!-------------------------------------------------------------------

!     --- FFT_SETUP ---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft_setup here]
subroutine fft_setup(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork)
   
   use trace
   implicit none
   
   _REAL_  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork

   _REAL_  alpha,beta,tmpy
   integer isign
   _REAL_  scale
   
   call trace_enter( 'fft_setup' )
   isign = 0
   scale = 1.d0
   call fft3d0rc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta )

   call trace_exit( 'fft_setup' )
   return
end subroutine fft_setup 

#ifdef MPI

!**************************************************************************

!            PAR_FFT_SETUP
!     Setup routine for parallel 3d fft.
!          Author M. Crowley
!     Only run by MASTER PE at present
!         then all is broadcast from startup() in parallel.f
!   Sets the values of the variables in ew_parallel.h
!      These have the parameters of the distributed work in the
!      ffts and the other ewald routines.
!**************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine par_fft_setup here]
subroutine par_fft_setup(n1,n2,n3,ldx,ldx2,ldx3, indz, &
      ntxyslab,ntxzslab,mxyslabs,mxzslabs, &
      nxyslab, nxzslab, mxystart,mxzstart,num_recip)

   use trace
   implicit none
   
   !  #include "ew_parallel.h"
   integer indz,ntxyslab,ntxzslab
   integer mxyslabs,mxzslabs
   integer nxyslab(0:*)
   integer nxzslab(0:*)
   integer mxystart(0:*)
   integer mxzstart(0:*)
   integer num_recip

#  include "parallel.h"

   integer n1, n2, n3, ldx, ldx2, ldx3, numtasks0 ,i,n3less
   integer n3all,n2all,n3left,n2left

   call trace_enter( 'par_fft_setup' )
   !This is good for pub 1d fft and SGI 1d fft:
   indz = 2*n2

   !-----------------------------------------------------------------------
   !      Calculate the work distribution:
   !         each PE should get approx. the same number of slabs in the
   !            x-y direction.
   !         n3  =  total number of slabs
   !       Slabs will be distributed in a non-circular fashion so that
   !               the shmem puts will be of blocks of memory, i.e.,
   !               there will be only one shmem_put per PE
   !         n3all = number of slabs each pe gets
   !         n3left= number of slabs left over
   !                 these are given to the first n3left PEs as
   !                 an increment to their slabs (n3all + 1)
   !         ntxyslab = total number of allocated spots in a slab
   !                    this is the number of complex values that has to
   !                    be passed per slab to each PE for 2D transforming
   !         nxyslab(ipe) = number of slabs ipe will have to do.
   !         mxyslabs = number of slabs this PE will have to do.
   !         mxzstart(ipe) = last slab of the previous pe.
   
   !       Same type of variables for the xz slabs.
   !         n2all, n2left, ntxzslab, mxzslabs, nxzslab(0:255), mxzstart(0:255)
   !----------------------------------------------------------------------
   !      The grid n1xn2xn3 is divided into
   !                   x-y slabs n1xn2xnxyslab(taskid)
   !                   x-z slabs n1xn3xnxzslab(taskid)
   !      each processor will get a slab with the same number of
   !         x-y layers in it unless the layers do not evenly divide
   !         among the processors. In that case, the first processors get
   !         extra each until all layers are gone.
   
   numtasks0 = numtasks
   if(num_recip /= 0)numtasks = num_recip
   ntxyslab = ldx * ldx2 *2
   ntxzslab = ldx * n3 *2


   n3all  = n3/numtasks
   n3left = mod(n3,numtasks)
   n2all  = n2/numtasks
   n2left = mod(n2,numtasks)
   n3less=numtasks-n3left

   mxystart(0) = 0
   if(n3left == 0)then
      nxyslab(0)  = n3all
   else
      nxyslab(0)  = n3all+1
   end if

   do i=1,n3left-1
      nxyslab(i) = n3all+1
   end do
   do i=1,n3left-1
      mxystart(i) = mxystart(i-1) + nxyslab(i-1)
   end do

   do i=max(n3left,1),numtasks-1
      nxyslab(i) = n3all
   end do
   do i=max(n3left,1),numtasks-1
      mxystart(i) = mxystart(i-1) + nxyslab(i-1)
   end do
   mxyslabs = nxyslab(mytaskid)

   mxzstart(0) = 0
   if(n2left == 0)then
      nxzslab(0)  = n2all
   else
      nxzslab(0)  = n2all+1
   end if

   do i=1,n2left-1
      nxzslab(i) = n2all+1
   end do
   do i=1,n2left-1
      mxzstart(i) = mxzstart(i-1) + nxzslab(i-1)
   end do

   do i=max(n2left,1),numtasks-1
      nxzslab(i) = n2all
   end do
   do i=max(n2left,1),numtasks-1
      mxzstart(i) = mxzstart(i-1) + nxzslab(i-1)
   end do
   mxzslabs = nxzslab(mytaskid)

   numtasks = numtasks0
   call trace_exit( 'par_fft_setup' )
   return
end subroutine par_fft_setup 

!****************************************************************
!                   XY_ZX_TRANSPOSE
!****************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine xy_zx_transpose here]
subroutine xy_zx_transpose(targ,src,ldx,n3,tmp,tmp1)

   use trace
   implicit none
   _REAL_  targ(*), src(*), tmp(*),tmp1(*)
   integer n3, ldx, ktask, num_recv, numval, numtrans, numtries
   integer i,j,k,ks,k0,k00,jjtask,jtask,ibuff,ifoo,rtask
#  include "parallel.h"
   include 'mpif.h'
   integer ierr
#  include "ew_parallel.h"
   integer ireq,isnd_stat(mpi_status_size)

   call trace_enter( 'xy_zx_transpose' )

   rtask = mytaskid
   do jjtask = mytaskid+1, mytaskid+numtasks-1
      jtask = mod(jjtask, numtasks)
      numval = 0
      if(jtask /= mytaskid) then
         do i = 0,nxzslab(jtask)-1
            k00 = (mxzstart(jtask) + i )* ldx*2
            do j = 0, ldx-1
               k0 = k00 + j*2
               do ks = 0,nxyslab(mytaskid)-1
                  k = k0 + ks * ntxyslab + 1
                  numval = numval+1
                  tmp(numval) = src(k)
                  numval = numval+1
                  tmp(numval) = src(k+1)
               end do
            end do
         end do

         call trace_mpi('mpi_*send', &
               numval,'MPI_DOUBLE_PRECISION',jtask)
#   ifdef MPI_BUFFER_SIZE
         call mpi_bsend( &
               tmp(1), numval, MPI_DOUBLE_PRECISION, &
               jtask, 11, &
               recip_comm, ierr )
#   else
         call mpi_isend( &
               tmp(1), numval, MPI_DOUBLE_PRECISION, &
               jtask, 11, &
               recip_comm, ireq, ierr  )
#   endif

         rtask = rtask - 1
         if ( rtask < 0 ) rtask = rtask + numtasks
         call xy_zx_trans_recv(targ,src,ldx,n3,tmp1,rtask)
      end if
#   ifndef MPI_BUFFER_SIZE
      call mpi_wait(ireq, isnd_stat,ierr)
#   endif
   end do  !  jjtask = mytaskid+1, mytaskid+numtasks-1

   !-----------------------------------------------------------------------

   call trace_exit( 'xy_zx_transpose' )
   return
end subroutine xy_zx_transpose 

!**************************************************************************
!               XY_ZX_TRANS_RECV
!**************************************************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine xy_zx_trans_recv here]
subroutine xy_zx_trans_recv(targ,src,ldx,n3,foo,ktask)

   use trace
   implicit none

   logical flag
#  include "parallel.h"
   include 'mpif.h'
   integer ierr
#  include "ew_parallel.h"

   _REAL_  targ(*), src(*), foo(*)
   integer n3, ldx, ktask, num_recv, numval
   integer i,j,ks,k0,k00,jtask, ibuff,ifoo
   integer istart, iend

   integer status(mpi_status_size)

   call trace_enter( 'xy_zx_trans_recv' )
   numval = 2*ldx*nxzslab(mytaskid)*nxyslab(ktask)

   call trace_mpi('mpi_recv',numval,'MPI_DOUBLE_PRECISION',ktask)
   call mpi_recv(foo(1), numval, MPI_DOUBLE_PRECISION, &
         ktask, 11, &
         recip_comm, status, ierr )

   numval = 0
   do i = 0,mxzslabs-1
      k00 = i*ntxzslab + mxystart(ktask)*2
      do j = 0, ldx-1
         k0 = k00 + j * n3*2
         do ks = 1,nxyslab(ktask)*2,2
            numval = numval+1
            targ(ks + k0) = foo(numval)
            numval = numval+1
            targ(ks + k0+1) = foo(numval)
         end do
      end do
   end do

   call trace_exit( 'xy_zx_trans_recv' )
   return
end subroutine xy_zx_trans_recv 

!****************************************************************
!     ZX_XY_TRANSPOSE
!****************************************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zx_xy_transpose here]
subroutine zx_xy_transpose(targ,src,ldx,n3,tmp,tmp1)

   use trace
   implicit none
   _REAL_  targ(*), src(*), tmp(*), tmp1(*)
   integer n3, ldx, ktask, num_recv, numval, numtrans, numtries
   integer i,j,k,ks,k0,k00,jjtask,jtask, ibuff,ifoo,rtask

#  include "parallel.h"
   include 'mpif.h'
   integer ierr
#  include "ew_parallel.h"
   integer ireq,isnd_stat(mpi_status_size)

   call trace_enter( 'zx_xy_transpose' )
   rtask = mytaskid
   do jjtask = mytaskid+1, mytaskid+numtasks-1
      jtask = mod(jjtask, numtasks)
      numval = 0
      if(jtask /= mytaskid) then
         do i = 0,nxyslab(jtask)-1
            k00 = (mxystart(jtask) + i)*2
            do j = 0, nxzslab(mytaskid)-1
               k0 = k00 + j*ntxzslab
               do ks = 0,ldx-1
                  k = k0 + ks * n3*2 + 1
                  numval = numval+1
                  tmp(numval) = src(k)
                  numval = numval+1
                  tmp(numval) = src(k+1)
               end do
            end do
         end do

         call trace_mpi('mpi_*send', &
               numval,'MPI_DOUBLE_PRECISION',jtask)
#   ifdef MPI_BUFFER_SIZE
         call mpi_bsend( &
               tmp(1), numval, MPI_DOUBLE_PRECISION, &
               jtask, 11, &
               recip_comm, ierr )
#   else
         call mpi_isend( &
               tmp(1), numval, MPI_DOUBLE_PRECISION, &
               jtask, 11, &
               recip_comm, ireq, ierr )
#   endif

         rtask = rtask - 1
         if ( rtask < 0 ) rtask = rtask + numtasks
         call zx_trans_recv(targ,src,ldx,n3,tmp1,rtask)
      end if
#   ifndef MPI_BUFFER_SIZE
      call mpi_wait(ireq, isnd_stat,ierr)
#   endif
   end do  !  jjtask = mytaskid+1, mytaskid+numtasks-1

   call trace_exit( 'zx_xy_transpose' )
   return
end subroutine zx_xy_transpose 

!****************************************************************
!     ZX_TRANS_RECV
!****************************************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zx_trans_recv here]
subroutine zx_trans_recv(targ,src,ldx,n3,foo,ktask)

   use trace
   implicit none

   logical flag

   _REAL_  targ(*), src(*), foo(*)
   integer n3, ldx, m1, m2, m3, ktask, num_recv, numval
   integer i,j,k,ks,k0,k00,jtask, ibuff,ifoo
   integer istart, iend

#  include "parallel.h"
   include 'mpif.h'
   integer ierr
#  include "ew_parallel.h"

   integer status(mpi_status_size)

   call trace_enter( 'zx_trans_recv' )
   numval = 2*ldx*nxyslab(mytaskid)*nxzslab(ktask)

   call trace_mpi('mpi_recv',numval,'MPI_DOUBLE_PRECISION',ktask)
   call mpi_recv(foo(1), numval, MPI_DOUBLE_PRECISION, &
         ktask, 11, &
         recip_comm, status, ierr )

   numval = 0
   do i = 0,mxyslabs-1
      k00 = i*ntxyslab + mxzstart(ktask)*ldx*2
      do j = 0,nxzslab(ktask)-1
         k0 = k00 + j*ldx*2
         do ks = 1,ldx*2,2
            numval = numval+1
            targ(k0+ks) = foo(numval)
            numval = numval+1
            targ(k0+ks+1) = foo(numval)
         end do
      end do
   end do

   call trace_exit( 'zx_trans_recv' )
   return
end subroutine zx_trans_recv 

!**************************************************************
!      subroutine fft2d
!   Complex version May 22 1995
!**************************************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft2d here]
subroutine fft2d(isign,n1, n2, scale, x, ldx, y, ldy, &
      table, work, isys)
   !**************************************************************
   !     1 PE non-distributed 2 dim. FFT
   !         calls a 1 dim FFT (CCFFT for the 1)
   
   !         Author: Michael F. Crowley
   !                 Pittsburgh Supercomputing Center
   !                 Oct 20, 1994
   !**************************************************************

   use trace
   use constants, only : one
   implicit none
   _REAL_  x(*), y(*)
   _REAL_  work(*), table(*)
   _REAL_  scale
   integer isys,isign,n1,n2,n3,ldx,ldy
   integer i, ibegw, idx, idy, j, jy

#  include "parallel.h"
#  include "ew_parallel.h"

   !---------------------------------------------------------------

   call trace_enter( 'fft2d' )
   scale = one
   if(isign == 0)then

      call cffti(n1,table)
      call cffti(n2,table(4*n1+15 + 1))

      call trace_exit( 'fft2d' )
      return
   end if

   if(isign == -1)then
      !---------------------------------------------------------------
      !     First the x direction, the data is already contiguous
      
      !-----------------------------------------------------------------------
      do i = 0,n2-1
         idx = i*ldx*2+1
         idy = i*ldy*2+1
         do j = 0, 2*n1 -2
            y(idy+j)=x(idx+j)
         end do
         if(isign == 1)then
            call cfftb(n1, y(idy), table)
         else
            call cfftf(n1, y(idy), table)
         end if
      end do
      !-----------------------------------------------------------------------
      !     Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !     then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !     area necessary for the machine specific fft.
      !     for pubfft, there is no work area used so ibegw=0
      
      ibegw=0
      do i = 1, 2*n1-1, 2
         do j = 1, 2*n2, 2
            jy = i  + (j-1)*ldy
            work(ibegw+j)   = y(jy)
            work(ibegw+j+1) = y(jy+1)
         end do

         if(isign == 1)then
            call cfftb(n2, work(ibegw+1), table(4*n1+16))
         else
            call cfftf(n2, work(ibegw+1), table(4*n1+16))
         end if

         do j = 1, n2*2, 2
            jy = i + (j-1)*ldy
            y(jy)   = work(ibegw+j)
            y(jy+1) = work(ibegw+j+1)
         end do
      end do
      
      !-----------------------------------------------------------------------
      !                  FORWARD 2D FFT
      !-----------------------------------------------------------------------
   else if( isign == 1) then
      !     Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !     then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !     area necessary for the machine specific fft.
      !     for pubfft, there is no work area used so ibegw=0
      
      ibegw=0
      do i = 1, 2*n1-1, 2
         do j = 1, 2*n2, 2
            jy = i  + (j-1)*ldy
            work(ibegw+j)   = y(jy)
            work(ibegw+j+1) = y(jy+1)
         end do

         if(isign == 1)then
            call cfftb(n2, work(ibegw+1), table(4*n1+16))
         else
            call cfftf(n2, work(ibegw+1), table(4*n1+16))
         end if

         do j = 1, n2*2, 2
            jy = i + (j-1)*ldy
            y(jy)   = work(ibegw+j)
            y(jy+1) = work(ibegw+j+1)
         end do
      end do
      !---------------------------------------------------------------
      !     Now the x direction, the data is already contiguous
      
      !-----------------------------------------------------------------------
      do i = 0,n2-1
         idx = i*ldx*2+1
         idy = i*ldy*2+1
         do j = 0, 2*n1 -2
            y(idy+j)=x(idx+j)
         end do
         if(isign == 1)then
            call cfftb(n1, y(idy), table)
         else
            call cfftf(n1, y(idy), table)
         end if

      end do
      
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
   else
      call sander_bomb("EWALD_FFT_2D", &
            "Unknown FFT flag","Bombs away")
   end if  ! (isign == -1)

   call trace_exit( 'fft2d' )
   return
end subroutine fft2d 
#endif /* MPI */

!   +----------------------------------------------------------------+
!   |**************************************************************  |
!   |   *****************************************************        |
!   |       ********************************************             |
!   |              REAL - COMPLEX FFT STUFF HERE                     |
!   |       ********************************************             |
!   |    ****************************************************        |
!   |**************************************************************  |
!   +----------------------------------------------------------------+


!**************************************************************
!         Author: Michael F. Crowley
!                 TSRI
!                 July 2000
!**************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft3d0rc here]
subroutine fft3d0rc(isign,n1,n2,n3, scale, &
      x, ldx,ldx2, table, work, tmpy, alpha, beta )

   use trace
   use constants, only : one
   implicit none

#  include "def_time.h"

   integer isign
   integer ldx,ldx2, n1, n2, n3
   _REAL_  x(*)
   _REAL_  table(*), work(*), scale
   _REAL_  tmpy(*),alpha(*),beta(*)

   integer i, k, k0, ks
   integer j, ja, ja00, jz, jz0, jz00, jj, jidx, jtask

#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#  include "ew_parallel.h"
#else
   integer mxyslabs,mxzslabs,ntxyslab,ntxzslab
   integer indz
#endif

   integer my_first_slab_xy
   integer my_first_slab_xz
   integer n1x

   
   !--------------------------------------------------------------
   !         startup: initialize tables for all three dimensions
   

   call trace_enter( 'fft3d0rc' )
#ifndef MPI
   mxyslabs=n3
   mxzslabs=n2
   ntxyslab = (ldx) * ldx2 *2
   ntxzslab = (ldx) * n3 *2
   my_first_slab_xy = 0
   my_first_slab_xz = 0
   indz = 2*n2
# endif
   n1x = n1/2
   scale = one

   if(isign == 0)then

      call fft2drc(0,n1,n2,one, x,ldx, table,work, &
            tmpy,alpha,beta)
      call cffti(n3,table(4*n1+15 + 4*n2+15 + 1))

      call trace_exit( 'fft3d0rc' )
      return
   end if

   !-----------------------------------------------------------------------
   ! each PE should do their 2D ffts now
   !-----------------------------------------------------------------------
   do j = 1, mxyslabs
      jj = (j-1)*ntxyslab+1
      call fft2drc(isign, n1, n2, one, &
            x(jj),ldx, table, work, tmpy, alpha, beta)
   end do

   !mfc ????  WHAT is this for???
   call zero_array(work(indz+1),n3*ldx*mxzslabs)
   !mfc ????  WHAT was that for???

#ifdef MPI
   my_first_slab_xy = mxystart(mytaskid)
   my_first_slab_xz = mxzstart(mytaskid)

   call timer_barrier( recip_comm )
   call timer_start(TIME_FFTCOMM)
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   if(numtasks > 1) &
         call xy_zx_transpose(work(indz+1),x,ldx,n3, &
         work(ind_tr_tmp),work(ind_tr_tmp1))
   !----------------------------------------------------------------------

   call timer_stop(TIME_FFTCOMM)

#endif /*   MPI   */

   do ks = 0,mxzslabs-1
      ja00 = (my_first_slab_xz+ks)*ldx*2
      jz00 = ks*ntxzslab + 2*my_first_slab_xy
      do j = 0,mxyslabs-1
         jz = jz00 + 2*j +1
         ja = ja00 + j*ntxyslab+1
         do i = 0,n1x
            work(indz+jz+i*n3*2) = x(ja+i*2)
            work(indz+jz+1+i*n3*2) = x(ja+i*2+1)
         end do
      end do
   end do

   !-----------------------------------------------------------------------
   !         END of TRANSPOSE
   !-----------------------------------------------------------------------
   !    Now do Z-FFTs
   !-----------------------------------------------------------------------
   do k = 0,mxzslabs-1
      k0 = k*ntxzslab
      do j = 0, n1x
         jidx=k0 + j*n3*2 +1
         call cfftf(n3,work(indz+jidx),table(4*n1+15 + 4*n2+15 +1))
      end do
   end do
   !*************************************************
   !     Leave it   DISTRIBUTED AS Z-X SLABS
   !*************************************************
   do k = 1,ntxzslab*mxzslabs
      x(k) = work(indz+k)
   end do

   call trace_exit( 'fft3d0rc' )
   return
end subroutine fft3d0rc 



!**************************************************************
!                   FFT3D_ZXYRC Complex-to-Real
!**************************************************************
!         Author: Michael F. Crowley
!                 Pittsburgh Supercomputing Center
!                 Oct 20, 1994
!**************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft3d_zxyrc here]
subroutine fft3d_zxyrc(isign,n1,n2,n3, scale, &
      x, ldx,ldx2, y, ldy, ldy2, &
      table, work,tmpy,alpha,beta )

   use trace
   use constants, only : one
   implicit none
   !****************************************************************
   
   !  isign = 1 for forward fft
   !        =-1 for reverse fft
   !        = 0 initializes table for all three dimensions
   !  n1, n2, n3   dimensions of the fft to be performed
   !  scale  data will be scaled by this number on return
   !         see ccfft3d from Cray for details of how to use this one
   !  x, ldx, ldx2  complex 3-d array
   !                the input array with declared dimensions ldx and ldx2
   !                in the first two dimensions
   !  y, ldy, ldy2  complex 3-d array output 3D array
   !  table  real size 2*(n1+n2+n3) (for Cray CCFFT)
   !  work   workspace size 4*( max(n1,n2,n3) ) (for Cray CCFFT)
   !                     + ldx * ldx2 * 2 * (n3/numtasks + 1)
   !                     + ldx * n3   * 2 * (n2/numtasks + 1)
   !  isys   use 0 (read cray docs on ccfft); currently unused.
   !*************************************************************

#  include "def_time.h"
   integer ldx,ldx2, ldy, ldy2, n1, n2, n3
   _REAL_  x(*), y(*)
   _REAL_  table(*), work(*), scale
   _REAL_  alpha(*),beta(*),tmpy(*)
   integer k, k0, ks
   integer ja, ja00, jz, jz0, jz00, jj, jidx, jtask
   integer  i, j, isign, ndim

#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#  include "ew_parallel.h"
#else
   integer mxyslabs,mxzslabs,ntxyslab,ntxzslab
   integer indz,mytaskid
#endif

   integer my_first_slab_xy
   integer my_first_slab_xz
   integer n1x, loop_limit

   !--------------------------------------------------------------
   !         startup: no initialization possible in this routine
   
   call trace_enter( 'fft3d_zxyrc' )
   if(isign == 0)then
      write(6,*)"fft3d_zxy ERROR, cannot do an isign of 0, I QUIT"
      stop
   end if

   if (isign /= 1 ) then
      write(6,*)'isign for 2nd fft should be 1'
      stop
   end if

#ifndef MPI
   mxyslabs=n3
   mxzslabs=n2
   ntxyslab = ldx * ldx2 *2
   ntxzslab = ldx * n3 *2
   my_first_slab_xy = 0
   my_first_slab_xz = 0
   indz = 2*n2
#endif
   n1x=n1/2

   !***********************************************************************
   !       DISTRIBUTED AS Z-X SLABS, put the data into z area of work
   !***********************************************************************
   loop_limit = ntxzslab*mxzslabs
   do k = 1,loop_limit
      work(indz+k) = x(k)
   end do
   !***********************************************************************
   !**** DO Z FFTs NOW **************************************************
   !***********************************************************************
   do k = 0,mxzslabs-1
      k0 = k*ntxzslab
      !           do j = 0, n1-1
      do j = 0, n1x
         jidx=k0 + j*n3*2 +1

         call cfftb(n3,work(indz+jidx),table(4*n1+15 + 4*n2+15 +1))
      end do
   end do
#ifdef MPI
   !***********************************************************************
   !*****  REDISTRIBUTE INTO XY SLABS *************************************
   !***********************************************************************
   ! Redistribute the data with strided shmem getsc
   
   ! split into two do loops so that the same PEs are not trying
   ! to read each other at the same time.
   
   ! Caution: MFC warns SRB that MPI_BARRIER calls may be
   ! needed before these loops with the shmem gets are entered.
   
   !***********************************************************************
   
   my_first_slab_xy = mxystart(mytaskid)
   my_first_slab_xz = mxzstart(mytaskid)

   call timer_barrier( recip_comm )
   call timer_start(TIME_FFTCOMM)
   
   if(numtasks > 1) &
         call zx_xy_transpose(x,work(indz+1),ldx,n3, &
         work(ind_tr_tmp),work(ind_tr_tmp1))
   
   call timer_stop(TIME_FFTCOMM)

#endif /*   MPI   */

   do ks = 0,mxzslabs-1
      ja00 = (my_first_slab_xz+ks)*ldx*2
      jz00 = ks*ntxzslab + 2*my_first_slab_xy
      do j = 0,mxyslabs-1
         jz = jz00 + 2*j +1
         ja = ja00 + j*ntxyslab+1
         do i = 0,n1x
            x(ja+i*2) = work(indz+jz+i*n3*2)
            x(ja+i*2+1) = work(indz+jz+1+i*n3*2)
         end do
      end do
   end do


   !-----------------------------------------------------------------------
   ! each PE should do their 2D ffts now
   !-----------------------------------------------------------------------
   
   do j = 1, mxyslabs
      jj = (j-1)*ntxyslab+1
      call fft2drc(isign, n1, n2, one, &
            x(jj),ldx, table, work, tmpy, alpha, beta)
   end do

   call trace_exit( 'fft3d_zxyrc' )
   return
end subroutine fft3d_zxyrc 

!**************************************************************
!**************************************************************
!**************************************************************

!************ 2D  FFT real-complex-real ***********************

!**************************************************************
!     pubfft implementation for REAL-to-COMPLEX fft
!**************************************************************
!      subroutine fft2drc
!   Complex version May 22 1995
!**************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft2drc here]
subroutine fft2drc(isign,n1, n2, scale, x, ldx, &
      table, work, tmpy, alpha, beta)
   
   !**************************************************************
   !     1 PE non-distributed 2 dim. FFT
   !         calls a 1 dim FFT (CCFFT for the 1)
   
   !         Author: Michael F. Crowley
   !                 Pittsburgh Supercomputing Center
   !                 Oct 20, 1994
   !**************************************************************
   
   use trace
   use constants, only : TWOPI, one, half, zero
   implicit none
   integer isign,n1,n2,ldx
   _REAL_  scale
   _REAL_  x(2, 0:ldx-1, 0:n2-1)
   _REAL_  work(*), table(*)
   _REAL_  tmpy(2,0:ldx-1)
   _REAL_  alpha(0:n1),beta(0:n1)

   integer i, idx, idy, j, jy
   integer n1rc,kr,kkr,ki,kki,idt,n1x
   integer j1,j2,j3,j4,j1rl,j1im,j2rl,j2im,j3rl,j3im,j4rl,j4im,jdx
   integer istar
   _REAL_  a,b,c,d,pi2n,theta

   !---------------------------------------------------------------

   call trace_enter( 'fft2drc' )
   pi2n=TWOPI/n1
   n1x=n1/2

   !=================================================
   !             initialize fft tables
   
   if(isign == 0)then
      if(mod(n1,2) /= 0)then
         write(6,*)" NEED factor 2 value for nfft1 for RC fft"
         call mexit(6,1)
      end if

      call cffti(n1x,table)
      call cffti(n2,table(4*n1+15 + 1))

      call trace_exit( 'fft2drc' )
      return
   end if

   !-----------------------------------------------------
   do i=0,n1x-1
      theta=pi2n*i
      alpha(i) = cos(theta)
      beta(i)  = sin(theta)
   end do
   !---------------------------------------------------------------
   !       Backward fft real to complex
   !---------------------------------------------------------------
   if(isign == -1)then
      !---------------------------------------------------------------
      !  First the x direction, the data is already contiguous
      
      !-----------------------------------------------------------------------
      do j = 0,n2-1
         do i = 0, n1x-1
            tmpy(1,i)=x(1,i,j)
            tmpy(2,i)=x(2,i,j)
         end do
         call cfftf(n1x, tmpy(1,0), table)
         do i = 1, n1x-1
            a =  half*(tmpy(1,i)+tmpy(1,n1x-i)) ! Real F even
            b =  half*(tmpy(2,i)-tmpy(2,n1x-i)) ! Imag F even
            c =  half*(tmpy(2,i)+tmpy(2,n1x-i)) ! Real F odd
            d = -half*(tmpy(1,i)-tmpy(1,n1x-i)) ! Imag F odd
            x(1,i,j) = a + alpha(i)*c + beta(i)*d
            x(2,i,j) = b + alpha(i)*d - beta(i)*c
         end do
         !--------------------------------------------
         !     DC and nyquist
         x(1,0,j)  =tmpy(1,0)+tmpy(2,0)
         x(2,0,j)  =zero
         x(1,n1x,j)=tmpy(1,0)-tmpy(2,0)
         x(2,n1x,j)=zero
      end do
      
      !-----------------------------------------------------------------------
      !     Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !     then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !     area necessary for the machine specific fft.
      !     for pubfft, there is no work area used so ibegw=0
      
      do j1 = 0, n1x
         i=2*j1+1
         istar=2*(n1-j1)
         do j2 = 0, n2-1
            j=2*j2+1
            work(j)   = x(1,j1,j2)
            work(j+1) = x(2,j1,j2)
         end do
         call cfftf(n2, work(1), table(4*n1+16))
         do j2 = 0, n2-1
            j=2*j2+1
            x(1,j1,j2) = work(j)
            x(2,j1,j2) = work(j+1)
         end do
      end do
      !---------------------------------------------------------------
   else
      !---------------------------------------------------------------
      !     Forward fft complex to real
      !---------------------------------------------------------------
      !-----------------------------------------------------------------------
      !  Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !           then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !           area necessary for the machine specific fft.
      !           for pubfft, there is no work area used so ibegw=0
      
      
      
      do i = 0, n1x
         do j = 0, n2-1

            work(2*j+1) = x(1,i,j)
            work(2*j+2) = x(2,i,j)
         end do
         call cfftb(n2, work(1), table(4*n1+16))
         do j = 0, n2-1
            x(1,i,j) = work(2*j+1)
            x(2,i,j) = work(2*j+2)
         end do
      end do

      do j = 0,n2-1
         do i = 1, n1x-1
            a =  (x(1,i,j)+x(1,n1x-i,j)) ! Real F even
            b =  (x(2,i,j)-x(2,n1x-i,j)) ! Imag F even
            c =  (x(2,i,j)+x(2,n1x-i,j)) ! F odd contrib
            d =  (x(1,i,j)-x(1,n1x-i,j)) ! F odd contrib
            tmpy(1,i) = a - alpha(i)*c - beta(i)*d
            tmpy(2,i) = b + alpha(i)*d - beta(i)*c
         end do
         tmpy(1,0) = (x(1,0,j)+x(1,n1x,j))
         tmpy(2,0) = (x(1,0,j)-x(1,n1x,j))
         call cfftb(n1x, tmpy(1,0), table)
         do i = 0, n1x-1
            x(1,i,j)=tmpy(1,i)
            x(2,i,j)=tmpy(2,i)
         end do
      end do
   end if  ! (isign == -1)

   call trace_exit( 'fft2drc' )
   return
end subroutine fft2drc 

#include "pubfft.F90"

