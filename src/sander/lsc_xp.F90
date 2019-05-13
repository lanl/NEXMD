!#include "../include/assert.fh"
#include "dprec.fh"
!
!------------------------------------------------------------------------------
subroutine lsc_xp(x,p)
!
  use lscivr_vars
  use constants, only : pi, hbar, kB
  implicit none
# include "extra.h"
#ifdef MPI
# include "parallel.h"
#endif
!
  _REAL_, intent(inout):: x(ndof_lsc), p(ndof_lsc)
  integer:: i, ierr, j, iout
  _REAL_ :: pxcom, pycom, pzcom, psqcom, psqnew, rat_psq, psum    
  _REAL_ :: ubeta
! quantum correction factor 
  _REAL_, allocatable :: qcor_lsc(:)
! an auxiliariary matrix for the Kubo-transformed corr. func.
  _REAL_, allocatable :: p_kubo(:)
  integer:: info, nrand
  _REAL_:: xrand
!
   allocate(p_kubo(ndof_lsc), stat = ierr)
   allocate(qcor_lsc(ndof_lsc), stat = ierr)
!

! determine the randome number, this is necessary because each run has to be independent

    read(file_rhoa_lsc,*)nrand
    if(nrand == 0)then
       call gauss(0.0d0, 1.0d0, xrand)  
       rewind(file_rhoa_lsc)
    elseif(nrand >= 10001)then
       call gauss(0.0d0, 1.0d0, xrand)
       rewind(file_rhoa_lsc)
       nrand = 0
    else      
     call gauss(0.0d0, 1.0d0, xrand)  
     do i = 1, nrand
          do j = 1, ndof_lsc
             call gauss(0.0d0, 1.0d0, xrand)
          end do
       enddo  
       rewind(file_rhoa_lsc)
    endif

! 2nd derivative of the mass-weighted potenital
   do ilsc = 1, ndof_lsc
      v2_lsc(ilsc,1:ndof_lsc) = v2_lsc(ilsc,1:ndof_lsc)/sqrt(mass_lsc(ilsc))
   end do
   do ilsc = 1, ndof_lsc
      v2_lsc(1:ndof_lsc,ilsc) = v2_lsc(1:ndof_lsc,ilsc)/sqrt(mass_lsc(ilsc))
   end do   



!
! diagonalize the mass-weighted Hesian matrix
   call dsyev( 'V', 'L', ndof_lsc, v2_lsc, ndof_lsc, v2_freq, & 
               work_lsc, 10*ndof_lsc, info )




! Now the array 'v2_ivr' stores the eigenvectors
   p_lsc = 0.0d0
   sigma_lsc = 0.0d0
   do ilsc = 1, ndof_lsc
   if ( v2_freq(ilsc) >= 1.0d-48 ) then
      omega_lsc = sqrt( v2_freq(ilsc) )
      alpha_lsc(ilsc) = omega_lsc/hbar*cosh(0.50d0*beta_lsc*hbar*omega_lsc) &
                    /sinh(0.50d0*beta_lsc*hbar*omega_lsc)
      ubeta = 0.50d0 * beta_lsc * hbar * omega_lsc  
      qcor_lsc(ilsc) = ubeta * cosh(ubeta) / sinh(ubeta)
   elseif( abs(v2_freq(ilsc)) .lt. 1.0d-48 )then 
      omega_lsc = sqrt( abs(v2_freq(ilsc)) )
      alpha_lsc(ilsc) = cosh(0.50d0*beta_lsc*hbar*omega_lsc)/(0.5d0*beta_lsc*hbar**2)     
      ubeta = 0.50d0 * beta_lsc * hbar * omega_lsc
      qcor_lsc(ilsc) = cosh(ubeta) / (1.0d0+ubeta**2/6.0d0) 
   elseif (v2_freq(ilsc) <= -1.0d-48)then
      omega_lsc = sqrt( - v2_freq(ilsc) )
      alpha_lsc(ilsc) = 2.0d0/(hbar**2*beta_lsc)/(0.50d0*beta_lsc*hbar*omega_lsc)    &
                        *sinh(0.50d0*beta_lsc*hbar*omega_lsc)/cosh(0.50d0*beta_lsc*hbar*omega_lsc)
      ubeta = 0.50d0 * beta_lsc * hbar * omega_lsc
      qcor_lsc(ilsc) = sinh(ubeta) / ( ubeta * cosh(ubeta) ) 
!
   endif
   end do


   do ilsc = 1, ndof_lsc
      if(alpha_lsc(ilsc) .gt. 0.0d0)then
         sigma_lsc = sqrt( 0.5d0 * hbar**2 * alpha_lsc(ilsc) )
         call gauss( 0.0d0, sigma_lsc, p_lsc(ilsc) )
      else
         p_lsc(ilsc) = 0.0d0
      end if
   end do
   nrand = nrand + 1
  
   psum = 0.0d0
   do ilsc = 1, natom_lsc, 3
         psum = psum + p_lsc(3*ilsc-2)**2/(hbar**2 * alpha_lsc(3*ilsc-2) )   &
                     + p_lsc(3*ilsc-1)**2/(hbar**2 * alpha_lsc(3*ilsc-1) )   &
                     + p_lsc(3*ilsc)**2/(hbar**2 * alpha_lsc(3*ilsc) )
   end do
   psum = psum /dble(ndof_lsc/3)
             

   do i = 1, ndof_lsc
      if(i==1) then
         p(i) = 0.0
         do j=1, ndof_lsc
            p(i) = p(i) + v2_lsc(1,j)*p_lsc(j)
         end do
      end if
      p(i) = sum( v2_lsc(i,:) * p_lsc(:) )* sqrt(mass_lsc(i))
   enddo

!  At this point, transfer the momentum p(:) 
!   to the velocity v(:)
   do i = 1, ndof_lsc 
       p(i) = p(i) / mass_lsc(i)
   end do 
!
#ifdef MPI
   if( master )then
   write(file_rhoa_lsc,*)nrand
   write(file_rhoa_lsc,'(10f8.3)') x
   write(file_rhoa_lsc,'(10f8.3)') p
   write(file_rhoa_lsc,*)mass_lsc
   do ilsc = 1, ndof_lsc
      if(v2_freq(ilsc) .le. 0.0d0)then
         write(file_rhoa_lsc,*) -beta_lsc*hbar*sqrt(abs(v2_freq(ilsc)))
      else
         write(file_rhoa_lsc,*) beta_lsc*hbar*sqrt(abs(v2_freq(ilsc)))
      endif
   enddo
   endif
#else
   write(file_rhoa_lsc,*)nrand
   write(file_rhoa_lsc,'(10f8.3)') x
   write(file_rhoa_lsc,'(10f8.3)') p
   write(file_rhoa_lsc,*)mass_lsc
   do ilsc = 1, ndof_lsc
      if(v2_freq(ilsc) .le. 0.0d0)then
         write(file_rhoa_lsc,*) -beta_lsc*hbar*sqrt(abs(v2_freq(ilsc)))
      else
         write(file_rhoa_lsc,*) beta_lsc*hbar*sqrt(abs(v2_freq(ilsc)))
      endif
   enddo
!

#endif

!  stop

! the sign to determine which type of the correlation function
! icorf_lsc = 0,       scalar function
!             1,       position <x(0)x(t)>
!             2,       velocity <v(0)v(t)>
!             3,       nonlinear function <f(0)f(t)>
!             4,       Kubo-transformed <v(0)v(t)>
! Also see lsc_init.f lscivr_vars.f
   if(icorf_lsc == 1)then 
       rho_a = x
   elseif(icorf_lsc == 2)then
       rho_a = p
   elseif(icorf_lsc == 3)then
! Write the quantum correction factor qcor_lsc
!   and the eigenvector matrix v2_lsc
       write(file_Tmat_lsc, *)  qcor_lsc
       write(file_Tmat_lsc, *)  v2_lsc
   elseif(icorf_lsc == 4)then
      p_kubo = 0.0d0
      do i = 1, ndof_lsc
         p_kubo(i) = p_lsc(i)/alpha_lsc(i)
      enddo
      p_kubo = matmul(v2_lsc, p_kubo)
      do i = 1, ndof_lsc
       p_kubo(i) = sqrt(mass_lsc(i))*p_kubo(i)
      end do
      p_kubo = p_kubo * 2.0d0/(hbar**2 * beta_lsc) /(2.0d0*kb)
!
       rho_a = 0.0d0
       rho_a = p_kubo
#ifdef MPI
       if(master)then
       do i = 1, ndof_lsc
          write(file_rhoa_lsc,'(f10.3, I6)') p_kubo(i), i
          write(file_pos_lsc,'(f10.3, I6)') p_kubo(i), i
          write(file_cor_lsc,'(f10.3, I6)') p_kubo(i), i
       end do
       endif
#else
       do i = 1, ndof_lsc
          write(file_rhoa_lsc,'(f10.3, I6)') p_kubo(i), i
          write(file_pos_lsc,'(f10.3, I6)') p_kubo(i), i
          write(file_cor_lsc,'(f10.3, I6)') p_kubo(i), i
       end do
#endif

   endif
!

   if(icorf_lsc == 0)then
     write(file_pos_lsc,* )    &
          rho_a
   elseif(icorf_lsc == 4)then

      
#ifdef MPI
   if(master)then
      write(file_rhoa_lsc,*)rho_a
      write(file_pos_lsc,*)rho_a
   endif
#else
      write(file_rhoa_lsc,*)rho_a
      write(file_pos_lsc,*)rho_a 
#endif

#ifdef MPI
   if(master)then
     write(file_cor_lsc,*)psum
     write(file_rhoa_lsc,*)'psum=',psum 
  endif
#else
     write(file_cor_lsc,*)psum
     write(file_rhoa_lsc,*)'psum=',psum
#endif

   endif

   deallocate(p_kubo, stat = ierr)
   deallocate(qcor_lsc, stat = ierr)
!   stop

   return

  end subroutine lsc_xp
