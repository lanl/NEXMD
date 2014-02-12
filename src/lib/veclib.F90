! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

! Several vendors have vectorized math libraries.
! To address portability and to not degrade readability of the code
! too much we use the MKL name of the vectorized routine in our code.
! This module interfaces to non-MKL libraries and defines the
! vectorized routines for those platforms without vector libraries.

#ifndef MKL
!  Implement vectorized routines for platforms without MKL.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized square-root
subroutine vdsqrt( n, x, y )
   
   implicit none
   integer  n
   _REAL_   x(n), y(n)

#ifdef MASSLIB   
   call vsqrt ( y, x, n )
#else
   y(1:n) = sqrt(x(1:n))
#endif
   
   return
end subroutine vdsqrt 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized exponential
subroutine vdexp( n, x, y )
   
   implicit none
   integer  n
   _REAL_   x(n), y(n)
   
#  ifdef MASSLIB
      call vexp ( y, x, n )
#  else
      y(1:n) = exp(x(1:n))
#  endif
   
   return
end subroutine vdexp 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized logarithm
subroutine vdln( n, x, y )
   
   implicit none
   integer  n
   _REAL_   x(n), y(n)
   
#  ifdef MASSLIB
      call vlog ( y, x, n )
#  else
      y(1:n) = log(x(1:n))
#  endif
   
   return
end subroutine vdln 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse square root
subroutine vdinvsqrt( n, x, y )
   implicit none
   integer  n
   _REAL_   x(n), y(n)
  
#  ifdef MASSLIB
      call vrsqrt ( y, x, n )
#  else
      y(1:n) = 1.d0/sqrt( x(1:n) )
#  endif

   return
end subroutine vdinvsqrt 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse
subroutine vdinv( n, x, y )
   implicit none
   integer  n
   _REAL_   x(n), y(n)
   
#  ifdef MASSLIB
      call vrec ( y, x, n )
#  else
      y(1:n) = 1.d0/x(1:n)
#  endif
   
   return
end subroutine vdinv 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Vectorized cosine
subroutine vdcos( n, x, y )
  implicit none
  integer n
  _REAL_  x(n), y(n)
# ifdef MASSLIB
    call vcos ( y, x, n )
# else
    y(1:n) = cos(x(1:n))
# endif

  return
end subroutine vdcos
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Vectorized tanh
subroutine vdtanh( n, x, y )
  implicit none
  integer n
  _REAL_  x(n), y(n)
# ifdef MASSLIB
    call vtanh ( y, x, n )
# else
    y(1:n) = tanh(x(1:n))
# endif

  return
end subroutine vdtanh
!--------------------------------------------------------------

#else

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dummy subroutine to avoid compiling an empty file
subroutine vd_dummy_to_avoid_empty_file( )
   
   implicit none
   return
end subroutine vd_dummy_to_avoid_empty_file 
!--------------------------------------------------------------

#endif
