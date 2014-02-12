#include "copyright.h"
#define _REAL_ double precision

module iim_use
  
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine matvec here]
subroutine matvec(m,n,a,x,y)
   implicit none
   integer m,n,i,j
   _REAL_  a(m,n), x(n),y(m)

   do i=1,m
      y(i) = 0.0
      do j=1,n
         y(i) = y(i) + a(i,j)*x(j)
      end do
   end do

   return
end subroutine matvec 


!-------matmat finds the product a*b of two matrices a and b
end module iim_use
