#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine decnvh here]
subroutine decnvh(asol,bsol,nphb,radhb,welhb)
   
   
   ! SUBROUTINE DECoNVolute H-bond parameters
   
   ! This subroutine takes the A and C h-bond coefficients for 10-12 h-bond
   ! interactions and extracts the values of epsilon (well depth) and
   ! r* (ideal interaction distance) corresponding to these sets of coefficients.
   
   ! Author: David A. Pearlman
   ! Date: 4/90
   
   ! INPUT:
   !   ASOL(I)  : The "A" coefficient, corresponding to the 12-power term
   !   BSOL(I)  : The "C" coefficient, corresponding to the 10-power term
   !   NPHB     : The total number of h-bond interaction types
   
   ! OUTPUT
   !   RADHB(I) : The optimal interaction distance for h-bond I.
   !   WELHB(I) : The maximum well depth of h-bond interaction I.
   
   implicit none
   integer:: i, nphb
   _REAL_ :: asol, bsol, coef1, coef2, radhb, small, welhb
   dimension asol(1),bsol(1),radhb(1),welhb(1)
   data small/1.0d-7/
   !     SAVE SMALL
   
   coef1 = 12.0d0/10.0d0
   coef2 = (5.0d0**5)/(6.0d0**6)
   do i = 1,nphb
      if (abs(asol(i)) >= small .and. abs(bsol(i)) > small) then
         radhb(i) = sqrt(coef1*(asol(i)/bsol(i)))
         welhb(i) = coef2*bsol(i)*(bsol(i)/asol(i))**5
      else
         radhb(i) = 0.0d0
         welhb(i) = 0.0d0
      end if
   end do
   
   return
end subroutine decnvh 
