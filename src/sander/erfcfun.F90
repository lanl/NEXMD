! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the complementary error function.
!-----------------------------------------------------------------------
!     --- ERFCFUN ---
!---------------------------------------------------------------------
!  The algorithm is ??

subroutine erfcfun(x,erfc)

   implicit none
   _REAL_, intent(in)  :: x
   _REAL_, intent(out) :: erfc
   _REAL_ :: absx, c, p, q, nonexperfc, erf

   absx=abs(x)
   if (x > 26.d0) then
      erfc = 0.d0

   else if (x < -5.5d0) then
      erfc = 2.0d0

   else if (absx <= 0.5d0) then
      c = x * x
      p=((-0.356098437018154d-1*c+0.699638348861914d1)*c+   &
            0.219792616182942d2)*c+0.242667955230532d3
      q=((c+0.150827976304078d2)*c+0.911649054045149d2)*c+  &
            0.215058875869861d3
      erf = x*p/q
      erfc = 1.d0-erf

   else if (absx < 4.d0) then
      c=absx
      p=((((((-0.136864857382717d-6*c+0.564195517478974d0)*c+    &
         0.721175825088309d1)*c+0.431622272220567d2)*c+    &
         0.152989285046940d3)*c+0.339320816734344d3)*c+    &
         0.451918953711873d3)*c+0.300459261020162d3
      q=((((((c+0.127827273196294d2)*c+0.770001529352295d2)*c+    &
         0.277585444743988d3)*c+0.638980264465631d3)*c+    &
         0.931354094850610d3)*c+0.790950925327898d3)*c+    &
         0.300459260956983d3
      if ( x > 0.d0 ) then
         nonexperfc = p/q
      else
         nonexperfc = 2.d0*exp(x*x) - p/q
      end if
      erfc = exp(-absx*absx)*nonexperfc

   else
      c=1.d0/(x*x)
      p=(((0.223192459734185d-1*c+0.278661308609648d0)*c+      &
         0.226956593539687d0)*c+0.494730910623251d-1)*c+      &
         0.299610707703542d-2
      q=(((c+0.198733201817135d1)*c+0.105167510706793d1)*c+      &
         0.191308926107830d0)*c+0.106209230528468d-1
      c=(-c*p/q + 0.564189583547756d0)/absx
      if( x > 0.d0 ) then
         nonexperfc = c
      else
         nonexperfc = 2.d0*exp(x*x) - c
      end if
      erfc = exp(-absx*absx)*nonexperfc
   end if
   return
end subroutine erfcfun 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ a function wrapper for the erfcfun subroutine
function erfc(x)
  implicit none
  _REAL_ :: erfc,x
  call erfcfun(x,erfc)
end function erfc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_ee_func here]
subroutine get_ee_func(x,switch,d_switch_dx,ee_type)
   use constants, only : zero, one, two, INVSQRTPI
   implicit none
   _REAL_ x,switch,d_switch_dx
   integer ee_type
   
   !     ---get switch function multiplying the Coulomb interaction 1/r
   !     r has been converted to x by x = dxdr*r for convenience
   
   if ( ee_type == 1 )then
      
      !          ---erfc function for ewald
      
      call erfcfun(x,switch)
      d_switch_dx = -two*exp(-x*x)*INVSQRTPI
      
   else if ( ee_type == 2 )then
      
      !          ---force shift cutoff
      
      if ( x < one )then
         switch = (one - x)**2
         d_switch_dx = -two*(one - x)
      else
         switch = zero
         d_switch_dx = zero
      end if
   end if
   return
end subroutine get_ee_func 


