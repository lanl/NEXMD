#include "copyright.h"
#include "dprec.fh"

module opnq_switching

  implicit none

contains

  subroutine switchoff(x,xlo,xmax,f,deriv)

    implicit none
!
! 1 if x <= xlo
! 0 if x >= xhi
!

    _REAL_,intent(in) :: x,xmax,xlo
    _REAL_,intent(out) :: f
    _REAL_,intent(out),optional :: deriv

    _REAL_ :: s,sf,ds,dsf

    if ( present( deriv ) ) then
       if ( x >= xmax ) then
          f = 0.0d0
          deriv = 0.0d0
       else if ( x <= xlo ) then
          f = 1.0d0
          deriv = 0.0d0
       else
          s  = (xmax-x)/(xmax-xlo)
          ds = -1.0d0 / (xmax-xlo)

          sf =    10.0d0 * s**3 &
               & -15.0d0 * s**4 &
               & + 6.0d0 * s**5

          dsf =    30.0d0 * s**2 &
               & - 60.0d0 * s**3 &
               & + 30.0d0 * s**4

          f = sf
          deriv = dsf * ds

       end if

    else
       if ( x >= xmax ) then
          f = 0.0d0
       else if ( x <= xlo ) then
          f = 1.0d0
       else
          s  = (xmax-x)/(xmax-xlo)

          sf =    10.0d0 * s**3 &
               & -15.0d0 * s**4 &
               & + 6.0d0 * s**5

          f = sf

       end if
    end if

  end subroutine switchoff

end module opnq_switching
