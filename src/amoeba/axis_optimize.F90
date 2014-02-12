!*******************************************************************************
!
! Module: axis_optimize_mod
!
! Description:  Module to optimize axis orientation in MPI runs, which has a
!               considerable impact on atom distribution for non-cubic systems.
!               The z axis is made the longest axis for internal calcs, but
!               this module provides a mechanism to correct the xyx coords for
!               output.  As it turns out, for single processor runs, slightly
!               worse performance is obtained, so we disable this function
!               for single processor code.
!              
!*******************************************************************************

module axis_optimize_mod

  implicit none

  ! The longest axis prior to flipping.  This axis will be flipped to the
  ! z axis using a righthand rule.  The values are the normal ordinals, ie.,
  ! x=1, y=2, z=3

  integer, save :: longest_axis

  ! If these ordinals are used for coordinates, velocities, etc., the axes
  ! will effectively be flipped back to the original values.

  integer, save :: axis_flipback_ords(3)

contains

!*******************************************************************************
!
! Subroutine:  setup_axis_opt
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_axis_opt(x, y, z)

  implicit none

! Formal arguments:

  double precision      :: x, y, z

! Local variables:

  double precision      :: max_len

  ! If there is no difference between the length of the axes, then set up
  ! the flip control variables to minimize the work required.

  if (x .eq. y .and. y .eq. z) then

    longest_axis = 3
    axis_flipback_ords(1) = 1
    axis_flipback_ords(2) = 2
    axis_flipback_ords(3) = 3

  else

    max_len = max(x, y, z)

    if (max_len .eq. x) then
    
      longest_axis = 1
      axis_flipback_ords(1) = 3
      axis_flipback_ords(2) = 1
      axis_flipback_ords(3) = 2

    else if (max_len .eq. y) then
    
      longest_axis = 2
      axis_flipback_ords(1) = 2
      axis_flipback_ords(2) = 3
      axis_flipback_ords(3) = 1

    else
    
      longest_axis = 3
      axis_flipback_ords(1) = 1
      axis_flipback_ords(2) = 2
      axis_flipback_ords(3) = 3

    end if

  end if

  return

end subroutine setup_axis_opt

!*******************************************************************************
!
! Subroutine:  axes_flip
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine axes_flip(x, y, z)

  implicit none

! Formal arguments:

  double precision      :: x, y, z

! Local variables:

  double precision      :: tmp

  if (longest_axis .eq. 1) then

    tmp = x
    x = y
    y = z
    z = tmp     ! tmp is x

  else if (longest_axis .eq. 2) then

    tmp = y
    y = x
    x = z
    z = tmp     ! tmp is y

  end if

  ! No flip needed if longest axis is z (ie., longest_axis .eq. 3)

  return

end subroutine axes_flip

!*******************************************************************************
!
! Subroutine:  int_axes_flip
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine int_axes_flip(x, y, z)

  implicit none

! Formal arguments:

  integer :: x, y, z

! Local variables:

  integer :: tmp

  if (longest_axis .eq. 1) then

    tmp = x
    x = y
    y = z
    z = tmp     ! tmp is x

  else if (longest_axis .eq. 2) then

    tmp = y
    y = x
    x = z
    z = tmp     ! tmp is y

  end if

  ! No flip needed if longest axis is z (ie., longest_axis .eq. 3)

  return

end subroutine int_axes_flip

end module axis_optimize_mod
