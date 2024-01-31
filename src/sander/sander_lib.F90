! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A collection of useful (but difficult-to-place) subroutines.
module sander_lib

    implicit none

contains

! These subroutines should have few-to-no dependencies, or we will quickly
! reach a point where we get cyclic dependencies that kills compilation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Checks if **s exist in a given inpcrd line and emit helpful messages
    subroutine check_inpcrd_overflow(line, periodic)

        implicit none

        ! Passed variables
        character(*), intent(in) :: line     ! line buffer that offends the parser
        logical, intent(in)      :: periodic ! Does our simulation have PBCs?

        ! Local variables
        integer      :: i ! counter

        ! My advice to you depends on whether you're running a simulation with
        ! PBCs or in implicit solvent/vacuum
        do i = 1, len_trim(line)
            if (line(i:i) .eq. '*' .and. periodic) then
                write (6, '(a)') '*s in the inpcrd file often indicate an overflow of &
                &the Fortran format used', &
                    'to store coordinates in the inpcrd/restart files. &
                &This often happens when', &
                    'coordinates are not wrapped into the center cell &
                &(when iwrap = 0) and some', &
                    'particles diffuse too far away. Try restarting from your &
                &last good restart', &
                    'file and setting iwrap=1 or using a NetCDF restart file &
                &format. See the', &
                    'Amber manual for details'
                ! Only print error message once. Bail out here.
                exit
            else if (line(i:i) .eq. '*' .and. .not. periodic) then
                write (6, '(a)') '*s in the inpcrd file often indicate an overflow of &
                &the Fortran format used', &
                    'to store coordinates in the inpcrd/restart files. &
                &This often happens when', &
                    'particles diffuse very far away from each other. Make &
                &sure you are removing', &
                    'center-of-mass translation (nscm /= 0) or check if you &
                &have multiple, mobile', &
                    'molecules that have diffused very far away from each other. &
                &This condition is', &
                    'highly unusual for non-periodic simulations.'
                ! Only print error message once. Bail out here.
                exit
            end if
        end do

        return

    end subroutine check_inpcrd_overflow

end module sander_lib
