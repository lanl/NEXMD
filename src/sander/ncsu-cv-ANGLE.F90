! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! ANGLE "subclass"
!

module ncsu_cv_ANGLE

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

!=============================================================================

contains

!=============================================================================

function colvar_value(cv, x) result(value)

#ifndef NCSU_DISABLE_ASSERT
   use ncsu_utils
   use ncsu_sander_proxy
#endif /* NCSU_DISABLE_ASSERT */

   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_colvar_math

   implicit none

   NCSU_REAL :: value

   type(colvar_t), intent(in) :: cv

   NCSU_REAL, intent(in) :: x(*)

#include "ncsu-mpi.h"

   integer :: a1, a2, a3

   ncsu_assert(cv%type.eq.COLVAR_ANGLE)

   ncsu_assert(associated(cv%i))
   ncsu_assert(size(cv%i).eq.3)

   ncsu_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   ncsu_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   ncsu_assert(cv%i(3) > 0 .and. cv%i(3) <= sander_natoms())
   a3 = 3*cv%i(3) - 2

   ncsu_assert(a1 /= a2 .and. a1 /= a3 .and. a2 /= a3)

#ifdef MPI
   if (sanderrank.eq.0) then
#endif /* MPI */
      value = angle(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2))
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

#ifndef NCSU_DISABLE_ASSERT
   use ncsu_utils
   use ncsu_sander_proxy
#endif /* NCSU_DISABLE_ASSERT */

   use ncsu_colvar_type
   use ncsu_colvar_math

   implicit none

   type(colvar_t), intent(in) :: cv

   NCSU_REAL, intent(in) :: x(*), fcv

   NCSU_REAL, intent(inout) :: f(*)

#  include "ncsu-mpi.h"

   NCSU_REAL :: d1(3), d2(3), d3(3)

   integer :: a1, a2, a3

   ncsu_assert(cv%type == COLVAR_ANGLE)

   ncsu_assert(associated(cv%i))
   ncsu_assert(size(cv%i) == 3)

   ncsu_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   ncsu_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   ncsu_assert(cv%i(3) > 0 .and. cv%i(3) <= sander_natoms())
   a3 = 3*cv%i(3) - 2

   ncsu_assert(a1 /= a2 .and. a1 /= a3 .and. a2 /= a3)

   NCSU_MASTER_ONLY_BEGIN

   call angle_d(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), d1, d2, d3)

   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
   f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3

   NCSU_MASTER_ONLY_END

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   NCSU_USE_AFAILED

   use ncsu_colvar_type
   use ncsu_colvar_utils

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NCSU_REAL,      intent(in)    :: amass(*)

   ncsu_assert(cv%type == COLVAR_ANGLE)
   call check_i(cv%i, cvno, 'ANGLE', 3)

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

#ifndef NCSU_DISABLE_ASSERT
   use ncsu_utils
   use ncsu_sander_proxy
#endif /* NCSU_DISABLE_ASSERT */

   use ncsu_colvar_type
   use ncsu_colvar_utils

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   ncsu_assert(is_master())
   ncsu_assert(cv%type == COLVAR_ANGLE)
   ncsu_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

end module ncsu_cv_ANGLE
