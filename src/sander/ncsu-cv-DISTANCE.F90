! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! DISTANCE "subclass"
!

module ncsu_cv_DISTANCE

!=============================================================================

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

private :: do_PBC

private :: PBC_r2f
private :: PBC_f2r

private :: PBC_distance
private :: PBC_distance_d

!=============================================================================

contains

!=============================================================================

logical function do_PBC(cv)

   NCSU_USE_AFAILED

   use ncsu_constants, only : ZERO
   use ncsu_colvar_type
   use ncsu_sander_proxy, only : sander_ntb

   implicit none

   type(colvar_t), intent(in) :: cv

   do_PBC = .false.
   if (0.ne.sander_ntb()) then
      ncsu_assert(associated(cv%r))
      if (.not.cv%r(1).gt.ZERO) then
         do_PBC = .false.
      else
         do_PBC = .true.
      end if
   end if

end function do_PBC

!=============================================================================

subroutine PBC_r2f(r, f)

   use nblist, only : recip

   implicit none

   NCSU_REAL, intent(in) :: r(3)
   NCSU_REAL, intent(out) :: f(3)

   f(1) = r(1)*recip(1,1) + r(2)*recip(2,1) + r(3)*recip(3,1)
   f(2) = r(1)*recip(1,2) + r(2)*recip(2,2) + r(3)*recip(3,2)
   f(3) = r(1)*recip(1,3) + r(2)*recip(2,3) + r(3)*recip(3,3)

   f(1) = f(1) - floor(f(1))
   f(2) = f(2) - floor(f(2))
   f(3) = f(3) - floor(f(3))

end subroutine PBC_r2f

!=============================================================================

subroutine PBC_f2r(f, r, t1, t2, t3)

   use nblist, only : ucell

   implicit none

   NCSU_REAL, intent(in) :: f(3)
   NCSU_REAL, intent(out) :: r(3)

   integer, intent(in) :: t1, t2, t3

   r(1) = &
      (f(1) + t1)*ucell(1,1) + (f(2) + t2)*ucell(1,2) + (f(3) + t3)*ucell(1,3)
   r(2) = &
      (f(1) + t1)*ucell(2,1) + (f(2) + t2)*ucell(2,2) + (f(3) + t3)*ucell(2,3)
   r(3) = &
      (f(1) + t1)*ucell(3,1) + (f(2) + t2)*ucell(3,2) + (f(3) + t3)*ucell(3,3)

end subroutine PBC_f2r

!=============================================================================

function PBC_distance(r1, r2, d0) result(value)

   NCSU_USE_AFAILED

   use ncsu_constants, only : ZERO
   use ncsu_colvar_math, only : distance

   implicit none

   NCSU_REAL :: value
   NCSU_REAL, intent(in) :: r1(*), r2(*), d0

   NCSU_REAL :: f1(3), f2(3), d(-1:1,-1:1,-1:1)
   NCSU_REAL :: x1(3), x2(3), d_min

   integer :: i, j, k

   ncsu_assert(d0.gt.ZERO)

   ! get the fractionals

   call PBC_r2f(r1(1:3), f1)
   call PBC_r2f(r2(1:3), f2)

   ! wrap r1 to the primary cell

   call PBC_f2r(f1, x1, 0, 0, 0)
   d_min = ZERO

   ! wrap r2 to the primary/neighboring cells

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)
            d(i,j,k) = distance(x1, x2)
            if (d(i,j,k).lt.d_min) &
               d_min = d(i,j,k)
         end do
      end do
   end do

   value = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            value = value + exp((d_min - d(i,j,k))/d0)
         end do
      end do
   end do

   value = d_min - d0*log(value)

end function PBC_distance

!=============================================================================

subroutine PBC_distance_d(r1, r2, d0, d1, d2)

   NCSU_USE_AFAILED

   use ncsu_constants, only : ZERO
   use ncsu_colvar_math, only : distance

   implicit none

   NCSU_REAL, intent(in) :: r1(*), r2(*), d0
   NCSU_REAL, intent(out) :: d1(*), d2(*)

   NCSU_REAL :: f1(3), f2(3), d(-1:1,-1:1,-1:1)
   NCSU_REAL :: x1(3), x2(3), d_min, sum, ttt

   integer :: i, j, k, l

   ncsu_assert(d0.gt.ZERO)

   ! get the fractionals

   call PBC_r2f(r1(1:3), f1)
   call PBC_r2f(r2(1:3), f2)

   ! wrap r1 to the primary cell

   call PBC_f2r(f1, x1, 0, 0, 0)

   ! find d_min
   d_min = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)
            d(i,j,k) = distance(x1, x2)
            if (d(i,j,k).lt.d_min) &
               d_min = d(i,j,k)
         end do
      end do
   end do

   do l = 1, 3
      d1(l) = ZERO
      d2(l) = ZERO
   end do

   sum = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)

            ttt = exp((d_min - d(i,j,k))/d0)
            sum = sum + ttt
            ttt = ttt/d(i,j,k)

            do l = 1, 3
               d1(l) = d1(l) + (x1(l) - x2(l))*ttt
               d2(l) = d2(l) - (x1(l) - x2(l))*ttt
            end do
         end do
      end do
   end do

   do l = 1, 3
      d1(l) = d1(l)/sum
      d2(l) = d2(l)/sum
   end do

end subroutine PBC_distance_d

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

#  include "ncsu-mpi.h"

   integer :: a1, a2

   ncsu_assert(cv%type == COLVAR_DISTANCE)

   ncsu_assert(associated(cv%i))
   ncsu_assert(size(cv%i) == 2)

   ncsu_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   ncsu_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   ncsu_assert(a1 /= a2)

#ifdef MPI
   if (sanderrank.eq.0) then
#endif /* MPI */
   if (do_PBC(cv)) then
      value = PBC_distance(x(a1:a1 + 2), x(a2:a2 + 2), cv%r(1))
   else
      value = distance(x(a1:a1 + 2), x(a2:a2 + 2))
   end if ! do_PBC(cv)
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

   NCSU_REAL :: d1(3), d2(3)

   integer :: a1, a2

   ncsu_assert(cv%type == COLVAR_DISTANCE)

   ncsu_assert(associated(cv%i))
   ncsu_assert(size(cv%i) == 2)

   ncsu_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   ncsu_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   ncsu_assert(a1 /= a2)

   NCSU_MASTER_ONLY_BEGIN

   if (do_PBC(cv)) then
      call PBC_distance_d(x(a1:a1 + 2), x(a2:a2 + 2), cv%r(1), d1, d2)
   else
      call distance_d(x(a1:a1 + 2), x(a2:a2 + 2), d1, d2)
   end if ! do_PBC(cv)

   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2

   NCSU_MASTER_ONLY_END

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_colvar_utils
   use ncsu_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NCSU_REAL,      intent(in)    :: amass(*)

#  include "ncsu-mpi.h"

   integer :: error

   ncsu_assert(cv%type == COLVAR_DISTANCE)
   call check_i(cv%i, cvno, 'DISTANCE', 2)

   if (0.eq.sander_ntb()) &
      return

   if (.not.associated(cv%r)) then
      allocate(cv%r(1), stat = error)
      if (error.ne.0) &
         NCSU_OUT_OF_MEMORY
      cv%r(1) = ONE
   end if

   if (size(cv%r).ne.1) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NCSU_ERROR, 'CV #', cvno, &
            ' (DISTANCE) : unexpected number of reals (at most 1 is possible)'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   use ncsu_utils
   use ncsu_constants
   use ncsu_sander_proxy

   use ncsu_colvar_type
   use ncsu_colvar_utils

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   ncsu_assert(is_master())
   ncsu_assert(cv%type == COLVAR_DISTANCE)
   ncsu_assert(associated(cv%i))

   if (0.ne.sander_ntb()) then
      ncsu_assert(associated(cv%r))
      ncsu_assert(size(cv%r).eq.1)
      if (cv%r(1).gt.ZERO) then
         write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1), 3)//',a)') &
         NCSU_INFO, '    (between the closest PBC images smoothed with d0 = ', &
         cv%r(1), ' A)'
      end if ! cv%r(1).gt.ZERO
   end if

   call print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

end module ncsu_cv_DISTANCE
