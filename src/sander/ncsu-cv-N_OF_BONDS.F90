! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! cv%r = (d0)
! cv%i = (i1, j1, i2, j2, ... ) ! list of pairs
!
! value = \sum_{pairs}\frac{1 - (r_i - r_j)^6/d0^6}{1 - (r_i - r_j)^12/d0^12}
!

module ncsu_cv_N_OF_BONDS

!=============================================================================

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

private :: f6_12_value
private :: f6_12_derivative

private :: partition

!=============================================================================

contains

!=============================================================================

NCSU_REAL pure function f6_12_value(x)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, intent(in) :: x

   f6_12_value = ONE/(ONE + x**3)

end function f6_12_value

!=============================================================================

NCSU_REAL pure function f6_12_derivative(x)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, parameter :: MINUS_THREE = -3.000000000000000000000D0 ! NCSU_TO_REAL(-3)

   NCSU_REAL, intent(in) :: x

   NCSU_REAL :: x2, tmp

   x2 = x*x
   tmp = ONE + x*x2

   f6_12_derivative = MINUS_THREE*x2/(tmp*tmp)

end function f6_12_derivative

!=============================================================================

subroutine partition(cv, first, last)

   NCSU_USE_AFAILED

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(out) :: first, last

#  ifdef MPI
#     include "ncsu-mpi.h"
      integer :: tmp
#  endif /* MPI */

   ncsu_assert(cv%type == COLVAR_N_OF_BONDS)

   ncsu_assert(associated(cv%i))
   ncsu_assert(mod(size(cv%i), 2).eq.0)

#  ifdef MPI
   tmp = (size(cv%i)/2)/sandersize
   if (tmp.gt.0) then
      if (sanderrank.eq.(sandersize - 1)) then
         first = 2*tmp*sanderrank + 1
         last = size(cv%i) - 1
      else
         first = 2*tmp*sanderrank + 1
         last = 2*(sanderrank + 1)*tmp - 1
      end if
   else
      if (sanderrank.eq.(sandersize - 1)) then
         first = 1
         last = size(cv%i) - 1
      else
         first = 1
         last = 0
      end if
   end if
#  else
   first = 1
   last = size(cv%i) - 1
#  endif /* MPI */

end subroutine partition

!=============================================================================

function colvar_value(cv, x) result(value)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   NCSU_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NCSU_REAL, intent(in) :: x(*)

#  ifdef MPI
#  include "ncsu-mpi.h"
   integer :: error
   NCSU_REAL :: accu
#  endif /* MPI */

   integer :: first, last
   integer :: i, i3, j3

   NCSU_REAL :: r2

   ncsu_assert(cv%type == COLVAR_N_OF_BONDS)

   ncsu_assert(associated(cv%i))
   ncsu_assert(associated(cv%r))

   ncsu_assert(mod(size(cv%i), 2).eq.0)

   value = ZERO

   call partition(cv, first, last)

   do i = first, last, 2

      i3 = 3*cv%i(i) - 2
      j3 = 3*cv%i(i + 1) - 2

      r2 = (x(i3 + 0) - x(j3 + 0))**2 &
         + (x(i3 + 1) - x(j3 + 1))**2 &
         + (x(i3 + 2) - x(j3 + 2))**2

      r2 = r2/cv%r(1)**2

      value = value + f6_12_value(r2)

   end do

#  ifdef MPI
   call mpi_reduce(value, accu, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, 0, commsander, error)
   ncsu_assert(error.eq.0)
   value = accu
#  endif /* MPI */

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   NCSU_REAL, intent(in) :: x(*), fcv
   NCSU_REAL, intent(inout) :: f(*)

   integer   :: first, last, i, i3, j3
   NCSU_REAL :: dx(3), r2, tmp

   ncsu_assert(cv%type == COLVAR_N_OF_BONDS)

   ncsu_assert(associated(cv%i))
   ncsu_assert(associated(cv%r))

   call partition(cv, first, last)

   do i = first, last, 2

      i3 = 3*cv%i(i) - 2
      j3 = 3*cv%i(i + 1) - 2

      dx(1) = x(i3 + 0) - x(j3 + 0)
      dx(2) = x(i3 + 1) - x(j3 + 1)
      dx(3) = x(i3 + 2) - x(j3 + 2)

      r2 = (dx(1)**2 + dx(2)**2 + dx(3)**2)/cv%r(1)**2
      tmp = (2*fcv/cv%r(1)**2)*f6_12_derivative(r2)

      dx(1) = tmp*dx(1)
      f(i3 + 0) = f(i3 + 0) + dx(1)
      f(j3 + 0) = f(j3 + 0) - dx(1)

      dx(2) = tmp*dx(2)
      f(i3 + 1) = f(i3 + 1) + dx(2)
      f(j3 + 1) = f(j3 + 1) - dx(2)

      dx(3) = tmp*dx(3)
      f(i3 + 2) = f(i3 + 2) + dx(3)
      f(j3 + 2) = f(j3 + 2) - dx(3)

   end do

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NCSU_REAL,      intent(in)    :: amass(*)

#  include "ncsu-mpi.h"

  integer :: a

   ncsu_assert(cv%type == COLVAR_N_OF_BONDS)

   if (.not.associated(cv%i)) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : no integers'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (size(cv%i).lt.2) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : too few integers'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.3

   if (mod(size(cv%i), 2).ne.0) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, &
            ' (N_OF_BONDS) : number of integers is odd'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! sep.eq.size(cv%i)

   do a = 1, size(cv%i) - 1, 2
      if (cv%i(a).lt.1.or.cv%i(a).gt.sander_natoms()) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
              (a)//',a,'//pfmt(cv%i(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (N_OF_BONDS) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if
      if (cv%i(a + 1).lt.1.or.cv%i(a + 1).gt.sander_natoms()) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
               (a + 1)//',a,'//pfmt(cv%i(a + 1))//',a,'//pfmt &
               (sander_natoms())//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (N_OF_BONDS) : integer #', a + 1, ' (', cv%i(a + 1), &
               ') is out of range [1, ', sander_natoms(), ']'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if
      if (cv%i(a).eq.cv%i(a + 1)) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
               (a)//',a,'//pfmt(a + 1)//',a,'//pfmt(cv%i(a))//',a/)') &
               NCSU_WARNING, 'CV #', cvno, &
               ' (N_OF_BONDS) : integers #', a, ' and ', a + 1, &
               ' are equal (', cv%i(a), ')'
         NCSU_MASTER_ONLY_END
      end if
   end do

   if (.not.associated(cv%r)) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : no reals found'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (size(cv%r).ne.1) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : number of reals is not 1'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (cv%r(1).le.ZERO) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : r(1).le.0.0D0 is .true.'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   use ncsu_utils
   use ncsu_colvar_type
   use ncsu_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: a
   character(4) :: aname

   ncsu_assert(cv%type == COLVAR_N_OF_BONDS)

   ncsu_assert(associated(cv%i))
   ncsu_assert(associated(cv%r))

   ncsu_assert(mod(size(cv%i), 2).eq.0)

   write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1), 3)//')') &
      NCSU_INFO, '    d0 = ', cv%r(1)
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, ' pairs = ('

   do a = 1, size(cv%i)

      ncsu_assert(cv%i(a).gt.0.and.cv%i(a).le.sander_natoms())
      aname = sander_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a.eq.size(cv%i)) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(a, 4).eq.0) then
         write (unit = lun, fmt = '(a,/a,10x)', advance = 'NO') ',', NCSU_INFO
      else if (mod(a, 2).eq.1) then
         write (unit = lun, fmt = '(a)', advance = 'NO') ' <=> '
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

   end do

end subroutine print_details

!=============================================================================

end module ncsu_cv_N_OF_BONDS
