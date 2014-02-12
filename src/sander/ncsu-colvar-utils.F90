#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar_utils

implicit none

private

!=============================================================================

public :: check_i
public :: print_i
public :: print_pca

public :: com_check_i
public :: com_print_i

public :: com_init_weights

public :: group_com
public :: group_com_d

!=============================================================================

contains

!=============================================================================

subroutine check_i(cvi, cvno, cvtype, expected_isize)

   use ncsu_utils
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: cvno
   character(*), intent(in) :: cvtype
   integer, optional, intent(in) :: expected_isize

   integer :: a, b

#include "ncsu-mpi.h"

   ncsu_assert(cvno > 0)

   if (.not.associated(cvi)) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : no integers found'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (present(expected_isize)) then
      if (size(cvi) /= expected_isize) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
            (size(cvi))//',a,'//pfmt(expected_isize)//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : unexpected number of integers (', &
               size(cvi), ' instead of ', expected_isize, ')'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if ! size(cvi) /= isize
   end if ! present(expected_isize)

   do a = 1, size(cvi)
      if (cvi(a) < 1 .or. cvi(a) > sander_natoms()) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a,'//pfmt &
               (cvi(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : integer #', a, ' (', cvi(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates

   do a = 1, size(cvi)
      do b = a + 1, size(cvi)
         if (cvi(a) == cvi(b)) then
            NCSU_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
                  (a)//',a,'//pfmt(b)//',a,'//pfmt(cvi(a))//',a/)') &
                  NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
                  ') : integers #', a, ' and #', b, ' are equal (', cvi(a), ')'
            NCSU_MASTER_ONLY_END
            call terminate()
         end if ! cvi(a) == cvi(b)
      end do
   end do

end subroutine check_i

!=============================================================================

subroutine print_i(cvi, lun)

   use ncsu_utils
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun

   integer :: a
   character(4) :: aname

   ncsu_assert(is_master())
   ncsu_assert(associated(cvi))
   ncsu_assert(size(cvi) > 0)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, '  atoms = ('

   do a = 1, size(cvi)

      ncsu_assert(cvi(a) > 0 .and. cvi(a) <= sander_natoms())
      aname = sander_atom_name(cvi(a))

      write (unit = lun, fmt = '('//pfmt(cvi(a))//',a,a,a)', advance = 'NO') &
         cvi(a), ' [', trim(aname), ']'

      if (a == size(cvi)) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(a, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,a)', advance = 'NO') &
            ',', NCSU_INFO, '          '
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

   end do

end subroutine print_i

!=============================================================================

subroutine com_check_i(cvi, cvno, cvtype, ngroups)

   use ncsu_utils
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: cvno
   character(*), intent(in) :: cvtype
   integer, intent(out) :: ngroups

   integer :: a, a0

#  include "ncsu-mpi.h"

   ncsu_assert(cvno > 0)

   if (.not.associated(cvi)) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : no integers found'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (.not.size(cvi).ge.3) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : too few integers'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   a0 = 1
   ngroups = 0

   do a = 1, size(cvi)
      if (cvi(a).eq.0) then
         ngroups = ngroups + 1
         if (a.eq.a0) then
            NCSU_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a/)') &
                  NCSU_ERROR, 'CV #', cvno, &
                  ' (', cvtype, ') : unexpected zero (integer #', a, ')'
            NCSU_MASTER_ONLY_END
            call terminate()
         end if ! a.eq.a0
         a0 = a + 1
      else if (cvi(a).lt.1.or.cvi(a).gt.sander_natoms()) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a,'//pfmt &
               (cvi(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : integer #', a, ' (', cvi(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if
   end do

   if (cvi(size(cvi)).gt.0) &
      ngroups = ngroups + 1

end subroutine com_check_i

!=============================================================================

subroutine com_print_i(cvi, lun)

   use ncsu_utils
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun

   integer :: a, g, c, ncvi
   character(4) :: aname

   ncsu_assert(is_master())
   ncsu_assert(associated(cvi))
   ncsu_assert(size(cvi).gt.0)

   c = 1
   g = 1

   ncvi = size(cvi)

   if (cvi(ncvi).eq.0) &
      ncvi = ncvi - 1

   write (unit = lun, fmt = '(a,a,i1,a)', advance = 'NO') &
      NCSU_INFO, '  group #', g, ' = ('

   do a = 1, ncvi
      if (cvi(a).eq.0) then
         g = g + 1
         c = 1

         write (unit = lun, fmt = '(a,a,i1,a)', advance = 'NO') &
            NCSU_INFO, '  group #', g, ' = ('
      else ! cvi(a).ne.0
         ncsu_assert(cvi(a).gt.0.and.cvi(a).le.sander_natoms())
         aname = sander_atom_name(cvi(a))

      write (unit = lun, fmt = '('//pfmt(cvi(a))//',a,a,a)', advance = 'NO') &
         cvi(a), ' [', trim(aname), ']'

         if (next_is_atom(a)) then
            if (mod(c, 5) == 0) then
               write (unit = lun, fmt = '(a,/a,a)', advance = 'NO') &
                  ',', NCSU_INFO, '              '
            else
               write (unit = lun, fmt = '(a)', advance = 'NO') ', '
            end if
         else
            write (unit = lun, fmt = '(a)') ')'
         end if ! next_is_atom(a)
         c = c + 1
      end if ! cvi(a).eq.0
   end do

contains

logical function next_is_atom(ii)

   implicit none

   integer, intent(in) :: ii

   next_is_atom = .false.
   if (ii.lt.size(cvi)) &
      next_is_atom = cvi(ii + 1).gt.0

end function next_is_atom

end subroutine com_print_i

!=============================================================================

! cvr holds atomic masses of different groups padded by zeros
! com_init_weights() divides the individual masses by total group mass
subroutine com_init_weights(cvr)

   use ncsu_utils
   use ncsu_constants

   implicit none

   NCSU_REAL, pointer :: cvr(:)

   integer :: a, a0, n
   NCSU_REAL :: mass

   ncsu_assert(associated(cvr))
   ncsu_assert(size(cvr).gt.2)

   mass = ZERO
   a0 = 1

   do a = 1, size(cvr)
      mass = mass + cvr(a)
      if (abs(cvr(a)).lt.1.0d-10.or.a.eq.size(cvr)) then
         ncsu_assert(a.ge.a0)
         ncsu_assert(mass.gt.ZERO)
         do n = a0, a
            cvr(n) = cvr(n)/mass
         end do
         a0 = a + 1
         mass = ZERO
      end if
   end do

end subroutine com_init_weights

!=============================================================================

! see ncsu-cv-COM_*.f
subroutine group_com(cv, x, pos, cm)

   use ncsu_constants, only : ZERO
   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv
   NCSU_REAL, intent(in) :: x(*)
   integer, intent(inout) :: pos
   NCSU_REAL, intent(out) :: cm(3)

   integer :: a

   cm = ZERO

   do while(pos.le.size(cv%i))
      a = cv%i(pos)
      if (a.gt.0) then
         a = 3*a - 2
         cm = cm + cv%r(pos)*x(a:a + 2)
      else
         exit
      end if
      pos = pos + 1
   end do

   pos = pos + 1 ! skip zero

end subroutine group_com

!=============================================================================

! see ncsu-cv-COM_*.f
subroutine group_com_d(cv, f, d, pos)

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv
   NCSU_REAL, intent(inout) :: f(*)
   NCSU_REAL, intent(in) :: d(3)
   integer, intent(inout) :: pos

   integer :: a

   do while(pos.le.size(cv%i))
      a = cv%i(pos)
      if (a.gt.0) then
         a = 3*a - 2
         f(a:a + 2) = f(a:a + 2) + cv%r(pos)*d(1:3)
      else
         exit
      end if
      pos = pos + 1
   end do

   pos = pos + 1 ! skip zero

end subroutine group_com_d

!=============================================================================

subroutine print_pca(cvi, lun)

   use ncsu_utils
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun


   ncsu_assert(is_master())
   ncsu_assert(associated(cvi))
   ncsu_assert(size(cvi) == 3)

   write (unit = lun, fmt = '(a,a,i5)') NCSU_INFO, '     solute part total atom number = ', cvi(1)
   write (unit = lun, fmt = '(a,a,i5)') NCSU_INFO, '        ref part total atom number = ', cvi(2)
   write (unit = lun, fmt = '(a,a,i5)') NCSU_INFO, '        pca part total atom number = ', cvi(3)   

end subroutine print_pca

!=============================================================================
end module ncsu_colvar_utils
