! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_umbrella

#ifndef NCSU_NO_NETCDF
use netcdf, only : nf90_double, nf90_float
#endif /* NCSU_NO_NETCDF */

use ncsu_constants, only : ZERO, ONE, TWO, THREE, FOUR

implicit none

!============================================================================+
!             + x +   P U B L I C    I N T E R F A C E   + x +               !
!============================================================================+

integer, public, parameter :: UMBRELLA_MIN_EXTENT   = 5
integer, public, parameter :: UMBRELLA_MAX_NEXTENTS = 4

type, public :: umbrella_t

   private

   integer :: nextents, extents(UMBRELLA_MAX_NEXTENTS)
   logical :: periodicity(UMBRELLA_MAX_NEXTENTS)

   NCSU_REAL :: origin(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: spacing(UMBRELLA_MAX_NEXTENTS)

   NCSU_REAL, pointer :: coeffs(:) ! row-major

#ifndef NCSU_DISABLE_ASSERT
   logical :: inited = .false.
#endif /* NCSU_DISABLE_ASSERT */
end type umbrella_t

public :: umbrella_init
public :: umbrella_fini

public :: umbrella_nextents
public :: umbrella_extent
public :: umbrella_origin
public :: umbrella_spacing
public :: umbrella_periodicity

public :: umbrella_hill
public :: umbrella_swap
public :: umbrella_transfer

public :: umbrella_eval_v
public :: umbrella_eval_vdv

#ifdef MPI
public :: umbrella_bcast
public :: umbrella_send_coeffs
public :: umbrella_recv_coeffs
#endif /* MPI */

#ifndef NCSU_NO_NETCDF
public :: umbrella_load
public :: umbrella_save
#  ifdef NCSU_REAL_IS_DOUBLE
integer, private, parameter :: coeffs_type = nf90_double
#  else
integer, private, parameter :: coeffs_type = nf90_float ! not tested
#  endif /* NCSU_REAL_IS_DOUBLE */
#endif /* NCSU_NO_NETCDF */

!============================================================================+
!                      + x +   D E T A I L S   + x +                         !
!============================================================================+

private :: m4_v
private :: m4_vdv

NCSU_REAL, private, parameter :: MINUS_ONE = -ONE
NCSU_REAL, private, parameter :: ONE_THIRD = ONE/THREE
NCSU_REAL, private, parameter :: MINUS_ONE_THIRD = -ONE_THIRD
NCSU_REAL, private, parameter :: ONE_SIXTH = ONE/6
NCSU_REAL, private, parameter :: TWO_THIRD = TWO/THREE
NCSU_REAL, private, parameter :: HALF = ONE/TWO
NCSU_REAL, private, parameter :: MINUS_HALF = -HALF

!=============================================================================

contains

!=============================================================================

subroutine umbrella_init(this, &
   nextents, extents, origin, spacing, periodicity)

   use ncsu_utils

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: nextents
   integer, intent(in) :: extents(*)

   NCSU_REAL, intent(in) :: origin(*)
   NCSU_REAL, intent(in) :: spacing(*)

   logical, intent(in) :: periodicity(*)

   integer :: n, ncoeffs, ierr, tmp

   ncsu_assert(.not.this%inited)

   ncsu_assert(nextents.ge.1)
   ncsu_assert(nextents.le.UMBRELLA_MAX_NEXTENTS)

   this%nextents = nextents

   ncoeffs = 1
   do n = 1, nextents
      ncsu_assert(extents(n).ge.UMBRELLA_MIN_EXTENT)
      ncsu_assert(spacing(n).gt.ZERO)

      this%origin(n) = origin(n)
      this%spacing(n) = spacing(n)

      this%extents(n) = extents(n)
      this%periodicity(n) = periodicity(n)

      tmp = ncoeffs
      ncoeffs = ncoeffs*extents(n)
      if (ncoeffs/extents(n).ne.tmp) &
         call fatal('overflow in umbrella_init()')
   end do

   allocate(this%coeffs(ncoeffs), stat = ierr)
   if (ierr /= 0) &
      NCSU_OUT_OF_MEMORY

   this%coeffs(1:ncoeffs) = ZERO

#ifndef NCSU_DISABLE_ASSERT
   this%inited = .true.
#endif /* NCSU_DISABLE_ASSERT */

end subroutine umbrella_init

!=============================================================================

subroutine umbrella_fini(this)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this

   ncsu_assert(this%inited)
   ncsu_assert(associated(this%coeffs))

   deallocate(this%coeffs)

#ifndef NCSU_DISABLE_ASSERT
   this%inited = .false.
#endif /* NCSU_DISABLE_ASSERT */

end subroutine umbrella_fini

!=============================================================================

integer function umbrella_nextents(this)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this

   ncsu_assert(this%inited)
   umbrella_nextents = this%nextents

end function umbrella_nextents

!=============================================================================

integer function umbrella_extent(this, n)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   ncsu_assert(this%inited)
   ncsu_assert(n.ge.1)
   ncsu_assert(n.le.this%nextents)

   umbrella_extent = this%extents(n)

end function umbrella_extent

!=============================================================================

NCSU_REAL function umbrella_origin(this, n)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   ncsu_assert(this%inited)
   ncsu_assert(n.ge.1)
   ncsu_assert(n.le.this%nextents)

   umbrella_origin = this%origin(n)

end function umbrella_origin

!=============================================================================

NCSU_REAL function umbrella_spacing(this, n)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   ncsu_assert(this%inited)
   ncsu_assert(n.ge.1)
   ncsu_assert(n.le.this%nextents)

   umbrella_spacing = this%spacing(n)

end function umbrella_spacing

!=============================================================================

logical function umbrella_periodicity(this, n)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   ncsu_assert(this%inited)
   ncsu_assert(n.ge.1)
   ncsu_assert(n.le.this%nextents)

   umbrella_periodicity = this%periodicity(n)

end function umbrella_periodicity

!=============================================================================

subroutine umbrella_hill(this, x, alt)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this

   NCSU_REAL, intent(in) :: x(*)
   NCSU_REAL, intent(in) :: alt

   NCSU_REAL :: xs, accum, hill_v(4, UMBRELLA_MAX_NEXTENTS)
   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   ncsu_assert(this%inited)
   ncsu_assert(this%nextents.ge.1)
   ncsu_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   ncsu_assert(associated(this%coeffs))

   do i = 1, this%nextents
      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         hill_v(j, i) = hill(xs - NCSU_TO_REAL(j + gc(i) - 2))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      accum = alt

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         accum = accum*hill_v(iand(p, 3) + 1, i)

         p = ishft(p, -2)

      end do ! loop over i

      this%coeffs(o) = this%coeffs(o) + accum

   end do outer

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

! "hill" centered at 0 (-2 < x < 2)
pure NCSU_REAL function hill(x)

   implicit none

   NCSU_REAL, intent(in) :: x

   NCSU_REAL, parameter :: SCALE = 48.000000000000000000000D0 / 41.000000000000000000000D0!NCSU_TO_REAL(48)/NCSU_TO_REAL(41)

   NCSU_REAL :: x2

   x2 = (HALF*x)**2

!   if (x2.lt.ONE) then
    hill = SCALE*(x2 - ONE)**2
!   else
!      hill = ZERO
!   end if

end function hill

!-----------------------------------------------------------------------------

end subroutine umbrella_hill

!=============================================================================

! sets 'this' from 'other' [does *NOT* interpolate rigorously]
subroutine umbrella_transfer(this, other)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this
   type(umbrella_t), intent(in) :: other

   NCSU_REAL :: pos(UMBRELLA_MAX_NEXTENTS)

   integer :: i, j, p, n, o, ncoeffs

   ncsu_assert(this%inited)
   ncsu_assert(other%inited)

   ncsu_assert(this%nextents.ge.1)
   ncsu_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)

   ncsu_assert(other%nextents.ge.1)
   ncsu_assert(other%nextents.le.UMBRELLA_MAX_NEXTENTS)

   ncsu_assert(this%nextents.eq.other%nextents)

   ncoeffs = 1
   do n = 1, this%nextents
      ncoeffs = ncoeffs*this%extents(n)
   end do

   ! this is *not* interpolation
   do n = 1, ncoeffs
      p = n
      o = 0
      do i = 1, this%nextents
         j = mod(p, this%extents(i))
         pos(i) = this%origin(i) + j*this%spacing(i)
         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if
         p = p/this%extents(i)
      end do
      this%coeffs(o) = umbrella_eval_v(other, pos)
   end do

end subroutine umbrella_transfer

!=============================================================================

! swaps without any checks
subroutine umbrella_swap(this, other)

   implicit none

   type(umbrella_t), intent(inout) :: this, other

   type(umbrella_t) :: tmp

   tmp%nextents = this%nextents
   this%nextents = other%nextents
   other%nextents = tmp%nextents

   tmp%extents = this%extents
   this%extents = other%extents
   other%extents = tmp%extents

   tmp%periodicity = this%periodicity
   this%periodicity = other%periodicity
   other%periodicity = tmp%periodicity

   tmp%origin = this%origin
   this%origin = other%origin
   other%origin = tmp%origin

   tmp%spacing = this%spacing
   this%spacing = other%spacing
   other%spacing = tmp%spacing

#ifndef NCSU_DISABLE_ASSERT
   tmp%inited = this%inited
   this%inited = other%inited
   other%inited = tmp%inited
#endif /* NCSU_DISABLE_ASSERT */

   tmp%coeffs => this%coeffs
   this%coeffs => other%coeffs
   other%coeffs => tmp%coeffs

end subroutine umbrella_swap

!=============================================================================

NCSU_REAL function umbrella_eval_v(this, x) result(v)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL :: xs, m4v_prod, m4v(4, UMBRELLA_MAX_NEXTENTS)
   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   ncsu_assert(this%inited)
   ncsu_assert(this%nextents.ge.1)
   ncsu_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   ncsu_assert(associated(this%coeffs))

   v = ZERO

   do i = 1, this%nextents
      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         m4v(j, i) = m4_v(xs - NCSU_TO_REAL(j + gc(i) - 2))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      m4v_prod = ONE

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         m4v_prod = m4v_prod*m4v(iand(p, 3) + 1, i)

         p = ishft(p, -2)

      end do ! loop over i

      v = v + this%coeffs(o)*m4v_prod

   end do outer

end function umbrella_eval_v

!=============================================================================

subroutine umbrella_eval_vdv(this, x, v, dv)

   NCSU_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this

   NCSU_REAL, intent(in)  :: x(*)
   NCSU_REAL, intent(out) :: v
   NCSU_REAL, intent(out) :: dv(*)

   NCSU_REAL :: xs, m4v_prod, m4dv_prod(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: m4v(4, UMBRELLA_MAX_NEXTENTS), m4dv(4, UMBRELLA_MAX_NEXTENTS)

   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   ncsu_assert(this%inited)
   ncsu_assert(this%nextents.ge.1)
   ncsu_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   ncsu_assert(associated(this%coeffs))

   v = ZERO

   do i = 1, this%nextents
      dv(i) = ZERO

      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         call m4_vdv(xs - NCSU_TO_REAL(j + gc(i) - 2), m4v(j, i), m4dv(j, i))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      m4v_prod = ONE

      do i = 1, this%nextents
         m4dv_prod(i) = ONE
      end do

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         m4v_prod = m4v_prod*m4v(iand(p, 3) + 1, i)

         do j = 1, i - 1
            m4dv_prod(j) = m4dv_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         m4dv_prod(i) = m4dv_prod(i)*m4dv(iand(p, 3) + 1, i)

         do j = i + 1, this%nextents
            m4dv_prod(j) = m4dv_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         p = ishft(p, -2)

      end do ! loop over i

      v = v + this%coeffs(o)*m4v_prod

      do i = 1, this%nextents
         dv(i) = dv(i) + this%coeffs(o)*m4dv_prod(i)
      end do

   end do outer

   do i = 1, this%nextents
      dv(i) = dv(i)/this%spacing(i)
   end do

end subroutine umbrella_eval_vdv

!=============================================================================

#ifdef MPI

subroutine umbrella_bcast(this, comm, root) ! inefficient && memory-leak prone

   use ncsu_utils

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: comm, root

#  include "ncsu-mpi.h"

   integer :: p(UMBRELLA_MAX_NEXTENTS)

   integer :: n, error, commrank, commsize, ncoeffs

   ncsu_assert(comm.ne.mpi_comm_null)

   call mpi_comm_rank(comm, commrank, error)
   ncsu_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   ncsu_assert(error.eq.0)

   ncsu_assert(root.ge.0)
   ncsu_assert(root.lt.commsize)

#  ifndef NCSU_DISABLE_ASSERT
   if (commrank.eq.root) then
      ncsu_assert(this%inited)
   else
      ncsu_assert(.not.this%inited)
   end if
#  endif /* NCSU_DISABLE_ASSERT */

   call mpi_bcast(this%nextents, 1, MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

   call mpi_bcast(this%extents, this%nextents, MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

   do n = 1, this%nextents
      if (this%periodicity(n)) then
         p(n) = 1
      else
         p(n) = 0
      end if
   end do

   call mpi_bcast(p, this%nextents, MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

   do n = 1, this%nextents
      this%periodicity(n) = (p(n).eq.1)
   end do

   call mpi_bcast(this%origin, this%nextents, MPI_DOUBLE_PRECISION, &
      root, comm, error)
   ncsu_assert(error.eq.0)

   call mpi_bcast(this%spacing, this%nextents, MPI_DOUBLE_PRECISION, &
      root, comm, error)
   ncsu_assert(error.eq.0)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_bcast()')
   end do

   if (commrank.ne.root) then
      allocate(this%coeffs(ncoeffs), stat = error)
      if (error.ne.0) &
         NCSU_OUT_OF_MEMORY
   end if

   call mpi_bcast(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, root, comm, error)
   ncsu_assert(error.eq.0)

#  ifndef NCSU_DISABLE_ASSERT
   this%inited = .true.
#  endif /* NCSU_DISABLE_ASSERT */

end subroutine umbrella_bcast

!=============================================================================

subroutine umbrella_send_coeffs(this, dst, comm) ! unsafe

   use ncsu_utils

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: dst, comm

#  include "ncsu-mpi.h"

   integer :: n, ncoeffs, error

   ncsu_assert(this%inited)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_send_coeffs()')
   end do

   call mpi_send(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, &
                 dst, this%nextents, comm, error)

   ncsu_assert(error.eq.0)

end subroutine umbrella_send_coeffs

!=============================================================================

subroutine umbrella_recv_coeffs(this, src, comm) ! unsafe

   use ncsu_utils

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: src, comm

#  include "ncsu-mpi.h"

   integer :: n, ncoeffs, error

   ncsu_assert(this%inited)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_recv_coeffs()')
   end do

   call mpi_recv(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, &
                 src, this%nextents, comm, MPI_STATUS_IGNORE, error)

   ncsu_assert(error.eq.0)

end subroutine umbrella_recv_coeffs

#endif /* MPI */

!=============================================================================

#ifndef NCSU_NO_NETCDF
subroutine umbrella_load(this, filename)

   use netcdf

   use ncsu_utils
   use ncsu_constants, only : ERR_UNIT
   use ncsu_sander_proxy

   implicit none

   type(umbrella_t), intent(inout) :: this
   character(*), intent(in) :: filename

   integer :: n, rc, setid, dimid, varid, ncoeffs, tmp, nextents, coeffs_len
   integer :: periodicity(UMBRELLA_MAX_NEXTENTS) ! netCDF doesn't like logical

   ncsu_assert(.not.this%inited)

   rc = nf90_open(filename, nf90_nowrite, setid)
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'nextents', nextents)
   call check_rc()

   if (nextents.lt.1.or.nextents.gt.UMBRELLA_MAX_NEXTENTS) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,'//pfmt(nextents)//',a/)') &
         NCSU_ERROR, 'umbrella_load(filename=''', trim(filename), &
         ''') : nextents is out of range (', nextents, ')'
      call terminate()
   end if

   this%nextents = nextents

   rc = nf90_get_att(setid, nf90_global, 'extents', this%extents(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'periodicity', periodicity(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'origin', this%origin(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'spacing', this%spacing(1:nextents))
   call check_rc()

   ncoeffs = 1
   do n = 1, nextents
      tmp = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.tmp) &
         call fatal('overflow in umbrella_load()')
      if (this%extents(n).lt.UMBRELLA_MIN_EXTENT) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a,i1,a,i1/)') NCSU_ERROR, &
            'umbrella_load(filename=''', trim(filename), &
            ''') : extents(', n, ').lt.', UMBRELLA_MIN_EXTENT
         call terminate()
      end if
      if (this%spacing(n).le.ZERO) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a,i1,a/)') NCSU_ERROR, &
            'umbrella_load(filename=''', trim(filename), &
            ''') : spacing(', n, ').le.0'
         call terminate()
      end if
      if (periodicity(n).ne.0) then
         this%periodicity(n) = .true.
      else
         this%periodicity(n) = .false.
      end if
   end do

   rc = nf90_inq_dimid(setid, 'row-major', dimid)
   call check_rc()

   rc = nf90_inquire_dimension(setid, dimid, len = coeffs_len)
   call check_rc()

   if (coeffs_len.ne.ncoeffs) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NCSU_ERROR, &
         'umbrella_load(filename=''', trim(filename), &
         ''') : number of ''coeffs'' is wrong'
      call terminate()
   end if

   rc = nf90_inq_varid(setid, 'coeffs', varid)
   call check_rc()

   allocate(this%coeffs(ncoeffs), stat = rc)
   if (rc.ne.0) &
      NCSU_OUT_OF_MEMORY

   rc = nf90_get_var(setid, varid, this%coeffs)
   call check_rc()

   rc = nf90_close(setid)
   call check_rc()

#ifndef NCSU_DISABLE_ASSERT
   this%inited = .true.
#endif /* NCSU_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine check_rc()

   implicit none

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') NCSU_ERROR, &
         'umbrella_load(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_rc

!-----------------------------------------------------------------------------

end subroutine umbrella_load

!=============================================================================

subroutine umbrella_save(this, filename)

   use netcdf

   use ncsu_utils
   use ncsu_sander_proxy

   implicit none

   type(umbrella_t), intent(in) :: this
   character(*), intent(in) :: filename

   integer :: n, rc, setid, dimid, varid, nextents, ncoeffs
   integer :: periodicity(UMBRELLA_MAX_NEXTENTS) ! netCDF doesn't like logical

   ncsu_assert(this%inited)
   ncsu_assert(this%nextents.ge.1)
   ncsu_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   ncsu_assert(associated(this%coeffs))

   rc = nf90_create(filename, cmode = nf90_clobber, ncid = setid)
   call check_rc()

   nextents = this%nextents

   ncoeffs = 1
   do n = 1, nextents
      if (this%periodicity(n)) then
         periodicity(n) = 1
      else
         periodicity(n) = 0
      end if
      ncoeffs = ncoeffs*this%extents(n)
   end do

   rc = nf90_put_att(setid, nf90_global, 'nextents', nextents)
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'extents', this%extents(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'periodicity', periodicity(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'origin', this%origin(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'spacing', this%spacing(1:nextents))
   call check_rc()

   rc = nf90_def_dim(setid, 'row-major', ncoeffs, dimid)
   call check_rc()

   rc = nf90_def_var(setid, 'coeffs', coeffs_type, dimid, varid)
   call check_rc()

   rc = nf90_enddef(setid)
   call check_rc()

   rc = nf90_put_var(setid, varid, this%coeffs)
   call check_rc()

   rc = nf90_close(setid)
   call check_rc()

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine check_rc()

   use ncsu_constants, only : ERR_UNIT

   implicit none

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') NCSU_ERROR, &
         'umbrella_save(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_rc

!-----------------------------------------------------------------------------

end subroutine umbrella_save
#endif /* NCSU_NO_NETCDF */

!=============================================================================

!
! evaluates value(f) of cubic B-spline centered
! at 0.0D0 (works correctly only for -2 < x < 2)
!

pure NCSU_REAL function m4_v(x)

   implicit none

   NCSU_REAL, intent(in)  :: x

   if (x.lt.MINUS_ONE) then
      m4_v = ONE_SIXTH*(x + TWO)**3
   else if (x.lt.ZERO) then
      m4_v = MINUS_HALF*x*x*(x + TWO) + TWO_THIRD
   else if (x.lt.ONE) then
      m4_v = MINUS_HALF*x*x*(TWO - x) + TWO_THIRD
   else
      m4_v = ONE_SIXTH*(TWO - x)**3
   end if

end function m4_v

!=============================================================================

!
! evaluates value(f) and derivative(df)
! of cubic B-spline centered at 0.0D0
! (works correctly only for -2 < x < 2)
!

subroutine m4_vdv(x, f, df)

   implicit none

   NCSU_REAL, intent(in)  :: x
   NCSU_REAL, intent(out) :: f, df

   NCSU_REAL :: x2

   if (x.lt.MINUS_ONE) then
      x2 = x + TWO
      df = HALF*x2*x2
       f = ONE_THIRD*x2*df
   else if (x.lt.ZERO) then
      x2 = MINUS_HALF*x
       f = x2*x*(x + TWO) + TWO_THIRD
      df = x2*(THREE*x + FOUR)
   else if (x.lt.ONE) then
      x2 = HALF*x
       f = x2*x*(x - TWO) + TWO_THIRD
      df = x2*(THREE*x - FOUR)
   else
      x2 = TWO - x
      df = MINUS_HALF*x2*x2
       f = MINUS_ONE_THIRD*x2*df
   end if

end subroutine m4_vdv

!=============================================================================

end module ncsu_umbrella
