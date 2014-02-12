!<compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_smd_hooks

use ncsu_constants, only : SL => STRING_LENGTH, SMD_OUTPUT_UNIT
use ncsu_colvar_type, only : colvar_base_t => colvar_t

implicit none

private

public :: on_sander_init
public :: on_sander_exit

public :: on_force

!- - - - - - - - - - - - - - - - P R I V A T E - - - - - - - - - - - - - - - -

integer, private, parameter :: smd_UNIT = SMD_OUTPUT_UNIT

character(*), private, parameter :: SECTION = 'ncsu_smd'
character(*), private, parameter :: DEFAULT_OUTPUT_FILE = 'ncsu-smd.txt'

integer, private, parameter :: DEFAULT_OUTPUT_FREQ = 50

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type, private :: colvar_t

   type(colvar_base_t) :: parent

   NCSU_REAL, pointer :: cvar_path(:) => null() ! master only
   NCSU_REAL, pointer :: harm_path(:) => null() ! master only

   NCSU_REAL :: inst, curr, harm

   integer :: harm_mode
   integer :: path_mode

end type colvar_t

integer, private, save :: ncolvars = 0 ! .gt.0 means "active"
type(colvar_t), private, allocatable, save :: cv(:)

NCSU_REAL, private, allocatable, save :: fcv_curr(:)

NCSU_REAL, private, save :: work

integer, private, save :: output_freq
integer, private, save :: mdstep ! = runmd.f::nstep + 1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, private, parameter :: MODE_LINES = 4123
integer, private, parameter :: MODE_SPLINE = 3112

private :: mode_eval
private :: mode_write
private :: mode_from_string

private :: lines, spline

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_sander_init(mdin_name, mdin_unit, amass, acrds)

   use ncsu_utils
   use ncsu_value
   use ncsu_cftree
   use ncsu_parser
   use ncsu_colvar
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   character(SL), intent(in) :: mdin_name
   integer, intent(in) :: mdin_unit

   NCSU_REAL, intent(in) :: amass(*)
   NCSU_REAL, intent(in) :: acrds(*)

   type(node_t), pointer :: root
   type(child_t), pointer :: child
   type(value_node_t), pointer :: alist, aiter

   logical :: found
   integer :: i, n, error
   character(SL) :: output_file, modestr

#  include "ncsu-mpi.h"

#  ifdef MPI
   ncsu_assert(multisander_rem().eq.0)
#  endif /* MPI */

   ncsu_assert(ncolvars.eq.0)

   NCSU_MASTER_ONLY_BEGIN

   root => parse_cf(mdin_name, SECTION, mdin_unit)
   if (.not.associated(root)) &
      goto 1

   if (sander_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   ncolvars = 0

   child => node_children(root)
   do while (associated(child))
      if (node_title(child%node) == 'variable') &
         ncolvars = ncolvars + 1
      child => child%next
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the '''//SECTION//''' section')

   allocate(cv(ncolvars), fcv_curr(ncolvars), stat = error)
   if (error /= 0) &
      NCSU_OUT_OF_MEMORY

   n = 1
   child => node_children(root)

   do while (associated(child))
      if (node_title(child%node) == 'variable') then
         ncsu_assert(n.le.ncolvars)
         call colvar_mdread(cv(n)%parent, child%node, n) 
         n = n + 1
      end if ! title == 'variable'
      child => child%next
   end do

   ! output
   found = node_lookup_string(root, 'output_file', output_file)
   if (.not.found) &
      output_file = DEFAULT_OUTPUT_FILE

   found = node_lookup_positive_integer(root, 'output_freq', output_freq)
   if (.not.found) &
      output_freq = DEFAULT_OUTPUT_FREQ

   output_freq = min(output_freq, sander_nstlim())
   output_freq = max(1, output_freq)

1  continue

   NCSU_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   ncsu_assert(error.eq.0)
#endif /* MPI */

   if (ncolvars.eq.0) &
      return

#ifdef MPI
   if (sanderrank.ne.0) then
      allocate(cv(ncolvars), fcv_curr(ncolvars), stat = error)
      if (error.ne.0) &
         NCSU_OUT_OF_MEMORY
   end if
#endif /* MPI */

   do n = 1, ncolvars
      call colvar_bootstrap(cv(n)%parent, n, amass)
      cv(n)%inst = colvar_value(cv(n)%parent, acrds)
      cv(n)%curr = ZERO
      cv(n)%harm = ZERO
      fcv_curr(n) = ZERO
   end do

   mdstep = 0
   work = ZERO

#ifdef MPI
   if (sanderrank.ne.0) &
      return
#endif /* MPI */

   ncsu_assert(associated(root))

   n = 1
   child => node_children(root)

   do while (associated(child))
      if (node_title(child%node) == 'variable') then
         ncsu_assert(n.le.ncolvars)

         ! path
         found = node_lookup_list(child%node, 'path', alist)
         if (.not.found) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
               NCSU_ERROR, 'there is no ''path'' defined for CV #', n
            call terminate()
         end if

         i = 0
         aiter => alist
         do while (associated(aiter))
            i = i + 1
            aiter => aiter%next
         end do

         if (i.lt.2) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//',a/)') &
              NCSU_ERROR, '''path'' for CV #', n, ' contains less than 2 values'
            call terminate()
         end if

         allocate(cv(n)%cvar_path(i), stat = error)
         if (error.ne.0) &
            NCSU_OUT_OF_MEMORY

         i = 0
         aiter => alist
         do while (associated(aiter))
            i = i + 1
            if (value_is_real(aiter%value)) then
               cv(n)%cvar_path(i) = value_get_real(aiter%value)
            else if (value_is_integer(aiter%value)) then
               cv(n)%cvar_path(i) = NCSU_TO_REAL(value_get_integer(aiter%value))
            else if (value_is_string(aiter%value)) then
               cv(n)%cvar_path(i) = cv(n)%inst
            else
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(n)//',a,'//pfmt(i)//',a/)') NCSU_ERROR, &
                  'CV #', n, ' : unexpected type of path(', i, ')'
               call terminate()
            end if
            aiter => aiter%next
         end do

         ! harm
         found = node_lookup_list(child%node, 'harm', alist)
         if (.not.found) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') NCSU_ERROR, &
               'there is no ''harm'' defined for CV #', n
            call terminate()
         end if

         i = 0
         aiter => alist
         do while (associated(aiter))
            i = i + 1
            aiter => aiter%next
         end do

         if (i.lt.1) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//',a/)') &
               NCSU_ERROR, '''harm'' for CV #', n, ' contains less than 1 value'
            call terminate()
         end if

         allocate(cv(n)%harm_path(i), stat = error)
         if (error.ne.0) &
            NCSU_OUT_OF_MEMORY

         i = 0
         aiter => alist
         do while (associated(aiter))
            i = i + 1
            if (value_is_real(aiter%value)) then
               cv(n)%harm_path(i) = value_get_real(aiter%value)
            else if (value_is_integer(aiter%value)) then
               cv(n)%harm_path(i) = NCSU_TO_REAL(value_get_integer(aiter%value))
            else
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(n)//',a,'//pfmt(i)//',a/)') NCSU_ERROR, &
                  'CV #', n, ' : unexpected type of harm(', i, ')'
            end if
            if (cv(n)%harm_path(i).lt.ZERO) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,'//pfmt(i)//',a/)') NCSU_ERROR, &
                  'CV #', n, ' : harm(', i, ') is negatve'
               call terminate()
            end if
            aiter => aiter%next
         end do

         found = node_lookup_string(child%node, 'harm_mode', modestr)
         if (.not.found) then
            cv(n)%harm_mode = MODE_SPLINE
         else
            cv(n)%harm_mode = mode_from_string(modestr)
            if (cv(n)%harm_mode.lt.0) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,a,a/)') NCSU_ERROR, &
                  'CV #', n, ' : unknown harm_mode = ''', trim(modestr), ''''
               call terminate()
            end if
         end if ! .not.found

         found = node_lookup_string(child%node, 'path_mode', modestr)
         if (.not.found) then
            cv(n)%path_mode = MODE_SPLINE
         else
            cv(n)%path_mode = mode_from_string(modestr)
            if (cv(n)%path_mode.lt.0) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,a,a/)') NCSU_ERROR, &
                  'CV #', n, ' : unknown path_mode = ''', trim(modestr), ''''
               call terminate()
            end if
         end if ! .not.found

         n = n + 1
      end if ! title == 'variable'
      child => child%next
   end do

   ! done with parsing

   call node_cleanup(root)
   deallocate(root)

   open (unit = smd_UNIT, file = output_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NCSU_ERROR, 'failed to open ''', trim(output_file), ''' for writing'
      call terminate()
   end if

   write (unit = smd_UNIT, fmt = '(a,/a,/a)') '#', &
      '# MD time (ps), CV, handle_position, spring_constant, work', '#'

   call flush_UNIT(smd_UNIT)

   ! print summary & we'r done

   write (unit = OUT_UNIT, fmt = '(a,a)') NCSU_INFO, &
      ' *  *  *  *  *  *  *  *  S T E E R E D  M.D.  *  *  *  *  *  *  *  *  *'

   write (unit = OUT_UNIT, fmt = '(a,/a,a,a)') NCSU_INFO, NCSU_INFO, &
      'output_file = ', trim(output_file)
   write (unit = OUT_UNIT, &
fmt = '(a,a,'//pfmt(output_freq)//',a,'//pfmt(output_freq*sander_timestep(), 4)//',a)') &
      NCSU_INFO, &
      'output_freq = ', output_freq, ' (', &
      output_freq*sander_timestep(), ' ps)'

   write (unit = OUT_UNIT, fmt = '(a)') NCSU_INFO
   do n = 1, ncolvars
      write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(n)//',/a)') &
         NCSU_INFO, 'CV #', n, NCSU_INFO
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NCSU_INFO, ' <> path = ('
      do i = 1, size(cv(n)%cvar_path)
         write (unit = OUT_UNIT, &
            fmt = '('//pfmt(cv(n)%cvar_path(i), 4)//')', advance = 'NO') &
            cv(n)%cvar_path(i)
         if (i.eq.size(cv(n)%cvar_path)) then
            write (unit = OUT_UNIT, fmt = '(a)') ')'
         else if (mod(i + 1, 5).eq.0) then
            write (unit = OUT_UNIT, fmt = '(a,/a,a)', advance = 'NO') &
               ',', NCSU_INFO, '     '
         else
            write (unit = OUT_UNIT, fmt = '(a)', advance = 'NO') ', '
         end if
      end do
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NCSU_INFO, ' <> path_mode = '
      call mode_write(cv(n)%path_mode, OUT_UNIT)

      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NCSU_INFO, ' <> harm = ('
      do i = 1, size(cv(n)%harm_path)
         write (unit = OUT_UNIT, &
            fmt = '('//pfmt(cv(n)%harm_path(i), 4)//')', advance = 'NO') &
            cv(n)%harm_path(i)
         if (i.eq.size(cv(n)%harm_path)) then
            write (unit = OUT_UNIT, fmt = '(a)') ')'
         else if (mod(i + 1, 5).eq.0) then
            write (unit = OUT_UNIT, fmt = '(a,/a,a)', advance = 'NO') &
               ',', NCSU_INFO, '     '
         else
            write (unit = OUT_UNIT, fmt = '(a)', advance = 'NO') ', '
         end if
      end do
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NCSU_INFO, ' <> harm_mode = '
      call mode_write(cv(n)%harm_mode, OUT_UNIT)
      write (unit = OUT_UNIT, fmt = '(a)') NCSU_INFO
      call colvar_print(cv(n)%parent, OUT_UNIT)
      write (unit = OUT_UNIT, fmt = '(a)') NCSU_INFO
   end do

   write (unit = OUT_UNIT, fmt = '(a,a/)') NCSU_INFO, &
      ' *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'

end subroutine on_sander_init

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_sander_exit()

   use ncsu_colvar, only : colvar_cleanup

   implicit none

   integer :: n

#  include "ncsu-mpi.h"

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n)%parent)
         NCSU_MASTER_ONLY_BEGIN
         deallocate(cv(n)%cvar_path, cv(n)%harm_path)
         NCSU_MASTER_ONLY_END
      end do

      deallocate(cv, fcv_curr)

      NCSU_MASTER_ONLY_BEGIN
      write (unit = smd_UNIT, fmt = '(a/,a,f16.10/,a)') &
         '#', '# <> total work done: ', work, '#'
      close (smd_UNIT)
      NCSU_MASTER_ONLY_END
   end if

   ncolvars = 0

end subroutine on_sander_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_force(x, f, virial)

   NCSU_USE_AFAILED

   use ncsu_colvar
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL, intent(inout) :: f(*)
   NCSU_REAL, intent(inout) :: virial(4)

   logical :: real_mdstep
   integer :: n

   NCSU_REAL :: position, dwork, f2, dcv, dharm

#  ifdef MPI
#     include "ncsu-mpi.h"
   integer :: error
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

   ncsu_assert(multisander_rem().eq.0)
   real_mdstep = (sander_init().eq.4)

   position = NCSU_TO_REAL(mdstep)/NCSU_TO_REAL(sander_nstlim())

   do n = 1, ncolvars
      cv(n)%inst = colvar_value(cv(n)%parent, x)
   end do

   NCSU_MASTER_ONLY_BEGIN
   dwork = ZERO
   do n = 1, ncolvars
      dharm = cv(n)%harm
      cv(n)%harm = mode_eval(cv(n)%harm_mode, cv(n)%harm_path, position)
      dharm = cv(n)%harm - dharm
      dcv = cv(n)%curr
      cv(n)%curr = mode_eval(cv(n)%path_mode, cv(n)%cvar_path, position)
      dcv = cv(n)%curr - dcv
      f2 = fcv_curr(n)
      fcv_curr(n) = cv(n)%harm*(cv(n)%curr - cv(n)%inst)
      f2 = (f2 + fcv_curr(n))/2
      dwork = dwork + f2*dcv + dharm*(cv(n)%curr - cv(n)%inst)**2/2
   end do
   NCSU_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(fcv_curr, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   ncsu_assert(error.eq.0)
#  endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(cv(n)%parent, x, fcv_curr(n), f)
   end do

   NCSU_MASTER_ONLY_BEGIN
   if (real_mdstep) then

      if (mdstep.gt.0) &
         work = work + dwork

      if (mod(mdstep, output_freq).eq.0) then
         write (unit = smd_UNIT, fmt = '(f12.4,1x)', advance = 'NO') &
            sander_mdtime()
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%inst
         end do
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%curr
         end do
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%harm
         end do
         write (unit = smd_UNIT, fmt = '(f16.8)') work
         call flush_UNIT(smd_UNIT)
      end if

      mdstep = mdstep + 1

   end if ! real_mdstep
   NCSU_MASTER_ONLY_END

end subroutine on_force

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function mode_eval(mode, path, t) result(y)

   NCSU_USE_AFAILED

   implicit none

   NCSU_REAL :: y
   integer, intent(in) :: mode
   NCSU_REAL, intent(in) :: path(:), t

   if (mode.eq.MODE_LINES) then
      y = lines(path, t)
   else if (mode.eq.MODE_SPLINE) then
      y = spline(path, t)
   else
      ncsu_assert_not_reached()
      y = 77.0D0
   end if

end function mode_eval

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function lines(path, t) result(y)

   NCSU_USE_AFAILED

   use ncsu_constants, only : ZERO, ONE

   implicit none

   NCSU_REAL :: y
   NCSU_REAL, intent(in) :: path(:), t

   integer :: npoints, n

   NCSU_REAL :: s

   npoints = size(path)

   if (npoints.eq.1) then
      y = path(1)
      return
   end if

   if (t.le.ZERO) then
      y = path(1)
   else if (t.ge.ONE) then
      y = path(npoints)
   else
      n = 1 + int(floor((npoints - 1)*t))
      ncsu_assert(n.ge.1)
      ncsu_assert(n.lt.npoints)

      ! 0 < s < 1 between path(n) and path(n + 1)
      s = (npoints - 1)*t - n + 1
      y = (ONE - s)*path(n) + s*path(n + 1)
   end if

end function lines

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function spline(path, t) result(y)

   NCSU_USE_AFAILED

   use ncsu_constants, only : ZERO, ONE, TWO

   implicit none

   NCSU_REAL :: y
   NCSU_REAL, intent(in) :: path(:), t

   integer :: npoints, n

   NCSU_REAL :: m1, m2
   NCSU_REAL :: s, s2, s3

   npoints = size(path)

   if (npoints.eq.1) then
      y = path(1)
      return
   end if

   if (t.le.ZERO) then
      y = path(1)
   else if (t.ge.ONE) then
      y = path(npoints)
   else
      n = 1 + int(floor((npoints - 1)*t))
      ncsu_assert(n.ge.1)
      ncsu_assert(n.lt.npoints)

      if (npoints.eq.2) then
         m1 = ZERO
         m2 = ZERO
      else if (n.eq.1) then
         m1 = ZERO
         m2 = (path(3) - path(1))/2
      else if ((n + 1).eq.npoints) then
         m1 = (path(n + 1) - path(n - 1))/2
         m2 = ZERO
      else
         m1 = (path(n + 1) - path(n - 1))/2
         m2 = (path(n + 2) - path(n))/2
      end if

      ! compute the value
      s = (npoints - 1)*t - n + 1

      s2 = s*s
      s3 = 2*s - NCSU_TO_REAL(3)

      y = (s2*s3 + ONE)*path(n) &
        + s*(s*(s - TWO) + ONE)*m1 &
        - s2*s3*path(n + 1) &
        + s2*(s - ONE)*m2
   end if

end function spline

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine mode_write(m, lun)

   implicit none

   integer, intent(in) :: m
   integer, intent(in) :: lun

   if (m.eq.MODE_LINES) then
      write (unit = lun, fmt = '(a)') 'LINES'
   else if (m.eq.MODE_SPLINE) then
      write (unit = lun, fmt = '(a)') 'SPLINE'
   else
      write (unit = lun, fmt = '(a)') 'UNKNOWN'
   end if

end subroutine mode_write

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function mode_from_string(str) result(mode)

    implicit none

    integer :: mode
    character(*), intent(in) :: str

    if (str.eq.'LINES') then
        mode = MODE_LINES
    else if (str.eq.'SPLINE') then
        mode = MODE_SPLINE
    else
        mode = -1
    end if

end function mode_from_string

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module ncsu_smd_hooks
