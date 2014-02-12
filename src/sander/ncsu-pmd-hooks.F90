!<compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_pmd_hooks

use ncsu_constants, only : SL => STRING_LENGTH, PMD_OUTPUT_UNIT
use ncsu_colvar_type, only : colvar_t

#ifdef MPI
use ncsu_constants, only : PMD_REMLOG_UNIT
#endif /* MPI */

implicit none

private

#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */

public :: on_multisander_exit

public :: on_sander_init
public :: on_sander_exit

public :: on_force

!- - - - - - - - - - - - - - - - P R I V A T E - - - - - - - - - - - - - - - -

integer, private, parameter :: pmd_UNIT = PMD_OUTPUT_UNIT

character(*), private, parameter :: SECTION = 'ncsu_pmd'
character(*), private, parameter :: DEFAULT_OUTPUT_FILE = 'ncsu-pmd'

integer, private, parameter :: DEFAULT_OUTPUT_FREQ = 50

#ifdef MPI
character(SL), private, parameter :: REMLOG_FILE = 'ncsu-pmd.log'
integer, private, parameter :: REMLOG_UNIT = PMD_REMLOG_UNIT
#endif /* MPI */

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, private, save :: ncolvars = 0 ! .gt.0 means "active"
type(colvar_t), private, allocatable, save :: cv(:)

NCSU_REAL, private, pointer, save :: anchor(:) => null() ! master
NCSU_REAL, private, pointer, save :: a_position(:) => null() ! master
NCSU_REAL, private, pointer, save :: a_strength(:) => null() ! master

NCSU_REAL, private, allocatable, save :: cv_inst(:)
NCSU_REAL, private, allocatable, save :: f_cv(:)

character(SL), private, save :: output_file
character(SL), private, save :: output_fmt

integer, private, save :: output_freq
integer, private, save :: mdstep ! = runmd.f::nstep + 1 (not zeroed on exchange)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI
NCSU_REAL, private, pointer, save :: o_anchor(:) => null() ! master
#endif /* MPI */

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI
subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   use ncsu_colvar, only : colvar_difference
   use ncsu_constants, only : ZERO

#  ifndef NCSU_DISABLE_ASSERT
   use ncsu_utils
   use ncsu_sander_proxy
#  endif /* NCSU_DISABLE_ASSERT */

   implicit none

   integer, intent(in) :: o_masterrank
   logical, intent(in) :: need_U_xx

   NCSU_REAL, intent(inout) :: U_mm, U_mo, U_om, U_oo

#  include "ncsu-mpi.h"

   NCSU_REAL :: U_o(2), U_m(2)

   integer :: n, error

   if (ncolvars.eq.0) &
      return

   ncsu_assert(multisander_rem().ne.0)
   ncsu_assert(sanderrank.eq.0) ! master
   ncsu_assert(commmaster.ne.mpi_comm_null)

   ! exchange cv_inst(:) with the partner [store partner values in f_cv]
   call mpi_sendrecv &
      (cv_inst, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
          f_cv, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)

   ! evaluate 'my' values
   U_m(1) = ZERO ! U_mm = U_m(x_m)
   U_m(2) = ZERO ! U_mo = U_m(x_o)

   do n = 1, ncolvars
      U_m(1) = U_m(1) &
       + a_strength(n)*colvar_difference(cv(n), cv_inst(n), a_position(n))**2/2
      U_m(2) = U_m(2) &
       + a_strength(n)*colvar_difference(cv(n), f_cv(n), a_position(n))**2/2
   end do

   ! get partner's U_m? (i.e., U_o? in this replica)
   call mpi_sendrecv &
      (U_m, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       U_o, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)

   if (need_U_xx) then
      U_mm = U_mm + U_m(1)
      U_mo = U_mo + U_m(2)
      U_om = U_om + U_o(2)
      U_oo = U_oo + U_o(1)
   end if

end subroutine on_delta

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_exchange(o_masterrank)

#  ifndef NCSU_DISABLE_ASSERT
   use ncsu_utils
   use ncsu_sander_proxy
#  endif /* NCSU_DISABLE_ASSERT */

   implicit none

   integer, intent(in) :: o_masterrank

#  include "ncsu-mpi.h"

   character(SL) :: o_output_file
   integer :: o_output_freq, error

   if (ncolvars.eq.0) &
      return

   ncsu_assert(multisander_rem().ne.0)
   ncsu_assert(sanderrank.eq.0) ! master
   ncsu_assert(commmaster.ne.mpi_comm_null) ! master

   ! slow & naive

   call mpi_sendrecv(output_file, SL, MPI_CHARACTER, o_masterrank, 5, &
                   o_output_file, SL, MPI_CHARACTER, o_masterrank, 5, &
                     commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)
   output_file = o_output_file

   call mpi_sendrecv(output_freq, 1, MPI_INTEGER, o_masterrank, 6, &
                   o_output_freq, 1, MPI_INTEGER, o_masterrank, 6, &
                     commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)
   output_freq = o_output_freq

   ncsu_assert(associated(anchor))
   ncsu_assert(associated(o_anchor))

   call mpi_sendrecv &
      (anchor, 2*ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, 7, &
     o_anchor, 2*ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, 7, &
       commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)

   anchor(1:2*ncolvars) = o_anchor(1:2*ncolvars)

end subroutine on_exchange
#endif /* MPI */

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_multisander_exit()

   use ncsu_colvar, only : colvar_cleanup
#  ifdef MPI
   use ncsu_sander_proxy, only : multisander_rem
#  endif /* MPI */

   implicit none

   integer :: n

#  include "ncsu-mpi.h"

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n))
      end do
      deallocate(cv, f_cv, cv_inst)
      NCSU_MASTER_ONLY_BEGIN
      nullify(a_position, a_strength)
      deallocate(anchor)
#     ifdef MPI
      if (multisander_rem().ne.0) &
         deallocate(o_anchor)
#     endif /* MPI */
      NCSU_MASTER_ONLY_END
   end if

   mdstep = 0
   ncolvars = 0

end subroutine on_multisander_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI
subroutine on_sander_init(mdin_name, mdin_unit, amass, rem_idx)
#else
subroutine on_sander_init(mdin_name, mdin_unit, amass)
#endif /* MPI */

   use ncsu_utils
   use ncsu_cftree
   use ncsu_parser
   use ncsu_colvar
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   character(SL), intent(in) :: mdin_name
   integer, intent(in) :: mdin_unit

   NCSU_REAL, intent(in) :: amass(*)

#  ifdef MPI
   integer, intent(in) :: rem_idx
#  endif /* MPI */

   type(node_t), pointer :: root
   type(child_t), pointer :: child

   logical :: found
   integer :: n, error

#  ifdef MPI
   logical, save :: first_time = .true.
   integer :: LOG_UNIT
#  endif /* MPI */

#  include "ncsu-mpi.h"

#  ifdef MPI
   ncsu_assert(first_time.or.multisander_rem().ne.0)

   if (.not.first_time) then
      ncsu_assert(multisander_rem().ne.0)
      NCSU_MASTER_ONLY_BEGIN
      if (ncolvars.gt.0) then
         ! re-open pmd_UNIT after exchange (closed in on_sander_exit())
         open (unit = pmd_UNIT, file = output_file, &
               iostat = error, form = 'FORMATTED', action = 'WRITE', &
               position = 'APPEND', status = 'OLD')

         if (error.ne.0) then
            write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NCSU_ERROR, &
               'could not open ''', trim(output_file), ''' file for writing'
            call terminate()
         end if
      end if ! ncolvars.gt.0
      NCSU_MASTER_ONLY_END
      return
   end if

   first_time = .false.

   NCSU_MASTER_ONLY_BEGIN
#  endif /* MPI */

   ncsu_assert(ncolvars.eq.0)
   root => parse_cf(mdin_name, SECTION, mdin_unit)

   if (.not.associated(root)) &
      goto 1

   ncolvars = 0

   child => node_children(root)
   do while (associated(child))
      if (node_title(child%node) == 'variable') &
         ncolvars = ncolvars + 1
      child => child%next
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the '''//SECTION//''' section')

   allocate(cv(ncolvars), anchor(2*ncolvars), &
      f_cv(ncolvars), cv_inst(ncolvars), stat = error)
   if (error /= 0) &
      NCSU_OUT_OF_MEMORY

   a_position => anchor(0*ncolvars + 1:1*ncolvars)
   a_strength => anchor(1*ncolvars + 1:2*ncolvars)

   n = 1
   child => node_children(root)

   do while (associated(child))
      if (node_title(child%node) == 'variable') then
         ncsu_assert(n.le.ncolvars)

         call colvar_mdread(cv(n), child%node, n)

         found = node_lookup_real(child%node, 'anchor_strength', a_strength(n))
         if (.not.found) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
               NCSU_ERROR, 'could not find ''anchor_strength'' for CV #', n
            call terminate()
         end if

         found = node_lookup_real(child%node, 'anchor_position', a_position(n))
         if (.not.found) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
               NCSU_ERROR, 'could not find ''anchor_position'' for CV #', n
            call terminate()
         end if

         n = n + 1
      end if ! title == 'variable'
      child => child%next
   end do

   ! output
   found = node_lookup_string(root, 'output_file', output_file)
   if (.not.found) then
#  ifdef MPI
      if (multisander_rem().ne.0) then
         write (unit = output_file, fmt = '(a,a,i3.3,a)') &
            DEFAULT_OUTPUT_FILE, '.', rem_idx, '.txt'
      else
#  endif /* MPI */
         output_file = DEFAULT_OUTPUT_FILE//'.txt'
#  ifdef MPI
      end if ! multisander_rem().ne.0
#  endif /* MPI */
   end if ! .not.found

   found = node_lookup_positive_integer(root, 'output_freq', output_freq)
   if (.not.found) &
      output_freq = DEFAULT_OUTPUT_FREQ

   output_freq = min(output_freq, sander_nstlim())
   output_freq = max(1, output_freq)

   call node_cleanup(root)
   deallocate(root)

1  continue

   NCSU_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   ncsu_assert(error.eq.0)
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

#  ifdef MPI
   if (sanderrank.ne.0) then
      allocate(cv(ncolvars), f_cv(ncolvars), cv_inst(ncolvars), stat = error)
      if (error.ne.0) &
         NCSU_OUT_OF_MEMORY
   end if
#  endif /* MPI */

   do n = 1, ncolvars
      call colvar_bootstrap(cv(n), n, amass)
   end do

   mdstep = 0

   NCSU_MASTER_ONLY_BEGIN

   open (unit = pmd_UNIT, file = output_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NCSU_ERROR, 'failed to open ''', trim(output_file), ''' for writing'
      call terminate()
   end if

   write (unit = pmd_UNIT, fmt = '(a,66(''=''))') '# = NCSU%PMD '
   do n = 1, ncolvars
      write (unit = pmd_UNIT, &
         fmt = '(a,'//pfmt(n)//',a,'//pfmt(a_position(n), 6)//',a,'//pfmt(a_strength(n), 6)//',a)') &
         '#   << anchor(', n, ') : position = ', a_position(n), &
         ', strength = ', a_strength(n), ' >>'
   end do
   write (unit = pmd_UNIT, &
      fmt = '(a,77(''-''),/a,'//pfmt(ncolvars)//',a,/a,77(''=''))') '# ', &
      '# MD time (ps), CV(1:', ncolvars, ')', '# '

   call flush_UNIT(pmd_UNIT)

   write (unit = output_fmt, fmt = '(a,'//pfmt(ncolvars)//',a)') &
      '(f12.4,', ncolvars, '(1x,f16.8))'

   ! print summary & we'r done

#  ifdef MPI
   if (multisander_rem().eq.0) then
      LOG_UNIT = OUT_UNIT ! write to MDOUT
   else
      allocate(o_anchor(2*ncolvars), stat = error)
      if (error.ne.0) &
         NCSU_OUT_OF_MEMORY

      LOG_UNIT = REMLOG_UNIT

      if (masterrank.eq.0) then
         open (unit = LOG_UNIT, file = REMLOG_FILE, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

         if (error.ne.0) then
            write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NCSU_ERROR, &
               'failed to open ''', trim(REMLOG_FILE), ''' for writing'
            call terminate()
         end if
         write (unit = LOG_UNIT, fmt = '(80(''*''))')
         close (unit = LOG_UNIT)
      end if ! masterrank.eq.0

      call cpus_enter(commmaster, 123456)

      open (unit = LOG_UNIT, file = REMLOG_FILE, &
            iostat = error, form = 'FORMATTED', action = 'WRITE', &
            position = 'APPEND', status = 'OLD')

      if (error.ne.0) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NCSU_ERROR, &
            'could not open ''', trim(REMLOG_FILE), ''' file for writing'
         call terminate()
      end if

      write (unit = LOG_UNIT, fmt = '(/23x,a,'//pfmt(masterrank)//',a,f7.3,a/)') &
        'REPLICA #', masterrank, ' (temp0 = ', sander_temp0(), ')'
      write (unit = LOG_UNIT, fmt = '(1x,a,a,a/)') &
         'MDIN = ''', trim(mdin_name), ''''

   end if ! multisander_rem().ne.0
#  else
#     define LOG_UNIT OUT_UNIT
#  endif

   write (unit = LOG_UNIT, fmt = '(a,a)') NCSU_INFO, &
      '~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ P I N N E D  M.D. ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~'

   write (unit = LOG_UNIT, fmt = '(a,/a,a,a)') NCSU_INFO, NCSU_INFO, &
      'output_file = ', trim(output_file)
   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(output_freq)//',a,'//pfmt(output_freq*sander_timestep(), 4)//',a)') &
        NCSU_INFO, 'output_freq = ', output_freq, ' (', &
        output_freq*sander_timestep(), ' ps)'

   write (unit = LOG_UNIT, fmt = '(a)') NCSU_INFO
   do n = 1, ncolvars
      write (unit = LOG_UNIT, &
         fmt = '(a,a,'//pfmt(n)//',a,'//pfmt(a_position(n), 6)//',a,'//pfmt(a_strength(n), 6)//',a)') &
         NCSU_INFO, 'CV #', n, ' << anchor : position = ', a_position(n), &
         ', strength = ', a_strength(n), ' >>'
      call colvar_print(cv(n), LOG_UNIT)
      write (unit = LOG_UNIT, fmt = '(a)') NCSU_INFO
   end do

   write (unit = LOG_UNIT, fmt = '(a,a/)') NCSU_INFO, &
      '~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~'

#  ifdef MPI
   if (multisander_rem().ne.0) then
      write (unit = LOG_UNIT, fmt = '(/80(''*''))')
      close (unit = LOG_UNIT)
      call cpus_leave(commmaster, 123456)
   end if
#  endif /* MPI */

   NCSU_MASTER_ONLY_END

end subroutine on_sander_init

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_sander_exit()

   use ncsu_utils, only : close_UNIT

   implicit none

#  include "ncsu-mpi.h"

   NCSU_MASTER_ONLY_BEGIN
   call close_UNIT(pmd_UNIT)
   NCSU_MASTER_ONLY_END

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

#  ifdef MPI
#     include "ncsu-mpi.h"
   integer :: error
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

   real_mdstep = (.not.multisander_initremd().and.sander_init().eq.4)

   do n = 1, ncolvars
      cv_inst(n) = colvar_value(cv(n), x)
   end do

   NCSU_MASTER_ONLY_BEGIN
   do n = 1, ncolvars
      f_cv(n) = &
         - a_strength(n)*colvar_difference(cv(n), cv_inst(n), a_position(n))
   end do
   NCSU_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(f_cv, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   ncsu_assert(error.eq.0)
#  endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(cv(n), x, f_cv(n), f)
   end do

   NCSU_MASTER_ONLY_BEGIN
   if (real_mdstep) then

      if (mod(mdstep, output_freq).eq.0) then
         write (unit = pmd_UNIT, fmt = output_fmt) &
            sander_mdtime(), cv_inst(1:ncolvars)
         call flush_UNIT(pmd_UNIT)
      end if

      mdstep = mdstep + 1

   end if ! real_mdstep
   NCSU_MASTER_ONLY_END

end subroutine on_force

end module ncsu_pmd_hooks
