!<compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! barf [ba:rf]  2. "He suggested using FORTRAN, and everybody barfed."
!    - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)
!

module ncsu_sander_hooks

use ncsu_constants, only : SL => STRING_LENGTH

implicit none

private

!-----------------------------------------------------------------------------
!                          T H E    H O O K S
!-----------------------------------------------------------------------------

! remd.f
#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */

! multisander.f
public :: on_multisander_exit

!-----------------------------------------------------------------------------

! sander.f
public :: on_sander_init
public :: on_sander_exit

! force.f
public :: on_force

! mdwrit.f
public :: on_mdwrit

! mdread.f
#ifdef MPI
public :: on_mdread1
#endif /* MPI */

! runmd.f
#ifdef NCSU_ENABLE_BBMD
public :: on_mdstep
#endif /* NCSU_ENABLE_BBMD */

!-----------------------------------------------------------------------------

#ifdef MPI
character(SL), private, save :: initial_mdin_name
NCSU_REAL, private, save :: initial_mdin_temp0
NCSU_REAL, private, save :: initial_mdin_sgft,initial_mdin_tempsg

private :: rem_preinit
#endif /* MPI */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   use ncsu_constants, only : ZERO

   use ncsu_pmd_hooks, only : pmd_on_delta => on_delta
   use ncsu_abmd_hooks, only : abmd_on_delta => on_delta

   implicit none

   integer, intent(in) :: o_masterrank
   logical, intent(in) :: need_U_xx

   NCSU_REAL, intent(out) :: U_mm, U_mo, U_om, U_oo

   U_mm = ZERO
   U_mo = ZERO
   U_om = ZERO
   U_oo = ZERO

   call pmd_on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)
   call abmd_on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

end subroutine on_delta

!-----------------------------------------------------------------------------

subroutine on_exchange(o_masterrank)

   use ncsu_pmd_hooks, only : pmd_on_exchange => on_exchange
   use ncsu_abmd_hooks, only : abmd_on_exchange => on_exchange

   implicit none

   integer, intent(in) :: o_masterrank

   call pmd_on_exchange(o_masterrank)
   call abmd_on_exchange(o_masterrank)

end subroutine on_exchange
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine on_multisander_exit()

   use ncsu_sander_proxy, only : proxy_finalize
   use ncsu_pmd_hooks, only : pmd_on_multisander_exit => on_multisander_exit
   use ncsu_abmd_hooks, only : abmd_on_multisander_exit => on_multisander_exit

   implicit none

#  include "ncsu-mpi.h"

   call pmd_on_multisander_exit()
   call abmd_on_multisander_exit()

   NCSU_MASTER_ONLY_BEGIN
   call proxy_finalize()
   NCSU_MASTER_ONLY_END

end subroutine on_multisander_exit

!-----------------------------------------------------------------------------

!
! 'on_sander_init()' is called at the point when
! MPI is initialized and MDIN/PRMTOP/INPCRD are loaded
!

subroutine on_sander_init(ih, amass, acrds, rem)

   use ncsu_constants
   use ncsu_sander_proxy

   use ncsu_smd_hooks, only : smd_on_sander_init => on_sander_init
   use ncsu_pmd_hooks, only : pmd_on_sander_init => on_sander_init
   use ncsu_abmd_hooks, only : abmd_on_sander_init => on_sander_init
#ifdef NCSU_ENABLE_BBMD
   use ncsu_bbmd_hooks, only : bbmd_on_sander_init => on_sander_init
#endif /* NCSU_ENABLE_BBMD */

   implicit none

   character(len = 4), intent(in) :: ih(*)

   NCSU_REAL, intent(in) :: amass(*)
   NCSU_REAL, intent(in) :: acrds(*)

   integer, intent(in) :: rem

#ifdef MPI
#  include "ncsu-mpi.h"
   logical, save :: first_call = .true.

   integer, save :: my_unit = 5, my_idx = -1
   character(STRING_LENGTH), save :: my_mdin = ''

   integer :: ios

   if (first_call) then
#endif /* MPI */
      call remember_rem(rem)

      NCSU_MASTER_ONLY_BEGIN
      call remember_atom_names(ih)
#ifdef MPI
      if (multisander_rem().lt.3 .and. multisander_rem().ne.0) then
         call rem_preinit(my_mdin, my_idx)
         if (my_mdin.ne.sander_mdin_name()) then
            my_unit = REM_MDIN_UNIT
            open (unit = my_unit, file = my_mdin, iostat = ios, &
                  form = 'FORMATTED', action = 'READ', status = 'OLD')
            if (ios > 0) then
               write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
                  NCSU_ERROR, 'failed to open ''', trim(my_mdin), &
                  ''' for reading'
               call terminate()
            end if
         else
            my_unit = 5
         end if ! my_mdin.ne.sander_mdin_name()
      else
         my_mdin = sander_mdin_name()
         my_unit = 5
      end if ! multisander_rem().ne.0
      NCSU_MASTER_ONLY_END
      first_call = .false.
   endif ! first_call
#endif /* MPI */

#ifdef MPI
   if (multisander_rem().eq.0) &
      call smd_on_sander_init(my_mdin, my_unit, amass, acrds)

   call pmd_on_sander_init(my_mdin, my_unit, amass, my_idx)
   call abmd_on_sander_init(my_mdin, my_unit, amass, my_idx)

#ifdef NCSU_ENABLE_BBMD
   if (multisander_rem().eq.0.and.multisander_numgroup().gt.1) &
      call bbmd_on_sander_init(my_mdin, my_unit, amass)
#endif /* NCSU_ENABLE_BBMD */

#else
   call smd_on_sander_init(sander_mdin_name(), 5, amass, acrds)
   call pmd_on_sander_init(sander_mdin_name(), 5, amass)
   call abmd_on_sander_init(sander_mdin_name(), 5, amass)
#endif /* MPI */

end subroutine on_sander_init

!-----------------------------------------------------------------------------

subroutine on_sander_exit()

   use ncsu_smd_hooks, only : smd_on_sander_exit => on_sander_exit
   use ncsu_pmd_hooks, only : pmd_on_sander_exit => on_sander_exit
   use ncsu_abmd_hooks, only : abmd_on_sander_exit => on_sander_exit
#ifdef NCSU_ENABLE_BBMD
   use ncsu_bbmd_hooks, only : bbmd_on_sander_exit => on_sander_exit
#endif /* NCSU_ENABLE_BBMD */

   implicit none

   call smd_on_sander_exit()
   call pmd_on_sander_exit()
   call abmd_on_sander_exit()

#ifdef NCSU_ENABLE_BBMD
   call bbmd_on_sander_exit()
#endif /* NCSU_ENABLE_BBMD */

end subroutine on_sander_exit

!-----------------------------------------------------------------------------

subroutine on_force(x, f, virial)

   use ncsu_smd_hooks, only : smd_on_force => on_force
   use ncsu_pmd_hooks, only : pmd_on_force => on_force
   use ncsu_abmd_hooks, only : abmd_on_force => on_force
#ifdef NCSU_ENABLE_BBMD
   use ncsu_bbmd_hooks, only : bbmd_on_force => on_force
#endif /* NCSU_ENABLE_BBMD */

   implicit none

   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL, intent(inout) :: f(*)
   NCSU_REAL, intent(inout) :: virial(4)

   call smd_on_force(x, f, virial)
   call pmd_on_force(x, f, virial)
   call abmd_on_force(x, f, virial)
#ifdef NCSU_ENABLE_BBMD
   call bbmd_on_force(x, f, virial)
#endif /* NCSU_ENABLE_BBMD */

end subroutine on_force

!-----------------------------------------------------------------------------

subroutine on_mdwrit()

   use ncsu_abmd_hooks, only : abmd_on_mdwrit => on_mdwrit
#ifdef NCSU_ENABLE_BBMD
   use ncsu_bbmd_hooks, only : bbmd_on_mdwrit => on_mdwrit
#endif /* NCSU_ENABLE_BBMD */

   implicit none

   call abmd_on_mdwrit()
#ifdef NCSU_ENABLE_BBMD
   call bbmd_on_mdwrit()
#endif /* NCSU_ENABLE_BBMD */

end subroutine on_mdwrit

!-----------------------------------------------------------------------------

#ifdef NCSU_ENABLE_BBMD
subroutine on_mdstep(eptot, v, ekmh)

   use ncsu_bbmd_hooks, only : bbmd_on_mdstep => on_mdstep

   implicit none

   NCSU_REAL, intent(in) :: eptot
   NCSU_REAL, intent(inout) :: v(*)
   NCSU_REAL, intent(inout) :: ekmh ! self-explaining, scientific variable

   call bbmd_on_mdstep(eptot, v, ekmh)

end subroutine on_mdstep
#endif /* NCSU_ENABLE_BBMD */

#ifdef MPI

! ><><><><><><><><><><><><><><><>< R E M D ><><><><><><><><><><><><><><><><><><

!
! remembers MDIN/temp0 [needed for proper REMD restart]
!

subroutine on_mdread1()

   use ncsu_sander_proxy

   implicit none

   initial_mdin_name  = sander_mdin_name()
   initial_mdin_temp0 = sander_temp0()
   initial_mdin_sgft = sander_sgft()
   initial_mdin_tempsg = sander_tempsg()

end subroutine on_mdread1

!
! gets proper MDIN filename
!

subroutine rem_preinit(my_mdin, my_idx)

   use ncsu_utils
   use ncsu_constants
   use ncsu_sander_proxy
   use sgld, only: trxsgld,sorttempsg,tempsglookup

   implicit none

   character(len = SL), intent(out) :: my_mdin
   integer, intent(out) :: my_idx

#  include "ncsu-mpi.h"

   NCSU_REAL, allocatable :: all_initial_temp0(:)
   NCSU_REAL, allocatable :: all_current_temp0(:)

   NCSU_REAL, allocatable :: all_initial_sgft(:)
   NCSU_REAL, allocatable :: all_current_sgft(:)
   NCSU_REAL, allocatable :: all_initial_tempsg(:)
   NCSU_REAL, allocatable :: all_current_tempsg(:)

   NCSU_REAL, parameter :: TINY = 0.00010000000000000000D0 ! NCSU_TO_REAL(0.0001)

   integer :: error, i, j, src_rank, dst_rank

   ncsu_assert(multisander_rem().ne.0)
   ncsu_assert(multisander_numgroup().gt.1)
   ncsu_assert(commmaster.ne.mpi_comm_null)

   ! get all the temperatures

   allocate(all_initial_temp0(multisander_numgroup()), &
            all_current_temp0(multisander_numgroup()), stat = error)
   if (error.ne.0) &
      NCSU_OUT_OF_MEMORY

   call mpi_allgather(initial_mdin_temp0, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_temp0, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
   ncsu_assert(error.eq.0)

   call mpi_allgather(sander_temp0(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_temp0, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
   ncsu_assert(error.eq.0)

   if(trxsgld)then
      allocate(all_initial_sgft(multisander_numgroup()), &
            all_current_sgft(multisander_numgroup()), &
            all_initial_tempsg(multisander_numgroup()), &
            all_current_tempsg(multisander_numgroup()), stat = error)
      if (error.ne.0) NCSU_OUT_OF_MEMORY

      call mpi_allgather(initial_mdin_sgft, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_sgft, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_allgather(sander_sgft(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_sgft, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      ncsu_assert(error.eq.0)
      call mpi_allgather(initial_mdin_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_allgather(sander_tempsg(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      ncsu_assert(error.eq.0)
   endif
   ! src_rank -- rank of the master that has MDIN <--> sander_temp0()
   ! dst_rank -- rank of the master that needs MDIN <--> mdin_temp0

   src_rank = mpi_proc_null
   dst_rank = mpi_proc_null

   do i = 1, multisander_numgroup()
      if(trxsgld)then
         do j = i + 1, multisander_numgroup()
            if (abs(all_initial_temp0(i) - all_initial_temp0(j)).lt.TINY .and. &
                abs(all_initial_sgft(i) - all_initial_sgft(j)).lt.TINY .and. &
                abs(all_initial_tempsg(i) - all_initial_tempsg(j)).lt.TINY ) &
               call fatal('same temp0, sgft, and tempsg in different replicas')
         end do
         if (abs(all_initial_temp0(i) - sander_temp0()).lt.TINY .and. &
             abs(all_initial_sgft(i) - sander_sgft()).lt.TINY .and. &
             abs(all_initial_tempsg(i) - sander_tempsg()).lt.TINY ) &
            src_rank = i - 1
         if (abs(all_current_temp0(i) - initial_mdin_temp0).lt.TINY .and. &
            abs(all_current_sgft(i) - initial_mdin_sgft).lt.TINY .and. &
             abs(all_current_tempsg(i) - initial_mdin_tempsg).lt.TINY ) &
            dst_rank = i - 1
      else
         do j = i + 1, multisander_numgroup()
            if (abs(all_initial_temp0(i) - all_initial_temp0(j)).lt.TINY) &
               call fatal('same temp0 in different replicas')
         end do
         if (abs(all_initial_temp0(i) - sander_temp0()).lt.TINY) &
            src_rank = i - 1
         if (abs(all_current_temp0(i) - initial_mdin_temp0).lt.TINY) &
            dst_rank = i - 1
      end if
   end do

   if (src_rank.eq.mpi_proc_null) then
      write (unit = ERR_UNIT, fmt = '(/a,a,f8.3/)') NCSU_ERROR, &
         'could not find MDIN with temp0 = ', sander_temp0()
      call terminate()
   end if

   if (dst_rank.eq.mpi_proc_null) then
      write (unit = ERR_UNIT, fmt = '(/a,a,f8.3/)') NCSU_ERROR, &
         'could not find replica that needs MDIN with temp0 = ', &
         initial_mdin_temp0
      call terminate()
   end if

   call mpi_sendrecv(initial_mdin_name, SL, MPI_CHARACTER, dst_rank, 3, &
                     my_mdin, SL, MPI_CHARACTER, src_rank, 3, &
                     commmaster, MPI_STATUS_IGNORE, error)
   ncsu_assert(error.eq.0)

   if(trxsgld)then
      ! Sort temperatures
      call sorttempsg(multisander_numgroup(),all_current_temp0,all_current_tempsg,all_current_sgft)
      ! Determine this replca's ID
      my_idx=tempsglookup(multisander_numgroup(),sander_temp0(),sander_tempsg(),sander_sgft(), &
                   all_current_temp0,all_current_tempsg,all_current_sgft)
      deallocate(all_current_sgft, all_initial_sgft, &
                     all_current_tempsg, all_initial_tempsg)
   else
   ! (bubble) sort the temperatures
   do i = 1, multisander_numgroup()
      do j = i + 1, multisander_numgroup()
         if (all_current_temp0(j).lt.all_current_temp0(i)) &
            call swap(all_current_temp0(i), all_current_temp0(j))
      end do
   end do

   my_idx = 0 ! my index in sorted T-table
   do i = 1, multisander_numgroup()
      if (abs(sander_temp0() - all_current_temp0(i)).lt.TINY) then
         my_idx = i
         exit
      end if
   end do
   endif
   
   ncsu_assert(my_idx.gt.0)

   deallocate(all_current_temp0, all_initial_temp0)

end subroutine rem_preinit

! ><><><><><><><><><><><><><><><><><><><><+><><><><><><><><><><><><><><><><><><

#endif /* MPI */

end module ncsu_sander_hooks
