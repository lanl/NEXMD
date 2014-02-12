#include "copyright.i"

!*******************************************************************************
!
! Module: parallel_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module parallel_dat_mod

#ifdef MPI
  use gbl_datatypes_mod
  implicit none
  include 'mpif.h'
#endif

  logical, save         :: master        ! Is this task master of pmemd_comm?

  integer, save         :: numgroups = 1 ! Num of pmemd replicas. Default to 1

#ifdef MPI
! MAJOR NOTE:  Some of the stuff here is basically unnecessary under
!              Generalized Born; due to large cutoffs and small atom counts
!              in that environment, spatial decomposition is somewhat
!              impractical.  We DO use some of the data structures; these are
!              noted by a "VGB" sidebar - meaning "valid for GB".

! Global variables necessary for MPI. These are NOT broadcast. All are VGB

  integer, save         :: mytaskid      ! rank in pmemd replica
  integer, save         :: numtasks      ! size of pmemd communicator
  integer, save         :: err_code_mpi  ! error code for MPI subroutine calls
  integer, save         :: notdone       ! used to synchronize efforts
  integer, save         :: worldrank     ! rank in MPI_COMM_WORLD
  integer, save         :: worldsize     ! size of MPI_COMM_WORLD
  integer, save         :: master_rank   ! rank in the pmemd_master_comm

  integer, save         :: repnum        ! replica number (essentially just a
                                         ! bcast of master_rank to pmemd_comm)

  logical, save         :: master_master ! is this rank 0 of pmemd_masters?

! Set up the different communicators. The global one is mpi_comm_world, each
! pmemd replica will have a pmemd_comm, and the masters will be linked via
! pmemd_master_comm. All are VGB

  integer, save         :: pmemd_comm         ! replica communicator
  integer, save         :: pmemd_comm_number  ! Numeric ID for replica communicator
  integer, save         :: pmemd_group        ! intra-replica group
  integer, save         :: pmemd_master_comm  ! comm linking replica masters

  ! Atom ownership data structures:

  integer, allocatable, save    :: gbl_atm_owner_map(:)                     !VGB
  integer, save                 :: my_atm_cnt                               !VGB
  integer, allocatable, save    :: gbl_my_atm_lst(:)                        !VGB

  ! Shared MPI data buffers:

  integer, save :: siz_dbl_mpi_bufs = 0                                     !VGB

  double precision, allocatable, save   :: dbl_mpi_send_buf(:)              !VGB
  double precision, allocatable, save   :: dbl_mpi_recv_buf(:)              !VGB

! MPI communications tag parameters, centralized here to insure we don't
! double assign a value...

  integer, parameter    :: dif_tag = 10
  integer, parameter    :: dc_tag = 11
  integer, parameter    :: gifd_tag = 12
  integer, parameter    :: xyzxt_tag = 13       ! used with slabs, xy->zx
  integer, parameter    :: zxxyt_tag = 14       ! used with slabs, zx->xy
  integer, parameter    :: xyyxt_tag = 15       ! used with blks, xy->yx
  integer, parameter    :: yzzyt_tag = 16       ! used with blks, yz->zy
  integer, parameter    :: yxxyt_tag = 17       ! used with blks, yx->xy
  integer, parameter    :: zyyzt_tag = 18       ! used with blks, zy->yz
  integer, parameter    :: remd_tag = 19

  ! Logging data, for detailed atom and image usage stats.  The first variable
  ! in each pair is the value last set in a step.  The second variable is a
  ! running average (done by at each step dividing the current value of the
  ! first variable by the total step count and adding it).

  double precision, save        :: log_owned_img_cnt = 0.d0
  double precision, save        :: log_owned_img_cnt_avg = 0.d0

  double precision, save        :: log_used_img_cnt = 0.d0
  double precision, save        :: log_used_img_cnt_avg = 0.d0

  double precision, save        :: log_owned_atm_cnt = 0.d0
  double precision, save        :: log_owned_atm_cnt_avg = 0.d0

  double precision, save        :: log_used_atm_cnt = 0.d0
  double precision, save        :: log_used_atm_cnt_avg = 0.d0

  double precision, save        :: log_used_atm_source_cnt = 0.d0
  double precision, save        :: log_used_atm_source_cnt_avg = 0.d0

  double precision, save        :: log_provided_atm_cnt = 0.d0
  double precision, save        :: log_provided_atm_cnt_avg = 0.d0

  double precision, save        :: log_provided_atm_sink_cnt = 0.d0
  double precision, save        :: log_provided_atm_sink_cnt_avg = 0.d0

  double precision, save        :: log_recip_nstep = 0.d0
  
  integer, save                 :: log_listbuild_call_ctr = 0

  integer, parameter            :: log_owner_user_stats_cnt = 7

contains

!*******************************************************************************
!
! Subroutine:  parallel_dat_setup
!
! Description:  Allocation of shared parallel data structures.
!*******************************************************************************

subroutine parallel_dat_setup(atm_cnt, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in)           :: atm_cnt
  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

! Allocate storage for structures associated with keeping track of the
! division of the atom and image workload in all parallel implementations.

  allocate(gbl_atm_owner_map(atm_cnt), &
           gbl_my_atm_lst(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_atm_owner_map) + &
                        size(gbl_my_atm_lst)

  gbl_atm_owner_map(:) = 0
  gbl_my_atm_lst(:) = 0

  my_atm_cnt = 0        ! until assigned...

  return

end subroutine parallel_dat_setup

!*******************************************************************************
!
! Subroutine:  set_minimum_mpi_bufs_size
!
! Description:  This routine grows dbl_mpi_send_buf and dbl_mpi_recv_buf as
!               required.  These buffers will always be the maximum value of
!               min_buf_size seen so far.  There is no shrinking, and the 
!               buffers are always the same size.
!
!*******************************************************************************

subroutine set_minimum_mpi_bufs_size(min_buf_size, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer       :: min_buf_size
  integer       :: num_reals

! Local variables:

  integer       :: alloc_failed

  if (min_buf_size .le. siz_dbl_mpi_bufs) return

  if (allocated(dbl_mpi_send_buf)) then
    num_reals = num_reals - siz_dbl_mpi_bufs
    deallocate(dbl_mpi_send_buf)
  end if

  allocate(dbl_mpi_send_buf(min_buf_size), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + min_buf_size
  
  if (allocated(dbl_mpi_recv_buf)) then
    num_reals = num_reals - siz_dbl_mpi_bufs
    deallocate(dbl_mpi_recv_buf)
  end if

  allocate(dbl_mpi_recv_buf(min_buf_size), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + min_buf_size

  siz_dbl_mpi_bufs = min_buf_size

  return

end subroutine set_minimum_mpi_bufs_size

!*******************************************************************************
!
! Subroutine: assign_group_suffix
!
! Description: Assigns the output file suffix for this pmemd replica
!
!*******************************************************************************

subroutine assign_group_suffix(outfile_suffix)

  implicit none

! Passed arguments

  character(*), intent(out) :: outfile_suffix

! We'll make the format 3 digits. After that, we'll just use the number (which
! means they won't appear ordered according to a filesystem, necessarily)

  if (master_rank .lt. 10) then
    write(outfile_suffix, '(a,i1)') '00', master_rank

  else if (master_rank .lt. 100) then
    write(outfile_suffix, '(a,i2)') '0', master_rank

  else if (master_rank .lt. 1000) then
    write(outfile_suffix, '(i3)') master_rank

  else if (master_rank .lt. 10000) then
    write(outfile_suffix, '(i4)') master_rank

  else if (master_rank .lt. 100000) then
    write(outfile_suffix, '(i5)') master_rank

  end if

end subroutine assign_group_suffix

#endif /* MPI */

end module parallel_dat_mod
