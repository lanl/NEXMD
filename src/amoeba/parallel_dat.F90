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
#ifdef USE_MPI_MODULE
  USE_MPI_MODULE
  implicit none
#else
  implicit none
  include 'mpif.h'
#endif
#endif

  logical, save         :: master       ! Is this task the mpi master task?

  ! Processor-specific atom count:

  integer, save         :: my_atm_cnt

#ifdef MPI

! MAJOR NOTE:  Much of the stuff here is basically unnecessary under
!              Generalized Born; due to large cutoffs and small atom counts
!              in that environment, spatial decomposition is somewhat
!              impractical.  We DO use some of the data structures; these are
!              noted by a "VGB" sidebar - meaning "valid for GB".

! Global variables necessary for MPI implementation. These are NOT broadcast.

  integer, save         :: mytaskid, numtasks, world_group, err_code_mpi   ! VGB
  integer, save         :: notdone                                         ! VGB
  logical, save         :: i_do_recip

  integer, save         :: my_send_atm_cnts_total = 0
  integer, save         :: my_send_atm_cnts_sums = 0

  integer, parameter    :: dif_tag = 10
  integer, parameter    :: dc_tag = 11
  integer, parameter    :: gifd_tag = 12
  integer, parameter    :: xyzxt_tag = 13
  integer, parameter    :: zxxyt_tag = 14

  ! Used for async mpi sends/recvs in fft and force code.

  integer, allocatable, save    :: gbl_taskmap(:)
  integer, allocatable, save    :: gbl_inv_taskmap(:)

  integer, allocatable, save    :: gbl_img_div_tbl(:)

  ! The following three arrays are for use in mpi routines defined in this
  ! module; the vectors they are used on (atm_frc, atm_crd, atm_vel) DO NOT
  ! contain data that is segregated by owning task (ie., these are not
  ! really division tables).  Instead, data must be marshalled into and out of
  ! temporary buffers according to gbl_atm_owner_map().

  integer, save         :: extra_used_atm_cnt

  integer, allocatable, save    :: gbl_atm_offsets(:)                       !VGB
  integer, allocatable, save    :: gbl_vec_offsets(:)                       !VGB
  integer, allocatable, save    :: gbl_vec_rcvcnts(:)                       !VGB

  integer, allocatable, save    :: gbl_atm_owner_map(:)                     !VGB

  integer, allocatable, save    :: gbl_my_atm_lst(:)                        !VGB

  integer, allocatable, save    :: gbl_send_atm_lst(:)
  integer, allocatable, save    :: gbl_send_atm_cnts(:)
  integer, allocatable, save    :: gbl_recv_atm_lsts(:,:)
  integer, allocatable, save    :: gbl_recv_atm_cnts(:)
  integer, allocatable, save    :: gbl_extra_used_atms(:)

  integer, save :: siz_dbl_mpi_bufs = 0                                     !VGB

  double precision, allocatable, save   :: dbl_mpi_send_buf(:)              !VGB
  double precision, allocatable, save   :: dbl_mpi_recv_buf(:)              !VGB

#endif /* MPI */

contains

! mexit ends up in this routine due to mpi implications.

!*******************************************************************************
!
! Subroutine:  mexit
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mexit(ifil, i)

  implicit none

! Formal arguments:

  integer       :: ifil
  integer       :: i

! Local variables:

#ifdef MPI
  integer       :: err_ret_code ! value returned by mpi calls 
                                ! (but these mpi calls don't return...)
#endif

  if (ifil .ne. 0) then
    close(unit = ifil)
  end if

#ifdef MPI

! i .gt. 0 implies an error condition, therefore we want to
! kill all the nodes.  This is accomplished with mpi_abort.  If
! it is not an error, exit gracefully with the mpi_finalize.
! NOTE: no mpi functions may be called after a call to mpi_finalize.

  call mpi_group_free(world_group, err_ret_code)
  if (i .eq. 0) then
    call mpi_finalize(err_ret_code)
  else
    call mpi_abort(mpi_comm_world, i, err_ret_code)
  end if

#endif

  if (i .eq. 0) then
    stop
  else
    stop 'PMEMD Terminated Abnormally!'
  end if

end subroutine mexit

end module parallel_dat_mod
