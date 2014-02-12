#include "copyright.i"

!*******************************************************************************
!
! Module:  loadbal_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module loadbal_mod

  implicit none

#ifdef MPI

  ! Load balancing timers and other data structures:

  integer, private, save        :: start_sec, start_usec
  integer, save                 :: elapsed_100usec_dirfrc = 0
  integer, save                 :: elapsed_100usec_recipfrc = 0
  integer, save                 :: elapsed_100usec_other = 0  ! bnd,ang,dihed,
                                                              ! other atm work.
  integer, save                 :: elapsed_100usec_listbld = 0 

  ! Some counters start at -1 to inhibit 1st pass redistribution.

  integer, private, save        :: img_redist_ctr = -1  
  integer, private, save        :: img_redist_trigger = 1 ! 1,2,4,8,16,32,64
  integer, parameter            :: img_redist_trigger_max = 64
  integer, parameter            :: img_redist_trigger_divisor = 4
  double precision, parameter   :: img_redist_time_tol = 0.01d0

  logical, save                 :: fft_slab_redist_enabled = .false.
  logical, save                 :: fft_slab_redist_needed = .false.
  integer, private, save        :: fft_slab_redist_ctr = -1
  integer, private, save        :: retry_fft_slab_redist = 0
  integer, private, save        :: retry_fft_slab_redist_step_ctr = 0
  integer, private, save        :: retry_fft_slab_redist_interval = 2000
  integer, parameter            :: retry_fft_slab_redist_incr = 2000

  logical, save                 :: atm_redist_needed = .false.
  logical, save                 :: force_atm_redist = .false.
  double precision, parameter   :: atm_redist_tol = 1.05d0
  integer, save, private        :: last_send_atm_cnts_total = 0
  logical, save, private        :: recalc_last_send_atm_cnts_total = .false.

  integer, save                 :: loadbal_step_ctr = 0
  integer, save                 :: last_loadbal_step_ctr = 0

  integer, allocatable, save    :: gbl_loadbal_node_dat(:,:)

  integer, save                 :: recip_numtasks = 0

contains

!*******************************************************************************
!
! Subroutine:  start_loadbal_timer
!
! Description: Used by cit code.
!              
!*******************************************************************************

subroutine start_loadbal_timer

  implicit none

  call get_wall_time(start_sec, start_usec)

  return

end subroutine start_loadbal_timer

!*******************************************************************************
!
! Subroutine:  update_loadbal_timer
!
! Description: Used by cit code. Note - Only can run one loadbal timer at a
!              time (ie., the adjustable and nonadjustable timer start/stop
!              cycles must not overlap (but they would never need to)).
!              
!*******************************************************************************

subroutine update_loadbal_timer(elapsed_100usec)

  implicit none

! Formal arguments:

  integer               :: elapsed_100usec

! Local variables:

  integer               :: stop_sec, stop_usec
  integer               :: next_start_sec, next_start_usec
  integer               :: this_exec_time

  call get_wall_time(stop_sec, stop_usec)

  next_start_sec = stop_sec
  next_start_usec = stop_usec
    
  ! For next round:

  if (stop_usec .lt. start_usec) then
    stop_usec = stop_usec + 1000000
    stop_sec = stop_sec - 1
  end if

  this_exec_time = (stop_sec - start_sec) * 10000 + &
                   (stop_usec - start_usec) / 100

  elapsed_100usec = elapsed_100usec + this_exec_time

  start_sec = next_start_sec
  start_usec = next_start_usec

  return

end subroutine update_loadbal_timer

!*******************************************************************************
!
! Subroutine:  do_load_balancing
!
! Description: The main entry into load balancing code; should be called at
!              the beginning of each force evaluation that has a new pairlist.
!              
!*******************************************************************************

subroutine do_load_balancing(new_list, atm_cnt)

  use mdin_ctrl_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  logical               :: new_list
  integer               :: atm_cnt

  loadbal_step_ctr = loadbal_step_ctr + 1

  if (.not. new_list) return

  if (fft_slab_redist_enabled) then
    if (.not. fft_slab_redist_needed) then
      if (retry_fft_slab_redist .ne. 0) then
        if (retry_fft_slab_redist_step_ctr .le. loadbal_step_ctr) then
          fft_slab_redist_needed = .true.
          fft_slab_redist_ctr = -1
          retry_fft_slab_redist_step_ctr = 0
          atm_redist_needed = .false.   ! Cancel outstanding atm redist, as
                                        ! with fft slab dist coming, it is
                                        ! useless.
        end if
      end if 
    end if
  end if

  if (fft_slab_redist_needed) then

    fft_slab_redist_ctr = fft_slab_redist_ctr + 1
    if (fft_slab_redist_ctr .gt. 0) call do_fft_slab_redistribution(atm_cnt)
    img_redist_ctr = 0
    img_redist_trigger = 1
    call update_pme_time(fft_slab_reassign_timer)

  else if (atm_redist_needed) then

    call do_atm_redistribution(atm_cnt, .true.)
    atm_redist_needed = .false.
    recalc_last_send_atm_cnts_total = .true.
    img_redist_ctr = 0
    img_redist_trigger = 1
    call update_pme_time(atm_reassign_timer)

  else

    img_redist_ctr = img_redist_ctr + 1
    if (img_redist_ctr .ge. 1 .and. img_redist_ctr .ge. &
        img_redist_trigger * nrespa / img_redist_trigger_divisor) then
      call do_img_redistribution
      img_redist_ctr = 0
      img_redist_trigger = min(img_redist_trigger * 2, img_redist_trigger_max)
      call update_pme_time(img_reassign_timer)
    end if

  end if

  return

end subroutine do_load_balancing

!*******************************************************************************
!
! Subroutine:  check_new_list_limit
!
! Description: Used by cit code. Note - Only can run one loadbal timer at a
!              time (ie., the adjustable and nonadjustable timer start/stop
!              cycles must not overlap (but they would never need to)).
!              
!*******************************************************************************

subroutine check_new_list_limit(new_list)

  implicit none

! Formal arguments:

  logical               :: new_list

! Local variables:
                                                                                  integer, save         :: last_new_list_cnt = 0
  integer, save         :: last_new_list_limit = 16
  integer, parameter    :: max_last_new_list_limit = 32

  if (new_list) then
    last_new_list_cnt = 0
  else
    last_new_list_cnt = last_new_list_cnt + 1
    if (last_new_list_cnt .ge. last_new_list_limit) then
      last_new_list_cnt = 0
      last_new_list_limit = last_new_list_limit + 1
      last_new_list_limit = min(last_new_list_limit, max_last_new_list_limit)
      new_list = .true.
    end if
  end if

  return

end subroutine check_new_list_limit

!*******************************************************************************
!
! Subroutine:  do_fft_slab_redistribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine do_fft_slab_redistribution(atm_cnt)

  use pme_fft_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt

! Local variables:

  double precision              :: avg_total_time
  double precision              :: recip_nodes_time
  double precision              :: avg_recip_nodes_time
  double precision, save        :: last_avg_recip_nodes_time = 0.d0
  double precision              :: recip_slab_time
  double precision              :: recip_workload
  double precision              :: total_time
  logical                       :: not_distributed
  integer                       :: i
  integer                       :: node
  integer                       :: num_ints, num_reals  ! values are not kept.
  integer                       :: use_recip_workload

  ! We place the following data in static storage to keep the Myrinet MPI
  ! implementation from going ape over mpi buffers on the stack (causes
  ! problems between mpi and ifc, resulting in huge stalls in processing).

  integer, save                 :: my_node_dat(5)
  integer, save                 :: fft_slab_dat(6)

  num_ints = 0
  num_reals = 0
  use_recip_workload = 0

  ! Gather times from the other processes to the master process.  The
  ! reference to gbl_loadbal_node_dat is only used in the master process.
  ! The dat collected from each node are the direct, reciprocal, other
  ! (basically associated with atom ownership) and listbuild work times.

  my_node_dat(1) = elapsed_100usec_dirfrc
  my_node_dat(2) = elapsed_100usec_recipfrc
  my_node_dat(3) = elapsed_100usec_other
  my_node_dat(4) = elapsed_100usec_listbld
  my_node_dat(5) = 0

  call mpi_gather(my_node_dat, 5, mpi_integer, gbl_loadbal_node_dat, 5, &
                  mpi_integer, 0, mpi_comm_world, err_code_mpi)

  if (master) then

    if (retry_fft_slab_redist .ne. 0) then

      retry_fft_slab_redist = 0
      recip_workload = fft_workload_estimate
      not_distributed = .true.
      avg_recip_nodes_time = 0.d0
      use_recip_workload = 1
      fft_slab_redist_ctr = 0   ! Allows triggering of next cycle which uses
                                ! a computed recip_workload.
    else

      total_time = 0.d0
      recip_nodes_time = 0.d0

      do i = 0, numtasks - 1
        if (xy_slab_cnt(i) .gt. 0) then
          recip_nodes_time = recip_nodes_time + dble(gbl_loadbal_node_dat(2, i))
        end if
        ! We intentionally omit list building time from the total time because
        ! it does not occur every cycle, is not completely variable with image
        ! count, and basically adds an undesirable element of randomness.
        total_time = total_time + &
                     dble(gbl_loadbal_node_dat(1, i) + &
                          gbl_loadbal_node_dat(2, i) + &
                          gbl_loadbal_node_dat(3, i))
      end do

      recip_slab_time = recip_nodes_time / dble(fft_z_dim)
      avg_total_time = total_time / dble(numtasks)
      avg_recip_nodes_time = recip_nodes_time / &
                             dble(recip_numtasks * &
                             (loadbal_step_ctr - last_loadbal_step_ctr + 1))

      recip_workload = recip_nodes_time / total_time

! BEGIN DBG
!     write(logfile,*)'DBG: avg_recip_nodes_time=', &
!                     avg_recip_nodes_time
!     write(logfile,*)'DBG: last_avg_recip_nodes_time=', &
!                     last_avg_recip_nodes_time
!     write(logfile,*)'DBG: recip_slab_time*max_xy_slabs=',&
!                     recip_slab_time*dble(max_xy_slab_cnt)
!     write(logfile,*)'DBG: avg_total_time=',avg_total_time
!     write(logfile,*)'DBG: elapsed steps=',&
!                     loadbal_step_ctr - last_loadbal_step_ctr + 1
! END DBG

      if (recip_slab_time * dble(max_xy_slab_cnt) .gt. avg_total_time) then

        if (fft_slab_redist_ctr .eq. 1) then
          not_distributed = .true.
          use_recip_workload = 1
        else
          if (avg_recip_nodes_time .ge. last_avg_recip_nodes_time .or. &
              max_xy_slab_cnt .eq. 1 .or. &
              max_zx_slab_cnt .eq. 1 .or. &
              recip_numtasks .eq. numtasks) then

            not_distributed = .false.
            
            ! There is either a performance degradation, or we have split the
            ! slabs up as far as we can.  In either case, we may have an
            ! interconnect performance degradation, and we set things up for
            ! a retry later.  The retry's are backed off linearly.

            retry_fft_slab_redist = 1
            retry_fft_slab_redist_step_ctr = loadbal_step_ctr + &
                                             retry_fft_slab_redist_interval
            retry_fft_slab_redist_interval = retry_fft_slab_redist_interval + &
                                             retry_fft_slab_redist_incr
          else
            not_distributed = .true.
          end if
        end if
      else
        not_distributed = .false.
      end if

    end if

    if (.not. not_distributed) then
      ! Negative value indicates balancing done; don't need to redistrib
      ! slabs, but DO adjust images and atoms using workload.
      recip_workload = -recip_workload

      ! BEGIN DBG - Checking out multiple fft slab redist cycles...
      ! Artificially trigger another cycle...
      !     retry_fft_slab_redist = 1
      !     retry_fft_slab_redist_step_ctr = loadbal_step_ctr + &
      !                                      retry_fft_slab_redist_interval
      !     retry_fft_slab_redist_interval = retry_fft_slab_redist_interval + &
      !                                      retry_fft_slab_redist_incr
      ! END DBG
    end if

    ! Integerize the reciprocal workload, and send the other fft slab
    ! distribution control variables to the slaves.

    fft_slab_dat(1) = int(recip_workload * 1000000.d0)
    fft_slab_dat(2) = retry_fft_slab_redist
    fft_slab_dat(3) = retry_fft_slab_redist_step_ctr
    fft_slab_dat(4) = retry_fft_slab_redist_interval
    fft_slab_dat(5) = fft_slab_redist_ctr
    fft_slab_dat(6) = use_recip_workload

    call mpi_bcast(fft_slab_dat, 6, mpi_integer, 0, mpi_comm_world, &
                   err_code_mpi)

    if (loadbal_verbose .gt. 1 .and. recip_workload .gt. 0.d0) then
      write(logfile, '(2x,a,f5.1,a)') 'Recip force workload estimate is ', &
                         recip_workload * 100.d0, ' percent'
    end if

    last_avg_recip_nodes_time = avg_recip_nodes_time

  else

    call mpi_bcast(fft_slab_dat, 6, mpi_integer, 0, mpi_comm_world, &
                   err_code_mpi)

    retry_fft_slab_redist = fft_slab_dat(2)
    retry_fft_slab_redist_step_ctr = fft_slab_dat(3)
    retry_fft_slab_redist_interval = fft_slab_dat(4)
    fft_slab_redist_ctr = fft_slab_dat(5)
    use_recip_workload = fft_slab_dat(6)

  end if

  ! Recip workload recalc'd in master also to insure exact same result as
  ! in slaves.

  recip_workload = dble(fft_slab_dat(1))/1000000.d0

  if (recip_workload .gt. 0.d0) then

    if (use_recip_workload .ne. 0) then
      recip_numtasks = ceiling(dble(numtasks) * recip_workload)
    else
      ! We must only execute this code if max_xy_slab_cnt .gt. 1, currently
      ! insured by .not. not_distributed --> recip_workload .lt. 0.d0
      recip_numtasks = fft_z_dim / (max_xy_slab_cnt - 1)
      if (mod(fft_z_dim, max_xy_slab_cnt - 1) .ne. 0) then
        recip_numtasks = recip_numtasks + 1
      end if
    end if

    call distribute_fft_slabs(recip_numtasks, fft_slab_redist_enabled, &
                              loadbal_step_ctr, num_ints, num_reals)
  else
    recip_workload = -recip_workload
    fft_slab_redist_needed = .false.
  end if

  if (recip_numtasks .eq. numtasks) then
    call divide_images_evenly(gbl_img_div_tbl, .false.)
  else
    call divide_images_recip_biased(recip_workload, gbl_img_div_tbl, .false.)
  end if
  call do_atm_redistribution(atm_cnt, .false.)
  force_atm_redist = .true.

  elapsed_100usec_dirfrc = 0
  elapsed_100usec_recipfrc = 0
  elapsed_100usec_other = 0
  elapsed_100usec_listbld = 0
  last_loadbal_step_ctr = loadbal_step_ctr

  return

end subroutine do_fft_slab_redistribution

!*******************************************************************************
!
! Subroutine:  do_img_redistribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine do_img_redistribution

  use img_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Local variables:

  integer               :: i
  integer               :: node
  integer               :: send_atm_cnts_total
  double precision      :: total_time
  double precision      :: recip_time
  logical, save         :: initial_img_div_logged = .false.
  logical               :: write_log

  ! We place the following data in static storage to keep the Myrinet MPI
  ! implementation from going ape over mpi buffers on the stack (causes
  ! problems between mpi and ifc, resulting in huge stalls in processing).

  integer, save :: my_node_dat(5)

  ! Gather times from the other processes to the master process.  The
  ! reference to gbl_loadbal_node_dat is only used in the master process.
  ! The dat collected from each node are the "direct force time" due to image
  ! nonbonded calcs, the "reciprocal force time" due to pme nonbonded
  ! reciprocal force calcs, the "other workload time", due to bond, angle,
  ! dihedral calc times and other atom-ownership based workload in the timespan
  ! that we attempt to loadbalance (there are other items of this sort on the
  ! side we basically can't simultaneously balance; bummer), list build time,
  ! and the send_atm_cnts total for each node, which can be used to determine
  ! when to do redistribution of the atom workload.  The reciprocal force
  ! component is artificially inflated by a factor of nrespa to keep
  ! pathological loadbalancing from occurring if respa is in effect.  For
  ! non-respa runs, this has no effect, as nrespa is 1.  For respa runs where
  ! all tasks have a slice of the reciprocal force workload, this has no effect.
  ! For respa runs where some tasks don't have a reciprocal force workload,
  ! this keeps the reciprocal force tasks from being given a load that they
  ! can't handle when they are actually doing reciprocal force calcs.  The
  ! list build times are currently included in the calcs as part of the
  ! image-associated workload; this probably works well in systems that are
  ! relatively near equilibrium with a fairly regular number of steps between
  ! listbuilds; we have excluded this in the past because it can have a
  ! destabilizing influence on timings.
 
  my_node_dat(1) = elapsed_100usec_dirfrc
  my_node_dat(2) = elapsed_100usec_recipfrc * nrespa
  my_node_dat(3) = elapsed_100usec_other
  my_node_dat(4) = elapsed_100usec_listbld
  my_node_dat(5) = my_send_atm_cnts_total / my_send_atm_cnts_sums

  call mpi_gather(my_node_dat, 5, mpi_integer, gbl_loadbal_node_dat, 5, &
                  mpi_integer, 0, mpi_comm_world, err_code_mpi)

  if (master) then

10 format(/, a, i7, a /)
20 format('  ', a, i4)
30 format(/, '  ', a)
40 format('    ', 8i9)

    if (loadbal_verbose .gt. 2) then

      write(logfile, 10) &
        'Image Distribution Calc Data at run step ', loadbal_step_ctr, ':'
      write(logfile, 20) 'Pairlist builds done between calls: ', &
                         my_send_atm_cnts_sums
      write(logfile, 30) &
        'Direct force time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(1, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Reciprocal force time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(2, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Bond-Angle-Dihedral force time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(3, i), i = 0, numtasks - 1)
#ifdef LISTBLD_USED
      write(logfile, 30) &
        'List build time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(4, i), i = 0, numtasks - 1)
#endif /* LISTBLD_USED */
      write(logfile, 30) &
        'Total force time (in units of 100 usec) for each task:'
      write(logfile, 40) &
       (gbl_loadbal_node_dat(1, i) + &
        gbl_loadbal_node_dat(2, i) + &
        gbl_loadbal_node_dat(3, i) + &
        gbl_loadbal_node_dat(4, i), i = 0, numtasks - 1)

    end if

    call calc_new_img_distribution(gbl_loadbal_node_dat, total_time)

    ! First turn off atom redistribution.

    gbl_img_div_tbl(numtasks + 1) = 0

    ! If image redistribution is stable (at max trigger), we check to see
    ! if atom redistribution is actually needed.  We also update the last
    ! count of atoms sent if an atom redistribution was recently done.

    if (img_redist_ctr .eq. &
        img_redist_trigger_max * nrespa / img_redist_trigger_divisor) then

      send_atm_cnts_total = 0

      do i = 0, numtasks - 1
        send_atm_cnts_total = send_atm_cnts_total + gbl_loadbal_node_dat(5,i)
      end do

!     write(0,*)'| DBG: send_atm_cnts_total = ', send_atm_cnts_total

      if (recalc_last_send_atm_cnts_total) then

        ! If the last time this was done, we did need an atom redistribution,
        ! save a new value for last_send_atm_cnts_total.  We weight the value
        ! by a factor of 9/10 toward its previous value, thus slowing the growth
        ! of the factor when the molecule configuration requires lots of atom
        ! sending, but still allowing for upward adjustment (prevents scenarios
        ! where a particularly favorable starting arrangement could result in
        ! constant atom redistribution).  This also "primes the pump" when we
        ! first execute this code (saves a value for last_send_atm_cnts_total).
        ! If the new value is lower than the old value, we just save it.

        if (last_send_atm_cnts_total .gt. 0) then
          if (last_send_atm_cnts_total .gt. send_atm_cnts_total) then
            last_send_atm_cnts_total = send_atm_cnts_total
          else
            last_send_atm_cnts_total = (9 * last_send_atm_cnts_total + &
                                        send_atm_cnts_total) / 10
          end if
        else
          last_send_atm_cnts_total = send_atm_cnts_total ! First save.
        end if

        recalc_last_send_atm_cnts_total = .false.

!       write(0,*)'| DBG: Saved send atom cnts total = ', &
!                 last_send_atm_cnts_total

      end if

      ! Check to see if an atom redistribution is needed.  If needed, this is
      ! delayed until the next cycle to allow for distribution of coordinates
      ! and velocities.  There is a debug facility.
      ! BUGBUG - We probably need to change the max trigger for the debug
      !          facility to actually be useful; otherwise too many cycles
      !          pass before the first atom redistribution.

      if (dbg_atom_redistribution .eq. 0) then
        if (force_atm_redist .or. dble(send_atm_cnts_total) .ge. &
            dble(last_send_atm_cnts_total) * atm_redist_tol) then
              gbl_img_div_tbl(numtasks + 1) = 1
        end if
      else
!       write(0,*)'| DBG: ATOM REDISTRIBUTION DEBUGGING ENABLED!'
        gbl_img_div_tbl(numtasks + 1) = 1
      end if

    end if

    call mpi_bcast(gbl_img_div_tbl, size(gbl_img_div_tbl), mpi_integer, 0, &
                     mpi_comm_world, err_code_mpi)

    write_log = .false.

    if (.not. initial_img_div_logged) then
      if (img_redist_ctr .eq. &
          img_redist_trigger_max / img_redist_trigger_divisor) then
        if (fft_slab_redist_enabled) then
          write_log = .true.
          initial_img_div_logged = .true.
        end if
      end if
    end if
    if (loadbal_verbose .gt. 0) then
      write_log = .true.
    end if

    if (write_log) then
        write(logfile, 110) &
          'Image Distribution at run step ', loadbal_step_ctr, ':'
        write(logfile, 120) 'Count of images assigned to each task:'
        write(logfile, 130) &
          (gbl_img_div_tbl(i + 1) - gbl_img_div_tbl(i), i = 0, numtasks - 1)
    end if

110 format(/, a, i7, a /)
120 format('  ', a)
130 format('    ', 8i9)

  else

    call mpi_bcast(gbl_img_div_tbl, size(gbl_img_div_tbl), mpi_integer, 0, &
                   mpi_comm_world, err_code_mpi)

  end if

  if (gbl_img_div_tbl(numtasks + 1) .ne. 0) then
    atm_redist_needed = .true.
    force_atm_redist = .false.
  end if

  my_img_lo = gbl_img_div_tbl(mytaskid) + 1
  my_img_hi = gbl_img_div_tbl(mytaskid + 1)

  elapsed_100usec_dirfrc = 0
  elapsed_100usec_recipfrc = 0
  elapsed_100usec_other = 0
  elapsed_100usec_listbld = 0
  my_send_atm_cnts_total = 0
  my_send_atm_cnts_sums = 0
  last_loadbal_step_ctr = loadbal_step_ctr

  return

end subroutine do_img_redistribution

!*******************************************************************************
!
! Subroutine:  calc_new_img_distribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine calc_new_img_distribution(node_dat, total_time)

  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: node_dat(5, 0:numtasks - 1) ! used by master
  double precision      :: total_time

! Local variables:

  integer               :: img_cnt
  integer               :: node
  integer               :: images_per_node(0:numtasks - 1) ! final result
  double precision      :: avg_time_per_img
  double precision      :: factor
  double precision      :: total_images
  double precision      :: target_time
  double precision      :: target_node_adj_time
  double precision      :: dirfrc_avg_time
  double precision      :: dirfrc_time_per_image
  double precision      :: dirfrc_total_time
  integer               :: max_total_time
  integer               :: min_total_time
  integer               :: dirfrc_limit         ! below this, use different
                                                ! algorithm for assigning
                                                ! images to node
  integer               :: images_limit         ! ditto...
  double precision      :: max_time_variation
  double precision      :: old_images_per_node(0:numtasks - 1)
  double precision      :: new_images_per_node(0:numtasks - 1)
  integer               :: total_times(0:numtasks - 1)

  ! Get total time in last run and average (target) time you want each 
  ! node to use.

  total_times(:) = node_dat(1, :) + &
                   node_dat(2, :) + &
                   node_dat(3, :) + &
                   node_dat(4, :)

  max_total_time = maxval(total_times)
  min_total_time = minval(total_times)

  total_time = 0.d0

  do node = 0, numtasks - 1
    total_time = total_time + dble(total_times(node))
  end do

  ! target_time is passed back for use in the caller as an average task
  ! time (for recip workload balancing)

  target_time = total_time / dble(numtasks)

  max_time_variation = dble(max_total_time - min_total_time) * &
                       dble(numtasks) / total_time

  ! We only readjust the image division if it is worthwhile.

! write(0,*)'| DBG: max img redist time variation = ', max_time_variation

  if (max_time_variation .gt. img_redist_time_tol) then

    ! Load balancing is handled differently on nodes with less than 5% the
    ! average workload because the standard algorithm does not work well at
    ! very low times / image counts (to say nothing of the possibility of
    ! dividing by 0...)

    dirfrc_total_time = 0.d0
    do node = 0, numtasks - 1
      dirfrc_total_time = dirfrc_total_time + dble(node_dat(1, node) + &
                                                   node_dat(4, node))
    end do

    dirfrc_avg_time = dirfrc_total_time / dble(numtasks)
    dirfrc_limit = ceiling(0.05d0 * dirfrc_avg_time)
    dirfrc_time_per_image = dirfrc_total_time / dble(natom)

    images_limit = ceiling(0.05d0 * dble(natom) / dble(numtasks))

    ! Fill in old images per node array:

    do node = 0, numtasks - 1
      old_images_per_node(node) = dble(gbl_img_div_tbl(node + 1) - &
                                       gbl_img_div_tbl(node))
    end do

    ! Calculate new image distributions. We correct for the nonadjustable 
    ! (recip and bnd-angle-dihed) times, make an estimate of processing time
    ! required per image based on the last cycle, and interpolated between old
    ! and new calculated values in order to damp the fluctuations, unless we
    ! are below limits for either the direct force time or image count for the
    ! node.  In that case, we use the average time required per image and
    ! the available time for direct force processing to estimate the correct
    ! number of images.

    total_images = 0.d0

    do node = 0, numtasks - 1

      target_node_adj_time = target_time - &
                             dble(node_dat(2, node) + node_dat(3, node))

      if (old_images_per_node(node) .ge. images_limit .and. &
          node_dat(1, node) + node_dat(4, node) .ge. dirfrc_limit) then
        factor = target_node_adj_time / &
                 dble(node_dat(1, node) + node_dat(4, node))
        if (factor .gt. 0.d0) then
          new_images_per_node(node) = factor * old_images_per_node(node)
        else
          new_images_per_node(node) = 0.d0
        end if
      else
        if (target_node_adj_time .gt. 0.d0) then
          new_images_per_node(node) = &
            target_node_adj_time / dirfrc_time_per_image
        else
          new_images_per_node(node) = 0.d0
        end if
      end if

      new_images_per_node(node) = (new_images_per_node(node) + &
                                   old_images_per_node(node)) / 2.d0
      total_images = total_images + new_images_per_node(node)
    end do

    ! Now basically scale the results to get the correct image total as an
    ! integer.  This can still be incorrect after all this scaling and
    ! adjustment, so we have to make a final adjustment pass on the integers,
    ! making sure the total is correct and none are 0.

    factor = dble(natom) / dble(total_images)
  
    img_cnt = 0

    do node = 0, numtasks - 1
      images_per_node(node) = nint(new_images_per_node(node) * factor)
      img_cnt = img_cnt + images_per_node(node)
    end do

    if (img_cnt .lt. natom) then

!     write(0,*)'DBG: Img redist - correcting low img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .lt. natom)

        images_per_node(node) = images_per_node(node) + 1

        img_cnt = img_cnt + 1

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    else if (img_cnt .gt. natom) then

!     write(0,*)'DBG: Img redist - correcting high img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .gt. natom)

        if (images_per_node(node) .gt. 1) then
          images_per_node(node) = images_per_node(node) - 1
          img_cnt = img_cnt - 1
        end if

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    end if

    gbl_img_div_tbl(0) = 0
    gbl_img_div_tbl(numtasks) = natom

    do node = 0, numtasks - 2
      gbl_img_div_tbl(node + 1) = gbl_img_div_tbl(node) + images_per_node(node)
    end do

  end if

  return

end subroutine calc_new_img_distribution

!*******************************************************************************
!
! Subroutine:  do_atm_redistribution
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine do_atm_redistribution(atm_cnt, write_log)

  use angles_mod
  use bonds_mod
  use cit_mod
  use dihedrals_mod
  use dynamics_mod
  use dynamics_dat_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use pmemd_lib_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use shake_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  logical               :: write_log
    
! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: num_ints, num_reals  ! values are not kept.
  integer               :: extra_used_atms(atm_cnt)
  integer               :: use_atm_map(atm_cnt)

  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)
  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)

! Divide atoms up among the processors.  The atom division is redone
! periodically under cit, and is either residue or molecule-based, with
! locality.  In other words, under cit a contiguous block of atoms owned by
! each process is a thing of the past.

! if (master) write(0,*)'| DBG: Doing atom redistribution.'

  num_ints = 0
  num_reals = 0

  extra_used_atms(:) = 0
  use_atm_map(:) = 0

  call get_fract_crds(atm_cnt, atm_crd, fraction)

  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  call divide_atoms(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst, write_log)

  ! BUGBUG - Must have different mechanism for Amoeba.  At present, though
  !          Amoeba does not do loadbalancing...

  call bonds_setup(num_ints, num_reals, use_atm_map)
  call angles_setup(num_ints, num_reals, use_atm_map)
  call dihedrals_setup(num_ints, num_reals, use_atm_map)

  extra_used_atm_cnt = 0

  do i = 1, atm_cnt
    if (use_atm_map(i) .ne. 0) then
      if (gbl_atm_owner_map(i) .ne. mytaskid) then
        extra_used_atm_cnt = extra_used_atm_cnt + 1
        extra_used_atms(extra_used_atm_cnt) = i
      end if
    end if
  end do

  if (extra_used_atm_cnt .gt. 0) then

    if (allocated(gbl_extra_used_atms)) then
      if (size(gbl_extra_used_atms) .lt. extra_used_atm_cnt) then
        deallocate(gbl_extra_used_atms)
        allocate(gbl_extra_used_atms(extra_used_atm_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
      end if
    else
      allocate(gbl_extra_used_atms(extra_used_atm_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if

    gbl_extra_used_atms(1:extra_used_atm_cnt) = &
      extra_used_atms(1:extra_used_atm_cnt)

  end if

  ! These allocations may need to be redone:

  if (size(gbl_recv_atm_lsts, 1) .lt. my_atm_cnt) then
    deallocate(gbl_recv_atm_lsts)
    allocate(gbl_recv_atm_lsts(my_atm_cnt, numtasks - 1), &
             stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
  end if

  gbl_recv_atm_lsts(:,:) = 0

#ifdef SLOW_NONBLOCKING_MPI
#else
  call set_minimum_mpi_bufs_size(max(3 * atm_cnt, &
                                 3 * my_atm_cnt * (numtasks-1)),&
                                 num_reals)
#endif

! Redo the shake per-process setup:

  call claim_my_fastwat_residues(num_ints, num_reals)
  call claim_my_nonfastwat_bonds(num_ints, num_reals)

! Recreate the entire molecule COM array.  You don't actually need the whole
! thing under mpi, but you have all the crds from which it is derived, and
! this code should execute very infrequently.

  if (ntp .gt. 0 .and. imin .eq. 0) then
    call get_all_mol_com(nspm, atm_crd, atm_mass, gbl_mol_atms, &
                         gbl_mol_mass_inv, gbl_mol_com)
  end if

! Clear the force array.  This avoids potential problems with non-zero data
! for atoms we don't own in the nmr routines and possibly elsewhere.

  atm_frc(:,:) = 0.d0

  return

end subroutine do_atm_redistribution

!*******************************************************************************
!
! Subroutine:  divide_atoms
!
! Description:  Set up atom workload division for parallel processing.
!              
!*******************************************************************************

subroutine divide_atoms(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst, &
                        write_log)

  use cit_mod
  use dynamics_dat_mod
  use pme_fft_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)
  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)
  logical               :: write_log


! Local variables:

  integer               :: atm_id, nxt_atm_id
  integer               :: i, j, k
  integer               :: frag_idx
  integer               :: frag_mol_idx
  integer               :: img_id, img_lo, img_hi
  integer               :: nxt_idx
  integer               :: mol_idx, res_idx
  integer               :: task_id 
  integer               :: added_atms_cnt
  integer, save         :: distrib_cnt = 0
  integer               :: target_cnt(0:numtasks - 1)

  integer               :: img_owner_map(atm_cnt)

  integer               :: atm_img_map(atm_cnt) ! NOTE local copy...

  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  distrib_cnt = distrib_cnt + 1 ! count of calls to this routine

  ! Initialize the atom owner map to values that will catch bugs:

  gbl_atm_owner_map(:) = -1

  ! Make the image owner map, as well as the target cnts:

  do task_id = 0, numtasks - 1
    img_lo = gbl_img_div_tbl(task_id) + 1
    img_hi = gbl_img_div_tbl(task_id + 1)
    img_owner_map(img_lo:img_hi) = task_id
    target_cnt(task_id) = img_hi - img_lo + 1
  end do

  img_hi = 0

  do k = 0, cit_tbl_z_dim - 1
    do j = 0, cit_tbl_y_dim - 1
      do i = 0, cit_tbl_x_dim - 1

        nxt_idx = crd_idx_lst_tbl(i, j, k)

        if (nxt_idx .ne. 0) then

          img_hi = img_hi + 1
          img_lo = img_hi

          do
            atm_img_map(atm_lst(nxt_idx)%idx) = img_hi
            nxt_idx = atm_lst(nxt_idx)%nxt
            if (nxt_idx .ne. 0) then
              img_hi = img_hi + 1
            else
              exit
            end if
          end do

          crd_idx_tbl(i, j, k)%img_lo = img_lo
          crd_idx_tbl(i, j, k)%img_hi = img_hi

          nxt_idx = crd_idx_lst_tbl(i, j, k)

        else
          crd_idx_tbl(i, j, k)%img_lo = 0
          crd_idx_tbl(i, j, k)%img_hi = -1
        end if

      end do
    end do
  end do

  ! Find the first atom in residues for CV and molecules for CP.  Map this atom
  ! to its image, find the image owner, and claim all the atoms in the residue
  ! or molecule for the owner.  Make my_atm_lst and my_mol_lst as appropriate.

  my_atm_cnt = 0
  my_mol_cnt = 0
  my_frag_mol_cnt = 0
  frag_mol_idx = 0
  gbl_vec_rcvcnts(:) = 0

#ifdef AMOEBA
  ! We can treat Amoeba as if it were a NVT run because it uses atom-based
  ! pressure scaling.

  if (ntp .eq. 0 .or. iamoeba .ne. 0) then  ! Constant volume, or amoeba.
#else
  if (ntp .eq. 0) then  ! Constant volume.
#endif /* AMOEBA */

    ! Note that single H residues are now disallowed in the input, since they
    ! serve no real purpose and are a pain to deal with.  Thus we don't need
    ! to check for them here.

    do res_idx = 1, nres

      atm_id = gbl_res_atms(res_idx)
      nxt_atm_id = gbl_res_atms(res_idx + 1)

      img_id = atm_img_map(atm_id)
      task_id = img_owner_map(img_id)

      if (task_id .ne. mytaskid) then

        gbl_vec_rcvcnts(task_id) = gbl_vec_rcvcnts(task_id) + &
                                   nxt_atm_id - atm_id

        do while (atm_id .lt. nxt_atm_id)
          gbl_atm_owner_map(atm_id) = task_id
          atm_id = atm_id + 1
        end do

      else

        do while (atm_id .lt. nxt_atm_id)
          gbl_atm_owner_map(atm_id) = task_id
          my_atm_cnt = my_atm_cnt + 1
          gbl_my_atm_lst(my_atm_cnt) = atm_id
          atm_id = atm_id + 1
        end do

      end if

    end do

    gbl_vec_rcvcnts(mytaskid) = my_atm_cnt

  else                  ! Constant pressure. Need molecule data.

    ! Constant pressure is a bear for a variety of reasons.  One big problem
    ! is uneven atom division caused by big molecules.  To counteract this
    ! effect, we have implemented a molecule fragment model, and we try to not
    ! exceed a target count by spilling over into other nodes when necessary.
    ! The target count is set to match the image count for a task. We do not
    ! overflow to tasks with a recip force workload if we have a choice.

!BEGIN DBG
!   if (master) write(0,*)'| DBG: Target node atm count =', target_cnt(:)
!END DBG

    if (recip_numtasks .eq. numtasks) then

      do mol_idx = 1, nspm
  
        atm_id = gbl_mol_atms(mol_idx)
        nxt_atm_id = gbl_mol_atms(mol_idx + 1)
        added_atms_cnt = nxt_atm_id - atm_id

        if (added_atms_cnt .gt. max_unfrag_mol_len) then

          frag_mol_idx = frag_mol_idx + 1

          do frag_idx = gbl_frag_mols(frag_mol_idx)%first_frag_idx, &
                        gbl_frag_mols(frag_mol_idx)%first_frag_idx + &
                        gbl_frag_mols(frag_mol_idx)%frag_cnt - 1

            atm_id = gbl_mol_frags(frag_idx)%first_atm_id
            added_atms_cnt = gbl_mol_frags(frag_idx)%atm_cnt
            nxt_atm_id = atm_id + added_atms_cnt
            img_id = atm_img_map(atm_id)
            task_id = img_owner_map(img_id)

            if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
                target_cnt(task_id)) then
  
              ! If there is no previous assignment to this task, go ahead and
              ! assign this; otherwise look for a more appropriate task.
  
              if (gbl_vec_rcvcnts(task_id) .ne. 0) then
                do i = 1, numtasks
                  task_id = task_id + 1
                  if (task_id .ge. numtasks) task_id = 0
                  if (gbl_vec_rcvcnts(task_id) .eq. 0 .or. &
                      gbl_vec_rcvcnts(task_id) + added_atms_cnt .le. &
                      target_cnt(task_id)) exit
                end do
              end if

            end if

            gbl_vec_rcvcnts(task_id) = gbl_vec_rcvcnts(task_id) + added_atms_cnt

            if (task_id .ne. mytaskid) then
              do while (atm_id .lt. nxt_atm_id)
                gbl_atm_owner_map(atm_id) = task_id
                atm_id = atm_id + 1
              end do
            else
              if (my_frag_mol_cnt .ne. 0) then
                if (gbl_my_frag_mol_lst(my_frag_mol_cnt) .ne. frag_mol_idx) then
                  my_frag_mol_cnt = my_frag_mol_cnt + 1
                  gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
                end if
              else
                my_frag_mol_cnt = my_frag_mol_cnt + 1
                gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
              end if
              do while (atm_id .lt. nxt_atm_id)
                gbl_atm_owner_map(atm_id) = task_id
                my_atm_cnt = my_atm_cnt + 1
                gbl_my_atm_lst(my_atm_cnt) = atm_id
                atm_id = atm_id + 1
              end do
            end if
            gbl_mol_frags(frag_idx)%owner = task_id

          end do

        else
  
          img_id = atm_img_map(atm_id)
          task_id = img_owner_map(img_id)
        
          if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
              target_cnt(task_id)) then
  
            ! If there is no previous assignment to this task, go ahead and
            ! assign this; otherwise look for a more appropriate task.
  
            if (gbl_vec_rcvcnts(task_id) .ne. 0) then
              do i = 1, numtasks
                task_id = task_id + 1
                if (task_id .ge. numtasks) task_id = 0
                if (gbl_vec_rcvcnts(task_id) .eq. 0 .or. &
                    gbl_vec_rcvcnts(task_id) + added_atms_cnt .le. &
                    target_cnt(task_id)) exit
              end do
            end if

          end if
            
          gbl_vec_rcvcnts(task_id) = gbl_vec_rcvcnts(task_id) + added_atms_cnt
  
          if (task_id .ne. mytaskid) then
            do while (atm_id .lt. nxt_atm_id)
              gbl_atm_owner_map(atm_id) = task_id
              atm_id = atm_id + 1
            end do
          else
            my_mol_cnt = my_mol_cnt + 1
            gbl_my_mol_lst(my_mol_cnt) = mol_idx
            do while (atm_id .lt. nxt_atm_id)
              gbl_atm_owner_map(atm_id) = task_id
              my_atm_cnt = my_atm_cnt + 1
              gbl_my_atm_lst(my_atm_cnt) = atm_id
              atm_id = atm_id + 1
            end do
          end if

        end if
  
      end do

    else

      do mol_idx = 1, nspm
  
        atm_id = gbl_mol_atms(mol_idx)
        nxt_atm_id = gbl_mol_atms(mol_idx + 1)
        added_atms_cnt = nxt_atm_id - atm_id

        if (added_atms_cnt .gt. max_unfrag_mol_len) then
  
          frag_mol_idx = frag_mol_idx + 1

          do frag_idx = gbl_frag_mols(frag_mol_idx)%first_frag_idx, &
                        gbl_frag_mols(frag_mol_idx)%first_frag_idx + &
                        gbl_frag_mols(frag_mol_idx)%frag_cnt - 1

            atm_id = gbl_mol_frags(frag_idx)%first_atm_id
            added_atms_cnt = gbl_mol_frags(frag_idx)%atm_cnt
            nxt_atm_id = atm_id + added_atms_cnt
            img_id = atm_img_map(atm_id)
            task_id = img_owner_map(img_id)

            if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
                target_cnt(task_id)) then
  
              ! Just go ahead and assign to the current task unless there is
              ! already some assignment or the current task has a recip
              ! workload.  Otherwise, look for a task that would not be
              ! overloaded.
  
              if (gbl_vec_rcvcnts(task_id) .ne. 0 .or. &
                  xy_slab_cnt(task_id) .ne. 0) then
              
                do i = 1, numtasks - 1
                  task_id = task_id + 1
                  if (task_id .ge. numtasks) task_id = 0
                  if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .le. &
                    target_cnt(task_id)) exit
                end do
  
                ! Okay, did we end up on a node where we overflowed anyway,
                ! and it also has a recip force workload?  If so, if there are
                ! any tasks without a recip force workload, overflow to them
                ! regardless of their target count.
  
                if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
                    target_cnt(task_id)) then
                  if (xy_slab_cnt(task_id) .gt. 0) then
                    task_id = img_owner_map(img_id)
                    do i = 1, numtasks
                      task_id = task_id + 1
                      if (task_id .ge. numtasks) task_id = 0
                      if (xy_slab_cnt(task_id) .le. 0) exit
                    end do
                  end if
                end if
  
              end if
            end if
            
            gbl_vec_rcvcnts(task_id) = gbl_vec_rcvcnts(task_id) + added_atms_cnt
  
            if (task_id .ne. mytaskid) then
              do while (atm_id .lt. nxt_atm_id)
                gbl_atm_owner_map(atm_id) = task_id
                atm_id = atm_id + 1
              end do
            else
              if (my_frag_mol_cnt .ne. 0) then
                if (gbl_my_frag_mol_lst(my_frag_mol_cnt) .ne. frag_mol_idx) then
                  my_frag_mol_cnt = my_frag_mol_cnt + 1
                  gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
                end if
              else
                my_frag_mol_cnt = my_frag_mol_cnt + 1
                gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
              end if
              do while (atm_id .lt. nxt_atm_id)
                gbl_atm_owner_map(atm_id) = task_id
                my_atm_cnt = my_atm_cnt + 1
                gbl_my_atm_lst(my_atm_cnt) = atm_id
                atm_id = atm_id + 1
              end do
            end if
            gbl_mol_frags(frag_idx)%owner = task_id
  
          end do

        else

          img_id = atm_img_map(atm_id)
          task_id = img_owner_map(img_id)

          if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
              target_cnt(task_id)) then
  
            ! Just go ahead and assign to the current task unless there is
            ! already some assignment or the current task has a recip workload.
            ! Otherwise, look for a task that would not be overloaded.
  
            if (gbl_vec_rcvcnts(task_id) .ne. 0 .or. &
                xy_slab_cnt(task_id) .ne. 0) then
              
              do i = 1, numtasks - 1
                task_id = task_id + 1
                if (task_id .ge. numtasks) task_id = 0
                if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .le. &
                    target_cnt(task_id)) exit
              end do
  
              ! Okay, did we end up on a node where we overflowed anyway, and it
              ! also has a recip force workload?  If so, if there are any tasks
              ! without a recip force workload, overflow to them regardless of
              ! their target count.
  
              if (gbl_vec_rcvcnts(task_id) + added_atms_cnt .gt. &
                  target_cnt(task_id)) then
                if (xy_slab_cnt(task_id) .gt. 0) then
                  task_id = img_owner_map(img_id)
                  do i = 1, numtasks
                    task_id = task_id + 1
                    if (task_id .ge. numtasks) task_id = 0
                    if (xy_slab_cnt(task_id) .le. 0) exit
                  end do
                end if
              end if
  
            end if
          end if
            
          gbl_vec_rcvcnts(task_id) = gbl_vec_rcvcnts(task_id) + added_atms_cnt
  
          if (task_id .ne. mytaskid) then
            do while (atm_id .lt. nxt_atm_id)
              gbl_atm_owner_map(atm_id) = task_id
              atm_id = atm_id + 1
            end do
          else
            my_mol_cnt = my_mol_cnt + 1
            gbl_my_mol_lst(my_mol_cnt) = mol_idx
            do while (atm_id .lt. nxt_atm_id)
              gbl_atm_owner_map(atm_id) = task_id
              my_atm_cnt = my_atm_cnt + 1
              gbl_my_atm_lst(my_atm_cnt) = atm_id
              atm_id = atm_id + 1
            end do
          end if

        end if
  
      end do

    end if

  end if

  gbl_atm_offsets(0) = 0

  do task_id = 0, numtasks - 1
    gbl_atm_offsets(task_id + 1) = gbl_atm_offsets(task_id) + &
                                   gbl_vec_rcvcnts(task_id)
  end do

!BEGIN DBG
! if (master) write(0,*)'| DBG: Node atom cnts =', gbl_vec_rcvcnts(:)

! do i = 1, atm_cnt
!   if (gbl_atm_owner_map(i) .lt. 0 .or. &
!       gbl_atm_owner_map(i) .ge. numtasks) &
!       write(0,*)'| DBG_ERR: atom ', i, 'assigned to bad task ', &
!                 gbl_atm_owner_map(i)
! end do
!END DBG

  do task_id = 0, numtasks - 1
    gbl_vec_rcvcnts(task_id) = 3 * gbl_vec_rcvcnts(task_id)
  end do

  gbl_vec_rcvcnts(numtasks) = 0

  gbl_vec_offsets(0) = 0

  do task_id = 0, numtasks - 1
    gbl_vec_offsets(task_id + 1) = gbl_vec_offsets(task_id) + &
                                   gbl_vec_rcvcnts(task_id)
  end do

  ! Set up for distributed molecule handling if needed. ALL tasks have to call
  ! this stuff, even if they don't participate in processing a particular, or
  ! any molecule.

  if (frag_mol_cnt .gt. 0) then
    call destroy_communicators()        ! only does something if comm's exist...
    call create_communicators()
  end if

  if (master) then
    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Atom Distribution No. ', distrib_cnt, ' at run step ', &
          loadbal_step_ctr, ':'
        write(logfile, 20) 'Count of atoms assigned to each task:'
        write(logfile, 30) (gbl_vec_rcvcnts(i)/3, i = 0, numtasks - 1)
    end if
  end if

10 format(/, a, i5, a, i7, a /)
20 format('  ', a)
30 format('    ', 8i9)

  return

end subroutine divide_atoms

!*******************************************************************************
!
! Subroutine:  divide_images_evenly
!
! Description:  Set up image boundaries for parallel processing.
!
!*******************************************************************************

subroutine divide_images_evenly(div_tbl, write_log)

  use img_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: div_tbl(0:)
  logical       :: write_log

! Local variables:

  double precision      :: fraction
  double precision      :: next_div
  integer               :: i
  
  ! The fraction is guaranteed .gt. 0 by an earlier check...

  fraction = dble(natom) / dble(numtasks)
  next_div = 0.d0
  div_tbl(0) = 0
  div_tbl(numtasks) = natom

  do i = 1, numtasks - 1
    next_div = next_div + fraction
    div_tbl(i) = dnint(next_div)
  end do

  my_img_lo = div_tbl(mytaskid) + 1
  my_img_hi = div_tbl(mytaskid + 1)

  if (master) then
    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Image Distribution at run step ', loadbal_step_ctr, ':'
        write(logfile, 20) 'Count of images assigned to each task:'
        write(logfile, 30) (div_tbl(i + 1) - div_tbl(i), i = 0, numtasks - 1)
    end if
  end if

10 format(/, a, i7, a /)
20 format('  ', a)
30 format('    ', 8i9)

  return

end subroutine divide_images_evenly

!*******************************************************************************
!
! Subroutine:  divide_images_recip_biased
!
! Description:  Set up image boundaries for parallel processing.
!
!*******************************************************************************

subroutine divide_images_recip_biased(recip_workload, div_tbl, write_log)

  use pme_fft_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: recip_workload
  integer               :: div_tbl(0:)
  logical               :: write_log

! Local variables:

  double precision              :: nonrecip_frac
  double precision              :: recip_frac
  double precision              :: next_div
  integer                       :: i

  ! If all tasks have a reciprocal workload, there is no bias in the
  ! assignment of workload.

  if (recip_numtasks .eq. numtasks) then
    call divide_images_evenly(div_tbl, write_log)
    return
  end if

  ! The direct workload is 1.d0.  We first apportion the direct workload to
  ! give all tasks an equal workload, based on the assumption that the input
  ! recip workload, expressed as a fraction of the direct workload, is correct.
  ! Then, if the direct workload assigned to the tasks handling a reciprocal
  ! workload is less than 10% of the total direct workload, we readjust the
  ! direct workload, distributing 10% of it to the tasks that also handle a
  ! reciprocal workload.  This adjustment will place a higher total workload
  ! on the reciprocal tasks, but it will be quickly adjusted down as
  ! appropriate.  This algorithm is used because direct force load balancing
  ! does not adapt quickly if the initially assigned direct force load is
  ! near 0.

  nonrecip_frac = (1.d0 + recip_workload) / dble(numtasks)
  recip_frac = nonrecip_frac - (recip_workload / dble(recip_numtasks))

  ! Division by 0 insured to not occur by above check for
  ! recip_numtasks .eq. numtasks.
  if (recip_frac * dble(recip_numtasks) / &
      (nonrecip_frac * dble(numtasks - recip_numtasks)) .lt. 0.1d0) then
    recip_frac = 0.1d0 / dble(recip_numtasks)
    nonrecip_frac = 0.9d0 / dble(numtasks - recip_numtasks)
  end if

  recip_frac = recip_frac * dble(natom)
  nonrecip_frac = nonrecip_frac * dble(natom)

  next_div = 0.d0
  div_tbl(0) = 0
  div_tbl(numtasks) = natom

  do i = 1, numtasks - 1
    if (xy_slab_cnt(i - 1) .eq. 0) then
      next_div = next_div + nonrecip_frac
    else
      next_div = next_div + recip_frac
    end if
    div_tbl(i) = dnint(next_div)
  end do

  my_img_lo = div_tbl(mytaskid) + 1
  my_img_hi = div_tbl(mytaskid + 1)

10 format(/, a, i7, a /)
20 format('  ', a, i4, a)
22 format('  ', a, f5.1, a)
30 format('  ', a)
40 format('    ', 8i9)

  if (master) then

    if (loadbal_verbose .gt. 2) then
      write(logfile, 10) &
        'Image Distribution Bias Data at run step ', loadbal_step_ctr, ':'
      write(logfile, 20) 'Recip force workload assigned to ', recip_numtasks, &
                         ' tasks'
      write(logfile, 22) 'Recip force workload estimate is ', &
            recip_workload * 100.d0, ' percent'
    end if

    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Image Distribution at run step ', loadbal_step_ctr, ':'
        write(logfile, 30) 'Count of images assigned to each task:'
        write(logfile, 40) (div_tbl(i + 1) - div_tbl(i), i = 0, numtasks - 1)
    end if

  end if


  return

end subroutine divide_images_recip_biased

#endif /* MPI */

end module loadbal_mod
