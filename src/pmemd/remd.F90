#include "copyright.i"

!*******************************************************************************
!
! Module: remd_mod
!
! Description: 
!
! Module for controlling replica exchange functionality. Many parts were adapted
! from sander's remd module for pmemd. The sander implementation was done by
! Daniel Roe, based on the original implementation by Guanglei Cui and Carlos
! Simmerling. 
!
! This is being designed such that exchanges can be done in an arbitrary number
! of dimensions. This is accomplished through the use of a 2-D "partners" array.
! The first dimension is the number of exchange dimensions that we have, and is
! of arbitrary size. The second dimension is the list of the partners that it 
! can exchange with (typically 2, one above and one below). MPI_Recv will ALWAYS
! be done by the *lower* partner (and it will then calculate the exchange prob).
! Note that the "partners" are indexed beginning from 1, whereas MPI ranks start
! from 0, so when using the partner array to indicate a replica rank, always
! subtract 1. Written by Jason Swails, 3/2011
!              
!*******************************************************************************

module remd_mod

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

  integer, save                      :: remd_method = 0 ! REMD method to use

  double precision, parameter :: BIG = 9999999.d0

#ifdef MPI
  
  ! Everything is private by default

! Variables
!
!  exchange_successes   : used to track acceptance ratio
!  exch_num             : number of exchanges we've done in each dimension
!  partners             : which replicas this one can attempt exchanges with
!                       : allocated like (remd_dimension, 2)
!  remd_dimension       : how many dimensions we will attempt to exchange in
!  num_from_min         : How many temps are less than this replica's temp
!  even_replica         : is our replica number "even"? Determines whether we
!  even_exchanges       : Does the even replica exchange or not?
!  crd_temp             : temporary coordinate array for hamiltonian exchanges
!  force_temp           : temporary force array for hamiltonian exchanges
!  temperatures         : array with temperatures for every replica
!                       : exchange or let our partner do it
!  total_left_fe        : calculated FEP free energy for exchanging to the left
!  total_right_fe       : calculated FEP free energy for exchanging to the right
!  num_right_exchg      : # of times we've exchanged to the right
!  num_left_exchg       : # of times we've exchanged to the left
!  remd_modwt           : Forces modwt to re-read the temperatures, etc., or the
!                         NMR restraint facility will overwrite temp0

! General for all REMD methods

  integer, allocatable, save         :: exchange_successes(:,:)
  integer, allocatable, save         :: exch_num(:) ! how many exchgs we've done
  
  integer, allocatable, save         :: partners(:,:)
  integer, save                      :: remd_dimension

  logical, allocatable, save         :: even_replica(:)
  logical, allocatable, save         :: even_exchanges(:)

  logical, save                      :: remd_modwt
  
! Specific to T-REMD

  integer, save                      :: num_from_min

  double precision, allocatable      :: temperatures(:)

! Specific to H-REMD

  double precision, allocatable      :: crd_temp(:,:)

  double precision, save             :: total_left_fe
  double precision, save             :: total_right_fe

  integer, save                      :: num_right_exchg
  integer, save                      :: num_left_exchg

contains

!*******************************************************************************
!
! Subroutine: remd_setup
!
! Description: Sets up the REMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine remd_setup(numexchg)

  use file_io_mod,       only : amopen
  use gbl_constants_mod
  use mdin_ctrl_dat_mod, only : ig, temp0
  use pmemd_lib_mod,     only : mexit

  implicit none

  ! Passed variables

  integer, intent(in) :: numexchg

  ! Local variables

  integer             :: alloc_failed

  ! Force an initial re-reading of the variables (doesn't hurt)

  remd_modwt = .true.

  ! Open remlog file and write some initial info

  if (master_master) then

    call amopen(remlog, remlog_name, owrite, 'F', 'W')

    write(remlog, '(a)')     '# Replica Exchange log file'
    write(remlog, '(a,i10)') '# numexchg is ', numexchg
    write(remlog, '(a)')     '# REMD filenames:'
    write(remlog, '(a,a)')   '#   remlog= ', trim(remlog_name)
    write(remlog, '(a,a)')   '#   remtype= ', trim(remtype_name)

  end if

  ! If a REMD exchange definition file was specified, open that and set up
  ! the partners. Otherwise, define them the way they'd be defined in sander

  if (len_trim(remd_partners_name) .gt. 0) then
    call set_partners_from_file
  else
  
    if (remd_method .eq. 1) then
      call temp_remd_setup
    else if (remd_method .eq. 3) then
      call h_remd_setup
    end if

  end if

end subroutine remd_setup

!*******************************************************************************
!
! Subroutine: slave_remd_setup
!
! Description: Performs necessary setup for slave nodes
!
!*******************************************************************************

subroutine slave_remd_setup
  
  use prmtop_dat_mod, only : natom

  implicit none

  integer :: alloc_failed

  if (remd_method .eq. 3) then

    allocate(crd_temp(3, natom), stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('slave_remd_setup', 'Error allocating crd_temp array')

  end if

end subroutine slave_remd_setup

!*******************************************************************************
!
! Subroutine: temp_remd_setup
!
! Description: This defines the partners array by sorting the temperatures. This
!              subroutine allocates partners(:,:), so don't call this with the
!              other partner subroutines!
!
!*******************************************************************************

subroutine temp_remd_setup

  use gbl_constants_mod, only : error_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

! Local variables

  double precision   :: my_temp       ! temperature of my replica

  integer            :: alloc_failed
  integer            :: i

  remd_dimension = 1 ! only exchange in T-space

  ! Allocate temperature table and exchange successes, and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            exch_num(remd_dimension),                      &
            temperatures(numgroups),                       &
            even_replica(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('temp_remd_setup', 'allocation error')

  exchange_successes(:,:) = 0
  exch_num(:) = 1
  even_exchanges(:) = .false.

  ! Set up the temperature table

  call collect_temps

  ! Write some T-REMD specific rem.log info

  if (master_master) &
    write(remlog, '(a)') '# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, &
                         &Success rate (i,i+1), ResStruct#'

  my_temp = temperatures(master_rank)

  ! First thing we want to do is allocate partners

  allocate(partners(1,2), stat = alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('temp_remd_setup', 'error allocating partner array')

  ! Check to make sure we have no duplicate temperatures. Then set the partners
  ! array via a call to set_temp_partners

  do i = 1, numgroups

    if (temperatures(i) .eq. my_temp .and. i .ne. master_rank) then
      write(mdout, '(a,a)') error_hdr, 'two temperatures are identical &
        &in temp_remd_setup'
    end if
    
  end do

  call set_temp_partners(1)

  return

end subroutine temp_remd_setup

!*******************************************************************************
!
! Subroutine: h_remd_setup
!
! Description: Sets up Hamiltonian-REMD in which coordinates are traded between
!              replicas.
!
!*******************************************************************************

subroutine h_remd_setup

  use prmtop_dat_mod,    only : natom

  implicit none

  integer :: alloc_failed

  remd_dimension = 1 ! only exchange in H-space
  num_right_exchg = 0
  num_left_exchg = 0
  total_left_fe = 0.d0
  total_right_fe = 0.d0

  ! First allocate necessary data structures and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            exch_num(remd_dimension),                      &
            even_replica(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            crd_temp(3, natom),                            &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) call alloc_error('h_remd_setup', 'allocation error')

  exchange_successes(:,:) = 0
  exch_num(:) = 1
  even_exchanges(:) = .false.

  ! Write the header line for the rem.log

  if (master_master) &
    write(remlog, '(a)') '# Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), left_fe,&
                        & right_fe, Success, Success rate (i,i+1)'

  ! The partners are just the adjacent neighbors here

  call set_adjacent_partners

end subroutine h_remd_setup

!*******************************************************************************
!
! Subroutine: alloc_error
!
! Description: Called in the case of an allocation error
!
!*******************************************************************************

subroutine alloc_error(routine, message)

  use gbl_constants_mod, only : error_hdr, extra_line_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

  ! Passed variables

  character(*), intent(in) :: routine, message

  write(mdout, '(a,a,a)') error_hdr, 'Error in ', routine
  write(mdout, '(a,a)') extra_line_hdr, message

  call mexit(mdout, 1)

end subroutine alloc_error

!*******************************************************************************
!
! Subroutine: set_partners_from_file
!
! Description: Parses remd_partners file to define the partners array
!
!*******************************************************************************

subroutine set_partners_from_file

  use pmemd_lib_mod, only : mexit

  implicit none

  ! This is not implemented yet

  write(mdout, '(a)') 'partner file is not implemented yet!'
  call mexit(mdout, 1)

  return

end subroutine set_partners_from_file

!*******************************************************************************
!
! Subroutine: set_temp_partners
!
! Description: When temperatures change, we need to reassign partners. Thus,
!              make this an easy subroutine to call. Should be cheaper than
!              collective communication (broadcasting the temperature table).
!              Thus, here we just shift our indices in the partners array
!              instead of collecting, sorting, etc. all temperatures each step
!
!*******************************************************************************

subroutine set_temp_partners(t_dim)

  use pmemd_lib_mod, only : mexit

  implicit none

! Passed variables

  integer, intent(in) :: t_dim ! which dimension in partners is Temp

! Local variables

  double precision    :: minimum_temp  ! Minimum temp in temperature array
  double precision    :: maximum_temp  ! Maximum temp in temperature array

  double precision    :: step_up       ! Lowest temp greater than my temp
  double precision    :: step_down     ! Highest temp lower than my temp

  double precision    :: my_temp       ! temperature of my replica

  integer             :: min_idx       ! index of the max temp
  integer             :: max_idx       ! index of the min temp
  integer             :: step_up_idx   ! index of the step down
  integer             :: step_down_idx ! index of the step up

  integer             :: i

  my_temp = temperatures(master_rank + 1)

  step_up = BIG
  step_down = -BIG
  step_up_idx = 0
  step_down_idx = 0

  maximum_temp = my_temp
  minimum_temp = my_temp
  min_idx = master_rank + 1
  max_idx = master_rank + 1

  num_from_min = 0

  do i = 1, numgroups

    if (temperatures(i) .gt. maximum_temp) then
      maximum_temp = temperatures(i)
      max_idx = i
    else if (temperatures(i) .lt. minimum_temp) then
      minimum_temp = temperatures(i)
      min_idx = i
    end if

    if (temperatures(i) .gt. my_temp .and. &
        temperatures(i) .lt. step_up) then

      step_up = temperatures(i)
      step_up_idx = i

    else if (temperatures(i) .lt. my_temp) then
      num_from_min = num_from_min + 1
      
      if (temperatures(i) .gt. step_down) then
        step_down = temperatures(i)
        step_down_idx = i
      end if

    end if

  end do

  ! Now assign the partners

  if (my_temp .eq. maximum_temp) then
    partners(t_dim,1) = step_down_idx
    partners(t_dim,2) = min_idx
  else if (my_temp .eq. minimum_temp) then
    partners(t_dim,1) = max_idx
    partners(t_dim,2) = step_up_idx
  else
    partners(t_dim,1) = step_down_idx
    partners(t_dim,2) = step_up_idx
  end if

  ! Now define if we are an "even" replica or not. The minimum temp will be
  ! 0, next highest will be 1, etc. This will allow us to distinguish b/w
  ! the replicas that should be doing the send or receive each step

  even_replica(1) = mod(num_from_min, 2) .eq. 0

  return

end subroutine set_temp_partners

!*******************************************************************************
!
! Subroutine: set_adjacent_partners
!
! Description: Assign partners according to the order they're listed in the
!              groupfile
!
!*******************************************************************************

subroutine set_adjacent_partners

  implicit none

  integer alloc_failed

  ! First we allocate partners and set it. Note that master_rank spans 0 -
  ! numgroups - 1. Therefore, this replica is indexed master_rank + 1. The
  ! 'previous' partner should be indexed master_rank and the 'next' partner 
  ! should be indexed master_rank + 2

  allocate(partners(1,2), stat=alloc_failed)

  if (master_rank .eq. 0) then
    partners(1,1) = numgroups
    partners(1,2) = master_rank + 2
  else if (master_rank .eq. numgroups - 1) then
    partners(1,1) = master_rank
    partners(1,2) = 1
  else
    partners(1,1) = master_rank
    partners(1,2) = master_rank + 2
  end if

  even_replica(1) = mod(master_rank, 2) .eq. 0

  return

end subroutine set_adjacent_partners

!*******************************************************************************
!
! Subroutine: rescale_velocities
!
! Description: rescale the velocities for the new temperature
!
!*******************************************************************************

subroutine rescale_velocities(atm_cnt, vel, vel_scale_factor)

  use mdin_ctrl_dat_mod, only : temp0

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: vel_scale_factor

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(a,f8.3,a,f8.3)') 'REMD: scaling velocities by ', &
      vel_scale_factor, ' to match new bath T ', temp0
#endif

#ifdef CUDA
  call gpu_scale_velocities(vel_scale_factor)
#else
  vel(:, :) = vel(:, :) * vel_scale_factor
#endif

end subroutine rescale_velocities

!*******************************************************************************
!
! Subroutine: bcast_remd_method
!
! Description: Broadcasts the replica exchange method we're using so everyone in
!              the world knows what it is
!
!*******************************************************************************

subroutine bcast_remd_method

  implicit none

  call mpi_bcast(remd_method, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)

  return

end subroutine bcast_remd_method

!*******************************************************************************
!
! Subroutine: collect_temps
!
! Description: Collects new temperatures in case they were changed by NMR
!              restraints
!
!*******************************************************************************

subroutine collect_temps

  use mdin_ctrl_dat_mod, only : temp0

  implicit none

  if (.not. master) return

  temperatures(:) = 0.d0

  call mpi_allgather(   temp0,     1, mpi_double_precision, &
                     temperatures, 1, mpi_double_precision, &
                     pmemd_master_comm, err_code_mpi)
  return

end subroutine collect_temps

!*******************************************************************************
!
! Subroutine: remd_cleanup
!
! Description: Cleans up after a REMD run. It closes the REMD files and 
!              deallocates arrays
!
!*******************************************************************************

subroutine remd_cleanup

  implicit none

! Return if no REMD was run

  if (remd_method .eq. 0) return

! If REMD was run, deallocate any of the allocated data structures

  if (allocated(exchange_successes)) deallocate(exchange_successes)

  if (allocated(partners)) deallocate(partners)

  if (allocated(crd_temp)) deallocate(crd_temp)

  if (allocated(temperatures)) deallocate(temperatures)

  if (allocated(exch_num)) deallocate(exch_num)

  if (allocated(even_replica)) deallocate(even_replica)

  if (allocated(even_exchanges)) deallocate(even_exchanges)

! Close any files that were opened
  
  close(remlog)

end subroutine remd_cleanup
#endif /* MPI */

end module remd_mod
