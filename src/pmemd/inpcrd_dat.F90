!*******************************************************************************
!
! Module: inpcrd_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module inpcrd_dat_mod

  use file_io_dat_mod
#ifdef MPI
  use remd_mod, only : remd_method
#endif

  implicit none

  ! atm_crd = the atom xyz coordinates as a 2d array.
  ! atm_frc = the atom xyz forces as a 2d array. Not read here, but
  !           kept here as part of the fundamental data.
  ! atm_vel = the atom velocity xyz vectors as a 2d array.

  double precision, allocatable, save   :: atm_crd(:,:)
  double precision, allocatable, save   :: atm_frc(:,:)
  double precision, allocatable, save   :: atm_vel(:,:)
  double precision, allocatable, save   :: atm_last_vel(:,:)

! Hide internal routines:

  private       alloc_inpcrd_mem

contains

!*******************************************************************************
!
! Subroutine:  init_inpcrd_dat
!
! Description: Read "old" style inpcrd, without section tags.  These files are
!              pre-amber 9 compliant, and compliant for non-amoeba amber 9 and
!              10.
!              
!*******************************************************************************

subroutine init_inpcrd_dat(num_ints, num_reals, inpcrd_natom, &
                           inpcrd_alpha, inpcrd_beta, inpcrd_gamma, &
                           inpcrd_box, tt, title)

  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use binrestart_mod, only: read_nc_restart, read_nc_restart_atoms, check_nc_restart

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     
  integer, intent(out)          :: inpcrd_natom
  double precision, intent(out) :: inpcrd_alpha
  double precision, intent(out) :: inpcrd_beta
  double precision, intent(out) :: inpcrd_gamma
  double precision, intent(out) :: inpcrd_box(3)
  double precision, intent(out) :: tt
  character(80),    intent(out) :: title

! Local variables:

  logical               :: formatted_input
  logical               :: netcdf_input
  integer               :: i, line_num
  character(80)         :: read_buf                ! For format checking
  integer               :: inpcrd_version          ! For pmemd 2.01 format
  integer               :: zero_natom              ! For pmemd 2.01 format
  logical               :: box_found
  logical               :: angles_found
  logical               :: velocities_found
  logical               :: velocities_needed
  double precision      :: local_temp0        ! for REMD restart

  angles_found = .false.
  box_found = .false.
  velocities_found = .false.
  local_temp0 = 0.d0
  velocities_needed = (ntx==4 .or. ntx==5 .or. ntx==6 .or. ntx==7)

! Open the inpcrd file:

  formatted_input =  (ntx .eq. 1 .or. ntx .eq. 5 .or. ntx .eq. 7)
  netcdf_input = check_nc_restart(inpcrd_name) 

  inpcrd_box(:) = 0.d0
  inpcrd_alpha = 0.d0
  inpcrd_beta = 0.d0
  inpcrd_gamma = 0.d0

  ! ---------- NetCDF restart ----------
  if (netcdf_input) then
    ! Get number of atoms for memory allocation
    call read_nc_restart_atoms(inpcrd_name, inpcrd_natom)
    ! Allocate memory
    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)
    ! Read the entire restart
    call read_nc_restart(inpcrd_name,title,inpcrd_natom,atm_crd,atm_vel,local_temp0,tt,&
                         inpcrd_box,inpcrd_alpha,inpcrd_beta,inpcrd_gamma,&
                         box_found,velocities_found)
#ifdef MPI
    ! If we're restarting a REMD trajectory, then replace the temp0 from the
    ! mdin file with the temperature from the inpcrd file
    if (remd_method .gt. 0 .and. local_temp0 .gt. 0 .and. irest .eq. 1) then
      write(mdout, '(a)') '| Overwriting temp0 from mdin with temp0 from &
                           &netcdf inpcrd file'
      temp0 = local_temp0
    end if
#endif
  ! ---------- Formatted (ASCII) restart ----------
  else if (formatted_input) then 
    call amopen(inpcrd, inpcrd_name, 'O', 'F', 'R')

    ! Sander 7/8 uses a mechanism of checking if the 6th char is blank to see if
    ! this is an old i5 format sander 6 file or a new i6 format sander 7 file.
    ! We also check for the versioned pmemd 2.01 format, with an initial natom
    ! of 0.  This was really the best approach, but unfortunately sander 7/8
    ! did not implement something like it.

    read(inpcrd, 9008) title
    read(inpcrd, '(a80)') read_buf
    rewind(inpcrd)
    read(inpcrd, 9008) title

    if (read_buf(6:6) .eq. ' ') then ! Sander 6 or pmemd 2.01 format...

       read(inpcrd, '(i5)', err=666) inpcrd_natom
       rewind(inpcrd)
       read(inpcrd, 9008) title

       if (inpcrd_natom .ne. 0) then
         read(inpcrd, 9018, err=666) inpcrd_natom, tt, local_temp0
       else
         read(inpcrd, '(2i5, i10, 4e15.7)') &
              zero_natom, inpcrd_version, inpcrd_natom, tt
         if (inpcrd_version .ne. 2) then
           write(mdout, '(a,a)') error_hdr, &
             'Unrecognized format for inpcrd/restrt file!'
           call mexit(6, 1)
         end if
       end if
    else if (read_buf(7:7) .eq. ' ') then ! Sander 7/8/9/10 large system format...
      read(inpcrd, 9019, err=666) inpcrd_natom, tt, local_temp0
    else if (read_buf(8:8) .eq. ' ') then ! Sander 11 - 1 mil+ format
      read(inpcrd, 9020, err=666) inpcrd_natom, tt, local_temp0
    else        ! assume amber 11 VERY large system format. 10 mil+
      read(inpcrd, 9021, err=666) inpcrd_natom, tt, local_temp0
    end if

#ifdef MPI
    ! If we're restarting a REMD trajectory, then replace the temp0 from the
    ! mdin file with the temperature from the inpcrd file

    if (remd_method .gt. 0 .and. local_temp0 .gt. 0 .and. irest .eq. 1) then
      write(mdout, '(a)') '| Overwriting temp0 from mdin with temp0 from &
                           &inpcrd file'
      temp0 = local_temp0
    end if
#endif
    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)

    ! We should always find coordinates:

    read(inpcrd, 9028, end=1000, err=1001) (atm_crd(1:3,i), i = 1, inpcrd_natom)

    ! We may or may not find velocities.  Thus we may inadvertently read
    ! box/angle info here if the velocities are missing. We clear the first
    ! 6 entries in atm_vel(:,:) to be able to pick up box/angle info if this
    ! happens.

    atm_vel(1:3, 1:2) = 0.d0
    read(inpcrd, 9028, end=700, err=1011) (atm_vel(1:3,i), i = 1, inpcrd_natom)
    velocities_found = .true.

700 continue
    if (.not. velocities_found) then 
      if (atm_vel(1,1) .ne. 0.d0) then
        box_found = .true.
        inpcrd_box(1:3) = atm_vel(1:3, 1)
        inpcrd_alpha = atm_vel(1, 2)
        inpcrd_beta = atm_vel(2, 2)
        inpcrd_gamma = atm_vel(3, 2)
      end if
      ! If no processing occurs:
      ! You have hit EOF without finding velocities or anything else.
      ! This is okay for nonperiodic simulations that also don't need
      ! velocities.
    else
      read(inpcrd, 9038, advance='no', end=900, err=1020) inpcrd_box(:)
      box_found = .true.
      ! If this read eof's, the angle values will remain 0.d0.
      read(inpcrd, 9038, end=900, err=1020) inpcrd_alpha, inpcrd_beta, &
                                            inpcrd_gamma
    end if

  ! ---------- Unformatted (Binary) Restart ----------      
  else  
    call amopen(inpcrd, inpcrd_name, 'O', 'U', 'R')

    read(inpcrd) title
    read(inpcrd) inpcrd_natom, tt
    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)
    ! We should always find coordinates:
    read(inpcrd, end=1000, err=1001) (atm_crd(1:3,i), i = 1, inpcrd_natom)
    ! We may or may not find velocities.  Thus we may inadvertently read
    ! box/angle info here if the velocities are missing. We clear the first
    ! 6 entries in atm_vel(:,:) to be able to pick up box/angle info if this
    ! happens.
    atm_vel(1:3, 1:2) = 0.d0
    read(inpcrd, end = 800, err = 800) (atm_vel(1:3,i), i = 1, inpcrd_natom)
    velocities_found = .true.

800 continue

    if (.not. velocities_found) then 
      if (atm_vel(1,1) .ne. 0.d0) then
        box_found = .true.
        inpcrd_box(1:3) = atm_vel(1:3, 1)
        inpcrd_alpha = atm_vel(1, 2)
        inpcrd_beta = atm_vel(2, 2)
        inpcrd_gamma = atm_vel(3, 2)
      end if
      ! If no processing occurs:
      ! You have hit EOF without finding velocities or anything else.
      ! This is okay for nonperiodic simulations that also don't need
      ! velocities.
    else
      ! If this read eof's, the angle values will remain 0.d0.
      read(inpcrd, end=900, err=900) inpcrd_box(1:3), &
                                     inpcrd_alpha, inpcrd_beta, inpcrd_gamma
      box_found = .true.
    end if

  endif

900 continue ! We branch here on an "okay" EOF

  if (velocities_needed) then
    if (.not. velocities_found) goto 1010       ! Velocity input error!
  else
    atm_vel(:,:) = 0.d0                         ! Velocities not used.
  end if

! Determine whether you got 0, 1 (beta), or 3 angles:

  if (inpcrd_alpha .ne. 0.d0) then
    if (inpcrd_beta .eq. 0.d0 .or. inpcrd_gamma .eq. 0.d0) then
      inpcrd_beta = inpcrd_alpha
      inpcrd_alpha = 90.d0
      inpcrd_gamma = 90.d0
    end if
    angles_found = .true.
  end if

  ! Only close inpcrd if not reading NetCDF
  if (.not. netcdf_input) close(inpcrd)

  if (.not. box_found .and. ntb .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'Box parameters not found in inpcrd file!'
    call mexit(6, 1)
  end if

  if (.not. angles_found) then
    inpcrd_alpha = 90.d0
    inpcrd_beta =  90.d0
    inpcrd_gamma = 90.d0
  end if

  return

666 continue

  write(mdout, '(3a)') error_hdr, 'Could not read second line of ', &
                       inpcrd_name
  write(mdout, '(2a)') extra_line_hdr, 'Bad inpcrd file!'

  call mexit(6, 1)

1000 continue

  write(mdout, '(a,a,a)') error_hdr, 'I could not find enough coordinates in ',&
                          inpcrd_name
  call mexit(6, 1)

1010 continue

  write(mdout, '(a,a,a)') error_hdr, 'I could not find enough velocities in ', &
                          inpcrd_name
  call mexit(6, 1)

1001 continue
  ! We hit here if we had an error reading floating points into the crd array

  ! Calculate which line the error occurred on. The first 2 lines are title,
  ! and info (NATOM, etc.). i-1 is the last successfully-read float, and there
  ! are 6 floats on each line, starting on line 3
  line_num = 3 + (i-1) / 2
  go to 1111

1011 continue
   ! We hit here if we had an error reading floating points into the vel array

   ! Calculate which line the error occurred on. Same as above for coordinates,
   ! but adjusted for the fact that all coordinates are present. There are 6 #s
   ! (2 atoms) on each line, with an odd atom getting its own line at the end
   line_num = inpcrd_natom / 2 + mod(inpcrd_natom, 2) + 3 + (i-1) / 2
   go to 1111

1111 continue

  write(mdout, '(2a,i5)') error_hdr, 'I could not understand line ', line_num

  ! Rewind the inpcrd file and read up to that line
  rewind(inpcrd)
  do i = 1, line_num
    read(inpcrd, '(a80)') read_buf
  end do
  
  write(mdout, '(a)') read_buf
  write(mdout, '()')

  call check_inpcrd_overflow(read_buf, ntb .gt. 0)

  call mexit(6, 1)

1020 continue

  write(mdout, '(a,a,a)') error_hdr, 'Could not read box lengths from ', &
                          inpcrd_name
  call mexit(6, 1)

9008 format(a80)
9018 format(i5, 5e15.7)
9019 format(i6, 5e15.7)
9020 format(i7, 5e15.7)
9021 format(i8, 5e15.7)
9028 format(6f12.7)
9038 format(3f12.7)

end subroutine init_inpcrd_dat

!*******************************************************************************
!
! Subroutine:  alloc_inpcrd_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)

  use mdin_ctrl_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     
  integer, intent(in)           :: inpcrd_natom

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_crd(3, inpcrd_natom), &
           atm_frc(3, inpcrd_natom), &
           atm_vel(3, inpcrd_natom), &
           atm_last_vel(3, inpcrd_natom), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_crd) + &
                          size(atm_frc) + &
                          size(atm_vel) + &
                          size(atm_last_vel)

  atm_last_vel(:,:) = 0.d0

#ifdef MPI
  atm_frc(:,:) = 0.d0           ! In case of nmr access to unowned atoms
#endif

  return

end subroutine alloc_inpcrd_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_inpcrd_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_inpcrd_dat(inpcrd_natom)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: inpcrd_natom

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)
  end if

  call mpi_bcast(atm_crd, 3 * inpcrd_natom, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  ! No need to broadcast atm_frc.

  if (ntx .eq. 1 .or. ntx .eq. 2) then
    atm_vel(:,:) = 0.d0 ! Will be initialized to random values later.
  else
    call mpi_bcast(atm_vel, 3 * inpcrd_natom, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  ! Make sure we broadcast temp0 if we read it from the inpcrd file.

  if (remd_method .gt. 1 .and. irest .eq. 1) &
    call mpi_bcast(temp0, 1, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)

  ! No need to broadcast atm_last_vel.

  return

end subroutine bcast_inpcrd_dat
#endif /* MPI */
end module inpcrd_dat_mod
