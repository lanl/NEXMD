!*******************************************************************************
!
! Module: inpcrd_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module inpcrd_dat_mod

  use file_io_dat_mod

  implicit none

  ! atm_crd = the atom xyz coordinates as a 2d array.
  ! atm_frc = the atom xyz forces as a 2d array. Not read here, but
  !           kept here as part of the fundamental data.
  ! atm_vel = the atom velocity xyz vectors as a 2d array.

  double precision, allocatable, save   :: atm_crd(:,:)
  double precision, allocatable, save   :: atm_frc(:,:)
  double precision, allocatable, save   :: atm_vel(:,:)
  double precision, allocatable, save   :: atm_last_vel(:,:)

#ifdef AMOEBA
  ! Used for amoeba with Beeman integrator:

  double precision, allocatable, save   :: acceleration(:,:)
  double precision, allocatable, save   :: old_acceleration(:,:)

  integer, save                         :: inpcrd_type = crdfile_type_undefined
#endif /* AMOEBA */

! Hide internal routines:

  private       alloc_inpcrd_mem

contains

!*******************************************************************************
!
! Subroutine:  init_old_type_inpcrd_dat
!
! Description: Read "old" style inpcrd, without section tags.  These files are
!              pre-amber 9 compliant, and compliant for amber 9 as long as the
!              Beeman integrator is not being used with the Amoeba forcefield.
!              
!*******************************************************************************

subroutine init_old_type_inpcrd_dat(num_ints, num_reals, inpcrd_natom, &
                                    inpcrd_alpha, inpcrd_beta, inpcrd_gamma, &
                                    inpcrd_box, tt, title)

  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
#ifdef AMOEBA
  use mdin_amoeba_dat_mod
#endif /* AMOEBA */
  use parallel_dat_mod

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

  character(1)          :: fform
  logical               :: formatted_input
  integer               :: i
  character             :: read_buf(80)            ! For format checking
  integer               :: inpcrd_version          ! For pmemd 2.01 format
  integer               :: zero_natom              ! For pmemd 2.01 format
  logical               :: box_found
  logical               :: angles_found
  logical               :: velocities_found
  logical               :: velocities_needed

  angles_found = .false.
  box_found = .false.
  velocities_found = .false.
  velocities_needed = (ntx .ge. 4)

#ifdef AMOEBA
! If using Amoeba with the Beeman integrator, an old style inpcrd cannot
! possibly provide needed accelerations:

  if (iamoeba .eq. 1) then
    if (beeman_integrator .eq. 1 .and. ntx .ge. 4) then
      write(mdout, '(a,a)') error_hdr, &
        'Old style inpcrd does not contain accelerations needed'
      write(mdout, '(a,a)') error_hdr, &
        'when restarting an Amoeba run with the Beeman integrator.'
      call mexit(mdout, 1)
    end if
  end if
#endif /* AMOEBA */

! Open the inpcrd file:

  formatted_input =  (ntx .eq. 1 .or. ntx .eq. 5 .or. ntx .eq. 7)

  inpcrd_box(:) = 0.d0
  inpcrd_alpha = 0.d0
  inpcrd_beta = 0.d0
  inpcrd_gamma = 0.d0

  if (formatted_input) then
    fform = 'F'
  else
    fform = 'U'
  end if

  call amopen(inpcrd, inpcrd_name, 'O', fform, 'R')

  if (formatted_input) then

    ! Sander 7/8 uses a mechanism of checking if the 6th char is blank to see if
    ! this is an old i5 format sander 6 file or a new i6 format sander 7 file.
    ! We also check for the versioned pmemd 2.01 format, with an initial natom
    ! of 0.  This was really the best approach, but unfortunately sander 7/8
    ! did not implement something like it.

    read(inpcrd, 9008) title
    read(inpcrd, '(80a1)') read_buf
    rewind(inpcrd)
    read(inpcrd, 9008) title

    if (read_buf(6) .eq. ' ') then ! Sander 6 or pmemd 2.01 format...

       read(inpcrd, '(i5)') inpcrd_natom
       rewind(inpcrd)
       read(inpcrd, 9008) title

       if (inpcrd_natom .ne. 0) then
         read(inpcrd, 9018) inpcrd_natom, tt
       else
         read(inpcrd, '(2i5, i10, 4e15.7)') &
              zero_natom, inpcrd_version, inpcrd_natom, tt
         if (inpcrd_version .ne. 2) then
           write(mdout, '(a,a)') error_hdr, &
             'Unrecognized format for inpcrd/restrt file!'
           call mexit(6, 1)
         end if
       end if

    else        ! It is probably sander 7/8 large system format...

      read(inpcrd, 9019) inpcrd_natom, tt

    end if

    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)

    ! We should always find coordinates:

    read(inpcrd, 9028, end=1000, err=1000) (atm_crd(1:3,i), i = 1, inpcrd_natom)

    ! We may or may not find velocities.  Thus we may inadvertently read
    ! box/angle info here if the velocities are missing. We clear the first
    ! 6 entries in atm_vel(:,:) to be able to pick up box/angle info if this
    ! happens.

    atm_vel(1:3, 1:2) = 0.d0
    read(inpcrd, 9028, end=700, err=1010) (atm_vel(1:3,i), i = 1, inpcrd_natom)
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
      
  else  ! Binary input...

    read(inpcrd) title
    read(inpcrd) inpcrd_natom, tt
    call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)
    ! We should always find coordinates:
    read(inpcrd, end=1000, err=1000) (atm_crd(1:3,i), i = 1, inpcrd_natom)
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

  close(inpcrd)

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

1000 continue

  write(mdout, '(a,a,a)') error_hdr, 'Could not read coords from ', &
                          inpcrd_name
  call mexit(6, 1)

1010 continue

  write(mdout, '(a,a,a)') error_hdr, 'Could not read velocities from ', &
                          inpcrd_name
  call mexit(6, 1)

1020 continue

  write(mdout, '(a,a,a)') error_hdr, 'Could not read box lengths from ', &
                          inpcrd_name
  call mexit(6, 1)

9008 format(a80)
9018 format(i5, 5e15.7)
9019 format(i6, 5e15.7)
9028 format(6f12.7)
9038 format(3f12.7)

end subroutine init_old_type_inpcrd_dat

#ifdef AMOEBA
!*******************************************************************************
!
! Subroutine:  init_new_type_inpcrd_dat
!
! Description:  This format is currently used for the Amoeba forcefield when the
!               Beeman integrator is also used, for Amber 9.  The scope of use
!               may be expanded in the future.
!              
!*******************************************************************************

subroutine init_new_type_inpcrd_dat(num_ints, num_reals, inpcrd_natom, &
                                    inpcrd_alpha, inpcrd_beta, inpcrd_gamma, &
                                    inpcrd_box, tt, title)

  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use mdin_amoeba_dat_mod
  use parallel_dat_mod
  use nextinpcrd_section_mod

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

  integer                       :: i, j
  integer                       :: array_cnt
  integer                       :: iok, ionerr
  logical                       :: velocities_needed
  logical                       :: accelerations_needed
  character(len = 80)           :: fmt, fmtin, dtype
  double precision, parameter   :: sander_to_tinker_time_convert = 20.455d0

  inpcrd_box(:) = 0.d0
  inpcrd_alpha = 0.d0
  inpcrd_beta = 0.d0
  inpcrd_gamma = 0.d0

  call amopen(inpcrd, inpcrd_name, 'O', 'F', 'R')

  call nxtsec_crd_reset()

  dtype = 'TITLE'
  fmtin = '(a)'
  ionerr = 0 ! fatal if missing

  call nxtsec_crd(inpcrd,  mdout,  ionerr, fmtin,  dtype,  fmt,  iok)

  read(inpcrd, fmt) title

  dtype = 'ATOMIC_COORDS_SIMULATION_TIME'
  fmtin = '(E16.8)'
  ionerr = 1 ! not fatal if missing ---but should be there

  call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

  if (iok .eq. 0) then
    read(inpcrd, fmt) tt
  else
    tt = 0.d0
  end if

  dtype = 'ATOMIC_COORDS_NUM_LIST'
  fmtin = '(I8)'
  ionerr = 0 ! fatal if missing

  call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

  read(inpcrd, fmt) inpcrd_natom

  ! We don't yet know natom, and will do consistency checks afterward.

  call alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)

  dtype = 'ATOMIC_COORDS_LIST'
  fmtin = '(5E16.8)'
  ionerr = 0 ! fatal if missing

  call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

  read(inpcrd, fmt)((atm_crd(i, j), i = 1, 3), j = 1, inpcrd_natom)

  if (ntx .ge. 4) then ! look for velocities, etc.

    dtype = 'ATOMIC_VELOCITIES_NUM_LIST'
    fmtin = '(I8)'
    ionerr = 0 ! fatal if missing

    call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

    read(inpcrd, fmt) array_cnt

    if (array_cnt .ne. inpcrd_natom) then
      write(mdout, '(a,a)') error_hdr, &
        'Atomic velocity count in amoeba inpcrd is inconsistent!'
      call mexit(mdout, 1)
    end if

    dtype = 'ATOMIC_VELOCITIES_LIST'
    fmtin = '(5E16.8)'
    ionerr = 0 ! fatal if missing

    call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

    read(inpcrd, fmt)((atm_vel(i, j), i = 1, 3), j = 1, array_cnt)

    atm_vel(:,:) = sander_to_tinker_time_convert * atm_vel(:,:)

  end if

  ! If using amoeba and the Beeman integrator, accelerations are needed if
  ! ntx >= 4. Otherwise we are setting random initial velocities and using
  ! 0 accels/old_accels.

  if (iamoeba .eq. 1 .and. ntx .ge. 4) then
    if (beeman_integrator .eq. 1) then

      dtype = 'ATOMIC_ACCELERATIONS_NUM_LIST'
      fmtin = '(I8)'
      ionerr = 0 ! fatal if missing

      call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)


      read(inpcrd, fmt) array_cnt

      if (array_cnt .ne. inpcrd_natom) then
        write(mdout, '(a,a)') error_hdr, &
          'Atomic acceleration count in amoeba inpcrd is inconsistent!'
        call mexit(mdout, 1)
      end if

      dtype = 'ATOMIC_ACCELERATIONS_LIST'
      fmtin = '(5E16.8)'
      ionerr = 0

      call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

      read(inpcrd, fmt)((acceleration(i, j), i = 1, 3), j = 1, array_cnt)

      acceleration(:,:) = acceleration(:,:) * sander_to_tinker_time_convert**2 

      dtype = 'OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
      fmtin = '(I8)'
      ionerr = 0 ! fatal if missing

      call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

      read(inpcrd, fmt) array_cnt

      if (array_cnt .ne. inpcrd_natom) then
        write(mdout, '(a,a)') error_hdr, &
          'Old Atomic acceleration count in amoeba inpcrd is inconsistent!'
        call mexit(mdout, 1)
      end if

      dtype = 'OLD_ATOMIC_ACCELERATIONS_LIST'
      fmtin = '(5E16.8)'
      ionerr = 0 ! fatal if missing

      call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

      read(inpcrd, fmt)((old_acceleration(i, j), i = 1, 3), j = 1, array_cnt)

      old_acceleration(:,:) = old_acceleration(:,:) * &
                              sander_to_tinker_time_convert**2 

    end if
  end if

  dtype = 'UNIT_CELL_PARAMETERS'
  fmtin = '(3e16.8)'
  
  if (ntb .ne. 0) then
    ionerr = 0 ! fatal if missing
  else
    ionerr = 1 ! not fatal if missing
  end if

  call nxtsec_crd(inpcrd, mdout, ionerr, fmtin, dtype, fmt, iok)

  if (iok .eq. 0) then
    read(inpcrd, fmt) inpcrd_box(1:3), inpcrd_alpha, inpcrd_beta, inpcrd_gamma
  else
    inpcrd_box(1:3) = 0.d0
    inpcrd_alpha = 0.d0
    inpcrd_beta = 0.d0
    inpcrd_gamma = 0.d0
  end if

  call nxtsec_crd_reset()

  close(inpcrd)

  return

end subroutine init_new_type_inpcrd_dat
#endif /* AMOEBA */

!*******************************************************************************
!
! Subroutine:  alloc_inpcrd_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_inpcrd_mem(num_ints, num_reals, inpcrd_natom)

  use mdin_ctrl_dat_mod
#ifdef AMOEBA
  use mdin_amoeba_dat_mod
#endif /* AMOEBA */
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

#ifdef AMOEBA
  if (iamoeba .eq. 1) then
    if (beeman_integrator .eq. 1) then

      allocate(acceleration(3, inpcrd_natom), &
               old_acceleration(3, inpcrd_natom), &
               stat = alloc_failed)

      if (alloc_failed .ne. 0) call setup_alloc_error
    
      num_reals = num_reals + size(acceleration) + &
                              size(old_acceleration)

      ! Accelerations are zero'd in case this is an initial run with
      ! no accelerations available; the old_crd() zeroing is just for
      ! insurance...

      acceleration(:,:) = 0.d0
      old_acceleration(:,:) = 0.d0

    end if
  end if
#endif /* AMOEBA */

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
#ifdef AMOEBA
  use mdin_amoeba_dat_mod
#endif /* AMOEBA */
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
                 mpi_comm_world, err_code_mpi)

  ! No need to broadcast atm_frc.

  if (ntx .eq. 1 .or. ntx .eq. 2) then
    atm_vel(:,:) = 0.d0 ! Will be initialized to random values later.
  else
    call mpi_bcast(atm_vel, 3 * inpcrd_natom, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  ! No need to broadcast atm_last_vel.

#ifdef AMOEBA
  ! For amoeba, Beeman integrator, restart run we need the accelerations.

  if (iamoeba .eq. 1 .and. ntx .ge. 4) then
    if (beeman_integrator .eq. 1) then
      call mpi_bcast(acceleration, 3 * inpcrd_natom, mpi_double_precision, &
                     0, mpi_comm_world, err_code_mpi)
      call mpi_bcast(old_acceleration, 3 * inpcrd_natom, mpi_double_precision, &
                     0, mpi_comm_world, err_code_mpi)
    end if
  end if
#endif /* AMOEBA */

  return

end subroutine bcast_inpcrd_dat
#endif /* MPI */
end module inpcrd_dat_mod
