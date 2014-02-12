#include "copyright.i"

!*******************************************************************************!
! Module:  bintraj_mod
!
! Description:  Module for generating binary trajectory-type output in NetCDF
!               format.  Originally developed by John Mongan, November 2005
!*******************************************************************************                                                                                
module bintraj_mod

  implicit none

  private

  ! Variables used for mdcrd:

  integer, save         :: mdcrd_ncid
  integer, save         :: coord_var_id
  integer, save         :: cell_length_var_id
  integer, save         :: cell_angle_var_id
  integer, save         :: mdcrd_time_var_id
  integer, save         :: mdcrd_veloc_var_id
  integer, save         :: mdcrd_frame

  ! Variables used for mdvel:

  integer, save         :: mdvel_ncid
  integer, save         :: mdvel_time_var_id
  integer, save         :: mdvel_veloc_var_id
  integer, save         :: mdvel_frame

  ! In this context, count of atoms dumped to file:
  
  integer, save         :: atom_wrt_cnt

  public        open_binary_files, &
                close_binary_files, &
                write_binary_crds, &
                write_binary_vels, &
                write_binary_cell_dat, &
                end_binary_frame
contains

!*******************************************************************************!
! Subroutine:  open_binary_files
!
! Description:  Open the coordinate, velocity, and energy files.
!
!*******************************************************************************

subroutine open_binary_files

#ifdef BINTRAJ
  use netcdf
#endif
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

#ifdef BINTRAJ

! Local variables:

  integer       :: cmode
  integer       :: frame_dim_id
  integer       :: spatial_dim_id
  integer       :: cell_spatial_dim_id
  integer       :: cell_angular_dim_id
  integer       :: atom_dim_id
  integer       :: label_dim_id

  integer       :: old_mode
  integer       :: spatial_var_id
  integer       :: cell_spatial_var_id
  integer       :: cell_angular_var_id
  integer       :: status

  integer, parameter    :: nc_label_len = 5

  ! Create netCDF file

  cmode = nf90_64bit_offset

  if (owrite .eq. 'N') cmode = ior(cmode, nf90_noclobber)

  if (ntwprt .gt. 0) then
    atom_wrt_cnt = ntwprt
  else
    atom_wrt_cnt = natom
  end if

  if (ntwx .gt. 0) then

    status = nf90_create(path=mdcrd_name, cmode=cmode, ncid=mdcrd_ncid)

    call checkerror(status, 'create file')

    if (status .ne. nf90_noerr) then
#ifndef DUMB
      write (0,*) 'Error on opening ', mdcrd_name
#endif
      call mexit(mdout,1)
    end if

    ! Define dimensions

    call checkerror(nf90_def_dim(mdcrd_ncid, 'frame', nf90_unlimited, &
                    frame_dim_id))

    call checkerror(nf90_def_dim(mdcrd_ncid, 'spatial', 3, spatial_dim_id))

    call checkerror(nf90_def_dim(mdcrd_ncid, 'atom', atom_wrt_cnt, atom_dim_id))

    call checkerror(nf90_def_dim(mdcrd_ncid, 'label', nc_label_len, &
                    label_dim_id))

    ! Set global attributes

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'title', &
                    prmtop_ititl), 'define title')

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'application', &
                    'AMBER'), 'define application')

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'program', &
                    'pmemd'), 'define program')

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'programVersion', &
                    '9.0'), 'define programVersion')

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'Conventions', &
                    'AMBER'), 'define Convention')

    call checkerror(nf90_put_att(mdcrd_ncid, nf90_global, 'ConventionVersion', &
                    '1.0'), 'define ConventionVersion')

    ! Define non-optional variables

    call checkerror(nf90_def_var(mdcrd_ncid, 'spatial', nf90_char, &
                    (/ spatial_dim_id /), spatial_var_id))
   
    call checkerror(nf90_def_var(mdcrd_ncid, 'time', nf90_float, &
                    (/ frame_dim_id /), mdcrd_time_var_id))

    call checkerror(nf90_put_att(mdcrd_ncid, mdcrd_time_var_id, 'units', &
                    'picosecond'), 'define time units')

    call checkerror(nf90_def_var(mdcrd_ncid, 'coordinates', nf90_float, &
                    (/ spatial_dim_id, atom_dim_id, frame_dim_id /), &
                    coord_var_id), 'define coordinates')

    call checkerror(nf90_put_att(mdcrd_ncid, coord_var_id, 'units', &
                    'angstrom'), 'define coordinates units')

    ! Define optional variables

    if (ntb .gt. 0 ) then

      ! Dimensions:

      call checkerror(nf90_def_dim(mdcrd_ncid, 'cell_spatial', 3, &
                      cell_spatial_dim_id), 'define cell spatial dim')

      call checkerror(nf90_def_dim(mdcrd_ncid, 'cell_angular', 3, &
                      cell_angular_dim_id), 'define cell angular dim')

      ! Label variables:

      call checkerror(nf90_def_var(mdcrd_ncid, 'cell_spatial', nf90_char, &
                      (/ cell_spatial_dim_id /), cell_spatial_var_id), &
                      'define cell spatial var')

      call checkerror(nf90_def_var(mdcrd_ncid, 'cell_angular', nf90_char, &
                      (/ label_dim_id, cell_angular_dim_id /), &
                      cell_angular_var_id), "define cell angular var")

      ! Data variables and attributes:

      call checkerror(nf90_def_var(mdcrd_ncid, 'cell_lengths', nf90_double, &
                      (/ cell_spatial_dim_id, frame_dim_id /), &
                      cell_length_var_id), 'define cell lengths')

      call checkerror(nf90_put_att(mdcrd_ncid, cell_length_var_id, 'units', &
                      'angstrom'), 'define cell length units')
  
      call checkerror(nf90_def_var(mdcrd_ncid, 'cell_angles', nf90_double, &
                      (/ cell_angular_dim_id, frame_dim_id /), &
                      cell_angle_var_id), 'define cell angles')

      call checkerror(nf90_put_att(mdcrd_ncid, cell_angle_var_id, 'units', &
                      'degree'), 'define cell angle units')

    end if

    if (ntwv .lt. 0 ) then

      call checkerror(nf90_def_var(mdcrd_ncid, 'velocities', nf90_float, &
                      (/ spatial_dim_id, atom_dim_id, frame_dim_id /), &
                      mdcrd_veloc_var_id), 'define velocities')

      call checkerror(nf90_put_att(mdcrd_ncid, mdcrd_veloc_var_id, 'units', &
                      'angstrom/picosecond'), 'define velocity units')

      call checkerror(nf90_put_att(mdcrd_ncid, mdcrd_veloc_var_id, &
                      'scale_factor', 20.455), 'define velocity scale_factor')

    end if

    ! Set NoFill and end definition mode

    call checkerror(nf90_set_fill(mdcrd_ncid, nf90_nofill, old_mode), &
                    'set no fill')

    call checkerror(nf90_enddef(mdcrd_ncid), 'end define')

    ! Fill dimension label variables

    call checkerror(nf90_put_var(mdcrd_ncid, spatial_var_id, &
                    (/ 'x', 'y', 'z' /), start = (/ 1 /), count = (/ 3 /)), &
                    'write spatial variable')

    if (ntb .gt. 0) then

      call checkerror(nf90_put_var(mdcrd_ncid, cell_spatial_var_id, &
                      (/ 'a','b','c' /), start=(/ 1 /), count=(/ 3 /)), &
                      "write spatial variable")

      call checkerror(nf90_put_var(mdcrd_ncid, cell_angular_var_id, &
                      (/ 'alpha','beta ','gamma' /), &
                      start=(/ 1, 1 /), count=(/ nc_label_len, 3 /)), &
                      "write spatial variable")

    end if

    mdcrd_frame = 1

  end if

  if (ntwv .gt. 0) then

    status = nf90_create(path=mdvel_name, cmode=cmode, ncid=mdvel_ncid)

    call checkerror(status, 'create file')

    if (status .ne. nf90_noerr) then
#ifndef DUMB
      write (0,*) 'Error on opening ', mdvel_name
#endif
      call mexit(mdout,1)
    end if

    ! Define dimensions

    call checkerror(nf90_def_dim(mdvel_ncid, 'frame', nf90_unlimited, &
                    frame_dim_id))

    call checkerror(nf90_def_dim(mdvel_ncid, 'spatial', 3, spatial_dim_id))

    call checkerror(nf90_def_dim(mdvel_ncid, 'atom', atom_wrt_cnt, atom_dim_id))

    call checkerror(nf90_def_dim(mdvel_ncid, 'label', nc_label_len, &
                    label_dim_id))

    ! Set global attributes

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'title', &
                    prmtop_ititl), 'define title')

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'application', &
                    'AMBER'), 'define application')

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'program', &
                    'pmemd'), 'define program')

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'programVersion', &
                    '9.0'), 'define programVersion')

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'Conventions', &
                    'AMBER'), 'define Convention')

    call checkerror(nf90_put_att(mdvel_ncid, nf90_global, 'ConventionVersion', &
                    '1.0'), 'define ConventionVersion')

    ! Define non-optional variables

    call checkerror(nf90_def_var(mdvel_ncid, 'spatial', nf90_char, &
                    (/ spatial_dim_id /), spatial_var_id))
   
    call checkerror(nf90_def_var(mdvel_ncid, 'time', nf90_float, &
                    (/ frame_dim_id /), mdvel_time_var_id))

    call checkerror(nf90_put_att(mdvel_ncid, mdvel_time_var_id, 'units', &
                    'picosecond'), 'define time units')

    call checkerror(nf90_def_var(mdvel_ncid, 'velocities', nf90_float, &
                    (/ spatial_dim_id, atom_dim_id, frame_dim_id /), &
                    mdvel_veloc_var_id), 'define velocities')

    call checkerror(nf90_put_att(mdvel_ncid, mdvel_veloc_var_id, 'units', &
                    'angstrom/picosecond'), 'define velocity units')

    call checkerror(nf90_put_att(mdvel_ncid, mdvel_veloc_var_id, &
                    'scale_factor', 20.455), 'define velocity scale_factor')

    ! Set NoFill and end definition mode

    call checkerror(nf90_set_fill(mdvel_ncid, nf90_nofill, old_mode), &
                    'set no fill')

    call checkerror(nf90_enddef(mdvel_ncid), 'end define')

    ! Fill dimension label variables

    call checkerror(nf90_put_var(mdvel_ncid, spatial_var_id, &
                    (/ 'x', 'y', 'z' /), start = (/ 1 /), count = (/ 3 /)), &
                    'write spatial variable')

    mdvel_frame = 1

  end if


#else

  call no_bintraj_err

#endif

  return
   
end subroutine open_binary_files


!*******************************************************************************!
! Subroutine:  close_binary_files
!
! Description:  Close the coordinate, velocity, and energy files.
!
!*******************************************************************************

subroutine close_binary_files

  use mdin_ctrl_dat_mod
#ifdef BINTRAJ
  use netcdf
#endif

  implicit none

#ifdef BINTRAJ

  if (ntwx .gt. 0) call checkerror(nf90_close(mdcrd_ncid))
  if (ntwv .gt. 0) call checkerror(nf90_close(mdvel_ncid))

#else

  call no_bintraj_err

#endif

  return

end subroutine close_binary_files

!*******************************************************************************!
! Subroutine:   write_binary_crds
!
! Description:  Emit coordinates to trajectory file mdcrd
!
!*******************************************************************************

subroutine write_binary_crds

#ifdef BINTRAJ
  use netcdf
#endif
  use axis_optimize_mod
  use inpcrd_dat_mod

  implicit none

#ifdef BINTRAJ

! Local variables:

  integer       :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then

    call checkerror(nf90_put_var(mdcrd_ncid, coord_var_id, atm_crd(:,:), &
                    start=(/ 1, 1, mdcrd_frame /), &
                    count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write atom coords')
  else

    call write_binary_crds_axis_flipped(ord1, ord2, ord3)

  end if
#else

  call no_bintraj_err

#endif

  return

end subroutine write_binary_crds

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   write_binary_crds_axis_flipped
!
! Description:  Emit coordinates to trajectory file mdcrd
!
!*******************************************************************************

subroutine write_binary_crds_axis_flipped(ord1, ord2, ord3)

  use netcdf
  use inpcrd_dat_mod

  implicit none


! Formal arguments:

  integer       :: ord1, ord2, ord3

! Local variables:

  integer       :: i
  real          :: buf(3, atom_wrt_cnt)

  do i = 1, atom_wrt_cnt
    buf(1, i) = atm_crd(ord1, i)
    buf(2, i) = atm_crd(ord2, i)
    buf(3, i) = atm_crd(ord3, i)
  end do

  call checkerror(nf90_put_var(mdcrd_ncid, coord_var_id, buf(:,:), &
                  start=(/ 1, 1, mdcrd_frame /), &
                  count=(/ 3, atom_wrt_cnt, 1 /)), &
                  'write atom coords')

  return

end subroutine write_binary_crds_axis_flipped
#endif /* BINTRAJ */

!*******************************************************************************!
! Subroutine:   write_binary_vels
!
! Description:  Emit velocities to trajectory or velocities file (mdrcrd or
!               mdvel, dependent on value of ntwv).
!
!*******************************************************************************

subroutine write_binary_vels

#ifdef BINTRAJ
  use netcdf
#endif
  use axis_optimize_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod

  implicit none

#ifdef BINTRAJ

! Local variables:

  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3

  if (ntwv .lt. 0) then
    ncid = mdcrd_ncid
    frame = mdcrd_frame
    var_id = mdcrd_veloc_var_id
  else
    ncid = mdvel_ncid
    frame = mdvel_frame
    var_id = mdvel_veloc_var_id
  end if

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then

    call checkerror(nf90_put_var(ncid, var_id, atm_vel(:,:), &
                    start=(/ 1, 1, frame /), &
                    count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write velocities')
  else

    call write_binary_vels_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3)

  end if

#else

  call no_bintraj_err

#endif

  return

end subroutine write_binary_vels

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   write_binary_vels_axis_flipped
!
! Description:  Emit velocities to trajectory or velocities file (mdrcrd or
!               mdvel, dependent on value of ntwv).
!
!*******************************************************************************

subroutine write_binary_vels_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3)

  use netcdf
  use inpcrd_dat_mod

  implicit none

! Formal arguments:

  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3

! Local variables:

  integer       :: i
  real          :: buf(3, atom_wrt_cnt)

  do i = 1, atom_wrt_cnt
    buf(1, i) = atm_vel(ord1, i)
    buf(2, i) = atm_vel(ord2, i)
    buf(3, i) = atm_vel(ord3, i)
  end do

  call checkerror(nf90_put_var(ncid, var_id, buf(:,:), &
                    start=(/ 1, 1, frame /), &
                    count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write velocities')
  return

end subroutine write_binary_vels_axis_flipped
#endif /* BINTRAJ */

!*******************************************************************************!
! Subroutine:   write_binary_cell_dat
!
! Description:  Emit cell length and angle data to mdcrd.
!
!*******************************************************************************

subroutine write_binary_cell_dat

#ifdef BINTRAJ
  use netcdf
#endif
  use axis_optimize_mod
  use file_io_dat_mod
  use pbc_mod

  implicit none

#ifdef BINTRAJ

! Local variables:

  integer       :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  call checkerror(nf90_put_var(mdcrd_ncid, cell_length_var_id, &
                  (/ pbc_box(ord1), pbc_box(ord2), pbc_box(ord3) /), &
                  start=(/ 1, mdcrd_frame /), count=(/ 3, 1 /)), &
                  'write cell lengths')

  call checkerror(nf90_put_var(mdcrd_ncid, cell_angle_var_id, &
                  (/ pbc_alpha, pbc_beta, pbc_gamma /), &
                  start=(/ 1, mdcrd_frame /), count=(/ 3, 1 /)), &
                  'write cell angles')
#else

  call no_bintraj_err

#endif

  return

end subroutine write_binary_cell_dat

!*******************************************************************************!
! Subroutine:   end_binary_frame
!
! Description:  Write scalar data and increment frame counter.
!
!*******************************************************************************

subroutine end_binary_frame(unit)

#ifdef BINTRAJ
  use netcdf
#endif
  use file_io_dat_mod
  use mdin_ctrl_dat_mod

  implicit none

  integer, intent(in)           :: unit

#ifdef BINTRAJ

  if (unit .eq. mdcrd) then

    call checkerror(nf90_put_var(mdcrd_ncid, mdcrd_time_var_id, (/ t /), &
                    start=(/ mdcrd_frame /), count=(/ 1 /)), 'write time')
   
    call checkerror(nf90_sync(mdcrd_ncid))

    mdcrd_frame = mdcrd_frame + 1

  else if (unit .eq. mdvel) then

    call checkerror(nf90_put_var(mdvel_ncid, mdvel_time_var_id, (/ t /), &
                    start=(/ mdvel_frame /), count=(/ 1 /)), 'write time')
   
    call checkerror(nf90_sync(mdvel_ncid))

    mdvel_frame = mdvel_frame + 1

  else

    write (mdout, *) 'Error: unhandled unit ', unit, &
                 ' selected for flush in bintraj'

  end if

#else

  call no_bintraj_err

#endif

  return

end subroutine end_binary_frame

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   checkerror
!
! Description:  Checks error return from netCDF routines.
!
!*******************************************************************************

subroutine checkerror(status, location)

  use netcdf
  use file_io_dat_mod

  implicit none

  integer, intent(in)                   :: status       ! netCDF return code
  character(*), optional, intent(in)    :: location     ! purpose of call
  
  if (status .ne. nf90_noerr) then

    write(mdout, *) 'NetCDF error: ', trim(nf90_strerror(status))

    if (present(location)) then
      write(mdout, *) '  at ', location
    end if

  end if

  return

end subroutine checkerror

#else

!*******************************************************************************!
! Subroutine:   no_bintraj_err
!
! Description:  Checks error return from netCDF routines.
!
!*******************************************************************************

subroutine no_bintraj_err

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

  write(mdout, *) &
    'No binary trajectory support in this build of pmemd!'
  write(mdout, *) &
    'To write binary trajectories, reconfigure using the bintraj option.'
  call mexit(mdout, 1)


end subroutine no_bintraj_err

#endif

end module bintraj_mod
