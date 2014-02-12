#include "copyright.i"

!*******************************************************************************
! Module:  binrestart_mod
!
! Description:  Module for generating Amber NetCDF restart files. Originally 
!               developed by Dan Roe, August 2011. Based on the Amber NetCDF 
!               trajectory format originally developed by John Mongan.
!*******************************************************************************

module binrestart_mod
   private

   integer, save :: atomDID, coordVID, velocityVID 
   integer, save :: cellAngleVID, cellLengthVID, spatialDID, labelDID
   integer, save :: cell_spatialDID, cell_angularDID, spatialVID, timeVID
   integer, save :: cell_spatialVID, cell_angularVID, TempVID

   ! If WRITE_NC_RESTART is called for the main restart file, this variable
   ! indicates whether the file has been set up.
   logical, save :: mainRestartSetup = .false.

   public write_nc_restart, &
          read_nc_restart, &
          read_nc_restart_atoms, &
          read_nc_refc, &
          check_nc_restart
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write Netcdf restart file.
!-------------------------------------------------------------------
!     --- WRITE_NC_RESTART ---
!-------------------------------------------------------------------
!     Write Netcdf Restart file with given filename and title. 
!     owrite indicates overwrite status (N is no overwrite), natom 
!     is the # atoms, ntb>0 indicates presence of box coords, isMain
!     indicates this is the main restart file. The main restart file
!     can only be created and set-up once, all others will always
!     be created and set-up. 
!     Coords and Velo are the coordinates and velocities, temp0 is 
!     the current  temperature and Time is the current simulation 
!     time. If imin is 0, this is for MD so velocities will be 
!     written, otherwise no velocity info will be saved.
subroutine write_nc_restart(filename,title,owrite,natom,ntb,isMain,Coords,Velo,temp0,Time,imin,&
                            box,alpha,beta,gamma)
#ifdef BINTRAJ
   use netcdf
   use axis_optimize_mod
#  ifdef MPI
   use remd_mod, only: remd_method
#  endif
#endif
   use file_io_dat_mod, only: mdout
   use pmemd_lib_mod, only: mexit
   implicit none
   ! Formal arguments:
   character(len=*), intent(in) :: filename
   character(len=*), intent(in)  :: title
   character, intent(in) :: owrite
   integer, intent(in) :: natom,ntb
   logical, intent(in) :: isMain 
   double precision, dimension(:, :), intent(in) :: Coords, Velo
   double precision, intent(in) :: temp0, Time
   integer, intent(in) :: imin
   double precision, dimension(3), intent(in) :: box
   double precision, intent(in) :: alpha, beta, gamma
#ifdef BINTRAJ
   ! Local variables:
   integer :: ncid
   integer :: cmode, err, oldMode
   integer :: ord1, ord2, ord3

   ! ------------------------- File Setup --------------------------
   ! If first call, create the file and set up all dimensions and vars
   if ( (isMain .eqv. .false.) .or. (mainRestartSetup .eqv. .false.) ) then
      ! owrite status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
      ! pmemd flag -O='R', -A='U', default='N'
      cmode = nf90_64bit_offset
      if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
      !if (owrite == 'U' .and. facc == 'A') cmode = ior(cmode, nf90_noclobber)
      err = nf90_create(path=filename,cmode=cmode,ncid=ncid)
      if (err == nf90_eexist) then
         write(mdout,*) 'Error: write_nc_restart(): File exists and -O not specified: ', filename 
         call mexit(mdout,1)
      endif
      call checkerror(err, "write_nc_restart(): Creating netcdf restart")
      if ( err /= nf90_noerr ) call mexit(mdout,1)
      ! Time variable
      call checkerror( nf90_def_var(ncid, "time", nf90_double, timeVID),"define timeVID" )
      call checkerror( nf90_put_att(ncid, timeVID, "units", "picosecond"), &
                        "define timeVID units")
      ! Spatial dimension and variable
      call checkerror( nf90_def_dim(ncid, "spatial", 3, spatialDID),"define spatialDID")
      call checkerror( nf90_def_var(ncid, "spatial", nf90_char, &
                         (/ spatialDID /), spatialVID), "define spatialVID")
      ! Atom dimension
      call checkerror( nf90_def_dim(ncid, "atom", natom, atomDID),"define atomDID")
      ! Coord variable
      call checkerror( nf90_def_var(ncid, "coordinates", nf90_double, &
                         (/ spatialDID, atomDID /), coordVID ), "define coordVID")
      call checkerror( nf90_put_att(ncid, coordVID, "units", "angstrom"), &
                        "define coordVID units")
      ! Velocity variable
      if (imin .eq. 0) then
         call checkerror( nf90_def_var(ncid, "velocities", nf90_double, &
                            (/ spatialDID, atomDID /), velocityVID), "define velocityVID")
         call checkerror( nf90_put_att(ncid, velocityVID, "units", "angstrom/picosecond"), &
                           "define velocityVID units")
         call checkerror( nf90_put_att(ncid, velocityVID, "scale_factor", 20.455), &
                           "define velocityVID scale factor")
      endif
      ! Box info
      if (ntb>0) then
         ! Cell Spatial
         call checkerror( nf90_def_dim(ncid, "cell_spatial", 3, cell_spatialDID), &
                            "define cell_spatialDID")
         call checkerror( nf90_def_var(ncid, "cell_spatial", nf90_char, &
                            (/ cell_spatialDID /), cell_spatialVID), "define cell_spatialVID")
         ! Cell angular
         call checkerror( nf90_def_dim(ncid, "label", 5, labelDID), "define labelDID")
         call checkerror( nf90_def_dim(ncid, "cell_angular", 3, cell_angularDID), &
                            "define cell_angularDID")
         call checkerror( nf90_def_var(ncid, "cell_angular", nf90_char, &
                            (/ labelDID, cell_angularDID /), cell_angularVID), &
                           "define cell_angularVID")
         ! Cell length
         call checkerror( nf90_def_var(ncid, "cell_lengths", nf90_double, &
                            (/ cell_spatialDID /), cellLengthVID), "define cellLengthVID")
         call checkerror( nf90_put_att(ncid, cellLengthVID, "units", "angstrom"), &
                           "define cellLengthVID units")
         ! Cell angle
         call checkerror( nf90_def_var(ncid, "cell_angles", nf90_double, &
                            (/ cell_angularDID /), cellAngleVID), "define cellAngleVID")
         call checkerror( nf90_put_att(ncid, cellAngleVID, "units", "degree"), &
                           "define cellAngleVID units")
      endif
#  ifdef MPI
      ! Replica Temperature
      if (remd_method>0) then
         call checkerror( nf90_def_var(ncid, "temp0", nf90_double, TempVID), &
                           "define TempVID")
         call checkerror( nf90_put_att(ncid, TempVID, "units", "kelvin"), &
                           "define TempVID units")
      endif
#  endif
      ! Global attributes: Title etc
      call checkerror( nf90_put_att(ncid, nf90_global, &
                        'title', title), "define title")
      call checkerror(nf90_put_att(ncid,nf90_global, "application", &
                        'AMBER'), "define application")
      call checkerror(nf90_put_att(ncid,nf90_global, "program", &
                        'pmemd'), "define program")
      call checkerror(nf90_put_att(ncid,nf90_global, "programVersion", &
                        '11.0'), "define programVersion")
      call checkerror(nf90_put_att(ncid,nf90_global, "Conventions", &
                        'AMBERRESTART'), "define Convention")
      call checkerror(nf90_put_att(ncid,nf90_global, "ConventionVersion", &
                        '1.0'), "define ConventionVersion")
      ! Set fill mode
      call checkerror( nf90_set_fill(ncid, nf90_nofill, oldMode), "Setting fill mode")
      ! End definitions
      call checkerror( nf90_enddef(ncid), "end define")
      ! Specify dimension labels
      call checkerror( nf90_put_var(ncid, spatialVID, &
         (/ 'x','y','z' /), start = (/ 1 /), count = (/ 3 /)), "write spatial variable")
      if (ntb>0) then
         call checkerror(nf90_put_var(ncid, cell_spatialVID, &
               (/ 'a','b','c' /), start = (/ 1 /), count = (/ 3 /)), &
               "write spatial variable")
         call checkerror(nf90_put_var(ncid, cell_angularVID, &
               (/ 'alpha','beta ','gamma' /), &
               start = (/ 1, 1 /), count = (/ 5, 3 /)), &
               "write spatial variable")
      end if
      ! If this is the main restart, indicate it has been setup
      if (isMain) mainRestartSetup=.true.

   ! If not the first call, just reopen the existing file 
   else
      cmode = nf90_64bit_offset
      err = nf90_open(path = filename, mode = nf90_write, ncid = ncid)
      if ( err /= nf90_noerr ) then
         call checkerror(err,'write_nc_restart()');
         write(mdout,*) 'Error: write_nc_restart(): Could not open restart ', filename
         call mexit(mdout,1)
      endif
   endif

   ! ------------------------- File Write --------------------------
   ord1 = axis_flipback_ords(1)
   ord2 = axis_flipback_ords(2)
   ord3 = axis_flipback_ords(3)

   ! Write time
   call checkerror(nf90_put_var(ncid, timeVID, Time), 'write time')
   ! Write coords
   call write_nc_coords(ncid,coordVID,natom,Coords,ord1,ord2,ord3)
   ! Write velocities
   if (imin .eq. 0) then
      call write_nc_coords(ncid,velocityVID,natom,Velo,ord1,ord2,ord3)
   endif
   ! Write box information
   if (ntb > 0) then
      call checkerror(nf90_put_var(ncid,cellLengthVID, &
              (/ box(ord1),box(ord2),box(ord3) /), &
              start = (/ 1 /), count = (/ 3 /) ), 'write cell lengths')
      call checkerror(nf90_put_var(ncid,cellAngleVID, &
              (/ alpha,beta,gamma /), &
              start = (/ 1 /), count = (/ 3 /) ), 'write cell angles')
   endif
#  ifdef MPI
   ! Write replica temperature
   if (remd_method>0) then
      call checkerror(nf90_put_var(ncid, TempVID, temp0), 'write temp0')
   endif
#  endif
   ! Close restart file       
   err = nf90_close(ncid)
   call checkerror(err, "Closing netcdf restart")
#else
   write(mdout,*) 'No binary trajectory support in this version'
   write(mdout,*) 'recompile using the -DBINTRAJ flag'
   call mexit(mdout,1)
#endif
end subroutine write_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read box information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_BOX ---
!-------------------------------------------------------------------
!     Read box information from the Netcdf Restart file with 
!     specified filename.
!     The box read is called from load_ewald_info() in ew_setup.f
!     and is separate from the coord/velocity read since the box
!     information is needed to set up certain ewald parameters.
subroutine read_nc_restart_box(filename,a,b,c,alpha,beta,gamma)
#ifdef BINTRAJ
   use netcdf
#endif
   use file_io_dat_mod, only: mdout
   use pmemd_lib_mod, only: mexit
   implicit none

   character(len=*), intent(in) :: filename
   double precision, intent(out) :: a,b,c,alpha,beta,gamma
#ifdef BINTRAJ
   double precision, dimension(3) :: box
   integer :: ncid, err
   
   ! Open file
   err = nf90_open(filename,nf90_nowrite,ncid)
   call checkerror(err, "read_nc_restart_box(): Could not open coordinate file.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! Cell length VID
   err = nf90_inq_varid(ncid, "cell_lengths", cellLengthVID)
   call checkerror(err, "read_nc_restart_box(): Getting cellLengthVID.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! Cell angle VID
   err = nf90_inq_varid(ncid, "cell_angles", cellAngleVID)
   call checkerror(err, "read_nc_restart_box(): Getting cellAngleVID.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! Cell lengths
   err = nf90_get_var(ncid, cellLengthVID, box(1:3), start = (/ 1 /), count = (/ 3 /))
   call checkerror(err, "read_nc_restart_box(): Getting cell lengths.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   a = box(1)
   b = box(2)
   c = box(3)
   ! Cell angles
   err = nf90_get_var(ncid, cellAngleVID, box(1:3), start = (/ 1 /), count = (/ 3 /))
   call checkerror(err, "read_nc_restart_box(): Getting cell angles.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   alpha = box(1)
   beta  = box(2)
   gamma = box(3)
   ! Close file
   err = nf90_close(ncid)
   write(mdout,'(a)') '| check_nc_restart_box: Box info found'
#else
   write(mdout,*) 'No binary trajectory support in this version'
   write(mdout,*) 'recompile using the -DBINTRAJ flag'
   call mexit(mdout,1)
#endif
end subroutine read_nc_restart_box

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Return true if file is amber netcdf restart.
!-------------------------------------------------------------------
!     --- CHECK_NC_RESTART ---
!-------------------------------------------------------------------
!     First check that this is a netcdf file, then check the 
!     Conventions attribute for AMBERRESTART
logical function check_nc_restart(filename) 
#ifdef BINTRAJ
   use netcdf
#endif
   implicit none

   character(len=*), intent(in) :: filename
#ifdef BINTRAJ
   integer :: err, ncid
   character(len=80) :: attribute

   check_nc_restart=.false.
   ! Open file
   err = nf90_open(filename,nf90_nowrite,ncid)
   if (err==nf90_noerr) then
      ! Get Conventions
      err = nf90_get_att(ncid, nf90_global, 'Conventions', attribute)
      if (err==nf90_noerr) then
         ! Check for AMBERRESTART 
         if ( attribute .eq. "AMBERRESTART" ) then
            check_nc_restart=.true.
         endif
      endif
      ! Close file if it was opened with no errors
      err = nf90_close(ncid)
   endif
#else
   check_nc_restart=.false.
#endif
end function check_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read number of atoms from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_ATOMS ---
!-------------------------------------------------------------------
!     This routine is necessary so PMEMD can allocate memory for 
!     input coordinates before reading them in.
subroutine read_nc_restart_atoms(filename,natom)
#ifdef BINTRAJ
   use netcdf
#endif
   use file_io_dat_mod, only: mdout
   use pmemd_lib_mod, only: mexit
   implicit none
   ! Formal arguments
   character(len=*), intent(in) :: filename
   integer, intent(out) :: natom
#ifdef BINTRAJ
   ! Local variables
   integer :: ncid, err
   character(len=80) :: attribute

   ! Open file
   err = nf90_open(filename,nf90_nowrite,ncid)
   call checkerror(err, "read_nc_restart_atoms(): Could not open coordinate file.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! Check conventions 
   call checkerror( nf90_get_att(ncid, nf90_global, 'Conventions', attribute), &
        "read_nc_restart_atoms(): Getting netcdf restart conventions.")
   if ( attribute .ne. "AMBERRESTART" ) then
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      write(mdout,'(a,a)') "ERROR: INPCRD has convention that is not AMBERRESTART: ",attribute
      write(mdout,'(a)')   "       Use of this file is NOT recommended."
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      call mexit(mdout,1)
   endif
   call checkerror( nf90_get_att(ncid, nf90_global, 'ConventionVersion', attribute), &
        "read_nc_restart_atoms(): Getting netcdf restart convention version.")
   if ( attribute .ne. "1.0" ) then
      write(mdout,'(a,a)') "WARNING: INPCRD has convention version that is not 1.0: ", attribute
   endif
   ! Atom dimension and number of atoms
   if ( GetDimInfo(ncid, "atom", atomDID, natom) ) call mexit(mdout,1)
   ! Close file
   call checkerror( nf90_close(ncid), "read_nc_restart_atoms(): Closing file.")
#else
   write(mdout,*) 'No binary trajectory support in this version'
   write(mdout,*) 'recompile using the -DBINTRAJ flag'
   call mexit(mdout,1)
#endif
end subroutine read_nc_restart_atoms

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file 
!     with specified filename. Title will be read in and set. 
!     natom is the number of atoms in the restart file.
!     Coords and Velo are the coordinates and velocities, temp0 is 
!     the temperature (if present) and Time is the time.
!     box_found is set to true if restart contains box coords.
!     velocities_found set to true if restart contains velocities.
subroutine read_nc_restart(filename,title,natom,Coords,Velo,temp0,Time,&
                           box,alpha,beta,gamma,box_found,velocities_found)
#ifdef BINTRAJ
   use netcdf
#endif
   use file_io_dat_mod, only: mdout
   use pmemd_lib_mod, only: mexit
   implicit none
   ! Formal Arguments
   character(len=*), intent(in) :: filename
   character(len=80), intent(out) :: title
   integer, intent(out) :: natom
   double precision, dimension(:, :), intent(out) :: Coords, Velo
   double precision, intent(out) :: temp0, Time
   double precision, dimension(3), intent(out) :: box
   double precision, intent(out) :: alpha, beta, gamma
   logical, intent(out) :: box_found, velocities_found
#ifdef BINTRAJ
   ! Local variables
   character(len=80) :: attribute
   integer :: ncid, err, ncatom, spatial
   double precision, dimension(3) :: angles

   angles(:) = 0.d0
   ! ---=== Open file
   err = nf90_open(filename,nf90_nowrite,ncid)
   call checkerror(err, "read_nc_restart(): Could not open coordinate file.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! ---=== Read global attributes
   err = nf90_get_att(ncid, nf90_global, 'title', title)
   call checkerror(err, "read_nc_restart(): Getting netcdf restart title.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   call checkerror( nf90_get_att(ncid, nf90_global, 'Conventions', attribute), &
        "read_nc_restart(): Getting netcdf restart conventions.")
   if ( attribute .ne. "AMBERRESTART" ) then
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      write(mdout,'(a,a)') "ERROR: INPCRD has convention that is not AMBERRESTART: ",attribute
      write(mdout,'(a)')   "       Use of this file is NOT recommended."
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      call mexit(mdout,1)
   endif
   call checkerror( nf90_get_att(ncid, nf90_global, 'ConventionVersion', attribute), &
        "read_nc_restart(): Getting netcdf restart convention version.")
   if ( attribute .ne. "1.0" ) then
      write(mdout,'(a,a)') "WARNING: INPCRD has convention version that is not 1.0: ", attribute
   endif
   ! ---=== Atom dimension and number of atoms
   if ( GetDimInfo(ncid, "atom", atomDID, ncatom) ) call mexit(mdout,1)
   natom = ncatom
   ! ---=== Coords Variable ID and units check
   ! NOTE: Move spatial above this?
   err = nf90_inq_varid(ncid, "coordinates", coordVID)
   call checkerror(err, "read_nc_restart(): Getting coordinates VID")
   if (err/=nf90_noerr) call mexit(mdout,1)
   call checkerror( nf90_get_att(ncid, coordVID, "units", attribute), &
        "read_nc_restart(): Getting coordinate units")
   if (attribute .ne. "angstrom") then
      write(mdout,'(a,a)') "WARNING: INPCRD has coordinate units not angstrom: ", attribute
   endif
   err = nf90_get_var(ncid, coordVID, Coords(1:3,1:ncatom), &
                      start = (/ 1, 1 /), count = (/ 3, ncatom /))
   call checkerror(err, "read_nc_restart(): Getting coordinates")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! ---=== Spatial Dimension and VID
   if ( GetDimInfo(ncid, "spatial", spatialDID, spatial) ) call mexit(mdout,1)
   if (spatial /= 3) then
      write(mdout,'(a,i6)') "ERROR: read_nc_restart(): expected 3 spatial dimensions, got ", spatial
      call mexit(mdout,1)
   endif
   err = nf90_inq_varid(ncid, "spatial", spatialVID)
   call checkerror(err, "read_nc_restart(): Getting spatial VID")
   if (err/=nf90_noerr) call mexit(mdout,1)
   ! ---=== Velocity Variable ID and units check
   !        Just check for the velocityVID to determine if velocities are
   !        present, and let the calling routine worry about ntx
   err = nf90_inq_varid(ncid, "velocities", velocityVID)
   if (err==nf90_noerr) then
      call checkerror(err, "read_nc_restart(): Getting velocities VID")
      if (err/=nf90_noerr) call mexit(mdout,1)
      call checkerror( nf90_get_att(ncid, velocityVID, "units", attribute), &
           "read_nc_restart(): Getting velocity units")
      if (attribute .ne. "angstrom/picosecond") then
        write(mdout,'(a,a)') "WARNING: INPCRD has velocity units not angstrom/picosecond: ", attribute
      endif
      ! NOTE: Check scale_factor?
      err = nf90_get_var(ncid, velocityVID, Velo(1:3,1:ncatom), &
                         start = (/ 1, 1 /), count = (/ 3, ncatom /))
      call checkerror(err, "read_nc_restart(): Getting velocities")
      if (err /= nf90_noerr) call mexit(mdout,1)
      velocities_found=.true.
   endif
   ! ---=== Restart Time Info
   err = nf90_inq_varid(ncid, "time", timeVID)
   call checkerror(err, "read_nc_restart(): Getting time VID")
   if (err/=nf90_noerr) call mexit(mdout,1)
   call checkerror( nf90_get_att(ncid, timeVID, "units", attribute), &
        "read_nc_restart(): Getting time units")
   if (attribute .ne. "picosecond") then
      write(mdout,'(a,a)') "WARNING: INPCRD has time units not picosecond: ", attribute
   endif
   err = nf90_get_var(ncid, timeVID, Time)
   call checkerror(err, "read_nc_restart(): Getting restart time")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! ---=== Box Info 
   !        Check for Cell length VID to determine if box info is present
   err = nf90_inq_varid(ncid, "cell_lengths", cellLengthVID)
   if (err == nf90_noerr) then 
     !call checkerror(err, "read_nc_restart_box(): Getting cellLengthVID.")
     !if (err /= nf90_noerr) call mexit(mdout,1)
     ! Cell angle VID
     err = nf90_inq_varid(ncid, "cell_angles", cellAngleVID)
     call checkerror(err, "read_nc_restart_box(): Getting cellAngleVID.")
     if (err /= nf90_noerr) call mexit(mdout,1)
     ! Cell lengths
     err = nf90_get_var(ncid, cellLengthVID, box(1:3), start = (/ 1 /), count = (/ 3 /))
     call checkerror(err, "read_nc_restart_box(): Getting cell lengths.")
     if (err /= nf90_noerr) call mexit(mdout,1)
     ! Cell angles
     err = nf90_get_var(ncid, cellAngleVID, angles(1:3), start = (/ 1 /), count = (/ 3 /))
     call checkerror(err, "read_nc_restart_box(): Getting cell angles.")
     if (err /= nf90_noerr) call mexit(mdout,1)
     alpha = angles(1)
     beta  = angles(2)
     gamma = angles(3)
     box_found=.true.
   end if
   ! ---=== Replica Temperature
   temp0=0.0d0
   err = nf90_inq_varid(ncid, "temp0", TempVID)
   if (err == nf90_noerr) then
      err = nf90_get_var(ncid, TempVID, temp0)
      call checkerror(err, "read_nc_restart(): Getting restart temperature")
      if (err /= nf90_noerr) temp0=0.0d0
   endif

   ! NOTE: TO BE ADDED
   !labelDID;
   !int cell_spatialDID, cell_angularDID;
   !int spatialVID, cell_spatialVID, cell_angularVID;
  
   ! ---=== Close file
   err = nf90_close(ncid) 
#else
   write(mdout,*) 'No binary trajectory support in this version'
   write(mdout,*) 'recompile using the -DBINTRAJ flag'
   call mexit(mdout,1)
#endif
end subroutine read_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REFC ---
!-------------------------------------------------------------------
!     Read coordinates from the Netcdf Restart file with specified 
!     filename. Title will be read in and set. 
!     natom is the expected number of atoms in the restart file.
!     Coords are the coordinates. 
subroutine read_nc_refc(filename,title,natom,Coords)
#ifdef BINTRAJ
   use netcdf
#endif
   use file_io_dat_mod, only: mdout
   use pmemd_lib_mod, only: mexit
   implicit none
   ! Formal Arguments
   character(len=*), intent(in) :: filename
   character(len=80), intent(out) :: title
   integer, intent(in) :: natom
   double precision, dimension(*), intent(out) :: Coords
#ifdef BINTRAJ
   ! Local variables
   character(len=80) :: attribute
   integer :: ncid, err, ncatom, nr3, spatial

   ! ---=== Open file
   err = nf90_open(filename,nf90_nowrite,ncid)
   call checkerror(err, "read_nc_restart(): Could not open coordinate file.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! ---=== Read global attributes
   err = nf90_get_att(ncid, nf90_global, 'title', title)
   call checkerror(err, "read_nc_restart(): Getting netcdf restart title.")
   if (err /= nf90_noerr) call mexit(mdout,1)
   call checkerror( nf90_get_att(ncid, nf90_global, 'Conventions', attribute), &
        "read_nc_restart(): Getting netcdf restart conventions.")
   if ( attribute .ne. "AMBERRESTART" ) then
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      write(mdout,'(a,a)') "ERROR: INPCRD has convention that is not AMBERRESTART: ",attribute
      write(mdout,'(a)')   "       Use of this file is NOT recommended."
      write(mdout,'(a)')   "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      call mexit(mdout,1)
   endif
   call checkerror( nf90_get_att(ncid, nf90_global, 'ConventionVersion', attribute), &
        "read_nc_restart(): Getting netcdf restart convention version.")
   if ( attribute .ne. "1.0" ) then
      write(mdout,'(a,a)') "WARNING: INPCRD has convention version that is not 1.0: ", attribute
   endif
   ! ---=== Atom dimension and number of atoms
   if ( GetDimInfo(ncid, "atom", atomDID, ncatom) ) call mexit(mdout,1)
   if (natom .ne. ncatom) then
     write(mdout, '(/2x,a)') 'FATAL: NATOM mismatch in NetCDF constraint coord and prmtop files'
     call mexit(6, 1)
   end if
   nr3 = ncatom * 3
   ! ---=== Coords Variable ID and units check
   ! NOTE: Move spatial above this?
   err = nf90_inq_varid(ncid, "coordinates", coordVID)
   call checkerror(err, "read_nc_restart(): Getting coordinates VID")
   if (err/=nf90_noerr) call mexit(mdout,1)
   call checkerror( nf90_get_att(ncid, coordVID, "units", attribute), &
        "read_nc_restart(): Getting coordinate units")
   if (attribute .ne. "angstrom") then
      write(mdout,'(a,a)') "WARNING: INPCRD has coordinate units not angstrom: ", attribute
   endif
   err = nf90_get_var(ncid, coordVID, Coords(1:nr3), &
                      start = (/ 1, 1 /), count = (/ 3, ncatom /))
   call checkerror(err, "read_nc_restart(): Getting coordinates")
   if (err /= nf90_noerr) call mexit(mdout,1)
   ! ---=== Spatial Dimension and VID
   if ( GetDimInfo(ncid, "spatial", spatialDID, spatial) ) call mexit(mdout,1)
   if (spatial /= 3) then
      write(mdout,'(a,i6)') "ERROR: read_nc_restart(): expected 3 spatial dimensions, got ", spatial
      call mexit(mdout,1)
   endif
   err = nf90_inq_varid(ncid, "spatial", spatialVID)
   call checkerror(err, "read_nc_restart(): Getting spatial VID")
   if (err/=nf90_noerr) call mexit(mdout,1)
   ! ---=== Close file
   err = nf90_close(ncid) 
#else
   write(mdout,*) 'No binary trajectory support in this version'
   write(mdout,*) 'recompile using the -DBINTRAJ flag'
   call mexit(mdout,1)
#endif
end subroutine read_nc_refc

! ======================== PRIVATE SUBROUTINES =========================
#ifdef BINTRAJ
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Get dimension ID and length of given attribute.
!-------------------------------------------------------------------
!     --- GetDimInfo ---
!-------------------------------------------------------------------
!     Given an attribute, set the dimension ID and length.
!     Return true if successful, false if not.
logical function GetDimInfo(ncid, attribute, dimID, length)
   use netcdf
   use file_io_dat_mod, only: mdout
   implicit none

   integer, intent(in) :: ncid
   character(len=*), intent(in) :: attribute
   integer, intent(out) :: dimID, length

   integer :: err

   GetDimInfo=.false.
   ! First get dimension ID
   err = nf90_inq_dimid(ncid, attribute, dimID)
   call checkerror(err, "GetDimInfo(): Getting dimension ID")
   if (err /= nf90_noerr) then
      write(mdout,'(2x,a,a)') 'On attribute ', attribute
      GetDimInfo=.true. 
      return
   endif
   ! Then get dimension length
   err = nf90_inquire_dimension(ncid, dimID, len = length)
   call checkerror(err, "GetDimInfo(): Getting dimension length")
   if (err /= nf90_noerr) then
      write(mdout,'(2x,a,a)') 'On attribute ', attribute
      GetDimInfo=.true.
      return
   endif
end function GetDimInfo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write given Netcdf coordinates.
!-------------------------------------------------------------------
!     --- WRITE_NC_COORDS ---
!-------------------------------------------------------------------
!     Write given coordinates to VID in given ncid. Based on the 
!     axis, write flipped if necessary.
subroutine write_nc_coords(ncid,VID,natom,arrayIn,ord1,ord2,ord3)
  use netcdf
  implicit none
  integer, intent(in) :: ncid, VID, natom
  double precision, dimension(3, natom), intent(in) :: arrayIn
  integer, intent(in) :: ord1, ord2, ord3
  ! Local variables for dealing with flipped axis
  integer       :: i
  real          :: buf(3, natom)

  ! Normal axis
  if (ord1 .eq. 1 .and. ord2 .eq. 2) then

    call checkerror(nf90_put_var(ncid, VID, arrayIn(:,:), &
                                 start=(/ 1, 1 /), &
                                 count=(/ 3, natom /)), &
                    'NetCDF write coords')
  ! Flipped axis
  else
    do i = 1, natom
      buf(1, i) = arrayIn(ord1, i)
      buf(2, i) = arrayIn(ord2, i)
      buf(3, i) = arrayIn(ord3, i)
    end do

    call checkerror(nf90_put_var(ncid, VID, arrayIn(:,:), &
                                 start=(/ 1, 1 /), &
                                 count=(/ 3, natom /)), &
                  'NetCDF write flipped coords')

  end if
end subroutine write_nc_coords

!*******************************************************************************!
! Subroutine:   checkerror
!
! Description:  Checks error return from netCDF routines. Put here to avoid
!               dependency on bintraj
!
!*******************************************************************************
subroutine checkerror(status, location)
  use netcdf
  use file_io_dat_mod, only: mdout
  
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
#endif

end module binrestart_mod
