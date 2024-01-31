! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!! Module for generating binary restart-type output in NetCDF format
!! Developed by Dan Roe
!! 2010-01-10

module binrestart
    private

    integer, save :: atomDID, coordVID, velocityVID
    integer, save :: cellAngleVID, cellLengthVID, spatialDID, labelDID
    integer, save :: cell_spatialDID, cell_angularDID, spatialVID, timeVID
    integer, save :: cell_spatialVID, cell_angularVID, TempVID
    logical, save :: hasV

    public write_nc_restart, read_nc_restart_box, read_nc_restart, check_nc_restart
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write Netcdf restart file.
!-------------------------------------------------------------------
!     --- WRITE_NC_RESTART ---
!-------------------------------------------------------------------
!     Write Netcdf Restart file with given filename and title.
!     owrite indicates overwrite status (N is no overwrite), natom
!     is the # atoms, ntb>0 indicates presence of box coords, first
!     indicates the file should be created and set-up, Coords and
!     Velo are the coordinates and velocities, temp0 is the current
!     temperature and Time is the current simulation time. If hasV
!     is false, no velocities will be written (for e.g. during min)
    subroutine write_nc_restart(filename, title, owrite, natom, ntb, first, Coords, Velo, temp0, Time, hasVin)
#ifdef BINTRAJ
        use netcdf
        use bintraj, only : checkNCerror
        use nblist, only : a, b, c, alpha, beta, gamma
#  ifdef MPI
        use remd, only : rem
#  endif
#endif
        implicit none

        character(len=*), intent(in) :: filename
        character(len=*), intent(in)  :: title
        character, intent(in) :: owrite
        integer, intent(in) :: natom, ntb
        logical, intent(in) :: first
        _REAL_, dimension(*), intent(in) :: Coords, Velo
        _REAL_, intent(in) :: temp0, Time
        logical, intent(in) :: hasVin
#ifdef BINTRAJ
        integer :: ncid, natom3
        integer :: cmode, err, oldMode

        ! If first call, create the file and set up all dimensions and vars
        if (first) then
            ! owrite status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
            ! sander flag -O='R', -A='U', default='N'
            cmode = nf90_64bit_offset
            if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
            !if (owrite == 'U' .and. facc == 'A') cmode = ior(cmode, nf90_noclobber)
            err = nf90_create(path=filename, cmode=cmode, ncid=ncid)
            if (err == nf90_eexist) then
                write (6, *) 'Error: write_nc_restart(): File exists and -O not specified: ', filename
                call mexit(6, 1)
            end if
            call checkNCerror(err, "write_nc_restart(): Creating netcdf restart")
            if (err /= nf90_noerr) call mexit(6, 1)
            ! Time variable
            call checkNCerror(nf90_def_var(ncid, "time", nf90_double, timeVID), "define timeVID")
            call checkNCerror(nf90_put_att(ncid, timeVID, "units", "picosecond"), &
                "define timeVID units")
            ! Spatial dimension and variable
            call checkNCerror(nf90_def_dim(ncid, "spatial", 3, spatialDID), "define spatialDID")
            call checkNCerror(nf90_def_var(ncid, "spatial", nf90_char, &
                (/spatialDID/), spatialVID), "define spatialVID")
            ! Atom dimension
            call checkNCerror(nf90_def_dim(ncid, "atom", natom, atomDID), "define atomDID")
            ! Coord variable
            call checkNCerror(nf90_def_var(ncid, "coordinates", nf90_double, &
                (/spatialDID, atomDID/), coordVID), "define coordVID")
            call checkNCerror(nf90_put_att(ncid, coordVID, "units", "angstrom"), &
                "define coordVID units")
            ! Velocity variable
            hasV = hasVin
            if (hasV) then
                call checkNCerror(nf90_def_var(ncid, "velocities", nf90_double, &
                    (/spatialDID, atomDID/), velocityVID), "define velocityVID")
                call checkNCerror(nf90_put_att(ncid, velocityVID, "units", "angstrom/picosecond"), &
                    "define velocityVID units")
                call checkNCerror(nf90_put_att(ncid, velocityVID, "scale_factor", 20.455), &
                    "define velocityVID scale factor")
            end if
            ! Box info
            if (ntb > 0) then
                ! Cell Spatial
                call checkNCerror(nf90_def_dim(ncid, "cell_spatial", 3, cell_spatialDID), &
                    "define cell_spatialDID")
                call checkNCerror(nf90_def_var(ncid, "cell_spatial", nf90_char, &
                    (/cell_spatialDID/), cell_spatialVID), "define cell_spatialVID")
                ! Cell angular
                call checkNCerror(nf90_def_dim(ncid, "label", 5, labelDID), "define labelDID")
                call checkNCerror(nf90_def_dim(ncid, "cell_angular", 3, cell_angularDID), &
                    "define cell_angularDID")
                call checkNCerror(nf90_def_var(ncid, "cell_angular", nf90_char, &
                    (/labelDID, cell_angularDID/), cell_angularVID), &
                    "define cell_angularVID")
                ! Cell length
                call checkNCerror(nf90_def_var(ncid, "cell_lengths", nf90_double, &
                    (/cell_spatialDID/), cellLengthVID), "define cellLengthVID")
                call checkNCerror(nf90_put_att(ncid, cellLengthVID, "units", "angstrom"), &
                    "define cellLengthVID units")
                ! Cell angle
                call checkNCerror(nf90_def_var(ncid, "cell_angles", nf90_double, &
                    (/cell_angularDID/), cellAngleVID), "define cellAngleVID")
                call checkNCerror(nf90_put_att(ncid, cellAngleVID, "units", "degree"), &
                    "define cellAngleVID units")
            end if
#  ifdef MPI
            ! Replica Temperature
            if (rem > 0) then
                call checkNCerror(nf90_def_var(ncid, "temp0", nf90_double, TempVID), &
                    "define TempVID")
                call checkNCerror(nf90_put_att(ncid, TempVID, "units", "kelvin"), &
                    "define TempVID units")
            end if
#  endif
            ! Global attributes: Title etc
            call checkNCerror(nf90_put_att(ncid, nf90_global, &
                'title', title), "define title")
            call checkNCerror(nf90_put_att(ncid, nf90_global, "application", &
                'AMBER'), "define application")
            call checkNCerror(nf90_put_att(ncid, nf90_global, "program", &
                'sander'), "define program")
            call checkNCerror(nf90_put_att(ncid, nf90_global, "programVersion", &
                '11.0'), "define programVersion")
            call checkNCerror(nf90_put_att(ncid, nf90_global, "Conventions", &
                'AMBERRESTART'), "define Convention")
            call checkNCerror(nf90_put_att(ncid, nf90_global, "ConventionVersion", &
                '1.0'), "define ConventionVersion")
            ! Set fill mode
            call checkNCerror(nf90_set_fill(ncid, nf90_nofill, oldMode), "Setting fill mode")
            ! End definitions
            call checkNCerror(nf90_enddef(ncid), "end define")
            ! Specify dimension labels
            call checkNCerror(nf90_put_var(ncid, spatialVID, &
                (/'x', 'y', 'z'/), start=(/1/), count=(/3/)), "write spatial variable")
            if (ntb > 0) then
                call checkNCerror(nf90_put_var(ncid, cell_spatialVID, &
                    (/'a', 'b', 'c'/), start=(/1/), count=(/3/)), &
                    "write spatial variable")
                call checkNCerror(nf90_put_var(ncid, cell_angularVID, &
                    (/'alpha', 'beta ', 'gamma'/), &
                    start=(/1, 1/), count=(/5, 3/)), &
                    "write spatial variable")
            end if

            ! If not the first call, just reopen the existing file
        else
            cmode = nf90_64bit_offset
            err = nf90_open(path=filename, mode=nf90_write, ncid=ncid)
            if (err /= nf90_noerr) then
                call checkNCerror(err);
                write (6, *) 'Error: write_nc_restart(): Could not open restart ', filename
                call mexit(6, 1)
            end if
        end if

        natom3 = natom*3
        ! Write time
        call checkNCerror(nf90_put_var(ncid, timeVID, Time), 'write time')
        ! Write coords
        call checkNCerror(nf90_put_var(ncid, coordVID, Coords(1:natom3), &
            start=(/1, 1/), count=(/3, natom/)), 'write atom coords')
        ! Write velocities
        if (hasV) then
            call checkNCerror(nf90_put_var(ncid, velocityVID, Velo(1:natom3), &
                start=(/1, 1/), count=(/3, natom/)), 'write velocities')
        end if
        ! Write box information
        if (ntb > 0) then
            call checkNCerror(nf90_put_var(ncid, cellLengthVID, &
                (/a, b, c/), start=(/1/), count=(/3/)), 'write cell lengths')
            call checkNCerror(nf90_put_var(ncid, cellAngleVID, &
                (/alpha, beta, gamma/), start=(/1/), count=(/3/)), &
                'write cell angles')
        end if
#  ifdef MPI
        ! Write replica temperature
        if (rem > 0) then
            call checkNCerror(nf90_put_var(ncid, TempVID, temp0), 'write temp0')
        end if
#  endif
        ! Close restart file
        err = nf90_close(ncid)
        call checkNCerror(err, "Closing netcdf restart")
#else
        write (6, *) 'No binary trajectory support in this version'
        write (6, *) 'recompile using the -DBINTRAJ flag'
        call mexit(6, 1)
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
    subroutine read_nc_restart_box(filename, a, b, c, alpha, beta, gamma)
#ifdef BINTRAJ
        use netcdf
        use bintraj, only : checkNCerror
#endif
        implicit none

        character(len=*), intent(in) :: filename
        _REAL_, intent(out) :: a, b, c, alpha, beta, gamma
#ifdef BINTRAJ
        _REAL_, dimension(3) :: box
        integer :: ncid, err

        ! Open file
        err = nf90_open(filename, nf90_nowrite, ncid)
        call checkNCerror(err, "read_nc_restart_box(): Could not open coordinate file.")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! Cell length VID
        err = nf90_inq_varid(ncid, "cell_lengths", cellLengthVID)
        call checkNCerror(err, "read_nc_restart_box(): Getting cellLengthVID.")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! Cell angle VID
        err = nf90_inq_varid(ncid, "cell_angles", cellAngleVID)
        call checkNCerror(err, "read_nc_restart_box(): Getting cellAngleVID.")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! Cell lengths
        err = nf90_get_var(ncid, cellLengthVID, box(1:3), start=(/1/), count=(/3/))
        call checkNCerror(err, "read_nc_restart_box(): Getting cell lengths.")
        if (err /= nf90_noerr) call mexit(6, 1)
        a = box(1)
        b = box(2)
        c = box(3)
        ! Cell angles
        err = nf90_get_var(ncid, cellAngleVID, box(1:3), start=(/1/), count=(/3/))
        call checkNCerror(err, "read_nc_restart_box(): Getting cell angles.")
        if (err /= nf90_noerr) call mexit(6, 1)
        alpha = box(1)
        beta = box(2)
        gamma = box(3)
        ! Close file
        err = nf90_close(ncid)
        write (6, '(a)') '| check_nc_restart_box: Box info found'
#else
        write (6, *) 'No binary trajectory support in this version'
        write (6, *) 'recompile using the -DBINTRAJ flag'
        call mexit(6, 1)
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

        check_nc_restart = .false.
        ! Open file
        err = nf90_open(filename, nf90_nowrite, ncid)
        if (err == nf90_noerr) then
            ! Get Conventions
            err = nf90_get_att(ncid, nf90_global, 'Conventions', attribute)
            if (err == nf90_noerr) then
                ! Check for AMBERRESTART
                if (attribute .eq. "AMBERRESTART") then
                    check_nc_restart = .true.
                end if
            end if
            ! Close file if it was opened with no errors
            err = nf90_close(ncid)
        end if
#else
        check_nc_restart = .false.
#endif
    end function check_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file
!     with specified filename. This is called from getcor.f. Title
!     will be read in and set.
!     ntx specifies whether coords and velo or just coords will be
!     read. ntx=1 means read coords only, and ntx=5 means read
!     coords and velocities.
!     parmatoms is the expected number of atoms in the restart.
!     Coords and Velo are the coordinates and
!     velocities, temp0 is the temperature (if present) and Time
!     is the time.
!     NOTE: Box info is not read here; it is obtained using
!     read_nc_restart_box in load_ewald_info.
    subroutine read_nc_restart(filename, title, ntx, parmatoms, Coords, Velo, temp0, Time)
#ifdef BINTRAJ
        use netcdf
        use bintraj, only : checkNCerror
#endif
        implicit none

        character(len=*), intent(in) :: filename
        character(len=80), intent(out) :: title
        integer, intent(in) :: ntx, parmatoms
        _REAL_, dimension(*), intent(out) :: Coords, Velo
        _REAL_, intent(out) :: temp0, Time
#ifdef BINTRAJ
        character(len=80) :: attribute
        integer :: ncid, err, ncatom, ncatom3, spatial

        ! ---=== Open file
        err = nf90_open(filename, nf90_nowrite, ncid)
        call checkNCerror(err, "read_nc_restart(): Could not open coordinate file.")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! ---=== Read global attributes
        err = nf90_get_att(ncid, nf90_global, 'title', title)
        call checkNCerror(err, "read_nc_restart(): Getting netcdf restart title.")
        if (err /= nf90_noerr) call mexit(6, 1)
        call checkNCerror(nf90_get_att(ncid, nf90_global, 'Conventions', attribute), &
            "read_nc_restart(): Getting netcdf restart conventions.")
        if (attribute .ne. "AMBERRESTART") then
            write (6, '(a)') "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
            write (6, '(a,a)') "ERROR: INPCRD has convention that is not AMBERRESTART: ", attribute
            write (6, '(a)') "       Use of this file is NOT recommended with ntx 8 or 9!"
            write (6, '(a)') "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
            call mexit(6, 1)
        end if
        call checkNCerror(nf90_get_att(ncid, nf90_global, 'ConventionVersion', attribute), &
            "read_nc_restart(): Getting netcdf restart convention version.")
        if (attribute .ne. "1.0") then
            write (6, '(a,a)') "WARNING: INPCRD has convention version that is not 1.0: ", attribute
        end if
        ! ---=== Atom dimension and number of atoms
        if (GetDimInfo(ncid, "atom", atomDID, ncatom)) call mexit(6, 1)
        if (ncatom /= parmatoms) then
            write (6, '(2x,a)') "FATAL: NATOM mismatch in coord and topology files."
            call mexit(6, 1)
        end if
        ncatom3 = ncatom*3
        ! Spatial Dimension and VID
        if (GetDimInfo(ncid, "spatial", spatialDID, spatial)) call mexit(6, 1)
        if (spatial /= 3) then
            write (6, '(a,i6)') "ERROR: read_nc_restart(): expected 3 spatial dimensions, got ", spatial
            call mexit(6, 1)
        end if
        err = nf90_inq_varid(ncid, "spatial", spatialVID)
        call checkNCerror(err, "read_nc_restart(): Getting spatial VID")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! ---=== Coords Variable ID and units check
        err = nf90_inq_varid(ncid, "coordinates", coordVID)
        call checkNCerror(err, "read_nc_restart(): Getting coordinates VID")
        if (err /= nf90_noerr) call mexit(6, 1)
        call checkNCerror(nf90_get_att(ncid, coordVID, "units", attribute), &
            "read_nc_restart(): Getting coordinate units")
        if (attribute .ne. "angstrom") then
            write (6, '(a,a)') "WARNING: INPCRD has coordinate units not angstrom: ", attribute
        end if
        err = nf90_get_var(ncid, coordVID, Coords(1:ncatom3), &
            start=(/1, 1/), count=(/3, ncatom/))
        call checkNCerror(err, "read_nc_restart(): Getting coordinates")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! ---=== Velocity Variable ID and units check
        !        ntx=1 No Velocity Read
        !        ntx=5 Read Velocity
        if (ntx == 5) then
            err = nf90_inq_varid(ncid, "velocities", velocityVID)
            call checkNCerror(err, "read_nc_restart(): Getting velocities VID")
            if (err /= nf90_noerr) call mexit(6, 1)
            call checkNCerror(nf90_get_att(ncid, velocityVID, "units", attribute), &
                "read_nc_restart(): Getting velocity units")
            if (attribute .ne. "angstrom/picosecond") then
                write (6, '(a,a)') "WARNING: INPCRD has velocity units not angstrom/picosecond: ", attribute
            end if
            ! NOTE: Check scale_factor?
            err = nf90_get_var(ncid, velocityVID, Velo(1:ncatom3), &
                start=(/1, 1/), count=(/3, ncatom/))
            call checkNCerror(err, "read_nc_restart(): Getting velocities")
            if (err /= nf90_noerr) call mexit(6, 1)
        end if
        ! ---=== Restart Time Info
        err = nf90_inq_varid(ncid, "time", timeVID)
        call checkNCerror(err, "read_nc_restart(): Getting time VID")
        if (err /= nf90_noerr) call mexit(6, 1)
        call checkNCerror(nf90_get_att(ncid, timeVID, "units", attribute), &
            "read_nc_restart(): Getting time units")
        if (attribute .ne. "picosecond") then
            write (6, '(a,a)') "WARNING: INPCRD has time units not picosecond: ", attribute
        end if
        err = nf90_get_var(ncid, timeVID, Time)
        call checkNCerror(err, "read_nc_restart(): Getting restart time")
        if (err /= nf90_noerr) call mexit(6, 1)
        ! ---=== Replica Temperature
        temp0 = 0.0d0
        err = nf90_inq_varid(ncid, "temp0", TempVID)
        if (err == nf90_noerr) then
            err = nf90_get_var(ncid, TempVID, temp0)
            call checkNCerror(err, "read_nc_restart(): Getting restart temperature")
            if (err /= nf90_noerr) temp0 = 0.0d0
        end if

        ! NOTE: TO BE ADDED
        !labelDID;
        !int cell_spatialDID, cell_angularDID;
        !int spatialVID, cell_spatialVID, cell_angularVID;

        ! ---=== Close file
        err = nf90_close(ncid)
#else
        write (6, *) 'No binary trajectory support in this version'
        write (6, *) 'recompile using the -DBINTRAJ flag'
        call mexit(6, 1)
#endif
    end subroutine read_nc_restart

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
        use bintraj, only : checkNCerror
        implicit none

        integer, intent(in) :: ncid
        character(len=*), intent(in) :: attribute
        integer, intent(out) :: dimID, length

        integer :: err

        GetDimInfo = .false.
        ! First get dimension ID
        err = nf90_inq_dimid(ncid, attribute, dimID)
        call checkNCerror(err, "GetDimInfo(): Getting dimension ID")
        if (err /= nf90_noerr) then
            write (6, '(2x,a,a)') 'On attribute ', attribute
            GetDimInfo = .true.
            return
        end if
        ! Then get dimension length
        err = nf90_inquire_dimension(ncid, dimID, len=length)
        call checkNCerror(err, "GetDimInfo(): Getting dimension length")
        if (err /= nf90_noerr) then
            write (6, '(2x,a,a)') 'On attribute ', attribute
            GetDimInfo = .true.
            return
        end if
    end function GetDimInfo

#endif

end module binrestart
