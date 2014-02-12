! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!! Module for generating binary trajectory-type output in NetCDF format
!! Developed by John Mongan <jmongan@mccammon.ucsd.edu>, November 2005

module bintraj
   private

   integer :: crd_ncid, vel_ncid, oldMode, cmode, atomCnt
   integer :: FrameDimID,SpatialDimID, AtomDimID, LabelDimID, Cell_spatialDimID, Cell_angularDimID
   integer :: CoordVarID, Cell_lengthVarID, Cell_angleVarID, SpatialVarID, VelocVarID, TempVarID
   integer :: Cell_spatialVarID, Cell_angularVarID, crd_TimeVarID, vel_TimeVarID
   integer :: crd_frame, vel_frame

   public open_binary_files, close_binary_files, write_binary_traj, end_binary_frame, checkNCerror
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open the coordinate, velocity file(s)
subroutine open_binary_files

#ifdef BINTRAJ
   use file_io_dat
   use netcdf
#  ifdef MPI
      use remd, only: rem
#  endif 

   implicit none

#  include "md.h"
#  include "memory.h"
#  include "box.h"
#  define NCLABELLEN 5

   logical :: crd_file_exists = .false., vel_file_exists = .false.

   if (ntwv == 0 .and. ntwx == 0) return

   atomCnt = natom
   if (ntwprt > 0) atomCnt = ntwprt
   
   ! Define coordinate variables
   if (ntwx > 0 ) then
      call open_nc_file(crd_ncid,mdcrd,owrite,facc,atomCnt,title, &
            crd_file_exists,crd_frame,crd_TimeVarID)
      if (crd_file_exists) then
         call checkNCerror(nf90_inq_varid(crd_ncid, "coordinates", CoordVarID), &
               "query coordinates")
      else
         call checkNCerror(nf90_def_var(crd_ncid, "coordinates", nf90_float, &
               (/ SpatialDimID, AtomDimID, FrameDimID /), CoordVarID), "define coordinates")
         call checkNCerror(nf90_put_att(crd_ncid, CoordVarID, "units", &
               "angstrom"), "define coordinates units")
      end if
#  ifdef MPI
      ! For remd define replica temperature
      if (rem>0 .and. rem/=3) then
         call checkNCerror(nf90_def_var(crd_ncid,"temp0",nf90_double, &
               (/ FrameDimID /), TempVarID), "define temp0")
         call checkNCerror(nf90_put_att(crd_ncid, TempVarID, "units", &
               "kelvin"), "define temp0 units")
      endif
#  endif
      ! Define unit cell data for periodic simulations
      if (ntb > 0 ) then
         if (crd_file_exists) then
            call checkNCerror(nf90_inq_varid(crd_ncid, "cell_lengths", Cell_lengthVarID), &
                  "query cell lengths")
            call checkNCerror(nf90_inq_varid(crd_ncid, "cell_angles",  Cell_angleVarID), &
                  "query cell angles")
         else
            ! Dimensions
            call checkNCerror(nf90_def_dim(crd_ncid,"cell_spatial", 3, Cell_spatialDimID), &
                  "define cell spatial dim")
            call checkNCerror(nf90_def_dim(crd_ncid,"cell_angular", 3, Cell_angularDimID), &
                  "define cell angular dim")
            ! Label variables
            call checkNCerror(nf90_def_var(crd_ncid, "cell_spatial", nf90_char, &
                  (/  Cell_spatialDimID /), Cell_spatialVarID),"define cell spatial var" )
            call checkNCerror(nf90_def_var(crd_ncid, "cell_angular", nf90_char, &
                  (/ LabelDimID, Cell_angularDimID /), Cell_angularVarID), "define cell angular var")
            ! Data variables and attributes
            call checkNCerror(nf90_def_var(crd_ncid, "cell_lengths", nf90_double, &
                  (/ Cell_spatialDimID, FrameDimID /), Cell_lengthVarID), "define cell lengths")
            call checkNCerror(nf90_put_att(crd_ncid, Cell_lengthVarID, "units", &
                  "angstrom"), "define cell length units")
            call checkNCerror(nf90_def_var(crd_ncid, "cell_angles", nf90_double, &
                  (/ Cell_angularDimID, FrameDimID /), Cell_angleVarID), "define cell angles")
            call checkNCerror(nf90_put_att(crd_ncid, Cell_angleVarID, "units", &
                  "degree"), "define cell angle units")
         end if
      end if
   end if

   ! Point vel_ncid at either a separate file or the combined (crd) traj file
   if (ntwv > 0) then
      call open_nc_file(vel_ncid,mdvel,owrite,facc,atomCnt,title, &
            vel_file_exists,vel_frame,vel_TimeVarID)
   else if (ntwv == -1) then
      vel_ncid = crd_ncid
      vel_file_exists = crd_file_exists
      vel_frame = crd_frame
   end if

   ! Define velocity variables in file pointed to by vel_ncid
   if (ntwv /= 0 ) then
      if (vel_file_exists) then
         call checkNCerror(nf90_inq_varid(vel_ncid, "velocities", VelocVarID), &
               "query velocities")
      else
         call checkNCerror(nf90_def_var(vel_ncid, "velocities", nf90_float, &
               (/ SpatialDimID, AtomDimID, FrameDimID /), VelocVarID), &
               "define velocities")
         call checkNCerror(nf90_put_att(vel_ncid, VelocVarID, "units", &
               "angstrom/picosecond"), "define velocity units")
         call checkNCerror(nf90_put_att(vel_ncid, VelocVarID, "scale_factor", &
               20.455), "define velocity scale_factor")
      end if
   end if


   ! Prepare files for data writing
   if (ntwx > 0 .and. .not. crd_file_exists) then 
      call end_nc_define(crd_ncid)
      ! Fill dimension label variables
      if (ntb > 0) then
         call checkNCerror(nf90_put_var(crd_ncid, Cell_spatialVarID, &
               (/ 'a','b','c' /), start = (/ 1 /), count = (/ 3 /)), &
               "write spatial variable")
         call checkNCerror(nf90_put_var(crd_ncid, Cell_angularVarID, &
               (/ 'alpha','beta ','gamma' /), &
               start = (/ 1, 1 /), count = (/ NCLABELLEN, 3 /)), &
               "write spatial variable")
      end if
   end if
   if (ntwv > 0 .and. .not. vel_file_exists) call end_nc_define(vel_ncid)
   
#else
   write(6,*) 'No binary trajectory support in this version'
   write(6,*) 'reconfigure using the -bintraj flag'
   call mexit(6,1)
#endif

end subroutine open_binary_files
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Close the binary coordinate, velocity file(s).
subroutine close_binary_files
#ifdef BINTRAJ
   use file_io_dat
   use netcdf

   implicit none

   if (ntwx > 0) call checkNCerror(nf90_close(crd_ncid),"close mdcrd")   
   if (ntwv > 0) call checkNCerror(nf90_close(vel_ncid),"close mdvel")

#else
   write(6,*) 'No binary trajectory support in this version'
   write(6,*) 'reconfigure using the -bintraj flag'
   call mexit(6,1)
#endif

end subroutine close_binary_files
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit coordinates or velocities, r(istart:n), to  netcdf.
#ifdef BINTRAJ
subroutine write_binary_traj(r,istart,n,unit)
   use constantph, only : target_ph
   use file_io_dat
   use netcdf
   use nblist, only: a,b,c,alpha,beta,gamma
#  ifdef MPI
      use remd, only: rem,mytargettemp
#  endif 
   
   implicit none
#  include "box.h"
   integer, intent(in) :: istart,n,unit
   _REAL_, intent(in) ::  r(n)

   select case (unit)
      case (MDCRD_UNIT)
         if (n == 3 .and. ntb > 0) then       ! Assume this is box (fails on one atom systems)
            call checkNCerror(nf90_put_var(crd_ncid,Cell_lengthVarID,(/ a,b,c/), &
                  start = (/ 1, crd_frame /), count = (/ 3, 1 /)), 'write cell lengths')
            call checkNCerror(nf90_put_var(crd_ncid,Cell_angleVarID,(/ alpha,beta,gamma /), &
                  start = (/ 1, crd_frame /), count = (/ 3, 1 /)), 'write cell angles')
         else
            call checkNCerror(nf90_put_var(crd_ncid,CoordVarID, r(istart:n), &
                  start = (/ 1, 1, crd_frame /),count = (/ 3, (n-istart+1)/3, 1 /)), &
                  'write atom coords')
#  ifdef MPI
            ! If this is a replica run write temp0
            if (rem>0.and.rem<3) then
               call checkNCerror(nf90_put_var(crd_ncid,TempVarID,mytargettemp, &
                     start = (/ crd_frame /) ), &
                     'write replica mytargettemp')
            else if (rem==4) then
               call checkNCerror(nf90_put_var(crd_ncid,TempVarID,target_ph, &
                     start = (/ crd_frame /) ), &
                     'write replica pH')
            endif
#  endif
         end if
      case (MDVEL_UNIT)
         call checkNCerror(nf90_put_var(vel_ncid,VelocVarID, r(istart:n), &
               start = (/ 1, 1, vel_frame /), count = (/ 3, (n-istart+1)/3, 1 /)), &
               'write velocities')
      case default
         write (6,*) 'Error: unhandled unit ',unit,' selected for output in bintraj'
   end select

end subroutine write_binary_traj

#else
subroutine write_binary_traj(r,istart,n,unit)
   integer, intent(in) :: istart,n,unit
   _REAL_, intent(in) ::  r(n)
   write(6,*) 'No binary trajectory support in this version'
   write(6,*) 'reconfigure using the -bintraj flag'
   call mexit(6,1)
end subroutine write_binary_traj
#endif

!-----------------------------------------------------------------------
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write scalar data and increment frame counter
subroutine end_binary_frame(unit)
#ifdef BINTRAJ
   use file_io_dat
   use netcdf
   implicit none
#  include "md.h"
   integer, intent(in) :: unit
 select case (unit)
      case (MDCRD_UNIT)
         call checkNCerror(nf90_put_var(crd_ncid, crd_TimeVarID, &
               (/ t /), start = (/ crd_frame /), count = (/ 1 /)), 'write time')
         
         call checkNCerror(nf90_sync(crd_ncid))

         crd_frame = crd_frame + 1
         ! Sync frame counters if combined trajectory
         if (ntwv == -1) vel_frame = vel_frame + 1
      case (MDVEL_UNIT)
         call checkNCerror(nf90_put_var(vel_ncid, vel_TimeVarID, &
               (/ t /), start = (/ vel_frame /), count = (/ 1 /)), 'write time')
         
         call checkNCerror(nf90_sync(vel_ncid))
         
         vel_frame = vel_frame +1
      case default
         write (6,*) 'Error: unhandled unit ',unit,' selected for end frame in bintraj'
   end select

#else
   integer, intent(in) :: unit
   write(6,*) 'No binary trajectory support in this version'
   write(6,*) 'reconfigure using the -bintraj flag'
   call mexit(6,1)
#endif

end subroutine end_binary_frame
!-----------------------------------------------------------------------


!! Private subroutines follow

#ifdef BINTRAJ
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open new or existing NetCDF file
subroutine open_nc_file(ncid,filename,owrite,facc,natom,title,exists,frame,TimeVarID)
   use netcdf
   use file_io_dat, only : MAX_FN_LEN

   implicit none

   integer, intent(out) :: ncid
   character(len=MAX_FN_LEN), intent(in) :: filename
   character(len=80), intent(in)  :: title
   character, intent(in) :: owrite, facc
   integer, intent(in) :: natom
   logical, intent(out) :: exists
   integer, intent(out) :: frame, TimeVarID

   integer :: cmode, status
   character(len=nf90_max_name) :: name

   ! Create netCDF file
   cmode = nf90_64bit_offset
   exists = .false.
   if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
   if (owrite == 'U' .and. facc == 'A') cmode = ior(cmode, nf90_noclobber)
   status = nf90_create(path=filename,cmode=cmode,ncid=ncid)
   frame = 1
   if (status == nf90_eexist) then
      exists = .true.
      status = nf90_open(path = filename, mode = nf90_write, ncid = ncid)
      if (status == nf90_noerr) then !else fall through to error handling below and exit
         call checkNCerror(nf90_inq_dimid(ncid, "frame",FrameDimID), "query Frame ID")
         call checkNCerror(nf90_Inquire_Dimension(ncid,FrameDimID,name,frame),"query current frame")
         frame=frame + 1
         call checkNCerror(nf90_inq_varid(ncid, "time", TimeVarID), &
               "query time ID")         
         return 
      end if
   end if
   call checkNCerror(status, "create file")
   if (status /= nf90_noerr) then
#ifndef DUMB
      write (0,*) 'Error on opening ', filename
#endif
      call mexit(6,1)
   end if

   ! Define dimensions
   call checkNCerror(nf90_def_dim(ncid,"frame", nf90_unlimited, FrameDimID))
   call checkNCerror(nf90_def_dim(ncid,"spatial", 3, SpatialDimID))
   call checkNCerror(nf90_def_dim(ncid,"atom", natom, AtomDimID))
   call checkNCerror(nf90_def_dim(ncid,"label", NCLABELLEN, LabelDimID))

   ! Set global attributes
   call checkNCerror(nf90_put_att(ncid,nf90_global, "title", &
         title), "define title")
   call checkNCerror(nf90_put_att(ncid,nf90_global, "application", &
         'AMBER'), "define application")
   call checkNCerror(nf90_put_att(ncid,nf90_global, "program", &
         'sander'), "define program")
   call checkNCerror(nf90_put_att(ncid,nf90_global, "programVersion", &
         '10.0'), "define programVersion")
   call checkNCerror(nf90_put_att(ncid,nf90_global, "Conventions", &
         'AMBER'), "define Convention")
   call checkNCerror(nf90_put_att(ncid,nf90_global, "ConventionVersion", &
         '1.0'), "define ConventionVersion")

   ! Define non-optional variables
   call checkNCerror(nf90_def_var(ncid, "spatial", nf90_char, &
         (/  SpatialDimID /), SpatialVarID))
   
   call checkNCerror(nf90_def_var(ncid, "time", nf90_float, &
         (/  FrameDimID /), TimeVarID))
   call checkNCerror(nf90_put_att(ncid, TimeVarID, "units", &
         "picosecond"), "define time units")

end subroutine open_nc_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ End define mode and write label variables for specified nc file
subroutine end_nc_define(ncid)
   use netcdf

   implicit none

   integer, intent(in) :: ncid
   
   ! Set NoFill and end definition mode
   call checkNCerror(nf90_set_fill(ncid, nf90_nofill, oldMode), "set no fill")
   call checkNCerror(nf90_enddef(ncid), "end define")

   ! Fill dimension label variables
   call checkNCerror(nf90_put_var(ncid, SpatialVarID, &
         (/ 'x','y','z' /), start = (/ 1 /), count = (/ 3 /)), "write spatial variable")

end subroutine end_nc_define
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check NetCDF error status and format error message if appropriate
subroutine checkNCerror(status, location)
  use netcdf

  implicit none

  integer, intent(in) :: status !Status returned by netCDF lib call
  character(*), optional, intent(in) :: location !purpose of call
  
  if(status /= nf90_noerr) then
     write (6,*) 'NetCDF error: ',trim(nf90_strerror(status))
     if (present(location)) then
        write (6,*) '  at ',location
     end if
  end if

end subroutine checkNCerror
#endif

end module bintraj
