#include "copyright.i"

!*******************************************************************************
!
! Module: pmemd_lib_mod
!
! Description: Contains useful wrappers for common tasks
!              
!*******************************************************************************

module pmemd_lib_mod

use file_io_dat_mod

  implicit none

#ifdef MPI
  ! This needs to be visible for adding groups; sigh...
  integer, save         :: world_group
#endif /* MPI */

contains

!*******************************************************************************
!
! Subroutine:  mexit
!
! Description: Exit wrapper routine; cleans up MPI stuff as it quits
!              
!*******************************************************************************

subroutine mexit(ifil, i)

#ifdef MPI
  implicit none
  include 'mpif.h'
#endif

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

  if (i .eq. 0) then
    call mpi_finalize(err_ret_code)
  else
    call mpi_abort(mpi_comm_world, i, err_ret_code)
  end if
#endif /* MPI */

  if (i .eq. 0) then
    stop
  else
    stop 'PMEMD Terminated Abnormally!'
  end if

end subroutine mexit

!*******************************************************************************
!
! Subroutine:  alloc_error
!
! Description: Error subroutine to call if an allocation failed
!              
!*******************************************************************************

subroutine alloc_error(routine, string1)

  implicit none
  
  character(*) routine, string1
    
  write(mdout, '(1x,2a)') &
    'FATAL dynamic memory allocation error in subroutine ', routine
  
  write(mdout, '(1x,a)') string1
  
  call mexit(6, 1)
  
end subroutine alloc_error

!*******************************************************************************
!
! Subroutine:  setup_alloc_error
!
! Description: Error subroutine to call if global allocation failed
!              
!*******************************************************************************

subroutine setup_alloc_error

  implicit none
  
  write(mdout, '(1x,2a)') 'FATAL global dynamic memory setup allocation error!'
 
  call mexit(6, 1)

end subroutine setup_alloc_error

!*******************************************************************************
!
! Subroutine: get_num_tokens
!
! Description: Determines how many whitespace-delimited words are in a string
!
!*******************************************************************************

subroutine get_num_tokens(string, token_num)

  implicit none

! Passed arguments

  character(*), intent(in) :: string

  integer, intent(out)     :: token_num

! Local variables

  integer :: string_loc  ! our location in the string
  integer :: iend        ! last non-whitespace character location

  string_loc = 1
  iend = len_trim(string)
  token_num = 0

  do while (string_loc .le. iend)

    if ( string(string_loc:string_loc) .le. ' ' ) then
      string_loc = string_loc + 1
    else

      do while ( string(string_loc:string_loc) .gt. ' ' )
        string_loc = string_loc + 1
      end do

      token_num = token_num + 1
    end if
  end do

end subroutine get_num_tokens

!*******************************************************************************
!
! Subroutine: get_token
!
! Description: Gets the num'th token in a string
!
!*******************************************************************************

subroutine get_token(string, num, token)

  implicit none

! Passed arguments
  
  character(*), intent(in)  :: string  ! The string to parse
  character(*), intent(out) :: token   ! The token to return

  integer, intent(in)       :: num     ! Which token to return

! Local variables

  integer   :: num_tokens
  integer   :: istart
  integer   :: iend
  integer   :: string_loc
  integer   :: token_count

  ! Uncomment the below chunk of code for a "safe" get_num_tokens at the 
  ! expense of calling get_num_tokens() each time a specific token is
  ! pulled from the string. When it's commented out, token will just be
  ! a blank string upon return

! call get_num_tokens(string, num_tokens)

! if (num .gt. num_tokens)
!   write(mdout, *) ' Error in get_token: Looking for more tokens than &
!                     &there are in string'
!   call mexit(6,1)
! end if

  ! Now get the num'th token

  token_count = 0
  istart = 1
  iend = len_trim(string)
  token = ' '

  do while (istart .le. iend)

    if (string(istart:istart) .le. ' ') then

      istart = istart + 1
      
    else

      do string_loc = istart, iend
        if ( string(string_loc:string_loc) .le. ' ' ) exit
      end do

      token_count = token_count + 1

      ! If this is the token we want, store it and return
      if ( token_count .eq. num ) then
        token = string(istart:string_loc-1)
        return
      end if

      istart = string_loc ! Move to the next token

    end if

  end do

end subroutine get_token

!*******************************************************************************
!
! Subroutine: check_inpcrd_overflow
!
! Description: Checks a line from the inpcrd file to see if it has *s, then
!              prints a helpful warning if it does
!
!*******************************************************************************

subroutine check_inpcrd_overflow(line, periodic)

   implicit none

   ! Passed variables
   character(*), intent(in) :: line     ! line buffer that offends the parser
   logical, intent(in)      :: periodic ! Does our simulation have PBCs?

   ! Local variables
   integer      :: i ! counter

   
   ! My advice to you depends on whether you're running a simulation with
   ! PBCs or in implicit solvent/vacuum
   do i = 1, len_trim(line)
      if (line(i:i) .eq. '*' .and. periodic) then
         write(6, '(a)') '*s in the inpcrd file often indicate an overflow of &
            &the Fortran format used', &
            'to store coordinates in the inpcrd/restart files. &
            &This often happens when', &
            'coordinates are not wrapped into the center cell &
            &(when iwrap = 0) and some', &
            'particles diffuse too far away. Try restarting from your &
            &last good restart', &
            'file and setting iwrap=1 or using a NetCDF restart file &
            &format. See the', &
            'Amber manual for details'
         ! Only print error message once. Bail out here.
         exit
      else if (line(i:i) .eq. '*' .and. .not. periodic) then
         write(6, '(a)') '*s in the inpcrd file often indicate an overflow of &
            &the Fortran format used', &
            'to store coordinates in the inpcrd/restart files. &
            &This often happens when', &
            'particles diffuse very far away from each other. Make &
            &sure you are removing', &
            'center-of-mass translation (nscm /= 0) or check if you &
            &have multiple, mobile', &
            'molecules that have diffused very far away from each other. &
            &This condition is', &
            'highly unusual for non-periodic simulations.'
         ! Only print error message once. Bail out here.
         exit
      end if
   end do

   return

end subroutine check_inpcrd_overflow

!*******************************************************************************
!
! Subroutine: basename
!
! Description: Operates like the GNU basename utility -- it strips off the
!              leading path of a file name
!
!*******************************************************************************

subroutine basename(filename)

   implicit none

   ! Passed variables
   character(*), intent(in out)  :: filename

   ! Local variables
   integer :: i, pos

   pos = 1
   do i = 1, len_trim(filename)
      if (filename(i:i) .eq. '/') pos = i + 1
   end do

   filename = filename(pos:len_trim(filename))

   return

end subroutine basename

!*******************************************************************************
!
! Function: upcase
!
! Description: Converts a string to uppercase
!
!*******************************************************************************

function upcase(string) result (upper)

   implicit none

   character(len=*), intent(in) :: string
   character(len=len_trim(string)) :: upper
   integer :: i, ic

   do i = 1, len_trim(string)
      ic = iachar(string(i:i))
      if (ic .gt. 96 .and. ic .lt. 123) ic = ic - 32
      upper(i:i) = achar(ic)
   end do

end function upcase

!*******************************************************************************
!
! Subroutine: get_atomic_number
!
! Description:
!  Determines the atomic number based on the mass and the first letter of the
!  atom name which has been read in from the topology file. The assumption here
!  is that an atom name matches the first letter of the element. If it doesn't,
!  then this subroutine will need to be modified. However, topology files
!  starting with Amber 12 should include an ATOMIC_NUMBER section in it, so this
!  routine shouldn't be required in too many cases. Same code as what's used in
!  sqm to do the same task.
!  
!*******************************************************************************

subroutine get_atomic_number(atom_name, atom_mass, atomic_number, errorFlag)
    
  implicit none
    
! Passed variables

  character(len=4), intent(in)   :: atom_name
  double precision, intent(in)   :: atom_mass
  logical, optional, intent(out) :: errorFlag
  integer,          intent(out)  :: atomic_number

! Local variables

  logical ::localErrorFlag = .false.

  ! Lanthanides are not supported.
  ! Actinides are not supported.
    
  if( upcase(atom_name(1:1)) .eq. 'A' ) then
    if(atom_mass > 24.0d0 .and. atom_mass <= 28.0d0) then
      atomic_number =  13 !Aluminium
    else if(atom_mass > 35.0d0 .and. atom_mass <= 40.0d0) then
      atomic_number =  18 !Argon
    else if(atom_mass > 73.0d0 .and. atom_mass <= 77.0d0) then
      atomic_number =  33 !Arsenic
    else if(atom_mass > 106.0d0 .and. atom_mass <= 109.0d0) then
      atomic_number =  47 !Silver
    else if(atom_mass > 195.0d0 .and. atom_mass <= 199.0d0) then
      atomic_number =  79 !Gold
    else if(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
      atomic_number =  85 !Astatine
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'B' ) then
    if(atom_mass > 8.0d0 .and. atom_mass <= 10.0d0) then
      atomic_number =  4 !Beryllium
    else if(atom_mass > 10.0d0 .and. atom_mass <= 12.0d0) then
      atomic_number =  5 !Boron
    else if(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
      atomic_number =  35 !Bromine
    else if(atom_mass > 135.0d0 .and. atom_mass <= 139.0d0) then
      atomic_number =  56 !Barium
    else if(atom_mass > 207.0d0 .and. atom_mass <= 211.0d0) then
      atomic_number =  83 !Bismuth
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'C' ) then
    if(atom_mass > 10.0d0 .and. atom_mass <= 14.0d0) then
      atomic_number =  6 !Carbon
    else if(atom_mass > 33.0d0 .and. atom_mass <= 37.0d0) then
      atomic_number =  17 !Chlorine
    else if(atom_mass > 38.0d0 .and. atom_mass <= 42.0d0) then
      atomic_number =  20 !Calcium
    else if(atom_mass > 50.0d0 .and. atom_mass <= 54.0d0) then
      atomic_number =  24 !Chromium
    else if(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
      atomic_number =  27 !Cobalt
    else if(atom_mass > 61.0d0 .and. atom_mass <= 65.0d0) then
      atomic_number =  29 !Copper
    else if(atom_mass > 110.0d0 .and. atom_mass <= 114.0d0) then
      atomic_number =  48 !Cadmium
    else if(atom_mass > 131.0d0 .and. atom_mass <= 135.0d0) then
      atomic_number =  55 !Cesium
    else
      localErrorFlag=.true.
    end if
  
  else if( upcase(atom_name(1:1)) .eq. 'D' ) then
    localErrorFlag=.true.

  else if( upcase(atom_name(1:1)) .eq. 'E' ) then
    localErrorFlag=.true.

  else if( upcase(atom_name(1:1)) .eq. 'F' ) then
    if(atom_mass > 17.0d0 .and. atom_mass <= 21.0d0) then
      atomic_number =  9 !Fluorine
    else if(atom_mass > 54.0d0 .and. atom_mass <= 58.0d0) then
      atomic_number =  26 !Iron
    else if(atom_mass > 218.0d0 .and. atom_mass <= 228.0d0) then
      atomic_number =  87 !Francium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'G' ) then
    if(atom_mass > 67.0d0 .and. atom_mass <= 71.0d0) then
      atomic_number =  31 !Gallium
    else if(atom_mass > 71.0d0 .and. atom_mass <= 75.0d0) then
      atomic_number =  32 !Germanium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'H' ) then
    if(atom_mass > 0.0d0 .and. atom_mass <= 2.0d0) then
      atomic_number =  1 !Hydrogen
    else if(atom_mass > 3.0d0 .and. atom_mass <= 5.0d0) then
      atomic_number =  2 !Helium
    else if(atom_mass > 176.0d0 .and. atom_mass <= 180.0d0) then
      atomic_number =  72 !Hafnium
    else if(atom_mass > 198.0d0 .and. atom_mass <= 202.0d0) then
      atomic_number =  80 !Mercury
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'I' ) then
    if(atom_mass > 112.0d0 .and. atom_mass <= 116.0d0) then
      atomic_number = 49 !Indium
    else if(atom_mass > 125.0d0 .and. atom_mass <= 129.0d0) then
      atomic_number =  53 !Iodine
    else if(atom_mass > 190.0d0 .and. atom_mass <= 194.0d0) then
      atomic_number =  77 !Iridium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'J' ) then
    localErrorFlag=.true.

  else if( upcase(atom_name(1:1)) .eq. 'K' ) then
    if(atom_mass > 37.0d0 .and. atom_mass <= 41.0d0) then
      atomic_number = 19 !Potassium
    else if(atom_mass > 77.0d0 .and. atom_mass <= 86.0d0) then
      atomic_number = 36 !Krypton
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'L') then
    if(atom_mass > 6.0d0 .and. atom_mass <= 8.0d0) then
      atomic_number = 3 !Lithium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'M' ) then
    if(atom_mass > 22.0d0 .and. atom_mass <= 26.0d0) then
      atomic_number = 12 !Magnesium
    else if(atom_mass > 53.0d0 .and. atom_mass <= 57.0d0) then
      atomic_number = 25 !Manganese
    else if(atom_mass > 94.0d0 .and. atom_mass <= 98.0d0) then
      atomic_number = 42 !Molybdenem
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'N') then
    if(atom_mass > 13.0d0 .and. atom_mass <= 15.0d0) then
      atomic_number = 7 !Nitrogen
    else if(atom_mass > 19.0d0 .and. atom_mass <= 22.0d0) then
      atomic_number = 10 !Neon
    else if(atom_mass > 22.1d0 .and. atom_mass <= 23.0d0) then
      atomic_number = 11 !Sodium
    else if(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
      atomic_number = 28 !Nickel
    else if(atom_mass > 95.0d0 .and. atom_mass <= 99.0d0) then
      atomic_number = 41 !Niobium
    else
        localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'O' ) then
    if(atom_mass > 14.0d0 .and. atom_mass <= 18.0d0) then
      atomic_number = 8 !Oxygen
    else if(atom_mass > 188.0d0 .and. atom_mass <= 192.0d0) then
      atomic_number = 76 !Osmium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'P' ) then
    if(atom_mass > 29.0d0 .and. atom_mass <= 33.0d0) then
      atomic_number = 15 !Phosphorus
    else if(atom_mass > 104.0d0 .and. atom_mass <= 108.0d0) then
      atomic_number = 46 !Palladium
    else if(atom_mass > 193.0d0 .and. atom_mass <= 197.0d0) then
      atomic_number = 78 !Platinum
    else if(atom_mass > 205.0d0 .and. atom_mass <= 208.0d0) then
      atomic_number = 82 !Lead
    else if(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
      atomic_number = 84 !Polonium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'Q' ) then
    localErrorFlag=.true.

  else if( upcase(atom_name(1:1)) .eq. 'R' ) then
    if(atom_mass > 84.0d0 .and. atom_mass <= 88.0d0) then
      atomic_number = 37 !Rubidium
    else if(atom_mass > 99.0d0 .and. atom_mass <= 102.0d0) then
      atomic_number = 44 !Ruthenium
    else if(atom_mass > 102.0d0 .and. atom_mass <= 105.0d0) then
      atomic_number = 45 !Rhodium
    else if(atom_mass > 184.0d0 .and. atom_mass <= 188.0d0) then
      atomic_number = 75 !Rhenium
    else if(atom_mass > 210.0d0 .and. atom_mass <= 222.5d0) then
      atomic_number = 86 !Radon
    else if(atom_mass > 223.0d0 .and. atom_mass <= 229.0d0) then
      atomic_number = 88 !Radium
    else
      localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'S' ) then
    if(atom_mass > 26.0d0 .and. atom_mass <= 30.0d0) then
       atomic_number = 14 !Silicon
    else if(atom_mass > 30.0d0 .and. atom_mass <= 34.0d0) then
       atomic_number = 16 !Sulphur
    else if(atom_mass > 43.0d0 .and. atom_mass <= 47.0d0) then
       atomic_number = 21 !Scandium
    else if(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
       atomic_number = 34 !Selenium
    else if(atom_mass > 86.0d0 .and. atom_mass <= 89.0d0) then
       atomic_number = 38 !Strontium
    else if(atom_mass > 116.0d0 .and. atom_mass <= 120.0d0) then
       atomic_number = 50 !Tin
    else if(atom_mass > 120.0d0 .and. atom_mass <= 124.0d0) then
       atomic_number = 51 !Antimony
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'T' ) then
    if(atom_mass > 46.0d0 .and. atom_mass <= 50.0d0) then
       atomic_number = 22 !Titanium
    else if(atom_mass > 96.0d0 .and. atom_mass <= 100.0d0) then
       atomic_number = 43 !Technetium
    else if(atom_mass > 125.0d0 .and. atom_mass <= 130.0d0) then
       atomic_number = 52 !Tellurium
    else if(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
       atomic_number = 73 !Tantalum
    else if(atom_mass > 201.0d0 .and. atom_mass <= 206.0d0) then
       atomic_number = 81 !Thallium
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'U' ) then
    write(6,*) 'Unable to correctly identify element ', atom_name
    call mexit(6,1)

  else if( upcase(atom_name(1:1)) .eq. 'V' ) then
    if(atom_mass > 49.0d0 .and. atom_mass <= 53.0d0) then
       atomic_number = 23 !Vanadium
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'W' ) then
    if(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
       atomic_number = 74 !Tungsten
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'X' ) then
    if (atom_mass > 123.0d0 .and. atom_mass < 136.0d0) then
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'Z' ) then
    if(atom_mass > 61.0d0 .and. atom_mass <= 69.0d0) then
       atomic_number = 30 !Zinc
    else if(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
       atomic_number = 40 !Zirconium
    else
       localErrorFlag=.true.
    end if

  else if( upcase(atom_name(1:1)) .eq. 'Z' ) then
    if(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
       atomic_number = 40 !Zirconium
    else
      localErrorFlag=.true.
    end if

  else
    localErrorFlag=.true.
  end if
  
  if (present(errorFlag)) then
    errorFlag = localErrorFlag
  else if (localErrorFlag) then
    write(6,*) 'Unable to correctly identify element ', atom_name
    call mexit(6, 1)
  end if

end subroutine get_atomic_number

#ifdef MPI
!*******************************************************************************
!
! Subroutine: bcast_logical
!
! Description: Broadcasts a single logical variable across a given communicator
!
!*******************************************************************************

subroutine bcast_logical(to_bcast, comm)

  implicit none
  include 'mpif.h'

! Passed variables
  logical, intent(in out) :: to_bcast
  integer, intent(in)     :: comm

! Local variables
  integer err_code_mpi ! Define a local copy here

  call mpi_bcast(to_bcast, 1, mpi_logical, 0, comm, err_code_mpi)

  return

end subroutine bcast_logical

#endif /* MPI */

end module pmemd_lib_mod
