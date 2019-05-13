#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! barf [ba:rf]  2. "He suggested using FORTRAN, and everybody barfed."
!    - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)
!

module ncsu_sander_proxy

use file_io_dat

implicit none

private

public :: sander_mdin_name
public :: sander_mdin_unit

public :: multisander_rem
public :: multisander_initremd
public :: multisander_numgroup

public :: multisander_numwatkeep

public :: terminate
public :: is_master

public :: proxy_finalize

public :: sander_imin
public :: sander_natoms
public :: sander_mdtime
public :: sander_sgft
public :: sander_tempsg
public :: sander_temp0
public :: sander_timestep
public :: sander_init
public :: sander_nstlim
public :: sander_ntp
public :: sander_ntb
public :: sander_nsolut 

#ifdef MPI
public :: set_sander_temp0
#endif

public :: sander_atom_name
public :: remember_atom_names

character(len = 4), private, pointer, save :: atom_names(:) => null()

public :: flush_UNIT

#ifndef NCSU_DISABLE_ASSERT
private :: afailed
#endif /* NCSU_DISABLE_ASSERT */

public :: remember_rem
public :: remember_initremd

integer, private, save :: saved_rem = -3212341
#ifdef MPI
logical, private, save :: saved_initremd = .true.
#else
logical, private, save :: saved_initremd = .false.
#endif /* MPI */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine terminate()
   implicit none
   call mexit(6, 1)
end subroutine terminate

!-----------------------------------------------------------------------------

character(len=MAX_FN_LEN) function sander_mdin_name()
   implicit none
#  include "files.h"
   sander_mdin_name = mdin
end function sander_mdin_name

!-----------------------------------------------------------------------------

pure integer function sander_mdin_unit()
   implicit none
   sander_mdin_unit = 5
end function sander_mdin_unit

!-----------------------------------------------------------------------------

pure integer function multisander_numgroup()
   implicit none
#  include "files.h"
   multisander_numgroup = numgroup
end function multisander_numgroup

!-----------------------------------------------------------------------------

NCSU_PURE_EXCEPT_ASSERT integer function multisander_rem()
   implicit none
   multisander_rem = saved_rem
end function multisander_rem

!-----------------------------------------------------------------------------

subroutine remember_rem(r)
   implicit none
   integer, intent(in) :: r
   saved_rem = r
end subroutine remember_rem

!-----------------------------------------------------------------------------

subroutine remember_initremd(i)
   implicit none
   logical, intent(in) :: i
   saved_initremd = i
end subroutine remember_initremd

!-----------------------------------------------------------------------------

pure logical function multisander_initremd()
   implicit none
   multisander_initremd = saved_initremd
end function multisander_initremd

!-----------------------------------------------------------------------------

pure integer function multisander_numwatkeep()
   implicit none
#  include "md.h"
   multisander_numwatkeep = numwatkeep
end function multisander_numwatkeep

!-----------------------------------------------------------------------------

!
! for use in ncsu_assert() & co [where unneeded indirection is acceptable]
!

logical function is_master()

   implicit none

#ifdef MPI
#  include "ncsu-mpi.h"
   ncsu_assert(commsander /= mpi_comm_null)
   is_master = (sanderrank == 0)
#else
   is_master = .true.
#endif /* MPI */

end function is_master

!-----------------------------------------------------------------------------

subroutine proxy_finalize()
   implicit none

   if (associated(atom_names)) &
      deallocate(atom_names)
end subroutine proxy_finalize

!-----------------------------------------------------------------------------

pure integer function sander_imin()
   implicit none
#include "md.h"
   sander_imin = imin
end function sander_imin

!-----------------------------------------------------------------------------
pure integer function sander_nsolut()
   implicit none
#include "md.h"
#include "memory.h"
   sander_nsolut = natom-(nres-ibgwat+1)*4
end function sander_nsolut
!-----------------------------------------------------------------------------

pure integer function sander_natoms()
   implicit none
#include "memory.h"
   sander_natoms = natom
end function sander_natoms

!-----------------------------------------------------------------------------

pure NCSU_REAL function sander_mdtime()
   implicit none
#include "md.h"
   sander_mdtime = t
end function sander_mdtime

!-----------------------------------------------------------------------------

pure NCSU_REAL function sander_sgft()
   use sgld, only:sgft
   implicit none
#include "md.h"
   sander_sgft = sgft
end function sander_sgft

!-----------------------------------------------------------------------------

pure NCSU_REAL function sander_tempsg()
   use sgld, only:tempsg
   implicit none
#include "md.h"
   sander_tempsg = tempsg
end function sander_tempsg

!-----------------------------------------------------------------------------

pure NCSU_REAL function sander_temp0()
   implicit none
#include "md.h"
   sander_temp0 = temp0
end function sander_temp0

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine set_sander_temp0(new_temp0)
   implicit none
   NCSU_REAL, intent(in) :: new_temp0
#include "md.h"
   temp0 = new_temp0
end subroutine set_sander_temp0
#endif /* MPI */

!-----------------------------------------------------------------------------

pure NCSU_REAL function sander_timestep()
   implicit none
#include "md.h"
   sander_timestep = dt
end function sander_timestep

!-----------------------------------------------------------------------------

pure integer function sander_init()
   implicit none
#include "md.h"
   sander_init = init
end function sander_init

!-----------------------------------------------------------------------------

pure integer function sander_nstlim()
   implicit none
#include "md.h"
   sander_nstlim = nstlim
end function sander_nstlim

!-----------------------------------------------------------------------------

pure integer function sander_ntp()
   implicit none
#include "md.h"
   sander_ntp = ntp
end function sander_ntp

!-----------------------------------------------------------------------------

pure integer function sander_ntb()
   implicit none
#include "box.h"
   sander_ntb = ntb
end function sander_ntb

!-----------------------------------------------------------------------------

character(len = 4) function sander_atom_name(n)

   implicit none

   integer, intent(in) :: n

#  include "memory.h"

   ncsu_assert(n > 0)
   ncsu_assert(n <= sander_natoms())
   ncsu_assert(associated(atom_names))

   sander_atom_name = atom_names(n)

end function sander_atom_name

!-----------------------------------------------------------------------------

subroutine remember_atom_names(ih)

   use ncsu_constants

   implicit none

   character(len = 4), intent(in) :: ih(*)

#include "memory.h"

   integer :: n, error

   if (associated(atom_names)) &
      deallocate(atom_names)

   ncsu_assert(natom > 0)

   allocate(atom_names(natom), stat = error)
   if (error /= 0) then
      write (unit = ERR_UNIT, fmt = '(a,a)') &
         NCSU_ERROR, 'out of memory in remember_atom_names()'
      call terminate()
   end if

   do n = 1, natom
      atom_names(n) = ih(m04 + n - 1)
   end do

end subroutine remember_atom_names

!-----------------------------------------------------------------------------

subroutine flush_UNIT(lun)

   implicit none

   integer, intent(in) :: lun

   call amflsh(lun)

end subroutine flush_UNIT

!-----------------------------------------------------------------------------

#ifndef NCSU_DISABLE_ASSERT
subroutine afailed(filename, lineno)

   use ncsu_constants, only : ERR_UNIT

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write (unit = ERR_UNIT, fmt = '(/a,a,a,i3,a/)') &
      NCSU_ERROR, filename, ':', lineno, ': ncsu_assert() failed'

   call terminate()

end subroutine afailed
#endif /* NCSU_DISABLE_ASSERT */


end module ncsu_sander_proxy
