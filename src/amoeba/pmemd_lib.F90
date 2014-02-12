#include "copyright.i"

!*******************************************************************************
!
! Module: pmemd_lib_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pmemd_lib_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  alloc_error
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_error(routine, string1)

  use parallel_dat_mod
  
  implicit none
  
  character(*) routine, string1
    
  write(mdout, '(1x,2a)') &
    'FATAL dynamic memory allocation error in subroutine ', routine
  
  write(mdout, '(1x,a)') string1
  
  call mexit(6, 1)
  
end subroutine alloc_error

!*******************************************************************************
!
! Subroutine:  
!
! Description: setup_alloc_error
!              
!*******************************************************************************

subroutine setup_alloc_error

  use parallel_dat_mod

  implicit none
  
  write(mdout, '(1x,2a)') 'FATAL global dynamic memory setup allocation error!'
 
  call mexit(6, 1)

end subroutine setup_alloc_error

end module pmemd_lib_mod
