#include "copyright.i"

!*******************************************************************************
!
! Module:  gbl_constants_mod
!
! Description: A central repository of constants that otherwise create
!              circular dependency headaches.
!              
!*******************************************************************************

module gbl_constants_mod

  implicit none

! Global constants:

  double precision, parameter   :: PI = 3.1415926535897932384626433832795d0

  ! big_int is largest int that fits in an i8 field, for max nscm, etc.
  integer, parameter    :: big_int = 99999999

  integer, parameter    :: RETIRED_INPUT_OPTION = -10301        ! from sander8
  integer, parameter    :: UNSUPPORTED_INPUT_OPTION = -10302

  character(11), parameter      :: info_hdr =       '| INFO:    '
  character(11), parameter      :: warn_hdr =       '| WARNING: '
  character(11), parameter      :: error_hdr =      '| ERROR:   '
  character(11), parameter      :: extra_line_hdr = '|          '
  character(15), parameter      :: prog_name =      'PMEMD.AMOEBA 12'
  character(41), parameter      :: use_sander = &
                                   '|           Please use sander 12 instead.'

end module gbl_constants_mod
