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
  double precision, parameter   :: DEG_TO_RAD = PI / 180.d0
  double precision, parameter   :: RAD_TO_DEG = 180.d0 / PI

  ! big_int is largest int that fits in an i8 field, for max nscm, etc.
  integer, parameter    :: big_int = 99999999

  integer, parameter    :: RETIRED_INPUT_OPTION = -10301        ! from sander8
  integer, parameter    :: UNSUPPORTED_INPUT_OPTION = -10302
  integer, parameter    :: NO_INPUT_VALUE = 12344321

  character(11), parameter      :: info_hdr =       '| INFO:    '
  character(11), parameter      :: warn_hdr =       '| WARNING: '
  character(11), parameter      :: error_hdr =      '| ERROR:   '
  character(11), parameter      :: extra_line_hdr = '|          '
  character(8), parameter       :: prog_name =      'PMEMD 12'
  character(12), parameter      :: VERSION = 'Version 12.0'
  character(41), parameter      :: use_sander = &
                                   '|           Please use sander 12 instead.'

end module gbl_constants_mod
