#include "dprec.fh"

!------------------------------------------------------------------------------
module cmd_vars

! internal variables used in CMD
!------------------------------------------------------------------------------

  implicit none

  save

  logical :: activate = .true.
  logical :: restart_cmd = .false.
  logical :: eq_cmd = .false.

  integer :: file_nrg_cmd = 2001
  integer :: file_pos_cmd  = 2002
  integer :: file_vel_cmd  = 2003

  integer :: nstep_cmd

  _REAL_ :: t_cmd, temp_cmd, etot_cmd, eke_cmd
  _REAL_ :: adiab_param = 1.d0
  _REAL_, allocatable :: lambda_nmode( : )
  _REAL_, allocatable :: omega_nmode( : )
  _REAL_, allocatable :: mass_nmode( : )
  _REAL_, allocatable :: nmode_to_cart( :, : )
  _REAL_, allocatable :: cart_to_nmode( :, : )
  _REAL_, allocatable :: fcart_to_fnmode( :, : )

end module
