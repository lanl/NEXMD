#include "copyright.i"

!*******************************************************************************
!
! Module:  gbl_datatypes_mod
!
! Description: A central repository of datatypes that otherwise create
!              circular dependency headaches.
!              
!*******************************************************************************

module gbl_datatypes_mod

  implicit none

! Global types:

! These types don't exist per se in prmtop; we convert a bunch of bond, angle
! and dihedral information in arrays into arrays of records of these types:

! The atm_[ijkl] values are now actual atom indices; they previously were
! indices into 1D arrays.

  type bond_rec
    integer             :: atm_i
    integer             :: atm_j
    integer             :: parm_idx
  end type bond_rec

  integer, parameter    :: bond_rec_ints = 3    ! don't use for allocation!

  type angle_rec
    integer             :: atm_i
    integer             :: atm_j
    integer             :: atm_k
    integer             :: parm_idx
  end type angle_rec

  integer, parameter    :: angle_rec_ints = 4    ! don't use for allocation!

  type dihed_rec
    integer             :: atm_i
    integer             :: atm_j
    integer             :: atm_k
    integer             :: atm_l
    integer             :: parm_idx
  end type dihed_rec

  integer, parameter    :: dihed_rec_ints = 5    ! don't use for allocation!

  type maskdata_rec
    integer             :: maskptr
    integer             :: nummask
  end type maskdata_rec

end module gbl_datatypes_mod
