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

  type(dihed_rec), parameter &
                        :: null_dihed_rec = dihed_rec(0, 0, 0, 0, 0)

  type listdata_rec
    integer             :: offset
    integer             :: cnt
  end type listdata_rec

  integer, parameter    :: listdata_rec_ints = 2    ! don't use for allocation!

  type atm_lst_rec
    integer             :: idx
    integer             :: nxt
  end type atm_lst_rec

  integer, parameter    :: atm_lst_rec_ints = 2    ! don't use for allocation!

!CHARMM Urey-Bradley Terms 
  type angle_ub_rec
    integer             :: atm_i
    integer             :: atm_k
    integer             :: parm_idx
  end type angle_ub_rec

  integer, parameter    :: angle_ub_rec_ints  =  3

  type dihed_imp_rec
    integer             :: atm_i
    integer             :: atm_j
    integer             :: atm_k
    integer             :: atm_l
    integer             :: parm_idx
  end type dihed_imp_rec

  integer, parameter    :: dihed_imp_rec_ints = 5    ! don't use for allocation!

  type cmap_rec
    integer             :: atm_i
    integer             :: atm_j
    integer             :: atm_k
    integer             :: atm_l
    integer             :: atm_m
    integer             :: parm_idx
  end type cmap_rec

  integer, parameter    :: cmap_rec_ints = 6    ! don't use for allocation!

  ! Get a 1 byte integer with integer(byte) :: foo.
  integer, parameter :: byte = selected_int_kind(2)
  ! Get a 2 byte integer with integer(short) :: foo.
  integer, parameter :: short = selected_int_kind(4)

end module gbl_datatypes_mod
