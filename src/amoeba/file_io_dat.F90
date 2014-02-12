!*******************************************************************************
!
! Module: file_io_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module file_io_dat_mod

  implicit none

! File-related data.  No need to broadcast this stuff, since only the master
! does i/o:

  integer, parameter            :: max_fn_len = 80
  character(max_fn_len), save   :: mdin_name
  character(max_fn_len), save   :: mdout_name
  character(max_fn_len), save   :: mdinfo_name
  character(max_fn_len), save   :: prmtop_name
  character(max_fn_len), save   :: inpcrd_name
  character(max_fn_len), save   :: refc_name
  character(max_fn_len), save   :: mdcrd_name
  character(max_fn_len), save   :: mdvel_name
  character(max_fn_len), save   :: mden_name
  character(max_fn_len), save   :: restrt_name
  character(max_fn_len), save   :: logfile_name
  character(max_fn_len), save   :: infile_suffix
  character(max_fn_len), save   :: outfile_suffix
  character, save       :: owrite

! Logical unit numbers for pmemd files.

  integer, parameter    :: mdin      =  5
  integer, parameter    :: mdout     =  6
  integer, parameter    :: mdinfo    =  7
  integer, parameter    :: prmtop    =  8
  integer, parameter    :: inpcrd    =  9
  integer, parameter    :: refc      = 10
  ! who was 11?
  integer, parameter    :: mdcrd     = 12
  integer, parameter    :: mdvel     = 13
  ! who was 14?
  integer, parameter    :: mden      = 15
  integer, parameter    :: restrt    = 16
  integer, parameter    :: restrt2   = 17
  integer, parameter    :: logfile   = 18

#ifdef AMOEBA
! Enumerations for crdfile_type:

  integer, parameter                    :: crdfile_type_undefined = 0
  integer, parameter                    :: crdfile_type_old = 1
  integer, parameter                    :: crdfile_type_new = 2
#endif /* AMOEBA */

end module file_io_dat_mod
