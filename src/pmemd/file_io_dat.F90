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
#ifdef MPI
  integer, parameter            :: grpfl_line_len = 4096
#endif

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
!AMD file for reweighting values
  character(max_fn_len), save   :: amdlog_name = 'amd.log'
#ifdef MPI
! Set file names for multipmemd and REMD. proc_map refers
! to the file describing how many processors to assign to
! each communicator
  character(max_fn_len), save   :: groupfile_name = ' '
  character(max_fn_len), save   :: remlog_name = 'rem.log'
  character(max_fn_len), save   :: remtype_name = 'rem.type' ! not supported yet
  character(max_fn_len), save   :: proc_map_name
  character(max_fn_len), save   :: remd_partners_name = ' '
  character(grpfl_line_len), save :: groupline_buffer = ' '
  character(grpfl_line_len), save :: proc_map_buffer = ' '
#endif
  character, save       :: owrite

! Logical unit numbers for pmemd files.

  integer, parameter    :: mdin          =  5
  integer, parameter    :: mdout         =  6
  integer, parameter    :: mdinfo        =  7
  integer, parameter    :: prmtop        =  8
  integer, parameter    :: inpcrd        =  9
  integer, parameter    :: refc          = 10
  ! who was 11?
  integer, parameter    :: mdcrd         = 12
  integer, parameter    :: mdvel         = 13
  ! who was 14?
  integer, parameter    :: mden          = 15
  integer, parameter    :: restrt        = 16
  integer, parameter    :: restrt2       = 17
  integer, parameter    :: logfile       = 18
#ifdef MPI
  integer, parameter    :: groupfile     = 19
  integer, parameter    :: remlog        = 20
  integer, parameter    :: proc_map      = 21
  integer, parameter    :: remd_partners = 22

! Maybe not quite the right place to put this, but it can't be in
! multipmemd_mod, since we'll get a cyclic dependency.

  logical, save         :: ng_nonsequential = .false.
#endif
  integer, parameter    :: amdlog        = 23

end module file_io_dat_mod
