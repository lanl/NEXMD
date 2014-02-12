#include "dprec.fh"
module qm2_extern_module
! ----------------------------------------------------------------
! Interface for QM and QM/MM MD with SANDER and external programs
!
! Currently supports:
! QM           with ADF and GAMESS
! QM and QM/MM with TeraChem, Gaussian and Orca
!
! Initial implementation for ADF by
! Matthew Clark and Prithvi Undavalli
! (SDSC LSSI summer highschool students)
! under supervision of
! Andreas Goetz and Ross Walker (SDSC)
! 
! Date: August 2010
!
! Extensions by Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! ----------------------------------------------------------------

  use qmmm_module, only: qmmm_nml
  use qm2_extern_adf_module, only: get_adf_forces
  use qm2_extern_gms_module, only: get_gms_forces
  use qm2_extern_tc_module, only:  get_tc_forces
  use qm2_extern_gau_module, only: get_gau_forces
  use qm2_extern_orc_module, only:  get_orc_forces
  use qm2_extern_nw_module, only:  get_nw_forces

  implicit none

  private
  public :: qm2_extern_get_qm_forces, qm2_extern_finalize
  
  contains

  subroutine qm2_extern_get_qm_forces(nstep, nqmatoms, qmcoords, nclatoms, clcoords, &
       escf, dxyzqm, dxyzcl)

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, CODATA08_AU_TO_DEBYE

    use pimd_vars, only: ipimd

    use full_pimd_vars, only: mybeadid
    use file_io_dat
#if defined(MPI)
    use remd, only : rem
#   include "parallel.h"
#endif /* MPI */
#   include "md.h" ! We are only using this for values of nstlim and maxcyc

    integer, intent(in) :: nstep
    integer, intent(in) :: nqmatoms             ! Number QM of atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM gradient
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM gradient

    ! List of supported external programs; will search for namelists in this order
    character(len=20), save :: extern_program
    character(len=3) :: id
    logical, save :: first_call = .true.
    logical, save :: do_gradient = .true.
    integer :: charge

    if(first_call) then
      ! first_call = .false. ! SET TO .false. BELOW (check_charge)
      ! If doing post-processing
      if( nstlim==0 .or. maxcyc==0 ) do_gradient = .false. 
      call select_program(extern_program)
      if (nclatoms > 0) then
        ! if electrostatic embedding in use
        if ( qmmm_nml%qmmm_int /= 5 ) then
           ! QM/MM not possible with ADF or GAMESS at present
           if ( (extern_program == 'adf') .or. (extern_program == 'gms') ) then
             call sander_bomb("qm2_extern_get_qm_forces","nquant /= natom", &
             trim(extern_program)//" does not support QM/MM")
           end if
        end if
      end if
      ! print constants that are used for conversion
      write (6,'(3(/,a),/)') ' Constants for unit conversion taken from', &
                             ' Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730', &
                             ' and using the thermochemical calorie (1 cal = 4.184 J):'
      write (6, '(a, es19.12)') ' A_TO_BOHRS  = ', CODATA08_A_TO_BOHRS
      write (6, '(a, es17.10)') ' AU_TO_KCAL  = ', CODATA08_AU_TO_KCAL
      write (6, '(a, es15.8)')  ' AU_TO_DEBYE = ', CODATA08_AU_TO_DEBYE
    endif

    ! Determine id
    id = ''
    if ( ipimd > 0 ) then
       ! Add number in case of parallel PIMD runs
       write (id,'(i3.3)') mybeadid
#if defined(MPI)
    else if ((rem > 0) .or. (qmmm_nml%vsolv > 1)) then
            ! Add rank for parallel REMD run
       write (id,'(i3.3)') masterrank
#endif /* MPI */
    end if

    ! Call chosen program
    select case (extern_program)
    case('adf')
      call get_adf_forces(do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        escf, dxyzqm, charge)
    case('gau')
      call get_gau_forces(do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)
    case('gms')
      call get_gms_forces(do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        escf, dxyzqm, charge)
    case('tc')
      call get_tc_forces( do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)
    case('orc')
      call get_orc_forces( do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)
    case('nw')
      call get_nw_forces( do_gradient, nstep, ntpr, id, nqmatoms, qmcoords,&
        nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)
    case default
      call sander_bomb("qm2_extern_get_qm_forces","External namelist not found", &
        "Please check your input.")
    end select

    ! We need to make sure that the charge specified for the external program
    ! is the same as the charge specified in the &qmmm namelist
    ! Otherwise the charge adjustment in qmmm_adjust_q() does the wrong thing
    if (first_call) then
      first_call = .false.
      ! This only matters for electronic embedding
      if ( (nclatoms > 0) .and. (qmmm_nml%qmmm_int /= 5) ) then
         call check_charge(extern_program, qmmm_nml%qmcharge, charge)
      end if
   end if

  end subroutine qm2_extern_get_qm_forces

  subroutine select_program(extern_program)

    integer :: i, ifind
    character(len=20) :: programs(6) = (/'adf', 'gms', 'tc ', 'gau', 'orc', 'nw '/)
    character(len=20), intent(out) :: extern_program

    ! Select which external program to use
    extern_program='none'
    do i=1, size(programs)
      rewind 5
      call nmlsrc(programs(i),5,ifind)
      if(ifind>0) then
        extern_program=programs(i)
        exit
      endif
    end do

 end subroutine select_program

 ! At present, both the &qmmm and &extern_program namelists
 ! need to specify the charge (to adjust MM charges correctly for charged QM systems)
 subroutine check_charge(extern_program, qmcharge, charge)

   implicit none
   character(len=20), intent(in) :: extern_program
   integer, intent(in) :: qmcharge, charge

   if (qmcharge /= charge) then
      call sander_bomb('check_charge() (qm2_extern_module)', &
           'Charge in &qmmm and &'//trim(extern_program)//' inconsistent', &
           'Please check your input and specify charge in both namelists.')
   end if
   
 end subroutine check_charge

 ! Used for tasks that must be run after the last step
 subroutine qm2_extern_finalize()
   use qm2_extern_tc_module, only: tc_finalize

   character(len=20) :: extern_program
   call select_program(extern_program)

   select case (extern_program)
     case('tc')
       call tc_finalize()
   end select

 end subroutine qm2_extern_finalize


end module qm2_extern_module
