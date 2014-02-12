module qmmm_qmtheorymodule
! ----------------------------------------------------------------------
! PURPOSE: Data type defining the theory we are running
! 
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : February 2010
! ----------------------------------------------------------------------
  
  implicit none

  private
  public :: qmTheoryType
  public :: Set, String
#ifdef MPI
  public :: Broadcast
#endif
  public :: CheckRetiredQmTheoryInputOption ! This is legacy

  type qmTheoryType
     logical PM3
     logical AM1
     logical AM1D           ! This is AM1/d (with d orbitals), *NOT* AM1-D (with dispersion correction)
     logical MNDO           
     logical MNDOD          ! This is MNDO/d
     logical PDDGPM3
     logical PDDGMNDO
     logical PM3CARB1
     logical PM3ZNB
     logical DFTB
     logical RM1
     logical PDDGPM3_08
     logical PM6
     logical DISPERSION     ! This is the D  correction (dispersion)
     logical DISPERSION_HYDROGENPLUS ! This is the DH+ correction (dispersion and hydrogen bond)
     logical PM3MAIS
     logical EXTERN         !External interface to ADF/GAMESS/Gaussian/TeraChem
  end type qmTheoryType
  
  interface Set
     module procedure SetQmTheoryType
  end interface

  interface String
     module procedure QmTheoryString
  end interface

#ifdef MPI
  interface Broadcast
     module procedure BroadcastQmTheoryType
  end interface
#endif

contains

  ! --------------------------
  ! Check retired input option
  ! --------------------------
  subroutine CheckRetiredQmTheoryInputOption(qmtheory)

    use constants, only : RETIRED_INPUT_OPTION
   implicit none
    integer, intent(in) :: qmtheory

    if ( qmtheory /= RETIRED_INPUT_OPTION ) then
       call sander_bomb('read_qmmm_nm_and_alloc','qmtheory is specified in the qmmm namelist.', &
                           'It is deprecated, please use qm_theory.')
    end if
    
  end subroutine CheckRetiredQmTheoryInputOption

  ! ---------------------------------------------------------
  ! Set the QM theory type
  ! NOTE: We default to PM3 if none is given in the input !!!
  ! ---------------------------------------------------------
  subroutine SetQmTheoryType(self, qm_theory)

    use UtilitiesModule, only : Upcase
    
    implicit none
    type (qmTheoryType), intent(out) :: self
    character(len=12), intent(in) :: qm_theory
    
    self%PM3        = .false. 
    self%AM1        = .false.
    self%AM1D       = .false.
    self%MNDO       = .false.
    self%MNDOD      = .false.
    self%PDDGPM3    = .false.
    self%PDDGMNDO   = .false.
    self%PM3CARB1   = .false.
    self%PM3ZNB     = .false.
    self%DFTB       = .false.
    self%RM1        = .false.
    self%PDDGPM3_08 = .false.
    self%PM6        = .false.
    self%DISPERSION = .false.
    self%DISPERSION_HYDROGENPLUS = .false.
    self%PM3MAIS    = .false.
    self%EXTERN     = .false.
    
    select case (Upcase(qm_theory))
    case ('', 'PM3')
       self%PM3 = .true.
    case ('AM1')
       self%AM1 = .true.
    case ('AM1_D*', 'AM1-D*')
       self%AM1 = .true.
       self%dispersion = .true.
    case ('AM1_DH+', 'AM1-DH+')
       self%AM1 = .true.
       self%DISPERSION_HYDROGENPLUS = .true.
    case ('AM1D', 'AM1_D', 'AM1/D')
       self%AM1D = .true.
    case ('MNDO')
       self%MNDO = .true.
    case ('MNDOD', 'MNDO_D', 'MNDO/D')
       self%MNDOD = .true.
    case ('PM3-PDDG', 'PM3PDDG', 'PM3_PDDG', 'PDDG-PM3', 'PDDGPM3', 'PDDG_PM3')
       self%PDDGPM3 = .true.
    case ('MNDO-PDDG', 'MNDOPDDG', 'MNDO_PDDG', 'PDDG-MNDO', 'PDDGMNDO', 'PDDG_MNDO')
       self%PDDGMNDO = .true.
    case ('PM3-CARB1', 'PM3CARB1', 'PM3_CARB1')
       self%PM3CARB1 = .true.
    case ('PM3ZNB', 'PM3-ZNB', 'PM3_ZNB', 'PM3/ZNB', 'ZNB')
       self%PM3ZNB = .true.
    case ('DFTB', 'SCCDFTB', 'SCC-DFTB', 'SCC_DFTB')
       self%DFTB = .true.
    case ('RM1')
       self%RM1 = .true.
    case ('PM3-PDDG08', 'PM3PDDG08', 'PM3_PDDG08', 'PDDG-PM308', 'PDDGPM308', 'PDDG_PM308', &
         'PM3-PDDG-08', 'PDDGPM3_08', 'PDDG_PM3_08', 'PM3-PDDG_08')
       self%PDDGPM3_08 = .true.
    case ('PM6')
       self%PM6 = .true.
    case ('PM6_D', 'PM6-D')
       self%PM6   = .true.
       self%DISPERSION = .true.
    case ('PM6_DH+', 'PM6-DH+')
       self%PM6       = .true.
       self%DISPERSION_HYDROGENPLUS = .true.
    case ('PM3-MAIS', 'PM3MAIS', 'PM3_MAIS', 'MAIS')
       self%PM3MAIS = .true.
    case ('EXTERN')
       self%EXTERN = .true.
    case default
       call sander_bomb('qmtheorymodule:SetQmTheoryType','Unknown method specified for qm_theory', &
            'Valid options are: PM3, AM1, RM1, MNDO, PM3-PDDG, PM3-PDDG_08, MNDO-PDDG, PM3-CARB1, '&
            //'PM3-ZNB, PM3-MAIS, PM6, DFTB and EXTERN (external)')
    end select

  end subroutine SetQmTheoryType

  ! ---------------------------------
  ! Return string with QM theory type
  ! ---------------------------------
  function QmTheoryString(self)

    implicit none
    type(qmTheoryType), intent(in) :: self
    character(len=12)   :: qmTheoryString
  
    if (self%PM3) then
       qmTheoryString = 'PM3'
    else if (self%AM1) then
       qmTheoryString = 'AM1'
    else if (self%AM1D) then
       qmTheoryString = 'AM1/D'
    else if (self%MNDO) then
       qmTheoryString = 'MNDO'
    else if (self%MNDOD) then
       qmTheoryString = 'MNDO/D'
    else if (self%PDDGPM3) then
       qmTheoryString = 'PDDG/PM3'
    else if (self%PDDGMNDO) then
       qmTheoryString = 'PDDG/MNDO'
    else if (self%PM3CARB1) then
       qmTheoryString = 'PM3/CARB1'
    else if (self%PM3ZNB) then
       qmTheoryString = 'PM3/ZNB'
    else if (self%DFTB) then
       qmTheoryString = 'DFTB'
    else if (self%RM1) then
       qmTheoryString = 'RM1'
    else if (self%PDDGPM3_08) then
       qmTheoryString = 'PDDG/PM3(08)'
    else if (self%PM6) then
       qmTheoryString = 'PM6'
    else if (self%PM3MAIS) then
       qmTheoryString = 'PM3-MAIS'
    else if (self%EXTERN) then
       qmTheoryString = 'EXTERN'
    end if

  end function QmTheoryString

#ifdef MPI
  subroutine BroadcastQmTheoryType(self)

    implicit none
#include "parallel.h"
    include 'mpif.h'

    type(qmTheoryType) :: self
    integer :: ier

    call mpi_bcast(self%PM3       , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%AM1       , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%AM1D      , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%MNDO      , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%MNDOD     , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PDDGPM3   , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PDDGMNDO  , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PM3CARB1  , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PM3ZNB    , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%DFTB      , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%RM1       , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PDDGPM3_08, 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PM6       , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%PM3MAIS   , 1, mpi_logical, 0, commsander, ier)
    call mpi_bcast(self%EXTERN    , 1, mpi_logical, 0, commsander, ier)

  end subroutine BroadcastQmTheoryType
#endif
   
end module qmmm_qmtheorymodule
