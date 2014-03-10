!i <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
module qmmm_nml_module
! ----------------------------------------------------------------------
! PURPOSE: Data type holding qmmm namelist data
! 
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : October 2010
!
! Based on previous data type and subroutines contained in qmmm_module
! as written by Ross Walker and Mike Crowley
!
! Methods:
!   new      : allocate
!   delete   : deallocate
!   broadcast: broadcast data
!              (includes allocation in case arrays do not exist on slaves)
!   print    : print content of all variables to stdout
!              (useful for debugging)
! ----------------------------------------------------------------------

  use qmmm_qmtheorymodule, only : qmTheoryType

  implicit none

  private
  
  public :: qmmm_nml_type
#ifdef MPI
  public :: broadcast
#endif
  public :: new, delete, print

  type qmmm_nml_type

     ! REMEMBER TO MODIFY qmmm_mpi_setup below if you modify this structure
     ! Contains QMMM namelist variables - all cpus should have this info available.

     ! Cutoff in angstroms for QM-MM electrostatics - default = same as MM-MM cutoff
     _REAL_ :: qmcut     

     ! Cutoff^2
     _REAL_ :: qmcut2    

     ! Distance in angstroms between QM atom and link atom
     ! A value of <0.0 means the link atom gets placed on the MM link pair atom's coordinates on every MD step.
     _REAL_ :: lnk_dis   

     ! SCF Convergence threshold for SCF routine
     ! Default = 1.0D-8, min (tightest conv) = 1.0D-12, max = 1.0D0
     _REAL_ :: scfconv

     ! SCF Convergence threshold for density matrix.
     ! If tight_p_conv = true then = scfconv else it is 0.5*sqrt(scf_conv)
     _REAL_ :: density_conv 
  
     !+TJG 1/26/2010
     ! Convergence criteria for maximum value in the error matrix FP-PF
     _REAL_ :: errconv 

     ! Maximum number of matrices used in a diis extrapolation
     integer :: ndiis_matrices      

     ! The initial number of diis tokens... the number of iterations
     ! that diis extrapolations will be attempted before giving up on diis
     integer :: ndiis_attempts      

     ! NOTE:
     ! One can turn off diis by 
     !  (1) ndiis_matrices=1 
     !    or 
     !  (2) ndiis_attempts = 0
     !-TJG 1/26/2010
  
     ! Criteria for the maximum change in the density matrix between two 
     ! SCF iterations before allowing pseudo diagonalisations.
     ! Default = 0.05.
     _REAL_  :: pseudo_diag_criteria           
     
     ! Atomic number of link atom.
     integer :: lnk_atomic_no 

     ! This defines how classical valence terms that cross the QM/MM boundary are dealt with.
     ! 1 = (Default) in this case any bond, angle or dihedral that involves at least one
     !    MM atom, including the MM link pair atom is included. This means the following
     !    where QM = QM atom, MM = MM atom, MML = MM link pair atom.
     !    Bonds = MM-MM, MM-MML, MML-QM
     !    Angles = MM-MM-MM, MM-MM-MML, MM-MML-QM, MML-QM-QM
     !    Dihedrals = MM-MM-MM-MM, MM-MM-MM-MML, MM-MM-MM-MML-QM, MM-MML-QM-QM, MML-QM-QM-QM
     ! 2 = Only include valence terms that include a full MM atom. I.e count the MM link
     !    pair atom as effectively being a QM atom. This therefore gives
     !    Bonds = MM-MM, MM-MML
     !    Angles = MM-MM-MM, MM-MM-MML, MM-MML-QM
     !    Dihedrals = MM-MM-MM-MM, MM-MM-MM-MML, MM-MM-MML-QM, MM-MML-QM-QM
     integer :: lnk_method 

     ! GB method used in QMMM (when igb /= 0 and /= 6)
     ! 0 = (default) leave QM charges at 0 effectively doing gas phase QM (EGB(MM-MM only). 
     ! 1 = use resp charges from prmtop for qmmm charges.
     ! 2 = use mulliken charges for QM atoms consistent with GB field, modified fock matrix.
     ! 3 = (FOR_DEBUG_ ONLY) use mulliken charges from gas phase QMMM calc, not consistent
     !     with GB field and so gradients will not be accurate.
     integer :: qmgb

     ! QM theory selected, see qmmm_qmtheorymodule
     type(qmTheoryType) :: qmtheory 

     ! Charge of the QM region in electron units - must be an integer charge. (Default = 0)
     integer :: qmcharge  

     ! Spin state - default = 1 (singlet). Current Valid values = 1 (alternate options not currently available).
     integer :: spin

     ! Controls amount of info printed about qm part of calc - (Default = 0)
     integer :: verbosity

     ! Maximum number of SCF cycles to conduct before assuming convergence has failed. (Default = 1000).
     integer :: itrmax

     ! Whether to shake qm atoms if ntc>1 (default = 1 - shake QM atoms) 0 = do not shake.
     integer :: qmshake

     ! Used for qmewald - Maximum K space vectors
     integer :: kmaxqx, kmaxqy, kmaxqz 

     ! Used for qmewald - Maximum K squared values for spherical cutoff in k space.
     integer :: ksqmaxq

     ! Used for qmewald - Kappa                    
     _REAL_ :: kappa

     ! Flag for doing an ewald sum for periodic QM-QM interactions - default = 1
     integer :: qm_ewald

     ! Flag controlling the way QM-MM interactions are handled:
     !   0 = mechanical embedding. No electrostatic QM-MM is calculated. Only interactions
     !       come from VDW and bonds, angles and dihedrals that cross the boundary.
     !   1 = Full QM-MM interaction is calculated by placing a 1S gaussian orbital on the MM
     !       atom and calculating the full multipole interaction (Default)
     !   2 = As 1 but also include extra Gaussian core-core terms if AM1, PM3, RM1, PM6 etc or derivative
     !       type hamiltonians are in used.
     !   3 = Along with qm_theory=PM3 will invoke the modified QM-MM core-charge interaction term
     integer :: qmmm_int

     ! Flag to turn on and off the QM/MM switching function in the code
     logical :: qmmm_switch

     ! Lower bound of the QM/MM switching function
     _REAL_ :: r_switch_lo

     ! Upper bound of the QM/MM switching function
     _REAL_ :: r_switch_hi

     ! Flag for whether to adjust the charges of certain MM atoms so that charge is conserved.
     !   0 = No charge adjustment
     !   1 = Missing charge is distributed over the nearest atom to each MM link pair.
     !   2 = Missing charge is distributed evenly over all MM atoms - excluding MM link pairs. (default)
     integer :: adjust_q
  
     ! Flag controlling which diagonalization routine to use when doing full diag in SCF.
     !   0 = Automatically pick fastest routine.
     !   1 = Use internal diagonalization routine. (default)
     !   2 = Use lapack dspev.
     !   3 = Use lapack dspevd.
     !   4 = Use lapack dspevx.
     !   5 = Use lapack dsyev.
     !   6 = Use lapack dsyevd.
     !   7 = Use lapack dsyevr. (not currently implemented)
     !   8 = Use lapack dsyevx. (not currently implemented)
     integer :: diag_routine

     ! Flag controlling the way in which the density matrix for MD step t+dt is predicted based on previous MD steps.
     !   0 = Use final converged density matrix from previous MD step. (default)
     !   1 = Use predictive algorithm based on Phys. Rev. Lett., 2006, 97, 123001
     !       Here Pguess(t) = 2Pconv(t-dt) - Pguess(t-2dt)
     !   2 = xl-bomd, see xlbomd module for more information! Added by JAB (LANL, 2014)
     integer :: density_predict

     ! Are we using variable solvent or adaptive QM/MM?
     ! vsolv = 0 (default, no variable solvent or adQMMM)
     ! vsolv = 1 (simple variable solvent, uses &vsolv namelist)
     ! vsolv = 2 (proper adaptive QMMM with variable solvent, uses both &vsolv and &adqmmm namelists)
     integer :: vsolv
  
     ! Flag controlling the way in which the fock matrix for MD step t+dt is predicted based on previous MD steps.
     !   0 = (Default) Do not attempt to predict the Fock matrix
     !   1 = Use predictive algorithm based on a modification of Chem. Phys. Lett., 2004, 386, 272
     !       by Ross Walker and Gustavo Seabra to only extrapolate the electronic part of
     !       the Fock matrix. Incompatible with density_predict > 0.
     integer :: fock_predict

     ! Prefactor to multiply each final stored fock matrix by when building the new predicted fock matrix.
     _REAL_ :: fockp_d1
     _REAL_ :: fockp_d2
     _REAL_ :: fockp_d3
     _REAL_ :: fockp_d4   
     
     !sah added for divcon/MOPAC
     integer :: idc       

     !added for DivPB
     integer :: divpb     

     ! Number of quantum atoms as specified in the input qmmm namelist
     ! Does not contain link atoms !!!
     ! The number of quantum atoms including link atoms is contained in qmmm_struct_type%nquant (qmmm_struct_module)
     integer :: nquant
  
     ! List of atom numbers of the qm atoms as numbered in prmtop (nquant)
     ! These are only the quantum atoms as specified in the input qmmm namelist
     ! Does not contain link atoms !!!
     ! A list of all quantum atoms including link atoms is contained in qmmm_struct_type (qmmm_struct_module)
     integer, pointer :: iqmatoms(:) => null()

#ifdef OPENMP
     ! Maximum openmp threads to use inside QMMM routines.
     ! If diag_routine /= 0 then this value will be used as the argument
     ! to omp_set_num_threads for all threaded QMMM functions.
     ! If diag_routine = 0 then the code will test the performance from
     ! 1 thread to this number to find the optimum value to use.
     integer :: qmmm_omp_max_threads
#endif
   
     ! GMS: Charge scaling factor for QM-FEP
     _REAL_  :: chg_lambda

     !! DFTB Options
     ! Maximum number of SCC iterations before resetting Broyden (default: 70 )
     integer :: dftb_maxiter    

     ! Use dispersion?  (default: 0 = false)
     integer :: dftb_disper

     ! DFTB Charge output: (default: 0 = Mulliken, 1 = CM3)
     integer :: dftb_chg

     ! Flag to control the printing of dipole. (default: 0 = not print, 1 = QM, 2 = QM + MM) 
     integer :: printdipole

     ! Electronic temperature, in Kelvins. Default: 0.0K
     _REAL_  :: dftb_telec

     ! Step size for automatic telec increse, for convergence. Default = 0.0K
     _REAL_  :: dftb_telec_step

     ! 3rd order SCC-DFTB (default: 0 == No third order)
     !     'PA': Do 3rd order, Proton Affinities parameterization
     !     'PR':               Phosphate reactions parameterization
     !     'READ': read the parameters from a user-specified file (TO IMPLEMENT)
     character(Len=256) :: dftb_3rd_order  
     !! End DFT Options
  
     ! Flag for whether QMMM is on (True) or off (False)
     logical :: ifqnt

     ! Flag for analytical vs (pseudo) numerical qm-qm derivates in qm2 routines. (Default = true, analytical).
     logical :: qmqm_analyt

! CML: This was added in Amber Tools 11 for eventual support of excited state analytical derivatives
! CML: 7/10/12

     ! Flag for analytical vs (pseudo) numerical qm-qm derivates for excited state. (Default = false, pseudo-numerical).
     logical :: qmqm_exc_analyt

     ! Flag to control convergence criteria for density matrix. If 0 SCF routine will converge energy
     ! to SCFCRT (def = 1*10^-8) and density to 0.05 * Sqrt(SCFCONV)
     ! If 1 then both energy and density will be converged to SCFCONV.
     logical :: tight_p_conv 

     ! Flag to control the printing of mulliken, cm1a and cm2a charges. If set to true (1) then
     ! charges are calculated and printed on every step. Otherwise no charges are printed.
     logical :: printcharges 
     logical :: printbondorders !Flag to control the printing of bond orders at the end of the calculation

     ! Default = .false. - add correction to peptide linkages?
     logical :: peptide_corr

     !Store qmqm 1-electron repulsion integrals in memory?
     logical :: qmqm_erep_incore
     
     ! Whether or not to allow pseudo diagonalisations in SCF when possible.
     logical :: allow_pseudo_diag

     ! Flag to store qmmm rij and related equations in memory - default = true.
     logical :: qmmmrij_incore

     ! if true then a crude pdb of the qm system will be written to file qmmm_region.pdb
     logical :: writepdb

     ! Flag for using PME instead of regular ewald for QM-MM interactions when qm_ewald>0.
     ! Default = true     
     logical :: qm_pme
  
  end type qmmm_nml_type

  interface new
     module procedure new_qmmm_nml_type
  end interface

  interface delete
     module procedure delete_qmmm_nml_type
  end interface

#ifdef MPI
  interface broadcast
     module procedure broadcast_qmmm_nml_type
  end interface
#endif

  interface print
     module procedure print_qmmm_nml_type
  end interface


contains


subroutine new_qmmm_nml_type(self)

  ! allocates iqmatoms array
  ! at this point self%nquant needs to be known !

  implicit none

  ! needs to be inout since it could already contain data which we don't want to destroy
  type(qmmm_nml_type), intent(inout) :: self

  integer :: nquant, ier

  nquant = self%nquant

  if ( .not. associated(self%iqmatoms) ) then
     allocate ( self%iqmatoms(nquant), stat=ier )
     REQUIRE(ier == 0)
  end if

end subroutine new_qmmm_nml_type


subroutine delete_qmmm_nml_type(self)

  implicit none

  type(qmmm_nml_type), intent(inout) :: self

  integer :: ier

  if ( associated(self%iqmatoms) ) then
     deallocate ( self%iqmatoms, stat = ier)
     REQUIRE(ier == 0)
  end if

end subroutine delete_qmmm_nml_type


#ifdef MPI
  subroutine broadcast_qmmm_nml_type(self)

    use qmmm_qmtheorymodule, only : qmTheoryType, Broadcast

    implicit none
#include "parallel.h"
    include 'mpif.h'

    type(qmmm_nml_type), intent(inout) :: self
    integer :: ier

     call broadcast(self%qmtheory)
     !I don't think we can assume the structure is linear in memory so we have to broadcast each seperately
     call mpi_bcast(self%qmcut,                1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
     call mpi_bcast(self%qmcut2,               1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
     call mpi_bcast(self%lnk_dis,              1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%scfconv,              1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%density_conv,         1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     !+TJG 01/26/2010
     call mpi_bcast(self%errconv,              1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%ndiis_matrices,       1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%ndiis_attempts,       1, mpi_integer,          0, commsander, ier) 
     !-TJG 01/26/2010
     call mpi_bcast(self%pseudo_diag_criteria, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%lnk_atomic_no,        1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%lnk_method,           1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qmgb,                 1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qmcharge,             1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%spin,                 1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%verbosity,            1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%itrmax,               1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%adjust_q,             1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%diag_routine,         1, mpi_integer,          0, commsander, ier) 
#ifdef OPENMP
     call mpi_bcast(self%qmmm_omp_max_threads, 1, mpi_integer,          0, commsander, ier) 
#endif
     call mpi_bcast(self%density_predict,      1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%fock_predict,         1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%fockp_d1,             1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%fockp_d2,             1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%fockp_d3,             1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%fockp_d4,             1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%vsolv,                1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%idc,                  1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qmshake,              1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qm_ewald,             1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%kmaxqx,               1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%kmaxqy,               1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%kmaxqz,               1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%ksqmaxq,              1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%kappa,                1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
     call mpi_bcast(self%qmmm_int,             1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qmqm_analyt,          1, mpi_logical,          0, commsander, ier)
     call mpi_bcast(self%qmqm_exc_analyt,      1, mpi_logical,          0, commsander, ier) ! CML 7/10/12
     call mpi_bcast(self%tight_p_conv,         1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%printcharges,         1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%printdipole,          1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%peptide_corr,         1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%allow_pseudo_diag,    1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%qmqm_erep_incore,     1, mpi_logical,          0, commsander, ier)
     call mpi_bcast(self%qmmmrij_incore,       1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%qm_pme,               1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%writepdb,             1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%nquant,               1, mpi_integer,          0, commsander, ier) 
     call mpi_bcast(self%qmmm_switch,          1, mpi_logical,          0, commsander, ier) 
     call mpi_bcast(self%r_switch_lo,          1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
     call mpi_bcast(self%r_switch_hi,          1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  !! GMS
     call mpi_bcast(self%chg_lambda ,      1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  !! DFTB options
     call mpi_bcast(self%dftb_maxiter    , 1, mpi_integer         , 0, commsander, ier) 
     call mpi_bcast(self%dftb_disper     , 1, mpi_integer         , 0, commsander, ier) 
     call mpi_bcast(self%dftb_chg        , 1, mpi_integer         , 0, commsander, ier) 
     call mpi_bcast(self%dftb_telec      , 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
     call mpi_bcast(self%dftb_telec_step , 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  !!

     ! arrays may not have been allocated on slaves yet, therefore try to create object
     call new(self)
     call mpi_bcast(self%iqmatoms, self%nquant, mpi_integer, 0, commsander, ier)

  end subroutine broadcast_qmmm_nml_type
#endif


  subroutine print_qmmm_nml_type(self)

    use qmmm_qmtheorymodule, only : String

    type(qmmm_nml_type), intent(in) :: self

    integer :: i, j, jstart, jend
    integer, parameter :: jstep = 10
    character(len=30) :: tmpstr

    write(6,'(/a)') '*** qmmm namelist'
    write(6,'(a,e10.2)') 'qmcut                      = ', self%qmcut
    write(6,'(a,e10.2)') 'qmcut2                     = ', self%qmcut2    
    write(6,'(a,e10.2)') 'lnk_dis                    = ', self%lnk_dis
    write(6,'(a,e10.2)') 'scfconv                    = ', self%scfconv
    write(6,'(a,e10.2)') 'density_conv               = ', self%density_conv 
    write(6,'(a,e10.2)') 'errconv                    = ', self%errconv 
    write(6,'(a,i5)')    'ndiis_matrices             = ', self%ndiis_matrices      
    write(6,'(a,i5)')    'ndiis_attempts             = ', self%ndiis_attempts      
    write(6,'(a,e10.2)') 'pseudo_diag_criteria       = ', self%pseudo_diag_criteria           
    write(6,'(a,i5)')    'lnk_atomic_no              = ', self%lnk_atomic_no
    write(6,'(a,i5)')    'lnk_method                 = ', self%lnk_method
    write(6,'(a,i5)')    'qmgb                       = ', self%qmgb
    write(6,'(2a)')      'qmtheory                   = ', String(self%qmtheory)
    write(6,'(a,i5)')    'qmcharge                   = ', self%qmcharge
    write(6,'(a,i5)')    'spin                       = ', self%spin
    write(6,'(a,i5)')    'verbosity                  = ', self%verbosity
    write(6,'(a,i5)')    'itrmax                     = ', self%itrmax
    write(6,'(a,i5)')    'qmshake                    = ', self%qmshake
    write(6,'(a,i5)')    'kmaxqx                     = ', self%kmaxqx
    write(6,'(a,i5)')    'kmaxqy                     = ', self%kmaxqy
    write(6,'(a,i5)')    'kmaxqz                     = ', self%kmaxqz
    write(6,'(a,i5)')    'ksqmaxq                    = ', self%ksqmaxq
    write(6,'(a,i5)')    'qm_ewald                   = ', self%qm_ewald
    write(6,'(a,i5)')    'qmmm_int                   = ', self%qmmm_int
    write(6,'(a,i5)')    'adjust_q                   = ', self%adjust_q
    write(6,'(a,i5)')    'diag_routine               = ', self%diag_routine
    write(6,'(a,i5)')    'density_predict            = ', self%density_predict
    write(6,'(a,i5)')    'vsolv                      = ', self%vsolv
    write(6,'(a,i5)')    'fock_predict               = ', self%fock_predict
    write(6,'(a,e10.2)') 'fockp_d1                   = ', self%fockp_d1
    write(6,'(a,e10.2)') 'fockp_d2                   = ', self%fockp_d2
    write(6,'(a,e10.2)') 'fockp_d3                   = ', self%fockp_d3
    write(6,'(a,e10.2)') 'fockp_d4                   = ', self%fockp_d4
    write(6,'(a,i5)')    'idc                        = ', self%idc
    write(6,'(a,i5)')    'divpb                      = ', self%divpb
    write(6,'(a,i5)')    'nquant (from iqmatoms)     = ', self%nquant
    jstart = 1
    do i = 1, size(self%iqmatoms) / jstep + 1
       if ( i == 1 ) then
          tmpstr = 'iqmatoms                   = '
       else
          tmpstr = '                           = '
       end if
       jend = min ( (jstart + jstep - 1), size(self%iqmatoms) )
       write(6,'(a,10(i5))') tmpstr, (self%iqmatoms(j), j = jstart, jend)
       jstart = jstart + jstep
    end do
#ifdef OPENMP
    write(6,'(a,i5)')    'qmmm_omp_max_threads       = ', self%qmmm_omp_max_threads
#endif
    write(6,'(a,e10.2)') 'chg_lambda                 = ', self%chg_lambda
    write(6,'(a,i5)')    'dftb_maxiter               = ', self%dftb_maxiter
    write(6,'(a,i5)')    'dftb_disper                = ', self%dftb_disper
    write(6,'(a,i5)')    'dftb_chg                   = ', self%dftb_chg
    write(6,'(a,e10.2)') 'dftb_telec                 = ', self%dftb_telec
    write(6,'(a,e10.2)') 'dftb_telec_step            = ', self%dftb_telec_step
    write(6,'(2a)')      'dftb_3rd_order             = ', trim(self%dftb_3rd_order)
    write(6,'(a,l)')     'qmqm_analyt                = ', self%qmqm_analyt
    write(6,'(a,l)')     'qmqm_exc_analyt            = ', self%qmqm_exc_analyt
    write(6,'(a,l)')     'tight_p_conv               = ', self%tight_p_conv
    write(6,'(a,l)')     'printcharges               = ', self%printcharges
    write(6,'(a,i3)')    'printdipole                = ', self%printdipole
    write(6,'(a,l)')     'peptide_corr               = ', self%peptide_corr
    write(6,'(a,l)')     'qmqm_erep_incore           = ', self%qmqm_erep_incore
    write(6,'(a,l)')     'allow_pseudo_diag          = ', self%allow_pseudo_diag
    write(6,'(a,l)')     'qmmmrij_incore             = ', self%qmmmrij_incore
    write(6,'(a,l)')     'writepdb                   = ', self%writepdb
    write(6,'(a,l)')     'qm_pme                     = ', self%qm_pme
    
  end subroutine print_qmmm_nml_type

end module qmmm_nml_module
