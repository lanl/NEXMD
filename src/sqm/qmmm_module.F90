! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains allocatable arrays
!and variables used for coupled potential
!qmmm calculations. If you need access to one
!of the public variables from module you should
!include it in your code with
!use qmmm_module
!
!Main Authors:
!     Ross Walker
!     Mike Crowley
!+++++++++++++++++++++++++++++++++++++++++++++

module qmmm_module

  use constants         , only : A2_TO_BOHRS2, one, zero
  use qmmm_nml_module   , only : qmmm_nml_type
  use qmmm_struct_module, only : qmmm_struct_type
  use qmmm_vsolv_module , only : qmmm_vsolv_type
  use qm2_params_module,  only : qm2_params_type
  
  implicit none

  private

  ! constants
  public :: ALPH_MM, AXIS_TOL, OVERLAP_CUTOFF, EXPONENTIAL_CUTOFF

  ! data types
  public :: qmmm_mpi_structure
  public :: qmmm_scratch_structure, qmmm_div_structure
  public :: qm2_structure, qm2_rij_eqns_structure
  public :: qm_ewald_structure, qm_gb_structure
  public :: qmmm_opnq_structure

  ! objects - these *should* not be public (meaning, globally accessible)
  ! do not use these in new subroutines but rather pass them to
  ! subroutines/functions via explicit interfaces
  ! these live here (in this module) only for historic reasons
  public :: qmmm_nml, qmmm_mpi
  public :: qmmm_scratch, qmmm_div, qmmm_vsolv, qmmm_opnq
  public :: qm2_struct, qm2_rij_eqns, qm2_params
  public :: qmewald, qm_gb

  ! functions and subroutines
  public :: validate_qm_atoms
  public :: qmsort
  public :: get_atomic_number
  public :: allocate_qmmm
  public :: deallocate_qmmm
#ifdef MPI
  public :: qmmm_mpi_setup
#endif

  
  ! -----------------
  ! CUTOFF PARAMETERS
  ! -----------------
  ! Exponential decay factor for the MM atom in the core-core interaction
  _REAL_, parameter :: ALPH_MM = 5.0d0
  
  _REAL_, parameter :: AXIS_TOL = 1.0d-8  !Tolerance at which to define a vector is along the axis.
  _REAL_, parameter :: OVERLAP_CUTOFF = 100.0d0*A2_TO_BOHRS2 !Distance^2 in bohrs at which to assume
                                                             !Gaussian overlap is zero.
  _REAL_, parameter :: EXPONENTIAL_CUTOFF = 30.0d0 !Value of x at which to assume Exp(-x) = zero.
                                                   !AWG: MNDO97, MOPAC2007, DCQTP use 25.0d0
  
  
  ! ----------
  ! DATA TYPES
  ! ----------
  type qm2_structure  !Variables that are specific to qm_routine=2 (qm2)

   !+TJG 01/26/2010   Not sure if this should be in qm
  
     ! The total density matrix
     _REAL_, dimension(:), pointer :: den_matrix => null()

     ! Old total density matrix from previous step
     ! Allocated in qm2_load_params on first call - deallocated by deallocate_qmmm
     _REAL_, dimension(:), pointer :: old_den_matrix => null()

     ! Used by qm2_cnvg as workspace, norbs
     _REAL_, dimension(:), pointer :: old2_density => null()
  
     ! These two guesses are only used when density_predict=1
     ! They contain Pguess(t-1) and Pguess(t-2).
     _REAL_, dimension(:), pointer :: md_den_mat_guess1 => null()
     _REAL_, dimension(:), pointer :: md_den_mat_guess2 => null()
  
     ! Final Fock matrices from previous MD steps. In the case of fock_predict=1 
     ! it contains the previous 4 MD step fock matrices. F4 = t-4, F3 = t-3, F2 = t-2, F1 = t-1.
     _REAL_, dimension(:), pointer :: fock_mat_final4 => null()
     _REAL_, dimension(:), pointer :: fock_mat_final3 => null()
     _REAL_, dimension(:), pointer :: fock_mat_final2 => null()
     _REAL_, dimension(:), pointer :: fock_mat_final1 => null() 
  
     ! Fock matrix
     _REAL_, dimension(:), pointer :: fock_matrix => null()
     
     ! Fock matrix displace plus
     _REAL_, dimension(:,:), pointer :: fock_matrix_dp => null()
     
     ! Fock matrix displace minus
     !_REAL_, dimension(:,:), pointer :: fock_matrix_dm => null()

     ! QM-MM electron repulsion integrals
     _REAL_, dimension(:,:), pointer :: qm_mm_e_repul => null()

     ! QM-QM 2-electron repulsion integrals. This is a big array, allocated by qm_mm and
     ! deallocated by deallocate_qmmm - it needs to be a total of n2el long
     _REAL_, dimension(:), pointer :: qm_qm_2e_repul => null()

     ! The 1-electron matrix - used by routines called from qm2_energy
     ! allocated in qm_mm on first call - deallocated by  deallocate_qmmm
     _REAL_, dimension(:), pointer :: hmatrix => null()

     ! QM-QM electron repulsion integrals
     ! This was originally written as a file to disk in the energy routines and then re-read in the derivative 
     ! routines. Now it is stored in memory. Allocated in qm_mm on first call - deallocated by deallocate_qmmm
     _REAL_, dimension(:,:), pointer :: qm_qm_e_repul => null()

     ! Used in qm2_fock2 routine - 16,nquant_nlink allocated in qm2_load_params deallocated in deallocate_qmmm
     _REAL_, dimension(:,:), pointer :: fock2_ptot2 => null()

     ! Eigen vectors during the SCF - allocated in qm2_load_params
     _REAL_, dimension(:,:), pointer :: eigen_vectors => null()

     ! Mulliken charges at each scf step if calc_mchg_scf is true.
     ! Otherwise only do it at the end of the scf.
     _REAL_, dimension(:), pointer :: scf_mchg => null()
  
     !+TJG 01/26/2010
     ! previous fock matrices for diis extrapolation, matsize x ndiis
     _REAL_, dimension(:,:), pointer :: diis_fock => null()

     ! previous error matrices for diis extrapolation, matsize x ndiis
     ! The error matrices are antisymmetric
     _REAL_, dimension(:,:), pointer :: diis_errmat => null()

     ! The diis extrapolation matrix.  The i,j element is the inner-product of error matrices
     ! (packed as vectors, e.g., the error vectors)
     ! The extra dimension enforces the constraint that the sum of extrapolation coefficients sum to unity
     ! ndiis+1 x ndiis+1
     _REAL_, dimension(:,:), pointer :: diis_mat => null()
     !-TJG 01/26/2010
  
     ! Size of the various packed symmetric matrices. (norbs(norbs+1)/2)
     integer :: matsize

     ! Number of 2 electron repulsion integrals, calculated by
     ! moldat = 50*nheavy(nheavy-1)+10*nheavy*nlight+(nlight*(nlight-1))/2
     integer :: n2el
                                             
     ! Total Number of atomic orbitals.
     integer :: norbs

     ! Number of doubly occupied orbitals
     integer :: nclosed

     ! Number of doubly occupied and singly occupied orbitals
     integer :: nopenclosed

     ! This was originally written as a file to disk in the energy routines and then re-read
     ! in the derivative routines. Now it is stored in memory. Allocated in qm_mm on first call
     ! or if pair list changes too much. See qm_mm_e_repul array above - deallocated by deallocate_qmmm
     integer :: qm_mm_e_repul_allocated

     ! Number of peptide linkages in QM region to apply MM correction to.
     integer :: n_peptide_links

     ! Identity of peptide linkages 4,n_peptide_linkages. 1 to 4 = H-N-C-O atom numbers.
     integer, dimension(:,:), pointer :: peptide_links

     ! If set to true the mulliken charges will be calculated on each SCF iteration.
     logical calc_mchg_scf

  end type qm2_structure  
  
  type  qm2_rij_eqns_structure !This structure is used to store RIJ info for each QM-QM pair and related equations
  !QM-MM                                           !equations. See array_locations.h for the first dimension
   _REAL_, dimension(:,:), pointer ::  qmmmrijdata !offsets.
   integer :: qmmmrij_allocated
  end type qm2_rij_eqns_structure
  
  
  !QMEwald specific structure
  type qm_ewald_structure
     _REAL_, dimension(:), pointer :: kvec !Array storing K vectors (totkq long
     _REAL_, dimension(:,:), pointer :: dkvec !Array storing derivative K vectors (3,totkq long
     _REAL_, dimension(:,:), pointer :: dmkv !used for calculating the ktable
     _REAL_, dimension(:,:,:), pointer :: ktable !Table for storing complex exp(ik,r[j])
                                                 !dimensions = 6,natom,totkq
                                                 !1,x,y = x_cos
                                                 !2,x,y = x_sin
                                                 !3,x,y = y_cos
                                                 !4,x,y = y_sin
                                                 !5,x,y = z_cos
                                                 !6,x,y = z_sin
     _REAL_, dimension(:,:,:), pointer :: qmktable !As Ktable but stores the qmatom copies in a linear 1->nquant fashion.
     _REAL_, dimension(:), pointer :: mmpot !Nquant long, stores the potential at each QM atom due to the MM field.
     _REAL_, dimension(:), pointer :: qmpot  !Nquant long, stores the self energy of the QM atoms to avoid double counting
     _REAL_, dimension(:), pointer :: coulpot  !Nquant long, stores the coulombic potential at each QM atom due to MM atoms.
     _REAL_, dimension(:,:), pointer :: d_ewald_mm !3,natom long stores gradients on MM atoms due to QM-MM Ewald field.
                                                   !Reciprocal forces.
     _REAL_ :: ewald_core !Ewald Potential with QM CORE charges - energy in eV.
     _REAL_ :: mm_recip_e !Reciprocal energy from MM atoms - qm_pme.
     _REAL_ :: kappa ! the EWald coefficient--it should be different from the one sander uses-- by TL
     integer :: totkq  !Total number of kspace vectors
     integer :: natom  !Same as sander's natom, copied here by qm_mm for convenience.
     logical :: ewald_startup !True if this is the very first MD step and we are doing qmewald.
  end type qm_ewald_structure
  
  type qm_gb_structure
     _REAL_, dimension(:), pointer :: qmqm_onefij !Stores the 1.0/fij equations for qm-qm pairs. Since these
                                               !values only depend on Rij and the effective radii 
                                               !they remain fixed during the SCF and so are calculated
                                               !once outside the SCF and then reused inside.
     _REAL_, dimension(:), pointer :: qmqm_kappafij !Stores exp(-kappa*fij) - These are calculated outside of the
                                                    !scf and reused inside. Only allocated and calculated if kappa/=0.0d0
                                                    !in other words saltcon /= 0.0d0.
     _REAL_, dimension(:), pointer :: gb_mmpot  !Nquant long, stores GB potential at each QM atom due to MM atoms.
     _REAL_, dimension(:), pointer :: gb_qmpot  !Nquant long, stores GB potential at each QM atom due to QM atoms.
     _REAL_ :: intdieli, extdieli !1.0d0/intdiel and extdiel respectively
     _REAL_ :: kappa    !Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
                        !T = 298.15, epsext=78.5, kappa = sqrt( 0.10806d0 * saltcon )
                        !scaled kappa by 0.73 to account(?) for lack of ion exlcusions. 
     _REAL_ :: mmcut2   !cut^2 in angstroms^2 from cntrl namelist
     _REAL_ :: one_Arad_beta !alpb_beta/Arad when alpb/=0
     integer, dimension(:,:), pointer :: qmqm_gb_list !1+nquant,nquant - list of qm atoms that interact with current qm atom.
                                                      !+1 because the first entry stores the number of interactions for QM atom y. 
     logical :: saltcon_on !True if saltcon /= 0.0d0
     logical :: alpb_on !True if alpb = 1
  end type qm_gb_structure
  
  type qmmm_mpi_structure
    integer :: commqmmm !Communications within a given set of QMMM threads potentially a subset of processors of commsander.
    integer :: numthreads !Number of threads in commqmmm.
    integer :: mytaskid   !Task id of this thread in commqmmm.
    integer :: openmp_numthreads !When openmp is in use this will contain how many threads would be spawned by the master.
    integer :: natom_start !Where this thread should start for 1,natom loops
    integer :: natom_end   !Where this thread should end for 1,natom loops
    integer :: nquant_nlink_start !Where this thread should start for 1,nquant_nlink loops
    integer :: nquant_nlink_end !Where this thread should end for 1,nquant_nlink loops
    integer :: totkq_count
    integer :: kvec_start !Where this thread starts and ends in 1 to totkq loops.
    integer :: kvec_end
    integer :: two_e_offset !Offset into two electron matrix for this thread
  
  !Below are for matrix type double loop load balancing
    integer ::                        nquant_nlink_istart !These three control load balancing
    integer ::                        nquant_nlink_iend
    integer :: nquant_nlink_loop_extent_begin !
    integer :: nquant_nlink_loop_extent_end !
    integer, dimension(:,:), pointer :: nquant_nlink_jrange !within a do i = 1, nquant_nlink
                                                          !           do j=1,i-1
                                                          !loop. Basically istart and iend define
                                                          !the extent of the outer loop and then
                                                          !(1,i) and (2,i) of jrange define the 
                                                          !limits of the inner loop for this processor.
    logical :: commqmmm_master !True if master thread of commqmmm
  
  end type qmmm_mpi_structure
  
#ifdef OPENMP
  type qmmm_openmp_structure
    integer :: diag_threads  !number of threads to use for diagonalization routines.
    integer :: pdiag_threads !number of threads to use for openmp diagonalization routines.
  end type qmmm_openmp_structure
#endif

  type qmmm_scratch_structure
  !Various scratch arrays used as part of QMMM - one should typically assume that upon leaving a routine
  !the contents of these arrays can be assumed to be junk.
    _REAL_, dimension(:), pointer  :: matsize_red_scratch !Allocated as qm2_struct%matsize when doing MPI during the 
                                                          !load parameters routine. ONLY ALLOCATED IF WE CAN'T
                                                          !DO MPI_IN_PLACE.
                                                          !+1 in size so we can pack extra energies on the end etc.
    _REAL_, dimension(:), pointer :: qm_pme_scratch !natom long scratch for qm_pme - only allocated if doing qm_pme
    _REAL_, dimension(:,:), pointer :: mat_diag_workspace !Matrix diagonalisation workspace - allocated in qm2_load_params_and_allocate
                                                          !norbs,6 for internal diag
                                                          !norbs,1 if lapack diag.
  !The Pseudo diagonalizer needs a total of 5 real scratch arrays and 1 integer scratch array.
  !Passed in Scratch Arrays - these are all allocated in qm2_load_params and only if allow_pseudo_diag is true.
  !ONLY ALLOCATED ON COMMQMMM MASTER THREAD
    _REAL_, dimension(:,:), pointer :: pdiag_scr_norbs_norbs     !(norbs,norbs)
    _REAL_, dimension(:,:), pointer :: pdiag_scr_noccupied_norbs !(noccupied,norbs)
    _REAL_, dimension(:), pointer :: pdiag_vectmp1               !(noccupied*(norbs-noccupied))
    _REAL_, dimension(:), pointer :: pdiag_vectmp2               !(noccupied*(norbs-noccupied))
    _REAL_, dimension(:), pointer :: pdiag_vectmp3               !(noccupied*(norbs-noccupied))
    integer, dimension(:,:), pointer :: pdiag_vecjs                !(2,noccupied*(norbs-noccupied))
  
   _REAL_, dimension(:), pointer :: lapack_dc_real_scr
   _REAL_, dimension(:), pointer :: lapack_dc_int_scr
  !END ONLY ALLOCATED ON COMMQMMM MASTER THREAD
  !Scratch Arrays - These arrays are available for any subroutine to use - they are allocated by allocate_qmmm
  !Each routine that uses one of these arrays for temporary storage should assume that the contents of such array
  !will be garbage once the routine is left or another routine is called.
   _REAL_, dimension(:), pointer :: qm_real_scratch !Real scratch array - 4*natom.
   integer, dimension(:), pointer :: qm_int_scratch !Integer scratch array - 3*natom.
   integer :: lapack_dc_real_scr_aloc !Number of reals allocated for lapack_dc_real_scr
   integer :: lapack_dc_int_scr_aloc !Number of ints allocated for lapack_dc_int_scr_aloc
   integer :: qm_mm_pairs_allocated !Size of expected qm_mm_pairs that scratch arrays and calc_rij_array was allocated to.
  end type qmmm_scratch_structure
  
  type qmmm_div_structure
     !Specific for the divcon version of QMMM
     integer :: ntotatm
     integer, dimension(:), pointer :: all_atom_numbers !atomic numbers of ALL atoms (MM and QM)
  end type qmmm_div_structure

    type qmmm_opnq_structure
        ! only variables to be set up in qm_mm startup are stored here
        ! the qmmm_opnq variabel serves as the common data connection between
        ! the OPNQ modules and other modules
        ! the "real" opnq parameters are stored in qm2_params
        logical::useOPNQ=.false.
        _REAL_::OPNQCorrection, vdWCorrection ! in EV
        logical::switching=.true.
        _REAL_::NB_cutoff  ! the non-bond cutoff used in the MM region--required for MM-correction in OPNQ
        _REAL_::switch_cutoff1  ! the distance where the opnq correction will begin being switched off        
        _REAL_::switch_cutoff2  ! the distance where the opnq correction will be totally zeroed
        integer, dimension(:), pointer::MM_atomType  ! the MM atom types for each atom
        logical, dimension(:), pointer::supported
        integer, dimension(:), pointer::atomic_number ! atomic number for each atom
        _REAL_, dimension(:), pointer::LJ_r, LJ_epsilon  ! the MM 6-12 parameters for each type
    end type qmmm_opnq_structure

  ! --------------
  ! GLOBAL OBJECTS
  ! --------------
  ! these should really not live here
  ! but need to be global for historic reasons - too much work to disentangle sander
  ! do *not* use these as globals in new subroutines!
  type(qmmm_nml_type)   , save :: qmmm_nml
  type(qm2_structure)   , save,target :: qm2_struct
  type(qm2_params_type)        :: qm2_params
  type(qm2_rij_eqns_structure) :: qm2_rij_eqns
  type(qm_ewald_structure)     :: qmewald
  type(qm_gb_structure)        :: qm_gb
  type(qmmm_mpi_structure)     :: qmmm_mpi
#ifdef OPENMP
  type(qmmm_openmp_structure)  :: qmmm_omp
#endif
  type(qmmm_scratch_structure) :: qmmm_scratch
  type(qmmm_vsolv_type) , save :: qmmm_vsolv
  type(qmmm_div_structure)     :: qmmm_div
  type(qmmm_opnq_structure),save::qmmm_opnq
  
contains
  
  
#ifdef MPI
  subroutine qmmm_mpi_setup(qmmm_struct, master, natom )

    ! QMMM specific mpi setup and broadcasts

    use qmmm_nml_module, only : broadcast
    use qmmm_struct_module, only : broadcast
    use qmmm_vsolv_module, only : broadcast

    implicit none
  
#  include "parallel.h"
     include 'mpif.h'

     logical, intent(in) :: master
     integer, intent(in) :: natom

     type(qmmm_struct_type), intent(inout) :: qmmm_struct
     integer :: mpi_division, i, istartend(2)
     integer :: loop_extent, loop_extent_begin, loop_extent_end
     integer :: j, jstart, jend
     integer :: ier=0, istatus=0

     !Note, this routine should be used for namelist variables only. It also broadcasts a few things
     !      that are static params from the prmtop such as the MM charges. It should NOT be used for
     !      any QM parameters since ALL threads call the qm2_load_params_and_allocate routine.

     call broadcast(qmmm_nml)

     call broadcast(qmmm_struct, qmmm_nml%qmmm_int, qmmm_nml%qmmm_switch)

     call broadcast(qmmm_vsolv)

#ifdef OPENMP
     ! for the moment diag_threads and pdiag_threads are just
     ! set to qmmm_omp_max_threads
     call mpi_bcast(qmmm_omp%diag_threads, 1, mpi_integer, 0, commsander, ier) 
     call mpi_bcast(qmmm_omp%pdiag_threads, 1, mpi_integer, 0, commsander, ier) 
#endif

     call mpi_bcast(qm2_struct%calc_mchg_scf,1, mpi_logical, 0, commsander, ier)
  
     call mpi_bcast(qm_gb%alpb_on,1,mpi_logical, 0, commsander, ier)
     call mpi_bcast(qm_gb%saltcon_on,1,mpi_logical, 0, commsander, ier)
     call mpi_bcast(qm_gb%kappa,1,MPI_DOUBLE_PRECISION, 0, commsander, ier)

     ! Now we do the arrays, note we have to allocate on all the non-master nodes
     ! Note: We don't need to broadcast the pair list as this is done in qm_mm routine.
     if ( .not. master ) then

        allocate (qmmm_scratch%qm_real_scratch(4*natom), stat=ier )
        REQUIRE(ier == 0)

        allocate (qmmm_scratch%qm_int_scratch(3*natom), stat=ier )
        REQUIRE(ier == 0)

        ! scf_mchg has not been allocated yet on non-master threads so allocate it.
        allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier )
        REQUIRE(ier == 0) !Deallocated in deallocate qmmm
  
     end if

     ! Set the master flag on all threads
     qmmm_mpi%commqmmm_master = master
  
     ! Setup the commqmmm communicator. For regular QMMM MD simulations this will
     ! be the same size as commsander and have the same members. For QMMM LES
     ! simulations this will contain only the current processor and be of numtasks=1
  
     ! In theory we could use this later to make a smaller commqmmm than commsander for efficiency reasons.
  
     ! A single commqmm that is the same on all threads that make up commsander
     qmmm_mpi%commqmmm = commsander
     qmmm_mpi%numthreads = sandersize
     qmmm_mpi%mytaskid = sanderrank
     qmmm_mpi%commqmmm_master = master

     !AWG: The following should not be required if external QM programs are in use
     !     But check below for variable solvent MPI setup
     if (qmmm_nml%qmtheory%EXTERN) return
  
     ! Divide up i=2,nquant_nlink
     !              j=1,i-1
     ! loops
  
     ! Matrix looks like this. E.g. for 4 cpus and 12 QM atoms - thread id's listed
     ! Total elements = 66
     ! cpus 0,1 and 2 do 17 elements each
     ! cpu 3 does 15 elements
     !  i 1 2 3 4 5 6 7 8 9 10 11 12
     ! j1   0 0 0 0 0 0 1 1 2  2  3
     !  2     0 0 0 0 0 1 1 2  2  3
     !  3       0 0 0 1 1 1 2  2  3
     !  4         0 0 1 1 1 2  2  3
     !  5           0 1 1 1 2  2  3
     !  6             1 1 1 2  2  3
     !  7               1 2 2  3  3
     !  8                 2 2  3  3
     !  9                   2  3  3
     ! 10                      3  3
     ! 11                         3
     ! 12 
  
     ! first of all work out out the total extent of the loop.
     ! allocate the memory for the jrange array
  
     allocate(qmmm_mpi%nquant_nlink_jrange(2,qmmm_struct%nquant_nlink),stat=ier)
     REQUIRE(ier==0)
     loop_extent = qmmm_struct%nquant_nlink*(qmmm_struct%nquant_nlink-1)/2
     mpi_division = (loop_extent+(qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
     loop_extent_end = min(mpi_division*(qmmm_mpi%mytaskid+1),loop_extent)
     loop_extent_begin = mpi_division*qmmm_mpi%mytaskid+1
     !loop_extent_begin = (istart-1)(istart-2)/2 + jstart
     !loop_extent_end = (iend-1)(iend-2)/2 + jend
     !s = 1+sqrt(1+8x)/2
     !i = int(s) - ROUNDED UP
     qmmm_mpi%nquant_nlink_istart = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_begin)))/2.0d0)
     qmmm_mpi%nquant_nlink_iend   = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_end)))/2.d0)
     qmmm_mpi%nquant_nlink_loop_extent_begin = loop_extent_begin
     qmmm_mpi%nquant_nlink_loop_extent_end = loop_extent_end
  
     !Now we need to work out what range of j values we do for each i we will be doing.
     !What value of j would, when coupled with our istart give us loop_extent_begin?
     ! j = loop_extent_begin -((-i-1)(i-2)/2)
  
     jstart = loop_extent_begin - ((qmmm_mpi%nquant_nlink_istart-1)*(qmmm_mpi%nquant_nlink_istart-2)/2)
     jend   = loop_extent_end - ((qmmm_mpi%nquant_nlink_iend-1)*(qmmm_mpi%nquant_nlink_iend-2)/2)
  
     do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
  
       if (i == qmmm_mpi%nquant_nlink_istart) then
         qmmm_mpi%nquant_nlink_jrange(1,i) = jstart
       else
         qmmm_mpi%nquant_nlink_jrange(1,i) = 1
       end if
  
       if (i == qmmm_mpi%nquant_nlink_iend) then
         qmmm_mpi%nquant_nlink_jrange(2,i) = jend
       else
         qmmm_mpi%nquant_nlink_jrange(2,i) = i-1
       end if
  
     end do
  
     !------- End Matrix Type Calc Load Balancing ------------
  
     ! Now divide up the atoms between threads
     ! We will divide up evenly as best we can. E.g. with 1603 atoms on 4 threads
     ! we would ideally do 1->401, 402->802, 803->1203, 1204->1603
     ! But this is quite difficult to manage. For the moment I will just divide
     ! up using a simple formula that ensures that most threads have an even spread.
     ! Thus threads 1 to 3 do 401 atoms and thread 4 does 400 atoms. This approach
     ! allows linear movement in memory.
     ! Note the last thread is responsible for ALL link atoms.
  
     ! NATOM division first.
     mpi_division = (natom + (qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
     qmmm_mpi%natom_end = min(mpi_division*(qmmm_mpi%mytaskid+1),natom)
     qmmm_mpi%natom_start = min(mpi_division*qmmm_mpi%mytaskid+1,natom+1)
  
     ! write info about atom division - Get each thread to send the master it's values
     ! This allows a sanity check.
     if (qmmm_mpi%commqmmm_master) then
        write (6,'(/a,i4,a)') '|QMMM: Running QMMM calculation in parallel mode on ',qmmm_mpi%numthreads,' threads.'
        write (6,'(a)') '|QMMM: All atom division among threads:'
        write (6,'(a)') '|QMMM:                  Start       End      Count'
        !Already know my own.
        write(6,'(a,i8,a,i8,a,i8,a)') &
              '|QMMM: Thread(   0): ',qmmm_mpi%natom_start,'->',qmmm_mpi%natom_end, &
                                    '  (',qmmm_mpi%natom_end-qmmm_mpi%natom_start+1,')'
        do i = 1, qmmm_mpi%numthreads-1
           call mpi_recv(istartend,2,mpi_integer,i,0,qmmm_mpi%commqmmm,istatus,ier)
           write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
               '|QMMM: Thread(',i,'): ',istartend(1),'->',istartend(2), &
                                    '  (',istartend(2)-istartend(1)+1,')'
       end do
     else
       ! Send a message to the master with our counts in.
       istartend(1) = qmmm_mpi%natom_start
       istartend(2) = qmmm_mpi%natom_end
       call mpi_send(istartend,2,mpi_integer,0,0,qmmm_mpi%commqmmm,ier)
     end if
  
     ! Nquant
     mpi_division = (qmmm_struct%nquant + (qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
     qmmm_mpi%nquant_nlink_end = min(mpi_division*(qmmm_mpi%mytaskid+1),qmmm_struct%nquant)
     qmmm_mpi%nquant_nlink_start = min(mpi_division*qmmm_mpi%mytaskid+1,qmmm_struct%nquant+1)
  
     ! Nquant_nlink - Last thread gets all the link atoms
     if (qmmm_mpi%mytaskid == qmmm_mpi%numthreads - 1) then
        ! Last thread
        qmmm_mpi%nquant_nlink_end = qmmm_mpi%nquant_nlink_end + qmmm_struct%nlink
     end if
  
     if (qmmm_mpi%commqmmm_master) then
       write (6,'(/a)') '|QMMM: Quantum atom + link atom division among threads:'
       write (6,'(a)') '|QMMM:                  Start       End      Count'
       ! Already know my own
       write(6,'(a,i8,a,i8,a,i8,a)') &
             '|QMMM: Thread(   0): ',qmmm_mpi%nquant_nlink_start,'->',qmmm_mpi%nquant_nlink_end, &
             '  (',qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1,')'
       do i = 1, qmmm_mpi%numthreads-1
         call mpi_recv(istartend,2,mpi_integer,i,0,qmmm_mpi%commqmmm,istatus,ier)
         write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
             '|QMMM: Thread(',i,'): ',istartend(1),'->',istartend(2), &
                                    '  (',istartend(2)-istartend(1)+1,')'
       end do
     else
       ! Send a message to the master with our counts in.
       istartend(1) = qmmm_mpi%nquant_nlink_start
       istartend(2) = qmmm_mpi%nquant_nlink_end
       call mpi_send(istartend,2,mpi_integer,0,0,qmmm_mpi%commqmmm,ier)
     end if
  
     ! Now we need to calculate the offset into the 2 electron matrix that
     ! this thread will have - this depends on the number of light and heavy
     ! atoms and so is calculated at the end of qm2_load_params_and_allocate.
  
  end subroutine qmmm_mpi_setup
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
  
  subroutine allocate_qmmm( qmmm_nml, qmmm_struct, natom )

     ! allocates space for qmmm variables and arrays that depend only on nquant or natom

     use qmmm_nml_module   , only : qmmm_nml_type, new
     use qmmm_struct_module, only : qmmm_struct_type, new
    
     implicit none

     type(qmmm_nml_type)   , intent(inout) :: qmmm_nml
     type(qmmm_struct_type), intent(inout) :: qmmm_struct
     integer, intent(in) :: natom

     integer :: ier
  
     !WARNING - nlink is NOT known when the master thread enters this routine (is zero). 
     !          Thus you cannot do allocations which depend on nlink in here. The exception
     !          is iqmatoms which first gets allocated with nlink = 0 on the master thread
     !          and then gets reallocated to include the link atoms later on.

     !iqmatoms and iqm_atomic_numbers are allocated here as nquant+nlink long but 
     !initially nlink is zero so it only gets allocated as nquant long on the master
     !thread. It is then resized when nlink are known about. When all other mpi threads
     !call this allocation routine nquant+nlink should include both qm and link atoms.
     call new(qmmm_nml)

     call new(qmmm_struct, qmmm_nml%qmmm_int, qmmm_nml%qmmm_switch)

     allocate (qmmm_scratch%qm_real_scratch(4*natom), stat=ier )
     REQUIRE(ier == 0)

     allocate (qmmm_scratch%qm_int_scratch(3*natom), stat=ier )
     REQUIRE(ier == 0)

  end subroutine allocate_qmmm
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
    
    use qmmm_nml_module   , only : qmmm_nml_type, delete
    use qmmm_struct_module, only : qmmm_struct_type, delete
    use qmmm_vsolv_module , only : qmmm_vsolv_type, delete
    use qm2_params_module , only : qm2_params_type, delete

    implicit none

    type(qmmm_nml_type)   , intent(inout) :: qmmm_nml
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(qmmm_vsolv_type) , intent(inout) :: qmmm_vsolv
    type(qm2_params_type), intent(inout) :: qm2_params

    integer :: ier=0

    !If this is a parallel run the non master threads will only have
    !allocated this memory if LES is on since otherwise QM calc is
    !currently only done on master thread.
    !Deallocate pointers 

    call delete(qmmm_struct, qmmm_nml%qmmm_int, qmmm_nml%idc, qmmm_nml%qmmm_switch)

    call delete(qmmm_nml)

    call delete(qmmm_vsolv)

    ! The rest of this subroutine still needs to be cleaned up

    !Deallocate qm_ewald memory if it is on.
    if ( qmmm_nml%qm_ewald>0 ) then
       deallocate ( qmewald%dmkv, stat = ier)
       REQUIRE(ier == 0)
       deallocate ( qmewald%dkvec, stat = ier)
       REQUIRE(ier == 0)
       deallocate ( qmewald%kvec, stat = ier)
       REQUIRE(ier == 0)
       if (.not. qmmm_nml%qm_pme) then
          deallocate ( qmewald%ktable, stat = ier )
          REQUIRE(ier == 0)
          deallocate ( qmewald%d_ewald_mm, stat = ier )
          REQUIRE(ier == 0)
       else
          deallocate ( qmmm_scratch%qm_pme_scratch, stat = ier )
          REQUIRE(ier == 0)
       end if
       deallocate ( qmewald%qmktable, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qmewald%mmpot, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qmewald%qmpot, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qmewald%coulpot, stat = ier )
       REQUIRE(ier == 0)       
    end if

    !Deallocate QM-GB arrays - only used with qmgb=2.
    if ( qmmm_nml%qmgb == 2 ) then
       deallocate ( qm_gb%qmqm_onefij, stat = ier )
       REQUIRE(ier == 0)
       if ( qm_gb%saltcon_on) then
          deallocate ( qm_gb%qmqm_kappafij, stat = ier )
          REQUIRE(ier == 0)
       end if
       deallocate ( qm_gb%gb_mmpot, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm_gb%gb_qmpot, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm_gb%qmqm_gb_list, stat = ier )
       REQUIRE(ier == 0)
    end if

    ! The rest was not allocated if we are using an external program
    if (qmmm_nml%qmtheory%EXTERN) return

    !MPI Specific deallocations
#ifdef MPI
    deallocate(qmmm_mpi%nquant_nlink_jrange, stat=ier )
    REQUIRE(ier == 0)
# ifndef USE_MPI_IN_PLACE
    deallocate(qmmm_scratch%matsize_red_scratch, stat=ier )
    REQUIRE(ier == 0)
# endif
#endif

    if (qmmm_nml%qmqm_erep_incore) then
       deallocate ( qm2_struct%qm_qm_e_repul, stat = ier )
       REQUIRE(ier == 0)
    end if
    if (qmmm_nml%idc==0) then
       if (qm2_struct%n_peptide_links>0) then
          deallocate ( qm2_struct%peptide_links, stat=ier )
          REQUIRE(ier == 0)
       end if
       if (.not. qmmm_nml%qmtheory%DFTB) then
          deallocate ( qm2_struct%eigen_vectors, stat=ier )
          REQUIRE(ier == 0)
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
             deallocate ( qmmm_scratch%mat_diag_workspace, stat=ier )
             REQUIRE(ier == 0)
          end if
       end if
       if (qmmm_scratch%lapack_dc_real_scr_aloc /= 0) then
          deallocate ( qmmm_scratch%lapack_dc_real_scr, stat=ier)
          REQUIRE(ier == 0)
       end if
       if (qmmm_scratch%lapack_dc_int_scr_aloc /= 0) then
          deallocate ( qmmm_scratch%lapack_dc_int_scr, stat=ier)
          REQUIRE(ier == 0)
       end if
       if (qmmm_nml%allow_pseudo_diag .and. qmmm_mpi%commqmmm_master) then
          deallocate(qmmm_scratch%pdiag_scr_norbs_norbs, &
               qmmm_scratch%pdiag_scr_noccupied_norbs, &
               qmmm_scratch%pdiag_vectmp1, &
               qmmm_scratch%pdiag_vectmp2, &
               qmmm_scratch%pdiag_vectmp3, &
               qmmm_scratch%pdiag_vecjs, &
               stat=ier)
          REQUIRE(ier == 0)
       end if

       if (qmmm_mpi%commqmmm_master .and. qmmm_nml%diag_routine==7 .and. (.not. qmmm_nml%allow_pseudo_diag)) then
          !we allocated pdiag_scr_norbs_norbs to store the unpacked matrix so make sure
          !we deallocate it.
          deallocate(qmmm_scratch%pdiag_scr_norbs_norbs, stat=ier)
          REQUIRE(ier == 0)
       end if
       deallocate ( qm2_struct%fock2_ptot2, stat=ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%hmatrix, stat = ier )
       REQUIRE(ier == 0)
       deallocate (qm2_struct%qm_qm_2e_repul, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%fock_matrix, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%old2_density, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%old_den_matrix, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%den_matrix, stat = ier )
       REQUIRE(ier == 0)



       if (qmmm_nml%density_predict == 1) then
          !We are using Niklasson et al density matrix prediction algorithm.
          deallocate ( qm2_struct%md_den_mat_guess1, stat = ier )
          REQUIRE(ier == 0)
          deallocate ( qm2_struct%md_den_mat_guess2, stat = ier )
          REQUIRE(ier == 0)
       end if

       if (qmmm_nml%fock_predict == 1) then
          !We are using Pulay et al matrix prediction algorithm.
          deallocate ( qm2_struct%fock_mat_final4, stat = ier )
          REQUIRE(ier == 0)
          deallocate ( qm2_struct%fock_mat_final3, stat = ier )
          REQUIRE(ier == 0)
          deallocate ( qm2_struct%fock_mat_final2, stat = ier )
          REQUIRE(ier == 0)
          deallocate ( qm2_struct%fock_mat_final1, stat = ier )
          REQUIRE(ier == 0)
       end if

       !+TJG 01/26/2010  Hm.  You should check association before deallocate()
       ier = 0
       IF ( ASSOCIATED( qm2_struct%diis_fock ) ) DEALLOCATE( qm2_struct%diis_fock , stat = ier )
       REQUIRE(ier == 0)
       IF ( ASSOCIATED( qm2_struct%diis_errmat ) ) DEALLOCATE( qm2_struct%diis_errmat , stat = ier )
       REQUIRE(ier == 0)
       IF ( ASSOCIATED( qm2_struct%diis_mat ) ) DEALLOCATE( qm2_struct%diis_mat , stat = ier )
       REQUIRE(ier == 0)
       !-TJG 01/26/2010

       call delete(qm2_params, qmmm_nml%qmtheory, qmmm_struct%PM3MMX_INTERFACE)

    end if !if idc==0

    if (qmmm_nml%qmmmrij_incore) then
       if (qmmm_mpi%commqmmm_master) then
          deallocate (qm2_rij_eqns%qmmmrijdata, stat = ier )
          REQUIRE(ier == 0)
       end if
    end if

    deallocate ( qmmm_scratch%qm_int_scratch, stat = ier )
    REQUIRE(ier == 0)

    deallocate ( qmmm_scratch%qm_real_scratch, stat = ier )
    REQUIRE(ier == 0)

    !Deallocate the scf Mulliken charge array.
    !Was allocated on all cpus.
    deallocate ( qm2_struct%scf_mchg, stat = ier )
    REQUIRE(ier == 0)

    ! GMS: Dellocate DFTB arrays
    if ((qmmm_nml%qmtheory%DFTB).and.(qmmm_mpi%commqmmm_master)) call qm2_dftb_deallocate

  end subroutine deallocate_qmmm
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !+Sorts the array iqmatoms into numerical order
  subroutine qmsort(qmmm_struct, iqmatoms)
  
    implicit none
    integer, intent(inout) :: iqmatoms(*)
    type(qmmm_struct_type), intent(in) :: qmmm_struct
  
    ! Local
    integer i,j,lcurrent
  
    !  sort array iqmatoms in ascending order
    !  sort only over nquant atoms don't sort the link atom
    !  MM link pair atoms on the end.
    
    do i = 1, qmmm_struct%nquant
       lcurrent = iqmatoms(i)
       do j = i+1,qmmm_struct%nquant
          if (lcurrent.gt.iqmatoms(j)) then
             iqmatoms(i) = iqmatoms(j)
             iqmatoms(j) = lcurrent
             lcurrent = iqmatoms(i)
          endif
       end do
    end do
  end subroutine qmsort
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ----------------------
  ! Identify atom elements
  ! ----------------------
  subroutine get_atomic_number(atom_name,atom_mass,atomic_number, errorFlag)
    !
    ! Assign atomic number based upon the first letter of the atom symbol
    ! which has been read in from the topology file and based upon the mass.
    ! The assumption here is that an atom name matches the first letter of the
    ! element. If it doesn't then this routine will need to be modified.
    
    use UtilitiesModule, only : Upcase

    implicit none
    
    character(len=4), intent(in)  :: atom_name
    _REAL_,           intent(in)  :: atom_mass
    logical, optional, intent(out) :: errorFlag
    integer,          intent(out) :: atomic_number
    
    logical::localErrorFlag
    localErrorFlag=.false.
    
    ! Lanthanides are not supported.
    ! Actinides are not supported.
    
    if( Upcase(atom_name(1:1)) .eq. 'A' ) then
       if(atom_mass > 24.0d0 .and. atom_mass <= 28.0d0) then
          atomic_number =  13 !Aluminium
       elseif(atom_mass > 35.0d0 .and. atom_mass <= 40.0d0) then
          atomic_number =  18 !Argon
       elseif(atom_mass > 73.0d0 .and. atom_mass <= 77.0d0) then
          atomic_number =  33 !Arsenic
       elseif(atom_mass > 106.0d0 .and. atom_mass <= 109.0d0) then
          atomic_number =  47 !Silver
       elseif(atom_mass > 195.0d0 .and. atom_mass <= 199.0d0) then
          atomic_number =  79 !Gold
       elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number =  85 !Astatine
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'B' ) then
       if(atom_mass > 8.0d0 .and. atom_mass <= 10.0d0) then
          atomic_number =  4 !Beryllium
       elseif(atom_mass > 10.0d0 .and. atom_mass <= 12.0d0) then
          atomic_number =  5 !Boron
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number =  35 !Bromine
       elseif(atom_mass > 135.0d0 .and. atom_mass <= 139.0d0) then
          atomic_number =  56 !Barium
       elseif(atom_mass > 207.0d0 .and. atom_mass <= 211.0d0) then
          atomic_number =  83 !Bismuth
       else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'C' ) then
       if(atom_mass > 10.0d0 .and. atom_mass <= 14.0d0) then
          atomic_number =  6 !Carbon
       elseif(atom_mass > 33.0d0 .and. atom_mass <= 37.0d0) then
          atomic_number =  17 !Chlorine
       elseif(atom_mass > 38.0d0 .and. atom_mass <= 42.0d0) then
          atomic_number =  20 !Calcium
       elseif(atom_mass > 50.0d0 .and. atom_mass <= 54.0d0) then
          atomic_number =  24 !Chromium
       elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number =  27 !Cobalt
       elseif(atom_mass > 61.0d0 .and. atom_mass <= 65.0d0) then
          atomic_number =  29 !Copper
       elseif(atom_mass > 110.0d0 .and. atom_mass <= 114.0d0) then
          atomic_number =  48 !Cadmium
       elseif(atom_mass > 131.0d0 .and. atom_mass <= 135.0d0) then
          atomic_number =  55 !Cesium
       else
          localErrorFlag=.true.
       endif
    
    elseif( Upcase(atom_name(1:1)) .eq. 'D' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'E' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'F' ) then
       if(atom_mass > 17.0d0 .and. atom_mass <= 21.0d0) then
          atomic_number =  9 !Fluorine
       elseif(atom_mass > 54.0d0 .and. atom_mass <= 58.0d0) then
          atomic_number =  26 !Iron
       elseif(atom_mass > 218.0d0 .and. atom_mass <= 228.0d0) then
          atomic_number =  87 !Francium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'G' ) then
       if(atom_mass > 67.0d0 .and. atom_mass <= 71.0d0) then
          atomic_number =  31 !Gallium
       elseif(atom_mass > 71.0d0 .and. atom_mass <= 75.0d0) then
          atomic_number =  32 !Germanium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'H' ) then
       if(atom_mass > 0.0d0 .and. atom_mass <= 2.0d0) then
          atomic_number =  1 !Hydrogen
       elseif(atom_mass > 3.0d0 .and. atom_mass <= 5.0d0) then
          atomic_number =  2 !Helium
       elseif(atom_mass > 176.0d0 .and. atom_mass <= 180.0d0) then
          atomic_number =  72 !Hafnium
       elseif(atom_mass > 198.0d0 .and. atom_mass <= 202.0d0) then
          atomic_number =  80 !Mercury
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'I' ) then
       if(atom_mass > 112.0d0 .and. atom_mass <= 116.0d0) then
          atomic_number = 49 !Indium
       elseif(atom_mass > 125.0d0 .and. atom_mass <= 129.0d0) then
          atomic_number =  53 !Iodine
       elseif(atom_mass > 190.0d0 .and. atom_mass <= 194.0d0) then
          atomic_number =  77 !Iridium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'J' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'K' ) then
       if(atom_mass > 37.0d0 .and. atom_mass <= 41.0d0) then
          atomic_number = 19 !Potassium
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 86.0d0) then
          atomic_number = 36 !Krypton
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'L') then
       if(atom_mass > 6.0d0 .and. atom_mass <= 8.0d0) then
          atomic_number = 3 !Lithium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'M' ) then
       if(atom_mass > 22.0d0 .and. atom_mass <= 26.0d0) then
          atomic_number = 12 !Magnesium
       elseif(atom_mass > 53.0d0 .and. atom_mass <= 57.0d0) then
          atomic_number = 25 !Manganese
       elseif(atom_mass > 94.0d0 .and. atom_mass <= 98.0d0) then
          atomic_number = 42 !Molybdenem
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'N') then
       if(atom_mass > 13.0d0 .and. atom_mass <= 15.0d0) then
          atomic_number = 7 !Nitrogen
       elseif(atom_mass > 19.0d0 .and. atom_mass <= 22.0d0) then
          atomic_number = 10 !Neon
       elseif(atom_mass > 22.1d0 .and. atom_mass <= 23.0d0) then
          atomic_number = 11 !Sodium
       elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number = 28 !Nickel
       elseif(atom_mass > 95.0d0 .and. atom_mass <= 99.0d0) then
          atomic_number = 41 !Niobium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'O' ) then
       if(atom_mass > 14.0d0 .and. atom_mass <= 18.0d0) then
          atomic_number = 8 !Oxygen
       elseif(atom_mass > 188.0d0 .and. atom_mass <= 192.0d0) then
          atomic_number = 76 !Osmium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'P' ) then
       if(atom_mass > 29.0d0 .and. atom_mass <= 33.0d0) then
          atomic_number = 15 !Phosphorus
       elseif(atom_mass > 104.0d0 .and. atom_mass <= 108.0d0) then
          atomic_number = 46 !Palladium
       elseif(atom_mass > 193.0d0 .and. atom_mass <= 197.0d0) then
          atomic_number = 78 !Platinum
       elseif(atom_mass > 205.0d0 .and. atom_mass <= 208.0d0) then
          atomic_number = 82 !Lead
       elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number = 84 !Polonium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'Q' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'R' ) then
       if(atom_mass > 84.0d0 .and. atom_mass <= 88.0d0) then
          atomic_number = 37 !Rubidium
       elseif(atom_mass > 99.0d0 .and. atom_mass <= 102.0d0) then
          atomic_number = 44 !Ruthenium
       elseif(atom_mass > 102.0d0 .and. atom_mass <= 105.0d0) then
          atomic_number = 45 !Rhodium
       elseif(atom_mass > 184.0d0 .and. atom_mass <= 188.0d0) then
          atomic_number = 75 !Rhenium
       elseif(atom_mass > 210.0d0 .and. atom_mass <= 222.5d0) then
          atomic_number = 86 !Radon
       elseif(atom_mass > 223.0d0 .and. atom_mass <= 229.0d0) then
          atomic_number = 88 !Radium
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'S' ) then
       if(atom_mass > 26.0d0 .and. atom_mass <= 30.0d0) then
          atomic_number = 14 !Silicon
       elseif(atom_mass > 30.0d0 .and. atom_mass <= 34.0d0) then
          atomic_number = 16 !Sulphur
       elseif(atom_mass > 43.0d0 .and. atom_mass <= 47.0d0) then
          atomic_number = 21 !Scandium
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number = 34 !Selenium
       elseif(atom_mass > 86.0d0 .and. atom_mass <= 89.0d0) then
          atomic_number = 38 !Strontium
       elseif(atom_mass > 116.0d0 .and. atom_mass <= 120.0d0) then
          atomic_number = 50 !Tin
       elseif(atom_mass > 120.0d0 .and. atom_mass <= 124.0d0) then
          atomic_number = 51 !Antimony
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'T' ) then
       if(atom_mass > 46.0d0 .and. atom_mass <= 50.0d0) then
          atomic_number = 22 !Titanium
       elseif(atom_mass > 96.0d0 .and. atom_mass <= 100.0d0) then
          atomic_number = 43 !Technetium
       elseif(atom_mass > 125.0d0 .and. atom_mass <= 130.0d0) then
          atomic_number = 52 !Tellurium
       elseif(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 73 !Tantalum
       elseif(atom_mass > 201.0d0 .and. atom_mass <= 206.0d0) then
          atomic_number = 81 !Thallium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'U' ) then
       write(6,*) 'Unable to correctly identify element ', atom_name
       call mexit(6,1)

    elseif( Upcase(atom_name(1:1)) .eq. 'V' ) then
       if(atom_mass > 49.0d0 .and. atom_mass <= 53.0d0) then
          atomic_number = 23 !Vanadium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'W' ) then
       if(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 74 !Tungsten
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'X' ) then
       if (atom_mass > 123.0d0 .and. atom_mass < 136.0d0) then
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'Z' ) then
       if(atom_mass > 61.0d0 .and. atom_mass <= 69.0d0) then
          atomic_number = 30 !Zinc
       elseif(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
          atomic_number = 40 !Zirconium
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'Z' ) then
       if(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
          atomic_number = 40 !Zirconium
       else
         localErrorFlag=.true.
       end if

    else
       localErrorFlag=.true.
    endif
    
    if (localErrorFlag) then
       if (present(errorFlag)) then
           ! not to force the program to terminate--left to the calling program
           ! what to do next---Taisung Lee (Rutgers, 2011)
           errorFlag=.true.
       else
           write(6,*) 'Unable to correctly identify element ', atom_name
           call mexit(6,1)
       end if          
    end if
  
  end subroutine get_atomic_number
  
  subroutine validate_qm_atoms(iqmatoms, nquant, natom)

    ! Validate range of QM atom numbers
    !
    ! Check the list of atoms numbers stored in iqmatoms for
    !
    ! 1) All are >= 1 .and. <= natoms
    ! 2) All are unique integer numbers
    !
    ! Written by Ross Walker, TSRI, 2004
  
    implicit none
    
    !Passed in
    integer, intent(in) :: nquant, natom
    integer, intent(in) :: iqmatoms(nquant)
    
    !Local
    integer :: icount1, icount2, iatom
    
    ! Sanity check 1, ensure nquant isn't bigger than natom (it can't be)
    if ((nquant < 1) .OR. (nquant > natom)) then
       write (6,'(" QM ATOM VALIDATION: nquant has a value of ",i8)') nquant
       write (6,'(" which is bigger than natom of ",i8,". Need 0 < nquant <= natom.")') natom
       call sander_bomb('validate_qm_atoms','nquant illegal', 'Need 0 < nquant <= natom')
    end if

    ! Check 2 - loop over nquant atoms and check it is > 1 and <= natom and it is unique
    do icount1=1,nquant
       !check atom number is legal
       iatom = iqmatoms(icount1)
       if ( (iatom > 0) .AND. (iatom <= natom)  ) then
          !QM atom ID is valid - check it is unique
          do icount2=(icount1+1),nquant
             if ( iatom == iqmatoms(icount2) ) then
                write (6,'(" QM ATOM VALIDATION: qm atom ",i8," is not unique.")') iatom
                call sander_bomb('validate_qm_atoms',&
                     'QM atoms specified with iqmatoms do not form a unique set.', &
                     'Require a unique list of qm atom numbers')
             end if
          end do
       else
          !it is not legal
          write (6,'(" QM ATOM VALIDATION: iqmatom ID number of ",i8," is not valid.")') iatom
          call sander_bomb('validate_qm_atoms','invalid QM atom ID', 'Need 0 < QM atom ID <= natom')
       end if
    end do
  
  end subroutine validate_qm_atoms
  
  !END SUBROUTINES
end module qmmm_module

