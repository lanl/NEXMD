!
!********************************************************************
!
!  Module for Davidson calculations
!
!  by Kirill A Velizhanin (kirill{at}lanl{dot}gov)
!
!********************************************************************
!
#include "dprec.fh"
#include "assert.fh"
!
   module qm2_davidson_module
   implicit none
!
!--------------------------------------------------------------------
!
!  Global variables/arrays/user-defined types and structures
!
!--------------------------------------------------------------------
!
   type qm2_davidson_structure_type

      logical :: calcxdens ! flag for calculating cross densities !JAKB

      logical :: has_been_run = .FALSE. ! Has the Davidson calculation been run at least once?

      integer::verbosity=0 ! verbosity of Davidson output

      ! General davidson control parameters
      integer::mdflag=0 ! default for no batching and no shifting
      integer::irflag=0 ! ???
      integer::idav=1 ! 1 for Tamm-Dancoff (CIS), 2 for RPA
      integer::dav_guess=1 ! use (1) or not (0) previous Davidson as guess
      integer::icount_M=100 ! max iterations within Davidson
      integer::iloop_M=50 ! max restarts within Davidson
      integer::iderivfl=0 ! analytic derivatives ???

   ! Controlling SCF quitting CML TEST 7/15/12
	  logical::minimization=.FALSE. ! T if a minimization process is to be done

      ! Lancsoz parameters
      integer::Mx=0 ! number of excited states to be computed 
      _REAL_::ftol=0.d0 !  Min tolerance (|emin-eold|)
      _REAL_::ftol0=1.d-8 !  Acceptance tol.(|emin-eold|)
      _REAL_::ftol1=1.d-8 ! Accept.tol.for residual norm

      ! Non-adiabatic dynamics
      _REAL_, allocatable :: dij(:)	! NAD coupling
      _REAL_, allocatable :: cmdqt(:,:) ! Coefficients
      ! Various sizes
      integer Nb ! basis size, total number of atomic orbitals
      integer Lt ! =Nb*(Nb+1)/2 - triagonal
      integer Np ! number of pairs of electrons, i.e., #occupied orbitals
      integer Nh ! number of virtual orbitals, Nh=Nb-Np
      integer Ncis ! size of cis matrix, Np*Nh
      integer Nrpa ! size of the RPA matrix, 2*Np*Nh
      integer Nqat ! number of "quantum" atoms

      ! Ground state Hartree-Fock (HF) data
      _REAL_,pointer::ehf(:) ! ehf(1:Nb), HF eigenvalues
      _REAL_,pointer::vhf(:,:) ! vhf(1:Nb,1:Nb), HF eigenvectors
      _REAL_,allocatable::vhf_old(:,:) ! old HF eigenvectors to preserve quirality
      _REAL_ Eground ! total energy of ground state

      ! Pointer to Coulomb matrix elements
      _REAL_,pointer::W(:)

      ! Transition energies - Davidson eigenvalues
      _REAL_,allocatable::e0(:)

      ! Total energies of excited states 
      ! i.e., ground state total energies+transition energies
      _REAL_,allocatable::Etot(:)

      ! Total energy of the required state
      ! excited (Mx>0), ground (Mx==0)
      _REAL_ Ereq
      integer struct_opt_state ! cml-test state to optimize structure

      ! Gradient (force) of required excited state
      _REAL_,allocatable::grad(:) 

      ! Transition density - Davidson eigenvectors
      _REAL_,allocatable::v0(:,:),v0_old(:,:)
      _REAL_,allocatable::v2(:,:) ! small space for cross-terms

      ! Various auxiliary data
      _REAL_,allocatable::ferr(:) ! resulting error of found modes 
      _REAL_,allocatable::rrwork(:) ! Lanczos vectors
      integer,allocatable::ix(:),jx(:),kx(:) ! index arrays
      _REAL_,allocatable::vexp1(:),vexp(:) ! Davidson expansion vectors
      _REAL_,allocatable::temp1(:),temp2(:)

      ! Spare space
      _REAL_,allocatable::eta(:),etas(:),eta_scratch(:),eta_scratch_2(:)
      _REAL_,allocatable::xi(:),xis(:),xi_scratch(:),xi_scratch_2(:)
      _REAL_,allocatable::eta_tz(:),xi_tz(:)
      _REAL_,allocatable::tz_scratch(:) ! replacement for rrwork in calc_rhotz()
      _REAL_,allocatable::nacr_scratch(:) ! replacement for rrwork in nacr functions

      ! Davidson Rayleigh matrices
      _REAL_,allocatable::ray1(:,:),ray2(:,:),ray1a(:,:)
      _REAL_,allocatable::ray(:,:)
      _REAL_,allocatable::raye(:),raye1(:)
      _REAL_,allocatable::rayv(:,:),rayvL(:,:),rayvR(:,:)

      ! Dimension of Krylov expansion in davidson
      integer::nd

      ! Batching
      integer Mj ! ???
      integer::idavfrac=10 ! Fraction of the davidson space for 1 batch
      _REAL_::fs=2.d0 ! Frequency shift for many states

      ! Discrimination values for expansion vectors in Davidson
      _REAL_::tresh=0.9
      _REAL_::tresh1=0.999999

      ! Density matrix information
      _REAL_,allocatable :: rhoTZ(:)
      _REAL_,allocatable :: rhoT(:)
      _REAL_,allocatable :: rhoLZ(:)

      ! Excited state gradient information
      _REAL_, allocatable :: dxyz(:) 

      ! COSMO parameters do not really belong here
      ! but it is easy to have them here for time being
      !_REAL_::ceps ! COSMO dielectric permittivity

   end type qm2_davidson_structure_type


   public :: allocate_davidson, deallocate_davidson
   
   type(qm2_davidson_structure_type),target,save :: qm2ds

   logical, private :: initialized = 0 ! initially zero, i.e. not initialized
   contains
!
!********************************************************************
!
!  Davidson initialization and allocation
!
!********************************************************************
!
   subroutine allocate_davidson()
   use qmmm_module,only:qm2_struct,qmmm_struct,qmmm_scratch,qmmm_nml

   implicit none

   integer i,j,ierr

   REQUIRE(qm2ds%Mx .GT. 0)

   if (initialized) return 


   initialized = .TRUE. ! allocated

!-----------INITIALIZE "CONSTANTS"--------------
   ! Nb - the basis size, i.e., total number of atomic orbitals
   qm2ds%Nb=qm2_struct%norbs
   qm2ds%Lt=qm2ds%Nb*(qm2ds%Nb+1)/2

   ! Np - number of pairs of electrons, i.e. number of occupied orbitals
   qm2ds%Np=qm2_struct%nclosed
   qm2ds%Nh=qm2ds%Nb-qm2ds%Np

   ! Size of the davidson CIS and RPA matrices
   qm2ds%Ncis=qm2ds%Np*qm2ds%Nh

   qm2ds%Nrpa=2*qm2ds%Ncis
   qm2ds%Mx=min(qm2ds%Ncis,qm2ds%Mx) ! Mx has to be smaller than Ncis
   ! Nqat - number of quantum atoms (nlink?)
   qm2ds%Nqat=qmmm_struct%nquant_nlink

   ! Maximum dimension of Krylov expansion in davidson
   qm2ds%nd=5*qm2ds%Mx ! empirical factor of 5
   ! apparently this is not enough
   qm2ds%nd=500 ! kav: testing
   
   ! Pointer to ground state Hartree-Fock eigenenergies, (1:Nb)
   qm2ds%ehf=>qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,1)

   ! Pointer to ground state Hartree-Fock eigenvectors
   qm2ds%vhf=>qm2_struct%eigen_vectors

   ! Pointer to all the Coulomb matrix elements, (:qm2_struct%n2el)
   qm2ds%W=>qm2_struct%qm_qm_2e_repul

   ! Verbosity
   qm2ds%verbosity=qmmm_nml%verbosity

!-----------------BEGIN ALLOCATION-----------------------
   ! Transition energies - Davidson eigenvalues
   allocate(qm2ds%e0(qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)

   ! Total energies of excited states
   allocate(qm2ds%Etot(qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)

   ! Gradient of required (last) excited state
   if(.not.allocated(qm2ds%grad)) then
      allocate(qm2ds%grad(3*qm2ds%Nqat),stat=ierr)
      REQUIRE(ierr==0)
   end if

   ! Transition densities - Davidson eigenvectors
   allocate(qm2ds%v0(qm2ds%Nrpa,qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%v0_old(qm2ds%Nrpa,qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%v2(qm2ds%Nb**2,qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)

   ! Various auxiliary arrays
   allocate(qm2ds%ferr(qm2ds%Mx),stat=ierr) ! convergence error of found modes 
   REQUIRE(ierr==0)
   allocate(qm2ds%rrwork(4*qm2ds%Nrpa),stat=ierr) ! Lanczos vectors
   REQUIRE(ierr==0)
   allocate(qm2ds%temp1(qm2ds%Nrpa),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%temp2(qm2ds%Nrpa),stat=ierr)
   REQUIRE(ierr==0)

   ! Spare space
   allocate(qm2ds%eta(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%xi(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%etas(qm2ds%Lt),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%xis(qm2ds%Lt),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%xi_tz(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%eta_tz(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%tz_scratch(2*qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%xi_scratch(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%cmdqt(2*qm2ds%Np*qm2ds%Nh, qm2ds%Mx),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%eta_scratch(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%xi_scratch_2(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%eta_scratch_2(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%nacr_scratch(qm2ds%Nb*(qm2ds%Nb+1)),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%dij(qmmm_struct%nquant_nlink*3),stat=ierr)
   REQUIRE(ierr==0)


   ! Davidson expansion vectors
   allocate(qm2ds%vexp1(qm2ds%Ncis*qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%vexp(qm2ds%Nrpa*qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)

   ! Davidson Rayleigh matrices
   allocate(qm2ds%ray1(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%ray2(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%ray1a(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%ray(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%raye(qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%raye1(qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%rayv(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%rayvL(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%rayvR(qm2ds%nd,qm2ds%nd),stat=ierr)
   REQUIRE(ierr==0)

   ! index arrays
   allocate(qm2ds%ix(qm2ds%Ncis),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%jx(qm2ds%Ncis),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%kx(qm2ds%Mx),stat=ierr)  
   REQUIRE(ierr==0)
 
   ! Density matrix information
   allocate(qm2ds%rhoTZ(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%rhoT(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)
   allocate(qm2ds%rhoLZ(qm2ds%Nb**2),stat=ierr)
   REQUIRE(ierr==0)

   ! Excited state gradient information
   allocate (qm2ds%dxyz(1:3*qm2ds%Nqat),stat=ierr)
   REQUIRE(ierr==0)


   return
   end subroutine allocate_davidson
!
!********************************************************************
!
!  Davidson deallocation
!
!********************************************************************
!
   subroutine deallocate_davidson()

	use qmmm_module, only: qmmm_struct

	integer :: ierr

   print*,'davidson deallocation'

   if(.NOT. initialized) then
      write(6,*) ' Davidson was never initialized. Exiting deallocation procedure'
      return
   end if
   
   ! "de-pointing"
   qm2ds%ehf=>null()
   qm2ds%W=>null()

	! Deallocation

   deallocate(qm2ds%e0, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%Etot, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%grad, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%v0, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%v0_old, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%v2, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ferr, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rrwork, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%temp1, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%temp2, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%eta, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%xi, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%etas, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%xis, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%xi_tz, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%eta_tz, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%tz_scratch, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%xi_scratch, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%cmdqt, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%eta_scratch, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%xi_scratch_2, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%eta_scratch_2, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%nacr_scratch, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%dij, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%vexp1, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%vexp, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ray1, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ray2, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ray1a, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ray, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%raye, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%raye1, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rayv, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rayvL, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rayvR, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%ix, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%jx, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%kx, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rhoTZ, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rhoT, stat=ierr)
   REQUIRE(ierr==0)
   deallocate(qm2ds%rhoLZ, stat=ierr)
   REQUIRE(ierr==0)


   return
   end subroutine deallocate_davidson

   end module qm2_davidson_module
!
