! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-
#include "assert.fh"
#include "dprec.fh"

module qm2_dftb_module

!use qmmm_module, only: qmmm_struct
!! 
!! (GMS) DFTB Module defines the arrays used in DFTB.
!! ==================================================
!!
!! Currently, the Common blocks used in DFTB were all 
!! converted into structures. Also, all arrays from DFTB 
!! are now defined here to ba allocated dinamically.
!!


!! Static paramters

   ! These parameters define the format of the integral table files, and
   ! should not change.
   integer, parameter :: MAXINT = 250       ! Maximum number of intervals in spline files
   integer, parameter :: MAXTAB = 700       ! Maximum number of integrals in SK files
   integer, parameter :: LDIM   = 9         ! Maximum number of orbitals per atom (1s + 3p + 5d)

   ! Parameters needed by the Broyden mixing, to be set in qm2_dftb_load_params. 
   integer :: MAX_BRD_ITER
   integer :: IMATSZ      

   ! Parameters that are supposed to be different for every run - filled in qm2_dftb_load_params
   integer :: NNDIM    ! (NNDIM  = 100) Maximum number of quantum atoms
   integer :: MDIM     ! (MDIM   = 400) Maximum number of orbitals
   integer :: MAXSIZ   ! (MAXSIZ = MDIM)

   integer :: NDIM     ! To be determined in load_params, as the actual size of H and S matrices
   logical :: do_scf   ! .false. -> do non-SCC DFTB. Default: .true.

!! Arrays
   _REAL_ , dimension(:,:), pointer :: dummy     ! Used as scratch space in slkode <qm2_dftb_slkide.f.
   _REAL_ , dimension(:),   pointer :: scr_space ! (3*MDIM), used in ewevge

!!-------------------------------------------------------
!! Old common blocks from DFTB, converted into structures
!!-------------------------------------------------------

   ! common /disper/ Edis, dispers
   ! used by:
   !         dispersion_egr       (qm2_dftb_dispersion_egr.f)
   !         dylcao               (qm2_dftb_dylcao.f)
   !         qm2_dftb_load_params (qm2_dftb_load_params.f)
   !         qm2_dftb_main        (qm2_dftb_main.f)
   !         eglcao               (qm2_dftb_eglcao.f)
   !         
   type disper_structure
      _REAL_  :: Edis    ! Dispersion energy
      ! DEPRECATED. Substituted by qmmm_nml%dftb_disper == 1(T) or 0(F, default):
      ! logical :: dispers ! Do dispersion? 
   end type disper_structure
   type (disper_structure) disper


   ! common /dispertmp/ A,B,C,r0,rv,C6,Rvdw
   ! Parameters for dispersion calculations.
   ! Used by:
   !         dispersion_egr       (qm2_dftb_dispersion_egr.f)
   !         dis_e                (qm2_dftb_dispersion_egr.f)
   !         dis_gr               (qm2_dftb_dispersion_egr.f)
   !         dispersion_params    (qm2_dftb_dispersion_params.f)
   !         dispersionread       (qm2_dftb_dispersionread.f)
   !
   type dispertmp_structure
      _REAL_, dimension(:,:), pointer :: C6,Rvdw
      _REAL_ :: A,B,C,r0,rv
   end type dispertmp_structure
   type (dispertmp_structure), target :: dispertmp

   !  Common /dispfile/ has the information from the dispersion parameters file.
   !  common /dispfile/ h1,h2,Ni0,scale,read_DISP_DOT_INP
   !  Used by:
   !          dispersion_params  (qm2_dftb_dispersion_params.f)
   !          dispersionread     (qm2_dftb_dispersionread.f)
   !
   type dispfile_structure
      _REAL_, dimension(:)  , pointer :: Ni0 ! Ni0(MAXTYP)
      _REAL_, dimension(:,:), pointer :: h1  ! h1(MAXTYP,4)
      _REAL_, dimension(:,:), pointer :: h2  ! h2(MAXTYP,4)

      _REAL_  :: scale
      logical :: read_DISP_DOT_INP

      ! The following are used only by dispersion_params <qm2_dftb_dispersion_params.f> :
      integer, dimension(:)  , pointer :: nei   ! (NNDIM)
      _REAL_ , dimension(:)  , pointer :: Ni    ! (NNDIM)
      _REAL_ , dimension(:)  , pointer :: hh1   ! (NNDIM)
      _REAL_ , dimension(:)  , pointer :: hh2   ! (NNDIM)
      _REAL_ , dimension(:,:), pointer :: R0vdw ! (MAXTYP,4)
      
   end type dispfile_structure
   type (dispfile_structure), target :: dispfile


   ! Third-order SCC-DFTB 
   type DFTB_3rd_order_structure
      _REAL_ :: Gaussian_D0
      _REAL_ :: Gaussian_G0
      _REAL_ :: Gaussian_Q0
      _REAL_, dimension(:), pointer :: Hubbard_deriv
      integer :: file_unit
      logical :: do_3rd_order ! Default: .false.
      logical :: debug_print  ! Default: .false.
      character(len=1024) :: file_name
   end type DFTB_3rd_order_structure
   type (DFTB_3rd_order_structure), target :: DFTB_3rd_order_str

   ! Common /ctrl/
   ! common /ctrl/ fmax,scftol,deltat,tatom,telec,wvscale, &
   !                maxcyc,mode, &
   !                chrr,scfhelp, &
   !                outfile
   ! Used by:
   !         dylcao                (qm2_dftb_dylcao.f)
   !         qm2_dftb_load_params  (qm2_dftb_load_params.f)
   !

   ! common /mol/ x,qmat,nn,ntype,atyp,atnames
   ! Used by:
   !         dylcao               (qm2_dftb_dylcao.f)
   !         qm2_dftb_load_params (qm2_dftb_load_params.f)
   !         qm2_dftb_main        (qm2_dftb_main.f)
   !
   type mol_structure
      _REAL_, dimension(:), pointer :: qmat
      character(len=2), dimension(:), pointer :: atyp ! (MAXTYP)
   end type mol_structure
   type (mol_structure) mol

   ! common /sktab/ skhtab,skstab,skself,sr,dimens
   ! Information from the S-K tables about
   ! Hamiltonian and Overlap integrals.
   ! Used by:
   !           dylcao          (qm2_dftb_dylcao.f)
   !           gettab          (qm2_dftb_gettab.f)
   !  FUNCTION skspar          (qm2_dftb_skpar.f)
   !  FUNCTION skhpar          (qm2_dftb_skpar.f)
   !
   type sktab_structure

      character(len=1024), dimension(:,:), pointer :: skfiles ! (MAXTYP,MAXTYP)
      character(len=2 ), dimension(:)  , pointer :: latyp   ! (MAXTYP)
      _REAL_ , dimension(:,:), pointer :: sr     ! Step widht
      integer, dimension(:,:), pointer :: dimens ! Number of distances stored in SK file
      _REAL_ , dimension(:,:), pointer :: skself ! d, p and s (self) energies.
      _REAL_ , dimension(:,:,:,:), pointer  :: skhtab ! Parameters for H
      _REAL_ , dimension(:,:,:,:), pointer  :: skstab ! Parameters for S

      _REAL_  :: slkcutoff                        ! cutoff
      _REAL_  :: slkcutoff2    
   end type sktab_structure
   type (sktab_structure), target :: sktab


   ! common /spltab/ coeff,xr,efkt,cutoff,numint
   ! Spline data for repulsive force
   ! Used by:
   !         gettab 
   !         FUNCTION repen   (qm2_dftb_repulsiv.f)
   !         FUNCTION grdrep  (qm2_dftb_repulsiv.f)
   type spltab_structure

      integer, dimension(:,:)    , pointer :: numint      ! Number of splines (intervals)
      _REAL_ , dimension(:,:)    , pointer :: cutoff      ! Spline cutoff
      _REAL_ , dimension(:,:,:)  , pointer :: efkt        ! Short-distance exponential parameters
      _REAL_ , dimension(:,:,:,:), pointer :: xr          ! Spline intervals (1=begin, 2=end)
      _REAL_ , dimension(:,:,:,:), pointer :: coeff       ! Spline coefficients

   end type spltab_structure
   type (spltab_structure) spltab


   ! common /spin/ espin
   ! Atomic spin polarization energy
   ! Used by:
   !         dylcao                 (qm2_dftb_dylcao.f)
   !         gettab                 (qm2_dftb_gettab.f)
   !_REAL_  :: espin(MAXTYP)
   _REAL_, dimension(:), pointer  :: espin ! Spin polarization energy

   ! Common /lmax/ --> Maximum number of orbitals
   ! common /lmax/ lmax
   ! Used by:
   !         dylcao                 (qm2_dftb_dylcao.f)
   !         eglcao                 (qm2_dftb_eglcao.f)
   !         qm2_dftb_load_params   (qm2_dftb_load_params.f)
   !         qm2_dftb_main          (qm2_dftb_main.f)
   !         FUNCTION skpar         (qm2_dftb_skpar.f)
   !         slkmatrices            (qm2_dftb_slkode.f)
   !         slkode                 (qm2_dftb_slkode.f)
   integer, dimension(:), pointer :: lmax    !! lmax(qmmm_struct%qm_ntypes): Max # of orbitals for type "ntype"

   ! common /machine/ dacc           Machine precision*4.0d0
   !                  log_racc       log machine precision
   ! Used by:
   !         dylcao
   !         eglcao
   !         fermi                   (qm2_dftb_fermi.f)
   _REAL_  :: dacc
   _REAL_  :: log_racc


   ! common /mcharge/ qzero, uhubb  :: Charge parameters
   ! Used by:
   !         dylcao
   !         eglcao
   !         externalchgrad          (qm2_dftb_externalchgrad.f)
   !         externalshift           (qm2_dftb_externalshift.f)
   !         gettab                  (qm2_dftb_gettab.f)
   !
   type mcharge_structure

      _REAL_, dimension(:), pointer :: qzero
      _REAL_, dimension(:), pointer :: uhubb
   end type mcharge_structure
   type (mcharge_structure) mcharge

   ! common /izp/ izp, nel, nbeweg    !!  More moleular info...
   ! Used by:
   !         dylcao
   !         eglcao
   !         qm2_dftb_load_params  (qm2_dftb_load_parmas.f)
   !         qm2_dftb_main         (qm2_dftb_main.f)
   !         slkmatrices           (qm2_dftb_slkode.f)
   !         outeigenvectors       (qm2_dftb_output.f)
   type izp_structure
      _REAL_  :: nel            !! CHARGE on the qm region (BAD NAMING!!)
      integer, dimension(:), pointer :: izp  !! (NNDIM) =>  atom types
   end type izp_structure
   type (izp_structure) izp_str

!!------------------------
!! Locally defined arrays
!!------------------------

   ! Arrays defined in eglcao
   ! These arrays are used in the solution of the 
   ! Kohn-Sham equations for DFTB
   type ks_dftb_structure


      integer, dimension(:)  , pointer :: ind      ! (NNDIM+1)     => Indices in the H and S matrices
      _REAL_ , dimension(:,:), pointer :: hgrad    ! (3,NNDIM)     => receipient for the gradients from each part.
      _REAL_ , dimension(:,:), pointer :: au       ! (LDIM,LDIM)   => H data from SK matrices
      _REAL_ , dimension(:,:), pointer :: bu       ! (LDIM,LDIM)   => S data from SK matrices
      _REAL_ , dimension(:,:), pointer :: auh      ! (LDIM,LDIM)   => H from numerical derivative step (grad)
      _REAL_ , dimension(:,:), pointer :: buh      ! (LDIM,LDIM)   => S from numerical derivative step (grad)
      _REAL_ , dimension(:,:), pointer :: a        ! (MDIM,MDIM)   => Upper triangle of H (for eigenvector solver)
      _REAL_ , dimension(:,:), pointer :: b        ! (MDIM,MDIM)   => Upper triangle of S (for eigenvector solver)
      _REAL_ , dimension(:,:), pointer :: gammamat ! (NNDIM,NNDIM) => Scratch space (?) for gammamatrix (in hamilshift)
      _REAL_ , dimension(:,:), pointer :: derivx   ! (NNDIM,NNDIM) => SCRATCH SPACE Used only in gammagrad call
      _REAL_ , dimension(:,:), pointer :: derivy   ! (NNDIM,NNDIM) => SCRATCH SPACE Used only in gammagrad call
      _REAL_ , dimension(:,:), pointer :: derivz   ! (NNDIM,NNDIM) => SCRATCH SPACE Used only in gammagrad call
      _REAL_ , dimension(:)  , pointer :: ev       ! (MDIM)        => Orbital energies (eigenvalues)
      _REAL_ , dimension(:)  , pointer :: occ      ! (MDIM)        => Occupation numbers
      _REAL_ , dimension(:)  , pointer :: qmulli   ! (MDIM)        => Electron population per atomic orbital
      _REAL_ , dimension(:)  , pointer :: qmold    ! (NNDIM)       => Old mulliken charges (for SCC)
      _REAL_ , dimension(:,:), pointer :: hamil    ! (MDIM,MDIM)   => hamiltonian matrix
      _REAL_ , dimension(:,:), pointer :: overl    ! (MDIM,MDIM)   => Overlap matrix
      _REAL_ , dimension(:,:), pointer :: density  ! (NDIM,NDIM)   => Density matrix
      _REAL_ , dimension(:)  , pointer :: shift    ! (NNDIM)       => Energy shift due to SCC (gamma)
      _REAL_ , dimension(:)  , pointer :: shiftE   ! (NNDIM)       => Energy shift due to external charges


      _REAL_ , dimension(:,:), pointer :: scr1     ! (MDIM,MDIM)   => Overlap matrix
      _REAL_ , dimension(:,:), pointer :: scr2     ! (MDIM,MDIM)   => Overlap matrix
      _REAL_ , dimension(:,:), pointer :: xtrans   ! (MDIM,MDIM)   => Overlap matrix

      
      _REAL_ :: dipol(3)                               ! Dipole vector

   end type ks_dftb_structure
   type (ks_dftb_structure), target :: ks_struct

   type broyden_strucure
      _REAL_, dimension(:)    , pointer :: f       ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: ui      ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: vti     ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: t1      ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: dumvi   ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: df      ! (maxsiz)
      _REAL_, dimension(:)    , pointer :: cm      ! (imatsz)
      _REAL_, dimension(:)    , pointer :: w       ! (imatsz)
      _REAL_, dimension(:,:)  , pointer :: a       ! (imatsz,imatsz)
      _REAL_, dimension(:,:)  , pointer :: b       ! (imatsz,imatsz)
      _REAL_, dimension(:,:)  , pointer :: d       ! (imatsz,imatsz)
      _REAL_, dimension(:,:)  , pointer :: unit31  ! (maxsiz,2)
      _REAL_, dimension(:,:,:), pointer :: unit32  ! (maxsiz,2,MAX_BRD_ITER)

      _REAL_, dimension(:)    , pointer :: td      ! (imatsz)
      _REAL_, dimension(:)    , pointer :: ad      ! (imatsz)
      _REAL_, dimension(:)    , pointer :: bd      ! (imatsz)

   end type broyden_strucure
   type (broyden_strucure), target :: brd_struct

   ! Structure for CM3 charges
   type cm3_structure
      integer :: num_params = 15                   ! Number of atom pairs available [JPC-A 2004, 108, 2545]
      _REAL_ , dimension(:)  , pointer :: qcm3     ! (natoms) Store cm3 charges
      _REAL_ , dimension(:,:), pointer :: c        ! (ntype, ntype), CM3 C parameters
      _REAL_ , dimension(:,:), pointer :: d        ! (ntype, ntype), CM3 D parameters
      _REAL_ , dimension(:,:), pointer :: t        ! (natoms,natoms), CM3 T matrix
      _REAL_ , dimension(:,:), pointer :: b        ! (natoms,natoms), Mayer bond order matrix
   end type cm3_structure
   type (cm3_structure), target, save :: cm3

   ! Structure for Fermi distributions
   type fermi_structure
      _REAL_ :: telec      ! Electronic Temperature
      _REAL_ :: telec_step ! Step size when increasing telec for convergence.
   end type fermi_structure
   type (fermi_structure), target, save :: fermi_str

end module qm2_dftb_module

