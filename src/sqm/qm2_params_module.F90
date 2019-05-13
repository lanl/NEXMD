#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
module qm2_params_module
! ----------------------------------------------------------------------
! PURPOSE: Data type holding semiempirical parameters and
!          pre-computed multipole expansion terms and integrals
! 
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : October 2011
! ----------------------------------------------------------------------
  
  use ElementOrbitalIndex, only :NumberElements

  implicit none

  private
  public :: qm2_params_type
  public :: new, delete

  type qm2_params_type

     ! Calculated in qm2_load_params - independent of structure so constant for a sander run.
     _REAL_ :: tot_heat_form

     integer, dimension(:), pointer :: sp_quantum_number ! qm_ntypes
     integer, dimension(:), pointer :: d_quantum_number  ! qm_ntypes  
   
     _REAL_, dimension(:), pointer :: gss ! qm_ntypes
     _REAL_, dimension(:), pointer :: hsp ! qm_ntypes
     _REAL_, dimension(:), pointer :: hpp ! qm_ntypes  

     _REAL_, dimension(:,:), pointer :: dd ! 6 x ntypes
     _REAL_, dimension(:,:), pointer :: po ! 9 x ntypes    
   
     ! The core charge on each atom as seen by the electrons
     _REAL_, dimension(:), pointer :: core_chg 

     ! Orbital electron kinetic energy integrals (2,nquant_nlink) (1=s, 2=p)
     _REAL_, dimension(:,:), pointer :: orb_elec_ke 

     _REAL_, dimension(:,:), pointer :: betasas !betas(ni)+betas(nj) qm_ntypes x qm_ntypes
     _REAL_, dimension(:,:), pointer :: betasap !betas(ni)+betap(nj) qm_ntypes x qm_ntypes
     _REAL_, dimension(:,:), pointer :: betasad !betas(ni)+betad(nj) qm_ntypes x qm_ntypes
     _REAL_, dimension(:,:), pointer :: betapap !betap(ni)+betap(nj) qm_ntypes x qm_ntypes
     _REAL_, dimension(:,:), pointer :: betapad !betap(ni)+betad(nj) qm_ntypes x qm_ntypes
     _REAL_, dimension(:,:), pointer :: betadad !betad(ni)+betad(nj) qm_ntypes x qm_ntypes 

     _REAL_, dimension(:), pointer :: GNN, rho_core      

     ! Pre-computed Slater-Condon parameters F0SD and G2SD
     ! Required for some elements for PM6
     _REAL_, dimension(:), pointer :: F0SD ! qm_ntypes long
     _REAL_, dimension(:), pointer :: G2SD ! qm_ntypes long

     ! PM3 / AM1 specific parameters for core-core repulsions. (also RM1, PM6 etc)     
     _REAL_, dimension(:,:), pointer :: FN1     
     _REAL_, dimension(:,:), pointer :: FN2
     _REAL_, dimension(:,:), pointer :: FN3


     !Coulomb and exchange one centre-two electron integral params.
     _REAL_, dimension(:,:), pointer :: onec2elec_params

     ! Parameters for the multipole expansion of the 2 centre 2 electron integerals.
     ! 9, nquant_nlink in order DD,QQ,AM,AD,AQ,AM2,AD2,AQ2
     _REAL_, dimension(:,:), pointer :: multip_2c_elec_params

     ! Exponents for core core repulsions
     _REAL_, dimension(:), pointer :: cc_exp_params

     ! Exponents for pairwise core core repulsions
     _REAL_, dimension(:,:), pointer :: pm6_alpab
     ! Coefficients for pairwise core core repulsions
     _REAL_, dimension(:,:), pointer :: pm6_xab
  
     !PM3-MAIS: Coefficients for pairwise core core repulsions.
     _REAL_, dimension(:,:,:), pointer :: pm3mais_alpab
     !PM3-MAIS: Exponents for pairwise core core repulsions.   
     _REAL_, dimension(:,:,:), pointer :: pm3mais_betab
     !PM3-MAIS: Coefficients for pairwise core core repulsions.
     _REAL_, dimension(:,:,:), pointer :: pm3mais_gamab

     ! The arrays for s and p functions are deallocated by qm2_setup_orb_exp as it is not needed after this is called
     ! S/P/D orbital expansion coefficients for the Slater orbital expansion
     _REAL_, dimension(:), pointer :: s_orb_exp_by_type
     _REAL_, dimension(:), pointer :: p_orb_exp_by_type
     _REAL_, dimension(:), pointer :: d_orb_exp_by_type

     ! S/P/D orbital expansion coefficients for the Slater orbital expansion
     _REAL_, dimension(:), pointer :: s_orb_exp_tail_by_type
     _REAL_, dimension(:), pointer :: p_orb_exp_tail_by_type
     _REAL_, dimension(:), pointer :: d_orb_exp_tail_by_type  

     !Arrays for PDDG Hamiltonians
     _REAL_, dimension(:), pointer :: pddge1, pddge2

     ! Two params used in PM3/MM* calculations.
     _REAL_, dimension(:,:), pointer :: scale_factor1_pm3mmx
     _REAL_, dimension(:,:), pointer :: scale_factor2_pm3mmx
     _REAL_, dimension(:), pointer :: rho_pm3mmx

     ! --------------------------------------------
     ! Arrays for pre-computed orbital interactions
     ! --------------------------------------------
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_sxs_over_sas !atom_orb_zz_s_x_s/atom_orb_zz_one_s_a_s
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_sxp_over_sap !atom_orb_zz_s_x_p/atom_orb_zz_one_s_a_p
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_sxd_over_sad !atom_orb_zz_s_x_d/atom_orb_zz_one_s_a_d
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_pxp_over_pap !atom_orb_zz_p_x_p/atom_orb_zz_one_p_a_p
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_pxd_over_pad !atom_orb_zz_p_x_d/atom_orb_zz_one_p_a_d 
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_dxd_over_dad !atom_orb_zz_d_x_d/atom_orb_zz_one_d_a_d 

     ! sqrt((two*sqrt(atom_orb_zz_s_x_s)*atom_orb_zz_one_s_a_s)**3*atom_orb_cc_s_x_s
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_ss_eqn 

     ! 2.0D0 * qm2_params%atom_orb_zz_s_by_type(K,qmitype)* SQRT(qm2_params%atom_orb_zz_p_by_type(L,qmjtype))*
     ! atom_orb_zz_one_s_a_p(k,l,qmitype,qmjtype)*atom_orb_sp_eqn
     ! Used in gover for Si-Pj overlap energy and -Sj-Pi overlap.
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_ovlp
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sd_ovlp, atom_orb_pd_ovlp

     ! -4.0D0*sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p
     ! *qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)*atom_orb_pp_eqn
     ! Used in gover for Pi-Pj overlap energy when i!=j  
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_inj
  
     ! -4.0D0*atom_orb_pp_eqn*sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p*
     ! qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)
     ! Used in gover for Pi-Pj overlap energy when i==j
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_ieqj1 

     ! 2.0D0*atom_orb_pp_eqn*sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p  
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_ieqj2
  
     ! -4.0D0*sqrt(atom_orb_zz_d_x_d)*atom_orb_zz_one_d_a_d
     ! *qm2_params%atom_orb_zz_dxd_over_dad(k,l,qmitype,qmjtype)*atom_orb_dd_eqn
     ! Used in gover for di-dj overlap energy when i!=j
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_dd_ovlp_inj 
  
     ! -4.0D0*atom_orb_dd_eqn*sqrt(atom_orb_zz_d_x_d)*atom_orb_zz_one_d_a_d*
     ! qm2_params%atom_orb_zz_dxd_over_dad(k,l,qmitype,qmjtype)
     ! Used in gover for di-dj overlap energy when i==j
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_dd_ovlp_ieqj1 
  
     ! 2.0D0*atom_orb_dd_eqn*sqrt(atom_orb_zz_d_x_d)*atom_orb_zz_one_d_a_d
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_dd_ovlp_ieqj2 
     
     ! -2.0D0*A2_TO_BOHRS2*qm2_params%atom_orb_ss_eqn_adb(i,j,qmitype,qmjtype)*
     ! qm2_params%atom_orb_zz_sxs_over_sas(i,j,qmitype,qmjtype)
     ! Used for S-S overlap in QM-QM derivatives
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_ss_eqn_adb 
  
     ! -four*A3_TO_BOHRS3*qm2_params%atom_orb_zz_sxp_over_sap(i,j,qmitype,qmjtype)**2* &
     ! (one/(SQRT(qm2_params%atom_orb_zz_p_by_type(J,qmjtype))))*atom_orb_sp_eqn
     ! Used for S-P overlap in QM-QM derivatives where P... /= axis
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xy  

     ! Used for S-P overlap in QM-QM derivatives where P... == axis
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xx1
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xx2 
  
     ! -four*A2_TO_BOHRS2*(ADB_array(inner_index)**2)*(one/(SQRT(atom_orb_zz_p_x_p)))*atom_orb_pp_eqn
     ! Used for P-P overlap in QM-QM derivatives where P... = P... /= axis
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_eqn_xxy1 

     ! eight*A2_TO_BOHRS2*A2_TO_BOHRS2*(ADB_array(inner_index)**3)*(one/(SQRT(atom_orb_zz_p_x_p)))*atom_orb_pp_eqn
     _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_eqn_xxy2
   
     ! ------------------------------------------------
     ! End arrays for pre-computed orbital interactions
     ! ------------------------------------------------
  
     ! Pre-computed PDDG parameters if PDDG Hamiltonian is in use.
     ! Ntypes*ntypes - Stores the pre-exponential part of the PDDG equation.
     _REAL_, dimension(:,:), pointer :: pddg_term1, pddg_term2, pddg_term3, pddg_term4

     ! Number of atomic orbitals on atom.
     integer, dimension(:), pointer :: natomic_orbs
     ! Locations of orbitals. 2,nquant_nlink. 1,x gives beginning of
     ! of orbitals on atom x. 2,x gives last orbital on atom x.
     integer, dimension(:,:), pointer :: orb_loc

     ! Lower half triangle indices (pascal's triangle)
     integer, dimension(:), pointer :: pascal_tri1 
     integer, dimension(:), pointer :: pascal_tri2

     ! Number of FNX terms (first dimension) that are not zero.
     integer, dimension(:), pointer :: NUM_FN
   
     ! OPNQ parameters
     logical, dimension(:), pointer :: qxd_supported
     _REAL_, dimension(:), pointer :: qxd_s, qxd_z0, qxd_zq, qxd_d0, qxd_dq, qxd_q0, qxd_qq, qxd_neff

     !save datas from MNDOCharegSeparation Module subroutines    
     integer :: DDP_qmtype_saved=-1
     _REAL_  :: DD_saved(6), PO_saved(9)
     logical :: DDP_initialized=.false. 
   
     _REAL_  :: REPD(52)=0.0D0
     logical :: OC2E_initialized=.false.
     integer :: OC2E_qmType_saved=-1

     !save from qm2_scf
     _REAL_  :: smallsum, abstol !Underflow limit for dspevr 
     integer :: itrmax_local

     !save from rotation module
    _REAL_   :: matrix_saved(15,45)
    _REAL_  :: xij_saved(3)=(/ 1.0D9, 1.0D9, 1.0D9 /) 
    
     !save from Vxi_pack
     logical :: vxi_first=.true.
     logical :: vxia_first=.true.
     character :: keywr*6    
 
  end type qm2_params_type

  interface new
     module procedure allocate_qm2_params_type
  end interface

  interface delete
     module procedure deallocate_qm2_params_type
  end interface

contains
  
  subroutine allocate_qm2_params_type(self, qm_ntypes, nquant_nlink, qmtheory, PM3MMX_INTERFACE, useOPNQ)

    use qmmm_qmtheorymodule, only : qmTheoryType
    use ElementOrbitalIndex, only : MaxAngularQuantumNumber
    use constants, only : zero

    implicit none

    type(qm2_params_type), intent(out) :: self
    integer, intent(in) :: qm_ntypes
    integer, intent(in) :: nquant_nlink
    type(qmTheoryType), intent(in) :: qmtheory
    logical, intent(in) :: PM3MMX_INTERFACE
    logical, intent(in) :: useOPNQ

    integer :: ier

    ! Memory allocation - Needs cleaning up for DFTB at some point.
    allocate ( self%sp_quantum_number(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
    allocate ( self%d_quantum_number(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)  
    allocate ( self%gss(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
    allocate ( self%hsp(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
    allocate ( self%hpp(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)                  
 
    allocate ( self%GNN(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
    self%GNN=1.0D0
 
    allocate ( self%rho_core(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
    self%rho_core = zero
      
    allocate ( self%dd(6, qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%po(9, qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 

    allocate ( self%F0SD(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    self%F0SD = zero
    allocate ( self%G2SD(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    self%G2SD = zero

    allocate ( self%core_chg(nquant_nlink), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%natomic_orbs(nquant_nlink), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%orb_loc(2,nquant_nlink), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%onec2elec_params(5,nquant_nlink), stat=ier )
    REQUIRE(ier == 0)
    allocate ( self%multip_2c_elec_params(6,nquant_nlink), stat=ier )
    REQUIRE(ier == 0)

    if (qmtheory%PM6) then
       !Allocation of pairwise core core terms (PM6)
       allocate ( self%pm6_alpab(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pm6_xab(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
    elseif (qmtheory%PM3MAIS) then
       !Allocation of pairwise core core terms (PM3MAIS)
       allocate ( self%pm3mais_alpab(qm_ntypes,qm_ntypes,3), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pm3mais_betab(qm_ntypes,qm_ntypes,3), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pm3mais_gamab(qm_ntypes,qm_ntypes,3), stat=ier )
       REQUIRE(ier == 0)
    end if
    allocate ( self%cc_exp_params(nquant_nlink), stat=ier )
    REQUIRE(ier == 0)

    allocate (self%orb_elec_ke(MaxAngularQuantumNumber+1,nquant_nlink), stat=ier )
    REQUIRE(ier == 0)
    allocate ( self%betasas(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%betasap(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%betasad(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%betapap(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%betapad(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate ( self%betadad(qm_ntypes,qm_ntypes), stat=ier )
    REQUIRE(ier == 0)
      
    ! OPNQ parameters
    if (useOPNQ) then
       allocate ( self%qxd_supported(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)      
       allocate ( self%qxd_s(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)   
       allocate ( self%qxd_z0(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%qxd_zq(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%qxd_d0(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)        
       allocate ( self%qxd_dq(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%qxd_q0(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)     
       allocate ( self%qxd_qq(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)  
       allocate ( self%qxd_neff(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0)                                                                
    end if
            
    if (qmtheory%AM1 .OR. qmtheory%AM1D .or.  &
         qmtheory%PM3 .OR. qmtheory%PDDGPM3 .OR. &
         qmtheory%PM3CARB1 .OR. qmtheory%RM1 .OR. qmtheory%PDDGPM3_08 .OR. &
         qmtheory%PM6 .OR. qmtheory%PM3ZNB .OR. qmtheory%PM3MAIS) then
       allocate ( self%FN1(4,qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%FN2(4,qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%FN3(4,qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
       allocate ( self%NUM_FN(qm_ntypes), stat=ier ) 
       REQUIRE(ier == 0) 
    end if

    if (qmtheory%PDDGPM3 .OR. qmtheory%PDDGMNDO .OR. qmtheory%PDDGPM3_08) then
       allocate ( self%pddge1(nquant_nlink), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pddge2(nquant_nlink), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pddg_term1(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pddg_term2(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pddg_term3(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate ( self%pddg_term4(qm_ntypes,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
    end if

    ! PM3/MM*
    if (PM3MMX_INTERFACE) then
       !should really be based on atom types.
       allocate (self%scale_factor1_pm3mmx(2,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate (self%scale_factor2_pm3mmx(2,qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
       allocate (self%rho_pm3mmx(qm_ntypes), stat=ier )
       REQUIRE(ier == 0)
    end if

    ! Note: s_orb_exp and p_orb_exp are only needed on the first call to fill 
    ! qm2_setup_orb_exp so after they are used in qm2_setup_orb_exp they are deallocated.
    ! change to type-based index--TL
    allocate (self%s_orb_exp_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate (self%p_orb_exp_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate (self%d_orb_exp_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)     
 
    allocate (self%s_orb_exp_tail_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate (self%p_orb_exp_tail_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0) 
    allocate (self%d_orb_exp_tail_by_type(qm_ntypes), stat=ier )
    REQUIRE(ier == 0)  
        
  end subroutine allocate_qm2_params_type


  ! Note: Part of the arrays get deallocated somewhere else
  subroutine deallocate_qm2_params_type(self, qmtheory, PM3MMX_INTERFACE)

    use qmmm_qmtheorymodule, only : qmTheoryType

    implicit none

    type(qm2_params_type), intent(inout) :: self
    type(qmTheoryType), intent(in) :: qmtheory
    logical, intent(in) :: PM3MMX_INTERFACE

    integer :: ier

    if (.not. qmtheory%DFTB) then
       
       deallocate (self%s_orb_exp_by_type, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%p_orb_exp_by_type, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%d_orb_exp_by_type, stat = ier)
       REQUIRE(ier == 0)
         
       deallocate (self%s_orb_exp_tail_by_type, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%p_orb_exp_tail_by_type, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%d_orb_exp_tail_by_type, stat = ier)
       REQUIRE(ier == 0)         
                  
       deallocate (self%atom_orb_pp_eqn_xxy2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_pp_eqn_xxy1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_sp_eqn_xx2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_sp_eqn_xx1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_sp_eqn_xy, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_ss_eqn_adb, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_pp_ovlp_ieqj2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_pp_ovlp_ieqj1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_pp_ovlp_inj, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_sp_ovlp, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_ss_eqn, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_zz_pxp_over_pap, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_zz_sxp_over_sap, stat = ier)
       REQUIRE(ier == 0)
       deallocate (self%atom_orb_zz_sxs_over_sas, stat = ier)
       REQUIRE(ier == 0)

       if (qmtheory%AM1 .OR. qmtheory%PM3 .OR. qmtheory%PDDGPM3 .OR. &
            qmtheory%PM3CARB1 .OR. qmtheory%RM1 .OR. qmtheory%PDDGPM3_08 .OR. &
            qmtheory%PM6 .OR. qmtheory%PM3ZNB .OR. qmtheory%PM3MAIS) then
          deallocate (self%NUM_FN, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%FN3, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%FN2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%FN1, stat = ier)
          REQUIRE(ier == 0)
       end if

       if (qmtheory%PDDGPM3 .OR. qmtheory%PDDGMNDO .OR. qmtheory%PDDGPM3_08) then
          deallocate (self%pddge1, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%pddge2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%pddg_term1, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%pddg_term2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%pddg_term3, stat = ier)
          REQUIRE(ier == 0)
          deallocate (self%pddg_term4, stat = ier)
          REQUIRE(ier == 0)
       end if

       ! PM3/MM* 
       if (PM3MMX_INTERFACE) then
          deallocate ( self%scale_factor1_pm3mmx, stat = ier)
          REQUIRE(ier == 0)
          deallocate ( self%scale_factor2_pm3mmx, stat = ier)
          REQUIRE(ier == 0)
          deallocate ( self%rho_pm3mmx, stat = ier)
          REQUIRE(ier == 0)
       end if

    end if

    if (qmtheory%PM6) then
       ! Deallocate pairwise core core parameters
       deallocate ( self%pm6_alpab, stat=ier )
       REQUIRE(ier == 0)
       deallocate ( self%pm6_xab, stat=ier )
       REQUIRE(ier == 0)
    end if

    ! PM3-MAIS
    if (qmtheory%PM3MAIS) then
       !Deallocate pairwise core core parameters
       deallocate ( self%pm3mais_alpab, stat=ier )
       REQUIRE(ier == 0)
       deallocate ( self%pm3mais_betab, stat=ier )
       REQUIRE(ier == 0)
       deallocate ( self%pm3mais_gamab, stat=ier )
       REQUIRE(ier == 0)
    end if

    ! Deallocate atom-specific core core parameters
    deallocate (self%cc_exp_params, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%multip_2c_elec_params, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%onec2elec_params, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%orb_elec_ke, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betasas, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betasap, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betasad, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betapap, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betapad, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%betadad, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%orb_loc, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%natomic_orbs, stat = ier)
    REQUIRE(ier == 0)
    deallocate (self%core_chg, stat = ier )
    REQUIRE(ier == 0)
    deallocate (self%sp_quantum_number, stat = ier )
    REQUIRE(ier == 0)  
    deallocate (self%d_quantum_number, stat = ier )
    REQUIRE(ier == 0)  
    deallocate (self%gss, stat = ier )
    REQUIRE(ier == 0) 
    deallocate (self%hsp, stat = ier )
    REQUIRE(ier == 0) 
    deallocate (self%hpp, stat = ier )
    REQUIRE(ier == 0) 
    deallocate (self%dd, stat = ier )
    REQUIRE(ier == 0)   
    deallocate (self%rho_core, stat = ier )
    REQUIRE(ier == 0)   
    deallocate (self%gnn, stat = ier )
    REQUIRE(ier == 0)   
    deallocate (self%po, stat = ier )
    REQUIRE(ier == 0)                                           
    deallocate (self%F0SD, stat = ier )
    REQUIRE(ier == 0)                                           
    deallocate (self%G2SD, stat = ier )
    REQUIRE(ier == 0)                                           
    deallocate (self%pascal_tri1, stat = ier )
    REQUIRE(ier == 0)
    deallocate (self%pascal_tri2, stat = ier )
    REQUIRE(ier == 0)
  
  end subroutine deallocate_qm2_params_type


end module qm2_params_module
