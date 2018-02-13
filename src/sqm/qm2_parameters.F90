! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
!***********************************************************
! Written by: Ross Walker (TSRI 2005)
! Updates by: Andreas Goetz (SDSC 2009)
! Converted into module by: Taisung Lee (Rutgers, 2011) 

! File contains the semi-empirical parameters for use in amber qmmm.

! This file is not intended to be used directly but instead, for speed, is designed
! to be loaded into a parameter array for each atom.
!***********************************************************
module QM2_parameters

  use ElementOrbitalIndex, only :NumberElements
  use qmmm_qmtheorymodule, only :qmTheoryType

  implicit none

  _REAL_, dimension(1:NumberElements) :: heat_of_form
  _REAL_, dimension(1:NumberElements) :: elec_eng_pddgpm3,elec_eng_pddgpm3_08, elec_eng_pddgmndo
  _REAL_, dimension(1:NumberElements) :: s_orb_exp_mndo, s_orb_exp_am1, s_orb_exp_pm3, s_orb_exp_pddgpm3, &
                                    s_orb_exp_pddgpm3_08, s_orb_exp_pddgmndo, s_orb_exp_rm1, &
                                    s_orb_exp_pm6, s_orb_exp_mndod, s_orb_exp_am1d
  _REAL_, dimension(1:NumberElements) :: p_orb_exp_mndo, p_orb_exp_am1, p_orb_exp_pm3, p_orb_exp_pddgpm3, &
                                    p_orb_exp_pddgpm3_08, p_orb_exp_pddgmndo, p_orb_exp_rm1, &
                                    p_orb_exp_pm6, p_orb_exp_mndod, p_orb_exp_am1d
  _REAL_, dimension(1:NumberElements) :: d_orb_exp_mndod, d_orb_exp_am1d, d_orb_exp_pm6 
  _REAL_, dimension(1:NumberElements) :: s_orb_exp_tail_mndod, p_orb_exp_tail_mndod, d_orb_exp_tail_mndod
  _REAL_, dimension(1:NumberElements) :: s_orb_exp_tail_am1d, p_orb_exp_tail_am1d, d_orb_exp_tail_am1d
  _REAL_, dimension(1:NumberElements) :: s_orb_exp_tail_pm6, p_orb_exp_tail_pm6, d_orb_exp_tail_pm6
  _REAL_, dimension(1:NumberElements) :: betas_mndo, betas_am1, betas_pm3, betas_pddgpm3, betas_pddgpm3_08, &
                                    betas_pddgmndo, betas_pm3carb1, betas_rm1, betas_pm6, betas_mndod, betas_am1d 
  _REAL_, dimension(1:NumberElements) :: betap_mndo, betap_am1, betap_pm3, betap_pddgpm3, betap_pddgpm3_08, &
                                    betap_pddgmndo, betap_pm3carb1, betap_rm1, betap_pm6, betap_mndod, betap_am1d
  _REAL_, dimension(1:NumberElements) ::betad_mndod, betad_am1d, betad_pm6  
  _REAL_, dimension(1:NumberElements) ::GNN_am1d
  _REAL_, dimension(1:NumberElements) :: rho_core_am1d, rho_core_mndod, rho_core_pm6
  _REAL_, dimension(1:4,NumberElements) :: FN1_am1, FN2_am1, FN3_am1
  _REAL_, dimension(1:4,NumberElements) :: FN1_am1d, FN2_am1d, FN3_am1d  
  _REAL_, dimension(1:4,NumberElements) :: FN1_rm1, FN2_rm1, FN3_rm1
  _REAL_, dimension(1:4,NumberElements) :: FN1_pm3, FN2_pm3, FN3_pm3
  _REAL_, dimension(1:4,NumberElements) :: FN1_pm6, FN2_pm6, FN3_pm6
  _REAL_, dimension(1:4,NumberElements) :: FN1_pddgpm3, FN2_pddgpm3, FN3_pddgpm3
  _REAL_, dimension(1:4,NumberElements) :: FN1_pddgpm3_08, FN2_pddgpm3_08, FN3_pddgpm3_08
  _REAL_, dimension(1:NumberElements) :: GSS_mndo, GSP_mndo, GPP_mndo, GP2_mndo, HSP_mndo
  _REAL_, dimension(1:NumberElements) :: GSS_mndod, GSP_mndod, GPP_mndod, GDD_mndod, GP2_mndod, HSP_mndod 
  _REAL_, dimension(1:NumberElements) :: GSS_am1, GSP_am1, GPP_am1, GP2_am1, HSP_am1
  _REAL_, dimension(1:NumberElements) :: GSS_am1d, GSP_am1d, GPP_am1d, GDD_am1d, GP2_am1d, HSP_am1d  
  _REAL_, dimension(1:NumberElements) :: GSS_rm1, GSP_rm1, GPP_rm1, GP2_rm1, HSP_rm1
  _REAL_, dimension(1:NumberElements) :: GSS_pm3, GSP_pm3, GPP_pm3, GP2_pm3, HSP_pm3
  _REAL_, dimension(1:NumberElements) :: GSS_pm6, GSP_pm6, GPP_pm6, GP2_pm6, HSP_pm6
  _REAL_, dimension(1:NumberElements) :: GSS_pddgpm3, GSP_pddgpm3, GPP_pddgpm3, GP2_pddgpm3, HSP_pddgpm3
  _REAL_, dimension(1:NumberElements) :: GSS_pddgpm3_08, GSP_pddgpm3_08, GPP_pddgpm3_08, GP2_pddgpm3_08, HSP_pddgpm3_08
  _REAL_, dimension(1:NumberElements) :: GSS_pddgmndo, GSP_pddgmndo, GPP_pddgmndo, GP2_pddgmndo, HSP_pddgmndo
  _REAL_, dimension(1:NumberElements) :: alp_mndo, alp_am1, alp_pm3, alp_pddgpm3, alp_pddgpm3_08, &
                                    alp_pddgmndo, alp_pm3carb1, alp_rm1, alp_pm6, alp_mndod, alp_am1d
! Note: alp_pm6 for PM6 does not contain PM6 parameters but PM3 parameters. For elements which are not parametrized
!       for PM3, the parameters have been guessed by Andreas Goetz (AWG). These parameters are only used for QM/MM
!       core/resp charge interactions.
!       PM6 uses a pair wise core core potential so we need 2D arrays containing both exponent and coefficient. 
  _REAL_, dimension(NumberElements,NumberElements) :: alpab_pm6 !Exponent
  _REAL_, dimension(NumberElements,NumberElements) :: xab_pm6 !Coefficient

! PM3-MAIS arrays. Implementation by the Paesani group
  _REAL_, dimension(NumberElements,NumberElements,1:3) :: alpab_pm3mais
  _REAL_, dimension(NumberElements,NumberElements,1:3) :: betab_pm3mais
  _REAL_, dimension(NumberElements,NumberElements,1:3) :: gamab_pm3mais

  _REAL_, dimension(1:NumberElements) :: USS_mndo, USS_mndod, USS_am1, USS_am1d, USS_pm3, USS_pddgpm3, USS_pddgpm3_08, &
                                    USS_pddgmndo, USS_pm3carb1, USS_rm1, USS_pm6
  _REAL_, dimension(1:NumberElements) :: UPP_mndo, UPP_mndod, UPP_am1, UPP_am1d, UPP_pm3, UPP_pddgpm3, UPP_pddgpm3_08, &
                                    UPP_pddgmndo, UPP_pm3carb1, UPP_rm1, UPP_pm6
  _REAL_, dimension(1:NumberElements) :: UDD_mndod, UDD_am1d, UDD_pm6                                    
  _REAL_, dimension(1:NumberElements) :: PDDGC1_pm3, PDDGC2_pm3, PDDGE1_pm3, PDDGE2_pm3
  _REAL_, dimension(1:NumberElements) :: PDDGC1_pm3_08, PDDGC2_pm3_08, PDDGE1_pm3_08, PDDGE2_pm3_08
  _REAL_, dimension(1:NumberElements) :: PDDGC1_mndo, PDDGC2_mndo, PDDGE1_mndo, PDDGE2_mndo
  _REAL_, dimension(1:2,NumberElements) :: scale_f1_PM3MMX, scale_f2_PM3MMX  ! PM3/MM*
  _REAL_, dimension(1:2,NumberElements) :: scale_f1_PM3MMX2, scale_f2_PM3MMX2  ! PM3/MM* 2nd version
  _REAL_, dimension(1:NumberElements) :: rho_PM3MMX2  ! PM3/MM* 2nd version

  ! Pre-computed Slater-Condon parameters for PM6 (some elements)
  _REAL_, dimension(1:NumberElements) :: F0SD_pm6, G2SD_pm6

  ! The atomic electronic energy (called EISOL in MOPAC and elec_eng in SQM) should be
  ! computed from one-center two-eletron parameters
  ! This is currently only implemented for atoms with s and p orbitals
  ! For d orbitals we currently use EISOL from MOPAC
  ! This should be changed
  _REAL_, dimension(1:NumberElements) :: EISOL_pm6

  integer, dimension(1:NumberElements) :: mndo_ref_index, mndod_ref_index, am1_ref_index, am1d_ref_index, pm3_ref_index, &
                                     pm3carb1_ref_index, rm1_ref_index, pm6_ref_index, pm3znb_ref_index
  integer, dimension(1:NumberElements) :: pddgpm3_ref_index, pddgpm3_08_ref_index, pddgmndo_ref_index
                                     !Index for printing parameter origin in qm_print_ref()
  integer, dimension(1:NumberElements) :: core_chg, ios, iop, nsshell
  integer, dimension(1:NumberElements) :: natomic_orbs
  integer, dimension(1:NumberElements) :: NUM_FN_am1, NUM_FN_am1d, NUM_FN_pm3, NUM_FN_pddgpm3, &
                                     NUM_FN_pddgpm3_08, NUM_FN_rm1, NUM_FN_pm6
                                     !Number of FN arrays that are not zero
                                     
! OPNQ parameters 
  logical, dimension(1:NumberElements) ::qxd_supported
  _REAL_, dimension(1:NumberElements) ::qxd_s, qxd_z0, qxd_zq, qxd_d0, qxd_dq, qxd_q0, qxd_qq, qxd_neff

  logical, dimension(1:NumberElements) :: element_supported_mndo, element_supported_mndod,  &
                                 element_supported_am1, element_supported_am1d,  &
                                 element_supported_pm3, element_supported_pddgpm3, &
                                 element_supported_pddgpm3_08, element_supported_pddgmndo, &
                                 element_supported_rm1, element_supported_pm6, &
                                 element_supported_pm3mais, &
                                 element_supported_pm3mmx,   & ! qmmm_int==3, MODIFIED PM3/MM* 
                                 element_supported_pm3mmx2,  & ! qmmm_int==4, MODIFIED PM3/MM* 2nd version
                                 element_supported_opnq  ! OPNQ
contains

subroutine InitializeParameter(qm2_struct,currentTheory)

  use ParameterReader, only: ParameterEntry, GetNumberParameterEntries, GetParameterEntry
#ifdef MPI
  use qmmm_module, only : qmmm_mpi,qm2_structure
#else 
  use qmmm_module, only : qm2_structure
#endif /* MPI */

    implicit none

#ifdef MPI
    include 'mpif.h'
    integer::ierr
#endif /* MPI */

  type (qm2_structure), intent(inout)::qm2_struct
  type (qmTheoryType), intent(in)::currentTheory
  type(ParameterEntry)::temp
  
! local variables

  integer :: atomic_number, n, i, j
  logical :: userDefinedVariable=.false.
  _REAL_ ::  tempIntegerInReal
  
  
! SECTION 1: Parameters common to all semi-empirical methods
!-----------------------------------------------------------
! core_chg = core charges of elements
! heat_of_form = heat of formation

! natomic_orbs = number of atomic orbitals
!***********************************************************************
!*                      VALENCE SHELLS ARE DEFINED AS                  *
!*  PQN   VALENCE SHELLS                                               *
!*                 P-GROUP              F-GROUP    TRANSITION METALS   *
!*   1       1S                                                        *
!*   2       2S 2P                                                     *
!*   3       3S 3P  OR  3S 3P 3D                                       *
!*   4       4S 4P                                    4S 4P 3D         *
!*   5       5S 5P                                    5S 5P 4D         *
!*   6       6S 6P                       6S 4F        6S 6P 5D         *
!*   7  NOT ASSIGNED YET  ****DO  NOT  USE****                         *
!***********************************************************************

!Enthalpies of formation of gaseous atoms are take from 'Annual reports,
!1974, 71B, P117'
!NOTE: for natomic_orbs only values of 1 and 4 are currently supported.

  core_chg(  1) = 1; natomic_orbs(  1) = 1; heat_of_form(  1) = 52.102D0 !H
  core_chg(  2) = 0; natomic_orbs(  2) = 1; heat_of_form(  2) =  0.000D0 !He

  core_chg(  3) = 1; natomic_orbs(  3) = 4; heat_of_form(  3) = 38.410D0 !Li
  core_chg(  4) = 2; natomic_orbs(  4) = 4; heat_of_form(  4) = 76.960D0 !Be
  core_chg(  5) = 3; natomic_orbs(  5) = 4; heat_of_form(  5) =135.700D0 !B
  core_chg(  6) = 4; natomic_orbs(  6) = 4; heat_of_form(  6) =170.890D0 !C
  core_chg(  7) = 5; natomic_orbs(  7) = 4; heat_of_form(  7) =113.000D0 !N
  core_chg(  8) = 6; natomic_orbs(  8) = 4; heat_of_form(  8) = 59.559D0 !O
  core_chg(  9) = 7; natomic_orbs(  9) = 4; heat_of_form(  9) = 18.890D0 !F
  core_chg( 10) = 0; natomic_orbs( 10) = 4; heat_of_form( 10) =  0.000D0 !Ne

  core_chg( 11) = 1; natomic_orbs( 11) = 4; heat_of_form( 11) = 25.850D0 !Na
  core_chg( 12) = 2; natomic_orbs( 12) = 4; heat_of_form( 12) = 35.000D0 !Mg
  core_chg( 13) = 3; natomic_orbs( 13) = 4; heat_of_form( 13) = 79.490D0 !Al
  core_chg( 14) = 4; natomic_orbs( 14) = 4; heat_of_form( 14) =108.390D0 !Si
  core_chg( 15) = 5; natomic_orbs( 15) = 4; heat_of_form( 15) = 75.570D0 !P
  core_chg( 16) = 6; natomic_orbs( 16) = 4; heat_of_form( 16) = 66.400D0 !S
  core_chg( 17) = 7; natomic_orbs( 17) = 4; heat_of_form( 17) = 28.990D0 !Cl
  core_chg( 18) = 0; natomic_orbs( 18) = 4; heat_of_form( 18) =  0.000D0 !Ar

  core_chg( 19) = 1; natomic_orbs( 19) = 4; heat_of_form( 19) = 21.420D0 !K
  core_chg( 20) = 2; natomic_orbs( 20) = 4; heat_of_form( 20) = 42.600D0 !Ca
  core_chg( 21) = 3; natomic_orbs( 21) = 9; heat_of_form( 21) = 90.300D0 !Sc
  core_chg( 22) = 4; natomic_orbs( 22) = 9; heat_of_form( 22) =112.300D0 !Ti
  core_chg( 23) = 5; natomic_orbs( 23) = 9; heat_of_form( 23) =122.900D0 !V
  core_chg( 24) = 6; natomic_orbs( 24) = 9; heat_of_form( 24) = 95.000D0 !Cr
  core_chg( 25) = 7; natomic_orbs( 25) = 9; heat_of_form( 25) = 67.700D0 !Mn
  core_chg( 26) = 8; natomic_orbs( 26) = 9; heat_of_form( 26) = 99.300D0 !Fe
  core_chg( 27) = 9; natomic_orbs( 27) = 9; heat_of_form( 27) =102.400D0 !Co
  core_chg( 28) =10; natomic_orbs( 28) = 9; heat_of_form( 28) =102.800D0 !Ni
  core_chg( 29) =11; natomic_orbs( 29) = 9; heat_of_form( 29) = 80.700D0 !Cu
  core_chg( 30) = 2; natomic_orbs( 30) = 4; heat_of_form( 30) = 31.170D0 !Zn
  core_chg( 31) = 3; natomic_orbs( 31) = 4; heat_of_form( 31) = 65.400D0 !Ga
  core_chg( 32) = 4; natomic_orbs( 32) = 4; heat_of_form( 32) = 89.500D0 !Ge
  core_chg( 33) = 5; natomic_orbs( 33) = 4; heat_of_form( 33) = 72.300D0 !As
  core_chg( 34) = 6; natomic_orbs( 34) = 4; heat_of_form( 34) = 54.300D0 !Se
  core_chg( 35) = 7; natomic_orbs( 35) = 4; heat_of_form( 35) = 26.740D0 !Br
  core_chg( 36) = 0; natomic_orbs( 36) = 4; heat_of_form( 36) =  0.000D0 !Kr

  core_chg( 37) = 1; natomic_orbs( 37) = 4; heat_of_form( 37) = 19.600D0 !Rb
  core_chg( 38) = 2; natomic_orbs( 38) = 4; heat_of_form( 38) = 39.100D0 !Sr
  core_chg( 39) = 3; natomic_orbs( 39) = 9; heat_of_form( 39) =101.500D0 !Y
  core_chg( 40) = 4; natomic_orbs( 40) = 9; heat_of_form( 40) =145.500D0 !Zr
  core_chg( 41) = 5; natomic_orbs( 41) = 9; heat_of_form( 41) =172.400D0 !Nb
  core_chg( 42) = 6; natomic_orbs( 42) = 9; heat_of_form( 42) =157.300D0 !Mo
  core_chg( 43) = 7; natomic_orbs( 43) = 9; heat_of_form( 43) =  0.000D0 !Tc
  core_chg( 44) = 8; natomic_orbs( 44) = 9; heat_of_form( 44) =155.500D0 !Ru
  core_chg( 45) = 9; natomic_orbs( 45) = 9; heat_of_form( 45) =133.000D0 !Rh
  core_chg( 46) =10; natomic_orbs( 46) = 9; heat_of_form( 46) = 90.000D0 !Pd
  core_chg( 47) =11; natomic_orbs( 47) = 9; heat_of_form( 47) = 68.100D0 !Ag
  core_chg( 48) = 2; natomic_orbs( 48) = 4; heat_of_form( 48) = 26.720D0 !Cd
  core_chg( 49) = 3; natomic_orbs( 49) = 4; heat_of_form( 49) = 58.000D0 !In
  core_chg( 50) = 4; natomic_orbs( 50) = 4; heat_of_form( 50) = 72.200D0 !Sn
  core_chg( 51) = 5; natomic_orbs( 51) = 4; heat_of_form( 51) = 63.200D0 !Sb
  core_chg( 52) = 6; natomic_orbs( 52) = 4; heat_of_form( 52) = 47.000D0 !Te
  core_chg( 53) = 7; natomic_orbs( 53) = 4; heat_of_form( 53) = 25.517D0 !I
  core_chg( 54) = 0; natomic_orbs( 54) = 4; heat_of_form( 54) =  0.000D0 !Xe

  core_chg( 55) = 1; natomic_orbs( 55) = 4; heat_of_form( 55) = 18.700D0 !Cs
  core_chg( 56) = 2; natomic_orbs( 56) = 4; heat_of_form( 56) = 42.500D0 !Ba
  core_chg( 57) = 3; natomic_orbs( 57) = 8; heat_of_form( 57) =  0.000D0 !La
  core_chg( 58) = 4; natomic_orbs( 58) = 8; heat_of_form( 58) =101.300D0 !Ce
  core_chg( 59) = 5; natomic_orbs( 59) = 8; heat_of_form( 59) =  0.000D0 !Pr
  core_chg( 60) = 6; natomic_orbs( 60) = 8; heat_of_form( 60) =  0.000D0 !Nd
  core_chg( 61) = 7; natomic_orbs( 61) = 8; heat_of_form( 61) =  0.000D0 !Pm
  core_chg( 62) = 8; natomic_orbs( 62) = 8; heat_of_form( 62) = 49.400D0 !Sm
  core_chg( 63) = 9; natomic_orbs( 63) = 8; heat_of_form( 63) =  0.000D0 !Eu
  core_chg( 64) =10; natomic_orbs( 64) = 8; heat_of_form( 64) =  0.000D0 !Gd
  core_chg( 65) =11; natomic_orbs( 65) = 8; heat_of_form( 65) =  0.000D0 !Tb
  core_chg( 66) =12; natomic_orbs( 66) = 8; heat_of_form( 66) =  0.000D0 !Dy
  core_chg( 67) =13; natomic_orbs( 67) = 8; heat_of_form( 67) =  0.000D0 !Ho
  core_chg( 68) =14; natomic_orbs( 68) = 8; heat_of_form( 68) = 75.800D0 !Er
  core_chg( 69) =15; natomic_orbs( 69) = 8; heat_of_form( 69) =  0.000D0 !Tm
  core_chg( 70) =16; natomic_orbs( 70) = 8; heat_of_form( 70) = 36.350D0 !Yb
  core_chg( 71) = 3; natomic_orbs( 71) = 9; heat_of_form( 71) =  0.000D0 !Lu
  core_chg( 72) = 4; natomic_orbs( 72) = 9; heat_of_form( 72) =148.000D0 !Hf
  core_chg( 73) = 5; natomic_orbs( 73) = 9; heat_of_form( 73) =186.900D0 !Ta
  core_chg( 74) = 6; natomic_orbs( 74) = 9; heat_of_form( 74) =203.100D0 !W
  core_chg( 75) = 7; natomic_orbs( 75) = 9; heat_of_form( 75) =185.000D0 !Re
  core_chg( 76) = 8; natomic_orbs( 76) = 9; heat_of_form( 76) =188.000D0 !Os
  core_chg( 77) = 9; natomic_orbs( 77) = 9; heat_of_form( 77) =160.000D0 !Ir
  core_chg( 78) =10; natomic_orbs( 78) = 9; heat_of_form( 78) =135.200D0 !Pt
  core_chg( 79) =11; natomic_orbs( 79) = 9; heat_of_form( 79) = 88.000D0 !Au
  core_chg( 80) = 2; natomic_orbs( 80) = 4; heat_of_form( 80) = 14.690D0 !Hg
  core_chg( 81) = 3; natomic_orbs( 81) = 4; heat_of_form( 81) = 43.550D0 !Tl
  core_chg( 82) = 4; natomic_orbs( 82) = 4; heat_of_form( 82) = 46.620D0 !Pb
  core_chg( 83) = 5; natomic_orbs( 83) = 4; heat_of_form( 83) = 50.100D0 !Bi
  core_chg( 84) = 6; natomic_orbs( 84) = 4; heat_of_form( 84) =  0.000D0 !Po
  core_chg( 85) = 7; natomic_orbs( 85) = 4; heat_of_form( 85) =  0.000D0 !At
  core_chg( 86) = 0; natomic_orbs( 86) = 4; heat_of_form( 86) =  0.000D0 !Rn

! d-orbital support for 3rd row elements (except Na) (MNDOD)
if (currentTheory%MNDOD) then 
  core_chg( 11) = 1; natomic_orbs( 11) = 4; heat_of_form( 11) = 25.850D0 !Na
  core_chg( 12) = 2; natomic_orbs( 12) = 4; heat_of_form( 12) = 35.000D0 !Mg
  core_chg( 13) = 3; natomic_orbs( 13) = 9; heat_of_form( 13) = 79.490D0 !Al
  core_chg( 14) = 4; natomic_orbs( 14) = 9; heat_of_form( 14) =108.390D0 !Si
  core_chg( 15) = 5; natomic_orbs( 15) = 9; heat_of_form( 15) = 75.570D0 !P 
  core_chg( 16) = 6; natomic_orbs( 16) = 9; heat_of_form( 16) = 66.400D0 !S
  core_chg( 17) = 7; natomic_orbs( 17) = 9; heat_of_form( 17) = 28.990D0 !Cl  
end if

if (currentTheory%AM1D) then 
  core_chg( 11) = 1; natomic_orbs( 11) = 4; heat_of_form( 11) = 25.850D0 !Na
  core_chg( 12) = 2; natomic_orbs( 12) = 9; heat_of_form( 12) = 35.000D0 !Mg
  core_chg( 13) = 3; natomic_orbs( 13) = 4; heat_of_form( 13) = 79.490D0 !Al
  core_chg( 14) = 4; natomic_orbs( 14) = 4; heat_of_form( 14) =108.390D0 !Si
  core_chg( 15) = 5; natomic_orbs( 15) = 9; heat_of_form( 15) = 75.570D0 !P 
  core_chg( 16) = 6; natomic_orbs( 16) = 4; heat_of_form( 16) = 66.400D0 !S
  core_chg( 17) = 7; natomic_orbs( 17) = 4; heat_of_form( 17) = 28.990D0 !Cl  
end if

! d-orbital support for 3rd row elements (PM6)
if (currentTheory%PM6) then 
  natomic_orbs( 13) = 9 !Al
  natomic_orbs( 14) = 9 !Si
  natomic_orbs( 15) = 9 !P
  natomic_orbs( 16) = 9 !S
  natomic_orbs( 17) = 9 !Cl
  natomic_orbs( 33) = 9 !As
  natomic_orbs( 35) = 9 !Br
  natomic_orbs( 51) = 9 !Sb
  natomic_orbs( 53) = 9 !I
end if

 !Originally eisol the electronic energy used to be stored as a parameter. As did
 !DD,QQ,AM,AD,AQ - Now we calculate these as they are actually derivatives of other
 !parameters.

 !For the electronic energy the calculation is:
 !elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
 !           + HSP*HSPC
 ! where gssc is the number of two-electron terms of type <SS|SS>
 !             = max(ios-1,0)
 !       gspc is the number of two-electron terms of type <SS|PP>
 !             = ios * iop
 !       gp2c is the number of two-electron terms of type <PP|PP> + o.5 of the number of HPP
 !            integrals which is not used but instead is replaced by 0.5(GPP-GP2)
 !             = (iop * (iop-1))/2 + 0.5*(min(iop,6-iop)*((min(iop,6-iop)-1))/2
 !       gppc is minus 0.5 x the number of HPP integrals.
 !             = -0.5*(min(iop,6-iop)*((min(iop,6-iop)-1))/2
 !       HSPC is the number of two electron terms of type <SP|SP>.
 !            S and P must have the same spin. If P is non-zero there are two S electrons.
 !             = -iop
 !
 !          H                                                               He
 !          Li Be                                            B  C  N  O  F  Ne
 !          Na Mg                                            Al Si P  S  Cl Ar
 !          K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
 !          Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
 !          Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
 !          Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
 ! IOS
 !       &/ 1,                                                                2, &!    2
 !       &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   10
 !       &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   18
 !       &  1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0, &!   36
 !       &  1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0, &!   54
 !       &  1, 2, 2, 5*0,3*2,6*2, 2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0   !   86

 ! IOP
 !            / 0 ,                                                           0, &!    2
 !           &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   10
 !           &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   18
 !           &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   36
 !           &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   54
 !           &  0, 0, 0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6   !   86

 ! so for PM3 carbon this would be:
 !
 ! USS = -47.2703200D0
 ! IOS = 2
 ! UPP = -36.2669180D0
 ! IOP = 2
 ! GSS = 11.2007080D0
 ! GSSC = 1
 ! GPP = 10.7962920D0
 ! GPPC = -0.5   (-0.5*(2*1)/2)
 ! GSP = 10.2650270D0
 ! GSPC = 4
 ! GP2 = 9.0425660D0
 ! GP2C = 1 + 0.5 = 1.5
 ! HSP = 2.2909800D0
 ! HSPC = -2
 !
 ! Therefore elec_eng_pm3(carbon) = -11.229917

 !FILL IOS ARRAY - Initial S orbital occupancies
   ios = (/ &
            &  1,                                                                2,   &!    2
            &  1, 2,                                              2, 2, 2, 2, 2, 0,   &!   10
            &  1, 2,                                              2, 2, 2, 2, 2, 0,   &!   18
            &  1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0,   &!   36
            &  1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0,   &!   54
            &  1, 2, 2, &
                        0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                                     2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0 /)

 !FILL IOP ARRAY - Initial P orbital occupancies
   iop = (/ &
            &  0 ,                                                           0,   &!    2
            &  0, 0,                                          1, 2, 3, 4, 5, 6,   &!   10
            &  0, 0,                                          1, 2, 3, 4, 5, 6,   &!   18
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6,   &!   36
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6,   &!   54
            &  0, 0, 0, &
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                                 0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6 /)

 ! For calculation of DD, QQ, AM, AD and AQ Arrays:
 !
 ! DD = ( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)
 !      ----------------------------------------------------------------------
 !        ( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0)
 !
 ! QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
 !
 ! AM = gss/AU_TO_EV
 !
 ! So for PM3 carbon this would be:
 !
 ! nsshell = 2
 ! GSS = 11.2007080D0
 ! s_orb_exp = 1.5650850D0
 ! p_orb_exp = 1.8423450D0
 ! DD = ( (4.0d0*1.5650850D0*1.8423450D0)**(2+0.5d0) ) * (2.0d0*2 + 1)
 !      --------------------------------------------------------------------
 !        ( (1.5650850D0+1.8423450D0)**(2.0d0*2 + 2.0d0) ) * (sqrt(3.0d0)
 !    = 0.8332396384d0
 ! QQ = sqrt((4.0d0*4+6.0d0*2+2.0d0)/20.0)/1.8423450D0 = 0.664775d0
 ! AD = gdd1 + df*((hsp(i)/ev-hsp1)/(hsp2 - hsp1)

 !FILL NSSHELL ARRAY - I believe that this is the S shell number - or essentially the row
 !                     number in the periodic table.
   nsshell = (/ &
            &  1 ,                                                           1,   &!    2
            &  2, 2,                                          2, 2, 2, 2, 2, 2,   &!   10
            &  3, 3,                                          3, 3, 3, 3, 3, 3,   &!   18
            &  4, 4, 4,          4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4,   &!   36
            &  5, 5, 5,          5, 5, 5, 5, 5, 5, 5, 5, 5,   5, 5, 5, 5, 5, 5,   &!   54
            &  6, 6, 6, &
                        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,   &
                                 6, 6, 6, 6, 6, 6, 6, 6, 6,   6, 6, 6, 6, 6, 6 /)


!-----------------------------------------------------------
!END OF SECTION 1 - COMMON PARAMETERS
!-----------------------------------------------------------


!------------------------------------------------------------------------------
!SECTION 2 - MNDO, AM1, PM3, RM1, PM6, PDDG/PM3 and PDDG/MNDO PARAMS BY ELEMENT
!------------------------------------------------------------------------------
!This section contains the MNDO, AM1, RM1, PM3, PM6 and PDDG parameter sets
!for each element that is supported as well as additional scaling terms for PM3/MM*
!Set element supported array to false and then set each entry true
!as and when we add the parameters.

  element_supported_mndo  = .false.
  element_supported_mndod  = .false.
  element_supported_am1  = .false.
  element_supported_am1d  = .false.  
  element_supported_rm1  = .false.
  element_supported_pm3  = .false.
  element_supported_pm3mmx  = .false.  ! PM3/MM*
  element_supported_pm3mmx2  = .false.  ! PM3/MM* 2nd version
  element_supported_pm3mais  = .false.
  element_supported_pm6  = .false.
  element_supported_pddgpm3  = .false.
  element_supported_pddgpm3_08  = .false.
  element_supported_pddgmndo  = .false.

!Parameter meanings
!   s_orb_exp, p_orb_exp - The Slater exponents of the basis functions
!   betas, betap - two centre, one electron core integral parameters.
!   FN1,FN2,FN3 - PM3 / AM1 / RM1 specific parameters for the core-core interactions
!    GSS ::= (SS,SS)
!    GPP ::= (PP,PP)
!    GSP ::= (SS,PP)
!    GP2 ::= (PP,P*P*)
!    HSP ::= (SP,SP)
!   GSS, GSP, GPP, GP2, HSP - Coulomb and exchange one centre-two electron integral params.
!   DD, QQ, AD, AM, AQ - parameters for multipole expansion of the two centre, two electron integrals.
!   ALP - Exponents for the core-core repulsion terms
!   USS, UPP - electron kinetic energy integral parameters.

! Initialise the arrays
  elec_eng_pddgmndo  = 0.0d0 !Only pddg routines treat this as an independent parameter
  elec_eng_pddgpm3  = 0.0d0  !in all others is calculated from the other parameters.
  elec_eng_pddgpm3_08  = 0.0d0
  s_orb_exp_mndo  = 0.0d0
  s_orb_exp_mndod  = 0.0d0
  s_orb_exp_am1  = 0.0d0
  s_orb_exp_am1d  = 0.0d0  
  s_orb_exp_rm1  = 0.0d0
  s_orb_exp_pm3  = 0.0d0
  s_orb_exp_pm6  = 0.0d0
  s_orb_exp_pddgmndo  = 0.0d0
  s_orb_exp_pddgpm3  = 0.0d0
  s_orb_exp_pddgpm3_08  = 0.0d0
  p_orb_exp_mndo  = 0.0d0
  p_orb_exp_mndod  = 0.0d0
  p_orb_exp_am1  = 0.0d0
  p_orb_exp_am1d  = 0.0d0  
  p_orb_exp_rm1  = 0.0d0
  p_orb_exp_pm3  = 0.0d0
  p_orb_exp_pm6  = 0.0d0
  p_orb_exp_pddgmndo  = 0.0d0
  p_orb_exp_pddgpm3  = 0.0d0
  p_orb_exp_pddgpm3_08  = 0.0d0
  d_orb_exp_mndod  = 0.0d0  
  d_orb_exp_am1d  = 0.0d0    
  d_orb_exp_pm6   = 0.0d0    
  s_orb_exp_tail_am1d  = 0.0d0 
  p_orb_exp_tail_am1d  = 0.0d0 
  d_orb_exp_tail_am1d  = 0.0d0   
  s_orb_exp_tail_mndod  = 0.0d0 
  p_orb_exp_tail_mndod  = 0.0d0 
  d_orb_exp_tail_mndod  = 0.0d0     
  s_orb_exp_tail_pm6 = 0.0d0 
  p_orb_exp_tail_pm6 = 0.0d0 
  d_orb_exp_tail_pm6 = 0.0d0     
  betas_mndo  = 0.0d0
  betas_mndod  = 0.0d0
  betas_am1  = 0.0d0
  betas_am1d  = 0.0d0  
  betas_rm1  = 0.0d0
  betas_pm3  = 0.0d0
  betas_pm6  = 0.0d0
  betas_pddgmndo  = 0.0d0
  betas_pddgpm3  = 0.0d0
  betas_pddgpm3_08  = 0.0d0
  betas_pm3carb1  = 0.0d0
  betap_mndo  = 0.0d0
  betap_mndod  = 0.0d0
  betap_am1  = 0.0d0
  betap_am1d  = 0.0d0  
  betap_rm1  = 0.0d0
  betap_pm3  = 0.0d0
  betap_pm6  = 0.0d0
  betap_pddgmndo  = 0.0d0
  betap_pddgpm3  = 0.0d0
  betap_pddgpm3_08  = 0.0d0
  betap_pm3carb1  = 0.0d0
  betad_mndod  = 0.0d0  
  betad_am1d  = 0.0d0    
  betad_pm6 = 0.0d0  
  GNN_am1d=1.0d0
  rho_core_am1d=0.0d0
  rho_core_mndod=0.0d0  
  rho_core_pm6=0.0d0  
  FN1_am1 = 0.0d0
  FN2_am1 = 0.0d0
  FN3_am1 = 0.0d0
  NUM_FN_am1 = 0
  FN1_am1d = 0.0d0
  FN2_am1d = 0.0d0
  FN3_am1d = 0.0d0
  NUM_FN_am1d = 0  
  FN1_rm1 = 0.0d0
  FN2_rm1 = 0.0d0
  FN3_rm1 = 0.0d0
  NUM_FN_rm1 = 0
  FN1_pm6 = 0.0d0
  FN2_pm6 = 0.0d0
  FN3_pm6 = 0.0d0
  NUM_FN_pm6 = 0
  FN1_pm3 = 0.0d0
  FN2_pm3 = 0.0d0
  FN3_pm3 = 0.0d0
  NUM_FN_pm3 = 0
  FN1_pddgpm3 = 0.0d0
  FN2_pddgpm3 = 0.0d0
  FN3_pddgpm3 = 0.0d0
  NUM_FN_pddgpm3 = 0
  FN1_pddgpm3_08 = 0.0d0
  FN2_pddgpm3_08 = 0.0d0
  FN3_pddgpm3_08 = 0.0d0
  NUM_FN_pddgpm3_08 = 0
  GSS_mndo = 0.0d0; GSP_mndo = 0.0d0; GPP_mndo = 0.0d0; GP2_mndo = 0.0d0; HSP_mndo = 0.0d0
  GSS_mndod = 0.0d0; GSP_mndod = 0.0d0; GPP_mndod = 0.0d0; GDD_mndod = 0.0d0; GP2_mndod = 0.0d0; HSP_mndod = 0.0d0
  GSS_am1 = 0.0d0; GSP_am1 = 0.0d0; GPP_am1 = 0.0d0; GP2_am1 = 0.0d0; HSP_am1 = 0.0d0
  GSS_am1d = 0.0d0; GSP_am1d = 0.0d0; GPP_am1d = 0.0d0; GDD_am1d = 0.0d0; GP2_am1d = 0.0d0; HSP_am1d = 0.0d0  
  GSS_rm1 = 0.0d0; GSP_rm1 = 0.0d0; GPP_rm1 = 0.0d0; GP2_rm1 = 0.0d0; HSP_rm1 = 0.0d0
  GSS_pm3 = 0.0d0; GSP_pm3 = 0.0d0; GPP_pm3 = 0.0d0; GP2_pm3 = 0.0d0; HSP_pm3 = 0.0d0
  GSS_pm6 = 0.0d0; GSP_pm6 = 0.0d0; GPP_pm6 = 0.0d0; GP2_pm6 = 0.0d0; HSP_pm6 = 0.0d0
  GSS_pddgmndo = 0.0d0; GSP_pddgmndo = 0.0d0; GPP_pddgmndo = 0.0d0; GP2_pddgmndo = 0.0d0; HSP_pddgmndo = 0.0d0
  GSS_pddgpm3 = 0.0d0; GSP_pddgpm3 = 0.0d0; GPP_pddgpm3 = 0.0d0; GP2_pddgpm3 = 0.0d0; HSP_pddgpm3 = 0.0d0
  GSS_pddgpm3_08 = 0.0d0; GSP_pddgpm3_08 = 0.0d0; GPP_pddgpm3_08 = 0.0d0; GP2_pddgpm3_08 = 0.0d0; HSP_pddgpm3_08 = 0.0d0
  alp_mndo = 0.0d0; alp_mndod = 0.0d0; alp_am1=0.0d0; alp_am1d=0.0d0; alp_rm1=0.0d0; alp_pm3 = 0.0d0; alp_pddgmndo = 0.0d0
  alp_pddgpm3 = 0.0d0; alp_pddgpm3_08 = 0.0d0; alp_pm3carb1 = 0.0d0; alp_pm6 = 0.0d0
  uss_mndo = 0.0d0; uss_mndod = 0.0d0; uss_am1=0.0d0; uss_am1d=0.0d0     
  uss_rm1=0.0d0; uss_pm3=0.0d0; uss_pm6=0.0d0; uss_pddgmndo = 0.0d0
  uss_pddgpm3 = 0.0d0; uss_pddgpm3_08 = 0.0d0; uss_pm3carb1 = 0.0d0
  upp_mndo = 0.0d0; upp_mndod = 0.0d0; upp_am1=0.0d0; upp_am1d=0.0d0   
  upp_rm1=0.0d0; upp_pm3=0.0d0; upp_pm6=0.0d0; upp_pddgmndo = 0.0d0
  upp_pddgpm3 = 0.0d0; upp_pddgpm3_08 = 0.0d0; upp_pm3carb1 = 0.0d0
  udd_mndod = 0.0d0; udd_am1d = 0.0d0; udd_pm6 = 0.0d0
  PDDGC1_pm3 = 0.0d0; PDDGC2_pm3 = 0.0d0; PDDGE1_pm3 = 0.0d0; PDDGE2_pm3 = 0.0d0
  PDDGC1_pm3_08 = 0.0d0; PDDGC2_pm3_08 = 0.0d0; PDDGE1_pm3_08 = 0.0d0; PDDGE2_pm3_08 = 0.0d0
  PDDGC1_mndo = 0.0d0; PDDGC2_mndo = 0.0d0; PDDGE1_mndo = 0.0d0; PDDGE2_mndo = 0.0d0
  
  F0SD_pm6 = 0.0d0
  G2SD_pm6 = 0.0d0

  ! To be removed
  EISOL_pm6 = 0.0d0

  mndo_ref_index = 0; mndod_ref_index = 0; am1_ref_index = 0; am1d_ref_index = 0 
  rm1_ref_index = 0; pm3_ref_index = 0; pm6_ref_index = 0
  pm3carb1_ref_index = 0
  pddgpm3_ref_index = 0; pddgpm3_08_ref_index = 0; pddgmndo_ref_index = 0

  alpab_pm6 = 0.0d0
  xab_pm6 = 0.0d0

! PM3-MAIS parameters.
  alpab_pm3mais = 0.0d0
  betab_pm3mais = 0.0d0
  gamab_pm3mais = 0.0d0

  scale_f1_pm3mmx= 0.0d0; scale_f2_pm3mmx = 0.0d0  ! PM3/MM* 
  scale_f1_pm3mmx2= 0.0d0; scale_f2_pm3mmx2 = 0.0d0  ! PM3/MM* 2nd version
  rho_pm3mmx2 = 0.0d0
  
! OPNQ 
    qxd_supported=.false.
    qxd_s=0.0d0;  qxd_z0=0.0d0;  qxd_zq=0.0d0;  qxd_d0=0.0d0;  qxd_dq=0.0d0
    qxd_q0=0.0d0;  qxd_qq=0.0d0;  qxd_neff=0.0d0

!------- EXPLANATION OF PARAMETER FIELDS ------------
!xxx_ref_index = a character array index for where to find the text to
!                print as the citation for this parameter set.
!element_supported_xxx = flag for whether parameters exist for this element.
!s_orb_exp_xxx = s-type Slater atomic orbital exponent (zeta-s)
!p_orb_exp_xxx = p-type Slater atomic orbital exponent (zeta-p)
!betas_xxx = s atomic orbital one-electron two-center resonance integral term.
!betap_xxx = p atomic orbital one-electron two-center resonance integral term.
!alp_xx = atom core-core repulsion term.
!FN1(1:4) = Gaussian multiplier for the ith Gaussian of atom. (ai)
!FN2(1:4) = Gaussian exponent multiplier for the ith Gaussian of atom. (bi)
!FN3(1:4) = radial center of the ith Gaussian of atom. (ci)
!NUM_FN_xxx = Number of i's.
!GSS_xxx = s-s atomic orbitals one-center two-electron repulsion integral.
!GSP_xxx = s-p atomic orbitals one-center two-electron replusion integral.
!GPP_xxx = p-p atomic orbitals one-center two-electron replusion integral.
!GP2_xxx = p-p' atomic orbitals one-center two-electron replusion integral.
!HSP_xxx = s-p atomic orbital one-center two-electron exchange integral.
!USS_xxx = s atomic orbital one-electron one-center integral.
!UPP_xxx = p atomic orbital one-electron one-center integral.
!----------------------------------------------------


 !-------------------
 !HYDROGEN
 !-------------------
 !Notes for elec eng: IOS = 1, IOP = 0, GSSC = 0, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0

   atomic_number = 1
  !MNDO
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (Index = 1)
     mndo_ref_index(atomic_number) = 1
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.3319670D0
     p_orb_exp_mndo(atomic_number) = 0.0d0
     betas_mndo(atomic_number) = -6.9890640D0
     betap_mndo(atomic_number) = 0.0d0
     GSS_mndo(atomic_number) = 12.848D00
     alp_mndo(atomic_number) = 2.5441341D0
     USS_mndo(atomic_number) = -11.9062760D0
     UPP_mndo(atomic_number) = 0.0d0

  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (Index = 1)
     mndod_ref_index(atomic_number) = 1
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.3319670D0
     p_orb_exp_mndod(atomic_number) = 0.0d0
     betas_mndod(atomic_number) = -6.9890640D0
     betap_mndod(atomic_number) = 0.0d0
     GSS_mndod(atomic_number) = 12.848D00
     alp_mndod(atomic_number) = 2.5441341D0
     USS_mndod(atomic_number) = -11.9062760D0
     UPP_mndod(atomic_number) = 0.0d0

  !AM1
   ! Reference: M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (Index = 17)
     am1_ref_index(atomic_number) = 17
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.1880780D0
     p_orb_exp_am1(atomic_number) = 0.0d0
     betas_am1(atomic_number) = -6.1737870D0
     betap_am1(atomic_number) = 0.0d0
     FN1_am1(1,atomic_number) = 0.1227960D0
     FN2_am1(1,atomic_number) = 5.0000000D0
     FN3_am1(1,atomic_number) = 1.2000000D0
     FN1_am1(2,atomic_number) = 0.0050900D0
     FN2_am1(2,atomic_number) = 5.0000000D0
     FN3_am1(2,atomic_number) = 1.8000000D0
     FN1_am1(3,atomic_number) =-0.0183360D0
     FN2_am1(3,atomic_number) = 2.0000000D0
     FN3_am1(3,atomic_number) = 2.1000000D0
     FN1_am1(4,atomic_number) = 0.0d0
     FN2_am1(4,atomic_number) = 0.0d0
     FN3_am1(4,atomic_number) = 0.0d0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 12.8480000D0
     alp_am1(atomic_number) = 2.8823240D0
     USS_am1(atomic_number) = -11.3964270D0
     UPP_am1(atomic_number) = 0.0d0

  !AM1D
   ! Reference: K. Nam, Q. Cui, J. Gao, D. York. J. CHEM. THEO. COMP., 3, 486, (2007) (102) AM1/d-PhoT
     am1d_ref_index(atomic_number) = 102
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.143846D0
     p_orb_exp_am1d(atomic_number) = 0.0d0
     betas_am1d(atomic_number) = -5.9111080D0
     betap_am1d(atomic_number) = 0.0d0
     FN1_am1d(1,atomic_number) = 0.106238D0
     FN2_am1d(1,atomic_number) = 5.735290D0
     FN3_am1d(1,atomic_number) = 1.261430D0
     FN1_am1d(2,atomic_number) = 0.004043D0
     FN2_am1d(2,atomic_number) = 7.080122D0
     FN3_am1d(2,atomic_number) = 2.084095D0
     FN1_am1d(3,atomic_number) =-0.002800D0
     FN2_am1d(3,atomic_number) = 0.739913D0
     FN3_am1d(3,atomic_number) = 3.649474D0
     FN1_am1d(4,atomic_number) = 0.0d0
     FN2_am1d(4,atomic_number) = 0.0d0
     FN3_am1d(4,atomic_number) = 0.0d0
     NUM_FN_am1d(atomic_number) = 3
     GSS_am1d(atomic_number) = 13.737453D0
     alp_am1d(atomic_number) = 2.884915D0
     USS_am1d(atomic_number) = -10.934610D0
     UPP_am1d(atomic_number) = 0.0d0


  !PM3
   ! Reference J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 0.9678070D0
     p_orb_exp_pm3(atomic_number) = 0.0d0
     betas_pm3(atomic_number) = -5.6265120D0
     betap_pm3(atomic_number) = 0.0d0
     FN1_pm3(1,atomic_number) = 1.1287500D0
     FN2_pm3(1,atomic_number) = 5.0962820D0
     FN3_pm3(1,atomic_number) = 1.5374650D0
     FN1_pm3(2,atomic_number) =-1.0603290D0
     FN2_pm3(2,atomic_number) = 6.0037880D0
     FN3_pm3(2,atomic_number) = 1.5701890D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 14.7942080D0
     alp_pm3(atomic_number) = 3.3563860D0
     USS_pm3(atomic_number) = -13.0733210D0
     UPP_pm3(atomic_number) = 0.0d0

  !PM3/MM* - A reformulated QM core-MM charge interface. Has improved performance particularly
  !          for hydrogen-bonded molecules. Current parameters are only available for H, C, N
  !          and O QM atoms. For H, seperate parameter sets are assigned to H-H and H-heavy pairs.
   ! Reference: WANG,Q.T. and BRYCE,R.A. (2009) JCTC, 5, p2206
     element_supported_pm3mmx(atomic_number) = .true.
     ! For H(QM)-H(MM)
     scale_f1_pm3mmx(1,atomic_number) = 3.4D0
     scale_f2_pm3mmx(1,atomic_number) = 3.6D0
     ! For H(QM)-Other MM
     scale_f1_pm3mmx(2,atomic_number) = 2.2D0
     scale_f2_pm3mmx(2,atomic_number) = 2.7D0

  !PM3/MM* 2nd version
     element_supported_pm3mmx2(atomic_number) = .true.
     rho_pm3mmx2(atomic_number) = -0.293D0
     ! For H(QM)-H(MM)
     scale_f1_pm3mmx2(1,atomic_number) = 3.147D0
     scale_f2_pm3mmx2(1,atomic_number) = 3.210D0
     ! For H(QM)-Other MM
     scale_f1_pm3mmx2(2,atomic_number) = 3.040D0
     scale_f2_pm3mmx2(2,atomic_number) = 3.748D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.268641d0
     p_orb_exp_pm6(atomic_number) = 0.0d0
     betas_pm6(atomic_number) = -8.352984d0
     betap_pm6(atomic_number) = 0.0d0
     FN1_pm6(1,atomic_number) = 0.024184d0
     FN2_pm6(1,atomic_number) = 3.055953d0
     FN3_pm6(1,atomic_number) = 1.786011d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 14.448686d0
     USS_pm6(atomic_number) = -11.246958d0
     UPP_pm6(atomic_number) = 0.0d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 3.3563860D0
     !For pairwise core core terms see section 3 below.

  !PDDG/PM3
   ! Reference Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgpm3_ref_index(atomic_number) = 29
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -13.120566198192D0
     s_orb_exp_pddgpm3(atomic_number) = 0.97278550084430D0
     p_orb_exp_pddgpm3(atomic_number) = 0.0d0
     betas_pddgpm3(atomic_number) = -6.1526542062173D0
     betap_pddgpm3(atomic_number) = 0.0d0
     FN1_pddgpm3(1,atomic_number) = 1.12224395962630D0
     FN2_pddgpm3(1,atomic_number) = 4.70779030777590D0
     FN3_pddgpm3(1,atomic_number) = 1.54709920873910D0
     FN1_pddgpm3(2,atomic_number) =-1.0697373657305D0
     FN2_pddgpm3(2,atomic_number) = 5.85799464741120D0
     FN3_pddgpm3(2,atomic_number) = 1.56789274832050D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 14.7942080D0
     alp_pddgpm3(atomic_number) = 3.38168610300700D0
     USS_pddgpm3(atomic_number) = -12.893272003385D0
     UPP_pddgpm3(atomic_number) = 0.0d0
     PDDGC1_pm3(atomic_number) = 0.05719290135800D0
     PDDGC2_pm3(atomic_number) = -0.0348228612590D0
     PDDGE1_pm3(atomic_number) = 0.66339504047230D0
     PDDGE2_pm3(atomic_number) = 1.08190071942210D0

  !PDDG/PM3_08
   ! Reference J. Tirado-Rives et al. J. CHEM. THEO. COMP., 4, 297, (2008) (Index 34)
     pddgpm3_08_ref_index(atomic_number) = 34
     element_supported_pddgpm3_08(atomic_number) = .true.
     elec_eng_pddgpm3_08(atomic_number) = -13.248091d0
     s_orb_exp_pddgpm3_08(atomic_number) = 0.988391d0
     p_orb_exp_pddgpm3_08(atomic_number) = 0.0d0
     betas_pddgpm3_08(atomic_number) = -6.162383d0
     betap_pddgpm3_08(atomic_number) = 0.0d0
     FN1_pddgpm3_08(1,atomic_number) = 1.127822d0
     FN2_pddgpm3_08(1,atomic_number) = 4.750023d0
     FN3_pddgpm3_08(1,atomic_number) = 1.549373d0
     FN1_pddgpm3_08(2,atomic_number) =-1.074605d0
     FN2_pddgpm3_08(2,atomic_number) = 5.870974d0
     FN3_pddgpm3_08(2,atomic_number) = 1.566692d0
     FN1_pddgpm3_08(3,atomic_number) = 0.0d0
     FN2_pddgpm3_08(3,atomic_number) = 0.0d0
     FN3_pddgpm3_08(3,atomic_number) = 0.0d0
     FN1_pddgpm3_08(4,atomic_number) = 0.0d0
     FN2_pddgpm3_08(4,atomic_number) = 0.0d0
     FN3_pddgpm3_08(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3_08(atomic_number) = 2
     GSS_pddgpm3_08(atomic_number) = 14.7942080D0
     alp_pddgpm3_08(atomic_number) = 3.34016d0
     USS_pddgpm3_08(atomic_number) = -13.043714D0
     UPP_pddgpm3_08(atomic_number) = 0.0d0
     PDDGC1_pm3_08(atomic_number) = 0.057812d0
     PDDGC2_pm3_08(atomic_number) =-0.035533d0
     PDDGE1_pm3_08(atomic_number) = 0.683017d0
     PDDGE2_pm3_08(atomic_number) = 1.113826d0

  !PDDG/MNDO
   ! Reference: Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (index 29)
     pddgmndo_ref_index(atomic_number) = 29
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -12.015955786557D0
     s_orb_exp_pddgmndo(atomic_number) = 1.32243115467370D0
     p_orb_exp_pddgmndo(atomic_number) = 0.0d0
     betas_pddgmndo(atomic_number) = -7.4935039195719D0
     betap_pddgmndo(atomic_number) = 0.0d0
     GSS_pddgmndo(atomic_number) = 12.848D00
     alp_pddgmndo(atomic_number) = 2.49181323064320D0
     USS_pddgmndo(atomic_number) = -11.724114276410D0
     UPP_pddgmndo(atomic_number) = 0.0d0
     PDDGC1_mndo(atomic_number) = -0.1088607444359D0
     PDDGC2_mndo(atomic_number) = -0.0247060666203D0
     PDDGE1_mndo(atomic_number) = 0.46072116172000D0
     PDDGE2_mndo(atomic_number) = 1.29873123436820D0

  !PM3CARB-1
  !Reference: McNamara, J.P., Muslim, A.M., Abdel-Aal, H., Wang, H., Mohr, M., Hillier, I.H., Bryce, R.A.,
  !           2004, Chem Phys Lett, 394, 429-436 (Index 28)
  !Note: PM3CARB-1 is not a complete parameter set, it is simply PM3 with a few parameters changed for
  !      Oxygen and Hydrogen. So when loading these parameters into the correct structure the code
  !      needs to load pm3 parameters except for those that are different for PM3CARB-1
    pm3carb1_ref_index(atomic_number) = 28
    USS_PM3CARB1(atomic_number) = -13.514849D0
    UPP_PM3CARB1(atomic_number) = 0.0d0
    betas_PM3CARB1(atomic_number) = -4.011786D0
    betap_PM3CARB1(atomic_number) = 0.0d0
    alp_PM3CARB1(atomic_number) = 2.753199D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 1.08267366D0
     p_orb_exp_rm1(atomic_number) = 0.0d0
     betas_rm1(atomic_number) = -5.76544469D0
     betap_rm1(atomic_number) = 0.0d0
     FN1_rm1(1,atomic_number) = 0.10288875D0
     FN2_rm1(1,atomic_number) = 5.90172268D0
     FN3_rm1(1,atomic_number) = 1.17501185D0
     FN1_rm1(2,atomic_number) = 0.06457449D0
     FN2_rm1(2,atomic_number) = 6.41785671D0
     FN3_rm1(2,atomic_number) = 1.93844484D0
     FN1_rm1(3,atomic_number) =-0.03567387D0
     FN2_rm1(3,atomic_number) = 2.80473127D0
     FN3_rm1(3,atomic_number) = 1.63655241D0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 3
     GSS_rm1(atomic_number) = 13.98321296D0
     alp_rm1(atomic_number) = 3.06835947D0
     USS_rm1(atomic_number) = -11.96067697D0
     UPP_rm1(atomic_number) = 0.0d0

 !-------------------
 !END HYDROGEN
 !-------------------

 !-------------------
 !HELIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0
     atomic_number = 2
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 3.313204d0
     p_orb_exp_pm6(atomic_number) = 3.657133d0
     betas_pm6(atomic_number) = -58.903774d0
     betap_pm6(atomic_number) = -37.039974d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  9.445299d0
     GSP_pm6(atomic_number) = 11.201419d0
     GPP_pm6(atomic_number) =  9.214548d0
     GP2_pm6(atomic_number) = 13.046115d0
     HSP_pm6(atomic_number) =  0.299954d0
     USS_pm6(atomic_number) = -31.770969d0
     UPP_pm6(atomic_number) =  -5.856382d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 3.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END HELIUM
 !-------------------

 !-------------------
 !LITHIUM
 !-------------------
 !Notes for elec eng: IOS = 1, IOP = 0, GSSC = 0, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 3
  !MNDO
   ! Reference:   TAKEN FROM MNDOC BY W.THIEL,QCPE NO.438, V. 2, P.63, (1982). (index = 2)
     mndo_ref_index(atomic_number) = 2
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 0.7023800D0
     p_orb_exp_mndo(atomic_number) = 0.7023800D0
     betas_mndo(atomic_number) = -1.3500400D0
     betap_mndo(atomic_number) = -1.3500400D0
     GSS_mndo(atomic_number) = 7.3000000D0
     GSP_mndo(atomic_number) = 5.4200000D0
     GPP_mndo(atomic_number) = 5.0000000D0
     GP2_mndo(atomic_number) = 4.5200000D0
     HSP_mndo(atomic_number) = 0.8300000D0
     alp_mndo(atomic_number) = 1.2501400D0
     USS_mndo(atomic_number) = -5.1280000D0
     UPP_mndo(atomic_number) = -2.7212000D0
  !MNDOD (same as MNDO)
   ! Reference:   TAKEN FROM MNDOC BY W.THIEL,QCPE NO.438, V. 2, P.63, (1982). (index = 2)
     mndod_ref_index(atomic_number) = 2
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 0.7023800D0
     p_orb_exp_mndod(atomic_number) = 0.7023800D0
     betas_mndod(atomic_number) = -1.3500400D0
     betap_mndod(atomic_number) = -1.3500400D0
     GSS_mndod(atomic_number) = 7.3000000D0
     GSP_mndod(atomic_number) = 5.4200000D0
     GPP_mndod(atomic_number) = 5.0000000D0
     GP2_mndod(atomic_number) = 4.5200000D0
     HSP_mndod(atomic_number) = 0.8300000D0
     alp_mndod(atomic_number) = 1.2501400D0
     USS_mndod(atomic_number) = -5.1280000D0
     UPP_mndod(atomic_number) = -2.7212000D0     
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 0.981041d0
     p_orb_exp_pm6(atomic_number) = 2.953445d0
     betas_pm6(atomic_number) = -2.283946d0
     betap_pm6(atomic_number) = -7.535573d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 11.035907d0
     GSP_pm6(atomic_number) = 19.998647d0
     GPP_pm6(atomic_number) = 11.543650d0
     GP2_pm6(atomic_number) =  9.059036d0
     HSP_pm6(atomic_number) =  1.641886d0
     USS_pm6(atomic_number) = -4.709912d0
     UPP_pm6(atomic_number) = -2.722581d0
     ! alp_pm6 is a guess by AWG; 1.25 is the MNDO value
     alp_pm6(atomic_number) = 1.25d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END LITHIUM
 !-------------------

 !-------------------
 !BERYLLIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0, GSSC = 1, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 4
  !MNDO
   ! Reference: M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978) (index=3)
     mndo_ref_index(atomic_number) = 3
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.0042100D0
     p_orb_exp_mndo(atomic_number) = 1.0042100D0
     betas_mndo(atomic_number) = -4.0170960D0
     betap_mndo(atomic_number) = -4.0170960D0
     GSS_mndo(atomic_number) = 9.00D00
     GSP_mndo(atomic_number) = 7.43D00
     GPP_mndo(atomic_number) = 6.97D00
     GP2_mndo(atomic_number) = 6.22D00
     HSP_mndo(atomic_number) = 1.28D00
     alp_mndo(atomic_number) = 1.6694340D0
     USS_mndo(atomic_number) = -16.6023780D0
     UPP_mndo(atomic_number) = -10.7037710D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978) (index=3)
     mndod_ref_index(atomic_number) = 3
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.0042100D0
     p_orb_exp_mndod(atomic_number) = 1.0042100D0
     betas_mndod(atomic_number) = -4.0170960D0
     betap_mndod(atomic_number) = -4.0170960D0
     GSS_mndod(atomic_number) = 9.00D00
     GSP_mndod(atomic_number) = 7.43D00
     GPP_mndod(atomic_number) = 6.97D00
     GP2_mndod(atomic_number) = 6.22D00
     HSP_mndod(atomic_number) = 1.28D00
     alp_mndod(atomic_number) = 1.6694340D0
     USS_mndod(atomic_number) = -16.6023780D0
     UPP_mndod(atomic_number) = -10.7037710D0
  !AM1
   ! Reference: None
     element_supported_am1(atomic_number) = .false.

  !AM1
   ! Reference: None
     element_supported_am1d(atomic_number) = .false.

  !PM3
   ! Reference:  J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 0.8774390D0
     p_orb_exp_pm3(atomic_number) = 1.5087550D0
     betas_pm3(atomic_number) = -3.9620530D0
     betap_pm3(atomic_number) = -2.7806840D0
     FN1_pm3(1,atomic_number) = 1.6315720D0
     FN2_pm3(1,atomic_number) = 2.6729620D0
     FN3_pm3(1,atomic_number) = 1.7916860D0
     FN1_pm3(2,atomic_number) =-2.1109590D0
     FN2_pm3(2,atomic_number) = 1.9685940D0
     FN3_pm3(2,atomic_number) = 1.7558710D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 9.0128510D0
     GSP_pm3(atomic_number) = 6.5761990D0
     GPP_pm3(atomic_number) = 6.0571820D0
     GP2_pm3(atomic_number) = 9.0052190D0
     HSP_pm3(atomic_number) = 0.5446790D0
     alp_pm3(atomic_number) = 1.5935360D0
     USS_pm3(atomic_number) = -17.2647520D0
     UPP_pm3(atomic_number) = -11.3042430D0
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.212539d0
     p_orb_exp_pm6(atomic_number) = 1.276487d0
     betas_pm6(atomic_number) = -3.199549d0
     betap_pm6(atomic_number) = -4.451920d0
     FN1_pm6(1,atomic_number) = 0.164180d0
     FN2_pm6(1,atomic_number) = 1.704828d0
     FN3_pm6(1,atomic_number) = 1.785591d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  7.552804d0
     GSP_pm6(atomic_number) = 10.203146d0
     GPP_pm6(atomic_number) = 12.862153d0
     GP2_pm6(atomic_number) = 13.602858d0
     HSP_pm6(atomic_number) =  1.501452d0
     USS_pm6(atomic_number) = -16.360315d0
     UPP_pm6(atomic_number) = -16.339216d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.5935360D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END BERYLLIUM
 !-------------------

 !-------------------
 !BORON
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 1, GSSC = 1, GPPC = 0, GSPC = 2, GP2C = 0, HSP = -1
   atomic_number = 5
  !MNDO
   ! Reference: M.J.S. DEWAR, M.L. MCKEE, J. AM. CHEM. SOC., 99, 5231, (1977) (index = 4)
     mndo_ref_index(atomic_number) = 4
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.5068010D0
     p_orb_exp_mndo(atomic_number) = 1.5068010D0
     betas_mndo(atomic_number) = -8.2520540D0
     betap_mndo(atomic_number) = -8.2520540D0
     GSS_mndo(atomic_number) = 10.59D00
     GSP_mndo(atomic_number) = 9.56D00
     GPP_mndo(atomic_number) = 8.86D00
     GP2_mndo(atomic_number) = 7.86D00
     HSP_mndo(atomic_number) = 1.81D00
     alp_mndo(atomic_number) = 2.1349930D0
     USS_mndo(atomic_number) = -34.5471300D0
     UPP_mndo(atomic_number) = -23.1216900D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, M.L. MCKEE, J. AM. CHEM. SOC., 99, 5231, (1977) (index = 4)
     mndod_ref_index(atomic_number) = 4
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.5068010D0
     p_orb_exp_mndod(atomic_number) = 1.5068010D0
     betas_mndod(atomic_number) = -8.2520540D0
     betap_mndod(atomic_number) = -8.2520540D0
     GSS_mndod(atomic_number) = 10.59D00
     GSP_mndod(atomic_number) = 9.56D00
     GPP_mndod(atomic_number) = 8.86D00
     GP2_mndod(atomic_number) = 7.86D00
     HSP_mndod(atomic_number) = 1.81D00
     alp_mndod(atomic_number) = 2.1349930D0
     USS_mndod(atomic_number) = -34.5471300D0
     UPP_mndod(atomic_number) = -23.1216900D0
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.634174d0
     p_orb_exp_pm6(atomic_number) = 1.479195d0
     betas_pm6(atomic_number) = -4.959706d0
     betap_pm6(atomic_number) = -4.656753d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 8.179341d0
     GSP_pm6(atomic_number) = 7.294021d0
     GPP_pm6(atomic_number) = 7.829395d0
     GP2_pm6(atomic_number) = 6.401072d0
     HSP_pm6(atomic_number) = 1.252845d0
     USS_pm6(atomic_number) = -25.967679d0
     UPP_pm6(atomic_number) = -19.115864d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 2.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END BORON
 !-------------------

 !-------------------
 !CARBON
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 2, GSSC = 1, GPPC = -0.5, GSPC = 4, GP2C = 1.5, HSP = -2
   atomic_number = 6
  !MNDO
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndo_ref_index(atomic_number) = 1
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.7875370D0
     p_orb_exp_mndo(atomic_number) = 1.7875370D0
     betas_mndo(atomic_number) = -18.9850440D0
     betap_mndo(atomic_number) = -7.9341220D0
     GSS_mndo(atomic_number) = 12.23D00
     GSP_mndo(atomic_number) = 11.47D00
     GPP_mndo(atomic_number) = 11.08D00
     GP2_mndo(atomic_number) = 9.84D00
     HSP_mndo(atomic_number) = 2.43D00
     alp_mndo(atomic_number) = 2.5463800D0
     USS_mndo(atomic_number) = -52.2797450D0
     UPP_mndo(atomic_number) = -39.2055580D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndod_ref_index(atomic_number) = 1
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.7875370D0
     p_orb_exp_mndod(atomic_number) = 1.7875370D0
     betas_mndod(atomic_number) = -18.9850440D0
     betap_mndod(atomic_number) = -7.9341220D0
     GSS_mndod(atomic_number) = 12.23D00
     GSP_mndod(atomic_number) = 11.47D00
     GPP_mndod(atomic_number) = 11.08D00
     GP2_mndod(atomic_number) = 9.84D00
     HSP_mndod(atomic_number) = 2.43D00
     alp_mndod(atomic_number) = 2.5463800D0
     USS_mndod(atomic_number) = -52.2797450D0
     UPP_mndod(atomic_number) = -39.2055580D0
  !AM1
   ! Reference:  M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (index = 17)
     am1_ref_index(atomic_number) = 17
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.8086650D0
     p_orb_exp_am1(atomic_number) = 1.6851160D0
     betas_am1(atomic_number) = -15.7157830D0
     betap_am1(atomic_number) = -7.7192830D0
     FN1_am1(1,atomic_number) = 0.0113550D0
     FN2_am1(1,atomic_number) = 5.0000000D0
     FN3_am1(1,atomic_number) = 1.6000000D0
     FN1_am1(2,atomic_number) = 0.0459240D0
     FN2_am1(2,atomic_number) = 5.0000000D0
     FN3_am1(2,atomic_number) = 1.8500000D0
     FN1_am1(3,atomic_number) =-0.0200610D0
     FN2_am1(3,atomic_number) = 5.0000000D0
     FN3_am1(3,atomic_number) = 2.0500000D0
     FN1_am1(4,atomic_number) =-0.0012600D0
     FN2_am1(4,atomic_number) = 5.0000000D0
     FN3_am1(4,atomic_number) = 2.6500000D0
     NUM_FN_am1(atomic_number) = 4
     GSS_am1(atomic_number) = 12.2300000D0
     GSP_am1(atomic_number) = 11.4700000D0
     GPP_am1(atomic_number) = 11.0800000D0
     GP2_am1(atomic_number) = 9.8400000D0
     HSP_am1(atomic_number) = 2.4300000D0
     alp_am1(atomic_number) = 2.6482740D0
     USS_am1(atomic_number) = -52.0286580D0
     UPP_am1(atomic_number) = -39.6142390D0
  !AM1D
   ! Reference:  M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (index = 17)
     am1d_ref_index(atomic_number) = 17
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.8086650D0
     p_orb_exp_am1d(atomic_number) = 1.6851160D0
     betas_am1d(atomic_number) = -15.7157830D0
     betap_am1d(atomic_number) = -7.7192830D0
     FN1_am1d(1,atomic_number) = 0.0113550D0
     FN2_am1d(1,atomic_number) = 5.0000000D0
     FN3_am1d(1,atomic_number) = 1.6000000D0
     FN1_am1d(2,atomic_number) = 0.0459240D0
     FN2_am1d(2,atomic_number) = 5.0000000D0
     FN3_am1d(2,atomic_number) = 1.8500000D0
     FN1_am1d(3,atomic_number) =-0.0200610D0
     FN2_am1d(3,atomic_number) = 5.0000000D0
     FN3_am1d(3,atomic_number) = 2.0500000D0
     FN1_am1d(4,atomic_number) =-0.0012600D0
     FN2_am1d(4,atomic_number) = 5.0000000D0
     FN3_am1d(4,atomic_number) = 2.6500000D0
     NUM_FN_am1d(atomic_number) = 4
     GSS_am1d(atomic_number) = 12.2300000D0
     GSP_am1d(atomic_number) = 11.4700000D0
     GPP_am1d(atomic_number) = 11.0800000D0
     GP2_am1d(atomic_number) = 9.8400000D0
     HSP_am1d(atomic_number) = 2.4300000D0
     alp_am1d(atomic_number) = 2.6482740D0
     USS_am1d(atomic_number) = -52.0286580D0
     UPP_am1d(atomic_number) = -39.6142390D0
  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.5650850D0
     p_orb_exp_pm3(atomic_number) = 1.8423450D0
     betas_pm3(atomic_number) = -11.9100150D0
     betap_pm3(atomic_number) = -9.8027550D0
     FN1_pm3(1,atomic_number) = 0.0501070D0
     FN2_pm3(1,atomic_number) = 6.0031650D0
     FN3_pm3(1,atomic_number) = 1.6422140D0
     FN1_pm3(2,atomic_number) = 0.0507330D0
     FN2_pm3(2,atomic_number) = 6.0029790D0
     FN3_pm3(2,atomic_number) = 0.8924880D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 11.2007080D0
     GSP_pm3(atomic_number) = 10.2650270D0
     GPP_pm3(atomic_number) = 10.7962920D0
     GP2_pm3(atomic_number) = 9.0425660D0
     HSP_pm3(atomic_number) = 2.2909800D0
     alp_pm3(atomic_number) = 2.7078070D0
     USS_pm3(atomic_number) = -47.2703200D0
     UPP_pm3(atomic_number) = -36.2669180D0

  !PM3/MM*
   ! Reference: WANG,Q.T. and BRYCE,R.A. (2009) JCTC, 5, p2206
     element_supported_pm3mmx(atomic_number) = .true.
     scale_f1_pm3mmx(1,atomic_number) = 3.4D0
     scale_f2_pm3mmx(1,atomic_number) = 3.9D0

  !PM3/MM* 2nd version
     element_supported_pm3mmx2(atomic_number) = .true.
     rho_pm3mmx2(atomic_number) = 0.931D0
     scale_f1_pm3mmx2(1,atomic_number) = 2.915D0
     scale_f2_pm3mmx2(1,atomic_number) = 3.577D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.047558d0
     p_orb_exp_pm6(atomic_number) = 1.702841d0
     betas_pm6(atomic_number) = -15.385236d0
     betap_pm6(atomic_number) = -7.471929d0
     FN1_pm6(1,atomic_number) = 0.046302d0
     FN2_pm6(1,atomic_number) = 2.100206d0
     FN3_pm6(1,atomic_number) = 1.333959d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 13.335519d0
     GSP_pm6(atomic_number) = 11.528134d0
     GPP_pm6(atomic_number) = 10.778326d0
     GP2_pm6(atomic_number) = 9.486212d0
     HSP_pm6(atomic_number) = 0.717322d0
     USS_pm6(atomic_number) = -51.089653d0
     UPP_pm6(atomic_number) = -39.937920d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 2.7078070D0
     !For pairwise core core terms see section 3 below.

  !PDDG/PM3
   ! Reference Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgpm3_ref_index(atomic_number) = 29
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -113.42824208974D0
     s_orb_exp_pddgpm3(atomic_number) = 1.56786358751710D0
     p_orb_exp_pddgpm3(atomic_number) = 1.84665852120070D0
     betas_pddgpm3(atomic_number) = -11.952818190434D0
     betap_pddgpm3(atomic_number) = -9.9224112120852D0
     FN1_pddgpm3(1,atomic_number) = 0.04890550330860D0
     FN2_pddgpm3(1,atomic_number) = 5.76533980799120D0
     FN3_pddgpm3(1,atomic_number) = 1.68223169651660D0
     FN1_pddgpm3(2,atomic_number) = 0.04769663311610D0
     FN2_pddgpm3(2,atomic_number) = 5.97372073873460D0
     FN3_pddgpm3(2,atomic_number) = 0.89440631619350D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 11.2007080D0
     GSP_pddgpm3(atomic_number) = 10.2650270D0
     GPP_pddgpm3(atomic_number) = 10.7962920D0
     GP2_pddgpm3(atomic_number) = 9.0425660D0
     HSP_pddgpm3(atomic_number) = 2.2909800D0
     alp_pddgpm3(atomic_number) = 2.72577212540530D0
     USS_pddgpm3(atomic_number) = -48.241240946951D0
     UPP_pddgpm3(atomic_number) = -36.461255999939D0
     PDDGC1_pm3(atomic_number) = -0.0007433618099D0
     PDDGC2_pm3(atomic_number) = 0.00098516072940D0
     PDDGE1_pm3(atomic_number) = 0.83691519687330D0
     PDDGE2_pm3(atomic_number) = 1.58523608520060D0

  !PDDG/PM3_08
   ! Reference J. Tirado-Rives et al. J. CHEM. THEO. COMP., 4, 297, (2008) (Index 34)
     pddgpm3_08_ref_index(atomic_number) = 34
     element_supported_pddgpm3_08(atomic_number) = .true.
     elec_eng_pddgpm3_08(atomic_number) = -112.7969d0
     s_orb_exp_pddgpm3_08(atomic_number) = 1.565931d0
     p_orb_exp_pddgpm3_08(atomic_number) = 1.840669d0
     betas_pddgpm3_08(atomic_number) = -11.76394d0
     betap_pddgpm3_08(atomic_number) = -9.883236d0
     FN1_pddgpm3_08(1,atomic_number) = 0.051192d0
     FN2_pddgpm3_08(1,atomic_number) = 5.762521d0
     FN3_pddgpm3_08(1,atomic_number) = 1.706747d0
     FN1_pddgpm3_08(2,atomic_number) = 0.0475d0
     FN2_pddgpm3_08(2,atomic_number) = 6.034004d0
     FN3_pddgpm3_08(2,atomic_number) = 0.932312d0
     FN1_pddgpm3_08(3,atomic_number) = 0.0d0
     FN2_pddgpm3_08(3,atomic_number) = 0.0d0
     FN3_pddgpm3_08(3,atomic_number) = 0.0d0
     FN1_pddgpm3_08(4,atomic_number) = 0.0d0
     FN2_pddgpm3_08(4,atomic_number) = 0.0d0
     FN3_pddgpm3_08(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3_08(atomic_number) = 2
     GSS_pddgpm3_08(atomic_number) = 11.2007080D0
     GSP_pddgpm3_08(atomic_number) = 10.2650270D0
     GPP_pddgpm3_08(atomic_number) = 10.7962920D0
     GP2_pddgpm3_08(atomic_number) = 9.0425660D0
     HSP_pddgpm3_08(atomic_number) = 2.2909800D0
     alp_pddgpm3_08(atomic_number) = 2.719791d0
     USS_pddgpm3_08(atomic_number) = -48.09596d0
     UPP_pddgpm3_08(atomic_number) = -36.38891d0
     PDDGC1_pm3_08(atomic_number) = -0.0007433618099D0
     PDDGC2_pm3_08(atomic_number) = 0.00098516072940D0
     PDDGE1_pm3_08(atomic_number) = 0.83691519687330D0
     PDDGE2_pm3_08(atomic_number) = 1.58523608520060D0

  !PDDG/MNDO
   ! Reference: Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgmndo_ref_index(atomic_number) = 29
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -123.86441152368D0
     s_orb_exp_pddgmndo(atomic_number) = 1.80981702301050D0
     p_orb_exp_pddgmndo(atomic_number) = 1.82500792388930D0
     betas_pddgmndo(atomic_number) = -18.841334137411D0
     betap_pddgmndo(atomic_number) = -7.9222341225346D0
     GSS_pddgmndo(atomic_number) = 12.23D00
     GSP_pddgmndo(atomic_number) = 11.47D00
     GPP_pddgmndo(atomic_number) = 11.08D00
     GP2_pddgmndo(atomic_number) = 9.84D00
     HSP_pddgmndo(atomic_number) = 2.43D00
     alp_pddgmndo(atomic_number) = 2.55552238806810D0
     USS_pddgmndo(atomic_number) = -53.837582488984D0
     UPP_pddgmndo(atomic_number) = -39.936408766823D0
     PDDGC1_mndo(atomic_number) = -0.0068893327627D0
     PDDGC2_mndo(atomic_number) = -0.0277514418977D0
     PDDGE1_mndo(atomic_number) = 1.19245557326430D0
     PDDGE2_mndo(atomic_number) = 1.32952163414800D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 1.85018803D0
     p_orb_exp_rm1(atomic_number) = 1.76830093D0
     betas_rm1(atomic_number) = -15.45932428D0
     betap_rm1(atomic_number) = -8.23608638D0
     FN1_rm1(1,atomic_number) = 0.07462271D0
     FN2_rm1(1,atomic_number) = 5.73921605D0
     FN3_rm1(1,atomic_number) = 1.04396983D0
     FN1_rm1(2,atomic_number) = 0.01177053D0
     FN2_rm1(2,atomic_number) = 6.92401726D0
     FN3_rm1(2,atomic_number) = 1.66159571D0
     FN1_rm1(3,atomic_number) = 0.03720662D0
     FN2_rm1(3,atomic_number) = 6.26158944D0
     FN3_rm1(3,atomic_number) = 1.63158721D0
     FN1_rm1(4,atomic_number) = -0.00270657D0
     FN2_rm1(4,atomic_number) = 9.00003735D0
     FN3_rm1(4,atomic_number) = 2.79557901D0
     NUM_FN_rm1(atomic_number) = 4
     GSS_rm1(atomic_number) = 13.05312440D0
     GSP_rm1(atomic_number) = 11.33479389D0
     GPP_rm1(atomic_number) = 10.95113739D0
     GP2_rm1(atomic_number) = 9.72395099D0
     HSP_rm1(atomic_number) = 1.55215133D0
     alp_rm1(atomic_number) = 2.79282078D0
     USS_rm1(atomic_number) = -51.72556032D0
     UPP_rm1(atomic_number) = -39.40728943D0

 !-------------------
 !END CARBON
 !-------------------

 !-------------------
 !NITROGEN
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 3, GSSC = 1, GPPC = -1.5, GSPC = 6, GP2C = 4.5, HSP = -3
   atomic_number = 7
  !MNDO
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndo_ref_index(atomic_number) = 1
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.2556140D0
     p_orb_exp_mndo(atomic_number) = 2.2556140D0
     betas_mndo(atomic_number) = -20.4957580D0
     betap_mndo(atomic_number) = -20.4957580D0
     GSS_mndo(atomic_number) = 13.59D00
     GSP_mndo(atomic_number) = 12.66D00
     GPP_mndo(atomic_number) = 12.98D00
     GP2_mndo(atomic_number) = 11.59D00
     HSP_mndo(atomic_number) = 3.14D00
     alp_mndo(atomic_number) = 2.8613420D0
     USS_mndo(atomic_number) = -71.9321220D0
     UPP_mndo(atomic_number) = -57.1723190D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndod_ref_index(atomic_number) = 1
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.2556140D0
     p_orb_exp_mndod(atomic_number) = 2.2556140D0
     betas_mndod(atomic_number) = -20.4957580D0
     betap_mndod(atomic_number) = -20.4957580D0
     GSS_mndod(atomic_number) = 13.59D00
     GSP_mndod(atomic_number) = 12.66D00
     GPP_mndod(atomic_number) = 12.98D00
     GP2_mndod(atomic_number) = 11.59D00
     HSP_mndod(atomic_number) = 3.14D00
     alp_mndod(atomic_number) = 2.8613420D0
     USS_mndod(atomic_number) = -71.9321220D0
     UPP_mndod(atomic_number) = -57.1723190D0
  !AM1
   ! Reference:  M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (index = 17)
     am1_ref_index(atomic_number) = 17
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 2.3154100D0
     p_orb_exp_am1(atomic_number) = 2.1579400D0
     betas_am1(atomic_number) = -20.2991100D0
     betap_am1(atomic_number) = -18.2386660D0
     FN1_am1(1,atomic_number) = 0.0252510D0
     FN2_am1(1,atomic_number) = 5.0000000D0
     FN3_am1(1,atomic_number) = 1.5000000D0
     FN1_am1(2,atomic_number) = 0.0289530D0
     FN2_am1(2,atomic_number) = 5.0000000D0
     FN3_am1(2,atomic_number) = 2.1000000D0
     FN1_am1(3,atomic_number) =-0.0058060D0
     FN2_am1(3,atomic_number) = 2.0000000D0
     FN3_am1(3,atomic_number) = 2.4000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 13.5900000D0
     GSP_am1(atomic_number) = 12.6600000D0
     GPP_am1(atomic_number) = 12.9800000D0
     GP2_am1(atomic_number) = 11.5900000D0
     HSP_am1(atomic_number) = 3.1400000D0
     alp_am1(atomic_number) = 2.9472860D0
     USS_am1(atomic_number) = -71.8600000D0
     UPP_am1(atomic_number) = -57.1675810D0

  !AM1D
   ! Reference:  M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (index = 17)
     am1d_ref_index(atomic_number) = 17
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 2.3154100D0
     p_orb_exp_am1d(atomic_number) = 2.1579400D0
     betas_am1d(atomic_number) = -20.2991100D0
     betap_am1d(atomic_number) = -18.2386660D0
     FN1_am1d(1,atomic_number) = 0.0252510D0
     FN2_am1d(1,atomic_number) = 5.0000000D0
     FN3_am1d(1,atomic_number) = 1.5000000D0
     FN1_am1d(2,atomic_number) = 0.0289530D0
     FN2_am1d(2,atomic_number) = 5.0000000D0
     FN3_am1d(2,atomic_number) = 2.1000000D0
     FN1_am1d(3,atomic_number) =-0.0058060D0
     FN2_am1d(3,atomic_number) = 2.0000000D0
     FN3_am1d(3,atomic_number) = 2.4000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 3
     GSS_am1d(atomic_number) = 13.5900000D0
     GSP_am1d(atomic_number) = 12.6600000D0
     GPP_am1d(atomic_number) = 12.9800000D0
     GP2_am1d(atomic_number) = 11.5900000D0
     HSP_am1d(atomic_number) = 3.1400000D0
     alp_am1d(atomic_number) = 2.9472860D0
     USS_am1d(atomic_number) = -71.8600000D0
     UPP_am1d(atomic_number) = -57.1675810D0

  !PM3
   ! Reference:  J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989). (index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.0280940D0
     p_orb_exp_pm3(atomic_number) = 2.3137280D0
     betas_pm3(atomic_number) = -14.0625210D0
     betap_pm3(atomic_number) = -20.0438480D0
     FN1_pm3(1,atomic_number) = 1.5016740D0
     FN2_pm3(1,atomic_number) = 5.9011480D0
     FN3_pm3(1,atomic_number) = 1.7107400D0
     FN1_pm3(2,atomic_number) = -1.5057720D0
     FN2_pm3(2,atomic_number) = 6.0046580D0
     FN3_pm3(2,atomic_number) = 1.7161490D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 11.9047870D0
     GSP_pm3(atomic_number) = 7.3485650D0
     GPP_pm3(atomic_number) = 11.7546720D0
     GP2_pm3(atomic_number) = 10.8072770D0
     HSP_pm3(atomic_number) = 1.1367130D0
     alp_pm3(atomic_number) = 2.8305450D0
     USS_pm3(atomic_number) = -49.3356720D0
     UPP_pm3(atomic_number) = -47.5097360D0

  !PM3/MM*
   ! Reference: WANG,Q.T. and BRYCE,R.A. (2009) JCTC, 5, p2206
     element_supported_pm3mmx(atomic_number) = .true.
     scale_f1_pm3mmx(1,atomic_number) = 2.9D0
     scale_f2_pm3mmx(1,atomic_number) = 3.4D0

  !PM3/MM* 2nd version
     element_supported_pm3mmx2(atomic_number) = .true.
     rho_pm3mmx2(atomic_number) = 0.395D0
     scale_f1_pm3mmx2(1,atomic_number) = 2.852D0
     scale_f2_pm3mmx2(1,atomic_number) = 3.997D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.380406d0
     p_orb_exp_pm6(atomic_number) = 1.999246d0
     betas_pm6(atomic_number) = -17.979377d0
     betap_pm6(atomic_number) = -15.055017d0
     FN1_pm6(1,atomic_number) = -0.001436d0
     FN2_pm6(1,atomic_number) = 0.495196d0
     FN3_pm6(1,atomic_number) = 1.704857d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 12.357026d0
     GSP_pm6(atomic_number) =  9.636190d0
     GPP_pm6(atomic_number) = 12.570756d0
     GP2_pm6(atomic_number) = 10.576425d0
     HSP_pm6(atomic_number) =  2.871545d0
     USS_pm6(atomic_number) = -57.784823d0
     UPP_pm6(atomic_number) = -49.893036d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 2.8305450D0
     !For pairwise core core terms see section 3 below.

  !PDDG/PM3
   ! Reference Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgpm3_ref_index(atomic_number) = 29
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -158.41620481951D0
     s_orb_exp_pddgpm3(atomic_number) = 2.03580684361910D0
     p_orb_exp_pddgpm3(atomic_number) = 2.32432725808280D0
     betas_pddgpm3(atomic_number) = -14.117229602371D0
     betap_pddgpm3(atomic_number) = -19.938508878969D0
     FN1_pddgpm3(1,atomic_number) = 1.51332030575080D0
     FN2_pddgpm3(1,atomic_number) = 5.90439402634500D0
     FN3_pddgpm3(1,atomic_number) = 1.72837621719040D0
     FN1_pddgpm3(2,atomic_number) = -1.5118916914302D0
     FN2_pddgpm3(2,atomic_number) = 6.03001440913320D0
     FN3_pddgpm3(2,atomic_number) = 1.73410826456840D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 11.9047870D0
     GSP_pddgpm3(atomic_number) = 7.3485650D0
     GPP_pddgpm3(atomic_number) = 11.7546720D0
     GP2_pddgpm3(atomic_number) = 10.8072770D0
     HSP_pddgpm3(atomic_number) = 1.1367130D0
     alp_pddgpm3(atomic_number) = 2.84912399973850D0
     USS_pddgpm3(atomic_number) = -49.454546358059D0
     UPP_pddgpm3(atomic_number) = -47.757406358412D0
     PDDGC1_pm3(atomic_number) = -0.0031600751673D0
     PDDGC2_pm3(atomic_number) = 0.01250092178130D0
     PDDGE1_pm3(atomic_number) = 1.00417177651930D0
     PDDGE2_pm3(atomic_number) = 1.51633618021020D0

  !PDDG/PM3_08
   ! Reference J. Tirado-Rives et al. J. CHEM. THEO. COMP., 4, 297, (2008) (Index 34)
     pddgpm3_08_ref_index(atomic_number) = 34
     element_supported_pddgpm3_08(atomic_number) = .true.
     elec_eng_pddgpm3_08(atomic_number) = -157.6938d0
     s_orb_exp_pddgpm3_08(atomic_number) = 2.026598d0
     p_orb_exp_pddgpm3_08(atomic_number) = 2.334183d0
     betas_pddgpm3_08(atomic_number) = -14.08164d0
     betap_pddgpm3_08(atomic_number) = -19.69538d0
     FN1_pddgpm3_08(1,atomic_number) = 1.508427d0
     FN2_pddgpm3_08(1,atomic_number) = 5.957281d0
     FN3_pddgpm3_08(1,atomic_number) = 1.72277d0
     FN1_pddgpm3_08(2,atomic_number) = -1.508203d0
     FN2_pddgpm3_08(2,atomic_number) = 6.025113d0
     FN3_pddgpm3_08(2,atomic_number) = 1.731257d0
     FN1_pddgpm3_08(3,atomic_number) = 0.0d0
     FN2_pddgpm3_08(3,atomic_number) = 0.0d0
     FN3_pddgpm3_08(3,atomic_number) = 0.0d0
     FN1_pddgpm3_08(4,atomic_number) = 0.0d0
     FN2_pddgpm3_08(4,atomic_number) = 0.0d0
     FN3_pddgpm3_08(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3_08(atomic_number) = 2
     GSS_pddgpm3_08(atomic_number) = 11.9047870D0
     GSP_pddgpm3_08(atomic_number) = 7.3485650D0
     GPP_pddgpm3_08(atomic_number) = 11.7546720D0
     GP2_pddgpm3_08(atomic_number) = 10.8072770D0
     HSP_pddgpm3_08(atomic_number) = 1.1367130D0
     alp_pddgpm3_08(atomic_number) = 2.846987d0
     USS_pddgpm3_08(atomic_number) = -49.42949d0
     UPP_pddgpm3_08(atomic_number) = -47.64097d0
     PDDGC1_pm3_08(atomic_number) = -0.003229d0
     PDDGC2_pm3_08(atomic_number) = 0.012714d0
     PDDGE1_pm3_08(atomic_number) = 1.007904d0
     PDDGE2_pm3_08(atomic_number) = 1.511671d0

  !PDDG/MNDO
   ! Reference: Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgmndo_ref_index(atomic_number) = 29
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -206.46662581723320D0
     s_orb_exp_pddgmndo(atomic_number) = 2.23142379586030D0
     p_orb_exp_pddgmndo(atomic_number) = 2.25345995688440D0
     betas_pddgmndo(atomic_number) = -20.37577411084280D0
     betap_pddgmndo(atomic_number) = -21.08537341050740D0
     GSS_pddgmndo(atomic_number) = 13.59D00
     GSP_pddgmndo(atomic_number) = 12.66D00
     GPP_pddgmndo(atomic_number) = 12.98D00
     GP2_pddgmndo(atomic_number) = 11.59D00
     HSP_pddgmndo(atomic_number) = 3.14D00
     alp_pddgmndo(atomic_number) = 2.84367788492060D0
     USS_pddgmndo(atomic_number) = -71.87189435530930D0
     UPP_pddgmndo(atomic_number) = -58.21661676886340D0
     PDDGC1_mndo(atomic_number) = 0.03502693409010D0
     PDDGC2_mndo(atomic_number) = -0.00172140634590D0
     PDDGE1_mndo(atomic_number) = 1.01162987147870D0
     PDDGE2_mndo(atomic_number) = 2.27842256966070D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 2.37447159D0
     p_orb_exp_rm1(atomic_number) = 1.97812569D0
     betas_rm1(atomic_number) = -20.87124548D0
     betap_rm1(atomic_number) = -16.67171853D0
     FN1_rm1(1,atomic_number) = 0.06073380D0
     FN2_rm1(1,atomic_number) = 4.58892946D0
     FN3_rm1(1,atomic_number) = 1.37873881D0
     FN1_rm1(2,atomic_number) = 0.02438558D0
     FN2_rm1(2,atomic_number) = 4.62730519D0
     FN3_rm1(2,atomic_number) = 2.08370698D0
     FN1_rm1(3,atomic_number) = -0.02283430D0
     FN2_rm1(3,atomic_number) = 2.05274659D0
     FN3_rm1(3,atomic_number) = 1.86763816D0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 3
     GSS_rm1(atomic_number) = 13.08736234D0
     GSP_rm1(atomic_number) = 13.21226834D0
     GPP_rm1(atomic_number) = 13.69924324D0
     GP2_rm1(atomic_number) = 11.94103953D0
     HSP_rm1(atomic_number) = 5.00000846D0
     alp_rm1(atomic_number) = 2.96422542D0
     USS_rm1(atomic_number) = -70.85123715D0
     UPP_rm1(atomic_number) = -57.97730920D0

 !-------------------
 !END NITROGEN
 !-------------------

 !-------------------
 !OXYGEN
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 4, GSSC = 1, GPPC = -0.5, GSPC = 8, GP2C = 6.5, HSP = -4
   atomic_number = 8
  !MNDO
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndo_ref_index(atomic_number) = 1
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.6999050D0
     p_orb_exp_mndo(atomic_number) = 2.6999050D0
     betas_mndo(atomic_number) = -32.6880820D0
     betap_mndo(atomic_number) = -32.6880820D0
     GSS_mndo(atomic_number) = 15.42D00
     GSP_mndo(atomic_number) = 14.48D00
     GPP_mndo(atomic_number) = 14.52D00
     GP2_mndo(atomic_number) = 12.98D00
     HSP_mndo(atomic_number) = 3.94D00
     alp_mndo(atomic_number) = 3.1606040D0
     USS_mndo(atomic_number) = -99.6443090D0
     UPP_mndo(atomic_number) = -77.7974720D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977) (index = 1)
     mndod_ref_index(atomic_number) = 1
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.6999050D0
     p_orb_exp_mndod(atomic_number) = 2.6999050D0
     betas_mndod(atomic_number) = -32.6880820D0
     betap_mndod(atomic_number) = -32.6880820D0
     GSS_mndod(atomic_number) = 15.42D00
     GSP_mndod(atomic_number) = 14.48D00
     GPP_mndod(atomic_number) = 14.52D00
     GP2_mndod(atomic_number) = 12.98D00
     HSP_mndod(atomic_number) = 3.94D00
     alp_mndod(atomic_number) = 3.1606040D0
     USS_mndod(atomic_number) = -99.6443090D0
     UPP_mndod(atomic_number) = -77.7974720D0
  !AM1
   ! Reference: M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985) (index = 17)
     am1_ref_index(atomic_number) = 17
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 3.1080320D0
     p_orb_exp_am1(atomic_number) = 2.5240390D0
     betas_am1(atomic_number) = -29.2727730D0
     betap_am1(atomic_number) = -29.2727730D0
     FN1_am1(1,atomic_number) = 0.2809620D0
     FN2_am1(1,atomic_number) = 5.0000000D0
     FN3_am1(1,atomic_number) = 0.8479180D0
     FN1_am1(2,atomic_number) = 0.0814300D0
     FN2_am1(2,atomic_number) = 7.0000000D0
     FN3_am1(2,atomic_number) = 1.4450710D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 2
     GSS_am1(atomic_number) = 15.4200000D0
     GSP_am1(atomic_number) = 14.4800000D0
     GPP_am1(atomic_number) = 14.5200000D0
     GP2_am1(atomic_number) = 12.9800000D0
     HSP_am1(atomic_number) = 3.9400000D0
     alp_am1(atomic_number) = 4.4553710D0
     USS_am1(atomic_number) = -97.8300000D0
     UPP_am1(atomic_number) = -78.2623800D0

  !AM1D
   ! Reference: K. Nam, Q. Cui, J. Gao, D. York. J. CHEM. THEO. COMP., 3, 486, (2007) (102) AM1/d-PhoT
     am1d_ref_index(atomic_number) = 102
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 3.057965D0
     p_orb_exp_am1d(atomic_number) = 2.515332D0
     betas_am1d(atomic_number) = -29.472306D0
     betap_am1d(atomic_number) = -28.515785D0
     FN1_am1d(1,atomic_number) = 0.288526D0
     FN2_am1d(1,atomic_number) = 4.883265D0
     FN3_am1d(1,atomic_number) = 0.850910D0
     FN1_am1d(2,atomic_number) = 0.061586D0
     FN2_am1d(2,atomic_number) = 4.435791D0
     FN3_am1d(2,atomic_number) = 1.353681D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 2
     GSS_am1d(atomic_number) = 14.234714D0
     GSP_am1d(atomic_number) = 14.539451D0
     GPP_am1d(atomic_number) = 14.454530D0
     GP2_am1d(atomic_number) = 12.942259D0
     HSP_am1d(atomic_number) = 4.339705D0
     alp_am1d(atomic_number) = 4.404417D0
     USS_am1d(atomic_number) = -96.760676D0
     UPP_am1d(atomic_number) = -78.776203D0

  !PM3
   ! Reference:  J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 3.7965440D0
     p_orb_exp_pm3(atomic_number) = 2.3894020D0
     betas_pm3(atomic_number) = -45.2026510D0
     betap_pm3(atomic_number) = -24.7525150D0
     FN1_pm3(1,atomic_number) = -1.1311280D0
     FN2_pm3(1,atomic_number) = 6.0024770D0
     FN3_pm3(1,atomic_number) = 1.6073110D0
     FN1_pm3(2,atomic_number) = 1.1378910D0
     FN2_pm3(2,atomic_number) = 5.9505120D0
     FN3_pm3(2,atomic_number) = 1.5983950D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 15.7557600D0
     GSP_pm3(atomic_number) = 10.6211600D0
     GPP_pm3(atomic_number) = 13.6540160D0
     GP2_pm3(atomic_number) = 12.4060950D0
     HSP_pm3(atomic_number) = 0.5938830D0
     alp_pm3(atomic_number) = 3.2171020D0
     USS_pm3(atomic_number) = -86.9930020D0
     UPP_pm3(atomic_number) = -71.8795800D0

  !PM3/MM*
   ! Reference: WANG,Q.T. and BRYCE,R.A. (2009) JCTC, 5, p2206
     element_supported_pm3mmx(atomic_number) = .true.
     scale_f1_pm3mmx(1,atomic_number) = 3.6D0
     scale_f2_pm3mmx(1,atomic_number) = 3.6D0

  !PM3/MM* 2nd version
     element_supported_pm3mmx2(atomic_number) = .true.
     rho_pm3mmx2(atomic_number) = 0.875D0
     scale_f1_pm3mmx2(1,atomic_number) = 2.682D0
     scale_f2_pm3mmx2(1,atomic_number) = 4.015D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 5.421751d0
     p_orb_exp_pm6(atomic_number) = 2.270960d0
     betas_pm6(atomic_number) = -65.635137d0
     betap_pm6(atomic_number) = -21.622604d0
     FN1_pm6(1,atomic_number) = -0.017771d0
     FN2_pm6(1,atomic_number) = 3.058310d0
     FN3_pm6(1,atomic_number) = 1.896435d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  11.304042d0
     GSP_pm6(atomic_number) =  15.807424d0
     GPP_pm6(atomic_number) =  13.618205d0
     GP2_pm6(atomic_number) =  10.332765d0
     HSP_pm6(atomic_number) =   5.010801d0
     USS_pm6(atomic_number) = -91.678761d0
     UPP_pm6(atomic_number) = -70.460949d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 3.2171020D0
     !For pairwise core core terms see section 3 below.


  !PDDG/PM3
   ! Reference Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgpm3_ref_index(atomic_number) = 29
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -292.18876564023D0
     s_orb_exp_pddgpm3(atomic_number) = 3.81456531095080D0
     p_orb_exp_pddgpm3(atomic_number) = 2.31801122165690D0
     betas_pddgpm3(atomic_number) = -44.874553472211D0
     betap_pddgpm3(atomic_number) = -24.601939339720D0
     FN1_pddgpm3(1,atomic_number) = -1.1384554300359D0
     FN2_pddgpm3(1,atomic_number) = 6.00004254473730D0
     FN3_pddgpm3(1,atomic_number) = 1.62236167639400D0
     FN1_pddgpm3(2,atomic_number) = 1.14600702743950D0
     FN2_pddgpm3(2,atomic_number) = 5.96349383486760D0
     FN3_pddgpm3(2,atomic_number) = 1.61478803799000D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 15.7557600D0
     GSP_pddgpm3(atomic_number) = 10.6211600D0
     GPP_pddgpm3(atomic_number) = 13.6540160D0
     GP2_pddgpm3(atomic_number) = 12.4060950D0
     HSP_pddgpm3(atomic_number) = 0.5938830D0
     alp_pddgpm3(atomic_number) = 3.22530882036500D0
     USS_pddgpm3(atomic_number) = -87.412505208248D0
     UPP_pddgpm3(atomic_number) = -72.183069806393D0
     PDDGC1_pm3(atomic_number) = -0.00099962677420D0
     PDDGC2_pm3(atomic_number) = -0.00152161350520D0
     PDDGE1_pm3(atomic_number) = 1.36068502987020D0
     PDDGE2_pm3(atomic_number) = 1.36640659538530D0

  !PDDG/PM3_08
   ! Reference J. Tirado-Rives et al. J. CHEM. THEO. COMP., 4, 297, (2008) (Index 34)
     pddgpm3_08_ref_index(atomic_number) = 34
     element_supported_pddgpm3_08(atomic_number) = .true.
     elec_eng_pddgpm3_08(atomic_number) = -294.7122d0
     s_orb_exp_pddgpm3_08(atomic_number) = 3.811544d0
     p_orb_exp_pddgpm3_08(atomic_number) = 2.302506d0
     betas_pddgpm3_08(atomic_number) = -44.6312d0
     betap_pddgpm3_08(atomic_number) = -24.71147d0
     FN1_pddgpm3_08(1,atomic_number) = -1.135968d0
     FN2_pddgpm3_08(1,atomic_number) = 5.988441d0
     FN3_pddgpm3_08(1,atomic_number) = 1.620971d0
     FN1_pddgpm3_08(2,atomic_number) = 1.146007d0
     FN2_pddgpm3_08(2,atomic_number) = 5.963494d0
     FN3_pddgpm3_08(2,atomic_number) = 1.614788d0
     FN1_pddgpm3_08(3,atomic_number) = 0.0d0
     FN2_pddgpm3_08(3,atomic_number) = 0.0d0
     FN3_pddgpm3_08(3,atomic_number) = 0.0d0
     FN1_pddgpm3_08(4,atomic_number) = 0.0d0
     FN2_pddgpm3_08(4,atomic_number) = 0.0d0
     FN3_pddgpm3_08(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3_08(atomic_number) = 2
     GSS_pddgpm3_08(atomic_number) = 15.7557600D0
     GSP_pddgpm3_08(atomic_number) = 10.6211600D0
     GPP_pddgpm3_08(atomic_number) = 13.6540160D0
     GP2_pddgpm3_08(atomic_number) = 12.4060950D0
     HSP_pddgpm3_08(atomic_number) = 0.5938830D0
     alp_pddgpm3_08(atomic_number) = 3.221091d0
     USS_pddgpm3_08(atomic_number) = -87.92097d0
     UPP_pddgpm3_08(atomic_number) = -72.4924d0
     PDDGC1_pm3_08(atomic_number) = -0.00099962677420D0
     PDDGC2_pm3_08(atomic_number) = -0.00152161350520D0
     PDDGE1_pm3_08(atomic_number) = 1.36068502987020D0
     PDDGE2_pm3_08(atomic_number) = 1.36640659538530D0

  !PDDG/MNDO
   ! Reference: Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 29)
     pddgmndo_ref_index(atomic_number) = 29
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -310.87974465179D0
     s_orb_exp_pddgmndo(atomic_number) = 2.56917199926090D0
     p_orb_exp_pddgmndo(atomic_number) = 2.69715151721790D0
     betas_pddgmndo(atomic_number) = -33.606335624658D0
     betap_pddgmndo(atomic_number) = -27.984442042827D0
     GSS_pddgmndo(atomic_number) = 15.42D00
     GSP_pddgmndo(atomic_number) = 14.48D00
     GPP_pddgmndo(atomic_number) = 14.52D00
     GP2_pddgmndo(atomic_number) = 12.98D00
     HSP_pddgmndo(atomic_number) = 3.94D00
     alp_pddgmndo(atomic_number) = 3.23884200872830D0
     USS_pddgmndo(atomic_number) = -97.884970179897D0
     UPP_pddgmndo(atomic_number) = -77.342673804072D0
     PDDGC1_mndo(atomic_number) = 0.08634413812890D0
     PDDGC2_mndo(atomic_number) = 0.03040342779910D0
     PDDGE1_mndo(atomic_number) = 0.72540784783600D0
     PDDGE2_mndo(atomic_number) = 0.70972848794410D0

  !PM3CARB-1
  !Reference: McNamara, J.P., Muslim, A.M., Abdel-Aal, H., Wang, H., Mohr, M., Hillier, I.H., Bryce, R.A.,
  !           2004, Chem Phys Lett, 394, 429-436 (Index 28)
  !Note: PM3CARB-1 is not a complete parameter set, it is simply PM3 with a few parameters changed for
  !      Oxygen and Hydrogen. So when loading these parameters into the correct structure the code
  !      needs to load pm3 parameters except for those that are different for PM3CARB-1
    pm3carb1_ref_index(atomic_number) = 28
    USS_PM3CARB1(atomic_number) = -90.938073D0
    UPP_PM3CARB1(atomic_number) = -76.932200D0
    betas_PM3CARB1(atomic_number) = -44.449581D0
    betap_PM3CARB1(atomic_number) = -35.343869D0
    alp_PM3CARB1(atomic_number) = 3.031867D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 3.17936914D0
     p_orb_exp_rm1(atomic_number) = 2.55361907D0
     betas_rm1(atomic_number) = -29.85101212D0
     betap_rm1(atomic_number) = -29.15101314D0
     FN1_rm1(1,atomic_number) = 0.23093552D0
     FN2_rm1(1,atomic_number) = 5.21828736D0
     FN3_rm1(1,atomic_number) = 0.90363555D0
     FN1_rm1(2,atomic_number) = 0.05859873D0
     FN2_rm1(2,atomic_number) = 7.42932932D0
     FN3_rm1(2,atomic_number) = 1.51754610D0
     FN1_rm1(3,atomic_number) = 0.0d0
     FN2_rm1(3,atomic_number) = 0.0d0
     FN3_rm1(3,atomic_number) = 0.0d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 2
     GSS_rm1(atomic_number) = 14.00242788D0
     GSP_rm1(atomic_number) = 14.95625043D0
     GPP_rm1(atomic_number) = 14.14515138D0
     GP2_rm1(atomic_number) = 12.70325497D0
     HSP_rm1(atomic_number) = 3.93217161D0
     alp_rm1(atomic_number) = 4.17196717D0
     USS_rm1(atomic_number) = -96.94948069D0
     UPP_rm1(atomic_number) = -77.89092978D0

 !-------------------
 !END OXYGEN
 !-------------------

 !-------------------
 !FLUORINE
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 5, GSSC = 1, GPPC = 0, GSPC = 10, GP2C = 10, HSP = -5
   atomic_number = 9
  !MNDO
   ! Reference: M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978) (index = 3)
     mndo_ref_index(atomic_number) = 3
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.8484870D0
     p_orb_exp_mndo(atomic_number) = 2.8484870D0
     betas_mndo(atomic_number) = -48.2904660D0
     betap_mndo(atomic_number) = -36.5085400D0
     GSS_mndo(atomic_number) = 16.92D00
     GSP_mndo(atomic_number) = 17.25D00
     GPP_mndo(atomic_number) = 16.71D00
     GP2_mndo(atomic_number) = 14.91D00
     HSP_mndo(atomic_number) = 4.83D00
     alp_mndo(atomic_number) = 3.4196606D0
     USS_mndo(atomic_number) = -131.0715480D0
     UPP_mndo(atomic_number) = -105.7821370D0
  !MNDOD (same as MNDO)
   ! Reference: M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978) (index = 3)
     mndod_ref_index(atomic_number) = 3
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.8484870D0
     p_orb_exp_mndod(atomic_number) = 2.8484870D0
     betas_mndod(atomic_number) = -48.2904660D0
     betap_mndod(atomic_number) = -36.5085400D0
     GSS_mndod(atomic_number) = 16.92D00
     GSP_mndod(atomic_number) = 17.25D00
     GPP_mndod(atomic_number) = 16.71D00
     GP2_mndod(atomic_number) = 14.91D00
     HSP_mndod(atomic_number) = 4.83D00
     alp_mndod(atomic_number) = 3.4196606D0
     USS_mndod(atomic_number) = -131.0715480D0
     UPP_mndod(atomic_number) = -105.7821370D0
  !AM1
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1_ref_index(atomic_number) = 18
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 3.7700820D0
     p_orb_exp_am1(atomic_number) = 2.4946700D0
     betas_am1(atomic_number) = -69.5902770D0
     betap_am1(atomic_number) = -27.9223600D0
     FN1_am1(1,atomic_number) = 0.2420790D0
     FN2_am1(1,atomic_number) = 4.8000000D0
     FN3_am1(1,atomic_number) = 0.9300000D0
     FN1_am1(2,atomic_number) = 0.0036070D0
     FN2_am1(2,atomic_number) = 4.6000000D0
     FN3_am1(2,atomic_number) = 1.6600000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 2
     GSS_am1(atomic_number) = 16.9200000D0
     GSP_am1(atomic_number) = 17.2500000D0
     GPP_am1(atomic_number) = 16.7100000D0
     GP2_am1(atomic_number) = 14.9100000D0
     HSP_am1(atomic_number) = 4.8300000D0
     alp_am1(atomic_number) = 5.5178000D0
     USS_am1(atomic_number) = -136.1055790D0
     UPP_am1(atomic_number) = -104.8898850D0

  !AM1D
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1d_ref_index(atomic_number) = 18
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 3.7700820D0
     p_orb_exp_am1d(atomic_number) = 2.4946700D0
     betas_am1d(atomic_number) = -69.5902770D0
     betap_am1d(atomic_number) = -27.9223600D0
     FN1_am1d(1,atomic_number) = 0.2420790D0
     FN2_am1d(1,atomic_number) = 4.8000000D0
     FN3_am1d(1,atomic_number) = 0.9300000D0
     FN1_am1d(2,atomic_number) = 0.0036070D0
     FN2_am1d(2,atomic_number) = 4.6000000D0
     FN3_am1d(2,atomic_number) = 1.6600000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 2
     GSS_am1d(atomic_number) = 16.9200000D0
     GSP_am1d(atomic_number) = 17.2500000D0
     GPP_am1d(atomic_number) = 16.7100000D0
     GP2_am1d(atomic_number) = 14.9100000D0
     HSP_am1d(atomic_number) = 4.8300000D0
     alp_am1d(atomic_number) = 5.5178000D0
     USS_am1d(atomic_number) = -136.1055790D0
     UPP_am1d(atomic_number) = -104.8898850D0
  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 4.7085550D0
     p_orb_exp_pm3(atomic_number) = 2.4911780D0
     betas_pm3(atomic_number) = -48.4059390D0
     betap_pm3(atomic_number) = -27.7446600D0
     FN1_pm3(1,atomic_number) = -0.0121660D0
     FN2_pm3(1,atomic_number) = 6.0235740D0
     FN3_pm3(1,atomic_number) = 1.8568590D0
     FN1_pm3(2,atomic_number) = -0.0028520D0
     FN2_pm3(2,atomic_number) = 6.0037170D0
     FN3_pm3(2,atomic_number) = 2.6361580D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 10.4966670D0
     GSP_pm3(atomic_number) = 16.0736890D0
     GPP_pm3(atomic_number) = 14.8172560D0
     GP2_pm3(atomic_number) = 14.4183930D0
     HSP_pm3(atomic_number) = 0.7277630D0
     alp_pm3(atomic_number) = 3.3589210D0
     USS_pm3(atomic_number) = -110.4353030D0
     UPP_pm3(atomic_number) = -105.6850470D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgpm3_ref_index(atomic_number) = 30
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -442.457133D0
     s_orb_exp_pddgpm3(atomic_number) = 5.538033D0
     p_orb_exp_pddgpm3(atomic_number) = 2.538066D0
     betas_pddgpm3(atomic_number) = -50.937301D0
     betap_pddgpm3(atomic_number) = -31.636976D0
     FN1_pddgpm3(1,atomic_number) = -0.008079D0
     FN2_pddgpm3(1,atomic_number) = 5.938969D0
     FN3_pddgpm3(1,atomic_number) = 1.863949D0
     FN1_pddgpm3(2,atomic_number) = -0.002659D0
     FN2_pddgpm3(2,atomic_number) = 5.925105D0
     FN3_pddgpm3(2,atomic_number) = 2.388864D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 10.4966670D0
     GSP_pddgpm3(atomic_number) = 16.0736890D0
     GPP_pddgpm3(atomic_number) = 14.8172560D0
     GP2_pddgpm3(atomic_number) = 14.4183930D0
     HSP_pddgpm3(atomic_number) = 0.7277630D0
     alp_pddgpm3(atomic_number) = 3.200571D0
     USS_pddgpm3(atomic_number) = -111.400432D0
     UPP_pddgpm3(atomic_number) = -106.395264D0
     PDDGC1_pm3(atomic_number) = -0.012866D0
     PDDGC2_pm3(atomic_number) = 0.007315D0
     PDDGE1_pm3(atomic_number) = 1.305681D0
     PDDGE2_pm3(atomic_number) = 1.842572D0

  !PDDG/MNDO
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgmndo_ref_index(atomic_number) = 30
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -488.703243D0
     s_orb_exp_pddgmndo(atomic_number) = 4.328519D0
     p_orb_exp_pddgmndo(atomic_number) = 2.905042D0
     betas_pddgmndo(atomic_number) = -67.827612D0
     betap_pddgmndo(atomic_number) = -40.924818D0
     GSS_pddgmndo(atomic_number) = 16.92D00
     GSP_pddgmndo(atomic_number) = 17.25D00
     GPP_pddgmndo(atomic_number) = 16.71D00
     GP2_pddgmndo(atomic_number) = 14.91D00
     HSP_pddgmndo(atomic_number) = 4.83D00
     alp_pddgmndo(atomic_number) = 3.322382D0
     USS_pddgmndo(atomic_number) = -134.220379D0
     UPP_pddgmndo(atomic_number) = -107.155961D0
     PDDGC1_mndo(atomic_number) = -0.011579D0
     PDDGC2_mndo(atomic_number) = -0.012943D0
     PDDGE1_mndo(atomic_number) = 0.834606D0
     PDDGE2_mndo(atomic_number) = 1.875603D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 4.40337913d0
     p_orb_exp_rm1(atomic_number) = 2.64841556d0
     betas_rm1(atomic_number) = -70.00000512d0
     betap_rm1(atomic_number) = -32.67982711d0
     FN1_rm1(1,atomic_number) = 0.40302025d0
     FN2_rm1(1,atomic_number) = 7.20441959d0
     FN3_rm1(1,atomic_number) = 0.81653013d0
     FN1_rm1(2,atomic_number) = 0.07085831d0
     FN2_rm1(2,atomic_number) = 9.00001562d0
     FN3_rm1(2,atomic_number) = 1.43802381d0
     FN1_rm1(3,atomic_number) = 0.0d0
     FN2_rm1(3,atomic_number) = 0.0d0
     FN3_rm1(3,atomic_number) = 0.0d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 2
     GSS_rm1(atomic_number) = 16.72091319d0
     GSP_rm1(atomic_number) = 16.76142629d0
     GPP_rm1(atomic_number) = 15.22581028d0
     GP2_rm1(atomic_number) = 14.86578679d0
     HSP_rm1(atomic_number) = 1.99766171d0
     alp_rm1(atomic_number) = 6.00000062d0
     USS_rm1(atomic_number) = -134.18369591d0
     UPP_rm1(atomic_number) = -107.84660920d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 6.043849D0
     p_orb_exp_pm6(atomic_number) = 2.906722D0
     betas_pm6(atomic_number) = -69.922593D0
     betap_pm6(atomic_number) = -30.448165D0
     FN1_pm6(1,atomic_number) = -0.010792D0
     FN2_pm6(1,atomic_number) =  6.004648D0
     FN3_pm6(1,atomic_number) =  1.847724D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 12.446818D0
     GSP_pm6(atomic_number) = 18.496082D0
     GPP_pm6(atomic_number) =  8.417366D0
     GP2_pm6(atomic_number) = 12.179816D0
     HSP_pm6(atomic_number) =  2.604382D0
     USS_pm6(atomic_number) = -140.225626D0
     UPP_pm6(atomic_number) =  -98.778044D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 3.3589210D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END FLUORINE
 !-------------------

 !-------------------
 !NEON
 !-------------------
     atomic_number = 10
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 6.000148d0
     p_orb_exp_pm6(atomic_number) = 3.834528d0
     betas_pm6(atomic_number) = -69.793475d0
     betap_pm6(atomic_number) = -33.261962d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  19.999574d0
     GSP_pm6(atomic_number) =  16.896951d0
     GPP_pm6(atomic_number) =   8.963560d0
     GP2_pm6(atomic_number) =  16.027799d0
     HSP_pm6(atomic_number) =   1.779280d0
     USS_pm6(atomic_number) =  -2.978729d0
     UPP_pm6(atomic_number) = -85.441118d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 4.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END NEON
 !-------------------

 !-------------------
 !SODIUM
 !-------------------
     atomic_number = 11
     
  !MNDO--copied from MNDOD for test--should be re-entered
  ! TL Work
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndo_ref_index(atomic_number) = 101
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 0.98750834D0
     p_orb_exp_mndo(atomic_number) = 0.89334983D0
     betas_mndo(atomic_number) = -1.08738166D0
     betap_mndo(atomic_number) = -0.48623935D0
     GSS_mndo(atomic_number) = 4.59444476D00
     GSP_mndo(atomic_number) = 4.14757400D00
     GPP_mndo(atomic_number) = 4.29919761D00
     GP2_mndo(atomic_number) = 3.79695732D00
     HSP_mndo(atomic_number) = 0.53440874D00
     alp_mndo(atomic_number) = 1.17010202D0
     USS_mndo(atomic_number) = -5.20100000D0
     UPP_mndo(atomic_number) = -2.71257317D0

         
   !MNDOD 
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 0.98750834D0
     p_orb_exp_mndod(atomic_number) = 0.89334983D0
     s_orb_exp_tail_mndod(atomic_number) = 0.65411258D0
     p_orb_exp_tail_mndod(atomic_number) = 0.56440874D0
     betas_mndod(atomic_number) = -1.08738166D0
     betap_mndod(atomic_number) = -0.48623935D0
     rho_core_mndod(atomic_number) = 1.53055325
     GSS_mndod(atomic_number) = 4.59444476D00
     GSP_mndod(atomic_number) = 4.14757400D00
     GPP_mndod(atomic_number) = 4.29919761D00
     GP2_mndod(atomic_number) = 3.79695732D00
     HSP_mndod(atomic_number) = 0.53440874D00
     alp_mndod(atomic_number) = 1.17010202D0
     USS_mndod(atomic_number) = -5.20100000D0
     UPP_mndod(atomic_number) = -2.71257317D0
  
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 0.686327D0
     p_orb_exp_pm6(atomic_number) = 0.950068D0
     betas_pm6(atomic_number) = 0.244853D0
     betap_pm6(atomic_number) = 0.491998D0
     FN1_pm6(1,atomic_number) = -1.026036D0
     FN2_pm6(1,atomic_number) =  2.014506D0
     FN3_pm6(1,atomic_number) =  1.271202D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  4.059972D0
     GSP_pm6(atomic_number) =  7.061183D0
     GPP_pm6(atomic_number) =  9.283540D0
     GP2_pm6(atomic_number) = 17.034978D0
     HSP_pm6(atomic_number) =  0.640715D0
     USS_pm6(atomic_number) = -4.537153D0
     UPP_pm6(atomic_number) = -2.433015D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.2d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END SODIUM
 !-------------------


 !-------------------
 !MAGNESIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0, GSSC = 1, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 12
  !MNDO
  ! TL Work
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
   ! mndo_ref_index(atomic_number) = 101
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.44890446D0
     p_orb_exp_mndo(atomic_number) = 0.95293002D0
     betas_mndo(atomic_number) = -1.89588355D0
     betap_mndo(atomic_number) = -2.14108943D0
     GSS_mndo(atomic_number) = 7.37513258D0
     GSP_mndo(atomic_number) = 6.88890741D0
     GPP_mndo(atomic_number) = 7.04795383D0
     GP2_mndo(atomic_number) = 6.22459871D0
     HSP_mndo(atomic_number) = 0.72673390D0
     alp_mndo(atomic_number) = 1.62146984D0
     USS_mndo(atomic_number) = -15.09700000D0
     UPP_mndo(atomic_number) = -10.65000000D0
  !MNDOD--sp_atom only for now
  ! TL Work
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.44890446D0
     p_orb_exp_mndod(atomic_number) = 0.95293002D0
     betas_mndod(atomic_number) = -1.89588355D0
     betap_mndod(atomic_number) = -2.14108943D0
     s_orb_exp_tail_mndod(atomic_number) = 1.05000000D0
     p_orb_exp_tail_mndod(atomic_number) = 0.92527190D0    
     rho_core_mndod(atomic_number) = 1.35077620
     GSS_mndod(atomic_number) = 7.37513258D0
     GSP_mndod(atomic_number) = 6.88890741D0
     GPP_mndod(atomic_number) = 7.04795383D0
     GP2_mndod(atomic_number) = 6.22459871D0
     HSP_mndod(atomic_number) = 0.72673390D0
     alp_mndod(atomic_number) = 1.62146984D0
     USS_mndod(atomic_number) = -15.09700000D0
     UPP_mndod(atomic_number) = -10.65000000D0 
     
!AM1
   ! Reference:  M.C. HUTTER, J.R. REIMERS, N.S. HUSH, J.PHYS.CHEM. (1998)      
     am1_ref_index(atomic_number) = 103
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) =  1.22339270D0
     p_orb_exp_am1(atomic_number) = 1.02030798D0
     betas_am1(atomic_number) = -1.25974355D0
     betap_am1(atomic_number) = -0.77836604D0
     FN1_am1(1,atomic_number) = 2.55017735D0
     FN2_am1(1,atomic_number) = 4.29397225D0
     FN3_am1(1,atomic_number) = 0.79989601D0
     FN1_am1(2,atomic_number) =-0.00565806D0
     FN2_am1(2,atomic_number) = 2.96053910D0
     FN3_am1(2,atomic_number) = 1.47499983D0
     FN1_am1(3,atomic_number) = -0.00610286D0
     FN2_am1(3,atomic_number) = 2.61416919D0
     FN3_am1(3,atomic_number) = 2.42604040D0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 7.50132277D0
     GSP_am1(atomic_number) = 6.34591536D0
     GPP_am1(atomic_number) = 4.77534467D0
     GP2_am1(atomic_number) = 4.34017279D0
     HSP_am1(atomic_number) = 0.48930466D0
     alp_am1(atomic_number) = 1.67049799D0
     USS_am1(atomic_number) =  -14.96959313D0
     UPP_am1(atomic_number) = -11.56229248D0

  !AM1D--
   ! Reference: P. Imhof et al, J. Chem.  Theo.  Comp., 2, 1050-1056 (2006) (index = 103
     am1d_ref_index(atomic_number) = 104
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.16850D0
     p_orb_exp_am1d(atomic_number) = 1.07072D0
     d_orb_exp_am1d(atomic_number) = 0.93469D0      
     betas_am1d(atomic_number) = -3.60785D0
     betap_am1d(atomic_number) = -2.07794D0
     betad_am1d(atomic_number) = -3.30858D0 
     s_orb_exp_tail_am1d(atomic_number) = 1.61862D0
     p_orb_exp_tail_am1d(atomic_number) = 1.48840D0
     d_orb_exp_tail_am1d(atomic_number) = 1.07347D0  
     GSS_am1d(atomic_number) = 14.645747D00 ! Need to be fixed
     GSP_am1d(atomic_number) = 7.48305D00 
     GPP_am1d(atomic_number) = 11.694918D00 ! Need to be fixed
     GDD_am1d(atomic_number) = 11.694918D00   ! Need to be fixed
     GP2_am1d(atomic_number) = 10.328696D00 !Need to be fixed
     HSP_am1d(atomic_number) = 0.67433D00
     alp_am1d(atomic_number) = 1.28263D0
     USS_am1d(atomic_number) = -16.63758D0
     UPP_am1d(atomic_number) = -11.97469D0
     UDD_am1d(atomic_number) = -10.90361D0
     FN1_am1d(1,atomic_number) = 1.84869D0
     FN2_am1d(1,atomic_number) = 4.22931D0
     FN3_am1d(1,atomic_number) = 0.66917D0
     FN1_am1d(2,atomic_number) = 0.03381D0
     FN2_am1d(2,atomic_number) = 3.57399D0
     FN3_am1d(2,atomic_number) = 2.33163D0
     FN1_am1d(3,atomic_number) = 0.02860D0
     FN2_am1d(3,atomic_number) = 2.27472D0
     FN3_am1d(3,atomic_number) = 2.89337D0
     NUM_FN_am1d(atomic_number) = 3
     rho_core_am1d(atomic_number) = 0.94048D0

  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 0.6985520D0
     p_orb_exp_pm3(atomic_number) = 1.4834530D0
     betas_pm3(atomic_number) = -2.0716910D0
     betap_pm3(atomic_number) = -0.5695810D0
     FN1_pm3(1,atomic_number) = 2.1170500D0
     FN2_pm3(1,atomic_number) = 6.0094770D0
     FN3_pm3(1,atomic_number) = 2.0844060D0
     FN1_pm3(2,atomic_number) = -2.5477670D0
     FN2_pm3(2,atomic_number) = 4.3953700D0
     FN3_pm3(2,atomic_number) = 2.0636740D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 6.6943000D0
     GSP_pm3(atomic_number) = 6.7939950D0
     GPP_pm3(atomic_number) = 6.9104460D0
     GP2_pm3(atomic_number) = 7.0908230D0
     HSP_pm3(atomic_number) = 0.5433000D0
     alp_pm3(atomic_number) = 1.3291470D0
     USS_pm3(atomic_number) = -14.6236880D0
     UPP_pm3(atomic_number) = -14.1734600D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.310830D0
     p_orb_exp_pm6(atomic_number) = 1.388897D0
     betas_pm6(atomic_number) = -9.604932D0
     betap_pm6(atomic_number) =  3.416908D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 7.115328D0
     GSP_pm6(atomic_number) = 3.253024D0
     GPP_pm6(atomic_number) = 4.737311D0
     GP2_pm6(atomic_number) = 8.428485D0
     HSP_pm6(atomic_number) = 0.877379D0
     USS_pm6(atomic_number) = -14.574226D0
     UPP_pm6(atomic_number) =  -7.583850D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.3291470D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END MAGNESIUM
 !-------------------

 !-------------------
 !ALUMINIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 1, GSSC = 1, GPPC = 0, GSPC = 2, GP2C = 0, HSP = -1
   atomic_number = 13
  !MNDO
   ! Reference: L.P. DAVIS, ET.AL.  J. COMP. CHEM., 2, 433, (1981) (index = 5)
     mndo_ref_index(atomic_number) = 5
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.4441610D0
     p_orb_exp_mndo(atomic_number) = 1.4441610D0
     betas_mndo(atomic_number) = -2.6702840D0
     betap_mndo(atomic_number) = -2.6702840D0
     GSS_mndo(atomic_number) = 8.09D00
     GSP_mndo(atomic_number) = 6.63D00
     GPP_mndo(atomic_number) = 5.98D00
     GP2_mndo(atomic_number) = 5.40D00
     HSP_mndo(atomic_number) = 0.70D00
     alp_mndo(atomic_number) = 1.8688394D0
     USS_mndo(atomic_number) = -23.8070970D0
     UPP_mndo(atomic_number) = -17.5198780D0

  !AM1
   ! Reference: M. J. S. Dewar, A. J. Holder, Organometallics, 9, 508-511 (1990). (index = 19)
     am1_ref_index(atomic_number) = 19
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.5165930D0
     p_orb_exp_am1(atomic_number) = 1.3063470D0
     betas_am1(atomic_number) = -3.8668220D0
     betap_am1(atomic_number) = -2.3171460D0
     FN1_am1(1,atomic_number) = 0.0900000D0
     FN2_am1(1,atomic_number) = 12.3924430D0
     FN3_am1(1,atomic_number) = 2.0503940D0
     FN1_am1(2,atomic_number) = 0.0000000D0
     FN2_am1(2,atomic_number) = 0.0000000D0
     FN3_am1(2,atomic_number) = 0.0000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 1
     GSS_am1(atomic_number) = 8.0900000D0
     GSP_am1(atomic_number) = 6.6300000D0
     GPP_am1(atomic_number) = 5.9800000D0
     GP2_am1(atomic_number) = 5.4000000D0
     HSP_am1(atomic_number) = 0.7000000D0
     alp_am1(atomic_number) = 1.9765860D0
     USS_am1(atomic_number) = -24.3535850D0
     UPP_am1(atomic_number) = -18.3636450D0

  !AM1D
   ! Reference: M. J. S. Dewar, A. J. Holder, Organometallics, 9, 508-511 (1990). (index = 19)
     am1d_ref_index(atomic_number) = 19
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.5165930D0
     p_orb_exp_am1d(atomic_number) = 1.3063470D0
     betas_am1d(atomic_number) = -3.8668220D0
     betap_am1d(atomic_number) = -2.3171460D0
     FN1_am1d(1,atomic_number) = 0.0900000D0
     FN2_am1d(1,atomic_number) = 12.3924430D0
     FN3_am1d(1,atomic_number) = 2.0503940D0
     FN1_am1d(2,atomic_number) = 0.0000000D0
     FN2_am1d(2,atomic_number) = 0.0000000D0
     FN3_am1d(2,atomic_number) = 0.0000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 1
     GSS_am1d(atomic_number) = 8.0900000D0
     GSP_am1d(atomic_number) = 6.6300000D0
     GPP_am1d(atomic_number) = 5.9800000D0
     GP2_am1d(atomic_number) = 5.4000000D0
     HSP_am1d(atomic_number) = 0.7000000D0
     alp_am1d(atomic_number) = 1.9765860D0
     USS_am1d(atomic_number) = -24.3535850D0
     UPP_am1d(atomic_number) = -18.3636450D0

  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.7028880D0
     p_orb_exp_pm3(atomic_number) = 1.0736290D0
     betas_pm3(atomic_number) = -0.5943010D0
     betap_pm3(atomic_number) = -0.9565500D0
     FN1_pm3(1,atomic_number) = -0.4730900D0
     FN2_pm3(1,atomic_number) = 1.9158250D0
     FN3_pm3(1,atomic_number) = 1.4517280D0
     FN1_pm3(2,atomic_number) = -0.1540510D0
     FN2_pm3(2,atomic_number) = 6.0050860D0
     FN3_pm3(2,atomic_number) = 2.5199970D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 5.7767370D0
     GSP_pm3(atomic_number) = 11.6598560D0
     GPP_pm3(atomic_number) = 6.3477900D0
     GP2_pm3(atomic_number) = 6.1210770D0
     HSP_pm3(atomic_number) = 4.0062450D0
     alp_pm3(atomic_number) = 1.5217030D0
     USS_pm3(atomic_number) = -24.8454040D0
     UPP_pm3(atomic_number) = -22.2641590D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.364264D0
     p_orb_exp_pm6(atomic_number) = 1.749102D0
     d_orb_exp_pm6(atomic_number) = 1.269384D0
     betas_pm6(atomic_number) = -18.375229D0
     betap_pm6(atomic_number) = -9.382700D0
     betad_pm6(atomic_number) = -20.840474D0
     s_orb_exp_tail_pm6(atomic_number) = 4.742341D0
     p_orb_exp_tail_pm6(atomic_number) = 4.669626D0
     d_orb_exp_tail_pm6(atomic_number) = 7.131138D0
     FN1_pm6(1,atomic_number) =  1.002222D0
     FN2_pm6(1,atomic_number) =  1.517400D0
     FN3_pm6(1,atomic_number) =  0.659101D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  6.652155D0
     GSP_pm6(atomic_number) =  7.459435D0
     GPP_pm6(atomic_number) =  7.668857D0
     GP2_pm6(atomic_number) =  6.673299D0
     HSP_pm6(atomic_number) =  0.435060D0
     USS_pm6(atomic_number) = -24.546778D0
     UPP_pm6(atomic_number) = -20.104434D0
     UDD_pm6(atomic_number) =   8.004394D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.5217030D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END ALUMINIUM
 !-------------------

 !-------------------
 !SILICON
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 2, GSSC = 1, GPPC = -0.5, GSPC = 4, GP2C = 1.5, HSP = -2
   atomic_number = 14
  !MNDO
   ! Reference: M.J.S.DEWAR, ET. AL. ORGANOMETALLICS  5, 375 (1986) (index = 6)
     mndo_ref_index(atomic_number) = 6
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.3159860D0
     p_orb_exp_mndo(atomic_number) = 1.7099430D0
     betas_mndo(atomic_number) = -9.0868040D0
     betap_mndo(atomic_number) = -1.0758270D0
     GSS_mndo(atomic_number) = 9.82D00
     GSP_mndo(atomic_number) = 8.36D00
     GPP_mndo(atomic_number) = 7.31D00
     GP2_mndo(atomic_number) = 6.54D00
     HSP_mndo(atomic_number) = 1.32D00
     alp_mndo(atomic_number) = 2.2053160D0
     USS_mndo(atomic_number) = -37.0375330D0
     UPP_mndo(atomic_number) = -27.7696780D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =   1.91565460D0
     p_orb_exp_mndod(atomic_number)      =   1.68161130D0
     d_orb_exp_mndod(atomic_number)      =   0.96677166D0   
     betas_mndod(atomic_number)          =  -8.21073420D0
     betap_mndod(atomic_number)          =  -4.88462030D0
     betad_mndod(atomic_number)          =  -2.60801150D0  
     s_orb_exp_tail_mndod(atomic_number) =   1.52929180D0
     p_orb_exp_tail_mndod(atomic_number) =   0.97628075D0
     d_orb_exp_tail_mndod(atomic_number) =   0.93816441D0
     GSS_mndod(atomic_number)            =  10.74164700D0
     GSP_mndod(atomic_number)            =   7.56066640D0
     GPP_mndod(atomic_number)            =   7.43649690D0
     GDD_mndod(atomic_number)            =   7.05875000D0                                        
     GP2_mndod(atomic_number)            =   6.56775150D0
     HSP_mndod(atomic_number)            =   0.87753880D0
     alp_mndod(atomic_number)            =   1.66006930D0
     USS_mndod(atomic_number)            = -36.05153000D0
     UPP_mndod(atomic_number)            = -27.53569100D0
     UDD_mndod(atomic_number)            = -14.67743900D0

  !AM1
   ! Reference: M.J.S.DEWAR, C. JIE, ORGANOMETALLICS, 6, 1486-1490 (1987). (index = 20)
     am1_ref_index(atomic_number) = 20
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.830697D0
     p_orb_exp_am1(atomic_number) = 1.2849530D0
     betas_am1(atomic_number) = -3.784852D0
     betap_am1(atomic_number) = -1.968123D0
     FN1_am1(1,atomic_number) = 0.25D0
     FN2_am1(1,atomic_number) = 9.000D0
     FN3_am1(1,atomic_number) = 0.911453D0
     FN1_am1(2,atomic_number) = 0.061513D0
     FN2_am1(2,atomic_number) = 5.00D0
     FN3_am1(2,atomic_number) = 1.995569D0
     FN1_am1(3,atomic_number) = 0.0207890D0
     FN2_am1(3,atomic_number) = 5.00D0
     FN3_am1(3,atomic_number) = 2.990610D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 9.8200000D0
     GSP_am1(atomic_number) = 8.3600000D0
     GPP_am1(atomic_number) = 7.3100000D0
     GP2_am1(atomic_number) = 6.5400000D0
     HSP_am1(atomic_number) = 1.3200000D0
     alp_am1(atomic_number) = 2.257816D0
     USS_am1(atomic_number) = -33.9536220D0
     UPP_am1(atomic_number) = -28.9347490D0

  !AM1D
   ! Reference: M.J.S.DEWAR, C. JIE, ORGANOMETALLICS, 6, 1486-1490 (1987). (index = 20)
     am1d_ref_index(atomic_number) = 20
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.830697D0
     p_orb_exp_am1d(atomic_number) = 1.2849530D0
     betas_am1d(atomic_number) = -3.784852D0
     betap_am1d(atomic_number) = -1.968123D0
     FN1_am1d(1,atomic_number) = 0.25D0
     FN2_am1d(1,atomic_number) = 9.000D0
     FN3_am1d(1,atomic_number) = 0.911453D0
     FN1_am1d(2,atomic_number) = 0.061513D0
     FN2_am1d(2,atomic_number) = 5.00D0
     FN3_am1d(2,atomic_number) = 1.995569D0
     FN1_am1d(3,atomic_number) = 0.0207890D0
     FN2_am1d(3,atomic_number) = 5.00D0
     FN3_am1d(3,atomic_number) = 2.990610D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 3
     GSS_am1d(atomic_number) = 9.8200000D0
     GSP_am1d(atomic_number) = 8.3600000D0
     GPP_am1d(atomic_number) = 7.3100000D0
     GP2_am1d(atomic_number) = 6.5400000D0
     HSP_am1d(atomic_number) = 1.3200000D0
     alp_am1d(atomic_number) = 2.257816D0
     USS_am1d(atomic_number) = -33.9536220D0
     UPP_am1d(atomic_number) = -28.9347490D0

  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.6350750D0
     p_orb_exp_pm3(atomic_number) = 1.3130880D0
     betas_pm3(atomic_number) = -2.8621450D0
     betap_pm3(atomic_number) = -3.9331480D0
     FN1_pm3(1,atomic_number) = -0.3906000D0
     FN2_pm3(1,atomic_number) = 6.0000540D0
     FN3_pm3(1,atomic_number) = 0.6322620D0
     FN1_pm3(2,atomic_number) = 0.0572590D0
     FN2_pm3(2,atomic_number) = 6.0071830D0
     FN3_pm3(2,atomic_number) = 2.0199870D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 5.0471960D0
     GSP_pm3(atomic_number) = 5.9490570D0
     GPP_pm3(atomic_number) = 6.7593670D0
     GP2_pm3(atomic_number) = 5.1612970D0
     HSP_pm3(atomic_number) = 0.9198320D0
     alp_pm3(atomic_number) = 2.1358090D0
     USS_pm3(atomic_number) = -26.7634830D0
     UPP_pm3(atomic_number) = -22.8136350D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Guimaraes and Jorgensen, 2005, J Chem Theory Comput., 1: 817-823
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -66.839000D0
     s_orb_exp_pddgpm3(atomic_number) = 1.586389D0
     p_orb_exp_pddgpm3(atomic_number) = 1.485958D0
     betas_pddgpm3(atomic_number) = -3.376445D0
     betap_pddgpm3(atomic_number) = -3.620969D0
     FN1_pddgpm3(1,atomic_number) = -0.071314D0
     FN2_pddgpm3(1,atomic_number) = 6.000000D0
     FN3_pddgpm3(1,atomic_number) = 0.237995D0
     FN1_pddgpm3(2,atomic_number) = 0.089451D0
     FN2_pddgpm3(2,atomic_number) = 6.000000D0
     FN3_pddgpm3(2,atomic_number) = 1.897728D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 5.0471960D0
     GSP_pddgpm3(atomic_number) = 5.9490570D0
     GPP_pddgpm3(atomic_number) = 6.7593670D0
     GP2_pddgpm3(atomic_number) = 5.1612970D0
     HSP_pddgpm3(atomic_number) = 0.9198320D0
     alp_pddgpm3(atomic_number) = 2.215157D0
     USS_pddgpm3(atomic_number) = -26.332522D0
     UPP_pddgpm3(atomic_number) = -22.602540D0
     PDDGC1_pm3(atomic_number) = -0.091928D0
     PDDGC2_pm3(atomic_number) = -0.040753D0
     PDDGE1_pm3(atomic_number) =  1.163190D0
     PDDGE2_pm3(atomic_number) = 2.190526D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.752741D0
     p_orb_exp_pm6(atomic_number) = 1.198413D0
     d_orb_exp_pm6(atomic_number) = 2.128593D0
     betas_pm6(atomic_number) = -8.686909D0
     betap_pm6(atomic_number) = -1.856482D0
     betad_pm6(atomic_number) = -6.360627D0
     s_orb_exp_tail_pm6(atomic_number) = 8.388111D0
     p_orb_exp_tail_pm6(atomic_number) = 1.843048D0
     d_orb_exp_tail_pm6(atomic_number) = 0.708600D0
     FN1_pm6(1,atomic_number) =  0.208571D0
     FN2_pm6(1,atomic_number) =  6.000483D0
     FN3_pm6(1,atomic_number) =  1.185245D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  5.194805D0
     GSP_pm6(atomic_number) =  5.090534D0
     GPP_pm6(atomic_number) =  5.185150D0
     GP2_pm6(atomic_number) =  4.769775D0
     HSP_pm6(atomic_number) =  1.425012D0
     USS_pm6(atomic_number) = -27.358058D0
     UPP_pm6(atomic_number) = -20.490578D0
     UDD_pm6(atomic_number) = -22.751900D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 2.1358090D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END SILICON
 !-------------------

 !-------------------
 !PHOSPHORUS
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 3, GSSC = 1, GPPC = -1.5, GSPC = 6, GP2C = 4.5, HSP = -3
   atomic_number = 15
  !MNDO
   ! Reference: M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA,J. AM. CHEM. SOC., 100 3607 1978 (index = 7)
     mndo_ref_index(atomic_number) = 7
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.1087200D0
     p_orb_exp_mndo(atomic_number) = 1.7858100D0
     betas_mndo(atomic_number) = -6.7916000D0
     betap_mndo(atomic_number) = -6.7916000D0
     GSS_mndo(atomic_number) = 11.56D00
     GSP_mndo(atomic_number) = 10.08D00
     GPP_mndo(atomic_number) = 8.64D00
     GP2_mndo(atomic_number) = 7.68D00
     HSP_mndo(atomic_number) = 1.92D00
     alp_mndo(atomic_number) = 2.4152800D0
     USS_mndo(atomic_number) = -56.1433600D0
     UPP_mndo(atomic_number) = -42.8510800D0

  !MNDOD
  ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.266463D0
     p_orb_exp_mndod(atomic_number) = 1.940015D0
     d_orb_exp_mndod(atomic_number) = 1.100109D0     
     betas_mndod(atomic_number) = -8.902104D0
     betap_mndod(atomic_number) = -9.386110D0
     betad_mndod(atomic_number) = -2.091701D0 
     s_orb_exp_tail_mndod(atomic_number) = 1.63437610D0
     p_orb_exp_tail_mndod(atomic_number) = 1.08291170D0
     d_orb_exp_tail_mndod(atomic_number) = 1.00651470D0  
     GSS_mndod(atomic_number) = 11.479753D00
     GSP_mndod(atomic_number) = 8.55756910D00
     GPP_mndod(atomic_number) = 8.2487228D00
     GDD_mndod(atomic_number) = 7.573017D00
     GP2_mndod(atomic_number) = 7.2850917D00
     HSP_mndod(atomic_number) = 2.1078044D00
     alp_mndod(atomic_number) = 1.852551D0
     USS_mndod(atomic_number) = -47.055529D0
     UPP_mndod(atomic_number) = -38.067059D0
     UDD_mndod(atomic_number) = -23.691597D0
  !AM1
   ! Reference: M.J.S.DEWAR, JIE, C, THEOCHEM, 187,1 (1989) (index = 21)
     am1_ref_index(atomic_number) = 21
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.9812800D0
     p_orb_exp_am1(atomic_number) = 1.8751500D0
     betas_am1(atomic_number) = -6.3537640D0
     betap_am1(atomic_number) = -6.5907090D0
     FN1_am1(1,atomic_number) = -0.0318270D0
     FN2_am1(1,atomic_number) = 6.0000000D0
     FN3_am1(1,atomic_number) = 1.4743230D0
     FN1_am1(2,atomic_number) = 0.0184700D0
     FN2_am1(2,atomic_number) = 7.0000000D0
     FN3_am1(2,atomic_number) = 1.7793540D0
     FN1_am1(3,atomic_number) = 0.0332900D0
     FN2_am1(3,atomic_number) = 9.0000000D0
     FN3_am1(3,atomic_number) = 3.0065760D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 11.5600050D0
     GSP_am1(atomic_number) = 5.2374490D0
     GPP_am1(atomic_number) = 7.8775890D0
     GP2_am1(atomic_number) = 7.3076480D0
     HSP_am1(atomic_number) = 0.7792380D0
     alp_am1(atomic_number) = 2.4553220D0
     USS_am1(atomic_number) = -42.0298630D0
     UPP_am1(atomic_number) = -34.0307090D0

  !AM1D--
   ! Reference: K. Nam, Q. Cui, J. Gao, D. York. J. CHEM. THEO. COMP., 3, 486, (2007) (102) AM1/d-PhoT
     am1d_ref_index(atomic_number) = 102
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.909168D0
     p_orb_exp_am1d(atomic_number) = 2.008466D0
     d_orb_exp_am1d(atomic_number) = 0.840667D0      
     betas_am1d(atomic_number) = -11.194791D0
     betap_am1d(atomic_number) = -11.985621D0
     betad_am1d(atomic_number) = -2.360095D0 
     s_orb_exp_tail_am1d(atomic_number) = 2.085120D0
     p_orb_exp_tail_am1d(atomic_number) = 1.535336D0
     d_orb_exp_tail_am1d(atomic_number) = 1.236266D0  
     GSS_am1d(atomic_number) = 14.645747D00
     GSP_am1d(atomic_number) = 5.689654D00
     GPP_am1d(atomic_number) = 11.694918D00
     GDD_am1d(atomic_number) = 11.694918D00  ! GDD is not used
     GP2_am1d(atomic_number) = 10.328696D00
     HSP_am1d(atomic_number) = 1.175115D00
     alp_am1d(atomic_number) = 1.883237D0
     USS_am1d(atomic_number) = -46.250810D0
     UPP_am1d(atomic_number) = -40.712918D0
     UDD_am1d(atomic_number) = -24.504161D0
     GNN_am1d(atomic_number) = 0.353722D0
     rho_core_am1d(atomic_number) = 1.185130D0
     FN1_am1d(1,atomic_number) = -0.344529D0
     FN2_am1d(1,atomic_number) = 3.034933D0
     FN3_am1d(1,atomic_number) = 1.134275D0
     FN1_am1d(2,atomic_number) = -0.021847D0
     FN2_am1d(2,atomic_number) = 1.684515D0
     FN3_am1d(2,atomic_number) = 2.716684D0
     FN1_am1d(3,atomic_number) = -0.036003D0
     FN2_am1d(3,atomic_number) = 5.243357D0
     FN3_am1d(3,atomic_number) = 1.924175D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 3

  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.0175630D0
     p_orb_exp_pm3(atomic_number) = 1.5047320D0
     betas_pm3(atomic_number) = -12.6158790D0
     betap_pm3(atomic_number) = -4.1600400D0
     FN1_pm3(1,atomic_number) = -0.6114210D0
     FN2_pm3(1,atomic_number) = 1.9972720D0
     FN3_pm3(1,atomic_number) = 0.7946240D0
     FN1_pm3(2,atomic_number) = -0.0939350D0
     FN2_pm3(2,atomic_number) = 1.9983600D0
     FN3_pm3(2,atomic_number) = 1.9106770D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 7.8016150D0
     GSP_pm3(atomic_number) = 5.1869490D0
     GPP_pm3(atomic_number) = 6.6184780D0
     GP2_pm3(atomic_number) = 6.0620020D0
     HSP_pm3(atomic_number) = 1.5428090D0
     alp_pm3(atomic_number) = 1.9405340D0
     USS_pm3(atomic_number) = -40.4130960D0
     UPP_pm3(atomic_number) = -29.5930520D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Guimaraes and Jorgensen, 2005, J Chem Theory Comput., 1: 817-823 (Index 31)
     pddgpm3_ref_index(atomic_number) = 31
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -117.212854D0
     s_orb_exp_pddgpm3(atomic_number) = 2.395882D0
     p_orb_exp_pddgpm3(atomic_number) = 1.742213D0
     betas_pddgpm3(atomic_number) = -12.676297D0
     betap_pddgpm3(atomic_number) = -7.093318D0
     FN1_pddgpm3(1,atomic_number) = -0.398055D0
     FN2_pddgpm3(1,atomic_number) = 1.997272D0
     FN3_pddgpm3(1,atomic_number) = 0.950073D0
     FN1_pddgpm3(2,atomic_number) = -0.079653D0
     FN2_pddgpm3(2,atomic_number) = 1.998360D0
     FN3_pddgpm3(2,atomic_number) = 2.336959D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 7.8016150D0
     GSP_pddgpm3(atomic_number) = 5.1869490D0
     GPP_pddgpm3(atomic_number) = 6.6184780D0
     GP2_pddgpm3(atomic_number) = 6.0620020D0
     HSP_pddgpm3(atomic_number) = 1.5428090D0
     alp_pddgpm3(atomic_number) = 2.005294D0
     USS_pddgpm3(atomic_number) = -37.882113D0
     UPP_pddgpm3(atomic_number) = -30.312979D0
     PDDGC1_pm3(atomic_number) = 0.462741D0
     PDDGC2_pm3(atomic_number) = -0.020444D0
     PDDGE1_pm3(atomic_number) = 0.714296D0
     PDDGE2_pm3(atomic_number) = 2.041209D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 2.12240118d0
     p_orb_exp_rm1(atomic_number) = 1.74327954d0
     betas_rm1(atomic_number) = -6.13514969d0
     betap_rm1(atomic_number) = -5.94442127d0
     FN1_rm1(1,atomic_number) = -0.41063467d0
     FN2_rm1(1,atomic_number) = 6.08752832d0
     FN3_rm1(1,atomic_number) = 1.31650261d0
     FN1_rm1(2,atomic_number) = -0.16299288d0
     FN2_rm1(2,atomic_number) = 7.09472602d0
     FN3_rm1(2,atomic_number) = 1.90721319d0
     FN1_rm1(3,atomic_number) = -0.04887125d0
     FN2_rm1(3,atomic_number) = 8.99979308d0
     FN3_rm1(3,atomic_number) = 2.658577780d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 3
     GSS_rm1(atomic_number) = 11.08059265d0
     GSP_rm1(atomic_number) = 5.68339201d0
     GPP_rm1(atomic_number) = 7.60417563d0
     GP2_rm1(atomic_number) = 7.40265182d0
     HSP_rm1(atomic_number) = 1.16181792d0
     alp_rm1(atomic_number) = 1.90993294d0
     USS_rm1(atomic_number) = -41.81533184d0
     UPP_rm1(atomic_number) = -34.38342529d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.158033D0
     p_orb_exp_pm6(atomic_number) = 1.805343D0
     d_orb_exp_pm6(atomic_number) = 1.230358D0
     betas_pm6(atomic_number) = -14.583780D0
     betap_pm6(atomic_number) = -11.744725D0
     betad_pm6(atomic_number) = -20.099893D0
     s_orb_exp_tail_pm6(atomic_number) = 6.042706D0
     p_orb_exp_tail_pm6(atomic_number) = 2.376473D0
     d_orb_exp_tail_pm6(atomic_number) = 7.147750D0
     FN1_pm6(1,atomic_number) = -0.034320D0
     FN2_pm6(1,atomic_number) =  6.001394D0
     FN3_pm6(1,atomic_number) =  2.296737D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  8.758856D0
     GSP_pm6(atomic_number) =  8.483679D0
     GPP_pm6(atomic_number) =  8.662754D0
     GP2_pm6(atomic_number) =  7.734264D0
     HSP_pm6(atomic_number) =  0.871681D0
     USS_pm6(atomic_number) = -48.729905D0
     UPP_pm6(atomic_number) = -40.354689D0
     UDD_pm6(atomic_number) = -7.349246D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.9405340D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END PHOSPHORUS
 !-------------------

 !-------------------
 !SULPHUR
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 4, GSSC = 1, GPPC = -0.5, GSPC = 8, GP2C = 6.5, HSP = -4
   atomic_number = 16
  !MNDO
   ! Reference: M.J.S.DEWAR, C.H. REYNOLDS, J. COMP. CHEM. 7, 140-143 (1986) (index = 8)
     mndo_ref_index(atomic_number) = 8
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.3129620D0
     p_orb_exp_mndo(atomic_number) = 2.0091460D0
     betas_mndo(atomic_number) = -10.7616700D0
     betap_mndo(atomic_number) = -10.1084330D0
     GSS_mndo(atomic_number) = 12.88D00
     GSP_mndo(atomic_number) = 11.26D00
     GPP_mndo(atomic_number) = 9.90D00
     GP2_mndo(atomic_number) = 8.83D00
     HSP_mndo(atomic_number) = 2.26D00
     alp_mndo(atomic_number) = 2.4780260D0
     USS_mndo(atomic_number) = -72.2422810D0
     UPP_mndo(atomic_number) = -56.9732070D0
  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.225851D0
     p_orb_exp_mndod(atomic_number) = 2.099706D0
     d_orb_exp_mndod(atomic_number) = 1.231472D0     
     betas_mndod(atomic_number) = -10.999545D0
     betap_mndod(atomic_number) = -12.215437D0
     betad_mndod(atomic_number) = -1.880669D0   
     s_orb_exp_tail_mndod(atomic_number) = 1.7363910D0
     p_orb_exp_tail_mndod(atomic_number) = 1.1211820D0
     d_orb_exp_tail_mndod(atomic_number) = 1.0508470D0 
     GSS_mndod(atomic_number) = 12.196301D00
     GSP_mndod(atomic_number) = 8.853901D00
     GPP_mndod(atomic_number) = 8.540233D00
     GDD_mndod(atomic_number) = 7.906571D00
     GP2_mndod(atomic_number) = 7.542547D00
     HSP_mndod(atomic_number) = 2.646352D00
     alp_mndod(atomic_number) = 2.023060D0
     USS_mndod(atomic_number) = -56.889130D0
     UPP_mndod(atomic_number) = -47.274746D0
     UDD_mndod(atomic_number) = -25.095118D0
  !AM1
   ! Reference: M.J.S.DEWAR, Y-C YUAN, INORGANIC CHEMISTRY 29 (19): 3881-3890 SEP 19 1990 (Index = 22)
     am1_ref_index(atomic_number) = 22
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 2.3665150D0
     p_orb_exp_am1(atomic_number) = 1.6672630D0
     betas_am1(atomic_number) = -3.9205660D0
     betap_am1(atomic_number) = -7.9052780D0
     FN1_am1(1,atomic_number) = -0.5091950D0
     FN2_am1(1,atomic_number) = 4.5936910D0
     FN3_am1(1,atomic_number) = 0.7706650D0
     FN1_am1(2,atomic_number) = -0.0118630D0
     FN2_am1(2,atomic_number) = 5.8657310D0
     FN3_am1(2,atomic_number) = 1.5033130D0
     FN1_am1(3,atomic_number) = 0.0123340D0
     FN2_am1(3,atomic_number) = 13.5573360D0
     FN3_am1(3,atomic_number) = 2.0091730D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 3
     GSS_am1(atomic_number) = 11.7863290D0
     GSP_am1(atomic_number) = 8.6631270D0
     GPP_am1(atomic_number) = 10.0393080D0
     GP2_am1(atomic_number) = 7.7816880D0
     HSP_am1(atomic_number) = 2.5321370D0
     alp_am1(atomic_number) = 2.4616480D0
     USS_am1(atomic_number) = -56.6940560D0
     UPP_am1(atomic_number) = -48.7170490D0

  !AM1D
   ! Reference: M.J.S.DEWAR, Y-C YUAN, INORGANIC CHEMISTRY 29 (19): 3881-3890 SEP 19 1990 (Index = 22)
     am1d_ref_index(atomic_number) = 22
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 2.3665150D0
     p_orb_exp_am1d(atomic_number) = 1.6672630D0
     betas_am1d(atomic_number) = -3.9205660D0
     betap_am1d(atomic_number) = -7.9052780D0
     FN1_am1d(1,atomic_number) = -0.5091950D0
     FN2_am1d(1,atomic_number) = 4.5936910D0
     FN3_am1d(1,atomic_number) = 0.7706650D0
     FN1_am1d(2,atomic_number) = -0.0118630D0
     FN2_am1d(2,atomic_number) = 5.8657310D0
     FN3_am1d(2,atomic_number) = 1.5033130D0
     FN1_am1d(3,atomic_number) = 0.0123340D0
     FN2_am1d(3,atomic_number) = 13.5573360D0
     FN3_am1d(3,atomic_number) = 2.0091730D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 3
     GSS_am1d(atomic_number) = 11.7863290D0
     GSP_am1d(atomic_number) = 8.6631270D0
     GPP_am1d(atomic_number) = 10.0393080D0
     GP2_am1d(atomic_number) = 7.7816880D0
     HSP_am1d(atomic_number) = 2.5321370D0
     alp_am1d(atomic_number) = 2.4616480D0
     USS_am1d(atomic_number) = -56.6940560D0
     UPP_am1d(atomic_number) = -48.7170490D0


  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.8911850D0
     p_orb_exp_pm3(atomic_number) = 1.6589720D0
     betas_pm3(atomic_number) = -8.8274650D0
     betap_pm3(atomic_number) = -8.0914150D0
     FN1_pm3(1,atomic_number) = -0.3991910D0
     FN2_pm3(1,atomic_number) = 6.0006690D0
     FN3_pm3(1,atomic_number) = 0.9621230D0
     FN1_pm3(2,atomic_number) = -0.0548990D0
     FN2_pm3(2,atomic_number) = 6.0018450D0
     FN3_pm3(2,atomic_number) = 1.5799440D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 8.9646670D0
     GSP_pm3(atomic_number) = 6.7859360D0
     GPP_pm3(atomic_number) = 9.9681640D0
     GP2_pm3(atomic_number) = 7.9702470D0
     HSP_pm3(atomic_number) = 4.0418360D0
     alp_pm3(atomic_number) = 2.2697060D0
     USS_pm3(atomic_number) = -49.8953710D0
     UPP_pm3(atomic_number) = -44.3925830D0

  !PM3/MM* 2nd version
     element_supported_pm3mmx2(atomic_number) = .true.
     rho_pm3mmx2(atomic_number) = 0.268D0
     scale_f1_pm3mmx2(1,atomic_number) = 3.289D0
     scale_f2_pm3mmx2(1,atomic_number) = 2.392D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Guimaraes and Jorgensen, 2005, J Chem Theory Comput., 1: 817-823 (Index 31)
     pddgpm3_ref_index(atomic_number) = 31
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -166.3365540000D0
     s_orb_exp_pddgpm3(atomic_number) = 1.0120020000D0
     p_orb_exp_pddgpm3(atomic_number) = 1.8769990000D0
     betas_pddgpm3(atomic_number) = -2.9539120000D0
     betap_pddgpm3(atomic_number) = -8.5077790000D0
     FN1_pddgpm3(1,atomic_number) = -0.330692000D0
     FN2_pddgpm3(1,atomic_number) = 6.000000D0
     FN3_pddgpm3(1,atomic_number) = 0.823837000D0
     FN1_pddgpm3(2,atomic_number) = 0.024171000D0
     FN2_pddgpm3(2,atomic_number) = 6.000000D0
     FN3_pddgpm3(2,atomic_number) = 2.017756000D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 8.9646670D0
     GSP_pddgpm3(atomic_number) = 6.7859360D0
     GPP_pddgpm3(atomic_number) = 9.9681640D0
     GP2_pddgpm3(atomic_number) = 7.9702470D0
     HSP_pddgpm3(atomic_number) = 4.0418360D0
     alp_pddgpm3(atomic_number) = 2.5397510000D0
     USS_pddgpm3(atomic_number) = -43.9063660000D0
     UPP_pddgpm3(atomic_number) = -43.4613480000D0
     PDDGC1_pm3(atomic_number) = 0.120434000D0
     PDDGC2_pm3(atomic_number) = -0.002663D0
     PDDGE1_pm3(atomic_number) = 0.672870D0
     PDDGE2_pm3(atomic_number) = 2.032340D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 2.13344308d0
     p_orb_exp_rm1(atomic_number) = 1.87460650d0
     betas_rm1(atomic_number) = -1.95910719d0
     betap_rm1(atomic_number) = -8.77430652
     FN1_rm1(1,atomic_number) = -0.74601055d0
     FN2_rm1(1,atomic_number) = 4.81038002d0
     FN3_rm1(1,atomic_number) = 0.59380129d0
     FN1_rm1(2,atomic_number) = -0.06519286d0
     FN2_rm1(2,atomic_number) = 7.20760864d0
     FN3_rm1(2,atomic_number) = 1.29492008d0
     FN1_rm1(3,atomic_number) = -0.00655977d0
     FN2_rm1(3,atomic_number) = 9.00000180d0
     FN3_rm1(3,atomic_number) = 1.80060151d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 3
     GSS_rm1(atomic_number) = 12.48828408d0
     GSP_rm1(atomic_number) = 8.56910574d0
     GPP_rm1(atomic_number) = 8.52301167d0
     GP2_rm1(atomic_number) = 7.66863296d0
     HSP_rm1(atomic_number) = 3.88978932d0
     alp_rm1(atomic_number) = 2.44015636d0
     USS_rm1(atomic_number) = -55.16775121d0
     UPP_rm1(atomic_number) = -46.52930422d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true. 
     s_orb_exp_pm6(atomic_number) = 2.192844D0
     p_orb_exp_pm6(atomic_number) = 1.841078D0
     d_orb_exp_pm6(atomic_number) = 3.109401D0
     betas_pm6(atomic_number) = -13.827440D0
     betap_pm6(atomic_number) = -7.664613D0
     betad_pm6(atomic_number) = -9.986172D0
     s_orb_exp_tail_pm6(atomic_number) = 0.479722D0 
     p_orb_exp_tail_pm6(atomic_number) = 1.015507D0
     d_orb_exp_tail_pm6(atomic_number) = 4.317470D0
     FN1_pm6(1,atomic_number) = -0.036928D0
     FN2_pm6(1,atomic_number) =  1.795067D0
     FN3_pm6(1,atomic_number) =  2.082618D0
     FN1_pm6(2,atomic_number) = 0.0d0 
     FN2_pm6(2,atomic_number) = 0.0d0 
     FN3_pm6(2,atomic_number) = 0.0d0 
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  9.170350D0
     GSP_pm6(atomic_number) =  5.944296D0
     GPP_pm6(atomic_number) =  8.165473D0
     GP2_pm6(atomic_number) =  7.301878D0
     HSP_pm6(atomic_number) =  5.005404D0
     USS_pm6(atomic_number) = -47.530706D0
     UPP_pm6(atomic_number) = -39.191045D0
     UDD_pm6(atomic_number) = -46.306944D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 2.2697060D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END SULPHUR
 !-------------------

 !-------------------
 !CHLORINE
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 5, GSSC = 1, GPPC = 0, GSPC = 10, GP2C = 10, HSP = -5
   atomic_number = 17
  !MNDO
   ! Reference:  M.J.S.DEWAR, H.S.RZEPA, J. COMP. CHEM., 4, 158, (1983) (index = 9)
     mndo_ref_index(atomic_number) = 9
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 3.7846450D0
     p_orb_exp_mndo(atomic_number) = 2.0362630D0
     betas_mndo(atomic_number) = -14.2623200D0
     betap_mndo(atomic_number) = -14.2623200D0
     GSS_mndo(atomic_number) = 15.03D00
     GSP_mndo(atomic_number) = 13.16D00
     GPP_mndo(atomic_number) = 11.30D00
     GP2_mndo(atomic_number) = 9.97D00
     HSP_mndo(atomic_number) = 2.42D00
     alp_mndo(atomic_number) = 2.5422010D0
     USS_mndo(atomic_number) = -100.2271660D0
     UPP_mndo(atomic_number) = -77.3786670D0
  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.561611D0
     p_orb_exp_mndod(atomic_number) = 2.389338D0
     d_orb_exp_mndod(atomic_number) = 1.251398D0     
     betas_mndod(atomic_number) = -6.037292D0
     betap_mndod(atomic_number) = -19.183386D0
     betad_mndod(atomic_number) = -1.877782D0      
     s_orb_exp_tail_mndod(atomic_number) = 1.88087547D0
     p_orb_exp_tail_mndod(atomic_number) = 1.18104227D0
     d_orb_exp_tail_mndod(atomic_number) = 1.14061555D0   
     GSS_mndod(atomic_number) = 13.211148D00
     GSP_mndod(atomic_number) = 9.419496D00
     GPP_mndod(atomic_number) = 8.996201D00
     GDD_mndod(atomic_number) = 8.581992D00
     GP2_mndod(atomic_number) = 7.945248D00
     HSP_mndod(atomic_number) = 3.081499D00
     alp_mndod(atomic_number) = 2.180300D0
     USS_mndod(atomic_number) = -69.622971D0
     UPP_mndod(atomic_number) = -59.100731D0
     UDD_mndod(atomic_number) = -36.674572D0
  !AM1
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1_ref_index(atomic_number) = 18
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 3.6313760D0
     p_orb_exp_am1(atomic_number) = 2.0767990D0
     betas_am1(atomic_number) = -24.5946700D0
     betap_am1(atomic_number) = -14.6372160D0
     FN1_am1(1,atomic_number) = 0.0942430D0
     FN2_am1(1,atomic_number) = 4.0000000D0
     FN3_am1(1,atomic_number) = 1.3000000D0
     FN1_am1(2,atomic_number) = 0.0271680D0
     FN2_am1(2,atomic_number) = 4.0000000D0
     FN3_am1(2,atomic_number) = 2.1000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 2
     GSS_am1(atomic_number) = 15.0300000D0
     GSP_am1(atomic_number) = 13.1600000D0
     GPP_am1(atomic_number) = 11.3000000D0
     GP2_am1(atomic_number) = 9.9700000D0
     HSP_am1(atomic_number) = 2.4200000D0
     alp_am1(atomic_number) = 2.9193680D0
     USS_am1(atomic_number) = -111.6139480D0
     UPP_am1(atomic_number) = -76.6401070D0

  !AM1D
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1d_ref_index(atomic_number) = 18
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 3.6313760D0
     p_orb_exp_am1d(atomic_number) = 2.0767990D0
     betas_am1d(atomic_number) = -24.5946700D0
     betap_am1d(atomic_number) = -14.6372160D0
     FN1_am1d(1,atomic_number) = 0.0942430D0
     FN2_am1d(1,atomic_number) = 4.0000000D0
     FN3_am1d(1,atomic_number) = 1.3000000D0
     FN1_am1d(2,atomic_number) = 0.0271680D0
     FN2_am1d(2,atomic_number) = 4.0000000D0
     FN3_am1d(2,atomic_number) = 2.1000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 2
     GSS_am1d(atomic_number) = 15.0300000D0
     GSP_am1d(atomic_number) = 13.1600000D0
     GPP_am1d(atomic_number) = 11.3000000D0
     GP2_am1d(atomic_number) = 9.9700000D0
     HSP_am1d(atomic_number) = 2.4200000D0
     alp_am1d(atomic_number) = 2.9193680D0
     USS_am1d(atomic_number) = -111.6139480D0
     UPP_am1d(atomic_number) = -76.6401070D0

  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.2462100D0
     p_orb_exp_pm3(atomic_number) = 2.1510100D0
     betas_pm3(atomic_number) = -27.5285600D0
     betap_pm3(atomic_number) = -11.5939220D0
     FN1_pm3(1,atomic_number) = -0.1715910D0
     FN2_pm3(1,atomic_number) = 6.0008020D0
     FN3_pm3(1,atomic_number) = 1.0875020D0
     FN1_pm3(2,atomic_number) = -0.0134580D0
     FN2_pm3(2,atomic_number) = 1.9666180D0
     FN3_pm3(2,atomic_number) = 2.2928910D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 16.0136010D0
     GSP_pm3(atomic_number) = 8.0481150D0
     GPP_pm3(atomic_number) = 7.5222150D0
     GP2_pm3(atomic_number) = 7.5041540D0
     HSP_pm3(atomic_number) = 3.4811530D0
     alp_pm3(atomic_number) = 2.5172960D0
     USS_pm3(atomic_number) = -100.6267470D0
     UPP_pm3(atomic_number) = -53.6143960D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgpm3_ref_index(atomic_number) = 30
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -305.715201D0
     s_orb_exp_pddgpm3(atomic_number) = 2.548268D0
     p_orb_exp_pddgpm3(atomic_number) = 2.284624D0
     betas_pddgpm3(atomic_number) = -26.913129D0
     betap_pddgpm3(atomic_number) = -14.991178D0
     FN1_pddgpm3(1,atomic_number) = -0.112222D0
     FN2_pddgpm3(1,atomic_number) = 5.963719D0
     FN3_pddgpm3(1,atomic_number) = 1.027719D0
     FN1_pddgpm3(2,atomic_number) = -0.013061D0
     FN2_pddgpm3(2,atomic_number) = 1.999556D0
     FN3_pddgpm3(2,atomic_number) = 2.286377D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 16.0136010D0
     GSP_pddgpm3(atomic_number) = 8.0481150D0
     GPP_pddgpm3(atomic_number) = 7.5222150D0
     GP2_pddgpm3(atomic_number) = 7.5041540D0
     HSP_pddgpm3(atomic_number) = 3.4811530D0
     alp_pddgpm3(atomic_number) = 2.497617D0
     USS_pddgpm3(atomic_number) = -95.094434D0
     UPP_pddgpm3(atomic_number) = -53.921651D0
     PDDGC1_pm3(atomic_number) = -0.016552D0
     PDDGC2_pm3(atomic_number) = -0.016646D0
     PDDGE1_pm3(atomic_number) = 1.727690D0
     PDDGE2_pm3(atomic_number) = 1.784655D0

  !PDDG/MNDO
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgmndo_ref_index(atomic_number) = 30
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -378.909727D0
     s_orb_exp_pddgmndo(atomic_number) = 4.212404D0
     p_orb_exp_pddgmndo(atomic_number) = 2.037647D0
     betas_pddgmndo(atomic_number) = -15.663317D0
     betap_pddgmndo(atomic_number) = -15.399331D0
     GSS_pddgmndo(atomic_number) = 15.03D0
     GSP_pddgmndo(atomic_number) = 13.16D0
     GPP_pddgmndo(atomic_number) = 11.30D0
     GP2_pddgmndo(atomic_number) = 9.97D0
     HSP_pddgmndo(atomic_number) = 2.42D0
     alp_pddgmndo(atomic_number) = 2.602846D0
     USS_pddgmndo(atomic_number) = -111.133653D0
     UPP_pddgmndo(atomic_number) = -78.062493D0
     PDDGC1_mndo(atomic_number) = -0.017119D0
     PDDGC2_mndo(atomic_number) = 0.005497D0
     PDDGE1_mndo(atomic_number) = 1.466335D0
     PDDGE2_mndo(atomic_number) = 2.236842D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 3.86491071d0
     p_orb_exp_rm1(atomic_number) = 1.89593144d0
     betas_rm1(atomic_number) = -19.92430432d0
     betap_rm1(atomic_number) = -11.52935197d0
     FN1_rm1(1,atomic_number) = 0.12947108d0
     FN2_rm1(1,atomic_number) = 2.97724424d0
     FN3_rm1(1,atomic_number) = 1.46749784d0
     FN1_rm1(2,atomic_number) = 0.00288899d0
     FN2_rm1(2,atomic_number) = 7.09827589d0
     FN3_rm1(2,atomic_number) = 2.50002723d0
     FN1_rm1(3,atomic_number) = 0.0d0
     FN2_rm1(3,atomic_number) = 0.0d0
     FN3_rm1(3,atomic_number) = 0.0d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 2
     GSS_rm1(atomic_number) = 15.36023105d0
     GSP_rm1(atomic_number) = 13.30671171d0
     GPP_rm1(atomic_number) = 12.56502640d0
     GP2_rm1(atomic_number) = 9.66397083d0
     HSP_rm1(atomic_number) = 1.76489897d0
     alp_rm1(atomic_number) = 3.69358828d0
     USS_rm1(atomic_number) = -118.47306918d0
     UPP_rm1(atomic_number) = -76.35330340d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.637050D0
     p_orb_exp_pm6(atomic_number) = 2.118146D0
     d_orb_exp_pm6(atomic_number) = 1.324033D0
     betas_pm6(atomic_number) = -2.367988D0
     betap_pm6(atomic_number) = -13.802139D0
     betad_pm6(atomic_number) = -4.037751D0
     s_orb_exp_tail_pm6(atomic_number) = 0.956297D0
     p_orb_exp_tail_pm6(atomic_number) = 2.464067D0
     d_orb_exp_tail_pm6(atomic_number) = 6.410325D0
     FN1_pm6(1,atomic_number) = -0.013213D0
     FN2_pm6(1,atomic_number) =  3.687022D0
     FN3_pm6(1,atomic_number) =  2.544635D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  11.142654D0
     GSP_pm6(atomic_number) =  7.487881D0
     GPP_pm6(atomic_number) =  9.551886D0
     GP2_pm6(atomic_number) =  8.128436D0
     HSP_pm6(atomic_number) =  5.004267D0
     USS_pm6(atomic_number) = -61.389930D0
     UPP_pm6(atomic_number) = -54.482801D0
     UDD_pm6(atomic_number) = -38.258155D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 2.5172960D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END CHLORINE
 !-------------------

 !-------------------
 !ARGON
 !-------------------
     atomic_number = 18
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 6.000272D0
     p_orb_exp_pm6(atomic_number) = 5.949170D0
     betas_pm6(atomic_number) =  -8.839842D0 
     betap_pm6(atomic_number) = -28.427303D0 
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 17.858776D0
     GSP_pm6(atomic_number) =  4.168451D0
     GPP_pm6(atomic_number) = 11.852500D0
     GP2_pm6(atomic_number) = 15.669543D0
     HSP_pm6(atomic_number) =  4.574549D0
     USS_pm6(atomic_number) =  -7.797931D0
     UPP_pm6(atomic_number) = -83.211487D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 3.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END ARGON
 !-------------------

 !-------------------
 !POTASSIUM
 !-------------------
     atomic_number = 19
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 6.000478D0
     p_orb_exp_pm6(atomic_number) = 1.127503D0
     betas_pm6(atomic_number) = -8.755195D0
     betap_pm6(atomic_number) = -1.788061D0
     FN1_pm6(1,atomic_number) = 0.157519D0
     FN2_pm6(1,atomic_number) = 6.000566D0
     FN3_pm6(1,atomic_number) = 2.047539D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =  3.369251D0
     GSP_pm6(atomic_number) =  6.129351D0
     GPP_pm6(atomic_number) =  0.999505D0
     GP2_pm6(atomic_number) = 18.999148D0
     HSP_pm6(atomic_number) =  0.300325D0
     USS_pm6(atomic_number) = -3.801108D0
     UPP_pm6(atomic_number) = -3.339656D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.2d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END POTASSIUM
 !-------------------

 !-------------------
 !CALCIUM
 !-------------------
     atomic_number = 20
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.528258D0
     p_orb_exp_pm6(atomic_number) = 2.060094D0
     betas_pm6(atomic_number) = -4.343881D0
     betap_pm6(atomic_number) = -1.296612D0
     FN1_pm6(1,atomic_number) = -0.025275D0
     FN2_pm6(1,atomic_number) =  0.500017D0
     FN3_pm6(1,atomic_number) =  2.329051D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 5.725773D0
     GSP_pm6(atomic_number) = 4.781065D0
     GPP_pm6(atomic_number) = 7.172103D0
     GP2_pm6(atomic_number) = 7.431876D0
     HSP_pm6(atomic_number) = 1.240572D0
     USS_pm6(atomic_number) = -10.770058D0
     UPP_pm6(atomic_number) =  -9.754177D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END CALCIUM
 !-------------------

 !-------------------
 !SCANDIUM
 !-------------------
     atomic_number = 21
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.402469d0
     p_orb_exp_pm6(atomic_number) = 1.345196d0
     d_orb_exp_pm6(atomic_number) = 1.859012d0
     betas_pm6(atomic_number) = -8.620944d0
     betap_pm6(atomic_number) =  3.075948d0
     betad_pm6(atomic_number) = -9.768661d0
     s_orb_exp_tail_pm6(atomic_number) = 0.848418d0
     p_orb_exp_tail_pm6(atomic_number) = 2.451729d0
     d_orb_exp_tail_pm6(atomic_number) = 0.789372d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   4.638215d0
     GSP_pm6(atomic_number) =   5.739164d0
     GPP_pm6(atomic_number) =  14.604872d0
     GP2_pm6(atomic_number) =  12.802595d0
     HSP_pm6(atomic_number) =   0.193835d0
     USS_pm6(atomic_number) = -15.544461d0
     UPP_pm6(atomic_number) = -18.646295d0
     UDD_pm6(atomic_number) = -16.069444d0
     F0SD_pm6(atomic_number) =  4.798313d0
     G2SD_pm6(atomic_number) =  5.380136d0
     rho_core_pm6(atomic_number) = 3.173734d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -33.99955186d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END SCANDIUM
 !-------------------

 !-------------------
 !TITANIUM
 !-------------------
     atomic_number = 22
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 5.324777d0
     p_orb_exp_pm6(atomic_number) = 1.164068d0
     d_orb_exp_pm6(atomic_number) = 1.418280d0
     betas_pm6(atomic_number) =  3.389142d0
     betap_pm6(atomic_number) = -3.355350d0
     betad_pm6(atomic_number) = -1.842829d0
     s_orb_exp_tail_pm6(atomic_number) = 1.045904d0
     p_orb_exp_tail_pm6(atomic_number) = 1.076844d0
     d_orb_exp_tail_pm6(atomic_number) = 0.717945d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 5.717851d0
     GSP_pm6(atomic_number) = 5.800015d0
     GPP_pm6(atomic_number) = 6.414726d0
     GP2_pm6(atomic_number) = 5.623133d0
     HSP_pm6(atomic_number) = 1.403732d0
     USS_pm6(atomic_number) = -25.507973d0
     UPP_pm6(atomic_number) = -17.260909d0
     UDD_pm6(atomic_number) = -23.809486d0
     F0SD_pm6(atomic_number) = 6.56056200d0
     G2SD_pm6(atomic_number) = 3.39623500d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -63.46031221d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END TITANIUM
 !-------------------

 !-------------------
 !VANADIUM
 !-------------------
     atomic_number = 23
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  1.974330d0
     p_orb_exp_pm6(atomic_number) =  1.063106d0
     d_orb_exp_pm6(atomic_number) =  1.394806d0
     betas_pm6(atomic_number) =  -1.211330d0
     betap_pm6(atomic_number) =   0.740746d0
     betad_pm6(atomic_number) =   3.153669d0
     s_orb_exp_tail_pm6(atomic_number) =  1.094426d0
     p_orb_exp_tail_pm6(atomic_number) =  0.755378d0
     d_orb_exp_tail_pm6(atomic_number) =  1.099367d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   5.983116d0
     GSP_pm6(atomic_number) =   4.736769d0
     GPP_pm6(atomic_number) =   4.499763d0
     GP2_pm6(atomic_number) =   3.944481d0
     HSP_pm6(atomic_number) =   0.901105d0
     USS_pm6(atomic_number) = -32.162276d0
     UPP_pm6(atomic_number) = -21.572501d0
     UDD_pm6(atomic_number) = -34.506245d0
     F0SD_pm6(atomic_number) =  6.810021d0
     G2SD_pm6(atomic_number) =  1.831407d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -100.61396315d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END VANADIUM
 !-------------------

 !-------------------
 !CHROMIUM
 !-------------------
     atomic_number = 24
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  3.283460d0
     p_orb_exp_pm6(atomic_number) =  1.029394d0
     d_orb_exp_pm6(atomic_number) =  1.623119d0
     betas_pm6(atomic_number) =  -5.122615d0
     betap_pm6(atomic_number) =   3.926711d0
     betad_pm6(atomic_number) =  -4.230550d0
     s_orb_exp_tail_pm6(atomic_number) =  1.619853d0
     p_orb_exp_tail_pm6(atomic_number) =  0.848266d0
     d_orb_exp_tail_pm6(atomic_number) =  1.405015d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   8.855572d0
     GSP_pm6(atomic_number) =   5.588631d0
     GPP_pm6(atomic_number) =   5.053094d0
     GP2_pm6(atomic_number) =   4.429530d0
     HSP_pm6(atomic_number) =   0.648039d0
     USS_pm6(atomic_number) = -34.864339d0
     UPP_pm6(atomic_number) = -26.978615d0
     UDD_pm6(atomic_number) = -54.431036d0
     F0SD_pm6(atomic_number) =  6.150136d0
     G2SD_pm6(atomic_number) =  2.000300d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -185.72482255d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END CHROMIUM
 !-------------------

 !-------------------
 !MANGANESE
 !-------------------
     atomic_number = 25
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  2.131680d0
     p_orb_exp_pm6(atomic_number) =  1.525880d0
     d_orb_exp_pm6(atomic_number) =  2.607800d0
     betas_pm6(atomic_number) =  -4.185290d0
     betap_pm6(atomic_number) =  -3.479630d0
     betad_pm6(atomic_number) = -13.473190d0 
     s_orb_exp_tail_pm6(atomic_number) =  1.132450d0
     p_orb_exp_tail_pm6(atomic_number) =  1.390740d0
     d_orb_exp_tail_pm6(atomic_number) =  0.962550d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   6.190990d0
     GSP_pm6(atomic_number) =   6.757427d0
     GPP_pm6(atomic_number) =   8.284594d0
     GP2_pm6(atomic_number) =   7.262255d0
     HSP_pm6(atomic_number) =   1.520518d0
     USS_pm6(atomic_number) = -51.460000d0
     UPP_pm6(atomic_number) = -37.543990d0
     UDD_pm6(atomic_number) = -47.655370d0
     F0SD_pm6(atomic_number) =  7.690920d0
     G2SD_pm6(atomic_number) =  1.105330d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -195.80157708d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END MANGANESE
 !-------------------

 !-------------------
 !IRON
 !-------------------
   atomic_number = 26
  !MNDO
   ! Reference: None
     element_supported_mndo(atomic_number) = .false.

  !AM1
   ! Reference: None
     element_supported_am1(atomic_number) = .false.

  !PM3
   ! Reference: J.P. McNamara et al., J. Mol. Gra. Mod., 24, 128, (2005) (Index 32)
   ! NOTE: The PM3 Iron parameters are a restricted set in that they require a pairwise
   !       definition for core-core interactions. Thus at present the only combinations
   !       available are Fe-H, Fe-C, Fe-N, Fe-O, Fe-S, Fe-Cl and Fe-Br. Thus you can only
   !       have one Fe atom in the QM region.
     pm3_ref_index(atomic_number) = 32
     element_supported_pm3(atomic_number) = .false.

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.479150d0
     p_orb_exp_pm6(atomic_number) = 6.002246d0
     d_orb_exp_pm6(atomic_number) = 1.080747d0
     betas_pm6(atomic_number) =  8.027621d0
     betap_pm6(atomic_number) = -1.125760d0
     betad_pm6(atomic_number) = -3.507531d0
     s_orb_exp_tail_pm6(atomic_number) = 1.459152d0
     p_orb_exp_tail_pm6(atomic_number) = 1.392614d0
     d_orb_exp_tail_pm6(atomic_number) = 2.161909d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =    7.977036d0
     GSP_pm6(atomic_number) =    7.786867d0
     GPP_pm6(atomic_number) =    8.295758d0
     GP2_pm6(atomic_number) =    7.272041d0
     HSP_pm6(atomic_number) =    1.880189d0
     USS_pm6(atomic_number) =  -70.515047d0
     UPP_pm6(atomic_number) =  -62.963069d0
     UDD_pm6(atomic_number) = -103.631790d0
     F0SD_pm6(atomic_number) =  9.300165d0
     G2SD_pm6(atomic_number) =  1.601345d0
     rho_core_pm6(atomic_number) = 1.272092d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -426.83526638d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END IRON
 !-------------------

 !-------------------
 !COBALT
 !-------------------
     atomic_number = 27
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  1.166613d0
     p_orb_exp_pm6(atomic_number) =  3.000000d0
     d_orb_exp_pm6(atomic_number) =  1.860218d0
     betas_pm6(atomic_number) =  -8.992062d0
     betap_pm6(atomic_number) =  -0.100000d0
     betad_pm6(atomic_number) =  -2.481509d0
     s_orb_exp_tail_pm6(atomic_number) =  0.519518d0
     p_orb_exp_tail_pm6(atomic_number) =  1.000000d0
     d_orb_exp_tail_pm6(atomic_number) =  0.352115d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   2.840152d0
     GSP_pm6(atomic_number) =   3.425933d0
     GPP_pm6(atomic_number) =   5.956968d0
     GP2_pm6(atomic_number) =   5.221864d0
     HSP_pm6(atomic_number) =   0.390087d0
     USS_pm6(atomic_number) = -21.039413d0
     UPP_pm6(atomic_number) =  10.000000d0
     UDD_pm6(atomic_number) = -28.068971d0
     F0SD_pm6(atomic_number) =  1.446283d0
     G2SD_pm6(atomic_number) =  1.680225d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -167.65660749d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END COBALT
 !-------------------

 !-------------------
 !NICKEL
 !-------------------
     atomic_number = 28
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.591828d0
     p_orb_exp_pm6(atomic_number) = 2.304739d0
     d_orb_exp_pm6(atomic_number) = 2.514761d0
     betas_pm6(atomic_number) = -9.151521d0
     betap_pm6(atomic_number) = -8.086696d0
     betad_pm6(atomic_number) = -8.655910d0
     s_orb_exp_tail_pm6(atomic_number) = 0.746470d0
     p_orb_exp_tail_pm6(atomic_number) = 0.753327d0
     d_orb_exp_tail_pm6(atomic_number) = 1.461345d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   4.080876d0
     GSP_pm6(atomic_number) =   4.099452d0
     GPP_pm6(atomic_number) =   4.487545d0
     GP2_pm6(atomic_number) =   3.933771d0
     HSP_pm6(atomic_number) =   0.993498d0
     USS_pm6(atomic_number) = -47.620247d0
     UPP_pm6(atomic_number) = -32.878408d0
     UDD_pm6(atomic_number) = -93.026395d0
     F0SD_pm6(atomic_number) =  4.651664d0
     G2SD_pm6(atomic_number) =  1.880502d0
     rho_core_pm6(atomic_number) = 1.586979d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -485.16555989d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END NICKEL
 !-------------------

 !-------------------
 !COPPER
 !-------------------
     atomic_number = 29
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  1.669096d0
     p_orb_exp_pm6(atomic_number) =  3.000000d0
     d_orb_exp_pm6(atomic_number) =  2.734990d0
     betas_pm6(atomic_number) =  -9.369508d0
     betap_pm6(atomic_number) =  -0.100000d0
     betad_pm6(atomic_number) = -16.982092d0 
     s_orb_exp_tail_pm6(atomic_number) =  1.899598d0
     p_orb_exp_tail_pm6(atomic_number) =  3.000000d0
     d_orb_exp_tail_pm6(atomic_number) =  1.484317d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   10.384910d0
     GSP_pm6(atomic_number) =   12.145361d0
     GPP_pm6(atomic_number) =   17.870905d0
     GP2_pm6(atomic_number) =   15.665592d0
     HSP_pm6(atomic_number) =    2.037394d0
     USS_pm6(atomic_number) =  -92.002205d0
     UPP_pm6(atomic_number) =   -1.000000d0
     UDD_pm6(atomic_number) = -110.442592d0
     F0SD_pm6(atomic_number) =  9.848807d0
     G2SD_pm6(atomic_number) =  9.847577d0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     ! EISOL should be calculated from the above parameters...
     EISOL_pm6(atomic_number) = -656.59528257d0 ! taken from MOPAC2009 output
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END COPPER
 !-------------------

 !-------------------
 !ZINC
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0, GSSC = 1, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 30
  !MNDO
   ! Reference: M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 5, 1494-1496 (1986) (index = 10)
     mndo_ref_index(atomic_number) = 10
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.0473590D0
     p_orb_exp_mndo(atomic_number) = 1.4609460D0
     betas_mndo(atomic_number) = -1.0000000D0
     betap_mndo(atomic_number) = -2.0000000D0
     GSS_mndo(atomic_number) = 11.8000000D0
     GSP_mndo(atomic_number) = 11.1820180D0
     GPP_mndo(atomic_number) = 13.3000000D0
     GP2_mndo(atomic_number) = 12.9305200D0
     HSP_mndo(atomic_number) = 0.4846060D0
     alp_mndo(atomic_number) = 1.5064570D0
     USS_mndo(atomic_number) = -20.8397160D0
     UPP_mndo(atomic_number) = -19.6252240D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =   1.73150352D0
     p_orb_exp_mndod(atomic_number)      =   1.39358305D0
     s_orb_exp_tail_mndod(atomic_number) =   1.56600000D0
     p_orb_exp_tail_mndod(atomic_number) =   0.86283981D0
     betas_mndod(atomic_number)          =  -5.01726076D0 
     betap_mndod(atomic_number)          =  -0.71205972D0
     rho_core_mndod(atomic_number)       =   1.58923393D0
     GSS_mndod(atomic_number)            =   8.56072836D0
     GSP_mndod(atomic_number)            =   7.49003598D0
     GPP_mndod(atomic_number)            =   5.13964830D0
     GP2_mndod(atomic_number)            =   4.50540309D0
     HSP_mndod(atomic_number)            =   0.53294610D0
     alp_mndod(atomic_number)            =   1.51763697D0
     USS_mndod(atomic_number)            = -18.02300000D0
     UPP_mndod(atomic_number)            = -12.24216585D0

  !AM1
   ! Reference: M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 7, 522-524 (1988) (Index 23)
     am1_ref_index(atomic_number) = 23
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.9542990D0
     p_orb_exp_am1(atomic_number) = 1.3723650D0
     betas_am1(atomic_number) = -1.9974290D0
     betap_am1(atomic_number) = -4.7581190D0
     FN1_am1(1,atomic_number) = 0.0000000D0
     FN2_am1(1,atomic_number) = 0.0000000D0
     FN3_am1(1,atomic_number) = 0.0000000D0
     FN1_am1(2,atomic_number) = 0.0000000D0
     FN2_am1(2,atomic_number) = 0.0000000D0
     FN3_am1(2,atomic_number) = 0.0000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 0
     GSS_am1(atomic_number) = 11.8000000D0
     GSP_am1(atomic_number) = 11.1820180D0
     GPP_am1(atomic_number) = 13.3000000D0
     GP2_am1(atomic_number) = 12.9305200D0
     HSP_am1(atomic_number) = 0.4846060D0
     alp_am1(atomic_number) = 1.4845630D0
     USS_am1(atomic_number) = -21.0400080D0
     UPP_am1(atomic_number) = -17.6555740D0
     
  !AM1D
   ! Reference: M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 7, 522-524 (1988) (Index 23)
     am1d_ref_index(atomic_number) = 23
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.9542990D0
     p_orb_exp_am1d(atomic_number) = 1.3723650D0
     betas_am1d(atomic_number) = -1.9974290D0
     betap_am1d(atomic_number) = -4.7581190D0
     FN1_am1d(1,atomic_number) = 0.0000000D0
     FN2_am1d(1,atomic_number) = 0.0000000D0
     FN3_am1d(1,atomic_number) = 0.0000000D0
     FN1_am1d(2,atomic_number) = 0.0000000D0
     FN2_am1d(2,atomic_number) = 0.0000000D0
     FN3_am1d(2,atomic_number) = 0.0000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 0
     GSS_am1d(atomic_number) = 11.8000000D0
     GSP_am1d(atomic_number) = 11.1820180D0
     GPP_am1d(atomic_number) = 13.3000000D0
     GP2_am1d(atomic_number) = 12.9305200D0
     HSP_am1d(atomic_number) = 0.4846060D0
     alp_am1d(atomic_number) = 1.4845630D0
     USS_am1d(atomic_number) = -21.0400080D0
     UPP_am1d(atomic_number) = -17.6555740D0
     

  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.8199890D0
     p_orb_exp_pm3(atomic_number) = 1.5069220D0
     betas_pm3(atomic_number) = -0.7155780D0
     betap_pm3(atomic_number) = -6.3518640D0
     FN1_pm3(1,atomic_number) = -0.1112340D0
     FN2_pm3(1,atomic_number) = 6.0014780D0
     FN3_pm3(1,atomic_number) = 1.5160320D0
     FN1_pm3(2,atomic_number) = -0.1323700D0
     FN2_pm3(2,atomic_number) = 1.9958390D0
     FN3_pm3(2,atomic_number) = 2.5196420D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 9.6771960D0
     GSP_pm3(atomic_number) = 7.7362040D0
     GPP_pm3(atomic_number) = 4.9801740D0
     GP2_pm3(atomic_number) = 4.6696560D0
     HSP_pm3(atomic_number) = 0.6004130D0
     alp_pm3(atomic_number) = 1.3501260D0
     USS_pm3(atomic_number) = -18.5321980D0
     UPP_pm3(atomic_number) = -11.0474090D0

!  !PM3 with ZnB parameters
   ! Reference: E.N. Brothers, D. Suarez, D. W. Deerfield II and K. Merz,
   !            J. Comp. Chem. 25, 1677, (2004)(Index 36)
   ! Note: PM3-ZnB is PM3 with special parameters for Zn for biological systems
   if (currentTheory%PM3ZNB) then
     pm3znb_ref_index(atomic_number) = 36
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.434189D0
     p_orb_exp_pm3(atomic_number) = 1.454582D0
     betas_pm3(atomic_number) = -0.570514D0
     betap_pm3(atomic_number) = -4.124587D0
     FN1_pm3(1,atomic_number) = -0.262170D0
     FN2_pm3(1,atomic_number) =  4.730939D0
     FN3_pm3(1,atomic_number) =  1.802245D0
     FN1_pm3(2,atomic_number) = -0.132917D0
     FN2_pm3(2,atomic_number) =  0.959929D0
     FN3_pm3(2,atomic_number) =  2.382463D0
     FN1_pm3(3,atomic_number) = 0.0D0      
     FN2_pm3(3,atomic_number) = 0.0D0      
     FN3_pm3(3,atomic_number) = 0.0D0      
     FN1_pm3(4,atomic_number) = 0.0D0      
     FN2_pm3(4,atomic_number) = 0.0D0      
     FN3_pm3(4,atomic_number) = 0.0D0      
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 8.496571D0  
     GSP_pm3(atomic_number) = 8.301945D0  
     GPP_pm3(atomic_number) = 6.485329D0  
     GP2_pm3(atomic_number) = 6.235893D0  
     HSP_pm3(atomic_number) = 0.802358D0  
     alp_pm3(atomic_number) = 1.360252D0  
     USS_pm3(atomic_number) = -16.974636D0
     UPP_pm3(atomic_number) =  -9.794156D0
   end if

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.512875D0
     p_orb_exp_pm6(atomic_number) = 1.789482D0
     betas_pm6(atomic_number) = -13.276583D0
     betap_pm6(atomic_number) =   1.479642D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  8.707424D0
     GSP_pm6(atomic_number) =  3.436116D0
     GPP_pm6(atomic_number) = 20.000041D0
     GP2_pm6(atomic_number) =  6.782785D0
     HSP_pm6(atomic_number) =  0.662036D0
     USS_pm6(atomic_number) = -18.040862D0
     UPP_pm6(atomic_number) =  -7.834895D0
     ! alp_pm6 is the same value as PM3
     alp_pm6(atomic_number) = 1.3501260D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END ZINC
 !-------------------

 !-------------------
 !GALLIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 1, GSSC = 1, GPPC = 0, GSPC = 2, GP2C = 0, HSP = -1
   atomic_number = 31
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.8470400D0
     p_orb_exp_pm3(atomic_number) = 0.8394110D0
     betas_pm3(atomic_number) = -4.9456180D0
     betap_pm3(atomic_number) = -0.4070530D0
     FN1_pm3(1,atomic_number) = -0.5601790D0
     FN2_pm3(1,atomic_number) = 5.6232730D0
     FN3_pm3(1,atomic_number) = 1.5317800D0
     FN1_pm3(2,atomic_number) = -0.2727310D0
     FN2_pm3(2,atomic_number) = 1.9918430D0
     FN3_pm3(2,atomic_number) = 2.1838640D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 8.4585540D0
     GSP_pm3(atomic_number) = 8.9256190D0
     GPP_pm3(atomic_number) = 5.0868550D0
     GP2_pm3(atomic_number) = 4.9830450D0
     HSP_pm3(atomic_number) = 2.0512600D0
     alp_pm3(atomic_number) = 1.6051150D0
     USS_pm3(atomic_number) = -29.8555930D0
     UPP_pm3(atomic_number) = -21.8753710D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.339067D0
     p_orb_exp_pm6(atomic_number) = 1.729592D0
     betas_pm6(atomic_number) = -10.808320D0
     betap_pm6(atomic_number) =  -4.185500D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 10.354885D0
     GSP_pm6(atomic_number) =  7.993674D0
     GPP_pm6(atomic_number) =  6.090184D0
     GP2_pm6(atomic_number) =  6.299226D0
     HSP_pm6(atomic_number) =  1.295974D0
     USS_pm6(atomic_number) = -30.600226D0
     UPP_pm6(atomic_number) = -21.032425D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.6051150D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END GALLIUM
 !-------------------

 !-------------------
 !GERMANIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 2, GSSC = 1, GPPC = -0.5, GSPC = 4, GP2C = 1.5, HSP = -2
   atomic_number = 32
  !MNDO
   ! Reference: M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,ORGANOMETALLICS 6 186-189, (1987) (index = 11)
     mndo_ref_index(atomic_number) = 11
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 1.2931800D0
     p_orb_exp_mndo(atomic_number) = 2.0205640D0
     betas_mndo(atomic_number) = -4.5164790D0
     betap_mndo(atomic_number) = -1.7555170D0
     GSS_mndo(atomic_number) = 9.8000000D0
     GSP_mndo(atomic_number) = 8.3000000D0
     GPP_mndo(atomic_number) = 7.3000000D0
     GP2_mndo(atomic_number) = 6.5000000D0
     HSP_mndo(atomic_number) = 1.3000000D0
     alp_mndo(atomic_number) = 1.9784980D0
     USS_mndo(atomic_number) = -33.9493670D0
     UPP_mndo(atomic_number) = -27.4251050D0

  !MNDOD
   ! Reference: M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,ORGANOMETALLICS 6 186-189, (1987) (index = 11)
     mndod_ref_index(atomic_number) = 11
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 1.2931800D0
     p_orb_exp_mndod(atomic_number) = 2.0205640D0
     betas_mndod(atomic_number) = -4.5164790D0
     betap_mndod(atomic_number) = -1.7555170D0
     GSS_mndod(atomic_number) = 9.8000000D0
     GSP_mndod(atomic_number) = 8.3000000D0
     GPP_mndod(atomic_number) = 7.3000000D0
     GP2_mndod(atomic_number) = 6.5000000D0
     HSP_mndod(atomic_number) = 1.3000000D0
     alp_mndod(atomic_number) = 1.9784980D0
     USS_mndod(atomic_number) = -33.9493670D0
     UPP_mndod(atomic_number) = -27.4251050D0

  !AM1
   ! Reference: M.J.S.Dewar and C.Jie, Organometallics, 8, 1544, (1989) (index = 24)
     am1_ref_index(atomic_number) = 24
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 1.2196310D0
     p_orb_exp_am1(atomic_number) = 1.9827940D0
     betas_am1(atomic_number) = -4.3566070D0
     betap_am1(atomic_number) = -0.9910910D0
     FN1_am1(1,atomic_number) = 0.0000000D0
     FN2_am1(1,atomic_number) = 0.0000000D0
     FN3_am1(1,atomic_number) = 0.0000000D0
     FN1_am1(2,atomic_number) = 0.0000000D0
     FN2_am1(2,atomic_number) = 0.0000000D0
     FN3_am1(2,atomic_number) = 0.0000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 0
     GSS_am1(atomic_number) = 10.1686050D0
     GSP_am1(atomic_number) = 8.1444730D0
     GPP_am1(atomic_number) = 6.6719020D0
     GP2_am1(atomic_number) = 6.2697060D0
     HSP_am1(atomic_number) = 0.9370930D0
     alp_am1(atomic_number) = 2.1364050D0
     USS_am1(atomic_number) = -34.1838890D0
     UPP_am1(atomic_number) = -28.6408110D0

  !AM1D
   ! Reference: M.J.S.Dewar and C.Jie, Organometallics, 8, 1544, (1989) (index = 24)
     am1d_ref_index(atomic_number) = 24
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 1.2196310D0
     p_orb_exp_am1d(atomic_number) = 1.9827940D0
     betas_am1d(atomic_number) = -4.3566070D0
     betap_am1d(atomic_number) = -0.9910910D0
     FN1_am1d(1,atomic_number) = 0.0000000D0
     FN2_am1d(1,atomic_number) = 0.0000000D0
     FN3_am1d(1,atomic_number) = 0.0000000D0
     FN1_am1d(2,atomic_number) = 0.0000000D0
     FN2_am1d(2,atomic_number) = 0.0000000D0
     FN3_am1d(2,atomic_number) = 0.0000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 0
     GSS_am1d(atomic_number) = 10.1686050D0
     GSP_am1d(atomic_number) = 8.1444730D0
     GPP_am1d(atomic_number) = 6.6719020D0
     GP2_am1d(atomic_number) = 6.2697060D0
     HSP_am1d(atomic_number) = 0.9370930D0
     alp_am1d(atomic_number) = 2.1364050D0
     USS_am1d(atomic_number) = -34.1838890D0
     UPP_am1d(atomic_number) = -28.6408110D0

  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.2373526D0
     p_orb_exp_pm3(atomic_number) = 1.5924319D0
     betas_pm3(atomic_number) = -5.3250024D0
     betap_pm3(atomic_number) = -2.2501567D0
     FN1_pm3(1,atomic_number) = 0.9631726D0
     FN2_pm3(1,atomic_number) = 6.0120134D0
     FN3_pm3(1,atomic_number) = 2.1633655D0
     FN1_pm3(2,atomic_number) =-0.9593891D0
     FN2_pm3(2,atomic_number) = 5.7491802D0
     FN3_pm3(2,atomic_number) = 2.1693724D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 5.3769635D0
     GSP_pm3(atomic_number) = 10.2095293D0
     GPP_pm3(atomic_number) = 7.6718647D0
     GP2_pm3(atomic_number) = 6.9242663D0
     HSP_pm3(atomic_number) = 1.3370204D0
     alp_pm3(atomic_number) = 1.9723370D0
     USS_pm3(atomic_number) = -35.4671955D0
     UPP_pm3(atomic_number) = -31.5863583D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.546073D0
     p_orb_exp_pm6(atomic_number) = 1.709130D0
     betas_pm6(atomic_number) = -14.854297D0
     betap_pm6(atomic_number) =  -2.591260D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 7.518301D0
     GSP_pm6(atomic_number) = 6.594443D0
     GPP_pm6(atomic_number) = 6.066801D0
     GP2_pm6(atomic_number) = 5.305947D0
     HSP_pm6(atomic_number) = 0.290742D0
     USS_pm6(atomic_number) = -32.747338D0
     UPP_pm6(atomic_number) = -24.709016D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.9723370D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END GERMANIUM
 !-------------------

 !-------------------
 !ARSENIC
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 3, GSSC = 1, GPPC = -1.5, GSPC = 6, GP2C = 4.5, HSP = -3
   atomic_number = 33
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.6361770D0
     p_orb_exp_pm3(atomic_number) = 1.7038890D0
     betas_pm3(atomic_number) = -8.2321650D0
     betap_pm3(atomic_number) = -5.0173860D0
     FN1_pm3(1,atomic_number) =-0.4600950D0
     FN2_pm3(1,atomic_number) = 1.9831150D0
     FN3_pm3(1,atomic_number) = 1.0867930D0
     FN1_pm3(2,atomic_number) =-0.0889960D0
     FN2_pm3(2,atomic_number) = 1.9929440D0
     FN3_pm3(2,atomic_number) = 2.1400580D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 8.7890010D0
     GSP_pm3(atomic_number) = 5.3979830D0
     GPP_pm3(atomic_number) = 8.2872500D0
     GP2_pm3(atomic_number) = 8.2103460D0
     HSP_pm3(atomic_number) = 1.9510340D0
     alp_pm3(atomic_number) = 1.7944770D0
     USS_pm3(atomic_number) = -38.5074240D0
     UPP_pm3(atomic_number) = -35.1524150D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.926171d0
     p_orb_exp_pm6(atomic_number) = 1.765191d0
     d_orb_exp_pm6(atomic_number) = 1.392142d0
     betas_pm6(atomic_number) =  -11.963725d0
     betap_pm6(atomic_number) =   -7.340073d0
     betad_pm6(atomic_number) =    3.753005d0
     s_orb_exp_tail_pm6(atomic_number) =  2.006543d0
     p_orb_exp_tail_pm6(atomic_number) =  3.316832d0
     d_orb_exp_tail_pm6(atomic_number) =  4.653440d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =    6.665030d0
     GSP_pm6(atomic_number) =    6.213867d0
     GPP_pm6(atomic_number) =    9.310836d0
     GP2_pm6(atomic_number) =    8.721542d0
     HSP_pm6(atomic_number) =    0.280662d0
     USS_pm6(atomic_number) =  -37.956965d0
     UPP_pm6(atomic_number) =  -38.453701d0
     UDD_pm6(atomic_number) =  -30.282658d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.7944770D0 
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END ARSENIC
 !-------------------

 !-------------------
 !SELENIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 4, GSSC = 1, GPPC = -0.5, GSPC = 8, GP2C = 6.5, HSP = -4
   atomic_number = 34
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.8280510D0
     p_orb_exp_pm3(atomic_number) = 1.7325360D0
     betas_pm3(atomic_number) = -6.1578220D0
     betap_pm3(atomic_number) = -5.4930390D0
     FN1_pm3(1,atomic_number) = 0.0478730D0
     FN2_pm3(1,atomic_number) = 6.0074000D0
     FN3_pm3(1,atomic_number) = 2.0817170D0
     FN1_pm3(2,atomic_number) = 0.1147200D0
     FN2_pm3(2,atomic_number) = 6.0086720D0
     FN3_pm3(2,atomic_number) = 1.5164230D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 7.4325910D0
     GSP_pm3(atomic_number) = 10.0604610D0
     GPP_pm3(atomic_number) = 9.5683260D0
     GP2_pm3(atomic_number) = 7.7242890D0
     HSP_pm3(atomic_number) = 4.0165580D0
     alp_pm3(atomic_number) = 3.0439570D0
     USS_pm3(atomic_number) = -55.3781350D0
     UPP_pm3(atomic_number) = -49.8230760D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.512366d0
     p_orb_exp_pm6(atomic_number) = 2.007576d0
     betas_pm6(atomic_number) =  2.636001d0
     betap_pm6(atomic_number) = -9.557700d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =   5.522356d0
     GSP_pm6(atomic_number) =   2.907562d0
     GPP_pm6(atomic_number) =   8.042391d0
     GP2_pm6(atomic_number) =   6.735106d0
     HSP_pm6(atomic_number) =   3.095789d0
     USS_pm6(atomic_number) = -32.671088d0
     UPP_pm6(atomic_number) = -32.522220d0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 3.0439570D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END SELENIUM
 !-------------------

 !-------------------
 !BROMINE
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 5, GSSC = 1, GPPC = 0, GSPC = 10, GP2C = 10, HSP = -5
   atomic_number = 35
  !MNDO
   ! Reference: M.J.S.DEWAR, E.F. HEALY, J. COMP. CHEM., 4, 542, (1983)  (index = 12)
     mndo_ref_index(atomic_number) = 12
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 3.8543019D0
     p_orb_exp_mndo(atomic_number) = 2.1992091D0
     betas_mndo(atomic_number) = -8.9171070D0
     betap_mndo(atomic_number) = -9.9437400D0
     GSS_mndo(atomic_number) = 15.03643948D0
     GSP_mndo(atomic_number) = 13.03468242D0
     GPP_mndo(atomic_number) = 11.27632539D0
     GP2_mndo(atomic_number) = 9.85442552D0
     HSP_mndo(atomic_number) = 2.45586832D0
     alp_mndo(atomic_number) = 2.4457051D0
     USS_mndo(atomic_number) = -99.9864405D0
     UPP_mndo(atomic_number) = -75.6713075D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =   2.59054101D0
     p_orb_exp_mndod(atomic_number)      =   2.33085649D0
     d_orb_exp_mndod(atomic_number)      =   1.35736115D0
     betas_mndod(atomic_number)          =  -8.31497607D0
     betap_mndod(atomic_number)          = -10.50704145D0
     betad_mndod(atomic_number)          =  -0.96259930D0
     s_orb_exp_tail_mndod(atomic_number) =   2.23581544D0
     p_orb_exp_tail_mndod(atomic_number) =   1.43292654D0
     d_orb_exp_tail_mndod(atomic_number) =   1.24257826D0
     GSS_mndod(atomic_number)            =  12.22235546D0
     GSP_mndod(atomic_number)            =   8.26372010D0
     GPP_mndod(atomic_number)            =   8.53546437D0
     GDD_mndod(atomic_number)            =   7.31095300D0
     GP2_mndod(atomic_number)            =   7.48216712D0
     HSP_mndod(atomic_number)            =   2.74952230D0
     alp_mndod(atomic_number)            =   2.09105000D0
     USS_mndod(atomic_number)            = -65.40277790D0
     UPP_mndod(atomic_number)            = -54.55375352D0
     UDD_mndod(atomic_number)            = -13.72809929D0

  !AM1
   ! Reference:  Br: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1_ref_index(atomic_number) = 18
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 3.0641330D0
     p_orb_exp_am1(atomic_number) = 2.0383330D0
     betas_am1(atomic_number) = -19.3998800D0
     betap_am1(atomic_number) = -8.9571950D0
     FN1_am1(1,atomic_number) = 0.0666850D0
     FN2_am1(1,atomic_number) = 4.0000000D0
     FN3_am1(1,atomic_number) = 1.5000000D0
     FN1_am1(2,atomic_number) = 0.0255680D0
     FN2_am1(2,atomic_number) = 4.0000000D0
     FN3_am1(2,atomic_number) = 2.3000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 2
     GSS_am1(atomic_number) = 15.0364395D0
     GSP_am1(atomic_number) = 13.0346824D0
     GPP_am1(atomic_number) = 11.2763254D0
     GP2_am1(atomic_number) = 9.8544255D0
     HSP_am1(atomic_number) = 2.4558683D0
     alp_am1(atomic_number) = 2.5765460D0
     USS_am1(atomic_number) = -104.6560630D0
     UPP_am1(atomic_number) = -74.9300520D0

  !AM1D
   ! Reference:  Br: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1d_ref_index(atomic_number) = 18
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 3.0641330D0
     p_orb_exp_am1d(atomic_number) = 2.0383330D0
     betas_am1d(atomic_number) = -19.3998800D0
     betap_am1d(atomic_number) = -8.9571950D0
     FN1_am1d(1,atomic_number) = 0.0666850D0
     FN2_am1d(1,atomic_number) = 4.0000000D0
     FN3_am1d(1,atomic_number) = 1.5000000D0
     FN1_am1d(2,atomic_number) = 0.0255680D0
     FN2_am1d(2,atomic_number) = 4.0000000D0
     FN3_am1d(2,atomic_number) = 2.3000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 2
     GSS_am1d(atomic_number) = 15.0364395D0
     GSP_am1d(atomic_number) = 13.0346824D0
     GPP_am1d(atomic_number) = 11.2763254D0
     GP2_am1d(atomic_number) = 9.8544255D0
     HSP_am1d(atomic_number) = 2.4558683D0
     alp_am1d(atomic_number) = 2.5765460D0
     USS_am1d(atomic_number) = -104.6560630D0
     UPP_am1d(atomic_number) = -74.9300520D0

  !PM3
   ! Reference:  J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989). (Index 26)
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 5.3484570D0
     p_orb_exp_pm3(atomic_number) = 2.1275900D0
     betas_pm3(atomic_number) = -31.1713420D0
     betap_pm3(atomic_number) = -6.8140130D0
     FN1_pm3(1,atomic_number) = 0.9604580D0
     FN2_pm3(1,atomic_number) = 5.9765080D0
     FN3_pm3(1,atomic_number) = 2.3216540D0
     FN1_pm3(2,atomic_number) =-0.9549160D0
     FN2_pm3(2,atomic_number) = 5.9447030D0
     FN3_pm3(2,atomic_number) = 2.3281420D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 15.9434250D0
     GSP_pm3(atomic_number) = 16.0616800D0
     GPP_pm3(atomic_number) = 8.2827630D0
     GP2_pm3(atomic_number) = 7.8168490D0
     HSP_pm3(atomic_number) = 0.5788690D0
     alp_pm3(atomic_number) = 2.5118420D0
     USS_pm3(atomic_number) = -116.6193110D0
     UPP_pm3(atomic_number) = -74.2271290D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgpm3_ref_index(atomic_number) = 30
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -351.013887D0
     s_orb_exp_pddgpm3(atomic_number) = 4.345079D0
     p_orb_exp_pddgpm3(atomic_number) = 2.190961D0
     betas_pddgpm3(atomic_number) = -21.538044D0
     betap_pddgpm3(atomic_number) = -8.524764D0
     FN1_pddgpm3(1,atomic_number) = 0.961362D0
     FN2_pddgpm3(1,atomic_number) = 6.013600D0
     FN3_pddgpm3(1,atomic_number) = 2.340445D0
     FN1_pddgpm3(2,atomic_number) = -0.948834D0
     FN2_pddgpm3(2,atomic_number) = 5.976329D0
     FN3_pddgpm3(2,atomic_number) = 2.348745D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 15.9434250D0
     GSP_pddgpm3(atomic_number) = 16.0616800D0
     GPP_pddgpm3(atomic_number) = 8.2827630D0
     GP2_pddgpm3(atomic_number) = 7.8168490D0
     HSP_pddgpm3(atomic_number) = 0.5788690D0
     alp_pddgpm3(atomic_number) = 2.424673D0
     USS_pddgpm3(atomic_number) = -115.841963D0
     UPP_pddgpm3(atomic_number) = -74.205146D0
     PDDGC1_pm3(atomic_number) = -0.013772D0
     PDDGC2_pm3(atomic_number) = 0.008849D0
     PDDGE1_pm3(atomic_number) = 1.852030D0
     PDDGE2_pm3(atomic_number) = 2.338958D0

  !PDDG/MNDO
   ! Reference: Repasky, Chandrasekhar and Jorgensen, 2002, J Comp Chem, 23: 1601-1622 (Index 30)
     pddgmndo_ref_index(atomic_number) = 30
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -349.564096D0
     s_orb_exp_pddgmndo(atomic_number) = 3.999975D0
     p_orb_exp_pddgmndo(atomic_number) = 2.245040D0
     betas_pddgmndo(atomic_number) = -7.054170D0
     betap_pddgmndo(atomic_number) = -10.221030D0
     GSS_pddgmndo(atomic_number) = 15.03643948D0
     GSP_pddgmndo(atomic_number) = 13.03468242D0
     GPP_pddgmndo(atomic_number) = 11.27632539D0
     GP2_pddgmndo(atomic_number) = 9.85442552D0
     HSP_pddgmndo(atomic_number) = 2.45586832D0
     alp_pddgmndo(atomic_number) = 2.414265D0
     USS_pddgmndo(atomic_number) = -100.637007D0
     UPP_pddgmndo(atomic_number) = -76.015735D0
     PDDGC1_mndo(atomic_number) = -0.017133D0
     PDDGC2_mndo(atomic_number) = -0.016964D0
     PDDGE1_mndo(atomic_number) = 2.201539D0
     PDDGE2_mndo(atomic_number) = 2.255764D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 5.73157215d0
     p_orb_exp_rm1(atomic_number) = 2.03147582d0
     betas_rm1(atomic_number) = -1.34139841d0
     betap_rm1(atomic_number) = -8.20225991d0
     FN1_rm1(1,atomic_number) = 0.98689937d0
     FN2_rm1(1,atomic_number) = 4.28484191d0
     FN3_rm1(1,atomic_number) = 2.00019696d0
     FN1_rm1(2,atomic_number) = -0.92731247d0
     FN2_rm1(2,atomic_number) = 4.54004910d0
     FN3_rm1(2,atomic_number) = 2.01617695d0
     FN1_rm1(3,atomic_number) = 0.0d0
     FN2_rm1(3,atomic_number) = 0.0d0
     FN3_rm1(3,atomic_number) = 0.0d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 2
     GSS_rm1(atomic_number) = 17.11563074d0
     GSP_rm1(atomic_number) = 15.62419251d0
     GPP_rm1(atomic_number) = 10.73546293d0
     GP2_rm1(atomic_number) = 8.86056199d0
     HSP_rm1(atomic_number) = 2.23512762d0
     alp_rm1(atomic_number) = 2.86710531d0
     USS_rm1(atomic_number) = -113.48398183d0
     UPP_rm1(atomic_number) = -76.18720023d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  4.670684d0
     p_orb_exp_pm6(atomic_number) =  2.035626d0
     d_orb_exp_pm6(atomic_number) =  1.521031d0
     betas_pm6(atomic_number) =  -32.131665d0
     betap_pm6(atomic_number) =   -9.514484d0
     betad_pm6(atomic_number) =   -9.839124d0
     s_orb_exp_tail_pm6(atomic_number) =  3.094777d0
     p_orb_exp_tail_pm6(atomic_number) =  3.065764d0
     d_orb_exp_tail_pm6(atomic_number) =  2.820003d0
     FN1_pm6(1,atomic_number) = -0.004996d0
     FN2_pm6(1,atomic_number) =  6.001292d0
     FN3_pm6(1,atomic_number) =  2.895153d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 7.616791d0
     GSP_pm6(atomic_number) = 5.010425d0
     GPP_pm6(atomic_number) = 9.649216d0
     GP2_pm6(atomic_number) = 8.343792d0
     HSP_pm6(atomic_number) = 4.996553d0
     USS_pm6(atomic_number) =  -45.834364d0
     UPP_pm6(atomic_number) =  -50.293675d0
     UDD_pm6(atomic_number) =    7.086738d0
     ! alp_pm6 is taken from PM3
     alp_pm6(atomic_number) = 2.5118420D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END BROMINE
 !-------------------

 !-------------------
 !KRYPTON
 !-------------------
     atomic_number = 36
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.312248D0
     p_orb_exp_pm6(atomic_number) = 4.491371D0
     betas_pm6(atomic_number) =  -2.727088D0
     betap_pm6(atomic_number) = -16.142951D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 19.999857D0
     GSP_pm6(atomic_number) =  1.175304D0
     GPP_pm6(atomic_number) =  9.174784D0
     GP2_pm6(atomic_number) = 14.926948D0
     HSP_pm6(atomic_number) =  0.299867D0
     USS_pm6(atomic_number) =   8.535384D0
     UPP_pm6(atomic_number) = -80.484321D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 3.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END KRYPTON
 !-------------------

 !-------------------
 !RUBIDIUM
 !-------------------
     atomic_number = 37
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 5.510145D0
     p_orb_exp_pm6(atomic_number) = 1.335170D0
     betas_pm6(atomic_number) = 9.998744D0
     betap_pm6(atomic_number) = 1.343004D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  6.680824D0
     GSP_pm6(atomic_number) = 20.001098D0
     GPP_pm6(atomic_number) =  5.068874D0
     GP2_pm6(atomic_number) =  2.747860D0
     HSP_pm6(atomic_number) =  3.602834D0
     USS_pm6(atomic_number) = -3.636505D0
     UPP_pm6(atomic_number) = -2.500671D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.2d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END RUBIDIUM
 !-------------------

 !-------------------
 !STRONTIUM
 !-------------------
     atomic_number = 38
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.197303D0
     p_orb_exp_pm6(atomic_number) = 1.730137D0
     betas_pm6(atomic_number) = -6.253108D0
     betap_pm6(atomic_number) = -9.844498D0
     FN1_pm6(1,atomic_number) = -0.012948D0
     FN2_pm6(1,atomic_number) =  6.000126D0
     FN3_pm6(1,atomic_number) =  3.011964D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 4.603664D0
     GSP_pm6(atomic_number) = 5.716069D0
     GPP_pm6(atomic_number) = 7.334620D0
     GP2_pm6(atomic_number) = 7.443088D0
     HSP_pm6(atomic_number) = 0.831527D0
     USS_pm6(atomic_number) = -10.427671D0
     UPP_pm6(atomic_number) =  -9.943751D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END STRONTIUM
 !-------------------

 !-------------------
 !CADMIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0, GSSC = 1, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 48
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.6793510D0
     p_orb_exp_pm3(atomic_number) = 2.0664120D0
     betas_pm3(atomic_number) = -8.5819440D0
     betap_pm3(atomic_number) = -0.6010340D0
     FN1_pm3(1,atomic_number) = 0.0d0
     FN2_pm3(1,atomic_number) = 0.0d0
     FN3_pm3(1,atomic_number) = 0.0d0
     FN1_pm3(2,atomic_number) = 0.0d0
     FN2_pm3(2,atomic_number) = 0.0d0
     FN3_pm3(2,atomic_number) = 0.0d0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 0
     GSS_pm3(atomic_number) = 9.2069600D0
     GSP_pm3(atomic_number) = 8.2315390D0
     GPP_pm3(atomic_number) = 4.9481040D0
     GP2_pm3(atomic_number) = 4.6696560D0
     HSP_pm3(atomic_number) = 1.6562340D0
     alp_pm3(atomic_number) = 1.5253820D0
     USS_pm3(atomic_number) = -15.8285840D0
     UPP_pm3(atomic_number) = 8.7497950D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =   1.74880559D0
     p_orb_exp_mndod(atomic_number)      =   1.56321473D0
     s_orb_exp_tail_mndod(atomic_number) =   1.76314840D0
     p_orb_exp_tail_mndod(atomic_number) =   1.52551900D0
     betas_mndod(atomic_number)          =  -2.77154379D0 
     betap_mndod(atomic_number)          =  -1.80565019D0
     rho_core_mndod(atomic_number)       =   1.72118577D0
     GSS_mndod(atomic_number)            =   7.90443438D0
     GSP_mndod(atomic_number)            =   7.51570687D0
     GPP_mndod(atomic_number)            =   7.47999993D0
     GP2_mndod(atomic_number)            =   6.51866416D0
     HSP_mndod(atomic_number)            =   0.63674441D0
     alp_mndod(atomic_number)            =   1.42461329D0
     USS_mndod(atomic_number)            = -16.96970000D0
     UPP_mndod(atomic_number)            = -12.40096476D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.384108D0
     p_orb_exp_pm6(atomic_number) = 1.957413D0
     betas_pm6(atomic_number) = -11.613183D0
     betap_pm6(atomic_number) =   1.663178D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  6.677284D0
     GSP_pm6(atomic_number) =  5.953373D0
     GPP_pm6(atomic_number) = 18.729843D0
     GP2_pm6(atomic_number) =  9.917452D0
     HSP_pm6(atomic_number) =  0.825192D0
     USS_pm6(atomic_number) = -14.645792D0
     UPP_pm6(atomic_number) =  -9.318664D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.5253820D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END CADMIUM
 !-------------------

 !-------------------
 !INDIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 1, GSSC = 1, GPPC = 0, GSPC = 2, GP2C = 0, HSP = -1
   atomic_number = 49
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.0161160D0
     p_orb_exp_pm3(atomic_number) = 1.4453500D0
     betas_pm3(atomic_number) = -2.9933190D0
     betap_pm3(atomic_number) = -1.8289080D0
     FN1_pm3(1,atomic_number) = -0.3431380D0
     FN2_pm3(1,atomic_number) = 1.9940340D0
     FN3_pm3(1,atomic_number) = 1.6255160D0
     FN1_pm3(2,atomic_number) = -0.1095320D0
     FN2_pm3(2,atomic_number) = 5.6832170D0
     FN3_pm3(2,atomic_number) = 2.8670090D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 6.5549000D0
     GSP_pm3(atomic_number) = 8.2298730D0
     GPP_pm3(atomic_number) = 6.2992690D0
     GP2_pm3(atomic_number) = 4.9842110D0
     HSP_pm3(atomic_number) = 2.6314610D0
     alp_pm3(atomic_number) = 1.4183850D0
     USS_pm3(atomic_number) = -26.1762050D0
     UPP_pm3(atomic_number) = -20.0058220D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.023087D0
     p_orb_exp_pm6(atomic_number) = 2.106618D0
     betas_pm6(atomic_number) = -1.982376D0
     betap_pm6(atomic_number) = -3.330294D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  9.906091D0
     GSP_pm6(atomic_number) = 10.520060D0
     GPP_pm6(atomic_number) =  4.826006D0
     GP2_pm6(atomic_number) =  7.906563D0
     HSP_pm6(atomic_number) =  3.500299D0
     USS_pm6(atomic_number) = -28.339246D0
     UPP_pm6(atomic_number) = -23.373875D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.4183850D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END INDIUM
 !-------------------

 !-------------------
 !TIN
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 2, GSSC = 1, GPPC = -0.5, GSPC = 4, GP2C = 1.5, HSP = -2
   atomic_number = 50
  !MNDO
   ! Reference: M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART, J.AM.CHEM.SOC.,106 6771 (1984) (index = 13)
     mndo_ref_index(atomic_number) = 13
     element_supported_mndo(atomic_number) = .false.
     s_orb_exp_mndo(atomic_number) = 2.0803800D0
     p_orb_exp_mndo(atomic_number) = 1.9371060D0
     betas_mndo(atomic_number) = -3.2351470D0
     betap_mndo(atomic_number) = -4.2904160D0
     GSS_mndo(atomic_number) = 9.8000000D0
     GSP_mndo(atomic_number) = 8.3000000D0
     GPP_mndo(atomic_number) = 7.3000000D0
     GP2_mndo(atomic_number) = 6.5000000D0
     HSP_mndo(atomic_number) = 1.3000000D0
     alp_mndo(atomic_number) = 1.8008140D0
     USS_mndo(atomic_number) = -40.8518020D0
     UPP_mndo(atomic_number) = -28.5602490D0

  !MNDOD (same as MNDO)
   ! Reference: M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART, J.AM.CHEM.SOC.,106 6771 (1984) (index = 13)
     mndod_ref_index(atomic_number) = 13
     element_supported_mndod(atomic_number) = .false.
     s_orb_exp_mndod(atomic_number) = 2.0803800D0
     p_orb_exp_mndod(atomic_number) = 1.9371060D0
     betas_mndod(atomic_number) = -3.2351470D0
     betap_mndod(atomic_number) = -4.2904160D0
     GSS_mndod(atomic_number) = 9.8000000D0
     GSP_mndod(atomic_number) = 8.3000000D0
     GPP_mndod(atomic_number) = 7.3000000D0
     GP2_mndod(atomic_number) = 6.5000000D0
     HSP_mndod(atomic_number) = 1.3000000D0
     alp_mndod(atomic_number) = 1.8008140D0
     USS_mndod(atomic_number) = -40.8518020D0
     UPP_mndod(atomic_number) = -28.5602490D0

  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.3733280D0
     p_orb_exp_pm3(atomic_number) = 1.6382330D0
     betas_pm3(atomic_number) = -2.7858020D0
     betap_pm3(atomic_number) = -2.0059990D0
     FN1_pm3(1,atomic_number) =-0.1503530D0
     FN2_pm3(1,atomic_number) = 6.0056940D0
     FN3_pm3(1,atomic_number) = 1.7046420D0
     FN1_pm3(2,atomic_number) =-0.0444170D0
     FN2_pm3(2,atomic_number) = 2.2573810D0
     FN3_pm3(2,atomic_number) = 2.4698690D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 10.1900330D0
     GSP_pm3(atomic_number) = 7.2353270D0
     GPP_pm3(atomic_number) = 5.6738100D0
     GP2_pm3(atomic_number) = 5.1822140D0
     HSP_pm3(atomic_number) = 1.0331570D0
     alp_pm3(atomic_number) = 1.6996500D0
     USS_pm3(atomic_number) = -34.5501920D0
     UPP_pm3(atomic_number) = -25.8944190D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.383941D0
     p_orb_exp_pm6(atomic_number) = 2.057908D0
     betas_pm6(atomic_number) = -8.621087D0
     betap_pm6(atomic_number) = -4.989752D0
     FN1_pm6(1,atomic_number) =  -1.004587D0
     FN2_pm6(1,atomic_number) =   4.706252D0
     FN3_pm6(1,atomic_number) =   1.180218D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 8.269655D0
     GSP_pm6(atomic_number) = 5.013349D0
     GPP_pm6(atomic_number) = 6.584874D0
     GP2_pm6(atomic_number) = 5.855159D0
     HSP_pm6(atomic_number) = 0.531212D0
     USS_pm6(atomic_number) = -29.888217D0
     UPP_pm6(atomic_number) = -22.156954D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.6996500D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END TIN
 !-------------------

 !-------------------
 !ANTIMONY
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 3, GSSC = 1, GPPC = -1.5, GSPC = 6, GP2C = 4.5, HSP = -3
   atomic_number = 51
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 2.3430390D0
     p_orb_exp_pm3(atomic_number) = 1.8999920D0
     betas_pm3(atomic_number) = -14.7942170D0
     betap_pm3(atomic_number) = -2.8179480D0
     FN1_pm3(1,atomic_number) = 3.0020280D0
     FN2_pm3(1,atomic_number) = 6.0053420D0
     FN3_pm3(1,atomic_number) = 0.8530600D0
     FN1_pm3(2,atomic_number) =-0.0188920D0
     FN2_pm3(2,atomic_number) = 6.0114780D0
     FN3_pm3(2,atomic_number) = 2.7933110D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 9.2382770D0
     GSP_pm3(atomic_number) = 5.2776800D0
     GPP_pm3(atomic_number) = 6.3500000D0
     GP2_pm3(atomic_number) = 6.2500000D0
     HSP_pm3(atomic_number) = 2.4244640D0
     alp_pm3(atomic_number) = 2.0343010D0
     USS_pm3(atomic_number) = -56.4321960D0
     UPP_pm3(atomic_number) = -29.4349540D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  2.391178d0
     p_orb_exp_pm6(atomic_number) =  1.773006d0
     d_orb_exp_pm6(atomic_number) =  2.465590d0
     betas_pm6(atomic_number) =  -7.472322d0
     betap_pm6(atomic_number) =  -5.940750d0
     betad_pm6(atomic_number) =  -3.979108d0
     s_orb_exp_tail_pm6(atomic_number) =  5.993591d0
     p_orb_exp_tail_pm6(atomic_number) =  6.145086d0
     d_orb_exp_tail_pm6(atomic_number) =  5.704031d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  10.588832d0
     GSP_pm6(atomic_number) =   7.310023d0
     GPP_pm6(atomic_number) =   9.281609d0
     GP2_pm6(atomic_number) =   8.954081d0
     HSP_pm6(atomic_number) =   0.779112d0
     USS_pm6(atomic_number) = -41.688879d0
     UPP_pm6(atomic_number) = -39.541180d0
     UDD_pm6(atomic_number) =  -6.581663d0
     ! alp_pm6 is taken from PM3
     alp_pm6(atomic_number) = 2.0343010D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END ANTIMONY
 !-------------------

 !-------------------
 !TELLURIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 4, GSSC = 1, GPPC = -0.5, GSPC = 8, GP2C = 6.5, HSP = -4
   atomic_number = 52
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 4.1654920D0
     p_orb_exp_pm3(atomic_number) = 1.6475550D0
     betas_pm3(atomic_number) = -2.6651460D0
     betap_pm3(atomic_number) = -3.8954300D0
     FN1_pm3(1,atomic_number) = 0.0333910D0
     FN2_pm3(1,atomic_number) = 5.9563790D0
     FN3_pm3(1,atomic_number) = 2.2775750D0
     FN1_pm3(2,atomic_number) =-1.9218670D0
     FN2_pm3(2,atomic_number) = 4.9732190D0
     FN3_pm3(2,atomic_number) = 0.5242430D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 10.2550730D0
     GSP_pm3(atomic_number) = 8.1691450D0
     GPP_pm3(atomic_number) = 7.7775920D0
     GP2_pm3(atomic_number) = 7.7551210D0
     HSP_pm3(atomic_number) = 3.7724620D0
     alp_pm3(atomic_number) = 2.4850190D0
     USS_pm3(atomic_number) = -44.9380360D0
     UPP_pm3(atomic_number) = -46.3140990D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  2.769862d0
     p_orb_exp_pm6(atomic_number) =  1.731319d0
     betas_pm6(atomic_number) =  -70.001062d0
     betap_pm6(atomic_number) =   -6.151642d0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =    7.030626d0
     GSP_pm6(atomic_number) =   12.601389d0
     GPP_pm6(atomic_number) =    7.883479d0
     GP2_pm6(atomic_number) =    6.973163d0
     HSP_pm6(atomic_number) =    5.000826d0
     USS_pm6(atomic_number) = -114.733316d0
     UPP_pm6(atomic_number) =  -50.096389d0
     ! alp_pm6 is taken from PM3
     alp_pm6(atomic_number) = 2.4850190D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END TELLURIUM
 !-------------------

 !-------------------
 !IODINE
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 5, GSSC = 1, GPPC = 0, GSPC = 10, GP2C = 10, HSP = -5
   atomic_number = 53
  !MNDO
   ! Reference: M.J.S.DEWAR, E.F. HEALY, J.J.P. STEWART, J.COMP.CHEM., 5,358,(1984) (index = 14)
     mndo_ref_index(atomic_number) = 14
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.2729610D0
     p_orb_exp_mndo(atomic_number) = 2.1694980D0
     betas_mndo(atomic_number) = -7.4144510D0
     betap_mndo(atomic_number) = -6.1967810D0
     GSS_mndo(atomic_number) = 15.04044855D0
     GSP_mndo(atomic_number) = 13.05655798D0
     GPP_mndo(atomic_number) = 11.14778369D0
     GP2_mndo(atomic_number) = 9.91409071D0
     HSP_mndo(atomic_number) = 2.45638202D0
     alp_mndo(atomic_number) = 2.2073200D0
     USS_mndo(atomic_number) = -100.0030538D0
     UPP_mndo(atomic_number) = -74.6114692D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =   2.75654324D0
     p_orb_exp_mndod(atomic_number)      =   2.25307954D0
     d_orb_exp_mndod(atomic_number)      =   1.50233509D0   
     betas_mndod(atomic_number)          = -10.69948666D0
     betap_mndod(atomic_number)          =  -4.94117776D0
     betad_mndod(atomic_number)          =  -2.35046098D0  
     s_orb_exp_tail_mndod(atomic_number) =   2.67241100D0
     p_orb_exp_tail_mndod(atomic_number) =   1.57229871D0
     d_orb_exp_tail_mndod(atomic_number) =   1.25884802D0
     GSS_mndod(atomic_number)            =  11.98078196D0
     GSP_mndod(atomic_number)            =   7.85590192D0
     GPP_mndod(atomic_number)            =   7.70937227D0
     GDD_mndod(atomic_number)            =   6.09729900D0
     GP2_mndod(atomic_number)            =   6.71855729D0
     HSP_mndod(atomic_number)            =   2.07147462D0
     alp_mndod(atomic_number)            =   1.90617441D0
     USS_mndod(atomic_number)            = -62.76535256D0
     UPP_mndod(atomic_number)            = -50.29211568D0
     UDD_mndod(atomic_number)            = -12.24830501D0

  !AM1
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1_ref_index(atomic_number) = 18
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 2.1028580D0
     p_orb_exp_am1(atomic_number) = 2.1611530D0
     betas_am1(atomic_number) = -8.4433270D0
     betap_am1(atomic_number) = -6.3234050D0
     FN1_am1(1,atomic_number) = 0.0043610D0
     FN2_am1(1,atomic_number) = 2.3000000D0
     FN3_am1(1,atomic_number) = 1.8000000D0
     FN1_am1(2,atomic_number) = 0.0157060D0
     FN2_am1(2,atomic_number) = 3.0000000D0
     FN3_am1(2,atomic_number) = 2.2400000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 2
     GSS_am1(atomic_number) = 15.0404486D0
     GSP_am1(atomic_number) = 13.0565580D0
     GPP_am1(atomic_number) = 11.1477837D0
     GP2_am1(atomic_number) = 9.9140907D0
     HSP_am1(atomic_number) = 2.4563820D0
     alp_am1(atomic_number) = 2.2994240D0
     USS_am1(atomic_number) = -103.5896630D0
     UPP_am1(atomic_number) = -74.4299970D0
     
  !AM1D
   ! Reference: M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988). (index = 18)
     am1d_ref_index(atomic_number) = 18
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 2.1028580D0
     p_orb_exp_am1d(atomic_number) = 2.1611530D0
     betas_am1d(atomic_number) = -8.4433270D0
     betap_am1d(atomic_number) = -6.3234050D0
     FN1_am1d(1,atomic_number) = 0.0043610D0
     FN2_am1d(1,atomic_number) = 2.3000000D0
     FN3_am1d(1,atomic_number) = 1.8000000D0
     FN1_am1d(2,atomic_number) = 0.0157060D0
     FN2_am1d(2,atomic_number) = 3.0000000D0
     FN3_am1d(2,atomic_number) = 2.2400000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 2
     GSS_am1d(atomic_number) = 15.0404486D0
     GSP_am1d(atomic_number) = 13.0565580D0
     GPP_am1d(atomic_number) = 11.1477837D0
     GP2_am1d(atomic_number) = 9.9140907D0
     HSP_am1d(atomic_number) = 2.4563820D0
     alp_am1d(atomic_number) = 2.2994240D0
     USS_am1d(atomic_number) = -103.5896630D0
     UPP_am1d(atomic_number) = -74.4299970D0     

  !PM3
   ! Reference: J. J. P. STEWART, J. COMP. CHEM.10, 209 (1989).
     pm3_ref_index(atomic_number) = 26
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 7.0010130D0
     p_orb_exp_pm3(atomic_number) = 2.4543540D0
     betas_pm3(atomic_number) = -14.4942340D0
     betap_pm3(atomic_number) = -5.8947030D0
     FN1_pm3(1,atomic_number) =-0.1314810D0
     FN2_pm3(1,atomic_number) = 5.2064170D0
     FN3_pm3(1,atomic_number) = 1.7488240D0
     FN1_pm3(2,atomic_number) =-0.0368970D0
     FN2_pm3(2,atomic_number) = 6.0101170D0
     FN3_pm3(2,atomic_number) = 2.7103730D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 13.6319430D0
     GSP_pm3(atomic_number) = 14.9904060D0
     GPP_pm3(atomic_number) = 7.2883300D0
     GP2_pm3(atomic_number) = 5.9664070D0
     HSP_pm3(atomic_number) = 2.6300350D0
     alp_pm3(atomic_number) = 1.9901850D0
     USS_pm3(atomic_number) = -96.4540370D0
     UPP_pm3(atomic_number) = -61.0915820D0

  !PDDG/PM3
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgpm3_ref_index(atomic_number) = 30
     element_supported_pddgpm3(atomic_number) = .true.
     elec_eng_pddgpm3(atomic_number) = -291.537869D0
     s_orb_exp_pddgpm3(atomic_number) = 5.062801D0
     p_orb_exp_pddgpm3(atomic_number) = 2.417757D0
     betas_pddgpm3(atomic_number) = -16.592621D0
     betap_pddgpm3(atomic_number) = -6.599816D0
     FN1_pddgpm3(1,atomic_number) = -0.136003D0
     FN2_pddgpm3(1,atomic_number) = 3.852912D0
     FN3_pddgpm3(1,atomic_number) = 1.697455D0
     FN1_pddgpm3(2,atomic_number) = -0.037287D0
     FN2_pddgpm3(2,atomic_number) = 5.229264D0
     FN3_pddgpm3(2,atomic_number) = 2.768669D0
     FN1_pddgpm3(3,atomic_number) = 0.0d0
     FN2_pddgpm3(3,atomic_number) = 0.0d0
     FN3_pddgpm3(3,atomic_number) = 0.0d0
     FN1_pddgpm3(4,atomic_number) = 0.0d0
     FN2_pddgpm3(4,atomic_number) = 0.0d0
     FN3_pddgpm3(4,atomic_number) = 0.0d0
     NUM_FN_pddgpm3(atomic_number) = 2
     GSS_pddgpm3(atomic_number) = 13.6319430D0
     GSP_pddgpm3(atomic_number) = 14.9904060D0
     GPP_pddgpm3(atomic_number) = 7.2883300D0
     GP2_pddgpm3(atomic_number) = 5.9664070D0
     HSP_pddgpm3(atomic_number) = 2.6300350D0
     alp_pddgpm3(atomic_number) = 1.978170D0
     USS_pddgpm3(atomic_number) = -97.664174D0
     UPP_pddgpm3(atomic_number) = -61.167137D0
     PDDGC1_pm3(atomic_number) = 0.012901D0
     PDDGC2_pm3(atomic_number) = -0.012825D0
     PDDGE1_pm3(atomic_number) = 1.994299D0
     PDDGE2_pm3(atomic_number) = 2.263417D0

  !PDDG/MNDO
   ! Reference Tubert-Brohman, Werneck, Repasky and Jorgensen, 2003, J Comp Chem, 25: 138-150 (Index 30)
     pddgpm3_ref_index(atomic_number) = 30
     element_supported_pddgmndo(atomic_number) = .true.
     elec_eng_pddgmndo(atomic_number) = -356.076398D0
     s_orb_exp_pddgmndo(atomic_number) = 2.718404D0
     p_orb_exp_pddgmndo(atomic_number) = 2.461813D0
     betas_pddgmndo(atomic_number) = -6.698375D0
     betap_pddgmndo(atomic_number) = -5.693814D0
     GSS_pddgmndo(atomic_number) = 15.04044855D0
     GSP_pddgmndo(atomic_number) = 13.05655798D0
     GPP_pddgmndo(atomic_number) = 11.14778369D0
     GP2_pddgmndo(atomic_number) = 9.91409071D0
     HSP_pddgmndo(atomic_number) = 2.45638202D0
     alp_pddgmndo(atomic_number) = 2.242446D0
     USS_pddgmndo(atomic_number) = -106.588422D0
     UPP_pddgmndo(atomic_number) = -75.282605D0
     PDDGC1_mndo(atomic_number) = 0.009616D0
     PDDGC2_mndo(atomic_number) = -0.007505D0
     PDDGE1_mndo(atomic_number) = 2.572332D0
     PDDGE2_mndo(atomic_number) = 2.936456D0

  !RM1
   ! Reference: G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)
     rm1_ref_index(atomic_number) = 33
     element_supported_rm1(atomic_number) = .true.
     s_orb_exp_rm1(atomic_number) = 2.53003753d0
     p_orb_exp_rm1(atomic_number) = 2.31738678d0
     betas_rm1(atomic_number) = -4.19316149d0
     betap_rm1(atomic_number) = -4.40038412d0
     FN1_rm1(1,atomic_number) = -0.08147724d0
     FN2_rm1(1,atomic_number) = 1.56065072d0
     FN3_rm1(1,atomic_number) = 2.00002063d0
     FN1_rm1(2,atomic_number) = 0.05913991d0
     FN2_rm1(2,atomic_number) = 5.76111270d0
     FN3_rm1(2,atomic_number) = 2.20488800d0
     FN1_rm1(3,atomic_number) = 0.0d0
     FN2_rm1(3,atomic_number) = 0.0d0
     FN3_rm1(3,atomic_number) = 0.0d0
     FN1_rm1(4,atomic_number) = 0.0d0
     FN2_rm1(4,atomic_number) = 0.0d0
     FN3_rm1(4,atomic_number) = 0.0d0
     NUM_FN_rm1(atomic_number) = 2
     GSS_rm1(atomic_number) = 19.99974131d0
     GSP_rm1(atomic_number) = 7.68957672d0
     GPP_rm1(atomic_number) = 7.30488343d0
     GP2_rm1(atomic_number) = 6.85424614d0
     HSP_rm1(atomic_number) = 1.41602940d0
     alp_rm1(atomic_number) = 2.14157092d0
     USS_rm1(atomic_number) = -74.89997837d0
     UPP_rm1(atomic_number) = -51.41023805d0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) =  4.498653d0
     p_orb_exp_pm6(atomic_number) =  1.917072d0
     d_orb_exp_pm6(atomic_number) =  1.875175d0
     betas_pm6(atomic_number) =  -30.522481d0
     betap_pm6(atomic_number) =   -5.942120d0
     betad_pm6(atomic_number) =   -7.676107d0
     s_orb_exp_tail_pm6(atomic_number) =  9.135244d0
     p_orb_exp_tail_pm6(atomic_number) =  6.888191d0
     d_orb_exp_tail_pm6(atomic_number) =  3.791523d0
     FN1_pm6(1,atomic_number) = -0.035519d0
     FN2_pm6(1,atomic_number) =  1.744389d0
     FN3_pm6(1,atomic_number) =  1.223844d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) =   7.234759d0
     GSP_pm6(atomic_number) =   9.154406d0
     GPP_pm6(atomic_number) =   9.877466d0
     GP2_pm6(atomic_number) =   8.035916d0
     HSP_pm6(atomic_number) =   5.004215d0
     USS_pm6(atomic_number) = -59.973232d0
     UPP_pm6(atomic_number) = -56.459835d0
     UDD_pm6(atomic_number) = -28.822603d0
     ! alp_pm6 is taken from PM3
     alp_pm6(atomic_number) = 1.9901850D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END IODINE
 !-------------------

 !-------------------
 ! XENON
 !-------------------
     atomic_number = 54
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.759787D0
     p_orb_exp_pm6(atomic_number) = 1.977446D0
     betas_pm6(atomic_number) =  -3.980622D0
     betap_pm6(atomic_number) = -38.822792D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 20.000252D0
     GSP_pm6(atomic_number) =  4.175902D0
     GPP_pm6(atomic_number) =  2.305787D0
     GP2_pm6(atomic_number) =  4.063220D0
     HSP_pm6(atomic_number) =  4.418843D0
     USS_pm6(atomic_number) =  -18.270227D0
     UPP_pm6(atomic_number) = -167.163063D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 3.0d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END XENON
 !-------------------

 !-------------------
 ! CESIUM
 !-------------------
     atomic_number = 55
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 5.956008D0
     p_orb_exp_pm6(atomic_number) = 1.619485D0
     betas_pm6(atomic_number) =  2.287838D0
     betap_pm6(atomic_number) = -5.908071D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  6.464751D0
     GSP_pm6(atomic_number) =  4.004501D0
     GPP_pm6(atomic_number) = 13.775390D0
     GP2_pm6(atomic_number) = 12.912537D0
     HSP_pm6(atomic_number) =  1.026928D0
     USS_pm6(atomic_number) = -3.748609D0
     UPP_pm6(atomic_number) = -2.348109D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.2d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END CESIUM
 !-------------------

 !-------------------
 ! BARIUM
 !-------------------
     atomic_number = 56
  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 1.395379D0
     p_orb_exp_pm6(atomic_number) = 1.430139D0
     betas_pm6(atomic_number) = 10.003125D0
     betap_pm6(atomic_number) = -6.335160D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 3.600823D0
     GSP_pm6(atomic_number) = 4.740579D0
     GPP_pm6(atomic_number) = 3.345166D0
     GP2_pm6(atomic_number) = 3.142783D0
     HSP_pm6(atomic_number) = 0.929429D0
     USS_pm6(atomic_number) = -9.306985D0
     UPP_pm6(atomic_number) = -8.826713D0
     ! alp_pm6 is a guess by AWG
     alp_pm6(atomic_number) = 1.3d0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END BARIUM
 !-------------------

 !-------------------
 !MERCURY
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 0, GSSC = 1, GPPC = 0, GSPC = 0, GP2C = 0, HSP = 0
   atomic_number = 80
  !MNDO
   ! Reference: M.J.S.DEWAR,  ET. AL. ORGANOMETALLICS 4, 1964, (1985) (index = 15)
     mndo_ref_index(atomic_number) = 15
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.2181840D0
     p_orb_exp_mndo(atomic_number) = 2.0650380D0
     betas_mndo(atomic_number) = -0.4045250D0
     betap_mndo(atomic_number) = -6.2066830D0
     GSS_mndo(atomic_number) = 10.8000000D0
     GSP_mndo(atomic_number) = 9.3000000D0
     GPP_mndo(atomic_number) = 14.3000000D0
     GP2_mndo(atomic_number) = 13.5000000D0
     HSP_mndo(atomic_number) = 1.3000000D0
     alp_mndo(atomic_number) = 1.3356410D0
     USS_mndo(atomic_number) = -19.8095740D0
     UPP_mndo(atomic_number) = -13.1025300D0

  !MNDOD
   ! Reference: W. Thiel and A. Voityuk,J. Phys. CHEM., 100 616-626 1996 (index = 101)
     mndod_ref_index(atomic_number) = 101
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number)      =    2.33310757D0
     p_orb_exp_mndod(atomic_number)      =    1.70831069D0
     s_orb_exp_tail_mndod(atomic_number) =    2.18600011D0
     p_orb_exp_tail_mndod(atomic_number) =    1.70500461D0
     betas_mndod(atomic_number)          =   -2.21872239D0 
     betap_mndod(atomic_number)          =   -2.90978573D0
     rho_core_mndod(atomic_number)       =    1.63607185D0
     GSS_mndod(atomic_number)            =    8.31564948D0
     GSP_mndod(atomic_number)            =    8.21217300D0
     GPP_mndod(atomic_number)            =    7.11525878D0
     GP2_mndod(atomic_number)            =    6.17124983D0
     HSP_mndod(atomic_number)            =    0.83594100D0
     alp_mndod(atomic_number)            =    1.38224172D0
     USS_mndod(atomic_number)            =  -18.81564903D0
     UPP_mndod(atomic_number)            =  -13.39711352D0

  !AM1
   ! Reference: M.J.S.Dewar and C.Jie, Organometallics 8, 1547, (1989)
     am1_ref_index(atomic_number) = 25
     element_supported_am1(atomic_number) = .true.
     s_orb_exp_am1(atomic_number) = 2.0364130D0
     p_orb_exp_am1(atomic_number) = 1.9557660D0
     betas_am1(atomic_number) = -0.9086570D0
     betap_am1(atomic_number) = -4.9093840D0
     FN1_am1(1,atomic_number) = 0.0000000D0
     FN2_am1(1,atomic_number) = 0.0000000D0
     FN3_am1(1,atomic_number) = 0.0000000D0
     FN1_am1(2,atomic_number) = 0.0000000D0
     FN2_am1(2,atomic_number) = 0.0000000D0
     FN3_am1(2,atomic_number) = 0.0000000D0
     FN1_am1(3,atomic_number) = 0.0000000D0
     FN2_am1(3,atomic_number) = 0.0000000D0
     FN3_am1(3,atomic_number) = 0.0000000D0
     FN1_am1(4,atomic_number) = 0.0000000D0
     FN2_am1(4,atomic_number) = 0.0000000D0
     FN3_am1(4,atomic_number) = 0.0000000D0
     NUM_FN_am1(atomic_number) = 0
     GSS_am1(atomic_number) = 10.8000000D0
     GSP_am1(atomic_number) = 9.3000000D0
     GPP_am1(atomic_number) = 14.3000000D0
     GP2_am1(atomic_number) = 13.5000000D0
     HSP_am1(atomic_number) = 1.3000000D0
     alp_am1(atomic_number) = 1.4847340D0
     USS_am1(atomic_number) = -19.9415780D0
     UPP_am1(atomic_number) = -11.1108700D0

  !AM1D
   ! Reference: M.J.S.Dewar and C.Jie, Organometallics 8, 1547, (1989)
     am1d_ref_index(atomic_number) = 25
     element_supported_am1d(atomic_number) = .true.
     s_orb_exp_am1d(atomic_number) = 2.0364130D0
     p_orb_exp_am1d(atomic_number) = 1.9557660D0
     betas_am1d(atomic_number) = -0.9086570D0
     betap_am1d(atomic_number) = -4.9093840D0
     FN1_am1d(1,atomic_number) = 0.0000000D0
     FN2_am1d(1,atomic_number) = 0.0000000D0
     FN3_am1d(1,atomic_number) = 0.0000000D0
     FN1_am1d(2,atomic_number) = 0.0000000D0
     FN2_am1d(2,atomic_number) = 0.0000000D0
     FN3_am1d(2,atomic_number) = 0.0000000D0
     FN1_am1d(3,atomic_number) = 0.0000000D0
     FN2_am1d(3,atomic_number) = 0.0000000D0
     FN3_am1d(3,atomic_number) = 0.0000000D0
     FN1_am1d(4,atomic_number) = 0.0000000D0
     FN2_am1d(4,atomic_number) = 0.0000000D0
     FN3_am1d(4,atomic_number) = 0.0000000D0
     NUM_FN_am1d(atomic_number) = 0
     GSS_am1d(atomic_number) = 10.8000000D0
     GSP_am1d(atomic_number) = 9.3000000D0
     GPP_am1d(atomic_number) = 14.3000000D0
     GP2_am1d(atomic_number) = 13.5000000D0
     HSP_am1d(atomic_number) = 1.3000000D0
     alp_am1d(atomic_number) = 1.4847340D0
     USS_am1d(atomic_number) = -19.9415780D0
     UPP_am1d(atomic_number) = -11.1108700D0



  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 1.4768850D0
     p_orb_exp_pm3(atomic_number) = 2.4799510D0
     betas_pm3(atomic_number) = -3.1013650D0
     betap_pm3(atomic_number) = -3.4640310D0
     FN1_pm3(1,atomic_number) = 1.0827200D0
     FN2_pm3(1,atomic_number) = 6.4965980D0
     FN3_pm3(1,atomic_number) = 1.1951460D0
     FN1_pm3(2,atomic_number) =-0.0965530D0
     FN2_pm3(2,atomic_number) = 3.9262810D0
     FN3_pm3(2,atomic_number) = 2.6271600D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 6.6247200D0
     GSP_pm3(atomic_number) = 10.6392970D0
     GPP_pm3(atomic_number) = 14.7092830D0
     GP2_pm3(atomic_number) = 16.0007400D0
     HSP_pm3(atomic_number) = 2.0363110D0
     alp_pm3(atomic_number) = 1.5293770D0
     USS_pm3(atomic_number) = -17.7622290D0
     UPP_pm3(atomic_number) = -18.3307510D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.104896D0
     p_orb_exp_pm6(atomic_number) = 1.516293D0
     betas_pm6(atomic_number) = -3.045239D0
     betap_pm6(atomic_number) = -5.693556D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  6.372822D0
     GSP_pm6(atomic_number) = 10.143176D0
     GPP_pm6(atomic_number) = 10.397393D0
     GP2_pm6(atomic_number) = 14.794056D0
     HSP_pm6(atomic_number) =  0.926128D0
     USS_pm6(atomic_number) = -17.608732D0
     UPP_pm6(atomic_number) = -18.369417D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.5293770D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END MERCURY
 !-------------------

 !-------------------
 !THALLIUM
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 1, GSSC = 1, GPPC = 0, GSPC = 2, GP2C = 0, HSP = -1
   atomic_number = 81
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 6.8679210D0
     p_orb_exp_pm3(atomic_number) = 1.9694450D0
     betas_pm3(atomic_number) = -1.0844950D0
     betap_pm3(atomic_number) = -7.9467990D0
     FN1_pm3(1,atomic_number) =-1.3613990D0
     FN2_pm3(1,atomic_number) = 3.5572260D0
     FN3_pm3(1,atomic_number) = 1.0928020D0
     FN1_pm3(2,atomic_number) =-0.0454010D0
     FN2_pm3(2,atomic_number) = 2.3069950D0
     FN3_pm3(2,atomic_number) = 2.9650290D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 10.4604120D0
     GSP_pm3(atomic_number) = 11.2238830D0
     GPP_pm3(atomic_number) = 4.9927850D0
     GP2_pm3(atomic_number) = 8.9627270D0
     HSP_pm3(atomic_number) = 2.5304060D0
     alp_pm3(atomic_number) = 1.3409510D0
     USS_pm3(atomic_number) = -30.0531700D0
     UPP_pm3(atomic_number) = -26.9206370D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 3.335883D0
     p_orb_exp_pm6(atomic_number) = 1.766141D0
     betas_pm6(atomic_number) = -7.230170D0
     betap_pm6(atomic_number) = -7.575544D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) =  5.015118D0
     GSP_pm6(atomic_number) = 13.932049D0
     GPP_pm6(atomic_number) = 10.495551D0
     GP2_pm6(atomic_number) = 10.526198D0
     HSP_pm6(atomic_number) =  0.293760D0
     USS_pm6(atomic_number) = -29.518621D0
     UPP_pm6(atomic_number) = -29.826907D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.3409510D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END THALLIUM
 !-------------------

 !-------------------
 !LEAD
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 2, GSSC = 1, GPPC = -0.5, GSPC = 4, GP2C = 1.5, HSP = -2
   atomic_number = 82
  !MNDO
   ! Reference: M.J.S.DEWAR, ET.AL ORGANOMETALLICS 4 1973-1980 (1985) (index = 16)
     mndo_ref_index(atomic_number) = 16
     element_supported_mndo(atomic_number) = .true.
     s_orb_exp_mndo(atomic_number) = 2.4982860D0
     p_orb_exp_mndo(atomic_number) = 2.0820710D0
     betas_mndo(atomic_number) = -8.0423870D0
     betap_mndo(atomic_number) = -3.0000000D0
     GSS_mndo(atomic_number) = 9.8000000D0
     GSP_mndo(atomic_number) = 8.3000000D0
     GPP_mndo(atomic_number) = 7.3000000D0
     GP2_mndo(atomic_number) = 6.5000000D0
     HSP_mndo(atomic_number) = 1.3000000D0
     alp_mndo(atomic_number) = 1.7283330D0
     USS_mndo(atomic_number) = -47.3196920D0
     UPP_mndo(atomic_number) = -28.8475600D0

  !MNDOD (same as MNDO)
   ! Reference: M.J.S.DEWAR, ET.AL ORGANOMETALLICS 4 1973-1980 (1985) (index = 16)
     mndod_ref_index(atomic_number) = 16
     element_supported_mndod(atomic_number) = .true.
     s_orb_exp_mndod(atomic_number) = 2.4982860D0
     p_orb_exp_mndod(atomic_number) = 2.0820710D0
     betas_mndod(atomic_number) = -8.0423870D0
     betap_mndod(atomic_number) = -3.0000000D0
     GSS_mndod(atomic_number) = 9.8000000D0
     GSP_mndod(atomic_number) = 8.3000000D0
     GPP_mndod(atomic_number) = 7.3000000D0
     GP2_mndod(atomic_number) = 6.5000000D0
     HSP_mndod(atomic_number) = 1.3000000D0
     alp_mndod(atomic_number) = 1.7283330D0
     USS_mndod(atomic_number) = -47.3196920D0
     UPP_mndod(atomic_number) = -28.8475600D0

  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 3.1412890D0
     p_orb_exp_pm3(atomic_number) = 1.8924180D0
     betas_pm3(atomic_number) = -6.1260240D0
     betap_pm3(atomic_number) = -1.3954300D0
     FN1_pm3(1,atomic_number) =-0.1225760D0
     FN2_pm3(1,atomic_number) = 6.0030620D0
     FN3_pm3(1,atomic_number) = 1.9015970D0
     FN1_pm3(2,atomic_number) =-0.0566480D0
     FN2_pm3(2,atomic_number) = 4.7437050D0
     FN3_pm3(2,atomic_number) = 2.8618790D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 7.0119920D0
     GSP_pm3(atomic_number) = 6.7937820D0
     GPP_pm3(atomic_number) = 5.1837800D0
     GP2_pm3(atomic_number) = 5.0456510D0
     HSP_pm3(atomic_number) = 1.5663020D0
     alp_pm3(atomic_number) = 1.6200450D0
     USS_pm3(atomic_number) = -30.3227560D0
     UPP_pm3(atomic_number) = -24.4258340D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 2.368901D0
     p_orb_exp_pm6(atomic_number) = 1.685246D0
     betas_pm6(atomic_number) = -8.323792D0
     betap_pm6(atomic_number) = -2.237891D0
     FN1_pm6(1,atomic_number) = -0.239463D0
     FN2_pm6(1,atomic_number) =  5.444338D0
     FN3_pm6(1,atomic_number) =  1.613682D0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 1
     GSS_pm6(atomic_number) = 5.254128D0
     GSP_pm6(atomic_number) = 7.061016D0
     GPP_pm6(atomic_number) = 6.818551D0
     GP2_pm6(atomic_number) = 5.603019D0
     HSP_pm6(atomic_number) = 1.018819D0
     USS_pm6(atomic_number) = -35.038145D0
     UPP_pm6(atomic_number) = -25.413401D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.6200450D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END LEAD
 !-------------------

 !-------------------
 !BISMUTH
 !-------------------
 !Notes for elec eng: IOS = 2, IOP = 3, GSSC = 1, GPPC = -1.5, GSPC = 6, GP2C = 4.5, HSP = -3
   atomic_number = 83
  !PM3
   ! Reference: J.J.P.STEWART, JCC, 3, 320 (1991) (Index 27)
     pm3_ref_index(atomic_number) = 27
     element_supported_pm3(atomic_number) = .true.
     s_orb_exp_pm3(atomic_number) = 4.9164510D0
     p_orb_exp_pm3(atomic_number) = 1.9349350D0
     betas_pm3(atomic_number) = -5.6072830D0
     betap_pm3(atomic_number) = -5.8001520D0
     FN1_pm3(1,atomic_number) = 2.5816930D0
     FN2_pm3(1,atomic_number) = 5.0940220D0
     FN3_pm3(1,atomic_number) = 0.4997870D0
     FN1_pm3(2,atomic_number) = 0.0603200D0
     FN2_pm3(2,atomic_number) = 6.0015380D0
     FN3_pm3(2,atomic_number) = 2.4279700D0
     FN1_pm3(3,atomic_number) = 0.0d0
     FN2_pm3(3,atomic_number) = 0.0d0
     FN3_pm3(3,atomic_number) = 0.0d0
     FN1_pm3(4,atomic_number) = 0.0d0
     FN2_pm3(4,atomic_number) = 0.0d0
     FN3_pm3(4,atomic_number) = 0.0d0
     NUM_FN_pm3(atomic_number) = 2
     GSS_pm3(atomic_number) = 4.9894800D0
     GSP_pm3(atomic_number) = 6.1033080D0
     GPP_pm3(atomic_number) = 8.6960070D0
     GP2_pm3(atomic_number) = 8.3354470D0
     HSP_pm3(atomic_number) = 0.5991220D0
     alp_pm3(atomic_number) = 1.8574310D0
     USS_pm3(atomic_number) = -33.4959380D0
     UPP_pm3(atomic_number) = -35.5210260D0

  !PM6
   ! Reference J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007). (Index 35)
     pm6_ref_index(atomic_number) = 35
     element_supported_pm6(atomic_number) = .true.
     s_orb_exp_pm6(atomic_number) = 3.702377D0
     p_orb_exp_pm6(atomic_number) = 1.872327D0
     betas_pm6(atomic_number) = -34.951578D0
     betap_pm6(atomic_number) =  -7.359060D0
     FN1_pm6(1,atomic_number) = 0.0d0
     FN2_pm6(1,atomic_number) = 0.0d0
     FN3_pm6(1,atomic_number) = 0.0d0
     FN1_pm6(2,atomic_number) = 0.0d0
     FN2_pm6(2,atomic_number) = 0.0d0
     FN3_pm6(2,atomic_number) = 0.0d0
     FN1_pm6(3,atomic_number) = 0.0d0
     FN2_pm6(3,atomic_number) = 0.0d0
     FN3_pm6(3,atomic_number) = 0.0d0
     FN1_pm6(4,atomic_number) = 0.0d0
     FN2_pm6(4,atomic_number) = 0.0d0
     FN3_pm6(4,atomic_number) = 0.0d0
     NUM_FN_pm6(atomic_number) = 0
     GSS_pm6(atomic_number) = 5.851803D0
     GSP_pm6(atomic_number) = 6.790583D0
     GPP_pm6(atomic_number) = 8.389442D0
     GP2_pm6(atomic_number) = 7.724219D0
     HSP_pm6(atomic_number) = 0.295606D0
     USS_pm6(atomic_number) = -42.409177D0
     UPP_pm6(atomic_number) = -36.393746D0
     ! alp_pm6 contains the PM3 parameter value
     alp_pm6(atomic_number) = 1.8574310D0
     !For pairwise core core terms see section 3 below.

 !-------------------
 !END BISMUTH
 !-------------------

!----------------------------------------------------------------------
!END OF SECTION 2 - MNDO, AM1, RM1, PM3, PM6 and PDDG PARAMS BY ELEMENT
!----------------------------------------------------------------------

!-----------------------------------------
!SECTION 3 - PM6 pair wise core core terms
!-----------------------------------------

!In PM6 the core core terms are dealt with on an element by element
!pairwise basis. There are two values a coefficient for the exponent
!and a coefficient for the multiplier. We use the same terminology
!as the PM6 paper.

! H
    alpab_pm6( 1, 1) =   3.540942d0  !    Hydrogen -     Hydrogen
    xab_pm6  ( 1, 1) =   2.243587d0  !    Hydrogen -     Hydrogen
! He
    alpab_pm6( 2, 1) =   2.989881d0  !      Helium -     Hydrogen
    xab_pm6  ( 2, 1) =   2.371199d0  !      Helium -     Hydrogen
    alpab_pm6( 2, 2) =   3.783559d0  !      Helium -       Helium
    xab_pm6  ( 2, 2) =   3.450900d0  !      Helium -       Helium
! Li
    alpab_pm6( 3, 1) =   2.136265d0  !     Lithium -     Hydrogen
    xab_pm6  ( 3, 1) =   2.191985d0  !     Lithium -     Hydrogen
    alpab_pm6( 3, 2) =   3.112403d0  !     Lithium -       Helium
    xab_pm6  ( 3, 2) =   9.273676d0  !     Lithium -       Helium
    alpab_pm6( 3, 3) =   4.714674d0  !     Lithium -      Lithium
    xab_pm6  ( 3, 3) =  16.116384d0  !     Lithium -      Lithium
! Be
    alpab_pm6( 4, 1) =   2.475418d0  !   Beryllium -     Hydrogen
    xab_pm6  ( 4, 1) =   2.562831d0  !   Beryllium -     Hydrogen
    alpab_pm6( 4, 2) =   3.306702d0  !   Beryllium -       Helium
    xab_pm6  ( 4, 2) =  12.544878d0  !   Beryllium -       Helium
    alpab_pm6( 4, 3) =   2.236728d0  !   Beryllium -      Lithium
    xab_pm6  ( 4, 3) =   3.287165d0  !   Beryllium -      Lithium
    alpab_pm6( 4, 4) =   1.499907d0  !   Beryllium -    Beryllium
    xab_pm6  ( 4, 4) =   0.238633d0  !   Beryllium -    Beryllium
! B
    alpab_pm6( 5, 1) =   2.615231d0  !       Boron -     Hydrogen
    xab_pm6  ( 5, 1) =   1.321394d0  !       Boron -     Hydrogen
    alpab_pm6( 5, 2) =   3.163140d0  !       Boron -       Helium
    xab_pm6  ( 5, 2) =   1.974170d0  !       Boron -       Helium
    alpab_pm6( 5, 3) =   3.759397d0  !       Boron -      Lithium
    xab_pm6  ( 5, 3) =   7.886018d0  !       Boron -      Lithium
    alpab_pm6( 5, 4) =   1.888998d0  !       Boron -    Beryllium
    xab_pm6  ( 5, 4) =   1.151792d0  !       Boron -    Beryllium
    alpab_pm6( 5, 5) =   3.318624d0  !       Boron -        Boron
    xab_pm6  ( 5, 5) =   3.593619d0  !       Boron -        Boron
! C
    alpab_pm6( 6, 1) =   1.027806d0  !      Carbon -     Hydrogen
    xab_pm6  ( 6, 1) =   0.216506d0  !      Carbon -     Hydrogen
    alpab_pm6( 6, 2) =   3.042705d0  !      Carbon -       Helium
    xab_pm6  ( 6, 2) =   3.213971d0  !      Carbon -       Helium
    alpab_pm6( 6, 3) =   3.241874d0  !      Carbon -      Lithium
    xab_pm6  ( 6, 3) =  16.180002d0  !      Carbon -      Lithium
    alpab_pm6( 6, 4) =   4.212882d0  !      Carbon -    Beryllium
    xab_pm6  ( 6, 4) =  25.035879d0  !      Carbon -    Beryllium
    alpab_pm6( 6, 5) =   2.919007d0  !      Carbon -        Boron
    xab_pm6  ( 6, 5) =   1.874859d0  !      Carbon -        Boron
    alpab_pm6( 6, 6) =   2.613713d0  !      Carbon -       Carbon
    xab_pm6  ( 6, 6) =   0.813510d0  !      Carbon -       Carbon
! N
    alpab_pm6( 7, 1) =   0.969406d0  !     Nitrogen -     Hydrogen
    xab_pm6  ( 7, 1) =   0.175506d0  !     Nitrogen -     Hydrogen
    alpab_pm6( 7, 2) =   2.814339d0  !     Nitrogen -       Helium
    xab_pm6  ( 7, 2) =   1.077861d0  !     Nitrogen -       Helium
    alpab_pm6( 7, 3) =   2.640623d0  !     Nitrogen -      Lithium
    xab_pm6  ( 7, 3) =   2.823403d0  !     Nitrogen -      Lithium
    alpab_pm6( 7, 4) =   2.580895d0  !     Nitrogen -    Beryllium
    xab_pm6  ( 7, 4) =   1.740605d0  !     Nitrogen -    Beryllium
    alpab_pm6( 7, 5) =   2.477004d0  !     Nitrogen -        Boron
    xab_pm6  ( 7, 5) =   0.952882d0  !     Nitrogen -        Boron
    alpab_pm6( 7, 6) =   2.686108d0  !     Nitrogen -       Carbon
    xab_pm6  ( 7, 6) =   0.859949d0  !     Nitrogen -       Carbon
    alpab_pm6( 7, 7) =   2.574502d0  !     Nitrogen -     Nitrogen
    xab_pm6  ( 7, 7) =   0.675313d0  !     Nitrogen -     Nitrogen
! O
    alpab_pm6( 8, 1) =   1.260942d0  !       Oxygen -     Hydrogen
    xab_pm6  ( 8, 1) =   0.192295d0  !       Oxygen -     Hydrogen
    alpab_pm6( 8, 2) =   3.653775d0  !       Oxygen -       Helium
    xab_pm6  ( 8, 2) =   6.684525d0  !       Oxygen -       Helium
    alpab_pm6( 8, 3) =   2.584442d0  !       Oxygen -      Lithium
    xab_pm6  ( 8, 3) =   1.968598d0  !       Oxygen -      Lithium
    alpab_pm6( 8, 4) =   3.051867d0  !       Oxygen -    Beryllium
    xab_pm6  ( 8, 4) =   3.218155d0  !       Oxygen -    Beryllium
    alpab_pm6( 8, 5) =   2.695351d0  !       Oxygen -        Boron
    xab_pm6  ( 8, 5) =   1.269801d0  !       Oxygen -        Boron
    alpab_pm6( 8, 6) =   2.889607d0  !       Oxygen -       Carbon
    xab_pm6  ( 8, 6) =   0.990211d0  !       Oxygen -       Carbon
    alpab_pm6( 8, 7) =   2.784292d0  !       Oxygen -     Nitrogen
    xab_pm6  ( 8, 7) =   0.764756d0  !       Oxygen -     Nitrogen
    alpab_pm6( 8, 8) =   2.623998d0  !       Oxygen -       Oxygen
    xab_pm6  ( 8, 8) =   0.535112d0  !       Oxygen -       Oxygen
! F
    alpab_pm6( 9, 1) =   3.136740d0  !     Fluorine -     Hydrogen
    xab_pm6  ( 9, 1) =   0.815802d0  !     Fluorine -     Hydrogen
    alpab_pm6( 9, 2) =   2.856543d0  !     Fluorine -       Helium
    xab_pm6  ( 9, 2) =   0.745107d0  !     Fluorine -       Helium
    alpab_pm6( 9, 3) =   3.043901d0  !     Fluorine -      Lithium
    xab_pm6  ( 9, 3) =   1.975985d0  !     Fluorine -      Lithium
    alpab_pm6( 9, 4) =   3.726923d0  !     Fluorine -    Beryllium
    xab_pm6  ( 9, 4) =   3.882993d0  !     Fluorine -    Beryllium
    alpab_pm6( 9, 5) =   2.823837d0  !     Fluorine -        Boron
    xab_pm6  ( 9, 5) =   0.862761d0  !     Fluorine -        Boron
    alpab_pm6( 9, 6) =   3.027600d0  !     Fluorine -       Carbon
    xab_pm6  ( 9, 6) =   0.732968d0  !     Fluorine -       Carbon
    alpab_pm6( 9, 7) =   2.856646d0  !     Fluorine -     Nitrogen
    xab_pm6  ( 9, 7) =   0.635854d0  !     Fluorine -     Nitrogen
    alpab_pm6( 9, 8) =   3.015444d0  !     Fluorine -       Oxygen
    xab_pm6  ( 9, 8) =   0.674251d0  !     Fluorine -       Oxygen
    alpab_pm6( 9, 9) =   3.175759d0  !     Fluorine -     Fluorine
    xab_pm6  ( 9, 9) =   0.681343d0  !     Fluorine -     Fluorine
! Ne
    alpab_pm6(10, 1) =   5.999680d0  !         Neon -     Hydrogen
    xab_pm6  (10, 1) =   5.535021d0  !         Neon -     Hydrogen
    alpab_pm6(10, 2) =   3.677758d0  !         Neon -       Helium
    xab_pm6  (10, 2) =   1.960924d0  !         Neon -       Helium
    alpab_pm6(10, 3) =   2.193666d0  !         Neon -      Lithium
    xab_pm6  (10, 3) =   0.704958d0  !         Neon -      Lithium
    alpab_pm6(10, 4) =   1.316588d0  !         Neon -    Beryllium
    xab_pm6  (10, 4) =   0.392628d0  !         Neon -    Beryllium
    alpab_pm6(10, 5) =   2.756190d0  !         Neon -        Boron
    xab_pm6  (10, 5) =   2.764140d0  !         Neon -        Boron
    alpab_pm6(10, 6) =   3.441188d0  !         Neon -       Carbon
    xab_pm6  (10, 6) =   5.468780d0  !         Neon -       Carbon
    alpab_pm6(10, 7) =   4.426370d0  !         Neon -     Nitrogen
    xab_pm6  (10, 7) =  29.999609d0  !         Neon -     Nitrogen
    alpab_pm6(10, 8) =   2.889587d0  !         Neon -       Oxygen
    xab_pm6  (10, 8) =   0.763899d0  !         Neon -       Oxygen
    alpab_pm6(10, 9) =   3.675611d0  !         Neon -     Fluorine
    xab_pm6  (10, 9) =   2.706754d0  !         Neon -     Fluorine
    alpab_pm6(10,10) =   3.974567d0  !         Neon -         Neon
    xab_pm6  (10,10) =   2.794830d0  !         Neon -         Neon
! Na
    alpab_pm6(11, 1) =     0.500326d0 !      Sodium -     Hydrogen
    xab_pm6  (11, 1) =     0.207831d0 !      Sodium -     Hydrogen
    alpab_pm6(11, 2) =     1.703029d0 !      Sodium -       Helium
    xab_pm6  (11, 2) =     4.282517d0 !      Sodium -       Helium
    alpab_pm6(11, 3) =     1.267299d0 !      Sodium -      Lithium
    xab_pm6  (11, 3) =     0.881482d0 !      Sodium -      Lithium
    alpab_pm6(11, 4) =     1.255480d0 !      Sodium -    Beryllium
    xab_pm6  (11, 4) =     3.121620d0 !      Sodium -    Beryllium
    alpab_pm6(11, 5) =     1.569961d0 !      Sodium -        Boron
    xab_pm6  (11, 5) =     3.188608d0 !      Sodium -        Boron
    alpab_pm6(11, 6) =     2.196050d0 !      Sodium -       Carbon
    xab_pm6  (11, 6) =     4.520429d0 !      Sodium -       Carbon
    alpab_pm6(11, 7) =     2.494384d0 !      Sodium -     Nitrogen
    xab_pm6  (11, 7) =     8.586387d0 !      Sodium -     Nitrogen
    alpab_pm6(11, 8) =     1.981449d0 !      Sodium -       Oxygen
    xab_pm6  (11, 8) =     3.270079d0 !      Sodium -       Oxygen
    alpab_pm6(11, 9) =     2.619551d0 !      Sodium -     Fluorine
    xab_pm6  (11, 9) =     7.047351d0 !      Sodium -     Fluorine
    alpab_pm6(11,10) =     1.774236d0 !      Sodium -         Neon
    xab_pm6  (11,10) =     1.343037d0 !      Sodium -         Neon
    alpab_pm6(11,11) =     0.446435d0 !      Sodium -       Sodium
    xab_pm6  (11,11) =     0.287137d0 !      Sodium -       Sodium
! Mg
    alpab_pm6(12, 1) =     2.651594d0 !   Magnesium -     Hydrogen
    xab_pm6  (12, 1) =     7.758237d0 !   Magnesium -     Hydrogen
    alpab_pm6(12, 2) =     2.210603d0 !   Magnesium -       Helium
    xab_pm6  (12, 2) =     3.725850d0 !   Magnesium -       Helium
    alpab_pm6(12, 3) =     1.184380d0 !   Magnesium -      Lithium
    xab_pm6  (12, 3) =     2.490250d0 !   Magnesium -      Lithium
    alpab_pm6(12, 4) =     1.557591d0 !   Magnesium -    Beryllium
    xab_pm6  (12, 4) =     2.066392d0 !   Magnesium -    Beryllium
    alpab_pm6(12, 5) =     2.527441d0 !   Magnesium -        Boron
    xab_pm6  (12, 5) =     6.146701d0 !   Magnesium -        Boron
    alpab_pm6(12, 6) =     3.040946d0 !   Magnesium -       Carbon
    xab_pm6  (12, 6) =    10.517690d0 !   Magnesium -       Carbon
    alpab_pm6(12, 7) =     2.079125d0 !   Magnesium -     Nitrogen
    xab_pm6  (12, 7) =     1.208075d0 !   Magnesium -     Nitrogen
    alpab_pm6(12, 8) =     2.251520d0 !   Magnesium -       Oxygen
    xab_pm6  (12, 8) =     1.535734d0 !   Magnesium -       Oxygen
    alpab_pm6(12, 9) =     3.362208d0 !   Magnesium -     Fluorine
    xab_pm6  (12, 9) =     5.859023d0 !   Magnesium -     Fluorine
    alpab_pm6(12,10) =     2.031676d0 !   Magnesium -         Neon
    xab_pm6  (12,10) =     1.214859d0 !   Magnesium -         Neon
    alpab_pm6(12,11) =     1.506773d0 !   Magnesium -       Sodium
    xab_pm6  (12,11) =     8.675619d0 !   Magnesium -       Sodium
    alpab_pm6(12,12) =     1.093573d0 !   Magnesium -    Magnesium
    xab_pm6  (12,12) =     0.465645d0 !   Magnesium -    Magnesium
! Al
    alpab_pm6(13, 1) =     2.025996d0 !   Aluminium -     Hydrogen
    xab_pm6  (13, 1) =     2.958379d0 !   Aluminium -     Hydrogen
    alpab_pm6(13, 2) =     2.255830d0 !   Aluminium -       Helium
    xab_pm6  (13, 2) =     2.701400d0 !   Aluminium -       Helium
    alpab_pm6(13, 3) =     1.581593d0 !   Aluminium -      Lithium
    xab_pm6  (13, 3) =     1.106819d0 !   Aluminium -      Lithium
    alpab_pm6(13, 4) =     1.938237d0 !   Aluminium -    Beryllium
    xab_pm6  (13, 4) =     5.037214d0 !   Aluminium -    Beryllium
    alpab_pm6(13, 5) =     2.059569d0 !   Aluminium -        Boron
    xab_pm6  (13, 5) =     2.741479d0 !   Aluminium -        Boron
    alpab_pm6(13, 6) =     2.267440d0 !   Aluminium -       Carbon
    xab_pm6  (13, 6) =     2.928056d0 !   Aluminium -       Carbon
    alpab_pm6(13, 7) =     2.009754d0 !   Aluminium -     Nitrogen
    xab_pm6  (13, 7) =     1.345202d0 !   Aluminium -     Nitrogen
    alpab_pm6(13, 8) =     2.498660d0 !   Aluminium -       Oxygen
    xab_pm6  (13, 8) =     2.131396d0 !   Aluminium -       Oxygen
    alpab_pm6(13, 9) =     3.084258d0 !   Aluminium -     Fluorine
    xab_pm6  (13, 9) =     1.975635d0 !   Aluminium -     Fluorine
    alpab_pm6(13,10) =     2.447869d0 !   Aluminium -         Neon
    xab_pm6  (13,10) =     1.709200d0 !   Aluminium -         Neon
    alpab_pm6(13,11) =     1.202871d0 !   Aluminium -       Sodium
    xab_pm6  (13,11) =     2.071847d0 !   Aluminium -       Sodium
    alpab_pm6(13,12) =     1.972530d0 !   Aluminium -    Magnesium
    xab_pm6  (13,12) =    13.472443d0 !   Aluminium -    Magnesium
    alpab_pm6(13,13) =     1.387714d0 !   Aluminium -    Aluminium
    xab_pm6  (13,13) =     2.139200d0 !   Aluminium -    Aluminium
!Si
    alpab_pm6(14, 1) =     1.896950d0 !     Silicon -     Hydrogen
    xab_pm6  (14, 1) =     0.924196d0 !     Silicon -     Hydrogen
    alpab_pm6(14, 2) =     2.040498d0 !     Silicon -       Helium
    xab_pm6  (14, 2) =     1.853583d0 !     Silicon -       Helium
    alpab_pm6(14, 3) =     1.789609d0 !     Silicon -      Lithium
    xab_pm6  (14, 3) =     3.090791d0 !     Silicon -      Lithium
    alpab_pm6(14, 4) =     1.263132d0 !     Silicon -    Beryllium
    xab_pm6  (14, 4) =     0.623433d0 !     Silicon -    Beryllium
    alpab_pm6(14, 5) =     1.982653d0 !     Silicon -        Boron
    xab_pm6  (14, 5) =     1.028287d0 !     Silicon -        Boron
    alpab_pm6(14, 6) =     1.984498d0 !     Silicon -       Carbon
    xab_pm6  (14, 6) =     0.785745d0 !     Silicon -       Carbon
    alpab_pm6(14, 7) =     1.818988d0 !     Silicon -     Nitrogen
    xab_pm6  (14, 7) =     0.592972d0 !     Silicon -     Nitrogen
    alpab_pm6(14, 8) =     1.923600d0 !     Silicon -       Oxygen
    xab_pm6  (14, 8) =     0.751095d0 !     Silicon -       Oxygen
    alpab_pm6(14, 9) =     2.131028d0 !     Silicon -     Fluorine
    xab_pm6  (14, 9) =     0.543516d0 !     Silicon -     Fluorine
    alpab_pm6(14,10) =     2.287784d0 !     Silicon -         Neon
    xab_pm6  (14,10) =    14.378676d0 !     Silicon -         Neon
    alpab_pm6(14,11) =     2.007615d0 !     Silicon -       Sodium
    xab_pm6  (14,11) =     9.237644d0 !     Silicon -       Sodium
    alpab_pm6(14,12) =     3.139749d0 !     Silicon -    Magnesium
    xab_pm6  (14,12) =    29.994520d0 !     Silicon -    Magnesium
    alpab_pm6(14,13) =     1.900000d0 !     Silicon -    Aluminium
    xab_pm6  (14,13) =     2.000000d0 !     Silicon -    Aluminium
    alpab_pm6(14,14) =     1.329000d0 !     Silicon -      Silicon
    xab_pm6  (14,14) =     0.273477d0 !     Silicon -      Silicon
!P  
    alpab_pm6(15, 1) =     1.926537d0 !  Phosphorus -     Hydrogen
    xab_pm6  (15, 1) =     1.234986d0 !  Phosphorus -     Hydrogen
    alpab_pm6(15, 2) =     2.093158d0 !  Phosphorus -       Helium
    xab_pm6  (15, 2) =     1.490218d0 !  Phosphorus -       Helium
    alpab_pm6(15, 3) =     1.394544d0 !  Phosphorus -      Lithium
    xab_pm6  (15, 3) =     1.122950d0 !  Phosphorus -      Lithium
    alpab_pm6(15, 4) =     1.800070d0 !  Phosphorus -    Beryllium
    xab_pm6  (15, 4) =     1.684831d0 !  Phosphorus -    Beryllium
    alpab_pm6(15, 5) =     1.923168d0 !  Phosphorus -        Boron
    xab_pm6  (15, 5) =     1.450886d0 !  Phosphorus -        Boron
    alpab_pm6(15, 6) =     1.994653d0 !  Phosphorus -       Carbon
    xab_pm6  (15, 6) =     0.979512d0 !  Phosphorus -       Carbon
    alpab_pm6(15, 7) =     2.147042d0 !  Phosphorus -     Nitrogen
    xab_pm6  (15, 7) =     0.972154d0 !  Phosphorus -     Nitrogen
    alpab_pm6(15, 8) =     2.220768d0 !  Phosphorus -       Oxygen
    xab_pm6  (15, 8) =     0.878705d0 !  Phosphorus -       Oxygen
    alpab_pm6(15, 9) =     2.234356d0 !  Phosphorus -     Fluorine
    xab_pm6  (15, 9) =     0.514575d0 !  Phosphorus -     Fluorine
    alpab_pm6(15,10) =     2.219036d0 !  Phosphorus -         Neon
    xab_pm6  (15,10) =     0.774954d0 !  Phosphorus -         Neon
    alpab_pm6(15,11) =     1.500320d0 !  Phosphorus -       Sodium
    xab_pm6  (15,11) =     2.837095d0 !  Phosphorus -       Sodium
    alpab_pm6(15,12) =     1.383773d0 !  Phosphorus -    Magnesium
    xab_pm6  (15,12) =     1.177881d0 !  Phosphorus -    Magnesium
    alpab_pm6(15,13) =     1.980727d0 !  Phosphorus -    Aluminium
    xab_pm6  (15,13) =     5.050816d0 !  Phosphorus -    Aluminium
    alpab_pm6(15,14) =     3.313466d0 !  Phosphorus -      Silicon
    xab_pm6  (15,14) =    13.239121d0 !  Phosphorus -      Silicon
    alpab_pm6(15,15) =     1.505792d0 !  Phosphorus -   Phosphorus
    xab_pm6  (15,15) =     0.902501d0 !  Phosphorus -   Phosphorus
! S 
    alpab_pm6(16, 1) =     2.215975d0 !      Sulfur -     Hydrogen
    xab_pm6  (16, 1) =     0.849712d0 !      Sulfur -     Hydrogen
    alpab_pm6(16, 2) =     1.959149d0 !      Sulfur -       Helium
    xab_pm6  (16, 2) =     0.437618d0 !      Sulfur -       Helium
    alpab_pm6(16, 3) =     2.294275d0 !      Sulfur -      Lithium
    xab_pm6  (16, 3) =     2.642502d0 !      Sulfur -      Lithium
    alpab_pm6(16, 4) =     2.781736d0 !      Sulfur -    Beryllium
    xab_pm6  (16, 4) =     3.791565d0 !      Sulfur -    Beryllium
    alpab_pm6(16, 5) =     2.403696d0 !      Sulfur -        Boron
    xab_pm6  (16, 5) =     1.125394d0 !      Sulfur -        Boron
    alpab_pm6(16, 6) =     2.210305d0 !      Sulfur -       Carbon
    xab_pm6  (16, 6) =     0.666849d0 !      Sulfur -       Carbon
    alpab_pm6(16, 7) =     2.289990d0 !      Sulfur -     Nitrogen
    xab_pm6  (16, 7) =     0.738710d0 !      Sulfur -     Nitrogen
    alpab_pm6(16, 8) =     2.383289d0 !      Sulfur -       Oxygen
    xab_pm6  (16, 8) =     0.747215d0 !      Sulfur -       Oxygen
    alpab_pm6(16, 9) =     2.187186d0 !      Sulfur -     Fluorine
    xab_pm6  (16, 9) =     0.375251d0 !      Sulfur -     Fluorine
    alpab_pm6(16,10) =     2.787058d0 !      Sulfur -         Neon
    xab_pm6  (16,10) =     3.296160d0 !      Sulfur -         Neon
    alpab_pm6(16,11) =     1.400850d0 !      Sulfur -       Sodium
    xab_pm6  (16,11) =     0.852434d0 !      Sulfur -       Sodium
    alpab_pm6(16,12) =     1.500163d0 !      Sulfur -    Magnesium
    xab_pm6  (16,12) =     0.500748d0 !      Sulfur -    Magnesium
    alpab_pm6(16,13) =     1.976705d0 !      Sulfur -    Aluminium
    xab_pm6  (16,13) =     2.347384d0 !      Sulfur -    Aluminium
    alpab_pm6(16,14) =     1.885916d0 !      Sulfur -      Silicon
    xab_pm6  (16,14) =     0.876658d0 !      Sulfur -      Silicon
    alpab_pm6(16,15) =     1.595325d0 !      Sulfur -   Phosphorus
    xab_pm6  (16,15) =     0.562266d0 !      Sulfur -   Phosphorus
    alpab_pm6(16,16) =     1.794556d0 !      Sulfur -       Sulfur
    xab_pm6  (16,16) =     0.473856d0 !      Sulfur -       Sulfur
! Cl
    alpab_pm6(17, 1) =     2.402886d0 !    Chlorine -     Hydrogen
    xab_pm6  (17, 1) =     0.754831d0 !    Chlorine -     Hydrogen
    alpab_pm6(17, 2) =     1.671677d0 !    Chlorine -       Helium
    xab_pm6  (17, 2) =     0.272964d0 !    Chlorine -       Helium
    alpab_pm6(17, 3) =     2.783001d0 !    Chlorine -      Lithium
    xab_pm6  (17, 3) =     4.227794d0 !    Chlorine -      Lithium
    alpab_pm6(17, 4) =     2.822676d0 !    Chlorine -    Beryllium
    xab_pm6  (17, 4) =     2.507275d0 !    Chlorine -    Beryllium
    alpab_pm6(17, 5) =     2.259323d0 !    Chlorine -        Boron
    xab_pm6  (17, 5) =     0.822129d0 !    Chlorine -        Boron
    alpab_pm6(17, 6) =     2.162197d0 !    Chlorine -       Carbon
    xab_pm6  (17, 6) =     0.515787d0 !    Chlorine -       Carbon
    alpab_pm6(17, 7) =     2.172134d0 !    Chlorine -     Nitrogen
    xab_pm6  (17, 7) =     0.520745d0 !    Chlorine -     Nitrogen
    alpab_pm6(17, 8) =     2.323236d0 !    Chlorine -       Oxygen
    xab_pm6  (17, 8) =     0.585510d0 !    Chlorine -       Oxygen
    alpab_pm6(17, 9) =     2.313270d0 !    Chlorine -     Fluorine
    xab_pm6  (17, 9) =     0.411124d0 !    Chlorine -     Fluorine
    alpab_pm6(17,10) =     1.703151d0 !    Chlorine -         Neon
    xab_pm6  (17,10) =     0.125133d0 !    Chlorine -         Neon
    alpab_pm6(17,11) =     1.816429d0 !    Chlorine -       Sodium
    xab_pm6  (17,11) =     1.357894d0 !    Chlorine -       Sodium
    alpab_pm6(17,12) =     2.391806d0 !    Chlorine -    Magnesium
    xab_pm6  (17,12) =     2.430856d0 !    Chlorine -    Magnesium
    alpab_pm6(17,13) =     2.125939d0 !    Chlorine -    Aluminium
    xab_pm6  (17,13) =     2.153451d0 !    Chlorine -    Aluminium
    alpab_pm6(17,14) =     1.684978d0 !    Chlorine -      Silicon
    xab_pm6  (17,14) =     0.513000d0 !    Chlorine -      Silicon
    alpab_pm6(17,15) =     1.468306d0 !    Chlorine -   Phosphorus
    xab_pm6  (17,15) =     0.352361d0 !    Chlorine -   Phosphorus
    alpab_pm6(17,16) =     1.715435d0 !    Chlorine -       Sulfur
    xab_pm6  (17,16) =     0.356971d0 !    Chlorine -       Sulfur
    alpab_pm6(17,17) =     1.823239d0 !    Chlorine -     Chlorine
    xab_pm6  (17,17) =     0.332919d0 !    Chlorine -     Chlorine

! Ar
    alpab_pm6(18, 1) =     4.056167d0 !       Argon -     Hydrogen
    xab_pm6  (18, 1) =     3.933445d0 !       Argon -     Hydrogen
    alpab_pm6(18, 2) =     2.716562d0 !       Argon -       Helium
    xab_pm6  (18, 2) =     1.177211d0 !       Argon -       Helium
    alpab_pm6(18, 3) =     3.122895d0 !       Argon -      Lithium
    xab_pm6  (18, 3) =     3.362910d0 !       Argon -      Lithium
    alpab_pm6(18, 4) =     3.044007d0 !       Argon -    Beryllium
    xab_pm6  (18, 4) =     2.755492d0 !       Argon -    Beryllium
    alpab_pm6(18, 5) =     2.415471d0 !       Argon -        Boron
    xab_pm6  (18, 5) =     1.931586d0 !       Argon -        Boron
    alpab_pm6(18, 6) =     1.471309d0 !       Argon -       Carbon
    xab_pm6  (18, 6) =     0.122309d0 !       Argon -       Carbon
    alpab_pm6(18, 7) =     2.326805d0 !       Argon -     Nitrogen
    xab_pm6  (18, 7) =     0.562581d0 !       Argon -     Nitrogen
    alpab_pm6(18, 8) =     2.240673d0 !       Argon -       Oxygen
    xab_pm6  (18, 8) =     0.355795d0 !       Argon -       Oxygen
    alpab_pm6(18, 9) =     3.920658d0 !       Argon -     Fluorine
    xab_pm6  (18, 9) =     9.269715d0 !       Argon -     Fluorine
    alpab_pm6(18,10) =     2.963747d0 !       Argon -         Neon
    xab_pm6  (18,10) =     1.304697d0 !       Argon -         Neon
    alpab_pm6(18,11) =     2.167677d0 !       Argon -       Sodium
    xab_pm6  (18,11) =     3.398138d0 !       Argon -       Sodium
    alpab_pm6(18,12) =     2.092664d0 !       Argon -    Magnesium
    xab_pm6  (18,12) =     1.970638d0 !       Argon -    Magnesium
    alpab_pm6(18,13) =     2.645165d0 !       Argon -    Aluminium
    xab_pm6  (18,13) =     1.852009d0 !       Argon -    Aluminium
    alpab_pm6(18,14) =     1.780350d0 !       Argon -      Silicon
    xab_pm6  (18,14) =     1.067890d0 !       Argon -      Silicon
    alpab_pm6(18,15) =     4.372516d0 !       Argon -   Phosphorus
    xab_pm6  (18,15) =     0.171014d0 !       Argon -   Phosphorus
    alpab_pm6(18,16) =     2.049398d0 !       Argon -       Sulfur
    xab_pm6  (18,16) =     0.653769d0 !       Argon -       Sulfur
    alpab_pm6(18,17) =     2.554449d0 !       Argon -     Chlorine
    xab_pm6  (18,17) =     2.256094d0 !       Argon -     Chlorine
    alpab_pm6(18,18) =     2.306432d0 !       Argon -        Argon
    xab_pm6  (18,18) =     0.972699d0 !       Argon -        Argon
! K
    alpab_pm6(19, 1) =     0.648173d0 !   Potassium -     Hydrogen
    xab_pm6  (19, 1) =     0.369340d0 !   Potassium -     Hydrogen
    alpab_pm6(19, 2) =     1.418501d0 !   Potassium -       Helium
    xab_pm6  (19, 2) =     2.895045d0 !   Potassium -       Helium
    alpab_pm6(19, 3) =     1.036487d0 !   Potassium -      Lithium
    xab_pm6  (19, 3) =     4.374567d0 !   Potassium -      Lithium
    alpab_pm6(19, 4) =     1.931888d0 !   Potassium -    Beryllium
    xab_pm6  (19, 4) =     6.732221d0 !   Potassium -    Beryllium
    alpab_pm6(19, 5) =     2.031768d0 !   Potassium -        Boron
    xab_pm6  (19, 5) =     8.900541d0 !   Potassium -        Boron
    alpab_pm6(19, 6) =     2.241757d0 !   Potassium -       Carbon
    xab_pm6  (19, 6) =    10.317987d0 !   Potassium -       Carbon
    alpab_pm6(19, 7) =     2.325859d0 !   Potassium -     Nitrogen
    xab_pm6  (19, 7) =     7.977707d0 !   Potassium -     Nitrogen
    alpab_pm6(19, 8) =     1.508571d0 !   Potassium -       Oxygen
    xab_pm6  (19, 8) =     1.012275d0 !   Potassium -       Oxygen
    alpab_pm6(19, 9) =     3.182817d0 !   Potassium -     Fluorine
    xab_pm6  (19, 9) =     6.592971d0 !   Potassium -     Fluorine
    alpab_pm6(19,10) =     1.138021d0 !   Potassium -         Neon
    xab_pm6  (19,10) =     0.233995d0 !   Potassium -         Neon
    alpab_pm6(19,11) =     0.884307d0 !   Potassium -       Sodium
    xab_pm6  (19,11) =     5.563027d0 !   Potassium -       Sodium
    alpab_pm6(19,12) =     0.884810d0 !   Potassium -    Magnesium
    xab_pm6  (19,12) =     3.290502d0 !   Potassium -    Magnesium
    alpab_pm6(19,13) =     1.976076d0 !   Potassium -    Aluminium
    xab_pm6  (19,13) =    29.944708d0 !   Potassium -    Aluminium
    alpab_pm6(19,14) =     1.675930d0 !   Potassium -      Silicon
    xab_pm6  (19,14) =     8.279200d0 !   Potassium -      Silicon
    alpab_pm6(19,15) =     1.443738d0 !   Potassium -   Phosphorus
    xab_pm6  (19,15) =     4.475384d0 !   Potassium -   Phosphorus
    alpab_pm6(19,16) =     2.512156d0 !   Potassium -       Sulfur
    xab_pm6  (19,16) =    29.528951d0 !   Potassium -       Sulfur
    alpab_pm6(19,17) =     1.622163d0 !   Potassium -     Chlorine
    xab_pm6  (19,17) =     1.231481d0 !   Potassium -     Chlorine
    alpab_pm6(19,18) =     2.302803d0 !   Potassium -        Argon
    xab_pm6  (19,18) =     9.710508d0 !   Potassium -        Argon
    alpab_pm6(19,19) =     1.435514d0 !   Potassium -    Potassium
    xab_pm6  (19,19) =     5.934329d0 !   Potassium -    Potassium
! Ca
    alpab_pm6(20, 1) =     2.141859d0 !     Calcium -     Hydrogen
    xab_pm6  (20, 1) =     7.728606d0 !     Calcium -     Hydrogen
    alpab_pm6(20, 2) =     1.719847d0 !     Calcium -       Helium
    xab_pm6  (20, 2) =     2.913852d0 !     Calcium -       Helium
    alpab_pm6(20, 5) =     1.700010d0 !     Calcium -        Boron
    xab_pm6  (20, 5) =     1.700010d0 !     Calcium -        Boron
    alpab_pm6(20, 6) =     1.035305d0 !     Calcium -       Carbon
    xab_pm6  (20, 6) =     0.148450d0 !     Calcium -       Carbon
    alpab_pm6(20, 7) =     2.386600d0 !     Calcium -     Nitrogen
    xab_pm6  (20, 7) =     2.988074d0 !     Calcium -     Nitrogen
    alpab_pm6(20, 8) =     3.263897d0 !     Calcium -       Oxygen
    xab_pm6  (20, 8) =    17.028946d0 !     Calcium -       Oxygen
    alpab_pm6(20, 9) =     2.645053d0 !     Calcium -     Fluorine
    xab_pm6  (20, 9) =     3.482821d0 !     Calcium -     Fluorine
    alpab_pm6(20,10) =     0.954530d0 !     Calcium -         Neon
    xab_pm6  (20,10) =     0.332586d0 !     Calcium -         Neon
    alpab_pm6(20,11) =     3.107104d0 !     Calcium -       Sodium
    xab_pm6  (20,11) =     9.657509d0 !     Calcium -       Sodium
    alpab_pm6(20,12) =     2.299800d0 !     Calcium -    Magnesium
    xab_pm6  (20,12) =     8.599800d0 !     Calcium -    Magnesium
    alpab_pm6(20,13) =     1.612565d0 !     Calcium -    Aluminium
    xab_pm6  (20,13) =     4.188555d0 !     Calcium -    Aluminium
    alpab_pm6(20,14) =     1.218788d0 !     Calcium -      Silicon
    xab_pm6  (20,14) =     0.336233d0 !     Calcium -      Silicon
    alpab_pm6(20,15) =     1.024142d0 !     Calcium -   Phosphorus
    xab_pm6  (20,15) =     0.410840d0 !     Calcium -   Phosphorus
    alpab_pm6(20,16) =     0.958171d0 !     Calcium -       Sulfur
    xab_pm6  (20,16) =     0.325739d0 !     Calcium -       Sulfur
    alpab_pm6(20,17) =     2.383391d0 !     Calcium -     Chlorine
    xab_pm6  (20,17) =     5.956144d0 !     Calcium -     Chlorine
    alpab_pm6(20,18) =     1.034881d0 !     Calcium -        Argon
    xab_pm6  (20,18) =     0.291072d0 !     Calcium -        Argon
    alpab_pm6(20,19) =     1.119200d0 !     Calcium -    Potassium
    xab_pm6  (20,19) =     1.240320d0 !     Calcium -    Potassium
    alpab_pm6(20,20) =     1.889674d0 !     Calcium -      Calcium
    xab_pm6  (20,20) =    30.003591d0 !     Calcium -      Calcium
! Sc
    alpab_pm6(21, 1) =     1.179485d0 !    Scandium -     Hydrogen
    xab_pm6  (21, 1) =     0.351199d0 !    Scandium -     Hydrogen
    alpab_pm6(21, 6) =     2.630490d0 !    Scandium -       Carbon
    xab_pm6  (21, 6) =     8.608052d0 !    Scandium -       Carbon
    alpab_pm6(21, 7) =     2.270004d0 !    Scandium -     Nitrogen
    xab_pm6  (21, 7) =     3.231881d0 !    Scandium -     Nitrogen
    alpab_pm6(21, 8) =     2.256516d0 !    Scandium -       Oxygen
    xab_pm6  (21, 8) =     3.058672d0 !    Scandium -       Oxygen
    alpab_pm6(21, 9) =     3.107985d0 !    Scandium -     Fluorine
    xab_pm6  (21, 9) =     7.252347d0 !    Scandium -     Fluorine
    alpab_pm6(21,13) =     1.003550d0 !    Scandium -    Aluminium
    xab_pm6  (21,13) =     0.500620d0 !    Scandium -    Aluminium
    alpab_pm6(21,14) =     2.016870d0 !    Scandium -      Silicon
    xab_pm6  (21,14) =     3.219070d0 !    Scandium -      Silicon
    alpab_pm6(21,15) =     0.868165d0 !    Scandium -   Phosphorus
    xab_pm6  (21,15) =     0.626749d0 !    Scandium -   Phosphorus
    alpab_pm6(21,16) =     0.422939d0 !    Scandium -       Sulfur
    xab_pm6  (21,16) =     0.211850d0 !    Scandium -       Sulfur
    alpab_pm6(21,17) =     2.141474d0 !    Scandium -     Chlorine
    xab_pm6  (21,17) =     2.996129d0 !    Scandium -     Chlorine
    alpab_pm6(21,21) =     1.132838d0 !    Scandium -     Scandium 
    xab_pm6  (21,21) =     2.598166d0 !    Scandium -     Scandium
! Ti
    alpab_pm6(22, 1) =     0.832669d0 !    Titanium -     Hydrogen
    xab_pm6  (22, 1) =     0.143722d0 !    Titanium -     Hydrogen
    alpab_pm6(22, 5) =     1.628710d0 !    Titanium -        Boron
    xab_pm6  (22, 5) =     0.649360d0 !    Titanium -        Boron
    alpab_pm6(22, 6) =     1.597973d0 !    Titanium -       Carbon
    xab_pm6  (22, 6) =     0.416706d0 !    Titanium -       Carbon
    alpab_pm6(22, 7) =     1.678686d0 !    Titanium -     Nitrogen
    xab_pm6  (22, 7) =     0.545461d0 !    Titanium -     Nitrogen
    alpab_pm6(22, 8) =     1.789118d0 !    Titanium -       Oxygen
    xab_pm6  (22, 8) =     0.799486d0 !    Titanium -       Oxygen
    alpab_pm6(22, 9) =     2.307087d0 !    Titanium -     Fluorine
    xab_pm6  (22, 9) =     1.085742d0 !    Titanium -     Fluorine
    alpab_pm6(22,12) =     1.911340d0 !    Titanium -    Magnesium
    xab_pm6  (22,12) =     4.330240d0 !    Titanium -    Magnesium
    alpab_pm6(22,13) =     1.369486d0 !    Titanium -    Aluminium
    xab_pm6  (22,13) =     2.091841d0 !    Titanium -    Aluminium
    alpab_pm6(22,14) =     2.856038d0 !    Titanium -      Silicon
    xab_pm6  (22,14) =     6.773815d0 !    Titanium -      Silicon
    alpab_pm6(22,15) =     2.151929d0 !    Titanium -   Phosphorus
    xab_pm6  (22,15) =     4.150500d0 !    Titanium -   Phosphorus
    alpab_pm6(22,16) =     1.846439d0 !    Titanium -       Sulfur
    xab_pm6  (22,16) =     0.943784d0 !    Titanium -       Sulfur
    alpab_pm6(22,17) =     1.461034d0 !    Titanium -     Chlorine
    xab_pm6  (22,17) =     0.333297d0 !    Titanium -     Chlorine
    alpab_pm6(22,20) =     2.000000d0 !    Titanium -      Calcium
    xab_pm6  (22,20) =     4.109141d0 !    Titanium -      Calcium
    alpab_pm6(22,22) =     2.648597d0 !    Titanium -     Titanium 
    xab_pm6  (22,22) =     2.000000d0 !    Titanium -     Titanium
! V
    alpab_pm6(23, 1) =     1.280133d0 !    Vanadium -     Hydrogen
    xab_pm6  (23, 1) =     0.105204d0 !    Vanadium -     Hydrogen
    alpab_pm6(23, 6) =     2.789855d0 !    Vanadium -       Carbon
    xab_pm6  (23, 6) =     1.938760d0 !    Vanadium -       Carbon
    alpab_pm6(23, 7) =     1.607540d0 !    Vanadium -     Nitrogen
    xab_pm6  (23, 7) =     0.276725d0 !    Vanadium -     Nitrogen
    alpab_pm6(23, 8) =     1.623973d0 !    Vanadium -       Oxygen
    xab_pm6  (23, 8) =     0.415312d0 !    Vanadium -       Oxygen
    alpab_pm6(23, 9) =     1.825160d0 !    Vanadium -     Fluorine
    xab_pm6  (23, 9) =     0.342815d0 !    Vanadium -     Fluorine
    alpab_pm6(23,11) =     2.551010d0 !    Vanadium -       Sodium
    xab_pm6  (23,11) =     8.276020d0 !    Vanadium -       Sodiium
    alpab_pm6(23,15) =     2.549154d0 !    Vanadium -   Phosphorus
    xab_pm6  (23,15) =     6.250624d0 !    Vanadium -   Phosphorus
    alpab_pm6(23,16) =     2.704124d0 !    Vanadium -       Sulfur
    xab_pm6  (23,16) =     2.035039d0 !    Vanadium -       Sulfur
    alpab_pm6(23,17) =     1.688529d0 !    Vanadium -     Chlorine
    xab_pm6  (23,17) =     0.243657d0 !    Vanadium -     Chlorine
    alpab_pm6(23,19) =     4.521360d0 !    Vanadium -    Potassium
    xab_pm6  (23,19) =     2.026590d0 !    Vanadium -    Potassium 
    alpab_pm6(23,23) =     4.832391d0 !    Vanadium -     Vanadium 
    xab_pm6  (23,23) =    10.779892d0 !    Vanadium -     Vanadium
! Cr
    alpab_pm6(24, 1) =     0.882661d0 !    Chromium -     Hydrogen
    xab_pm6  (24, 1) =     0.044469d0 !    Chromium -     Hydrogen
    alpab_pm6(24, 6) =     3.656754d0 !    Chromium -       Carbon
    xab_pm6  (24, 6) =     6.110187d0 !    Chromium -       Carbon
    alpab_pm6(24, 7) =     3.029186d0 !    Chromium -     Nitrogen
    xab_pm6  (24, 7) =     1.920324d0 !    Chromium -     Nitrogen
    alpab_pm6(24, 8) =     2.500000d0 !    Chromium -       Oxygen
    xab_pm6  (24, 8) =     1.055511d0 !    Chromium -       Oxygen
    alpab_pm6(24, 9) =     2.716521d0 !    Chromium -     Fluorine
    xab_pm6  (24, 9) =     0.737607d0 !    Chromium -     Fluorine
    alpab_pm6(24,11) =     2.295056d0 !    Chromium -       Sodium
    xab_pm6  (24,11) =     8.364274d0 !    Chromium -       Sodiium
    alpab_pm6(24,14) =     1.860760d0 !    Chromium -      Silicon
    xab_pm6  (24,14) =     1.029110d0 !    Chromium -      Silicon
    alpab_pm6(24,15) =     1.695383d0 !    Chromium -   Phosphorus
    xab_pm6  (24,15) =     0.600177d0 !    Chromium -   Phosphorus
    alpab_pm6(24,16) =     2.260978d0 !    Chromium -       Sulfur
    xab_pm6  (24,16) =     0.550334d0 !    Chromium -       Sulfur
    alpab_pm6(24,17) =     2.152618d0 !    Chromium -     Chlorine
    xab_pm6  (24,17) =     0.369073d0 !    Chromium -     Chlorine
    alpab_pm6(24,19) =     2.000000d0 !    Chromium -    Potassium
    xab_pm6  (24,19) =     2.000000d0 !    Chromium -    Potassium 
    alpab_pm6(24,24) =     4.655419d0 !    Chromium -     Chromium 
    xab_pm6  (24,24) =    10.318607d0 !    Chromium -     Chromium
! Mn
    alpab_pm6(25, 1) =     2.309940d0 !   Manganese -     Hydrogen
    xab_pm6  (25, 1) =     1.269210d0 !   Manganese -     Hydrogen
    alpab_pm6(25, 6) =     3.000750d0 !   Manganese -       Carbon
    xab_pm6  (25, 6) =     2.583110d0 !   Manganese -       Carbon
    alpab_pm6(25, 7) =     2.921470d0 !   Manganese -     Nitrogen
    xab_pm6  (25, 7) =     1.956750d0 !   Manganese -     Nitrogen
    alpab_pm6(25, 8) =     2.577540d0 !   Manganese -       Oxygen
    xab_pm6  (25, 8) =     1.285620d0 !   Manganese -       Oxygen
    alpab_pm6(25, 9) =     2.791950d0 !   Manganese -     Fluorine
    xab_pm6  (25, 9) =     1.113070d0 !   Manganese -     Fluorine
    alpab_pm6(25,13) =     1.768360d0 !   Manganese -     Aluminum
    xab_pm6  (25,13) =     1.040790d0 !   Manganese -     Aluminum
    alpab_pm6(25,14) =     1.937959d0 !   Manganese -      Silicon
    xab_pm6  (25,14) =     0.950580d0 !   Manganese -      Silicon
    alpab_pm6(25,15) =     1.947020d0 !   Manganese -   Phosphorus
    xab_pm6  (25,15) =     1.130320d0 !   Manganese -   Phosphorus
    alpab_pm6(25,16) =     2.482510d0 !   Manganese -       Sulfur
    xab_pm6  (25,16) =     1.612650d0 !   Manganese -       Sulfur
    alpab_pm6(25,17) =     1.657010d0 !   Manganese -     Chlorine
    xab_pm6  (25,17) =     0.201850d0 !   Manganese -     Chlorine
    alpab_pm6(25,20) =     1.491440d0 !   Manganese -      Calcium
    xab_pm6  (25,20) =     0.620180d0 !   Manganese -      Calcium 
    alpab_pm6(25,25) =     2.665420d0 !   Manganese -    Manganese
    xab_pm6  (25,25) =     2.460040d0 !   Manganese -    Manganese
! Fe
    alpab_pm6(26, 1) =     0.854488d0 !        Iron -     Hydrogen
    xab_pm6  (26, 1) =     0.025195d0 !        Iron -     Hydrogen
    alpab_pm6(26, 6) =     3.991343d0 !        Iron -       Carbon
    xab_pm6  (26, 6) =     0.366835d0 !        Iron -       Carbon
    alpab_pm6(26, 7) =     2.500486d0 !        Iron -     Nitrogen
    xab_pm6  (26, 7) =     0.155342d0 !        Iron -     Nitrogen
    alpab_pm6(26, 8) =     1.726313d0 !        Iron -       Oxygen
    xab_pm6  (26, 8) =     0.136422d0 !        Iron -       Oxygen
    alpab_pm6(26, 9) =     4.294707d0 !        Iron -     Fluorine
    xab_pm6  (26, 9) =     3.657350d0 !        Iron -     Fluorine
    alpab_pm6(26,15) =     2.567534d0 !        Iron -   Phosphorus
    xab_pm6  (26,15) =     0.431291d0 !        Iron -   Phosphorus
    alpab_pm6(26,16) =     0.988991d0 !        Iron -       Sulfur
    xab_pm6  (26,16) =     0.033478d0 !        Iron -       Sulfur
    alpab_pm6(26,17) =     1.229793d0 !        Iron -     Chlorine
    xab_pm6  (26,17) =     0.019473d0 !        Iron -     Chlorine
    alpab_pm6(26,19) =     2.000000d0 !        Iron -    Potassium
    xab_pm6  (26,19) =     6.000000d0 !        Iron -    Potassium
    alpab_pm6(26,26) =     2.720785d0 !        Iron -         Iron
    xab_pm6  (26,26) =     1.846890d0 !        Iron -         Iron
! Co
    alpab_pm6(27, 1) =     2.966518d0 !      Cobalt -     Hydrogen
    xab_pm6  (27, 1) =     2.472465d0 !      Cobalt -     Hydrogen
    alpab_pm6(27, 6) =     3.716233d0 !      Cobalt -       Carbon
    xab_pm6  (27, 6) =     2.123930d0 !      Cobalt -       Carbon
    alpab_pm6(27, 7) =     3.618638d0 !      Cobalt -     Nitrogen
    xab_pm6  (27, 7) =     2.653836d0 !      Cobalt -     Nitrogen
    alpab_pm6(27, 8) =     3.726911d0 !      Cobalt -       Oxygen
    xab_pm6  (27, 8) =     5.252022d0 !      Cobalt -       Oxygen
    alpab_pm6(27, 9) =     3.956347d0 !      Cobalt -     Fluorine
    xab_pm6  (27, 9) =     4.585030d0 !      Cobalt -     Fluorine
    alpab_pm6(27,14) =     2.469805d0 !      Cobalt -      Silicon
    xab_pm6  (27,14) =     1.090240d0 !      Cobalt -      Silicon
    alpab_pm6(27,15) =     1.152505d0 !      Cobalt -   Phosphorus
    xab_pm6  (27,15) =     0.105936d0 !      Cobalt -   Phosphorus
    alpab_pm6(27,16) =     2.429255d0 !      Cobalt -       Sulfur
    xab_pm6  (27,16) =     0.436707d0 !      Cobalt -       Sulfur
    alpab_pm6(27,17) =     3.217497d0 !      Cobalt -     Chlorine
    xab_pm6  (27,17) =     1.033414d0 !      Cobalt -     Chlorine
    alpab_pm6(27,27) =     3.288166d0 !      Cobalt -       Cobalt
    xab_pm6  (27,27) =     3.919618d0 !      Cobalt -       Cobalt
! Ni
    alpab_pm6(28, 1) =     2.635280d0 !      Nickel -     Hydrogen
    xab_pm6  (28, 1) =     1.763124d0 !      Nickel -     Hydrogen
    alpab_pm6(28, 6) =     4.285513d0 !      Nickel -       Carbon
    xab_pm6  (28, 6) =     7.133324d0 !      Nickel -       Carbon
    alpab_pm6(28, 7) =     3.845215d0 !      Nickel -     Nitrogen
    xab_pm6  (28, 7) =     4.286800d0 !      Nickel -     Nitrogen
    alpab_pm6(28, 8) =     2.937232d0 !      Nickel -       Oxygen
    xab_pm6  (28, 8) =     0.885942d0 !      Nickel -       Oxygen
    alpab_pm6(28, 9) =     3.440241d0 !      Nickel -     Fluorine
    xab_pm6  (28, 9) =     1.088208d0 !      Nickel -     Fluorine
    alpab_pm6(28,14) =     2.068881d0 !      Nickel -      Silicon
    xab_pm6  (28,14) =     0.938646d0 !      Nickel -      Silicon
    alpab_pm6(28,15) =     3.260283d0 !      Nickel -   Phosphorus
    xab_pm6  (28,15) =     5.059727d0 !      Nickel -   Phosphorus
    alpab_pm6(28,16) =     2.002752d0 !      Nickel -       Sulfur
    xab_pm6  (28,16) =     0.274852d0 !      Nickel -       Sulfur
    alpab_pm6(28,17) =     2.200512d0 !      Nickel -     Chlorine
    xab_pm6  (28,17) =     0.202313d0 !      Nickel -     Chlorine
    alpab_pm6(28,28) =     1.097960d0 !      Nickel -       Cobalt
    xab_pm6  (28,28) =     0.035474d0 !      Nickel -       Cobalt
! Cu
    alpab_pm6(29, 1) =     2.335359d0 !      Copper -     Hydrogen
    xab_pm6  (29, 1) =     0.603591d0 !      Copper -     Hydrogen
    alpab_pm6(29, 6) =     4.638773d0 !      Copper -       Carbon
    xab_pm6  (29, 6) =     7.067794d0 !      Copper -       Carbon
    alpab_pm6(29, 7) =     4.214337d0 !      Copper -     Nitrogen
    xab_pm6  (29, 7) =     3.228667d0 !      Copper -     Nitrogen
    alpab_pm6(29, 8) =     3.959951d0 !      Copper -       Oxygen
    xab_pm6  (29, 8) =     2.000000d0 !      Copper -       Oxygen
    alpab_pm6(29, 9) =     4.478832d0 !      Copper -     Fluorine
    xab_pm6  (29, 9) =     1.282108d0 !      Copper -     Fluorine
    alpab_pm6(29,15) =     0.210640d0 !      Copper -   Phosphorus
    xab_pm6  (29,15) =     0.020126d0 !      Copper -   Phosphorus
    alpab_pm6(29,16) =     0.273112d0 !      Copper -       Sulfur
    xab_pm6  (29,16) =     0.005248d0 !      Copper -       Sulfur
    alpab_pm6(29,17) =     2.776531d0 !      Copper -     Chlorine
    xab_pm6  (29,17) =     0.139065d0 !      Copper -     Chlorine
    alpab_pm6(29,29) =     3.616846d0 !      Copper -       Cobalt
    xab_pm6  (29,29) =     5.184376d0 !      Copper -       Cobalt
!Zn
    alpab_pm6(30, 1) =     1.987891d0 !        Zinc -     Hydrogen
    xab_pm6  (30, 1) =     3.109193d0 !        Zinc -     Hydrogen
    alpab_pm6(30, 6) =     1.802327d0 !        Zinc -       Carbon
    xab_pm6  (30, 6) =     0.991465d0 !        Zinc -       Carbon
    alpab_pm6(30, 7) =     1.844579d0 !        Zinc -     Nitrogen
    xab_pm6  (30, 7) =     0.952476d0 !        Zinc -     Nitrogen
    alpab_pm6(30, 8) =     2.335054d0 !        Zinc -       Oxygen
    xab_pm6  (30, 8) =     2.265313d0 !        Zinc -       Oxygen
    alpab_pm6(30, 9) =     2.410021d0 !        Zinc -     Fluorine
    xab_pm6  (30, 9) =     1.225545d0 !        Zinc -     Fluorine
    alpab_pm6(30,14) =     1.832058d0 !        Zinc -      Silicon
    xab_pm6  (30,14) =     3.783905d0 !        Zinc -      Silicon
    alpab_pm6(30,15) =     1.220480d0 !        Zinc -   Phosphorus
    xab_pm6  (30,15) =     0.581530d0 !        Zinc -   Phosphorus
    alpab_pm6(30,16) =     1.455000d0 !        Zinc -       Sulfur
    xab_pm6  (30,16) =     0.648000d0 !        Zinc -       Sulfur
    alpab_pm6(30,17) =     1.625176d0 !        Zinc -     Chlorine
    xab_pm6  (30,17) =     0.721351d0 !        Zinc -     Chlorine
    alpab_pm6(30,20) =     1.119180d0 !        Zinc -      Calcium
    xab_pm6  (30,20) =     1.240290d0 !        Zinc -      Calcium
    alpab_pm6(30,30) =     0.929000d0 !        Zinc -         Zinc
    xab_pm6  (30,30) =     0.465000d0 !        Zinc -         Zinc
!Ga
    alpab_pm6(31, 1) =     1.847350d0 !     Gallium -     Hydrogen
    xab_pm6  (31, 1) =     1.386652d0 !     Gallium -     Hydrogen
    alpab_pm6(31, 6) =     2.325410d0 !     Gallium -       Carbon
    xab_pm6  (31, 6) =     1.962990d0 !     Gallium -       Carbon
    alpab_pm6(31, 7) =     2.121820d0 !     Gallium -     Nitrogen
    xab_pm6  (31, 7) =     1.188338d0 !     Gallium -     Nitrogen
    alpab_pm6(31, 8) =     2.348347d0 !     Gallium -       Oxygen
    xab_pm6  (31, 8) =     1.523644d0 !     Gallium -       Oxygen
    alpab_pm6(31, 9) =     2.679869d0 !     Gallium -     Fluorine
    xab_pm6  (31, 9) =     1.416942d0 !     Gallium -     Fluorine
    alpab_pm6(31,14) =     1.913780d0 !     Gallium -      Silicon
    xab_pm6  (31,14) =     1.002290d0 !     Gallium -      Silicon
    alpab_pm6(31,15) =     2.979650d0 !     Gallium -   Phosphorus
    xab_pm6  (31,15) =     0.500000d0 !     Gallium -   Phosphorus
    alpab_pm6(31,16) =     2.232108d0 !     Gallium -       Sulfur
    xab_pm6  (31,16) =     2.456284d0 !     Gallium -       Sulfur
    alpab_pm6(31,17) =     2.024710d0 !     Gallium -     Chlorine
    xab_pm6  (31,17) =     1.186661d0 !     Gallium -     Chlorine
    alpab_pm6(31,31) =     1.334643d0 !     Gallium -      Gallium
    xab_pm6  (31,31) =     1.198394d0 !     Gallium -      Gallium
! Ge
    alpab_pm6(32, 1) =     2.206793d0 !   Germanium -     Hydrogen
    xab_pm6  (32, 1) =     1.733226d0 !   Germanium -     Hydrogen
    alpab_pm6(32, 6) =     2.257469d0 !   Germanium -       Carbon
    xab_pm6  (32, 6) =     1.297510d0 !   Germanium -       Carbon
    alpab_pm6(32, 7) =     1.988226d0 !   Germanium -     Nitrogen
    xab_pm6  (32, 7) =     0.637506d0 !   Germanium -     Nitrogen
    alpab_pm6(32, 8) =     2.139413d0 !   Germanium -       Oxygen
    xab_pm6  (32, 8) =     0.826964d0 !   Germanium -       Oxygen
    alpab_pm6(32, 9) =     2.384777d0 !   Germanium -     Fluorine
    xab_pm6  (32, 9) =     0.651977d0 !   Germanium -     Fluorine
    alpab_pm6(32,14) =     0.299721d0 !   Germanium -      Silicon
    xab_pm6  (32,14) =     0.178680d0 !   Germanium -      Silicon
    alpab_pm6(32,15) =     2.469291d0 !   Germanium -   Phosphorus
    xab_pm6  (32,15) =     5.616349d0 !   Germanium -   Phosphorus
    alpab_pm6(32,16) =     2.024588d0 !   Germanium -       Sulfur
    xab_pm6  (32,16) =     1.160957d0 !   Germanium -       Sulfur
    alpab_pm6(32,17) =     1.771228d0 !   Germanium -     Chlorine
    xab_pm6  (32,17) =     0.545239d0 !   Germanium -     Chlorine
    alpab_pm6(32,25) =     2.382834d0 !   Germanium -    Manganese
    xab_pm6  (32,25) =     2.255151d0 !   Germanium -    Manganese
    alpab_pm6(32,27) =     2.852610d0 !   Germanium -       Cobalt
    xab_pm6  (32,27) =     2.151850d0 !   Germanium -       Cobalt
    alpab_pm6(32,32) =     2.019000d0 !   Germanium -    Germanium
    xab_pm6  (32,32) =     3.023000d0 !   Germanium -    Germanium
! As
    alpab_pm6(33, 1) =     1.993527d0 !     Arsenic -     Hydrogen
    xab_pm6  (33, 1) =     1.090589d0 !     Arsenic -     Hydrogen
    alpab_pm6(33, 6) =     1.855069d0 !     Arsenic -       Carbon
    xab_pm6  (33, 6) =     0.579098d0 !     Arsenic -       Carbon
    alpab_pm6(33, 7) =     1.496543d0 !     Arsenic -     Nitrogen
    xab_pm6  (33, 7) =     0.273337d0 !     Arsenic -     Nitrogen
    alpab_pm6(33, 8) =     2.003950d0 !     Arsenic -       Oxygen
    xab_pm6  (33, 8) =     0.701614d0 !     Arsenic -       Oxygen
    alpab_pm6(33, 9) =     2.012583d0 !     Arsenic -     Fluorine
    xab_pm6  (33, 9) =     0.402628d0 !     Arsenic -     Fluorine
    alpab_pm6(33,13) =     1.152786d0 !     Arsenic -     Aluminum
    xab_pm6  (33,13) =     1.003580d0 !     Arsenic -     Aluminum
    alpab_pm6(33,14) =     1.915600d0 !     Arsenic -      Silicon
    xab_pm6  (33,14) =     1.430706d0 !     Arsenic -      Silicon
    alpab_pm6(33,16) =     1.954368d0 !     Arsenic -       Sulfur
    xab_pm6  (33,16) =     1.033784d0 !     Arsenic -       Sulfur
    alpab_pm6(33,17) =     1.691070d0 !     Arsenic -     Chlorine
    xab_pm6  (33,17) =     0.454433d0 !     Arsenic -     Chlorine
    alpab_pm6(33,22) =     1.932911d0 !     Arsenic -     Titanium
    xab_pm6  (33,22) =     1.581317d0 !     Arsenic -     Titanium
    alpab_pm6(33,27) =     3.368140d0 !     Arsenic -       Cobalt
    xab_pm6  (33,27) =     1.675240d0 !     Arsenic -       Cobalt
    alpab_pm6(33,30) =     1.459130d0 !     Arsenic -         Zinc
    xab_pm6  (33,30) =     3.156571d0 !     Arsenic -         Zinc
    alpab_pm6(33,31) =     1.730977d0 !     Arsenic -      Gallium
    xab_pm6  (33,31) =     1.686298d0 !     Arsenic -      Gallium
    alpab_pm6(33,33) =     1.588264d0 !     Arsenic -      Arsenic
    xab_pm6  (33,33) =     0.737307d0 !     Arsenic -      Arsenic
! Se
    alpab_pm6(34, 1) =     2.035068d0 !    Selenium -     Hydrogen
    xab_pm6  (34, 1) =     0.847998d0 !    Selenium -     Hydrogen
    alpab_pm6(34, 6) =     2.387118d0 !    Selenium -       Carbon
    xab_pm6  (34, 6) =     1.114787d0 !    Selenium -       Carbon
    alpab_pm6(34, 7) =     1.937764d0 !    Selenium -     Nitrogen
    xab_pm6  (34, 7) =     0.482840d0 !    Selenium -     Nitrogen
    alpab_pm6(34, 8) =     2.484263d0 !    Selenium -       Oxygen
    xab_pm6  (34, 8) =     0.955161d0 !    Selenium -       Oxygen
    alpab_pm6(34, 9) =     2.302180d0 !    Selenium -     Fluorine
    xab_pm6  (34, 9) =     0.444806d0 !    Selenium -     Fluorine
    alpab_pm6(34,14) =     1.529817d0 !    Selenium -      Silicon
    xab_pm6  (34,14) =     0.518227d0 !    Selenium -      Silicon
    alpab_pm6(34,15) =     1.048183d0 !    Selenium -   Phosphorus
    xab_pm6  (34,15) =     0.292052d0 !    Selenium -   Phosphorus
    alpab_pm6(34,16) =     1.479606d0 !    Selenium -       Sulfur
    xab_pm6  (34,16) =     0.391721d0 !    Selenium -       Sulfur
    alpab_pm6(34,17) =     2.128861d0 !    Selenium -     Chlorine
    xab_pm6  (34,17) =     0.981067d0 !    Selenium -     Chlorine
    alpab_pm6(34,25) =     2.648038d0 !    Selenium -    Manganese
    xab_pm6  (34,25) =     2.180720d0 !    Selenium -    Manganese
    alpab_pm6(34,27) =     2.523450d0 !    Selenium -       Cobalt
    xab_pm6  (34,27) =     2.202410d0 !    Selenium -       Cobalt
    alpab_pm6(34,30) =     1.186242d0 !    Selenium -         Zinc
    xab_pm6  (34,30) =     0.511594d0 !    Selenium -         Zinc
    alpab_pm6(34,32) =     2.669057d0 !    Selenium -    Germanium
    xab_pm6  (34,32) =     5.872051d0 !    Selenium -    Germanium
    alpab_pm6(34,33) =     1.665280d0 !    Selenium -      Arsenic
    xab_pm6  (34,33) =     0.711261d0 !    Selenium -      Arsenic
    alpab_pm6(34,34) =     1.795894d0 !    Selenium -     Selenium
    xab_pm6  (34,34) =     0.821823d0 !    Selenium -     Selenium
! Br
    alpab_pm6(35, 1) =     2.192803d0 !     Bromine -     Hydrogen
    xab_pm6  (35, 1) =     0.850378d0 !     Bromine -     Hydrogen
    alpab_pm6(35, 2) =     2.128275d0 !     Bromine -       Helium
    xab_pm6  (35, 2) =     1.062043d0 !     Bromine -       Helium
    alpab_pm6(35, 3) =     2.074441d0 !     Bromine -      Lithium
    xab_pm6  (35, 3) =     1.858866d0 !     Bromine -      Lithium
    alpab_pm6(35, 4) =     2.367146d0 !     Bromine -    Beryllium
    xab_pm6  (35, 4) =     1.940933d0 !     Bromine -    Beryllium
    alpab_pm6(35, 5) =     2.307890d0 !     Bromine -        Boron
    xab_pm6  (35, 5) =     1.226420d0 !     Bromine -        Boron
    alpab_pm6(35, 6) =     2.015086d0 !     Bromine -       Carbon
    xab_pm6  (35, 6) =     0.570686d0 !     Bromine -       Carbon
    alpab_pm6(35, 7) =     4.224901d0 !     Bromine -     Nitrogen
    xab_pm6  (35, 7) =    30.000133d0 !     Bromine -     Nitrogen
    alpab_pm6(35, 8) =     2.283046d0 !     Bromine -       Oxygen
    xab_pm6  (35, 8) =     0.706584d0 !     Bromine -       Oxygen
    alpab_pm6(35, 9) =     2.031765d0 !     Bromine -     Fluorine
    xab_pm6  (35, 9) =     0.293500d0 !     Bromine -     Fluorine
    alpab_pm6(35,10) =     2.464172d0 !     Bromine -         Neon
    xab_pm6  (35,10) =     1.006159d0 !     Bromine -         Neon
    alpab_pm6(35,11) =     1.622218d0 !     Bromine -       Sodium
    xab_pm6  (35,11) =     1.752937d0 !     Bromine -       Sodium
    alpab_pm6(35,12) =     2.195697d0 !     Bromine -    Magnesium
    xab_pm6  (35,12) =     2.916280d0 !     Bromine -    Magnesium
    alpab_pm6(35,13) =     1.894141d0 !     Bromine -     Aluminum
    xab_pm6  (35,13) =     2.357130d0 !     Bromine -     Aluminum
    alpab_pm6(35,14) =     1.570825d0 !     Bromine -      Silicon
    xab_pm6  (35,14) =     0.589511d0 !     Bromine -      Silicon
    alpab_pm6(35,15) =     1.402139d0 !     Bromine -   Phosphorus
    xab_pm6  (35,15) =     0.456521d0 !     Bromine -   Phosphorus
    alpab_pm6(35,16) =     1.509874d0 !     Bromine -       Sulfur
    xab_pm6  (35,16) =     0.286688d0 !     Bromine -       Sulfur
    alpab_pm6(35,17) =     1.710331d0 !     Bromine -     Chlorine
    xab_pm6  (35,17) =     0.389238d0 !     Bromine -     Chlorine
    alpab_pm6(35,18) =     2.450801d0 !     Bromine -        Argon
    xab_pm6  (35,18) =     3.262668d0 !     Bromine -        Argon
    alpab_pm6(35,19) =     1.616093d0 !     Bromine -    Potassium
    xab_pm6  (35,19) =     3.322795d0 !     Bromine -    Potassium
    alpab_pm6(35,20) =     2.078405d0 !     Bromine -      Calcium
    xab_pm6  (35,20) =     4.052910d0 !     Bromine -      Calcium
    alpab_pm6(35,21) =     1.793486d0 !     Bromine -     Scandium
    xab_pm6  (35,21) =     2.098251d0 !     Bromine -     Scandium
    alpab_pm6(35,22) =     1.674847d0 !     Bromine -     Titanium
    xab_pm6  (35,22) =     0.883434d0 !     Bromine -     Titanium
    alpab_pm6(35,23) =     1.902904d0 !     Bromine -     Vanadium
    xab_pm6  (35,24) =     0.612698d0 !     Bromine -     Vanadium
    alpab_pm6(35,24) =     1.566028d0 !     Bromine -     Chromium
    xab_pm6  (35,24) =     0.217853d0 !     Bromine -     Chromium
    alpab_pm6(35,25) =     2.283820d0 !     Bromine -    Manganese
    xab_pm6  (35,25) =     1.183580d0 !     Bromine -    Manganese
    alpab_pm6(35,26) =     3.641782d0 !     Bromine -         Iron
    xab_pm6  (35,26) =     6.061921d0 !     Bromine -         Iron
    alpab_pm6(35,27) =     2.632688d0 !     Bromine -       Cobalt
    xab_pm6  (35,27) =     0.425148d0 !     Bromine -       Cobalt
    alpab_pm6(35,28) =     2.772136d0 !     Bromine -       Nickel
    xab_pm6  (35,28) =     0.632145d0 !     Bromine -       Nickel
    alpab_pm6(35,29) =     5.826407d0 !     Bromine -       Copper
    xab_pm6  (35,29) =     0.768517d0 !     Bromine -       Copper
    alpab_pm6(35,30) =     1.416120d0 !     Bromine -         Zinc
    xab_pm6  (35,30) =     0.747027d0 !     Bromine -         Zinc
    alpab_pm6(35,31) =     1.819105d0 !     Bromine -      Gallium
    xab_pm6  (35,31) =     1.261036d0 !     Bromine -      Gallium
    alpab_pm6(35,32) =     1.602366d0 !     Bromine -    Germanium
    xab_pm6  (35,32) =     0.627737d0 !     Bromine -    Germanium
    alpab_pm6(35,33) =     1.520170d0 !     Bromine -      Arsenic
    xab_pm6  (35,33) =     0.514153d0 !     Bromine -      Arsenic
    alpab_pm6(35,34) =     1.483713d0 !     Bromine -     Selenium
    xab_pm6  (35,34) =     0.319342d0 !     Bromine -     Selenium
    alpab_pm6(35,35) =     1.758146d0 !     Bromine -      Bromine
    xab_pm6  (35,35) =     0.615308d0 !     Bromine -      Bromine
! Kr
    alpab_pm6(36, 1) =     3.770453d0 !     Krypton -     Hydrogen
    xab_pm6  (36, 1) =     5.125897d0 !     Krypton -     Hydrogen
    alpab_pm6(36, 2) =     1.996943d0 !     Krypton -       Helium
    xab_pm6  (36, 2) =     0.627701d0 !     Krypton -       Helium
    alpab_pm6(36, 3) =     3.314562d0 !     Krypton -      Lithium
    xab_pm6  (36, 3) =     8.758697d0 !     Krypton -      Lithium
    alpab_pm6(36, 4) =     3.253048d0 !     Krypton -    Beryllium
    xab_pm6  (36, 4) =    10.237796d0 !     Krypton -    Beryllium
    alpab_pm6(36, 5) =     2.363169d0 !     Krypton -        Boron
    xab_pm6  (36, 5) =     2.946781d0 !     Krypton -        Boron
    alpab_pm6(36, 6) =     2.076738d0 !     Krypton -       Carbon
    xab_pm6  (36, 6) =     0.652623d0 !     Krypton -       Carbon
    alpab_pm6(36, 7) =     1.644052d0 !     Krypton -     Nitrogen
    xab_pm6  (36, 7) =     0.199606d0 !     Krypton -     Nitrogen
    alpab_pm6(36, 8) =     0.292300d0 !     Krypton -       Oxygen
    xab_pm6  (36, 8) =     0.006733d0 !     Krypton -       Oxygen
    alpab_pm6(36, 9) =     3.452321d0 !     Krypton -     Fluorine
    xab_pm6  (36, 9) =     4.134407d0 !     Krypton -     Fluorine
    alpab_pm6(36,10) =     2.813679d0 !     Krypton -         Neon
    xab_pm6  (36,10) =     1.433722d0 !     Krypton -         Neon
    alpab_pm6(36,11) =     2.480598d0 !     Krypton -       Sodium
    xab_pm6  (36,11) =     8.354448d0 !     Krypton -       Sodium
    alpab_pm6(36,12) =     1.391487d0 !     Krypton -    Magnesium
    xab_pm6  (36,12) =     0.888436d0 !     Krypton -    Magnesium
    alpab_pm6(36,13) =     2.467131d0 !     Krypton -    Aluminium
    xab_pm6  (36,13) =     5.091716d0 !     Krypton -    Aluminium
    alpab_pm6(36,14) =     1.764100d0 !     Krypton -      Silicon
    xab_pm6  (36,14) =     0.554250d0 !     Krypton -      Silicon
    alpab_pm6(36,17) =     1.884974d0 !     Krypton -     Chlorine
    xab_pm6  (36,17) =     0.520217d0 !     Krypton -     Chlorine
    alpab_pm6(36,18) =     1.995125d0 !     Krypton -        Argon
    xab_pm6  (36,18) =     0.554874d0 !     Krypton -        Argon
    alpab_pm6(36,19) =     2.182487d0 !     Krypton -    Potassium
    xab_pm6  (36,19) =     8.609782d0 !     Krypton -    Potassium
    alpab_pm6(36,20) =     1.305197d0 !     Krypton -      Calcium
    xab_pm6  (36,20) =     0.878891d0 !     Krypton -      Calcium
    alpab_pm6(36,35) =     1.529006d0 !     Krypton -      Bromine
    xab_pm6  (36,35) =     0.308098d0 !     Krypton -      Bromine
    alpab_pm6(36,36) =     1.135319d0 !     Krypton -      Krypton
    xab_pm6  (36,36) =     0.052099d0 !     Krypton -      Krypton
! Rb
    alpab_pm6(37, 1) =     2.443556d0 !    Rubidium -     Hydrogen
    xab_pm6  (37, 1) =    29.861632d0 !    Rubidium -     Hydrogen
    alpab_pm6(37, 2) =     1.270741d0 !    Rubidium -       Helium
    xab_pm6  (37, 2) =     1.862585d0 !    Rubidium -       Helium
    alpab_pm6(37, 5) =     5.532239d0 !    Rubidium -        Boron
    xab_pm6  (37, 5) =     9.040493d0 !    Rubidium -        Boron
    alpab_pm6(37, 6) =     2.765830d0 !    Rubidium -       Carbon
    xab_pm6  (37, 6) =    29.974031d0 !    Rubidium -       Carbon
    alpab_pm6(37, 7) =     0.761047d0 !    Rubidium -     Nitrogen
    xab_pm6  (37, 7) =     0.024636d0 !    Rubidium -     Nitrogen
    alpab_pm6(37, 8) =     1.334908d0 !    Rubidium -       Oxygen
    xab_pm6  (37, 8) =     1.125350d0 !    Rubidium -       Oxygen
    alpab_pm6(37, 9) =     3.638122d0 !    Rubidium -     Fluorine
    xab_pm6  (37, 9) =    28.815278d0 !    Rubidium -     Fluorine
    alpab_pm6(37,10) =     2.267591d0 !    Rubidium -         Neon
    xab_pm6  (37,10) =     7.736563d0 !    Rubidium -         Neon
    alpab_pm6(37,13) =     0.798774d0 !    Rubidium -    Aluminium
    xab_pm6  (37,13) =     2.992457d0 !    Rubidium -    Aluminium
    alpab_pm6(37,16) =     1.303184d0 !    Rubidium -       Sulfur
    xab_pm6  (37,16) =     0.964411d0 !    Rubidium -       Sulfur
    alpab_pm6(37,17) =     2.274411d0 !    Rubidium -     Chlorine
    xab_pm6  (37,17) =    10.384486d0 !    Rubidium -     Chlorine
    alpab_pm6(37,18) =     2.510977d0 !    Rubidium -        Argon
    xab_pm6  (37,18) =    18.433329d0 !    Rubidium -        Argon
    alpab_pm6(37,35) =     1.797766d0 !    Rubidium -      Bromine
    xab_pm6  (37,35) =     5.176214d0 !    Rubidium -      Bromine
    alpab_pm6(37,36) =     2.268753d0 !    Rubidium -      Krypton
    xab_pm6  (37,36) =    15.307503d0 !    Rubidium -      Krypton
    alpab_pm6(37,37) =     1.180818d0 !    Rubidium -     Rubidium
    xab_pm6  (37,37) =    20.147610d0 !    Rubidium -     Rubidium
! Sr
    alpab_pm6(38, 1) =     2.105914d0 !   Strontium -     Hydrogen
    xab_pm6  (38, 1) =    12.973316d0 !   Strontium -     Hydrogen
    alpab_pm6(38, 6) =     1.986688d0 !   Strontium -       Carbon
    xab_pm6  (38, 6) =     6.654657d0 !   Strontium -       Carbon
    alpab_pm6(38, 7) =     2.183629d0 !   Strontium -     Nitrogen
    xab_pm6  (38, 7) =     6.853866d0 !   Strontium -     Nitrogen
    alpab_pm6(38, 8) =     2.138399d0 !   Strontium -       Oxygen
    xab_pm6  (38, 8) =     3.561396d0 !   Strontium -       Oxygen
    alpab_pm6(38, 9) =     3.050666d0 !   Strontium -     Fluorine
    xab_pm6  (38, 9) =    10.971705d0 !   Strontium -     Fluorine
    alpab_pm6(38,14) =     2.969780d0 !   Strontium -      Silicon
    xab_pm6  (38,14) =     2.764750d0 !   Strontium -      Silicon
    alpab_pm6(38,15) =     2.789150d0 !   Strontium -   Phosphorus
    xab_pm6  (38,15) =     2.552100d0 !   Strontium -   Phosphorus
    alpab_pm6(38,16) =     1.598106d0 !   Strontium -       Sulfur
    xab_pm6  (38,16) =     3.129603d0 !   Strontium -       Sulfur
    alpab_pm6(38,17) =     1.854190d0 !   Strontium -     Chlorine
    xab_pm6  (38,17) =     3.783955d0 !   Strontium -     Chlorine
    alpab_pm6(38,22) =     2.880030d0 !   Strontium -     Titanium
    xab_pm6  (38,22) =     2.817250d0 !   Strontium -     Titanium
    alpab_pm6(38,35) =     1.524316d0 !   Strontium -      Bromine
    xab_pm6  (38,35) =     2.766567d0 !   Strontium -      Bromine
    alpab_pm6(38,38) =     1.000040d0 !   Strontium -    Strontium
    xab_pm6  (38,38) =     5.372120d0 !   Strontium -    Strontium
! Cd
    alpab_pm6(48, 1) =     2.628748d0 !     Cadmium -     Hydrogen
    xab_pm6  (48, 1) =    11.914201d0 !     Cadmium -     Hydrogen
    alpab_pm6(48, 6) =     1.425678d0 !     Cadmium -       Carbon
    xab_pm6  (48, 6) =     0.603441d0 !     Cadmium -       Carbon
    alpab_pm6(48, 7) =     0.970423d0 !     Cadmium -     Nitrogen
    xab_pm6  (48, 7) =     0.180663d0 !     Cadmium -     Nitrogen
    alpab_pm6(48, 8) =     1.696673d0 !     Cadmium -       Oxygen
    xab_pm6  (48, 8) =     0.926146d0 !     Cadmium -       Oxygen
    alpab_pm6(48, 9) =     2.312135d0 !     Cadmium -     Fluorine
    xab_pm6  (48, 9) =     1.353665d0 !     Cadmium -     Fluorine
    alpab_pm6(48,14) =     1.371225d0 !     Cadmium -      Silicon
    xab_pm6  (48,14) =     2.253346d0 !     Cadmium -      Silicon
    alpab_pm6(48,16) =     1.182202d0 !     Cadmium -       Sulfur
    xab_pm6  (48,16) =     0.361389d0 !     Cadmium -       Sulfur
    alpab_pm6(48,17) =     0.943547d0 !     Cadmium -     Chlorine
    xab_pm6  (48,17) =     0.140424d0 !     Cadmium -     Chlorine
    alpab_pm6(48,35) =     1.001451d0 !     Cadmium -      Bromine
    xab_pm6  (48,35) =     0.272267d0 !     Cadmium -      Bromine
    alpab_pm6(48,48) =     1.564044d0 !     Cadmium -      Cadmium
    xab_pm6  (48,48) =    18.617999d0 !     Cadmium -      Cadmium
! In
    alpab_pm6(49, 1) =     3.064144d0 !      Indium -     Hydrogen
    xab_pm6  (49, 1) =    14.975293d0 !      Indium -     Hydrogen
    alpab_pm6(49, 6) =     2.189272d0 !      Indium -       Carbon
    xab_pm6  (49, 6) =     2.187385d0 !      Indium -       Carbon
    alpab_pm6(49, 7) =     2.469868d0 !      Indium -     Nitrogen
    xab_pm6  (49, 7) =     3.369993d0 !      Indium -     Nitrogen
    alpab_pm6(49, 8) =     2.662095d0 !      Indium -       Oxygen
    xab_pm6  (49, 8) =     4.128583d0 !      Indium -       Oxygen
    alpab_pm6(49, 9) =     2.948797d0 !      Indium -     Fluorine
    xab_pm6  (49, 9) =     3.701016d0 !      Indium -     Fluorine
    alpab_pm6(49,16) =     2.542131d0 !      Indium -       Sulfur
    xab_pm6  (49,16) =     6.341105d0 !      Indium -       Sulfur
    alpab_pm6(49,17) =     2.233405d0 !      Indium -     Chlorine
    xab_pm6  (49,17) =     2.388552d0 !      Indium -     Chlorine
    alpab_pm6(49,31) =     1.628870d0 !      Indium -      Gallium
    xab_pm6  (49,31) =     2.421987d0 !      Indium -      Gallium
    alpab_pm6(49,33) =     2.299552d0 !      Indium -      Arsenic
    xab_pm6  (49,33) =     6.208350d0 !      Indium -      Arsenic
    alpab_pm6(49,34) =     1.906572d0 !      Indium -     Selenium
    xab_pm6  (49,34) =     2.319323d0 !      Indium -     Selenium
    alpab_pm6(49,35) =     2.257957d0 !      Indium -      Bromine
    xab_pm6  (49,35) =     3.728598d0 !      Indium -      Bromine
    alpab_pm6(49,49) =     2.073241d0 !      Indium -       Indium
    xab_pm6  (49,49) =     8.063491d0 !      Indium -       Indium
! Sn
    alpab_pm6(50, 1) =     2.648910d0 !         Tin -     Hydrogen
    xab_pm6  (50, 1) =     6.535162d0 !         Tin -     Hydrogen
    alpab_pm6(50, 6) =     2.440538d0 !         Tin -       Carbon
    xab_pm6  (50, 6) =     3.374355d0 !         Tin -       Carbon
    alpab_pm6(50, 7) =     2.085589d0 !         Tin -     Nitrogen
    xab_pm6  (50, 7) =     1.391900d0 !         Tin -     Nitrogen
    alpab_pm6(50, 8) =     2.727260d0 !         Tin -       Oxygen
    xab_pm6  (50, 8) =     4.374017d0 !         Tin -       Oxygen
    alpab_pm6(50, 9) =     3.724286d0 !         Tin -     Fluorine
    xab_pm6  (50, 9) =    18.598664d0 !         Tin -     Fluorine
    alpab_pm6(50,16) =     2.131542d0 !         Tin -       Sulfur
    xab_pm6  (50,16) =     2.314870d0 !         Tin -       Sulfur
    alpab_pm6(50,17) =     1.771522d0 !         Tin -     Chlorine
    xab_pm6  (50,17) =     0.807782d0 !         Tin -     Chlorine
    alpab_pm6(50,32) =     2.524633d0 !         Tin -    Germanium
    xab_pm6  (50,32) =    12.343411d0 !         Tin -    Germanium
    alpab_pm6(50,34) =     2.127377d0 !         Tin -     Selenium
    xab_pm6  (50,34) =     3.061885d0 !         Tin -     Selenium
    alpab_pm6(50,35) =     1.535089d0 !         Tin -      Bromine
    xab_pm6  (50,35) =     0.668798d0 !         Tin -      Bromine
    alpab_pm6(50,50) =     0.921000d0 !         Tin -          Tin
    xab_pm6  (50,50) =     0.287000d0 !         Tin -          Tin
! Sb
    alpab_pm6(51, 1) =     1.571272d0 !    Antimony -     Hydrogen
    xab_pm6  (51, 1) =     0.795343d0 !    Antimony -     Hydrogen
    alpab_pm6(51, 6) =     1.696206d0 !    Antimony -       Carbon
    xab_pm6  (51, 6) =     0.579212d0 !    Antimony -       Carbon
    alpab_pm6(51, 7) =     0.676115d0 !    Antimony -     Nitrogen
    xab_pm6  (51, 7) =     0.082065d0 !    Antimony -     Nitrogen
    alpab_pm6(51, 8) =     1.846384d0 !    Antimony -       Oxygen
    xab_pm6  (51, 8) =     0.634234d0 !    Antimony -       Oxygen
    alpab_pm6(51, 9) =     2.182922d0 !    Antimony -     Fluorine
    xab_pm6  (51, 9) =     0.650277d0 !    Antimony -     Fluorine
    alpab_pm6(51,13) =     1.422641d0 !    Antimony -     Aluminum
    xab_pm6  (51,13) =     1.616690d0 !    Antimony -     Aluminum
    alpab_pm6(51,14) =     2.686590d0 !    Antimony -      Silicon
    xab_pm6  (51,14) =     8.713749d0 !    Antimony -      Silicon
    alpab_pm6(51,16) =     1.418837d0 !    Antimony -       Sulfur
    xab_pm6  (51,16) =     0.396969d0 !    Antimony -       Sulfur
    alpab_pm6(51,17) =     1.117287d0 !    Antimony -     Chlorine
    xab_pm6  (51,17) =     0.156475d0 !    Antimony -     Chlorine
    alpab_pm6(51,25) =     2.400320d0 !    Antimony -    Manganese
    xab_pm6  (51,25) =     2.236710d0 !    Antimony -    Manganese
    alpab_pm6(51,27) =     2.204630d0 !    Antimony -       Cobalt
    xab_pm6  (51,27) =     2.276050d0 !    Antimony -       Cobalt
    alpab_pm6(51,35) =     1.063916d0 !    Antimony -      Bromine
    xab_pm6  (51,35) =     0.198044d0 !    Antimony -      Bromine
    alpab_pm6(51,43) =     2.204850d0 !    Antimony -   Technetium
    xab_pm6  (51,43) =     2.276260d0 !    Antimony -   Technetium
    alpab_pm6(51,44) =     2.204350d0 !    Antimony -    Ruthenium
    xab_pm6  (51,44) =     2.275760d0 !    Antimony -    Ruthenium
    alpab_pm6(51,45) =     2.204930d0 !    Antimony -      Rhodium
    xab_pm6  (51,45) =     2.276340d0 !    Antimony -      Rhodium
    alpab_pm6(51,49) =     2.141933d0 !    Antimony -       Indium
    xab_pm6  (51,49) =     6.660801d0 !    Antimony -       Indium
    alpab_pm6(51,51) =     1.348535d0 !    Antimony -     Antimony
    xab_pm6  (51,51) =     0.724885d0 !    Antimony -     Antimony
! Te
    alpab_pm6(52, 1) =     2.039130d0 !   Tellurium -     Hydrogen
    xab_pm6  (52, 1) =     1.807679d0 !   Tellurium -     Hydrogen
    alpab_pm6(52, 6) =     1.992816d0 !   Tellurium -       Carbon
    xab_pm6  (52, 6) =     0.970494d0 !   Tellurium -       Carbon
    alpab_pm6(52, 7) =     1.722269d0 !   Tellurium -     Nitrogen
    xab_pm6  (52, 7) =     0.358593d0 !   Tellurium -     Nitrogen
    alpab_pm6(52, 8) =     1.853064d0 !   Tellurium -       Oxygen
    xab_pm6  (52, 8) =     0.382926d0 !   Tellurium -       Oxygen
    alpab_pm6(52, 9) =     1.998576d0 !   Tellurium -     Fluorine
    xab_pm6  (52, 9) =     2.106812d0 !   Tellurium -     Fluorine
    alpab_pm6(52,13) =     1.387541d0 !   Tellurium -     Aluminum
    xab_pm6  (52,13) =     2.106812d0 !   Tellurium -     Aluminum
    alpab_pm6(52,15) =     1.453718d0 !   Tellurium -   Phosphorus
    xab_pm6  (52,15) =     1.109289d0 !   Tellurium -   Phosphorus
    alpab_pm6(52,16) =     1.830170d0 !   Tellurium -       Sulfur
    xab_pm6  (52,16) =     0.943925d0 !   Tellurium -       Sulfur
    alpab_pm6(52,17) =     1.300260d0 !   Tellurium -     Chlorine
    xab_pm6  (52,17) =     0.285478d0 !   Tellurium -     Chlorine
    alpab_pm6(52,30) =     1.218929d0 !   Tellurium -         Zinc
    xab_pm6  (52,30) =     1.756070d0 !   Tellurium -         Zinc
    alpab_pm6(52,32) =     2.342372d0 !   Tellurium -    Germanium
    xab_pm6  (52,32) =     7.019049d0 !   Tellurium -    Germanium
    alpab_pm6(52,33) =     1.189253d0 !   Tellurium -      Arsenic
    xab_pm6  (52,33) =     0.685774d0 !   Tellurium -      Arsenic
    alpab_pm6(52,34) =     1.566008d0 !   Tellurium -     Selenium
    xab_pm6  (52,34) =     1.187826d0 !   Tellurium -     Selenium
    alpab_pm6(52,35) =     1.250940d0 !   Tellurium -      Bromine
    xab_pm6  (52,35) =     0.394202d0 !   Tellurium -      Bromine
    alpab_pm6(52,48) =     1.307262d0 !   Tellurium -      Cadmium
    xab_pm6  (52,48) =     1.085919d0 !   Tellurium -      Cadmium
    alpab_pm6(52,49) =     1.540988d0 !   Tellurium -       Indium
    xab_pm6  (52,49) =     2.039582d0 !   Tellurium -       Indium
    alpab_pm6(52,50) =     1.763941d0 !   Tellurium -          Tin
    xab_pm6  (52,50) =     2.951976d0 !   Tellurium -          Tin
    alpab_pm6(52,52) =     1.164978d0 !   Tellurium -    Tellurium
    xab_pm6  (52,52) =     0.642486d0 !   Tellurium -    Tellurium
! Te
    alpab_pm6(53, 1) =     2.139913d0 !      Iodine -     Hydrogen
    xab_pm6  (53, 1) =     0.981898d0 !      Iodine -     Hydrogen
    alpab_pm6(53, 2) =     2.172984d0 !      Iodine -       Helium
    xab_pm6  (53, 2) =     1.630721d0 !      Iodine -       Helium
    alpab_pm6(53, 3) =     2.121251d0 !      Iodine -      Lithium
    xab_pm6  (53, 3) =     1.630721d0 !      Iodine -      Lithium
    alpab_pm6(53, 4) =     2.288023d0 !      Iodine -    Beryllium
    xab_pm6  (53, 4) =     2.351898d0 !      Iodine -    Beryllium
    alpab_pm6(53, 5) =     2.667605d0 !      Iodine -        Boron
    xab_pm6  (53, 5) =     3.161385d0 !      Iodine -        Boron
    alpab_pm6(53, 6) =     2.068710d0 !      Iodine -       Carbon
    xab_pm6  (53, 6) =     0.810156d0 !      Iodine -       Carbon
    alpab_pm6(53, 7) =     1.677518d0 !      Iodine -     Nitrogen
    xab_pm6  (53, 7) =     0.264903d0 !      Iodine -     Nitrogen
    alpab_pm6(53, 8) =     2.288919d0 !      Iodine -       Oxygen
    xab_pm6  (53, 8) =     0.866204d0 !      Iodine -       Oxygen
    alpab_pm6(53, 9) =     2.203580d0 !      Iodine -     Fluorine
    xab_pm6  (53, 9) =     0.392425d0 !      Iodine -     Fluorine
    alpab_pm6(53,10) =     2.414415d0 !      Iodine -         Neon
    xab_pm6  (53,10) =     1.503568d0 !      Iodine -         Neon
    alpab_pm6(53,11) =     1.403090d0 !      Iodine -       Sodium
    xab_pm6  (53,11) =     1.986112d0 !      Iodine -       Sodium
    alpab_pm6(53,12) =     2.045137d0 !      Iodine -    Magnesium
    xab_pm6  (53,12) =     3.276914d0 !      Iodine -    Magnesium
    alpab_pm6(53,13) =     1.816068d0 !      Iodine -     Aluminum
    xab_pm6  (53,13) =     2.929080d0 !      Iodine -     Aluminum
    alpab_pm6(53,14) =     1.559579d0 !      Iodine -      Silicon
    xab_pm6  (53,14) =     0.700299d0 !      Iodine -      Silicon
    alpab_pm6(53,15) =     2.131593d0 !      Iodine -   Phosphorus
    xab_pm6  (53,15) =     3.047207d0 !      Iodine -   Phosphorus
    alpab_pm6(53,16) =     1.855110d0 !      Iodine -       Sulfur
    xab_pm6  (53,16) =     0.709929d0 !      Iodine -       Sulfur
    alpab_pm6(53,17) =     1.574161d0 !      Iodine -     Chlorine
    xab_pm6  (53,17) =     0.310474d0 !      Iodine -     Chlorine
    alpab_pm6(53,18) =     1.576587d0 !      Iodine -        Argon
    xab_pm6  (53,18) =     0.305367d0 !      Iodine -        Argon
    alpab_pm6(53,19) =     1.539714d0 !      Iodine -    Potassium
    xab_pm6  (53,19) =     4.824353d0 !      Iodine -    Potassium
    alpab_pm6(53,20) =     2.196490d0 !      Iodine -      Calcium
    xab_pm6  (53,20) =     7.689921d0 !      Iodine -      Calcium
    alpab_pm6(53,21) =     1.814884d0 !      Iodine -     Scandium
    xab_pm6  (53,21) =     3.114282d0 !      Iodine -     Scandium
    alpab_pm6(53,22) =     1.933469d0 !      Iodine -     Titanium
    xab_pm6  (53,22) =     2.426747d0 !      Iodine -     Titanium
    alpab_pm6(53,23) =     2.683520d0 !      Iodine -     Vanadium
    xab_pm6  (53,24) =     6.198112d0 !      Iodine -     Vanadium
    alpab_pm6(53,24) =     2.634224d0 !      Iodine -     Chromium
    xab_pm6  (53,24) =     2.598590d0 !      Iodine -     Chromium
    alpab_pm6(53,25) =     2.266600d0 !      Iodine -    Manganese
    xab_pm6  (53,25) =     1.193410d0 !      Iodine -    Manganese
    alpab_pm6(53,26) =     1.912829d0 !      Iodine -         Iron
    xab_pm6  (53,26) =     0.532622d0 !      Iodine -         Iron
    alpab_pm6(53,27) =     3.235204d0 !      Iodine -       Cobalt
    xab_pm6  (53,27) =     1.105239d0 !      Iodine -       Cobalt
    alpab_pm6(53,28) =     1.085343d0 !      Iodine -       Nickel
    xab_pm6  (53,28) =     0.017459d0 !      Iodine -       Nickel
    alpab_pm6(53,29) =     0.834305d0 !      Iodine -       Copper
    xab_pm6  (53,29) =     0.006781d0 !      Iodine -       Copper
    alpab_pm6(53,30) =     1.394762d0 !      Iodine -         Zinc
    xab_pm6  (53,30) =     0.976607d0 !      Iodine -         Zinc
    alpab_pm6(53,31) =     1.671729d0 !      Iodine -      Gallium
    xab_pm6  (53,31) =     1.252168d0 !      Iodine -      Gallium
    alpab_pm6(53,32) =     1.817425d0 !      Iodine -    Germanium
    xab_pm6  (53,32) =     1.323267d0 !      Iodine -    Germanium
    alpab_pm6(53,33) =     1.245262d0 !      Iodine -      Arsenic
    xab_pm6  (53,33) =     0.310824d0 !      Iodine -      Arsenic
    alpab_pm6(53,35) =     1.579376d0 !      Iodine -      Bromine
    xab_pm6  (53,35) =     0.483054d0 !      Iodine -      Bromine
    alpab_pm6(53,36) =     1.238574d0 !      Iodine -      Krypton
    xab_pm6  (53,36) =     0.201136d0 !      Iodine -      Krypton
    alpab_pm6(53,37) =     1.432675d0 !      Iodine -     Rubidium
    xab_pm6  (53,37) =     4.092446d0 !      Iodine -     Rubidium
    alpab_pm6(53,38) =     1.262042d0 !      Iodine -    Strontium
    xab_pm6  (53,38) =     2.103941d0 !      Iodine -    Strontium
    alpab_pm6(53,39) =     1.279110d0 !      Iodine -      Yttrium
    xab_pm6  (53,39) =     1.021402d0 !      Iodine -      Yttrium
    alpab_pm6(53,40) =     1.995182d0 !      Iodine -    Zirconium
    xab_pm6  (53,40) =     4.513943d0 !      Iodine -    Zirconium 
    alpab_pm6(53,41) =     1.967251d0 !      Iodine -      Niobium
    xab_pm6  (53,41) =     2.399298d0 !      Iodine -      Niobium
    alpab_pm6(53,42) =     0.948461d0 !      Iodine -   Molybdenum 
    xab_pm6  (53,42) =     0.124695d0 !      Iodine -   Molybdenum
    alpab_pm6(53,43) =     1.292312d0 !      Iodine -   Technetium
    xab_pm6  (53,43) =     0.110594d0 !      Iodine -   Technetium 
    alpab_pm6(53,44) =     3.953203d0 !      Iodine -    Ruthenium
    xab_pm6  (53,44) =     7.837710d0 !      Iodine -    Ruthenium
    alpab_pm6(53,45) =     3.708170d0 !      Iodine -      Rhodium
    xab_pm6  (53,45) =     2.357944d0 !      Iodine -      Rhodium
    alpab_pm6(53,46) =     5.144544d0 !      Iodine -    Palladium
    xab_pm6  (53,46) =     3.522017d0 !      Iodine -    Palladium
    alpab_pm6(53,47) =     2.593161d0 !      Iodine -       Silver
    xab_pm6  (53,47) =     0.048904d0 !      Iodine -       Silver
    alpab_pm6(53,48) =     0.996238d0 !      Iodine -      Cadmium
    xab_pm6  (53,48) =     0.396784d0 !      Iodine -      Cadmium
    alpab_pm6(53,49) =     2.351758d0 !      Iodine -       Indium
    xab_pm6  (53,49) =     5.947821d0 !      Iodine -       Indium
    alpab_pm6(53,50) =     1.855633d0 !      Iodine -          Tin
    xab_pm6  (53,50) =     1.783163d0 !      Iodine -          Tin
    alpab_pm6(53,51) =     1.155315d0 !      Iodine -     Antimony
    xab_pm6  (53,51) =     0.318190d0 !      Iodine -     Antimony
    alpab_pm6(53,52) =     1.493951d0 !      Iodine -    Tellurium
    xab_pm6  (53,52) =     1.101116d0 !      Iodine -    Tellurium
    alpab_pm6(53,53) =     1.519925d0 !      Iodine -       Iodine
    xab_pm6  (53,53) =     0.510542d0 !      Iodine -       Iodine
! Xe
    alpab_pm6(54, 1) =     1.356861d0 !       Xenon -     Hydrogen
    xab_pm6  (54, 1) =     0.701016d0 !       Xenon -     Hydrogen
    alpab_pm6(54, 2) =     2.497832d0 !       Xenon -       Helium
    xab_pm6  (54, 2) =     2.599471d0 !       Xenon -       Helium
    alpab_pm6(54, 3) =     2.466895d0 !       Xenon -      Lithium
    xab_pm6  (54, 3) =     4.582081d0 !       Xenon -      Lithium
    alpab_pm6(54, 4) =     6.000003d0 !       Xenon -    Beryllium
    xab_pm6  (54, 4) =     0.660525d0 !       Xenon -    Beryllium
    alpab_pm6(54, 5) =     5.051957d0 !       Xenon -        Boron
    xab_pm6  (54, 5) =     1.100612d0 !       Xenon -        Boron
    alpab_pm6(54, 6) =     1.704440d0 !       Xenon -       Carbon
    xab_pm6  (54, 6) =     0.826727d0 !       Xenon -       Carbon
    alpab_pm6(54, 7) =     1.932952d0 !       Xenon -     Nitrogen
    xab_pm6  (54, 7) =     0.925624d0 !       Xenon -     Nitrogen
    alpab_pm6(54, 8) =     0.839233d0 !       Xenon -       Oxygen
    xab_pm6  (54, 8) =     0.035356d0 !       Xenon -       Oxygen
    alpab_pm6(54, 9) =     1.128812d0 !       Xenon -     Fluorine
    xab_pm6  (54, 9) =     0.065011d0 !       Xenon -     Fluorine
    alpab_pm6(54,10) =     1.330202d0 !       Xenon -         Neon
    xab_pm6  (54,10) =     0.293862d0 !       Xenon -         Neon
    alpab_pm6(54,11) =     2.103003d0 !       Xenon -       Sodium
    xab_pm6  (54,11) =     8.368204d0 !       Xenon -       Sodium
    alpab_pm6(54,12) =     2.698414d0 !       Xenon -    Magnesium
    xab_pm6  (54,12) =     9.723572d0 !       Xenon -    Magnesium
    alpab_pm6(54,13) =     2.412039d0 !       Xenon -    Aluminium
    xab_pm6  (54,13) =     7.404465d0 !       Xenon -    Aluminium
    alpab_pm6(54,14) =     3.087060d0 !       Xenon -      Silicon
    xab_pm6  (54,14) =    16.092000d0 !       Xenon -      Silicon
    alpab_pm6(54,17) =     1.546396d0 !       Xenon -     Chlorine
    xab_pm6  (54,17) =     0.463758d0 !       Xenon -     Chlorine
    alpab_pm6(54,18) =     0.591520d0 !       Xenon -        Argon
    xab_pm6  (54,18) =     0.049266d0 !       Xenon -        Argon
    alpab_pm6(54,19) =     1.171250d0 !       Xenon -    Potassium
    xab_pm6  (54,19) =     1.224889d0 !       Xenon -    Potassium
    alpab_pm6(54,20) =     1.510653d0 !       Xenon -      Calcium
    xab_pm6  (54,20) =     1.717121d0 !       Xenon -      Calcium
    alpab_pm6(54,35) =     1.439618d0 !       Xenon -      Bromine
    xab_pm6  (54,35) =     0.475116d0 !       Xenon -      Bromine
    alpab_pm6(54,36) =     0.551561d0 !       Xenon -      Krypton
    xab_pm6  (54,36) =     0.049793d0 !       Xenon -      Krypton
    alpab_pm6(54,37) =     1.087823d0 !       Xenon -     Rubidium
    xab_pm6  (54,37) =     0.974965d0 !       Xenon -     Rubidium
    alpab_pm6(54,53) =     0.799155d0 !       Xenon -       Iodine
    xab_pm6  (54,53) =     0.112090d0 !       Xenon -       Iodine
    alpab_pm6(54,54) =     1.244762d0 !       Xenon -        Xenon
    xab_pm6  (54,54) =     0.344474d0 !       Xenon -        Xenon
! Cs
    alpab_pm6(55, 1) =     0.264882d0 !      Cesium -     Hydrogen
    xab_pm6  (55, 1) =     0.096901d0 !      Cesium -     Hydrogen
    alpab_pm6(55, 5) =     1.487110d0 !      Cesium -        Boron
    xab_pm6  (55, 5) =    10.392610d0 !      Cesium -        Boron
    alpab_pm6(55, 6) =     2.147104d0 !      Cesium -       Carbon
    xab_pm6  (55, 6) =    24.514623d0 !      Cesium -       Carbon
    alpab_pm6(55, 7) =     2.446532d0 !      Cesium -     Nitrogen
    xab_pm6  (55, 7) =    29.711077d0 !      Cesium -     Nitrogen
    alpab_pm6(55, 8) =     2.085139d0 !      Cesium -       Oxygen
    xab_pm6  (55, 8) =     8.176843d0 !      Cesium -       Oxygen
    alpab_pm6(55, 9) =     2.834100d0 !      Cesium -     Fluorine
    xab_pm6  (55, 9) =    22.233416d0 !      Cesium -     Fluorine
    alpab_pm6(55,15) =     2.924953d0 !      Cesium -   Phosphorus
    xab_pm6  (55,15) =     0.506512d0 !      Cesium -   Phosphorus
    alpab_pm6(55,16) =     0.289412d0 !      Cesium -       Sulfur
    xab_pm6  (55,16) =     0.091743d0 !      Cesium -       Sulfur
    alpab_pm6(55,17) =     1.673663d0 !      Cesium -     Chlorine
    xab_pm6  (55,17) =     4.531965d0 !      Cesium -     Chlorine
    alpab_pm6(55,35) =     1.167189d0 !      Cesium -      Bromine
    xab_pm6  (55,35) =     1.658427d0 !      Cesium -      Bromine
    alpab_pm6(55,53) =     0.919562d0 !      Cesium -       Iodine
    xab_pm6  (55,53) =     1.072178d0 !      Cesium -       Iodine
    alpab_pm6(55,55) =     1.170843d0 !      Cesium -       Cesium
    xab_pm6  (55,55) =    25.320055d0 !      Cesium -       Cesium
! Ba
    alpab_pm6(56, 1) =     6.000135d0 !      Barium -     Hydrogen
    xab_pm6  (56, 1) =     2.040004d0 !      Barium -     Hydrogen
    alpab_pm6(56, 6) =     0.770626d0 !      Barium -       Carbon
    xab_pm6  (56, 6) =     0.119793d0 !      Barium -       Carbon
    alpab_pm6(56, 7) =     1.148233d0 !      Barium -     Nitrogen
    xab_pm6  (56, 7) =     0.207934d0 !      Barium -     Nitrogen
    alpab_pm6(56, 8) =     1.283018d0 !      Barium -       Oxygen
    xab_pm6  (56, 8) =     0.348945d0 !      Barium -       Oxygen
    alpab_pm6(56, 9) =     3.000618d0 !      Barium -     Fluorine
    xab_pm6  (56, 9) =     5.575255d0 !      Barium -     Fluorine
    alpab_pm6(56,13) =     2.105924d0 !      Barium -    Aluminium
    xab_pm6  (56,13) =     9.539099d0 !      Barium -    Aluminium
    alpab_pm6(56,14) =     1.240420d0 !      Barium -      Silicon
    xab_pm6  (56,14) =     1.212660d0 !      Barium -      Silicon
    alpab_pm6(56,16) =     0.705188d0 !      Barium -       Sulfur
    xab_pm6  (56,16) =     0.215386d0 !      Barium -       Sulfur
    alpab_pm6(56,17) =     1.071044d0 !      Barium -     Chlorine
    xab_pm6  (56,17) =     0.160177d0 !      Barium -     Chlorine
    alpab_pm6(56,22) =     2.176040d0 !      Barium -     Titanium
    xab_pm6  (56,22) =     9.493530d0 !      Barium -     Titanium
    alpab_pm6(56,35) =     1.190346d0 !      Barium -      Bromine
    xab_pm6  (56,35) =     0.828794d0 !      Barium -      Bromine
    alpab_pm6(56,53) =     0.982528d0 !      Barium -       Iodine
    xab_pm6  (56,53) =     0.835597d0 !      Barium -       Iodine
    alpab_pm6(56,56) =     0.339269d0 !      Barium -       Barium
    xab_pm6  (56,56) =     0.356186d0 !      Barium -       Barium
! Hg
    alpab_pm6(80, 1) =     1.136587d0 !     Mercury -     Hydrogen
    xab_pm6  (80, 1) =     0.799399d0 !     Mercury -     Hydrogen
    alpab_pm6(80, 6) =     0.795816d0 !     Mercury -       Carbon
    xab_pm6  (80, 6) =     0.147128d0 !     Mercury -       Carbon
    alpab_pm6(80, 7) =     0.332152d0 !     Mercury -     Nitrogen
    xab_pm6  (80, 7) =     0.050240d0 !     Mercury -     Nitrogen
    alpab_pm6(80, 8) =     1.052145d0 !     Mercury -       Oxygen
    xab_pm6  (80, 8) =     0.240720d0 !     Mercury -       Oxygen
    alpab_pm6(80, 9) =     1.240572d0 !     Mercury -     Fluorine
    xab_pm6  (80, 9) =     0.113827d0 !     Mercury -     Fluorine
    alpab_pm6(80,14) =     2.770860d0 !     Mercury -      Silicon
    xab_pm6  (80,14) =     3.680740d0 !     Mercury -      Silicon
    alpab_pm6(80,15) =     0.608604d0 !     Mercury -   Phosphorus
    xab_pm6  (80,15) =     0.214951d0 !     Mercury -   Phosphorus
    alpab_pm6(80,16) =     1.041682d0 !     Mercury -       Sulfur
    xab_pm6  (80,16) =     0.347383d0 !     Mercury -       Sulfur
    alpab_pm6(80,17) =     0.430731d0 !     Mercury -     Chlorine
    xab_pm6  (80,17) =     0.053660d0 !     Mercury -     Chlorine
    alpab_pm6(80,22) =     3.414630d0 !     Mercury -     Titanium
    xab_pm6  (80,22) =     2.957200d0 !     Mercury -     Titanium
    alpab_pm6(80,35) =     0.638717d0 !     Mercury -      Bromine
    xab_pm6  (80,35) =     0.172363d0 !     Mercury -      Bromine
    alpab_pm6(80,52) =     0.291500d0 !     Mercury -    Tellurium
    xab_pm6  (80,52) =     0.212732d0 !     Mercury -    Tellurium
    alpab_pm6(80,53) =     0.758162d0 !     Mercury -       Iodine
    xab_pm6  (80,53) =     0.342058d0 !     Mercury -       Iodine
    alpab_pm6(80,80) =     0.474413d0 !     Mercury -      Mercury
    xab_pm6  (80,80) =     0.423276d0 !     Mercury -      Mercury
! Tl
    alpab_pm6(81, 1) =     0.673658d0 !    Thallium -     Hydrogen
    xab_pm6  (81, 1) =     0.138205d0 !    Thallium -     Hydrogen
    alpab_pm6(81, 5) =     1.528347d0 !    Thallium -        Boron
    xab_pm6  (81, 5) =    10.504338d0 !    Thallium -        Boron
    alpab_pm6(81, 6) =     1.390345d0 !    Thallium -       Carbon
    xab_pm6  (81, 6) =     0.582895d0 !    Thallium -       Carbon
    alpab_pm6(81, 7) =     0.982335d0 !    Thallium -     Nitrogen
    xab_pm6  (81, 7) =     0.158812d0 !    Thallium -     Nitrogen
    alpab_pm6(81, 8) =     1.550068d0 !    Thallium -       Oxygen
    xab_pm6  (81, 8) =     0.636906d0 !    Thallium -       Oxygen
    alpab_pm6(81, 9) =     1.469516d0 !    Thallium -     Fluorine
    xab_pm6  (81, 9) =     0.226166d0 !    Thallium -     Fluorine
    alpab_pm6(81,16) =     0.994851d0 !    Thallium -       Sulfur
    xab_pm6  (81,16) =     0.303426d0 !    Thallium -       Sulfur
    alpab_pm6(81,17) =     0.846193d0 !    Thallium -     Chlorine
    xab_pm6  (81,17) =     0.162037d0 !    Thallium -     Chlorine
    alpab_pm6(81,35) =     0.874419d0 !    Thallium -      Bromine
    xab_pm6  (81,35) =     0.296836d0 !    Thallium -      Bromine
    alpab_pm6(81,53) =     0.902012d0 !    Thallium -       Iodine
    xab_pm6  (81,53) =     0.430033d0 !    Thallium -       Iodine
    alpab_pm6(81,81) =     1.191684d0 !    Thallium -     Thallium
    xab_pm6  (81,81) =     9.535127d0 !    Thallium -     Thallium
! Pb
    alpab_pm6(82, 1) =     1.522676d0 !        Lead -     Hydrogen
    xab_pm6  (82, 1) =     0.840096d0 !        Lead -     Hydrogen
    alpab_pm6(82, 3) =     1.001810d0 !        Lead -      Lithium
    xab_pm6  (82, 3) =     1.285064d0 !        Lead -      Lithium
    alpab_pm6(82, 5) =     0.911197d0 !        Lead -        Boron
    xab_pm6  (82, 5) =     1.138157d0 !        Lead -        Boron
    alpab_pm6(82, 6) =     1.525593d0 !        Lead -       Carbon
    xab_pm6  (82, 6) =     0.404656d0 !        Lead -       Carbon
    alpab_pm6(82, 7) =     1.317394d0 !        Lead -     Nitrogen
    xab_pm6  (82, 7) =     0.335787d0 !        Lead -     Nitrogen
    alpab_pm6(82, 8) =     1.763210d0 !        Lead -       Oxygen
    xab_pm6  (82, 8) =     0.782506d0 !        Lead -       Oxygen
    alpab_pm6(82, 9) =     3.288902d0 !        Lead -     Fluorine
    xab_pm6  (82, 9) =     8.368562d0 !        Lead -     Fluorine
    alpab_pm6(82,15) =     4.516800d0 !        Lead -   Phosphorus
    xab_pm6  (82,15) =     5.033200d0 !        Lead -   Phosphorus
    alpab_pm6(82,16) =     1.027519d0 !        Lead -       Sulfur
    xab_pm6  (82,16) =     0.175150d0 !        Lead -       Sulfur
    alpab_pm6(82,17) =     1.094123d0 !        Lead -     Chlorine
    xab_pm6  (82,17) =     0.164814d0 !        Lead -     Chlorine
    alpab_pm6(82,23) =     1.500000d0 !        Lead -     Vanadium
    xab_pm6  (82,23) =     1.000000d0 !        Lead -     Vanadium
    alpab_pm6(82,24) =     1.860760d0 !        Lead -     Chromium
    xab_pm6  (82,24) =     1.029110d0 !        Lead -     Chromium
    alpab_pm6(82,30) =     1.500000d0 !        Lead -         Zinc
    xab_pm6  (82,30) =     1.000000d0 !        Lead -         Zinc
    alpab_pm6(82,34) =     2.000000d0 !        Lead -     Selenium
    xab_pm6  (82,34) =     0.111195d0 !        Lead -     Selenium
    alpab_pm6(82,35) =     0.865550d0 !        Lead -      Bromine
    xab_pm6  (82,35) =     0.148229d0 !        Lead -      Bromine
    alpab_pm6(82,41) =     1.500000d0 !        Lead -      Niobium
    xab_pm6  (82,41) =     1.000000d0 !        Lead -      Niobium
    alpab_pm6(82,42) =     2.000000d0 !        Lead -   Molybdenum
    xab_pm6  (82,42) =     5.000000d0 !        Lead -   Molybdenum
    alpab_pm6(82,52) =     1.002559d0 !        Lead -    Tellurium
    xab_pm6  (82,52) =     0.809042d0 !        Lead -    Tellurium
    alpab_pm6(82,53) =     0.983474d0 !        Lead -       Iodine
    xab_pm6  (82,53) =     0.267426d0 !        Lead -       Iodine
    alpab_pm6(82,82) =     1.881764d0 !        Lead -         Lead
    xab_pm6  (82,82) =     2.362343d0 !        Lead -         Lead
! Bi
    alpab_pm6(83, 1) =     1.679905d0 !     Bismuth -     Hydrogen
    xab_pm6  (83, 1) =     1.397462d0 !     Bismuth -     Hydrogen
    alpab_pm6(83, 3) =     0.340140d0 !     Bismuth -      Lithium
    xab_pm6  (83, 3) =     0.695320d0 !     Bismuth -      Lithium
    alpab_pm6(83, 6) =     1.534025d0 !     Bismuth -       Carbon
    xab_pm6  (83, 6) =     0.576179d0 !     Bismuth -       Carbon
    alpab_pm6(83, 7) =     1.143876d0 !     Bismuth -     Nitrogen
    xab_pm6  (83, 7) =     0.152738d0 !     Bismuth -     Nitrogen
    alpab_pm6(83, 8) =     1.553297d0 !     Bismuth -       Oxygen
    xab_pm6  (83, 8) =     0.333042d0 !     Bismuth -       Oxygen
    alpab_pm6(83, 9) =     2.355400d0 !     Bismuth -     Fluorine
    xab_pm6  (83, 9) =     1.035324d0 !     Bismuth -     Fluorine
    alpab_pm6(83,16) =     1.466879d0 !     Bismuth -       Sulfur
    xab_pm6  (83,16) =     0.620997d0 !     Bismuth -       Sulfur
    alpab_pm6(83,17) =     1.272975d0 !     Bismuth -     Chlorine
    xab_pm6  (83,17) =     0.326871d0 !     Bismuth -     Chlorine
    alpab_pm6(83,34) =     1.344746d0 !     Bismuth -     Selenium
    xab_pm6  (83,34) =     0.651208d0 !     Bismuth -     Selenium
    alpab_pm6(83,35) =     1.146233d0 !     Bismuth -      Bromine
    xab_pm6  (83,35) =     0.381170d0 !     Bismuth -      Bromine
    alpab_pm6(83,53) =     1.302171d0 !     Bismuth -       Iodine
    xab_pm6  (83,53) =     0.862377d0 !     Bismuth -       Iodine
    alpab_pm6(83,83) =     1.074064d0 !     Bismuth -      Bismuth
    xab_pm6  (83,83) =     1.168214d0 !     Bismuth -      Bismuth

!---------------------------------------------
!END SECTION 3 - PM6 pair wise core core terms
!---------------------------------------------


!---------------------------------------------
!SECTION 4 - PM3MAIS  core core function terms
!---------------------------------------------
! Units: alpha(eV), beta(A**-2),gamma(A)
! unless otherwise stated, all parameters from Chem Phys Lett 330 (2000) 118
    element_supported_pm3mais(1)  = .true. ! H
    element_supported_pm3mais(8)  = .true. ! O
    element_supported_pm3mais(17) = .true. ! Cl
! H-H
    alpab_pm3mais( 1, 1, 1) = 0.009851d0
    alpab_pm3mais( 1, 1, 2) = 0.248658d0
    alpab_pm3mais( 1, 1, 3) = 0.195242d0
    betab_pm3mais( 1, 1, 1) = 1.376611d0
    betab_pm3mais( 1, 1, 2) = 0.812525d0
    betab_pm3mais( 1, 1, 3) = 10.759699d0
    gamab_pm3mais( 1, 1, 1) = 3.297249d0
    gamab_pm3mais( 1, 1, 2) = 1.170370d0
    gamab_pm3mais( 1, 1, 3) = 1.177800d0
! O-H for HCl dissociation from Theor Chem Acc 118 (2007) 425
    alpab_pm3mais( 8, 1, 1) = 1.0380054d0
    alpab_pm3mais( 8, 1, 2) = -1.162334d0
    alpab_pm3mais( 8, 1, 3) = -1.084972d0
    betab_pm3mais( 8, 1, 1) = 2.017471d0
    betab_pm3mais( 8, 1, 2) = 1.851292d0
    betab_pm3mais( 8, 1, 3) = 1.221409d0
    gamab_pm3mais( 8, 1, 1) = 2.184459d0
    gamab_pm3mais( 8, 1, 2) = 2.147501d0
    gamab_pm3mais( 8, 1, 3) = 0.280095d0
! O-H for protonated water from Chem Phys Lett, 2000,330,118
!      alpab_pm3mais( 8, 1, 1) = -1.947697d0
!      alpab_pm3mais( 8, 1, 2) = 0.415113d0
!      alpab_pm3mais( 8, 1, 3) = 0.804783d0
!      betab_pm3mais( 8, 1, 1) = 0.685841d0
!      betab_pm3mais( 8, 1, 2) = 4.586889d0
!      betab_pm3mais( 8, 1, 3) = 1.232103d0
!      gamab_pm3mais( 8, 1, 1) = 0.762936d0
!      gamab_pm3mais( 8, 1, 2) = 1.160914d0
!      gamab_pm3mais( 8, 1, 3) = 1.368463d0
! O-O
    alpab_pm3mais( 8, 8, 1) = -184.976043d0
    alpab_pm3mais( 8, 8, 2) = 188.095991d0
    alpab_pm3mais( 8, 8, 3) = 2.029969d0
    betab_pm3mais( 8, 8, 1) = 1.164736d0
    betab_pm3mais( 8, 8, 2) = 1.143302d0
    betab_pm3mais( 8, 8, 3) = 4.528650d0
    gamab_pm3mais( 8, 8, 1) = 0.641029d0
    gamab_pm3mais( 8, 8, 2) = 0.640013d0
    gamab_pm3mais( 8, 8, 3) = 1.365028d0
! Cl-H
    alpab_pm3mais( 17, 1, 1) = -6.731307d0
    alpab_pm3mais( 17, 1, 2) = 6.678218d0
    alpab_pm3mais( 17, 1, 3) = -0.250426d0
    betab_pm3mais( 17, 1, 1) = 5.429595d0
    betab_pm3mais( 17, 1, 2) = 5.892119d0
    betab_pm3mais( 17, 1, 3) = 0.434120d0
    gamab_pm3mais( 17, 1, 1) = 1.151939d0
    gamab_pm3mais( 17, 1, 2) = 1.176673d0
    gamab_pm3mais( 17, 1, 3) = 0.445505d0
! Cl-O
    alpab_pm3mais( 17, 8, 1) = 0.483247d0
    alpab_pm3mais( 17, 8, 2) = -0.042831d0
    alpab_pm3mais( 17, 8, 3) = -1.795162d0
    betab_pm3mais( 17, 8, 1) = 3.653024d0
    betab_pm3mais( 17, 8, 2) = 1.320588d0
    betab_pm3mais( 17, 8, 3) = 2.111719d0
    gamab_pm3mais( 17, 8, 1) = 2.404136d0
    gamab_pm3mais( 17, 8, 2) = 3.989375d0
    gamab_pm3mais( 17, 8, 3) = 1.052972d0
!FPG  end

!-----------------------------------------
!END SECTION 4 - PM3MAIS  core core function terms
!-----------------------------------------

!-----------------------------------------------
! OPNQ parameters
! Cl and F are from T. Giese and D. York, J. Chem. Phys, v127, p194101 (2007)
! Others: best fit to Parm99/10
!
! By Taisung Lee (Rutgers, 2011)
!-----------------------------------------------

        atomic_number=1  ! Parm99,  HC          1.4870  0.0157             OPLS
        qxd_supported(atomic_number)=.true.
        qxd_s(atomic_number)=5.242300D0
        qxd_z0(atomic_number)=3.078201D0
        qxd_zq(atomic_number)=0.0D0 
        qxd_d0(atomic_number)=1.745835D0
        qxd_dq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_q0(atomic_number)=0.0D0
        qxd_qq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_neff(atomic_number)= 0.824D0

        atomic_number=6  ! Parm99, C*           1.9080  0.0860 
        qxd_supported(atomic_number)=.true.
        qxd_s(atomic_number)=21.308085D0
        qxd_z0(atomic_number)=2.486688D0
        qxd_zq(atomic_number)=0.0D0 
        qxd_d0(atomic_number)=8.554845D0
        qxd_dq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_q0(atomic_number)=0.0D0
        qxd_qq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_neff(atomic_number)=2.657D0


        atomic_number=7  ! Parm99, N (OPLS)  N           1.8240  0.1700   
        qxd_supported(atomic_number)=.true.
        qxd_s(atomic_number)=24.833864D0
        qxd_z0(atomic_number)=2.535777D0
        qxd_zq(atomic_number)=0.0D0 
        qxd_d0(atomic_number)=10.759963D0
        qxd_dq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_q0(atomic_number)=0.0D0
        qxd_qq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_neff(atomic_number)=3.187D0
        
        
        atomic_number=8  ! Parm99, OW (TIP3P water) ! OW          1.7683  0.1520   
        qxd_supported(atomic_number)=.true.
        qxd_s(atomic_number)=21.730536D0
        qxd_z0(atomic_number)=2.599015D0
        qxd_zq(atomic_number)=0.0D0 
        qxd_d0(atomic_number)=8.526693D0
        qxd_dq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_q0(atomic_number)=0.0D0
        qxd_qq(atomic_number)=-3.0*0.D0  ! -3*B in the paper
        qxd_neff(atomic_number)=3.663D0
 

        atomic_number=9  ! F
        qxd_supported(atomic_number)=.true.
        qxd_s(atomic_number)=7.2D0
        qxd_z0(atomic_number)=2.75D0
        qxd_zq(atomic_number)=-0.37D0 
        qxd_d0(atomic_number)=3.759D0
        qxd_dq(atomic_number)=-3.0*0.22286D0  ! -3*B in the paper
        qxd_q0(atomic_number)=0.0D0
        qxd_qq(atomic_number)=-3.0*0.22286D0  ! -3*B in the paper
        qxd_neff(atomic_number)=4.086D0
     
        
        atomic_number=17  ! Cl
        qxd_supported(atomic_number)=.true.        
        qxd_s(atomic_number)=16.1D0
        qxd_z0(atomic_number)=2.15D0
        qxd_zq(atomic_number)=-0.26D0 
        qxd_d0(atomic_number)=14.71D0
        qxd_dq(atomic_number)=-3.0*0.30045D0  ! -3*B in the paper
        qxd_q0(atomic_number)=81.69D0
        qxd_qq(atomic_number)=-3.0*0.30045D0  ! -3*B in the paper
        qxd_neff(atomic_number)=5.551D0        
!-----------------------------------------------
!Load User-defined parameter file
!By Taisung Lee (Rutgers, 2011)
!-----------------------------------------------


    if(qm2_struct%ParameterFileExisting) userDefinedVariable=.true.
#ifdef MPI
       call mpi_bcast(userDefinedVariable,1,MPI_LOGICAL, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */

    if (userDefinedVariable) then 
    
        n=GetNumberParameterEntries(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting)
#ifdef MPI
            call mpi_bcast(n,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */

        !!  the read-in section
        !!  need to repeat for every type of hamitonian--kind of stupid...

        if (currentTheory%MNDO) then
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_mndo(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_mndo(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_mndo(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_mndo(atomic_number) ,temp)
                call SetUpUserParameter('GSS', GSS_mndo(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_mndo(atomic_number) ,temp)
                call SetUpUserParameter('GSP', GSP_mndo(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_mndo(atomic_number) ,temp)
                call SetUpUserParameter('HSP', HSP_mndo(atomic_number) ,temp)
                call SetUpUserParameter('ZS', s_orb_exp_mndo(atomic_number) ,temp)
                call SetUpUserParameter('ZP', p_orb_exp_mndo(atomic_number) ,temp)
                call SetUpUserParameter('ALP', alp_mndo(atomic_number) ,temp)
             end do
        end if



        if (currentTheory%PM3) then
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_pm3(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_pm3(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_pm3(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_pm3(atomic_number) ,temp)
                call SetUpUserParameter('GSS', GSS_pm3(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_pm3(atomic_number) ,temp)
                call SetUpUserParameter('GSP', GSP_pm3(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_pm3(atomic_number) ,temp)
                call SetUpUserParameter('HSP', HSP_pm3(atomic_number) ,temp)
                call SetUpUserParameter('ZS', s_orb_exp_pm3(atomic_number) ,temp)
                call SetUpUserParameter('ZP', p_orb_exp_pm3(atomic_number) ,temp)
                call SetUpUserParameter('ALP', alp_pm3(atomic_number) ,temp)
                call SetUpUserParameter('FN11', FN1_pm3(1,atomic_number) ,temp)
                call SetUpUserParameter('FN21', FN2_pm3(1,atomic_number) ,temp)
                call SetUpUserParameter('FN31', FN3_pm3(1,atomic_number) ,temp)
                call SetUpUserParameter('FN12', FN1_pm3(2,atomic_number) ,temp)
                call SetUpUserParameter('FN22', FN2_pm3(2,atomic_number) ,temp)
                call SetUpUserParameter('FN32', FN3_pm3(2,atomic_number) ,temp)
                call SetUpUserParameter('FN13', FN1_pm3(3,atomic_number) ,temp)
                call SetUpUserParameter('FN23', FN2_pm3(3,atomic_number) ,temp)
                call SetUpUserParameter('FN33', FN3_pm3(3,atomic_number) ,temp)
                call SetUpUserParameter('FN14', FN1_pm3(4,atomic_number) ,temp)
                call SetUpUserParameter('FN24', FN2_pm3(4,atomic_number) ,temp)
                call SetUpUserParameter('FN34', FN3_pm3(4,atomic_number) ,temp)
                tempIntegerInReal=NUM_FN_am1(atomic_number)     
                call SetUpUserParameter('NUM_FN', tempIntegerInReal ,temp)
                NUM_FN_pm3(atomic_number)=int(tempIntegerInReal+0.25)
            end do
        end if



        if (currentTheory%AM1) then 
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_am1(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_am1(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_am1(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_am1(atomic_number) ,temp)
                call SetUpUserParameter('GSS', GSS_am1(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_am1(atomic_number) ,temp)
                call SetUpUserParameter('GSP', GSP_am1(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_am1(atomic_number) ,temp)  
                call SetUpUserParameter('HSP', HSP_am1(atomic_number) ,temp)  
                call SetUpUserParameter('ZS', s_orb_exp_am1(atomic_number) ,temp)  
                call SetUpUserParameter('ZP', p_orb_exp_am1(atomic_number) ,temp) 
                call SetUpUserParameter('ALP', alp_am1(atomic_number) ,temp)
                call SetUpUserParameter('FN11', FN1_am1(1,atomic_number) ,temp)
                call SetUpUserParameter('FN21', FN2_am1(1,atomic_number) ,temp)
                call SetUpUserParameter('FN31', FN3_am1(1,atomic_number) ,temp)
                call SetUpUserParameter('FN12', FN1_am1(2,atomic_number) ,temp)
                call SetUpUserParameter('FN22', FN2_am1(2,atomic_number) ,temp)
                call SetUpUserParameter('FN32', FN3_am1(2,atomic_number) ,temp)
                call SetUpUserParameter('FN13', FN1_am1(3,atomic_number) ,temp)
                call SetUpUserParameter('FN23', FN2_am1(3,atomic_number) ,temp)
                call SetUpUserParameter('FN33', FN3_am1(3,atomic_number) ,temp)
                call SetUpUserParameter('FN14', FN1_am1(4,atomic_number) ,temp)
                call SetUpUserParameter('FN24', FN2_am1(4,atomic_number) ,temp)
                call SetUpUserParameter('FN34', FN3_am1(4,atomic_number) ,temp)
                tempIntegerInReal=NUM_FN_am1(atomic_number)     
                call SetUpUserParameter('NUM_FN', tempIntegerInReal ,temp)
                NUM_FN_am1(atomic_number)=int(tempIntegerInReal+0.25)
            end do
        end if
                
        if (currentTheory%MNDOD) then 
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_mndod(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_mndod(atomic_number) ,temp)
                call SetUpUserParameter('UDD', UDD_mndod(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_mndod(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_mndod(atomic_number) ,temp)
                call SetUpUserParameter('BETAD', betad_mndod(atomic_number) ,temp)
                call SetUpUserParameter('GSS', GSS_mndod(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_mndod(atomic_number) ,temp)
                call SetUpUserParameter('GDD', GDD_mndod(atomic_number) ,temp)
                call SetUpUserParameter('GSP', GSP_mndod(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_mndod(atomic_number) ,temp)
                call SetUpUserParameter('HSP', HSP_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ZS', s_orb_exp_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ZP', p_orb_exp_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ZD', d_orb_exp_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ALP', alp_mndod(atomic_number) ,temp)
                call SetUpUserParameter('POCORD', rho_core_mndod(atomic_number) ,temp) 
                call SetUpUserParameter('ZST', s_orb_exp_tail_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ZPT', p_orb_exp_tail_mndod(atomic_number) ,temp)
                call SetUpUserParameter('ZDT', d_orb_exp_tail_mndod(atomic_number) ,temp)
             end do
        end if


        if (currentTheory%AM1D) then 
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_am1d(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_am1d(atomic_number) ,temp)
                call SetUpUserParameter('UDD', UDD_am1d(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_am1d(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_am1d(atomic_number) ,temp)
                call SetUpUserParameter('BETAD', betad_am1d(atomic_number) ,temp)                                
                call SetUpUserParameter('GSS', GSS_am1d(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_am1d(atomic_number) ,temp)
                call SetUpUserParameter('GDD', GDD_am1d(atomic_number) ,temp)   
                call SetUpUserParameter('GSP', GSP_am1d(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_am1d(atomic_number) ,temp)  
                call SetUpUserParameter('HSP', HSP_am1d(atomic_number) ,temp)  
                call SetUpUserParameter('ZS', s_orb_exp_am1d(atomic_number) ,temp)  
                call SetUpUserParameter('ZP', p_orb_exp_am1d(atomic_number) ,temp) 
                call SetUpUserParameter('ZD', d_orb_exp_am1d(atomic_number) ,temp) 
                call SetUpUserParameter('ALP', alp_am1d(atomic_number) ,temp)
                call SetUpUserParameter('ZST', s_orb_exp_tail_am1d(atomic_number) ,temp)
                call SetUpUserParameter('ZPT', p_orb_exp_tail_am1d(atomic_number) ,temp)
                call SetUpUserParameter('ZDT', d_orb_exp_tail_am1d(atomic_number) ,temp) 
                call SetUpUserParameter('GNN', GNN_am1d(atomic_number) ,temp)
                call SetUpUserParameter('POCORD', rho_core_am1d(atomic_number) ,temp)                
                call SetUpUserParameter('FN11', FN1_am1d(1,atomic_number) ,temp)
                call SetUpUserParameter('FN21', FN2_am1d(1,atomic_number) ,temp)
                call SetUpUserParameter('FN31', FN3_am1d(1,atomic_number) ,temp)               
                call SetUpUserParameter('FN12', FN1_am1d(2,atomic_number) ,temp)
                call SetUpUserParameter('FN22', FN2_am1d(2,atomic_number) ,temp)
                call SetUpUserParameter('FN32', FN3_am1d(2,atomic_number) ,temp)
                call SetUpUserParameter('FN13', FN1_am1d(3,atomic_number) ,temp)
                call SetUpUserParameter('FN23', FN2_am1d(3,atomic_number) ,temp)
                call SetUpUserParameter('FN33', FN3_am1d(3,atomic_number) ,temp)
                call SetUpUserParameter('FN14', FN1_am1d(4,atomic_number) ,temp)
                call SetUpUserParameter('FN24', FN2_am1d(4,atomic_number) ,temp)
                call SetUpUserParameter('FN34', FN3_am1d(4,atomic_number) ,temp)
                tempIntegerInReal=NUM_FN_am1(atomic_number)     
                call SetUpUserParameter('NUM_FN', tempIntegerInReal ,temp)
                NUM_FN_am1d(atomic_number)=int(tempIntegerInReal+0.25)
            end do
        
        end if

        if (currentTheory%PM6) then
            do i=1, n
                temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
                atomic_number=temp%atomicNumber
                call SetUpUserParameter('USS', USS_pm6(atomic_number) ,temp)
                call SetUpUserParameter('UPP', UPP_pm6(atomic_number) ,temp)
                call SetUpUserParameter('UDD', UDD_pm6(atomic_number) ,temp)
                call SetUpUserParameter('BETAS', betas_pm6(atomic_number) ,temp)
                call SetUpUserParameter('BETAP', betap_pm6(atomic_number) ,temp)
                call SetUpUserParameter('BETAD', betad_pm6(atomic_number) ,temp)                    
                call SetUpUserParameter('GSS', GSS_pm6(atomic_number) ,temp)
                call SetUpUserParameter('GPP', GPP_pm6(atomic_number) ,temp)
                call SetUpUserParameter('GSP', GSP_pm6(atomic_number) ,temp)
                call SetUpUserParameter('GP2', GP2_pm6(atomic_number) ,temp)
                call SetUpUserParameter('HSP', HSP_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZS', s_orb_exp_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZP', p_orb_exp_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZD', d_orb_exp_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ALP', alp_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZST', s_orb_exp_tail_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZPT', p_orb_exp_tail_pm6(atomic_number) ,temp)
                call SetUpUserParameter('ZDT', d_orb_exp_tail_pm6(atomic_number) ,temp)
                call SetUpUserParameter('POCORD', rho_core_pm6(atomic_number) ,temp)

                call SetUpUserParameter('FN11', FN1_pm6(1,atomic_number) ,temp)
                call SetUpUserParameter('FN21', FN2_pm6(1,atomic_number) ,temp)
                call SetUpUserParameter('FN31', FN3_pm6(1,atomic_number) ,temp)
                call SetUpUserParameter('FN12', FN1_pm6(2,atomic_number) ,temp)
                call SetUpUserParameter('FN22', FN2_pm6(2,atomic_number) ,temp)
                call SetUpUserParameter('FN32', FN3_pm6(2,atomic_number) ,temp)
                call SetUpUserParameter('FN13', FN1_pm6(3,atomic_number) ,temp)
                call SetUpUserParameter('FN23', FN2_pm6(3,atomic_number) ,temp)
                call SetUpUserParameter('FN33', FN3_pm6(3,atomic_number) ,temp)
                call SetUpUserParameter('FN14', FN1_pm6(4,atomic_number) ,temp)
                call SetUpUserParameter('FN24', FN2_pm6(4,atomic_number) ,temp)
                call SetUpUserParameter('FN34', FN3_pm6(4,atomic_number) ,temp)
                tempIntegerInReal=NUM_FN_pm6(atomic_number)
                call SetUpUserParameter('NUM_FN', tempIntegerInReal ,temp)
                NUM_FN_pm6(atomic_number)=int(tempIntegerInReal+0.25)
            end do

        end if

        
        ! OPNQ (QXD) correction for vdW.  The parameters are hamitonian-independent
        do i=1, n
            temp=GetParameterEntry(qm2_struct%parameterEntries,qm2_struct%ParameterFileExisting,i)
#ifdef MPI
                call mpi_bcast(temp%name,8,MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%value,1,AMBER_MPI_REAL, 0, qmmm_mpi%commqmmm, ierr)
                call mpi_bcast(temp%atomicNumber,1,MPI_INTEGER, 0, qmmm_mpi%commqmmm, ierr)
#endif /* MPI */
            atomic_number=temp%atomicNumber
            call SetUpUserParameter('QXD_S', qxd_s(atomic_number) ,temp)
            call SetUpUserParameter('QXD_Z0', qxd_z0(atomic_number) ,temp)
            call SetUpUserParameter('QXD_ZQ', qxd_zq(atomic_number) ,temp) 
            call SetUpUserParameter('QXD_D0', qxd_d0(atomic_number) ,temp)
            call SetUpUserParameter('QXD_DQ', qxd_dq(atomic_number) ,temp)
            call SetUpUserParameter('QXD_Q0', qxd_q0(atomic_number) ,temp)
            call SetUpUserParameter('QXD_QQ', qxd_qq(atomic_number) ,temp)  
            call SetUpUserParameter('QXD_NEFF', qxd_neff(atomic_number) ,temp)
                                            
         end do
   
    end if  !(ParameterFileExisting)


end subroutine InitializeParameter


subroutine SetUpUserParameter(name, target, temp)

    use ParameterReader, only: ParameterEntry
    use UtilitiesModule, only: Upcase

    implicit none

    character(*), intent(in)::name
    _REAL_,intent(inout)::target
    type(ParameterEntry), intent(in)::temp


    if (Upcase(temp%name)==name) then
      target=temp%value
    end if 
      

end subroutine SetUpUserParameter 

end module QM2_parameters
