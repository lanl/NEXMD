! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!This file contains several routines for printing information
!about QMMM simulations.

subroutine qm2_print_info
! ----------------------------------------------------------------------
! PURPOSE: Print information on QM/MM calculation setup in sander / sqm
!
! Driver routine, uses other subroutines in from qm_print_info
! Part of this has been dissected from subroutine qm2_load_params_and_allocate
!
! Authors: Andreas W. Goetz, Ross C. Walker (SDSC)
! Date   : February 2010
! ----------------------------------------------------------------------

  use qmmm_module, only : qmmm_nml, &
                          qmmm_struct, &
                          qm2_struct
  use qmmm_qmtheorymodule, only : String
  use dh_correction_module, only : dh_correction_info
  ! Include the file containing all of the parameter constants.
  ! (required for the paper reference indices)
  ! This module should be cleaned up at some point!
  use QM2_parameters

  implicit none

  integer :: i
  integer :: reference_index

  ! -------------------------
  ! Print references (papers)
  ! -------------------------
  call qm_print_ref(.true.,1, 1, qmmm_nml%qmtheory%DFTB)

  ! -------------------------------------------------------
  ! Information on spin state and number of occupied levels
  ! -------------------------------------------------------
  write(6,'()')
  if (qmmm_nml%spin==1) then
     write (6,'(''QMMM: SINGLET STATE CALCULATION'')')
  else if (qmmm_nml%spin==2) then
     write (6,'(''QMMM: DOUBLET STATE CALCULATION'')')
  else if (qmmm_nml%spin==3) then
     write (6,'(''QMMM: TRIPLET STATE CALCULATION'')')
  else if (qmmm_nml%spin==4) then
     write (6,'(''QMMM: QUARTET STATE CALCULATION'')')
  else if (qmmm_nml%spin==5) then
     write (6,'(''QMMM: QUINTET STATE CALCULATION'')')
  else if (qmmm_nml%spin==6) then
     write (6,'(''QMMM: SEXTET STATE CALCULATION'')')
  end if

  if (qmmm_nml%spin == 1 ) then
     write(6,'(''QMMM: RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS ='',I3)') qm2_struct%nclosed
  else
     ! This will have to be updated one day - we should support UHF runs, as well
     write(6,'(''QMMM: ROHF CALCULATION'')')
     write(6,'(''QMMM: THERE ARE'',I3,'' DOUBLY FILLED LEVELS'')') qm2_struct%nclosed
     write(6,'(''QMMM: AND'')')
     write(6,'(''QMMM:          '',i3,'' SINGLY OCCUPIED LEVELS'')') qm2_struct%nopenclosed - qm2_struct%nclosed
  end if

  ! --------------------------------------------------------
  ! Print references for semiempirical parameter sets in use
  ! --------------------------------------------------------
  if (.not. qmmm_nml%qmtheory%DFTB) then
     write(6,'(/,a)') '| QMMM: *** Selected Hamiltonian *** '
     write(6,'(2a)')  '| QMMM: ', String(qmmm_nml%qmtheory)
     if (qmmm_nml%qmtheory%PM6) then
        write(6,'(a)') '| QMMM: J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007)'
        write(6,'(a)') '| QMMM: and unpublished corrections'
     end if
     call dh_correction_info(qmmm_nml%qmtheory)
     write(6,'(/,"| QMMM: *** Parameter sets in use ***")')
     do i = 1, qmmm_struct%qm_ntypes
        if (qmmm_nml%qmtheory%MNDO) then
           reference_index = mndo_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%MNDOD) then
           reference_index = mndod_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%AM1) then
           reference_index = am1_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%AM1D) then
           reference_index = am1d_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%PM3 .OR. qmmm_nml%qmtheory%PM3CARB1 &
             .OR. qmmm_nml%qmtheory%PM3ZNB .OR. qmmm_nml%qmtheory%PM3MAIS) then
           reference_index = pm3_ref_index(qmmm_struct%qm_type_id(i))
           if (qmmm_nml%qmtheory%PM3CARB1 .and. &
                (qmmm_struct%qm_type_id(i) == 1 .or. qmmm_struct%qm_type_id(i) == 8) ) then
              reference_index = pm3carb1_ref_index(qmmm_struct%qm_type_id(i))
           else if (qmmm_nml%qmtheory%PM3ZNB .and. &
                (qmmm_struct%qm_type_id(i) == 30 ) ) then
              reference_index = pm3znb_ref_index(qmmm_struct%qm_type_id(i))
           end if
        else if (qmmm_nml%qmtheory%PM6) then
           reference_index = pm6_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%PDDGPM3) then
           reference_index = pddgpm3_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%PDDGPM3_08) then
           reference_index = pddgpm3_08_ref_index(qmmm_struct%qm_type_id(i))
        elseif (qmmm_nml%qmtheory%PDDGMNDO) then
           reference_index = pddgmndo_ref_index(qmmm_struct%qm_type_id(i))
        else if (qmmm_nml%qmtheory%RM1) then
           reference_index = rm1_ref_index(qmmm_struct%qm_type_id(i))
        endif
        call qm_print_ref(.false., reference_index, qmmm_struct%qm_type_id(i), qmmm_nml%qmtheory%DFTB)
     end do
     if (qmmm_nml%qmtheory%PM3MAIS) then
        write(6,'(a)') '| QMMM: MAIS Ref: M.I. BERNAL-URUCHURTU et al. CPL 330, 118 (2000)'
        write(6,'(a)') '| QMMM:   for OH: O.I. ARILLO-FLORES et al. TCAcc 118, 425 (2007)'
     end if
  end if

  ! ---------------------------------------------
  ! Print the reference for the PM3/MM* interface
  ! ---------------------------------------------
  !if (qmmm_nml%qmmm_int == 3 .AND. qmmm_nml%qmtheory%PM3) then
  if (qmmm_struct%PM3MMX_INTERFACE) then
     write (6,'(/,"| QMMM: *** PM3/MM* (WITH MODIFIED QM-MM INTERFACE) APPLIED ***")')
     write (6,'("| QMMM: Ref: Q.T.WANG and R.A.BRYCE, JCTC, 5, 2206, (2009)")') 
  end if

  ! ---------------------------------
  ! Information on peptide correction
  ! ---------------------------------
  if (qmmm_nml%peptide_corr) then
     write (6,'(''QMMM: MOLECULAR MECHANICS CORRECTION APPLIED TO PEPTIDE LINKAGES'')')
     write (6,'(''QMMM: '',i5,'' PEPTIDE LINKAGES HAVE BEEN FOUND:'')') qm2_struct%n_peptide_links
     do i = 1, qm2_struct%n_peptide_links
        write(6,'(''QMMM:    '',i4,'' - '',i4,'' - '',i4,'' - '',i4)') qm2_struct%peptide_links(1,i), &
             qm2_struct%peptide_links(2,i), qm2_struct%peptide_links(3,i), qm2_struct%peptide_links(4,i)
     end do
  end if

#ifdef SQM
  ! ---------------------------------------
  ! Information on SCF convergence criteria
  ! ---------------------------------------
  write(6,'(/,"| QMMM: *** SCF convergence criteria ***")')
  write(6,'("| QMMM: Energy change                : ",D8.1," kcal/mol")') qmmm_nml%scfconv
  write(6,'("| QMMM: Error matrix |FP-PF|         : ",D8.1," au")')       qmmm_nml%errconv
  write(6,'("| QMMM: Density matrix change        : ",D8.1)')             qmmm_nml%density_conv
  write(6,'("| QMMM: Maximum number of SCF cycles : ",i8)')               qmmm_nml%itrmax
#endif

end subroutine qm2_print_info

!------------------------------------------------------------------------------
subroutine qm_print_coords(nstep,print_coor)
!This routine prints the qm region coordinates.
!Including link atoms
!
!Written by Ross Walker, TSRI, 2004
!
!================================================

  use ElementOrbitalIndex, only : elementSymbol
  use qmmm_module, only : qmmm_struct, qmmm_nml
  use file_io_dat, only : MAX_FN_LEN
  implicit none

!Passed in
      integer, intent(in) :: nstep
      logical, intent(in) :: print_coor

!Local
      integer i,j
      character(len=MAX_FN_LEN) :: qmmm_pdb_filename

!When doing nearest solvents this will get called every time an identify
!nearest solvents is done.

      if (print_coor .or. qmmm_nml%verbosity > 1) then
        write(6,'(/,''  QMMM: QM Region Cartesian Coordinates (*=link atom) '')')
        write(6,'(''  QMMM: QM_NO.   MM_NO.'',2X,''ATOM'',9X,''X'',9X,''Y'',9X,''Z'')')
        do I=1,qmmm_struct%nquant
           WRITE(6,'("  QMMM:",I6,2X,I7,6X,A2,3X,3F10.4)') I, &
                 qmmm_struct%iqmatoms(i), &
                 elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
                 (qmmm_struct%qm_coords(J,I),J=1,3)
        end do
        do i=qmmm_struct%nquant+1,qmmm_struct%nquant_nlink
           WRITE(6,'("  QMMM:",I6,2X,7X,5X,"*",A2,3X,3F10.4)') I, &
                 elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
                 (qmmm_struct%qm_coords(J,I),J=1,3)
        end do
      end if

      if (qmmm_nml%writepdb) then
        if (qmmm_nml%vsolv > 0) then
           write(qmmm_pdb_filename,'(I8)') nstep
           qmmm_pdb_filename = 'qmmm_region.pdb.' // adjustl(qmmm_pdb_filename)
           write(6,'(a,a)') 'QMMM: Writing QM coordinates to PDB file: ',qmmm_pdb_filename
           call qm_write_pdb(qmmm_pdb_filename)
        else
          !We write a pdb of the coordinates
          call qm_write_pdb('qmmm_region.pdb')
        end if
      end if

end subroutine qm_print_coords

!------------------------------------------------------------------------------

subroutine qm_write_pdb(filename)
! This routine will write a crude pdb of the coordinates representing
! the QM atoms and link atoms. This allows the user to visually check that
! the selected QM region is what they expected.

! Written by Ross Walker (TSRI, 2005)

  use ElementOrbitalIndex, only : elementSymbol
  use qmmm_module, only : qmmm_struct
  implicit none

  character (len=*) :: filename
  
!Local
  integer :: i,j

  call amopen(23,trim(filename),'R','F','W')
  write(23,'(a)') 'REMARK'
  do i=1,qmmm_struct%nquant_nlink
    write(23,60) i,elementSymbol(qmmm_struct%iqm_atomic_numbers(i)),1, &
                 (qmmm_struct%qm_coords(J,I),J=1,3)
  end do
  write(23,'(a)') 'END'

  close(unit=23) 

  60 FORMAT('ATOM',2X,I5,1X,A4,1X,'QMR ',1X, I4,4X,3F8.3)

end subroutine qm_write_pdb

!------------------------------------------------------------------------------

subroutine qm_print_dyn_mem(natom,npairs)
!This routine prints a summary of the dynamic
!memory allocated for use in the QM calculation.
!Written by Ross Walker, TSRI, 2004
!Assumes _REAL_ is double precision = 8 bytes.
!================================================

  use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, qm2_rij_eqns, &
                          qmewald, qm_gb, qmmm_mpi, qmmm_scratch
  implicit none

#include "qm2_array_locations.h"

!Passed in
  integer, intent(in) :: natom, npairs

!Local
  integer total_memory, element_memory
  _REAL_, parameter :: bytes_to_mb = 1.0D0/(1024D0*1024D0)
  integer bytes_per_int, bytes_per_real, bytes_per_logical
  total_memory = 0
  bytes_per_int = bit_size(element_memory)/8 
  bytes_per_real = 8     !Assume size of _REAL_ is 8 bytes
  bytes_per_logical = 1  !Assume size of logical is a single byte

  write(6,'(/"| QMMM: Estimated QM Dynamic Memory Usage (per thread)")')
  write(6,'("| QMMM: ---------------------------------------------------")')

  element_memory = size(qmmm_struct%qm_atom_type)*bytes_per_real &
                  +size(qmmm_struct%qm_type_id)*bytes_per_real
  write(6,'("| QMMM:              QM Atom Type Info : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%qm_resp_charges)*bytes_per_real
  write(6,'("| QMMM:         QM RESP Charge Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%iqmatoms)*bytes_per_int
  write(6,'("| QMMM:            QM Atom Number List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = 0
  if (qmmm_struct%nlink > 0) then  
    element_memory = size(qmmm_struct%link_pairs)*bytes_per_int
    element_memory = element_memory + &
                                 size(qmmm_struct%mm_link_pair_resp_charges) * bytes_per_real
    element_memory = element_memory + &
                                 size(qmmm_struct%mm_link_pair_saved_coords) * bytes_per_real
  end if
  write(6,'("| QMMM:                Link Atom Pairs : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if(qmmm_nml%peptide_corr .and. qm2_struct%n_peptide_links>0) then
    element_memory = size(qm2_struct%peptide_links)*bytes_per_int
    write(6,'("| QMMM:       Peptide Linkage Identity : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if

  element_memory = size(qmmm_struct%iqm_atomic_numbers)*bytes_per_int
  write(6,'("| QMMM:          QM Atomic Number List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  !Pair list is always allocated as at least 1 even if natom = nquant.
  element_memory = size(qmmm_struct%qm_mm_pair_list)*bytes_per_int
  write(6,'("| QMMM:                QM-MM Pair List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = (size(qmmm_struct%atom_mask)+size(qmmm_struct%mm_link_mask))*bytes_per_logical
  write(6,'("| QMMM:                   QM Atom Mask : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%qm_coords)*bytes_per_real
!include the qm_xcrd(4,natom) array here. This is declared
!in qm_mm.f
  element_memory = element_memory+4*natom*bytes_per_real
  write(6,'("| QMMM:           QM Coordinate Arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%scaled_mm_charges)*bytes_per_real
  write(6,'("| QMMM:         Scaled MM Charge Array : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  !SCF Mulliken Charge array
  element_memory = size(qm2_struct%scf_mchg)*bytes_per_real
  write(6,'("| QMMM:    SCF Mulliken Charge Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory  
  
  !qm ewald arrays
  if (qmmm_nml%qm_ewald>0) then
    element_memory = size(qmewald%kvec)+size(qmewald%dkvec)+size(qmewald%dmkv)
    if (.not. qmmm_nml%qm_pme) element_memory= element_memory + size(qmewald%ktable) + &
                                               size(qmewald%d_ewald_mm)
    element_memory = element_memory + size(qmewald%qmktable) + &
                     size(qmewald%mmpot) + size(qmewald%qmpot)
    element_memory = element_memory * bytes_per_real
    write(6,'("| QMMM:                QM Ewald Arrays : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if

  !QM GB arrays
  if (qmmm_nml%qmgb==2) then
    element_memory = size(qm_gb%qmqm_onefij)+size(qm_gb%gb_mmpot)+size(qm_gb%gb_qmpot)
    if (qm_gb%saltcon_on) element_memory = element_memory + size(qm_gb%qmqm_kappafij)
    element_memory = element_memory*bytes_per_real
    element_memory = element_memory + size(qm_gb%qmqm_gb_list)*bytes_per_int
    write(6,'("| QMMM:                   QM GB Arrays : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if
   
!Include the local force matrices, dxyzqm and dxyzcl here as well
  element_memory = 3*natom*bytes_per_real
  if (qmmm_nml%qmmm_int /= 0 ) then
    element_memory =element_memory+3*qmmm_struct%nquant_nlink*bytes_per_real
  end if
  write(6,'("| QMMM:                QM Force Arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%den_matrix)*bytes_per_real
  write(6,'("| QMMM:                 Density Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = (size(qm2_struct%old_den_matrix)+size(qm2_struct%old2_density))* &
                   bytes_per_real
  if (qmmm_nml%density_predict == 1) then
    element_memory = element_memory + (( size(qm2_struct%md_den_mat_guess1) + &
                                        size(qm2_struct%md_den_mat_guess2) ) * bytes_per_real)
  end if
  write(6,'("| QMMM:          Density Matrix Copies : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = qmmm_struct%nquant_nlink*16*bytes_per_real
  write(6,'("| QMMM: Fock2 Density Matrix Workspace : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%fock_matrix)*bytes_per_real
  write(6,'("| QMMM:                    Fock Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory
  if (qmmm_nml%fock_predict == 1) then
    element_memory = ( (size(qm2_struct%fock_mat_final1) + &
                                        size(qm2_struct%fock_mat_final2) + &
                                        size(qm2_struct%fock_mat_final3) + &
                                        size(qm2_struct%fock_mat_final4)) &
                                        * bytes_per_real)
    write(6,'("| QMMM:             Fock Matrix Copies : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if

  if (qmmm_nml%qmtheory%DFTB) then
    element_memory = 0
  else
    element_memory = size(qm2_struct%eigen_vectors)*bytes_per_real
  end if
  write(6,'("| QMMM:           Eigen Vector Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmqm_erep_incore) then
     element_memory = size(qm2_struct%qm_qm_e_repul)*bytes_per_real
  else
     element_memory = 0
  end if
  write(6,'("| QMMM: QM-QM Elec Repulsion Integrals : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%qm_qm_2e_repul)*bytes_per_real
  write(6,'("| QMMM:  QM 2-Elec Repulsion Integrals : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%hmatrix)*bytes_per_real
  write(6,'("| QMMM:              1-Electron Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_params%core_chg) + size(qm2_params%betasas) + &
                      size(qm2_params%orb_elec_ke) + &
                      size(qm2_params%onec2elec_params) + &
                      size(qm2_params%multip_2c_elec_params) + &
                      size(qm2_params%betasap) + size(qm2_params%betapap)
  if (.not. qmmm_nml%qmtheory%DFTB) then
    element_memory = element_memory + size(qm2_params%atom_orb_zz_sxs_over_sas) + &
                      size(qm2_params%atom_orb_zz_sxp_over_sap) + &
                      size(qm2_params%atom_orb_zz_pxp_over_pap) + &
                      size(qm2_params%atom_orb_ss_eqn) + &
                      size(qm2_params%atom_orb_sp_ovlp) + &
                      size(qm2_params%atom_orb_pp_ovlp_inj) + &
                      size(qm2_params%atom_orb_pp_ovlp_ieqj1) + &
                      size(qm2_params%atom_orb_pp_ovlp_ieqj2) + &
                      size(qm2_params%atom_orb_ss_eqn_adb) + &
                      size(qm2_params%atom_orb_sp_eqn_xy)+ &
                      size(qm2_params%atom_orb_sp_eqn_xx1) + &
                      size(qm2_params%atom_orb_sp_eqn_xx2) + &
                      size(qm2_params%atom_orb_pp_eqn_xxy1) + &
                      size(qm2_params%atom_orb_pp_eqn_xxy2)
  end if

  if (qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PDDGMNDO .OR. qmmm_nml%qmtheory%PDDGPM3_08) then
    element_memory = element_memory + size(qm2_params%pddge1) + size(qm2_params%pddge2)+ &
                     size(qm2_params%pddg_term1) +size(qm2_params%pddg_term2) + &
                     size(qm2_params%pddg_term3)+ size(qm2_params%pddg_term4)
  end if
  if (qmmm_nml%qmtheory%AM1 .OR. qmmm_nml%qmtheory%PM3 .OR. qmmm_nml%qmtheory%PDDGPM3 .OR. &
      qmmm_nml%qmtheory%PM3CARB1 .OR. qmmm_nml%qmtheory%RM1 .OR. qmmm_nml%qmtheory%PDDGPM3_08 .OR. &
      qmmm_nml%qmtheory%PM6 .OR. qmmm_nml%qmtheory%PM3ZNB) then
    element_memory = element_memory + size(qm2_params%FN1) + size(qm2_params%FN2) + &
                     size(qm2_params%FN3)
  end if
  if (qmmm_nml%qmtheory%PM6) then
    element_memory = element_memory + size(qm2_params%pm6_alpab) + size(qm2_params%pm6_xab)
  else
    element_memory = element_memory + size(qm2_params%cc_exp_params)
  end if
  element_memory = element_memory * bytes_per_real
  write(6,'("| QMMM:       _REAL_ parameter storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_params%natomic_orbs)+size(qm2_params%orb_loc)+ &
                   size(qm2_params%pascal_tri1)+size(qm2_params%pascal_tri2) 
  element_memory = element_memory*bytes_per_int
  write(6,'("| QMMM:      integer parameter storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmmmrij_incore) then
    element_memory = min(npairs*(qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)+npairs, &
                                 (qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)*(natom-qmmm_struct%nquant_nlink))
    element_memory = QMMMNORIJ*element_memory*bytes_per_real
  else
    element_memory = 0
  end if
  write(6,'("| QMMM:         QM-MM RIJ Eqns storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory
  
  if (qmmm_nml%qmtheory%DFTB) then
    element_memory = 0
  else 
    !qmmm_scratch%mat_diag_workspace
    element_memory=size(qmmm_scratch%mat_diag_workspace)*bytes_per_real
  end if

  element_memory = element_memory+size(qmmm_scratch%qm_real_scratch)*bytes_per_real
#ifdef MPI
# ifndef USE_MPI_IN_PLACE
  !Add in the reduction matsize temp array - only allocated if we can't do MPI_IN_PLACE
  element_memory = element_memory + size(qmmm_scratch%matsize_red_scratch)*bytes_per_real
# endif
#endif
  !Add in any scratch space for divide and conquer lapack diagonaliser
  element_memory = element_memory + qmmm_scratch%lapack_dc_real_scr_aloc*bytes_per_real
  !Add in scratch space for the vectmp by npairs arrays - these are always
  !allocated even if there are no mm atoms in which case they have size 1.
  if (qmmm_nml%allow_pseudo_diag) then
    element_memory = element_memory + size(qmmm_scratch%pdiag_scr_norbs_norbs)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_scr_noccupied_norbs)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp1)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp2)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp3)*bytes_per_real
  end if
  write(6,'("| QMMM:          _REAL_ Scratch arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_scratch%qm_int_scratch)*bytes_per_int
#ifdef MPI
  !Add in the jrange arrays to the scratch count.
  element_memory = element_memory + size(qmmm_mpi%nquant_nlink_jrange)*bytes_per_int
#endif
  !Add in any scratch space for divide and conquer lapack diagonaliser
  element_memory = element_memory + qmmm_scratch%lapack_dc_int_scr_aloc*bytes_per_int
  if (qmmm_nml%allow_pseudo_diag) &
    element_memory = element_memory + size(qmmm_scratch%pdiag_vecjs)*bytes_per_int
  write(6,'("| QMMM:         Integer Scratch arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  write(6,'("| QMMM: ---------------------------------------------------")')
  write(6,'("| QMMM:        Total Dynamic Memory Usage: ",f10.3," Mb")') &
            total_memory * bytes_to_mb

end subroutine qm_print_dyn_mem

!------------------------------------------------------------------------------

subroutine qm_print_ref(amber_papers, ref_index, atomic_number, isDFTB)
!This routine prints the reference corresponding to ref_index
!The idea here is that in the parameter list a paper index is
!given to each parameter set and the relevant reference can
!be printed at startup for each element in the calculation.
!
!Written by Ross Walker, TSRI, 2005
!
!================================================

  use ElementOrbitalIndex, only : elementSymbol
  implicit none

!Passed in
      logical, intent(in) :: amber_papers
      integer, intent(in) :: ref_index, atomic_number
      logical, intent(in) :: isDFTB

      if (amber_papers) then
        !First we print info about the QM/MM implementation.
         write(6,'()')
         write(6,'("| QMMM: Citation for AMBER QMMM Run:")')
         write(6,'("| QMMM: R.C. Walker, M.F. Crowley and D.A. Case, J. COMP. CHEM. 29:1019, 2008")')

         !Next if we are doing DFTB also print the DFTB paper.
         if (isDFTB) then
           write(6,'()')
           write(6,'("| QMMM: DFTB Calculation - Additional citation for AMBER DFTB QMMM Run:")')
           write(6,'("| QMMM:   Seabra, G.M., Walker, R.C. et al., J. PHYS. CHEM. A., 111, 5655, (2007)")')
           write(6,'()')
         endif
      else
        !Reference 1 - Most MNDO parameters
        !M.J.S. DEWAR, W. THIEL, JACS., 99, 4899, (1977)
        if (ref_index == 1) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 99, 4899, (1977)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 2) then
          write(6,'("| QMMM: ",A2,": TAKEN FROM MNDOC BY W.THIEL, QCPE 438, 2, p63, (1982)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 3) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 100, 777, (1978)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 4) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 99, 5231, (1977)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 5) then
          write(6,'("| QMMM: ",A2,": L.P.DAVIS et al. JCC, 2, 433, (1981)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 6) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS 5, 375 (1986)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 7) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 100, 3607, (1978)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 8) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 7, 140, (1986)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 9) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 4, 158, (1983)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 10) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 5, 1494, (1986)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 11) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 6, 186, (1987)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 12) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 4, 542, (1983)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 13) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS,106, 6771, (1984)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 14) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 5,358, (1984)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 15) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 4, 1964, (1985)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 16) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR, et al. ORGANOMETALLICS, 4, 1973, (1985)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 17) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 107, 3902, (1985)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 18) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. THEOCHEM, 180, 1, (1988)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 19) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 9, 508, (1990)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 20) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 6, 1486, (1987)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 21) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. THEOCHEM, 187,1, (1989)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 22) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. INORG. CHEM., 29, 3881, (1990)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 23) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 7, 522, (1988)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 24) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 8, 1544, (1989)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 25) then
          write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 8, 1547, (1989)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 26) then
          write(6,'("| QMMM: ",A2,": J.J.P.STEWART, JCC, 10, 209 (1989)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 27) then
          write(6,'("| QMMM: ",A2,": J.J.P.STEWART, JCC, 12, 320 (1991)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 28) then
          write(6,'("| QMMM: ",A2,": J.P.MCNAMARA et al. CHEM. PHYS. LETT., 394, 429, (2004)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 29) then
          write(6,'("| QMMM: ",A2,": REPASKY et al. JCC, 23, 1601, (2002)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 30) then
          write(6,'("| QMMM: ",A2,": TUBERT-BROHMAN et al. JCC, 25, 138, (2003)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 31) then
          write(6,'("| QMMM: ",A2,": TUBERT-BROHMAN et al. JCTC, 1, 817, (2005)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 32) then
          write(6,'("| QMMM: ",A2,": J.P.MCNAMARA et al. J. MOL. GRA. MOD., 24, 128, (2005)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 33) then
          write(6,'("| QMMM: ",A2,": G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 34) then
          write(6,'("| QMMM: ",A2,": J. Tirado-Rives et al. J. CHEM. THEO. COMP., 4, 297, (2008)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 35) then
          write(6,'("| QMMM: ",A2,": J.J.P. Stewart, J. Mol. Mod., 13, 1173, (2007)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 36) then
          write(6,'("| QMMM: ",A2,": E.N. Brothers, D. Suarez, D.W. Deerfield II, K. Merz, J. Comp. Chem., 25, 1677, (2004)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 101) then  ! MNDO/d paper
          write(6,'("| QMMM: ",A2,": W. Thiel, A. Voityuk, J. Phys. Chem., 100, 616, (1996)")') &
                        elementSymbol(atomic_number)
        else if (ref_index == 102) then  ! AM1/d-PhoT
          write(6,'("| QMMM: ",A2,": K. Nam, Q. Cui, J. Gao, D. York. J. CHEM. THEO. COMP., 3, 486, (2007)")') &
                        elementSymbol(atomic_number) 
        else if (ref_index == 103) then  ! Magnesum AM1
          write(6,'("| QMMM: ",A2,": M.C. HUTTER, J.R. REIMERS, N.S. HUSH, J.PHYS.CHEM. (1998)")') &
                        elementSymbol(atomic_number) 
        else if (ref_index == 104) then  ! Magnesum AM1d
          write(6,'("| QMMM: ",A2,": P. Imhof et al, J. Chem.  Theo.  Comp., 2, 1050-1056 (2006)")') &
                        elementSymbol(atomic_number)                        
        else
          !UNKNOWN REFERENCE
          write(6,'(''| QMMM: '',A2,'': REFERENCE UNKNOWN.'')') elementSymbol(atomic_number)
        endif
      endif

end subroutine qm_print_ref

!------------------------------------------------------------------------------
#ifdef OPENMP
subroutine qm_print_omp_info()

   use qmmm_module, only : qmmm_nml

   implicit none

   integer :: omp_get_num_threads !omp_get_num_threads is a function

!if OpenMP is going to be in use (to do the diagonalization) then the master thread
!when it gets to doing the diagonalisation will spawn shared memory threads to match
!the number of threads specified with the environment variable OMP_NUM_THREADS.
!Here we will print an informational message.
   call omp_set_num_threads(qmmm_nml%qmmm_omp_max_threads)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER

!Set the number of openmp threads for QMMM to the value specified by qmmm_nml%qmmm_omp_max_threads
!and check that it sticks.
   write(6,'(/a)') '| QMMM(OpenMP): SMP Multithreading is in use for matrix diagonalization routines.'
   write(6,'(a,i4,a/)') '| QMMM(OpenMP): ',omp_get_num_threads(),' threads will be spawned by the master MPI thread as needed.'
!return it to 1.
!$OMP END MASTER
!$OMP END PARALLEL
   call omp_set_num_threads(1)

   return

end subroutine qm_print_omp_info
#endif
