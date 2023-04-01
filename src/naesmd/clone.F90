#include "assert.fh"
#include "dprec.fh"
module clone_module
    use communism
    use qm2_davidson_module, only      : qm2_davidson_structure_type
    use qmmm_struct_module, only       : qmmm_struct_type
    use qmmm_module, only : qm2_structure, qm_ewald_structure, &
         qm2_rij_eqns_structure, qm_gb_structure, qmmm_mpi_structure, &
         qmmm_opnq_structure, qmmm_scratch_structure, qmmm_div_structure, &
         qmmm_openmp_structure 
    use qmmm_vsolv_module , only : qmmm_vsolv_type
    use naesmd_module
    use md_module, only            : md_structure
    use xlbomd_module, only            : xlbomd_structure
    use naesmd_constants
    use cosmo_C, only : cosmo_C_structure          
    use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type
    use AIMC_type_module, only : AIMC_type
implicit none
private
public :: clone_sim

interface clone
        module procedure clone_real_0
        module procedure clone_real_1
        module procedure clone_real_2
        module procedure clone_complex_0
        module procedure clone_complex_1
        module procedure clone_complex_2
        module procedure clone_int_0
        module procedure clone_int_1
        module procedure clone_int_2
        module procedure clone_character_0
        module procedure clone_logical
        module procedure clone_real_pointer_1
        module procedure clone_real_pointer_2
        module procedure clone_real_pointer_3
        module procedure clone_real_pointer_4
        module procedure clone_int_pointer_1
        module procedure clone_int_pointer_2
        module procedure clone_logical_pointer_1
        module procedure clone_logical_pointer_2
end interface

contains

subroutine clone_real_0(old,new)
_REAL_, intent(in) :: old
_REAL_, intent(out) :: new
    new=old
end subroutine

subroutine clone_real_1(old,new)
_REAL_, intent(in),allocatable :: old(:)
_REAL_, intent(out), allocatable :: new(:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_real_2(old,new)
_REAL_, intent(in),allocatable :: old(:,:)
_REAL_, intent(out), allocatable :: new(:,:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine

subroutine clone_real_pointer_1(old,new)
_REAL_, dimension(:), intent(in),pointer :: old
_REAL_, dimension(:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_real_pointer_2(old,new)
_REAL_, dimension(:,:), intent(in),pointer :: old
_REAL_, dimension(:,:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine

subroutine clone_real_pointer_3(old,new)
_REAL_, dimension(:,:,:), intent(in),pointer :: old
_REAL_, dimension(:,:,:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:,:,:)=old(:,:,:)
    endif
end subroutine

subroutine clone_real_pointer_4(old,new)
_REAL_, dimension(:,:,:,:), intent(in),pointer :: old
_REAL_, dimension(:,:,:,:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:,:,:,:)=old(:,:,:,:)
    endif
end subroutine

subroutine clone_int_pointer_1(old,new)
integer, dimension(:), intent(in),pointer :: old
integer, dimension(:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_int_pointer_2(old,new)
integer, dimension(:,:), intent(in),pointer :: old
integer, dimension(:,:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine

subroutine clone_logical_pointer_1(old,new)
logical, dimension(:), intent(in),pointer :: old
logical, dimension(:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_logical_pointer_2(old,new)
logical, dimension(:,:), intent(in),pointer :: old
logical, dimension(:,:), intent(out), pointer :: new
    if(associated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine

subroutine clone_complex_0(old,new)
complex, intent(in) :: old
complex, intent(out) :: new
    new=old
end subroutine

subroutine clone_complex_1(old,new)
complex, intent(in),allocatable :: old(:)
complex, intent(out), allocatable :: new(:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_complex_2(old,new)
complex, intent(in),allocatable :: old(:,:)
complex, intent(out), allocatable :: new(:,:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine
subroutine clone_int_0(old,new)
integer, intent(in) :: old
integer, intent(out) :: new
    new=old
end subroutine

subroutine clone_int_1(old,new)
integer, intent(in), allocatable :: old(:)
integer, intent(out), allocatable :: new(:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:)=old(:)
    endif
end subroutine

subroutine clone_int_2(old,new)
integer, intent(in), allocatable :: old(:,:)
integer, intent(out), allocatable :: new(:,:)
    if(allocated(old)) then
        allocate(new,SOURCE=old)
        new(:,:)=old(:,:)
    endif
end subroutine

subroutine clone_character_0(old,new)
character(*), intent(in) :: old
character(*), intent(out) :: new
    new=trim(old)
end subroutine

!subroutine clone_character_1(old,new)
!character(:), intent(in), allocatable :: old
!character(:), intent(out), allocatable :: new(:)
!    if(allocated(old)) then
!        allocate(new,SOURCE=old)
!        new(:)=trim(old(:))
!    endif
!end subroutine

subroutine clone_logical(old,new)
logical, intent(in) :: old
logical, intent(out) :: new
    new=old
end subroutine

subroutine clone_sim(sim1, sim2)
!sim 2 is new sim 1 is old
    type(simulation_t),pointer::sim1
    type(simulation_t),pointer::sim2
    
    integer :: i, j

        print *, "cloning sim"     
    sim2%Nsim=sim1%Nsim
    sim2%Na=sim1%Na
    sim2%nbasis=sim1%nbasis
    sim2%Ncharge=sim1%Ncharge
    sim2%ibo=sim1%ibo
    sim2%excN=sim1%excN
    sim2%dotrivial=sim1%dotrivial
    sim2%cohertype=sim1%cohertype
    sim2%lprint=sim1%lprint
    sim2%itime1=sim1%itime1
    sim2%outfile_1=sim1%outfile_1  
    sim2%outfile_2=sim1%outfile_2  
    sim2%outfile_3=sim1%outfile_3
    sim2%outfile_4=sim1%outfile_4
    sim2%outfile_5=sim1%outfile_5
    sim2%outfile_6=sim1%outfile_6
    sim2%outfile_7=sim1%outfile_7
    sim2%outfile_8=sim1%outfile_8
    sim2%outfile_9=sim1%outfile_9
    sim2%outfile_10=sim1%outfile_10
    sim2%outfile_11=sim1%outfile_11
    sim2%outfile_13=sim1%outfile_13
    sim2%outfile_14=sim1%outfile_14
    sim2%outfile_15=sim1%outfile_15
    sim2%outfile_16=sim1%outfile_16
    sim2%outfile_17=sim1%outfile_17
    sim2%outfile_18=sim1%outfile_18
    sim2%outfile_19=sim1%outfile_19
    sim2%outfile_20=sim1%outfile_20
    sim2%outfile_21=sim1%outfile_21
    sim2%outfile_22=sim1%outfile_22
    sim2%outfile_25=sim1%outfile_25
    sim2%outfile_27=sim1%outfile_27
    sim2%outfile_28=sim1%outfile_28
    sim2%outfile_29=sim1%outfile_29
    sim2%constcoherE0=sim1%constcoherE0
    sim2%constcoherC=sim1%constcoherC
    sim2%time_sqm_took=sim1%time_sqm_took
    sim2%time_davidson_took=sim1%time_davidson_took
    sim2%time_deriv_took=sim1%time_deriv_took
    sim2%time_nact_took=sim1%time_nact_took
    !Allocate then copy 
    if(.not.sim2%forces_allocated) then
        allocate(sim2%deriv_forces(3*sim2%Na))
        allocate(sim2%deriv_forces_state(sim2%excN,3*sim2%Na))
        allocate(sim2%coords(sim2%Na*3))
        allocate(sim2%input_line(size(sim1%input_line)))
    endif
    sim2%coords(:)=sim1%coords(:)
    sim2%deriv_forces(:)=sim1%deriv_forces(:)
    sim2%deriv_forces_state(:,:)=sim1%deriv_forces_state(:,:)
    sim2%input_line(:)=sim1%input_line(:)

    print *, "cloning naesmd data"
    sim2%naesmd=sim1%naesmd 
    print *, "cloning md data"     
    sim2%md=sim1%md 
    print *, "cloning qmmm data"     
    sim2%qmmm=sim1%qmmm
        !Pointers need to be deep copied to make independent versions (should have used allocatable for these instead)
        call clone(sim1%qmmm%qm_resp_charges, sim2%qmmm%qm_resp_charges)
        call clone(sim1%qmmm%qm_resp_charge_sum , sim2%qmmm%qm_resp_charge_sum )
        call clone(sim1%qmmm%mm_link_pair_resp_charges , sim2%qmmm%mm_link_pair_resp_charges)
        call clone(sim1%qmmm%mm_link_pair_saved_coords , sim2%qmmm%mm_link_pair_saved_coords)
        call clone(sim1%qmmm%qm_coords , sim2%qmmm%qm_coords )
        call clone(sim1%qmmm%scaled_mm_charges , sim2%qmmm%scaled_mm_charges )
        call clone(sim1%qmmm%dxyzqm, sim2%qmmm%dxyzqm )
        call clone(sim1%qmmm%dxyzcl , sim2%qmmm%dxyzcl )
        call clone(sim1%qmmm%qm_xcrd , sim2%qmmm%qm_xcrd )
        call clone(sim1%qmmm%switched_mmpot , sim2%qmmm%switched_mmpot )
        call clone(sim1%qmmm%qm_atom_type, sim2%qmmm%qm_atom_type )
        call clone(sim1%qmmm%link_pairs , sim2%qmmm%link_pairs )
        call clone(sim1%qmmm%iqmatoms , sim2%qmmm%iqmatoms )
        call clone(sim1%qmmm%iqm_atomic_numbers , sim2%qmmm%iqm_atomic_numbers )
        call clone(sim1%qmmm%qm_mm_pair_list , sim2%qmmm%qm_mm_pair_list )
        call clone(sim1%qmmm%qm_mm_pair_atom_numbers, sim2%qmmm%qm_mm_pair_atom_numbers )
        call clone(sim1%qmmm%atom_mask , sim2%qmmm%atom_mask )
        call clone(sim1%qmmm%mm_link_mask , sim2%qmmm%mm_link_mask )
    print *, "cloning qm2 data"    
    sim2%qm2=sim1%qm2
        call clone(sim1%qm2%den_matrix , sim2%qm2%den_matrix )
        call clone(sim1%qm2%old_den_matrix , sim2%qm2%old_den_matrix )
        call clone(sim1%qm2%old2_density , sim2%qm2%old2_density )
        call clone(sim1%qm2%md_den_mat_guess1 , sim2%qm2%md_den_mat_guess1 )
        call clone(sim1%qm2%md_den_mat_guess2 , sim2%qm2%md_den_mat_guess2 )
        call clone(sim1%qm2%fock_mat_final4 , sim2%qm2%fock_mat_final4 )
        call clone(sim1%qm2%fock_mat_final3 , sim2%qm2%fock_mat_final3 )
        call clone(sim1%qm2%fock_mat_final2 , sim2%qm2%fock_mat_final2 )
        call clone(sim1%qm2%fock_mat_final1 , sim2%qm2%fock_mat_final1 )
        call clone(sim1%qm2%fock_matrix , sim2%qm2%fock_matrix )
        call clone(sim1%qm2%fock_matrix_dp , sim2%qm2%fock_matrix_dp )
        !call clone(sim1%qm2%fock_matrix_dm , sim2%qm2%fock_matrix_dm )
        call clone(sim1%qm2%qm_mm_e_repul , sim2%qm2%qm_mm_e_repul )
        call clone(sim1%qm2%qm_qm_2e_repul , sim2%qm2%qm_qm_2e_repul )
        call clone(sim1%qm2%hmatrix , sim2%qm2%hmatrix )
        call clone(sim1%qm2%qm_qm_e_repul , sim2%qm2%qm_qm_e_repul )
        call clone(sim1%qm2%fock2_ptot2 , sim2%qm2%fock2_ptot2 )
        call clone(sim1%qm2%eigen_vectors , sim2%qm2%eigen_vectors )
        call clone(sim1%qm2%scf_mchg , sim2%qm2%scf_mchg )
        call clone(sim1%qm2%diis_fock , sim2%qm2%diis_fock )
        call clone(sim1%qm2%diis_errmat , sim2%qm2%diis_errmat )
        call clone(sim1%qm2%diis_mat , sim2%qm2%diis_mat )
        call clone(sim1%qm2%peptide_links , sim2%qm2%peptide_links )
        print *, "cloning xlbomd data"     
   sim2%xlbomd=sim1%xlbomd
        call clone(sim1%xlbomd%coef , sim2%xlbomd%coef )
        call clone(sim1%xlbomd%phi , sim2%xlbomd%phi )
        if(associated(sim1%xlbomd%phi)) then
        do i=1, sim2%xlbomd%K+2
        do j=1, sim2%xlbomd%K+1
                if(associated(sim1%xlbomd%phi_point(i)%guess,sim1%xlbomd%phi(1:sim1%xlbomd%mat_size,j))) then
                sim2%xlbomd%phi_point(i)%guess=>sim2%xlbomd%phi(1:sim1%xlbomd%mat_size,j)
                exit
                endif
        enddo
        enddo
        endif
        print *, "cloning cosmo data"
    sim2%cosmo=sim1%cosmo
        call clone(sim1%cosmo%v_solvent_difdens, sim2%cosmo%v_solvent_difdens )
        if(associated(sim1%cosmo%v_solvent_difdens)) then
        if(associated(sim1%cosmo%v_solvent_xi, sim1%cosmo%v_solvent_difdens)) then 
                sim2%cosmo%v_solvent_xi=>sim2%cosmo%v_solvent_difdens       
        endif 
        endif
        print *, "cloning qparams data"     
!!    call clone_qparams(sim1%qparams, sim2%qparams)
        sim2%qparams=sim1%qparams
!        call clone(sim1%qparams%sp_quantum_number , sim2%qparams%sp_quantum_number )
!        call clone(sim1%qparams%d_quantum_number , sim2%qparams%d_quantum_number )
!        call clone(sim1%qparams%gss , sim2%qparams%gss )
!        call clone(sim1%qparams%hsp , sim2%qparams%hsp )
!        call clone(sim1%qparams%hpp , sim2%qparams%hpp )
!        call clone(sim1%qparams%dd , sim2%qparams%dd )
!        call clone(sim1%qparams%po , sim2%qparams%po )
!        call clone(sim1%qparams%core_chg , sim2%qparams%core_chg )
!        call clone(sim1%qparams%orb_elec_ke , sim2%qparams%orb_elec_ke )
!        call clone(sim1%qparams%betasas , sim2%qparams%betasas )
!        call clone(sim1%qparams%betasap , sim2%qparams%betasap )
!        call clone(sim1%qparams%betasad , sim2%qparams%betasad )
!        call clone(sim1%qparams%betapap , sim2%qparams%betapap )
!        call clone(sim1%qparams%betapad , sim2%qparams%betapad )
!        call clone(sim1%qparams%betadad , sim2%qparams%betadad )
!        call clone(sim1%qparams%GNN , sim2%qparams%GNN )
!        call clone(sim1%qparams%rho_core , sim2%qparams%rho_core )
!        call clone(sim1%qparams%F0SD , sim2%qparams%F0SD )
!        call clone(sim1%qparams%G2SD , sim2%qparams%G2SD )
!        call clone(sim1%qparams%FN1 , sim2%qparams%FN1 )
!        call clone(sim1%qparams%FN2 , sim2%qparams%FN2 )
!        call clone(sim1%qparams%FN3 , sim2%qparams%FN3 )
!        call clone(sim1%qparams%onec2elec_params , sim2%qparams%onec2elec_params )
!        call clone(sim1%qparams%multip_2c_elec_params , sim2%qparams%multip_2c_elec_params )
!        call clone(sim1%qparams%cc_exp_params , sim2%qparams%cc_exp_params )
!        call clone(sim1%qparams%pm6_alpab , sim2%qparams%pm6_alpab )
!        call clone(sim1%qparams%pm6_xab , sim2%qparams%pm6_xab )
!        call clone(sim1%qparams%pm3mais_alpab , sim2%qparams%pm3mais_alpab )
!        call clone(sim1%qparams%pm3mais_betab , sim2%qparams%pm3mais_betab )
!        call clone(sim1%qparams%pm3mais_gamab , sim2%qparams%pm3mais_gamab )
!        call clone(sim1%qparams%s_orb_exp_by_type , sim2%qparams%s_orb_exp_by_type )
!        call clone(sim1%qparams%p_orb_exp_by_type , sim2%qparams%p_orb_exp_by_type )
!        call clone(sim1%qparams%d_orb_exp_by_type , sim2%qparams%d_orb_exp_by_type )
!        call clone(sim1%qparams%s_orb_exp_tail_by_type , sim2%qparams%s_orb_exp_tail_by_type )
!        call clone(sim1%qparams%p_orb_exp_tail_by_type , sim2%qparams%p_orb_exp_tail_by_type )
!        call clone(sim1%qparams%d_orb_exp_tail_by_type , sim2%qparams%d_orb_exp_tail_by_type )
!        call clone(sim1%qparams%pddge1 , sim2%qparams%pddge1 )
!        call clone(sim1%qparams%pddge2 , sim2%qparams%pddge2 )
!        call clone(sim1%qparams%scale_factor1_pm3mmx , sim2%qparams%scale_factor1_pm3mmx )
!        call clone(sim1%qparams%scale_factor2_pm3mmx , sim2%qparams%scale_factor2_pm3mmx )
!        call clone(sim1%qparams%rho_pm3mmx , sim2%qparams%rho_pm3mmx )
!        call clone(sim1%qparams%atom_orb_zz_sxs_over_sas , sim2%qparams%atom_orb_zz_sxs_over_sas )
!        call clone(sim1%qparams%atom_orb_zz_sxp_over_sap , sim2%qparams%atom_orb_zz_sxp_over_sap )
!        call clone(sim1%qparams%atom_orb_zz_sxd_over_sad , sim2%qparams%atom_orb_zz_sxd_over_sad )
!        call clone(sim1%qparams%atom_orb_zz_pxp_over_pap , sim2%qparams%atom_orb_zz_pxp_over_pap )
!        call clone(sim1%qparams%atom_orb_zz_pxd_over_pad , sim2%qparams%atom_orb_zz_pxd_over_pad )
!        call clone(sim1%qparams%atom_orb_zz_dxd_over_dad , sim2%qparams%atom_orb_zz_dxd_over_dad )
!        call clone(sim1%qparams%atom_orb_ss_eqn , sim2%qparams%atom_orb_ss_eqn )
!        call clone(sim1%qparams%atom_orb_sp_ovlp , sim2%qparams%atom_orb_sp_ovlp )
!        call clone(sim1%qparams%atom_orb_sd_ovlp , sim2%qparams%atom_orb_sd_ovlp )
!        call clone(sim1%qparams%atom_orb_pd_ovlp , sim2%qparams%atom_orb_pd_ovlp )
!        call clone(sim1%qparams%atom_orb_pp_ovlp_inj , sim2%qparams%atom_orb_pp_ovlp_inj )
!        call clone(sim1%qparams%atom_orb_pp_ovlp_ieqj1 , sim2%qparams%atom_orb_pp_ovlp_ieqj1 )
!        call clone(sim1%qparams%atom_orb_pp_ovlp_ieqj2 , sim2%qparams%atom_orb_pp_ovlp_ieqj2 )
!        call clone(sim1%qparams%atom_orb_dd_ovlp_inj , sim2%qparams%atom_orb_dd_ovlp_inj )
!        call clone(sim1%qparams%atom_orb_dd_ovlp_ieqj1 , sim2%qparams%atom_orb_dd_ovlp_ieqj1 )
!        call clone(sim1%qparams%atom_orb_dd_ovlp_ieqj2 , sim2%qparams%atom_orb_dd_ovlp_ieqj2 )
!        call clone(sim1%qparams%atom_orb_ss_eqn_adb , sim2%qparams%atom_orb_ss_eqn_adb )
!        call clone(sim1%qparams%atom_orb_sp_eqn_xy , sim2%qparams%atom_orb_sp_eqn_xy )
!        call clone(sim1%qparams%atom_orb_sp_eqn_xx1 , sim2%qparams%atom_orb_sp_eqn_xx1 )
!        call clone(sim1%qparams%atom_orb_sp_eqn_xx2 , sim2%qparams%atom_orb_sp_eqn_xx2 )
!        call clone(sim1%qparams%atom_orb_pp_eqn_xxy1 , sim2%qparams%atom_orb_pp_eqn_xxy1 )
!        call clone(sim1%qparams%atom_orb_pp_eqn_xxy2 , sim2%qparams%atom_orb_pp_eqn_xxy2 )
!        call clone(sim1%qparams%pddg_term1 , sim2%qparams%pddg_term1 )
!        call clone(sim1%qparams%pddg_term2 , sim2%qparams%pddg_term2 )
!        call clone(sim1%qparams%pddg_term3 , sim2%qparams%pddg_term3 )
!        call clone(sim1%qparams%pddg_term4 , sim2%qparams%pddg_term4 )
!        call clone(sim1%qparams%natomic_orbs , sim2%qparams%natomic_orbs )
!        call clone(sim1%qparams%orb_loc , sim2%qparams%orb_loc )
!        call clone(sim1%qparams%pascal_tri1 , sim2%qparams%pascal_tri1 )
!        call clone(sim1%qparams%pascal_tri2 , sim2%qparams%pascal_tri2 )
!        call clone(sim1%qparams%NUM_FN , sim2%qparams%NUM_FN )
!        call clone(sim1%qparams%qxd_supported , sim2%qparams%qxd_supported )
!        call clone(sim1%qparams%qxd_s , sim2%qparams%qxd_s )
!        call clone(sim1%qparams%qxd_z0 , sim2%qparams%qxd_z0 )
!        call clone(sim1%qparams%qxd_zq , sim2%qparams%qxd_zq )
!        call clone(sim1%qparams%qxd_d0 , sim2%qparams%qxd_d0 )
!        call clone(sim1%qparams%qxd_dq , sim2%qparams%qxd_dq )
!        call clone(sim1%qparams%qxd_q0, sim2%qparams%qxd_q0 )
!        call clone(sim1%qparams%qxd_qq , sim2%qparams%qxd_qq )
!        call clone(sim1%qparams%qxd_neff , sim2%qparams%qxd_neff )
    print *, "cloning qnml data"     
    sim2%qnml= sim1%qnml
        call clone(sim1%qnml%iqmatoms , sim2%qnml%iqmatoms)
    print *, "cloning rij data"     
    sim2%rij= sim1%rij
        call clone(sim1%rij%qmmmrijdata , sim2%rij%qmmmrijdata)
    print *, "cloning gb data"     
    sim2%gb= sim1%gb
        call clone(sim1%gb%qmqm_onefij , sim2%gb%qmqm_onefij)
        call clone(sim1%gb%qmqm_kappafij , sim2%gb%qmqm_kappafij)
        call clone(sim1%gb%gb_mmpot , sim2%gb%gb_mmpot)
        call clone(sim1%gb%gb_qmpot , sim2%gb%gb_qmpot)
        call clone(sim1%gb%qmqm_gb_list , sim2%gb%qmqm_gb_list)
    print *, "cloning qmpi data"     
    if(associated(sim1%qmpi)) then
    sim2%qmpi= sim1%qmpi
        call clone(sim1%qmpi%nquant_nlink_jrange , sim2%qmpi%nquant_nlink_jrange)
    endif
    print *, "cloning opnq data"     
        call clone(sim1%opnq%MM_atomType , sim2%opnq%MM_atomType )
        call clone(sim1%opnq%supported , sim2%opnq%supported )
        call clone(sim1%opnq%atomic_number , sim2%opnq%atomic_number )
        call clone(sim1%opnq%LJ_r , sim2%opnq%LJ_r )
        call clone(sim1%opnq%LJ_epsilon , sim2%opnq%LJ_epsilon )
    print *, "cloning div data"     
    sim2%div= sim1%div
        call clone(sim1%div%all_atom_numbers , sim2%div%all_atom_numbers )
    print *, "cloning vsolv data"     
    sim2%vsolv= sim1%vsolv
        call clone(sim1%vsolv%fixed_iqmatoms , sim2%vsolv%fixed_iqmatoms )
        call clone(sim1%vsolv%solvent_pointers , sim2%vsolv%solvent_pointers )
        call clone(sim1%vsolv%nearest_solvent_pointers , sim2%vsolv%nearest_solvent_pointers )
        call clone(sim1%vsolv%nearest_solvent_distances , sim2%vsolv%nearest_solvent_distances )
        call clone(sim1%vsolv%iibh , sim2%vsolv%iibh )
        call clone(sim1%vsolv%ijbh , sim2%vsolv%ijbh )
        call clone(sim1%vsolv%icbh , sim2%vsolv%icbh )
        call clone(sim1%vsolv%iiba , sim2%vsolv%iiba )
        call clone(sim1%vsolv%ijba , sim2%vsolv%ijba )
        call clone(sim1%vsolv%icba , sim2%vsolv%icba )
        call clone(sim1%vsolv%iith , sim2%vsolv%iith )
        call clone(sim1%vsolv%ijth , sim2%vsolv%ijth )
        call clone(sim1%vsolv%ikth , sim2%vsolv%ikth )
        call clone(sim1%vsolv%icth , sim2%vsolv%icth )
        call clone(sim1%vsolv%iita , sim2%vsolv%iita )
        call clone(sim1%vsolv%ijta , sim2%vsolv%ijta )
        call clone(sim1%vsolv%ikta , sim2%vsolv%ikta )
        call clone(sim1%vsolv%icta , sim2%vsolv%icta )
        call clone(sim1%vsolv%iiph , sim2%vsolv%iiph )
        call clone(sim1%vsolv%ijph , sim2%vsolv%ijph )
        call clone(sim1%vsolv%ikph , sim2%vsolv%ikph )
        call clone(sim1%vsolv%ilph , sim2%vsolv%ilph )
        call clone(sim1%vsolv%icph , sim2%vsolv%icph )
        call clone(sim1%vsolv%iipa , sim2%vsolv%iipa )
        call clone(sim1%vsolv%ijpa , sim2%vsolv%ijpa )
        call clone(sim1%vsolv%ikpa , sim2%vsolv%ikpa )
        call clone(sim1%vsolv%ilpa , sim2%vsolv%ilpa )
        call clone(sim1%vsolv%icpa , sim2%vsolv%icpa )
    print *, "cloning scratch data"     
    sim2%scratch=sim1%scratch
        call clone(sim1%scratch%matsize_red_scratch , sim2%scratch%matsize_red_scratch )
        call clone(sim1%scratch%qm_pme_scratch , sim2%scratch%qm_pme_scratch )
        call clone(sim1%scratch%mat_diag_workspace , sim2%scratch%mat_diag_workspace )
        call clone(sim1%scratch%pdiag_scr_norbs_norbs , sim2%scratch%pdiag_scr_norbs_norbs )
        call clone(sim1%scratch%pdiag_scr_noccupied_norbs , sim2%scratch%pdiag_scr_noccupied_norbs )
        call clone(sim1%scratch%pdiag_vectmp1 , sim2%scratch%pdiag_vectmp1 )
        call clone(sim1%scratch%pdiag_vectmp2 , sim2%scratch%pdiag_vectmp2 )
        call clone(sim1%scratch%pdiag_vectmp3 , sim2%scratch%pdiag_vectmp3 )
        call clone(sim1%scratch%pdiag_vecjs , sim2%scratch%pdiag_vecjs )
        call clone(sim1%scratch%lapack_dc_real_scr , sim2%scratch%lapack_dc_real_scr )
        call clone(sim1%scratch%lapack_dc_int_scr , sim2%scratch%lapack_dc_int_scr )
        call clone(sim1%scratch%qm_real_scratch , sim2%scratch%qm_real_scratch )
        call clone(sim1%scratch%qm_int_scratch , sim2%scratch%qm_int_scratch )
    print *, "cloning dav data"     
    sim2%dav=sim1%dav
        if(associated(sim1%scratch%mat_diag_workspace)) then
        if(associated(sim1%dav%ehf, sim1%scratch%mat_diag_workspace(1:sim1%qm2%norbs,1))) then 
                sim2%dav%ehf=>sim2%scratch%mat_diag_workspace(1:sim2%qm2%norbs,1)
        else if(associated(sim1%dav%ehf)) then
                print *, "Why is this pointing to 8"
        endif
        endif
        if(associated(sim1%qm2%eigen_vectors)) then
        if(associated(sim1%dav%vhf, sim1%qm2%eigen_vectors)) then 
                sim2%dav%vhf=>sim2%qm2%eigen_vectors
        else if(associated(sim1%dav%vhf)) then
                print *, "Why is this pointing to 9"
        endif
        endif
        if(associated(sim1%qm2%qm_qm_2e_repul)) then
        if(associated(sim1%dav%W, sim1%qm2%qm_qm_2e_repul)) then 
                sim2%dav%W=>sim2%qm2%qm_qm_2e_repul
        else if(associated(sim1%dav%W)) then
                print *, "Why is this pointing to 10"
        endif
        endif
    print *, "cloning ewald data"     
    sim2%ewald=sim1%ewald
        call clone(sim1%ewald%kvec , sim2%ewald%kvec )
        call clone(sim1%ewald%dkvec , sim2%ewald%dkvec )
        call clone(sim1%ewald%dmkv , sim2%ewald%dmkv )
        call clone(sim1%ewald%ktable , sim2%ewald%ktable )
        call clone(sim1%ewald%qmktable , sim2%ewald%qmktable )
        call clone(sim1%ewald%mmpot , sim2%ewald%mmpot )
        call clone(sim1%ewald%qmpot , sim2%ewald%qmpot )
        call clone(sim1%ewald%coulpot , sim2%ewald%coulpot )
        call clone(sim1%ewald%d_ewald_mm , sim2%ewald%d_ewald_mm )
    if(associated(sim1%qomp)) then
    print *, "cloning qomp data"     
    sim2%qomp= sim1%qomp
    endif

    print *, "cloning aimc data"     
    sim2%aimc=sim1%aimc
        call clone(sim1%aimc%FM,sim2%aimc%FM)
        call clone(sim1%aimc%FE,sim2%aimc%FE)
        call clone(sim1%aimc%Fmax,sim2%aimc%Fmax)
    return
end subroutine

end module clone_module
