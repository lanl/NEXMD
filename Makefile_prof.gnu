
FFLAG = -O3 -pg
#-mcmodel=medium
LDFLAGS = $(FFLAG)
LINK =  -llapack -lblas

FC = gfortran
CC = gcc

NAESMDDIR = naesmd
SQMDIR = sqm
PMEMDDIR = pmemd
PBSADIR = pbsa
SANDERDIR = sander
DIRSFF = sff
DIRPUBPME = pubpme
AMOEBADIR = amoeba
DIRDRIVERSRC = driver_src
LIBDIR = lib

OBJDIR = obj
SRCDIR = src
MODDIR = mod
LIB = lib 

SQMINC = -I$(MODDIR)/$(SQMDIR)/
AMOEBAINC = -I$(MODDIR)/$(AMOEBADIR)/

INC= -Iinc/ -Iinc/$(SANDERDIR) -Iinc/$(PBSADIR)/ -Imod/$(PBSADIR)/ -Iinc/$(LIBDIR)/ -Iinc/$(PMEMDDIR)/ -Iinc/$(NAESMDDIR)/old/ -Iinc/$(SQMDIR) -Imod/$(SQMDIR)/ -Imod/$(NAESMDDIR)

DIRECTORIES= $(MODDIR)/$(SQMDIR) $(MODDIR)/$(NAESMDDIR) $(MODDIR)/$(AMEOBADIR) $(MODDIR)/$(PBSADIR) inc/$(SQMDIR) inc/$(SANDERDIR) inc/$(PBSADIR) inc/$(LIBDIR) inc/$(PMEMDDIR) inc/$(NAESMDDIR) inc/$(NAESMDDIR)/old $(OBJDIR)/$(SQMDIR) $(OBJDIR)/$(NAESMDDIR) $(OBJDIR)/$(SANDERDIR) $(OBJDIR)/$(PMEMDIR) $(OBJDIR)/$(PBSADIR) $(OBJDIR)/$(AMOEBADIR) $(OBJDIR)/$(DIRSFF) $(OBJDIR)/$(DIRSFF)/$(DIRPUBPME) $(OBJDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC) $(OBJDIR)/$(LIBDIR)

$(shell   mkdir -p $(DIRECTORIES)) 

OBJSQM = \
	$(OBJDIR)/$(SQMDIR)/assert.o \
	$(OBJDIR)/$(SQMDIR)/mexit.o \
	$(OBJDIR)/$(SQMDIR)/findmask.o \
	$(OBJDIR)/$(SQMDIR)/constants.o \
	$(OBJDIR)/$(SQMDIR)/utilitiesModule.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_qmtheorymodule.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_nml_module.o \
	$(OBJDIR)/$(SQMDIR)/xlbomd_module.o \
	$(OBJDIR)/$(SQMDIR)/elementOrbitalIndex.o \
	$(OBJDIR)/$(SQMDIR)/parameterReader.o \
	$(OBJDIR)/$(SQMDIR)/rotation.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_struct_module.o \
	$(OBJDIR)/$(SQMDIR)/file_io_dat.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_vsolv_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_params_module.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_module.o \
	$(OBJDIR)/$(SQMDIR)/qm_gb.o \
	$(OBJDIR)/$(SQMDIR)/slater_overlap.o \
	$(OBJDIR)/$(SQMDIR)/MNDOChargeSeparation.o \
	$(OBJDIR)/$(SQMDIR)/qm2_fock_d.o \
	$(OBJDIR)/$(SQMDIR)/dh_correction_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_davidson_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_pm6_hof_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_get_qm_forces.o \
	$(OBJDIR)/$(SQMDIR)/dcart1.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_read_cm3.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_3rd_order.o \
	$(OBJDIR)/$(SQMDIR)/nmlsrc.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_dispersionread.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_gettab.o \
	$(OBJDIR)/$(SQMDIR)/qm2_core_core_repulsion.o \
	$(OBJDIR)/$(SQMDIR)/qm2_parameters.o \
	$(OBJDIR)/$(SQMDIR)/qm2_repp_d.o \
	$(OBJDIR)/$(SQMDIR)/qm2_rotate_qmqm.o \
	$(OBJDIR)/$(SQMDIR)/qm2_repp.o \
	$(OBJDIR)/$(SQMDIR)/qm_assign_atom_types.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_get_qm_forces.o \
	$(OBJDIR)/$(SQMDIR)/cosmo_C.o \
	$(OBJDIR)/$(SQMDIR)/qm2_print_energy.o \
	$(OBJDIR)/$(SQMDIR)/amopen.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_load_params.o \
	$(OBJDIR)/$(SQMDIR)/qm2_fock.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_slkode.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_self.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_slktrafo.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_skpar.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_gamma.o \
	$(OBJDIR)/$(SQMDIR)/timer_dummy.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_repulsiv.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_dispersion_egr.o \
	$(OBJDIR)/$(SQMDIR)/qm2_h1elec_d.o \
	$(OBJDIR)/$(SQMDIR)/qm2_h1elec.o \
	$(OBJDIR)/$(SQMDIR)/qm2_core_core_repulsion_dxyz.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dihed.o \
	$(OBJDIR)/$(SQMDIR)/qm2_calc_dipole.o \
	$(OBJDIR)/$(SQMDIR)/qm2_iterator_mod.o \
	$(OBJDIR)/$(SQMDIR)/qm2_diagonalizer_module.o \
	$(OBJDIR)/$(SQMDIR)/opnq_Erep.o \
	$(OBJDIR)/$(SQMDIR)/opnq_Edisp.o \
	$(OBJDIR)/$(SQMDIR)/opnq_Evdw.o \
	$(OBJDIR)/$(SQMDIR)/opnq_SwitchMod.o \
	$(OBJDIR)/$(SQMDIR)/opnq.o \
	$(OBJDIR)/$(SQMDIR)/qm2_scf.o \
	$(OBJDIR)/$(SQMDIR)/cosmo.o \
	$(OBJDIR)/$(SQMDIR)/qm2_fock_predict.o \
	$(OBJDIR)/$(SQMDIR)/qm2_calc_charges.o \
	$(OBJDIR)/$(SQMDIR)/qm_link_atoms.o \
	$(OBJDIR)/$(SQMDIR)/qm2_get_qmmm_forces.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_get_qmmm_forces.o \
	$(OBJDIR)/$(SQMDIR)/qm2_energy.o \
	$(OBJDIR)/$(SQMDIR)/qm2_calc_rij_and_eqns.o \
	$(OBJDIR)/$(SQMDIR)/qm_print_info.o \
	$(OBJDIR)/$(SQMDIR)/qm2_load_params_and_allocate.o \
	$(OBJDIR)/$(SQMDIR)/qm_zero_charges.o \
	$(OBJDIR)/$(SQMDIR)/qm2_print_charges.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_cm3.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_mulliken.o \
	$(OBJDIR)/$(SQMDIR)/qm2_setup_orb_exp.o \
	$(OBJDIR)/$(SQMDIR)/qm2_identify_peptide_links.o \
	$(OBJDIR)/$(SQMDIR)/qm2_allocate_e_repul.o \
	$(OBJDIR)/$(SQMDIR)/qm2_hcore_qmmm.o \
	$(OBJDIR)/$(SQMDIR)/qm2_hcore_qmqm.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_energy.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_scf.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_gb.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_broyden.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_fermi.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_ewevge.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_shift.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_gammamat.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_externalshift.o \
	$(OBJDIR)/$(SQMDIR)/qm2_dftb_dispersion_params.o \
	$(OBJDIR)/$(SQMDIR)/qm2_smallest_number.o \
	$(OBJDIR)/$(SQMDIR)/dipole.o \
	$(OBJDIR)/$(SQMDIR)/nacr.o \
	$(OBJDIR)/$(SQMDIR)/dcart2.o \
	$(OBJDIR)/$(SQMDIR)/qm2_get_exc_forces.o \
	$(OBJDIR)/$(SQMDIR)/calc_rhotz.o \
	$(OBJDIR)/$(SQMDIR)/md_util.o \
        $(OBJDIR)/$(SQMDIR)/md_util_new.o \
	$(OBJDIR)/$(SQMDIR)/util.o \
	$(OBJDIR)/$(SQMDIR)/davidson.o \
	$(OBJDIR)/$(SQMDIR)/liouville.o \
	$(OBJDIR)/$(SQMDIR)/timing.o \
	$(OBJDIR)/$(SQMDIR)/qm2_print_bondorders.o \
	$(OBJDIR)/$(SQMDIR)/printNM.o
OBJNAESMD = \
	$(OBJDIR)/$(NAESMDDIR)/apc.o \
        $(OBJDIR)/$(NAESMDDIR)/naesmd_constants.o\
        $(OBJDIR)/$(NAESMDDIR)/naesmd_space_module.o \
        $(OBJDIR)/$(NAESMDDIR)/additional-subroutines.o \
        $(OBJDIR)/$(NAESMDDIR)/dcart_xpm_module.o \
        $(OBJDIR)/$(NAESMDDIR)/fcn.o \
        $(OBJDIR)/$(NAESMDDIR)/quantum-prop.o \
        $(OBJDIR)/$(NAESMDDIR)/response.o \
	$(OBJDIR)/$(NAESMDDIR)/random.o \
	$(OBJDIR)/$(NAESMDDIR)/langevin-temperature.o \
	$(OBJDIR)/$(NAESMDDIR)/communism.o \
        $(OBJDIR)/$(NAESMDDIR)/statespecific.o \
	$(OBJDIR)/$(NAESMDDIR)/nacT_analytic.o \
	$(OBJDIR)/$(NAESMDDIR)/fewest-switches.o \
	$(OBJDIR)/$(NAESMDDIR)/coherence.o \
	$(OBJDIR)/$(NAESMDDIR)/quantum-prop-add.o \
	$(OBJDIR)/$(NAESMDDIR)/writeoutput.o \
	$(OBJDIR)/$(NAESMDDIR)/verlet.o \
	$(OBJDIR)/$(NAESMDDIR)/cadiab.o \
	$(OBJDIR)/$(NAESMDDIR)/main.o \
        $(OBJDIR)/$(NAESMDDIR)/sqm_subs.o \
	$(OBJDIR)/$(NAESMDDIR)/deriv.o \
	$(OBJDIR)/$(NAESMDDIR)/buildM.o \
	$(OBJDIR)/$(NAESMDDIR)/liouv_new.o
OBJPBSA = \
	$(OBJDIR)/$(PBSADIR)/timer.o \
	$(OBJDIR)/$(PBSADIR)/sa_driver.o \
	$(OBJDIR)/$(PBSADIR)/decomp.o \
	$(OBJDIR)/$(PBSADIR)/pb_force.o \
        $(OBJDIR)/$(PBSADIR)/pb_write.o \
	$(OBJDIR)/$(PBSADIR)/svbksb.o \
	$(OBJDIR)/$(PBSADIR)/svdcmp.o \
	$(OBJDIR)/$(PBSADIR)/pythag.o \
	$(OBJDIR)/$(PBSADIR)/gen_dx_file.o \
	$(OBJDIR)/$(PBSADIR)/pb_augdrv.o \
	$(OBJDIR)/$(PBSADIR)/gmresX.o \
	$(OBJDIR)/$(PBSADIR)/aug_iccg.o \
	$(OBJDIR)/$(PBSADIR)/matvec3.o \
	$(OBJDIR)/$(PBSADIR)/iimod.o \
	$(OBJDIR)/$(PBSADIR)/interpX.o \
	$(OBJDIR)/$(PBSADIR)/dsvdc.o \
	$(OBJDIR)/$(PBSADIR)/pb_fftsolv.o \
	$(OBJDIR)/$(PBSADIR)/irre32.o \
	$(OBJDIR)/$(PBSADIR)/problem.o \
	$(OBJDIR)/$(PBSADIR)/regular.o \
	$(OBJDIR)/$(PBSADIR)/qld.o \
	$(OBJDIR)/$(PBSADIR)/miniop.o \
	$(OBJDIR)/$(PBSADIR)/jumps.o \
	$(OBJDIR)/$(PBSADIR)/irre31.o \
	$(OBJDIR)/$(PBSADIR)/wint.o \
	$(OBJDIR)/$(PBSADIR)/coed6.o \
	$(OBJDIR)/$(PBSADIR)/qint.o \
	$(OBJDIR)/$(PBSADIR)/coed20.o \
	$(OBJDIR)/$(PBSADIR)/prodis.o \
	$(OBJDIR)/$(PBSADIR)/indexg.o \
	$(OBJDIR)/$(PBSADIR)/pb_iimdrv.o \
	$(OBJDIR)/$(PBSADIR)/IIM.o \
	$(OBJDIR)/$(PBSADIR)/bicg.o \
	$(OBJDIR)/$(PBSADIR)/dslubc.o \
	$(OBJDIR)/$(PBSADIR)/dbcg.o \
	$(OBJDIR)/$(PBSADIR)/isdbcg.o \
	$(OBJDIR)/$(PBSADIR)/d1mach.o \
	$(OBJDIR)/$(PBSADIR)/xermsg.o \
	$(OBJDIR)/$(PBSADIR)/xerhlt.o \
	$(OBJDIR)/$(PBSADIR)/xersve.o \
	$(OBJDIR)/$(PBSADIR)/i1mach.o \
	$(OBJDIR)/$(PBSADIR)/xgetua.o \
	$(OBJDIR)/$(PBSADIR)/j4save.o \
	$(OBJDIR)/$(PBSADIR)/xerprn.o \
	$(OBJDIR)/$(PBSADIR)/fdump.o \
	$(OBJDIR)/$(PBSADIR)/xercnt.o \
	$(OBJDIR)/$(PBSADIR)/dsluti.o \
	$(OBJDIR)/$(PBSADIR)/dslui4.o \
	$(OBJDIR)/$(PBSADIR)/dslui.o \
	$(OBJDIR)/$(PBSADIR)/dslui2.o \
	$(OBJDIR)/$(PBSADIR)/dsmtv.o \
	$(OBJDIR)/$(PBSADIR)/dsmv.o \
	$(OBJDIR)/$(PBSADIR)/dsilus.o \
	$(OBJDIR)/$(PBSADIR)/dchkw.o \
	$(OBJDIR)/$(PBSADIR)/ds2y.o \
	$(OBJDIR)/$(PBSADIR)/qs2i1d.o \
	$(OBJDIR)/$(PBSADIR)/gmres.o \
	$(OBJDIR)/$(PBSADIR)/dslugm.o \
	$(OBJDIR)/$(PBSADIR)/dgmres.o \
	$(OBJDIR)/$(PBSADIR)/dpigmr.o \
	$(OBJDIR)/$(PBSADIR)/dhels.o \
	$(OBJDIR)/$(PBSADIR)/drlcal.o \
	$(OBJDIR)/$(PBSADIR)/isdgmr.o \
	$(OBJDIR)/$(PBSADIR)/dheqr.o \
	$(OBJDIR)/$(PBSADIR)/dxlcal.o \
	$(OBJDIR)/$(PBSADIR)/dorth.o \
	$(OBJDIR)/$(PBSADIR)/curv.o \
	$(OBJDIR)/$(PBSADIR)/GrToPr.o \
	$(OBJDIR)/$(PBSADIR)/transf.o \
	$(OBJDIR)/$(PBSADIR)/pb_lsolver.o \
	$(OBJDIR)/$(PBSADIR)/pb_nlsolver.o \
	$(OBJDIR)/$(PBSADIR)/pb_fddrv.o \
	$(OBJDIR)/$(PBSADIR)/pb_exmol.o \
	$(OBJDIR)/$(PBSADIR)/membrane.o \
	$(OBJDIR)/$(PBSADIR)/pb_direct.o \
	$(OBJDIR)/$(PBSADIR)/pb_mpfrc.o \
	$(OBJDIR)/$(PBSADIR)/pb_list.o \
	$(OBJDIR)/$(PBSADIR)/periodic_cg/random.o \
	$(OBJDIR)/$(PBSADIR)/np_force.o

OBJSANDER = \
	$(OBJDIR)/$(SANDERDIR)/amoeba_mdin.o \
	$(OBJDIR)/$(SANDERDIR)/trace.o \
	$(OBJDIR)/$(SANDERDIR)/nonbond_list.o \
	$(OBJDIR)/$(SANDERDIR)/binrestart.o \
	$(OBJDIR)/$(SANDERDIR)/ew_box.o \
	$(OBJDIR)/$(SANDERDIR)/qm_ewald.o \
	$(OBJDIR)/$(SANDERDIR)/erfcfun.o \
	$(OBJDIR)/$(SANDERDIR)/stack.o \
	$(OBJDIR)/$(SANDERDIR)/state.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_runmd.o \
	$(OBJDIR)/$(SANDERDIR)/ew_bspline.o \
	$(OBJDIR)/$(SANDERDIR)/sander_lib.o \
	$(OBJDIR)/$(SANDERDIR)/spatial_fft.o \
	$(OBJDIR)/$(SANDERDIR)/ew_setup.o \
	$(OBJDIR)/$(SANDERDIR)/ew_recip_reg.o \
	$(OBJDIR)/$(SANDERDIR)/ew_fft.o \
	$(OBJDIR)/$(SANDERDIR)/parms.o \
	$(OBJDIR)/$(SANDERDIR)/softcore.o \
	$(OBJDIR)/$(SANDERDIR)/decomp.o \
	$(OBJDIR)/$(SANDERDIR)/icosasurf.o \
	$(OBJDIR)/$(SANDERDIR)/linear_response.o \
	$(OBJDIR)/$(SANDERDIR)/locmem.o \
	$(OBJDIR)/$(SANDERDIR)/crg_reloc.o \
	$(OBJDIR)/$(SANDERDIR)/charmm.o \
	$(OBJDIR)/$(SANDERDIR)/ew_recip.o \
	$(OBJDIR)/$(SANDERDIR)/debug.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-constants.o \
	$(OBJDIR)/$(SANDERDIR)/sglds.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-sander-proxy.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-colvar-type.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-utils.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-colvar-math.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-colvar-utils.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-ANGLE.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-TORSION.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-DISTANCE.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-rmsd.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-MULTI_RMSD.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-R_OF_GYRATION.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-HANDEDNESS.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-N_OF_BONDS.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-N_OF_STRUCTURES.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-LCOD.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-COS_OF_DIHEDRAL.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-COM_ANGLE.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-COM_TORSION.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-COM_DISTANCE.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cv-PCA.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-value.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-cftree.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-read-pca.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-colvar.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-lexer.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-parser.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-pmd-hooks.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-umbrella.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-abmd-hooks.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-smd-hooks.o \
	$(OBJDIR)/$(SANDERDIR)/ncsu-sander-hooks.o \
	$(OBJDIR)/$(SANDERDIR)/egb.o \
	$(OBJDIR)/$(SANDERDIR)/pimd_vars.o \
	$(OBJDIR)/$(SANDERDIR)/relax_mat.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_multipoles.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_induced.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_vdw.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_adjust.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_recip.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_valence.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_direct.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_self.o \
	$(OBJDIR)/$(SANDERDIR)/amoeba_interface.o \
	$(OBJDIR)/$(SANDERDIR)/amd.o \
	$(OBJDIR)/$(SANDERDIR)/ips.o \
	$(OBJDIR)/$(SANDERDIR)/memory_module.o \
	$(OBJDIR)/$(SANDERDIR)/emap.o \
	$(OBJDIR)/$(SANDERDIR)/xref.o \
	$(OBJDIR)/$(SANDERDIR)/force.o \
	$(OBJDIR)/$(SANDERDIR)/pimd_force.o \
	$(OBJDIR)/$(SANDERDIR)/ene.o \
	$(OBJDIR)/$(SANDERDIR)/nmr.o \
	$(OBJDIR)/$(SANDERDIR)/nmrcal.o \
	$(OBJDIR)/$(SANDERDIR)/set.o \
	$(OBJDIR)/$(SANDERDIR)/decnvh.o \
	$(OBJDIR)/$(SANDERDIR)/align.o \
	$(OBJDIR)/$(SANDERDIR)/csa.o \
	$(OBJDIR)/$(SANDERDIR)/pcshift.o \
	$(OBJDIR)/$(SANDERDIR)/pearsn.o \
	$(OBJDIR)/$(SANDERDIR)/cshf.o \
	$(OBJDIR)/$(SANDERDIR)/multitmd.o \
	$(OBJDIR)/$(SANDERDIR)/mtmdcall.o \
	$(OBJDIR)/$(SANDERDIR)/ew_dipole_recip.o \
	$(OBJDIR)/$(SANDERDIR)/spatial_recip.o \
	$(OBJDIR)/$(SANDERDIR)/ew_force.o \
	$(OBJDIR)/$(SANDERDIR)/ew_handle_dips.o \
	$(OBJDIR)/$(SANDERDIR)/extra_pts.o \
	$(OBJDIR)/$(SANDERDIR)/short_ene.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_adf_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_gms_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_tc_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_gau_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_orc_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_nw_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_extern_module.o \
	$(OBJDIR)/$(SANDERDIR)/qm_mm.o \
	$(OBJDIR)/$(SANDERDIR)/qm2_read_adf_results.o \
	$(OBJDIR)/$(SANDERDIR)/KFReader.o \
	$(OBJDIR)/$(SANDERDIR)/ArrayList.o \
	$(OBJDIR)/$(SANDERDIR)/molecule.o \
	$(OBJDIR)/$(SANDERDIR)/qmmm_vsolv.o \
	$(OBJDIR)/$(SANDERDIR)/iwrap2.o \
	$(OBJDIR)/$(SANDERDIR)/rgroup.o \
	$(OBJDIR)/$(SANDERDIR)/bintraj.o \
	$(OBJDIR)/$(SANDERDIR)/cmd_vars.o \
	$(OBJDIR)/$(SANDERDIR)/dynlib.o \
	$(OBJDIR)/$(SANDERDIR)/matinv.o \
	$(OBJDIR)/$(SANDERDIR)/constantph.o \
	$(OBJDIR)/$(SANDERDIR)/new_time.o \
	$(OBJDIR)/$(SANDERDIR)/lscivr_vars.o \
	$(OBJDIR)/$(SANDERDIR)/nose_hoover.o \
	$(OBJDIR)/$(SANDERDIR)/fastwt.o \
	$(OBJDIR)/$(SANDERDIR)/nose_hoover_vars.o \
	$(OBJDIR)/$(SANDERDIR)/runmd.o \
	$(OBJDIR)/$(SANDERDIR)/prn_dipoles.o \
	$(OBJDIR)/$(SANDERDIR)/prn_qmmm_dipole.o \
	$(OBJDIR)/$(SANDERDIR)/mdwrit.o \
	$(OBJDIR)/$(SANDERDIR)/cmd_matrix_a1st.o \
	$(OBJDIR)/$(SANDERDIR)/cmd_matrix.o \
	$(OBJDIR)/$(SANDERDIR)/shake.o \
	$(OBJDIR)/$(SANDERDIR)/pimd_init.o \
	$(OBJDIR)/$(SANDERDIR)/lsc_xp.o \
	$(OBJDIR)/$(SANDERDIR)/lsc_init.o \
	$(OBJDIR)/$(SANDERDIR)/nose_hoover_init.o \
	$(OBJDIR)/$(SANDERDIR)/degcnt.o 
OBJPMEMD = \
	$(OBJDIR)/$(PMEMDDIR)/veclib.o \
	$(OBJDIR)/$(PMEMDDIR)/file_io_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/gbl_constants.o \
	$(OBJDIR)/$(PMEMDDIR)/parallel_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/pmemd_lib.o \
	$(OBJDIR)/$(PMEMDDIR)/gbl_datatypes.o \
	$(OBJDIR)/$(PMEMDDIR)/mol_list.o \
	$(OBJDIR)/$(PMEMDDIR)/dynamics_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/prfs.o \
	$(OBJDIR)/$(PMEMDDIR)/random.o \
	$(OBJDIR)/$(PMEMDDIR)/file_io.o \
	$(OBJDIR)/$(PMEMDDIR)/mdin_ctrl_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/dynamics.o \
	$(OBJDIR)/$(PMEMDDIR)/nextprmtop_section.o \
	$(OBJDIR)/$(PMEMDDIR)/charmm.o \
	$(OBJDIR)/$(PMEMDDIR)/prmtop_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/pbc.o \
	$(OBJDIR)/$(PMEMDDIR)/axis_optimize.o \
	$(OBJDIR)/$(PMEMDDIR)/binrestart.o \
	$(OBJDIR)/$(PMEMDDIR)/inpcrd_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/bintraj.o \
	$(OBJDIR)/$(PMEMDDIR)/remd.o \
	$(OBJDIR)/$(PMEMDDIR)/timers.o \
	$(OBJDIR)/$(PMEMDDIR)/state_info.o \
	$(OBJDIR)/$(PMEMDDIR)/fft1d.o \
	$(OBJDIR)/$(PMEMDDIR)/mdin_ewald_dat.o \
	$(OBJDIR)/$(PMEMDDIR)/nmr_lib.o \
	$(OBJDIR)/$(PMEMDDIR)/nmr_calls.o \
	$(OBJDIR)/$(PMEMDDIR)/runfiles.o \
	$(OBJDIR)/$(PMEMDDIR)/pmemd_clib.o \
	$(OBJDIR)/$(PMEMDDIR)/erfcfun.o 

OBJSFF = \
	$(OBJDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC)/utility.o \
	$(OBJDIR)/$(DIRSFF)/xminC.o
OBJSFF = \
	$(OBJDIR)/$(DIRSFF)/xminC.o 

OBJAMOEBA = \
	$(OBJDIR)/$(AMOEBADIR)/nextprmtop_section.o \
	$(OBJDIR)/$(AMOEBADIR)/state_info.o \
	$(OBJDIR)/$(AMOEBADIR)/file_io_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/parallel_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/pmemd_lib.o \
	$(OBJDIR)/$(AMOEBADIR)/gbl_constants.o \
	$(OBJDIR)/$(AMOEBADIR)/fft1d.o \
	$(OBJDIR)/$(AMOEBADIR)/axis_optimize.o \
	$(OBJDIR)/$(AMOEBADIR)/file_io.o \
	$(OBJDIR)/$(AMOEBADIR)/mdin_ctrl_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/mdin_ewald_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/gbl_datatypes.o \
	$(OBJDIR)/$(AMOEBADIR)/prmtop_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/nmr_lib.o \
	$(OBJDIR)/$(AMOEBADIR)/nmr_calls.o \
	$(OBJDIR)/$(AMOEBADIR)/dynamics_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/random.o \
	$(OBJDIR)/$(AMOEBADIR)/dynamics.o \
	$(OBJDIR)/$(AMOEBADIR)/pbc.o \
	$(OBJDIR)/$(AMOEBADIR)/inpcrd_dat.o \
	$(OBJDIR)/$(AMOEBADIR)/bintraj.o \
	$(OBJDIR)/$(AMOEBADIR)/loadbal.o \
	$(OBJDIR)/$(AMOEBADIR)/runfiles.o \
	$(OBJDIR)/$(AMOEBADIR)/erfcfun.o \
	$(OBJDIR)/$(AMOEBADIR)/pmemd_clib.o

OBJLIB = \
	$(OBJDIR)/$(LIBDIR)/rfree.o \
	$(OBJDIR)/$(LIBDIR)/nxtsec.o \
	$(OBJDIR)/$(LIBDIR)/veclib.o \
	$(OBJDIR)/$(LIBDIR)/sys.o \
	$(OBJDIR)/$(LIBDIR)/wallclock.o \
	$(OBJDIR)/$(LIBDIR)/random.o

$(OBJDIR)/$(LIBDIR)/%.o: $(SRCDIR)/$(LIBDIR)/%.F90
	$(FC) $(INC) $(FFLAG) -o $@ -c $<

$(OBJDIR)/$(LIBDIR)/%.o: $(SRCDIR)/$(LIBDIR)/%.F
	$(FC) $(INC) $(FFLAG) -o $@ -c $<

$(OBJDIR)/$(AMOEBADIR)/%.o: $(SRCDIR)/$(AMOEBADIR)/%.F90
	$(FC) $(INC) $(FFLAG) -J$(MODDIR)/$(AMOEBADIR)/  -o $@ -c $<

$(OBJDIR)/$(AMOEBADIR)/%.o: $(SRCDIR)/$(AMOEBADIR)/%.c
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC)/%.o: $(SRCDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC)/%.f
	$(FC) $(INC) $(FFLAG) -o $@ -c $<


$(OBJDIR)/$(DIRSFF)/%.o: $(SRCDIR)/$(DIRSFF)/%.c
	$(CC) $(INC) -DSQM -o $@ -c $<

$(OBJDIR)/$(LIBDIR)/%.o: $(SRCDIR)/$(LIBDIR)/%.c
	$(CC) $(INC) $(FFLAG) -o $@ -c $<

$(OBJDIR)/$(SANDERDIR)/%.o: $(SRCDIR)/$(SANDERDIR)/%.F90
	$(FC) $(INC) $(FFLAG) $(SQMINC) -J$(MODDIR)/$(SANDERDIR)/ -o $@ -c $<

$(OBJDIR)/$(SANDERDIR)/%.o: $(SRCDIR)/$(SANDERDIR)/%.c
	$(CC) $(INC) $(FFLAG) -o $@ -c $<

$(OBJDIR)/$(SQMDIR)/%.o: $(SRCDIR)/$(SQMDIR)/%.F90
	$(FC) $(INC) $(FFLAG) -DSQM -J$(MODDIR)/$(SQMDIR)/ -o $@ -c $< 

$(OBJDIR)/$(PBSADIR)/%.o: $(SRCDIR)/$(PBSADIR)/%.F90
	$(FC) $(INC) $(FFLAG) -J$(MODDIR)/$(PBSADIR)/ -o $@ -c $<


$(OBJDIR)/$(PBSADIR)/periodic_cg/%.o: $(SRCDIR)/$(PBSADIR)/periodic_cg/%.F90
	$(FC) $(INC) $(FFLAG) -o $@ -c $<


$(OBJDIR)/$(NAESMDDIR)/%.o: $(SRCDIR)/$(NAESMDDIR)/%.F90
	$(FC) $(INC) $(FFLAG) $(SQMINC) -J$(MODDIR)/$(NAESMDDIR)/  -o $@ -c $<  

$(OBJDIR)/$(NAESMDDIR)/old/%.o: $(SRCDIR)/$(NAESMDDIR)/old/%.F90
	$(FC) $(INC) $(FFLAG) -o $@ -c $< 

$(OBJDIR)/$(NAESMDDIR)/old/%.o: $(SRCDIR)/$(NAESMDDIR)/old/%.f
	$(FC) $(INC) $(FFLAG) -o $@ -c $<


$(OBJDIR)/$(PMEMDDIR)/%.o: $(SRCDIR)/$(PMEMDDIR)/%.F90
	$(FC) $(INC) $(FFLAG) -J$(MODDIR)/$(PMEMDDIR) -o $@ -c $<

$(OBJDIR)/$(PMEMDDIR)/%.o: $(SRCDIR)/$(PMEMDDIR)/%.c
	$(CC) $(INC) $(FFLAG) -o $@ -c $<




sqmceonaesmd.exe: $(OBJSQM) $(OBJLIB) $(OBJNAESMD) $(OBJSFF)
	$(FC) $(LDFLAGS) -o sqmceonaesmd.exe $(OBJNAESMD) $(OBJSQM) $(OBJLIB) $(OBJSFF) -L$(LIB) $(LINK) 
		
clean :
	rm -f ob*/*.o obj/*/*.o  mod/*.mod *.mod mod/*/*.mod rm lib/*.a
