# Flag for gfortran compiler
FFLAG_GNU = -ffree-line-length-512
# Special options for gfortran version -- >10 need to allow argument mismatch
GCC_GTEQ_10 := $(shell expr `gfortran -dumpfullversion -dumpversion | cut -f1 -d.` \>= 10)
ifeq "$(GCC_GTEQ_10)" "1"
	FFLAG_GNU += -mcmodel=small -fallow-argument-mismatch
else
	FFLAG_GNU += -mcmodel=medium
endif

MODCMD = -module
NAESMDDIR = naesmd
SQMDIR = sqm
SANDERDIR = sander
DIRSFF = sff
DIRDRIVERSRC = driver_src
LIBDIR = lib

OBJDIR = obj
SRCDIR = src
MODDIR = mod

SQMINC = -I$(MODDIR)/$(SQMDIR)/
AMOEBAINC = -I$(MODDIR)/$(AMOEBADIR)/
MODOPT = -J

INC= -Iinc/ -Iinc/$(SANDERDIR) -Iinc/$(LIBDIR)/ -Iinc/$(SQMDIR) -Imod/$(SQMDIR)/ -Imod/$(NAESMDDIR)

DIRECTORIES= $(MODDIR)/$(SQMDIR) $(MODDIR)/$(NAESMDDIR) inc/$(SQMDIR) inc/$(SANDERDIR) inc/$(LIBDIR) $(OBJDIR)/$(SQMDIR) $(OBJDIR)/$(NAESMDDIR) $(OBJDIR)/$(SANDERDIR) $(OBJDIR)/$(DIRSFF) $(OBJDIR)/$(DIRSFF)/$(DIRPUBPME) $(OBJDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC) $(OBJDIR)/$(LIBDIR)

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
	$(OBJDIR)/$(SQMDIR)/qmmm_struct_module.o \
	$(OBJDIR)/$(SQMDIR)/file_io_dat.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_vsolv_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_params_module.o \
	$(OBJDIR)/$(SQMDIR)/rotation.o \
	$(OBJDIR)/$(SQMDIR)/qmmm_module.o \
	$(OBJDIR)/$(SQMDIR)/qm_gb.o \
	$(OBJDIR)/$(SQMDIR)/slater_overlap.o \
	$(OBJDIR)/$(SQMDIR)/MNDOChargeSeparation.o \
	$(OBJDIR)/$(SQMDIR)/qm2_fock_d.o \
	$(OBJDIR)/$(SQMDIR)/dh_correction_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_davidson_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_pm6_hof_module.o \
	$(OBJDIR)/$(SQMDIR)/qm2_get_qm_forces.o \
	$(OBJDIR)/$(NAESMDDIR)/dcart_xpm_module.o \
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
	$(OBJDIR)/$(SQMDIR)/printNM.o \
        $(OBJDIR)/$(SQMDIR)/polarizab0ab_new.o
OBJNAESMD = \
        $(OBJDIR)/$(NAESMDDIR)/naesmd_module.o\
	$(OBJDIR)/$(NAESMDDIR)/apc.o \
        $(OBJDIR)/$(NAESMDDIR)/md_module.o\
        $(OBJDIR)/$(NAESMDDIR)/naesmd_constants.o\
        $(OBJDIR)/$(NAESMDDIR)/aimc_type.o \
        $(OBJDIR)/$(NAESMDDIR)/communism.o \
        $(OBJDIR)/$(NAESMDDIR)/additional-subroutines.o \
	$(OBJDIR)/$(NAESMDDIR)/random.o \
	$(OBJDIR)/$(NAESMDDIR)/langevin-temperature.o \
        $(OBJDIR)/$(NAESMDDIR)/statespecific.o \
        $(OBJDIR)/$(NAESMDDIR)/nacr.o \
	$(OBJDIR)/$(NAESMDDIR)/nacT_analytic.o \
	$(OBJDIR)/$(NAESMDDIR)/fewest-switches.o \
	$(OBJDIR)/$(NAESMDDIR)/coherence.o \
	$(OBJDIR)/$(NAESMDDIR)/cadiab.o \
	$(OBJDIR)/$(NAESMDDIR)/quantum-prop-add.o \
	$(OBJDIR)/$(NAESMDDIR)/writeoutput.o \
	$(OBJDIR)/$(NAESMDDIR)/freezing.o\
	$(OBJDIR)/$(NAESMDDIR)/verlet.o \
        $(OBJDIR)/$(NAESMDDIR)/sqm_subs.o \
	$(OBJDIR)/$(NAESMDDIR)/deriv.o \
	$(OBJDIR)/$(NAESMDDIR)/nuclear_normal_modes.o \
	$(OBJDIR)/$(NAESMDDIR)/buildM.o \
	$(OBJDIR)/$(NAESMDDIR)/liouv_new.o \
        $(OBJDIR)/$(NAESMDDIR)/rescaleveloc.o\
        $(OBJDIR)/$(NAESMDDIR)/clone.o\
        $(OBJDIR)/$(NAESMDDIR)/aimc.o \
	$(OBJDIR)/$(NAESMDDIR)/check_files.o \
	$(OBJDIR)/$(NAESMDDIR)/mce_prop.o \
        $(OBJDIR)/$(NAESMDDIR)/dropout_module.o \
	$(OBJDIR)/$(NAESMDDIR)/main.o 

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

OBJSFF = \
	$(OBJDIR)/$(DIRSFF)/xminC.o

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

$(OBJDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC)/%.o: $(SRCDIR)/$(DIRSFF)/$(DIRPUBPME)/$(DIRDRIVERSRC)/%.f
	$(FC) $(INC) $(FFLAG) -o $@ -c $<

$(OBJDIR)/$(DIRSFF)/%.o: $(SRCDIR)/$(DIRSFF)/%.c
	$(CC) $(INC) $(CFLAG) -DSQM -o $@ -c $<

$(OBJDIR)/$(LIBDIR)/%.o: $(SRCDIR)/$(LIBDIR)/%.c
	$(CC) $(INC) $(CFLAG) -o $@ -c $<

$(OBJDIR)/$(SANDERDIR)/%.o: $(SRCDIR)/$(SANDERDIR)/%.F90
	$(FC) $(INC) $(FFLAG) $(SQMINC) ${MODOPT}$(MODDIR)/$(SANDERDIR)/ -o $@ -c $<

$(OBJDIR)/$(SANDERDIR)/%.o: $(SRCDIR)/$(SANDERDIR)/%.c
	$(CC) $(INC) $(CFLAG) -o $@ -c $<

$(OBJDIR)/$(SQMDIR)/%.o: $(SRCDIR)/$(SQMDIR)/%.F90
	$(FC) $(INC) $(FFLAG) -DSQM ${MODOPT}$(MODDIR)/$(SQMDIR)/ -o $@ -c $< 

$(OBJDIR)/$(NAESMDDIR)/%.o: $(SRCDIR)/$(NAESMDDIR)/%.F90
	$(FC) $(INC) $(FFLAG) $(SQMINC) ${MODOPT}$(MODDIR)/$(NAESMDDIR)/  -o $@ -c $<  

$(OBJDIR)/$(NAESMDDIR)/old/%.o: $(SRCDIR)/$(NAESMDDIR)/old/%.F90
	$(FC) $(INC) $(FFLAG) -o $@ -c $< 

$(OBJDIR)/$(NAESMDDIR)/old/%.o: $(SRCDIR)/$(NAESMDDIR)/old/%.f
	$(FC) $(INC) $(FFLAG) -o $@ -c $<


ic_mkl: FC = ifort
ic_mkl: CC = icc
ic_mkl: MODOPT = -module 
ic_mkl: LINALG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
ic_mkl: FFLAG= -O2 -I${MKLROOT}/include
ic_mkl: CFLAG= -O2 -I${MKLROOT}/include -DMKL_LP64
ic_mkl: LDFLAGS = $(FFLAG)
ic_mkl: nexmd.exe

ic:   FC = ifort
ic:   CC = icc
ic:   MODOPT = -module 
ic:   LINALG = -llapack -lblas
ic:   FFLAG = -O2 -mcmodel=medium
ic:   CFLAG = -O2 -mcmodel=medium
ic:   LDFLAGS = $(FFLAG)
ic:   nexmd.exe

debug_ic: FC = ifort
debug_ic: CC = icc
debug_ic: MODOPT = -module 
debug_ic: LINALG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
debug_ic: FFLAG= -g -I${MKLROOT}/include
debug_ic: CFLAG= -g -I${MKLROOT}/include -DMKL_LP64
debug_ic: LDFLAGS = $(FFLAG)
debug_ic: nexmd.exe

### archive performance_ic options used in previous work
#performance_ic: FC = ifort
#performance_ic: CC = icc
#performance_ic: MODOPT = -module 
#performance_ic: LINALG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
#performance_ic: FFLAG= -O2 -g -I${MKLROOT}/include
#performance_ic: CFLAG= -O2 -g -I${MKLROOT}/include -DMKL_LP64
#performance_ic: nexmd.exe

gnu:  FC = gfortran
gnu:  CC = gcc
gnu:  LINALG = -llapack -lblas
gnu:  FFLAG = -O2 $(FFLAG_GNU) -ffpe-summary='none'
gnu:  CFLAG = -O2 
gnu:  LDFLAGS = $(FFLAG)
gnu:  nexmd.exe

gnu_mkl: FC = gfortran
gnu_mkl: CC = gcc
gnu_mkl: LINALG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
gnu_mkl: FFLAG= -O2 $(FFLAG_GNU) -I${MKLROOT}/include
gnu_mkl: CFLAG= -O2 -I${MKLROOT}/include -DMKL_LP64
gnu_mkl: LDFLAGS = $(FFLAG)
gnu_mkl: nexmd.exe

gnu_debug:  FC = gfortran
gnu_debug:  CC = gcc
gnu_debug:  LINALG = -llapack -lblas
gnu_debug:  FFLAG = -g $(FFLAG_GNU) -ffpe-trap=invalid,zero,overflow,underflow
gnu_debug:  CFLAG = -g 
gnu_debug:  LDFLAGS = $(FFLAG)
gnu_debug:  nexmd.exe

### custom field
custom:  FC = 
custom:  CC = 
custom:  LINALG = 
custom:  FFLAG = 
custom:  CFLAG =
custom:  LDFLAGS = $(FFLAG)
custom:  nexmd.exe

LINK =  $(LINALG)

nexmd.exe: $(OBJSQM) $(OBJLIB) $(OBJNAESMD) $(OBJSFF)
	$(FC) $(LDFLAGS) -o nexmd.exe $(OBJNAESMD) $(OBJSQM) $(OBJLIB) $(OBJSFF) $(LINK) 
		
clean :
	rm -f ob*/*.o obj/*/*.o  mod/*.mod *.mod mod/*/*.mod rm lib/*.a
