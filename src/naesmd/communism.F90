#include "dprec.fh"
#include "assert.fh"
          
module communism
    use qm2_davidson_module, only      : qm2_davidson_structure_type
    use qmmm_struct_module, only       : qmmm_struct_type
    use qmmm_module, only : qm2_structure, qm_ewald_structure, &
         qm2_rij_eqns_structure, qm_gb_structure, qmmm_mpi_structure, &
         qmmm_opnq_structure, qmmm_scratch_structure, qmmm_div_structure, &
         qmmm_openmp_structure 
    use qmmm_vsolv_module , only : qmmm_vsolv_type
    use naesmd_module, only            : naesmd_structure, realp_t
    use md_module, only            : md_structure
    use xlbomd_module, only            : xlbomd_structure
    use naesmd_constants
    use cosmo_C, only : cosmo_C_structure          
    use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type
    use AIMC_type_module, only: AIMC_type
    implicit none
          
    type:: simulation_t
        integer                        :: Nsim
        _REAL_, dimension(:), pointer  :: coords
        integer                        :: Na, nbasis
        integer                        :: Ncharge
        integer                        :: ibo
        integer                        :: id
        integer                        :: excN ! number of excited states
        integer                        :: dotrivial
        integer                        :: cohertype
        integer                        :: lprint
        integer                        :: itime1
        integer :: outfile_1,outfile_2,outfile_3,outfile_4
        integer :: outfile_5,outfile_6,outfile_7,outfile_8
        integer :: outfile_9,outfile_10,outfile_11
        integer :: outfile_13,outfile_14,outfile_15,outfile_16
        integer :: outfile_17,outfile_18,outfile_19,outfile_20
        integer :: outfile_21,outfile_22
        integer :: outfile_25,outfile_27, outfile_28
        integer :: outfile_29
        character*150, allocatable, dimension(:) :: input_line !input lines for restart
        _REAL_                         :: escf
        _REAL_, dimension(:), pointer  :: deriv_forces ! forces as they come out
        _REAL_, dimension(:,:), pointer  :: deriv_forces_state ! forces as they come out
        logical :: forces_allocated=.false.
        ! of deriv (eV/A) for each electronic state
        _REAL_:: constcoherE0,constcoherC

        ! Various cumulative times (initial values=0.d0)
        _REAL_::time_sqm_took=0.d0
        _REAL_::time_davidson_took=0.d0
        _REAL_::time_deriv_took=0.d0 ! adiabatic derivatives (i.e., forces)
        _REAL_::time_nact_took=0.d0 ! non-adiabatic derivatives (nact)
          
        type(naesmd_structure),pointer            :: naesmd
        type(md_structure),pointer                :: md
        type(qm2_davidson_structure_type),pointer :: dav
        type(qmmm_struct_type),pointer            :: qmmm
        type(qm2_structure),pointer               :: qm2
        type(xlbomd_structure),pointer            :: xlbomd
        type(cosmo_C_structure),pointer           :: cosmo
        type(qm2_params_type),pointer :: qparams
        type(qmmm_nml_type),pointer :: qnml
        type(qm2_rij_eqns_structure),pointer :: rij
        type(qm_gb_structure),pointer :: gb
        type(qmmm_mpi_structure),pointer :: qmpi
        type(qmmm_opnq_structure),pointer :: opnq
        type(qmmm_div_structure),pointer :: div
        type(qmmm_vsolv_type),pointer :: vsolv
        type(qmmm_scratch_structure),pointer:: scratch
        type(qm_ewald_structure),pointer :: ewald
        type(qmmm_openmp_structure),pointer  :: qomp 
        type(AIMC_type),pointer  :: aimc
    end type simulation_t
    
    type sim_pointer_array
        type(simulation_t),pointer :: sim
    end type sim_pointer_array

    type :: MCE
        double complex, dimension(:), allocatable :: D !nuclear coeffients for AIMC
        _REAL_, dimension(:,:,:,:), allocatable :: sE !electronic overlaps 
        double complex, dimension(:,:), allocatable :: Heff, Heff_old !Effective Hamiltonian
        logical, dimension(:), allocatable :: cloned !Flag for clones
        _REAL_, dimension(:), allocatable :: pop !Electronic populations
        !Outputs:
        integer :: outfile_0 !for debugging
        integer :: outfile_1 !for populations
        integer :: outfile_2 !for nuclear coefficients
        integer :: outfile_3 !for electronic overlaps
        integer :: outfile_4 !for last Heff (for restart)
        integer :: outfile_5 !for last Pha (for restart)
        integer :: outfile_6 !for dropped trajectories (for restart)
    end type MCE
          
contains
          
    subroutine init0_simulation(a,Nsim)
        type(simulation_t), pointer  :: a
        integer :: Nsim
        a%Na       = 0
        a%Ncharge  = 0
        a%excN     = 0
        a%dav      => null()
        a%qmmm     => null()
        a%qm2      => null()
        a%naesmd      => null()
        a%md      => null()
        a%cosmo      => null()
        a%xlbomd      => null()
        a%qparams => null() 
        a%qnml => null() 
        a%rij => null() 
        a%gb => null() 
        a%qmpi => null() 
        a%opnq => null()   
        a%div => null()   
        a%vsolv => null() 
        a%scratch => null()  
        a%ewald => null()
        a%aimc => null()
        a%Nsim = Nsim
#ifdef OPENMP 
        a%qomp => null()
#endif  
        allocate(a%naesmd)
    end subroutine

    subroutine setp_simulation(a,qmmm_struct,qm2ds,qm2_struct,naesmd_struct,md_struct, &
                cosmo_c_struct, xlbomd_struct, qparams_struct, &
                qnml_struct, rij_struct, gb_struct, qmpi_struct, &
                opnq_struct, div_struct, vsolv_struct, scratch_struct, &
                ewald_struct, aimc_struct, qomp_struct)
        type(simulation_t), pointer  :: a
        type(qmmm_struct_type),target :: qmmm_struct
        type(qm2_davidson_structure_type),target :: qm2ds
        type(qm2_structure),target :: qm2_struct
        type(xlbomd_structure),target :: xlbomd_struct
        type(cosmo_C_structure),target :: cosmo_c_struct
        type(naesmd_structure),target :: naesmd_struct
        type(md_structure),target :: md_struct
        type(qm2_params_type),target :: qparams_struct
        type(qmmm_nml_type),target :: qnml_struct
        type(qm2_rij_eqns_structure),target:: rij_struct
        type(qm_gb_structure),target :: gb_struct
        type(qmmm_mpi_structure),target :: qmpi_struct
        type(qmmm_opnq_structure),target :: opnq_struct
        type(qmmm_div_structure), target :: div_struct
        type(qmmm_vsolv_type), target :: vsolv_struct
        type(qmmm_scratch_structure), target:: scratch_struct
        type(qm_ewald_structure), target :: ewald_struct 
        type(AIMC_type), target :: aimc_struct
        type(qmmm_openmp_structure),optional,target :: qomp_struct
        a%qmmm     => qmmm_struct
        a%dav     => qm2ds
        a%qm2     => qm2_struct
        a%naesmd  => naesmd_struct
        a%md  => md_struct
        a%cosmo => cosmo_c_struct
        a%xlbomd => xlbomd_struct
        a%qparams => qparams_struct
        a%qnml => qnml_struct
        a%rij => rij_struct
        a%gb => gb_struct
        a%qmpi => qmpi_struct 
        a%opnq => opnq_struct  
        a%div => div_struct  
        a%vsolv => vsolv_struct  
        a%scratch => scratch_struct  
        a%ewald => ewald_struct
        a%aimc => aimc_struct
        if(present(qomp_struct)) a%qomp => qomp_struct 
    end subroutine
          
    !
    !********************************************************************
    !
    !  Energy conversion
    !
    !********************************************************************
    !
    subroutine dav2naesmd_Omega(sim) 
        implicit none
        type(simulation_t),pointer::sim
        integer::mn
        if (sim%dav%Mx>0) then
            mn=min(size(sim%dav%e0),size(sim%naesmd%Omega))
            sim%naesmd%Omega(1:mn) = sim%dav%e0(1:mn)/feVmdqt;
        endif
        sim%naesmd%E0=sim%dav%Eground/feVmdqt
        return
    end subroutine
    !
    !********************************************************************
    !
    subroutine qmmm2naesmd_r(sim)
        type(simulation_t), pointer :: sim
        integer :: mn
        mn = min(size(sim%naesmd%rx), size(sim%qmmm%qm_coords(1,:)))
        sim%naesmd%rx(:mn) = sim%qmmm%qm_coords(1,:mn) / convl
        sim%naesmd%ry(:mn) = sim%qmmm%qm_coords(2,:mn) / convl
        sim%naesmd%rz(:mn) = sim%qmmm%qm_coords(3,:mn) / convl
        return
    end subroutine
          
    subroutine naesmd2qmmm_r(sim)
        type(simulation_t), pointer :: sim
        integer :: mn
        mn = min(size(sim%naesmd%rx), size(sim%qmmm%qm_coords(1,:)))
        sim%qmmm%qm_coords(1,:mn) = sim%naesmd%rx(:mn) * convl
        sim%qmmm%qm_coords(2,:mn) = sim%naesmd%ry(:mn) * convl
        sim%qmmm%qm_coords(3,:mn) = sim%naesmd%rz(:mn) * convl
        return
    end subroutine
          
    subroutine xyz2qmmm_r(sim, rx, ry, rz)
        type(simulation_t), pointer :: sim
        _REAL_, intent(in) :: rx(:), ry(:), rz(:)
        integer :: mn
        mn = min(size(rx), size(sim%qmmm%qm_coords(1,:)))
        sim%qmmm%qm_coords(1,:mn) = rx(:mn) * convl
        sim%qmmm%qm_coords(2,:mn) = ry(:mn) * convl
        sim%qmmm%qm_coords(3,:mn) = rz(:mn) * convl
        return
    end subroutine
          
    !
    !********************************************************************
    !
    !  Conversion from forces to acceleration
    !  The source is forces in eV/A
    !  the result is the acceleration on atomic units
    !
    !********************************************************************
    !
    subroutine deriv2naesmd_accel(sim)
        implicit none
        type(simulation_t),pointer::sim
        integer i
          
        do i=1,sim%Na
            sim%naesmd%ax(i)=sim%deriv_forces(i*3-2) &
                *convl/feVmdqt/sim%naesmd%massmdqt(i)
          
            sim%naesmd%ay(i)=sim%deriv_forces(i*3-1) &
                *convl/feVmdqt/sim%naesmd%massmdqt(i)
          
            sim%naesmd%az(i)=sim%deriv_forces(i*3) &
                *convl/feVmdqt/sim%naesmd%massmdqt(i)
        end do
        return
    end subroutine
          
    subroutine qmmm2coords_r(sim)
        type(simulation_t), pointer :: sim
        integer :: i
          
        do i=1,sim%Na
            sim%coords(i*3-2) = sim%qmmm%qm_coords(1,i)
            sim%coords(i*3-1) = sim%qmmm%qm_coords(2,i)
            sim%coords(i*3  ) = sim%qmmm%qm_coords(3,i)
        enddo
          
    end subroutine
          
    !
    !********************************************************************
    !
    !
    !********************************************************************
    !
    subroutine dav2cmdqt(sim,cmdqt)
        implicit none
        type(simulation_t), pointer::sim
        _REAL_,intent(inout)::cmdqt(:,:)
        integer i,j
        do i=1,sim%dav%Ncis
            do j=1,sim%excN
                cmdqt(i,j)=sim%dav%cmdqt(i,j)
            end do
        end do
        return
    end subroutine
          
    !********************************************************************
    !
    !  Doing sqm (ground state) and davidson (excited states)
    !  update, i.e., recalculating energies
    !
    !  Furthermore: it checks the sign conservation (quarility)
    !  of the molecular orbitals and the transition density matrix
    !
    !********************************************************************
    !
    subroutine do_sqm_and_davidson(sim,printTdipole,rx,ry,rz,r,statelimit)
        implicit none
          
        type(simulation_t),pointer::sim
        _REAL_,intent(in),optional::rx(:),ry(:),rz(:)
        _REAL_ born_radii(sim%Na), one_born_radii(sim%Na)
        _REAL_ intdiel, extdiel, Arad
        _REAL_ f
        type(realp_t),intent(in),optional::r(3)
        integer i,j
        _REAL_ ddot
        integer :: quir_cmdqt=0
        integer,optional::statelimit
        integer printTdipole
          
          
        if(present(rx).and.present(ry).and.present(rz)) then
            call xyz2qmmm_r(sim, rx,ry,rz)
        else if(present(r)) then
            call xyz2qmmm_r(sim, r(1)%p, r(2)%p, r(3)%p)
        else
            call naesmd2qmmm_r(sim)
        end if
          
        call qmmm2coords_r(sim)
        if(sim%excN>0) then
            sim%dav%Mx=sim%excN+2
        else
            sim%dav%Mx=sim%excN
        endif
        if(present(statelimit)) then !reduces the calculated number of states
            sim%dav%Mx=statelimit
        endif

        ! CML Includes call to Davidson within sqm_energy() 7/16/12
        if(sim%naesmd%tfemto.eq.0.0d0) then
            call sqm_energy(sim%ewald,sim%opnq, sim%rij,sim%qmpi,sim%gb,sim%scratch,sim%qnml,sim%qparams, sim%xlbomd,sim%cosmo,&
                sim%qm2,sim%dav,sim%qmmm,sim%Na,sim%coords,sim%escf,born_radii, &
                one_born_radii,intdiel,extdiel,Arad,sim%qm2%scf_mchg, &
                sim%time_sqm_took,sim%time_davidson_took,printTdipole,sim%outfile_28,sim%naesmd%tfemto) !The use of sim here is a hack right now and could be fixed
        elseif(sim%naesmd%icontw.ne.sim%naesmd%nstepw) then
            call sqm_energy(sim%ewald,sim%opnq, sim%rij,sim%qmpi,sim%gb,sim%scratch,sim%qnml,sim%qparams, sim%xlbomd,sim%cosmo,&
                sim%qm2,sim%dav,sim%qmmm,sim%Na,sim%coords,sim%escf,born_radii, &
                one_born_radii,intdiel,extdiel,Arad,sim%qm2%scf_mchg, &
                sim%time_sqm_took,sim%time_davidson_took,0,sim%outfile_28,sim%naesmd%tfemto) !The use of sim here is a hack right now and could be fixed
        else
            call sqm_energy(sim%ewald,sim%opnq, sim%rij,sim%qmpi,sim%gb,sim%scratch,sim%qnml,sim%qparams, sim%xlbomd,sim%cosmo,&
                sim%qm2,sim%dav,sim%qmmm,sim%Na,sim%coords,sim%escf,born_radii, &
                one_born_radii,intdiel,extdiel,Arad,sim%qm2%scf_mchg, &
                sim%time_sqm_took,sim%time_davidson_took,printTdipole,sim%outfile_28,sim%naesmd%tfemto) !The use of sim here is a hack right now and could be fixed
        endif
        ! ground state energy
        sim%dav%Eground=sim%dav%Eground
        sim%nbasis=sim%dav%Nb

        if (.not.allocated(sim%dav%vhf_old)) allocate(sim%dav%vhf_old(sim%dav%Nb,sim%dav%Nb));
        sim%dav%vhf_old(1:sim%dav%Nb,1:sim%dav%Nb)=sim%dav%vhf(1:sim%dav%Nb,1:sim%dav%Nb);        ! updating hartree-fock vectors
        ! Checking and restoring the sign of the transition density matrix (cmdqt)
        do j=1,sim%excN
          
            f=ddot(sim%dav%Ncis, &
                sim%dav%cmdqt(1,j),1,sim%dav%v0(1,sim%dav%kx(j)),1)
          
            if(f.lt.0.d0) then
                do i=1,sim%dav%Ncis
                    sim%dav%cmdqt(i,j)=-sim%dav%v0(i,sim%dav%kx(j))
                    quir_cmdqt=quir_cmdqt+1
                end do
            else
                do i=1,sim%dav%Ncis
                    sim%dav%cmdqt(i,j)=sim%dav%v0(i,sim%dav%kx(j))
                end do
            end if
        end do
        if((quir_cmdqt>0).and.(sim%dav%verbosity>0)) then ! quirality was restored at least once above
            write(6,*) ' Restoring quirality in cmdqt'
        end if
        return
    end subroutine do_sqm_and_davidson
    !
    !********************************************************************
    !
    !  Updated ground and excited state energies
    !  Furthermore:
    !  a) Updates vgs - energy of the ground state (if present)
    !  b) Update vmdqt - eneriges of excited states (if present)
    !  c) Update cmdqt - expansion coefficients (if present)
    !
    !********************************************************************
    !
    subroutine do_sqm_davidson_update(sim,printTdipole,cmdqt,vmdqt,vgs,rx,ry,rz,r,statelimit)
        use naesmd_module, only : realp_t
        
        implicit none
        type(simulation_t),pointer::sim
        integer printTdipole
        _REAL_,intent(inout),optional::cmdqt(:,:),vmdqt(:),vgs
        _REAL_,intent(in),optional::rx(:),ry(:),rz(:)
        type(realp_t),intent(in),optional::r(3)
        integer i
        integer,optional::statelimit

          
        sim%dav%mdflag=2
          
        call do_sqm_and_davidson(sim,printTdipole,rx,ry,rz,r,statelimit)
        call dav2naesmd_Omega(sim) ! energy conversion
        if(present(vgs)) vgs=sim%naesmd%E0
          
        if(present(vmdqt)) then
            do i=1,sim%excN
                vmdqt(i)=(sim%naesmd%Omega(i)+sim%naesmd%E0)
            end do
        end if
          
        if(present(cmdqt)) then
            call dav2cmdqt(sim,cmdqt)
        end if
        return
    end subroutine do_sqm_davidson_update
    !
    !********************************************************************
    !
    subroutine nacR_analytic_wrap(sim, ihop, icheck, dij)
        type(simulation_t), pointer :: sim
        integer, intent(in)         :: ihop,icheck
        _REAL_,intent(inout), optional ::  dij(:)
        integer :: mn
          
        call qmmm2coords_r(sim)
        call nacR_analytic(sim%qparams,sim%qnml,sim%qmpi,sim%rij,sim%qm2,sim%dav,sim%qmmm, sim%coords,ihop,icheck)
        if(present(dij))  then
            mn = min(size(dij), size(sim%dav%dij))
            dij(:mn) = -sim%dav%dij(:mn)
        endif
    end subroutine
    !
    !********************************************************************
    !
    !  Initialization of sqm
    !
    !********************************************************************
    !
    subroutine init_sqm(sim,fdes_in,fdes_out)
        use qmmm_module,only: deallocate_qmmm
          
        use constants,only:KCAL_TO_EV,EV_TO_KCAL
        use file_io_dat,only:MAX_FN_LEN
        use qm2_davidson_module ! CML 7/11/12
 
        implicit none
          
        type(simulation_t),pointer::sim
          
        integer,intent(in)::fdes_in,fdes_out
        character(len=8) atnam(sim%Na)
        integer ier, xmin_iter
        integer ntpr
        _REAL_ excharge(sim%qmmm%qm_mm_pairs*4)
        integer chgatnum(sim%Na)
        integer ncharge
        integer :: igb, maxcyc, do_nm
        _REAL_ :: deltaX
        _REAL_  :: grms_tol
        logical :: master=.true.

        ! ==== Initialise first_call flags for QMMM ====
        sim%qmmm%qm_mm_first_call = .true.
        sim%qmmm%fock_first_call = .true.
        sim%qmmm%fock2_2atm_first_call = .true.
        sim%qmmm%qm2_allocate_e_repul_first_call = .true.
        sim%qmmm%qm2_calc_rij_eqns_first_call = .true.
        sim%qmmm%qm2_scf_first_call = .true.
        sim%qmmm%zero_link_charges_first_call = .true.
        sim%qmmm%adj_mm_link_pair_crd_first_call = .true.
        ncharge=0;
        igb = 0
         
        if(sim%Nsim.eq.1) then 
             sim%dav%modes_b='modes.b'
             sim%dav%modes_unit=10
             sim%dav%ee_b='ee.b'
             sim%dav%ee_unit=15
             sim%dav%normalmodesao='normalmodes-ao.out'
             sim%dav%normmodesao_unit=76
             sim%dav%normalmodesmo='normalmodes-mo.out'
             sim%dav%normmodesmo_unit=77
             sim%dav%normalmodescf='normalmodes-cfit.out'
             sim%dav%normmodescf_unit=78
             sim%dav%muab_unit=13
             sim%dav%muab_out='muab.out'
             sim%dav%ceo_unit=14
             sim%dav%ceo_out='ceo.out'
        else 
             write (sim%dav%modes_b, "(A16,I4.4,A2)") "modes_", sim%id , ".b"
             sim%dav%modes_unit=10*10000+sim%id
             write (sim%dav%ee_b, "(A3,I4.4,A2)") "ee_", sim%id , ".b"
             sim%dav%ee_unit=15*10000+sim%id
             write (sim%dav%normalmodesao, "(A15,I4.4,A4)") "normalmodes-ao_", sim%id , ".out"
             sim%dav%normmodesao_unit=76*10000+sim%id
             write (sim%dav%normalmodesao, "(A15,I4.4,A4)") "normalmodes-mo_", sim%id , ".out"
             sim%dav%normmodesmo_unit=77*10000+sim%id
             write (sim%dav%normalmodesao, "(A17,I4.4,A4)") "normalmodes-cfit_", sim%id , ".out"
             sim%dav%normmodescf_unit=78*10000+sim%id
             write (sim%dav%muab_out, "(A5,I4.4,A4)") "muab_", sim%id , ".out"
             sim%dav%muab_unit=13*10000+sim%id
             write (sim%dav%ceo_out, "(A4,I4.4,A4)") "ceo_", sim%id , ".out"
             sim%dav%ceo_unit=14*10000+sim%id
        endif
 
        call sqm_read_and_alloc(sim%qnml,sim%scratch,sim%div,sim%opnq,sim%vsolv,sim%xlbomd,sim%cosmo,&
            sim%qm2,sim%dav,sim%qmmm,fdes_in, fdes_out, &
            sim%Na,igb,atnam,sim%naesmd%atomtype,maxcyc, &
            grms_tol,ntpr, ncharge,excharge,chgatnum, &
            sim%excN,sim%dav%struct_opt_state,sim%dav%idav,sim%dav%dav_guess, &
            sim%dav%ftol,sim%dav%ftol0,sim%dav%ftol1,sim%dav%icount_M,sim%naesmd%nstep,do_nm,deltaX)
          
        call qm_assign_atom_types(sim%qmmm)
          
        ! Set default QMMM MPI parameters - for single cpu operation.
        ! These will get overwritten by qmmm_mpi_setup if MPI is on.
        ! qmmm_mpi%master = master
        sim%qmpi%commqmmm_master = master
        sim%qmpi%numthreads = 1
        sim%qmpi%mytaskid = 0
        sim%qmpi%natom_start = 1
        sim%qmpi%natom_end = sim%naesmd%natom
        sim%qmpi%nquant_nlink_start = 1
        sim%qmpi%nquant_nlink_end = sim%qmmm%nquant_nlink
        call allocate_qmgb(sim%gb,sim%qmmm%nquant_nlink)
          
        allocate( sim%qmmm%dxyzqm(3, sim%qmmm%nquant_nlink), stat = ier )
        REQUIRE(ier == 0)
          
        allocate ( sim%qm2%scf_mchg(sim%qmmm%nquant_nlink), stat = ier )
        REQUIRE(ier == 0)
          
        allocate ( sim%qmmm%qm_coords(3,sim%qmmm%nquant_nlink), stat=ier )
        REQUIRE(ier == 0)
          
        !Check if we are doing Geometry Optimization or Dynamics or NM and act accordingly
        if (maxcyc.gt.0.and.do_nm.gt.0) then
            write(6,*) "You must run geometry optimization (maxcyc > 0) or calculate nuclear &
                normal modes (do_nm > 0). Running both simultaneously is not possible."    
            STOP 0;
        elseif (sim%naesmd%nstep.gt.0.and.do_nm.gt.0) then
            write(6,*) "You must run dynamics (n_class_steps > 0) or calculate nuclear normal &
                modes (do_nm > 0). Running both simultaneously is not possible."
            STOP 0;
        elseif (do_nm.gt.0) then
            sim%dav%minimization = .FALSE.
            sim%Ncharge = ncharge
            sim%dav%verbosity=0 !don't print output from scf calculations
            call hessian_calc(sim,deltaX)
        endif
        if (maxcyc < 1) then
            sim%dav%minimization = .FALSE.
            sim%Ncharge  = ncharge
            !sim%qm2      => qm2_struct
            call naesmd2qmmm_r(sim)
            !sim%dav%Mx = sim%excN
            !sim%dav%mdflag=2
            if ((sim%cosmo%solvent_model.eq.4).or.(sim%cosmo%solvent_model.eq.5)) then !solvent model that loops over ground sim%naesmd%state
                call calc_cosmo_4(sim)
            else
                call do_sqm_davidson_update(sim,sim%naesmd%printTdipole,sim%naesmd%cmdqt,sim%naesmd%vmdqt,sim%naesmd%vgs)
            endif
        else if (sim%naesmd%nstep<1) then ! sim%naesmd%nstep - number of classical steps for dynamics
          
            sim%dav%minimization = .TRUE.
            sim%Ncharge  = ncharge
            !sim%qm2 => qm2_struct
          
            if (sim%qnml%verbosity<5) then
                sim%dav%verbosity=0 !don't print output from scf calculations
            endif
            call naesmd2qmmm_r(sim)
            call xmin(sim%naesmd%natom, sim%qmmm%qm_coords, xmin_iter, maxcyc, grms_tol, ntpr, sim, sim%qnml%verbosity)
        else
            write(6,*) "You must run dynamics (n_class_steps > 0) or geometry optimization &
                (maxcyc > 0). Running both simultaneously is not possible."            
            STOP 0;
        end if
        return
    end subroutine init_sqm

          
          
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ Driver routine for XMIN minimization.
    !-----------------------------------------------------------------------
#include "xmin.h" 
    !includes return_flag definitions for select_case statement shared with XminC
    subroutine xmin(natom, x, xmin_iter, maxiter, grms_tol, ntpr, sim_pass, verbosity)
          
        use qm2_davidson_module
        use constants, only : zero, AU_TO_EV
        use iso_c_binding
    !    use Cosmo_C, only : solvent_model

        implicit none
   
        ! ------ External functions -----------------
        !_REAL_   xminC
        !external xminC
          
        !Fortran 90 interface with xminc by JAB
        interface
            real(c_double) function xminc( xyz_min, xmin_method_code, maxiter, grms_tol, &
                natom, lbfgs_memory_depth, mvpm_code, &
                x, escf, fg, grms, xmin_iter, xmin_time, &
                xmin_verbosity, ls_method, ls_maxiter, ls_iter, ls_maxatmov, &
                beta_armijo, c_armijo, mu_armijo, ftol_wolfe, gtol_wolfe,  &
                return_flag, status_flag ) bind(c)
                use, intrinsic                  ::iso_c_binding
                integer(c_int)                  ::xyz_min
                integer(c_int)                  ::xmin_method_code
                integer(c_int)                  ::maxiter
                real(c_double)                  ::grms_tol
                integer(c_int)                  ::natom
                integer(c_int)                  ::lbfgs_memory_depth
                integer(c_int)                  ::mvpm_code
                real(c_double)                  ::x(natom*3)
                real(c_double)                  ::escf
                real(c_double)                  ::fg(natom*3)
                real(c_double)                  ::grms
                integer(c_int)                  ::xmin_iter
                real(c_double)                  ::xmin_time
                integer(c_int)                  ::xmin_verbosity
                integer(c_int)                  ::ls_method
                integer(c_int)                  ::ls_maxiter
                integer(c_int)                  ::ls_iter
                real(c_double)                  ::ls_maxatmov
                real(c_double)                  ::beta_armijo
                real(c_double)                  ::c_armijo
                real(c_double)                  ::mu_armijo
                real(c_double)                  ::ftol_wolfe
                real(c_double)                  ::gtol_wolfe
                integer(c_int)                  ::return_flag
                integer(c_int)                  ::status_flag
            end function xminc
        end interface
        type(simulation_t),pointer::sim_pass
        integer, intent(in) :: natom, maxiter, ntpr
        _REAL_,  intent(inout) :: x(natom*3) !coordinates
        integer, intent(inout) :: xmin_iter
        integer, intent(in) :: verbosity
        _REAL_,  intent(in)  :: grms_tol
          
        ! ------ local variables --------------------
        _REAL_  ::fg(natom*3) !derivatives (should it be forces?)
        _REAL_  ::energy
        _REAL_  :: grms            = ZERO
        logical :: is_error
        logical :: is_xmin_done
        _REAL_  :: minimum_energy
        integer :: n_force_calls
        integer :: return_flag
        integer :: status_flag
        _REAL_  :: xmin_time
          
        ! The depth of the LBFGS memory for XMIN's LBFGS minimization or TNCG
        ! preconditioning.
        ! The value 0 turns off preconditioning in TNCG minimization.
        integer  :: lbfgs_memory_depth = 3
          
        ! XMIN's finite difference Hv matrix-vector product method.
        ! Renaming of Kolossvary's numdiff.
        integer  :: mvpm_code = 1
          
        ! XMIN minimization method.
        integer,          parameter :: XMIN_METHOD_PRCG_CODE  = 1
        integer,          parameter :: XMIN_METHOD_LBFGS_CODE = 2
        integer,          parameter :: XMIN_METHOD_TNCG_CODE  = 3
        integer       :: xmin_method_code = 3
          
        ! Verbosity of the internal status output from the XMIN package:
        ! 0 = none, 1 = minimization details, 2 = minimization and
        ! line search details plus CG details in TNCG.
        integer,      parameter :: MAXIMUM_XMIN_VERBOSITY = 2
        integer,      parameter :: MINIMUM_XMIN_VERBOSITY = 0
        integer                 :: xmin_verbosity = 0
        integer i
          
        integer ::  ls_method, ls_maxiter, ls_iter, xyz_min
        _REAL_  :: ls_maxatmov, beta_armijo, c_armijo, mu_armijo, ftol_wolfe, &
            gtol_wolfe
          
        logical :: first
        ! ------ External Functions -----------------
        _REAL_ ddot
100     format(I5,'     ',3F12.6) !for printing final coordinates
        first = .true.
        status_flag = 0
        xyz_min = 1
        ls_method = 2
        ls_maxiter = 20
        ls_maxatmov = 0.2
        beta_armijo = 0.5
        c_armijo = 0.4
        mu_armijo = 1.0
        ftol_wolfe = 0.0001
        gtol_wolfe = 0.9
          
        if (verbosity > 1) open(500,file='coords_opt.xyz',status='replace')

        ! keep xminC() from thinking that atoms are frozen:
        fg(1:3*natom) = 1.d0
          
        n_force_calls = 0
        is_error = .false.
        is_xmin_done = .false.
        if(sim_pass%qmmm%state_of_interest>0) then
            energy=(sim_pass%naesmd%Omega(sim_pass%qmmm%state_of_interest)+sim_pass%naesmd%E0)*AU_TO_EV
        else
            energy=sim_pass%naesmd%E0*AU_TO_EV
        endif
        do while ( .not. is_xmin_done .and. .not. is_error )
            minimum_energy = xminc( xyz_min, xmin_method_code, maxiter, grms_tol, &
                natom, lbfgs_memory_depth, mvpm_code, &
                x, energy, fg, grms, xmin_iter, xmin_time, &
                xmin_verbosity, ls_method, ls_maxiter, ls_iter, ls_maxatmov, &
                beta_armijo, c_armijo, mu_armijo, ftol_wolfe, gtol_wolfe,  &
                return_flag, status_flag )
            select case ( return_flag )
                case ( DONE )
                    ! Finished minimization.
                    is_xmin_done = .true.
                    is_error = status_flag < 0
                    if(maxiter==xmin_iter) then
                        write(6,*) 'Maximum number of iterations reached without convergence. Calculation failed...'
                        stop
                    else
                        write(6,'(a)') '  ... geometry converged!'
                        write(6,*)
                        write(6,*) 'Final Structure'
                        !call qm_print_coords(sim_pass%qnml,0,.true.)
                        do i=1,natom !NAESMD format
                            write(6,100) sim_pass%qmmm%iqm_atomic_numbers(i),x(3*(i-1)+1),x(3*(i-1)+2),x(3*(i-1)+3)
                        end do
          
                        if ( sim_pass%qnml%printbondorders ) then
                            write(6,*) ''
                            write(6,*) 'Bond Orders'
                            call qm2_print_bondorders(sim_pass%qparams,sim_pass%qm2,sim_pass%qmmm)
                        end if
                    endif
                case ( CALCENRG, CALCGRAD, CALCBOTH )
                    ! Normal Amber control of NB list updates.
                    is_error = status_flag < 0
                    if ( .not. is_error ) then
                        call qmmm2naesmd_r(sim_pass)
                        if((sim_pass%cosmo%solvent_model.eq.4).or.(sim_pass%cosmo%solvent_model.eq.5)) then
                            call calc_cosmo_4(sim_pass)
                        else
                            call do_sqm_davidson_update(sim_pass,0)
                        endif
                        call deriv(sim_pass,sim_pass%qmmm%state_of_interest)
                        fg(1:3*natom)=-sim_pass%deriv_forces(1:3*natom)
                        if(sim_pass%qmmm%state_of_interest>0) then
                            energy=(sim_pass%naesmd%Omega(sim_pass%qmmm%state_of_interest)+sim_pass%naesmd%E0)*AU_TO_EV
                        else
                            energy=sim_pass%naesmd%E0*AU_TO_EV
                        endif
                        n_force_calls = n_force_calls + 1
                    end if
                case ( CALCENRG_NEWNBL, CALCGRAD_NEWNBL, CALCBOTH_NEWNBL )
                    ! Coerce a NB list update.
                    is_error = status_flag < 0
                    if ( .not. is_error ) then
                        if( n_force_calls>0 .and. mod(xmin_iter,ntpr)==0 ) then
                            grms = sqrt(ddot(3*natom,fg,1,fg,1)/dble(3*natom))
                            if (first) then
                                write(6,'(a)') '       iter       energy            rms gradient'
                                write(6,'(a)') '       ----    --------------     -----------------'
                                first = .false.
                            end if
                            write(6,'(a,i5,f16.6,a,f16.6,a)') 'xmin ', xmin_iter, energy,' eV', grms, ' eV/A'
                        end if

                        if (verbosity > 1) then
                        !write intermediate steps
                        write(500,*) natom
                        write(500,*) 'optimization step ', xmin_iter
                          do i=1,natom !NAESMD format
                            write(500,100) sim_pass%qmmm%iqm_atomic_numbers(i),x(3*(i-1)+1),x(3*(i-1)+2),x(3*(i-1)+3)
                          end do
                        endif

                        call qmmm2naesmd_r(sim_pass);
                        if((sim_pass%cosmo%solvent_model.eq.4).or.(sim_pass%cosmo%solvent_model.eq.5)) then
                            call calc_cosmo_4(sim_pass)
                        else
                            call do_sqm_davidson_update(sim_pass,0)
                        endif
                        call deriv(sim_pass,sim_pass%qmmm%state_of_interest)
                        if(sim_pass%qmmm%state_of_interest>0) then
                            energy=(sim_pass%naesmd%Omega(sim_pass%qmmm%state_of_interest)+sim_pass%naesmd%E0)*AU_TO_EV
                        else
                            energy=sim_pass%naesmd%E0*AU_TO_EV
                        endif
                        fg(1:3*natom)=-sim_pass%deriv_forces(1:3*natom)
                        n_force_calls = n_force_calls + 1
                    end if
                case ( CALCENRG_OLDNBL, CALCGRAD_OLDNBL, CALCBOTH_OLDNBL )
                    ! Prevent a NB list update.
                    is_error = status_flag < 0
                    if ( .not. is_error ) then
                        call qmmm2naesmd_r(sim_pass);
                        call do_sqm_davidson_update(sim_pass,0);
                        if((sim_pass%cosmo%solvent_model.eq.4).or.(sim_pass%cosmo%solvent_model.eq.5)) then
                            call calc_cosmo_4(sim_pass)
                        else
                            call do_sqm_davidson_update(sim_pass,0)
                        endif
                        call deriv(sim_pass,sim_pass%qmmm%state_of_interest)
                        if(sim_pass%qmmm%state_of_interest>0) then
                            energy=(sim_pass%naesmd%Omega(sim_pass%qmmm%state_of_interest)+sim_pass%naesmd%E0)*AU_TO_EV
                        else
                            energy=sim_pass%naesmd%E0*AU_TO_EV
                        endif
                        fg(1:3*natom)=-sim_pass%deriv_forces(1:3*natom);
                        n_force_calls = n_force_calls + 1
                    end if
                case default
                    ! error from XMIN or the return_flag is corrupted.
                    is_error = status_flag < 0
                    ASSERT( is_error )
            end select
        end do
          
        if (verbosity > 1) close(500)
        if ( is_error ) then
            write(6,'(a,i4)') '  XMIN ERROR: Status is ', status_flag
            call mexit(6,1)
        end if
          
        return
    end subroutine xmin
          
end module communism
