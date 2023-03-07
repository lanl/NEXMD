#include "assert.fh"
#include "dprec.fh"

#ifdef MPI
   include 'mpif.h'
#endif

#ifdef MPI
   include 'mpif.h'
#endif
    
subroutine sqm_energy(qmewald, qmmm_opnq, qm2_rij_eqns,qmmm_mpi,qm_gb, qmmm_scratch, qmmm_nml, qm2_params, &
    xlbomd_struct,cosmo_c_struct, qm2_struct,qm2ds,qmmm_struct, natom,coords,escf,born_radii,one_born_radii, &
    intdiel,extdiel,Arad,scf_mchg, &
    time_sqm_took,time_davidson_took,printTdipole,outfile_28,tfemto)
    !
    !     Argument list variables:
    !
    !     coords(natom*3)                 - Cartesian coordinates for all atoms.
    !                                       Amber array
    !     natom                           - Total number of REAL atoms.
    !     qmmm_struct%nquant              - Number of REAL quantum atoms as specified in mdin.
    !     iqmatoms(qmmm_struct%nquant)
    !                                     - Atom numbers for quantum atoms link atoms given values of -1
    !     qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant) - Atomic numbers for qm atoms.
    !     qmmm_struct%nlink               - Number of link atoms.
    !     escf                            - Heat of formation from QM.
    !     qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink)
    !                                     - Cartesian coordinates of quantum atoms.
    !                                       (Extracted from coords by qm_extract_coords)
    !
    !    Locally defined arrays:
    !    born_radii(1->natom)      - Effective GB radii - only used when doing qm with gb (and qm_gb==2)
    !                                Calculated via an initial call to egb.
    !    one_born_radii(1->natom)  - 1.0d0/born_radii(i)
    !    scf_mchg                  - nquant long, gets filled with the mulliken charges during scf.

    use qmmm_module,only:qm2_structure,qmmm_scratch_structure, qm2_rij_eqns_structure, &
                         qm_gb_structure, qmmm_mpi_structure, qmmm_opnq_structure, &
                         qm_ewald_structure
    use constants,only:EV_TO_KCAL,KCAL_TO_EV,zero,one,alpb_alpha
    use qm2_davidson_module
    use cosmo_C,only:cosmo_C_structure !,cosmo_c_struct%EF,cosmo_c_struct%nps,cosmo_c_struct%tri_2D,cosmo_c_struct%potential_type,cosmo_c_struct%solvent_model,cosmo_c_struct%onsagE;
    use xlbomd_module,only:init_xlbomd, xlbomd_structure
    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_nml_module   , only : qmmm_nml_type
    use qm2_params_module , only : qm2_params_type
    
    implicit none

    ! Passed in

    integer,intent(in)::natom
   type(qm_ewald_structure),intent(inout) :: qmewald
   type(qmmm_opnq_structure),intent(inout) :: qmmm_opnq
   type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
    type(qm2_rij_eqns_structure),intent(inout)  :: qm2_rij_eqns
    type(qmmm_mpi_structure),intent(inout)      :: qmmm_mpi
    type(qm_gb_structure), intent(inout)        :: qm_gb
    type(qmmm_nml_type)   , intent(inout) :: qmmm_nml
    type(qm2_params_type), intent(inout) :: qm2_params
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(xlbomd_structure),intent(inout) :: xlbomd_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    _REAL_,intent(inout)::coords(natom*3) !Amber array - adjusted for link atoms
    _REAL_,intent(out)::escf
    _REAL_,intent(in)::born_radii(natom),one_born_radii(natom)
    _REAL_,intent(in)::intdiel,extdiel,Arad
    _REAL_,intent(inout)::scf_mchg(qmmm_struct%nquant_nlink)
    _REAL_,intent(inout)::time_sqm_took,time_davidson_took
    _REAL_ ,allocatable, dimension(:,:) :: density_matrix_unpacked
    integer, intent(in) :: printTdipole,outfile_28
    _REAL_, intent(in) :: tfemto

    ! Locals
    _REAL_,dimension(2,3)::bxbnd
    _REAL_::mulliken_charge,total_mulliken_charge,total_energy
    _REAL_::alpb_beta
    _REAL_::scaled_mm_charges(2)

    integer::ier=0
    integer i,j,offset,qm_no,i3,T2DS;

    ! Locals for link atoms
    _REAL_::forcemod(3)
    integer::lnk_no,mm_no
    _REAL_ ctest
    integer quir_ev,quir_cmdqt,l
    _REAL_ t_start,t_finish ! to monitor execution time

    !=============================================================================
    !                   START OF QMMM SETUP: allocate list memory
    !=============================================================================
    if (qmmm_struct%qm_mm_first_call) then
    end if ! ---- first call end if ----------

    !  If this is the first call to the routine, do some initial allocation
    !  that has not been done elsewhere.
    if (qmmm_struct%qm_mm_first_call) then

        !  allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant_nlink), stat=ier )
        !  Stores the REAL and link atom qm coordinates
        !  REQUIRE(ier == 0)

        !Allocation for QM_GB (qmgb==2)
        if (qmmm_nml%qmgb==2) then
            !Calculate dielectric factor
            if (qm_gb%alpb_on) then
                alpb_beta=alpb_alpha*(intdiel/extdiel)
                qm_gb%intdieli=one/(intdiel*(one+alpb_beta))
                qm_gb%extdieli=one/(extdiel*(one+alpb_beta))
                qm_gb%one_Arad_beta=alpb_beta/Arad
            else
                qm_gb%intdieli=1.d0/intdiel
                qm_gb%extdieli=1.d0/extdiel
            end if

            qm_gb%mmcut2=999.d0
        end if
    end if ! ---- first call end if ----------

    !=============================================================================
    !                   START OF REST OF QMMM SETUP
    !=============================================================================
    if(qmmm_struct%qm_mm_first_call) then
        if (qmmm_mpi%commqmmm_master) then
            write(6,'(/80(1H-)/''  QM CALCULATION INFO'',/80(1H-))')
        end if
        call qm2_load_params_and_allocate(qm2_params, qmmm_nml, qmmm_mpi, qmmm_opnq, qmmm_scratch, &
                                          xlbomd_struct,qm2_struct,qmmm_struct) !Load the parameters
           ! Also does a lot of memory allocation and pre-calculates all
           ! the STO-6G orbital expansions.

        if(qmmm_mpi%commqmmm_master) then

            call qm_print_coords(qmmm_nml,qmmm_struct,0,.true.)

            !Finally print the result header that was skipped in sander.
            write(6,'(/80(1H-)/''  RESULTS'',/80(1H-)/)')
        end if

        if(qm2ds%Mx>0) then
            call allocate_davidson(qmmm_scratch,qmmm_nml,qm2_struct, qm2ds, qmmm_struct) ! Davidson allocation
            if(qm2ds%dav_guess.gt.1) call init_xlbomd(xlbomd_struct,qm2ds%Nb**2*qm2ds%Mx)
        endif

        !cosmo_c_struct%ceps-Dielectric Permittivity from COSMO module
        if((cosmo_c_struct%solvent_model.gt.0).or.(cosmo_c_struct%EF.gt.0)) then !
            call cosini(qm2_params, qmmm_nml, cosmo_c_struct,qm2_struct,qmmm_struct)
        end if

    end if !if (qmmm_struct%qm_mm_first_call)

    !======================END OF QMMM SETUP ======================================

    ! Calculate RIJ and many related equations here. Necessary memory allocation
    ! is done inside the routine.
    ! Parallel

    call qm2_calc_rij_and_eqns(qm2_params,qmmm_nml,qmmm_mpi,qm2_rij_eqns, &
        qmmm_struct, qmmm_struct%qm_coords,qmmm_struct%nquant_nlink, &
        qmmm_struct%qm_xcrd,natom,qmmm_struct%qm_mm_pairs)
       !and store them in memory to save time later.

    !============================
    ! Constructing/updating COSMO cavity
    !============================
    !cosmo_c_struct%ceps-Dielectric Permittivity from COSMO module
    if((cosmo_c_struct%solvent_model.gt.0).and.(cosmo_c_struct%potential_type.eq.3)) then ! non-default non-vacuum permittivity
        call coscav(qm2_params, qmmm_nml, cosmo_c_struct,qmmm_struct) ! constructing COSMO cavity
        call mkbmat(qm2_params,qmmm_nml,cosmo_c_struct,qmmm_struct) ! constructing B matrix
    end if

    !============================
    !  Calculate SCF Energy (ground state)
    !============================
    call cpu_time(t_start)
    call qm2_energy(qm2_params,qmmm_nml,qmewald,qm2_rij_eqns, qm_gb, qmmm_mpi, qmmm_opnq, qmmm_scratch, &
        xlbomd_struct,cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct, escf,scf_mchg,natom,born_radii,one_born_radii, &
        coords,scaled_mm_charges)
    call cpu_time(t_finish)
    ! Incrementing the cumulative sqm time
    time_sqm_took=time_sqm_took+t_finish-t_start

    !Total energy of the ground state
    qm2ds%Eground=qmmm_struct%elec_eng+qmmm_struct%enuclr_qmmm &
        +qmmm_struct%enuclr_qmqm

    !=============================
    !  Print Mulliken Charges for the ground state
    !=============================

    if(qmmm_nml%printcharges.and.qmmm_mpi%commqmmm_master) then
        allocate(density_matrix_unpacked(qm2_struct%norbs,qm2_struct%norbs))
        call unpacking(qm2_struct%norbs,qm2_struct%den_matrix,density_matrix_unpacked,'s')
        do i=1,qmmm_struct%nquant_nlink
               call qm2_calc_mulliken(qm2_params,qm2_struct,i,scf_mchg(i),density_matrix_unpacked);
        end do
        deallocate(density_matrix_unpacked);
        write (6,*)
        write (6,'("QMMM: Mulliken Charges")')
        call qm2_print_charges(qmmm_nml,qmmm_mpi, qmmm_struct, 0,1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
            scf_mchg,qmmm_struct%iqm_atomic_numbers)
        write(6,*)
    end if

    ! preserving the quirality of the HF orbitals
    ! important for non-adiabatic dynamics
    if(.not.allocated(qm2_struct%ev_old)) then ! first call
        allocate(qm2_struct%ev_old(qm2_struct%norbs,qm2_struct%norbs))
        qm2_struct%ev_old(:,:)=0.d0
    end if
      
    quir_ev=0
    do j=1,qm2_struct%norbs
        ctest=0.d0

        do i=1,qm2_struct%norbs
            ctest=ctest+qm2_struct%ev_old(i,j)*qm2_struct%eigen_vectors(i,j)
        end do

        if(ctest<0.d0) then
            do l=1,qm2ds%Nb
                qm2_struct%eigen_vectors(l,j)=-qm2_struct%eigen_vectors(l,j)
                quir_ev=quir_ev+1
            end do
        end if
    end do

    if((quir_ev>0).and.(qm2ds%verbosity>0)) then ! quirality was restored at least once above
        print*,' Restoring quirality in HF orbitals'
    end if
    flush(6)

    ! updating
    qm2_struct%ev_old(:,:)=qm2_struct%eigen_vectors(:,:)
      
    !===================================
    !  Calculate Excited State Energies (Davidson)
    !===================================

    if(qm2ds%Mx>0) then
        call cpu_time(t_start)
        call dav_wrap(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct,printTdipole,outfile_28,tfemto)
        call cpu_time(t_finish)
        ! Incrementing cumulative Davidson execution time
        time_davidson_took=time_davidson_took+t_finish-t_start
    end if

    !=============================
    !   Print Dipole Charges
    !=============================

    qmmm_struct%qm_mm_first_call=.false.
    if ((qmmm_nml%printdipole>0).or.(qmmm_nml%printcharges)) then
        if((qmmm_nml%printdipole<3).and.(qm2ds%Mx>0)) then
            call qm2_calc_molecular_dipole_in_excited_state(qm2_params,qmmm_nml,qmmm_mpi,&
                 cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct)
            
        else if (qmmm_nml%printdipole>2) then
            call qm2_calc_dipole(qmmm_nml,qm2_params,qm2_struct,qm2ds, qmmm_struct,coords)
        end if
    end if

    ! Print some extra informatiom about energy contributions
    ! (This is really only required in sander since we print the energies anyways
    !  but kept here for historical reasons)

    if(qmmm_mpi%commqmmm_master) then
        call qm2_print_energy(qmmm_scratch,cosmo_c_struct,qm2_struct, qmmm_nml%verbosity,qmmm_nml%qmtheory,escf,qmmm_struct)
    end if
    return
end subroutine sqm_energy
!
!********************************************************************
!======================END OF QM_MM ======================================

!-------------------------------------------------
!     --- FLOAT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
    implicit none
    _REAL_ param,lo,hi
    character(len=*)string

    if ( param < lo .or. param > hi )then
        write(6,59)
        write(6,60)string,param
        write(6,61)
        write(6,62)lo,hi
        write(6,63)
        call mexit(6,1)
    end if
59  format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
60  format(1x,'parameter ',a,' has value ',e12.5)
61  format(1x,'This is outside the legal range')
62  format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
63  format(1x,'Check ew_legal.h')
    return
end subroutine float_legal_range

!-------------------------------------------------
!     --- INT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
    implicit none
    integer param,lo,hi
    character(len=*)string

    if ( param < lo .or. param > hi )then
        write(6,59)
        write(6,60)string,param
        write(6,61)
        write(6,62)lo,hi
        write(6,63)
        call mexit(6,1)
    end if
59  format(/,1x,'PARAMETER RANGE CHECKING: ')
60  format(1x,'parameter ',a,' has value ',i8)
61  format(1x,'This is outside the legal range')
62  format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
63  format(1x,'The limits may be adjustable; search in the .h files ')
    return
end subroutine int_legal_range

!-------------------------------------------------
!     --- SANDER_BOMB ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
    implicit none
    character(len=*) routine,string1,string2

    write(6, '(1x,2a)') &
        'SANDER BOMB in subroutine ', routine
    write(6, '(1x,a)') string1
    write(6, '(1x,a)') string2
    call mexit(6,1)
end subroutine sander_bomb
!-------------------------------------------------

subroutine getsqmx(natom,x,atnam,atnum,ncharge,excharge,chgnam,chgatnum)
       
    !     --- reads initial coords,

    implicit none
    _REAL_ x(*)
    integer i,j,i3,lun
    integer natom,ier,atnum(*)
    character(len=8) atnam(*)
    character(len=80) line
    ! test-local
    _REAL_ excharge(*)
    integer chgatnum(*)
    character(len=8) chgnam(*)
    integer ncharge
    integer ia, ic, ihead, iend
    logical mdin_qm_atom
    logical mdin_external_charge

    lun = 5
    mdin_qm_atom = .false.
    mdin_external_charge = .false.
    ncharge = 0

    ! check header names
    ihead=0
    iend=0
    do i=1,999
        read(lun,'(a)',end=10) line
        if (line(1:1) == "#") then
            if (line(1:80) == "#EXCHARGES") then
                mdin_external_charge = .true.
                ihead = ihead + 1
            else if (line(1:80) == "#END") then
                iend = iend + 1
            else
                write(0,*) 'Unrecognized header name'
                write(0,*) line(1:80)
                call mexit(6,1)
            end if
        end if
    end do

10  if (iend < ihead) then
        write(0,*) 'Missing "#END" termination sign, exit program'
        call mexit(6,1)
    end if

    rewind(lun)

    !  skip over the &qmmm namelist at the beginning:
    ! CML The following commented lines are the original SQM. This assumes the qmmm namelist
    ! CML is 20 lines or less which means when we add parameters, the program complains.
    ! CML The procedure underneath runs the search until we find the terminating string (' /')
    ! CML is found or the end of the file is found.


    !---------MODIFIED BY CML----------
    do
        read(5,'(a)',IOSTAT=i) line
        ! if we find the end of nlist or the EOF, then exit
        if( line(1:2) == " /" .OR. i == -1 ) exit
    end do

    if (i == -1) then   ! if you reach EOF, then quit
        write(0,*) 'Error in finding end of qmmm namelist'
        call mexit(6,1)
    end if
    !---------END MODIFICATION---------

       
    ! reading QM atoms
11  i3=0
    ia=0
    do i=1,999
        read(lun,'(a)',end=12) line
        if (line(1:80) /= "") then
            if (line(1:80) /= "#EXCHARGES") then
                ia = ia + 1
                read(line,*,err=15) atnum(ia),atnam(ia),x(i3+1),x(i3+2),x(i3+3)
                i3 = i3 + 3
            else
                go to 12
            end if
        end if
    end do

12  natom = ia

    ! reading external charges
    if (mdin_external_charge) then
        i3=0
        ic=0
        do i=1,999
            read(lun,'(a)',end=14) line
            if (line(1:80) /= "") then
                if (line(1:80) /= "#END") then
                    ic = ic + 1
                    read(line,*,err=16) chgatnum(ic), chgnam(ic),excharge(i3+1:i3+4)
                    i3 = i3 + 4
                else
                    go to 13
                end if
            end if
        end do
13      ncharge = ic

        write(6,'(/80(1H-)/''  EXTERNAL CHARGES FOUND IN INPUT'',/80(1H-))')
        write(6,'(2x,"QMMM: External Charge Info")')
        write(6,'(2x,"QMMM:",1x,"ATOMIC",3x,"NAME",8x,"X",9x,"Y",9X,"Z",8X,"CHARGE")')

        i3=0
        do i=1,ncharge
            write(6,'(2x,"QMMM:",3x,i2,6x,a6,4f10.4)') chgatnum(i), chgnam(i), excharge(i3+1:i3+4)
            i3=i3+4
        end do
    end if

    return

14  write(0,*) 'The termination sign "#END" is missing'
    call mexit(6,1)

15  write(0,*) 'Error in reading QM atoms'
    call mexit(6,1)

16  write(0,*) 'Error in reading external charges'
    call mexit(6,1)

end subroutine getsqmx


!+++++++++++++++++++++++++++++++++++++++++++++
!This subroutine reads the QMMM namelist
!and also calls allocation routines
!for QMMM based on natom.
!
!Author:
!     Ross Walker (SDSC)
!+++++++++++++++++++++++++++++++++++++++++++++


!
!********************************************************************
!
!  Reads the qmmm namelist and calls the qmmm memory allocation routines
!  modified by CML 7/10/12
!  modified by KGB 7/24/12 (removed all blocks that are #ifndef SQM ) 
!
!********************************************************************
!
subroutine sqm_read_and_alloc(qmmm_nml, qmmm_scratch, qmmm_div, qmmm_opnq, qmmm_vsolv, xlbomd_struct,cosmo_c_struct,qm2_struct,&
    qm2ds,qmmm_struct,fdes_in,fdes_out,natom_inout,igb,atnam, &
    atnum,maxcyc,grms_tol,ntpr,ncharge_in, &
    excharge,chgatnum, &
    excNin,struct_opt_state,exst_method,dav_guess, & ! Excited state
    ftol,ftol0,ftol1,dav_maxcyc,nstep,do_nm,deltaX)

    use findmask
    use constants, only : RETIRED_INPUT_OPTION
    use qmmm_module, only : qm2_structure, &
        validate_qm_atoms, qmsort, &
        allocate_qmmm, get_atomic_number, qmmm_div_structure, &
        qmmm_opnq_structure, qmmm_scratch_structure
    use qmmm_vsolv_module, only : read_vsolv_nml, qmmm_vsolv_type
    use qmmm_qmtheorymodule
    use ElementOrbitalIndex, only : numberElements
    use ParameterReader, only : ReadParameterFile
    use cosmo_C, only: cosmo_C_structure
    use xlbomd_module, only: xlbomd_structure 
    use qm2_davidson_module
    use naesmd_constants
    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_nml_module   , only : qmmm_nml_type

    !XL-BOMD parameters
    implicit none
   type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qmmm_opnq_structure),intent(inout) :: qmmm_opnq
   type(qmmm_div_structure),intent(inout) :: qmmm_div
   type(qmmm_vsolv_type),intent(inout) :: qmmm_vsolv
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(xlbomd_structure),intent(inout) :: xlbomd_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct

    !STATIC MEMORY
    integer :: max_quantum_atoms  !Needed in read qmmm namelist since
    ! namelists cannot contain pointers
    parameter(max_quantum_atoms=10000)
    !END STATIC MEMORY

    !Passed in
    integer :: igb        !Value of igb from cntrl namelist
    integer use_pme, ntb
    integer, intent(in) :: fdes_in, fdes_out ! file descriptor / unit  to read/write
    integer, intent(inout) :: natom_inout
    character(len=8), intent(in) :: atnam(*)
    integer, intent(in)  :: atnum(*)
    integer, intent(out) :: maxcyc, ntpr
    _REAL_, intent(out)  :: grms_tol
    _REAL_, intent(in) :: excharge(*)
    integer, intent(in) :: chgatnum(*)
    integer, intent(in) :: ncharge_in
    integer, intent(in) :: nstep
    integer, intent(out) :: do_nm !Flag for calculating nuclear normal modes [0]
    _REAL_, intent(out) :: deltaX !Displacement for Hessian calculation, A [1.0d-4]

    ! Input parameters for excited state calculations
  
    integer, intent(in)::excNin;
    integer :: excN ! number of excited states to calculate
    integer :: printtd ! Flag to print transition densities
    integer, intent(out) :: struct_opt_state ! CML exc st num to optimize structure
    integer, intent(out)::exst_method ! method used to calc excited states 1=CIS,2=RPA
    integer,intent(out)::dav_guess ! use (1) or not (0) guess for davidson and (2) for XL-BOXMD
    integer,intent(out)::dav_maxcyc !Maximum cycles for davidson, negative for fixed number of cycles
    _REAL_,intent(out)::ftol !   Min tolerance (|emin-eold|)
    _REAL_,intent(out)::ftol0!  Acceptance tol.(|emin-eold|)
    _REAL_,intent(out)::ftol1 ! Accept.tol.for residual norm
  


    !local
    _REAL_ ::qmcut ! local copied to qmmm_nml%qmcut - specified cutoff
    ! to use for QM-MM electrostatics.

                          ! Default = same as regular MM cutoff.
    _REAL_ ::lnk_dis     ! Distance from the QM atom to place link atom.
                          !A value of <0.0 means the link atom gets placed
                         ! on the MM link pair atom's coordinates
                          !on every MD step.
    _REAL_ :: scfconv     ! local copied to qmmm_nml%scfconv - Convergence criteria for SCF routine. Default = 1.0D-8.
                          ! Minimum (tightest criteria) = 1.0D-16

    !+TJG 01/26/2010
    _REAL_ :: errconv      ! Convergence criteria for maximum value in the error matrix FP-PF
    integer :: ndiis_matrices ! Maximum number of matrices used in a diis extrapolation
    integer :: ndiis_attempts ! The initial number of diis tokens... the number of iterations
                              ! that diis extrapolations will be attempted before giving up on diis
    !-TJG 01/26/2010

    _REAL_ :: pseudo_diag_criteria !Criteria - maximum change in density matrix between successive SCF iterations
                                   !in which to allow pseudo diagonalisation. Default = 0.05.
    integer :: lnk_atomic_no !Atomic number of link atom
    integer :: lnk_method !controls how QM-MM valence terms are dealt with.
    integer :: qmgb       ! local copied to qmmm_nml%qmgb - flag for type of GB do with QM region
    integer :: qmtheory   ! deprecated flag for level of theory to use for QM region
    integer :: qmcharge   ! local copied to qmmm_nml%qmcharge - value of charge on QM system
    integer :: spin       ! local copied to qmmm_nml%spin - spin state of system
    integer :: i,j        ! temporary counter for local loops
    integer :: ifind
    integer :: qmqmdx     ! local copied to qmmm_nml%qmqm_analyt - 1 = analytical, 2 = numerical QM-QM derivatives in qm2
    integer :: qmqmdx_exc     ! local copied to qmmm_nml%qmqm_exc_analyt - 1 = analytical, 2 = numerical QM-QM derivatives in excited state procedures CML 7/10/12
    integer :: verbosity  ! local copied to qmmm_nml%verbosity - Controls amount of info about QM part of calc that is printed (0=Def)
    integer :: tight_p_conv ! local copied to qmmm_nml%tight_p_conv - Controls convergence of density matrix. 0 (Def) = 0.05*sqrt(SCFCRT)
                            ! 1 = Converged both Energy and Density to SCFCONV
    integer :: printcharges !Local copied to qmmm_nml%printcharges as a logical. 1 = true - print mulliken and cm1a and cm2a charges
                            !on every step. 0 = false = don't print charges. Default = 0 (.false.)
    integer :: printdipole  !Local copied to qmmm_nml%printdipole as an integer 1 = QM dipole moment, 2 = QM + MM dipole moment, (0=Def)
    integer :: peptide_corr !Local copied to the logical qmmm_nml%peptide_corr
                            !Add MM correction to peptide linkages 0 = No (Default), 1 = Yes.
    integer :: itrmax       !Local copied to qmmm_nml%itrmax - Maximum number of scf cycles to run
                            !before assuming convergence has failed (default = 1000)
    integer :: printbondorders !Local copied to qmmm_nml%printbondorders as a logical.
                               ! 1 = true - print bondorders at the end of the
                               ! calculation. 0 = false = dont print bondorders. Default = 0 (.false.)
    integer :: qmshake      !Local copied to qmmm_nml%qmshake - shake QM atoms if ntc>1?
    integer :: qmmmrij_incore !Flag to store rij between qm-mm pairs and related equations in memory.
                              !1 (default) = store in memory. 0 = calc on fly.
    integer :: qmqm_erep_incore !Flag to store QM-QM 1 electron repulsion integrals in memory or to calculate
                                !them on the fly. Only available with QM-QM analytical derivatives.
                                !1 (default) = store in memory. 0 = calc on fly.
    integer :: pseudo_diag      !Whether to allow pseudo diagonalisations to be done when possible in SCF.
                                !0 (default) = Always do full diagonalisations.
                                !1 = do pseudo diagonalisations when possible.
    integer :: qm_ewald          !0 (default) do only regular QM-MM interaction in periodic calculations.
                                 !1           do ewald based periodic QM-MM interactions.
                                 !2           do ewald based periodic QM-MM but with the QM image charges
                                 !            fixed at the previous steps mulliken charges during the SCF.
    integer :: qm_pme            !0 use regular Ewald for doing QM-MM interactions when qm_ewald>0.
                                 !1 (default) use PME to do the reciprocal sum.
    integer :: kmaxqx, kmaxqy, kmaxqz !Maximum K space vectors
    integer :: ksqmaxq !Maximum K squared values for spherical cutoff in k space.
    _REAL_ :: kappa  ! the ewald coefficient for QM region ewald calculations
    integer :: writepdb
    integer :: qmmm_int !QM-MM interaction method
    integer :: adjust_q
    integer :: diag_routine !Controls diagonalization routine to use in SCF.
#ifdef OPENMP
   integer :: qmmm_omp_max_threads !Maximum number of openmp threads to use for parallel QMMM routines
#endif
    integer :: density_predict !Controls prediction of density matrix for next SCF step.
    integer :: fock_predict !Controls prediction of Fock matrix for next SCF step.
    integer :: vsolv ! = 0 by default, = 1 for simple vsolv QM/MM, = 2 for adaptive QM/MM with vsolv

    _REAL_ :: fockp_d1 !prefactor for fock matrix prediction.
    _REAL_ :: fockp_d2 !prefactor for fock matrix prediction.
    _REAL_ :: fockp_d3 !prefactor for fock matrix prediction.
    _REAL_ :: fockp_d4 !prefactor for fock matrix prediction.
    logical :: mdin_qmmm=.false.

    integer :: idc
    integer :: divpb
    !! (GMS)
    _REAL_  :: chg_lambda       ! Charge scaling factor for free energy calculation
    !! DFTB options
    integer :: dftb_maxiter     ! Max # of iterations before resetting Broyden (default: 70 ) ==> qmmm_nml%dftb_maxiter
    integer :: dftb_disper      ! Use dispersion?  (default: 0 = false) ==> qmmm_nml%dftb_disper
    integer :: dftb_chg         ! DFTB CM3 charges (default: 0 = Mulliken, 1 = CM3) ==> qmmm_nml%dftb_chg
    _REAL_  :: dftb_telec       ! Electronic temperature, in Kelvins. (Default = 0.0K) ==> qmmm_nml%dftb_telec
    _REAL_  :: dftb_telec_step  ! Telec step size for convergence accelerator (Default = 0.0K) ==> qmmm_nml%dftb_telec_step
    character(Len=256) :: dftb_3rd_order  ! 3rd order SCC-DFTB (default: 'NONE'== No third order)
                                        !     'PA' == Do 3rd order, Proton Affinities parameterization
                                        !     'PR' ==               Phosphate reactions parameterization
                                        !     'READ' == read the parameters from a user-specified file (TO IMPLEMENT)
    _REAL_ :: r_switch_lo    !Lower bound of the QM/MM switching function
    _REAL_ :: r_switch_hi    !Upper bound of the QM/MM switching function
    integer :: qmmm_switch   !0           Turn off QM/MM switching function
                             !1           Turn on QM/MM switching function

    !#include "memory.h"

    !Apparently you can't use a pointer in a namelist :-( Therefore
    !we need a local scratch array that will be big enough that
    !the iqmatoms list never exceeds it
    integer :: iqmatoms( max_quantum_atoms )
    integer :: natom;

    character(len=1024) :: qmmask
    character(len=12) :: qm_theory
         !Options=PM3,AM1,MNDO,PDDG-PM3,PM3PDDG,PDDG-MNDO,PDDGMNDO,
         !        PM3-CARB1,PM3CARB1,DFTB,SCC-DFTB,RM1,PM6,PM3-ZnB,PM3-MAIS
         !        EXTERNAL (for external programs like ADF/GAMESS/TeraChem)
    integer, dimension(:), pointer :: isqm
    integer :: ier=0
    character(len=80) :: parameter_file
    logical :: qxd
    logical :: calcxdens
    real(8) :: ceps, cosmo_scf_ftol,cosmo_scf_maxcyc,index_of_refraction,onsager_radius
    real(8) :: Ex, Ey, Ez, linmixparam  
    integer :: nspa, solvent_model, potential_type, EF
    logical :: doZ
    integer :: xlbomd_flag,K
    _REAL_  :: dt2w2,xlalpha   

    namelist /qmmm/ qmcut, iqmatoms,qmmask,qmgb,qm_theory, qmtheory, &
        qmcharge, qmqmdx, qmqmdx_exc, verbosity, tight_p_conv, scfconv, & ! CML 7/10/12
        errconv,ndiis_matrices,ndiis_attempts, &  !+TJG 01/26/2010
        parameter_file, qxd,&
        printcharges, printdipole, peptide_corr, itrmax, qmshake, &
        qmqm_erep_incore, qmmmrij_incore, &
        lnk_dis, lnk_atomic_no, lnk_method, spin, pseudo_diag,   &
        pseudo_diag_criteria, &
        qm_ewald, qm_pme, kmaxqx, kmaxqy, kmaxqz, ksqmaxq, kappa, &
        writepdb, qmmm_int, adjust_q, diag_routine, &
        density_predict, fock_predict, &
        fockp_d1, fockp_d2, fockp_d3, fockp_d4, idc, divpb, &
        dftb_maxiter, dftb_disper, dftb_3rd_order, dftb_chg, &
        dftb_telec, dftb_telec_step, printbondorders, &
        qmmm_switch, r_switch_lo, r_switch_hi, &
#ifdef OPENMP
                   qmmm_omp_max_threads, &
#endif
        !Geometry optimization
        maxcyc, ntpr, grms_tol, struct_opt_state, &
        !Excited state and davidson
        exst_method,dav_guess,dav_maxcyc,ftol0,ftol1,printtd, &
        !Solvent Model and E-Field parameters
        ceps, chg_lambda, vsolv, nspa, solvent_model, potential_type, cosmo_scf_ftol, cosmo_scf_maxcyc,&
        doZ,index_of_refraction,onsager_radius,EF,Ex,Ey,Ez,linmixparam, &
        !XL-BOMD parameters
        xlbomd_flag,K,dt2w2,xlalpha, &
        !Cross densities
        calcxdens, &
        !Nuclear normal modes
        do_nm,deltaX

    !the input value of excNin MUST NOT be changed.
    excN=excNin;
    ! Setup defaults
    qmcut = 9999.d0
    use_pme = 0
    ntb = 0
    maxcyc = 0
    grms_tol = 0.02
    ntpr=10

    ! Default parameters of excited state calculations
    excN = 0 ! Default is to not run the Davidson procedure at all
    struct_opt_state = 0 ! Optimize the ground state by default !JAB I think this parameter is deprecated now
    printtd=0 ! Do not print tds by default
    exst_method=1 ! CI singles
    dav_guess=0 ! use previous davidson result as guess
    ftol=0.d-5   ! Min tolerance (|emin-eold|)
    ftol0=1.d-5 !  Acceptance tol.(|emin-eold|)
    ftol1=1.d-4 ! Accept.tol.for residual norm
    dav_maxcyc=100 ! Maximum cycles for davidson diagonalization
   
    ! Default for COSM and Solvent
    ceps=1.d0 ! no dielectric screening by default
    index_of_refraction=1.d0
    onsager_radius=4.0
    solvent_model=0 !no solvent model
    potential_type=1 !COSMO default
    cosmo_scf_ftol=1.d-3
    cosmo_scf_maxcyc=300
    doZ=.true.
    nspa=42
    EF=0;    !Electric Field Flag (0:none,1:GS+ES,2:ES only)
    Ex=0.0;    !Electric field vectors a.u.
    Ey=0.0;
    Ez=0.0;
    linmixparam=1.0;

    !Other parameters
    lnk_dis=1.09d0  !Methyl C-H distance
    lnk_atomic_no=1 !Hydrogen
    lnk_method=1 !treat MMLink as being MM atom.
    qmgb = 2 !Gets set to zero if igb==6 or igb==0.
    qm_theory = ''
    qmtheory = RETIRED_INPUT_OPTION
    qmcharge = 0
    spin = 1
    qmqmdx = 1 !Depricated in NAESMD/SQM integration
    qmqmdx_exc=2 ! CML 7/10/12 Numerical excited state derivatives by default !Depricated
    verbosity = -1
    parameter_file=''
    qxd=.false.

    ! Defaults for XL-BOMD
    xlbomd_flag = 0
    K = 5
    dt2w2 = 1.82
    xlalpha = 18d-3

    !Cross densities
    calcxdens=.false.

    !Default for nuclear normal modes
    do_nm = 0  !Flag for calculating nuclear normal modes
    deltaX = 1.0d-4 !Displacement for Hessian calculation, A [1.0d-4]

    !+TJG 01/26/2010

    ! defaults for stand-alone (geom. opt.) code:
    !   dac: for now, just use the same values as in the non-stand-alone
    !        code, but these should be updated once we figure out the best
    !        values for geometry optimization
    tight_p_conv = 0
    !We want inputs to be in eV, SCF algorithm works in kcal/mol
    scfconv = 1.0D-6
    errconv = 1.0d-1
    ndiis_matrices = 6
    ndiis_attempts = 0

    !-TJG 01/26/2010
    printcharges = 0
    printbondorders = 0
    printdipole = -1
    peptide_corr = 0
    itrmax = 1000
    qmshake = 1
    qmmask=''
    iqmatoms(1:max_quantum_atoms) = 0
    qmmmrij_incore = 1
    qmqm_erep_incore = 1
    pseudo_diag = 1
    pseudo_diag_criteria = 0.05d0
    qm_ewald=1 !Default is to do QMEwald, with varying charges, if ntb=0 or use_pme=0 then this will get turned off
    qm_pme = 1 !use pme for QM-MM
    kmaxqx=8; kmaxqy=8; kmaxqz=8    !Maximum K space vectors
    kappa=-1.0
    ksqmaxq=100 !Maximum K squared values for spherical cutoff in k space.
    writepdb = 0 !Set to 1 to write a pdb on the first step with just the QM region in it.
    qmmm_int = 1 !Default, do full interaction without extra Gaussian terms for PM3 / AM1 etc.
    adjust_q = 2 !Default adjust q over all atoms.
    diag_routine = 1 !Use default internal diagonalizer.
#ifdef OPENMP
   qmmm_omp_max_threads = 1 !Use just 1 openmp thread by default.
#endif
    density_predict = 0 !Use density matrix from previous MD step.
    fock_predict = 0 !Do not attempt to predict the Fock matrix.
    fockp_d1 = 2.4d0
    fockp_d2 = -1.2d0
    fockp_d3 = -0.8d0
    fockp_d4 = 0.6d0
    idc = 0
    divpb = 0
    vsolv = 0 ! by default do not use simple vsolv QM/MM or adaptive QM/MM based on vsolv
    qmmm_switch = 0              !Use QM/MM switching function
    r_switch_hi = qmcut          !Set the default value to be equal to qmcut
    r_switch_lo = r_switch_hi - 2.0D0  !Set the default value to be 2 Angstrom shorter than r_switch_hi

    !DFTB
    dftb_maxiter     = 70
    dftb_disper      = 0
    dftb_chg         = 0
    dftb_telec       = 0.0d0
    dftb_telec_step  = 0.0d0
    chg_lambda  = 1.0d0
    dftb_3rd_order   = 'NONE'

    !Read qmmm namelist
    rewind fdes_in

    call nmlsrc('qmmm',fdes_in,ifind)
    if (ifind /= 0) mdin_qmmm=.true.

    if(nstep>0) then
        write(6,*) "nstep successful"
    endif

    !Read qmmm namelist
    rewind fdes_in
    if ( mdin_qmmm ) then
        read(fdes_in,nml=qmmm)
    else
        write(fdes_out, '(1x,a,/)') 'Could not find qmmm namelist'
        call mexit(fdes_out,1)
    endif
    cosmo_c_struct%ceps=ceps
    cosmo_c_struct%nspa=nspa
    cosmo_c_struct%solvent_model=solvent_model
    cosmo_c_struct%potential_type=potential_type
    cosmo_c_struct%cosmo_scf_ftol=cosmo_scf_ftol
    cosmo_c_struct%cosmo_scf_maxcyc=cosmo_scf_maxcyc
    cosmo_c_struct%doZ=doZ
    cosmo_c_struct%index_of_refraction=index_of_refraction
    cosmo_c_struct%onsager_radius=onsager_radius
    cosmo_c_struct%EF=EF
    cosmo_c_struct%Ex=Ex
    cosmo_c_struct%Ey=Ey
    cosmo_c_struct%Ez=Ez
    cosmo_c_struct%linmixparam=linmixparam
    xlbomd_struct%xlbomd_flag=xlbomd_flag 
    xlbomd_struct%xlalpha=xlalpha 
    xlbomd_struct%dt2w2=dt2w2 
    xlbomd_struct%K=K 
    !Begin qmmm input checks
    if(printbondorders>0) then
        write(6,*) "Print bonding orders under development. Do not use."
        stop
    endif
    
    if(density_predict>0) then
        write(6,*) "density_prediction under development. Do not use."
        stop
    endif
    
    if(dav_guess>1) then
        write(6,*) "XL-BOMD under development. dav_guess can be 1 or 0. Do not use."
        stop
    endif
    
    if(abs(cosmo_c_struct%index_of_refraction-100)>1d-1) then
        write(6,*) "Index of refraction has unverified effects. Do not use."
        stop
    endif
    !End qmmm input checks

    !convert scfconv to kcal/mol from input eV
    scfconv=scfconv*evkcal
    
    !Set Solvent_model
    !The input format is: (0) None, (1) Linear response, (2) Vertical excitation, ! or (3) State-specific [0]
    !Going forward the code treats: (0) None, (1) LR, (2) Nonequilibrium SS, (3) Same as last with Xi, (4) Equilibrium SS, (5) Same as last with Xi [0]
    if(cosmo_c_struct%solvent_model==3) then
        cosmo_c_struct%solvent_model=4
    else if(cosmo_c_struct%solvent_model>4) then
        write(fdes_out,*) 'Invalid cosmo_c_struct%solvent_model'
        call mexit(fdes_out,1)
    endif
    
    !Set cosmo_c_struct%potential_type
    !The input format is:  (1) COSMO or (2) Onsager [1]
    !Going forward the code treats:  (2) Onsager (3) COSMO, (4) Testing (0) Normal Correlation
    if(cosmo_c_struct%potential_type==1) then
        cosmo_c_struct%potential_type=3
    elseif(cosmo_c_struct%potential_type==2) then
        cosmo_c_struct%potential_type=2
    endif

    !set verbosity if not set at read
    if(verbosity==-1) then
        if(nstep>0 .OR. maxcyc>0) then
            verbosity=1
        else
            verbosity=5
        endif
    endif
    
    !et print dipoles if not set at read
    if(printdipole==-1) then
        if(nstep>0 .OR. maxcyc>0) then
            printdipole=1
        else
            printdipole=2
        endif
    endif

    !AWG NEW
    call CheckRetiredQmTheoryInputOption(qmtheory)
    call set(qmmm_nml%qmtheory, qm_theory)
    !AWG END NEW
   
    !  Read-in the user-defined parameter file
    !  TL (Rutgers, 2011)
    call ReadParameterFile(parameter_file, qm2_struct%parameterEntries, &
	qm2_struct%ParameterFileExisting)
    ! turn on OPNQ if necessary 
    qmmm_opnq%useOPNQ=qxd

    qm2ds%calcxdens=calcxdens !JAKB for calculating cross density

    ! Disable EXTERN in SQM since
    ! it does not make sense for SQM to be calling the external ADF interface.
    if (qmmm_nml%qmtheory%EXTERN) then
        call sander_bomb('read_qmmm_namelist','External interface is not supported in SQM.', &
            '(qm_theory = ''EXTERN'')')
    end if
    if (ncharge_in > 0) then ! we have external charge

        if (maxcyc > 0) then
            ! external charge calculation is not supported in gradient minimization
            call sander_bomb('read_qmmm_namelist','maxcyc > 0 but external charge found.', &
                'external charge calculation is for single point energy calculations only.')
        end if
        natom = natom_inout + ncharge_in
        qmmm_struct%nquant = natom_inout
        qmmm_struct%nlink  = 0
        qmmm_struct%nquant_nlink = natom_inout
        qmmm_struct%qm_mm_pairs = ncharge_in
        do i=1,natom_inout
            iqmatoms(i) = i
        end do
        !update natom_inout with the number of external charges
        natom_inout = natom
    else
        natom = natom_inout
        qmmm_struct%nquant = natom
        qmmm_struct%nlink  = 0
        qmmm_struct%nquant_nlink = natom
        do i=1,natom
            iqmatoms(i) = i
        end do
    end if
    ! Test to see if QM atom selection is legal.
    call validate_qm_atoms(iqmatoms,qmmm_struct%nquant,natom)

    !  check we don't bust our statically allocated max_quantum_atoms
    call int_legal_range('QMMM: (number of quantum atoms) ', &
        qmmm_struct%nquant, 1, max_quantum_atoms )

    call qmsort(qmmm_struct, iqmatoms) !ensure the list of qm atoms is sorted numerically

    ! --- Variable QM solvent region - has to be very early here because we
    !     will be changing nquant and iqmatoms.  ---
    qmmm_nml%vsolv = vsolv
    if (qmmm_nml%vsolv > 0) then
        write(fdes_out,*) 'SQM does not support the use of nearest_qm_solvent.'
        call mexit(6,1)

    end if
    ! --- End Variable QM water region ---
    call float_legal_range('QMMM: (Force Criterion for Geometry Optimization) ', grms_tol,0.0D0,1.0D0)
    call int_legal_range('QMMM: (Print results every ntpr cycles ) ', ntpr,0,99999999)
    call float_legal_range('QMMM: (Acceptance tol.(|emin-eold|)) ', ftol0,1.0D-16,1.0D0)
    call int_legal_range('QMMM: (Solvent Model ) ', solvent_model,0,3)
    call int_legal_range('QMMM: (Potential Type) ', potential_type,1,2)
    call float_legal_range('QMMM: (Onsager Radius) ', onsager_radius,1.0D0, 1.0D9)
    call float_legal_range('QMMM: (Linear mixing parameter for vertical excitation) ', linmixparam, 0.0D0,1.0D3)
    call float_legal_range('QMMM: (Vertical excitation or state-specific SCF tolerance) ', cosmo_scf_ftol,1.0D-16,1.0D0)

    call float_legal_range('QMMM: (QM-MM Cutoff) ', qmcut,0.0D0,1.0D30)
    call int_legal_range('QMMM: (variable solvent VSOLV) ', vsolv,0,3)
    call int_legal_range('QMMM: (QM GB Method) ', qmgb,0,3)
    call int_legal_range('QMMM: (QM-QM Derivatives,Depricated in NAESMD) ', qmqmdx,1,2)
    call int_legal_range('QMMM: (Excited State Derivatives,Depricated in NAESMD) ', qmqmdx_exc,1,2) ! CML 7/10/12
    call int_legal_range('QMMM: (Verbosity) ', verbosity,0,5)
    call int_legal_range('QMMM: (Max SCF Iterations) ', itrmax,-10000000,10000000)
    call int_legal_range('QMMM: (Shake on QM atoms) ', qmshake,0,1)
    call int_legal_range('QMMM: (Density Matrix Convergence) ', tight_p_conv,0,1)
    call float_legal_range('QMMM: (SCF Convergence) ', scfconv,1.0D-16,1.0D0)
    !+TJG 01/26/2010
    call float_legal_range('QMMM: (Error Matrix Convergence) ', errconv,1.0D-16,1.0D0)
    call int_legal_range('QMMM: (Max num matrices in DIIS extrapolation) ', ndiis_matrices,1,20)
    call int_legal_range('QMMM: (Max num of DIIS attempts) ', ndiis_attempts,0,1000)
    !-TJG 01/26/2010
    call int_legal_range('QMMM: (PRINT CHARGES) ', printcharges,0,1)
    call int_legal_range('QMMM: (PRINT BONDORDERS) ',printbondorders,0,1)
    call int_legal_range('QMMM: (PRINT QM/Dipole) ', printdipole,0,2)
    call int_legal_range('QMMM: (Spin State) ', spin,1,1)
    !RCW: Currently limit spin state to singlets only since the code for spin>1 does not exist / work at present.
    !     WARNING - IF WE LATER ALLOW SPIN>1 qm2_densit will need updating.
    call int_legal_range('QMMM: (Peptide Correction) ',peptide_corr,0,1)
    call int_legal_range('QMMM: (QM-MM RIJ in Core) ',qmmmrij_incore,0,1)
    call int_legal_range('QMMM: (QM-QM E-Rep in Core) ',qmqm_erep_incore,0,1)
    call int_legal_range('QMMM: (Link Atomic Number) ',lnk_atomic_no,1,numberElements)
    call int_legal_range('QMMM: (QM-MM Link Method) ',lnk_method,1,2)
    call int_legal_range('QMMM: (Pseudo Diag) ',pseudo_diag,0,1)
    call int_legal_range('QMMM: (QM Ewald) ',qm_ewald,0,2)
    call int_legal_range('QMMM: (QM PME) ',qm_pme,0,1)
    call int_legal_range('QMMM: (QM Ewald kmaxqx) ',kmaxqx,1,99999999)
    call int_legal_range('QMMM: (QM Ewald kmaxqy) ',kmaxqy,1,99999999)
    call int_legal_range('QMMM: (QM Ewald kmaxqz) ',kmaxqz,1,99999999)
    call int_legal_range('QMMM: (QM Ewald ksqmaxq) ',ksqmaxq,1,kmaxqx*kmaxqy*kmaxqz)
    call int_legal_range('QMMM: (QM-MM qmmm_int) ',qmmm_int,0,5)
    call int_legal_range('QMMM: (QM-MM adjust_q) ',adjust_q,0,2)
    call int_legal_range('QMMM: (QM-MM diag_routine) ',diag_routine,0,7)
#ifdef OPENMP
   call int_legal_range('QMMM: (QM-MM qmmm_omp_max_threads) ',qmmm_omp_max_threads,1,32)
#endif
    call int_legal_range('QMMM: (QM-MM density_predict) ',density_predict,0,2)
    call int_legal_range('QMMM: (QM-MM fock_predict) ',fock_predict,0,1)
    call float_legal_range('QMMM: (Pseudo Diag Criteria) ',pseudo_diag_criteria,1.0D-12,1.0D0)
    call int_legal_range('QMMM: (QM-MM qmmm_switch) ',qmmm_switch,0,1)
    call float_legal_range('QMMM: (QM-MM r_switch_lo) ',r_switch_lo,0.0D0,1.0D30)
    call float_legal_range('QMMM: (QM-MM r_switch_hi) ',r_switch_hi,0.0D0,1.0D30)
    if (lnk_dis>0.0d0) then
        !if lnk_dis is less than 0.0d0 then the link atom is just placed on top of
        !the MM link pair atom.
        call float_legal_range('QMMM: (Link Atom Distance) ',lnk_dis,0.7D0,4.0D0)
    endif

    !! GMS
    call float_legal_range('QMMM: (QM-MM chg_lambda)'     , chg_lambda     , 0.0D0 , 1.0D0  )
    qmmm_nml%chg_lambda  = chg_lambda
    !! DFTB
    call int_legal_range(  'QMMM: (QM-MM dftb_maxiter ) ' , dftb_maxiter   , 1     , 10000  )
    call int_legal_range(  'QMMM: (QM-MM dftb_disper) '   , dftb_disper    , 0     , 1      )
    call int_legal_range(  'QMMM: (QM-MM dftb_chg   ) '   , dftb_chg       , 0     , 1      )
    call float_legal_range('QMMM: (QM-MM dftb_telec)'     , dftb_telec     , 0.0D0 , 1.0D4  )

    ! Checking sanity of excited state parameters
    call int_legal_range(  'QMMM: (QM-MM Number of Excited States) ',excN, 0,1)
    call int_legal_range(  'QMMM: (QM-MM State to Optimize) ',struct_opt_state,0,excN)
    call int_legal_range('QMMM: (QM-MM Excited State Method) ',exst_method,1,2)
    call int_legal_range('QMMM: (QM-MM use (1) or not (0) Davidson guess or (2) for XL-BOXMD) ', &
        dav_guess,0,3)
    call int_legal_range('QMMM: (QM-MM use -10-500 for max cycles of Davidson, negative for fixed number) ', dav_maxcyc,-10,500)

    ! Checking COSMO parameters
    call float_legal_range('QMMM: (COSMO eps)',cosmo_c_struct%ceps,1.0d0,1.0d6)

    if (dftb_3rd_order /= 'NONE') then
        call check_dftb_3rd_order(dftb_3rd_order)
        qmmm_nml%dftb_3rd_order = dftb_3rd_order
    endif

    qmmm_nml%dftb_maxiter   = dftb_maxiter
    qmmm_nml%dftb_disper      = dftb_disper
    qmmm_nml%dftb_chg         = dftb_chg
    qmmm_nml%dftb_telec       = dftb_telec
    qmmm_nml%dftb_telec_step  = dftb_telec_step

    if (dftb_chg > 0) printcharges=1

    qmmm_nml%qmcut = qmcut
    qmmm_nml%qmcut2 = qmcut*qmcut
    qmmm_nml%lnk_dis = lnk_dis
    qmmm_nml%lnk_atomic_no = lnk_atomic_no
    qmmm_nml%lnk_method = lnk_method
    !DFTB Limitations - current things not supported in DFTB
    !These are silent limitations that are non fatal - just to
    !avoid problems with the default. Fatal errors are handled
    !later in this routine.
    if (qmmm_nml%qmtheory%DFTB) then
        qmmmrij_incore=0
        qmqm_erep_incore=0
        pseudo_diag=0
        qmqmdx=1
        qmqmdx_exc=1 ! CML 7/10/12
    end if

    !Divcon limitations - current things not supported in sander.DIVCON
    !These are silent changes to remove defaults. Fatal errors are handled
    !later in this routine.
    if (idc /= 0) then
        qmmmrij_incore=0
        qmqm_erep_incore=0
        pseudo_diag=0
    end if

    !qmgb values:
    ! 0 - do GB but leave QM charges as zero and add nothing to Fock matrix. This is like a vacuum
    !     QM molecule in a solvated MM system.
    ! 1 - do GB using the prmtop fixed resp charges for the GB calculation.
    ! 2 - do GB using Mulliken charges that are consistent with the GB field by modifying the fock
    !     matrix at every SCF step. (default)
    ! 3 - do GB using QM gas phase Mulliken charges - This is really a debugging option since the charges
    !     will not be consistent with the GB field since the fock matrix is not modified. This similarly
    !     means that the gradients will not be accurate. A warning will be printed at every QM call if this
    !     option is selected.

    !Make sure igb in &cntrl namelist is compatible with qmgb setting.
    if (igb==0 .or. igb==6) then
        !no qmgb available
        qmgb = 0
    end if
    !Print warning about qmgb being for debugging only.
    if (qmgb==3) then
        write(fdes_out,*) "QMMM: ------------------------------ WARNING --------------------------------"
        write(fdes_out,*) "QMMM: qmgb = 3 is designed for debugging purposes only. It gives GB"
        write(fdes_out,*) "QMMM:          energies based on gas phase QM Mulliken charges. These charges"
        write(fdes_out,*) "QMMM:          are NOT consistent with the GB field felt by the QM region and"
        write(fdes_out,*) "QMMM:          so any gradients calculated using this approach will"
        write(fdes_out,*) "QMMM:          NOT BE ACCURATE."
        write(fdes_out,*) "QMMM:          This option is really designed for:"
        write(fdes_out,*) "QMMM:                SINGLE POINT ENERGY EVALUATIONS ONLY"
        write(fdes_out,*) "QMMM: ------------------------------ WARNING --------------------------------"
    end if
    qmmm_nml%qmgb = qmgb

    qmmm_struct%AM1_OR_PM3 = (qmmm_nml%qmtheory%AM1 .or. qmmm_nml%qmtheory%AM1D .or. qmmm_nml%qmtheory%PM3 &
        .OR. qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PM3CARB1 .OR. &
        qmmm_nml%qmtheory%RM1 .OR. qmmm_nml%qmtheory%PDDGPM3_08 .OR. qmmm_nml%qmtheory%PM3ZNB)
    qmmm_struct%PDDG_IN_USE = (qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PDDGMNDO  &
        .OR. qmmm_nml%qmtheory%PDDGPM3_08 )
    if (qmmm_int == 3 .OR. qmmm_int == 4) then
        if (qmmm_nml%qmtheory%PM3) then
            qmmm_struct%PM3MMX_INTERFACE = .true.
        else
            qmmm_struct%PM3MMX_INTERFACE = .false.
            call sander_bomb('read_qmmm_namelist','qmmm_int == 3/4 (Modified QM-MM interface) but qm_theory /= PM3.', &
                'qmmm_int == 3/4 is only supported along with PM3 Hamiltonian.')
        end if
    else
        qmmm_struct%PM3MMX_INTERFACE = .false.
    end if
    qmmm_nml%qmcharge = qmcharge
    qmmm_nml%spin = spin
    qmmm_nml%verbosity = verbosity
    qmmm_nml%itrmax = itrmax
    qmmm_nml%qmshake = qmshake
    qmmm_nml%pseudo_diag_criteria = pseudo_diag_criteria

    ! Analytical or numerical integral derivatives for semiempirical methods:
    ! Note: for d orbitals only numerical derivatives are available
    !       Thus switch to numerical derivatives for MNDO/d and AM1/d
    !       For PM6 qmmm_nml%qmqm_analyt will be adjusted after parameters are
    !       read in qm2_load_params_and_allocate since we can do analytical
    !       derivatives if we don't have an element with d orbitals
    if (qmqmdx /= 1 .or. qmmm_nml%qmtheory%MNDOD .or. qmmm_nml%qmtheory%AM1D) then
        ! Do numerical QM-QM derivatives in qm2
        qmmm_nml%qmqm_analyt = .false.
    else
        !Do analytical QM-QM dericatives in qm2
        qmmm_nml%qmqm_analyt = .true.
    end if

    ! CML Added for excited state support 7/10/12
    ! CML False should be used for now because there is no support for analytical
    ! CML excited state derivatives

    if (qmqmdx_exc /= 1 .or. qmmm_nml%qmtheory%MNDOD .or. qmmm_nml%qmtheory%AM1D) then
        ! Do numerical QM-QM derivatives in qm2
        qmmm_nml%qmqm_exc_analyt = .false.
    else
        !Do analytical QM-QM dericatives in qm2
        qmmm_nml%qmqm_exc_analyt = .true.
    end if

    if (tight_p_conv /= 1) then
        ! Loose density matrix convergence (0.05*sqrt(SCFCRT))
        qmmm_nml%tight_p_conv = .false.
    else
        ! Tight density matrix convergence (SCFCRT)
        qmmm_nml%tight_p_conv = .true.
    end if
    !Write a warning about excessively tight convergence requests.
    if ( scfconv < 1.0D-12 ) then
        write(fdes_out,'(" QMMM: WARNING - SCF Conv = ",G8.2)') scfconv
        write(fdes_out,*) "QMMM:           There is a risk of convergence problems when the"
        write(fdes_out,*) "QMMM:           requested convergence is less that 1.0D-12 kcal/mol."
    end if
    qmmm_nml%scfconv = scfconv

    !How tight do we want the density convergence?
    if (qmmm_nml%tight_p_conv) then
        qmmm_nml%density_conv = qmmm_nml%scfconv
    else
        qmmm_nml%density_conv = 0.05D0 * sqrt(qmmm_nml%scfconv)
    end if

    !+TJG 01/26/2010
    qmmm_nml%errconv = errconv
    qmmm_nml%ndiis_matrices = ndiis_matrices
    qmmm_nml%ndiis_attempts = ndiis_attempts
    !-TJG 01/26/2010

    if ( printcharges /= 1) then
        qmmm_nml%printcharges=.false.
    else
        qmmm_nml%printcharges=.true.
    end if

    qmmm_nml%printdipole=printdipole
    qmmm_nml%printtd=printtd

    if ( printbondorders /= 1) then
        qmmm_nml%printbondorders=.false.
    else
        qmmm_nml%printbondorders=.true.
    end if
    if ( peptide_corr == 0) then
        qmmm_nml%peptide_corr = .false.
    else
        qmmm_nml%peptide_corr =  .true.
    end if
    if ( qmmmrij_incore == 0 .or. qmmm_int==0 .or. qmmm_int == 5 ) then
        qmmm_nml%qmmmrij_incore = .false.
    else
        qmmm_nml%qmmmrij_incore = .true. !Only available with qmmm_int>1 and qmmm_int /= 5
    end if
    if ( qmqm_erep_incore == 0 .or. qmqmdx == 2 ) then
        qmmm_nml%qmqm_erep_incore = .false.
    else
        !Only available with analytical derivatives.
        qmmm_nml%qmqm_erep_incore = .true.
    end if
    if ( pseudo_diag == 1 ) then
        qmmm_nml%allow_pseudo_diag = .true.
    else
        qmmm_nml%allow_pseudo_diag = .false.
    end if
    qmmm_nml%qm_ewald = qm_ewald
    qmmm_nml%ksqmaxq = ksqmaxq
    qmmm_nml%kmaxqx = kmaxqx
    qmmm_nml%kmaxqy = kmaxqy
    qmmm_nml%kmaxqz = kmaxqz
    qmmm_nml%kappa = kappa
    !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
    !have put in the namelist and set the value to false.
    if (ntb==0 .or. use_pme==0) then
        qmmm_nml%qm_ewald = 0
        qmmm_nml%qm_pme = .false.
    end if
    if (qmmm_nml%qm_ewald>0 .and. qm_pme>0) then
        qmmm_nml%qm_pme=.true.
    else
        qmmm_nml%qm_pme=.false.
    end if

    if ( writepdb == 0 ) then
        qmmm_nml%writepdb=.false.
    else
        qmmm_nml%writepdb=.true.
    end if

    qmmm_nml%adjust_q = adjust_q

    qmmm_nml%qmmm_int = qmmm_int

    if (qmmm_nml%qmmm_int == 5) then
        ! Mechanical embedding
        ! Do not use QM and QM/MM Ewald or PME
        write(6,'(a)') 'QMMM: Mechanical embedding in use'
        if ( qmmm_nml%qm_pme ) then
            write(6,'(a)') 'QMMM: WARNING'
            write(6,'(a)') 'QMMM: Switching off QM PME'
            qmmm_nml%qm_pme = .false.
        end if
        if ( qmmm_nml%qm_ewald > 0) then
            write(6,'(a)') 'QMMM: WARNING'
            write(6,'(a)') 'QMMM: Switching off QM Ewald'
            qmmm_nml%qm_ewald = 0
        end if
        ! Prevent adjust_q, there is nothing to adjust
        qmmm_nml%adjust_q = 0
        ! Set qmcut as zero
        qmmm_nml%qmcut = 0.1d0
        qmmm_nml%qmcut2 = qmmm_nml%qmcut * qmmm_nml%qmcut
    end if

    qmmm_nml%idc = idc
    qmmm_nml%divpb = divpb

    if (qmmm_switch == 1) then
        qmmm_nml%qmmm_switch = .true.
    else
        qmmm_nml%qmmm_switch = .false.
    end if
    qmmm_nml%r_switch_lo = r_switch_lo
    qmmm_nml%r_switch_hi = r_switch_hi
    if ( (qmmm_nml%r_switch_hi - qmmm_nml%r_switch_lo) < 0.0D0 ) then
        call sander_bomb('read_qmmm_namelist', &
            & 'r_switch_hi is smaller than r_switch_lo!', &
            & 'please try a different set.')
    end if

    !Setup some specific calculation flags that depend on namelist variables.
    !Need to make sure these get copied to other threads in an MPI run.

    !Will we be calculating the Mulliken charges on every SCF iteration?
    !Default is no. Will be set to true in a bit if certain options, such as qm_ewald
    !require it.
    qm2_struct%calc_mchg_scf = .false.

    !DFTB Calculates Mulliken charges anyway so we might as well store them in the correct place.
    if (qmmm_nml%qmtheory%DFTB) qm2_struct%calc_mchg_scf = .true.
    !At this point we know nquant and natom so we can allocate our arrays that depend on nquant or natom
    !Note if this is a LES run qmmm_struct%nquant is nqaunt
    !Note non master mpi threads need to call this allocation routine manually themselves.
    qmmm_nml%nquant = qmmm_struct%nquant
    qmmm_struct%natom = natom
    call allocate_qmmm(qmmm_scratch, qmmm_nml, qmmm_struct, natom )
    qmmm_struct%iqm_atomic_numbers(1:qmmm_struct%nquant_nlink) = atnum(1:qmmm_struct%nquant_nlink)
    qmmm_nml%iqmatoms(1:qmmm_struct%nquant_nlink) = iqmatoms(1:qmmm_struct%nquant_nlink)
    if (ncharge_in > 0) then
        qmmm_struct%qm_xcrd = 0.0D0
        j=0
        do i=1,qmmm_struct%qm_mm_pairs
            qmmm_struct%qm_xcrd(1,i) = excharge(j+1)
            qmmm_struct%qm_xcrd(2,i) = excharge(j+2)
            qmmm_struct%qm_xcrd(3,i) = excharge(j+3)
            qmmm_struct%qm_xcrd(4,i) = excharge(j+4)
            j=j+4
        end do

        if (qmmm_struct%PM3MMX_INTERFACE) then
            allocate(qmmm_struct%qm_mm_pair_atom_numbers(qmmm_struct%qm_mm_pairs), stat=ier)
            REQUIRE(ier == 0)
            do i = 1, qmmm_struct%qm_mm_pairs
                qmmm_struct%qm_mm_pair_atom_numbers(i) = chgatnum(i)
            end do
        end if
    end if

    ! From now on the code uses qmmm_struct%iqmatoms
    ! which will be extended to contain link atom info
    qmmm_struct%iqmatoms(:) = qmmm_nml%iqmatoms(:)

    ! Now we have a list of atom numbers for QM atoms we can build a true false (natom long) list
    ! specifying what the quantum atoms are. Useful for doing quick .OR. operations against other
    ! lists.
    qmmm_struct%atom_mask = .false. !Note, sets entire natom long array to false

    do i = 1, qmmm_struct%nquant
        qmmm_struct%atom_mask(qmmm_nml%iqmatoms(i)) = .true.
    end do

    qmmm_nml%diag_routine = diag_routine
#ifdef OPENMP
   qmmm_nml%qmmm_omp_max_threads = qmmm_omp_max_threads

    !For the time being the number of threads to use for diag and pdiag
    !routines is set to max_threads - later this will be optimized if
    !diag_routine=0.
   qmmm_omp%diag_threads = qmmm_omp_max_threads
   qmmm_omp%pdiag_threads = qmmm_omp_max_threads
#endif
    qmmm_nml%density_predict = density_predict
    qmmm_nml%fock_predict = fock_predict
    qmmm_nml%fockp_d1 = fockp_d1
    qmmm_nml%fockp_d2 = fockp_d2
    qmmm_nml%fockp_d3 = fockp_d3
    qmmm_nml%fockp_d4 = fockp_d4

    ! --- CHECK FOR LIMITATIONS ---

    !--- Mechanical embedding is not supported for GB at the moment. ---
    if ( (igb > 0) .and. (qmmm_nml%qmmm_int == 5) ) then
        call sander_bomb('read_qmmm_nm_and_alloc','Mechanical embedding currently not supported with GB models.', &
            'Cannot have igb > 0 and qmmm_int = 5')
    end if

    !--- You cannot mix Fock prediction with density prediction. ---
    if (qmmm_nml%fock_predict > 0 .and. qmmm_nml%density_predict > 0) then
        call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix and Density matrix prediction are mutually exclusive.', &
            'Cannot have fock_predict > 0 and density_predict > 0')
    end if

    !--- For Fock prediction the 4 pre-factors must sum to 1.0d0 ---
    if (qmmm_nml%fock_predict == 1) then
        if (abs(1.0d0-(qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)) > 1.0d-6) then
            write(6,*) 'QMMM: Failure, fockp_d1 to d4 must sum to 1.0d0 - current sum is', &
                (qmmm_nml%fockp_d1 + qmmm_nml%fockp_d2 + qmmm_nml%fockp_d3 + qmmm_nml%fockp_d4)
            call sander_bomb('read_qmmm_nm_and_alloc','Fock matrix prediction coefficients do not sum to 1.0.', &
                'adjust fockp_d1 to fockp_d4 so that they sum to 1.0.')
        end if
    end if

    !--- You cannot use variable solvent with GB calculations.
    if (qmmm_nml%qmgb > 0 .and. qmmm_nml%vsolv > 0) then
        call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent and qmgb are mutually exclusive.', &
            'Cannot have nearest_qm_solvent > 0 and qmgb > 0')
    end if

    !--- DIVCON LIMITATIONS ---
    !Divcon currently only works with gas phase simulations.
    if (qmmm_nml%idc>0) then
        if (ntb /= 0) then
            !This covers qmewald, qm_pme as well.
            call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but periodic boundaries are in use.', &
                'Periodic boundaries are currently only supported when idc == 0')
        end if
        if (qmmm_nml%qmgb/=0) then
            call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but qmgb/=0. QMMM GB is currently', &
                'only supported when idc == 0')
        end if
        if (qmmm_nml%peptide_corr) then
            call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but peptide_corr /= 0.', &
                'Peptide correction for divcon is handled by the divcon.in file. Set peptide_corr = 0 to proceed.')
        end if
        if (qmmm_nml%printcharges) then
            call sander_bomb('read_qmmm_nm_and_alloc','idc /= 0 (Divcon on) but printcharges /= 0.', &
                'Printcharges is not available with idc > 0.')
        end if
        if (qmmm_nml%vsolv > 0) then
            call sander_bomb('read_qmmm_nm_and_alloc','Nearest QM solvent is not available with idc>0.', &
                'Cannot have nearest_qm_solvent > 0 and idc > 0')
        end if
#ifdef MPI
        !No support for parallel Divcon
    write(6,*) 'Divcon capability (idc>0) can only run in serial mode for now'
    call mexit(6,1)
#endif
        !Write a warning about qm_theory being ignored with Divcon.
        write(6,'("|QMMM: WARNING DIVCON IN USE")')
        write(6,'("|QMMM: qm_theory IS IGNORED WHEN USING DIVCON - QM HAMILTONIAN MUST BE SELECTED")')
        write(6,'("|QMMM: IN DIVCON.IN FILE.")')
    end if
    !--- END DIVCON LIMITATIONS ---

    !--- DFTB LIMITATIONS ---
    if (qmmm_nml%qmtheory%DFTB ) then
        if (qmmm_nml%peptide_corr) then
            call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=DFTB but peptide_corr /= 0.', &
                'Peptide correction is not available, or required for  DFTB. Set peptide_corr = 0 to proceed.')
        end if
    end if
    !--- END DFTB LIMITATIONS ---

    !--- PM6 LIMITATIONS ---
    if (qmmm_nml%qmtheory%PM6 ) then
        if (qmmm_nml%peptide_corr) then
            ! AWG: Do not use peptide correction with PM6 since it has not been parametrized for PM6
            call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=PM6 but peptide_corr /= 0.', &
                'Peptide correction is not available for PM6. Set peptide_corr = 0 to proceed.')
        end if
    end if
    !--- END DFTB LIMITATIONS ---

    !--- EXTERNAL INTERFACE LIMITATIONS ---
    ! ADF/GAMESS support works through an external interface. Currently there
    ! are a number of limitations.
    if (qmmm_nml%qmtheory%EXTERN) then
        ! 1) PME and EWALD are not supported with EXTERN.
        if (qmmm_nml%qm_ewald /= 0) then
            call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=EXTERN but qm_ewald /= 0.', &
                'The external interface does not currently support EWALD or PME.')
        end if
        ! 2) GB is not currently supported with EXTERN.
        if (qmmm_nml%qmgb /= 0) then
            call sander_bomb('read_qmmm_nm_and_alloc','qm_theory=EXTERN but qmgb /= 0.', &
                'The external interface does not currently support Generalized Born.')
        end if
    end if
    !--- END EXTERNAL INTERFACE LIMITATIONS ---
 

    ! --- END CHECK FOR LIMITATIONS ---

    return
end subroutine sqm_read_and_alloc

!  following stub routine to avoid changing qmmm_module.f for sqm:
#ifdef MPI
subroutine qmmm_vsolv_mpi_setup
   return
end subroutine qmmm_vsolv_mpi_setup
#endif
! end module
