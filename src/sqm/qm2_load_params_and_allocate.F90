! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

subroutine qm2_load_params_and_allocate(qm2_params, qmmm_nml, qmmm_mpi, qmmm_opnq, qmmm_scratch, &
                                        xlbomd_struct,qm2_struct, qmmm_struct)


    ! Written by: Ross Walker (TSRI, 2005)
    ! Updates by: Andreas Goetz (SDSC, 2009, 2010, 2011)
    ! Updates by: Taisung Lee (Rutgers, 2011)

    ! This routine should be called before running any qm2 routine calculations.
    ! It is responsible for filling the parameter arrays with the designated
    ! parameters for the method chosen.

    ! All parameters are loaded into qm2_params structure.

    ! Note, allocation is done in this routine for all pointers in the qm2_params structure.
    ! Deallocation is done by deallocate qmmm routine.

    ! Parameter definitions:
    !     Common params for all methods:
    !          core(nquant_nlink) - The core charge on each atom
    !  natomic_orbs(nquant_nlink) - Number of atomic orbitals on each atom.
    !  orb_loc(2,nquant_nlink) - (1,x) = position of first orbital on atom x. 2,x = last orbital on x.
    !  heat_of_form(nquant_nlink) - Gas Phase heat of formation for atom.

    use constants, only : half, EV_TO_KCAL, AU_TO_EV
    use ElementOrbitalIndex, only: NumberElements, &
        MaxValenceOrbitals, SPPrincipalQuantumNumber, DPrincipalQuantumNumber
            
    use QM2_parameters
    use qm2_params_module, only : new
    use qmmm_module, only : qm2_structure, &
        qmmm_mpi_structure, qmmm_scratch_structure, qmmm_opnq_structure
    use MNDOChargeSeparation, only : GetDDAndPho
    use qm2_diagonalizer_module, only : qm2_diagonalizer_setup
    use xlbomd_module, only : xlbomd_structure                                
    use qmmm_struct_module, only : qmmm_struct_type
    use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type

    implicit none
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(xlbomd_structure),intent(inout) :: xlbomd_struct
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qmmm_opnq_structure),intent(inout) :: qmmm_opnq
   type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch

    !Locals
    _REAL_ :: pdiag_guess1, pdiag_guess2, pddg_zaf, pddg_zbf
    _REAL_ :: gssc, gspc, gppc, gp2c, hspc, elec_eng
    _REAL_ :: exponent_temp1, base_temp1, exponent_temp2, base_temp2
    _REAL_ :: HSP1_temp, HSP2_temp, DD1_temp, DD2_temp, DD_diff, DD3_temp, hpp
    _REAL_ :: temp
    integer :: iostmp, ioptmp, int_dummy
    integer :: i, j, k, iat, jat, itmp, iqm_atomic, n_atomic_orb, first_orb, last_orb
    integer :: ns_atoms, nsp_atoms, nspd_atoms, nelectrons, nopen
    integer :: ier=0

    _REAL_::DD(6), PO(9)

#ifdef MPI
   include 'mpif.h'
    integer :: jstart, jend, ia, ib, ja, jb, status
    integer :: loop_extent_begin, loop_extent_end, loop_extent, loop_count
    integer, dimension(:), allocatable :: gather_array !Allocated and deallocated in this routine.
#endif

    logical :: test

    !  Initialize the parameter module
    call InitializeParameter(qm2_struct,qmmm_nml%qmtheory)
    call new(qm2_params, qmmm_struct%qm_ntypes, qmmm_struct%nquant_nlink, &
        qmmm_nml%qmtheory, qmmm_struct%PM3MMX_INTERFACE, qmmm_opnq%useOPNQ)
    !Zero the total heat of formation before calculating it.
    qm2_params%tot_heat_form = 0.0D0
    ! We start by loading in the parameters that are common to all the semi-empirical methods.
    !---------- COMMON PARAMS ---------------
    qm2_struct%norbs=0
    ns_atoms=0
    nsp_atoms=0
    nspd_atoms=0
    nelectrons=-qmmm_nml%qmcharge
    !
    !  define the parameters for each type--most of parameters should be in this way
    !  TL Work
    do i=1,qmmm_struct%qm_ntypes
        qm2_params%sp_quantum_number(i)=SPPrincipalQuantumNumber(qmmm_struct%qm_type_id(i))
        qm2_params%d_quantum_number(i)=DPrincipalQuantumNumber(qmmm_struct%qm_type_id(i))
    end do

    do i=1,qmmm_struct%nquant_nlink
        iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
        qm2_params%core_chg(i)=dble(core_chg(iqm_atomic))
        
        nelectrons=nelectrons+core_chg(iqm_atomic)
         
        n_atomic_orb=natomic_orbs(iqm_atomic)
         
        if (n_atomic_orb==1) ns_atoms=ns_atoms+1
        if (n_atomic_orb==4) nsp_atoms=nsp_atoms+1
        if (n_atomic_orb==9) nspd_atoms=nspd_atoms+1
        ! Check we don't bust any static arrays
        ! DFTB is independent of this and checks are done in qm2_dftb_load_params
        if ( .not. qmmm_nml%qmtheory%DFTB ) then
            if (n_atomic_orb > MaxValenceOrbitals) then
                write (6,*) 'n_atomic_orb of ', n_atomic_orb, ' exceeds MaxValenceOrbitals of MaxValenceOrbitals'
                call sander_bomb('qm2_load_params_and_allocate', &
                    'exceeded max', &
                    'Check qmmm_module.f and parameters.h')
            end if
        end if
        qm2_params%natomic_orbs(i)=n_atomic_orb
        qm2_params%orb_loc(1,i)=qm2_struct%norbs+1
        qm2_params%orb_loc(2,i)=qm2_struct%norbs+n_atomic_orb
        qm2_struct%norbs = qm2_struct%norbs+n_atomic_orb
        qm2_params%tot_heat_form=qm2_params%tot_heat_form+heat_of_form(iqm_atomic)
    end do !i=1,qmmm_struct%nquant_nlink
    !Work out how many 2 electron integrals there will be

      
    !! TL_Work
    qm2_struct%n2el= max(ns_atoms*(ns_atoms-1)/2,0) &
        +10*nsp_atoms*ns_atoms + 45 *nspd_atoms*ns_atoms &
        +max(100*nsp_atoms*(nsp_atoms-1)/2,0) + 450*nspd_atoms*nsp_atoms &
        +max(2025*nspd_atoms*(nspd_atoms-1)/2,0)

    !QMMM e-repul memory depends on QM-MM pair list size so is
    !allocated later on and checked on every call.
    call qm2_allocate_qmqm_e_repul(qmmm_nml, qmmm_mpi, qm2_struct,qmmm_struct, qm2_struct%n2el)

    !Protect DUMB users from STUPID errors
    if (nelectrons > 2*qm2_struct%norbs) then
        if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
            write(6,'(''QMMM: ERROR-number of electrons: '',i5,'' is more'')') nelectrons
            write(6,'(''QMMM: than 2xnorbs of: '',i5)') qm2_struct%norbs
            write(6,'(''QMMM: Check qmcharge in qmmm namelist and rerun'')')
            write(6,'(''QMMM: the calculation.'')')
        end if                                                                                                    
        call mexit(6,1)
    end if
    !Now we know the number of electrons work out how many closed and open shells there are.
    if(qmmm_nml%spin==1 .OR. qmmm_nml%spin==3 .OR.qmmm_nml%spin==5)then
        !        Make sure we have an even number of electrons
        if((nelectrons/2)*2 /= nelectrons) then
            if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
                write(6,'(''QMMM: System specified with odd number of electrons ('',i5,'')'')') nelectrons
                write(6,'(''QMMM: but odd spin ('',i3,''). You most likely have the charge of'')') qmmm_nml%spin
                write(6,'(''QMMM: QM region (qmcharge) set incorrectly. Correct error and re-run calculation.'')')
            end if
            call mexit(6,1)
        end if
    else if (qmmm_nml%spin==2 .OR. qmmm_nml%spin==4 .OR. qmmm_nml%spin==6) then
        !        Make sure we have an odd number of electrons.
        !        Note spins other than 1 are not currently allowed because UHF is not implemented.
        if((nelectrons/2)*2 == nelectrons) then
            if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
                write(6,'(''QMMM: System specified with even number of electrons ('',i5,'')'')') nelectrons
                write(6,'(''QMMM: but even spin ('',i3,''). You most likely have the charge of'')') qmmm_nml%spin
                write(6,'(''QMMM: QM region (qmcharge) set incorrectly. Correct error and re-run calculation.'')')
            end if
            call mexit(6,1)
        end if
    end if

    if (qmmm_nml%spin == 1) then
        nopen = 0
    else if (qmmm_nml%spin == 2) then
        nopen = 1
    else if (qmmm_nml%spin == 3) then
        nopen = 2
    else if (qmmm_nml%spin == 4) then
        nopen = 3
    else if (qmmm_nml%spin == 5) then
        nopen = 4
    else if (qmmm_nml%spin == 6) then
        nopen = 5
    end if
    qm2_struct%nclosed = nelectrons/2
    if( nopen > 0 ) then
        qm2_struct%nclosed = qm2_struct%nclosed - nopen/2
        if ((qm2_struct%nclosed + nopen) > qm2_struct%norbs) then
            if (qmmm_mpi%commqmmm_master .and. qmmm_struct%qm_mm_first_call) then
                write(6,'(''QMMM: Number of doubly filled ('',i3,'') plus'')') qm2_struct%nclosed
                write(6,'(''QMMM: number of partly filled ('',i3,'') levels'')') nopen
                write(6,'(''QMMM: is greater than the total number of orbitals ('',i5,'').'')') qm2_struct%norbs
                write(6,'(''QMMM: Fix problem and re-run calculation.'')')
            end if
            call mexit(6,1)
        end if
    end if
    qm2_struct%nopenclosed=nopen+qm2_struct%nclosed

    !Allocate things that depend on norbs:

    allocate (qm2_params%pascal_tri1(qm2_struct%norbs), stat=ier )
    REQUIRE(ier == 0)
    allocate (qm2_params%pascal_tri2(qm2_struct%norbs), stat=ier )
    REQUIRE(ier == 0)

    qm2_struct%matsize = ishft(qm2_struct%norbs*(qm2_struct%norbs+1),-1) !ishift(x,-1) = integer divide by 2

    allocate ( qm2_struct%den_matrix(qm2_struct%matsize), stat=ier )
    REQUIRE ( ier == 0 )

    allocate ( qm2_struct%old_den_matrix(qm2_struct%matsize), stat=ier )
    REQUIRE ( ier == 0 )
    !zero the entire density matrix on the first call
    qm2_struct%den_matrix = 0.0d0; qm2_struct%old_den_matrix = 0.0d0;

    allocate ( qm2_struct%old2_density(qm2_struct%norbs), stat=ier )
    REQUIRE ( ier == 0 ) !Used by qm2_cnvg as workspace.


    if (qmmm_nml%fock_predict == 1) then
        !We are using Pulay et al. Fock matrix prediction.
        allocate ( qm2_struct%fock_mat_final1(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
        allocate ( qm2_struct%fock_mat_final2(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
        allocate ( qm2_struct%fock_mat_final3(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
        allocate ( qm2_struct%fock_mat_final4(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
    end if

    allocate ( qm2_struct%fock_matrix(qm2_struct%matsize), stat=ier )
    REQUIRE ( ier == 0 )
    allocate ( qm2_struct%hmatrix(qm2_struct%matsize), stat=ier )
    REQUIRE ( ier == 0 )
    
    !BTN 08/10/2017 allocate Fock matricies
    allocate ( qm2_struct%fock_matrix_dp(qmmm_struct%nquant_nlink*(qmmm_struct%nquant_nlink-1)/2, &
        MaxValenceOrbitals*(MaxValenceOrbitals*2+1)), stat=ier )
    REQUIRE ( ier == 0 )

    !+TJG 01/26/2010
    allocate ( qm2_struct%diis_fock(qm2_struct%matsize,qmmm_nml%ndiis_matrices), stat=ier )
    REQUIRE ( ier == 0 )
    qm2_struct%diis_fock = 0.d0
    ! this could be done in packed form because the error matrix is an antisymmetric matrix,
    ! but I'd rather go for clarity than for speed. -TJG
    allocate ( qm2_struct%diis_errmat(qm2_struct%norbs*qm2_struct%norbs,qmmm_nml%ndiis_matrices), stat=ier )
    REQUIRE ( ier == 0 )
    qm2_struct%diis_errmat = 0.d0
    allocate ( qm2_struct%diis_mat(qmmm_nml%ndiis_matrices+1,qmmm_nml%ndiis_matrices+1), stat=ier )
    REQUIRE ( ier == 0 )
    qm2_struct%diis_mat = 0.d0
    !-TJG 01/26/2010

#ifdef MPI
# ifndef USE_MPI_IN_PLACE
    !Allocate a temporary array for doing the reduce of P and F etc. Only needed if
    !we can't do things using MPI_IN_PLACE
      allocate ( qmmm_scratch%matsize_red_scratch(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
# endif
#endif
    allocate (qm2_struct%fock2_ptot2(MaxValenceOrbitals**2,qmmm_struct%nquant_nlink),stat=ier)
    REQUIRE(ier==0)

    !set up array of lower half triangle indicies (Pascal's triangle)
    do I=1,qm2_struct%norbs
        qm2_params%pascal_tri1(I)=ishft((I*(I-1)),-1)
        qm2_params%pascal_tri2(I)=qm2_params%pascal_tri1(I)+I
    end do

    !Fill the diagonal of the density matrix with the first guess:

    pdiag_guess1=dble(qmmm_nml%qmcharge)/(qm2_struct%norbs+1.D-10)
    do j=1, qm2_struct%norbs
        qm2_struct%den_matrix(qm2_params%pascal_tri2(j))=-pdiag_guess1
        qm2_struct%old_den_matrix(qm2_params%pascal_tri2(j))=-pdiag_guess1
    end do
      
    do i=1,qmmm_struct%nquant_nlink
        first_orb=qm2_params%orb_loc(1,i)
        last_orb=qm2_params%orb_loc(2,i)
         
        k=4  ! default number of orbitals to be populated
 
        ! S only
        if ((qm2_params%core_chg(i).le.2) .and. &
            (qm2_params%natomic_orbs(i)==1)) then
            k=1
        end if
 
        ! when the atom has d electrons
        if ((qm2_params%core_chg(i).gt.8) .and. &
            (qm2_params%natomic_orbs(i).gt.4)) then
            k= qm2_params%natomic_orbs(i)
        end if
 
        pdiag_guess2=dble(qm2_params%core_chg(i))/dble(k)
        do j=first_orb,first_orb+k-1
            qm2_struct%den_matrix(qm2_params%pascal_tri2(j))= pdiag_guess2 + &
                qm2_struct%den_matrix(qm2_params%pascal_tri2(j))
            qm2_struct%old_den_matrix(qm2_params%pascal_tri2(j))= pdiag_guess2 +&
                qm2_struct%old_den_matrix(qm2_params%pascal_tri2(j))
        end do
    end do

      
    !!            qxd_s, qxd_z0, qxd_zq, qxd_d0, qxd_dq, qxd_q0, qxd_qq, qxd_neff
    !----------------------------------------
    !
    ! OPNQ parameter loding
    !
    !--------------------------------------
         
    if (qmmm_opnq%useOPNQ) then   
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%qxd_supported(i) = qxd_supported(qmmm_struct%qm_type_id(i))
            if (qm2_params%qxd_supported(i)) then
                qm2_params%qxd_s(i) = qxd_s(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_z0(i) = qxd_z0(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_zq(i) = qxd_zq(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_d0(i) = qxd_d0(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_dq(i) = qxd_dq(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_q0(i) = qxd_q0(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_qq(i) = qxd_qq(qmmm_struct%qm_type_id(i))
                qm2_params%qxd_neff(i) = qxd_neff(qmmm_struct%qm_type_id(i))
            end if
        end do
    end if ! (qmmm_opnq%useOPNQ)
    !----------------------------------------
    !
    ! Now we fill up the data depending on the method we are using
    !
    !--------------------------------------
    !           MNDO PARAMS               *
    !--------------------------------------
    if (qmmm_nml%qmtheory%MNDO) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in MNDO
            if (.NOT. element_supported_mndo(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no MNDO parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate', &
                    'UNSUPPORTED ELEMENT', &
                    'QM MNDO NOT AVAILABLE FOR THIS ATOM')
            end if

          
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_mndo(iqm_atomic)*iostmp + upp_mndo(iqm_atomic)*ioptmp &
                + gss_mndo(iqm_atomic)*gssc + gsp_mndo(iqm_atomic)*gspc &
                + gpp_mndo(iqm_atomic)*gppc + gp2_mndo(iqm_atomic)*gp2c &
                + hsp_mndo(iqm_atomic)*hspc

            if (p_orb_exp_mndo(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_mndo(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_mndo(iqm_atomic)*p_orb_exp_mndo(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_mndo(iqm_atomic) + p_orb_exp_mndo(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_mndo(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_mndo(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_mndo(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                !AD
                dd1_temp = (HSP_mndo(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_mndo(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_mndo(iqm_atomic)-gp2_mndo(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if
      
            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_mndo(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_mndo(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_mndo(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_mndo(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_mndo(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_mndo(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_mndo(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_mndo(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_mndo(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_mndo(qmmm_struct%qm_type_id(i))
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_mndo(qmmm_struct%qm_type_id(i))+betas_mndo(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_mndo(qmmm_struct%qm_type_id(i))+betap_mndo(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_mndo(qmmm_struct%qm_type_id(i))+betap_mndo(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !       end  MNDO PARAMS              *
    !--------------------------------------

    !----------------------------------------
    !
    ! Now we fill up the data depending on the method we are using
    !
    !--------------------------------------
    !           MNDOD PARAMS               *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%MNDOD) then
      
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in MNDOD
            if (.NOT. element_supported_MNDOD(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no MNDOD parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate', &
                    'UNSUPPORTED ELEMENT', &
                    'QM MNDOD NOT AVAILABLE FOR THIS ATOM')
            end if

          
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_MNDOD(iqm_atomic)*iostmp + upp_MNDOD(iqm_atomic)*ioptmp &
                + gss_MNDOD(iqm_atomic)*gssc + gsp_MNDOD(iqm_atomic)*gspc &
                + gpp_MNDOD(iqm_atomic)*gppc + gp2_MNDOD(iqm_atomic)*gp2c &
                + hsp_MNDOD(iqm_atomic)*hspc



            if(p_orb_exp_MNDOD(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_MNDOD(iqm_atomic) .ne. 0.0d0) then

                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_MNDOD(iqm_atomic)*p_orb_exp_MNDOD(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_MNDOD(iqm_atomic) + p_orb_exp_MNDOD(iqm_atomic)

                qm2_params%multip_2c_elec_params(1,i)=((base_temp1**exponent_temp1) &
                    *(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    /((base_temp2**exponent_temp2) * sqrt(3.0d0))

                qm2_params%multip_2c_elec_params(2,i)=&
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_MNDOD(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_MNDOD(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_MNDOD(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                !AD
                dd1_temp = (HSP_MNDOD(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_MNDOD(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_MNDOD(iqm_atomic)-gp2_MNDOD(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if
      
            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_MNDOD(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_MNDOD(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_MNDOD(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_MNDOD(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_MNDOD(iqm_atomic)
     
            qm2_params%cc_exp_params(i) = alp_MNDOD(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_MNDOD(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_MNDOD(iqm_atomic)
            qm2_params%orb_elec_ke(3,i) = udd_MNDOD(iqm_atomic)
        end do
        
        ! move the zeta loading to here so that the derived dd and rho_0 can be calculated
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%d_orb_exp_by_type(i) = d_orb_exp_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%s_orb_exp_tail_by_type(i) = s_orb_exp_tail_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_tail_by_type(i) = p_orb_exp_tail_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%d_orb_exp_tail_by_type(i) = d_orb_exp_tail_MNDOD(qmmm_struct%qm_type_id(i))
           
            qm2_params%gss(i) = gss_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%hsp(i) = hsp_MNDOD(qmmm_struct%qm_type_id(i))
            qm2_params%hpp(i) = (gpp_MNDOD(qmmm_struct%qm_type_id(i))-gp2_MNDOD(qmmm_struct%qm_type_id(i)))*half
           
            DD=0.0D0
            PO=0.0D0
            call GetDDAndPho(qm2_params,qmmm_struct,i, DD, PO)

            qm2_params%dd(1:6,i)=DD
            qm2_params%po(1:9,i)=PO
 
            temp=rho_core_mndod(qmmm_struct%qm_type_id(i))
            if (abs(temp) > 1.0d-5 ) then
                qm2_params%po(9,i)=temp
            end if
                              
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_MNDOD(qmmm_struct%qm_type_id(i))+betas_MNDOD(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_MNDOD(qmmm_struct%qm_type_id(i))+betap_MNDOD(qmmm_struct%qm_type_id(j))
                qm2_params%betasad(i,j) = betas_MNDOD(qmmm_struct%qm_type_id(i))+betad_MNDOD(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_MNDOD(qmmm_struct%qm_type_id(i))+betap_MNDOD(qmmm_struct%qm_type_id(j))
                qm2_params%betapad(i,j) = betap_MNDOD(qmmm_struct%qm_type_id(i))+betad_MNDOD(qmmm_struct%qm_type_id(j))
                qm2_params%betadad(i,j) = betad_MNDOD(qmmm_struct%qm_type_id(i))+betad_MNDOD(qmmm_struct%qm_type_id(j))
            end do
        end do
   
    !--------------------------------------
    !       end  MNDOD PARAMS              *
    !--------------------------------------

    !--------------------------------------
    !           AM1 PARAMS                *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%AM1) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in AM1
            if (.NOT. element_supported_am1(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no AM1 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate', &
                    'UNSUPPORTED ELEMENT', &
                    'QM AM1 NOT AVAILABLE FOR THIS ATOM')
            end if

            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_am1(iqm_atomic)*iostmp + upp_am1(iqm_atomic)*ioptmp &
                + gss_am1(iqm_atomic)*gssc + gsp_am1(iqm_atomic)*gspc &
                + gpp_am1(iqm_atomic)*gppc + gp2_am1(iqm_atomic)*gp2c &
                + hsp_am1(iqm_atomic)*hspc

            if (p_orb_exp_am1(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_am1(iqm_atomic) .ne. 0.0d0) then

                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_am1(iqm_atomic)*p_orb_exp_am1(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_am1(iqm_atomic) + p_orb_exp_am1(iqm_atomic)

                qm2_params%multip_2c_elec_params(1,i)=((base_temp1**exponent_temp1) &
                    *(2.0d0*nsshell(iqm_atomic)+1.0d0)) &
                    /((base_temp2**exponent_temp2) * sqrt(3.0d0))

                qm2_params%multip_2c_elec_params(2,i)= &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_am1(iqm_atomic)

                ! expression above sometimes (e.g., for hydrogen) gives Inf
                ! because p_orb_exp_am1 vanishes !!!???
                ! to correct for that: kav
                if(p_orb_exp_am1(iqm_atomic)==0.d0) then
                    qm2_params%multip_2c_elec_params(2,i)=0.d0
                end if
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if

            if (GSS_am1(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_am1(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_am1(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_am1(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_am1(iqm_atomic)-gp2_am1(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            ! Next do the electronic energy - add it to the total heat of formation energy.
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_am1(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_am1(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_am1(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_am1(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_am1(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_am1(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_am1(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_am1(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            ! get the Slater orbital expansion coefficients
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_am1(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_am1(qmmm_struct%qm_type_id(i))
            qm2_params%NUM_FN(i) = NUM_FN_am1(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_am1(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_am1(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_am1(j,qmmm_struct%qm_type_id(i))
            end do
        end do

        do i=1,qmmm_struct%qm_ntypes
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_am1(qmmm_struct%qm_type_id(i))+betas_am1(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_am1(qmmm_struct%qm_type_id(i))+betap_am1(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_am1(qmmm_struct%qm_type_id(i))+betap_am1(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !        end  AM1 PARAMS              *
    !--------------------------------------

    !--------------------------------------
    !           AM1D PARAMS               *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%AM1D) then
      
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in AM1D
            if (.NOT. element_supported_AM1D(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no AM1D parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM AM1D NOT AVAILABLE FOR THIS ATOM')
            end if

          
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_AM1D(iqm_atomic)*iostmp + upp_AM1D(iqm_atomic)*ioptmp &
                + gss_AM1D(iqm_atomic)*gssc + gsp_AM1D(iqm_atomic)*gspc &
                + gpp_AM1D(iqm_atomic)*gppc + gp2_AM1D(iqm_atomic)*gp2c &
                + hsp_AM1D(iqm_atomic)*hspc



            if (p_orb_exp_AM1D(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_AM1D(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_AM1D(iqm_atomic)*p_orb_exp_AM1D(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_AM1D(iqm_atomic) + p_orb_exp_AM1D(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_AM1D(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_AM1D(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_AM1D(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                !AD
                dd1_temp = (HSP_AM1D(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_AM1D(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_AM1D(iqm_atomic)-gp2_AM1D(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if
      
            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_AM1D(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_AM1D(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_AM1D(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_AM1D(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_AM1D(iqm_atomic)
     
            qm2_params%cc_exp_params(i) = alp_AM1D(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_AM1D(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_AM1D(iqm_atomic)
            qm2_params%orb_elec_ke(3,i) = udd_AM1D(iqm_atomic)
        end do
        
        ! move the zeta loading to here so that the derived dd and rho_0 can be calculated
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%d_orb_exp_by_type(i) = d_orb_exp_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%s_orb_exp_tail_by_type(i) = s_orb_exp_tail_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_tail_by_type(i) = p_orb_exp_tail_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%d_orb_exp_tail_by_type(i) = d_orb_exp_tail_AM1D(qmmm_struct%qm_type_id(i))
           
            qm2_params%gss(i) = gss_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%hsp(i) = hsp_AM1D(qmmm_struct%qm_type_id(i))
            qm2_params%hpp(i) = (gpp_AM1D(qmmm_struct%qm_type_id(i))-gp2_AM1D(qmmm_struct%qm_type_id(i)))*half
           
            qm2_params%GNN(i) = GNN_AM1D(qmmm_struct%qm_type_id(i))
           
            DD=0.0D0
            PO=0.0D0
            call GetDDAndPho(qm2_params,qmmm_struct,i, DD, PO)
           
            qm2_params%dd(1:6,i)=DD
            qm2_params%po(1:9,i)=PO
           
            temp=rho_core_am1d(qmmm_struct%qm_type_id(i))
            if (abs(temp) > 1.0d-5 ) then
                qm2_params%po(9,i)=temp
            end if
                    
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_AM1D(qmmm_struct%qm_type_id(i))+betas_AM1D(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_AM1D(qmmm_struct%qm_type_id(i))+betap_AM1D(qmmm_struct%qm_type_id(j))
                qm2_params%betasad(i,j) = betas_AM1D(qmmm_struct%qm_type_id(i))+betad_AM1D(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_AM1D(qmmm_struct%qm_type_id(i))+betap_AM1D(qmmm_struct%qm_type_id(j))
                qm2_params%betapad(i,j) = betap_AM1D(qmmm_struct%qm_type_id(i))+betad_AM1D(qmmm_struct%qm_type_id(j))
                qm2_params%betadad(i,j) = betad_AM1D(qmmm_struct%qm_type_id(i))+betad_AM1D(qmmm_struct%qm_type_id(j))
            end do
          
            qm2_params%NUM_FN(i) = NUM_FN_am1d(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_am1d(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_am1d(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_am1d(j,qmmm_struct%qm_type_id(i))
            end do

        end do
        
    !--------------------------------------
    !       end  AM1D PARAMS              *
    !--------------------------------------

    !-----------------------------------------------
    !     PM3, PM3CARB1, PM3ZNB and PM3MAIS PARAMS *
    !-----------------------------------------------
    else if (qmmm_nml%qmtheory%PM3 .OR. qmmm_nml%qmtheory%PM3CARB1  &
        .OR. qmmm_nml%qmtheory%PM3ZNB .OR. qmmm_nml%qmtheory%PM3MAIS) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in PM3
            if (.NOT. element_supported_pm3(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no PM3 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PM3 NOT AVAILABLE FOR THIS ATOM')
            end if
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_pm3(iqm_atomic)*iostmp + upp_pm3(iqm_atomic)*ioptmp &
                + gss_pm3(iqm_atomic)*gssc + gsp_pm3(iqm_atomic)*gspc &
                + gpp_pm3(iqm_atomic)*gppc + gp2_pm3(iqm_atomic)*gp2c &
                + hsp_pm3(iqm_atomic)*hspc

            if (p_orb_exp_pm3(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_pm3(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_pm3(iqm_atomic)*p_orb_exp_pm3(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_pm3(iqm_atomic) + p_orb_exp_pm3(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pm3(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_pm3(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pm3(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_pm3(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_pm3(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_pm3(iqm_atomic)-gp2_pm3(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pm3(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_pm3(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pm3(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pm3(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pm3(iqm_atomic)
            if (qmmm_nml%qmtheory%PM3CARB1 .AND. (iqm_atomic==1 .OR. iqm_atomic==8)) then
                !Load the PM3CARB1 versions of O and H params in place of the default PM3 params
                elec_eng = uss_pm3carb1(iqm_atomic)*iostmp + upp_pm3carb1(iqm_atomic)*ioptmp &
                    + gss_pm3(iqm_atomic)*gssc + gsp_pm3(iqm_atomic)*gspc &
                    + gpp_pm3(iqm_atomic)*gppc + gp2_pm3(iqm_atomic)*gp2c &
                    + hsp_pm3(iqm_atomic)*hspc
                qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
                qm2_params%cc_exp_params(i) = alp_pm3carb1(iqm_atomic)
                qm2_params%orb_elec_ke(1,i) = uss_pm3carb1(iqm_atomic)
                qm2_params%orb_elec_ke(2,i) = upp_pm3carb1(iqm_atomic)
            else
                qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
                qm2_params%cc_exp_params(i) = alp_pm3(iqm_atomic)
                qm2_params%orb_elec_ke(1,i) = uss_pm3(iqm_atomic)
                qm2_params%orb_elec_ke(2,i) = upp_pm3(iqm_atomic)
            end if
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            ! get the Slater orbital expansion coefficients
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pm3(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pm3(qmmm_struct%qm_type_id(i))
            qm2_params%NUM_FN(i) = NUM_FN_pm3(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_pm3(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_pm3(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_pm3(j,qmmm_struct%qm_type_id(i))
            end do
        end do

        !RCW: PM3CARB1 update - we need to make sure we use the correct parameters here for hydrogen and oxygen atoms
        !Simple solution is to replace the data in the betas_pm3 and betap_pm3 arrays with the PM3CARB1 values
        !This avoids having to use complex if statements in the loop below.
        if (qmmm_nml%qmtheory%PM3CARB1) then
            !Replace PM3 betas and betap params for O and H with PM3CARB1 values BEFORE they are copied
            !into the working array.
            betas_pm3(1) = betas_pm3carb1(1)
            betas_pm3(8) = betas_pm3carb1(8)
            betap_pm3(1) = betap_pm3carb1(1)
            betap_pm3(8) = betap_pm3carb1(8)
        end if
        do i=1,qmmm_struct%qm_ntypes
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_pm3(qmmm_struct%qm_type_id(i))+betas_pm3(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_pm3(qmmm_struct%qm_type_id(i))+betap_pm3(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_pm3(qmmm_struct%qm_type_id(i))+betap_pm3(qmmm_struct%qm_type_id(j))
            end do
        end do

        ! PM3/MM*
        if (qmmm_nml%qmmm_int == 3) then
            ! Should ultimately be based on types here to save memory.
            !do i = 1,qmmm_struct%nquant_nlink
            do i = 1,qmmm_struct%qm_ntypes
                iqm_atomic=qmmm_struct%qm_type_id(i)
                ! Check parameter availability for PM3/MM*
                ! Current params are for QM H, C, N and O only.
                if (.NOT. element_supported_pm3mmx(iqm_atomic)) then
                    write(6,'("QMMM: Atom with atomic number ",i4,".")') iqm_atomic
                    write(6,'("QMMM: There are no PM3/MM* parameters for this element. Sorry.")')
                    call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PM3/MM* NOT AVAILABLE FOR THIS ATOM')
                end if
          
                ! Load params - Should ultimately be done on atom types.
                qm2_params%scale_factor1_pm3mmx(1,i) = scale_f1_pm3mmx(1,iqm_atomic)
                qm2_params%scale_factor2_pm3mmx(1,i) = scale_f2_pm3mmx(1,iqm_atomic)

                qm2_params%scale_factor1_pm3mmx(2,i) = scale_f1_pm3mmx(2,iqm_atomic)
                qm2_params%scale_factor2_pm3mmx(2,i) = scale_f2_pm3mmx(2,iqm_atomic)

                qm2_params%rho_pm3mmx(i) = 0.0D0

            end do
        else if (qmmm_nml%qmmm_int == 4) then
            do i = 1,qmmm_struct%qm_ntypes
                iqm_atomic=qmmm_struct%qm_type_id(i)
                ! Check parameter availability for PM3/MM* 2nd version
                ! Current params are for QM H, C, N, O and S only.
                if (.NOT. element_supported_pm3mmx2(iqm_atomic)) then
                    write(6,'("QMMM: Atom with atomic number ",i4,".")') iqm_atomic
                    write(6,'("QMMM: There are no PM3/MMX2 parameters for this element. Sorry.")')
                    call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM PM3/MMX2 NOT AVAILABLE FOR THIS ATOM')
                end if
          
                ! Load params - Should ultimately be done on atom types.
                qm2_params%scale_factor1_pm3mmx(1,i) = scale_f1_pm3mmx2(1,iqm_atomic)
                qm2_params%scale_factor2_pm3mmx(1,i) = scale_f2_pm3mmx2(1,iqm_atomic)

                qm2_params%scale_factor1_pm3mmx(2,i) = scale_f1_pm3mmx2(2,iqm_atomic)
                qm2_params%scale_factor2_pm3mmx(2,i) = scale_f2_pm3mmx2(2,iqm_atomic)
                qm2_params%rho_pm3mmx(i) = rho_pm3mmx2(iqm_atomic)
            end do
        end if

        !--------------------------------------
        !        PM3MAIS PARAMS               *
        !--------------------------------------
        if (qmmm_nml%qmtheory%PM3MAIS) then
            do i = 1, qmmm_struct%qm_ntypes
                do j = 1, i
                    iat = qmmm_struct%qm_type_id(i)
                    jat = qmmm_struct%qm_type_id(j)
                    if ( .not. element_supported_pm3mais(iat) ) then
                        write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                        write(6,'("QMMM: There are no PM3-MAIS parameters for atomic number ",i4,". Sorry.")') iat
                        call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT',&
                            'QM PM3-MAIS NOT AVAILABLE FOR THIS ATOM')
                    end if
                    ! --------------------------------
                    ! PM3MAIS pairwise core core terms
                    ! --------------------------------
                    if (iat < jat) then
                        itmp = iat
                        iat = jat
                        jat = itmp
                    end if
                    do k = 1, 3
                        qm2_params%pm3mais_alpab(i,j,k) = alpab_pm3mais(iat,jat,k)
                        qm2_params%pm3mais_betab(i,j,k) = betab_pm3mais(iat,jat,k)
                        qm2_params%pm3mais_gamab(i,j,k) = gamab_pm3mais(iat,jat,k)
                        qm2_params%pm3mais_alpab(j,i,k) = qm2_params%pm3mais_alpab(i,j,k)
                        qm2_params%pm3mais_betab(j,i,k) = qm2_params%pm3mais_betab(i,j,k)
                        qm2_params%pm3mais_gamab(j,i,k) = qm2_params%pm3mais_gamab(i,j,k)
                    end do
                end do
            end do
        end if
    !--------------------------------------
    !    END PM3MAIS PARAMS               *
    !--------------------------------------


    !------------------------------------------------
    ! END PM3, PM3CARB1, PM3ZnB and PM3-MAIS PARAMS *
    !------------------------------------------------

    !--------------------------------------
    !      PM6 PARAMS                     *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%PM6) then

        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)

            ! Check that parameters exist for this element in PM6
            if (.NOT. element_supported_pm6(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no PM6 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PM6 NOT AVAILABLE FOR THIS ATOM')
            end if
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
            !            + HSP*HSPC
            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_pm6(iqm_atomic)*iostmp + upp_pm6(iqm_atomic)*ioptmp &
                + gss_pm6(iqm_atomic)*gssc + gsp_pm6(iqm_atomic)*gspc &
                + gpp_pm6(iqm_atomic)*gppc + gp2_pm6(iqm_atomic)*gp2c &
                + hsp_pm6(iqm_atomic)*hspc

            if ( abs(EISOL_pm6(iqm_atomic)) > 1.0d-09 ) then
                elec_eng = EISOL_pm6(iqm_atomic)
            end if

            if (p_orb_exp_pm6(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_pm6(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_pm6(iqm_atomic)*p_orb_exp_pm6(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_pm6(iqm_atomic) + p_orb_exp_pm6(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pm6(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_pm6(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pm6(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_pm6(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_pm6(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_pm6(iqm_atomic)-gp2_pm6(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pm6(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_pm6(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pm6(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pm6(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pm6(iqm_atomic)
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%cc_exp_params(i) = alp_pm6(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pm6(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pm6(iqm_atomic)
            qm2_params%orb_elec_ke(3,i) = udd_pm6(iqm_atomic)

        end do

        do i = 1, qmmm_struct%qm_ntypes
            iat = qmmm_struct%qm_type_id(i)
            ! Load pre-computed Slater-Condon parameters
            qm2_params%F0SD(i) = F0SD_pm6(iat)
            qm2_params%G2SD(i) = G2SD_pm6(iat)
        end do
       
        ! Precompute some parameters to save time later
        do i=1,qmmm_struct%qm_ntypes

            iat = qmmm_struct%qm_type_id(i)

            qm2_params%NUM_FN(i) = NUM_FN_pm6(iat)

            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pm6(iat)
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pm6(iat)
            qm2_params%d_orb_exp_by_type(i) = d_orb_exp_pm6(iat)
            qm2_params%s_orb_exp_tail_by_type(i) = s_orb_exp_tail_pm6(iat)
            qm2_params%p_orb_exp_tail_by_type(i) = p_orb_exp_tail_pm6(iat)
            qm2_params%d_orb_exp_tail_by_type(i) = d_orb_exp_tail_pm6(iat)
           
            qm2_params%gss(i) = gss_pm6(iat)
            qm2_params%hsp(i) = hsp_pm6(iat)
            qm2_params%hpp(i) = (gpp_pm6(iat)-gp2_pm6(iat))*half
           
            DD=0.0D0
            PO=0.0D0
            call GetDDAndPho(qm2_params,qmmm_struct,i, DD, PO)

            qm2_params%dd(1:6,i)=DD
            qm2_params%po(1:9,i)=PO

            temp=rho_core_pm6(qmmm_struct%qm_type_id(i))
            if (abs(temp) > 1.0d-5 ) then
                qm2_params%po(9,i)=temp
            end if
                              
            do j=1,4
                qm2_params%FN1(j,i) = FN1_pm6(j,iat)
                qm2_params%FN2(j,i) = FN2_pm6(j,iat)
                qm2_params%FN3(j,i) = FN3_pm6(j,iat)
            end do
        end do

        do i=1,qmmm_struct%qm_ntypes
            iat = qmmm_struct%qm_type_id(i)
            do j = 1,qmmm_struct%qm_ntypes
                jat = qmmm_struct%qm_type_id(j)
                qm2_params%betasas(i,j) = betas_pm6(iat)+betas_pm6(jat)
                qm2_params%betasap(i,j) = betas_pm6(iat)+betap_pm6(jat)
                qm2_params%betasad(i,j) = betas_pm6(iat)+betad_pm6(jat)
                qm2_params%betapap(i,j) = betap_pm6(iat)+betap_pm6(jat)
                qm2_params%betapad(i,j) = betap_pm6(iat)+betad_pm6(jat)
                qm2_params%betadad(i,j) = betad_pm6(iat)+betad_pm6(jat)
            end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
            do j = 1, i
                iat = qmmm_struct%qm_type_id(i)
                jat = qmmm_struct%qm_type_id(j)
                ! ----------------------------
                ! PM6 pairwise core core terms
                ! ----------------------------
                if (iat < jat) then
                    itmp = iat
                    iat = jat
                    jat = itmp
                end if
                qm2_params%pm6_alpab(i,j) = alpab_pm6(iat,jat)
                qm2_params%pm6_xab(i,j)   = xab_pm6(iat,jat)
                qm2_params%pm6_alpab(j,i) = qm2_params%pm6_alpab(i,j)
                qm2_params%pm6_xab(j,i)   = qm2_params%pm6_xab(i,j)
            end do
        end do

    !--------------------------------------
    !    END PM6 PARAMS                   *
    !--------------------------------------

    !--------------------------------------
    !         PDDG/PM3 PARAMS             *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%PDDGPM3) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in PM3
            if (.NOT. element_supported_pddgpm3(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no PDDG-PM3 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PDDG-PM3 NOT AVAILABLE FOR THIS ATOM')
            end if
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            ! In PM3/PDDG the electronic energy is treated as an independent parameter
            ! and so we have to have it explicitly defined in the parameter file and
            ! cannot just calculate it.

            if (p_orb_exp_pddgpm3(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_pddgpm3(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_pddgpm3(iqm_atomic)*p_orb_exp_pddgpm3(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_pddgpm3(iqm_atomic) + p_orb_exp_pddgpm3(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pddgpm3(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_pddgpm3(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pddgpm3(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_pddgpm3(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_pddgpm3(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_pddgpm3(iqm_atomic)-gp2_pddgpm3(iqm_atomic))
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng_pddgpm3(iqm_atomic)*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pddgpm3(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_pddgpm3(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pddgpm3(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pddgpm3(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pddgpm3(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_pddgpm3(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pddgpm3(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pddgpm3(iqm_atomic)
            qm2_params%pddge1(i) = pddge1_pm3(iqm_atomic)
            qm2_params%pddge2(i) = pddge2_pm3(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pddgpm3(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pddgpm3(qmmm_struct%qm_type_id(i))
            qm2_params%NUM_FN(i) = NUM_FN_pddgpm3(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_pddgpm3(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_pddgpm3(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_pddgpm3(j,qmmm_struct%qm_type_id(i))
            end do
        end do

        do i=1,qmmm_struct%qm_ntypes
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_pddgpm3(qmmm_struct%qm_type_id(i)) &
                    +betas_pddgpm3(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_pddgpm3(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgpm3(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_pddgpm3(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgpm3(qmmm_struct%qm_type_id(j))
                pddg_zaf = dble(core_chg(qmmm_struct%qm_type_id(i))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                pddg_zbf = dble(core_chg(qmmm_struct%qm_type_id(j))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                qm2_params%pddg_term1(i,j) = &
                    pddg_zaf*pddgc1_pm3(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_pm3(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term2(i,j) = &
                    pddg_zaf*pddgc1_pm3(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_pm3(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term3(i,j) = &
                    pddg_zaf*pddgc2_pm3(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_pm3(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term4(i,j) = &
                    pddg_zaf*pddgc2_pm3(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_pm3(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !          END PDDG/PM3 PARAMS        *
    !--------------------------------------
    !--------------------------------------
    !     PDDG/PM3 PARAMS 2008 Variant    *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%PDDGPM3_08) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in PM3
            if (.NOT. element_supported_pddgpm3_08(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no PDDG-PM3_08 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PDDG-PM3_08 NOT AVAILABLE FOR THIS ATOM')
            end if
            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            ! In PM3/PDDG_08 the electronic energy is treated as an independent parameter
            ! and so we have to have it explicitly defined in the parameter file and
            ! cannot just calculate it.

            if (p_orb_exp_pddgpm3_08(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_pddgpm3_08(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_pddgpm3_08(iqm_atomic)*p_orb_exp_pddgpm3_08(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_pddgpm3_08(iqm_atomic) + p_orb_exp_pddgpm3_08(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pddgpm3_08(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_pddgpm3_08(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pddgpm3_08(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_pddgpm3_08(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_pddgpm3_08(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_pddgpm3_08(iqm_atomic)-gp2_pddgpm3_08(iqm_atomic))
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng_pddgpm3_08(iqm_atomic)*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pddgpm3_08(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_pddgpm3_08(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pddgpm3_08(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pddgpm3_08(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pddgpm3_08(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_pddgpm3_08(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pddgpm3_08(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pddgpm3_08(iqm_atomic)
            qm2_params%pddge1(i) = pddge1_pm3_08(iqm_atomic)
            qm2_params%pddge2(i) = pddge2_pm3_08(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pddgpm3_08(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pddgpm3_08(qmmm_struct%qm_type_id(i))
            qm2_params%NUM_FN(i) = NUM_FN_pddgpm3_08(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_pddgpm3_08(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_pddgpm3_08(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_pddgpm3_08(j,qmmm_struct%qm_type_id(i))
            end do
        end do

        do i=1,qmmm_struct%qm_ntypes
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_pddgpm3_08(qmmm_struct%qm_type_id(i)) &
                    +betas_pddgpm3_08(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_pddgpm3_08(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgpm3_08(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_pddgpm3_08(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgpm3_08(qmmm_struct%qm_type_id(j))
                pddg_zaf = dble(core_chg(qmmm_struct%qm_type_id(i))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                pddg_zbf = dble(core_chg(qmmm_struct%qm_type_id(j))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                qm2_params%pddg_term1(i,j) = &
                    pddg_zaf*pddgc1_pm3_08(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_pm3_08(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term2(i,j) = &
                    pddg_zaf*pddgc1_pm3_08(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_pm3_08(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term3(i,j) = &
                    pddg_zaf*pddgc2_pm3_08(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_pm3_08(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term4(i,j) = &
                    pddg_zaf*pddgc2_pm3_08(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_pm3_08(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !   END PDDG/PM3 PARAMS 2008 Variant  *
    !--------------------------------------
    !--------------------------------------
    !           PDDG/MNDO PARAMS          *
    !--------------------------------------
    elseif (qmmm_nml%qmtheory%PDDGMNDO) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in MNDO
            if (.NOT. element_supported_pddgmndo(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no PDDG-MNDO parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM PDDG-MNDO NOT AVAILABLE FOR THIS ATOM')
            end if

            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            !1) Electronic Energy (EISOL)
            ! In PM3/PDDG the electronic energy is treated as an independent parameter
            ! and so we have to have it explicitly defined in the parameter file and
            ! cannot just calculate it.

            if (p_orb_exp_pddgmndo(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_pddgmndo(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_pddgmndo(iqm_atomic)*p_orb_exp_pddgmndo(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_pddgmndo(iqm_atomic) + p_orb_exp_pddgmndo(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pddgmndo(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_pddgmndo(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pddgmndo(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_pddgmndo(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_pddgmndo(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_pddgmndo(iqm_atomic)-gp2_pddgmndo(iqm_atomic))
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            ! Next do the electronic energy - add it to the total heat of formation energy.
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng_pddgmndo(iqm_atomic)*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pddgmndo(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_pddgmndo(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pddgmndo(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pddgmndo(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pddgmndo(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_pddgmndo(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pddgmndo(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pddgmndo(iqm_atomic)
            qm2_params%pddge1(i) = pddge1_mndo(iqm_atomic)
            qm2_params%pddge2(i) = pddge2_mndo(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pddgmndo(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pddgmndo(qmmm_struct%qm_type_id(i))
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = &
                    betas_pddgmndo(qmmm_struct%qm_type_id(i)) &
                    +betas_pddgmndo(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = &
                    betas_pddgmndo(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgmndo(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = &
                    betap_pddgmndo(qmmm_struct%qm_type_id(i)) &
                    +betap_pddgmndo(qmmm_struct%qm_type_id(j))
                pddg_zaf = dble(core_chg(qmmm_struct%qm_type_id(i))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                pddg_zbf = dble(core_chg(qmmm_struct%qm_type_id(j))) &
                    /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
                    dble(core_chg(qmmm_struct%qm_type_id(j))))
                qm2_params%pddg_term1(i,j) = &
                    pddg_zaf*pddgc1_mndo(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_mndo(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term2(i,j) = &
                    pddg_zaf*pddgc1_mndo(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_mndo(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term3(i,j) = &
                    pddg_zaf*pddgc2_mndo(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc1_mndo(qmmm_struct%qm_type_id(j))
                qm2_params%pddg_term4(i,j) = &
                    pddg_zaf*pddgc2_mndo(qmmm_struct%qm_type_id(i)) &
                    +pddg_zbf*pddgc2_mndo(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !     END PDDG/MNDO PARAMS            *
    !--------------------------------------

    !--------------------------------------
    !           RM1 PARAMS                *
    !--------------------------------------
    else if (qmmm_nml%qmtheory%RM1) then
        do i = 1,qmmm_struct%nquant_nlink
            iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
            ! Check that parameters exist for this element in RM1
            if (.NOT. element_supported_rm1(iqm_atomic)) then
                write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
                write(6,'("QMMM: There are no RM1 parameters for this element. Sorry.")')
                call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED ELEMENT','QM RM1 NOT AVAILABLE FOR THIS ATOM')
            end if

            !----------------------------------------
            ! Calculate parameters that are actually
            ! derived from other parameters.
            !----------------------------------------

            iostmp = ios(iqm_atomic)
            ioptmp = iop(iqm_atomic)
            gssc = dble(max(iostmp-1,0))
            gspc = dble(iostmp * ioptmp)
            gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
            hspc = dble(-ioptmp)
            elec_eng = uss_rm1(iqm_atomic)*iostmp + upp_rm1(iqm_atomic)*ioptmp &
                + gss_rm1(iqm_atomic)*gssc + gsp_rm1(iqm_atomic)*gspc &
                + gpp_rm1(iqm_atomic)*gppc + gp2_rm1(iqm_atomic)*gp2c &
                + hsp_rm1(iqm_atomic)*hspc

            if (p_orb_exp_rm1(iqm_atomic) .ne. 0.0d0 .or. &
                s_orb_exp_rm1(iqm_atomic) .ne. 0.0d0) then
                exponent_temp1 = nsshell(iqm_atomic)+0.5d0
                base_temp1 = 4.0d0*s_orb_exp_rm1(iqm_atomic)*p_orb_exp_rm1(iqm_atomic)
                exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
                base_temp2 = s_orb_exp_rm1(iqm_atomic) + p_orb_exp_rm1(iqm_atomic)
                qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                    / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
                qm2_params%multip_2c_elec_params(2,i) = &
                    sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_rm1(iqm_atomic)
            else
                qm2_params%multip_2c_elec_params(1,i)= 0.0d0
                qm2_params%multip_2c_elec_params(2,i)= 0.0d0
            end if
            if (GSS_rm1(iqm_atomic) .ne. 0.0d0 ) then
                qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_rm1(iqm_atomic) !AM
            else
                qm2_params%multip_2c_elec_params(3,i) = 0.0d0
            end if
            ! Calculation of AD and AQ
            if (iqm_atomic == 1) then
                qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
                qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
            else
                dd1_temp = (HSP_rm1(iqm_atomic) &
                    /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
                    hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(HSP_rm1(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
                !end AD
                !AQ
                hpp = 0.5D0*(gpp_rm1(iqm_atomic)-gp2_rm1(iqm_atomic))
                hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                     !match mopac results for Chlorine and potentially other elements.
                dd1_temp = (16.0d0*hpp &
                    /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
                dd2_temp = dd1_temp + 0.04d0
                do j = 1, 5
                    dd_diff = dd2_temp - dd1_temp
                    hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd1_temp**2)
                    hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                        + 1.0d0/dd2_temp**2)
                    if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
                    dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                        / (hsp2_temp - hsp1_temp)
                    dd1_temp = dd2_temp
                    dd2_temp = dd3_temp
                end do
                qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
              !end AQ
            end if

            !----------------------------------------
            ! End calculation of derived parameters.
            !----------------------------------------

            ! Next do the electronic energy - add it to the total heat of formation energy.
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_rm1(iqm_atomic)
            qm2_params%onec2elec_params(2,i) = GSP_rm1(iqm_atomic)
            qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_rm1(iqm_atomic)
            qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_rm1(iqm_atomic)
            qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_rm1(iqm_atomic)
            qm2_params%cc_exp_params(i) = alp_rm1(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_rm1(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_rm1(iqm_atomic)
        end do

        ! Precompute some parameters to save time later
        ! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
        ! size of the arrays. While this does not save us much memory it greatly increases
        ! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
            qm2_params%s_orb_exp_by_type(i) = s_orb_exp_rm1(qmmm_struct%qm_type_id(i))
            qm2_params%p_orb_exp_by_type(i) = p_orb_exp_rm1(qmmm_struct%qm_type_id(i))
            qm2_params%NUM_FN(i) = NUM_FN_rm1(qmmm_struct%qm_type_id(i))
            do j=1,4
                qm2_params%FN1(j,i) = FN1_rm1(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN2(j,i) = FN2_rm1(j,qmmm_struct%qm_type_id(i))
                qm2_params%FN3(j,i) = FN3_rm1(j,qmmm_struct%qm_type_id(i))
            end do
        end do

        do i=1,qmmm_struct%qm_ntypes
            do j = 1,qmmm_struct%qm_ntypes
                qm2_params%betasas(i,j) = betas_rm1(qmmm_struct%qm_type_id(i))+betas_rm1(qmmm_struct%qm_type_id(j))
                qm2_params%betasap(i,j) = betas_rm1(qmmm_struct%qm_type_id(i))+betap_rm1(qmmm_struct%qm_type_id(j))
                qm2_params%betapap(i,j) = betap_rm1(qmmm_struct%qm_type_id(i))+betap_rm1(qmmm_struct%qm_type_id(j))
            end do
        end do

    !--------------------------------------
    !        END  RM1 PARAMS              *
    !--------------------------------------

    else if (qmmm_nml%qmtheory%DFTB) then
    ! AWG: Load the DFTB parameters AFTER calling the info printing routine
    ! AWG: since qm2_dftb_load_params also prints information
    ! AWG: Otherwise the tests will not pass because the output gets mixed up
    ! AWG call qm2_dftb_load_params
    continue
else if (qmmm_nml%qmtheory%EXTERN) then
! AWG: We are using the external interface and don't need 99% of the stuff
! AWG: that is allocated / precomputed here. This needs to be cleaned up!
continue
else
    !UNKNOWN method - should never actually get this far but might as well call
    !sander bomb just in case.
    write (6,'("QMMM ERROR: Semiempirical method is not supported.")')
    call sander_bomb('qm2_load_params_and_allocate','UNSUPPORTED METHOD (qm_theory)', &
        'SELECTED LEVEL OF THEORY IS NOT AVAILABLE - PLEASE CHECK YOUR INPUT FILE')
end if

if (.not. qmmm_nml%qmtheory%DFTB) then
    ! ------------------------------------------------------
    ! Now see if user wants an MM peptide torsion correction
    ! ------------------------------------------------------
    qm2_struct%n_peptide_links = 0

    ! Finally setup the STO-6G orbital expansions and allocate the memory required.
    ! Setup the STO-6G orbital expansions and pre-calculate as many overlaps by type
    ! as we can and store these in memory. This will help a lot with speed in the
    ! energy and derivative code.
    call qm2_setup_orb_exp(qm2_params, qmmm_struct)

end if

if (qmmm_mpi%commqmmm_master) then
    ! ------------------------------
    ! Print information for this run
    ! ------------------------------
    ! This should not be called from here but we have to take apart and
    ! modularize the code first and/or adjust the DFTB test outputs
    ! because qm2_dftb_load_params (called below) prints stuff as well...
    if ( .not. qmmm_nml%qmtheory%EXTERN ) then
        call qm2_print_info(qmmm_nml,qm2_struct,qmmm_struct)
    end if
end if

if (qmmm_nml%qmtheory%DFTB) then
    call qm2_dftb_load_params(qmmm_mpi,qmmm_nml,qm2_struct, qmmm_struct)
end if

!In Parallel calculate the offset into the two electron array for each thread.
!This depends on the number of light-light, light-heavy and heavy-heavy interactions
!that this thread will do.

!Simulate what my loop would be and work out what my ending offset would be.
#ifdef MPI
      loop_count = 0
      do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
        jstart = qmmm_mpi%nquant_nlink_jrange(1,i)
        jend = qmmm_mpi%nquant_nlink_jrange(2,i)
        ia = qm2_params%orb_loc(1,i)
        ib = qm2_params%orb_loc(2,i)
        do j = jstart, jend
          ja = qm2_params%orb_loc(1,j)
          jb = qm2_params%orb_loc(2,j)
          if (ib /= ia .and. ja/=jb) then
      !Heavy - Heavy
            loop_count = loop_count+100
          elseif (ia /= ib) then
      !Light - Heavy
            loop_count = loop_count+10
          elseif (ja /= jb) then
      !Heavy - Light
            loop_count = loop_count+10
          else
      !Light - Light
            loop_count = loop_count+1
          endif
        end do
      end do

!At the end of this loop loop_count should be the starting value for the next cpu.
!However, we need to add the offset from all the other cpus to this.
!Each thread in turn passes it's total so far to the next thread.
!This total becomes that threads offset
      allocate(gather_array(qmmm_mpi%numthreads),stat=ier)
      REQUIRE(ier==0)
      call mpi_allgather(loop_count,1,mpi_integer,gather_array,1,mpi_integer,qmmm_mpi%commqmmm,ier)

!Now, our starting offset is the sum of all the previous thread's loop count
      qmmm_mpi%two_e_offset = 0
      do i = 1, qmmm_mpi%mytaskid
        qmmm_mpi%two_e_offset = qmmm_mpi%two_e_offset + gather_array(i)
      end do
      deallocate(gather_array,stat=ier)
      REQUIRE(ier==0)
#else
qmmm_mpi%two_e_offset = 0
#endif


! ------------------------------------------------
! Choose diagonalizer and allocate required memory
! ------------------------------------------------
if (.not. ( qmmm_nml%qmtheory%DFTB .or. qmmm_nml%qmtheory%EXTERN ) ) then
    call qm2_diagonalizer_setup(qm2_params,qmmm_mpi,qmmm_nml,qmmm_nml%diag_routine, qmmm_nml%allow_pseudo_diag, &
        qmmm_nml%verbosity, &
        qmmm_mpi%commqmmm_master, qmmm_mpi%commqmmm, &
        qm2_struct, qmmm_scratch)
end if

! ----------------------------------------------------------------
! Adjust whether numerical integral derivatives are required
! Check also whether we are running in parallel and quit if d
! orbitals are present since the d orbital code is not parallel
! ----------------------------------------------------------------
! Analytical integral derivatives are not available for d orbitals
! for MNDO type Hamiltonians
if ( (.not. qmmm_nml%qmtheory%DFTB) .and. (.not. qmmm_nml%qmtheory%EXTERN) ) then

    test = .false.
    do i=1,qmmm_struct%qm_ntypes
        iat = qmmm_struct%qm_type_id(i)
        if ( natomic_orbs(iat) > 4 ) then
            test = .true.
        end if
    end do

    if ( test .and. qmmm_nml%qmqm_analyt ) then
        write (6,'(/,a)') '| QMMM: *** Integral code information ***'
        write (6,'(a)')   '| QMMM: One or more elements are using d orbitals.'
        write (6,'(a)')   '| QMMM: Analytical derivatives for d orbitals are not supported.'
        write (6,'(a)')   '| QMMM: Forcing numerical integral derivatives.'
        qmmm_nml%qmqm_analyt = .false.
    end if

    if ( test .and. (qmmm_mpi%numthreads > 1) ) then
        call sander_bomb('sander()', &
            'Using d orbitals but the d orbital code is not parallelized.', &
            'Please run in serial.')
    end if

end if

end subroutine qm2_load_params_and_allocate
