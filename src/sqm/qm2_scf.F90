! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "def_time.h"
#include "assert.fh"
subroutine qm2_scf(qmmm_opnq, qm2_params, qmewald, qmmm_nml, qm_gb, qmmm_mpi, qmmm_scratch, &
	xlbomd_struct,cosmo_c_struct,qm2_struct, qm2ds, qmmm_struct, &
	fock_matrix, hmatrix, W, escf, den_matrix, scf_mchg, num_qmmm_calls)
    !---------------------------------------------------------------------
    ! This is the main SCF routine.
    ! Written by Ross Walker (TSRI, 2005)
    !
    ! This routine generates a self consistent field
    ! and returns the energy in the variable ESCF.
    !
    ! The main arrays used are:
    !     fock_matrix - Starts off containing the one-electron matrix
    !                   and is used to hold the Fock matrix. Units = Ev/Bohr
    !     hmatrix - Only ever contains the one electron matrix.
    !     W - Only ever contains the two electron matrix.
    !     P - Only ever contains the total density matrix of the
    !         current SCF step.
    !
    ! qm2_struct%matsize = size of packed triangle (NORBS*(NORBS+1)/2)
    ! qm2_struct%n2el = Number of 2 electron integrals
    !                    50*nheavy(nheavy-1)+10*nheavy*nlight+(nlight*(nlight-1))/2
    !
    ! References:
    ! PSEUDO DIAGONALISATION : "FAST SEMIEMPIRICAL CALCULATIONS",
    !                           STEWART. J.J.P., CSASZAR, P., PULAY, P.,
    !                           J. COMP. CHEM.,3, 227, (1982)
    !---------------------------------------------------------------------

    use qmmm_module, only : qm2_structure,  qm_ewald_structure,  &
        qm_gb_structure, qmmm_mpi_structure, qmmm_scratch_structure, qmmm_opnq_structure
    use constants, only : EV_TO_KCAL, zero, two
    use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type

    use qm2_iterator_mod, only : scf_iterator_value     ! this is an integer function
    use qm2_iterator_mod, only : diis_iterator_value    ! this is an integer function
    use qm2_iterator_mod, only : remaining_diis_tokens  ! this is an integer function

    use qm2_davidson_module
    use xlbomd_module, only : xlbomd_structure !K
    use qmmm_struct_module, only : qmmm_struct_type
    use cosmo_C, only : cosmo_C_structure

    implicit none
    type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq 
    type(qm_ewald_structure),intent(inout) :: qmewald
    type(qm_gb_structure),intent(inout) :: qm_gb
    type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
    type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
    type(qm2_params_type),intent(inout) :: qm2_params
    type(qmmm_nml_type),intent(inout) :: qmmm_nml
 
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(xlbomd_structure),intent(inout) :: xlbomd_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds


#ifdef MPI
   include 'mpif.h'
   integer :: ier
#endif

    !Passed in
    _REAL_, intent(inout) :: fock_matrix(qm2_struct%matsize)
    _REAL_, intent(in) :: hmatrix(qm2_struct%matsize)
    _REAL_, intent(in) :: W(qm2_struct%n2el)
    _REAL_, intent(inout) :: den_matrix(qm2_struct%matsize)
    _REAL_, intent(inout) :: escf
    _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)
    integer, intent(in) :: num_qmmm_calls

    !Local
    _REAL_ eold !SCF energy on previous step
    _REAL_ energy_diff, density_diff !Difference in energy and density from previous step.a
    _REAL_ small!, qm2_params%smallsum, qm2_params%abstol !Precision limits of this machine.
    _REAL_ smallest_energy_diff(2) !Smallest energy diff found so far (1) and density diff for this step(2)
    _REAL_ scf_energy !Computed in parts on different cpus and then all reduced. Only master returns
                      !this value in escf.
    integer sm_energy_diff_step_number
    integer scf_iteration !Number of scf iteration that have been done.
    integer i,m,n,o,p,q
    logical converged, first_iteration
    logical doing_pseudo_diag !Initially false, set to true if we are doing a pseudo diagonalisation.
    logical allow_pseudo_diag !Initially set to false. Set to true after at least 2 SCF iterations have been done 
                              !and allowed by the namelist option. Set back to full on last step to force a full
                              !diagonalisation.
    logical pseudo_converged !Flag used to indicate when doing pseudo diagonalisations has converged the SCF.
                                  !Once this is set to try the value of allow_pseudo_diag will not be set to true
                                  !again. This is used to force the last few SCF steps to be full diagonalisations.
    logical fock_predict_active !Set to true if we are attempting to predict the fock matrix on this QMMM call.

    character(len=1) :: fock_extrap_flag

    ! DIIS
    ! While DIIS will accelerate SCF convergence it will do so at a very high cost
    ! since Fock matrix builds and diagonalizations for semiempirical methods are
    ! cheaper than the overhead generated by the matrix multiplications required for DIIS
    ! Thus switch DIIS only on if we have difficulty in converging the SCF
    ! switch DIIS on after SCF step 100 if user turns on DIIS by setting ndiis_attempts > 0
    integer, parameter :: switch_on_diis = 100 ! don't set this to values lower than 3!!!
    ! switch DIIS on as a last resort if SCF did not converge until iteration 800
    integer, parameter :: switch_on_diis_last_resort = 800
    logical diis_is_off, errmat_is_converged, errmat_is_on
    logical dont_turn_off_diis
    _REAL_ :: errval
    _REAL_ ,allocatable, dimension(:,:) :: density_matrix_unpacked; 
    _REAL_ ,allocatable, dimension(:,:) :: fock_matrix_unpacked;

    ! current_scf_errval is a function
    _REAL_ :: current_scf_errval

    !qm2_Helect is a function
    _REAL_ qm2_HELECT

    !dlamch is a function
    _REAL_ dlamch

#ifdef MPI
    _REAL_ tmp_recv(2)
    logical :: conv_controls(4)
# ifndef USE_MPI_IN_PLACE
    _REAL_ tmp_send(2)
# endif
#endif

    integer lapack_info

    !Initialisation on every call
    converged = .false.
    first_iteration = .true.
    doing_pseudo_diag = .false.
    allow_pseudo_diag = .false.
    pseudo_converged = .false.
    eold = zero
    energy_diff=huge(energy_diff)
    density_diff = huge(density_diff)
    smallest_energy_diff(1) = huge(smallest_energy_diff(1))
    smallest_energy_diff(2) = huge(smallest_energy_diff(2))
    sm_energy_diff_step_number = 0

    ! DIIS IS TURNED OFF BY Dcosmo_c_struct%EFAULT
    ! ONLY TURN ON IF CONVERGENCE GOES NUTTY
    !    i = remaining_diis_tokens( qmmm_nml%ndiis_attempts )
    i = remaining_diis_tokens(qmmm_nml, 0 )
    dont_turn_off_diis = .FALSE.
    diis_is_off = qmmm_nml%ndiis_attempts < 1 .OR. qmmm_nml%ndiis_matrices < 2
    errmat_is_on = qmmm_nml%errconv < 1.0D-02
    errval = -1.0D0 ! initialize errval to unrealistic value

    fock_extrap_flag = " "


    if (qmmm_nml%verbosity > 2 .and. qmmm_mpi%commqmmm_master) then
        write(6,'("QMMM: ")')
        write(6,'("QMMM: SCF Convergence Information")')
        if (qmmm_nml%allow_pseudo_diag) write(6,'("QMMM: (*) = Pseudo Diagonalisation")')
        if ( remaining_diis_tokens(qmmm_nml) > 1 .AND. qmmm_nml%ndiis_matrices > 1 ) THEN
            write(6,'("QMMM: (D) = Used Fock extrapolation via DIIS")')
            write(6,'("QMMM: (M) = Used Fock extrapolation via 50/50 mixing")')
        END IF
        write(6,'("QMMM: Cycle      Energy       ")',ADVANCE='NO')
        write(6,'("   dE             dP              |FP-PF|")')
    end if

    if (qmmm_nml%fock_predict == 1 .and. num_qmmm_calls > 4)  then
        fock_predict_active = .true.
    else
        fock_predict_active = .false.
    end if


    if ( fock_predict_active ) then
        CALL qm2_diag(qm2_params,qmmm_mpi,qmmm_scratch, qmmm_nml, &
            & qm2_struct,qm2_struct%matsize, fock_matrix, &
            & allow_pseudo_diag, density_diff, doing_pseudo_diag, &
            & qm2_params%smallsum, qm2_params%abstol )
    end if


    ! MAIN SCF LOOP
    do_scf: do scf_iteration=1,abs(qmmm_nml%itrmax)
        ! set the current scf iterator to the current value
        i = scf_iterator_value(qmmm_nml, scf_iteration )

        ! construct a trial density
        if ((.NOT. first_iteration) .or. fock_predict_active) then
            CALL qm2_densmat(qm2_params, qmmm_mpi,qmmm_nml,qm2_struct, scf_iteration, qm2_struct%matsize, den_matrix, density_diff)
        end if ! if (.NOT. first_iteration)

         !Calculate the Mulliken charges for the current density matrix if we require
         !them on every SCF step. Save the results in the scf_mchg array
         !We also need to do this if we are doing QMewald with the image charges fixed but
         !it is our first MD step.
         !For the moment let all threads calculate mulliken charges as this is
         !probably quicker than doing a reduce.

        if (qm2_struct%calc_mchg_scf .or. qmewald%ewald_startup) then
            write(6,*)'Calculating Mulliken Charges'
            allocate(density_matrix_unpacked(qm2_struct%norbs,qm2_struct%norbs));
            call unpacking(qm2_struct%norbs,qm2_struct%den_matrix,density_matrix_unpacked,'s');
            do i=1,qmmm_struct%nquant_nlink
                call qm2_calc_mulliken(qm2_params,qm2_struct,i,scf_mchg(i),density_matrix_unpacked);
            end do
            deallocate(density_matrix_unpacked);
        end if

        CALL qm2_cpt_fock_and_energy(qmmm_mpi, qmmm_scratch, qm_gb, qmewald, &
            & qm2_params, qmmm_nml, qmmm_opnq, &
            & cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct, &
            & SIZE(fock_matrix), fock_matrix, hmatrix, den_matrix, &
            & SIZE(W), W, SIZE(scf_mchg), scf_mchg, density_diff )

        call timer_start(TIME_QMMMENERGYSCFELEC)
        scf_energy = qmmm_struct%elec_eng*EV_TO_KCAL
        energy_diff = scf_energy - eold
        if (qmmm_mpi%commqmmm_master) then
            ! do this only if we need it, this seriously affects performance
            if ( .not. errmat_is_on ) then
                ! check if we shall switch calculation of the error matrix on for diis
                if ( ( ( .not. diis_is_off ) .and. (scf_iteration+qmmm_nml%ndiis_matrices >= switch_on_diis) ) &
                    .or.(scf_iteration+qmmm_nml%ndiis_matrices >= switch_on_diis_last_resort) ) then
                    errmat_is_on = .true.
                end if
            end if
            if ( errmat_is_on ) then
                CALL pack_diis(qmmm_mpi,qmmm_nml,qm2_struct,SIZE(fock_matrix),fock_matrix,den_matrix)
                errval = current_scf_errval(qmmm_nml,qm2_struct)
            end if
        end if
      
        if (qmmm_mpi%commqmmm_master) then
            if (abs(energy_diff) < abs(smallest_energy_diff(1))) then
                !Keep track of the smallest difference we have found. Useful for
                !when we fail to converge - we can tell the user how close we ever got.
                !Ignore changes that are zero.
                if (energy_diff/=zero)  smallest_energy_diff(1) = energy_diff
                if (density_diff/=zero) smallest_energy_diff(2) = density_diff
                if (density_diff/=zero .and. energy_diff/=zero) sm_energy_diff_step_number = scf_iteration
            end if
            !If verbosity is >2 then print some info about this SCF step.
            if (qmmm_nml%verbosity > 2 .and. qmmm_mpi%commqmmm_master) then
                if (doing_pseudo_diag) then
                    write(6,'("QMMM: ",i5," *",A,4G15.7)') scf_iteration, fock_extrap_flag, &
                        & scf_energy, energy_diff, density_diff, errval
                else
                    write(6,'("QMMM: ",i5,"  ",A,4G15.7)') scf_iteration, fock_extrap_flag, &
                        & scf_energy, energy_diff, density_diff, errval
                end if
                if (qmmm_nml%verbosity > 4) then
                    !also print info in KJ/mol
                    write(6,'("QMMM: KJ/mol  ",3G15.7)') scf_energy*4.184d0,energy_diff*4.184d0,density_diff*4.184d0
                end if
            end if !(qmmm_nml%verbosity > 2)
           !End step 4
        
        
        end if
        call timer_stop(TIME_QMMMENERGYSCFELEC)
        if (  qmmm_mpi%commqmmm_master ) then
            ! Step 5 - Check if we have converged.
            !    check energy difference is less than qmmm_nml%scfconv
            !    check density difference is either zero or less than density_conv
            !    Make sure we have done at least 2 SCF iterations to avoid quiting
            !         due to some fluke.
            if (scf_iteration>2) then
                ! Since we have now done more than 2 iterations we can allow the
                ! pseudo diagonalisation as long as it is allowed from the
                ! namelist pseudo_diag option.
                ! Also, if doing pseudo diagonalisations has allowed us to converge
                ! and we are just doing a few more scf steps with full pseudo
                ! diagonalisation then make sure we don't turn pseudo diagonalisation
                ! back on by accident.
           
                if (qmmm_nml%allow_pseudo_diag .AND. (.not. pseudo_converged)) allow_pseudo_diag = .true.
           

           

                ! HMM. DIIS Isn't doing too well once we get close to convergence.
                ! Turn it off in this case... let's just rely on the density matrix extrapolation
                if ( abs(energy_diff) < qmmm_nml%scfconv .AND. .NOT. dont_turn_off_diis ) i = remaining_diis_tokens(qmmm_nml,0)


                ! the SCF is going crazy.  Force diis on as a last resort
                ! IF THE USER REQUESTED US TO!! The default value of
                ! qmmm_nml%ndiis_attempts is 0, so this requires a form
                ! of human intervention
                if ( scf_iteration == switch_on_diis ) then
                    i = remaining_diis_tokens(qmmm_nml, qmmm_nml%ndiis_attempts )
                    dont_turn_off_diis = .TRUE.
              
                    ! Dont allow anymore pseudo diagaonlizations
                    allow_pseudo_diag = .FALSE.
                    doing_pseudo_diag = .FALSE.
                    ! This is a lie, but we don't want it to turn back on
                    pseudo_converged = .TRUE.
                end if


                ! The SCF is STILL going crazy.  Maybe the user didn't
                ! turn the (default OFF) diis parameters on.
                ! I guess we should turn diis ON for the user.
                ! I mean, we're at iteration #800.  What else can we do?
                if ( scf_iteration == switch_on_diis_last_resort ) then
                    i = remaining_diis_tokens(qmmm_nml, 100 )
                    dont_turn_off_diis = .TRUE.
              
                    ! Dont allow anymore pseudo diagaonlizations
                    allow_pseudo_diag = .FALSE.
                    doing_pseudo_diag = .FALSE.
                    ! This is a lie, but we don't want it to turn back on
                    pseudo_converged = .TRUE.
                end if


                ! I empirically observed that pseudo_diag can only get to errval ~1e-5 or so
                ! Further minimization of the pseudo_diag energy does not cause the error matrix to vanish any further.
                ! DIIS tries to nullify the error matrix; therefore, if we reach this point
                ! and still have some diis tokens, then start doing full diagonalizations.
           

                ! Is the error matrix converged? Set to true if we don't check for the error matrix
                errmat_is_converged = errval < qmmm_nml%errconv .OR. ( .not. errmat_is_on )

           
                if (abs(energy_diff) < qmmm_nml%scfconv .AND. &
                    & (density_diff < qmmm_nml%density_conv .OR. density_diff == zero) ) then
              
                    ! We have (pseudo_diag) converged. However, if the last SCF step was a pseudo
                    ! diagonalisation we need to do a few more loops with full
                    ! diagonalisation.  Otherwise the forces will not be accurate enough.
              
                    if (doing_pseudo_diag) then
                        doing_pseudo_diag=.false.
                        allow_pseudo_diag=.false.
                        pseudo_converged=.true. !Stops any more pseudo diags being done for the rest of the SCF loop.
                    else
                        ! pseudo diag may have converged, but the error matrix may not be small enough.
                        converged = errmat_is_converged !.OR. abs(energy_diff) < 0.01d0 * qmmm_nml%scfconv
                    end if
                end if
            end if !scf_iteration>2

#ifdef MPI
        conv_controls(1) = doing_pseudo_diag
        conv_controls(2) = allow_pseudo_diag
        conv_controls(3) = pseudo_converged
        conv_controls(4) = converged
#endif
        end if

#ifdef MPI
     CALL mpi_bcast(conv_controls,4, MPI_LOGICAL,  0, qmmm_mpi%commqmmm, ier)
     doing_pseudo_diag = conv_controls(1)
     allow_pseudo_diag = conv_controls(2)
     pseudo_converged  = conv_controls(3)
     converged         = conv_controls(4)
#endif

        if ( converged ) exit do_scf

         !End of step 5
        eold = scf_energy !Copy the energy that this step got us for use next time.


        ! We have decided that we have not converged
        ! We must diagonalize a fock matrix so that the next iteration can form a new guess
        ! density matrix, provided a set of trial orbitals
        if (qmmm_mpi%commqmmm_master .and. (.not. diis_is_off) &
            .and. ( scf_iteration >= switch_on_diis ) ) then

            if (qmmm_nml%verbosity>4 .and. qmmm_mpi%commqmmm_master) then
                write(6,'("QMMM: DIIS extrapolation is performed")')
            end if

            CALL diis_extrap(qmmm_nml, qm2_struct, qm2_struct%matsize, fock_matrix , fock_extrap_flag )
        end if
     
        CALL qm2_diag(qm2_params,qmmm_mpi,qmmm_scratch, qmmm_nml, &
            & qm2_struct,qm2_struct%matsize, fock_matrix, &
            & allow_pseudo_diag, density_diff, doing_pseudo_diag, &
            & qm2_params%smallsum, qm2_params%abstol )
     
     
        first_iteration = .false.
    end do do_scf ! do scf_iteration, qmmm_nml%itrmax
    ! END MAIN SCF LOOP




    !We get to this point in the code for three reasons.
    ! 1) because we exceeded the maximum number of scf iterations.
    ! 2) becuase we performed the fixed (required) number of iterations
    !     chosen by itrmax<=0
    ! 2) because we converged and so broke out of the above loop.
    ! Check which condition it is based on the value of the logical converged.
    ! Only the master thread returns the scf energy to force.
    if (qmmm_mpi%commqmmm_master) then
        !Here this is just the electronic energy.
        escf = scf_energy
    else
        escf = zero
    end if

    if (.NOT. converged .and. qmmm_nml%itrmax>0.and. qmmm_mpi%commqmmm_master) then
        !Convergence failed. Print a warning message and return with
        !the current unconverged density. This will mean the forces
        !are not accurate but in an MD simulation they may be good
        !enough to allow the user to get out of a potentially bad geometry.
#ifdef SQM
    write(6,'(/,''QMMM: ERROR!'')')
#else
        write(6,'(/,''QMMM: WARNING!'')')
#endif
        write(6,'(''QMMM: Unable to achieve self consistency to the tolerances specified'')')
        write(6,'(''QMMM: No convergence in SCF after '',i6,'' steps.'')') qmmm_nml%itrmax
#ifndef SQM
        write(6,'(''QMMM: Job will continue with unconverged SCF. Warning energies'')')
        write(6,'(''QMMM: and forces for this step will not be accurate.'')')
#endif
        write(6,'(''QMMM: E = '',E12.4,'' DeltaE = '',E12.4,'' DeltaP = '',E12.4)') escf, energy_diff, density_diff
        write(6,'(''QMMM: Smallest DeltaE = '',E12.4,'' DeltaP = '',E12.4,'' Step = '',i6,/)') &
            smallest_energy_diff(1), smallest_energy_diff(2), sm_energy_diff_step_number
#ifdef SQM
    ! CML TEST
	write(6,*)  'Failed geometry'
	do i = 1,qmmm_struct%nquant_nlink
		write(6,*)  qmmm_struct%qm_coords(1,i), qmmm_struct%qm_coords(2,i), qmmm_struct%qm_coords(3,i)
	end do
    ! END CML TEST
    call mexit(6,1)
#endif
    end if

    if (qmmm_nml%verbosity > 0 .and. converged .and. qmmm_mpi%commqmmm_master) then
        write(6,'("QMMM: SCF Converged to ",G10.4," in: ",i5," Cycles ")') &
            qmmm_nml%scfconv,scf_iteration
    end if

    if (qmmm_nml%verbosity > 0 .and. qmmm_nml%itrmax<=0 &
        .and. qmmm_mpi%commqmmm_master) then

        write(6,'("QMMM: SCF got down to ",G10.4," in: ",i5," required cycles ")') &
            abs(energy_diff),scf_iteration
    end if

    ! If we are trying to do Fock matrix prediction based on an extrapolation of previous steps then
    ! we need to store the final fock matrix.

    if (qmmm_nml%fock_predict == 1) call qm2_fock_store(qm2_struct,qm2_struct%matsize, fock_matrix, hmatrix)


    return
end subroutine qm2_scf


subroutine qm2_densit(qm2_params, eigen_vecs,norbs,ndubl, den_matrix,matsize)
    !***********************************************************************
    !
    !   Computes the density matrix from the eigen vector matrix and the
    !   MO occupancy.
    !
    !   eigen_vecs = square eigen vector matrix, of size norbs BY norbs
    !                eigen vectors are stored in top left corner.
    !        norbs = no. orbitals
    !        ndubl = number of doubly occupied M.O.S.
    !
    !    Not currently used:
    !        nsingl= no. double+single occupied M.O.S.
    !
    !   returns den_matrix  = density matrix
    !
    ! Routine and optimization / Open MP by Ross Walker (SDSC, 2007)
    !
    !***********************************************************************
    !

  use qm2_params_module,  only : qm2_params_type
    implicit none

    ! Passed in
    type(qm2_params_type),intent(inout) :: qm2_params
    integer, intent(in) :: norbs, ndubl, matsize
    _REAL_, intent(in) :: eigen_vecs(norbs,norbs)
    _REAL_, intent(out) :: den_matrix(matsize)

    ! Local
    integer unoc_start,i,j,k,l
    _REAL_ eigen_veci


    unoc_start = ndubl + 1


#ifdef OPENMP_DEN
    !$OMP PARALLEL &
    !$OMP Dcosmo_c_struct%EFAULT(PRIVATE) &
    !$OMP SHARED(den_matrix,matsize,unoc_start,norbs,qm2_params,eigen_vecs)
    ! den_matrix can be shared since no two threads will do the same value of L.
    !$OMP DO SCHEDULE(static)
#endif
    do i = 1,matsize
        den_matrix(i)=0.0d0
    end do
#ifdef OPENMP_DEN
    !$OMP END DO
#endif

    do k=unoc_start,norbs
#ifdef OPENMP_DEN
        !$OMP DO SCHEDULE(guided)
        do i=1,norbs
          L=qm2_params%pascal_tri1(i)
#else
        L = 0
        do i=1,norbs
#endif
            eigen_veci = eigen_vecs(i,k)
            do j=1,i
                L=L+1
                den_matrix(L)=den_matrix(L)-2.0d0*eigen_veci*eigen_vecs(j,k)
            end do
        end do
#ifdef OPENMP_DEN
    !$OMP END DO
#endif
    end do

#ifdef OPENMP_DEN
    !$OMP END PARALLEL
#endif

    do i=1,norbs
        L=qm2_params%pascal_tri2(i)
        den_matrix(L)=2.0d0+den_matrix(L)
    end do

    
end subroutine qm2_densit

!-----------------------------------------------------------------------

subroutine qm2_cnvg(qm2_params,den_matrix, old_den_matrix, old2_density,norbs, scf_iteration, density_diff)
    !***********************************************************************
    !
    !  CNVG IS A TWO-POINT INTERPOLATION ROUTINE FOR SPEEDING CONVERGENCE
    !       OF THE DENSITY MATRIX.
    !
    ! ON OUTPUT den_matrix   = NEW DENSITY MATRIX
    !           old2_density  = DIAGONAL OF OLD DENSITY MATRIX
    !           density_diff     = LARGEST DIFFERENCE BETWEEN OLD AND NEW DENSITY
    !                    MATRIX DIAGONAL ELEMENTS
    !***********************************************************************
  use qm2_params_module,  only : qm2_params_type
    implicit none

    !Passed in
    type(qm2_params_type),intent(inout) :: qm2_params
    integer, intent(in) :: norbs, scf_iteration
    _REAL_, intent(inout) :: old2_density(norbs), old_den_matrix(*), den_matrix(*)
    _REAL_, intent(out) :: density_diff

    !Local
    _REAL_ faca, damp, facb, facb_temp, fac, sum0, current_den_sum, sum2, sum3
    _REAL_ den_element_diff, den_matrix_element, old_den_matrix_element, old2_density_element
    integer iminus,k,i,j,ie

    if (scf_iteration > 4) then
        DAMP=0.05D0
    else
        DAMP=1.0d10
    endif

    density_diff=0.0D00
    FACA=0.0D00
    FACB=0.0D00
    FAC=0.0D00
    current_den_sum=0.0D0

    !      RETURN

#ifdef OPENMP_CNVG
    !$OMP PARALLEL &
    !$OMP Dcosmo_c_struct%EFAULT(PRIVATE) &
    !$OMP SHARED(norbs, qm2_params, scf_iteration, den_matrix, old_den_matrix, old2_density, current_den_sum, density_diff, FACA, FACB, FAC, DAMP, SUM2, SUM0, SUM3)
#endif
    if (MOD(scf_iteration-1,3) /= 0) then
#ifdef OPENMP_CNVG
        !$OMP DO SCHEDULE(static) REDUCTION(+:current_den_sum) REDUCTION(max:density_diff)
        do i =1, norbs
           k=qm2_params%pascal_tri2(i)
#else
        K=0
        do I=1,norbs
            K=K+I !pascal_tri2
#endif
            den_matrix_element=den_matrix(K)
            old_den_matrix_element=old_den_matrix(K)
            old2_density(I)=old_den_matrix_element
            density_diff=max(ABS(den_matrix_element-old_den_matrix_element),density_diff)
            old_den_matrix(K)=den_matrix_element
            current_den_sum=current_den_sum+den_matrix_element
        end do
#ifdef OPENMP_CNVG
    !$OMP END DO
#endif
    else
#ifdef OPENMP_CNVG
        !$OMP DO SCHEDULE(static) REDUCTION(+:current_den_sum, faca, facb) REDUCTION(max:density_diff)
        do i =1, norbs
           k=qm2_params%pascal_tri2(i)
#else
        K=0
        do I=1,norbs
            K=K+I !pascal_tri2
#endif
            den_matrix_element=den_matrix(K)
            old_den_matrix_element=old_den_matrix(K)

            current_den_sum=current_den_sum+den_matrix_element
            den_element_diff=ABS(den_matrix_element-old_den_matrix_element)
            density_diff=max(den_element_diff,density_diff)
  
            FACA=FACA+den_element_diff*den_element_diff
            facb_temp = den_matrix_element-2.D00*old_den_matrix_element+old2_density(I)
            FACB=FACB+facb_temp*facb_temp

            old2_density(I)=old_den_matrix_element
            old_den_matrix(K)=den_matrix_element
        end do
#ifdef OPENMP_CNVG
    !$OMP END DO
#endif
    end if
#ifdef OPENMP_CNVG
    !$OMP SINGLE
#endif

    if (FACB > 0.0D00) then
        IF (FACA < (100.D00*FACB)) FAC=SQRT(FACA/FACB)
    end if

    !do i=1 case of loop below
    IF(ABS(old_den_matrix(1)-old2_density(1)) > DAMP) THEN
        old_den_matrix(1)=old2_density(1)+SIGN(DAMP,old_den_matrix(1)-old2_density(1))
    ELSE
        old_den_matrix(1)=old_den_matrix(1)+FAC*(old_den_matrix(1)-old2_density(1))
    ENDIF

    old_den_matrix(1)=MIN(2.0D0,MAX(old_den_matrix(1),0.D0))
    SUM2=old_den_matrix(1)
    den_matrix(1)=old_den_matrix(1)

#ifdef OPENMP_CNVG
    !$OMP END SINGLE
    !$OMP DO SCHEDULE(guided) REDUCTION(+:sum2)
      do i=2,norbs
        ie=qm2_params%pascal_tri1(i)
#else
    IE=1
    do I=2,norbs
#endif
        iminus=I-1
        do J=1,iminus
            IE=IE+1
            den_matrix_element=den_matrix(IE)
            old_den_matrix(IE)=den_matrix_element+FAC*(den_matrix_element-old_den_matrix(IE))
            den_matrix(IE)=old_den_matrix(IE)
        end do
        IE=IE+1
        old2_density_element = old2_density(I)
        old_den_matrix_element = old_den_matrix(IE)
        IF(ABS(old_den_matrix_element-old2_density_element) > DAMP) THEN
            old_den_matrix_element=old2_density_element+SIGN(DAMP,old_den_matrix_element-old2_density_element)
        ELSE
            old_den_matrix_element=old_den_matrix_element+FAC*(old_den_matrix_element-old2_density_element)
        ENDIF
        old_den_matrix_element=MIN(2.0D0,MAX(old_den_matrix_element,0.D0))
        old_den_matrix(IE)=old_den_matrix_element
        den_matrix(IE)=old_den_matrix_element
        SUM2=SUM2+old_den_matrix_element
    end do
#ifdef OPENMP_CNVG
    !$OMP END DO
    !$OMP END PARALLEL
#endif
    !   RE-NORMALIZE IF ANY DENSITY MATRIX ELEMENTS HAVE BEEN TRUNCATED

    SUM0=current_den_sum
    do
        IF(SUM2 > 1.D-3)THEN
            SUM3=current_den_sum/SUM2
        ELSE
            SUM3=0.D0
        ENDIF
        current_den_sum=SUM0
        IF(SUM2 <= 1.D-3 .OR. ABS(SUM3-1.D0) <= 1.D-5) exit
        SUM2=0.D0
        J=0
        do I=1,norbs
            J=J+I
            !  1.0d-20 is in case any occupancy is exactly zero
            old_den_matrix_element=max(old_den_matrix(J)*SUM3+1.D-20, 0.0d0)
                                                                               
            !  SET UP RENORMALIZATION OVER PARTLY OCCUPIED M.O.'S ONLY.  FULL M.O.'S
            !  CAN'T BE FILLED ANY MORE
                                                                               
            IF(old_den_matrix_element > 2.0D0)THEN
                old_den_matrix_element=2.0D0
                current_den_sum=current_den_sum-2.0D0
            ELSE
                SUM2=SUM2+old_den_matrix_element
            ENDIF
            old_den_matrix(J)=old_den_matrix_element
            den_matrix(J)=old_den_matrix_element
        end do
    end do

end subroutine qm2_cnvg

subroutine qm2_cnvg_simple(qm2_params,den_matrix,old_den_matrix,old2_density,norbs, &
    scf_iteration, density_diff)
    !
    !***********************************************************************
    !
    !  Simple two-point interpolation
    !
    ! ON OUTPUT den_matrix   = NEW DENSITY MATRIX
    !           old2_density  = DIAGONAL OF OLD DENSITY MATRIX
    !           density_diff     = LARGEST DIFFERENCE BETWEEN OLD AND NEW DENSITY
    !                    MATRIX DIAGONAL ELEMENTS
    !***********************************************************************
  use qm2_params_module,  only : qm2_params_type
    implicit none

    !Passed in
    type(qm2_params_type),intent(inout) :: qm2_params
    integer, intent(in) :: norbs, scf_iteration
    _REAL_, intent(inout) :: old2_density(norbs), old_den_matrix(*), den_matrix(*)
    _REAL_, intent(out) :: density_diff

    !Local
    _REAL_ faca, damp, facb, facb_temp, fac, sum0, current_den_sum, sum2, sum3
    _REAL_ den_element_diff, den_matrix_element, old_den_matrix_element, old2_density_element
    integer iminus,k,i,j,ie
    integer nbs2

    density_diff=0.0D00

#ifdef OPENMP_CNVG
    !$OMP DO SCHEDULE(static) REDUCTION(+:current_den_sum) REDUCTION(max:density_diff)
        do i =1, norbs
           k=qm2_params%pascal_tri2(i)
#else
    K=0
    do I=1,norbs
        K=K+I !pascal_tri2
#endif
        den_matrix_element=den_matrix(K)
        old_den_matrix_element=old_den_matrix(K)
        density_diff=max(ABS(den_matrix(K)-old_den_matrix(K)),density_diff)

    end do

    ! simple two-point mixing
    nbs2=norbs*(norbs+1)/2
    den_matrix(1:nbs2)=0.5*den_matrix(1:nbs2) &
        +0.5*old_den_matrix(1:nbs2) ! simple interpolation

    old_den_matrix(1:nbs2)=den_matrix(1:nbs2) ! updating old density matrix


end subroutine qm2_cnvg_simple

!---------------------------------------------------------------------------

!---------------------------------------------------------------------------

subroutine qm2_mat_diag(qm2_params,a,n,m,v,e,w1,w2,w3,w4,w5)

    use qm2_params_module,  only : qm2_params_type
    use constants, only : zero, one, half
    implicit none

    !Passed in
    type(qm2_params_type),intent(inout) :: qm2_params
    integer, intent(in) :: m,n
    _REAL_, intent(inout) :: a(*) !Note original matrix will be corrupted by this routine.
    _REAL_, intent(out) ::  e(n), v(n,m)
    _REAL_, intent(out) ::  w1(n), w2(n), w3(n), w4(n), w5(n)

    !**********************************************************************
    !
    ! qm2_mat_diag is a diagonalisation routine, based on:
    ! Yoshitaka Beppu of Nagoya University, Japan.
    !       For details see 'Computers & Chemistry' vol.6 1982. page 000.
    !
    ! This version has been modernised and optimised by
    ! Ross Walker and David Case (TSRI 2005)
    !
    ! on input    a       = matrix to be diagonalised (packed canonical)
    !             n       = size of matrix to be diagonalised.
    !             m       = number of eigenvectors needed.
    !             e       = array of size at least n
    !             v       = array of size at least nmax*m
    !
    ! on output   e       = eigenvalues
    !             v       = eigenvectors in array of size nmax*m
    !
    !***********************************************************************

    !Local
    _REAL_ ssum, r, s, h, summ, ff, ww, hinv, rinv
    _REAL_ c, gersch, z, sinv, ee, ra, t, rn, vn, del, u
    integer nm1, nm2, irank, jrank, krank, i, j, k, l, kp1, ip1
    integer im1, ii, kk, kpiv, ig, itere, ll

    ! eps and eps3 are machine-precision dependent
    _REAL_, parameter :: eps=1.d-16
    _REAL_, parameter :: eps3=1.d-30

    ! Householder transformation

    nm1=n-1
    nm2=n-2
    krank=0
    do k=1,nm2
        kp1=k+1
        krank=krank+k
        w2(k)=a(krank)
        ssum=zero
        jrank=krank
        do j=kp1,n
            w2(j)=a(jrank+k)
            jrank=jrank+j
            ssum=w2(j)*w2(j)+ssum
        end do
        s=sign(sqrt(ssum),w2(kp1))
        w1(k)=-s
        w2(kp1)=w2(kp1)+s
        a(k+krank)=w2(kp1)
        h=w2(kp1)*s
        if(abs(h) < eps3) then
            a(krank) = h
        else
            hinv=one/h
            summ=zero
            irank=krank
            do i=kp1,nm1 !kp1->n-1
                ssum=zero
                do j=kp1,i
                    ssum=ssum+a(j+irank)*w2(j)
                end do
                ip1=i+1
                jrank=qm2_params%pascal_tri2(i)+i
                do j=ip1,n
                    ssum=ssum+a(jrank)*w2(j)
                    jrank=jrank+j
                end do
                w1(i)=ssum*hinv
                irank=irank+i
                summ=w1(i)*w2(i)+summ
            end do
            !Now do the case for i=n
            ssum=zero
            do j=kp1,n
                ssum=ssum+a(j+irank)*w2(j)
            end do
            w1(n)=ssum*hinv
            summ=w1(n)*w2(n)+summ

            u=summ*half*hinv
            jrank=krank
            do j=kp1,n
                w1(j)=w2(j)*u-w1(j)
                do i=kp1,j
                    a(i+jrank)=w1(i)*w2(j)+w1(j)*w2(i)+a(i+jrank)
                end do
                jrank=jrank+j
            end do
            a(krank)=h
        end if !if (abs(h) < eps3)
    end do

    w2(nm1)=a( qm2_params%pascal_tri2(nm1) )
    w2(n)=a( qm2_params%pascal_tri2(n))
    w1(nm1)=a( nm1+qm2_params%pascal_tri1(n) )
    w1(n)=zero
    gersch=abs(w2(1))+abs(w1(1))
    do i=1,nm1
        gersch=max(abs(w2(i+1))+abs(w1(i))+abs(w1(i+1)),gersch)
    end do
    del=eps*gersch
    w3(1:n)=w1(1:n)
    e(1:n)=w2(1:n)
    v(1:n,m)=e(1:n)
    if(abs(del) >= eps3) then

        ! QR-method with origin shift

        k=n
        do
            l=k
            do
                if(abs(w3(l-1)) < del) exit
                l=l-1
                if(l == 1)  exit
            end do
            if(l /= k) then
                ww=(e(k-1)+e(k))*half
                r=e(k)-ww
                z=sign(sqrt(w3(k-1)**2+r*r),r)+ww
                ee=e(l)-z
                e(l)=ee
                ff=w3(l)
                rinv=one/sqrt(ee*ee+ff*ff)
                c=e(l)*rinv
                s=w3(l)*rinv
                ww=e(l+1)-z
                e(l)=(ff*c+ww*s)*s+ee+z
                e(l+1)=c*ww-s*ff
                do j=l+1,k-1
                    rinv=one/sqrt(e(j)*e(j)+w3(j)*w3(j))
                    r = one/rinv
                    w3(j-1)=s*r
                    ee=e(j)*c
                    ff=w3(j)*c
                    c=e(j)*rinv
                    s=w3(j)*rinv
                    ww=e(j+1)-z
                    e(j)=(ff*c+ww*s)*s+ee+z
                    e(j+1)=c*ww-s*ff
                end do
                w3(k-1)=e(k)*s
                e(k)=e(k)*c+z
                cycle
            end if
            k=k-1
            if(k == 1) exit
        end do

        ! At this point the array 'e' contains the un-ordered eigenvalues
        ! Straight selection sort of eigenvalues:

        j=n
        do
            l=1
            ii=1
            ll=1
            do i=2,j
                if( (e(i)-e(l)) >= zero ) then
                    l=i
                else
                    ii=i
                    ll=l
                end if
            end do
            !           if(ii /= ll) then
            ww=e(ll)
            e(ll)=e(ii)
            e(ii)=ww
            !           end if
            j=ii-1
            if(j < 2) exit
        end do
    end if

    ! Ordering of eigenvalues is complete.

    ! Inverse-iteration for eigenvectors:

    rn=zero
    ra=eps*0.6180339887485d0
                 ![0.618... is the fibonacci number (-1+sqrt(5))/2.]
    ig=1
    do i=1,m
        im1=i-1
        do j=1,n
            w3(j)=zero
            w4(j)=w1(j)
            w5(j)=v(j,m)-e(i)
            rn=rn+ra
            if(rn.ge.eps) rn=rn-eps
            v(j,i)=rn
        end do
        do j=1,nm1
            if(abs(w5(j)) < abs(w1(j))) then
                w2(j)=-w5(j)/w1(j)
                w5(j)=w1(j)
                t=w5(j+1)
                w5(j+1)=w4(j)
                w4(j)=t
                w3(j)=w4(j+1)
                if(abs(w3(j)).lt.eps3) w3(j)=del
                w4(j+1)=zero
            else
                if(abs(w5(j)).lt.eps3) w5(j)=del
                w2(j)=-w1(j)/w5(j)
            end if
            w4(j+1)=w3(j)*w2(j)+w4(j+1)
            w5(j+1)=w4(j)*w2(j)+w5(j+1)
        end do
        if(abs(w5(n)) < eps3) w5(n)=del
        !itere=1
        v(n,i)=v(n,i)/w5(n)
        v(nm1,i)=(v(nm1,i)-v(n,i)*w4(nm1))/w5(nm1)
        vn=max(abs(v(n,i)),abs(v(nm1,i)),1.0d-20)
        do k=nm2,1,-1 !n-2,1,-1
            v(k,i)=(v(k,i)-v(k+1,i)*w4(k)-v(k+2,i)*w3(k))/w5(k)
            vn=max(abs(v(k,i)),vn,1.d-20)
        end do
        s=1.0d-5/vn
        v(1:n,i)=v(1:n,i)*s

        do itere=2,5
            do j=1,nm1
                if(abs(w3(j)) >= eps3) then
                    t=v(j,i)
                    v(j,i)=v(j+1,i)
                    v(j+1,i)=t
                    v(j+1,i)=v(j,i)*w2(j)+v(j+1,i)
                end if
            end do
            v(n,i)=v(n,i)/w5(n)
            v(nm1,i)=(v(nm1,i)-v(n,i)*w4(nm1))/w5(nm1)
            vn=max(abs(v(n,i)),abs(v(nm1,i)),1.0d-20)
            do k=nm2,1,-1 !n-2,1,-1
                v(k,i)=(v(k,i)-v(k+1,i)*w4(k)-v(k+2,i)*w3(k))/w5(k)
                vn=max(abs(v(k,i)),vn,1.d-20)
            end do
            s=1.0d-5/vn
            v(1:n,i)=v(1:n,i)*s
            if(vn>1.0d0) exit
        end do

        ! Transformation of eigenvectors

        krank=ishft(nm2*(n+1),-1)
        kpiv=ishft(nm2*nm1,-1)
        do k=nm2,1,-1 !n-2,1,-1
            kp1=k+1
            if(abs(a(kpiv)) > eps3) then
                ssum=zero
                do kk=kp1,n
                    ssum=ssum+a(krank)*v(kk,i)
                    krank=krank+kk
                end do
                s=-ssum/a(kpiv)
                do kk=n,kp1,-1
                    krank=krank-kk
                    v(kk,i)=a(krank)*s+v(kk,i)
                end do
            end if
            kpiv=kpiv-k
            krank=krank-kp1
        end do
        do j=ig,im1 !i-1
            if( abs(e(j)-e(i)) < 0.05d0 ) exit
        end do
        ig=j
        if(ig /= i) then
            ! re-orthogonalisation:

            do k=ig,im1
                ssum=zero
                do j = 1,n
                    ssum=v(j,k)*v(j,i)+ssum
                end do
                s=-ssum
                v(1:n,i)=v(1:n,k)*s+v(1:n,i)
            end do
        end if

        ! Normalisation:

        ssum=1.d-24
        do j=1,n
            ssum=ssum+v(j,i)**2
        end do
        sinv=one/sqrt(ssum)
        v(1:n,i)=v(1:n,i)*sinv
    end do

end subroutine qm2_mat_diag

!---------------------------------------------------------------------------

subroutine qm2_pseudo_diag(qmmm_mpi, qm2_struct,matrix,vectors,noccupied,eigen,norbs,smallsum, &
    matrix_workspace, scratch_matrix, vectmp1, vectmp2, vectmp3, vecjs)

    !------------------------------------------------
    !  "FAST" but "APPROXIMATE" matrix diagonalisation.
    !  See: Stewart, J.J.P., Csaszar, P., Pulay, P., J. Comp. Chem. 3:227,1982
    !
    !  Written by Ross Walker (TSRI, 2005)
    !  Vector Version by Ross Walker (SDSC, 2006)
    !  OpenMP Parallel Version by Ross Walker (SDSC, 2006)
    !
    !------------------------------------------------

    !  qm2_pseudo_diag is a "FAST" but "approximate" matrix diagonalisation routine.
    !  The vectors generated by this routine are more likely to be able to
    !  block-diagonalise the Fock matrix over the molecular orbitals than the
    !  starting vectors. It is based on the work by Stewart et al. and must be
    !  considered to be a pseudo diagonaliser because:
    !
    ! 1) Eigenvectors are not generated. Only the occupied-virtual intersection is
    !    diagonalised, not the secular determinant.
    !
    ! 2) When rotation is used to eliminate elements of the secular determinant the
    !    remainder is assumed to remain unchanged. Any elements that are created
    !    are thus ignored.
    !
    ! 3) The rotation required to eliminate those elements considered significant
    !    is approximated to using the eigenvalues of the exact diagonalisation
    !    throughout the rest of the iterative procedure. In other words the
    !    errors on each pseudo-step are propogated to the next SCF step. The errors
    !    here are assumed to be smaller than the actual error in the SCF at this
    !    point in the SCF.
    !

    !Variables:
    !
    !  MATRIX - On input should contain the matrix to be diagonalised in the
    !           form of a packed lower half triangle.
    !
#ifdef OPENMP
  use qmmm_module, only : qm2_structure, qmmm_mpi_structrue, qmmm_omp
#else
    use qmmm_module, only : qm2_structure, qmmm_mpi_structure
#endif
    implicit none

    ! Passed in
    type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
    type(qm2_structure),intent(inout) :: qm2_struct
    _REAL_, intent(in) :: matrix(qm2_struct%matsize)
    integer, intent(in) :: noccupied, norbs
    _REAL_, intent(inout) :: vectors(norbs,norbs)
    _REAL_, intent(in) :: eigen(norbs)
    _REAL_, intent(in) :: smallsum

    !Passed in Scratch Arrays
    _REAL_, intent(out) :: matrix_workspace(norbs,norbs)
    _REAL_, intent(out) :: scratch_matrix(noccupied,norbs)
    _REAL_, intent(out) :: vectmp1(noccupied*(norbs-noccupied))
    _REAL_, intent(out) :: vectmp2(noccupied*(norbs-noccupied))
    _REAL_, intent(out) :: vectmp3(noccupied*(norbs-noccupied))
    integer, intent(out) :: vecjs(2,noccupied*(norbs-noccupied))

    ! Local
    _REAL_ sum1,a,b,c,d,alpha,beta
    _REAL_ eigeni
    integer i,j,k,j1,k2,kk,m,lumo
    integer :: ii, jj, veccount

#ifdef OPENMP
  call omp_set_num_threads(qmmm_omp%pdiag_threads)
#endif

    lumo=noccupied+1
#ifdef OPENMP
    !$OMP PARALLEL &
    !$OMP Dcosmo_c_struct%EFAULT(PRIVATE) &
    !$OMP SHARED(lumo, norbs, noccupied, matrix, vectors, scratch_matrix, matrix_workspace, veccount, eigeni, eigen, c,d, smallsum, vectmp1,vectmp2,vectmp3,vecjs)
    !workspace can be shared for OMP since no two threads should do the same value of i.
    !$OMP DO SCHEDULE(guided)
#endif
    do i = lumo,norbs
        kk=0
        do j=1,norbs
            sum1 = 0.0d0
            do k=1,j
                kk=kk+1
                sum1=sum1+matrix(kk)*vectors(k,i)
            end do
            j1 = j+1
            k2 = kk
            do k=j1,norbs
                k2=k2+k-1
                sum1=sum1+matrix(k2)*vectors(k,i)
            end do !k=j1,norbs
            matrix_workspace(j,i) = sum1
        end do !j=1,norbs
    end do
#ifdef OPENMP
    !$OMP END DO

    !$OMP DO SCHEDULE(static)
#endif
    do i = lumo,norbs
        do j=1,noccupied
            scratch_matrix(j,i) = sum(matrix_workspace(1:norbs,i)*vectors(1:norbs,j))
        end do !j=1,noccupied
    end do !i=lumo,norbs
#ifdef OPENMP
    !$OMP END DO
    !$OMP SINGLE
#endif
    !Vectored Version by Ross Walker
    veccount=0
    do i=lumo,norbs
        eigeni=eigen(i)
        do j=1,noccupied
            C=scratch_matrix(j,i) !Filling this array in a seperate loop above is faster in most cases
            D=eigen(j)-eigeni     !than doing everything here as a single loop because of dependencies.
            !Check the machine precision for whether to do a 2x2 rotation
            if (abs(C) >= (smallsum*abs(D))) then
                veccount = veccount+1
                vectmp1(veccount) = c
                vectmp2(veccount) = d
                vecjs(1,veccount) = i
                vecjs(2,veccount) = j
            end if
        end do
    end do
#ifdef OPENMP
    !$OMP END SINGLE
    !$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) = 4.0d0*vectmp1(i)*vectmp1(i) + vectmp2(i)*vectmp2(i)
  end do 
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) = 1.0d0/sqrt(vectmp3(i))
  end do
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) =  0.5d0 + abs(0.5d0*vectmp2(i)*vectmp3(i))
  end do
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = 1.0d0-vectmp3(i)
  end do
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = sqrt(vectmp2(i))
  end do
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = -sign(vectmp2(i),vectmp1(i))
  end do
    !$OMP END DO
    !$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp3(i) = sqrt(vectmp3(i))
  end do
    !$OMP END DO
    !$OMP END PARALLEL
#else
    vectmp3(1:veccount) = 4.0d0*vectmp1(1:veccount)*vectmp1(1:veccount)+ &
        vectmp2(1:veccount)*vectmp2(1:veccount)
    call vdinvsqrt(veccount, vectmp3, vectmp3) !E=1.0D0/sqrt(4.0D0*c*c+d*d)
    vectmp3(1:veccount) = 0.5d0 + abs(0.5d0*vectmp2(1:veccount)*vectmp3(1:veccount)) !E=0.5D0+abs(0.5d0*D*E)
    vectmp2(1:veccount) = 1.0d0-vectmp3(1:veccount)
    call vdsqrt(veccount, vectmp2, vectmp2)
    vectmp2(1:veccount) = -sign(vectmp2(1:veccount),vectmp1(1:veccount)) !beta = -sign(sqrt(1.0d0-E),c)
    call vdsqrt(veccount, vectmp3, vectmp3) !alpha = sqrt(E)
#endif

    !-------------------------------------------------------------------
    ! Now we have done the squaring we can do a crude 2 by 2 rotation,
    ! which we will assume eliminates the significant elements.
    !-------------------------------------------------------------------

    do i=1, veccount
        ii = vecjs(1,i)
        jj = vecjs(2,i)
        beta = vectmp2(i)
        alpha = vectmp3(i)
        do m=1, norbs
            a=vectors(m,jj)  !jj is not part of the set of ii.
            b=vectors(m,ii)
            vectors(m,jj)=alpha*a+beta*b
            vectors(m,ii)=alpha*b-beta*a
        end do !m=1,norb
    end do !i=1,veccount

#ifdef OPENMP
  call omp_set_num_threads(1)
#endif

end subroutine qm2_pseudo_diag
 
function qm2_HELECT(nminus,den_matrix,hmatrix,F)
    !***********************************************************************
    !
    !    SUBROUTINE CALCULATES THE ELECTRONIC ENERGY OF THE SYSTEM IN EV.
    !
    !    ON ENTRY Nminus = NUMBER OF ATOMIC ORBITALS - 1
    !             den_matrix = DENSITY MATRIX, PACKED, LOWER TRIANGLE.
    !             hmatrix = ONE-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
    !             F = TWO-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
    !    ON EXIT
    !        qm2_HELECT = ELECTRONIC ENERGY.
    !
    !***********************************************************************
    implicit none

    _REAL_, intent(in) :: den_matrix(*), hmatrix(*), F(*)
    integer, intent(in) :: nminus

    _REAL_ ed, qm2_helect
    _REAL_ fmat;
    integer k,i,j

    ED=0.0D00
    qm2_helect=0.0D00
    fmat = 0.0D00;
    K=0
    do I=1,nminus
        K=K+1
        ED=ED+den_matrix(K)*0.5D0*(hmatrix(K)+F(K))
        fmat=fmat+F(K)
        do J=1,I
            K=K+1
            qm2_helect=qm2_helect+den_matrix(K)*(hmatrix(K)+F(K))
            fmat=fmat+2*F(K);
        end do
    end do

    K=K+1
    ED=ED+den_matrix(K)*0.5D0*(hmatrix(K)+F(K))
    fmat=fmat+F(K)
    qm2_helect=qm2_helect+ED
end function qm2_helect

subroutine qm2_full_diagonalize(qm2_params,qmmm_scratch,qm2_struct,diag_routine,matrix,matrix_dimension,eigen_vectors,abstol)
    !******************************
    !This routine is a central driver routine for doing the diagonalization,
    !essentially it can be called with the matrix to be diagonalized and the
    !type of diagonalization to be done. Note it assumes the matrix and
    !everything else is laid out correctly in memory and the relevant
    !timers have been started or stopped.
    !******************************

#ifdef OPENMP
  use qmmm_module, only : qmmm_scratch_structure, qm2_structure, qmmm_omp
#else
    use qmmm_module, only : qmmm_scratch_structure, qm2_structure
#endif
    use qm2_params_module,  only : qm2_params_type

    implicit none

    !Passed in
    type(qm2_params_type), intent(inout) :: qm2_params 
    type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
    type(qm2_structure),intent(inout) :: qm2_struct
    integer, intent(inout) :: diag_routine !Controls the diagonalization method to be used.
                                        !1=built in diagonalizer
                                        !2=dspev
                                        !3=dspevd
                                        !4=dspevx
                                        !5=dsyev
                                        !6=dsyevd

    _REAL_, intent(inout) :: matrix(*)  !The matrix to be diagonalized
    integer, intent(in) :: matrix_dimension !The dimension of the matrix (1D).
    _REAL_, intent(out) :: eigen_vectors(matrix_dimension,matrix_dimension) !Eigen vectors
    _REAL_, intent(in) :: abstol !Machine underflow limit for dspevr

    !Local
    integer :: ierr, eigen_value_count

#ifdef OPENMP
  call omp_set_num_threads(qmmm_omp%diag_threads)
#endif

    if (diag_routine == 1) then
        ! Note: on our test machine - Intel P4 and Altix this qm2_mat_diag
        ! routine is generally faster than the LAPACK diagonaliser for up
        ! to around 90 QM atoms.

        ! For systems larger than 90 QM atoms the LAPACK diagonaliser
        ! may be quicker but thorough testing on your specific machine,
        ! compiler and LAPACK library combination is necessary
        ! to find the crossover point.

        !Dimension 1 of mat_diag_workspace contains the eignvalues
        call qm2_mat_diag(qm2_params,matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,1), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,2), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,3), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,4), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,5), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,6) )

    else if (diag_routine == 2) then
        ! use dspev
        ! expects mat diag workspace of (norbs,4) - norbs,1 used for eigenvalues, norbs,2->4
        ! is used for scratch space.
        call dspev('V','U',matrix_dimension,matrix, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            eigen_vectors, matrix_dimension, &
            qmmm_scratch%mat_diag_workspace(1,2), &
            ierr)
        ! Check ierr returned zero (fall back to internal diagonalizer below)
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dspev failed to converge and returned ",i8)') ierr
        end if

    else if (diag_routine == 3) then
        !use a divide and conquer lapack diagonaliser
        call dspevd('V','U',matrix_dimension,matrix, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            eigen_vectors, matrix_dimension, &
            qmmm_scratch%lapack_dc_real_scr, &
            qmmm_scratch%lapack_dc_real_scr_aloc, &
            qmmm_scratch%lapack_dc_int_scr,&
            qmmm_scratch%lapack_dc_int_scr_aloc, ierr)
        ! Check ierr returned zero (fall back to internal diagonalizer below)
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dspevd failed to converge and returned ",i8)') ierr
        end if

    else if (diag_routine == 4) then
        ! Use dspevx
        call dspevx('V','A','U',matrix_dimension,matrix, &
            0,0,0,0,abstol, &
            eigen_value_count, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            eigen_vectors, matrix_dimension, &
            qmmm_scratch%lapack_dc_real_scr, &
            qmmm_scratch%lapack_dc_int_scr(matrix_dimension+1), &
            qmmm_scratch%lapack_dc_int_scr, ierr )
        ! Check ierr returned zero (fall back to internal diagonalizer below)
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dspevx failed to converge and returned ",i8)') ierr
        end if

    else if (diag_routine == 5) then
        ! use dsyev - this takes a matrix that is not in packed storage so
        ! since everything is currently done in packed storage we need to unpack the matrix
        ! first. We will unpack the matrix into the eigen vectors array since dsyev overwrites
        ! the matrix with the eigen_vectors anyway.
        call qm2_unpack_matrix(matrix, eigen_vectors, matrix_dimension)
        call dsyev('V','U',matrix_dimension,eigen_vectors,matrix_dimension, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            qmmm_scratch%lapack_dc_real_scr, &
            qmmm_scratch%lapack_dc_real_scr_aloc, &
            ierr)
        !Check ierr returned zero.
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dsyev failed to converge and returned ",i8)') ierr
        end if

    else if (diag_routine == 6) then
        ! use dsyevd - this takes a matrix that is not in packed storage so
        ! since everything is currently done in packed storage we need to unpack the matrix
        ! first. We will unpack the matrix into the eigen vectors array since dsyevd overwrites
        ! the matrix with the eigen_vectors anyway.
        call qm2_unpack_matrix(matrix, eigen_vectors, matrix_dimension)
        call dsyevd('V','U',matrix_dimension,eigen_vectors, matrix_dimension, &
            qmmm_scratch%mat_diag_workspace(1,1),&
            qmmm_scratch%lapack_dc_real_scr, &
            qmmm_scratch%lapack_dc_real_scr_aloc, &
            qmmm_scratch%lapack_dc_int_scr, &
            qmmm_scratch%lapack_dc_int_scr_aloc, &
            ierr)
        !Check ierr returned zero.
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dsyevd failed to converge and returned ",i8)') ierr
        end if

    else if (diag_routine == 7) then
        ! use dsyevr - this takes a matrix that is not in packed storage so we first need to unpack it.
        !              we can't use the eigenvectors array here since dsyevr expects that in a different
        !              memory space. Thus we will use pdiag_scr_norbs_norbs array.
        call qm2_unpack_matrix(matrix, qmmm_scratch%pdiag_scr_norbs_norbs, matrix_dimension)
        call dsyevr('V','A','U',matrix_dimension,qmmm_scratch%pdiag_scr_norbs_norbs,matrix_dimension, &
            0.0d0, 0.0d0, 0, 0, abstol, eigen_value_count, qmmm_scratch%mat_diag_workspace(1,1), &
            eigen_vectors, matrix_dimension, qmmm_scratch%lapack_dc_int_scr, &
            qmmm_scratch%lapack_dc_real_scr, &
            qmmm_scratch%lapack_dc_real_scr_aloc, &
            qmmm_scratch%lapack_dc_int_scr(2*matrix_dimension+1), &
            (qmmm_scratch%lapack_dc_int_scr_aloc-2*matrix_dimension), ierr)
        ! Check ierr returned zero.
        if (ierr /= 0) then
            write(6,'("| QMMM: ERROR dsyevr failed to converge and returned ",i8)') ierr
        end if

    else
        ! Unknown diag method
        call sander_bomb('qm2_full_diagonalize','UNKNOWN DIAGONALIZATION ROUTINE','RUN TERMINATED')

    end if

    if ( (diag_routine > 1) .and. (ierr /= 0) ) then
        ! Diagonalization with one of the lapack routines (diag_routine > 1) failed
        ! Fall back to using the built in diagonalizer (diag_routine = 1)
        write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
        diag_routine = 1
        deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
        REQUIRE(ierr==0)
        allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
        REQUIRE(ierr==0)
        call qm2_mat_diag(qm2_params,matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,1), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,2), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,3), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,4), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,5), &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,6) )
    end if

#ifdef OPENMP
  call omp_set_num_threads(1)
#endif

end subroutine qm2_full_diagonalize

subroutine qm2_unpack_matrix(packed_matrix, unpacked_matrix, matrix_dimension)
    !takes a packed upper triangle in packed_matrix
    !and unpacks it into unpacked_matrix

    implicit none

    integer, intent(in) :: matrix_dimension !the dimension of the matrix
    _REAL_, intent(in) :: packed_matrix(*)
    _REAL_, intent(out) :: unpacked_matrix(matrix_dimension,matrix_dimension)

    !local
    integer :: i, j, jcount

    jcount = 1

    do j = 1, matrix_dimension
        do i = 1, j
            unpacked_matrix(i,j) = packed_matrix(jcount+I-1)
        end do
        jcount = jcount + j
    end do

end subroutine qm2_unpack_matrix


subroutine qm2_pack_matrix(packed_matrix, unpacked_matrix, matrix_dimension)
    !takes an unpacked matrix and packs it as an upper triangle

    implicit none

    integer, intent(in) :: matrix_dimension !the dimension of the matrix
    _REAL_, intent(out) :: packed_matrix(*)
    _REAL_, intent(in) :: unpacked_matrix(matrix_dimension,matrix_dimension)

    !local
    integer :: i, j, jcount

    jcount = 1

    do j = 1, matrix_dimension
        do i = 1, j
            packed_matrix(jcount+I-1) = unpacked_matrix(i,j)
        end do
        jcount = jcount + j
    end do

end subroutine qm2_pack_matrix




SUBROUTINE qm2_diag(qm2_params,qmmm_mpi,qmmm_scratch, qmmm_nml,  qm2_struct, n, fock_matrix,    &
    & allow_pseudo_diag, density_diff, doing_pseudo_diag, & ! control
    & smallsum, abstol                 ) ! machine precision limits

    USE qmmm_module, ONLY : qmmm_mpi_structure
    USE qmmm_module, ONLY : qmmm_scratch_structure  
    use qmmm_nml_module   , only : qmmm_nml_type
    USE qmmm_module, ONLY : qm2_structure
    use qm2_params_module,  only : qm2_params_type
    IMPLICIT NONE
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
    type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
    type(qmmm_nml_type),intent(inout) :: qmmm_nml
    type(qm2_structure),intent(inout) :: qm2_struct
    INTEGER,INTENT(IN)   :: n
    integer :: o,p,q
    _REAL_,INTENT(INOUT) :: fock_matrix(n)
    LOGICAL,INTENT(IN) :: allow_pseudo_diag
    LOGICAL,INTENT(INOUT) :: doing_pseudo_diag
    _REAL_,INTENT(IN) :: density_diff
    _REAL_,INTENT(IN) :: smallsum
    _REAL_,INTENT(IN) :: abstol



#ifdef MPI
  include 'mpif.h'
  integer :: ier
#endif

    ! Step 1 - We haven't converged so we need to get a better Fock matrix.
    !         Diagonalise the RHF secular determinant
    ! Diagonalise the RHF secular determinant
    ! We have two options here. If we are allowed, and the density matrix
    ! fluctuations are small enough we should do a pseudo diagonalisation
    ! instead of a full one.


    if (allow_pseudo_diag .AND. density_diff <= qmmm_nml%pseudo_diag_criteria) then

        !We can do a pseudo diagonalisation.
        doing_pseudo_diag = .true. !Marker used to tell SCF routine to do a full diagonalisation before quitting.

        !Dimension 1 of mat_diag_workspace contains the eignvalues
        !Dimension 2-6 is used as scratch space.

#ifdef MPI
        !only master does the diagonalisation for time being.
     if (qmmm_mpi%commqmmm_master) then
#endif
        call timer_start(TIME_QMMMENERGYSCFPSEUDO)
        !OPENMP PARALLEL
        call qm2_pseudo_diag(qmmm_mpi,qm2_struct,fock_matrix,qm2_struct%eigen_vectors,qm2_struct%nopenclosed, &
            qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,1),qm2_struct%norbs,smallsum, &
            qmmm_scratch%pdiag_scr_norbs_norbs,qmmm_scratch%pdiag_scr_noccupied_norbs, &
            qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2,qmmm_scratch%pdiag_vectmp3, &
            qmmm_scratch%pdiag_vecjs)
        call timer_stop(TIME_QMMMENERGYSCFPSEUDO)
#ifdef MPI
     end if
#endif
    else
#ifdef MPI
        !only master does the diagonalisation for time being.
     if (qmmm_mpi%commqmmm_master) then
#endif
        !Do a full diagonalisation
        call timer_start(TIME_QMMMENERGYSCFDIAG)
        call qm2_full_diagonalize(qm2_params,qmmm_scratch,qm2_struct, qmmm_nml%diag_routine,fock_matrix, &
            qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)

        
        call timer_stop(TIME_QMMMENERGYSCFDIAG)
#ifdef MPI
     end if
#endif
    end if !pseudo diag.
  !End of step 1
        

END SUBROUTINE qm2_diag

!
!********************************************************************
!
!  Construction of Fock matrix
!
!********************************************************************
!

SUBROUTINE qm2_cpt_fock_and_energy(qmmm_mpi, qmmm_scratch, qm_gb, qmewald, qm2_params, qmmm_nml, qmmm_opnq, &
    & cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct, nfock, fock_matrix, hmatrix, den_matrix, &
    & nW, W, nchg, scf_mchg, density_diff)

    USE qmmm_module, ONLY : qmmm_mpi_structure  
    USE qmmm_module, ONLY : qmmm_scratch_structure 
    USE qmmm_module, ONLY : qm2_structure
    USE qmmm_module, ONLY : qm_gb_structure  
    USE qmmm_module, ONLY : qm_ewald_structure  
    use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type
    USE qmmm_module, ONLY : qmmm_opnq_structure
    USE constants, ONLY : zero
    ! module subroutines
    use qm2_fock_d, only : qm2_fock1_d, qm2_fock2_d
    use opnq, only : Opnq_fock
    use qm2_davidson_module
    use cosmo_C , only : cosmo_C_structure !cosmo_c_struct%ceps, cosmo_c_struct%SOLVENT_MODEL, cosmo_c_struct%potential_type, cosmo_c_struct%onsagE, cosmo_c_struct%EF, cosmo_c_struct%rhotzpacked_k
    use qmmm_struct_module, only : qmmm_struct_type

    IMPLICIT NONE
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qmmm_scratch_structure), intent(inout) :: qmmm_scratch
    type(qm_gb_structure), intent(inout) :: qm_gb
    type(qm_ewald_structure), intent(inout) :: qmewald
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qmmm_opnq_structure), intent(inout) :: qmmm_opnq
    type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
    type(qm2_structure),intent(inout) :: qm2_struct
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    INTEGER,INTENT(IN) :: nfock
    _REAL_ ,INTENT(OUT) :: fock_matrix(nfock)
    _REAL_ ,INTENT(IN) :: hmatrix(nfock)
    _REAL_ ,INTENT(IN) :: den_matrix(nfock)
    INTEGER,INTENT(IN) :: nW
    _REAL_ ,INTENT(IN) :: W(nW)
    INTEGER,INTENT(IN) :: nchg
    _REAL_ ,INTENT(IN) :: scf_mchg(nchg)
    _REAL_ ,INTENT(INOUT) :: density_diff
    integer	p,n,m,l;
    _REAL_ ,allocatable :: temp_op(:)
    !qm2_Helect is a function
    _REAL_ qm2_HELECT
    _REAL_ ,allocatable, dimension(:,:) :: density_matrix_unpacked;
    _REAL_ ,allocatable, dimension(:,:) :: fock_matrix_unpacked;
    integer i,j,k
#ifdef MPI
  _REAL_ tmp_recv(2)
# ifndef USE_MPI_IN_PLACE
  _REAL_ tmp_send(2)
# endif
#endif

#ifdef MPI
  include 'mpif.h'
  integer :: ier
#endif


    ! Step 3 - Build the fock matrix !Step 1 if this is our first SCF iteration.
    call timer_start(TIME_QMMMENERGYSCFFOCK)

    ! Copy the one electron matrix into the FOCK matrix
    fock_matrix(1:qm2_struct%matsize)=hmatrix(1:qm2_struct%matsize) !hmatrix added here AND in HELECT energy calculation? JAB

    ! Add the two electron part of the FOCK matrix
    call qm2_fock1_d(qmmm_mpi, qm2_params, qmmm_struct, fock_matrix, den_matrix) !Once center two electron
    call qm2_fock2_d(qm2_params, qmmm_mpi, qm2_struct, qmmm_struct, fock_matrix, den_matrix, W) !Two center two electron

    ! Add Solvent Model or Electric Field
    ! SOLVENT MODEL BLOCK !!JAB
    if ((cosmo_c_struct%solvent_model.gt.0).and.(cosmo_c_struct%solvent_model.ne.10)) then !IF USING SOLVENT MODEL
        allocate(temp_op(nfock))
        if (((cosmo_c_struct%solvent_model.eq.4).or.(cosmo_c_struct%solvent_model.eq.5)) &
		&.and.(.not.qmmm_struct%qm_mm_first_call)) then !Use the excited state density matrix
            if (cosmo_c_struct%potential_type.eq.3) then !USE COSMO
                call addfck(qm2_params, qmmm_nml,cosmo_c_struct,  fock_matrix,den_matrix+cosmo_c_struct%rhotzpacked_k);
                if(cosmo_c_struct%solvent_model.eq.5) then !Variational term
                    qm2ds%tz_scratch=0.d0; temp_op=0.d0; qm2ds%eta=0.d0
                    call addfck(qm2_params, qmmm_nml,cosmo_c_struct, temp_op,den_matrix);
                    call unpacking(qm2ds%nb,temp_op,qm2ds%tz_scratch(1),'s')
                    call calc_xicommutator(cosmo_c_struct,qm2ds,qm2ds%tz_scratch(1))
                    call packing(qm2ds%nb,qm2ds%tz_scratch(1),qm2ds%eta,'s')
                    fock_matrix=fock_matrix+qm2ds%eta(1:nfock)
                    call packing(qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%eta,'u')
                    fock_matrix=fock_matrix+qm2ds%eta(1:nfock)
                endif
            else if (cosmo_c_struct%potential_type.eq.2) then !USE ONSAGER
                qm2ds%tz_scratch=0.d0; temp_op=0.d0; qm2ds%eta=0.d0
                call rcnfld_fock(qm2_params,qmmm_nml,cosmo_c_struct, qm2_struct,qmmm_struct, &
			fock_matrix,den_matrix+cosmo_c_struct%rhotzpacked_k,qm2_struct%norbs);
                call rcnfld_fock(qm2_params,qmmm_nml,cosmo_c_struct, qm2_struct,qmmm_struct,temp_op,den_matrix,qm2_struct%norbs);
                call unpacking(qm2ds%nb,temp_op,qm2ds%tz_scratch(1),'s')
                call calc_xicommutator(cosmo_c_struct,qm2ds,qm2ds%tz_scratch(1))
                call packing(qm2ds%nb,qm2ds%tz_scratch(1),qm2ds%eta,'s')
                fock_matrix=fock_matrix+qm2ds%eta(1:qm2ds%Nb*(qm2ds%Nb+1)/2)
            endif
        else !Use the ground state density matrix
            if (cosmo_c_struct%potential_type.eq.3) then !USE COSMOO
                call addfck(qm2_params, qmmm_nml,cosmo_c_struct, fock_matrix,den_matrix);
            else if (cosmo_c_struct%potential_type.eq.2) then !USE ONSAGER
                call rcnfld_fock(qm2_params,qmmm_nml,cosmo_c_struct, qm2_struct,qmmm_struct, &
                                 fock_matrix,den_matrix,qm2_struct%norbs);
            endif
        endif
        deallocate(temp_op)
    endif
    if (cosmo_c_struct%EF.eq.1) then !USE CONSTANT ELECTRIC FIELD
        allocate(temp_op(nfock))
        temp_op=0.d0
        call efield_fock(qm2_params,qmmm_nml,cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct,temp_op,qm2_struct%norbs);
        fock_matrix=fock_matrix+2.d0*temp_op
        deallocate(temp_op)
    end if
     !!  END SOLVENT MODEL BLOCK

    !AWGTODO: check performance and use the old sp-orbital routines qm2_fock2 and qm2_fock1
    !         if they are faster and d orbitals are not used
    !  !Parallel
    !  call qm2_fock2(fock_matrix,den_matrix,W,qm2_params%orb_loc)
    !  !Parallel
    !  call qm2_fock1(fock_matrix,den_matrix)

#ifndef SQM
    if ( qmmm_nml%qm_ewald>0 ) then
        write(6,*)'Ewald correction'
        call timer_start(TIME_QMMMENERGYSCFFOCKEWALD)
        !Add the qm_ewald contributions to the diagonal elements of
        !the fock matrix. This requires the Mulliken Charges
        !Note we skip this if we are keeping the image charges fixed during the SCF but
        !only after the first MD step has been done.
        if (qmmm_nml%qm_ewald==1 .OR. qmewald%ewald_startup) then
            !Parallel
            call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg, &
                & qmmm_struct%qm_coords,qmewald%kvec)
        end if

        call qm_ewald_add_fock(fock_matrix, qmewald%qmpot, qmewald%mmpot)
        call timer_stop(TIME_QMMMENERGYSCFFOCKEWALD)

    end if
#endif
  
    ! Add OPNQ correction, Taisung Lee (Rutgers, 2011)
    ! the opnq correction of fock matrix should NOT be in hamiltonian when
    ! calculating total electronic energy
    if (qmmm_opnq%useOPNQ) then
        write(6,*)'OPNQ Correction'
        call Opnq_fock(qmmm_opnq, qm2_params,qm2_struct,qmmm_struct, fock_matrix, den_matrix)
    end if

    call timer_stop(TIME_QMMMENERGYSCFFOCK)
    !End step 3

    !Step 4 - Calculate the energy in eV
    ! The fock matrix is dE/dPij, and it is only coincidence
    ! that one can compute the energy as is done here for
    ! HF and semiempirical methods.  ...not very general.
    call timer_start(TIME_QMMMENERGYSCFELEC)
    qmmm_struct%elec_eng = qm2_HELECT(qm2_struct%NORBS-1,den_matrix,hmatrix,fock_matrix);

#ifndef SQM
    if (qmmm_nml%qm_ewald>0) call  qm_ewald_correct_ee(qmmm_struct%elec_eng, &
        & qmewald%mmpot, den_matrix)
#endif


#ifdef MPI
    !Open MP Diagonalisation - only master has density diff - reduce it along with energy
    !Reduce the elec_eng and density diff to all threads so they can check for convergence
    ! Um. (1) This doesn't diagonalize anything.
    !     (2) Why doesn't the master check and then broadcast a signal?
# ifdef USE_MPI_IN_PLACE
  tmp_recv(1) = qmmm_struct%elec_eng
  tmp_recv(2) = density_diff
  call mpi_allreduce(MPI_IN_PLACE,tmp_recv,2, &
       & MPI_DOUBLE_PRECISION,mpi_sum,qmmm_mpi%commqmmm,ier)
  qmmm_struct%elec_eng=tmp_recv(1)
  density_diff=tmp_recv(2)
# else
  tmp_send(1) = qmmm_struct%elec_eng
  tmp_send(2) = density_diff
  call mpi_allreduce(tmp_send,tmp_recv,2, &
       & MPI_DOUBLE_PRECISION,mpi_sum,qmmm_mpi%commqmmm,ier)
  qmmm_struct%elec_eng=tmp_recv(1)
  density_diff=tmp_recv(2)
# endif
#endif

    call timer_stop(TIME_QMMMENERGYSCFELEC)


    !Step 6 - add any extra terms to the SCF for which we only need density convergence.
    ! Fij = dE/dPij
    ! If you modify the fock matrix here, then you should add the corresponding
    ! energy onto the scf energy.  The energy is the variational quantity, afterall.

    !Add in GB terms if qmgb=2
    call timer_start(TIME_QMMMENERGYSCFFOCK)
    if ( qmmm_nml%qmgb == 2 ) then
        call timer_start(TIME_QMMMENERGYSCFFOCKGB)
        !Step 1 - calculate the potential at QM atoms due to QM atoms with current Mulliken charges
        !Parallel
        qm_gb%gb_qmpot(1:qmmm_struct%nquant_nlink)=zero
        call qmgb_calc_qm_pot(qm_gb, qmmm_mpi, qm_gb%gb_qmpot,qm_gb%qmqm_onefij,scf_mchg)
        !Step 2 - Add the mm potential and the qm potential to the fock matrix.
        call qmgb_add_fock(qm2_params, qmmm_struct%nquant_nlink, fock_matrix, qm_gb%gb_mmpot, qm_gb%gb_qmpot)
        call timer_stop(TIME_QMMMENERGYSCFFOCKGB)
    end if
    call timer_stop(TIME_QMMMENERGYSCFFOCK)
  !End step 6


  ! COPY EVERYTHING TO THE MASTER PROCESSOR


#ifdef MPI
  !Reduce the fock matrix to master thread to do diagonalisation
  call timer_start(TIME_QMMMENERGYSCFFOCKRED)
# ifdef USE_MPI_IN_PLACE
  if (qmmm_mpi%commqmmm_master) then
     call mpi_reduce(MPI_IN_PLACE,qm2_struct%fock_matrix,qm2_struct%matsize, &
          MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
  else
     call mpi_reduce(qm2_struct%fock_matrix,0,qm2_struct%matsize, &
          MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
  end if
# else
  call mpi_reduce(qm2_struct%fock_matrix,qmmm_scratch%matsize_red_scratch,qm2_struct%matsize, &
       MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
  if (qmmm_mpi%commqmmm_master) &
       qm2_struct%fock_matrix(1:qm2_struct%matsize)=qmmm_scratch%matsize_red_scratch(1:qm2_struct%matsize)
# endif
  call timer_stop(TIME_QMMMENERGYSCFFOCKRED)
#endif
END SUBROUTINE qm2_cpt_fock_and_energy



SUBROUTINE qm2_densmat(qm2_params, qmmm_mpi,qmmm_nml,qm2_struct, scf_iteration,n,den_matrix,density_diff)

    USE qmmm_module, ONLY : qmmm_mpi_structure, qm2_structure
    USE qm2_iterator_mod, ONLY : scf_iterator_value
    use qmmm_nml_module   , only : qmmm_nml_type
    use qm2_params_module,  only : qm2_params_type
#ifdef MPI
  USE constants, only : zero
#endif

    IMPLICIT NONE
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qm2_params_type), intent(inout) :: qm2_params 
    type(qm2_structure),intent(inout) :: qm2_struct
    INTEGER,INTENT(IN) :: scf_iteration
    INTEGER,INTENT(IN) :: n
    _REAL_,INTENT(OUT) :: den_matrix(n)
    _REAL_,INTENT(OUT) :: density_diff

#ifdef MPI
  include 'mpif.h'
  integer :: ier
#endif


    !Step 2 - Calculate the density matrix:
#ifdef MPI
    !only master calculates density matrix for time being.
  if (qmmm_mpi%commqmmm_master) then
#endif

    call timer_start(TIME_QMMMENERGYSCFDEN)
    !OPENMP PARALLEL
    call qm2_densit(qm2_params, qm2_struct%eigen_vectors,qm2_struct%norbs,qm2_struct%nclosed,den_matrix,qm2_struct%matsize)

    !OPENMP PARALLEL

    call qm2_cnvg(qm2_params,den_matrix, qm2_struct%old_den_matrix, &
        qm2_struct%old2_density,qm2_struct%norbs,scf_iteration,density_diff)

    !==============
    !===============

    call timer_stop(TIME_QMMMENERGYSCFDEN)

#ifdef MPI
  else
     density_diff=zero
  end if
  call timer_start(TIME_QMMMENERGYSCFDENBCAST)
  !only master calculated the density matrix
  call mpi_bcast(qm2_struct%den_matrix, qm2_struct%matsize, MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
  call timer_stop(TIME_QMMMENERGYSCFDENBCAST)
#endif
  !End of step 2

END SUBROUTINE qm2_densmat



SUBROUTINE pack_diis(qmmm_mpi,qmmm_nml,qm2_struct,npf,pf,pp)

    ! assumes nddo
    ! s = unit
    ! x = unit

    USE qmmm_module, ONLY : qm2_structure,qmmm_mpi_structure
    use qmmm_nml_module   , only : qmmm_nml_type

    use qm2_iterator_mod, only : diis_iterator_value

    IMPLICIT NONE
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qm2_structure),intent(inout) :: qm2_struct
    integer,intent(in) :: npf  ! number of packed fock elements
    _REAL_,intent(in) :: pf(npf) ! packed fock
    _REAL_,intent(in) :: pp(npf) ! packed density

    _REAL_,ALLOCATABLE :: f(:,:),p(:,:)

    integer :: nf ! number of fock elements
    integer :: idiis
    integer :: i,j,k


    nf    = qm2_struct%norbs
    idiis = diis_iterator_value(qmmm_nml)

    if (qmmm_nml%verbosity > 4 .and. qmmm_mpi%commqmmm_master) then
        write(6,*) 'pack_diis is performed'
    end if

    ! don't do complicated or obfuscating optimizations until you know it's right!
    ! ...that's why this looks so stupid!

    ALLOCATE( f(nf,nf) )
    f=0.d0
    call qm2_unpack_matrix(pf, f, nf)

    ALLOCATE( p(nf,nf) )
    p=0.d0
    call qm2_unpack_matrix(pp, p, nf)

    DO i=1,nf
        DO j=1,i
            f(i,j) = f(j,i)
        END DO
    END DO
    DO i=1,nf
        DO j=1,i
            p(i,j) = p(j,i)
        END DO
    END DO


    f = MATMUL(f,p)

    k=0
    DO j=1,nf
        DO i=1,nf
            k=k+1
            qm2_struct%diis_errmat(k,idiis) = f(i,j) - f(j,i)
        END DO
    END DO


    qm2_struct%diis_fock(:,idiis) = pf

    DEALLOCATE( f )
    DEALLOCATE( p )

END SUBROUTINE pack_diis



FUNCTION current_scf_errval(qmmm_nml, qm2_struct) RESULT(err)
  
    USE qmmm_module, ONLY : qm2_structure
    USE constants, ONLY : AU_TO_EV
    use qm2_iterator_mod, only : diis_iterator_value
    use qmmm_nml_module   , only : qmmm_nml_type

    IMPLICIT NONE
    type(qmmm_nml_type), intent(inout) :: qmmm_nml 
    type(qm2_structure),intent(inout) :: qm2_struct

    _REAL_ :: err

    integer :: i,idiis


    idiis = diis_iterator_value(qmmm_nml)

    err = 0.d0
    DO i=1,SIZE( qm2_struct%diis_errmat , 1 )
        IF ( err < ABS( qm2_struct%diis_errmat(i,idiis) ) ) err = ABS( qm2_struct%diis_errmat(i,idiis) )
    END DO

    err = err / AU_TO_EV

END FUNCTION current_scf_errval







SUBROUTINE diis_extrap(qmmm_nml, qm2_struct, npf,pf,extrap_flag)

    USE qmmm_module, ONLY : qm2_structure
    use qmmm_nml_module   , only : qmmm_nml_type

    use qm2_iterator_mod, only : scf_iterator_value
    use qm2_iterator_mod, only : diis_iterator_value
    use qm2_iterator_mod, only : diis_iterator_prev_value
    use qm2_iterator_mod, only : remaining_diis_tokens

    IMPLICIT NONE
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qm2_structure),intent(inout) :: qm2_struct
    integer,intent(in) :: npf ! number of packed fock elements
    _REAL_,intent(out) :: pf(npf) ! packed fock matrix
    character(len=1),intent(out) :: extrap_flag

    _REAL_  :: current_scf_errval

    logical :: HAVE_DIIS_DATA
    logical :: ERR_IS_SMALL
    logical :: HAVE_DIIS_TOKENS
    logical :: EXTRAPOLATE
    logical :: EXTRAPOLATION_ERROR

    _REAL_ :: v( 1 + qmmm_nml%ndiis_matrices )
    _REAL_ :: w( 1 + qmmm_nml%ndiis_matrices )
    _REAL_ :: A( 1 + qmmm_nml%ndiis_matrices , 1 + qmmm_nml%ndiis_matrices )

    integer :: scf_iter
    integer :: diis_idx
    integer :: max_diis
    integer :: i,j,n

    extrap_flag = " "

    ! diis decision
    ! should we attempt an extrapolation?

    scf_iter  = scf_iterator_value(qmmm_nml)
    diis_idx  = diis_iterator_value(qmmm_nml)
    max_diis  = qmmm_nml%ndiis_matrices

    HAVE_DIIS_DATA   = scf_iter > 1
    ERR_IS_SMALL     = current_scf_errval(qmmm_nml,qm2_struct) < 0.5d0
    HAVE_DIIS_TOKENS = remaining_diis_tokens(qmmm_nml) > 0
    EXTRAPOLATE      = HAVE_DIIS_DATA .AND. ERR_IS_SMALL .AND. HAVE_DIIS_TOKENS
    IF ( .NOT. EXTRAPOLATE ) THEN
        pf = qm2_struct%diis_fock( : , diis_idx )
        RETURN
    END IF

    ! We will attempt an extrapolation
    ! decrement our tokens

    i = remaining_diis_tokens(qmmm_nml, remaining_diis_tokens(qmmm_nml) - 1 )


    ! construct the diis matrix

    n = MIN( scf_iter , max_diis )

    qm2_struct%diis_mat(1,:) = -1.d0
    qm2_struct%diis_mat(:,1) = -1.d0
    qm2_struct%diis_mat(1,1) =  0.d0
    DO j=1,n
        DO i=1,j
            qm2_struct%diis_mat(i+1,j+1) = DOT_PRODUCT( qm2_struct%diis_errmat(:,i), qm2_struct%diis_errmat(:,j) )
            qm2_struct%diis_mat(j+1,i+1) = qm2_struct%diis_mat(i+1,j+1)
        END DO
    END DO

    v    =  0.d0
    v(1) = -1.d0

    w = 0.d0
    A = 0.d0
  
    CALL SvdInvert_SymMat(n+1,qm2_struct%diis_mat(1:n+1,1:n+1),A(1:n+1,1:n+1),1.D-16)

    w(1:n+1) = MATMUL(A(1:n+1,1:n+1),v(1:n+1))
    ! NOTE: w(1) is a lagrange multiplier


    EXTRAPOLATION_ERROR =  MAXVAL(ABS(w(2:n+1))) >= 1000.d0
    IF ( EXTRAPOLATION_ERROR ) THEN
        !
        ! DIIS YIELDED A NONSENSICAL LINEAR COMBINATION
        ! DO A SIMPLE 50/50 MIXING
        !


        i = diis_iterator_prev_value(qmmm_nml)
        pf = 0.5d0 * ( qm2_struct%diis_fock(:,diis_idx) + qm2_struct%diis_fock(:,i) )

        extrap_flag = "M"

    ELSE


        pf = 0.d0
        DO i=1,n
            pf = pf + w(i+1) * qm2_struct%diis_fock(:,i)
        END DO

        extrap_flag = "D"

    END IF


END SUBROUTINE diis_extrap



SUBROUTINE SvdInvert_SymMat(n,a,ainv,thresh) ! THRESH,ERR
    
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: n
    _REAL_,INTENT(IN) :: a(n,n)
    _REAL_,INTENT(OUT) :: ainv(n,n)
    _REAL_,INTENT(IN) :: thresh

    _REAL_ :: my_thresh
    INTEGER :: my_err


    _REAL_ :: S(n)
    _REAL_ :: U(n,n)
    _REAL_ :: VT(n,n)
    _REAL_ :: twork
    _REAL_,ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK

    INTEGER :: i,j

    my_err = 0
    my_thresh = thresh

    ainv = a
    LWORK = -1
    CALL DGESVD("A","A",n,n, ainv,n, S,U,n, VT,n, twork,LWORK,my_err)
  
    LWORK = NINT(twork)
    ALLOCATE( WORK(LWORK) )
    WORK = 0.d0
  
    CALL DGESVD("A","A",n,n, ainv,n, S,U,n, VT,n, WORK,LWORK,my_err)
  
    IF ( my_err /= 0 ) THEN
        IF ( my_err < 0 ) THEN
            WRITE(6,'(A)')"Svd_SymmMat_Invert: Arg ",-my_err," is illegal"
            WRITE(6,'(A)')"This is a programming error."
            STOP
        ELSE IF ( my_err > 0 ) THEN
            WRITE(6,'(A)')"Svd_SymmMat_Invert: DBDSDC did not converge"
            WRITE(6,'(A)')"Total failure.  Not sure what to do"
            STOP
        END IF
    END IF
  
    DO i=1,n
        IF ( ABS(S(i)) <= my_thresh ) THEN
            S(i) = 0.d0
        ELSE
            S(i) = 1.d0 / S(i)
        END IF
    END DO
  
    DO j=1,n
        DO i=1,n
            VT(i,j) = S(i) * VT(i,j)
        END DO
    END DO
    CALL DGEMM( "N","N",N,N,N, 1.d0,  U,N,  VT,N,  0.d0,  ainv,N )
  
END SUBROUTINE SvdInvert_SymMat
