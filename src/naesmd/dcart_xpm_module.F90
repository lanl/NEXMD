#include "dprec.fh"

module dcart_xpm_module
    implicit none
contains

    ! KGB copied with modifications from CML's version of dcart1
    ! + / - coordinates are passed instead being calculated with a fixed delta
    ! This is intended to be used in nacT.

    ! CML Keep an eye on qm2_get_qm_forces() in future versions to make sure we
    ! CML implement any changes made to them so we can make the same changes in
    ! CML this subroutine, as well as DCART1() and DCART2(), since they are based
    ! CML on the same subroutine. 7/13/12

    function dcart1_xpm(qmmm_nml, qm2_params, qmmm_mpi, qm2_struct, qm2ds, qmmm_struct, ex_dm, xyz_in_p, xyz_in_m) 
        !Current code maintained by: Ross Walker (TSRI 2004)

        !This routine calculates the derivatives of the energy for QM-QM
        !interactions.
        !
        ! xyz_in - QM Atom Coordinates
        ! dxyzqm  - Returned with the forces in for each QM atom.

        use constants          , only : EV_TO_KCAL
        use ElementOrbitalIndex, only: MaxValenceOrbitals
        use qmmm_module        , only : qm2_structure, qmmm_mpi_structure
        use qm2_pm6_hof_module
        use dh_correction_module, only : dh_correction_grad
        use qm2_davidson_module ! CML 7/13/12
        use qmmm_struct_module, only : qmmm_struct_type
	use qm2_params_module,  only : qm2_params_type
	use qmmm_nml_module   , only : qmmm_nml_type

     
        implicit none
        type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
        type(qm2_params_type),intent(inout) :: qm2_params
        type(qmmm_nml_type),intent(inout) :: qmmm_nml
        type(qm2_structure),intent(inout) :: qm2_struct
        type(qm2_davidson_structure_type), intent(inout) :: qm2ds
 
         _REAL_ :: dcart1_xpm
        !Passed in
         type(qmmm_struct_type), intent(in) :: qmmm_struct
        !_REAL_, intent(inout) :: dxyzqm(3,qmmm_struct%nquant_nlink)
        _REAL_, intent(in) :: xyz_in_p(3,qmmm_struct%nquant_nlink), xyz_in_m(3,qmmm_struct%nquant_nlink) !, xyz_in(3,qmmm_struct%nquant_nlink)
        !_REAL_, pointer    :: xyz_in(:,:)
        ! CML Just in case we don't update coords in qmmm_struct 7/13/12
        _REAL_, intent(in) :: ex_dm(qm2ds%Nb*(qm2ds%Nb+1)/2) ! CML 7/13/12
        !Local


        _REAL_ e_repul(22) !Used when qmqm_erep_incore = false
        _REAL_ pair_force(3)
        integer loop_count !Keeps track of number of times through nquant * (nquant-1)/2 loop
        _REAL_ ptzsum(MaxValenceOrbitals**2*3) ! CML for excited state DM 7/13/12
        _REAL_ xyz_qmi(3), xyz_qmj(3), vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
        integer natqmi, natqmj, qmitype, qmjtype
        integer ii, iif, iil, jj, jjf, jjl, ij
        integer i,j,k,l
        integer n_atomic_orbi, n_atomic_orbj
        integer jstart, jend
        _REAL_ aa,ee,deriv,angle,refh,heat,sum
        _REAL_ DENER  !BTN 08/10/17 to store energy derivative
        _REAL_ corei, corej
        _REAL_ betasas, betasap, betapas, betapap
        _REAL_ bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j
        _REAL_ qqi, qqi2, qqj, qqj2, ddi,ddj
        _REAL_ htype, fqmii(3)
        _REAL_, target :: F(MaxValenceOrbitals*(MaxValenceOrbitals*2+1)) !BTN 10/08/2017 place to store fock matrix
        integer natom
          
        !#define change 1.D-4
        !#define halfChange 5.D-5
        !!one/change = 10000
        !#define onechange 10000
        !#define delAdj 1.0D-8
        !#define TWOONEdelAdj 50000000
        dcart1_xpm = 0.d0
        

           !xyz_in => qmmm_struct%qm_coords

        if (qmmm_nml%qmqm_exc_analyt) then !We do analytical derivatives
            ! CML as of right now, fully analytical derivatives for the excited
            ! CML state are not possible. So, we must exit.
            write(6,*)  'QMMM: Analytical excited state derivatives are not implemented! Exiting.'
            call mexit(6,1)
        else !We will do (pseudo numerical derivatives)
            !************** PSEUDO NUMERICAL DERIVATIVES **************
#ifdef MPI
          do ii = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
             jstart =  qmmm_mpi%nquant_nlink_jrange(1,ii)
             jend = qmmm_mpi%nquant_nlink_jrange(2,ii)
#else
            do II=2,qmmm_struct%nquant_nlink
                jstart = 1
                jend = ii-1
#endif
                !Loop over all pairs of quantum atoms
                iif=qm2_params%orb_loc(1,II)
                iil=qm2_params%orb_loc(2,II)
                qmitype = qmmm_struct%qm_atom_type(ii)
                natqmi=qmmm_struct%iqm_atomic_numbers(II)
                do JJ=jstart,jend !jj=1,ii-1
                    !  FORM DIATOMIC MATRICES
                    jjf=qm2_params%orb_loc(1,JJ)
                    jjl=qm2_params%orb_loc(2,JJ)
                    !   GET FIRST ATOM
                    qmjtype = qmmm_struct%qm_atom_type(jj)
                    natqmj=qmmm_struct%iqm_atomic_numbers(JJ)
                    IJ=0
                    do I=jjf,jjl
                        K=qm2_params%pascal_tri1(i)+jjf-1
                        do J=jjf,I
                            IJ=IJ+1
                            K=K+1
                            ptzsum(IJ)=ex_dm(K)
                        end do
                    end do
                    ! GET SECOND ATOM FIRST ATOM INTERSECTION
                    do I=iif,iil
                        L=qm2_params%pascal_tri1(i)
                        K=L+jjf-1
                        do J=jjf,jjl
                            IJ=IJ+1
                            K=K+1
                            ptzsum(IJ)=ex_dm(K)
                        end do
                        K=L+iif-1
                        do L=iif,I
                            K=K+1
                            IJ=IJ+1
                            ptzsum(IJ)=ex_dm(K)
                        end do
                    end do

!                    DENER=qm2_helect1(iil-iif+jjl-jjf+1,ptzsum,qm2_struct%fock_matrix_dm(qm2_params%pascal_tri1(ii-1)+jj,:))   ! CML 7/13/12 BTN 08/10/2017
!                    AA=DENER*2.d0
!                    
!                    DENER=qm2_helect1(iil-iif+jjl-jjf+1,ptzsum,qm2_struct%fock_matrix_dp(qm2_params%pascal_tri1(ii-1)+jj,:))   ! CML 7/13/12 BTN 08/10/2017
!                    EE=DENER*2.d0
!                             
!                    DERIV=(EE-AA)*EV_TO_KCAL
                    DENER=qm2_helect1(iil-iif+jjl-jjf+1,ptzsum,qm2_struct%fock_matrix_dp(qm2_params%pascal_tri1(ii-1)+jj,:))   ! CML 7/13/12 BTN 08/10/2017
                    
                    DERIV=DENER*2.0*EV_TO_KCAL

                    
                    !write(6,*) "Deriv is"
                    !write(6,*) DERIV
                    
                    dcart1_xpm = dcart1_xpm + DERIV
                  !dxyzqm(K,II)=dxyzqm(K,II)-DERIV
                  !dxyzqm(K,JJ)=dxyzqm(K,JJ)+DERIV

                end do
            end do
        !************** end PSEUDO NUMERICAL DERIVATIVES **************
        end if
        return
    end function dcart1_xpm

    subroutine qm2_dhc1(qm2_rij_eqns, qm2_params, qmmm_nml, qm2_struct, qmmm_struct, &
        P,iqm, jqm,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi, & ! CML 7/13/12
        natqmj, iif, iil, jjf, jjl, F)
        !***********************************************************************
        !
        ! Ross Walker (SDSC, 2006) : Do 'Pseudo' Numerical Derivatives for QM
        !
        ! d-orbital extension,  (Taisung Lee, 2011)
        !
        !***********************************************************************

        use constants          , only: ONE, A_TO_BOHRS, A2_TO_BOHRS2, EV_TO_KCAL
        use ElementOrbitalIndex, only: MaxValenceOrbitals,MaxValenceDimension
        use qmmm_module        , only: OVERLAP_CUTOFF, qm2_structure, qm2_rij_eqns_structure
        use Rotation           , only: GetRotationMatrix, Rotate2Center2Electron, RotateCore
        use qm2_fock_d         , only: W2Fock_atompair
        use qmmm_struct_module, only : qmmm_struct_type
        use qm2_params_module,  only : qm2_params_type
        use qmmm_nml_module   , only : qmmm_nml_type
 
        implicit none

        !Passed in
        type(qm2_params_type),intent(inout) :: qm2_params 
        type(qmmm_nml_type),intent(inout) :: qmmm_nml
        type(qm2_structure),intent(inout) :: qm2_struct
        type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
        type(qmmm_struct_type), intent(inout) :: qmmm_struct
        _REAL_ P(:)
        _REAL_, intent(in)  :: xyz_qmi(3),xyz_qmj(3)
        integer, intent(in) :: iqm, jqm, natqmi, natqmj, qmitype, qmjtype
        integer, intent(in) :: iif, iil, jjf, jjl
        _REAL_, intent(out) :: F(:)
        !_REAL_, intent(out), dimension(:), pointer :: F

        ! Local
        integer :: n_atomic_orbi, n_atomic_orbj
        integer :: i_dimension, j_dimension
        integer :: i,j,k,j1,jj,i1, linear, i2, ii, j2
        integer :: firstIndexAO_i, lastIndexAO_i, firstIndexAO_j, lastIndexAO_j

        _REAL_ :: H(MaxValenceOrbitals*(MaxValenceOrbitals*2+1))
        _REAL_ :: SHMAT(MaxValenceOrbitals,MaxValenceOrbitals)
        _REAL_ :: W(MaxValenceDimension**2)
        _REAL_ :: enuclr, ee
        _REAL_ :: r2, rij, r2InAu, rijInAu, oneOverRij
          
        _REAL_ :: RI(22), CORE(10,2)
        _REAL_ :: rotationMatrix(15,45)
        _REAL_, allocatable:: WW(:,:)
        integer ::orb_loc(2,2),KR
        logical::hasDOrbital
        !qm2_Helect1 is a function
        !_REAL_ qm2_helect1    ! CML 7/13/12 BTN:8/7/17 Why would you force reference to external function when function is internal? 
        if (iif < jjf) then
            i=iif-1
            j=jjf-iil+iif-2
        else
            i=iif-jjl+jjf-2
            j=jjf-1
        end if
          
        firstIndexAO_i=iif-i
        lastIndexAO_i=iil-i
        firstIndexAO_j=jjf-j
        lastIndexAO_j=jjl-j
          
        n_atomic_orbi = lastIndexAO_i-firstIndexAO_i+1
        n_atomic_orbj = lastIndexAO_j-firstIndexAO_j+1
          
        linear=qm2_params%pascal_tri2(n_atomic_orbi+n_atomic_orbj)
        F(1:linear)=0.0D0
        H(1:linear)=0.0D0

        ! RCW: Caution, i and j reversed here.

        r2 = sum((xyz_qmj-xyz_qmi)**2)
        rij=sqrt(r2)
        rijInAu=rij*A_TO_BOHRS
        oneOverRij=one/rij
        r2InAu=r2*A2_TO_BOHRS2
            
        if (r2InAu < OVERLAP_CUTOFF) then
          
            if ((n_atomic_orbi.lt.9) .and. (n_atomic_orbj.lt.9)) then  ! SP case
                SHMAT=0.0d0
                call qm2_h1elec(qm2_params,r2InAu,xyz_qmi(1),                                  &
                    xyz_qmj(1),n_atomic_orbi, n_atomic_orbj, SHMAT,           &
                    qmitype,qmjtype)
                
                I2=0
                do I1=firstIndexAO_i,lastIndexAO_i
                    II=qm2_params%pascal_tri1(i1)+firstIndexAO_j-1
                    I2=I2+1
                    J2=0
                    JJ=MIN(I1,lastIndexAO_j)
                    do J1=firstIndexAO_j,JJ
                        II=II+1                                                       
                        J2=J2+1                                                       
                        H(II)=H(II)+SHMAT(I2,J2) 
                    end do
                end do
                
            else  ! for atoms with d orbitals
            
                call qm2_h1elec_d(qm2_params,qm2_struct,r2InAu,xyz_qmi(1:3), xyz_qmj(1:3),  &
                    n_atomic_orbi,n_atomic_orbj,                &
                    firstIndexAO_i, firstIndexAO_j, qmitype, qmjtype,  &
                    n_atomic_orbi+n_atomic_orbj, H)

            end if ! ((n_atomic_orbi.lt.9) .and. (n_atomic_orbj.lt.9))
        end if !(R2 < OVERLAP_CUTOFF)

        KR=1
        hasDOrbital=((n_atomic_orbi.ge.9) .or. (n_atomic_orbj.ge.9))
        call GetRotationMatrix(xyz_qmj-xyz_qmi, rotationMatrix, hasDOrbital)
        call qm2_rotate_qmqm(qmmm_nml, qm2_params, qm2_rij_eqns, qm2_struct,qmmm_struct, &
            -1,iqm,jqm,natqmi,natqmj,xyz_qmi,xyz_qmj,            &
            W(KR),KR, RI, core)

        if (hasDOrbital) then   ! spd case

            i_dimension=n_atomic_orbi*(n_atomic_orbi+1)/2
            j_dimension=n_atomic_orbj*(n_atomic_orbj+1)/2
           
            allocate(ww(j_dimension, i_dimension))
            WW=0.0D0
            ! calculate the 2-center integrals and core-core interaction integrals
            call qm2_repp_d(qm2_params,qmmm_struct, qmitype,qmjtype,rijInAu,RI,CORE,WW,i_dimension,j_dimension,1)

     
            ! put 2-center 2-electron integrals to the linearized matrix W

            k=0
            do ii=1,i_dimension
                do jj=1, j_dimension
                    k=k+1
                    W(k)=WW(jj,ii)
                end do
            end do
                       
            deallocate(ww)
            call Rotate2Center2Electron(W, i_dimension, j_dimension,rotationMatrix)
            
        end if ! ((n_atomic_orbi.ge.9) .and. (n_atomic_orbj.ge.9))
       
        ! calculate the core-core contribution to the H matrix
        ii=qm2_params%pascal_tri2(firstIndexAO_i)
        jj=qm2_params%pascal_tri2(firstIndexAO_j)   
          
        call RotateCore(firstIndexAO_i,firstIndexAO_j,              &
            n_atomic_orbi,n_atomic_orbj,  &
            ii,jj,core,rotationMatrix,H)
       
        call qm2_core_core_repulsion(qmmm_nml, qm2_params, qmmm_struct, iqm, jqm, rij, oneOverRij, RI, enuclr)
            
        ! put what we have now to the Fock matrix
        F(1:linear)=H(1:linear)
           
        ! 2-center 2-electron contribution to the Fock matrix      
        !write(6,*)'W',shape(W)
        !write(6,*)'F',shape(F)
        !write(6,*)'P',shape(P)
        !write(6,*)'n',n_atomic_orbj,n_atomic_orbi,firstIndexAO_J,firstIndexAO_i
        call W2Fock_atompair(W, F, P, n_atomic_orbj, n_atomic_orbi,  &
            firstIndexAO_j, firstIndexAO_i, qmmm_struct%W2Fock_atompair_initialized, qmmm_struct%w_index)
        !write(6,*)'here2.1'
        call W2Fock_atompair(W, F, P, n_atomic_orbi, n_atomic_orbj,  &
            firstIndexAO_i, firstIndexAO_j, qmmm_struct%W2Fock_atompair_initialized, qmmm_struct%w_index)
        !write(6,*)'here2.2' 
        


    end subroutine qm2_dhc1

    function qm2_helect1(nminus,den_matrix,F)
        !***********************************************************************
        !
        !    SUBROUTINE CALCULATES THE ELECTRONIC ENERGY OF THE SYSTEM IN EV.
        !
        !    ON ENTRY Nminus = NUMBER OF ATOMIC ORBITALS - 1
        !             den_matrix = DENSITY MATRIX, PACKED, LOWER TRIANGLE.
        !             hmatrix = ONE-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
        !             F = TWO-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
        !    ON EXIT
        !        qm2_HELECT1 = ELECTRONIC ENERGY.
        !
        !***********************************************************************
        implicit none

        _REAL_, intent(in) :: den_matrix(*), F(*)
        integer, intent(in) :: nminus

        _REAL_ ed, qm2_helect1
        integer k,i,j

        ED=0.0D00
        qm2_helect1=0.0D00
        K=0
        do I=1,nminus
            K=K+1
            ED=ED+den_matrix(K)*0.5D0*F(K)
            do J=1,I
                K=K+1
                qm2_helect1=qm2_helect1+den_matrix(K)*F(K)
            end do
        end do

        K=K+1
        ED=ED+den_matrix(K)*0.5D0*F(K)

        qm2_helect1=qm2_helect1+ED


    end function qm2_helect1


end module


