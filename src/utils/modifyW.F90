
subroutine modifyW(choice, Wreduce, Nmolec1)

!***********************************************************************
! Josiah A. Bjorgaard NDSU 2013
! Scale certain repulsion integrals (W) by Wreduce
! choice=0 ==> all integrals
! chioce=1 ==> integrals for n<=Nmolec1 and m>Nmolec1 (Scale Intermolecular)
! choice=2 ==> integrals for n<=Nmolec1 and m<=Nmolec1 (Scale Intramolecular)
!***********************************************************************

    use qmmm_module, only : qmmm_struct, qm2_struct, qm2_params
    use qm2_davidson_module
    real :: integrals(qm2_struct%n2el)
    real, intent(in) :: Wreduce
    integer, intent(in) :: Nmolec1, choice
    integer :: kr, ki, jstart, jend, i, j, k, &
        j_dimension, i_dimension, ii, jj

    integrals = qm2ds%W

!IF CHOICE 1 or 2 then LOOP OVER ATOMS
    if (choice > 0) then

        if (choice == 1) then
            write (6, *) 'Scaling interamolecular (inner) 2e integrals by', Wreduce
        elseif (choice == 2) then
            write (6, *) 'Scaling intramolecular (outer) 2e integrals by', Wreduce
            qm2_params%onec2elec_params = qm2_params%onec2elec_params*Wreduce
        elseif (choice == 3) then
            write (6, *) 'Scaling all 2e integrals by', Wreduce
            qm2_params%onec2elec_params = qm2_params%onec2elec_params*Wreduce
        end if

        kr = 1

        do I = 2, qmmm_struct%nquant_nlink
            jstart = 1
            jend = i - 1

            do J = jstart, jend

                !SCALE W FOR INTERMOLECULAR OR INTRAMOLECULAR TERMS
                !THIS MACHINERY COMES OUT OF qm2_hcore_qmqm
                i_dimension = qm2_params%natomic_orbs(i)*(qm2_params%natomic_orbs(i) + 1)/2
                j_dimension = qm2_params%natomic_orbs(j)*(qm2_params%natomic_orbs(j) + 1)/2
                ki = i_dimension*j_dimension

                k = 0
                do ii = 1, i_dimension
                    do jj = 1, j_dimension
                        if (( & !Only scale if it's a selected integral
                            (choice == 2) .AND. ((I <= Nmolec1) .OR. (J > Nmolec1)) &
                            ) .OR. ( &
                            (choice == 1) .AND. ((J <= Nmolec1) .AND. (I > Nmolec1))) &
                            .OR. (choice == 3)) then
                            integrals(kr + k) = integrals(kr + k)*Wreduce
                        end if
                        k = k + 1
                    end do
                end do
                kr = kr + ki !increment values
            end do  ! J=1,iminus
        end do !  I=1,qmmm_struct%nquant_nlink

        !OTHERWISE SCALE ALL INTEGRALS
    elseif (choice == 0) then !!JAB
        write (6, *) 'Reducing all 2e integrals'
        integrals = integrals*Wreduce !!JAB
        qm2_params%onec2elec_params = qm2_params%onec2elec_params*Wreduce
    end if !!JAB

    qm2ds%W = integrals

end subroutine modifyW

subroutine averageE()

!***********************************************************************
! Josiah A. Bjorgaard NDSU 2013
! Average and/or scale orbital energies for dimer calculations
! choice=0 ==> average each set of two energies
! chioce=1 ==> scale energies not associated with H(1),H(2),L(1),L(2) where 1
! and 2 refer to the molecule
!***********************************************************************

    use qm2_davidson_module
    integer :: p

    do p = 1, qm2ds%Nb/2
        qm2ds%ehf(2*p) = (qm2ds%ehf(2*p) + qm2ds%ehf(2*p - 1))/2
        qm2ds%ehf(2*p - 1) = qm2ds%ehf(2*p)

    end do

end subroutine averageE

subroutine average_eigen_vectors()

!***********************************************************************
! Josiah A. Bjorgaard NDSU 2013
! Average eigenvector coefficients for dimer calculations
! choice=0 ==> average each set of two energies
! chioce=1 ==> scale energies not associated with H(1),H(2),L(1),L(2) where 1
! and 2 refer to the molecule
!***********************************************************************

    use qmmm_module, only : qmmm_struct, qm2_struct, qm2_params
    use qm2_davidson_module
    !real, intent(in) :: Escale
    integer :: p, p2, h

    do h = 1, qm2_struct%Norbs/2
        do p = 1, qm2_struct%Norbs/2
            p2 = qm2_struct%Norbs/2 + p
            qm2_struct%eigen_vectors(p, 2*h) = 1!(qm2_struct%eigen_vectors(p,2*h)+qm2_struct%eigen_vectors(p2,2*h-1))/2
            qm2_struct%eigen_vectors(p2, 2*h - 1) = 0!qm2_struct%eigen_vectors(p,2*h)
            qm2_struct%eigen_vectors(p2, 2*h) = (qm2_struct%eigen_vectors(p2, 2*h) + qm2_struct%eigen_vectors(p, 2*h - 1))/2
            qm2_struct%eigen_vectors(p, 2*h - 1) = 0!qm2_struct%eigen_vectors(p2,2*h)
        end do
    end do

end subroutine average_eigen_vectors
