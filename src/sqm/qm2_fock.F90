! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_fock1(qmmm_mpi, qm2_params, qmmm_struct, F, PTOT)
! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
!
! Current routine streamlined and optimised by Ross Walker (TSRI, 2005)
!
! *********************************************************************
    use qmmm_module, only : qmmm_mpi_structure
    use qm2_params_module, only : qm2_params_type
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    _REAL_, intent(inout) :: F(*)
    _REAL_, intent(in) :: PTOT(*)

    integer ii, ia, ib, ka, l, m, j, iminus, iplus, icc

    _REAL_ ptpop, GSSII, GSPII, GPPII, GP2II, HSPII
    _REAL_ GSPIIHSPII, GSPIIHSPIIPTK, HSPIIGSPII, GP2IIGPPII
    _REAL_ PTOTKA, PTOTL, PTPOPTL, GPPIIGP2II

    do II = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
        GSSII = qm2_params%onec2elec_params(1, II)
        IA = qm2_params%orb_loc(1, II)
        KA = qm2_params%pascal_tri2(IA)
        PTOTKA = PTOT(KA)
        if (qm2_params%natomic_orbs(ii) == 1) then
            F(KA) = F(KA) + PTOTKA*GSSII
        else
            ! P Orbitals
            GSPII = qm2_params%onec2elec_params(2, II)
            GPPII = qm2_params%onec2elec_params(3, II)
            GP2II = qm2_params%onec2elec_params(4, II)
            HSPII = qm2_params%onec2elec_params(5, II)
            GSPIIHSPII = GSPII - HSPII
            GSPIIHSPIIPTK = GSPIIHSPII*PTOTKA
            HSPIIGSPII = 6.0d0*HSPII - GSPII
            GP2IIGPPII = GP2II - 0.5d0*GPPII
            GPPIIGP2II = 1.5d0*GPPII - GP2II

!
!        The above code assumes that p-orbitals are the last orbitals--which
!        is certanily not true for d-orbitals
!
!        modified by Taisung Lee with the following
!

            IB = qm2_params%orb_loc(1, II) + 1  ! beginning of the p orbital
            PTPOP = PTOT(qm2_params%pascal_tri2(IB)) + &
                PTOT(qm2_params%pascal_tri2(IB + 1)) + PTOT(qm2_params%pascal_tri2(IB + 2))

            F(KA) = F(KA) + PTOTKA*GSSII + PTPOP*GSPIIHSPII
            IPLUS = IA + 1
            L = KA
            do J = IB, IB + 2
                M = L + IA
                L = L + J
                PTOTL = PTOT(L)
                PTPOPTL = PTPOP - PTOTL
                F(L) = F(L) + GSPIIHSPIIPTK + PTOTL*GPPII + PTPOPTL*GP2IIGPPII
                F(M) = F(M) + 0.5d0*PTOT(M)*HSPIIGSPII
            end do

            do J = IB, IB + 2
                ICC = J + 1
                do L = ICC, IB + 2
                    M = qm2_params%pascal_tri1(L) + J
                    F(M) = F(M) + PTOT(M)*GPPIIGP2II
                end do
            end do
        end if
    end do

end subroutine qm2_fock1

subroutine qm2_fock2(qmmm_mpi, qm2_params, qm2_struct, qmmm_struct, F, PTOT, W, orb_loc)
!***********************************************************************
!
! FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
! MATRIX
! ON INPUT  PTOT = TOTAL DENSITY MATRIX.
!           W    = TWO-ELECTRON INTEGRAL MATRIX.
!
!  ON OUTPUT F   = PARTIAL FOCK MATRIX
!***********************************************************************
    use qmmm_module, only : qm2_structure, qmmm_mpi_structure
    use qm2_params_module, only : qm2_params_type
    use qmmm_struct_module, only : qmmm_struct_type
    implicit none

    type(qm2_params_type), intent(inout) :: qm2_params
    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qm2_structure), intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    _REAL_, intent(inout) :: F(*)
    _REAL_, intent(in) :: ptot(*)
    _REAL_, intent(in) :: W(qm2_struct%n2el)
    integer, intent(in) :: orb_loc(2, qmmm_struct%nquant_nlink)

!Local
    integer w_index(qmmm_struct%nquant_nlink, qmmm_struct%nquant_nlink)
    _REAL_ PK(16), PJA(16), PJB(16)
    integer m, i, j, ij, ji, k, l, kl, lk, kk, ii, ia, ib, jk, kj, jj, ja, jb
    integer i1, ll, j1, jstart, jend
    _REAL_ sumdia, sumoff, sum, wkk

    _REAL_::Ftest(10000)

    if (qmmm_struct%fock_first_call) then
        qmmm_struct%fock_first_call = .false.

        !   SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
        !   INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
        !   INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
        !               KINDEX ARE THE INDICES OF THE K-INTEGRALS

        M = 0
        do I = 1, 4
            do J = 1, 4
                IJ = MIN(I, J)
                JI = I + J - IJ
                do K = 1, 4
                    do L = 1, 4
                        M = M + 1
                        KL = MIN(K, L)
                        LK = K + L - KL
                        qm2_struct%fock2_JINDEX(M) = (qm2_params%pascal_tri1(JI) + IJ)*10 &
                            + qm2_params%pascal_tri1(LK) + KL - 10
                    end do
                end do
            end do
        end do
        !  end OF INITIALIZATION
    end if

    !Do the diagonal cases that we don't do below
    !All threads have to do this bit - not in parallel :-(

    do ii = 1, qmmm_struct%nquant_nlink
        IA = orb_loc(1, ii)
        IB = orb_loc(2, ii)
        M = 0
        do J = IA, IB
            do K = IA, IB
                M = M + 1
                JK = MIN(J, K)
                KJ = K + J - JK
                JK = JK + qm2_params%pascal_tri1(KJ)
                qm2_struct%fock2_PTOT2(M, ii) = PTOT(JK)
            end do
        end do
    end do

#ifdef MPI
    KK = qmmm_mpi%two_e_offset
    do ii = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
        jstart = qmmm_mpi%nquant_nlink_jrange(1, ii)
        jend = qmmm_mpi%nquant_nlink_jrange(2, ii)
#else
    KK = 0
    do ii = 2, qmmm_struct%nquant_nlink
        jstart = 1
        jend = ii - 1
#endif
        IA = orb_loc(1, II)
        IB = orb_loc(2, II)
        do JJ = jstart, jend
            JA = orb_loc(1, JJ)
            JB = orb_loc(2, JJ)
            if (IB /= IA .AND. JA /= JB) then ! SP-ATOM  - SP-ATOM
                !   EXTRACT COULOMB TERMS
                PJA(1:16) = qm2_struct%fock2_PTOT2(1:16, ii)
                PJB(1:16) = qm2_struct%fock2_PTOT2(1:16, jj)
                !  COULOMB TERMS
                call qm2_jab(qm2_params, IA, JA, PJA, PJB, W(KK + 1), F)
                !  EXCHANGE TERMS
                !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN
                !      DENSITY MATRIX
                L = 0
                do I = IA, IB
                    I1 = qm2_params%pascal_tri1(I) + JA
                    do J = I1, I1 + 3
                        L = L + 1
                        PK(L) = PTOT(J)*0.5D0
                    end do
                end do
                call qm2_kab(qm2_params, IA, JA, PK, W(KK + 1), F)
                KK = KK + 100
            elseif (IA /= IB) then ! S-ATOM  - SP-ATOM
                !   COULOMB TERMS
                SUMDIA = 0.D0
                SUMOFF = 0.D0
                LL = qm2_params%pascal_tri2(JA)
                K = 0
                do I = 0, 3
                    J1 = qm2_params%pascal_tri1(IA + I) + IA - 1
                    do J = 0, I - 1
                        K = K + 1
                        J1 = J1 + 1
                        WKK = W(KK + K)
                        F(J1) = F(J1) + PTOT(LL)*WKK
                        SUMOFF = SUMOFF + PTOT(J1)*WKK
                    end do
                    J1 = J1 + 1
                    K = K + 1
                    WKK = W(KK + K)
                    F(J1) = F(J1) + PTOT(LL)*WKK
                    SUMDIA = SUMDIA + PTOT(J1)*WKK
                end do
                F(LL) = F(LL) + SUMOFF*2.D0 + SUMDIA
                !  EXCHANGE TERMS
                !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE
                !     SPIN DENSITY MATRIX
                K = 0
                do I = IA, IB
                    I1 = qm2_params%pascal_tri1(I) + JA
                    SUM = 0.D0
                    do J = IA, IB
                        K = K + 1
                        J1 = qm2_params%pascal_tri1(J) + JA
                        SUM = SUM + PTOT(J1)*0.5D0*W(KK + qm2_struct%fock2_JINDEX(K))
                    end do
                    F(I1) = F(I1) - SUM
                end do
                KK = KK + 10
            elseif (JA /= JB) then ! SP-ATOM - S-ATOM
                !   COULOMB TERMS
                SUMDIA = 0.D0
                SUMOFF = 0.D0
                LL = qm2_params%pascal_tri2(IA)
                K = 0
                do I = 0, 3
                    J1 = qm2_params%pascal_tri1(JA + I) + JA - 1
                    do J = 0, I - 1
                        K = K + 1
                        J1 = J1 + 1
                        WKK = W(KK + K)
                        F(J1) = F(J1) + PTOT(LL)*WKK
                        SUMOFF = SUMOFF + PTOT(J1)*WKK
                    end do
                    J1 = J1 + 1
                    K = K + 1
                    WKK = W(KK + K)
                    F(J1) = F(J1) + PTOT(LL)*WKK
                    SUMDIA = SUMDIA + PTOT(J1)*WKK
                end do
                F(LL) = F(LL) + SUMOFF*2.D0 + SUMDIA
                !  EXCHANGE TERMS
                !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE
                !     SPIN DENSITY MATRIX
                K = qm2_params%pascal_tri1(IA) + JA
                J = 0
                do I = K, K + 3
                    SUM = 0.D0
                    do L = K, K + 3
                        J = J + 1
                        SUM = SUM + PTOT(L)*0.5D0*W(KK + qm2_struct%fock2_JINDEX(J))
                    end do
                    F(I) = F(I) - SUM
                end do
                KK = KK + 10
            else !  S-ATOM - S-ATOM
                I1 = qm2_params%pascal_tri2(IA)
                J1 = qm2_params%pascal_tri2(JA)
                IJ = I1 + JA - IA
                WKK = W(KK + 1)
                F(I1) = F(I1) + PTOT(J1)*WKK
                F(J1) = F(J1) + PTOT(I1)*WKK
                F(IJ) = F(IJ) - PTOT(IJ)*0.5D0*WKK
                KK = KK + 1
            end if
        end do !JJ=1,(II-1)
    end do !II=1,nquant_nlink

end subroutine qm2_fock2

subroutine qm2_jab(qm2_params, IA, JA, PJA, PJB, W, F)
    use qm2_params_module, only : qm2_params_type
    implicit none

    type(qm2_params_type), intent(inout) :: qm2_params
    _REAL_, intent(in) :: PJA(16), PJB(16), W(100)
    _REAL_, intent(inout) :: F(*)

    integer i, i5, i6, iia, ija, ia, ja, ioff, joff

    _REAL_ suma(10), sumb(10)
    SUMA(1) = &
        +PJA(1)*W(1) + PJA(2)*W(11) + PJA(3)*W(31) + PJA(4)*W(61) &
        + PJA(5)*W(11) + PJA(6)*W(21) + PJA(7)*W(41) + PJA(8)*W(71) &
        + PJA(9)*W(31) + PJA(10)*W(41) + PJA(11)*W(51) + PJA(12)*W(81) &
        + PJA(13)*W(61) + PJA(14)*W(71) + PJA(15)*W(81) + PJA(16)*W(91)
    SUMA(2) = &
        +PJA(1)*W(2) + PJA(2)*W(12) + PJA(3)*W(32) + PJA(4)*W(62) &
        + PJA(5)*W(12) + PJA(6)*W(22) + PJA(7)*W(42) + PJA(8)*W(72) &
        + PJA(9)*W(32) + PJA(10)*W(42) + PJA(11)*W(52) + PJA(12)*W(82) &
        + PJA(13)*W(62) + PJA(14)*W(72) + PJA(15)*W(82) + PJA(16)*W(92)
    SUMA(3) = &
        +PJA(1)*W(3) + PJA(2)*W(13) + PJA(3)*W(33) + PJA(4)*W(63) &
        + PJA(5)*W(13) + PJA(6)*W(23) + PJA(7)*W(43) + PJA(8)*W(73) &
        + PJA(9)*W(33) + PJA(10)*W(43) + PJA(11)*W(53) + PJA(12)*W(83) &
        + PJA(13)*W(63) + PJA(14)*W(73) + PJA(15)*W(83) + PJA(16)*W(93)
    SUMA(4) = &
        +PJA(1)*W(4) + PJA(2)*W(14) + PJA(3)*W(34) + PJA(4)*W(64) &
        + PJA(5)*W(14) + PJA(6)*W(24) + PJA(7)*W(44) + PJA(8)*W(74) &
        + PJA(9)*W(34) + PJA(10)*W(44) + PJA(11)*W(54) + PJA(12)*W(84) &
        + PJA(13)*W(64) + PJA(14)*W(74) + PJA(15)*W(84) + PJA(16)*W(94)
    SUMA(5) = &
        +PJA(1)*W(5) + PJA(2)*W(15) + PJA(3)*W(35) + PJA(4)*W(65) &
        + PJA(5)*W(15) + PJA(6)*W(25) + PJA(7)*W(45) + PJA(8)*W(75) &
        + PJA(9)*W(35) + PJA(10)*W(45) + PJA(11)*W(55) + PJA(12)*W(85) &
        + PJA(13)*W(65) + PJA(14)*W(75) + PJA(15)*W(85) + PJA(16)*W(95)
    SUMA(6) = &
        +PJA(1)*W(6) + PJA(2)*W(16) + PJA(3)*W(36) + PJA(4)*W(66) &
        + PJA(5)*W(16) + PJA(6)*W(26) + PJA(7)*W(46) + PJA(8)*W(76) &
        + PJA(9)*W(36) + PJA(10)*W(46) + PJA(11)*W(56) + PJA(12)*W(86) &
        + PJA(13)*W(66) + PJA(14)*W(76) + PJA(15)*W(86) + PJA(16)*W(96)
    SUMA(7) = &
        +PJA(1)*W(7) + PJA(2)*W(17) + PJA(3)*W(37) + PJA(4)*W(67) &
        + PJA(5)*W(17) + PJA(6)*W(27) + PJA(7)*W(47) + PJA(8)*W(77) &
        + PJA(9)*W(37) + PJA(10)*W(47) + PJA(11)*W(57) + PJA(12)*W(87) &
        + PJA(13)*W(67) + PJA(14)*W(77) + PJA(15)*W(87) + PJA(16)*W(97)
    SUMA(8) = &
        +PJA(1)*W(8) + PJA(2)*W(18) + PJA(3)*W(38) + PJA(4)*W(68) &
        + PJA(5)*W(18) + PJA(6)*W(28) + PJA(7)*W(48) + PJA(8)*W(78) &
        + PJA(9)*W(38) + PJA(10)*W(48) + PJA(11)*W(58) + PJA(12)*W(88) &
        + PJA(13)*W(68) + PJA(14)*W(78) + PJA(15)*W(88) + PJA(16)*W(98)
    SUMA(9) = &
        +PJA(1)*W(9) + PJA(2)*W(19) + PJA(3)*W(39) + PJA(4)*W(69) &
        + PJA(5)*W(19) + PJA(6)*W(29) + PJA(7)*W(49) + PJA(8)*W(79) &
        + PJA(9)*W(39) + PJA(10)*W(49) + PJA(11)*W(59) + PJA(12)*W(89) &
        + PJA(13)*W(69) + PJA(14)*W(79) + PJA(15)*W(89) + PJA(16)*W(99)
    SUMA(10) = &
        +PJA(1)*W(10) + PJA(2)*W(20) + PJA(3)*W(40) + PJA(4)*W(70) &
        + PJA(5)*W(20) + PJA(6)*W(30) + PJA(7)*W(50) + PJA(8)*W(80) &
        + PJA(9)*W(40) + PJA(10)*W(50) + PJA(11)*W(60) + PJA(12)*W(90) &
        + PJA(13)*W(70) + PJA(14)*W(80) + PJA(15)*W(90) + PJA(16)*W(100)
    SUMB(1) = &
        +PJB(1)*W(1) + PJB(2)*W(2) + PJB(3)*W(4) + PJB(4)*W(7) &
        + PJB(5)*W(2) + PJB(6)*W(3) + PJB(7)*W(5) + PJB(8)*W(8) &
        + PJB(9)*W(4) + PJB(10)*W(5) + PJB(11)*W(6) + PJB(12)*W(9) &
        + PJB(13)*W(7) + PJB(14)*W(8) + PJB(15)*W(9) + PJB(16)*W(10)
    SUMB(2) = &
        +PJB(1)*W(11) + PJB(2)*W(12) + PJB(3)*W(14) + PJB(4)*W(17) &
        + PJB(5)*W(12) + PJB(6)*W(13) + PJB(7)*W(15) + PJB(8)*W(18) &
        + PJB(9)*W(14) + PJB(10)*W(15) + PJB(11)*W(16) + PJB(12)*W(19) &
        + PJB(13)*W(17) + PJB(14)*W(18) + PJB(15)*W(19) + PJB(16)*W(20)
    SUMB(3) = &
        +PJB(1)*W(21) + PJB(2)*W(22) + PJB(3)*W(24) + PJB(4)*W(27) &
        + PJB(5)*W(22) + PJB(6)*W(23) + PJB(7)*W(25) + PJB(8)*W(28) &
        + PJB(9)*W(24) + PJB(10)*W(25) + PJB(11)*W(26) + PJB(12)*W(29) &
        + PJB(13)*W(27) + PJB(14)*W(28) + PJB(15)*W(29) + PJB(16)*W(30)
    SUMB(4) = &
        +PJB(1)*W(31) + PJB(2)*W(32) + PJB(3)*W(34) + PJB(4)*W(37) &
        + PJB(5)*W(32) + PJB(6)*W(33) + PJB(7)*W(35) + PJB(8)*W(38) &
        + PJB(9)*W(34) + PJB(10)*W(35) + PJB(11)*W(36) + PJB(12)*W(39) &
        + PJB(13)*W(37) + PJB(14)*W(38) + PJB(15)*W(39) + PJB(16)*W(40)
    SUMB(5) = &
        +PJB(1)*W(41) + PJB(2)*W(42) + PJB(3)*W(44) + PJB(4)*W(47) &
        + PJB(5)*W(42) + PJB(6)*W(43) + PJB(7)*W(45) + PJB(8)*W(48) &
        + PJB(9)*W(44) + PJB(10)*W(45) + PJB(11)*W(46) + PJB(12)*W(49) &
        + PJB(13)*W(47) + PJB(14)*W(48) + PJB(15)*W(49) + PJB(16)*W(50)
    SUMB(6) = &
        +PJB(1)*W(51) + PJB(2)*W(52) + PJB(3)*W(54) + PJB(4)*W(57) &
        + PJB(5)*W(52) + PJB(6)*W(53) + PJB(7)*W(55) + PJB(8)*W(58) &
        + PJB(9)*W(54) + PJB(10)*W(55) + PJB(11)*W(56) + PJB(12)*W(59) &
        + PJB(13)*W(57) + PJB(14)*W(58) + PJB(15)*W(59) + PJB(16)*W(60)
    SUMB(7) = &
        +PJB(1)*W(61) + PJB(2)*W(62) + PJB(3)*W(64) + PJB(4)*W(67) &
        + PJB(5)*W(62) + PJB(6)*W(63) + PJB(7)*W(65) + PJB(8)*W(68) &
        + PJB(9)*W(64) + PJB(10)*W(65) + PJB(11)*W(66) + PJB(12)*W(69) &
        + PJB(13)*W(67) + PJB(14)*W(68) + PJB(15)*W(69) + PJB(16)*W(70)
    SUMB(8) = &
        +PJB(1)*W(71) + PJB(2)*W(72) + PJB(3)*W(74) + PJB(4)*W(77) &
        + PJB(5)*W(72) + PJB(6)*W(73) + PJB(7)*W(75) + PJB(8)*W(78) &
        + PJB(9)*W(74) + PJB(10)*W(75) + PJB(11)*W(76) + PJB(12)*W(79) &
        + PJB(13)*W(77) + PJB(14)*W(78) + PJB(15)*W(79) + PJB(16)*W(80)
    SUMB(9) = &
        +PJB(1)*W(81) + PJB(2)*W(82) + PJB(3)*W(84) + PJB(4)*W(87) &
        + PJB(5)*W(82) + PJB(6)*W(83) + PJB(7)*W(85) + PJB(8)*W(88) &
        + PJB(9)*W(84) + PJB(10)*W(85) + PJB(11)*W(86) + PJB(12)*W(89) &
        + PJB(13)*W(87) + PJB(14)*W(88) + PJB(15)*W(89) + PJB(16)*W(90)
    SUMB(10) = &
        +PJB(1)*W(91) + PJB(2)*W(92) + PJB(3)*W(94) + PJB(4)*W(97) &
        + PJB(5)*W(92) + PJB(6)*W(93) + PJB(7)*W(95) + PJB(8)*W(98) &
        + PJB(9)*W(94) + PJB(10)*W(95) + PJB(11)*W(96) + PJB(12)*W(99) &
        + PJB(13)*W(97) + PJB(14)*W(98) + PJB(15)*W(99) + PJB(16)*W(100)
    I = 0
    do I5 = 1, 4
        IIA = IA + I5 - 1
        IJA = JA + I5 - 1
        IOFF = qm2_params%pascal_tri1(IIA) + IA - 1
        JOFF = qm2_params%pascal_tri1(IJA) + JA - 1
        do I6 = 1, I5
            IOFF = IOFF + 1
            JOFF = JOFF + 1
            I = I + 1
            F(IOFF) = F(IOFF) + SUMB(I)
            F(JOFF) = F(JOFF) + SUMA(I)
        end do
    end do

end subroutine qm2_jab

subroutine qm2_kab(qm2_params, IA, JA, PK, W, F)
    use qm2_params_module, only : qm2_params_type
    implicit none

    type(qm2_params_type), intent(inout) :: qm2_params

    _REAL_, intent(in) :: PK(*), W(100)
    _REAL_, intent(inout) :: F(*)

    integer m, ia, j, ja, j1, j2, j3

    _REAL_ SUM(16)
    SUM(1) = &
        +PK(1)*W(1) + PK(2)*W(2) + PK(3)*W(4) + PK(4)*W(7) &
        + PK(5)*W(11) + PK(6)*W(12) + PK(7)*W(14) + PK(8)*W(17) &
        + PK(9)*W(31) + PK(10)*W(32) + PK(11)*W(34) + PK(12)*W(37) &
        + PK(13)*W(61) + PK(14)*W(62) + PK(15)*W(64) + PK(16)*W(67)
    SUM(2) = &
        +PK(1)*W(2) + PK(2)*W(3) + PK(3)*W(5) + PK(4)*W(8) &
        + PK(5)*W(12) + PK(6)*W(13) + PK(7)*W(15) + PK(8)*W(18) &
        + PK(9)*W(32) + PK(10)*W(33) + PK(11)*W(35) + PK(12)*W(38) &
        + PK(13)*W(62) + PK(14)*W(63) + PK(15)*W(65) + PK(16)*W(68)
    SUM(3) = &
        +PK(1)*W(4) + PK(2)*W(5) + PK(3)*W(6) + PK(4)*W(9) &
        + PK(5)*W(14) + PK(6)*W(15) + PK(7)*W(16) + PK(8)*W(19) &
        + PK(9)*W(34) + PK(10)*W(35) + PK(11)*W(36) + PK(12)*W(39) &
        + PK(13)*W(64) + PK(14)*W(65) + PK(15)*W(66) + PK(16)*W(69)
    SUM(4) = &
        +PK(1)*W(7) + PK(2)*W(8) + PK(3)*W(9) + PK(4)*W(10) &
        + PK(5)*W(17) + PK(6)*W(18) + PK(7)*W(19) + PK(8)*W(20) &
        + PK(9)*W(37) + PK(10)*W(38) + PK(11)*W(39) + PK(12)*W(40) &
        + PK(13)*W(67) + PK(14)*W(68) + PK(15)*W(69) + PK(16)*W(70)
    SUM(5) = &
        +PK(1)*W(11) + PK(2)*W(12) + PK(3)*W(14) + PK(4)*W(17) &
        + PK(5)*W(21) + PK(6)*W(22) + PK(7)*W(24) + PK(8)*W(27) &
        + PK(9)*W(41) + PK(10)*W(42) + PK(11)*W(44) + PK(12)*W(47) &
        + PK(13)*W(71) + PK(14)*W(72) + PK(15)*W(74) + PK(16)*W(77)
    SUM(6) = &
        +PK(1)*W(12) + PK(2)*W(13) + PK(3)*W(15) + PK(4)*W(18) &
        + PK(5)*W(22) + PK(6)*W(23) + PK(7)*W(25) + PK(8)*W(28) &
        + PK(9)*W(42) + PK(10)*W(43) + PK(11)*W(45) + PK(12)*W(48) &
        + PK(13)*W(72) + PK(14)*W(73) + PK(15)*W(75) + PK(16)*W(78)
    SUM(7) = &
        +PK(1)*W(14) + PK(2)*W(15) + PK(3)*W(16) + PK(4)*W(19) &
        + PK(5)*W(24) + PK(6)*W(25) + PK(7)*W(26) + PK(8)*W(29) &
        + PK(9)*W(44) + PK(10)*W(45) + PK(11)*W(46) + PK(12)*W(49) &
        + PK(13)*W(74) + PK(14)*W(75) + PK(15)*W(76) + PK(16)*W(79)
    SUM(8) = &
        +PK(1)*W(17) + PK(2)*W(18) + PK(3)*W(19) + PK(4)*W(20) &
        + PK(5)*W(27) + PK(6)*W(28) + PK(7)*W(29) + PK(8)*W(30) &
        + PK(9)*W(47) + PK(10)*W(48) + PK(11)*W(49) + PK(12)*W(50) &
        + PK(13)*W(77) + PK(14)*W(78) + PK(15)*W(79) + PK(16)*W(80)
    SUM(9) = &
        +PK(1)*W(31) + PK(2)*W(32) + PK(3)*W(34) + PK(4)*W(37) &
        + PK(5)*W(41) + PK(6)*W(42) + PK(7)*W(44) + PK(8)*W(47) &
        + PK(9)*W(51) + PK(10)*W(52) + PK(11)*W(54) + PK(12)*W(57) &
        + PK(13)*W(81) + PK(14)*W(82) + PK(15)*W(84) + PK(16)*W(87)
    SUM(10) = &
        +PK(1)*W(32) + PK(2)*W(33) + PK(3)*W(35) + PK(4)*W(38) &
        + PK(5)*W(42) + PK(6)*W(43) + PK(7)*W(45) + PK(8)*W(48) &
        + PK(9)*W(52) + PK(10)*W(53) + PK(11)*W(55) + PK(12)*W(58) &
        + PK(13)*W(82) + PK(14)*W(83) + PK(15)*W(85) + PK(16)*W(88)
    SUM(11) = &
        +PK(1)*W(34) + PK(2)*W(35) + PK(3)*W(36) + PK(4)*W(39) &
        + PK(5)*W(44) + PK(6)*W(45) + PK(7)*W(46) + PK(8)*W(49) &
        + PK(9)*W(54) + PK(10)*W(55) + PK(11)*W(56) + PK(12)*W(59) &
        + PK(13)*W(84) + PK(14)*W(85) + PK(15)*W(86) + PK(16)*W(89)
    SUM(12) = &
        +PK(1)*W(37) + PK(2)*W(38) + PK(3)*W(39) + PK(4)*W(40) &
        + PK(5)*W(47) + PK(6)*W(48) + PK(7)*W(49) + PK(8)*W(50) &
        + PK(9)*W(57) + PK(10)*W(58) + PK(11)*W(59) + PK(12)*W(60) &
        + PK(13)*W(87) + PK(14)*W(88) + PK(15)*W(89) + PK(16)*W(90)
    SUM(13) = &
        +PK(1)*W(61) + PK(2)*W(62) + PK(3)*W(64) + PK(4)*W(67) &
        + PK(5)*W(71) + PK(6)*W(72) + PK(7)*W(74) + PK(8)*W(77) &
        + PK(9)*W(81) + PK(10)*W(82) + PK(11)*W(84) + PK(12)*W(87) &
        + PK(13)*W(91) + PK(14)*W(92) + PK(15)*W(94) + PK(16)*W(97)
    SUM(14) = &
        +PK(1)*W(62) + PK(2)*W(63) + PK(3)*W(65) + PK(4)*W(68) &
        + PK(5)*W(72) + PK(6)*W(73) + PK(7)*W(75) + PK(8)*W(78) &
        + PK(9)*W(82) + PK(10)*W(83) + PK(11)*W(85) + PK(12)*W(88) &
        + PK(13)*W(92) + PK(14)*W(93) + PK(15)*W(95) + PK(16)*W(98)
    SUM(15) = &
        +PK(1)*W(64) + PK(2)*W(65) + PK(3)*W(66) + PK(4)*W(69) &
        + PK(5)*W(74) + PK(6)*W(75) + PK(7)*W(76) + PK(8)*W(79) &
        + PK(9)*W(84) + PK(10)*W(85) + PK(11)*W(86) + PK(12)*W(89) &
        + PK(13)*W(94) + PK(14)*W(95) + PK(15)*W(96) + PK(16)*W(99)
    SUM(16) = &
        +PK(1)*W(67) + PK(2)*W(68) + PK(3)*W(69) + PK(4)*W(70) &
        + PK(5)*W(77) + PK(6)*W(78) + PK(7)*W(79) + PK(8)*W(80) &
        + PK(9)*W(87) + PK(10)*W(88) + PK(11)*W(89) + PK(12)*W(90) &
        + PK(13)*W(97) + PK(14)*W(98) + PK(15)*W(99) + PK(16)*W(100)

    if (IA .GT. JA) then
        M = 0
        do J1 = IA, IA + 3
            J = qm2_params%pascal_tri1(j1)
            do J2 = JA, JA + 3
                M = M + 1
                J3 = J + J2
                F(J3) = F(J3) - SUM(M)
            end do
        end do
    else !  IA IS LESS THAN JA, THEREFORE USE OTHER HALF OF TRIANGLE
        M = 0
        do J1 = IA, IA + 3
            do J2 = JA, JA + 3
                M = M + 1
                J3 = qm2_params%pascal_tri1(j2) + j1
                F(J3) = F(J3) - SUM(M)
            end do
        end do
    end if

end subroutine qm2_kab

subroutine qm2_fock2_2atm(qm2_params, qm2_struct, qmmm_struct, F, PTOT, W, orb_loc)
!***********************************************************************
!
! This subroutine is a repetition of qm2_fock2 but for the explicit case
! where there are only 2 atoms. It is mainly used for doing (pseudo)
! numerical QM-QM derivatives where you rebuild a modified 2e-2c fock
! matrix for each pair of atoms.
!
! Written by Ross Walker (SDSC, 2006)
!***********************************************************************

    use ElementOrbitalIndex, only : MaxValenceOrbitals, MaxValenceDimension
    use qm2_params_module, only : qm2_params_type
    use qmmm_module, only : qm2_structure
    use qmmm_struct_module, only : qmmm_struct_type
    implicit none

    type(qm2_params_type), intent(inout) :: qm2_params
    type(qm2_structure), intent(inout) :: qm2_struct
    type(qmmm_struct_type), intent(inout) :: qmmm_struct

! dimension 36 = max of 8*(8+1)/2 = 4 orbs with 4 orbs - S,3P with S,3P
! modified to d-orbital 18*(18+1)/2
! --TL
    integer, parameter :: TwoMaxValence = MaxValenceOrbitals*2*(MaxValenceOrbitals*2 + 1)/2
    _REAL_, intent(inout) :: F(TwoMaxValence)
    _REAL_, intent(in) :: ptot(TwoMaxValence)
    _REAL_, intent(in) :: W(MaxValenceDimension**2)
    integer, intent(in) :: orb_loc(2, 2)

!Local
    _REAL_ PK(MaxValenceOrbitals**2), fock2_ptot2_1(MaxValenceOrbitals**2), fock2_ptot2_2(MaxValenceOrbitals**2)
    integer m, i, j, ij, ji, k, l, kl, lk, ia, ib, jk, kj, ja, jb
    integer i1, ll, j1
    _REAL_ sumdia, sumoff, sum, wkk

    if (qmmm_struct%fock2_2atm_first_call) then
        qmmm_struct%fock2_2atm_first_call = .false.

        !   SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
        !   INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
        !   INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
        !               KINDEX ARE THE INDICES OF THE K-INTEGRALS

        M = 0
        do I = 1, 4
            do J = 1, 4
                IJ = MIN(I, J)
                JI = I + J - IJ
                do K = 1, 4
                    do L = 1, 4
                        M = M + 1
                        KL = MIN(K, L)
                        LK = K + L - KL
                        qm2_struct%fock_2atm_JINDEX(M) = (qm2_params%pascal_tri1(JI) + IJ)*10 &
                            + qm2_params%pascal_tri1(LK) + KL - 10
                    end do
                end do
            end do
        end do
        !  end OF INITIALIZATION
    end if
!RCW: Two loops below should be okay to factor into main loops and avoid excessive
!     array copies. We can do this at some point for speed.
!ii=1
    IA = orb_loc(1, 1)
    IB = orb_loc(2, 1)
    M = 0
    do J = IA, IB
        do K = IA, IB
            M = M + 1
            JK = MIN(J, K)
            KJ = K + J - JK
            JK = JK + qm2_params%pascal_tri1(KJ)
            fock2_PTOT2_1(M) = PTOT(JK)
        end do
    end do
!ii=2
    IA = orb_loc(1, 2)
    IB = orb_loc(2, 2)
    M = 0
    do J = IA, IB
        do K = IA, IB
            M = M + 1
            JK = MIN(J, K)
            KJ = K + J - JK
            JK = JK + qm2_params%pascal_tri1(KJ)
            fock2_PTOT2_2(M) = PTOT(JK)
        end do
    end do
    JA = orb_loc(1, 1)
    JB = orb_loc(2, 1)
    IA = orb_loc(1, 2)
    IB = orb_loc(2, 2)
    if (IB /= IA .AND. JA /= JB) then ! SP-ATOM  - SP-ATOM
        !   EXTRACT COULOMB TERMS
        !  COULOMB TERMS
        call qm2_jab(qm2_params, IA, JA, fock2_ptot2_2, fock2_ptot2_1, W, F)
        !  EXCHANGE TERMS
        !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN
        !      DENSITY MATRIX
        L = 0
        do I = IA, IB
            I1 = qm2_params%pascal_tri1(I) + JA
            do J = I1, I1 + 3
                L = L + 1
                PK(L) = PTOT(J)*0.5D0
            end do
        end do
        call qm2_kab(qm2_params, IA, JA, PK, W, F)
    elseif (IA /= IB) then ! S-ATOM  - SP-ATOM
        !   COULOMB TERMS
        SUMDIA = 0.D0
        SUMOFF = 0.D0
        LL = qm2_params%pascal_tri2(JA)
        K = 0
        do I = 0, 3
            J1 = qm2_params%pascal_tri1(IA + I) + IA - 1
            do J = 0, I - 1
                K = K + 1
                J1 = J1 + 1
                WKK = W(K)
                F(J1) = F(J1) + PTOT(LL)*WKK
                SUMOFF = SUMOFF + PTOT(J1)*WKK
            end do
            J1 = J1 + 1
            K = K + 1
            WKK = W(K)
            F(J1) = F(J1) + PTOT(LL)*WKK
            SUMDIA = SUMDIA + PTOT(J1)*WKK
        end do
        F(LL) = F(LL) + SUMOFF*2.D0 + SUMDIA
        !  EXCHANGE TERMS
        !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE
        !     SPIN DENSITY MATRIX
        K = 0
        do I = IA, IB
            I1 = qm2_params%pascal_tri1(I) + JA
            SUM = 0.D0
            do J = IA, IB
                K = K + 1
                J1 = qm2_params%pascal_tri1(J) + JA
                SUM = SUM + PTOT(J1)*0.5D0*W(qm2_struct%fock_2atm_JINDEX(K))
            end do
            F(I1) = F(I1) - SUM
        end do
    elseif (JA /= JB) then ! SP-ATOM - S-ATOM
        !   COULOMB TERMS
        SUMDIA = 0.D0
        SUMOFF = 0.D0
        LL = qm2_params%pascal_tri2(IA)
        K = 0
        do I = 0, 3
            J1 = qm2_params%pascal_tri1(JA + I) + JA - 1
            do J = 0, I - 1
                K = K + 1
                J1 = J1 + 1
                WKK = W(K)
                F(J1) = F(J1) + PTOT(LL)*WKK
                SUMOFF = SUMOFF + PTOT(J1)*WKK
            end do
            J1 = J1 + 1
            K = K + 1
            WKK = W(K)
            F(J1) = F(J1) + PTOT(LL)*WKK
            SUMDIA = SUMDIA + PTOT(J1)*WKK
        end do
        F(LL) = F(LL) + SUMOFF*2.D0 + SUMDIA
        !  EXCHANGE TERMS
        !  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE
        !     SPIN DENSITY MATRIX
        K = qm2_params%pascal_tri1(IA) + JA
        J = 0
        do I = K, K + 3
            SUM = 0.D0
            do L = K, K + 3
                J = J + 1
                SUM = SUM + PTOT(L)*0.5D0*W(qm2_struct%fock_2atm_JINDEX(J))
            end do
            F(I) = F(I) - SUM
        end do
    else !  S-ATOM - S-ATOM
        I1 = qm2_params%pascal_tri2(IA)
        J1 = qm2_params%pascal_tri2(JA)
        IJ = I1 + JA - IA
        WKK = W(1)
        F(I1) = F(I1) + PTOT(J1)*WKK
        F(J1) = F(J1) + PTOT(I1)*WKK
        F(IJ) = F(IJ) - PTOT(IJ)*0.5D0*WKK
    end if

end subroutine qm2_fock2_2atm

! CML Added to SQM12
subroutine qm2_fock1_skew(qm2_params, qmmm_mpi, qmmm_struct, F, PTOT)
!
! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
!
! Current routine streamlined and optimised by Ross Walker (TSRI, 2005)
!
!  Differs from qm2_fock1:
!
!  qm2_fock1 implicitly assumed symmetric ptot, since it is done
!  assuming the "diagonal" density matrix, not the transition one,
!  which is intrinsically non-symmetric. To deal with that, in Davidson
!  xi - transition density matrix, is split into symmetric and antisymmetric
!  (skew-symmetric) parts.
!  Symmetric part is fed to qm2_fock1. Anti-symmetric part is fed
!  to modified qm2_fock1, i.e., qm2_fock1_skew. qm2_fock2 is
!  general enougn and is not to be modified.
!
!  Notice also that some of the one-center integrals are redifined
!  here (for some reason ???) as compared to original AM1 formulation,
!  for example. Redifining coefficients are given in
!  qm2_load_params_and_allocate.f (around :445).
!
! *********************************************************************
!
    use qm2_params_module, only : qm2_params_type
    use qmmm_module, only : qmmm_mpi_structure
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

    type(qmmm_mpi_structure), intent(inout) :: qmmm_mpi
    type(qm2_params_type), intent(inout) :: qm2_params
    type(qmmm_struct_type), intent(inout) :: qmmm_struct

    _REAL_, intent(inout) :: F(*)
    _REAL_, intent(in) :: PTOT(*)

    integer ii, ia, ib, ka, l, m, j, iminus, iplus, icc
    _REAL_ ptpop, GSSII, GSPII, GPPII, GP2II, HSPII
    _REAL_ GSPIIHSPII, GSPIIHSPIIPTK, HSPIIGSPII, GP2IIGPPII
    _REAL_ PTOTKA, PTOTL, PTPOPTL, GPPIIGP2II

    ! Going through all the quantum atoms
    do II = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end

        GSSII = qm2_params%onec2elec_params(1, II)
        IA = qm2_params%orb_loc(1, II)
        KA = qm2_params%pascal_tri2(IA)
        PTOTKA = PTOT(KA)

        if (qm2_params%natomic_orbs(ii) == 1) then  ! Hydrogen
        else
            ! P Orbitals
            GSPII = qm2_params%onec2elec_params(2, II)
            GPPII = qm2_params%onec2elec_params(3, II)
            GP2II = qm2_params%onec2elec_params(4, II)
            HSPII = qm2_params%onec2elec_params(5, II)
            GSPIIHSPII = GSPII - HSPII
            GSPIIHSPIIPTK = GSPIIHSPII*PTOTKA
            HSPIIGSPII = 6.0d0*HSPII - GSPII
            GP2IIGPPII = GP2II - 0.5d0*GPPII
            GPPIIGP2II = 1.5d0*GPPII - GP2II

            IB = qm2_params%orb_loc(2, II)
            PTPOP = PTOT(qm2_params%pascal_tri2(IB)) &
                + PTOT(qm2_params%pascal_tri2(IB - 1)) + PTOT(qm2_params%pascal_tri2(IB - 2))

            IPLUS = IA + 1
            L = KA
            do J = IPLUS, IB
                M = L + IA
                L = L + J
                PTOTL = PTOT(L)
                PTPOPTL = PTPOP - PTOTL

                F(M) = F(M) + 0.5d0*PTOT(M)*(2.0d0*HSPII - GSPII) ! for antisymmetric
            end do

            IMINUS = IB - 1
            do J = IPLUS, IMINUS
                ICC = J + 1
                do L = ICC, IB
                    M = qm2_params%pascal_tri1(L) + J
                    F(M) = F(M) + PTOT(M)*(0.5*GPPII - 0.6*GP2II) ! for antisymmetric
                end do
            end do
        end if
    end do

    return
end subroutine qm2_fock1_skew
