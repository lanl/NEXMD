! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
#define VACUUM3 603

module genborn

    _REAL_, private, dimension(:), allocatable :: r2x, rjx, vectmp1, vectmp2, &
        vectmp3, vectmp4, vectmp5, sumdeijda, psi
    integer, private, dimension(:), allocatable :: temp_jj, k_vals, j_vals, neckidx
    logical, private, dimension(:), allocatable :: skipv

#ifdef LES
    _REAL_, private, dimension(:), allocatable :: scalefac, lesscalefac, rbornlong, rix, vtemp7, vthi2
    logical, private, dimension(:), allocatable :: longskipv, spreadfrc
    integer, private, dimension(:), allocatable :: nradii, iridx, jridx
#endif

    public allocate_gb, deallocate_gb, egb, igb7_init

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for egb()
    subroutine allocate_gb(natom, ncopy)

        implicit none
        integer, intent(in) :: natom, ncopy
        integer ier

#ifdef LES

        allocate (r2x(ncopy*natom), rjx(ncopy*natom), vectmp1(ncopy*natom), &
            vectmp2(ncopy*natom), &
            vectmp3(ncopy*natom), vectmp4(ncopy*natom), vectmp5(ncopy*natom), &
            sumdeijda(ncopy*natom), psi(ncopy*natom), temp_jj(ncopy*natom), &
            skipv(0:natom), k_vals(natom), j_vals(natom), neckidx(natom), &
            scalefac(natom*ncopy), spreadfrc(natom*ncopy), vthi2(ncopy), &
            longskipv(natom*ncopy), nradii(natom), lesscalefac(natom*ncopy), &
            rbornlong(natom*ncopy), &
            rix(natom*ncopy), iridx(natom*ncopy), jridx(natom*ncopy), &
            vtemp7(ncopy), stat=ier)
#else
        allocate (r2x(natom), rjx(natom), vectmp1(natom), vectmp2(natom), &
            vectmp3(natom), vectmp4(natom), vectmp5(natom), &
            sumdeijda(natom), psi(natom), temp_jj(natom), &
            skipv(0:natom), k_vals(natom), j_vals(natom), neckidx(natom), &
            stat=ier)
#endif

        REQUIRE(ier == 0)
        return

    end subroutine allocate_gb

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for egb()
    subroutine deallocate_gb()

        implicit none
        integer ier

        ! assume that if r2x is allocated then all are allocated
        if (allocated(r2x)) then
            deallocate (skipv, neckidx, j_vals, k_vals, temp_jj, psi, sumdeijda, &
                vectmp5, vectmp4, vectmp3, vectmp2, vectmp1, rjx, r2x, stat=ier)
            REQUIRE(ier == 0)
#ifdef LES
            deallocate (longskipv, scalefac, lesscalefac, nradii, rbornlong, rix, &
                iridx, jridx, vtemp7, vthi2, spreadfrc, stat=ier)
            REQUIRE(ier == 0)
#endif

        else
            REQUIRE(.false.)  ! cannot deallocate un-allocated array
        end if
        return

    end subroutine deallocate_gb

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize table of indexes for GB neck lookup table.
    subroutine igb7_init(natom, ncopy, rborn)

        implicit none
        integer, intent(in) :: natom, ncopy

        integer i, j, k

        _REAL_, intent(in) :: rborn(*)

        do i = 1, natom
            neckidx(i) = nint((rborn(i) - 1.0d0)*20d0)
            if (neckidx(i) < 0 .or. neckidx(i) > 20) then
                write (6, *) "Atom ", i, " has radius ", rborn(i), &
                    " outside of allowed range"
                write (6, *) "of 1.0 to 2.0 angstroms for igb=7. Regenerate &
                &prmtop file with bondi radii."
                call mexit(6, 1)
            end if
        end do

    end subroutine igb7_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ handles generalized Born functionality, plus reg. nonbon, plus surface area
    subroutine egb(x, f, rborn, fs, reff, onereff, charge, iac, ico, numex, &
        natex, dcharge, cut, ntypes, natom, natbel, &
        epol, eelt, evdw, esurf, dvdl, vdwrad, ineighbor, p1, p2, p3, p4, &
        rbmax, rbmin, rbave, rbfluct, ncopy, gbvalpha, gbvbeta, gbvgamma)
        ! Hai Nguyen add gbvalpha,gbvbeta,gbvgamma arrays for GB
        ! put all gbalpha, gbbeta, gbgamma for each atom in these 3 arrays
        ! (look at mdread.f for details)

        !--------------------------------------------------------------------------

        !     Compute nonbonded interactions with a generalized Born model,
        !     getting the "effective" Born radii via the approximate pairwise method
        !     Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
        !     (1996).  Aside from the scaling of the radii, this is the same
        !     approach developed in Schaefer and Froemmel, JMB 216:1045 (1990).

        !     The input coordinates are in the "x" array, and the forces in "f"
        !     get updated; energy components are returned in "epol", "eelt" and
        !     "evdw".

        !     Input parameters for the generalized Born model are "rborn(i)", the
        !     intrinsic dielectric radius of atom "i", and "fs(i)", which is
        !     set (in routine mdread) to (rborn(i) - offset)*si.

        !     Input parameters for the "gas-phase" electrostatic energies are
        !     the charges, in the "charge()" array.

        !     Input parameters for the van der Waals terms are "cn1()" and "cn2()",
        !     containing LJ 12-6 parameters, and "asol" and "bsol" containing
        !     LJ 12-10 parameters.  (The latter are not used in 1994 are more
        !     forcefields.)  The "iac" and "ico" arrays are used to point into
        !     these matrices of coefficients.

        !     The "numex" and "natex" arrays are used to find "excluded" pairs of
        !     atoms, for which gas-phase electrostatics and LJ terms are skipped;
        !     note that GB terms are computed for all pairs of atoms.

        !     If gbsa=1, then an approximate surface-area dependent term is
        !     computed, with the resulting energy placed into "esurf".  The
        !     algorithm is from J. Weiser, P.S. Shenkin, and W.C. Still,
        !     "Approximate atomic sufraces from linear combinations of pariwise
        !     overlaps (LCPO)", J. Computat. Chem. 20:217 (1999).

        !     The code also supports a multiple-time-step facility:

        !       pairs closer than sqrt(cut_inner) are evaluated every nrespai steps;
        !         "   between sqrt(cut_inner) and sqrt(cut) are evaluated
        !                        every nrespa steps
        !         "   beyond sqrt(cut) are ignored

        !       the forces arising from the derivatives of the GB terms with respect
        !          to the effective Born radii are evaulated every nrespa steps

        !       the surface-area dependent term is evaluated every nrespa steps

        !       the effective radii are only updated every nrespai steps

        !     (Be careful with the above: what seems to work is dt=0.001,
        !     nrespai=2, nrespa=4; anything beyond this seems dangerous.)

        !     Written 1999-2000, primarily by D.A. Case, with help from C. Brooks,
        !       T. Simonson, R. Sinkovits  and V. Tsui.  The LCPO implementation
        !       was written by V. Tsui.

        !     Vectorization and optimization 1999-2000, primarily by C. P. Sosa,
        !       T. Hewitt, and D. A. Case.  Work presented at CUG Fall of 2000.
        !--------------------------------------------------------------------------

        use icosasurf, only : icosa_init, icosa_sphere_approx
        use decomp, only : decsasa, decpair
        use qmmm_module, only : qmmm_nml, qmmm_struct, qm2_struct, qm_gb
        use parms, only : cn1, cn2, asol, bsol
        use constants, only : zero, one, two, three, four, five, six, seven, &
            eight, nine, ten, eleven, twelve, half, third, &
            fourth, eighth, pi, fourpi, alpb_alpha, &
            AMBER_ELECTROSTATIC
        use crg_reloc, only : ifcr, cr_add_dcdr_factor

#ifdef LES
        use pimd_vars, only : ipimd, nrg_all
#  ifdef MPI
        use remd, only : rem
#  endif
#endif
        implicit none

#ifdef MPI
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
        include 'mpif.h'
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#endif
#  include "md.h"
#  include "def_time.h"
#if defined(LES)
#  include "les.h"
        _REAL_ lfaci, templfac, tempscale, tmp, addi, addj
        logical first, spread
        integer ijcnum, idx1, idx2, istrt, jstrt, icnum, jcnum, tmpstop, icopy
        integer kstrt, kend
        _REAL_ nrg_vdw_tmp, nrg_egb_tmp, nrg_ele_tmp
#endif

        _REAL_ deel
        logical onstep, onstepi, oncpstep, doeel, dovdw
        _REAL_ x, f, rborn, fs, reff, onereff, charge, cut, &
            epol, eelt, evdw, esurf, vdwrad, p1, p2, p3, p4, &
            dcharge, rbmax, rbmin, rbave, rbfluct
        _REAL_ totsasa, extdieli, intdieli, &
            lastxj, lastyj, lastzj, xi, yi, zi, ri, four_ri, ri1i, xij, yij, zij, &
            dij1i, r2, dij, sj, sj2, frespa, si, sumaij, sumajk, &
            sumaijajk, sumdaijddijdxi, sumdaijddijdyi, sumdaijddijdzi, &
            sumdaijddijdxiajk, sumdaijddijdyiajk, sumdaijddijdziajk, &
            xj, yj, zj, rij, tmpaij, aij, daijddij, daijddijdxj, daijddijdyj, &
            daijddijdzj, sumajk2, sumdajkddjkdxj, sumdajkddjkdyj, &
            sumdajkddjkdzj, p3p4aij, xk, yk, zk, rjk2, djk1i, rjk, vdw2dif, &
            tmpajk, ajk, dajkddjk, dajkddjkdxj, dajkddjkdyj, dajkddjkdzj, &
            daidxj, daidyj, daidzj, ai, daidxi, daidyi, daidzi, qi, dumx, &
            dumy, dumz, de, rj, temp1, fgbi, rinv, r2inv, qiqj, fgbk, expmkf, &
            dl, e, temp4, temp5, temp6, eel, r6inv, f6, f12, &
            dedx, dedy, dedz, qi2h, temp7, dij2i, datmp, dij3i, &
            qid2h, dvdl, thi, thi2, self_e, reff_j
        _REAL_ :: alpb_beta, one_Arad_beta
        _REAL_ :: gbvalpha, gbvbeta, gbvgamma !Hai Nguyen: add gbvalpha,gbvbeta,gbvgamma

#ifdef HAS_10_12
        _REAL_ :: r10inv, f10
#endif

        integer count, count2, icount, ineighbor(*), max_count, iminus
        integer iac, ico, numex, natex, ntypes, natom, natbel, ncopy
        integer i, j, k, kk1, maxi, num_j_vals, jjj, count2_fin, num_k_vals, &
            iexcl, iaci, jexcl, jexcl_last, jjv, ic, kk
        integer j3
        integer loop_count
#ifdef MPI
        integer mpistart, ierr
#endif
        _REAL_ f_x, f_y, f_z, f_xi, f_yi, f_zi
        _REAL_ dumbo, tmpsd, rborn_i, psi_i
        _REAL_ gba_i, gbb_i, gbg_i ! Hai Nguyen: temporary variables for GB (calculating GB force)

        ! Variables for QMMM specific loops
        integer qm_temp_count
!Locals for link atoms
        _REAL_ :: forcemod(3)
        integer :: lnk_no, mm_no, qm_no

        ! variables needed for icosa surface area calculation
        integer ineighborpt

        ! variables needed for smooth integration cutoff in Reff:
        _REAL_ rgbmax2, rgbmax1i, rgbmax2i, rgbmaxpsmax2
        !     _REAL_  datmp2

        _REAL_ onekappa !1/kappa

#include "gbneck.h"

#ifdef LES
        dimension x(*), f(*), rborn(*), charge(*), iac(*), &
            ico(*), numex(*), natex(*), fs(*), reff(ncopy*natom), onereff(ncopy*natom), &
            vdwrad(*), p1(*), p2(*), p3(*), p4(*), &
            dcharge(*), rbmax(*), rbmin(*), rbave(*), rbfluct(*)
#else
        dimension x(*), f(*), rborn(*), charge(*), iac(*), &
            ico(*), numex(*), natex(*), fs(*), reff(natom), onereff(natom), &
            vdwrad(*), p1(*), p2(*), p3(*), p4(*), &
            dcharge(*), rbmax(*), rbmin(*), rbave(*), rbfluct(*), gbvalpha(*), gbvbeta(*), gbvgamma(*) ! Hai Nguyen
#endif

        !   FGB taylor coefficients follow
        !   from A to H :
        !   1/3 , 2/5 , 3/7 , 4/9 , 5/11
        !   4/3 , 12/5 , 24/7 , 40/9 , 60/11

        _REAL_ ta
        _REAL_ tb
        _REAL_ tc
        _REAL_ td
        _REAL_ tdd
        _REAL_ te
        _REAL_ tf
        _REAL_ tg
        _REAL_ th
        _REAL_ thh
        parameter(ta=third)
        parameter(tb=two/five)
        parameter(tc=three/seven)
        parameter(td=four/nine)
        parameter(tdd=five/eleven)
        parameter(te=four/three)
        parameter(tf=twelve/five)
        parameter(tg=three*eight/seven)
        parameter(th=five*eight/nine)
        parameter(thh=three*four*five/eleven)

!============ QMMM ==========================================================
!  If ifqnt == True and igb /=0 and igb /= 6 then we will do EGB with QMMM.

!  In this case the charge array should currently be zero for the QM atoms.
!  This corresponds to qmgb = 0 where we do only EGB(MM-MM) -- skipping
!  QM-MM and QM-QM interactions.

!  If qmgb==1 then we need to copy back the qm resp charges from
!  qm_resp_charges, so we can do EGB on everything. When done we need to
!  ensure we zero the QM charge array again.

!  Before doing the non-bonded electrostatics:
        if (qmmm_nml%ifqnt) then
            if (qmmm_nml%qmgb == 1) then

                call qmmm_restore_mm_charges(qmmm_nml, qmmm_struct%nquant, qmmm_struct%qm_resp_charges, charge, &
                    qmmm_struct%scaled_mm_charges, qmmm_struct%iqmatoms, &
                    qmmm_nml%chg_lambda, qmmm_struct%nlink, qmmm_struct%link_pairs, &
                    qmmm_struct%mm_link_pair_resp_charges, &
                    .false.)

                ! Unlike ewald where we don't explicitly skip QM-MM interactions,
                ! in GB we do explicitly skip them so it doesn't matter that the
                ! QM charges are not set to zero in the charge array.

            else if (qmmm_nml%qmgb > 1) then

                ! In this case we need to fill the charge array with the Mulliken
                ! charges from the SCF.  Note if igb==2 we skip the radii calculation
                ! below as this has already been done before the call to qm_mm.
                do qm_temp_count = 1, qmmm_struct%nquant_nlink
                    charge(qmmm_struct%iqmatoms(qm_temp_count)) = &
                        qm2_struct%scf_mchg(qm_temp_count)*AMBER_ELECTROSTATIC
                end do
            end if

            !We need to replace the MM link pair coordinates with
            !the MM link atom coordinates.
            call adj_mm_link_pair_crd(qmmm_nml, qmmm_struct, x)

            !We also need to zero the nlink part of dxyzqm so we can accumulate link atom forces in this routine.
            qmmm_struct%dxyzqm(1:3, qmmm_struct%nquant + 1:qmmm_struct%nquant_nlink) = zero

        end if
!============ END QMMM =======================================================

        epol = zero
        eelt = zero
        evdw = zero
        esurf = zero
        totsasa = zero
#ifdef LES
        elesp = zero
#endif
        onstep = mod(irespa, nrespa) == 0
        onstepi = mod(irespa, nrespai) == 0
        if (.not. onstepi) return

        oncpstep = icnstph /= 0 .and. mod(irespa, ntcnstph) == 0

        if (alpb == 1) then
!     Sigalov Onufriev ALPB (epsilon-dependent GB):
            alpb_beta = alpb_alpha*(intdiel/extdiel)
            extdieli = one/(extdiel*(one + alpb_beta))
            intdieli = one/(intdiel*(one + alpb_beta))
            one_Arad_beta = alpb_beta/Arad
            if (kappa /= zero) onekappa = one/kappa
        else
            !  Standard Still's GB - alpb=0
            extdieli = one/extdiel
            intdieli = one/intdiel
            one_Arad_beta = zero
        end if

        ! Smooth "cut-off" in calculating GB effective radii.
        ! Implementd by Andreas Svrcek-Seiler and Alexey Onufriev.
        ! The integration over solute is performed up to rgbmax and includes
        ! parts of spheres; that is an atom is not just "in" or "out", as
        ! with standard non-bonded cut.  As a result, calclated effective
        ! radii are less than rgbmax. This saves time, and there is no
        ! discontinuity in dReff/drij.

        ! Only the case rgbmax > 5*max(sij) = 5*fsmax ~ 9A is handled; this is
        ! enforced in mdread().  Smaller values would not make much physical
        ! sense anyway.

        rgbmax2 = rgbmax*rgbmax
        rgbmax1i = one/rgbmax
        rgbmax2i = rgbmax1i*rgbmax1i
        rgbmaxpsmax2 = (rgbmax + fsmax)**2

#ifdef LES
! initialize some things for GB+LES

! one over number of LES copies
! GB+LES only works with 1 LES region so we can't have multiple
! copy numbers like we can with PME

        lfac = float(ncopy)
        lfaci = one/(lfac)

#endif

        !---------------------------------------------------------------------------
        !      Step 1: loop over pairs of atoms to compute the effective Born radii.
        !---------------------------------------------------------------------------
        !
        ! The effective Born radii are now calculated via a call at the
        ! beginning of force.
        !
        !---------------------------------------------------------------------------

        iexcl = 1 !moved to outside the index loop from original location in "step 2"

#ifdef MPI
        do i = 1, mytaskid
            iexcl = iexcl + numex(i)
        end do
        mpistart = mytaskid + 1
#endif

        maxi = natom
        if (natbel > 0) maxi = natbel

        !--------------------------------------------------------------------------
        !
        !     Step 2: loop over all pairs of atoms, computing the gas-phase
        !             electrostatic energies, the LJ terms, and the off-diagonal
        !             GB terms.  Also accumulate the derivatives of these off-
        !             diagonal terms with respect to the inverse effective radii,
        !             sumdeijda(k) will hold  sum over i,j>i ( deij/dak ),  where
        !             "ak" is the inverse of the effective radius for atom "k".
        !
        !             Update the forces with the negative derivatives of the
        !             gas-phase terms, plus the derivatives of the explicit
        !             distance dependence in Fgb, i.e. the derivatives of the
        !             GB energy terms assuming that the effective radii are
        !             constant.
        !
        !--------------------------------------------------------------------------

#ifdef LES
        sumdeijda(1:natom*ncopy) = zero
#else
        sumdeijda(1:natom) = zero
#endif

        call timer_start(TIME_GBFRC)

        !       Note: this code assumes that the belly atoms are the first natbel
        !             atoms...this is checked in mdread.

#ifdef MPI
        do i = mpistart, maxi, numtasks
#else

        do i = 1, maxi
#endif

#ifdef LES
            lestmp = nlesty*(lestyp(i) - 1)
#endif

            xi = x(3*i - 2)
            yi = x(3*i - 1)
            zi = x(3*i)
            qi = charge(i)
            ri = reff(i)
            four_ri = four*reff(i)
            iaci = ntypes*(iac(i) - 1)
            jexcl = iexcl
            jexcl_last = iexcl + numex(i) - 1
            dumx = zero
            dumy = zero
            dumz = zero

#if defined(LES)
            icnum = cnum(i)
            nrg_vdw_tmp = zero
            nrg_ele_tmp = zero
            nrg_egb_tmp = zero
#endif

            !         -- check the exclusion list for eel and vdw:

            do k = i + 1, natom
                skipv(k) = .false.
            end do
            do jjv = jexcl, jexcl_last
                skipv(natex(jjv)) = .true.
            end do

            ! QMMM:  We have 2 lists here:
            !        atom_mask which is natom long and .true. for each qm atom
            !        skipv which is natom long and is .true. if atom v should
            !             be skipped for this atom

            !Step 1 - is current atom i a QM atom - NOT link atoms?

            if (qmmm_nml%ifqnt) then
                if (qmmm_struct%atom_mask(i)) then  !yes it is
                    !step 2 - exclude all other QM atoms
                    do qm_temp_count = 1, qmmm_struct%nquant
                        skipv(qmmm_struct%iqmatoms(qm_temp_count)) = .true.
                        !The VDW between the link atom and the MM atoms is not
                        !calculated. Instead MM link pair atom with Real QM atoms is calculated.
                    end do
                    !Is atom i an MM link pair atom?
                else if (qmmm_struct%mm_link_mask(i)) then !yes it is
                    !Exclude all interactions (eel and VDW) with this
                    !atom. We will do the MM link pair VDW terms manually
                    !at the end of this routine.
                    skipv(i + 1:natom) = .true. !We can't just cycle here because we want the GB interactions.
                end if
            end if

            ! END QMMM

            icount = 0
            do j = i + 1, natom

#ifdef LES
! skip pairs in different copies
                if (icnum .ne. 0 .and. cnum(j) .ne. 0 .and. icnum .ne. cnum(j)) cycle
#endif

                xij = xi - x(3*j - 2)
                yij = yi - x(3*j - 1)
                zij = zi - x(3*j)
                r2 = xij*xij + yij*yij + zij*zij

                if (r2 <= cut .and. (onstep .or. r2 <= cut_inner)) then
#ifdef LES
! set up the loop over (possible) multiple effective radii for i and j

                    ijcnum = icnum + cnum(j)

! tempscale is for normalizing the epol when we average over pairs of reff for i and j
! default to only a single reff for each.
! use a temp value because this may be assigned to several icount values if we have a loop
! also default to no scaling factor for LES:LES pairs
! tempscale is not one only for non-les:non-les pairs with multiple radii
! and only used for epol.
! templfac is not one only for les:les pairs and used for epol, eel and evdw
!
                    tempscale = one
                    templfac = one

                    spread = .false.

                    if (ijcnum == 0) then

! both atoms are non-LES, see if we need to average

                        if (nradii(i) .gt. 1 .or. nradii(j) .gt. 1) then

! neither are LES, average over all ncopy pairs of reff
! do NOT average over ncopy*ncopy pairs - only do pairs
! with the same index.
! do average even if only 1 has multiple reff, nothing to
! be saved by using only the one reff for 1 atom if the other
! has multiple reff

                            idx1 = 1
                            idx2 = ncopy
                            tempscale = lfaci
                        else

! no need to average, use first reff for both (since all are average reff)
! here is where we need to deal with the sumdeijda problem
! if we use only 1 radius pair, we still need to divide the
! interaction among the radii for the calculation of derivatives

                            idx1 = 1
                            idx2 = 1

! set a flag saying whether the force for this interaction needs
! to be spread across all of the sumdeijda for the atom
! even if we use 1 reff, we need the force to arise from all atoms
! or else the forces on the LES atoms will not be balanced
! we only do this for non-LES:non-LES pairs

                            spread = .true.
                        end if
                    else

! at least 1 is LES, no averaging, only use the cnum reff
! for LES:LES, use the only reff that these have (cnum)
! for LES:non-LES, use the reff for the non-LES atom that goes with this LES copy (cnum)
!
                        if (icnum .gt. 0) then

! i is LES (j may be too), use i's cnum (same as j's if j is LES or
! else we would have skipped pair)

                            idx1 = icnum
                            idx2 = icnum

! set les scaling fac for ele, vdw if both are LES

                            if (cnum(i) .eq. cnum(j)) templfac = lfac
                        else

! i is not LES, j must be, so use j's cnum

                            idx1 = cnum(j)
                            idx2 = cnum(j)
                        end if
                    end if

                    istrt = ncopy*(i - 1)
                    jstrt = ncopy*(j - 1)

                    first = .true.
                    do k = idx1, idx2

! non-LES uses just ri (reff(i)) in this loop, doesn't need to be vector

                        rix(icount + 1) = reff(istrt + k)

                        reff_j = reff(jstrt + k)

! set longskipv, which serves the purpose of telling if an atom is
! excluded and also telling if this is not the first in the loop over

! pairs of reff for an atom (only calculate nonbonds for the first pair of reff)
! longskipv is only set for the icount atoms! (not all j like skipv)

                        if ((.not. first) .or. skipv(j)) then
                            longskipv(icount + 1) = .true.
                        else
                            longskipv(icount + 1) = .false.
                        end if

! for LES we need to know which reff this is when we calculate the sumdeija
! so set a pointer here to tell us which reff this i,j pair are using
! this points directly into correct spot in reff or sumdeijda

                        iridx(icount + 1) = istrt + k
                        jridx(icount + 1) = jstrt + k

! set flag to tell whether sumdeijda needs to be dividing among components for non-LES atoms
! temporary value was set above

                        spreadfrc(icount + 1) = spread

#else
! set for non-LES case so we can refer to reff_j not reff(j)
! this is needed with LES since each atom j has multiple reff

                    reff_j = reff(j)
#endif

                    icount = icount + 1
                    temp_jj(icount) = j
                    r2x(icount) = r2

! carlos changed this, we are STILL INSIDE LOOP OVER K for multiple LES reff_j
! so we need to use reff_j not reff(j) (since LES may have multiple reff_j for atom j)

                    rjx(icount) = reff_j

#ifdef LES
                    ! set scaling factor here where we know if we are doing
                    ! averaging over multiple radii. in the loops below we don't
                    ! have that info anymore without checking nradii
                    ! we couldn't set these above since we didn't have an
                    ! icount value yet

                    scalefac(icount) = tempscale
                    lesscalefac(icount) = templfac

                    ! set first to F so that only the first of these i,j
                    ! entries will have nonbonds

                    first = .false.

                end do ! LOOP OVER K, pairs of reff for i and j
#endif
            end if !r2 <= cut
        end do

        if (igb /= 6) then

#ifdef LES
            ! rix needs to be used instead of ri since i doesn't have just 1 reff
            vectmp1(1:icount) = four*rix(1:icount)*rjx(1:icount)
#else
            vectmp1(1:icount) = four_ri*rjx(1:icount)
#endif

            call vdinv(icount, vectmp1, vectmp1) !Invert things
            vectmp1(1:icount) = -r2x(1:icount)*vectmp1(1:icount)
            call vdexp(icount, vectmp1, vectmp1)
            ! [ends up with Exp(-rij^2/[4*ai*aj])]
#ifdef LES
            ! ri is not the same for all of the icount!
            vectmp3(1:icount) = r2x(1:icount) + &
                rjx(1:icount)*rix(1:icount)*vectmp1(1:icount) !ends up with fij
#else
            vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount)*ri*vectmp1(1:icount)
            ! [ends up with fij]
#endif

            ! CARLOS: LES
            ! we now have two sets of vectors: one is 1 to last atom, the other
            ! is 1 to icount
            ! the atom one is for vdw/ele and for excl, the other is for gb offdiag
            ! which do we loop over? unless we have a pointer we need to loop
            ! over the bigger one, have a pointer j in it, and pull excl out
            ! of that. but how to know to skip the nonbonds for all after the
            ! first in a loop over a pair? maybe we should have the excl loop
            ! done for icount, not atoms.
            ! SKIPV PART OK, DO THE DISTANCE PART SO WE DON'T HAVE TO INV ALL RIJ
            ! FOR NOW IT'S OK, JUST SLOW

            call vdinvsqrt(icount, vectmp3, vectmp2) !vectmp2 = 1/fij

            if (kappa /= zero) then
                call vdinv(icount, vectmp2, vectmp3)
                vectmp3(1:icount) = -kappa*vectmp3(1:icount)
                call vdexp(icount, vectmp3, vectmp4) !exp(-kappa*fij)
            end if

        end if
        call vdinvsqrt(icount, r2x, vectmp5) !1/rij

        !distance to mm_link pairs not to QM link atoms.

        !---- Start first outer loop ----
!!! !dir$ ivdep

        do k = 1, icount
            j = temp_jj(k)

            xij = xi - x(3*j - 2)
            yij = yi - x(3*j - 1)
            zij = zi - x(3*j)
            r2 = r2x(k)
            qiqj = qi*charge(j)
            if (igb /= 6) then

                if (kappa == zero) then
                    fgbk = zero
                    expmkf = extdieli
                else
                    expmkf = vectmp4(k)*extdieli

                    fgbk = vectmp3(k)*expmkf !-kappa*fij*exp(-kappa*fij)/Eout
                    if (alpb == 1) then ! Sigalov Onufriev ALPB:
                        fgbk = fgbk + (fgbk*one_Arad_beta*(-vectmp3(k)*onekappa))

                        ! Note: -vectmp3(k)*onekappa = fij
                    end if
                end if
                dl = intdieli - expmkf

                fgbi = vectmp2(k)  !1/fij

#ifdef LES
                ! epol gets scaled with scalefac
                ! the epol gets scaled for multiple reff pairs for i/j, but
                ! ele is NOT scaled since the extras are skipped using longskipv)
                ! the epol for LES-LES pairs gets scaled with lfac, and ele does too
                ! so we can't mix the lfac and the reff pair scaling in a
                ! single array; also can't pre-scale qiji since it is sometimes
                ! used with scaling and sometimes not.

                temp1 = -dl*(fgbi + one_Arad_beta)*scalefac(k)*lesscalefac(k)
#else
                temp1 = -dl*(fgbi + one_Arad_beta)
#endif
                e = qiqj*temp1
                if (ifcr /= 0) then
                    call cr_add_dcdr_factor(i, temp1*charge(j))
                    call cr_add_dcdr_factor(j, temp1*qi)
                end if
                epol = epol + e

#ifdef LES
                ! icnum is independant to j loop, so we cache the nrg here
                nrg_egb_tmp = nrg_egb_tmp + e
#endif

                if (idecomp == 1 .or. idecomp == 2) then
                    call decpair(1, i, j, e)
                else if (idecomp == 3 .or. idecomp == 4) then
                    call decpair(-1, i, j, e)
                end if
#ifdef MPI
#  ifdef LES
                if (rem == 2) then
                    if (cnum(i) > 0 .or. cnum(j) > 0) then
                        elesp = elesp + e
                    end if
                end if
#  endif
#endif
                if (oncpstep) then
                    dvdl = dvdl - (dl*fgbi*dcharge(i)*dcharge(j) + e)
                end if
                temp4 = fgbi*fgbi*fgbi !1/fij^3

                !   [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]

#ifdef LES
                temp6 = -qiqj*temp4*(dl + fgbk)*scalefac(k)*lesscalefac(k)
#else
                temp6 = -qiqj*temp4*(dl + fgbk)
#endif

                temp1 = vectmp1(k) !Exp(-rij^2/[4*ai*aj])
                de = temp6*(one - fourth*temp1)

                rj = rjx(k)
#ifdef LES
                ! modified to use rix vector entry instead of ri, since i
                !has multiple reff
                temp5 = half*temp1*temp6*(rix(k)*rj + fourth*r2)
#else
                temp5 = half*temp1*temp6*(ri*rj + fourth*r2)
#endif

#ifdef LES
                ! use iridx and jridx to get the pointers into the longer
                ! sumdeijda arrays.  Need this for LES since we need to know
                ! which of i's reff this j uses

                ! check to see if the forces need to be spread across multiple
                ! reff for this atom pair. This would happen if they are both
                !non-LES and only a single reff pair was used

                if (spreadfrc(k)) then

                    ! distribute the sumdeijda over multiple reff for these atoms
                    ! iridx and jridx will be pointing to the first reff for each atom

                    addi = rix(k)*temp5*lfaci
                    addj = rj*temp5*lfaci

                    do icopy = 0, ncopy - 1
                        sumdeijda(iridx(k) + icopy) = sumdeijda(iridx(k) + icopy) + addi
                        sumdeijda(jridx(k) + icopy) = sumdeijda(jridx(k) + icopy) + addj
                    end do

                else

                    ! only use this for this particular set of reff

                    sumdeijda(iridx(k)) = sumdeijda(iridx(k)) + rix(k)*temp5
                    sumdeijda(jridx(k)) = sumdeijda(jridx(k)) + rj*temp5
                end if

#else
                sumdeijda(i) = sumdeijda(i) + ri*temp5
                sumdeijda(j) = sumdeijda(j) + rj*temp5
#endif

                !    -- skip exclusions for remaining terms:
            else
                de = zero
            end if !if( igb/=6 )

#ifdef LES
            ! LES uses longskipv, which is for the icount atoms, and
            ! is T for the excluded atoms as well as the icount pairs
            ! that are not the first in the loop over reff pairs for i,j
            ! note that the index for longskipv is not the same as for skipv
            ! since we need multiple values for each i,j pair

            if (.not. longskipv(k)) then
#else
            if (.not. skipv(j)) then
#endif

                !   -- gas-phase Coulomb energy:

                ! Note: With QMMM we don't need to explicitly exclude QM-MM
                !       interactions here since the QM charges should
                !       be zero.  Unless qmgb/=0 in which case they are not
                !       zero so we need to exclude them.
                ! We also need to exclude the QM link atom interactions.

                doeel = .true.
                dovdw = .true.
                !  check if i or j is a quantum atom:
                if (qmmm_nml%ifqnt) then
                    if (qmmm_struct%atom_mask(i) .OR. qmmm_struct%atom_mask(j)) then
                        doeel = .false.
                    end if

                    ! if i or j is an MM link pair atom we also need to exclude it.
                    ! We also need to skip the VDW term if either is a mm link
                    ! pair atom
                    if (qmmm_struct%mm_link_mask(i) .or. qmmm_struct%mm_link_mask(j)) then
                        doeel = .false.
                        dovdw = .false.
                    end if
                end if

                !we can use the cached values.
                rinv = vectmp5(k) !1/rij
                r2inv = rinv*rinv

                if (doeel) then
#ifdef LES
                    temp1 = intdieli*rinv*lesscalefac(k)
#else
                    temp1 = intdieli*rinv
#endif
                    eel = qiqj*temp1
                    if (ifcr /= 0) then
                        call cr_add_dcdr_factor(i, temp1*charge(j))
                        call cr_add_dcdr_factor(j, temp1*qi)
                    end if
                    deel = eel*r2inv
                else
                    eel = zero
                    deel = zero
                end if

                eelt = eelt + eel
#ifdef LES
                nrg_ele_tmp = nrg_ele_tmp + eel
#endif
                if (idecomp == 1 .or. idecomp == 2) then
                    call decpair(2, i, j, eel)
                else if (idecomp == 3 .or. idecomp == 4) then
                    call decpair(-2, i, j, eel)
                end if
                de = de + deel
                if (oncpstep) then
                    dvdl = dvdl + (intdieli*rinv*dcharge(i) &
                        *dcharge(j) - eel)

                end if

                !    -- van der Waals energy:

                ic = ico(iaci + iac(j))
                if (ic > 0 .and. dovdw) then
                    !                                    6-12 potential:
                    r6inv = r2inv*r2inv*r2inv
#ifdef LES
                    ! scale LES:LES interactions
                    f6 = cn2(ic)*r6inv*lesscalefac(k)
                    f12 = cn1(ic)*(r6inv*r6inv)*lesscalefac(k)
#else
                    f6 = cn2(ic)*r6inv
                    f12 = cn1(ic)*(r6inv*r6inv)
#endif /*LES*/
                    evdw = evdw + (f12 - f6)
                    if (idecomp == 1 .or. idecomp == 2) then
                        call decpair(3, i, j, f12 - f6)
                    else if (idecomp == 3 .or. idecomp == 4) then
                        call decpair(-3, i, j, f12 - f6)
                    end if
                    de = de + (twelve*f12 - six*f6)*r2inv

#ifdef HAS_10_12

                    !    ---The following could be commented out if the Cornell
                    !       et al. force field was always used, since then all hbond
                    !       terms are zero.

                else
                    !                                    10-12 potential:
                    r10inv = r2inv*r2inv*r2inv*r2inv*r2inv
                    f10 = bsol(-ic)*r10inv
                    f12 = asol(-ic)*r10inv*r2inv
                    evdw = evdw + f12 - f10
                    if (idecomp == 3 .or. idecomp == 2) then
                        call decpair(1, i, j, f12 - f10)
                    else if (idecomp == 3 .or. idecomp == 4) then
                        call decpair(-1, i, j, f12 - f10)
                    end if
                    de = de + (twelve*f12 - ten*f10)*r2inv
#endif

                end if  ! ( ic > 0 )
#ifdef MPI
#  ifdef LES
                if (rem == 2) then
                    if (cnum(i) > 0 .or. cnum(j) > 0) then
                        elesp = elesp + eel + f12 - f6
                    end if
                end if
#  endif
#endif
            end if  ! ( .not. skipv(j) )

            !    -- derivatives:
            if (onstep .and. r2 > cut_inner) then
                de = de*nrespa
            else
                de = de*nrespai
            end if
            dedx = de*xij
            dedy = de*yij
            dedz = de*zij
            dumx = dumx + dedx
            dumy = dumy + dedy
            dumz = dumz + dedz

            if (qmmm_nml%ifqnt) then
                !Currently if j is a link atom then
                !dedxy,y,z are force on the link atom, not on the MM link pair.
                !we have to use the chain rule to put the forces back onto the MM link
                !pair atom.
                !For the moment we just accumulate the forces on the end of dxyzqm
                !this should have been zeroed at the beginning of this routine.
                if (qmmm_struct%mm_link_mask(j)) then
                    !j is a link atom, instead of adding it to the main force
                    !array, accumulate it in dxyzqm
                    !Find the link id for this j
                    do qm_temp_count = 1, qmmm_struct%nlink
                        if (qmmm_struct%link_pairs(1, qm_temp_count) == j) exit
                    end do
                    qm_temp_count = qm_temp_count + qmmm_struct%nquant
                    qmmm_struct%dxyzqm(1, qm_temp_count) = qmmm_struct%dxyzqm(1, qm_temp_count) + dedx
                    qmmm_struct%dxyzqm(2, qm_temp_count) = qmmm_struct%dxyzqm(2, qm_temp_count) + dedy
                    qmmm_struct%dxyzqm(3, qm_temp_count) = qmmm_struct%dxyzqm(3, qm_temp_count) + dedz
                else
                    !Not a link atom, can just add to main force array.
                    f(3*j - 2) = f(3*j - 2) - dedx
                    f(3*j - 1) = f(3*j - 1) - dedy
                    f(3*j) = f(3*j) - dedz
                end if
            else
                f(3*j - 2) = f(3*j - 2) - dedx
                f(3*j - 1) = f(3*j - 1) - dedy
                f(3*j) = f(3*j) - dedz
            end if
        end do !k=1,icount

#ifdef LES
        if (ipimd > 0) then
            nrg_all(icnum) = nrg_all(icnum) + nrg_vdw_tmp
            nrg_all(icnum) = nrg_all(icnum) + nrg_ele_tmp
            nrg_all(icnum) = nrg_all(icnum) + nrg_egb_tmp
        end if
#endif

        !---- End first outer loop ----

        if (qmmm_nml%ifqnt) then
            if (qmmm_struct%mm_link_mask(i)) then
                do qm_temp_count = 1, qmmm_struct%nlink
                    if (qmmm_struct%link_pairs(1, qm_temp_count) == i) exit
                end do
                qm_temp_count = qm_temp_count + qmmm_struct%nquant
                qmmm_struct%dxyzqm(1, qm_temp_count) = qmmm_struct%dxyzqm(1, qm_temp_count) - dumx
                qmmm_struct%dxyzqm(2, qm_temp_count) = qmmm_struct%dxyzqm(2, qm_temp_count) - dumy
                qmmm_struct%dxyzqm(3, qm_temp_count) = qmmm_struct%dxyzqm(3, qm_temp_count) - dumz
            else
                !Not a link atom, can just add to main force array.
                f(3*i - 2) = f(3*i - 2) + dumx
                f(3*i - 1) = f(3*i - 1) + dumy
                f(3*i) = f(3*i) + dumz
            end if
        else
            f(3*i - 2) = f(3*i - 2) + dumx
            f(3*i - 1) = f(3*i - 1) + dumy
            f(3*i) = f(3*i) + dumz
        end if
#ifdef MPI
        do k = i, (min(i + numtasks - 1, natom))
            iexcl = iexcl + numex(k)
        end do
#else
        iexcl = iexcl + numex(i)
#endif
    end do  !  i=1,maxi
    call timer_stop(TIME_GBFRC)

    if (igb == 6) goto VACUUM3

    !--------------------------------------------------------------------------
    !
    !    Step 3:  Finally, do the reduction over the sumdeijda terms:, adding
    !             into the forces those terms that involve derivatives of
    !             the GB terms (including the diagonal or "self" terms) with
    !             respect to the effective radii.  This is done by computing
    !             the vector dai/dxj, and using the chain rule with the
    !             previously-computed sumdeijda vector.
    !
    !             Also, compute a surface-area dependent term if igbsa=1
    !
    !             Do these terms only at "nrespa" multiple-time step intervals;
    !             (when igb=2 or 5, one may need to do this at every step)
    !
    !--------------------------------------------------------------------------

    if (onstep) then
        count = 0
        frespa = nrespa
#ifdef MPI

        !       -- first, collect all the sumdeijda terms:

        call timer_start(TIME_GBRADDIST)
        if (numtasks > 1) then
#ifdef LES
            k = natom*ncopy
#else
            k = natom
#endif
            ! carlos changed this to use k as set above, not natom
#ifdef USE_MPI_IN_PLACE
            call mpi_allreduce(MPI_IN_PLACE, sumdeijda, k, &
                MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#else
            call mpi_allreduce(sumdeijda, vectmp1, k, &
                MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
            sumdeijda(1:k) = vectmp1(1:k)
#endif
        end if
        call timer_stop(TIME_GBRADDIST)
#endif
        call timer_start(TIME_GBRAD2)

        !  -- diagonal epol term, plus off-diag derivs wrt alpha == reff^-1:

#ifdef MPI

        do i = mpistart, maxi, numtasks
#else
        do i = 1, maxi
#endif

            f_xi = zero
            f_yi = zero
            f_zi = zero
            qi = charge(i)
#ifdef LES
            icnum = cnum(i)
            nrg_egb_tmp = zero
#endif

#ifdef LES
            ! here we need to loop over possible multiple effective radii for i
            ! and average the energy/force from each reff

            self_e = zero

            icnum = cnum(i)
            istrt = ncopy*(i - 1)
            if (icnum .eq. 0) then

                ! non-LES atom, average over multiple reff by dividng by lfac

                do k = 1, ncopy
                    expmkf = exp(-kappa*reff(istrt + k))*extdieli
                    dl = intdieli - expmkf
                    qi2h = half*qi*qi
                    qid2h = qi2h*dl
                    temp1 = (onereff(istrt + k) + one_Arad_beta)*lfaci
                    self_e = self_e + qid2h*temp1
                    if (ifcr /= 0) &
                        call cr_add_dcdr_factor(i, -qi*temp1*dl)

                    ! add in contribution with scaling factor for average over
                    ! reff : DECREASE using lfaci

                    vtemp7(k) = -sumdeijda(istrt + k) &
                        + qid2h*lfaci - kappa*qi2h*lfaci*expmkf*reff(istrt + k)* &
                        (one + one_Arad_beta*reff(i))
                end do
            else

                ! LES atom: charge is scaled by N copies, so q^2 is scaled N^2.
                ! We need to scale UP by ncopy;  doesn't matter what copy # it
                ! is for scaling factor, but copy # determines reff index

                expmkf = exp(-kappa*reff(istrt + icnum))*extdieli
                dl = intdieli - expmkf
                qi2h = half*qi*qi
                qid2h = qi2h*dl
                temp1 = (onereff(istrt + icnum) + one_Arad_beta)*lfac
                self_e = qid2h*temp1
                if (ifcr /= 0) &
                    call cr_add_dcdr_factor(i, -qi*temp1*dl)

                ! add in contribution with scaling factor for q^2 :
                ! INCREASE using lfac

                vtemp7(icnum) = -sumdeijda(istrt + icnum) &
                    + qid2h*lfac - kappa*qi2h*lfac*expmkf*reff(istrt + icnum)* &
                    (one + one_Arad_beta*reff(i))
            end if

#else /*not LES*/

            expmkf = exp(-kappa*reff(i))*extdieli
            dl = intdieli - expmkf
            qi2h = half*qi*qi
            qid2h = qi2h*dl

            temp1 = (onereff(i) + one_Arad_beta)
            self_e = qid2h*temp1
            if (ifcr /= 0) &
                call cr_add_dcdr_factor(i, -qi*temp1*dl)

            temp7 = -sumdeijda(i) + qid2h - &
                kappa*qi2h*expmkf*reff(i)*(one + one_Arad_beta*reff(i))
#endif /*LES*/

            !Ross Walker
            epol = epol - self_e

#ifdef LES
            nrg_egb_tmp = nrg_egb_tmp - self_e
#endif
            if (idecomp == 1 .or. idecomp == 2) then
                call decpair(1, i, i, -qid2h*onereff(i))
            else if (idecomp == 3 .or. idecomp == 4) then
                call decpair(-1, i, i, -qid2h*onereff(i))
            end if
#ifdef MPI
#  ifdef LES
            ! local REMD
            if (rem == 2) then
                if (cnum(i) > 0) then
                    elesp = elesp - qid2h*onereff(i)
                    if (ifcr /= 0) &
                        call cr_add_dcdr_factor(i, -qi*dl*onereff(i))
                end if
            end if
#  endif
#endif
            if (oncpstep) then
                dvdl = dvdl - (half*dl*dcharge(i)*dcharge(i) - qid2h)*onereff(i)
            end if

            ! temp7 = ... + qiqj/2*[1/Ein - exp[-kappa*effbornrad]/Eout]
            !         -kappa*qiqj/2*exp[-kappa*effbornrad]/Eout*effbornrad
            ! temp7 without the -sumdeijda part is the diagonal gradient.

            ! carlos: moved temp7 calculation up into the LES region above

            xi = x(3*i - 2)
            yi = x(3*i - 1)
            zi = x(3*i)
            ri = rborn(i) - offset
            ri1i = one/ri
            iaci = ntypes*(iac(i) - 1)

            if (igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) then

                !  --- new onufriev: we have to later scale values by a
                !      alpha,beta,gamma -dependent factor:

                rborn_i = rborn(i)
                ri = rborn_i - offset

#ifdef LES
                if (icnum == 0) then
                    if (nradii(i) == 1) then

                        ! get single reff and calc like non-les code, set flag
                        ! for distributing sumdei again?

                        spread = .true.
                        psi_i = psi((i - 1)*ncopy + 1)
                        thi = tanh((gbalpha + gbgamma*psi_i*psi_i - gbbeta*psi_i)*psi_i)
                        thi2 = (gbalpha + three*gbgamma*psi_i*psi_i &
                            - two*gbbeta*psi_i)*(one - thi*thi)*ri/rborn_i
                    else

                        ! loop over multiple psi values, store vector of thi2 values

                        do k = 1, ncopy
                            psi_i = psi((i - 1)*ncopy + k)
                            thi = tanh((gbalpha + gbgamma*psi_i*psi_i - gbbeta*psi_i)*psi_i)

                            ! key change, thi2 is a vector due to multiple psi and
                            ! needing the thi2 below

                            vthi2(k) = (gbalpha + three*gbgamma*psi_i*psi_i &
                                - two*gbbeta*psi_i)*(one - thi*thi)*ri/rborn_i
                            spread = .false.
                        end do
                    end if

                else !cnum>0, LES atom

                    ! use the icnum psi value

                    spread = .false.
                    psi_i = psi((i - 1)*ncopy + icnum)
                    thi = tanh((gbalpha + gbgamma*psi_i*psi_i - gbbeta*psi_i)*psi_i)
                    thi2 = (gbalpha + three*gbgamma*psi_i*psi_i &
                        - two*gbbeta*psi_i)*(one - thi*thi)*ri/rborn_i
                end if
#else

                psi_i = psi(i)
                gba_i = gbvalpha(i) ! Hai Nguyen: take gbalpha_i,gbbeta_i,gbgamma_i
                gbb_i = gbvbeta(i)
                gbg_i = gbvgamma(i)
                thi = tanh((gba_i + gbg_i*psi_i*psi_i - gbb_i*psi_i)*psi_i)
                thi2 = (gba_i + three*gbg_i*psi_i*psi_i - &
                    two*gbb_i*psi_i)*(one - thi*thi)*ri/rborn_i
#endif
            end if

            icount = 0
            do j = 1, natom
                if (i /= j) then
#ifdef LES
                    jcnum = cnum(j)

                    ! skip if both LES but not same copy
                    if ((icnum .ne. 0 .and. jcnum .ne. 0) .and. icnum .ne. jcnum) cycle

#endif

                    xij = xi - x(3*j - 2)
                    yij = yi - x(3*j - 1)
                    zij = zi - x(3*j)
                    r2 = xij*xij + yij*yij + zij*zij
                    if (r2 <= rgbmaxpsmax2) then
                        ! pairlist contains only atoms within rgbmax + safety margin
                        icount = icount + 1
                        temp_jj(icount) = j
                        r2x(icount) = r2
                    end if ! r2 <= rgbmaxpsmax2
                end if !i/=j
            end do
            call vdinvsqrt(icount, r2x, vectmp1)

            kk1 = 0
            do k = 1, icount
                j = temp_jj(k)
                r2 = r2x(k)
                sj = fs(j)

                dij1i = vectmp1(k)
                dij = r2*dij1i
                sj2 = sj*sj

                if (dij <= four*sj) then
                    kk1 = kk1 + 1
                    vectmp3(kk1) = dij + sj
                    if (dij > ri + sj) then
                        vectmp2(kk1) = r2 - sj2
                        vectmp4(kk1) = dij - sj
                    else if (dij > abs(ri - sj)) then
                        vectmp2(kk1) = dij + sj
                        vectmp4(kk1) = ri
                    else if (ri < sj) then
                        vectmp2(kk1) = r2 - sj2
                        vectmp4(kk1) = sj - dij
                    else
                        vectmp2(kk1) = one
                        vectmp4(kk1) = one
                    end if
                end if !dij <= four*sj
            end do

            call vdinv(kk1, vectmp2, vectmp2)
            call vdinv(kk1, vectmp3, vectmp3)
            vectmp4(1:kk1) = vectmp4(1:kk1)*vectmp3(1:kk1)
            call vdln(kk1, vectmp4, vectmp4)

            kk1 = 0
            do k = 1, icount
                j = temp_jj(k)
                j3 = 3*j
                r2 = r2x(k)
                xij = xi - x(j3 - 2)
                yij = yi - x(j3 - 1)
                zij = zi - x(j3)

                dij1i = vectmp1(k)
                dij = r2*dij1i
                sj = fs(j)
                if (dij <= rgbmax + sj) then
                    sj2 = sj*sj

                    !           datmp will hold (1/r)(dai/dr):

                    dij2i = dij1i*dij1i
                    dij3i = dij2i*dij1i

                    if (dij > rgbmax - sj) then

                        temp1 = 1.0d0/(dij - sj)
                        datmp = eighth*dij3i*((r2 + sj2)* &
                            (temp1*temp1 - rgbmax2i) - two*log(rgbmax*temp1))

                    else if (dij > four*sj) then

                        tmpsd = sj2*dij2i
                        dumbo = te + tmpsd*(tf + tmpsd*(tg + tmpsd*(th + tmpsd*thh)))
                        datmp = tmpsd*sj*dij2i*dij2i*dumbo

                        !     ---check accuracy of above Taylor series:

                    else if (dij > ri + sj) then

                        kk1 = kk1 + 1
                        datmp = vectmp2(kk1)*sj*(-half*dij2i + vectmp2(kk1)) + &
                            fourth*dij3i*vectmp4(kk1)

                    else if (dij > abs(ri - sj)) then
                        kk1 = kk1 + 1
                        datmp = -fourth*(-half*(r2 - ri*ri + sj2)*dij3i*ri1i*ri1i &
                            + dij1i*vectmp2(kk1)*(vectmp2(kk1) - dij1i) &
                            - dij3i*vectmp4(kk1))

                    else if (ri < sj) then
                        kk1 = kk1 + 1
                        datmp = -half*(sj*dij2i*vectmp2(kk1) &
                            - two*sj*vectmp2(kk1)*vectmp2(kk1) &
                            - half*dij3i*vectmp4(kk1))

                    else
                        kk1 = kk1 + 1
                        datmp = zero
                    end if  ! ( dij > 4.d0*sj )

                    if ((igb == 7 .or. igb == 8) .and. dij < rborn(i) + rborn(j) + GBNECKCUT) then

                        ! (no changes needed for LES)

                        ! Derivative of neck with respect to dij is:
                        !                     5
                        !              9 mdist
                        !   (2 mdist + --------) neckMaxVal gbneckscale
                        !                 5
                        ! -(------------------------)
                        !                        6
                        !             2   3 mdist  2
                        !   (1 + mdist  + --------)
                        !                    10

                        mdist = dij - neckMaxPos(neckidx(i), neckidx(j))
                        mdist2 = mdist*mdist
                        mdist3 = mdist2*mdist
                        mdist5 = mdist2*mdist3
                        mdist6 = mdist3*mdist3

                        ! temp1 will be divisor of above fraction * dij
                        ! (datmp is deriv * 1/r)

                        temp1 = 1 + mdist2 + (three/ten)*mdist6
                        temp1 = temp1*temp1*dij

                        ! (Note "+" means subtracting derivative, since above
                        !     expression has leading "-")

                        datmp = datmp + ((2*mdist + (nine/five)*mdist5) &
                            *neckMaxVal(neckidx(i), neckidx(j))*gbneckscale)/temp1
                    end if

#ifdef LES
                    ! loop over pairs of radii
                    ! SHOULD WE USE NRADII HERE?

                    jcnum = cnum(j)
                    if (icnum == 0) then
                        ! i is not a LES copy
                        if (jcnum == 0) then
                            ! j not LES either, calculate force using all pairs of reff
                            idx1 = 1
                            idx2 = ncopy
                        else
                            ! j is LES, add force using only one of the radii for
                            ! i and j (the one for j's cnum)
                            idx1 = jcnum
                            idx2 = jcnum
                        end if
                    else
                        ! i is LES, either j is not LES or is in same LES copy as i,
                        ! so add to force using i's cnum
                        idx1 = icnum
                        idx2 = icnum
                    end if

                    ! save datmp since the loop modifies it
                    tmp = datmp
                    do icopy = idx1, idx2

                        ! note that vtemp7 takes the place of temp7 used for non-LES
                        !simulations
                        datmp = -tmp*frespa*vtemp7(icopy)

                        if (igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) then

                            ! for nradii>1, we have to use the thi2 vector (vtemp7),
                            ! others have only 1 thi2 value

                            if (nradii(i) .gt. 1) then
                                datmp = datmp*vthi2(icopy)
                            else
                                datmp = datmp*thi2
                            end if
                        end if

                        f_x = xij*datmp
                        f_y = yij*datmp
                        f_z = zij*datmp
                        f(j3 - 2) = f(j3 - 2) + f_x
                        f(j3 - 1) = f(j3 - 1) + f_y
                        f(j3) = f(j3) + f_z
                        f_xi = f_xi - f_x
                        f_yi = f_yi - f_y
                        f_zi = f_zi - f_z
                    end do
#else /* not LES */

                    datmp = -datmp*frespa*temp7

                    if (igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) datmp = datmp*thi2
                    f_x = xij*datmp
                    f_y = yij*datmp
                    f_z = zij*datmp
                    if (qmmm_nml%ifqnt) then
                        if (qmmm_struct%mm_link_mask(j)) then
                            do qm_temp_count = 1, qmmm_struct%nlink
                                if (qmmm_struct%link_pairs(1, qm_temp_count) == j) exit
                            end do
                            qm_temp_count = qm_temp_count + qmmm_struct%nquant
                            qmmm_struct%dxyzqm(1, qm_temp_count) = qmmm_struct%dxyzqm(1, qm_temp_count) - f_x
                            qmmm_struct%dxyzqm(2, qm_temp_count) = qmmm_struct%dxyzqm(2, qm_temp_count) - f_y
                            qmmm_struct%dxyzqm(3, qm_temp_count) = qmmm_struct%dxyzqm(3, qm_temp_count) - f_z
                        else
                            !Not a link atom, can just add to main force array.
                            f(j3 - 2) = f(j3 - 2) + f_x
                            f(j3 - 1) = f(j3 - 1) + f_y
                            f(j3) = f(j3) + f_z
                        end if
                    else
                        f(j3 - 2) = f(j3 - 2) + f_x
                        f(j3 - 1) = f(j3 - 1) + f_y
                        f(j3) = f(j3) + f_z
                    end if

                    f_xi = f_xi - f_x
                    f_yi = f_yi - f_y
                    f_zi = f_zi - f_z

#endif /*LES*/

                end if
            end do

            if (qmmm_nml%ifqnt) then
                if (qmmm_struct%mm_link_mask(i)) then
                    do qm_temp_count = 1, qmmm_struct%nlink
                        if (qmmm_struct%link_pairs(1, qm_temp_count) == i) exit
                    end do
                    qm_temp_count = qm_temp_count + qmmm_struct%nquant
                    qmmm_struct%dxyzqm(1, qm_temp_count) = qmmm_struct%dxyzqm(1, qm_temp_count) - f_xi
                    qmmm_struct%dxyzqm(2, qm_temp_count) = qmmm_struct%dxyzqm(2, qm_temp_count) - f_yi
                    qmmm_struct%dxyzqm(3, qm_temp_count) = qmmm_struct%dxyzqm(3, qm_temp_count) - f_zi
                else
                    !Not a link atom, can just add to main force array.
                    f(3*i - 2) = f(3*i - 2) + f_xi
                    f(3*i - 1) = f(3*i - 1) + f_yi
                    f(3*i) = f(3*i) + f_zi
                end if
            else
                f(3*i - 2) = f(3*i - 2) + f_xi
                f(3*i - 1) = f(3*i - 1) + f_yi
                f(3*i) = f(3*i) + f_zi
            end if

            !  --- Define neighbor list ineighbor for calc of LCPO areas ---

            if (gbsa > 0) then
                do k = 1, icount
                    j = temp_jj(k)
                    dij = sqrt(r2x(k))
                    if ((vdwrad(i) + vdwrad(j)) > dij) then

                        ! Consider all atoms for icosa, only non-H's for LCPO:

                        if ((gbsa == 2) .or. &
                            (vdwrad(i) > 2.5) .and. (vdwrad(j) > 2.5)) then
                            count = count + 1
                            ineighbor(count) = j
                        end if
                    end if
                end do
                count = count + 1
                ineighbor(count) = 0
            end if

            if (gbsa == 2) then

                !  --- Calc surface area with icosasurf algo;
                !  does not provide forces ==> only for single point calculations

#ifdef MPI
                if (i == mpistart) then
                    ineighborpt = 1
                    call icosa_init(2, 3, zero)
                end if

#else
                if (i == 1) then
                    ineighborpt = 1
                    call icosa_init(2, 3, zero)
                end if
#endif

                totsasa = totsasa + icosa_sphere_approx(i, x, &
                    vdwrad, ineighborpt, ineighbor, idecomp)

            end if  !  ( gbsa == 2 )

#ifdef LES
            if (ipimd > 0) nrg_all(icnum) = nrg_all(icnum) + nrg_egb_tmp
#endif

        end do   ! end loop over atom i

        call timer_stop(TIME_GBRAD2)

        if (gbsa == 1) then

            !       --- calculate surface area by computing LCPO (Still) over
            !           all atoms ---

#        include "gbsa.h"

        end if  !  ( gbsa == 1 )
        esurf = surften*totsasa
    end if  !  i=mytaskid+1,maxi,numtasks

    VACUUM3 &
    !======== QMMM ==========
        if(qmmm_nml%ifqnt) then
    call timer_start(TIME_QMMM)
    if (qmmm_nml%qmgb /= 0) then

        ! We filled the main charge array with charges for QM atoms; make sure
        ! we zero it again so the 1-4's which are done outside of EGB will be
        ! correctly skipped for QM-MM on the next step. Note here we don't have the routine
        ! save the charges again since depending on the GB option they may have been filled
        ! with either RESP charges or mulliken charges.
        call qm_zero_charges(qmmm_nml, qmmm_struct, charge, qmmm_struct%scaled_mm_charges, .false.)
        if (qmmm_struct%nlink > 0) then
            !Don't save the charges here since they could be the resp charges.
            call qm_zero_mm_link_pair_main_chg(qmmm_nml, qmmm_struct, qmmm_struct%nlink, qmmm_struct%link_pairs, charge, &
                qmmm_struct%scaled_mm_charges, .false.)
        end if
    end if

    call timer_start(TIME_QMMMCOLLATEF)
    !We need to restore the MM link pair coordinates and then
    !Use the chain rule to put the link pair forces back onto the
    !QM and MM link pairs.
    call rst_mm_link_pair_crd(qmmm_nml, qmmm_struct, x)
    do i = 1, qmmm_struct%nlink
        mm_no = 3*qmmm_struct%link_pairs(1, i) - 2  !location of atom in x array
        lnk_no = qmmm_struct%link_pairs(2, i) !Nquant number of QM atom bound to link atom
        qm_no = 3*qmmm_struct%iqmatoms(lnk_no) - 2
        !Note this routine uses the flink in the form -flink.
        call distribute_lnk_f(forcemod, qmmm_struct%dxyzqm(1, qmmm_struct%nquant + i), x(mm_no), &
            x(qm_no), qmmm_nml%lnk_dis)

        !NOTE: forces are reversed in QM calc with respect to amber force array
        !so we subtract forcemod from MM atom and add it to QM atom.
        j = (qmmm_struct%link_pairs(1, i) - 1)*3 !Natom number of MM link pair.
        !MM atom's new force = FMM(x,y,z) - FORCEMOD(x,y,z)
        f(j + 1) = f(j + 1) - forcemod(1)
        f(j + 2) = f(j + 2) - forcemod(2)
        f(j + 3) = f(j + 3) - forcemod(3)

        j = (qmmm_struct%iqmatoms(lnk_no) - 1)*3
        !QM atom's new force = FQM(x,y,z) - Flink(x,y,z) + FORCEMOD(x,y,z)
        !Note QM forces should be subtracted from sander F array to leave total force.
        f(j + 1) = f(j + 1) - qmmm_struct%dxyzqm(1, qmmm_struct%nquant + i) + forcemod(1)
        f(j + 2) = f(j + 2) - qmmm_struct%dxyzqm(2, qmmm_struct%nquant + i) + forcemod(2)
        f(j + 3) = f(j + 3) - qmmm_struct%dxyzqm(3, qmmm_struct%nquant + i) + forcemod(3)
    end do
    call timer_stop_start(TIME_QMMMCOLLATEF, TIME_QMMMENERGY)
    call timer_start(TIME_QMMMGBENERGY)
    !Finally we need to calculate the VDW terms
    !for the MM link pair atoms as they were skipped above.

!------------ START MM LINK PAIR VDW -------------------
    do i = 1, qmmm_struct%nlink
        skipv(1:natom) = .false.
        jexcl = 1
        mm_no = qmmm_struct%link_pairs(1, i)
        iminus = mm_no - 1
        do j = 1, iminus
            !As we loop over each atom check if this atom is supposed to exclude this MM
            !link pair atom.
            jexcl_last = jexcl + numex(j) - 1
            do jjv = jexcl, jexcl_last
                iexcl = natex(jjv)
                !Should this atom exclude the MM atom?
                if (iexcl == mm_no) skipv(j) = .true.
            end do
            jexcl = jexcl + numex(j)
        end do
        jexcl_last = jexcl + numex(mm_no) - 1
        do jjv = jexcl, jexcl_last
            skipv(natex(jjv)) = .true.
        end do
        !We need to avoid double counting so we should exclude ourself and
        !all link atoms of lower number
        !Note this assumes the link atom list is numerically sorted...
        do k = 1, i !To i ensures we skip the self interaction.
            skipv(qmmm_struct%link_pairs(1, k)) = .true.
        end do
        xi = x(3*mm_no - 2)
        yi = x(3*mm_no - 1)
        zi = x(3*mm_no)

        iaci = ntypes*(iac(mm_no) - 1)
#ifdef LES
        icnum = cnum(mm_no)
        nrg_vdw_tmp = zero
#endif
        dumx = zero
        dumy = zero
        dumz = zero
        de = zero
#ifdef MPI
        do j = mpistart, maxi, numtasks
#else
        do j = 1, maxi !1 to natom
#endif
            !Do all atoms that are not excluded and
            !are within the cutoff for this i. Skip other
            !MM interactions that are less than i to avoid double
            !counting.
            if (.not. skipv(j)) then
                xij = xi - x(3*j - 2)
                yij = yi - x(3*j - 1)
                zij = zi - x(3*j)
                r2 = xij*xij + yij*yij + zij*zij
                if (r2 <= cut) then
                    r2inv = one/r2
                    ic = ico(iaci + iac(j))
                    if (ic > 0) then
                        !6-12 potential:
                        r6inv = r2inv*r2inv*r2inv
                        f6 = cn2(ic)*r6inv
                        f12 = cn1(ic)*(r6inv*r6inv)
                        evdw = evdw + (f12 - f6)

                        if (idecomp == 1 .or. idecomp == 2) then
                            call decpair(3, i, j, f12 - f6)
                        else if (idecomp == 3 .or. idecomp == 4) then
                            call decpair(-3, i, j, f12 - f6)
                        end if
                        de = (twelve*f12 - six*f6)*r2inv
#ifdef HAS_10_12
                    else
                        !10-12 potential:
                        r10inv = r2inv*r2inv*r2inv*r2inv*r2inv
                        f10 = bsol(-ic)*r10inv
                        f12 = asol(-ic)*r10inv*r2inv
                        evdw = evdw + f12 - f10
                        if (idecomp == 3 .or. idecomp == 2) then
                            call decpair(1, i, j, f12 - f10)
                        else if (idecomp == 3 .or. idecomp == 4) then
                            call decpair(-1, i, j, f12 - f10)
                        end if
                        de = (twelve*f12 - ten*f10)*r2inv
#endif
                    end if
                    de = de*nrespai
                    dedx = de*xij
                    dedy = de*yij
                    dedz = de*zij
                    dumx = dumx + dedx
                    dumy = dumy + dedy
                    dumz = dumz + dedz
                    f(3*j - 2) = f(3*j - 2) - dedx
                    f(3*j - 1) = f(3*j - 1) - dedy
                    f(3*j) = f(3*j) - dedz
                end if
            end if
        end do

#ifdef LES
        if (ipimd > 0) nrg_all(icnum) = nrg_all(icnum) + nrg_vdw_tmp
#endif
        f(3*mm_no - 2) = f(3*mm_no - 2) + dumx
        f(3*mm_no - 1) = f(3*mm_no - 1) + dumy
        f(3*mm_no) = f(3*mm_no) + dumz
    end do
    call timer_stop(TIME_QMMMGBENERGY)
    call timer_stop(TIME_QMMMENERGY)
!------------ END MM LINK PAIR VDW ---------------------
    call timer_stop(TIME_QMMM)
end if
!======== END QMMM ==========

return
end subroutine egb

subroutine egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
    rborn, offset, rbornstat, &
    rbave, rbfluct, rbmax, rbmin, gbneckscale, ncopy, rdt, &
    gbvalpha, gbvbeta, gbvgamma, &  ! Hai Nguyen: parameters for GB,
 ! x(l2404),x(l2405),x(l2406)
    gbalpha, gbbeta, gbgamma &      !Hai Nguyen: keep for sander.LES
#ifdef MPI
    , mpistart &
#endif
    )

!Hai Nguyen: don't have igb = 8 for sander.LES. Thus, instead of using gbavalpha,
!gbvbeta and gbvgamma array, I still keep original gbalpha, gbbeta, gbgamma
!(igb1,2,5,7) for LES part to get correct force

! Calculates effective GB radii and puts result in reff;
! onereff contains 1.0d0/reff.

    use constants, only : zero, eighth, fourth, third, half, &
        one, two, three, four, five, seven, eight, nine, &
        eleven, twelve, thirtieth

    implicit none

#include "gbneck.h"

#ifdef MPI
# include "parallel.h"
    include 'mpif.h'
    integer, intent(in) :: mpistart
#endif
    integer :: ierr

    integer, intent(in) :: igb, natom, rbornstat, ncopy
    _REAL_, intent(in) :: x(3*natom), fs(natom), rdt

#ifdef LES
# include "les.h"
    integer istrt, iend, k1
    integer idx1, idx2, jstrt, tmpstop
! copy number for atoms i and j
    integer icnum, jcnum
    _REAL_ rmin, rmax
    _REAL_, intent(out) :: reff(ncopy*natom), onereff(ncopy*natom)
#else
    _REAL_, intent(out) :: reff(natom), onereff(natom)
#endif

    _REAL_, intent(in) :: fsmax, rgbmax, gbneckscale
    _REAL_, intent(in) :: rborn(natom), offset, gbalpha, gbbeta, gbgamma

    !Hai Nguyen: add gbvalpha(natom), gbvbeta(natom), gbvgamma(natom) arrays
    _REAL_, intent(in) :: gbvalpha(natom), gbvbeta(natom), gbvgamma(natom)
    _REAL_, intent(out) :: rbave(natom), rbfluct(natom), rbmax(natom), rbmin(natom)

! Local:
    integer :: i, icount, iplus, kk1, kk2, k, j, vecend
    _REAL_ :: xi, yi, zi, xij, yij, zij, r2, reff_i, ri, ri1i, rj1i, si, si2, sj, sj2
    _REAL_ :: dij1i, dij2i, dij, rj, uij, tmpsd, dumbo, theta
    _REAL_ :: rgbmax1i, rgbmax2i, rgbmaxpsmax2, rborn_i, psi_i
    _REAL_ :: temp3, temp4

    !   FGB taylor coefficients follow
    !   from A to D :
    !   1/3 , 2/5 , 3/7 , 4/9 , 5/11
    _REAL_ ta
    _REAL_ tb
    _REAL_ tc
    _REAL_ td
    _REAL_ tdd
    parameter(ta=third)
    parameter(tb=two/five)
    parameter(tc=three/seven)
    parameter(td=four/nine)
    parameter(tdd=five/eleven)

    rgbmax1i = one/rgbmax
    rgbmax2i = rgbmax1i*rgbmax1i
    rgbmaxpsmax2 = (rgbmax + fsmax)**2

#ifdef LES
! need to expand range of onereff that are intialized
    onereff(1:natom*ncopy) = zero
! initialize rbornlong, see egb()
    k = 0
    do i = 1, natom
        do j = 1, ncopy
            k = k + 1
            rbornlong(k) = rborn(i)
        end do
    end do
#else
    onereff(1:natom) = zero
#endif

#ifdef MPI
    do i = mpistart, natom, numtasks
#else
    do i = 1, natom
#endif
        xi = x(3*i - 2)
        yi = x(3*i - 1)
        zi = x(3*i)

#ifdef LES
! copy # for atom i
        icnum = cnum(i)

! due to the multiple reff we will not use reff_i but loop over array directly
! pointers to entries in expanded reff array for atom i
! reff and onereff can't use atom index as pointer

        istrt = ncopy*(i - 1)
#else
        reff_i = onereff(i)
#endif

#ifdef LES
! these next variables are not affected by LES (rborn, fs)
#endif

        ri = rborn(i) - offset
        ri1i = one/ri
        si = fs(i)
        si2 = si*si

        !  Here, reff_i will sum the contributions to the inverse effective
        !  radius from all of the atoms surrounding atom "i"; later the
        !  inverse of its own intrinsic radius will be added in

        icount = 0
        iplus = i + 1

        do j = iplus, natom

#ifdef LES
            jcnum = cnum(j)
            ! LES atoms in different copies do not descreen each other
            if ((icnum .ne. 0 .and. jcnum .ne. 0) .and. icnum .ne. jcnum) cycle
#endif

            xij = xi - x(3*j - 2)
            yij = yi - x(3*j - 1)
            zij = zi - x(3*j)
            r2 = xij*xij + yij*yij + zij*zij

            if (r2 <= rgbmaxpsmax2) then
                icount = icount + 1
                temp_jj(icount) = j
                r2x(icount) = r2
            end if

        end do

        call vdinvsqrt(icount, r2x, vectmp1)

        kk1 = 0
        kk2 = 0
!!!      !dir$ ivdep
        do k = 1, icount

            j = temp_jj(k)
            r2 = r2x(k)
            sj = fs(j)

!        don't fill the remaining vectmp arrays if atoms don't see each other:
            dij1i = vectmp1(k) !1/sqrt(r^2)
            dij = r2*dij1i != rij
            if (dij <= rgbmax + si .or. dij <= rgbmax + sj) then
                rj = rborn(j) - offset

                if (dij <= four*sj) then
                    kk1 = kk1 + 1
                    vectmp2(kk1) = dij + sj
                    if (dij > ri + sj) then
                        vectmp4(kk1) = dij - sj
                    else if (dij > abs(ri - sj)) then
                        vectmp4(kk1) = ri
                    else if (ri < sj) then
                        vectmp4(kk1) = sj - dij
                    else
                        vectmp4(kk1) = one
                    end if
                end if

                if (dij <= four*si) then
                    kk2 = kk2 + 1
                    vectmp3(kk2) = dij + si
                    if (dij > rj + si) then
                        vectmp5(kk2) = dij - si
                    else if (dij > abs(rj - si)) then
                        vectmp5(kk2) = rj
                    else if (rj < si) then
                        vectmp5(kk2) = si - dij
                    else
                        vectmp5(kk2) = one
                    end if
                end if
            end if
        end do

#ifdef LES
! these arrays for vector calls are not changed by LES
#endif

        call vdinv(kk1, vectmp2, vectmp2)
        call vdinv(kk2, vectmp3, vectmp3)
        vectmp4(1:kk1) = vectmp2(1:kk1)*vectmp4(1:kk1)
        vectmp5(1:kk2) = vectmp3(1:kk2)*vectmp5(1:kk2)
        call vdln(kk1, vectmp4, vectmp4)
        call vdln(kk2, vectmp5, vectmp5)

        kk1 = 0
        kk2 = 0
        do k = 1, icount

            j = temp_jj(k)
#ifdef LES
            jstrt = ncopy*(j - 1)
            jcnum = cnum(j)
#endif

            r2 = r2x(k)

            rj = rborn(j) - offset
            rj1i = one/rj
            sj = fs(j)

            sj2 = sj*sj

            xij = xi - x(3*j - 2)
            yij = yi - x(3*j - 1)
            zij = zi - x(3*j)

            dij1i = vectmp1(k)
            dij = r2*dij1i

            temp3 = zero
            temp4 = zero

            if (dij <= rgbmax + sj) then

                if ((dij > rgbmax - sj)) then
                    uij = 1.0d0/(dij - sj)

! carlos: store descreening contrib in temp3, this makes it easier to
! apply to multiple atoms for LES

                    temp3 = temp3 - eighth*dij1i*(one + two*dij*uij + &
                        rgbmax2i*(r2 - four*rgbmax*dij - sj2) + &
                        two*log((dij - sj)*rgbmax1i))

                else if (dij > four*sj) then

                    dij2i = dij1i*dij1i
                    tmpsd = sj2*dij2i
                    dumbo = ta + tmpsd*(tb + tmpsd*(tc + tmpsd*(td + tmpsd*tdd)))

                    temp3 = temp3 - tmpsd*sj*dij2i*dumbo

                    !     ---following are from the Appendix of Schaefer and Froemmel,
                    !        J. Mol. Biol. 216:1045-1066, 1990, divided by (4*Pi):

                else if (dij > ri + sj) then

                    kk1 = kk1 + 1
                    temp3 = temp3 - half*(sj/(r2 - sj2) + half*dij1i*vectmp4(kk1))

                    !-----------------------------------------------------------------

                else if (dij > abs(ri - sj)) then

                    kk1 = kk1 + 1
                    theta = half*ri1i*dij1i*(r2 + ri*ri - sj2)
                    temp3 = temp3 - fourth*(ri1i*(two - theta) &
                        - vectmp2(kk1) + dij1i*vectmp4(kk1))

                    !-----------------------------------------------------------------

                else if (ri < sj) then
                    kk1 = kk1 + 1
                    temp3 = temp3 - half*sj/(r2 - sj2) + ri1i &
                        + fourth*dij1i*vectmp4(kk1)

                    !-----------------------------------------------------------------

                else
                    kk1 = kk1 + 1
                end if  ! ( dij > 4.d0*sj )

                if ((igb == 7 .or. igb == 8) .and. dij < rborn(i) + rborn(j) + GBNECKCUT) then
                    mdist = dij - neckMaxPos(neckidx(i), neckidx(j))
                    mdist2 = mdist*mdist
                    mdist3 = mdist2*mdist
                    mdist6 = mdist3*mdist3
                    neck = neckMaxVal(neckidx(i), neckidx(j))/ &
                        (one + mdist2 + 0.3d0*mdist6)
                    temp3 = temp3 - gbneckscale*neck
                end if

            end if

            !   --- Now the same thing, but swap i and j and use temp4 for descreening contrib:

            if (dij <= rgbmax + si) then

                if (dij > rgbmax - si) then
                    uij = 1.0d0/(dij - si)
                    temp4 = temp4 - eighth*dij1i*(one + two*dij*uij + &
                        rgbmax2i*(r2 - four*rgbmax*dij - si2) + &
                        two*log((dij - si)*rgbmax1i))
                else if (dij > four*si) then
                    dij2i = dij1i*dij1i
                    tmpsd = si2*dij2i
                    dumbo = ta + tmpsd*(tb + tmpsd*(tc + tmpsd*(td + tmpsd*tdd)))
                    temp4 = temp4 - tmpsd*si*dij2i*dumbo
                else if (dij > rj + si) then
                    kk2 = kk2 + 1
                    temp4 = temp4 - half*(si/(r2 - si2) + &
                        half*dij1i*vectmp5(kk2))
                    !-----------------------------------------------------------------
                else if (dij > abs(rj - si)) then
                    kk2 = kk2 + 1
                    theta = half*rj1i*dij1i*(r2 + rj*rj - si2)
                    temp4 = temp4 - fourth*(rj1i*(two - theta) &
                        - vectmp3(kk2) + dij1i*vectmp5(kk2))
                    !-----------------------------------------------------------------
                else if (rj < si) then
                    kk2 = kk2 + 1
                    temp4 = temp4 - half*si/(r2 - si2) + rj1i &
                        + fourth*dij1i*vectmp5(kk2)
                    !-----------------------------------------------------------------
                else
                    kk2 = kk2 + 1
                end if  ! ( dij > 4.d0*si )

                if ((igb == 7 .or. igb == 8) .and. dij < rborn(j) + rborn(i) + GBNECKCUT) then
                    mdist = dij - neckMaxPos(neckidx(j), neckidx(i))
                    mdist2 = mdist*mdist
                    mdist3 = mdist2*mdist
                    mdist6 = mdist3*mdist3
                    neck = neckMaxVal(neckidx(j), neckidx(i))/ &
                        (one + mdist2 + 0.3d0*mdist6)
                    temp4 = temp4 - gbneckscale*neck
                end if

            end if
! now add the calculated descreening component to onereff

#ifdef LES

            ! recall that reff and onereff do not use atom index as pointer
            ! but instead we need istrt and jstrt

            if (icnum == 0) then

                ! i is not a LES copy

                if (jcnum == 0) then

                    ! j not LES either, add temp3 to all copies of i's reff

                    do k1 = 1, ncopy

                        onereff(istrt + k1) = onereff(istrt + k1) + temp3
                        onereff(jstrt + k1) = onereff(jstrt + k1) + temp4
                    end do
                else

                    ! j is LES, add only to one of the radii for i and j
                    !  (the one for j's cnum)

                    onereff(istrt + jcnum) = onereff(istrt + jcnum) + temp3
                    onereff(jstrt + jcnum) = onereff(jstrt + jcnum) + temp4
                end if

            else

                ! i is LES, either j is not LES or is in same LES copy as i,
                ! so add to reff for i's cnum

                onereff(istrt + icnum) = onereff(istrt + icnum) + temp3
                onereff(jstrt + icnum) = onereff(jstrt + icnum) + temp4
            end if

#else /* not LES */

            reff_i = reff_i + temp3
            onereff(j) = onereff(j) + temp4
#endif

        end do                    !  k = 1, icount

        ! we are ending the do-i-loop, reassign the scalar to the original array:

#ifdef LES
        ! LES used onereff() directly so we don't need to reassign onereff(i)
        ! based on reff_i
#else
        onereff(i) = reff_i
#endif

    end do  !  i = 1,natom

#ifdef MPI
    call timer_stop_start(TIME_GBRAD1, TIME_GBRADDIST)

    !       collect the (inverse) effective radii from other nodes:

    if (numtasks > 1) then
#  ifdef LES
        ! LES has more reff
#   ifdef USE_MPI_IN_PLACE
        call mpi_allreduce(MPI_IN_PLACE, onereff, ncopy*natom, &
            MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#   else
        call mpi_allreduce(onereff, vectmp1, ncopy*natom, &
            MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
        onereff(1:ncopy*natom) = vectmp1(1:ncopy*natom)
#   endif
#  else
#   ifdef USE_MPI_IN_PLACE
        call mpi_allreduce(MPI_IN_PLACE, onereff, natom, &
            MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#   else
        call mpi_allreduce(onereff, vectmp1, natom, &
            MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
        onereff(1:natom) = vectmp1(1:natom)
#   endif
#  endif

    end if
    call timer_stop_start(TIME_GBRADDIST, TIME_GBRAD1)
#endif

! set end variable since LES arrays are longer than normal
! this will make the changes for LES code simpler to understand and modify

#ifdef LES
    vecend = natom*ncopy
#else
    vecend = natom
#endif

    if (igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) then

        ! --- apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:

#ifdef LES
! use rbornlong for these vectors, not rborn
        vectmp1(1:vecend) = rbornlong(1:vecend)
#else
        vectmp1(1:vecend) = rborn(1:vecend)
#endif
        vectmp2(1:vecend) = vectmp1(1:vecend) - offset
        call vdinv(vecend, vectmp1, vectmp1) !1.0d0/rborn

! use of onereff here means that we need changes for LES to accomodate multiple onereff per atom

        psi(1:vecend) = -vectmp2(1:vecend)*onereff(1:vecend)
        call vdinv(vecend, vectmp2, vectmp2) !1.0d0/(rborn-offset)

#ifdef LES
        !Hai Nguyen: Using gbalpha, gbbeta, gbgamma
        !Not use igb=8 for LES
        vectmp3(1:vecend) = ((gbalpha + gbgamma*psi(1:vecend)* &
            psi(1:vecend) - gbbeta*psi(1:vecend))*psi(1:vecend))
#else
        !if not LES: use gbvalpha, gbvgamma, gbvbeta arrays instead
        vectmp3(1:vecend) = ((gbvalpha(1:vecend) + gbvgamma(1:vecend)*psi(1:vecend)* &
            psi(1:vecend) - gbvbeta(1:vecend)*psi(1:vecend))*psi(1:vecend))
#endif
        call vdtanh(vecend, vectmp3, vectmp3)

        onereff(1:vecend) = vectmp2(1:vecend) - &
            (vectmp3(1:vecend)*vectmp1(1:vecend))
        do j = 1, vecend
            if (onereff(j) < zero) onereff(j) = thirtieth
        end do
        call vdinv(vecend, onereff, reff)
    else

        !       "standard" GB, including the "diagonal" term here:

#ifdef LES
! use rbornlong for these vectors, not rborn
        vectmp1(1:vecend) = rbornlong(1:vecend) - offset
#else
        vectmp1(1:vecend) = rborn(1:vecend) - offset
#endif

        call vdinv(vecend, vectmp1, vectmp1)
        onereff(1:vecend) = onereff(1:vecend) + vectmp1(1:vecend)
        call vdinv(vecend, onereff, reff)
    end if

#ifdef LES

    ! use RDT - determine if any of the non-LES atoms
    ! need multiple reff to be used

    do k = 1, natom

        ! loop over the multiple effective radii for the atom

        if (cnum(k) .eq. 0) then

            ! non-LES atoms have multiple radii, compare range to rdt

            if (rdt .gt. 0.) then

                istrt = ncopy*(k - 1) + 1
                iend = ncopy*(k)
                rmax = -10.
                rmin = 999.
                do i = istrt, iend
                    if (rmax <= reff(i)) rmax = reff(i)
                    if (rmin >= reff(i)) rmin = reff(i)
                end do
                if ((rmax - rmin) .gt. rdt) then
                    nradii(k) = ncopy
                else
                    nradii(k) = one
                end if
            else
                nradii(k) = ncopy
            end if
        else
            ! LES atom
            nradii(k) = one
        end if
    end do
#endif

    if (rbornstat == 1) then

        do k = 1, natom
#ifdef LES
            ! loop over the multiple effective radii for the atom

            if (cnum(i) .eq. 0) then

                ! non-LES atoms have multiple radii

                istrt = ncopy*(i - 1) + 1
                iend = ncopy*(i)
            else

                ! LES atoms only have one reff

                istrt = ncopy*(i - 1) + 1
                iend = ncopy*(i - 1) + 1
            end if
            do i = istrt, iend
#else
            i = k
#endif
            rbave(i) = rbave(i) + reff(i)
            rbfluct(i) = rbfluct(i) + reff(i)*reff(i)
            if (rbmax(i) <= reff(i)) rbmax(i) = reff(i)
            if (rbmin(i) >= reff(i)) rbmin(i) = reff(i)
#ifdef LES
        end do
#endif
    end do
end if

return

end subroutine egb_calc_radii

end module genborn
