! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

!Note: The routines contained here were updated in August 2008 by RCW
!      to add support for multiple SCEE and SCNB numbers. This is done
!      by having an array of SCEE and SCNB scale factors based on dihedral
!      type.

! (TD based on my understanding of caldwell's ideas)

!   At this time, all atom names beginning with EP or LP are
!   considered extra points.

!   For vacuum calculations (i.e. not periodic boundary conditions)
!   set use_pme=0
!   (this turns off all but direct sum of ewald)
!   set eedmeth = 4
!   (this uses straight coulomb and not erfc as in ewald)
!   use a big cutoff to get all pairs

!    dac note: in sander7, setting ntb=igb=0 automatically sets up a large
!      box and turns the above two parameters on.

!   If you want the lone pairs or extra points treated just like regular
!   atoms, set noextra = 1
!   (this disables the pattern to find them, so all are normal atoms)
!   The DEFAULT is noextra = 0

!   If you want bond angle and dihedral forces as usual for lone pairs
!   (i.e. as if they are amber atoms)
!   set frameon=0

!   If frameon is set to 1, (DEFAULT) the bonds, angles and dihedral interactions
!   involving the lone pairs/extra points are removed
!   except for constraints added during parm. The lone pairs are kept
!   in ideal geometry relative to local atoms, and resulting torques are
!   transferred to these atoms.

!   If chngmask=1 (DEFAULT), new 1-1, 1-2, 1-3 and 1-4 interactions are
!   calculated. An extra point belonging to an atom has a 1-1 interaction with it,
!   and participates in any 1-2, 1-3 or 1-4 interaction that atom has.

!   For example, suppose (excusing the geometry)
!   C1,C2,C3,C4 form a dihedral and each has 1 extra point attached as below

!           C1------C2------C3---------C4
!           |        |       |         |
!           |        |       |         |
!          Ep1      Ep2     Ep3       Ep4

!   The 1-4 interactions include  C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4

!   To see a printout of all 1-1, 1-2, 1-3 and 1-4 interactions
!   set verbose=1
!   These interactions are masked out of nonbonds. Thus the amber mask list is
!   rebuilt from these 1-1, 1-2, 1-3 and 1-4 pairs. I don't mask pairs that
!   aren't in the union of these.

!   A separate list of 1-4 nonbonds is then compiled. This list does not agree
!   in general with the above 1-4, since a 1-4 could also be a 1-3 if its
!   in a ring. I use the rules in EPHI to see who is included:

!   Here is that code

!             DO 700 JN = 1,MAXLEN
!               I3 = IP(JN+IST)
!               K3T = KP(JN+IST)
!               L3T = LP(JN+IST)
!               IC0 = ICP(JN+IST)
!               IDUMI = ISIGN(1,K3T)
!               IDUML = ISIGN(1,L3T)
!               KDIV = (2+IDUMI+IDUML)/4
!               L3 = IABS(L3T)
!               FMULN = FLOAT(KDIV)
!   C
!               II = (I3+3)/3
!               JJ = (L3+3)/3
!               IA1 = IAC(II)
!               IA2 = IAC(JJ)
!               IBIG = MAX0(IA1,IA2)
!               ISML = MIN0(IA1,IA2)
!               IC = IBIG*(IBIG-1)/2+ISML
!   C
!   C             ----- CALCULATE THE 14-EEL ENERGY -----
!   C
!               R2 = FMULN/CT(JN)
!               R1 = SQRT(R2)
!       ...........

!   so I include a pair in the 1-4 list if kdiv is > 0
!   this is decided at startup. This decision logic is applied to the parent
!   atoms, and if they are included, so are extra points attached:

!   That is, in the above situation, if C1 and C4 pass the test I include
!   C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4. I don't test the dihedrals
!   involving the extra points since the decision is based solely on
!   parent atoms.

!   The list of 1-4 nonbonds is also spit out if verbose=1.

!   To scale 1-4 charge-dipole and dipole-dipole interactions the same as
!   1-4 charge-charge (i.e. divided by scee)
!    set scaldip=1   (DEFAULT)
!   If scaldip=0 the 1-4 charge-dipole and dipole-dipole interactions
!   are treated the same as other dipolar interactions (i.e. divided by 1)

!-----------------------------------------------------------------------------
! Using the bond list (the one not involving hydrogens), find the number
! of neighbors (heavy atom, hydrogens and extra points) attached to each
! atom. For example if atom i is heavy, numnghbr(1,i) is the number of heavy
! atoms attached to atom i, while numnghbr(2,i) is the number of
! hydrogens attached, and numnghbr(3,i) is the number of
! extra points attached to i. The identities of neighbors are
! packed back to back in nghbrlst. the attyp array holds the
! atom types, usded to distinguish extra points or lone pairs from
! regular atoms.
!-----------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine init_extra_pts1 here]
subroutine init_extra_pts1()

    !     --- just figure out how much space is needed for extra_pts,
    !         and set up correct pointers.  We need to do this first for
    !         the new dynamic memory allocation scheme
    implicit none
#  include "memory.h"
#  include "extra_pts.h"
    integer max14

#ifndef LES
    call allocate_frames(numextra, ifrtyp, iatcen, inumep, &
        iepfr, ifrst, imid, ithrd, leploc)
    max14 = 12*(nphih + nphia + ndper)
    call allocate_14nb(inb_14, max14)
#endif

    return
end subroutine init_extra_pts1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine init_extra_pts here]
subroutine init_extra_pts(ibh, jbh, icbh, &
    ib, jb, icb, ith, jth, kth, icth, &
    it, jt, kt, ict, iph, jph, kph, lph, icph, &
    ipa, jpa, kpa, lpa, icpa, &
    isymbl, ix, x, iblo, inb, &
    nspm, nsp, tma, tmass, tmassinv, amass, amassinv, req)
    implicit none
    integer ibh(*), jbh(*), icbh(*), ib(*), jb(*), icb(*), &
        ith(*), jth(*), kth(*), icth(*), &
        it(*), jt(*), kt(*), ict(*), &
        iph(*), jph(*), kph(*), lph(*), icph(*), &
        ipa(*), jpa(*), kpa(*), lpa(*), icpa(*), &
        ix(*), iblo(*), inb(*)
    character(len=4) isymbl(*)
    integer iz
    _REAL_ x(*)
    integer nspm, nsp(*)
    _REAL_ tma(*), tmass, tmassinv, amass(*), amassinv(*), req(*)

#  include "extra_pts.h"
#  include "ew_cntrl.h"
#  include "memory.h"

    integer, dimension(:), allocatable :: epbtyp
    integer, dimension(:), allocatable :: nghbrs
    integer, dimension(:), allocatable :: hnghbrs
    integer, dimension(:), allocatable :: enghbrs
    integer, dimension(:), allocatable :: numnghbr
    integer, dimension(:), allocatable :: epowner
    integer, dimension(:), allocatable :: offset
    integer, dimension(:), allocatable :: test
    integer, dimension(:), allocatable :: i11
    integer, dimension(:), allocatable :: i12
    integer, dimension(:), allocatable :: i13
    integer, dimension(:), allocatable :: i14
    integer, dimension(:), allocatable :: nb_14_list
    integer, dimension(:), allocatable :: s3
    integer, dimension(:), allocatable :: s4

    character(len=4) ep, blank
    integer n, numextra_test
    integer num11, num12, num13, num14, max11, max12, max13, max14
    integer maxa, ier

    ep = 'EP  '
    blank = '    '

    numextra_test = 0
    do n = 1, natom
        if (isymbl(n) == ep) then
            numextra_test = numextra_test + 1
        end if
    end do
    if (numextra_test /= numextra) then
        write (6, *) 'Error in numextra_test'
        call mexit(6, 1)
    end if
    if (numextra == 0) frameon = 0
#ifdef LES
    if (numextra /= 0) then
        write (6, *) 'LES requires numextra=0'
        call mexit(6, 1)
    end if
    return
#endif

    max11 = natom + numextra
    max12 = 3*(nbonh + nbona + nbper)
    max13 = 3*(ntheth + ntheta + ngper)
    max14 = 12*(nphih + nphia + ndper)
    maxa = max(max11, max12, max13, max14)
    allocate (s3(maxa), stat=ier)
    REQUIRE(ier == 0)
    allocate (s4(maxa), stat=ier)
    REQUIRE(ier == 0)

    allocate (epbtyp(5*natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (nghbrs(5*natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (hnghbrs(5*natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (enghbrs(5*natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (numnghbr(3*natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (epowner(natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (offset(natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (test(natom), stat=ier)
    REQUIRE(ier == 0)
    allocate (i11(2*max11), stat=ier)
    REQUIRE(ier == 0)
    allocate (i12(2*max12), stat=ier)
    REQUIRE(ier == 0)
    allocate (i13(2*max13), stat=ier)
    REQUIRE(ier == 0)
    allocate (i14(2*max14), stat=ier)
    REQUIRE(ier == 0)
    allocate (nb_14_list(max14), stat=ier)
    REQUIRE(ier == 0)

    call get_nghbrs(ibh, jbh, ib, jb, icb, nbonh, nbona + nbper, &
        natom, isymbl, ep, &
        nghbrs, hnghbrs, enghbrs, &
        numnghbr, epowner, epbtyp)
    call define_frames(natom, isymbl, &
        nghbrs, hnghbrs, enghbrs, &
        numnghbr, &
        ix(ifrtyp), ix(iatcen), ix(inumep), ix(iepfr), &
        ix(ifrst), ix(imid), ix(ithrd), x(leploc), numfr, &
        epbtyp, req, verbose)

    call fill_bonded(max11, max12, max13, max14, &
        num11, num12, num13, num14, &
        i11, i12, i13, i14, &
        enghbrs, &
        numnghbr, epowner, natom, &
        ibh, jbh, ib, jb, ith, kth, it, kt, iph, lph, ipa, lpa, &
        nbonh, nbona, nbper, ntheth, ntheta, ngper, nphih, nphia, ndper, &
        offset, test, s3, verbose)

    if (chngmask == 1) then
        call redo_masked(natom, iblo, inb, nnb, &
            num11, num12, num13, num14, &
            i11, i12, i13, i14, &
            offset, test)
    end if

    call build_14nb(nb_14_list, numnb14, max14, &
        iph, jph, kph, lph, icph, ipa, jpa, kpa, lpa, icpa, &
        nphih, nphia, ndper, &
        epowner, numnghbr, enghbrs, &
        natom, offset, test, s3, s4, &
        chngmask, verbose)
    call copy_14nb(nb_14_list, ix(inb_14), numnb14)

    if (frameon == 1) then

        !      ---zero out mass and massinv for extra points:

        call fix_masses(natom, epowner, &
            nspm, nsp, tma, tmass, tmassinv, amass, amassinv)

        !      ---now remove bonds etc involving extra points:

        iz = 0
        call trim_bonds(ibh, jbh, icbh, nbonh, iz, epowner)
        call trim_bonds(ib, jb, icb, nbona, nbper, epowner)
        call trim_theta(ith, jth, kth, icth, ntheth, iz, epowner)
        call trim_theta(it, jt, kt, ict, ntheta, ngper, epowner)
        call trim_phi(iph, jph, kph, lph, icph, nphih, iz, epowner)
        call trim_phi(ipa, jpa, kpa, lpa, icpa, nphia, ndper, epowner)
    end if

    deallocate (s3)
    deallocate (s4)
    deallocate (epbtyp)
    deallocate (nghbrs)
    deallocate (hnghbrs)
    deallocate (enghbrs)
    deallocate (numnghbr)
    deallocate (epowner)
    deallocate (offset)
    deallocate (test)
    deallocate (i11)
    deallocate (i12)
    deallocate (i13)
    deallocate (i14)
    deallocate (nb_14_list)

    return
end subroutine init_extra_pts
!---------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+      GET_NGHBRS
!---------------------------------------------------------------
!+   find atoms involved in each center with an EP
!
subroutine get_nghbrs(ibh, jbh, ib, jb, icb, numbh, numb, &
    natom, isymbl, ep, &
    nghbrs, hnghbrs, enghbrs, numnghbr, epowner, epbtyp)
    implicit none
    integer ibh(*), jbh(*), ib(*), jb(*), icb(*), numbh, numb, natom
    character(len=4) isymbl(*), ep
    integer nghbrs(5, natom), hnghbrs(5, natom), enghbrs(5, natom), &
        numnghbr(3, natom), epowner(natom), epbtyp(5, natom)
#  include "extra_pts.h"
    integer n, ii, jj

    epowner = 0
    numnghbr = 0
    nghbrs = 0
    hnghbrs = 0
    enghbrs = 0
    epbtyp = 0

    do n = 1, numb
        ii = (ib(n) + 3)/3
        jj = (jb(n) + 3)/3
        if (isymbl(ii) == ep) then
            numnghbr(3, jj) = numnghbr(3, jj) + 1
            enghbrs(numnghbr(3, jj), jj) = ii
            epowner(ii) = jj
            epbtyp(numnghbr(3, jj), jj) = icb(n)
        else if (isymbl(jj) == ep) then
            numnghbr(3, ii) = numnghbr(3, ii) + 1
            enghbrs(numnghbr(3, ii), ii) = jj
            epowner(jj) = ii
            epbtyp(numnghbr(3, ii), ii) = icb(n)
        else
            numnghbr(1, ii) = numnghbr(1, ii) + 1
            numnghbr(1, jj) = numnghbr(1, jj) + 1
            nghbrs(numnghbr(1, ii), ii) = jj
            nghbrs(numnghbr(1, jj), jj) = ii
        end if
    end do
    do n = 1, numbh
        ii = (ibh(n) + 3)/3
        jj = (jbh(n) + 3)/3
        numnghbr(2, ii) = numnghbr(2, ii) + 1
        numnghbr(2, jj) = numnghbr(2, jj) + 1
        hnghbrs(numnghbr(2, ii), ii) = jj
        hnghbrs(numnghbr(2, jj), jj) = ii
    end do
    return
end subroutine get_nghbrs
!---------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+      DEFINE_FRAMES
!-----------------------------------------------------------------------
!+   fix the positions of EP in the local frame/coord depending
!+   on the kind of frame and atom types.
!
subroutine define_frames(natom, isymbl, &
    nghbrs, hnghbrs, enghbrs, numnghbr, &
    frtype, atcenter, numep, epframe, first, middle, third, &
    eplocal, numfr, epbtyp, req, verbose)
    use constants, only : DEG_TO_RAD
    implicit none
    character(len=4) isymbl(*)
    integer natom
    integer nghbrs(5, *), hnghbrs(5, *), enghbrs(5, *), &
        numnghbr(3, *), epbtyp(5, *)
    integer frtype(*), atcenter(*), numep(*), epframe(2, *), &
        first(*), middle(*), third(*), numfr, verbose
    _REAL_ eplocal(3, 2, *), req(*)
    integer k, m, n, l

    _REAL_ tetcos, tetsin, angle, scos, ssin
    character(len=4) sulf, sulfh

    sulf = 'S   '
    sulfh = 'SH  '

    !     --- get half-angle for tetrahedral

    angle = 54.735d0
    angle = angle*DEG_TO_RAD
    tetcos = cos(angle)
    tetsin = sin(angle)

    !     --- get cos,sin for 60

    scos = 0.5d0
    ssin = sqrt(1.d0 - scos*scos)

    numfr = 0
    do n = 1, natom

        if (numnghbr(3, n) > 0) then

            if (numnghbr(1, n) + numnghbr(2, n) > 2) then
                write (6, *) 'EXTRA_PTS: too many nghbrs!!'
                call mexit(6, 1)
            end if
            numfr = numfr + 1
            atcenter(numfr) = n
            numep(numfr) = numnghbr(3, n)
            do k = 1, 2
                epframe(k, numfr) = 0
            end do
            do k = 1, numep(numfr)
                epframe(k, numfr) = enghbrs(k, n)
            end do

            if (numnghbr(1, n) == 0 .and. numnghbr(2, n) == 2 &
                .and. numnghbr(3, n) == 1) then

                !----- TIP4P water: ---------------------------------
                !   (temporarily assign frtype to 3,
                !                     will set to 1 in next section)

                frtype(numfr) = 3
                first(numfr) = hnghbrs(1, n)
                middle(numfr) = n
                third(numfr) = hnghbrs(2, n)

            else if (numnghbr(1, n) == 0 .and. numnghbr(2, n) == 2 &
                .and. numnghbr(3, n) == 2) then

                !---  TIP5P water: -------------------------------

                frtype(numfr) = 1
                first(numfr) = hnghbrs(1, n)
                middle(numfr) = n
                third(numfr) = hnghbrs(2, n)

            else if (numnghbr(1, n) > 1) then

                !--- "ordinary" type of frame ----------------------
                !      defined by two other heavy atoms:

                frtype(numfr) = 1
                first(numfr) = nghbrs(1, n)
                middle(numfr) = n
                third(numfr) = nghbrs(2, n)

            else if (numnghbr(1, n) == 1 .and. numnghbr(2, n) == 1) then

                !--- frame defined by one heavy atom and one hydrogen: -------

                frtype(numfr) = 1
                first(numfr) = nghbrs(1, n)
                middle(numfr) = n
                third(numfr) = hnghbrs(1, n)

            else if (numnghbr(1, n) == 1 .and. numnghbr(2, n) == 0) then

                !--- Assume this is CARBONYL oxygen. -----------------------
                !   (Need to use midpoints of other two bonds of the carbon
                !   for first and third in orient force, thus (mis)use
                !   first middle third and atcenter to store 4 atoms for
                !   this special case.)

                frtype(numfr) = 2
                m = nghbrs(1, n)
                middle(numfr) = m
                if (numnghbr(1, m) /= 3 .or. numnghbr(2, m) > 0) then
                    write (6, *) 'EXTRA_PTS: frtype 2 Should not be here'
                    write (6, *) n, m, numnghbr(1, m), numnghbr(2, m)
                    call mexit(6, 1)
                end if

                ! numnghbr(1,m) = 3. One is n (the oxygen) and other 2 are
                ! other bonding partners of carbon
                first(numfr) = 0
                third(numfr) = 0
                k = 1
                do while ((k < 4) .and. (first(numfr) == 0))
                    if (nghbrs(k, m) /= n) then
                        first(numfr) = nghbrs(k, m)
                    end if
                    k = k + 1
                end do
                k = 1
                do while ((k < 4) .and. (third(numfr) == 0))
                    if (nghbrs(k, m) /= n .and. &
                        nghbrs(k, m) /= first(numfr)) then
                        third(numfr) = nghbrs(k, m)
                    end if
                    k = k + 1
                end do
                if ((first(numfr) == 0) .or. (third(numfr) == 0)) then
                    write (6, *) 'EXTRA_PTS: cannot find first or third frame point '
                    write (6, *) 'define: ', n, numnghbr(1, n), &
                        numnghbr(2, n), numnghbr(3, n), first(numfr), third(numfr)
                    call mexit(6, 1)
                end if
            else
                write (6, *) 'EXTRA_PTS: unexpected numnghbr array: '
                write (6, *) 'define: ', n, numnghbr(1, n), &
                    numnghbr(2, n), numnghbr(3, n)
                call mexit(6, 1)
            end if
        end if
    end do

    !--- get the local coords ---------------------------------------

    eplocal(1:3, 1:2, 1:numfr) = 0.d0
    do n = 1, numfr
        l = atcenter(n)
        if (frtype(n) == 1 .or. frtype(n) == 3) then

            !  z axis along symmetry axis of second atom opposite
            !      bisector of first,third;
            !  x axis along the diff vector third minus first
            !  y axis is cross product

            if (numep(n) == 1) then

                ! extra point is along the z-direction: positive for ordinary
                ! lone pair, negative for TIP4P water extra point:

                eplocal(3, 1, n) = req(epbtyp(1, l))
                if (frtype(n) == 3) then
                    eplocal(3, 1, n) = -req(epbtyp(1, l))
                    frtype(n) = 1
                end if

            else if (numep(n) == 2) then

                ! extra points are are in the z,y plane, tetrahedrally
                ! (unless middle atom is sulfur, in which case they
                ! are opposite along y):

                m = middle(n)
                if (isymbl(m) == sulf .or. isymbl(m) == sulfh) then
                    eplocal(2, 1, n) = req(epbtyp(1, l))
                    eplocal(2, 2, n) = -req(epbtyp(2, l))
                else
                    eplocal(3, 1, n) = tetcos*req(epbtyp(1, l))
                    eplocal(2, 1, n) = tetsin*req(epbtyp(1, l))
                    eplocal(3, 2, n) = tetcos*req(epbtyp(2, l))
                    eplocal(2, 2, n) = -tetsin*req(epbtyp(2, l))
                end if

            else
                write (6, *) 'EXTRA_PTS: unexpected numep value: ', numep(n)
                call mexit(6, 1)
            end if  ! ( numep(n) == 1 )

        else if (frtype(n) == 2) then

            ! z axis is along bond from middle to atcenter
            ! x axis in plane of atcenter and midpoints of first,middle
            !    and middle,third

            if (numep(n) == 1) then
                eplocal(3, 1, n) = req(epbtyp(1, l))
            else if (numep(n) == 2) then
                eplocal(3, 1, n) = scos*req(epbtyp(1, l))
                eplocal(1, 1, n) = ssin*req(epbtyp(1, l))
                eplocal(3, 2, n) = scos*req(epbtyp(2, l))
                eplocal(1, 2, n) = -ssin*req(epbtyp(2, l))
            else
                write (6, *) 'EXTRA_PTS: unexpected numep value: ', numep(n)
                call mexit(6, 1)
            end if
        else
            write (6, *) 'EXTRA_PTS: unexpected frtype value: ', frtype(n)
            call mexit(6, 1)
        end if
    end do

    if (verbose > 3) then
        write (6, *) 'frames:'
        do n = 1, numfr
            write (6, 666) n, atcenter(n), &
                isymbl(atcenter(n)), numep(n), &
                epframe(1, n), epframe(2, n), frtype(n), &
                first(n), middle(n), third(n)
            write (6, 667) eplocal(1, 1, n), eplocal(2, 1, n), eplocal(3, 1, n), &
                eplocal(1, 2, n), eplocal(2, 2, n), eplocal(3, 2, n)
        end do
    end if

666 format(1x, 2i7, 1x, a4, 7i7)
667 format(1x, 6(1x, f10.4))

    return
end subroutine define_frames

!---------------------------------------------------------------
!     LOCAL TO GLOBAL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calls do_local_to_global() to
!+ Put EP in position in world coord system based on the
!  position of the frame and the local coordinates.
!
subroutine local_to_global(crd, x, ix)
    implicit none

    _REAL_ crd(3, *), x(*)
    integer ix(*)
#  include "extra_pts.h"
#  include "ew_cntrl.h"
    !     integer iproc,numtasks

    if (frameon == 0 .or. numfr == 0) return
    call do_local_global(crd, ix(ifrtyp), ix(iatcen), ix(inumep), &
        ix(iepfr), ix(ifrst), ix(imid), ix(ithrd), x(leploc), numfr)
    return
end subroutine local_to_global
!---------------------------------------------------------------

!----------------------------------------------------------------------
!     DO_LOCAL_TO_GLOBAL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Put EP in position in world coord system based on the
!  position of the frame and the local coordinates.
!
subroutine do_local_global(crd, frtype, atcenter, numep, &
    epframe, first, middle, third, eplocal, numfr)
    !        iproc,numtasks)
    implicit none
    _REAL_ crd(3, *), eplocal(3, 2, *)
    integer frtype(*), atcenter(*), numep(*), epframe(2, *), numfr
    integer first(*), middle(*), third(*)
    !     integer iproc,numtasks

    !     Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985)
    !       plus weird modification for carbonyl to symmetrize

    _REAL_ uvec(3), vvec(3), ave(3), diff(3), &
        usiz, vsiz, asiz, dsiz, f(3, 3)
    _REAL_ a(3), b(3), c(3)
    integer j, k, m, n

    if (numfr == 0) return
    do n = 1, numfr
        if (frtype(n) == 1) then
            do m = 1, 3
                a(m) = crd(m, first(n))
                b(m) = crd(m, middle(n))
                c(m) = crd(m, third(n))
            end do
        else if (frtype(n) == 2) then
            do m = 1, 3
                a(m) = 0.5d0*(crd(m, first(n)) + crd(m, middle(n)))
                b(m) = crd(m, atcenter(n))
                c(m) = 0.5d0*(crd(m, third(n)) + crd(m, middle(n)))
            end do
        end if

        !       z-axis along symmmetry axis of b midway between
        !         unit vector to a and unit vector to c; points opposite

        usiz = 0.d0
        vsiz = 0.d0
        do m = 1, 3
            uvec(m) = a(m) - b(m)
            usiz = usiz + uvec(m)*uvec(m)
            vvec(m) = c(m) - b(m)
            vsiz = vsiz + vvec(m)*vvec(m)
        end do
        usiz = sqrt(usiz)
        vsiz = sqrt(vsiz)
        asiz = 0.d0
        dsiz = 0.d0
        do m = 1, 3
            uvec(m) = uvec(m)/usiz
            vvec(m) = vvec(m)/vsiz
            ave(m) = (uvec(m) + vvec(m))/2.d0
            asiz = asiz + ave(m)*ave(m)
            diff(m) = (vvec(m) - uvec(m))/2.d0
            dsiz = dsiz + diff(m)*diff(m)
        end do
        asiz = sqrt(asiz)
        dsiz = sqrt(dsiz)
        do m = 1, 3
            f(m, 3) = -ave(m)/asiz
            f(m, 1) = diff(m)/dsiz
        end do
        f(1, 2) = f(2, 3)*f(3, 1) - f(3, 3)*f(2, 1)
        f(2, 2) = f(3, 3)*f(1, 1) - f(1, 3)*f(3, 1)
        f(3, 2) = f(1, 3)*f(2, 1) - f(2, 3)*f(1, 1)
        do k = 1, numep(n)
            j = epframe(k, n)
            do m = 1, 3
                crd(m, j) = crd(m, atcenter(n)) + eplocal(1, k, n)*f(m, 1) + &
                    eplocal(2, k, n)*f(m, 2) + eplocal(3, k, n)*f(m, 3)
            end do
        end do
    end do  !  n = 1,numfr
    return
end subroutine do_local_global
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    ORIENT_FRC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calls do_orient_frc() to
!+ transfer forces from EP to main atoms in frame
!
subroutine orient_frc(crd, frc, framevir, ix)
    implicit none

    _REAL_ crd(3, *), frc(3, *)
    _REAL_ framevir(3, 3)
    integer ix(*)
#  include "extra_pts.h"
#  include "ew_cntrl.h"

    if (frameon == 0) return
    call do_orient_frc(crd, frc, framevir, &
        ix(ifrtyp), ix(iatcen), ix(inumep), &
        ix(iepfr), ix(ifrst), ix(imid), ix(ithrd), numfr)
    return
end subroutine orient_frc
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   DO_ ORIENT_FRC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ transfer forces from EP to main atoms in frame
!
subroutine do_orient_frc(crd, frc, framevir, frtype, atcenter, numep, &
    epframe, first, middle, third, numfr)
    implicit none
    _REAL_ crd(3, *), frc(3, *)
    _REAL_ framevir(3, 3)
    integer frtype(*), atcenter(*), numep(*), epframe(2, *), numfr
    integer first(*), middle(*), third(*)

    ! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985)

    _REAL_ u(3), v(3), w(3)
    _REAL_ up(3), vp(3), diff(3)
    _REAL_ usiz, vsiz, wsiz, upsiz, vpsiz
    _REAL_ dotdu, dotdv, dphidu, dphidv, dphidw
    _REAL_ c, s, uvdis, vudis, du(3), dv(3)
    _REAL_ force(3), torque(3), rel(3)
    integer i, j, k, l, j1, j2, j3, j4, m, n
    _REAL_ ap(3), bp(3), cp(3)

    ! motions of the frame can be described in terms of rotations about the
    ! unit vectors u and v from middle to first, third respecively
    ! and rotation about the unit cross product w

    if (numfr == 0) return

    framevir = 0.d0
    do n = 1, numfr
        force = 0.d0
        torque = 0.d0
        i = atcenter(n)
        do k = 1, numep(n)
            j = epframe(k, n)
            force(1:3) = force(1:3) + frc(1:3, j)
            rel(1:3) = crd(1:3, j) - crd(1:3, i)

#ifndef noVIRIAL
            !         ---get transferred force component of virial

            do m = 1, 3
                do l = 1, 3
                    framevir(l, m) = framevir(l, m) + frc(l, j)*rel(m)
                end do
            end do
#endif

            !---torque is rel x frc

            torque(1) = torque(1) + rel(2)*frc(3, j) - rel(3)*frc(2, j)
            torque(2) = torque(2) + rel(3)*frc(1, j) - rel(1)*frc(3, j)
            torque(3) = torque(3) + rel(1)*frc(2, j) - rel(2)*frc(1, j)
            frc(1, j) = 0.d0
            frc(2, j) = 0.d0
            frc(3, j) = 0.d0
        end do

        if (frtype(n) == 1) then
            do m = 1, 3
                ap(m) = crd(m, first(n))
                bp(m) = crd(m, middle(n))
                cp(m) = crd(m, third(n))
            end do
        else if (frtype(n) == 2) then
            do m = 1, 3
                ap(m) = 0.5d0*(crd(m, first(n)) + crd(m, middle(n)))
                bp(m) = crd(m, atcenter(n))
                cp(m) = 0.5d0*(crd(m, third(n)) + crd(m, middle(n)))
            end do
        end if
        usiz = 0.d0
        vsiz = 0.d0
        do m = 1, 3
            u(m) = ap(m) - bp(m)
            usiz = usiz + u(m)*u(m)
            v(m) = cp(m) - bp(m)
            vsiz = vsiz + v(m)*v(m)
        end do
        usiz = sqrt(usiz)
        vsiz = sqrt(vsiz)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
        wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
        dotdu = 0.d0
        dotdv = 0.d0
        do m = 1, 3
            u(m) = u(m)/usiz
            v(m) = v(m)/vsiz
            w(m) = w(m)/wsiz
            diff(m) = v(m) - u(m)
            dotdu = dotdu + u(m)*diff(m)
            dotdv = dotdv + v(m)*diff(m)
        end do

        !       ---get perps to u,v to get direction of motion of u or v
        !          due to rotation about the cross product vector w

        upsiz = 0.d0
        vpsiz = 0.d0
        do m = 1, 3
            up(m) = diff(m) - dotdu*u(m)
            vp(m) = diff(m) - dotdv*v(m)
            upsiz = upsiz + up(m)*up(m)
            vpsiz = vpsiz + vp(m)*vp(m)
        end do
        upsiz = sqrt(upsiz)
        vpsiz = sqrt(vpsiz)
        do m = 1, 3
            up(m) = up(m)/upsiz
            vp(m) = vp(m)/vpsiz
        end do

        !       ---negative of dot product of torque with unit vectors
        !          along u,v and w.  give result of infinitesmal rotation
        !          along these vectors, i.e. dphi/dtheta = dot product

        dphidu = -(torque(1)*u(1) + torque(2)*u(2) + torque(3)*u(3))
        dphidv = -(torque(1)*v(1) + torque(2)*v(2) + torque(3)*v(3))
        dphidw = -(torque(1)*w(1) + torque(2)*w(2) + torque(3)*w(3))

        !       ---get projected distances between vectors

        c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
        s = sqrt(1.d0 - c*c)
        uvdis = usiz*s
        vudis = vsiz*s

        !---------------------------------------------------------------------
        ! frame formed by bisector of u,v,  its perp, and w
        ! movement of u by dz out of plane -> rotation about v of -dz/uvdis
        ! since positive rotation about v move u in negative dir. wrt w
        ! dphi/dz = dphi/dtheta dtheta/dz = -dotvt /uvdis
        ! movement of v by dz out of plane -> rotation about u of dz/vudis
        ! movement of u by dy along up -> rotation about w of 1/2 dy/usiz
        ! since bisector only rotates 1/2 as much as u or v in isolation
        ! movement of v by dy along vperp -> rotation about w of 1/2 dy/vsiz
        ! movement of u by dx along u doesn't change frame
        ! movement of v by dx along v doesn't change frame
        ! So... du_du = 0, du_dw = -dotvt/uvdis, du_dup = dotwt/(2.d0*usiz)
        ! So... dv_dv = 0, dv_dw = dotut/vudis, du_dup = dotwt/(2.d0*usiz)
        !---------------------------------------------------------------------

        if (frtype(n) == 1) then
            j1 = first(n)
            j2 = middle(n)
            j3 = third(n)

            do m = 1, 3
                du(m) = -w(m)*dphidv/uvdis + up(m)*dphidw/(2.d0*usiz)
                dv(m) = w(m)*dphidu/vudis + vp(m)*dphidw/(2.d0*vsiz)
                frc(m, j1) = frc(m, j1) - du(m)
                frc(m, j3) = frc(m, j3) - dv(m)
                frc(m, j2) = frc(m, j2) + dv(m) + du(m) + force(m)
            end do

#ifndef noVIRIAL
            !         ---get torque contribution to virial

            do m = 1, 3
                do l = 1, 3
                    framevir(l, m) = framevir(l, m) + du(l)*(ap(m) - bp(m)) &
                        + dv(l)*(cp(m) - bp(m))
                end do
            end do
#endif

        else if (frtype(n) == 2) then

            !       ---need to transfer forces from midpoints to atoms

            j1 = first(n)
            j2 = middle(n)
            j3 = third(n)
            j4 = atcenter(n)
            do m = 1, 3
                du(m) = -w(m)*dphidv/uvdis + up(m)*dphidw/(2.d0*usiz)
                dv(m) = w(m)*dphidu/vudis + vp(m)*dphidw/(2.d0*vsiz)
                frc(m, j1) = frc(m, j1) - 0.5d0*du(m)
                frc(m, j3) = frc(m, j3) - 0.5d0*dv(m)
                frc(m, j2) = frc(m, j2) - 0.5d0*(du(m) + dv(m))
                frc(m, j4) = frc(m, j4) + dv(m) + du(m) + force(m)
            end do

#ifndef noVIRIAL
            !         ---get torque contribution to virial

            do m = 1, 3
                do l = 1, 3
                    framevir(l, m) = framevir(l, m) + du(l)*(ap(m) - bp(m)) &
                        + dv(l)*(cp(m) - bp(m))
                end do
            end do
#endif

        end if  ! ( frtype(n) == 1 )

        !---------------------------------------------------------------------
        ! OTHER TYPE FRAME; NOT SEEN YET
        ! frame formed by  u, its perp, and w
        ! movement of v in plane doesn't change frame
        ! movement of u by dz out of plane -> rotation about v of -dz/uvdis
        ! since positive rotation about v move u in negative dir. wrt w
        ! dphi/dz = dphi/dtheta dtheta/dz = -dotvt /uvdis
        ! movement of v by dz out of plane -> rotation about u of dz/vudis
        ! movement of u by dy along up -> rotation about w of dy/usiz
        ! since frame rotates as much as u in isolation
        !---------------------------------------------------------------------

199     continue
    end do
    return
end subroutine do_orient_frc
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    ZERO_EXTRA_PNTS_VEC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calls do_zero_extra_pnts_vec() to
!+ zero the force or velocity vectors associated with extra points.
!
subroutine zero_extra_pnts_vec(vec, ix)
    implicit none

    _REAL_ vec(3, *)
    integer ix(*)
#  include "extra_pts.h"
#  include "ew_cntrl.h"

    if (frameon == 0) return
    call do_zero_extra_pnts_vec(vec, ix(inumep), ix(iepfr), numfr)
    return
end subroutine zero_extra_pnts_vec
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   DO_ZERO_EXTRA_PNTS_VEC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Zero the force or velocity vectors associated with extra points.
!
subroutine do_zero_extra_pnts_vec(vec, numep, epframe, numfr)
    implicit none
    _REAL_ vec(3, *)
    integer numep(*), epframe(2, *), numfr

    integer j, k, n

    if (numfr == 0) return

    do n = 1, numfr
        do k = 1, numep(n)
            j = epframe(k, n)
            vec(:, j) = 0.d0
        end do
    end do  !  n = 1,numfr

    return

end subroutine do_zero_extra_pnts_vec
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!       FILL_BONDED
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine fill_bonded(max11, max12, max13, max14, &
    num11, num12, num13, num14, &
    list11, list12, list13, list14, &
    enghbrs, numnghbr, epowner, natom, &
    ibh, jbh, ib, jb, ith, kth, it, kt, iph, lph, ip, lp, &
    nbonh, nbona, nbper, ntheth, ntheta, ngper, nphih, nphia, ndper, &
    scr1, scr2, scr3, verbose)
    implicit none
    integer max11, max12, max13, max14, num11, num12, num13, num14
    integer list11(2, max11), list12(2, max12), &
        list13(2, max13), list14(2, max14)
    integer scr1(*), scr2(*), scr3(*)
    integer enghbrs(5, *), &
        numnghbr(3, *), epowner(*), natom
    integer ibh(*), jbh(*), ib(*), jb(*), &
        ith(*), kth(*), it(*), kt(*), iph(*), lph(*), ip(*), lp(*), &
        nbonh, nbona, ntheth, ntheta, nphih, nphia, verbose, &
        nbper, ngper, ndper
    integer n, ifail
    num11 = 0
    num12 = 0
    num13 = 0
    num14 = 0
    do n = 1, natom
        if (numnghbr(3, n) == 1) then
            num11 = num11 + 1
            if (num11 > max11) goto 100
            list11(1, num11) = n
            list11(2, num11) = enghbrs(1, n)
        else if (numnghbr(3, n) == 2) then
            if (num11 + 3 > max11) goto 100
            num11 = num11 + 1
            if (num11 > max11) goto 100
            list11(1, num11) = n
            list11(2, num11) = enghbrs(1, n)
            num11 = num11 + 1
            if (num11 > max11) goto 100
            list11(1, num11) = n
            list11(2, num11) = enghbrs(2, n)
            num11 = num11 + 1
            if (num11 > max11) goto 100
            list11(1, num11) = enghbrs(1, n)
            list11(2, num11) = enghbrs(2, n)
        else if (numnghbr(3, n) == 3) then
            goto 500
        end if
    end do
    !   --- bonds ---------
    call do_pairs(ibh, jbh, nbonh, list12, num12, max12, &
        epowner, numnghbr, enghbrs, ifail)
    if (ifail == 1) goto 200
    call do_pairs(ib, jb, nbona + nbper, list12, num12, max12, &
        epowner, numnghbr, enghbrs, ifail)
    call sort_pairs(list12, num12, natom, scr1, scr2, scr3)
    if (ifail == 1) goto 200
    !   --- angles --------
    call do_pairs(ith, kth, ntheth, list13, num13, max13, &
        epowner, numnghbr, enghbrs, ifail)
    if (ifail == 1) goto 300
    call do_pairs(it, kt, ntheta + ngper, list13, num13, max13, &
        epowner, numnghbr, enghbrs, ifail)
    if (ifail == 1) goto 300
    call sort_pairs(list13, num13, natom, scr1, scr2, scr3)
    !   --- dihedrals -----
    call do_pairs(iph, lph, nphih, list14, num14, max14, &
        epowner, numnghbr, enghbrs, ifail)
    if (ifail == 1) goto 400
    call do_pairs(ip, lp, nphia + ndper, list14, num14, max14, &
        epowner, numnghbr, enghbrs, ifail)
    if (ifail == 1) goto 400
    call sort_pairs(list14, num14, natom, scr1, scr2, scr3)

    if (verbose > 0) &
        write (6, '(a,4i6)') '| EXTRA PTS fill_bonded: num11-14 = ', &
        num11, num12, num13, num14
    if (verbose > 3) then
        write (6, *) '$$$$$$$$$$$$$$$$$ 1-1 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
        do n = 1, num11
            write (6, 666) n, list11(1, n), list11(2, n)
        end do
        write (6, *) '$$$$$$$$$$$$$$$$$ 1-2 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
        do n = 1, num12
            write (6, 666) n, list12(1, n), list12(2, n)
        end do
        write (6, *) '$$$$$$$$$$$$$$$$$ 1-3 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
        do n = 1, num13
            write (6, 666) n, list13(1, n), list13(2, n)
        end do
        write (6, *) '$$$$$$$$$$$$$$$$$ 1-4 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
        do n = 1, num14
            write (6, 666) n, list14(1, n), list14(2, n)
        end do
        write (6, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
666     format(1x, i5, ':', 3x, 'i,j = ', i5, 1x, i5)
    end if
    return
100 write (6, *) 'fill_bonded: max11 exceeded!!'
    call mexit(6, 1)
200 write (6, *) 'fill_bonded: max12 exceeded!!'
    call mexit(6, 1)
300 write (6, *) 'fill_bonded: max13 exceeded!!'
    call mexit(6, 1)
400 write (6, *) 'fill_bonded: max14 exceeded!!'
    call mexit(6, 1)
500 write (6, *) 'fill_bonded: should not be here!'
end subroutine fill_bonded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sort_pairs here]
subroutine sort_pairs(list, num, natom, scr1, scr2, scr3)
    implicit none
    integer list(2, *), num, natom, scr1(*), scr2(*), scr3(*)
    integer i, j, k, m, n, ntot
    do n = 1, num
        i = list(1, n)
        j = list(2, n)
        if (i < j) then
            list(1, n) = i
            list(2, n) = j
        else
            list(1, n) = j
            list(2, n) = i
        end if
    end do

    !     ---now get rid of duplicates

    !     ---first pass

    do n = 1, natom
        scr1(n) = 0
    end do
    do n = 1, num
        i = list(1, n)
        scr1(i) = scr1(i) + 1
    end do
    scr2(1) = 0
    do n = 2, natom
        scr2(n) = scr2(n - 1) + scr1(n - 1)
        scr1(n - 1) = 0
    end do
    scr1(natom) = 0

    !     ---second pass

    do n = 1, num
        i = list(1, n)
        j = list(2, n)
        scr1(i) = scr1(i) + 1
        scr3(scr1(i) + scr2(i)) = j
    end do
    do n = 1, natom
        scr2(n) = 0
    end do

    !     ---now trim them

    ntot = 0
    k = 0
    do n = 1, natom
        do m = 1, scr1(n)
            j = scr3(ntot + m)
            if (scr2(j) /= n) then
                k = k + 1
                list(1, k) = n
                list(2, k) = j
                scr2(j) = n
            end if
        end do
        ntot = ntot + scr1(n)
    end do
    num = k
    return
end subroutine sort_pairs
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sort_pairs here]
subroutine sort_pairs_14_nb(nb_14_list, num, natom, scr1, scr2, scr3, scr4)
    implicit none
    integer nb_14_list(3, *), num, natom, scr1(*), scr2(*), scr3(*), scr4(*)
    integer i, j, k, m, n, ntot, ic0

    !Indexes of nb_14_list = 1,j,ic0  (ic0=icp(n))
    !First pass, make sure the first index contains the lowest number. Ignore
    !index 3 here since this won't change anything.

    do n = 1, num
        i = nb_14_list(1, n)
        j = nb_14_list(2, n)
        if (i < j) then
            nb_14_list(1, n) = i
            nb_14_list(2, n) = j
        else
            nb_14_list(1, n) = j
            nb_14_list(2, n) = i
        end if
    end do

    !     ---now get rid of duplicates

    !     ---first pass

    scr1(1:natom) = 0
    do n = 1, num
        i = nb_14_list(1, n)
        scr1(i) = scr1(i) + 1
    end do
    scr2(1) = 0
    do n = 2, natom
        scr2(n) = scr2(n - 1) + scr1(n - 1)
        scr1(n - 1) = 0
    end do
    scr1(natom) = 0

    !     ---second pass

    do n = 1, num
        i = nb_14_list(1, n)
        j = nb_14_list(2, n)
        ic0 = nb_14_list(3, n)
        scr1(i) = scr1(i) + 1
        scr3(scr1(i) + scr2(i)) = j
        scr4(scr1(i) + scr2(i)) = ic0
    end do
    scr2(1:natom) = 0

    !     ---now trim them

    ntot = 0
    k = 0
    do n = 1, natom
        do m = 1, scr1(n)
            j = scr3(ntot + m)
            ic0 = scr4(ntot + m)
            if (scr2(j) /= n) then
                k = k + 1
                nb_14_list(1, k) = n
                nb_14_list(2, k) = j
                nb_14_list(3, k) = ic0
                scr2(j) = n
            end if
        end do
        ntot = ntot + scr1(n)
    end do
    num = k
    return
end subroutine sort_pairs_14_nb
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine redo_masked here]
subroutine redo_masked(natom, iblo, inb, nnb, &
    num11, num12, num13, num14, &
    list11, list12, list13, list14, offset, test)
    implicit none
    integer natom, iblo(*), inb(*), nnb, &
        num11, num12, num13, num14, &
        list11(2, *), list12(2, *), list13(2, *), list14(2, *), &
        offset(*), test(*)

    !     ---build the mask list from list11-14. Make sure no duplication

    integer j, n, m, ntot
    do n = 1, natom
        iblo(n) = 0
        offset(n) = 0
        test(n) = 0
    end do

    !     ---PASS 1 fill iblo

    call add_one_list_iblo(iblo, list11, num11)
    call add_one_list_iblo(iblo, list12, num12)
    call add_one_list_iblo(iblo, list13, num13)
    call add_one_list_iblo(iblo, list14, num14)

    !     ---check totals while finding offsets, resetting iblo

    ntot = 0
    do n = 1, natom
        offset(n) = ntot
        ntot = ntot + iblo(n)
        iblo(n) = 0
    end do
    if (ntot > 2*nnb) then
        write (6, *) 'EXTRA POINTS: nnb too small! '
        write (6, *) 'nnb,ntot = ', nnb, ntot
        call mexit(6, 1)
    end if

    !     ---PASS 2 fill inb, redo iblo

    call add_one_list_inb(iblo, inb, offset, list11, num11)
    call add_one_list_inb(iblo, inb, offset, list12, num12)
    call add_one_list_inb(iblo, inb, offset, list13, num13)
    call add_one_list_inb(iblo, inb, offset, list14, num14)

    !     ---PASS 3 filter inb, remove duplicate entries

    do n = 1, natom - 1
        do m = 1, iblo(n)
            j = inb(offset(n) + m)
            if (test(j) /= n) then
                test(j) = n
            else
                inb(offset(n) + m) = 0
            end if
        end do
    end do
    return
end subroutine redo_masked
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine add_one_list_iblo here]
subroutine add_one_list_iblo(iblo, list, num)
    implicit none
    integer iblo(*), list(2, *), num
    integer n, i
    do n = 1, num
        i = list(1, n)
        iblo(i) = iblo(i) + 1
    end do
    return
end subroutine add_one_list_iblo
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine add_one_list_inb here]
subroutine add_one_list_inb(iblo, inb, offset, list, num)
    implicit none
    integer iblo(*), inb(*), offset(*), list(2, *), num
    integer n, i, j, m
    do n = 1, num
        i = list(1, n)
        j = list(2, n)
        m = offset(i)
        iblo(i) = iblo(i) + 1
        inb(m + iblo(i)) = j
    end do
    return
end subroutine add_one_list_inb
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine build_14nb here]
subroutine build_14nb(nb_14_list, numnb14, maxnb14, &
    iph, jph, kph, lph, icph, ipa, jpa, kpa, lpa, icpa, &
    nphih, nphia, ndper, &
    epowner, numnghbr, enghbrs, &
    natom, scr1, scr2, scr3, scr4, &
    chngmask, verbose)
    implicit none

    integer nb_14_list(3, *), numnb14, maxnb14, &
        iph(*), jph(*), kph(*), lph(*), icph(*), &
        ipa(*), jpa(*), kpa(*), lpa(*), icpa(*), &
        nphih, nphia, ndper, chngmask, verbose
    integer enghbrs(5, *), numnghbr(3, *), epowner(*)
    integer natom, scr1(*), scr2(*), scr3(*), scr4(*)
    integer ifail, n
    numnb14 = 0
    call do_14pairs(iph, jph, kph, lph, icph, &
        nphih, nb_14_list, numnb14, maxnb14, &
        epowner, numnghbr, enghbrs, ifail, chngmask)
    if (ifail == 1) then
        write (6, *) 'exceeded maxnb14 in build14: check extra_pts.h'
        call mexit(6, 1)
    end if
    call do_14pairs(ipa, jpa, kpa, lpa, icpa, &
        nphia + ndper, nb_14_list, numnb14, maxnb14, &
        epowner, numnghbr, enghbrs, ifail, chngmask)
    if (ifail == 1) then
        write (6, *) 'exceeded maxnb14 in build14: check extra_pts.h'
        call mexit(6, 1)
    end if
    call sort_pairs_14_nb(nb_14_list, numnb14, natom, scr1, scr2, scr3, scr4)
    if (verbose > 0) write (6, '(a,i6)') &
        '| EXTRA_PTS, build_14: num of 14 terms = ', numnb14
    if (verbose > 3) then
        write (6, *) '$$$$$$$$$$$$$$$$$$$$$$$  1-4 nb list $$$$$$$$$$'
        do n = 1, numnb14
            write (6, 666) n, nb_14_list(1, n), nb_14_list(2, n)
        end do
666     format(1x, i5, ':', 3x, 'i,j = ', i5, 1x, i5)
        write (6, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
    end if
    return
end subroutine build_14nb

!---------------------------------------------------------------------
!       GET_14_CG
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine get_14_cg(charge, crd, frc, &
    iac, ico, ntypes, cn1, cn2, &
    ee14, enb14, one_scee, one_scnb, e14vir, ix, mytaskid, numtasks, eedmeth)
    implicit none
    _REAL_ charge(*), crd(3, *), frc(3, *), &
        cn1(*), cn2(*), e14vir(3, 3)
    integer iac(*), ico(*), ntypes, ix(*)
    _REAL_ ee14, enb14, one_scee(*), one_scnb(*)
    integer mytaskid, numtasks, eedmeth
#  include "extra_pts.h"
    call do_14_cg(charge, crd, frc, &
        iac, ico, ntypes, cn1, cn2, &
        ee14, enb14, one_scee, one_scnb, e14vir, ix(inb_14), numnb14, &
        mytaskid, numtasks, eedmeth)
    return
end subroutine get_14_cg

!--------------------------------------------------------------
!       DO_14_CG
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine do_14_cg(charge, crd, frc, &
    iac, ico, ntypes, cn1, cn2, &
    ee14, enb14, one_scee, one_scnb, e14vir, nb_14_list, numnb14, &
    mytaskid, numtasks, eedmeth)
    use constants, only : zero, one
    use decomp, only : decpr, decpair
    use crg_reloc, only : ifcr, cropt, cr_charge, cr_add_dcdr_factor
    use file_io_dat
#ifdef MPI /* SOFT CORE */
    use softcore, only : ifsc, nsc, sc_ener, oneweight
#endif
    implicit none
    _REAL_ charge(*), crd(3, *), frc(3, *), &
        cn1(*), cn2(*), e14vir(3, 3)
    integer iac(*), ico(*), ntypes, nb_14_list(3, *), numnb14
    _REAL_ ee14, enb14, one_scee(*), one_scnb(*)
    integer mytaskid, numtasks, eedmeth
#  include "flocntrl.h"
#  include "md.h"

    integer n, i, j, ic, ic0, ia1, ia2, ibig, isml
    _REAL_ dx, dy, dz, r2, r2inv, rinv, scnb0, scee0
    _REAL_ g, f6, f12, r6, df
    _REAL_ cgi, cgj
    ee14 = zero
    enb14 = zero
#ifdef MPI /* SOFT CORE */
    sc_ener(4) = 0.0d0
    sc_ener(5) = 0.0d0
#endif

    if (do_14 == 0) return

    e14vir = zero

    if (eedmeth == 5) then
        do n = mytaskid + 1, numnb14, numtasks
            i = nb_14_list(1, n)
            j = nb_14_list(2, n)
            ic0 = nb_14_list(3, n)
            scee0 = one_scee(ic0)
            scnb0 = one_scnb(ic0)
            dx = crd(1, j) - crd(1, i)
            dy = crd(2, j) - crd(2, i)
            dz = crd(3, j) - crd(3, i)
            r2 = dx*dx + dy*dy + dz*dz
            rinv = sqrt(1.d0/r2)
            r2inv = rinv*rinv
            r6 = r2inv*r2inv*r2inv
            if (ifcr /= 0 .and. cropt == 0) then
                cgi = cr_charge(i)
                cgj = cr_charge(j)
            else
                cgi = charge(i)
                cgj = charge(j)
            end if
            g = cgi*cgj*r2inv
            ee14 = ee14 + scee0*g
            !  always use the 6-12 parameters, even if 10-12 are available:
            ia1 = iac(i)
            ia2 = iac(j)
            ibig = max0(ia1, ia2)
            isml = min0(ia1, ia2)
            ic = ibig*(ibig - 1)/2 + isml
            f6 = cn2(ic)*r6
            f12 = cn1(ic)*(r6*r6)
            enb14 = enb14 + scnb0*(f12 - f6)
            df = (2.d0*scee0*g + scnb0*(12.d0*f12 - 6.d0*f6))*r2inv
            if (ifcr /= 0 .and. cropt /= 0) then
                call cr_add_dcdr_factor(i, cgj*r2inv*scee0)
                call cr_add_dcdr_factor(j, cgi*r2inv*scee0)
            end if
            ! -- ti decomp
            if (decpr .and. idecomp /= 0) then
                if (idecomp == 1) then
                    call decpair(4, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(4, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                elseif (idecomp == 2) then
                    call decpair(2, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(3, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                end if
            end if
            frc(1, j) = frc(1, j) + df*dx
            frc(2, j) = frc(2, j) + df*dy
            frc(3, j) = frc(3, j) + df*dz
            frc(1, i) = frc(1, i) - df*dx
            frc(2, i) = frc(2, i) - df*dy
            frc(3, i) = frc(3, i) - df*dz
#ifndef noVIRIAL
            e14vir(1, 1) = e14vir(1, 1) - df*dx*dx
            e14vir(1, 2) = e14vir(1, 2) - df*dx*dy
            e14vir(1, 3) = e14vir(1, 3) - df*dx*dz
            e14vir(2, 1) = e14vir(2, 1) - df*dy*dx
            e14vir(2, 2) = e14vir(2, 2) - df*dy*dy
            e14vir(2, 3) = e14vir(2, 3) - df*dy*dz
            e14vir(3, 1) = e14vir(3, 1) - df*dz*dx
            e14vir(3, 2) = e14vir(3, 2) - df*dz*dy
            e14vir(3, 3) = e14vir(3, 3) - df*dz*dz
#endif
        end do  !  n = 1,numnb14
    else
        do n = mytaskid + 1, numnb14, numtasks
            i = nb_14_list(1, n)
            j = nb_14_list(2, n)
            ic0 = nb_14_list(3, n)
            scee0 = one_scee(ic0)
            scnb0 = one_scnb(ic0)
            dx = crd(1, j) - crd(1, i)
            dy = crd(2, j) - crd(2, i)
            dz = crd(3, j) - crd(3, i)
            r2 = dx*dx + dy*dy + dz*dz
            rinv = sqrt(1.d0/r2)
            r2inv = rinv*rinv
            r6 = r2inv*r2inv*r2inv
            if (ifcr /= 0 .and. cropt == 0) then
                cgi = cr_charge(i)
                cgj = cr_charge(j)
            else
                cgi = charge(i)
                cgj = charge(j)
            end if
            g = cgi*cgj*rinv
            ee14 = ee14 + scee0*g
            !  always use the 6-12 parameters, even if 10-12 are available:
            ia1 = iac(i)
            ia2 = iac(j)
            ibig = max0(ia1, ia2)
            isml = min0(ia1, ia2)
            ic = ibig*(ibig - 1)/2 + isml
            f6 = cn2(ic)*r6
            f12 = cn1(ic)*(r6*r6)
            enb14 = enb14 + scnb0*(f12 - f6)
            df = (scee0*g + scnb0*(12.d0*f12 - 6.d0*f6))*r2inv
            if (ifcr /= 0 .and. cropt /= 0) then
                call cr_add_dcdr_factor(i, cgj*rinv*scee0)
                call cr_add_dcdr_factor(j, cgi*rinv*scee0)
            end if
#ifdef MPI /* SOFT CORE */
            ! For dual-topology softcore runs, 1-4 interactions involving sc atoms are modified here
            if (ifsc /= 0) then
                ! Check if a softcore atom is involved in this interaction
                if (nsc(i) == 1 .or. nsc(j) == 1) then
                    ! This interactions need to
                    ! a) get their energies removed from the 1-4 energies
                    ! b) get their force scaled up by 1/weight
                    ! This means reversing the work done above, but keeps the code simple
                    sc_ener(4) = sc_ener(4) + (scnb0*(f12 - f6))
                    enb14 = enb14 - (scnb0*(f12 - f6))
                    sc_ener(5) = sc_ener(5) + (scee0*g)
                    ee14 = ee14 - (scee0*g)
                    df = df*oneweight
                    if (ifcr /= 0 .and. cropt /= 0) then
                        call cr_add_dcdr_factor(i, -cgj*rinv*scee0)
                        call cr_add_dcdr_factor(j, -cgi*rinv*scee0)
                    end if
                end if
            end if

            ! -- ti decomp
            if (ifsc == 0) then
                if (decpr .and. idecomp == 1) then
                    call decpair(4, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(4, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                else if (decpr .and. idecomp == 2) then
                    call decpair(2, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(3, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                end if
            else if (nsc(i) /= 1 .and. nsc(j) /= 1) then
                if (decpr .and. idecomp == 1) then
                    call decpair(4, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(4, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                else if (decpr .and. idecomp == 2) then
                    call decpair(2, i, j, scee0*g/(nstlim/ntpr))
                    call decpair(3, i, j, scnb0*(f12 - f6)/(nstlim/ntpr))
                end if
            end if
#endif
            frc(1, j) = frc(1, j) + df*dx
            frc(2, j) = frc(2, j) + df*dy
            frc(3, j) = frc(3, j) + df*dz
            frc(1, i) = frc(1, i) - df*dx
            frc(2, i) = frc(2, i) - df*dy
            frc(3, i) = frc(3, i) - df*dz
#ifndef noVIRIAL
            e14vir(1, 1) = e14vir(1, 1) - df*dx*dx
            e14vir(1, 2) = e14vir(1, 2) - df*dx*dy
            e14vir(1, 3) = e14vir(1, 3) - df*dx*dz
            e14vir(2, 1) = e14vir(2, 1) - df*dy*dx
            e14vir(2, 2) = e14vir(2, 2) - df*dy*dy
            e14vir(2, 3) = e14vir(2, 3) - df*dy*dz
            e14vir(3, 1) = e14vir(3, 1) - df*dz*dx
            e14vir(3, 2) = e14vir(3, 2) - df*dz*dy
            e14vir(3, 3) = e14vir(3, 3) - df*dz*dz
#endif
        end do
    end if
    return
end subroutine do_14_cg
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_14_dipole here]
subroutine get_14_dipole(charge, crd, frc, &
    dipole, field, iac, ico, ntypes, cn1, cn2, &
 ! M-WJ, WJM, YD
    pol, dampfactor, pol2, &
 !
    ee14, enb14, epol14, one_scee, one_scnb, e14vir, ix, &
    mytaskid, numtasks)
    implicit none
    _REAL_ charge(*), crd(3, *), frc(3, *), &
        cn1(*), cn2(*), e14vir(3, 3)
! M-WJ, WJM, YD
    _REAL_ dipole(3, *), field(3, *), one_scee(*), one_scnb(*), pol(*), dipdamp
    _REAL_ dampfactor(*), pol2(*)
!
    integer iac(*), ico(*), ntypes, ix(*)
    _REAL_ ee14, enb14, epol14
    integer mytaskid, numtasks
#  include "extra_pts.h"
#  include "ew_cntrl.h"

    call do_14_dipole(charge, crd, frc, dipole, field, &
        iac, ico, ntypes, cn1, cn2, &
    ! M-WJ, WJM, YD
        pol, mpoltype, dampfactor, pol2, &
    !
        ee14, enb14, epol14, one_scee, one_scnb, &
        e14vir, ix(inb_14), numnb14, mytaskid, numtasks, scaldip)
    return
end subroutine get_14_dipole
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_14_dipole here]
subroutine do_14_dipole(charge, crd, frc, dipole, field, &
    iac, ico, ntypes, cn1, cn2, &
 ! M-WJ, WJM, YD
    pol, mpoltype, dampfactor, pol2, &
 !
    ee14, enb14, epol14, one_scee, one_scnb, &
    e14vir, nb_14_list, numnb14, mytaskid, numtasks, scaldip)
! M-WJ
    use constants, only : zero, one, two, three, four, five, six, twelve, third, half
!
    implicit none
    _REAL_ charge(*), crd(3, *), frc(3, *), &
        cn1(*), cn2(*), e14vir(3, 3)
! M-WJ, WJM, YD
    _REAL_ dipole(3, *), field(3, *), one_scee(*), one_scnb(*), pol(*), dipdamp
    _REAL_ dampfactor(*), pol2(*)
    integer mpoltype
!
    integer iac(*), ico(*), ntypes, nb_14_list(3, *), numnb14
    _REAL_ ee14, enb14, epol14
    integer mytaskid, numtasks, scaldip
#  include "flocntrl.h"

    integer n, i, j, ic, ic0, ia1, ia2, ibig, isml
    _REAL_ dx, dy, dz, r2, r2inv
    _REAL_ rinv
    _REAL_ scee0, scnb0, scdd0, f6, f12, r6, df
    _REAL_ term, term0, term1, termi, termj, termq, cgp
    _REAL_ dphii_dx, dphii_dy, dphii_dz, &
        dphij_dx, dphij_dy, dphij_dz
    _REAL_ c0, c1, c2, c3
    _REAL_ dfx, dfy, dfz
    _REAL_ cgi, cgj, dotir, dotjr, dotp, dotij
! M-WJ
    _REAL_ u3, au3, exp_au3, lambda3, lambda5, lambda7, c1_o
    _REAL_ u, au, a2u2, a3u3, a4u4, exp_au, v, v3, v4, termi_o, termj_o
!
    ee14 = 0.d0
    enb14 = 0.d0
    epol14 = 0.d0
    if (do_14 == 0) return

    scdd0 = 1.d0

    do j = 1, 3
        do i = 1, 3
            e14vir(i, j) = 0.d0
        end do
    end do
    do n = mytaskid + 1, numnb14, numtasks
        i = nb_14_list(1, n)
        j = nb_14_list(2, n)
        ic0 = nb_14_list(3, n)
        scee0 = one_scee(ic0)
        scnb0 = one_scnb(ic0)
        if (scaldip == 1) then
            scdd0 = one_scee(1)
        end if
        cgi = charge(i)
        cgj = charge(j)
        cgp = cgi*cgj
        dx = crd(1, j) - crd(1, i)
        dy = crd(2, j) - crd(2, i)
        dz = crd(3, j) - crd(3, i)
        r2 = dx*dx + dy*dy + dz*dz
        r2inv = 1.d0/r2
        rinv = sqrt(r2inv)

! M-WJ: calculate Thole damping function
!  Modified by WJM
        if (pol(i) == 0 .or. pol(j) == 0) then
            lambda3 = one
            lambda5 = one
            lambda7 = one
        else if (mpoltype == 1) then
            lambda3 = one
            lambda5 = one
            lambda7 = one
        else if (mpoltype == 2) then
            au3 = r2/rinv/(pol2(i)*pol2(j))
            exp_au3 = exp(-au3)
            lambda3 = one - exp_au3
            lambda5 = one - (one + au3)*exp_au3
            lambda7 = one - (one + au3 + three/five*au3*au3)*exp_au3
        else if (mpoltype == 3) then
            au = one/(pol2(i)*pol2(j))/rinv
            exp_au = exp(-au)
            a2u2 = au*au
            a3u3 = a2u2*au
            a4u4 = a3u3*au
            lambda3 = one - (a2u2/two + au + one)*exp_au
            lambda5 = lambda3 - a3u3/six*exp_au
            lambda7 = lambda5 - a4u4/six/five*exp_au
        else if (mpoltype == 4) then
            v = one/(pol2(i)*pol2(j))/rinv
            if (v < one) then
                v3 = v*v*v
                v4 = v3*v
                lambda3 = four*v3 - three*v4
                lambda5 = v4
                lambda7 = v4/five
            else
                lambda3 = one
                lambda5 = one
                lambda7 = one
            end if
        end if
!
        c0 = rinv
        c1 = c0*r2inv
        c2 = 3.d0*c1*r2inv
        c3 = 5.d0*c2*r2inv
! M-WJ damping B
        c1_o = c1
        c1 = c1*lambda3
        c2 = c2*lambda5
        c3 = c3*lambda7
!
        dotjr = dipole(1, j)*dx + dipole(2, j)*dy + dipole(3, j)*dz
        dotir = dipole(1, i)*dx + dipole(2, i)*dy + dipole(3, i)*dz
        dotp = dotjr*dotir
        dotij = dipole(1, i)*dipole(1, j) + dipole(2, i)*dipole(2, j) + &
            dipole(3, i)*dipole(3, j)
        term = cgj*dotir - cgi*dotjr + dotij
        term0 = cgp*c0 + term*c1 - dotp*c2
! M-WJ
!     termq = cgp*c1
        termq = cgp*c1_o
!
        term1 = term*c2 - dotp*c3
        termi = cgi*c1 + dotir*c2
        termj = cgj*c1 - dotjr*c2
! M-WJ
        termi_o = cgi*c1_o + dotir*c2
        termj_o = cgj*c1_o - dotjr*c2
!
        ee14 = ee14 + scee0*cgp*c0
        epol14 = epol14 + scdd0*(term*c1 - dotp*c2)
        dfx = (scee0*termq + scdd0*term1)*dx + &
            scdd0*(termi*dipole(1, j) - termj*dipole(1, i))
        dfy = (scee0*termq + scdd0*term1)*dy + &
            scdd0*(termi*dipole(2, j) - termj*dipole(2, i))
        dfz = (scee0*termq + scdd0*term1)*dz + &
            scdd0*(termi*dipole(3, j) - termj*dipole(3, i))
        r6 = r2inv*r2inv*r2inv
        !  always use the 6-12 parameters, even if 10-12 are available:
        ia1 = iac(i)
        ia2 = iac(j)
        ibig = max0(ia1, ia2)
        isml = min0(ia1, ia2)
        ic = ibig*(ibig - 1)/2 + isml
        f6 = cn2(ic)*r6
        f12 = cn1(ic)*(r6*r6)
        enb14 = enb14 + scnb0*(f12 - f6)
        df = (12.d0*f12 - 6.d0*f6)*r2inv
        dfx = dfx + scnb0*df*dx
        dfy = dfy + scnb0*df*dy
        dfz = dfz + scnb0*df*dz
        frc(1, j) = frc(1, j) + dfx
        frc(2, j) = frc(2, j) + dfy
        frc(3, j) = frc(3, j) + dfz
        frc(1, i) = frc(1, i) - dfx
        frc(2, i) = frc(2, i) - dfy
        frc(3, i) = frc(3, i) - dfz
! M-WJ
        dphii_dx = termj_o*dx + c1*dipole(1, j)
        dphii_dy = termj_o*dy + c1*dipole(2, j)
        dphii_dz = termj_o*dz + c1*dipole(3, j)
        dphij_dx = -termi_o*dx + c1*dipole(1, i)
        dphij_dy = -termi_o*dy + c1*dipole(2, i)
        dphij_dz = -termi_o*dz + c1*dipole(3, i)
!
        field(1, i) = field(1, i) - scdd0*dphii_dx
        field(2, i) = field(2, i) - scdd0*dphii_dy
        field(3, i) = field(3, i) - scdd0*dphii_dz
        field(1, j) = field(1, j) - scdd0*dphij_dx
        field(2, j) = field(2, j) - scdd0*dphij_dy
        field(3, j) = field(3, j) - scdd0*dphij_dz
#ifndef noVIRIAL
        e14vir(1, 1) = e14vir(1, 1) - dfx*dx
        e14vir(1, 2) = e14vir(1, 2) - dfx*dy
        e14vir(1, 3) = e14vir(1, 3) - dfx*dz
        e14vir(2, 1) = e14vir(2, 1) - dfy*dx
        e14vir(2, 2) = e14vir(2, 2) - dfy*dy
        e14vir(2, 3) = e14vir(2, 3) - dfy*dz
        e14vir(3, 1) = e14vir(3, 1) - dfz*dx
        e14vir(3, 2) = e14vir(3, 2) - dfz*dy
        e14vir(3, 3) = e14vir(3, 3) - dfz*dz
#endif
    end do  !  n = 1,numnb14
    return
end subroutine do_14_dipole
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fix_masses here]
subroutine fix_masses(natom, epowner, &
    nspm, nsp, tma, tmass, tmassinv, amass, amassinv)
    implicit none
    integer natom, epowner(*), nspm, nsp(*)
    _REAL_ tma(*), tmass, tmassinv, &
        amass(*), amassinv(*)
    integer n, k, l

    !     ---zero out mass and inverse masses for extra points;
    !        first fix amass,amassinv

    do n = 1, natom
        if (epowner(n) /= 0) then
            amass(n) = 0.d0
            amassinv(n) = 0.d0
        end if
    end do

    !     ---now redo tmass and tma

    tmass = 0.d0
    l = 0
    do k = 1, nspm
        tma(k) = 0.d0
        do n = 1, nsp(k)
            l = l + 1
            tma(k) = tma(k) + amass(l)
        end do
        tmass = tmass + tma(k)
    end do
    tmassinv = 1.d0/tmass
    return
end subroutine fix_masses
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trim_bonds here]
subroutine trim_bonds(ib, jb, icb, nbona, nbper, epowner)
    implicit none
    integer ib(*), jb(*), icb(*), nbona, nbper, epowner(*)
    integer ii, jj, n, m, nbper_new

    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_bonds: num bonds BEFORE trim =', &
        nbona, nbper

    m = 0
    do n = 1, nbona
        ii = (ib(n) + 3)/3
        jj = (jb(n) + 3)/3

        !     ---only keep if neither is extra

        if (epowner(ii) == 0 .and. epowner(jj) == 0) then
            m = m + 1
            ib(m) = ib(n)
            jb(m) = jb(n)
            icb(m) = icb(n)
        end if
    end do
    nbona = m

    if (nbper > 0) then
        do n = nbona + 1, nbona + nbper
            ii = (ib(n) + 3)/3
            jj = (jb(n) + 3)/3

            !     ---only keep if neither is extra

            if (epowner(ii) == 0 .and. epowner(jj) == 0) then
                m = m + 1
                ib(m) = ib(n)
                jb(m) = jb(n)
                icb(m) = icb(n)
                icb(m + nbper) = icb(n + nbper)
            end if
        end do
        nbper_new = m - nbona

        do n = 1, nbper - nbper_new
            icb(nbona + nbper_new + n) = icb(nbona + nbper + n)
        end do
        nbper = nbper_new
    end if

    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_bonds: num bonds AFTER  trim =', &
        nbona, nbper
    return
end subroutine trim_bonds
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trim_theta here]
subroutine trim_theta(it, jt, kt, ict, ntheta, ngper, epowner)
    implicit none
    integer it(*), jt(*), kt(*), ict(*), ntheta, ngper, epowner(*)
    integer ii, kk, n, m, ngper_new
    m = 0
    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_theta: num angle BEFORE trim =', &
        ntheta, ngper

    do n = 1, ntheta
        ii = (it(n) + 3)/3
        kk = (kt(n) + 3)/3

        !       ---only keep if neither is extra

        if (epowner(ii) == 0 .and. epowner(kk) == 0) then
            m = m + 1
            it(m) = it(n)
            jt(m) = jt(n)
            kt(m) = kt(n)
            ict(m) = ict(n)
        end if
    end do
    ntheta = m

    if (ngper > 0) then
        do n = ntheta + 1, ntheta + ngper
            ii = (it(n) + 3)/3
            kk = (kt(n) + 3)/3

            !       ---only keep if neither is extra

            if (epowner(ii) == 0 .and. epowner(kk) == 0) then
                m = m + 1
                it(m) = it(n)
                jt(m) = jt(n)
                kt(m) = kt(n)
                ict(m) = ict(n)
                ict(m + ngper) = ict(n + ngper)
            end if
        end do
        ngper_new = m - ntheta
        do n = 1, ngper - ngper_new
            ict(ntheta + ngper_new + n) = ict(ntheta + ngper + n)
        end do
        ngper = ngper_new
    end if

    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_theta: num angle AFTER  trim =', &
        ntheta, ngper
    return
end subroutine trim_theta
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trim_phi here]
subroutine trim_phi(ip, jp, kp, lp, icp, nphi, ndper, epowner)
    implicit none
    integer ip(*), jp(*), kp(*), lp(*), icp(*), nphi, ndper, epowner(*)
    integer ii, ll, n, m, ndper_new
    m = 0
    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_phi:  num diheds BEFORE trim =', &
        nphi, ndper

    do n = 1, nphi
        ii = (ip(n) + 3)/3
        ll = (iabs(lp(n)) + 3)/3

        !       ---only keep if neither is extra

        if (epowner(ii) == 0 .and. epowner(ll) == 0) then
            m = m + 1
            ip(m) = ip(n)
            jp(m) = jp(n)
            kp(m) = kp(n)
            lp(m) = lp(n)
            icp(m) = icp(n)
        end if
    end do
    nphi = m

    if (ndper > 0) then
        do n = nphi + 1, nphi + ndper
            ii = (ip(n) + 3)/3
            ll = (iabs(lp(n)) + 3)/3

            !       ---only keep if neither is extra

            if (epowner(ii) == 0 .and. epowner(ll) == 0) then
                m = m + 1
                ip(m) = ip(n)
                jp(m) = jp(n)
                kp(m) = kp(n)
                lp(m) = lp(n)
                icp(m) = icp(n)
                icp(m + ndper) = icp(n + ndper)
            end if
        end do
        ndper_new = m - nphi

        do n = 1, ndper - ndper_new
            icp(nphi + ndper_new + n) = icp(nphi + ndper + n)
        end do
        ndper = ndper_new
    end if

    write (6, '(a,2i6)') &
        '|      EXTRA_PTS, trim_phi:  num diheds AFTER  trim =', &
        nphi, ndper
    return
end subroutine trim_phi
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fix_degree_count here]
subroutine fix_degree_count(degrees)
    implicit none
    _REAL_ degrees
#  include "extra_pts.h"
    degrees = degrees - 3.d0*numextra
    return
end subroutine fix_degree_count
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine allocate_14nb here]
subroutine allocate_14nb(inb_14, numnb14)
    implicit none
    integer inb_14, numnb14
#  include "memory.h"
    integer i_ptr

    i_ptr = lasti
    call adj_mem_ptr(i_ptr, inb_14, 3*numnb14)
!3*numnb14 here since the array contains, ii, jj and ic0
    lasti = i_ptr
    return
end subroutine allocate_14nb
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine copy_14nb here]
subroutine copy_14nb(from14, to14, numnb14)
    implicit none
    integer from14(3, *), to14(3, *), numnb14
    integer n
    do n = 1, numnb14
        to14(1, n) = from14(1, n)
        to14(2, n) = from14(2, n)
        to14(3, n) = from14(3, n)
    end do
    return
end subroutine copy_14nb
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine allocate_frames here]
subroutine allocate_frames(numextra, ifrtyp, iatcen, inumep, &
    iepfr, ifrst, imid, ithrd, leploc)
    implicit none
    integer numextra
    integer ifrtyp, iatcen, inumep, iepfr, ifrst, imid, ithrd, leploc
#  include "memory.h"
    integer i_ptr, r_ptr, maxfr

    !     ---numextra is upper limit to maxfr

    maxfr = numextra
    i_ptr = lasti
    call adj_mem_ptr(i_ptr, ifrtyp, maxfr)
    call adj_mem_ptr(i_ptr, iatcen, maxfr)
    call adj_mem_ptr(i_ptr, inumep, maxfr)
    call adj_mem_ptr(i_ptr, iepfr, 2*maxfr)
    call adj_mem_ptr(i_ptr, ifrst, maxfr)
    call adj_mem_ptr(i_ptr, imid, maxfr)
    call adj_mem_ptr(i_ptr, ithrd, maxfr)
    lasti = i_ptr
    r_ptr = lastr
    call adj_mem_ptr(r_ptr, leploc, 3*2*maxfr)
    lastr = r_ptr
    return
end subroutine allocate_frames

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_pairs here]
subroutine do_pairs(ib, jb, npairs, list, num, maxp, &
    epowner, numnghbr, enghbrs, ifail)
    implicit none
    integer ib(*), jb(*), npairs, list(2, *), num, maxp
    integer enghbrs(5, *), numnghbr(3, *), epowner(*)
    integer ifail
    integer i, j, k, l, ii, jj, n
    ifail = 0
    do n = 1, npairs

        !       ---sometimes second index negative (improper dihedrals)

        ii = (ib(n) + 3)/3
        jj = (iabs(jb(n)) + 3)/3

        !       ---check neither is extra. count this bond and also extra attached

        if (epowner(ii) == 0 .and. epowner(jj) == 0) then
            k = numnghbr(3, ii)
            l = numnghbr(3, jj)
            if (num + 1 + k + l + k*l > maxp) goto 100
            num = num + 1
            list(1, num) = ii
            list(2, num) = jj
            do i = 1, k
                num = num + 1
                list(1, num) = enghbrs(i, ii)
                list(2, num) = jj
            end do
            do j = 1, l
                num = num + 1
                list(1, num) = ii
                list(2, num) = enghbrs(j, jj)
            end do
            do i = 1, k
                do j = 1, l
                    num = num + 1
                    list(1, num) = enghbrs(i, ii)
                    list(2, num) = enghbrs(j, jj)
                end do
            end do
        end if
    end do  !  n = 1,npairs
    return
100 ifail = 1
    return
end subroutine do_pairs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_14pairs here]
subroutine do_14pairs(ip, jp, kp, lp, icp, &
    nphi, nb_14_list, num, maxnb14, &
    epowner, numnghbr, enghbrs, ifail, chngmask)
    use qmmm_module, only : qmmm_nml, qmmm_struct
    implicit none
    integer ip(*), jp(*), kp(*), lp(*), icp(*)
    integer nphi, nb_14_list(3, *), num, maxnb14, ifail
    integer enghbrs(5, *), numnghbr(3, *), epowner(*), chngmask
    integer i, j, k, l, ii, ll, ic0, n, ni, nl, i3, l3
    logical qmmm_skip
    ifail = 0
    do n = 1, nphi
        i3 = ip(n)
        l3 = iabs(lp(n))
        ic0 = icp(n)
        k = kp(n)
        l = lp(n)
        ii = (i3 + 3)/3
        ll = (l3 + 3)/3
        if (qmmm_nml%ifqnt) then
            if (qmmm_struct%atom_mask(ii) .and. qmmm_struct%atom_mask(ll)) then !both true - Skip QM-QM 1-4 VDW
                qmmm_skip = .true.
            else
                qmmm_skip = .false.
            end if
        else
            qmmm_skip = .false.
        end if
        if ((k >= 0) .and. &
            (l >= 0) .and. &
            (.not. qmmm_skip) &
            ) then
            if (chngmask == 0) then
                if (num + 1 > maxnb14) then
                    ifail = 1
                    return
                end if
                num = num + 1
                nb_14_list(1, num) = ii
                nb_14_list(2, num) = ll
                nb_14_list(3, num) = ic0 !needed to be able to look up scnb and scee from arrays.
            else
                !       ---check neither is extra. count this bond and also
                !          extra attached
                if (epowner(ii) == 0 .and. epowner(ll) == 0) then
                    ni = numnghbr(3, ii)
                    nl = numnghbr(3, ll)
                    if (num + 1 + ni + nl + ni*nl > maxnb14) then
                        ifail = 1
                        return
                    end if
                    num = num + 1
                    nb_14_list(1, num) = ii
                    nb_14_list(2, num) = ll
                    nb_14_list(3, num) = ic0 !needed to be able to look up scnb and scee from arrays.
                    do i = 1, ni
                        num = num + 1
                        nb_14_list(1, num) = enghbrs(i, ii)
                        nb_14_list(2, num) = ll
                        nb_14_list(3, num) = ic0  !needed to be able to look up scnb and scee from arrays.
                    end do
                    do j = 1, nl
                        num = num + 1
                        nb_14_list(1, num) = ii
                        nb_14_list(2, num) = enghbrs(j, ll)
                        nb_14_list(3, num) = ic0  !needed to be able to look up scnb and scee from arrays.
                    end do
                    do i = 1, ni
                        do j = 1, nl
                            num = num + 1
                            nb_14_list(1, num) = enghbrs(i, ii)
                            nb_14_list(2, num) = enghbrs(j, ll)
                            nb_14_list(3, num) = ic0  !needed to be able to look up scnb and scee from arrays.
                        end do
                    end do
                end if
            end if
        end if
    end do
    return
end subroutine do_14pairs
