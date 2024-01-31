#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Parse and initialize data in cards and groups from the mdin
subroutine rgroup(natom, natc, nres, ngrp, ipres, lbres, igraph, isymbl, &
    itree, igroup, jgroup, indx, irespw, npdec, &
    weit, xc, konst, dotgtmd, belly, idecomp, nfu, writeout)

    implicit none
#  include "rgroup.h"
    integer:: i, i1, i2, idecomp, iend, ifld, igroup, igrp, iii, &
        indx, ipres, irespw, isear, isrch, istart, itime, ivar, izero, j, &
        jfld, jgroup, jgrp, k, l, lfind, lsign, natc, natmg, natom, nf, &
        nfu, ngrp, npdec, nres
    _REAL_ :: fvar, weit, wt, xc

    !     Mods for Rev A by gls:
    !     - cpp selectable precision
    !     - wildcard functionality fixed. (see findp)
    !     - format change to print less garbage
    !     - unreferenced statement labels removed
    !     - changed title write from 12 to 19 A4. (not 20 to avoid wraps)
    !     - changed to read from unit nf -- dac 12/90
    !     - made # of fields in a line an adjustable parameter, rather than
    !       hardwiring it to 7. This parameter is found in rgroup.h and must
    !       be included here and rfree.f. This parameter will now serve as
    !       the limit for how many fields are read in (rather than 80 as it was)
    !       so ifld and jfld no longer have to end with a 0. -- JMS 11/2010

    logical konst, misc, belly, flres, frres, dotgtmd
    character(len=4) lbres, igraph, isymbl, itree

    ! Ross Walker modifications for dipole printing - flag to control writing
    ! of group data to output file
    ! if (writeout) then  - wrapper placed around all write(6,*) statements
    logical writeout

    !     ----- READ IN GROUPS EACH ATOM IS IN -----

    !           WILL KEEP READING IN GROUPS OF CARDS UNTIL A BLANK CARD IS
    !           READ
    !           ALL ATOMS IN THIS CARD GROUP WILL BE PUT IN THE SAME GROUP
    !           THE ONLY EXCEPTION IS A RES CARD WHICH STARTS WITH A -
    !           RES NO. THE RESIDUES IN THIS GROUP WILL ALL BE PUT IN INDIV.
    !           GROUPS ANY OTHER GROUPS IN THIS SECTION MUST ALSO START
    !           WITH A  I.E. THEY MUST ALSO BE INDIV. GROUPS
    !           ANY TIME A "FIND" CARD IS READ, ONLY THOSE ATOMS THAT MEET
    !           THE SPECIFICATIONS ON AT LEAST 1 OF THE FIND CARDS
    !           WILL BE INCLUDED IN THE GROUP
    !           THIS OPTION MAY BE ENDED BY READING IN "NFIND",TERMINATING
    !           THE INPUT FOR THE GROUPS, OR BY READING IN ANOTHER "FIND"

    integer, parameter :: max_finds = 100

    character(len=4), dimension(1:max_finds) :: jgraph, jresnm, jsymbl, jtree
    common/propf/jgraph, jresnm, jsymbl, jtree, isrch
    character(len=4) title1(20)

    dimension ipres(*), igroup(*), jgroup(*), indx(*), irespw(*), &
        weit(*), igrp(MAX_RES + 1), jgrp(MAX_RES + 1)
    dimension lbres(*), igraph(*), isymbl(*), itree(*), xc(*)
    character(len=4) ihol, katn, ifind, nfind, iiend, ksear, ires, ilres, irres, ktypg
    ! JMS I think only ifld and ivar need to be extended to MAX_RES_FIELDS...
    dimension ifld(MAX_RES_FIELDS), jfld(20), ihol(20), &
        ivar(MAX_RES_FIELDS), fvar(20)
    character(len=4) igrap(40), isymb(40), imolc(40), iwild, jwild

    ! Can't use data statement since MAX_RES_FIELDS is variable now... JMS
    do i = 1, MAX_RES_FIELDS
        ifld(i) = 2
        if (i <= 2 .or. i == MAX_RES_FIELDS) ifld(i) = 0
    end do

    data jfld/4*1, 16*0/
    data katn/'ATOM'/
    data ifind/'FIND'/
    data nfind/'NFIN'/
    data iiend/'END '/
    data ksear/'SEAR'/
    data ires/'RES '/
    data ilres/'LRES'/
    data irres/'RRES'/
    data isear/40/
    data igrap/"C   ", "O   ", "N   ", "H   ", "CA  ", "HA  ", "HA2 ", "HA3 ", & ! proteins
        "N   ", "H1  ", "H2  ", "H3  ", "HA  ", "O   ", "OXT ", &        ! proteins N/C-terminal
        "P   ", "O1P ", "O2P ", "O5' ", "C5' ", "H5'1", "H5'2", &
        "C4' ", "H4' ", "O4' ", "C1' ", "H1' ", "C3' ", "H3' ", &
        "C2' ", "H2'1", "H2'2", "H2'1", "O2' ", "HO'2", "O3' ", &        ! DNA/RNA
        "H5T ", "O5' ", "O3' ", "H3T "/                               ! DNA/RNA 5'/3'-terminal
    data isymb/"C   ", "O   ", "N   ", "H   ", "CT  ", "H1  ", "H1  ", "H1  ", & ! proteins
        "N3  ", "H   ", "H   ", "H   ", "HP  ", "O2  ", "O2  ", &        ! proteins N/C-terminal
        "P   ", "O2  ", "O2  ", "OS  ", "CT  ", "H1  ", "H1  ", &
        "CT  ", "H1  ", "OS  ", "CT  ", "H2  ", "CT  ", "H1  ", &
        "CT  ", "HC  ", "HC  ", "H1  ", "OH  ", "HO  ", "OS  ", &        ! DNA/RNA
        "HO  ", "OH  ", "OH  ", "HO  "/                               ! DNA/RNA 5'/3'-terminal
    data imolc/"prot", "prot", "prot", "prot", "prot", "prot", "prot", "prot", & ! proteins
        "prot", "prot", "prot", "prot", "prot", "prot", "prot", &        ! proteins N/C-terminal
        "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", &
        "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", &
        "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", "nucs", &        ! DNA/RNA
        "nucs", "nucs", "nucs", "nucs"/                               ! DNA/RNA 5'/3'-terminal
    data iwild/'*   '/
    data jwild/'*   '/

    !     ----- RES CARD LISTS RESIDUE GROUPS IGRP(I) TO JGRP(I) ----

    !           IF 1ST VALUE IS NEGATIVE, THEN EVERY RESIDUE FROM IGRP(I)
    !           TO JGRP(I) IS CONSIDERED A SEPARATE GROUP
    !           IF 2ND VALUE   = 0, THEN SUBGROUP CONSISTS OF 1 RESIDUE
    !           ATOM CARD READS IN GROUPS OF ATOMS
    !           IF 2ND ATOM IN EACH PAIR  = 0 , THEN SUBGROUP CONSISTS OF
    !           1 ATOM RES AND ATOM CARDS MAY BE READ IN ANY ORDER
    !           END INPUT WITH AN "END " CARD

    !           ZERO NGRP BEFORE CALLING THIS ROUTINE
    !           ROUTINE WILL RETURN WITH THE NUMBER OF THE LAST GROUP READ

    ngrp = 0
    npdec = 0
    flres = .false.
    frres = .false.
    nf = nfu
    if (nf <= 0) nf = 5

    !     ----- INITIALISE THE GROUP ARRAY -----

    do i = 1, natom
        igroup(i) = 0
        if (konst) weit(i) = 0.0d0
    end do

22  continue
    itime = 0
    lsign = 0
    izero = 0
    isrch = 0
    natmg = 0

    !       ----- READ DIFFERENT GROUPS -----

    read (nf, 9208) (title1(k), k=1, 19)
    if (title1(1) == iiend) then
        if (writeout) write (6, '(4x,a)') '----- END OF GROUP READ -----'
        goto 900
    end if
    ngrp = ngrp + 1
    if (writeout) then
        write (6, '(4x,a,i5,a)') '----- READING GROUP ', ngrp, '; TITLE:'
        write (6, 9218) (title1(k), k=1, 19)
    end if

    !       ----- IF CONSTRAINED GROUPS READ THE WEIGHT FOR EACH GROUP -----
    !       ----- NOTE NO WEIGHT FOR TARGETED MD -----

    if (konst) then
        ifld(1) = 3
        ifld(2) = 0
        call rfree(ifld, ihol, ivar, fvar, nf, 6)
        wt = fvar(1)
        if (writeout) write (6, 9018) ngrp, wt
    end if

10  continue
    !       ----- READ THE GROUP CARDS -----

    ifld(1) = 1
    ifld(2) = 2
    call rfree(ifld, ihol, ivar, fvar, nf, 6)
    ktypg = ihol(1)
    k = 1
    do i = 1, MAX_RES ! extend to MAX_RES, JMS
        igrp(i) = ivar(k)
        jgrp(i) = ivar(k + 1)
        k = k + 2
    end do
    if (ktypg == iiend) goto 16
    if (idecomp > 0 .and. &
        (ktypg /= ires .and. ktypg /= ilres .and. ktypg /= irres)) &
        goto 36
    if (idecomp > 0 .and. ktypg == ilres) flres = .true.
    if (idecomp > 0 .and. ktypg == irres) frres = .true.
    if (ktypg == nfind) then
        if (writeout) write (6, 199)
        isrch = 0
        goto 10
    end if
    if (ktypg == ifind) then

        !         ----- FIND OPTION ... READ THE ATOM SPECIFICS -----

        if (writeout) write (6, 200)
        do iii = 1, max_finds
            call rfree(jfld, ihol, ivar, fvar, nf, 6)
            jgraph(iii) = ihol(1)
            jsymbl(iii) = ihol(2)
            jtree(iii) = ihol(3)
            jresnm(iii) = ihol(4)
            if (jgraph(iii) == ksear) exit
            if (writeout) write (6, 202) jgraph(iii), jsymbl(iii), jtree(iii), jresnm(iii)
        end do

        isrch = iii - 1
        if (isrch > max_finds) then
            if (writeout) write (6, 66) isrch
        end if
        !         ----- NOW READ IN RES AND ATOMS TO BE SEARCHED -----

        goto 10
    end if
    itime = itime + 1

    if (idecomp > 0 .and. &
        (ktypg == ilres .or. ktypg == irres)) then

        !         ----- CHECK LRES or RRES CARD -----

        do i = 1, MAX_RES ! Extend this to MAX_RES, JMS
            !           --- Sign does not play a role ---
            i1 = igrp(i)
            if (i1 == 0) goto 10
            if (i1 < 0) i1 = -i1
            if (i1 > nres) goto 36
            i2 = jgrp(i)
            if (i2 < 0) i2 = -i2
            if (i2 > nres) i2 = nres
            do j = i1, i2
                if (idecomp > 2) then
                    npdec = npdec + 1
                    indx(j) = npdec
                    irespw(npdec) = j
                end if
                istart = ipres(j)
                iend = ipres(j + 1) - 1
                do k = istart, iend
                    !               --- Find backbone atoms ---
                    do l = 1, isear
                        jgraph(l) = igrap(l)
                        jresnm(l) = imolc(l)
                        jsymbl(l) = isymb(l)
                        jtree(l) = jwild
                    end do
                    isrch = isear
                    call findp(k, j, lfind, nres, ipres, lbres, isymbl, itree, igraph)
                    if (lfind == 1) then
                        if (ktypg == irres) then
                            jgroup(k) = nres + j
                        else !Ligand atom
                            jgroup(k) = -nres - j
                        end if
                    else
                        !                 --- Store as sidechain atoms ---
                        if (ktypg == irres) then
                            jgroup(k) = j
                        else !Ligand atom
                            jgroup(k) = -j
                        end if
                    end if
                end do
            end do
        end do
        ngrp = ngrp - 1
        goto 10
    end if

    if (ktypg /= katn) then

        !         ----- CHECK RES CARD -----

        !         ----- 1ST GROUP OF 1ST CARD MUST BE - IF ANY - NUMBERS ARE
        !               FOUND -----

        if (itime == 1 .and. igrp(1) < 0) lsign = 1
        do i = 1, MAX_RES ! Extend this to MAX_RES, JMS
            i1 = igrp(i)
            if (i1 == 0) goto 10
            i2 = jgrp(i)
            if (i2 > nres) i2 = nres
            if (i1 > 0 .and. lsign == 1) goto 36
            if (i1 < 0 .and. lsign == 0) goto 36
            if (i1 < 0) i1 = -i1
            if (i2 <= 0) i2 = i1
            if (lsign == 0) then
                if (writeout) write (6, 14) ngrp, i1, i2
            end if
            do j = i1, i2
                istart = ipres(j)
                iend = ipres(j + 1) - 1
                do k = istart, iend
                    if (isrch > 0) &
                        call findp(k, j, lfind, nres, ipres, lbres, isymbl, itree, igraph)
                    if (isrch > 0 .and. lfind == 0) cycle
                    igroup(k) = ngrp
                    if (konst) weit(k) = wt
                    natmg = natmg + 1
                end do
                if (lsign == 1) then
                    if (writeout) write (6, 46) ngrp, j
                end if
                if (lsign == 1) ngrp = ngrp + 1
            end do
        end do
        goto 10
    end if

    !       ----- ATOM TYPE CONSTRAINTS -----

    if (lsign == 1) goto 36
    if (writeout) write (6, 51) ngrp
    do i = 1, MAX_RES ! Extend this to MAX_RES, JMS
        i1 = igrp(i)
        if (i1 < 0) goto 36
        if (i1 == 0) goto 10
        i2 = jgrp(i)
        if (i2 > natom) i2 = natom
        if (i2 <= 0) i2 = i1
        if (writeout) write (6, 52) i1, i2
        do j = i1, i2
            if (isrch > 0) &
                call findp(j, izero, lfind, nres, ipres, lbres, isymbl, itree, igraph)
            if (isrch > 0 .and. lfind == 0) cycle
            natmg = natmg + 1
            igroup(j) = ngrp
            if (konst) weit(j) = wt
        end do
    end do
    goto 10

16  continue
    if (lsign == 1) ngrp = ngrp - 1
    if (itime == 0) ngrp = ngrp - 1
    if (writeout) write (6, 222) natmg
    goto 22

36  continue
    if (writeout) write (6, 127) ktypg, (igrp(i), jgrp(i), i=1, 7)
    goto 10

    !     ----- ALL GROUPS ARE READ RETURN -----

900 continue
    if (konst .or. dotgtmd) then

        !       ----- GATHER ALL THE CONSTRAINED ATOMS TOGETHER -----

        natc = 0
        do i = 1, natom
            if (igroup(i) <= 0) cycle
            natc = natc + 1
            igroup(natc) = i

            ! WEIT will not be used for targeted MD

            weit(natc) = weit(i)
        end do

    else if (idecomp > 0) then

        !       ----- Special treatment for energy decomposition -----

        if (.not. flres .and. .not. frres) then
            !         --- Assign all atoms to "Protein" ---
            !         --- Set all residues to be printed ---
            do j = 1, nres
                if (idecomp > 2) then
                    npdec = npdec + 1
                    indx(j) = npdec
                    irespw(npdec) = j
                end if
                istart = ipres(j)
                iend = ipres(j + 1) - 1
                do k = istart, iend
                    igroup(k) = 1
                    !              --- Find backbone atoms ---
                    do l = 1, isear
                        jgraph(l) = igrap(l)
                        jresnm(l) = imolc(l)
                        jsymbl(l) = isymb(l)
                        jtree(l) = jwild
                    end do
                    isrch = isear
                    call findp(k, j, lfind, nres, ipres, lbres, isymbl, itree, igraph)
                    if (lfind == 1) then
                        jgroup(k) = nres + j
                    else
                        !                --- Store as sidechain atoms ---
                        jgroup(k) = j
                    end if
                end do
            end do
        end if

        !       --- Print assignment of atoms ---
        if (writeout) then
            do j = 1, nres
                istart = ipres(j)
                iend = ipres(j + 1) - 1
                do k = istart, iend
                    write (6, '(a,i5,a,i5,a,i5,i5)') &
                        'Atom ', k, ' (', j, ') : ', jgroup(k), igroup(k)
                end do
            end do
        end if

    else if (.not. belly) then

        !       ----- PUT THE ATOMS WHICH ARE NOT IN THE DEFINED GROUPS
        !             AS THE LAST GROUP -----

        misc = .false.
        do i = 1, natom
            if (igroup(i) /= 0) cycle
            misc = .true.
            igroup(i) = ngrp + 1
        end do
        if (misc) ngrp = ngrp + 1
    end if

    ! for targeted MD, support only 1 group for now

    if (dotgtmd .and. ngrp > 1) then
        write (6, 9200)
        write (6, 9205)
        call mexit(6, 1)
9200    format("ERROR IN TARGETED MD GROUP INPUT (ITGTMD=1)")
9205    format("ONLY 1 GROUP ALLOWED")
    end if

199 format(6x, 'END OF ATOM SPECIFICATION',/)
200 format(6x, 'ALL ATOMS THAT MEET 1 OF THE FOLLOWING', &
        ' SPECIFICATIONS WILL BE INCLUDED IN GROUP BELOW',/)
202 format(6x, 'GRAPH NAME  = ', a4, 2x, 'SYMBOL  = ', a2, 4x, &
        'TREE SYMBOL  = ', a1, 5x, 'RESIDUE TYPE  = ', a4,/)
66  format(6x, '**** NUMBER OF FIND CARDS  = ', i5, 2x, &
        'IS TOO BIG ******',/)
14  format(' GRP', i5, ' RES', i5, ' TO ', i5)
46  format(6x, 'GROUP', i5, 3x, 'CONSISTS OF RESIDUE', i5,/)
51  format(6x, 'GROUP', i5, 3x, 'CONSISTS OF ATOMS -',/)
52  format(35x, i5, 2x, 'TO', i5)
222 format(6x, 'Number of atoms in this group  = ', i5)
127 format(6x, '***PROBLEMS WITH GROUP', a4, 14i5, '*******',/)
9018 format(/5x, 'GROUP', i5, ' HAS HARMONIC CONSTRAINTS', f12.5)
9208 format(20a4)
9218 format(1x, 20a4)
    return

contains

    subroutine findp(iatom, ires, iloc, nres, ipres, lbres, isymbl, itree, igraph)

        implicit none
        integer:: iatom, iloc, ipres, ires, n, nres
        character(len=4) lbres, isymbl, itree, igraph, iwild, jwild

        !     Rev A mod (G. Seibel, Apr 89)
        !     Changed iblank, jblank to iwild, jwild, to give the wildcard
        !     functionality promised by the doc.
        !     isymbl() dimensioned. (this was long-standing bug)

        dimension ipres(*), lbres(*), itree(*), igraph(*), isymbl(*)
        character(len=92)  :: ipeps = 'ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP '// &
            'ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL '
        character(len=145) :: inucs = 'DA  DC  DG  DT  DA5 DC5 DG5 DT5 DA3 DC3 DG3 DT3 '// &
            'RA  RC  RG  RU  RA5 RC5 RG5 RU5 RA3 RC3 RG3 RU3 '// &
            ' A   C   G   U   A5  C5  G5  U5  A3  C3  G3  U3 '

        !     ----- CHECKS IF A GIVEN ATOM HAS CERTAIN CHARACTERISTICS -----

        iwild = '*   '
        jwild = '*   '
        iloc = 0
        if (ires == 0) call findrs(iatom, ires, nres, ipres)
        do n = 1, isrch
            if ((idecomp == 0) .and. (jresnm(n) /= iwild) .and. (jresnm(n) /= lbres(ires))) cycle
            if ((jgraph(n) /= iwild) .and. (jgraph(n) /= igraph(iatom))) cycle
            if ((jtree(n) /= jwild) .and. (jtree(n) /= itree(iatom))) cycle
            if ((jsymbl(n) /= jwild) .and. (jsymbl(n) /= isymbl(iatom))) cycle
            if (idecomp > 0) then
                if ((jresnm(n) .eq. 'prot') .and. (index(ipeps, lbres(ires)) == 0)) cycle
                if ((jresnm(n) .eq. 'nucs') .and. (index(inucs, lbres(ires)) == 0)) cycle
            end if
            iloc = 1
            exit
        end do
        return
    end subroutine findp
!-----------------------------------------------------------------------
end subroutine rgroup
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine findrs here]
subroutine findrs(numa, ires, nres, ipres)
    implicit none
    integer:: i, im, ipres, ires, nres, numa
    dimension ipres(*)

    if (numa <= 0) goto 11
    if (numa >= ipres(nres)) ires = nres
    if (numa >= ipres(nres)) return
    im = nres - 1
    do i = 1, im
        if (numa >= ipres(i) .and. numa < ipres(i + 1)) goto 12
    end do
11  write (6, 100) numa
    goto 200
12  continue
    ires = i
    return
200 continue
100 format(/2x, 'PROBLEMS FINDING RESIDUE OF ATOM ', i5)
    call mexit(6, 1)
end subroutine findrs
