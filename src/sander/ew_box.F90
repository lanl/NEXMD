! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

!     -----------------------------------------------------------------
!     All of the particle mesh Ewald code was written and contributed
!     by Tom Darden from the National Institute of Environmental Health
!     Sciences division of the NIH.  Originally written with a modified
!     version of AMBER 3A, the code was updated during the summer of 1994
!     to be compatible with AMBER 4.1.
!     -----------------------------------------------------------------

!     The routines herein are used in the particle mesh Ewald code
!     specifically for pressure scaling, volume calculation, imaging,
!     and other things related to the manipulating the periodic box

!---------------------------------------------------------------------
!     --- EW_PSCALE ---
!     ...pressure scaling routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ew_pscale here]
subroutine ew_pscale(natom, x, amass, nummols, molsiz, npscal)
    use nblist, only : oldrecip, ucell
    implicit none
    _REAL_ x(3, *), amass(*)
    integer nummols, molsiz(*), natom, npscal

    integer n
    _REAL_ f1, f2, f3

    integer imol, iatom, i, num, isave
    _REAL_ mass, massmol, xmol, ymol, zmol, fm1, fm2, fm3, &
        xmolnu, ymolnu, zmolnu

    if (npscal == 1) then

        !     ...apply CENTER OF MOLECULE BASED PRESSURE SCALING

        i = 0
        do imol = 1, nummols
            massmol = 0.d0
            xmol = 0.d0
            ymol = 0.d0
            zmol = 0.d0
            num = molsiz(imol)
            isave = i

            !   ...get c.o.m. of molecule - TODO: use precalced masses

            do iatom = 1, num
                i = i + 1
                mass = amass(i)
                massmol = massmol + mass
                xmol = xmol + mass*x(1, i)
                ymol = ymol + mass*x(2, i)
                zmol = zmol + mass*x(3, i)
            end do
            xmol = xmol/massmol
            ymol = ymol/massmol
            zmol = zmol/massmol

            !   ...now get fracs for c.o.m. using old cell params

            fm1 = xmol*oldrecip(1, 1) + ymol*oldrecip(2, 1) + &
                zmol*oldrecip(3, 1)
            fm2 = xmol*oldrecip(1, 2) + ymol*oldrecip(2, 2) + &
                zmol*oldrecip(3, 2)
            fm3 = xmol*oldrecip(1, 3) + ymol*oldrecip(2, 3) + &
                zmol*oldrecip(3, 3)

            !   ...use these with new cell params to get new c.o.m. cartesians

            xmolnu = fm1*ucell(1, 1) + fm2*ucell(1, 2) + fm3*ucell(1, 3)
            ymolnu = fm1*ucell(2, 1) + fm2*ucell(2, 2) + fm3*ucell(2, 3)
            zmolnu = fm1*ucell(3, 1) + fm2*ucell(3, 2) + fm3*ucell(3, 3)
            i = isave

            !   ...now rigidly translate molecule

            do iatom = 1, num
                i = i + 1
                x(1, i) = x(1, i) + xmolnu - xmol
                x(2, i) = x(2, i) + ymolnu - ymol
                x(3, i) = x(3, i) + zmolnu - zmol
            end do
        end do  !  imol = 1,nummols

    else

        !     ...apply ATOMIC BASED PRESSURE SCALING

        do n = 1, natom

            !     ...get fractional coordinates with old cell parameters

            f1 = x(1, n)*oldrecip(1, 1) + x(2, n)*oldrecip(2, 1) + &
                x(3, n)*oldrecip(3, 1)
            f2 = x(1, n)*oldrecip(1, 2) + x(2, n)*oldrecip(2, 2) + &
                x(3, n)*oldrecip(3, 2)
            f3 = x(1, n)*oldrecip(1, 3) + x(2, n)*oldrecip(2, 3) + &
                x(3, n)*oldrecip(3, 3)

            !     ...use these with the new cell parameters to get new
            !     cartesian coordinates

            x(1, n) = f1*ucell(1, 1) + f2*ucell(1, 2) + f3*ucell(1, 3)
            x(2, n) = f1*ucell(2, 1) + f2*ucell(2, 2) + f3*ucell(2, 3)
            x(3, n) = f1*ucell(3, 1) + f2*ucell(3, 2) + f3*ucell(3, 3)
        end do
    end if  ! ( npscal == 1 )
    return
end subroutine ew_pscale
!---------------------------------------------------------------------

!     --- DOT ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dot here]
subroutine dot(v1, v2, result)
    _REAL_ v1(3), v2(3), result
    result = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    return
end subroutine dot
!---------------------------------------------------------------------

!     --- CROSS ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cross here]
subroutine cross(v1, v2, v12)

    !    v12 is cross product of v1 and v2

    _REAL_ v1(3), v2(3), v12(3)
    v12(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v12(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v12(3) = v1(1)*v2(2) - v1(2)*v2(1)
    return
end subroutine cross
!---------------------------------------------------------------------

!     --- FILL_UCELL ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_ucell here]
subroutine fill_ucell(ta, tb, tc, talpha, tbeta, tgamma)
    use nblist, only : a, b, c, alpha, beta, gamma, ucell, recip, dirlng, reclng, &
        sphere, volume
    implicit none
#  include "ew_cntrl.h"
#  include "box.h"
    _REAL_ ta, tb, tc, talpha, tbeta, tgamma
    a = ta
    b = tb
    c = tc
    box(1) = a
    box(2) = b
    box(3) = c
    alpha = talpha
    beta = tbeta
    gamma = tgamma
    call get_ucell(a, b, c, alpha, beta, gamma, &
        ucell, recip, dirlng, reclng, sphere, volume)
    return
end subroutine fill_ucell
!---------------------------------------------------------------------

!     --- GET_UCELL ---

!     ...this routine produces the direct and reciprocal lattice
!     vectors from the unit cell edge lengths and angles which are
!     passed to it.  It is assumed that the 1st vector (length a)
!     lies along the cartesian x-axis the 2nd vector (length b) is
!     in the x-y plane with positive y, and that the direct lattice
!     vectors are a non-degenerate right handed system.  Thus the 3rd
!     vector has positive z component.  Alpha is the angle (in degrees)
!     between 2nd and 3rd vectors, beta is the angle (in degrees)
!     between 1st and 3rd vectors, and gamma is the angle (in degrees)
!     between 1st and 2nd vectors.  The lengths of the direct lattice
!     vectors are given by dirlng(1), dirlng(2) and dirlng(3),
!     whereas the reciprocal lengths of the reciprocal vectors are
!     reclng(1),reclng(2)  and reclng(3).

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_ucell here]
subroutine get_ucell(a, b, c, alpha, beta, gamma, &
    ucell, recip, dirlng, reclng, sphere, volume)
    use constants, only : DEG_TO_RAD
    implicit none
#  include "extra.h"
    _REAL_ a, b, c, alpha, beta, gamma
    _REAL_ sphere, volume
    _REAL_ ucell(3, 3), recip(3, 3), dirlng(3), reclng(3)
    _REAL_ u23(3), u31(3), u12(3)
    _REAL_ result, distance, onevolume
    integer i, j

    ucell(1, 1) = a
    ucell(2, 1) = 0.d0
    ucell(3, 1) = 0.d0
    ucell(1, 2) = b*cos(DEG_TO_RAD*gamma)
    ucell(2, 2) = b*sin(DEG_TO_RAD*gamma)
    ucell(3, 2) = 0.d0
    ucell(1, 3) = c*cos(DEG_TO_RAD*beta)
    ucell(2, 3) = &
        (b*c*cos(DEG_TO_RAD*alpha) - ucell(1, 3)*ucell(1, 2))/ucell(2, 2)
    ucell(3, 3) = sqrt(c*c - ucell(1, 3)*ucell(1, 3) - &
        ucell(2, 3)*ucell(2, 3))
    dirlng(1) = a
    dirlng(2) = b
    dirlng(3) = c

    !  now get reciprocal vectors

    call cross(ucell(1, 2), ucell(1, 3), u23)
    call cross(ucell(1, 3), ucell(1, 1), u31)
    call cross(ucell(1, 1), ucell(1, 2), u12)
    call dot(ucell(1, 1), u23, volume)
    onevolume = 1.0d0/volume
    do j = 1, 3
        recip(j, 1) = u23(j)*onevolume
        recip(j, 2) = u31(j)*onevolume
        recip(j, 3) = u12(j)*onevolume
    end do

    reclng(1) = 1.d0/sqrt(recip(1, 1)*recip(1, 1) + &
        recip(2, 1)*recip(2, 1) + &
        recip(3, 1)*recip(3, 1))
    reclng(2) = 1.d0/sqrt(recip(1, 2)*recip(1, 2) + &
        recip(2, 2)*recip(2, 2) + &
        recip(3, 2)*recip(3, 2))
    reclng(3) = 1.d0/sqrt(recip(1, 3)*recip(1, 3) + &
        recip(2, 3)*recip(2, 3) + &
        recip(3, 3)*recip(3, 3))

    ! interfacial distances given by dot of direct,recip
    ! sphere is radius of largest sphere inscribed in unit cell
    ! the minimum image cutoff must be less than or equal to this

    sphere = a + b + c
    do i = 1, 3
        call dot(recip(1, i), ucell(1, i), result)
        distance = result*reclng(i)
        if (distance < sphere) sphere = distance
    end do
    sphere = 0.5d0*sphere
    write (6, '(a,f9.3)') &
        '|Largest sphere to fit in unit cell has radius = ', sphere
    return
end subroutine get_ucell
!---------------------------------------------------------------------

!    ---- NEW_UCELL_GEN ------

! this is used by debug in checking the virial tensor
! eventually need this is or similar to do full 3x3 scaling

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine new_ucell_gen here]
subroutine new_ucell_gen(anew)
    use nblist, only : ucell, recip, olducell, oldrecip, dirlng, reclng, &
        volume
    implicit none
    _REAL_ anew(3, 3)

#  include "ew_cntrl.h"
    _REAL_ u23(3), u31(3), u12(3)
    integer i, j
    do j = 1, 3
        do i = 1, 3
            olducell(i, j) = ucell(i, j)
            oldrecip(i, j) = recip(i, j)
        end do
    end do
    do j = 1, 3
        do i = 1, 3
            ucell(i, j) = anew(i, j)
        end do
    end do
    !  now get reciprocal vectors
    dirlng(1) = sqrt(ucell(1, 1)*ucell(1, 1) + &
        ucell(2, 1)*ucell(2, 1) + ucell(3, 1)*ucell(3, 1))
    dirlng(2) = sqrt(ucell(1, 2)*ucell(1, 2) + &
        ucell(2, 2)*ucell(2, 2) + ucell(3, 2)*ucell(3, 2))
    dirlng(3) = sqrt(ucell(1, 3)*ucell(1, 3) + &
        ucell(2, 3)*ucell(2, 3) + ucell(3, 3)*ucell(3, 3))

    call cross(ucell(1, 2), ucell(1, 3), u23)
    call cross(ucell(1, 3), ucell(1, 1), u31)
    call cross(ucell(1, 1), ucell(1, 2), u12)
    call dot(ucell(1, 1), u23, volume)
    do j = 1, 3
        recip(j, 1) = u23(j)/volume
        recip(j, 2) = u31(j)/volume
        recip(j, 3) = u12(j)/volume
    end do
    reclng(1) = 1.d0/sqrt(recip(1, 1)*recip(1, 1) + &
        recip(2, 1)*recip(2, 1) + &
        recip(3, 1)*recip(3, 1))
    reclng(2) = 1.d0/sqrt(recip(1, 2)*recip(1, 2) + &
        recip(2, 2)*recip(2, 2) + &
        recip(3, 2)*recip(3, 2))
    reclng(3) = 1.d0/sqrt(recip(1, 3)*recip(1, 3) + &
        recip(2, 3)*recip(2, 3) + &
        recip(3, 3)*recip(3, 3))

    return
end subroutine new_ucell_gen
!---------------------------------------------------------------------

!     --- REDO_UCELL ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine redo_ucell here]
subroutine redo_ucell(factor)
    use nblist, only : ucell, recip, olducell, oldrecip
    implicit none
    _REAL_ factor(3)

#  include "ew_cntrl.h"
    integer i, j

    do j = 1, 3
        do i = 1, 3
            olducell(i, j) = ucell(i, j)
            oldrecip(i, j) = recip(i, j)
        end do
    end do
    call update_ucell_iso(factor, verbose)
    return
end subroutine redo_ucell
!---------------------------------------------------------------------

!     --- UPDATE_UCELL_ISO ---

!     ...scales unit cell uniformly by factor;
!     alpha,beta,gamma unchanged...

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine update_ucell_iso here]
subroutine update_ucell_iso(factor, verbose)
    use nblist, only : a, b, c, alpha, beta, gamma, ucell, recip, dirlng, reclng, &
        sphere, volume, olducell, oldrecip, cutlist
    implicit none
    _REAL_ factor(3)
    integer verbose
    _REAL_ facinv(3), result, distance
    integer i, j
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
#  include "extra.h"
#  include "box.h"

    a = a*factor(1)
    b = b*factor(2)
    c = c*factor(3)
    box(1) = a
    box(2) = b
    box(3) = c
    volume = volume*factor(1)*factor(2)*factor(3)

    if (verbose == 1 .and. master) then
        write (6, 99) a, b, c, volume
    end if
99  format(1x, 'a,b,c,volume now equal to ', 4f12.3)
    do i = 1, 3
        facinv(i) = 1.d0/factor(i)
        dirlng(i) = factor(i)*dirlng(i)
        reclng(i) = factor(i)*reclng(i)
    end do
    do j = 1, 3
        do i = 1, 3
            ucell(i, j) = factor(i)*ucell(i, j)
            recip(i, j) = facinv(i)*recip(i, j)
        end do
    end do
    sphere = a + b + c
    do i = 1, 3
        call dot(recip(1, i), ucell(1, i), result)
        distance = result*reclng(i)
        if (distance < sphere) sphere = distance
    end do
    sphere = 0.5d0*sphere
    ! check unit cell versus cutlist
    if (cutlist > sphere) then
        if (master) then
            write (6, *) 'Cutoff list exceeds largest sphere in unit cell!!'
            write (6, *) 'Big problems with imaging!!'
            write (6, *) 'a,b,c = ', a, b, c
            write (6, *) 'alpha,beta,gamma = ', alpha, beta, gamma
            write (6, *) 'cutlist,sphere = ', cutlist, sphere
        end if
        call mexit(6, 1)
    end if
    return
end subroutine update_ucell_iso
!---------------------------------------------------------------------

!     --- WRAP_MOLECULES ---

!     ...wrap the molecules/coordinates across the periodic box.
!     Geometric center of each molecule is checked to see if
!     it is within the unit cell or not...

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wrap_molecules here]
subroutine wrap_molecules(nummols, molsiz, crd)
    use nblist, only : ucell, recip
    implicit none
#ifdef MPI
#  include "parallel.h"
#endif

    integer nummols, molsiz(*)
    _REAL_ crd(3, *)

    integer i, j, lo, hi
    _REAL_ tran(3), f1, f2, f3, g1, g2, g3

    lo = 1
    do i = 1, nummols
        hi = lo + molsiz(i) - 1
        f1 = 0.0d0
        f2 = 0.0d0
        f3 = 0.0d0
        do j = lo, hi
            f1 = f1 + crd(1, j)*recip(1, 1) + crd(2, j)*recip(2, 1) + &
                crd(3, j)*recip(3, 1)
            f2 = f2 + crd(1, j)*recip(1, 2) + crd(2, j)*recip(2, 2) + &
                crd(3, j)*recip(3, 2)
            f3 = f3 + crd(1, j)*recip(1, 3) + crd(2, j)*recip(2, 3) + &
                crd(3, j)*recip(3, 3)
        end do
#if 0
        f1 = f1/molsiz(i)
        f2 = f2/molsiz(i)
        f3 = f3/molsiz(i)

        g1 = f1
        if (f1 < 0.d0) g1 = f1 + 1.d0
        if (f1 >= 1.d0) g1 = f1 - 1.d0
        g2 = f2
        if (f2 < 0.d0) g2 = f2 + 1.d0
        if (f2 >= 1.d0) g2 = f2 - 1.d0
        g3 = f3
        if (f3 < 0.d0) g3 = f3 + 1.d0
        if (f3 >= 1.d0) g3 = f3 - 1.d0
#else
        f1 = f1/molsiz(i) - 0.5d0
        f2 = f2/molsiz(i) - 0.5d0
        f3 = f3/molsiz(i) - 0.5d0

        g1 = f1 - anint(f1)
        g2 = f2 - anint(f2)
        g3 = f3 - anint(f3)

#endif
        if (f1 /= g1 .or. f2 /= g2 .or. f3 /= g3) then
            tran(1) = (g1 - f1)*ucell(1, 1) + (g2 - f2)*ucell(1, 2) + &
                (g3 - f3)*ucell(1, 3)
            tran(2) = (g1 - f1)*ucell(2, 1) + (g2 - f2)*ucell(2, 2) + &
                (g3 - f3)*ucell(2, 3)
            tran(3) = (g1 - f1)*ucell(3, 1) + (g2 - f2)*ucell(3, 2) + &
                (g3 - f3)*ucell(3, 3)
            do j = lo, hi
                crd(1, j) = crd(1, j) + tran(1)
                crd(2, j) = crd(2, j) + tran(2)
                crd(3, j) = crd(3, j) + tran(3)
            end do
#ifdef MPI
            ! DAN ROE: Since in the stripwat routine (hybrid REMD, multisander.f)
            !  all threads call wrap_molecules, make it so only the masters write.
            if (i == 1 .and. sanderrank == 0) write (6, '(a,3f15.5)') &
                'wrapping first mol.:', tran(1), tran(2), tran(3)
#else
            if (i == 1) write (6, '(a,3f15.5)') 'wrapping first mol.:', &
                tran(1), tran(2), tran(3)
#endif
        end if
        lo = hi + 1
    end do  !  i = 1,nummols
    return
end subroutine wrap_molecules
!---------------------------------------------------------------------

!     --- WRAP_TO ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wrap_to here]
subroutine wrap_to(nspm, nsp, xx, box)
    use constants, only : one, fourth, pi
    implicit none
    integer nsp(*), nspm
    _REAL_ xx(3, *), box(3)

    _REAL_ cx, cy, cz, x0, y0, z0, x, y, z, xt, yt, zt
    _REAL_ phi, cos1, sin1, cos2, sin2
    _REAL_ tobox, tobinv, facecoord
    _REAL_ t11, t12, t13, t21, t22, t23, t31, t32, t33
    integer m, i, m0
    _REAL_ sign
    _REAL_ anint

    !      The trunf. oct. has:
    !        * center at (0,0,0) where the corner of the triclinic cell is.
    !        * one hex face with the x axis for a normal
    !                       face is perpendicular to the x axis and
    !                       1/2 box(1) away from the origin.
    !        * one hex face perp to xy plane, second edge vector of the
    !                       triclinic cell is its normal, 109 degrees from
    !                       the x axis in xy plane (-x,+y quadrant)

    !      Approach to reconstruct the t.o. is to rotate the coordinates
    !         to put the hex faces in the (+-1,+-1,+-1) normal directions, and
    !         the square (diamond) faces perp to xyz axes. This is 3 rotations:
    !         We did +45 around z, +(90-tetra/2) around y, +90 around x to get
    !            the to oriented for the triclinic cell
    !         Now we do the opposite: -90(x), -(90-tetra/2)(y), -45(z)
    !            to reproduce the original orientation, map coords into
    !            a t.o. centered at origin, then do the rotations again
    !            to put it back in the orientation that matches the restrt.

    facecoord = box(1)/(2.d0*sqrt(3.d0))
    tobox = 2.d0*box(1)/sqrt(3.d0)
    tobinv = 1.d0/tobox
    phi = fourth*PI
    cos1 = cos(phi)
    sin1 = sin(phi)
    cos2 = sqrt(2.d0)/sqrt(3.d0)
    sin2 = 1.d0/sqrt(3.d0)
    t11 = cos2*cos1
    t12 = -cos2*sin1
    t13 = -sin2
    t21 = -sin2*cos1
    t22 = sin2*sin1
    t23 = -cos2
    t31 = sin1
    t32 = cos1
    t33 = 0

    m0 = 1
    do m = 1, nspm
        !     ------ calculate center of geometry of molecule----------------------
        cx = 0.d0
        cy = 0.d0
        cz = 0.d0
        do i = m0, m0 + nsp(m) - 1
            cx = xx(1, i) + cx
            cy = xx(2, i) + cy
            cz = xx(3, i) + cz
        end do
        cx = cx/nsp(m)
        cy = cy/nsp(m)
        cz = cz/nsp(m)
        !     ------Rotate---------------------------------------------
        x0 = cx*t11 + cy*t21 + cz*t31
        y0 = cx*t12 + cy*t22 + cz*t32
        z0 = cx*t13 + cy*t23 + cz*t33
        !     ------First map into cube of size 2*box(1)/sqrt(3) ------
        xt = anint(x0*tobinv)
        x = x0 - xt*tobox
        yt = anint(y0*tobinv)
        y = y0 - yt*tobox
        zt = anint(z0*tobinv)
        z = z0 - zt*tobox
        !     ------ wrap molecules external to diag faces -----
        xt = abs(x)
        yt = abs(y)
        zt = abs(z)
        if (xt + yt + zt > 3.d0*facecoord) then
            x = x - 2.d0*facecoord*sign(one, x)
            y = y - 2.d0*facecoord*sign(one, y)
            z = z - 2.d0*facecoord*sign(one, z)
        end if
        !     ------get the translation in the rotated space for this
        !     ------   molecules c-o-geom -----------------------------
        xt = x - x0
        yt = y - y0
        zt = z - z0
        !     ------Rotate---------------------------------------------
        cx = xt*t11 + yt*t12 + zt*t13
        cy = xt*t21 + yt*t22 + zt*t23
        cz = xt*t31 + yt*t32 + zt*t33
        !     ------Now move the molecule------------------------
        do i = m0, m0 + nsp(m) - 1
            xx(1, i) = xx(1, i) + cx
            xx(2, i) = xx(2, i) + cy
            xx(3, i) = xx(3, i) + cz
        end do
        m0 = m0 + nsp(m)
    end do  !  m=1,nspm
    return
end subroutine wrap_to

!---------------------------------------------------------------------

!      GET_POSITION

!    find the center of a set of atoms, and the extent

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_position here]
subroutine get_position(n, x, xc, yc, zc, e, verbose)
    implicit none

    !   Input:  n  number of atoms
    !           x  coordinates
    !   Output: xc,yc,zc  coordinates geometric center
    !           e   xmin,ymin,zmin,xmax,ymax,zmax

#  include "extra.h"
    integer n, i, j, verbose
    _REAL_ x(3, n), xc, yc, zc, e(3, 2)

    xc = 0.d0
    yc = 0.d0
    zc = 0.d0
    do i = 1, 3
        e(i, 1) = x(i, 1)
        e(i, 2) = x(i, 1)
    end do
    do i = 1, n
        do j = 1, 3
            e(j, 1) = min(e(j, 1), x(j, i))
            e(j, 2) = max(e(j, 2), x(j, i))
        end do
    end do
    xc = (e(1, 2) + e(1, 1))*0.5d0
    yc = (e(2, 2) + e(2, 1))*0.5d0
    zc = (e(3, 2) + e(3, 1))*0.5d0

    if (verbose > 0 .and. master) then
        write (6, *) "*********** Molecule dimensions ****************"
        write (6, '(a,3f8.3)') "   x min max center   ", e(1, 1), e(1, 2), xc
        write (6, '(a,3f8.3)') "   y min max center   ", e(2, 1), e(2, 2), yc
        write (6, '(a,3f8.3)') "   z min max center   ", e(3, 1), e(3, 2), zc
    end if
    return
end subroutine get_position

!---------------------------------------------------------------------

!      RE_POSITION

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine re_position here]
subroutine re_position(n, ntr, x, xr, &
    xc, yc, zc, x0, y0, z0, &
    vec, mv_flag, verbose)

    !    move the center of a set of atoms

    implicit none
    !   Input: n  number of atoms
    !          x  coordinates
    !          xr reference coordinates
    !          xc,yc,zc current center
    !          x0,y0,z0 new center

    !  Output: x,xr  moved coordinates

#  include "extra.h"
    integer n, ntr, i, j, verbose
    _REAL_ x(3, n), xr(3, n), vec(3), xc, yc, zc, x0, y0, z0
    _REAL_ xd, yd, zd
    logical mv_flag

    xd = x0 - xc
    yd = y0 - yc
    zd = z0 - zc

    if (master) &
        write (6, '(a,3f10.6)') "| RE_POSITION Moving by ", xd, yd, zd

    vec(1) = vec(1) + xd
    vec(2) = vec(2) + yd
    vec(3) = vec(3) + zd
    do i = 1, n
        x(1, i) = xd + x(1, i)
        x(2, i) = yd + x(2, i)
        x(3, i) = zd + x(3, i)
    end do
    if (ntr > 0) then
        do i = 1, n
            xr(1, i) = xd + xr(1, i)
            xr(2, i) = yd + xr(2, i)
            xr(3, i) = zd + xr(3, i)
        end do
    end if
    if (verbose > 0 .and. master) then
        write (6, *) "*********** Coords moved ****************"
        write (6, *) "delta x,y,z ", xd, yd, zd
    end if
    mv_flag = .true.
    return
end subroutine re_position

!==========================================================
!        NONPER_BOX

!       Resize the box for the size of the system

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nonper_box here]
subroutine nonper_box(extents, extraboxdim, a, b, c)
    implicit none

    _REAL_ extents(3, 2), a, b, c, extraboxdim

    a = extents(1, 2) - extents(1, 1) + 2*extraboxdim
    b = extents(2, 2) - extents(2, 1) + 2*extraboxdim
    c = extents(3, 2) - extents(3, 1) + 2*extraboxdim
    return
end subroutine nonper_box
!==============================================================================
!                FIRSTBOX

!       Make an initial determination of the box sizer
!         for nonperiodic systems so that
!         correct list grid can be made in ew_setup calls
!         before allocating memory for that grid.

!      Nonperiodic box is facilitated as follows:
!         Read_ewald determines if nonperiodic from ntb and igb
!                    calls this routine.
!             Firstbox determines box size, and whether it is
!                   a no-cutoff system.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine firstbox here]
subroutine firstbox()
    use nblist, only : cutoffnb, a, b, c
    use file_io_dat
    use binrestart, only : check_nc_restart, read_nc_restart
    implicit none
#  include "md.h"
#  include "ew_cntrl.h"
#  include "extra.h"
#  include "box.h"

#ifdef MPI
#  include "parallel.h"
#else
    integer mytaskid, numtasks
    parameter(mytaskid=0, numtasks=1)
#endif

    logical form, oldcrdtype

    _REAL_ xcen, ycen, zcen, extents(3, 2), extent, centertest
    ! For Netcdf restart, not efficient to just read 6 coords at a time.
    _REAL_, dimension(:), allocatable :: NCcoords

    integer i, j
    integer nread, nleft, natom
    _REAL_ x(6)

    integer ier, lun

    character(len=256) line_test

    ! -------------------------------------------------------------------------
    !   --------- Still to be improved:
    !     Restraints: when ntr=1, the extents of the restraints should
    !                 also be checked for moving the imagcrds and box size
    ! -------------------------------------------------------------------------
    !  ---- For now, we only need coordinates
    !  ----   Check if ntx agrees with a nonper box with formatted inpcrd file

    form = (ntx == 1 .or. ntx == 5 .or. ntx == 7) !.or. ntx == 9)
    if (.not. form) then
        write (6, '(a)') &
            " Cannot read unformatted coord for nonperiodic system"
        write (6, '(a)') &
            " For ntb=0, ntx must be 1 or 5"
        call mexit(6, 1)
    end if
    ! --------------------------------------------------------------------------
    ! ----- All is OK to read the coord file
    ! -----   We are only interested in the coordinates right now and will throw
    ! -----   them away upon return (when memory is dynamically allocated)
    ! -----   or overwrite them in getcor in sander().
    ! ----- extent(1-3,1) will have min X, Y, Z
    ! ----- extent(1-3,2) will have max X, Y, Z

    ! ---------- NetCDF Restart
    if (check_nc_restart(inpcrd)) then
        ! Allocate space to read in all coords
        allocate (NCcoords(natom*3), stat=i)
        REQUIRE(i == 0)
        ! Since calling this with ntx=1, velocities will not be written to
        ! so OK to pass NCcoords in twice. Also using xcen and ycen as dummy
        ! args for temperature and time.
        call read_nc_restart(inpcrd, line_test, 1, natom, NCcoords, NCcoords, xcen, ycen)
        ! Search for max and min X, Y, Z
        extents(1, 1) = NCcoords(1) ! MinX
        extents(2, 1) = NCcoords(2) ! MinY
        extents(3, 1) = NCcoords(3) ! MinZ
        extents(1, 2) = NCcoords(1) ! MaxX
        extents(2, 2) = NCcoords(2) ! MaxY
        extents(3, 2) = NCcoords(3) ! MaxZ
        j = 4
        do i = 2, natom
            if (NCcoords(j) .lt. extents(1, 1)) extents(1, 1) = NCcoords(j)
            if (NCcoords(j) .gt. extents(1, 2)) extents(1, 2) = NCcoords(j)
            j = j + 1
            if (NCcoords(j) .lt. extents(2, 1)) extents(2, 1) = NCcoords(j)
            if (NCcoords(j) .gt. extents(2, 2)) extents(2, 2) = NCcoords(j)
            j = j + 1
            if (NCcoords(j) .lt. extents(3, 1)) extents(3, 1) = NCcoords(j)
            if (NCcoords(j) .gt. extents(3, 2)) extents(3, 2) = NCcoords(j)
            j = j + 1
        end do
        ! Deallocate coordinates
        deallocate (NCcoords, stat=i)
        REQUIRE(i == 0)
        ! ---------- Formatted (ASCII) Restart
    else
        lun = 9
        call amopen(lun, inpcrd, 'O', 'F', 'R')
        read (lun, *) !skip the title

        read (lun, '(a80)') line_test
        if (line_test(6:6) == ' ') then ! this is an old, i5 file
            read (line_test, '(i5)') natom
        elseif (line_test(7:7) == ' ') then ! sander 7/8/9/10 large system format...
            read (line_test, '(i6)') natom
        elseif (line_test(8:8) == ' ') then ! Sander 11 - 1 mil+ format
            read (line_test, '(i7)') natom
        else                   ! assume amber 11 VERY large system format. 10 mil+
            read (line_test, '(i8)') natom
        end if

        nread = natom/2
        nleft = mod(natom, 1)

        read (lun, 9028) (x(i), i=1, 6)
        do i = 1, 3
            extents(i, 1) = x(i)
            extents(i, 2) = x(i)
        end do
        do j = 1, 3
            extents(j, 1) = min(extents(j, 1), x(j + 3))
            extents(j, 2) = max(extents(j, 2), x(j + 3))
        end do
        do i = 2, nread
            read (lun, 9028) (x(j), j=1, 6)
            do j = 1, 3
                extents(j, 1) = min(extents(j, 1), x(j))
                extents(j, 2) = max(extents(j, 2), x(j))
                extents(j, 1) = min(extents(j, 1), x(j + 3))
                extents(j, 2) = max(extents(j, 2), x(j + 3))
            end do
        end do
        if (nleft == 1) then
            read (lun, 9028) (x(i), i=1, 3)
            do j = 1, 3
                extents(j, 1) = min(extents(j, 1), x(j))
                extents(j, 2) = max(extents(j, 2), x(j))
            end do
        end if

9028    format(6f12.7)
        close (lun)
    end if
    ! --------------------------------------------------------------------------
    ! -----

    extent = max(extents(1, 2) - extents(1, 1), &
        extents(2, 2) - extents(2, 1), &
        extents(3, 2) - extents(3, 1))
    nocutoff = (extent < cutoffnb)
    if (master) then
        if (nocutoff) then
            write (6, '("|  ",a)') &
                "*** cutoff > system size, list only builds once"
        else
            if (verbose >= 1) write (6, '("|  ",a)') &
                "*** cutoff < system size, heuristic list update"
        end if
    end if
    call nonper_box(extents, extraboxdim, a, b, c)

    !   --- calculate offset to put stuff in the middle of the box
    !       when determining fractional coordinates. (map_coords in ew_force.f)

    xbox0 = a*0.5d0 - (extents(1, 2) + extents(1, 1))*0.5d0
    ybox0 = b*0.5d0 - (extents(2, 2) + extents(2, 1))*0.5d0
    zbox0 = c*0.5d0 - (extents(3, 2) + extents(3, 1))*0.5d0

    return
end subroutine firstbox
