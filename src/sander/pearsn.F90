#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pearsn here]
subroutine pearsn(x, y, n, r, prob, z, rms, xrfact, xrfac6, nuse)

    ! Subroutine PEARSoN linear correlation coefficient

    !   adapted from Numerical Recipes, p. 487
    !   ---gets linear correlation coefficient
    !      also calculates root-mean-square error of y vs. x,
    !      and xrfact == sum |y(i)-x(i)|/sum[x(i)],
    !      and xrfac6 == sum [y(i)**(-1/6) - |x(i)|**(-1/6)]/
    !                          sum[|x(i)|**(-1/6)]

    use constants, only : zero, one, half, INVSQRT2, SIXTH
    implicit none

    integer n
    _REAL_ x(n), y(n)
    _REAL_ r, prob, z, rms, xrfact, xrfac6
    integer nuse

    _REAL_ tiny
    parameter(tiny=1.d-20)

    _REAL_ ax, ay
    _REAL_ den
    _REAL_ erfcc
    integer j
    _REAL_ rnum, rden
    _REAL_ rnum6, rden6
    _REAL_ sxx, sxy, syy
    _REAL_ x6, y6
    _REAL_ xt, yt

    ax = zero
    ay = zero
    rnum = zero
    rden = zero
    rnum6 = zero
    rden6 = zero
    rms = zero
    nuse = 0
    xrfact = zero
    xrfac6 = zero
    r = zero
    if (n <= 2) return
    do j = 1, n
        nuse = nuse + 1
        ax = ax + x(j)
        ay = ay + y(j)
        rnum = rnum + abs(y(j) - x(j))
        rden = rden + x(j)
        rms = rms + (y(j) - x(j))**2
        if (abs(y(j)) > 1.0e-5) then
            y6 = abs(y(j))**(-SIXTH)
        else
            y6 = zero
        end if
        if (abs(x(j)) > tiny) then
            x6 = abs(x(j))**(-SIXTH)
        else
            x6 = zero
        end if
        rnum6 = rnum6 + abs(y6 - x6)
        rden6 = rden6 + x6
    end do
    ax = ax/nuse
    ay = ay/nuse
    if (rden == zero) then
        write (6, *) 'rden is bad in pearsn'
        return
    else
        xrfact = rnum/rden
    end if
    if (rden6 == zero) then
        write (6, *) 'rden6 is bad in pearsn'
        return
    else
        xrfac6 = rnum6/rden6
    end if
    rms = sqrt(rms/nuse)
    sxx = zero
    sxy = zero
    syy = zero
    do j = 1, n
        xt = x(j) - ax
        yt = y(j) - ay
        sxx = sxx + xt**2
        syy = syy + yt**2
        sxy = sxy + xt*yt
    end do
    den = sxx*syy
    if (den /= zero) then
        r = sxy/sqrt(sxx*syy)
    else
        write (6, *) 'Bad denominator in pearson'
        r = zero
        return
    end if
    z = half*log(((one + r) + tiny)/((one - r) + tiny))
    prob = erfcc(abs(r*sqrt(nuse - one))*INVSQRT2)
    return
end subroutine pearsn

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function erfcc here]
function erfcc(x)

    !  Subroutine ERror FunCtion, Complmentary:

    use constants, only : zero, one, two, half
    implicit none

    _REAL_ erfcc
    _REAL_ x

    _REAL_ t
    _REAL_ z

    z = abs(x)
    t = one/(one + half*z)
    erfcc = t*exp(-z*z - 1.26551223d0 + t*(1.00002368d0 + t*(.37409196d0 + &
        t*(.09678418d0 + t*(-.18628806d0 + t*(.27886807d0 + t*(-1.13520398d0 + &
        t*(1.48851587d0 + t*(-.82215223d0 + t*.17087277d0)))))))))
    if (x < zero) erfcc = two - erfcc
    return
end function erfcc
