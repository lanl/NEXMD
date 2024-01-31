#include "dprec.fh"

!---------------------------------------------------
!    copy triangular matrix to square matrix
!---------------------------------------------------
subroutine unpacking_test(n, ini, fin, su)
    implicit none
    integer :: n, i, j
    character*1 su
    _REAL_ ini(n*(n + 1)/2), fin(n, n)

    if (su .eq. 's') then
        do i = 1, n
            j = i*(i - 1)/2
            fin(i, 1:i) = ini(j:j + i)
            fin(1:i, i) = ini(j:j + i)
        end do
    elseif (su .eq. 'u') then
        do i = 1, n
            j = i*(i - 1)/2
            fin(i, i:n) = ini(j:j + i)
            fin(i:n, i) = -ini(j:j + i)
            fin(i, i) = 0.0
        end do
    else
        write (6, *) 'Unrecognized flag to unpacking'
    end if
    return

!---------------------------------------------------
!       This subrotine copy square matrix to triangular matrix
!---------------------------------------------------
  entry packing_test(n, fin, ini, su)

    if (su .eq. 's') then
        do i = 1, n
            do j = 1, i
                ini(i*(i - 1)/2 + j) = 0.5*(fin(j, i) + fin(i, j))
            end do
        end do
    elseif (su .eq. 'u') then
        do i = 1, n
            do j = 1, i
                ini(i*(i - 1)/2 + j) = 0.5*(fin(j, i) - fin(i, j))
            end do
            ini(i*(i - 1)/2 + i) = 0.0
        end do
    else
        write (6, *) 'Unrecognized flag to packing'
    end if

    return
end subroutine unpacking_test
