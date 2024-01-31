! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
!
!  Some commonly used routines
!
! Author: Taisung Lee
!         Andreas W. Goetz
! Date  : April 2011
!
module UtilitiesModule

    implicit none

    private
    public :: Cross
    public :: Upcase
    public :: print

    interface print
        module procedure print_integer_array
        module procedure print_real_array
        module procedure print_real_array_packed
    end interface print

contains

    subroutine Cross(v1, v2, v12)

        ! v12 is cross product of v1 and v2
        implicit none

        _REAL_, intent(in)  ::  v1(3), v2(3)
        _REAL_, intent(out) :: v12(3)

        v12(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v12(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v12(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end subroutine Cross

    function Upcase(string) result(upper)

        implicit none

        character(len=*), intent(in) :: string
        character(len=len(string)) :: upper
        integer :: i, ic

        do i = 1, len(string)
            ic = iachar(string(i:i))
            if (ic > 96 .and. ic < 123) then
                ic = ic - 32
            end if
            upper(i:i) = achar(ic)
        end do

    end function Upcase

    subroutine print_integer_array(name, array)

        implicit none

        character(len=*), intent(in) :: name
        integer, dimension(:), intent(in) :: array

        integer :: i, j, jstart, jend
        integer :: jstep = 10

        write (6, '(a)') trim(name)//':'
        jstart = 1
        do i = 1, size(array)/jstep + 1
            jend = min((jstart + jstep - 1), size(array))
            write (6, '(10(i8))') (array(j), j=jstart, jend)
            jstart = jstart + jstep
        end do

    end subroutine print_integer_array

    subroutine print_real_array(name, array)

        implicit none

        character(len=*), intent(in) :: name
        _REAL_, dimension(:), intent(in) :: array

        integer :: i, j, jstart, jend
        integer :: jstep = 8

        write (6, '(a)') trim(name)//':'
        jstart = 1
        do i = 1, size(array)/jstep + 1
            jend = min((jstart + jstep - 1), size(array))
            write (6, '(8(f9.4,x))') (array(j), j=jstart, jend)
            jstart = jstart + jstep
        end do

    end subroutine print_real_array

    subroutine print_real_array_packed(name, array, packed)

        implicit none

        character(len=*), intent(in) :: name
        _REAL_, dimension(:), intent(in) :: array
        logical, intent(in) :: packed

        integer :: i, j, jstart, jend, dim, ipos
        integer :: jstep = 8

        if (.not. packed) then
            call print_real_array(name, array)
        else

            write (6, '(a)') trim(name)//':'
            dim = size(array)
            ! determine size of array
            do i = 1, 1000
                j = i*(i + 1)/2
                if (j == dim) then
                    dim = i
                    exit
                end if
            end do

            jstart = 1

            do
                do i = jstart, dim
                    ipos = (i - 1)*i/2
                    jend = min((jstart + jstep - 1), i)
                    write (6, '(8(f9.4,x))') (array(ipos + j), j=jstart, jend)
                end do
                write (6, '(/)')
                call flush (6)
                jstart = jstart + jstep
                if (jstart > dim) exit
            end do

        end if

    end subroutine print_real_array_packed

end module UtilitiesModule
