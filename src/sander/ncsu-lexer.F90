#include "ncsu-utils.h"

module ncsu_lexer

    implicit none

    private

    character(len=*), public, parameter :: DIGITS = '0123456789'
    character(len=*), public, parameter :: LETTERS = &
        'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer, public, parameter :: TOKEN_ERROR_END_OF_LINE = -1
    integer, public, parameter :: TOKEN_ERROR_UNTERMINATED = -2
    integer, public, parameter :: TOKEN_ERROR_UNRECOGNIZED = -3
    integer, public, parameter :: TOKEN_ERROR_BAD_REAL = -4

    integer, public, parameter :: TOKEN_WORD = 1
    integer, public, parameter :: TOKEN_REAL = 2
    integer, public, parameter :: TOKEN_INTEGER = 3
    integer, public, parameter :: TOKEN_STRING = 4
    integer, public, parameter :: TOKEN_COMMA = 5
    integer, public, parameter :: TOKEN_EQUAL_SIGN = 6
    integer, public, parameter :: TOKEN_LEFT_PARENTHESIS = 7
    integer, public, parameter :: TOKEN_RIGHT_PARENTHESIS = 8

    public :: token

contains

!
! returns a typecode (TOKEN_*); 'first' is set to zero on error
!

    integer function token(str, first, last)

        use ncsu_utils

        implicit none

        character(len=*), intent(in) :: str
        integer, intent(inout) :: first
        integer, intent(out) :: last

        integer :: nchars

        nchars = len_trim(str)

        !
        ! find the first non-blank character
        !

        if (first < 1) &
            first = 1

        do while (first <= nchars)
            if (str(first:first) /= ' ') &
                exit
            first = first + 1
        end do

        if (first > nchars) then
            token = TOKEN_ERROR_END_OF_LINE
            return
        end if

        !
        ! handle comments
        !

        if (str(first:first) == '#' .or. str(first:first) == '!') then
            token = TOKEN_ERROR_END_OF_LINE
            return
        end if

        !
        ! find the last character of the token
        !

        last = first

        if (str(first:first) == ',') then
            token = TOKEN_COMMA
        else if (str(first:first) == '=') then
            token = TOKEN_EQUAL_SIGN
        else if (str(first:first) == '(') then
            token = TOKEN_LEFT_PARENTHESIS
        else if (str(first:first) == ')') then
            token = TOKEN_RIGHT_PARENTHESIS
        else if (str(first:first) == '''') then
            last = scan(str(first + 1:), '''')
            if (last == 0) then
                token = TOKEN_ERROR_UNTERMINATED
            else
                last = last + first
                token = TOKEN_STRING
                ncsu_assert(last <= nchars)
            end if
        else if (scan(LETTERS//'_', str(first:first)) /= 0) then
            do while (last < nchars)
                last = last + 1
                if (scan(LETTERS//DIGITS//'_', str(last:last)) == 0) &
                    exit
            end do
            if (scan(LETTERS//DIGITS//'_', str(last:last)) == 0) &
                last = last - 1
            token = TOKEN_WORD
            ncsu_assert(last <= nchars)
        else
            token = TOKEN_INTEGER
            if (last == nchars) then
                if (scan(DIGITS, str(last:last)) == 0) &
                    token = TOKEN_ERROR_UNRECOGNIZED
            else
                ncsu_assert(last < nchars)
                if (str(last:last) == '-' .or. str(last:last) == '+') &
                    last = last + 1

                if (scan(DIGITS, str(last:last)) == 0) then
                    token = TOKEN_ERROR_UNRECOGNIZED
                else
                    do while (last < nchars)
                        last = last + 1
                        if (scan(DIGITS, str(last:last)) == 0) then
                            if (str(last:last) /= '.') then
                                last = last - 1
                                exit
                            end if
                            if (token == TOKEN_REAL) then
                                token = TOKEN_ERROR_BAD_REAL
                                return
                            end if
                            last = last + 1
                            if (last > nchars) then
                                token = TOKEN_ERROR_BAD_REAL
                                return
                            end if
                            if (scan(DIGITS, str(last:last)) == 0) then
                                token = TOKEN_ERROR_BAD_REAL
                                return
                            end if
                            token = TOKEN_REAL
                        end if
                    end do
                end if
            end if
        end if

    end function token

end module ncsu_lexer
