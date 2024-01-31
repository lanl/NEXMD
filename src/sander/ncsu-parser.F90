#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! BNF grammar of what parse_cf() parses (more or less)
!
! letter := "a"|"b"|"c"| ... |"z"|"A"|"B"|"C"| ... |"Z"
!
! digit := "0"|"1"|"2"| ... |"9"
!
! word_start := letter | "_"
!
! word_body := letter | letter word_body | digit word_body | "_" word_body
!
! word := word_start | word_start word_body
!
! integer := "-" digits_list | "+" digits_list | digits_list
!
! digits_list := digit | digit digits_list
!
! real :=  digits_list "." digits_list
!    | "+" digits_list "." digits_list
!    | "-" digits_list "." digits_list
!
! string := "'" whatever "'"
!
! key := word
!
! value := integer | real | string | word | list
!
! list_items := value | value, list_items
!
! list := "(" list_items ")"
!
! key_value_pair := key "=" value
!
! key_value_pairs := key_value_pair | key_value_pair key_value_pairs
!
! section := word "end" word | word section_content "end" word
!
! sections := section | section sections
!
! section_content := key_value_pairs | key_value_pairs sections
!      | sections key_value_pairs | key_value_pairs sections key_value_pairs
!
! CONTROL_FILE := NULL | section
!
! FORTRAN restrictions apply:
!  * input lines are truncated to 95 characters (it's FORTRAN 95)
!  * neither keys nor scalar (!list) values can exceed one line (95 chars)
!

module ncsu_parser

    use ncsu_cftree, only : node_t

    implicit none

    private

    type, private                  :: stack_node_t
        type(node_t), pointer :: node
        type(stack_node_t), pointer :: prev
    end type stack_node_t

    type, private :: stack_t
        type(stack_node_t), pointer :: top => null()
    end type stack_t

    private :: stack_pop
    private :: stack_push
    private :: stack_empty

!
! returns pointer to the root-node (or null()); caller
! must call node_cleanup() and deallocate() on it after use
!

    public :: parse_cf

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

    logical function stack_empty(stack)

        implicit none

        type(stack_t), intent(in) :: stack

        stack_empty = .not. associated(stack%top)

    end function stack_empty

!-----------------------------------------------------------------------------

    subroutine stack_push(stack, nodeptr)

        use ncsu_utils

        implicit none

        type(stack_t), intent(inout) :: stack
        type(node_t), pointer :: nodeptr

        integer                     :: error
        type(stack_node_t), pointer :: new

        allocate (new, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        new%node => nodeptr
        new%prev => stack%top

        stack%top => new

    end subroutine stack_push

!-----------------------------------------------------------------------------

    subroutine stack_pop(stack)

        use ncsu_utils

        implicit none

        type(stack_t), intent(inout) :: stack

        type(stack_node_t), pointer :: pray

        ncsu_assert(associated(stack%top))

        pray => stack%top
        stack%top => pray%prev

        deallocate (pray)

    end subroutine stack_pop

!-----------------------------------------------------------------------------

    function parse_cf(filename, rootname, lun) result(root)

        use ncsu_utils
        use ncsu_lexer
        use ncsu_value
        use ncsu_cftree
        use ncsu_constants
        use ncsu_sander_proxy

        implicit none

        type(node_t), pointer :: root
        type(node_t), pointer :: curr

        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: rootname

        integer, optional, intent(in) :: lun

        integer :: ios, lineno
        integer :: tb, te, ttype
        integer :: error
        integer :: un

        character(len=STRING_LENGTH) :: line
        character(len=STRING_LENGTH) :: prev

        type(stack_t) :: sections

        ncsu_assert(len(rootname) <= STRING_LENGTH)

        root => null()

        if (present(lun)) then
            un = lun
            rewind (unit=un, iostat=ios)

            if (ios > 0) then
                write (unit=ERR_UNIT, fmt='(/a,a,i3/)') &
                    NCSU_ERROR, 'failed to rewind UNIT ''', un
                call bail_out()
            end if
        else
            un = 56 ! FORTRAN is _exceptionally_ powerful language indeed
            open (unit=un, file=filename, iostat=ios, &
                form='FORMATTED', action='READ', status='OLD')

            if (ios > 0) then
                write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') &
                    NCSU_ERROR, 'failed to open ''', trim(filename), ''' for reading'
                call bail_out()
            end if
        end if

        !
        ! scan for 'rootname' section
        !

        lineno = 0
        do
            call next_line()
            if (ios < 0) exit

            ncsu_assert(ios == 0)

            if (lineno > 1) then ! skip first line
                tb = 0
                ttype = token(line, tb, te)
                if (ttype == TOKEN_WORD) then
                    if (line(tb:te) == rootname) &
                        exit
                end if
            end if
        end do

        if (ios < 0) then
            return
        end if

        !
        ! there is 'rootname' section ... parse it out
        !

        ncsu_assert(ios == 0)
        ncsu_assert(stack_empty(sections))

        do
            ncsu_assert(ttype == TOKEN_WORD)
            prev = line(tb:te) ! prev is always of type WORD
            call next_token_not_eof()
            ncsu_assert(ttype > 0)

            select case (ttype)
              case (TOKEN_EQUAL_SIGN)
                call parse_key_value(key=prev)

              case (TOKEN_WORD)
                if (prev == 'end') then
                    ncsu_assert(associated(curr))
                    ncsu_assert(.not. stack_empty(sections))

                    prev = node_title(curr)
                    if (prev /= line(tb:te)) &
                        call parse_error('got ''end '//line(tb:te)// &
                        ''' instead of ''end '//trim(prev)//'''')

                    call stack_pop(sections)

                    if (stack_empty(sections)) then
                        curr => null()
                        exit
                    else
                        curr => sections%top%node
                    end if

                    call next_token_not_eof()
                    if (ttype /= TOKEN_WORD) &
                        call parse_error()
                else
                    if (.not. associated(root)) then
                        ncsu_assert(stack_empty(sections))
                        allocate (root, stat=error)
                        if (error /= 0) &
                            NCSU_OUT_OF_MEMORY
                        curr => root
                    else
                        ncsu_assert(.not. stack_empty(sections))
                        curr => node_create_child(curr)
                    end if

                    call node_set_title(curr, prev)
                    call stack_push(sections, curr)
                end if

              case default
                call parse_error()
            end select
        end do

        ncsu_assert(stack_empty(sections))

        if (.not. present(lun)) &
            close (unit=un)

    contains

!-----------------------------------------------------------------------------

        subroutine bail_out()

            implicit none

            do while (.not. stack_empty(sections))
                call stack_pop(sections)
            end do

            if (associated(root)) then
                call node_cleanup(root)
                deallocate (root)
            end if

            call terminate()

        end subroutine bail_out

!-----------------------------------------------------------------------------

        subroutine print_errloc()

            implicit none

            write (unit=ERR_UNIT, &
                fmt='(/a,a,a,a,'//pfmt(lineno)//',a,'//pfmt(tb)//',a)', advance='NO') &
                NCSU_ERROR, '''', trim(filename), ''', line:', lineno, &
                ', col:', tb, ' : '

        end subroutine print_errloc

!-----------------------------------------------------------------------------

        subroutine premature_eof(msg)

            implicit none

            character(len=*), intent(in), optional :: msg

            write (unit=ERR_UNIT, fmt='(/a,a,a,a,'//pfmt(lineno)//')', advance='NO') &
                NCSU_ERROR, '''', trim(filename), &
                ''' : premature end of file after line ', lineno

            if (present(msg)) &
                write (unit=ERR_UNIT, fmt='(a,a)', advance='NO') ' - ', trim(msg)

            write (unit=ERR_UNIT, fmt='(a/)') ''

            call bail_out()

        end subroutine premature_eof

!-----------------------------------------------------------------------------

        subroutine parse_error(msg)

            implicit none

            character(len=*), optional, intent(in) :: msg

            call print_errloc()

            write (unit=ERR_UNIT, fmt='(a)', advance='NO') 'parse error'
            if (present(msg)) &
                write (unit=ERR_UNIT, fmt='(a,a)', advance='NO') ' - ', msg

            write (unit=ERR_UNIT, fmt='(a/)') ''
            call bail_out()

        end subroutine parse_error

!-----------------------------------------------------------------------------

        subroutine next_line()

            implicit none

            read (unit=un, fmt='(a)', iostat=ios) line

            if (ios == 0) &
                lineno = lineno + 1

            if (ios > 0) then
                write (unit=ERR_UNIT, fmt='(/a,a,a,a,'//pfmt(lineno)//'/)') &
                    NCSU_ERROR, 'read() from ''', &
                    trim(filename), ''' failed after line ', lineno
                call bail_out()
            end if

        end subroutine next_line

!-----------------------------------------------------------------------------

        subroutine next_token()

            implicit none

            tb = te + 1
            do
                ttype = token(line, tb, te)
                if (ttype == TOKEN_ERROR_END_OF_LINE) then
                    tb = 0
                    call next_line()
                    if (ios < 0) exit
                else
                    exit
                end if
            end do

            if (ttype < 0 .and. ios >= 0) then
                call print_errloc()
                select case (ttype)
                  case (TOKEN_ERROR_UNTERMINATED)
                    write (unit=ERR_UNIT, fmt='(a/)') &
                        'unterminated character constant'
                  case (TOKEN_ERROR_UNRECOGNIZED)
                    write (unit=ERR_UNIT, fmt='(a/)') 'unexpected character'
                  case (TOKEN_ERROR_BAD_REAL)
                    write (unit=ERR_UNIT, fmt='(a/)') 'ill-formatted real number'
                  case default
                    write (unit=ERR_UNIT, fmt='(a,i4,a/)') &
                        'must not be here ('//__FILE__//':', __LINE__, ')'
                end select
                call bail_out()
            end if

        end subroutine next_token

!-----------------------------------------------------------------------------

        subroutine next_token_not_eof(msg)

            implicit none

            character(len=*), intent(in), optional :: msg

            call next_token()

            if (ios < 0) &
                call premature_eof(msg)

        end subroutine next_token_not_eof

!-----------------------------------------------------------------------------

        subroutine parse_key_value(key)

            implicit none

            character(len=*), intent(in) :: key

            type(value_t) :: list

            logical :: original

            ncsu_assert(associated(curr))

            call next_token_not_eof()
            ncsu_assert(ttype > 0)

            if (ttype == TOKEN_WORD) then
                original = node_insert(curr, key, line(tb:te))
            else if (ttype == TOKEN_REAL) then
                original = node_insert(curr, key, token_to_real())
            else if (ttype == TOKEN_INTEGER) then
                original = node_insert(curr, key, token_to_integer())
            else if (ttype == TOKEN_STRING) then
                original = node_insert(curr, key, line(tb + 1:te - 1))
            else if (ttype == TOKEN_LEFT_PARENTHESIS) then
                call parse_list(list)
                original = node_insert(curr, key, list)
                call value_cleanup(list)
            else
                ncsu_assert_not_reached()
                original = .false. ! for g95
            end if

            if (.not. original) &
                write (unit=ERR_UNIT, &
                fmt='(/a,a,a,a,'//pfmt(lineno)//',a,'//pfmt(tb)//',a,a,a/)') &
                NCSU_WARNING, '''', trim(filename), &
                ''', line:', lineno, ', col:', tb, &
                ' : key ''', trim(key), ''' has been already seen'

            call next_token_not_eof('unterminated section(s)')
            ncsu_assert(ttype > 0)

            if (ttype /= TOKEN_WORD) &
                call parse_error()

        end subroutine parse_key_value

!-----------------------------------------------------------------------------

        integer function token_to_integer()
            implicit none
            read (unit=line(tb:te), fmt=*) token_to_integer
        end function token_to_integer

!-----------------------------------------------------------------------------

        NCSU_REAL function token_to_real()
        implicit none
        read (unit=line(tb:te), fmt=*) token_to_real
    end function parse_cf

!-----------------------------------------------------------------------------

    recursive subroutine parse_list(accum) ! recursion could be avoided here

        implicit none

        type(value_t), intent(inout) :: accum

        type(value_t) :: list

        do
            call next_token_not_eof('unterminated list')
            ncsu_assert(ttype > 0)

            if (ttype == TOKEN_WORD) then
                call value_append(accum, line(tb:te))
            else if (ttype == TOKEN_REAL) then
                call value_append(accum, token_to_real())
            else if (ttype == TOKEN_INTEGER) then
                call value_append(accum, token_to_integer())
            else if (ttype == TOKEN_STRING) then
                call value_append(accum, line(tb + 1:te - 1))
            else if (ttype == TOKEN_COMMA) then
                call parse_error('unexpected comma')
            else if (ttype == TOKEN_EQUAL_SIGN) then
                call parse_error('unexpected equal sign')
            else if (ttype == TOKEN_LEFT_PARENTHESIS) then
                call parse_list(list)
                call value_append(accum, list)
                call value_cleanup(list)
            else if (ttype == TOKEN_RIGHT_PARENTHESIS) then
                call parse_error('unexpected right parenthesis')
            else
                ncsu_assert_not_reached()
                continue
            end if

            call next_token_not_eof('unterminated list')
            ncsu_assert(ttype > 0)

            if (ttype == TOKEN_COMMA) then
                continue
            else if (ttype == TOKEN_RIGHT_PARENTHESIS) then
                return
            else
                call parse_error('expected either comma or right parenthesis&
                &, got '''//line(tb:te)//''' instead')
            end if
        end do

    end subroutine parse_list

!-----------------------------------------------------------------------------

end module ncsu_parser

!-----------------------------------------------------------------------------

end module ncsu_parser
