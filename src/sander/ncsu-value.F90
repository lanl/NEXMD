#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! sort of C/C++ union
!

module ncsu_value

    use ncsu_constants, only : STRING_LENGTH

    implicit none

    private

    type, public :: value_t
        private
        integer, pointer :: iptr => null()
        NCSU_REAL, pointer :: rptr => null()
        character(len=STRING_LENGTH), pointer :: sptr => null()
        type(value_node_t), pointer :: head => null()
        type(value_node_t), pointer :: tail => null() ! %next is null()
    end type value_t

    private :: append_node

    type, public :: value_node_t
        type(value_t)               :: value
        type(value_node_t), pointer :: next => null()
    end type value_node_t

    private :: allocate_node

    interface assignment(=)
        module procedure assign_V_to_V, assign_I_to_V, assign_R_to_V, assign_S_to_V
    end interface assignment

    private :: assign_V_to_V, assign_I_to_V, assign_R_to_V, assign_S_to_V
    private :: assign_L_to_V

    public :: assignment(=)

!
! value must be cleaned when no longer needed (memory will be leaked otherwise)
!

    public :: value_cleanup

    public :: value_is_empty

    public :: value_is_real
    public :: value_is_list
    public :: value_is_string
    public :: value_is_integer

    public :: value_get_real
    public :: value_get_list
    public :: value_get_string
    public :: value_get_integer

    public :: value_print

    interface value_append
        module procedure append_I_to_V, append_R_to_V, append_S_to_V, append_V_to_V
    end interface value_append

    private :: append_I_to_V, append_R_to_V, append_S_to_V, append_V_to_V

    public :: value_append

#ifndef NCSU_DISABLE_ASSERT
    private :: value_ok
#endif /* NCSU_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

    function allocate_node() result(node)

        use ncsu_utils

        implicit none

        type(value_node_t), pointer :: node

        integer :: error

        allocate (node, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        node%next => null()

    end function allocate_node

!-----------------------------------------------------------------------------

    subroutine append_node(value, node)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: value
        type(value_node_t), pointer :: node

        ncsu_assert(associated(node))

        !
        ! if it's a scalar value, convert it to the list
        !

        if (associated(value%iptr) .or. &
            associated(value%rptr) .or. &
            associated(value%sptr)) then

            ncsu_assert(.not. associated(value%tail))
            ncsu_assert(.not. associated(value%head))

            value%head => allocate_node()
            value%tail => value%head

            value%head%value%iptr => value%iptr
            value%head%value%rptr => value%rptr
            value%head%value%sptr => value%sptr
            value%head%next => null()

            nullify (value%iptr, value%rptr, value%sptr)

        end if

        !
        ! actually append
        !

        if (associated(value%tail)) then

            ncsu_assert(associated(value%head))
            ncsu_assert(.not. associated(value%iptr))
            ncsu_assert(.not. associated(value%rptr))
            ncsu_assert(.not. associated(value%sptr))

            node%next => null()
            value%tail%next => node
            value%tail => node

        else

            ncsu_assert(value_is_empty(value))

            value%head => node
            value%tail => node

        end if

    end subroutine append_node

!-----------------------------------------------------------------------------

    recursive subroutine assign_V_to_V(var, expr)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: var
        type(value_t), intent(in)    :: expr

        logical :: self_assignment

        !
        ! check for self-assignment
        !

        self_assignment = .false.

        if (associated(var%iptr) .and. associated(expr%iptr)) then
            if (associated(var%iptr, expr%iptr)) &
                self_assignment = .true.
        end if

        if (associated(var%rptr) .and. associated(expr%rptr)) then
            if (associated(var%rptr, expr%rptr)) &
                self_assignment = .true.
        end if

        if (associated(var%sptr) .and. associated(expr%sptr)) then
            if (associated(var%sptr, expr%sptr)) &
                self_assignment = .true.
        end if

        if (associated(var%tail) .and. associated(expr%tail)) then
            if (associated(var%tail, expr%tail)) then
                ncsu_assert(associated(var%head, target=expr%head))
                self_assignment = .true.
            end if
        end if

        if (self_assignment) &
            return

        !
        ! not a self-assignment; dispatch appropriately
        !

        if (associated(expr%iptr)) then
            call assign_I_to_V(var, expr%iptr)
        else if (associated(expr%rptr)) then
            call assign_R_to_V(var, expr%rptr)
        else if (associated(expr%sptr)) then
            call assign_S_to_V(var, expr%sptr)
        else if (associated(expr%head)) then
            call assign_L_to_V(var, expr%head)
        else
            ncsu_assert(value_is_empty(expr))
            call value_cleanup(var)
        end if

    end subroutine assign_V_to_V

!-----------------------------------------------------------------------------

    subroutine assign_I_to_V(var, expr)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: var
        integer, intent(in)    :: expr

        integer :: error

        call value_cleanup(var)

        allocate (var%iptr, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        var%iptr = expr

    end subroutine assign_I_to_V

!-----------------------------------------------------------------------------

    subroutine assign_R_to_V(var, expr)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: var
        NCSU_REAL, intent(in)    :: expr

        integer :: error

        call value_cleanup(var)

        allocate (var%rptr, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        var%rptr = expr

    end subroutine assign_R_to_V

!-----------------------------------------------------------------------------

    subroutine assign_S_to_V(var, expr)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: var
        character(len=*), intent(in)    :: expr

        integer :: error, b, e

        call value_cleanup(var)

        allocate (var%sptr, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        e = len_trim(expr)
        b = verify(expr, ' ')

        ncsu_assert(e - b < STRING_LENGTH)

        var%sptr = expr(b:e)

    end subroutine assign_S_to_V

!----------------------------------------------------------------------------

    recursive subroutine assign_L_to_V(var, head)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: var
        type(value_node_t), pointer :: head

        type(value_node_t), pointer :: curr

        call value_cleanup(var)

        curr => head
        do while (associated(curr))
            call append_V_to_V(var, curr%value)
            curr => curr%next
        end do

    end subroutine assign_L_to_V

!----------------------------------------------------------------------------

    recursive subroutine append_V_to_V(value, other)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: value
        type(value_t), intent(in)    :: other

        type(value_node_t), pointer :: node

        ncsu_assert(value_ok(value))
        ncsu_assert(value_ok(other))

        node => allocate_node()
        call assign_V_to_V(node%value, other)
        call append_node(value, node)

    end subroutine append_V_to_V

!----------------------------------------------------------------------------

    subroutine append_I_to_V(value, x)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: value
        integer, intent(in)    :: x

        type(value_node_t), pointer :: node

        ncsu_assert(value_ok(value))

        node => allocate_node()
        call assign_I_to_V(node%value, x)
        call append_node(value, node)

    end subroutine append_I_to_V

!----------------------------------------------------------------------------

    subroutine append_R_to_V(value, x)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: value
        NCSU_REAL, intent(in)    :: x

        type(value_node_t), pointer :: node

        ncsu_assert(value_ok(value))

        node => allocate_node()
        call assign_R_to_V(node%value, x)
        call append_node(value, node)

    end subroutine append_R_to_V

!----------------------------------------------------------------------------

    subroutine append_S_to_V(value, x)

        use ncsu_utils

        implicit none

        type(value_t), intent(inout) :: value
        character(len=*), intent(in)    :: x

        type(value_node_t), pointer :: node

        ncsu_assert(value_ok(value))

        node => allocate_node()
        call assign_S_to_V(node%value, x)
        call append_node(value, node)

    end subroutine append_S_to_V

!----------------------------------------------------------------------------

    integer function value_get_integer(value)

        use ncsu_utils

        implicit none

        type(value_t), intent(in) :: value

        ncsu_assert(value_ok(value))
        ncsu_assert(associated(value%iptr))
        value_get_integer = value%iptr

    end function value_get_integer

!----------------------------------------------------------------------------

    NCSU_REAL function value_get_real(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    ncsu_assert(associated(value%rptr))
    value_get_real = value%rptr

end module ncsu_value

!----------------------------------------------------------------------------

function value_get_string(value) result(string)

    use ncsu_utils

    implicit none

    character(len=STRING_LENGTH) :: string
    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    ncsu_assert(associated(value%sptr))
    string = value%sptr

end function value_get_string

!----------------------------------------------------------------------------

function value_get_list(value) result(head)

    use ncsu_utils

    implicit none

    type(value_node_t), pointer :: head
    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    ncsu_assert(associated(value%head))
    head => value%head

end function value_get_list

!----------------------------------------------------------------------------

recursive subroutine value_cleanup(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(inout) :: value

    type(value_node_t), pointer :: pray

    ncsu_assert(value_ok(value))

    if (associated(value%iptr)) then
        deallocate (value%iptr)
    else if (associated(value%rptr)) then
        deallocate (value%rptr)
    else if (associated(value%sptr)) then
        deallocate (value%sptr)
    else if (associated(value%head)) then
        do while (associated(value%head))
            pray => value%head
            value%head => pray%next
            call value_cleanup(pray%value)
            deallocate (pray)
        end do
        value%tail => null()
    end if

end subroutine value_cleanup

!----------------------------------------------------------------------------

logical function value_is_empty(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))

    value_is_empty = ((.not. associated(value%iptr)) &
        .and. (.not. associated(value%rptr)) &
        .and. (.not. associated(value%sptr)) &
        .and. (.not. associated(value%head)))

end function value_is_empty

!----------------------------------------------------------------------------

#ifndef NCSU_DISABLE_ASSERT
logical function value_ok(value)

    implicit none

    type(value_t), intent(in) :: value

    if (associated(value%iptr)) then
        value_ok = ((.not. associated(value%rptr)) &
            .and. (.not. associated(value%sptr)) &
            .and. (.not. associated(value%head)) &
            .and. (.not. associated(value%tail)))
    else if (associated(value%rptr)) then
        value_ok = ((.not. associated(value%iptr)) &
            .and. (.not. associated(value%sptr)) &
            .and. (.not. associated(value%head)) &
            .and. (.not. associated(value%tail)))
    else if (associated(value%sptr)) then
        value_ok = ((.not. associated(value%iptr)) &
            .and. (.not. associated(value%rptr)) &
            .and. (.not. associated(value%head)) &
            .and. (.not. associated(value%tail)))
    else if (associated(value%tail)) then
        value_ok = ((.not. associated(value%iptr)) &
            .and. (.not. associated(value%rptr)) &
            .and. (.not. associated(value%sptr)) &
            .and. (associated(value%head)))
    else if (associated(value%head)) then
        value_ok = ((.not. associated(value%iptr)) &
            .and. (.not. associated(value%rptr)) &
            .and. (.not. associated(value%sptr)) &
            .and. (associated(value%tail)))
    else
        value_ok = .true.
    end if

end function value_ok
#endif /* NCSU_DISABLE_ASSERT */

!----------------------------------------------------------------------------

logical function value_is_integer(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    value_is_integer = associated(value%iptr)

end function value_is_integer

!----------------------------------------------------------------------------

logical function value_is_real(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    value_is_real = associated(value%rptr)

end function value_is_real

!----------------------------------------------------------------------------

logical function value_is_string(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    value_is_string = associated(value%sptr)

end function value_is_string

!----------------------------------------------------------------------------

logical function value_is_list(value)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value

    ncsu_assert(value_ok(value))
    value_is_list = associated(value%head)

end function value_is_list

!----------------------------------------------------------------------------

subroutine value_print(value, lun)

    use ncsu_utils

    implicit none

    type(value_t), intent(in) :: value
    integer, intent(in) :: lun

    ncsu_assert(value_ok(value))

    call dispatch(value)

contains

    recursive subroutine dispatch(val)

        implicit none

        type(value_t), intent(in) :: val

        type(value_node_t), pointer :: curr

        if (associated(val%iptr)) then
            write (unit=lun, fmt='('//pfmt(val%iptr)//')', advance='NO') val%iptr
        else if (associated(val%rptr)) then
            write (unit=lun, fmt='('//pfmt(val%rptr, 4)//')', advance='NO') val%rptr
        else if (associated(val%sptr)) then
            write (unit=lun, fmt='(a)', advance='NO') trim(val%sptr)
        else if (associated(val%tail)) then
            curr => val%head
            write (unit=lun, fmt='(a)', advance='NO') '('
            do while (associated(curr%next))
                call dispatch(curr%value)
                write (unit=lun, fmt='(a)', advance='NO') ', '
                curr => curr%next
            end do
            call dispatch(curr%value)
            write (unit=lun, fmt='(a)', advance='NO') ')'
        else
            write (unit=lun, fmt='(a)', advance='NO') '<EMPTY>'
        end if

    end subroutine dispatch

end subroutine value_print

end module ncsu_value
