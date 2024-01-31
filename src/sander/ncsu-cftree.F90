#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! data structure to represent the control-file
!    [looks rather quite ugly in FORTRAN]
!

module ncsu_cftree

    use ncsu_value, only : value_t
    use ncsu_constants, only : STRING_LENGTH

    implicit none

    private

    type, private :: bucket_node_t
        character(len=STRING_LENGTH) :: key
        type(value_t)                  :: value
        type(bucket_node_t), pointer   :: next
    end type bucket_node_t

    type, private :: bucket_t
        type(bucket_node_t), pointer :: tail => null()
    end type bucket_t

    private :: bucket_cleanup

    type, public :: node_t
        private
        character(len=STRING_LENGTH) :: title
        type(bucket_t)                 :: buckets(11)
        type(child_t), pointer         :: head => null()
        type(child_t), pointer         :: tail => null()
    end type node_t

    type, public :: child_t
        type(node_t) node
        type(child_t), pointer :: next => null()
    end type child_t

    public :: node_cleanup

    interface node_insert
        module procedure insert_V, insert_I, insert_R, insert_S
    end interface node_insert

    private :: insert_V, insert_I, insert_R, insert_S, insert_bucket_node

    public :: node_title
    public :: node_set_title

    public :: node_insert
    public :: node_lookup

    public :: node_keys
    public :: node_children

    public :: node_create_child

#ifdef NCSU_ENABLE_NODE_PRINT
    public :: node_print
#endif /* NCSU_ENABLE_NODE_PRINT */

    private :: hash
    private :: find_bucket_node

    public :: node_lookup_list
    public :: node_lookup_real
    public :: node_lookup_string
    public :: node_lookup_logical
    public :: node_lookup_integer

    public :: node_lookup_positive_real
    public :: node_lookup_positive_integer

    private :: type_mismatch

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

    pure integer function hash(string)

        implicit none

        character(len=*), intent(in) :: string

        integer :: n, nchars

        nchars = len_trim(string)
        hash = 5381

        do n = 1, nchars
            hash = 33*hash + ichar(string(n:n))
        end do

        hash = abs(hash)

    end function hash

!-----------------------------------------------------------------------------

    subroutine bucket_cleanup(b)

        use ncsu_value

        implicit none

        type(bucket_t), intent(inout) :: b
        type(bucket_node_t), pointer  :: pray

        do while (associated(b%tail))
            pray => b%tail
            b%tail => pray%next
            call value_cleanup(pray%value)
            deallocate (pray)
        end do

    end subroutine bucket_cleanup

!-----------------------------------------------------------------------------

    recursive subroutine node_cleanup(node)

        implicit none

        type(node_t), intent(inout) :: node

        integer                :: n
        type(child_t), pointer :: child

        do n = 1, size(node%buckets)
            call bucket_cleanup(node%buckets(n))
        end do

        do while (associated(node%head))
            child => node%head
            node%head => child%next
            call node_cleanup(child%node)
            deallocate (child)
        end do

        node%tail => null()

    end subroutine node_cleanup

!-----------------------------------------------------------------------------

    function node_title(node) result(title)

        implicit none

        type(node_t), intent(in) :: node
        character(len=STRING_LENGTH)  :: title

        title = node%title

    end function node_title

!-----------------------------------------------------------------------------

    subroutine node_set_title(node, title)

        use ncsu_utils

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: title

        ncsu_assert(len(title) <= STRING_LENGTH)

        node%title = title

    end subroutine node_set_title

!-----------------------------------------------------------------------------

    logical function insert_bucket_node(node, key, bn)

        use ncsu_utils

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: key
        type(bucket_node_t), pointer      :: bn

        integer :: error, idx

        ncsu_assert(len(key) <= STRING_LENGTH)

        bn => find_bucket_node(node, key)

        if (.not. associated(bn)) then

            allocate (bn, stat=error)
            if (error /= 0) &
                NCSU_OUT_OF_MEMORY

            idx = 1 + mod(hash(key), size(node%buckets))
            bn%next => node%buckets(idx)%tail
            node%buckets(idx)%tail => bn

            bn%key = key
            insert_bucket_node = .true.

        else

            ncsu_assert(bn%key == key)
            insert_bucket_node = .false.

        end if

    end function insert_bucket_node

!-----------------------------------------------------------------------------

    logical function insert_I(node, key, ival)

        use ncsu_utils
        use ncsu_value

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: key
        integer, intent(in)    :: ival

        type(bucket_node_t), pointer :: bn

        insert_I = insert_bucket_node(node, key, bn)

        ncsu_assert(associated(bn))
        ncsu_assert(bn%key == key)

        bn%value = ival

    end function insert_I

!-----------------------------------------------------------------------------

    logical function insert_R(node, key, rval)

        use ncsu_utils
        use ncsu_value

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: key
        NCSU_REAL, intent(in)    :: rval

        type(bucket_node_t), pointer :: bn

        insert_R = insert_bucket_node(node, key, bn)

        ncsu_assert(associated(bn))
        ncsu_assert(bn%key == key)

        bn%value = rval

    end function insert_R

!-----------------------------------------------------------------------------

    logical function insert_S(node, key, sval)

        use ncsu_utils
        use ncsu_value

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: sval

        type(bucket_node_t), pointer :: bn

        insert_S = insert_bucket_node(node, key, bn)

        ncsu_assert(associated(bn))
        ncsu_assert(bn%key == key)

        bn%value = sval

    end function insert_S

!-----------------------------------------------------------------------------

    logical function insert_V(node, key, val)

        use ncsu_utils
        use ncsu_value

        implicit none

        type(node_t), intent(inout) :: node
        character(len=*), intent(in)    :: key
        type(value_t), intent(in)    :: val

        type(bucket_node_t), pointer :: bn

        insert_V = insert_bucket_node(node, key, bn)

        ncsu_assert(associated(bn))
        ncsu_assert(bn%key == key)

        bn%value = val

    end function insert_V

!-----------------------------------------------------------------------------

    function node_lookup(node, key) result(vptr)

        use ncsu_utils

        implicit none

        type(value_t), pointer :: vptr

        type(node_t), intent(in) :: node
        character(len=*), intent(in) :: key

        type(bucket_node_t), pointer :: bn

        ncsu_assert(len(key) <= STRING_LENGTH)

        bn => find_bucket_node(node, key)

        if (associated(bn)) then
            vptr => bn%value
        else
            vptr => null()
        end if

    end function node_lookup

!-----------------------------------------------------------------------------

!
! returns newly allocated array of strings; caller must deallocate()
!

    function node_keys(node) result(keys)

        use ncsu_utils

        implicit none

        type(node_t), intent(in) :: node
        character(len=STRING_LENGTH), pointer :: keys(:)

        integer :: error, n, nkeys
        type(bucket_node_t), pointer :: bn

        ! count the keys

        nkeys = 0
        do n = 1, size(node%buckets)
            bn => node%buckets(n)%tail
            do while (associated(bn))
                nkeys = nkeys + 1
                bn => bn%next
            end do
        end do

        allocate (keys(nkeys), stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        ! populate array

        do n = 1, size(node%buckets)
            bn => node%buckets(n)%tail
            do while (associated(bn))
                ncsu_assert(nkeys > 0)
                keys(nkeys) = bn%key
                nkeys = nkeys - 1
                bn => bn%next
            end do
        end do

    end function node_keys

!-----------------------------------------------------------------------------

    function node_children(node) result(children)

        implicit none

        type(child_t), pointer :: children
        type(node_t), intent(in) :: node

        children => node%head

    end function node_children

!-----------------------------------------------------------------------------

    function node_create_child(node) result(child)

        use ncsu_utils

        implicit none

        type(node_t), pointer :: child
        type(node_t), intent(inout) :: node

        type(child_t), pointer :: new

        integer :: error

        allocate (new, stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        if (.not. associated(node%head)) then

            ncsu_assert(.not. associated(node%tail))

            node%head => new
            node%tail => new

        else

            ncsu_assert(associated(node%tail))

            node%tail%next => new
            node%tail => new

        end if

        child => new%node

    end function node_create_child

!-----------------------------------------------------------------------------

    function find_bucket_node(node, key) result(bn)

        use ncsu_utils

        implicit none

        type(bucket_node_t), pointer :: bn
        type(node_t), intent(in) :: node
        character(len=*), intent(in) :: key

        integer :: idx

        idx = 1 + mod(hash(key), size(node%buckets))

        ncsu_assert(idx > 0)
        ncsu_assert(idx <= size(node%buckets))

        bn => node%buckets(idx)%tail

        do while (associated(bn))
            if (bn%key == key) &
                exit
            bn => bn%next
        end do

    end function find_bucket_node

!-----------------------------------------------------------------------------

#ifdef NCSU_ENABLE_NODE_PRINT
    subroutine node_print(node, lun)

        use ncsu_value

        implicit none

        type(node_t), intent(in) :: node
        integer, intent(in) :: lun

        call recursive_print(node, 0)

    contains

        recursive subroutine recursive_print(anode, depth)

            implicit none

            type(node_t), intent(in) :: anode
            integer, intent(in) :: depth

            integer                      :: n
            type(bucket_node_t), pointer :: bn
            type(child_t), pointer :: child

            call indent(3*depth)
            write (unit=lun, fmt='(a)') trim(anode%title)

            do n = 1, size(anode%buckets)
                bn => anode%buckets(n)%tail
                do while (associated(bn))
                    call indent(3*(depth + 1))
                    write (unit=lun, fmt='(a,a)', advance='NO') trim(bn%key), ' = '
                    call value_print(bn%value, lun)
                    write (unit=lun, fmt='(a)') ''
                    bn => bn%next
                end do
            end do

            child => node_children(anode)
            do while (associated(child))
                call recursive_print(child%node, depth + 1)
                child => child%next
            end do

            call indent(3*depth)
            write (unit=lun, fmt='(a,a)') 'end ', trim(anode%title)

        end subroutine recursive_print

        subroutine indent(n)

            implicit none

            integer, intent(in) :: n

            integer :: i

            do i = 1, n
                write (unit=lun, fmt='(a)', advance='NO') ' '
            end do

        end subroutine indent

    end subroutine node_print
#endif /* NCSU_ENABLE_NODE_PRINT */

!-----------------------------------------------------------------------------

    subroutine type_mismatch(key, value, expected)

        use ncsu_utils
        use ncsu_value
        use ncsu_constants
        use ncsu_sander_proxy

        implicit none

        character(len=*), intent(in) :: key
        type(value_t), pointer :: value
        character(len=*), intent(in) :: expected

        ncsu_assert(associated(value))

        write (unit=ERR_UNIT, fmt='(/a,a,a,a,a,a)', advance='NO') &
            NCSU_ERROR, 'expected ', expected, ' value for key ''', &
            key, ''', got '''

        call value_print(value, ERR_UNIT)
        write (unit=ERR_UNIT, fmt='(a/)') ''' instead'

        call terminate()

    end subroutine type_mismatch

!-----------------------------------------------------------------------------

    function node_lookup_integer(node, key, val) result(ok)

        use ncsu_value

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        integer, intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (.not. value_is_integer(value)) &
                call type_mismatch(key, value, 'integer')
            val = value_get_integer(value)
        end if

    end function node_lookup_integer

!-----------------------------------------------------------------------------

    function node_lookup_real(node, key, val) result(ok)

        use ncsu_value

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        NCSU_REAL, intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (.not. value_is_real(value)) &
                call type_mismatch(key, value, 'real number')
            val = value_get_real(value)
        end if

    end function node_lookup_real

!-----------------------------------------------------------------------------

    function node_lookup_positive_real(node, key, val) result(ok)

        use ncsu_value
        use ncsu_constants

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        NCSU_REAL, intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (value_is_real(value)) then
                val = value_get_real(value)
                if (val > zero) &
                    return
            end if
            call type_mismatch(key, value, 'positive real number')
        end if

    end function node_lookup_positive_real

!-----------------------------------------------------------------------------

    function node_lookup_positive_integer(node, key, val) result(ok)

        use ncsu_value

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        integer, intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (value_is_integer(value)) then
                val = value_get_integer(value)
                if (val > 0) &
                    return
            end if
            call type_mismatch(key, value, 'positive integer')
        end if

    end function node_lookup_positive_integer

!-----------------------------------------------------------------------------

    function node_lookup_string(node, key, val) result(ok)

        use ncsu_value
        use ncsu_constants

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        character(len=STRING_LENGTH), intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (.not. value_is_string(value)) &
                call type_mismatch(key, value, 'string')
            val = value_get_string(value)
        end if

    end function node_lookup_string

!-----------------------------------------------------------------------------

    function node_lookup_list(node, key, val) result(ok)

        use ncsu_value

        implicit none

        logical :: ok

        type(node_t), intent(in) :: node
        character(len=*), intent(in) :: key
        type(value_node_t), pointer :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (.not. value_is_list(value)) &
                call type_mismatch(key, value, 'list')
            val => value_get_list(value)
        end if

    end function node_lookup_list

!-----------------------------------------------------------------------------

    function node_lookup_logical(node, key, val) result(ok)

        use ncsu_value

        implicit none

        logical :: ok

        type(node_t), intent(in)  :: node
        character(len=*), intent(in)  :: key
        logical, intent(out) :: val

        type(value_t), pointer :: value

        value => node_lookup(node, key)
        ok = associated(value)
        if (ok) then
            if (value_is_string(value)) then
                if (value_get_string(value) == 'on') then
                    val = .true.
                    return
                else if (value_get_string(value) == 'off') then
                    val = .false.
                    return
                else
                    continue
                end if
            end if
            call type_mismatch(key, value, 'logical (on/off)')
        end if

    end function node_lookup_logical

end module ncsu_cftree
