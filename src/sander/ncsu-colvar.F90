! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar

    implicit none

    private

!=============================================================================

    public :: colvar_value
    public :: colvar_force

    public :: colvar_difference

    public :: colvar_is_periodic

    public :: colvar_has_min
    public :: colvar_min

    public :: colvar_has_max
    public :: colvar_max

    public :: colvar_print
    public :: colvar_cleanup

    public :: colvar_mdread
    public :: colvar_bootstrap

!=============================================================================

contains

!=============================================================================

! the value is needed only on master
    function colvar_value(cv, x) result(value)

        NCSU_USE_AFAILED

        use ncsu_colvar_type

        use ncsu_cv_ANGLE, only : v_ANGLE => colvar_value
        use ncsu_cv_TORSION, only : v_TORSION => colvar_value
        use ncsu_cv_DISTANCE, only : v_DISTANCE => colvar_value
        use ncsu_cv_MULTI_RMSD, only : v_MULTI_RMSD => colvar_value
        use ncsu_cv_R_OF_GYRATION, only : v_R_OF_GYRATION => colvar_value
        use ncsu_cv_HANDEDNESS, only : v_HANDEDNESS => colvar_value
        use ncsu_cv_N_OF_BONDS, only : v_N_OF_BONDS => colvar_value
        use ncsu_cv_N_OF_STRUCTURES, only : v_N_OF_STRUCTURES => colvar_value
        use ncsu_cv_LCOD, only : v_LCOD => colvar_value
        use ncsu_cv_COS_OF_DIHEDRAL, only : v_COS_OF_DIHEDRAL => colvar_value
        use ncsu_cv_COM_ANGLE, only : v_COM_ANGLE => colvar_value
        use ncsu_cv_COM_TORSION, only : v_COM_TORSION => colvar_value
        use ncsu_cv_COM_DISTANCE, only : v_COM_DISTANCE => colvar_value
        use ncsu_cv_PCA, only : v_PCA => colvar_value
        implicit none

        NCSU_REAL :: value

        type(colvar_t) :: cv ! mutable
        NCSU_REAL, intent(in) :: x(*)

        select case (cv%type)
          case (COLVAR_ANGLE)
            value = v_ANGLE(cv, x)
          case (COLVAR_TORSION)
            value = v_TORSION(cv, x)
          case (COLVAR_DISTANCE)
            value = v_DISTANCE(cv, x)
          case (COLVAR_MULTI_RMSD)
            value = v_MULTI_RMSD(cv, x)
          case (COLVAR_R_OF_GYRATION)
            value = v_R_OF_GYRATION(cv, x)
          case (COLVAR_HANDEDNESS)
            value = v_HANDEDNESS(cv, x)
          case (COLVAR_N_OF_BONDS)
            value = v_N_OF_BONDS(cv, x)
          case (COLVAR_N_OF_STRUCTURES)
            value = v_N_OF_STRUCTURES(cv, x)
          case (COLVAR_LCOD)
            value = v_LCOD(cv, x)
          case (COLVAR_COS_OF_DIHEDRAL)
            value = v_COS_OF_DIHEDRAL(cv, x)
          case (COLVAR_COM_ANGLE)
            value = v_COM_ANGLE(cv, x)
          case (COLVAR_COM_TORSION)
            value = v_COM_TORSION(cv, x)
          case (COLVAR_COM_DISTANCE)
            value = v_COM_DISTANCE(cv, x)
          case (COLVAR_PCA)
            value = v_PCA(cv, x)
          case default
            ncsu_assert_not_reached()
            value = NCSU_TO_REAL(0)
        end select

    end function colvar_value

!=============================================================================

    subroutine colvar_force(cv, x, fcv, f)

        NCSU_USE_AFAILED

        use ncsu_colvar_type

        use ncsu_cv_ANGLE, only : f_ANGLE => colvar_force
        use ncsu_cv_TORSION, only : f_TORSION => colvar_force
        use ncsu_cv_DISTANCE, only : f_DISTANCE => colvar_force
        use ncsu_cv_MULTI_RMSD, only : f_MULTI_RMSD => colvar_force
        use ncsu_cv_R_OF_GYRATION, only : f_R_OF_GYRATION => colvar_force
        use ncsu_cv_HANDEDNESS, only : f_HANDEDNESS => colvar_force
        use ncsu_cv_N_OF_BONDS, only : f_N_OF_BONDS => colvar_force
        use ncsu_cv_N_OF_STRUCTURES, only : f_N_OF_STRUCTURES => colvar_force
        use ncsu_cv_LCOD, only : f_LCOD => colvar_force
        use ncsu_cv_COS_OF_DIHEDRAL, only : f_COS_OF_DIHEDRAL => colvar_force
        use ncsu_cv_COM_ANGLE, only : f_COM_ANGLE => colvar_force
        use ncsu_cv_COM_TORSION, only : f_COM_TORSION => colvar_force
        use ncsu_cv_COM_DISTANCE, only : f_COM_DISTANCE => colvar_force
        use ncsu_cv_PCA, only : f_PCA => colvar_force

        implicit none

        type(colvar_t) :: cv ! mutable

        NCSU_REAL, intent(in) :: x(*), fcv
        NCSU_REAL, intent(inout) :: f(*)

        select case (cv%type)
          case (COLVAR_ANGLE)
            call f_ANGLE(cv, x, fcv, f)
          case (COLVAR_TORSION)
            call f_TORSION(cv, x, fcv, f)
          case (COLVAR_DISTANCE)
            call f_DISTANCE(cv, x, fcv, f)
          case (COLVAR_MULTI_RMSD)
            call f_MULTI_RMSD(cv, x, fcv, f)
          case (COLVAR_R_OF_GYRATION)
            call f_R_OF_GYRATION(cv, x, fcv, f)
          case (COLVAR_HANDEDNESS)
            call f_HANDEDNESS(cv, x, fcv, f)
          case (COLVAR_N_OF_BONDS)
            call f_N_OF_BONDS(cv, x, fcv, f)
          case (COLVAR_N_OF_STRUCTURES)
            call f_N_OF_STRUCTURES(cv, x, fcv, f)
          case (COLVAR_LCOD)
            call f_LCOD(cv, x, fcv, f)
          case (COLVAR_COS_OF_DIHEDRAL)
            call f_COS_OF_DIHEDRAL(cv, x, fcv, f)
          case (COLVAR_COM_ANGLE)
            call f_COM_ANGLE(cv, x, fcv, f)
          case (COLVAR_COM_TORSION)
            call f_COM_TORSION(cv, x, fcv, f)
          case (COLVAR_COM_DISTANCE)
            call f_COM_DISTANCE(cv, x, fcv, f)
          case (COLVAR_PCA)
            call f_PCA(cv, x, fcv, f)
          case default
            ncsu_assert_not_reached()
        end select

    end subroutine colvar_force

!=============================================================================

    function colvar_difference(cv, v1, v2) result(diff)

        NCSU_USE_AFAILED

        use ncsu_colvar_type
        use ncsu_constants, only : ZERO
        use constants, only : PI

        implicit none

        NCSU_REAL :: diff
        type(colvar_t), intent(in) :: cv
        NCSU_REAL, intent(in) :: v1, v2

        NCSU_REAL :: t1, t2

        t1 = fix_value(v1)
        t2 = fix_value(v2)

        diff = t1 - t2

        if (cv%type .eq. COLVAR_TORSION .or. cv%type .eq. COLVAR_COM_TORSION) then
            ncsu_assert(-PI .le. t1 .and. t1 .le. PI)
            ncsu_assert(-PI .le. t2 .and. t2 .le. PI)
            if (diff .gt. PI) then
                diff = diff - PI - PI
            else if (diff .lt. -PI) then
                diff = diff + PI + PI
            end if
        end if

    contains

        function fix_value(v) result(t)

            implicit none

            NCSU_REAL :: t
            NCSU_REAL, intent(in) :: v

            if (cv%type .eq. COLVAR_ANGLE .or. cv%type .eq. COLVAR_COM_ANGLE) then
                if (ZERO .le. v .and. v .le. PI) then
                    t = v
                else
                    t = acos(cos(v))
                end if
            else if (cv%type .eq. COLVAR_TORSION .or. cv%type .eq. COLVAR_COM_TORSION) then
                if (-PI .le. v .and. v .le. PI) then
                    t = v
                else
                    t = atan2(sin(v), cos(v))
                end if
            else
                t = v
            end if

        end function fix_value

    end function colvar_difference

!=============================================================================

    logical function colvar_is_periodic(cv)

        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(in) :: cv

        if (cv%type .eq. COLVAR_TORSION .or. cv%type .eq. COLVAR_COM_TORSION) then
            colvar_is_periodic = .true.
        else
            colvar_is_periodic = .false.
        end if

    end function colvar_is_periodic

!=============================================================================

    logical function colvar_has_min(cv)

        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(in) :: cv

        select case (cv%type)
          case (COLVAR_ANGLE:COLVAR_R_OF_GYRATION)
            colvar_has_min = .true.
          case (COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
            colvar_has_min = .true.
          case (COLVAR_COM_ANGLE:COLVAR_COM_DISTANCE)
            colvar_has_min = .true.
          case default
            colvar_has_min = .false.
        end select

    end function colvar_has_min

!=============================================================================

    NCSU_REAL function colvar_min(cv)

    NCSU_USE_AFAILED

    use constants, only : PI
    use ncsu_constants, only : zero
    use ncsu_colvar_type

    implicit none

    type(colvar_t), intent(in) :: cv

    ncsu_assert(colvar_has_min(cv))

    select case (cv%type)
      case (COLVAR_ANGLE)
        colvar_min = ZERO
      case (COLVAR_TORSION)
        colvar_min = -PI
      case (COLVAR_DISTANCE:COLVAR_R_OF_GYRATION)
        colvar_min = ZERO
      case (COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
        colvar_min = ZERO
      case (COLVAR_COM_ANGLE)
        colvar_min = ZERO
      case (COLVAR_COM_TORSION)
        colvar_min = -PI
      case (COLVAR_COM_DISTANCE)
        colvar_min = ZERO
      case default
        ncsu_assert_not_reached()
        colvar_min = ZERO
    end select

end module ncsu_colvar

!=============================================================================

logical function colvar_has_max(cv)

    use ncsu_colvar_type

    implicit none

    type(colvar_t), intent(in) :: cv

    select case (cv%type)
      case (COLVAR_ANGLE:COLVAR_TORSION)
        colvar_has_max = .true.
      case (COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
        colvar_has_max = .true.
      case default
        colvar_has_max = .false.
    end select

end function colvar_has_max

!=============================================================================

NCSU_REAL function colvar_max(cv)

NCSU_USE_AFAILED

use constants, only : PI
use ncsu_constants, only : ZERO
use ncsu_colvar_type

implicit none

type(colvar_t), intent(in) :: cv

ncsu_assert(colvar_has_max(cv))

select case (cv%type)
  case (COLVAR_ANGLE:COLVAR_TORSION)
    colvar_max = PI
  case (COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
    colvar_max = PI
  case default
    ncsu_assert_not_reached()
    colvar_max = ZERO
end select

end function colvar_max

!=============================================================================

subroutine colvar_print(cv, lun)

    use ncsu_utils
    use ncsu_constants
    use ncsu_colvar_type
    use ncsu_sander_proxy

    use ncsu_cv_ANGLE, only : p_ANGLE => print_details
    use ncsu_cv_TORSION, only : p_TORSION => print_details
    use ncsu_cv_DISTANCE, only : p_DISTANCE => print_details
    use ncsu_cv_MULTI_RMSD, only : p_MULTI_RMSD => print_details
    use ncsu_cv_R_OF_GYRATION, only : p_R_OF_GYRATION => print_details
    use ncsu_cv_HANDEDNESS, only : p_HANDEDNESS => print_details
    use ncsu_cv_N_OF_BONDS, only : p_N_OF_BONDS => print_details
    use ncsu_cv_N_OF_STRUCTURES, only : p_N_OF_STRUCTURES => print_details
    use ncsu_cv_LCOD, only : p_LCOD => print_details
    use ncsu_cv_COS_OF_DIHEDRAL, only : p_COS_OF_DIHEDRAL => print_details
    use ncsu_cv_COM_ANGLE, only : p_COM_ANGLE => print_details
    use ncsu_cv_COM_TORSION, only : p_COM_TORSION => print_details
    use ncsu_cv_COM_DISTANCE, only : p_COM_DISTANCE => print_details
    use ncsu_cv_PCA, only : p_PCA => print_details

    implicit none

    type(colvar_t), intent(in) :: cv
    integer, intent(in) :: lun

    write (unit=lun, fmt='(a,a)', advance='NO') NCSU_INFO, '  type = '''

    select case (cv%type)
      case (COLVAR_ANGLE)
        write (unit=lun, fmt='(a)') 'ANGLE'''
        call p_ANGLE(cv, lun)
      case (COLVAR_TORSION)
        write (unit=lun, fmt='(a)') 'TORSION'''
        call p_TORSION(cv, lun)
      case (COLVAR_DISTANCE)
        write (unit=lun, fmt='(a)') 'DISTANCE'''
        call p_DISTANCE(cv, lun)
      case (COLVAR_MULTI_RMSD)
        write (unit=lun, fmt='(a)') 'MULTI_RMSD'''
        call p_MULTI_RMSD(cv, lun)
      case (COLVAR_R_OF_GYRATION)
        write (unit=lun, fmt='(a)') 'R_OF_GYRATION'''
        call p_R_OF_GYRATION(cv, lun)
      case (COLVAR_HANDEDNESS)
        write (unit=lun, fmt='(a)') 'HANDEDNESS'''
        call p_HANDEDNESS(cv, lun)
      case (COLVAR_N_OF_BONDS)
        write (unit=lun, fmt='(a)') 'N_OF_BONDS'''
        call p_N_OF_BONDS(cv, lun)
      case (COLVAR_N_OF_STRUCTURES)
        write (unit=lun, fmt='(a)') 'N_OF_STRUCTURES'''
        call p_N_OF_STRUCTURES(cv, lun)
      case (COLVAR_LCOD)
        write (unit=lun, fmt='(a)') &
            'LCOD'' (Linear Combination Of Distances)'
        call p_LCOD(cv, lun)
      case (COLVAR_COS_OF_DIHEDRAL)
        write (unit=lun, fmt='(a)') 'COS_OF_DIHEDRAL'''
        call p_COS_OF_DIHEDRAL(cv, lun)
      case (COLVAR_COM_ANGLE)
        write (unit=lun, fmt='(a)') 'COM_ANGLE'''
        call p_COM_ANGLE(cv, lun)
      case (COLVAR_COM_TORSION)
        write (unit=lun, fmt='(a)') 'COM_TORSION'''
        call p_COM_TORSION(cv, lun)
      case (COLVAR_COM_DISTANCE)
        write (unit=lun, fmt='(a)') 'COM_DISTANCE'''
        call p_COM_DISTANCE(cv, lun)
      case (COLVAR_PCA)
        write (unit=lun, fmt='(a)') 'PRINCIPAL COMPONENT'''
        call p_PCA(cv, lun)
      case default
        ncsu_assert_not_reached()
        continue
    end select

end subroutine colvar_print

!=============================================================================

subroutine colvar_cleanup(cv)

    use ncsu_colvar_type

    use ncsu_cv_MULTI_RMSD, only : c_MULTI_RMSD => colvar_cleanup
    use ncsu_cv_R_OF_GYRATION, only : c_R_OF_GYRATION => colvar_cleanup
    use ncsu_cv_N_OF_STRUCTURES, only : c_N_OF_STRUCTURES => colvar_cleanup
    use ncsu_cv_PCA, only : c_PCA => colvar_cleanup

    implicit none

    type(colvar_t), intent(inout) :: cv

    select case (cv%type)
      case (COLVAR_MULTI_RMSD)
        call c_MULTI_RMSD(cv)
      case (COLVAR_R_OF_GYRATION)
        call c_R_OF_GYRATION(cv)
      case (COLVAR_N_OF_STRUCTURES)
        call c_N_OF_STRUCTURES(cv)
      case (COLVAR_PCA)
        call c_PCA(cv)
      case default
        continue
    end select

    if (associated(cv%i)) &
        deallocate (cv%i)

    if (associated(cv%r)) &
        deallocate (cv%r)

    if (associated(cv%avgcrd)) &
        deallocate (cv%avgcrd)

    if (associated(cv%evec)) &
        deallocate (cv%evec)

    if (associated(cv%state_ref)) &
        deallocate (cv%state_ref)

    if (associated(cv%state_pca)) &
        deallocate (cv%state_pca)

    if (associated(cv%ipca_to_i)) &
        deallocate (cv%ipca_to_i)

    cv%type = -1

end subroutine colvar_cleanup

!=============================================================================

subroutine colvar_mdread(cv, node, cvno)

    use ncsu_utils
    use ncsu_value
    use ncsu_cftree
    use ncsu_constants
    use ncsu_colvar_type
    use ncsu_sander_proxy
    use ncsu_read_pca

    implicit none

    type(colvar_t), intent(inout) :: cv
    type(node_t), intent(in)    :: node
    integer, intent(in)    :: cvno

    integer :: n, error
    integer :: found, crdsize
    logical :: found2
    integer :: nsolut

    type(value_node_t), pointer :: alist, aiter

    character(len=STRING_LENGTH) :: type

    ! declare three strings: ref_file avg_file evec_file
    character(len=STRING_LENGTH) :: ref_file
    character(len=STRING_LENGTH) :: avg_file
    character(len=STRING_LENGTH) :: evec_file
    character(len=STRING_LENGTH) :: index_file
    integer :: first, last

    ncsu_assert(is_master())

    ncsu_assert(.not. associated(cv%i))
    ncsu_assert(.not. associated(cv%r))
    ncsu_assert(.not. associated(cv%avgcrd))
    ncsu_assert(.not. associated(cv%evec))
    ncsu_assert(.not. associated(cv%state_ref))
    ncsu_assert(.not. associated(cv%state_pca))

    ncsu_assert(node_title(node) == 'variable')

    !
    ! type
    !

    if (.not. node_lookup_string(node, 'type', type)) then
        write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//'/)') &
            NCSU_ERROR, 'type is not specified for CV #', cvno
        call terminate()
    end if

    if (type == 'ANGLE') then
        cv%type = COLVAR_ANGLE
    else if (type == 'TORSION') then
        cv%type = COLVAR_TORSION
    else if (type == 'DISTANCE') then
        cv%type = COLVAR_DISTANCE
    else if (type == 'MULTI_RMSD') then
        cv%type = COLVAR_MULTI_RMSD
    else if (type == 'R_OF_GYRATION') then
        cv%type = COLVAR_R_OF_GYRATION
    else if (type == 'HANDEDNESS') then
        cv%type = COLVAR_HANDEDNESS
    else if (type == 'N_OF_BONDS') then
        cv%type = COLVAR_N_OF_BONDS
    else if (type == 'N_OF_STRUCTURES') then
        cv%type = COLVAR_N_OF_STRUCTURES
    else if (type == 'LCOD') then
        cv%type = COLVAR_LCOD
    else if (type == 'COS_OF_DIHEDRAL') then
        cv%type = COLVAR_COS_OF_DIHEDRAL
    else if (type == 'COM_ANGLE') then
        cv%type = COLVAR_COM_ANGLE
    else if (type == 'COM_TORSION') then
        cv%type = COLVAR_COM_TORSION
    else if (type == 'COM_DISTANCE') then
        cv%type = COLVAR_COM_DISTANCE
    else if (type == 'PCA') then
        cv%type = COLVAR_PCA
        nsolut = sander_nsolut()
    else
        write (unit=ERR_UNIT, fmt='(/a,a,a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV type ''', trim(type), &
            ''' is not supported so far (CV #', cvno, ')'
        call terminate()
    end if

    !
    ! cv%i
    !

    if (node_lookup_list(node, 'i', alist)) then

        n = 0
        aiter => alist

        do while (associated(aiter))
            n = n + 1
            if (.not. value_is_integer(aiter%value)) then
                write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                    NCSU_ERROR, 'CV #', cvno, ' : unexpected &
                &(not an integer) element of ''i'' list'
                call terminate()
            end if
            aiter => aiter%next
        end do

        if (n > 0) then
            allocate (cv%i(n), stat=error)
            if (error /= 0) &
                NCSU_OUT_OF_MEMORY

            n = 0
            aiter => alist
            do while (associated(aiter))
                n = n + 1
                cv%i(n) = value_get_integer(aiter%value)
                aiter => aiter%next
            end do
        end if

    end if ! node_lookup_list(vnode, 'i', alist))

    !
    ! cv%r
    !

    if (type /= 'PCA') then
        if (node_lookup_list(node, 'r', alist)) then

            n = 0
            aiter => alist

            do while (associated(aiter))
                n = n + 1
                if (.not. value_is_real(aiter%value)) then
                    write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                        NCSU_ERROR, 'CV #', cvno, ' : unexpected &
                    &(not a real number) element of ''r'' list'
                    call terminate()
                end if
                aiter => aiter%next
            end do

            if (n > 0) then
                allocate (cv%r(n), stat=error)
                if (error /= 0) &
                    NCSU_OUT_OF_MEMORY

                n = 0
                aiter => alist
                do while (associated(aiter))
                    n = n + 1
                    cv%r(n) = value_get_real(aiter%value)
                    aiter => aiter%next
                end do
            end if

        end if ! node_lookup_list(vnode, 'r', alist))
    end if

    ! read refcrd from "refcrd"
    ! store in cv% r
    if (type == 'PCA') then

        found2 = node_lookup_string(node, 'refcrd', ref_file)
        if (.not. found2) then
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(n)//'/)') &
                NCSU_ERROR, 'could not find ''reference coordinates'' for CV #', n
            call terminate()
        else
            ! allocate and read reference crd file into CV
            if (cv%i(1) > 0) then
                allocate (cv%r(cv%i(1)*3), stat=error)
                if (error /= 0) &
                    NCSU_OUT_OF_MEMORY
                call read_refcrd(cv, ref_file)
                write (unit=OUT_UNIT, fmt='(a,a,a,a)', advance='NO') NCSU_INFO, &
                    'ref_file = ', trim(ref_file), ' ('
                write (unit=OUT_UNIT, fmt='(a)') 'loaded)'
            end if
        end if

        ! cv % evec
        found2 = node_lookup_string(node, 'evec', evec_file)
        if (.not. found2) then
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(n)//'/)') &
                NCSU_ERROR, 'could not find ''reference coordinates'' for CV #', n
            call terminate()
        else
            ! allocate and read evec crd into CV
            if (cv%i(3) > 0) then
                allocate (cv%evec(cv%i(3)*3), stat=error)
                if (error /= 0) &
                    NCSU_OUT_OF_MEMORY
                call read_evec(cv, evec_file)
                write (unit=OUT_UNIT, fmt='(a,a,a,a)', advance='NO') NCSU_INFO, &
                    'evec_file = ', trim(evec_file), ' ('
                write (unit=OUT_UNIT, fmt='(a)') 'loaded)'
            end if
        end if

        ! cv % avgcrd
        found2 = node_lookup_string(node, 'avgcrd', avg_file)
        if (.not. found2) then
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(n)//'/)') &
                NCSU_ERROR, 'could not find ''average coordinates'' for CV #', n
            call terminate()
        else

            ! allocate and read average crd into CV
            if (cv%i(3) > 0) then
                allocate (cv%avgcrd(cv%i(3)*3), stat=error)
                if (error /= 0) &
                    NCSU_OUT_OF_MEMORY

                call read_avgcrd(cv, avg_file)
                write (unit=OUT_UNIT, fmt='(a,a,a,a)', advance='NO') NCSU_INFO, &
                    'avg_file = ', trim(avg_file), ' ('
                write (unit=OUT_UNIT, fmt='(a)') 'loaded)'
            end if
        end if

        ! cv % index
        found2 = node_lookup_string(node, 'index', index_file)
        if (.not. found2) then
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(n)//'/)') &
                NCSU_ERROR, 'could not find ''index for ref and pca part'' for CV #', n
            call terminate()
        else
            ! allocate and read average crd into CV
            if (cv%i(1) > 0 .and. cv%i(3) > 0) then
                allocate (cv%state_ref(cv%i(1)), cv%state_pca(cv%i(1)), cv%ipca_to_i(cv%i(3)), stat=error)
                if (error /= 0) &
                    NCSU_OUT_OF_MEMORY
                call read_index(cv, index_file)
                write (unit=OUT_UNIT, fmt='(a,a,a,a)', advance='NO') NCSU_INFO, &
                    'index_file = ', trim(index_file), ' ('
                write (unit=OUT_UNIT, fmt='(a)') 'loaded)'
            end if
        end if

    end if

end subroutine colvar_mdread

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

    use ncsu_utils
    use ncsu_colvar_type

    use ncsu_cv_ANGLE, only : b_ANGLE => colvar_bootstrap
    use ncsu_cv_TORSION, only : b_TORSION => colvar_bootstrap
    use ncsu_cv_DISTANCE, only : b_DISTANCE => colvar_bootstrap
    use ncsu_cv_MULTI_RMSD, only : b_MULTI_RMSD => colvar_bootstrap
    use ncsu_cv_R_OF_GYRATION, only : b_R_OF_GYRATION => colvar_bootstrap
    use ncsu_cv_HANDEDNESS, only : b_HANDEDNESS => colvar_bootstrap
    use ncsu_cv_N_OF_BONDS, only : b_N_OF_BONDS => colvar_bootstrap
    use ncsu_cv_N_OF_STRUCTURES, only : b_N_OF_STRUCTURES => colvar_bootstrap
    use ncsu_cv_LCOD, only : b_LCOD => colvar_bootstrap
    use ncsu_cv_COS_OF_DIHEDRAL, only : b_COS_OF_DIHEDRAL => colvar_bootstrap
    use ncsu_cv_COM_ANGLE, only : b_COM_ANGLE => colvar_bootstrap
    use ncsu_cv_COM_TORSION, only : b_COM_TORSION => colvar_bootstrap
    use ncsu_cv_COM_DISTANCE, only : b_COM_DISTANCE => colvar_bootstrap
    use ncsu_cv_PCA, only : b_PCA => colvar_bootstrap

    implicit none

    type(colvar_t), intent(inout) :: cv
    integer, intent(in)    :: cvno
    NCSU_REAL, intent(in)    :: amass(*)

#ifdef MPI
    integer :: bcastdata(8), ierr
#include "ncsu-mpi.h"

    !
    ! bcast type/i/r first
    !

#ifndef NCSU_DISABLE_ASSERT
    if (sanderrank == 0) then
        ncsu_assert(cv%type > 0)
    else
        ncsu_assert(.not. associated(cv%i))
        ncsu_assert(.not. associated(cv%r))
        if (cv%type == COLVAR_PCA) then
            ncsu_assert(.not. associated(cv%avgcrd))
            ncsu_assert(.not. associated(cv%evec))
            ncsu_assert(.not. associated(cv%state_ref))
            ncsu_assert(.not. associated(cv%state_pca))
            ncsu_assert(.not. associated(cv%ipca_to_i))
        end if
    end if
#endif /* NCSU_DISABLE_ASSERT */

    if (sanderrank == 0) then
        bcastdata(1) = cv%type

        bcastdata(2) = 0
        if (associated(cv%i)) &
            bcastdata(2) = size(cv%i)

        bcastdata(3) = 0
        if (associated(cv%r)) &
            bcastdata(3) = size(cv%r)

        bcastdata(4) = 0
        bcastdata(5) = 0
        bcastdata(6) = 0
        bcastdata(7) = 0
        bcastdata(8) = 0

        if (cv%type == COLVAR_PCA) then

            if (associated(cv%avgcrd)) &
                bcastdata(4) = size(cv%avgcrd)

            if (associated(cv%evec)) &
                bcastdata(5) = size(cv%evec)

            if (associated(cv%state_ref)) &
                bcastdata(6) = size(cv%state_ref)

            if (associated(cv%state_pca)) &
                bcastdata(7) = size(cv%state_pca)

            if (associated(cv%ipca_to_i)) &
                bcastdata(8) = size(cv%ipca_to_i)

        end if

    end if ! sanderrank == 0

    call mpi_bcast(bcastdata, size(bcastdata), MPI_INTEGER, 0, commsander, ierr)
    ncsu_assert(ierr == 0)

    if (sanderrank /= 0) &
        cv%type = bcastdata(1)

    !
    ! cv%i
    !

    if (bcastdata(2) > 0) then
        if (.not. associated(cv%i)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%i(bcastdata(2)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%i)

        call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, 0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%i)
    end if ! bcastdata(2) > 0

    !
    ! cv%r
    !

    if (bcastdata(3) > 0) then
        if (.not. associated(cv%r)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%r(bcastdata(3)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%r)

        call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%r)
    end if

    !
    ! cv % avgcrd
    !
    if (bcastdata(4) > 0) then
        if (.not. associated(cv%avgcrd)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%avgcrd(bcastdata(4)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%avgcrd)

        call mpi_bcast(cv%avgcrd, bcastdata(4), MPI_DOUBLE_PRECISION, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%avgcrd)
    end if ! bcastdata(4) > 0

    !
    ! cv % evec
    !
    if (bcastdata(5) > 0) then
        if (.not. associated(cv%evec)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%evec(bcastdata(5)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%evec)

        call mpi_bcast(cv%evec, bcastdata(5), MPI_DOUBLE_PRECISION, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%evec)
    end if ! bcastdata(5) > 0

    !
    ! cv % state_ref
    !
    if (bcastdata(6) > 0) then
        if (.not. associated(cv%state_ref)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%state_ref(bcastdata(6)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%state_ref)

        call mpi_bcast(cv%state_ref, bcastdata(6), MPI_INTEGER, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%state_ref)
    end if ! bcastdata(6) > 0

    !
    ! cv % state_pca
    !
    if (bcastdata(7) > 0) then
        if (.not. associated(cv%state_pca)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%state_pca(bcastdata(7)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%state_pca)

        call mpi_bcast(cv%state_pca, bcastdata(7), MPI_INTEGER, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%state_pca)
    end if ! bcastdata(7) > 0

    !
    ! cv % ipca_to_i()
    !
    if (bcastdata(8) > 0) then
        if (.not. associated(cv%ipca_to_i)) then
            ncsu_assert(sanderrank > 0)
            allocate (cv%ipca_to_i(bcastdata(8)), stat=ierr)
            if (ierr /= 0) &
                NCSU_OUT_OF_MEMORY
        end if ! .not. associated(cv%state_pca)

        call mpi_bcast(cv%ipca_to_i, bcastdata(8), MPI_INTEGER, &
            0, commsander, ierr)
        ncsu_assert(ierr == 0)
    else
        nullify (cv%ipca_to_i)
    end if ! bcastdata(8) > 0

#endif /* MPI */

    !
    ! dispatch according to the type
    !

    select case (cv%type)
      case (COLVAR_ANGLE)
        call b_ANGLE(cv, cvno, amass)
      case (COLVAR_TORSION)
        call b_TORSION(cv, cvno, amass)
      case (COLVAR_DISTANCE)
        call b_DISTANCE(cv, cvno, amass)
      case (COLVAR_MULTI_RMSD)
        call b_MULTI_RMSD(cv, cvno, amass)
      case (COLVAR_R_OF_GYRATION)
        call b_R_OF_GYRATION(cv, cvno, amass)
      case (COLVAR_HANDEDNESS)
        call b_HANDEDNESS(cv, cvno, amass)
      case (COLVAR_N_OF_BONDS)
        call b_N_OF_BONDS(cv, cvno, amass)
      case (COLVAR_N_OF_STRUCTURES)
        call b_N_OF_STRUCTURES(cv, cvno, amass)
      case (COLVAR_LCOD)
        call b_LCOD(cv, cvno, amass)
      case (COLVAR_COS_OF_DIHEDRAL)
        call b_COS_OF_DIHEDRAL(cv, cvno, amass)
      case (COLVAR_COM_ANGLE)
        call b_COM_ANGLE(cv, cvno, amass)
      case (COLVAR_COM_TORSION)
        call b_COM_TORSION(cv, cvno, amass)
      case (COLVAR_COM_DISTANCE)
        call b_COM_DISTANCE(cv, cvno, amass)
      case (COLVAR_PCA)
        call b_PCA(cv, cvno, amass)
      case default
        ncsu_assert_not_reached()
        continue
    end select

end subroutine colvar_bootstrap

!=============================================================================

end module ncsu_colvar
