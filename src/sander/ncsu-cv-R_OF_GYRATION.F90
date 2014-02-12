! <compile=optimized>
#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! cv%i = (i1, ..., iN) -- list of participating atoms
!

module ncsu_cv_R_OF_GYRATION

!=============================================================================

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

public :: colvar_cleanup

type, private :: priv_t
   NCSU_REAL :: Rg, cm(3) ! value & center of mass
   NCSU_REAL, pointer :: weights(:) ! m_i/total_mass
#  include "ncsu-cv-priv.type"
end type priv_t

#include "ncsu-cv-priv.decl"

!=============================================================================

contains

!=============================================================================

#include "ncsu-cv-priv.impl"

!=============================================================================

function colvar_value(cv, x) result(value)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   NCSU_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NCSU_REAL, intent(in) :: x(*)

   integer :: natoms, a, a3, i
   NCSU_REAL :: tmp

   type(priv_t), pointer :: priv

   ncsu_assert(cv%type == COLVAR_R_OF_GYRATION)
   ncsu_assert(associated(cv%i))

   natoms = size(cv%i)
   ncsu_assert(natoms > 2)

   priv => get_priv(cv)

   priv%Rg = ZERO
   priv%cm(1:3) = ZERO

   do a = 1, natoms
      a3 = 3*(cv%i(a) - 1)
      tmp = ZERO
      do i = 1, 3
         tmp = tmp + x(a3 + i)**2
         priv%cm(i) = priv%cm(i) + priv%weights(a)*x(a3 + i)
      end do
      priv%Rg = priv%Rg + priv%weights(a)*tmp
   end do

   priv%Rg = sqrt(priv%Rg - priv%cm(1)**2 - priv%cm(2)**2 - priv%cm(3)**2)
   value = priv%Rg

end function colvar_value

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine colvar_force(cv, x, fcv, f)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   NCSU_REAL, intent(in) :: x(*), fcv
   NCSU_REAL, intent(inout) :: f(*)

#ifdef MPI
#  include "ncsu-mpi.h"
   integer :: a_first, a_last
#endif

   integer :: natoms, a, a3
   NCSU_REAL :: tmp

   type(priv_t), pointer :: priv

   ncsu_assert(cv%type == COLVAR_R_OF_GYRATION)
   ncsu_assert(associated(cv%i))

   natoms = size(cv%i)
   ncsu_assert(natoms > 2)

   priv => get_priv(cv)
   ncsu_assert(priv%Rg > ZERO)

   tmp = fcv/priv%Rg

#ifdef MPI
   a = natoms/sandersize
   if (a.gt.0) then
      if (sanderrank.ne.(sandersize - 1)) then
         a_first = 1 + sanderrank*a
         a_last = (sanderrank + 1)*a
      else
         a_first = 1 + sanderrank*a
         a_last = natoms
      end if
   else
      if (sanderrank.eq.0) then
         a_first = 1
         a_last = natoms
      else
         a_first = 1
         a_last = 0
      end if
   end if
   do a = a_first, a_last
#else
   do a = 1, natoms
#endif /* MPI */
      a3 = 3*cv%i(a)

      f(a3 - 2:a3) = f(a3 - 2:a3) &
         + tmp*priv%weights(a)*(x(a3 - 2:a3) - priv%cm(1:3))
   end do

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_colvar_utils
   use ncsu_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NCSU_REAL,      intent(in)    :: amass(*)

   integer   :: natoms, a, error
   NCSU_REAL :: total_mass

   type(priv_t), pointer :: priv

#  include "ncsu-mpi.h"

   ncsu_assert(cv%type == COLVAR_R_OF_GYRATION)

   natoms = size(cv%i)

   call check_i(cv%i, cvno, 'R_OF_GYRATION')
   if (.not. natoms > 2) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NCSU_ERROR, 'CV #', cvno, &
            ' (R_OF_GYRATION) : too few integers'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if

   total_mass = ZERO

   do a = 1, natoms
      total_mass = total_mass + amass(cv%i(a))
   end do

   ncsu_assert(total_mass > ZERO)

   priv => new_priv(cv)

   allocate(priv%weights(natoms), stat = error)
   if (error.ne.0) &
      NCSU_OUT_OF_MEMORY

   do a = 1, natoms
      priv%weights(a) = amass(cv%i(a))/total_mass
   end do

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   NCSU_USE_AFAILED

   use ncsu_colvar_type
   use ncsu_colvar_utils

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   ncsu_assert(cv%type == COLVAR_R_OF_GYRATION)
   ncsu_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

subroutine colvar_cleanup(cv)

   NCSU_USE_AFAILED

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   ncsu_assert(cv%type.eq.COLVAR_R_OF_GYRATION)

   priv => get_priv(cv)
   ncsu_assert(associated(priv))

   deallocate(priv%weights)
   call del_priv(cv)

end subroutine colvar_cleanup

!=============================================================================

end module ncsu_cv_R_OF_GYRATION
