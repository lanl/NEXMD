! <compile=optimized>
#include "ncsu-utils.h"
#include "ncsu-config.h"
! By Sishi Tang and Lin Fu, 
! June 17, 2011 
!
! input:
!
! cv%i = (a1, a2)
!
!     (a1, a2 : first and last atom for PCA projection) 
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, b1x, ...)
!
!        (reference coordinates)  
! cv%avgcrd = (b1x, b1y, b1z, b2x, b2y, b2z,... ) 
! cv%evec = (v1x, v1y, ... ) 
! 
! value = ((fitcrd-avgcrd)*transpose(evec)) 
!
! Note: for v2, all solute atoms are fitted, and the projection is only calculated 
!       for selected atoms from a1 - a2. 
! dimension of arrays (clarification): 
!                   evec:   3*n_solut 
!                   avgcrd: 3*n_solut 
!                   refcrd: 3*n_solut 
!                   cm_crd: 3*n_solut 
!                  fit_crd: 3*natoms   

! corrected one:    evec:   3*size(cv%i)
!                   avgcrd: 3*size(cv%i)
!                   refcrd: 3*size(cv%j) read from both j and ref.crd
!                   cm_crd: 3*size(cv%j)
!                  fit_crd: 3*size(cv%i)

module ncsu_cv_PCA

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
   
   NCSU_REAL :: value, total_mass 
   NCSU_REAL :: ref_cm(3), quaternion(4), U(3,3) ! COM of ref, Q, and rot matrix U  
   NCSU_REAL, pointer :: mass(:) => null()
   NCSU_REAL, pointer :: ref_crd(:) => null() ! translated ref  
   NCSU_REAL, pointer :: cm_crd(:)  => null()  ! translated orig 
!  NCSU_REAL, pointer :: fit_crd => null()  ! do we need this ?
    
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
   use ncsu_rmsd
   use ncsu_sander_proxy
   implicit none

   NCSU_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NCSU_REAL, intent(in) :: x(*)
   NCSU_REAL :: x_cm(3), lambda 
   NCSU_REAL, pointer :: fit_crd(:) => null() 
   integer :: natoms, a, a3, i1, i2, i3, n, i, error, nsolut 

   type(priv_t), pointer :: priv

#  ifdef MPI
#     include "ncsu-mpi.h"
#  endif /* MPI */

   ncsu_assert(cv%type == COLVAR_PCA)
   ncsu_assert(associated(cv%i))

!  i1: nsolut
!  i2: nref
!  i3: npca

   i1 = cv%i(1) 
   i2 = cv%i(2)
   i3 = cv%i(3)
  
!   nsolut = sander_nsolut() 
 
!   natoms = i2 - i1 + 1
!   ncsu_assert(natoms > 1)

    natoms = i1
    ncsu_assert(natoms > 1)
    
!    ncsu_assert(i1 == nsolut)

!   allocate(fit_crd(3*natoms), stat = error)   
   allocate(fit_crd(3*cv%i(3)), stat = error)
   
   if (error /= 0) &
      NCSU_OUT_OF_MEMORY

   priv => get_priv(cv)

   priv%value = ZERO
   x_cm = ZERO
   
   ! compute COM of moving atoms for all solutes 
   n = 1 
   do a = 1, cv%i(1)
   	  if(cv%state_ref(a) == 0) cycle
      a3 = 3*a 
!      x_cm = x_cm + priv%mass(n)*x(a3 - 2:a3)  
      x_cm = x_cm + priv%mass(a)*x(a3 - 2:a3)  

      n = n + 1 
   end do

   ncsu_assert(i2 == n-1)

   x_cm = x_cm/priv%total_mass 
   
!   write(*,*) "x_cm=, worldrank=", x_cm, worldrank
   
   ! computer translated moving atoms for all solutes     
   n = 1
   do a = 1, cv%i(1)
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(a - 1) + i) - x_cm(i)
      end do
      n = n + 1
   end do

   ! Calculate quaternion 
   call rmsd_q1(natoms, cv%state_ref, priv%mass, priv%cm_crd,& 
               priv%ref_crd, lambda, priv%quaternion)
  
   ! Calculate transposed rotation matrix U 
   call rmsd_q3u(priv%quaternion, priv%U) 
   
!   write(*,*) "priv%U=, worldrank=", priv%U, worldrank
 
   ! Calculate fitted crd (from i1-i2 ONLY) wrt to original ref_crd     
   ! best fit X for ref_crd : U*cm_crd + ref_cm 
   ! then subject <X> to get zero mean 
   ! make it mass weighted 
   ! avgcrd includes i1:i2  
   ! cm_crd includes 1:nsolut
   ! fit_crd includes i1:i2  
   n = 1 
!   do a = i1, i2
    do a = 1, cv%i(1)
    	if(cv%state_pca(a) == 0) cycle
      a3 = 3*a
!      fit_crd(3*n-2:3*n) = (matmul(priv%U, priv%cm_crd(a3-2:a3)) &
!                        + priv%ref_cm - cv%avgcrd(a3-2:a3))  & 
!                        * sqrt(priv%mass(a)) 

       fit_crd(3*n-2:3*n) = (matmul(priv%U, priv%cm_crd(a3-2:a3)) &
                         + priv%ref_cm - cv%avgcrd(3*n-2:3*n))  & 
                         * sqrt(priv%mass(a)) 
      n = n + 1 

   end do

   ncsu_assert(i3 == n-1)


   ! Calculate projection
   ! P = XT  
!   priv%value = dot_product(fit_crd,cv%evec(i1*3-2:i2*3)) 
   priv%value = dot_product(fit_crd,cv%evec(1:3*i3))
   value = priv%value  
   

   deallocate(fit_crd)

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

   integer :: natoms, npca, a, a3, i1, i2, i3  

   type(priv_t), pointer :: priv

   ncsu_assert(cv%type == COLVAR_PCA)
   ncsu_assert(associated(cv%i))

   i1 = cv%i(1) 
   i2 = cv%i(2)
   i3 = cv%i(3)
!   natoms = i2 - i1 + 1
!   natoms = i1
   npca = i3
   
!   ncsu_assert(natoms > 1)
   ncsu_assert(npca > 1)

   priv => get_priv(cv)
   !ncsu_assert(priv%value > ZERO)
   ncsu_assert(priv%total_mass > ZERO)

#ifdef MPI
!   a = natoms/sandersize
   a = npca/sandersize
   if (a.gt.0) then
      if (sanderrank.ne.(sandersize - 1)) then
         a_first = 1 + sanderrank*a
         a_last = (sanderrank + 1)*a
      else
         a_first = 1 + sanderrank*a
         ! a_last = natoms
         a_last = npca
         
      end if
   else
      if (sanderrank.eq.0) then
         a_first = 1
         ! a_last = natoms
         a_last = npca
      else
         a_first = 1
         a_last = 0
      end if
   end if
   do a = a_first, a_last
#else
!   do a = 1, natoms
   do a = 1, npca
#endif /* MPI */
!      a3 = 3*(i1 + a - 1)
       a3 = 3*cv%ipca_to_i(a)
!      if(state_pca(a) == 0) cycle
      ! dc/dx = transpose(T)*U 
!      f(a3 - 2:a3) = f(a3 - 2:a3) + fcv * matmul(cv%evec(a3 - 2:a3), priv%U)*sqrt(priv%mass(i1 + a - 1))! now evec has the same index as force 
       f(a3 - 2:a3) = f(a3 - 2:a3) + fcv * matmul(cv%evec(3*a - 2 : 3*a), priv%U)*sqrt(priv%mass(cv%ipca_to_i(a)))
       
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

   integer   :: a, b, error, crdsize, nsolut, natoms  
   NCSU_REAL :: total_mass

   type(priv_t), pointer :: priv

#  include "ncsu-mpi.h"

   ncsu_assert(cv%type == COLVAR_PCA)

!   natoms = cv%i(2) - cv%i(1) + 1 
!   nsolut = sander_nsolut()
!   ncsu_assert(nsolut == cv%i(1))

!  make sure that there are only three integers 
   call check_i(cv%i, cvno, 'PCA', 3)
   if (.not. cv%i(3) >= 2) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NCSU_ERROR, 'CV #', cvno, &
            ' (PCA) : too few integers'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if

   priv => new_priv(cv)

   allocate(priv%mass(cv%i(1)), priv%ref_crd(3*cv%i(1)), &
            priv%cm_crd(3*cv%i(1)), stat = error)      
   

   if (error.ne.0) &
      NCSU_OUT_OF_MEMORY

   priv%ref_cm = ZERO
   priv%total_mass = ZERO
  
   total_mass = ZERO

   ! calculate total mass and refcm 
   do a = 1, cv%i(1)
   	  priv%mass(a) = amass(a)
   	  if(cv%state_ref(a) == 0 ) cycle
      total_mass = total_mass + priv%mass(a)
      priv%ref_cm = priv%ref_cm + priv%mass(a)*cv%r(3*a-2:3*a)
   end do

   priv%total_mass = total_mass 
   priv%ref_cm = priv%ref_cm/total_mass

   ! translate reference coordinates to CM frame
   ! ref_crd and r should be the same size 
   do a = 1, cv%i(1)
     do b = 1, 3
        priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b ) - priv%ref_cm(b)
     end do
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

   ncsu_assert(cv%type == COLVAR_PCA)
   ncsu_assert(associated(cv%i))

!  call print_i(cv%i, lun)
   call print_pca(cv%i, lun)

end subroutine print_details

!=============================================================================

subroutine colvar_cleanup(cv)

   NCSU_USE_AFAILED

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   ncsu_assert(cv%type.eq.COLVAR_PCA)

   priv => get_priv(cv)
   ncsu_assert(associated(priv))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd) 

   call del_priv(cv)

end subroutine colvar_cleanup

!=============================================================================

end module ncsu_cv_PCA
