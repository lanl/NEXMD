#include "copyright.h"
#  define _REAL_ double precision

module aug_solver

   implicit none

#  include "pb_constants.h"

   integer l_xm, l_ym, l_zm, l_xmym, l_xmymzm
   _REAL_ l_fmiccg,l_wsor
   _REAL_ l_norm, l_inorm, l_accept, l_epsout, l_pbkappa, l_h
   integer l_itn,l_maxitn
   integer mg_nlevel,ncyc_before,ncyc_after

!  All
   _REAL_,allocatable :: l_am1(:), l_am2(:), l_am3(:), l_ad(:)
   _REAL_,allocatable :: l_bv(:)
   _REAL_,allocatable :: l_pv(:), l_tv(:), l_zv(:), l_rd(:)

contains

!===========================================================================

subroutine init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)

   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_wsor,p_h

   l_xm = nx
   l_ym = ny
   l_zm = nz
   l_xmym = nxny
   l_xmymzm = nxnynz
   l_maxitn = p_maxitn
   l_accept = p_accept

!  ICCG
   l_fmiccg = p_fmiccg
!  MG
   mg_nlevel = 4
   ncyc_before = 10
   ncyc_after = 10
   l_pbkappa = p_pbkappa
   l_epsout = p_epsout
   l_h       = p_h
!  SOR
   l_wsor = p_wsor

end subroutine

!===========================================================================

subroutine allocate_array(solvopt)

   implicit none
   integer solvopt

   integer l,m,n,i

      allocate(l_ad(1:l_xmymzm+l_xmym),l_am1(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_am2(1-l_xmym:l_xmymzm+l_xmym), l_am3(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_rd(1-l_xmym:l_xmymzm),l_bv(1-l_xmym:l_xmymzm))
      allocate(l_tv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_pv(1-l_xmym:l_xmymzm+l_xmym))

end subroutine allocate_array

!===========================================================================

subroutine deallocate_array(solvopt)

   implicit none
   integer solvopt

      deallocate(l_ad,l_am1)
      deallocate(l_am2,l_am3)
      deallocate(l_rd,l_bv)
      deallocate(l_tv,l_zv)
      deallocate(l_pv)

end subroutine deallocate_array

!==============================================================================

!===========================================================================

subroutine init_array( solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs )
   
   implicit none

   integer solvopt
   _REAL_ epsx,epsy,epsz
   _REAL_ p_bv(1:l_xmymzm),p_iv(1:l_xmymzm)
   _REAL_ p_xs(1-l_xmym:l_xmymzm+l_xmym)

   integer lxmym,l,m,n,i
   _REAL_,allocatable :: lepsx(:), lepsy(:), lepsz(:) 
   _REAL_ lfactor

      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      l_am1(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am2(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am3(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO

      call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)

      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      l_ad (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_bv (1-l_xmym:0) = ZERO
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_tv = ZERO
      l_zv = ZERO
      l_rd = ONE
contains

subroutine feedepsintoam(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ am1(1:xm,1:ym,1:zm)
    _REAL_ am2(1:xm,1:ym,1:zm)
    _REAL_ am3(1:xm,1:ym,1:zm)
    _REAL_ eps1
    _REAL_ eps2
    _REAL_ eps3
    am1(1:xm,1:ym,1:zm) = eps1
    am2(1:xm,1:ym,1:zm) = eps2
    am3(1:xm,1:ym,1:zm) = eps3
end subroutine feedepsintoam

subroutine feedepsintoad(xm, ym, zm, ad, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ ad(1:xm,1:ym,1:zm)
    _REAL_ eps1
    _REAL_ eps2
    _REAL_ eps3
    ad(1:xm,1:ym,1:zm) =                    eps1
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3
end subroutine feedepsintoad

end subroutine init_array

!===========================================================================

subroutine pb_iccg(phi,xs)

!  use poisson_boltzmann, only: level
   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   ! initialization

   do i = 1, l_xmymzm
      l_rd(i) = ONE/( l_ad(i) - &
         l_am1(i-1   )*(       l_am1(i-1   )+l_fmiccg*l_am2(i-1   )+l_fmiccg*l_am3(i-1   ))*l_rd(i-1   ) - &
         l_am2(i-l_xm  )*(l_fmiccg*l_am1(i-l_xm  )+       l_am2(i-l_xm  )+l_fmiccg*l_am3(i-l_xm  ))*l_rd(i-l_xm  ) - &
         l_am3(i-l_xmym)*(l_fmiccg*l_am1(i-l_xmym)+l_fmiccg*l_am2(i-l_xmym)+       l_am3(i-l_xmym))*l_rd(i-l_xmym) )
   end do

   do i = 1, l_xmymzm
      l_ad(i) = l_ad(i)*l_rd(i)
      l_rd(i) = sqrt(l_rd(i))
      l_bv(i) = l_bv(i)*l_rd(i)
      l_am1(i-1   ) = l_am1(i-1   )*l_rd(i)*l_rd(i-1   )
      l_am2(i-l_xm  ) = l_am2(i-l_xm  )*l_rd(i)*l_rd(i-l_xm  )
      l_am3(i-l_xmym) = l_am3(i-l_xmym)*l_rd(i)*l_rd(i-l_xmym)
   end do

   l_inorm = ZERO
   do i = 1, l_xmymzm
      l_ad(i) = l_ad(i) - TWO

      l_bv(i) = l_bv(i) + l_am1(i-1   )*l_bv(i-1   ) &
                    + l_am2(i-l_xm  )*l_bv(i-l_xm  ) &
                    + l_am3(i-l_xmym)*l_bv(i-l_xmym)
      l_inorm = l_inorm + abs(l_bv(i))
   end do

   do i = l_xmymzm, 1, -1
      l_tv(i) = xs(i) + l_am1(i     )*l_tv(i+1     ) &
                      + l_am2(i     )*l_tv(i+l_xm  ) &
                      + l_am3(i     )*l_tv(i+l_xmym)
   end do
   do i = 1, l_xmymzm
      l_zv(i) = xs(i) + l_ad (i       )*l_tv(i       ) &
                      + l_am1(i-1     )*l_zv(i-1     ) &
                      + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                      + l_am3(i-l_xmym)*l_zv(i-l_xmym)
   end do
   bdotb1 = ZERO
   do i = l_xmymzm, 1, -1
      l_zv(i) = l_zv(i) + l_tv(i)
      l_bv(i) = l_bv(i) - l_zv(i)

      ! iteration 0.

      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)

      ! first step of the matrix vector multiplication, see below

      l_tv(i) = l_pv(i) + l_am1(i     )*l_tv(i+1     ) &
                        + l_am2(i     )*l_tv(i+l_xm  ) &
                        + l_am3(i     )*l_tv(i+l_xmym)
   end do

   l_itn = 0
   uconvg = .true.

   ! the main loop of iccg solver

   do while ( uconvg )

      ! second and third steps of the matrix vector multiplication

      pdotz = ZERO
      do i = 1, l_xmymzm+l_xmym
         l_zv(i) = l_pv(i) + l_ad (i     )*l_tv(i     ) &
                           + l_am1(i-1   )*l_zv(i-1   ) &
                           + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                           + l_am3(i-l_xmym)*l_zv(i-l_xmym)

         j = i - l_xmym
         l_zv(j) = l_zv(j) + l_tv(j)

         pdotz = pdotz + l_pv(j)*l_zv(j)
      end do
      alpha = bdotb1/pdotz

      l_norm = ZERO
      bdotb2 = ZERO
      l_itn = l_itn + 1
      do i = 1, l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  + abs(l_bv(i))

         bdotb2= bdotb2+ l_bv(i)*l_bv(i)
      end do

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6, *) 'PB warning in pb_miccg(): CG l_maxitn exceeded!'
         end if

      else

         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication

         do i = l_xmymzm, 1, -1
            l_pv(i) = l_bv(i) + beta*l_pv(i)

            l_tv(i) = l_pv(i) + l_am1(i)*l_tv(i+1   ) &
                          + l_am2(i)*l_tv(i+l_xm  ) &
                          + l_am3(i)*l_tv(i+l_xmym)
         end do
      end if
   end do  !  while ( uconvg ), end of the main iccg loop

   ! back scaling of the solution

   do i = l_xmymzm, 1, -1
      l_tv(i)  = xs(i) + l_am1(i)*l_tv(i+1     ) &
                       + l_am2(i)*l_tv(i+l_xm  ) &
                       + l_am3(i)*l_tv(i+l_xmym)

      phi(i) = l_tv(i)*l_rd(i)
   end do

end subroutine pb_iccg

subroutine pb_cg(phi,xs)

    implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

!  compute b - A * x(0) and save it in r(0)
!  p(0) = r(0)
!
!  iteration 0:
!  compute <r(0),r(0)>
!
   l_itn = 0
   bdotb1 = ZERO
   do i = 1,l_xmymzm
      l_bv(i)  = l_bv(i)  + l_am3(i-l_xmym)*xs(i-l_xmym) &
                      + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                      + l_am1(i-1   )*xs(i-1   ) &
                      - l_ad(i)      *xs(i     ) &
                      + l_am1(i     )*xs(i+1   ) &
                      + l_am2(i     )*xs(i+l_xm  ) &
                      + l_am3(i     )*xs(i+l_xmym)
      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)
   end do
!
! the main loop of the CG solver
!
   uconvg = .true.
   do while ( uconvg )
!
! iteration i:
!
!   compute Ap(i) = A * p(i)
!   compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
!
      pdotz = ZERO
      do i = 1,l_xmymzm
         l_zv(i) = - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                 - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                 - l_am1(i-1   )*l_pv(i-1   ) &
                 - l_am1(i     )*l_pv(i+1   ) &
                 - l_am2(i     )*l_pv(i+l_xm  ) &
                 - l_am3(i     )*l_pv(i+l_xmym) &
                 + l_ad(i)      *l_pv(i)  
         pdotz = pdotz + l_pv(i)*l_zv(i)
      end do
!
! iteration i+1:
!
      l_itn = l_itn + 1
!
!   update x(i+1) = x(i) + alpha(i) p(i)
!          r(i+1) = r(i) - alpha(i) Ap(i)
!
      alpha  = bdotb1/pdotz
      l_norm   = ZERO
      bdotb2 = ZERO
      do i = 1,l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  +  abs(l_bv(i))
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
         bdotb2 = bdotb2 + l_bv(i)*l_bv(i)
      end do
!     write(6, *)  'itn & norm ',l_itn, l_norm
!
!   check convergence
!
      if ( l_itn .ge. l_maxitn .or. l_norm .le. l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: CG l_maxitn exceeded!'
         endif

      else
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two
!
         beta   = bdotb2/bdotb1
         bdotb1 = bdotb2
!
!   update p(i+1) = r(i+1) + beta(i) p(i)
!
         l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
      endif
   enddo
!
! end of the main CG loop
!
!
   phi(1:l_xmymzm) = xs(1:l_xmymzm)


end subroutine pb_cg

subroutine pb_setupper( l_am1, l_am2, l_am3 )

   implicit none

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)

   integer i,j,k

   do j = 1, l_ym; do k = 1, l_zm
      l_am1(l_xm,j,k) = ZERO
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      l_am2(i,l_ym,k) = ZERO
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      l_am3(i,j,l_zm) = ZERO
   end do; end do

end subroutine pb_setupper

end module aug_solver
