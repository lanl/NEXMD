#include "copyright.h"
#  define _REAL_ double precision

module pb_lsolver

   implicit none

#  include "pb_constants.h"

   integer l_xm, l_ym, l_zm, l_xmym, l_xmymzm
   _REAL_ l_fmiccg,l_wsor
   _REAL_ l_norm, l_inorm, l_accept, l_epsout, l_pbkappa, l_h
   integer l_itn,l_maxitn,l_bcopt
   integer mg_nlevel,ncyc_before,ncyc_after

!  All
   _REAL_,allocatable :: l_am1(:), l_am2(:), l_am3(:), l_ad(:)
   _REAL_,allocatable :: l_am4(:), l_am5(:), l_am6(:)
   _REAL_,allocatable :: l_bv(:)
   _REAL_,allocatable :: l_pv(:), l_tv(:), l_zv(:), l_rd(:)
!  MG
   integer,allocatable ::  mg_index(:), mg_index_ext(:),mg_x_idx(:),mg_size(:,:)
   _REAL_,allocatable ::  mg_onorm(:)
   _REAL_,allocatable ::  l_rv(:), l_iv(:), l_bz(:), l_xv(:)

contains

!===========================================================================

subroutine init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)

   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_wsor,p_h

   l_xm = nx
   l_ym = ny
   l_zm = nz
   l_xmym = nxny
   l_xmymzm = nxnynz
   l_maxitn = p_maxitn
   l_bcopt = p_bcopt
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

   select case (solvopt)
   case (1)
      allocate(l_ad(1:l_xmymzm+l_xmym),l_am1(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_am2(1-l_xmym:l_xmymzm+l_xmym), l_am3(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_rd(1-l_xmym:l_xmymzm),l_bv(1-l_xmym:l_xmymzm))
      allocate(l_tv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_pv(1-l_xmym:l_xmymzm+l_xmym))
   case (3)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_pv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1:l_xmymzm))
      if ( l_bcopt == 10 ) then
         allocate(l_am4(1-l_xmym:l_xmymzm),l_am5(1-l_xmym:l_xmymzm),l_am6(1-l_xmym:l_xmymzm))
      end if
   case (4)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_zv(1:l_xmymzm))
   case (5)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_am4(1-l_xmym:l_xmymzm),l_am5(1-l_xmym:l_xmymzm),l_am6(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_pv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1:l_xmymzm))
   case (2)
      allocate ( mg_index(1:mg_nlevel+1), mg_index_ext(1:mg_nlevel+1))
      allocate ( mg_x_idx(1:mg_nlevel+1),mg_size(1:3,1:mg_nlevel) )
      allocate ( mg_onorm(1:mg_nlevel) )

      mg_index_ext(1) = 1
      mg_index(1) = 1
      mg_x_idx(1) = 1
      mg_size(1,1) = l_xm
      mg_size(2,1) = l_ym
      mg_size(3,1) = l_zm
      m = l_xmymzm
      l = m + l_xmym
      n = l + l_xmym
      allocate( l_zv(1:m) )
      do i = 2, mg_nlevel
         mg_index_ext(i) = 1 + l
         mg_index(i) = 1 + m
         mg_x_idx(i) = 1 + n
         mg_size(1:3,i) = ( mg_size(1:3,i-1) - 1 ) / 2
         m = m + mg_size(1,i) * mg_size(2,i) * mg_size(3,i)
         l = l + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + mg_size(1,i) * mg_size(2,i)
         n = n + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + 2 * mg_size(1,i) * mg_size(2,i)
      end do
      mg_index_ext(i) = 1 + l
      mg_index(i) = 1 + m
      mg_x_idx(i) = 1 + n
      allocate ( l_ad(1:m), l_bv(1:m), l_rv(1:m), l_iv(1:m), l_bz(1:m) )
      allocate ( l_am1(1:l), l_am2(1:l), l_am3(1:l) )
      allocate ( l_xv(1:n) )
   end select

end subroutine allocate_array

!===========================================================================

subroutine deallocate_array(solvopt)

   implicit none
   integer solvopt

   select case (solvopt)
   case (1)
      deallocate(l_ad,l_am1)
      deallocate(l_am2,l_am3)
      deallocate(l_rd,l_bv)
      deallocate(l_tv,l_zv)
      deallocate(l_pv)
   case (3)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_bv,l_pv,l_zv)
      if ( l_bcopt == 10 ) then
         deallocate(l_am4,l_am5,l_am6)
      end if
   case (4)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_bv,l_zv)
   case (5)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_am4,l_am5,l_am6)
      deallocate(l_bv,l_pv,l_zv)
   case (2)
      deallocate( mg_index, mg_index_ext )
      deallocate( mg_x_idx,mg_size )
      deallocate( mg_onorm )
      deallocate( l_zv )
      deallocate( l_ad, l_bv, l_rv, l_iv, l_bz )
      deallocate( l_am1, l_am2, l_am3 )
      deallocate( l_xv )
   end select

end subroutine deallocate_array

!==============================================================================

!===========================================================================

subroutine init_array( solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs )
   
   implicit none

   integer solvopt
   _REAL_ epsx(*),epsy(*),epsz(*)
   _REAL_ p_bv(1:l_xmymzm),p_iv(1:l_xmymzm)
   _REAL_ p_xs(1-l_xmym:l_xmymzm+l_xmym)

   integer lxmym,l,m,n,i
   _REAL_,allocatable :: lepsx(:), lepsy(:), lepsz(:) 
   _REAL_ lfactor

   select case (solvopt)
   case (1)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      l_am1(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am2(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am3(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
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
   case (3)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      if ( l_bcopt == 10 ) then
         call pb_setupper2(l_am1(1),l_am2(1),l_am3(1), &
                           l_am4(1),l_am5(1),l_am6(1))
      else
         call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      end if
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
   case (4)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
   case (5)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      !call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      if ( l_bcopt == 10 ) call pb_setupper2(l_am1(1),l_am2(1),l_am3(1), &
                                           l_am4(1),l_am5(1),l_am6(1))
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
   case (2)
      l_ad  = ZERO
      l_am1 = ZERO
      l_am2 = ZERO
      l_am3 = ZERO
      l_iv  = ZERO
      l_bv  = ZERO
      l_rv  = ZERO
      l_zv  = ZERO
      l_xv  = ZERO
      l_xv(1+l_xmym:l_xmymzm+l_xmym) = p_xs(1:l_xmymzm)
      l_bv(1:l_xmymzm) = p_bv(1:l_xmymzm)

      m = 0
      do i = 1, mg_nlevel
         m = m + (mg_size(1,i)+1) * (mg_size(2,i)+1) * (mg_size(3,i)+1)
      end do
      allocate ( lepsx(1:m), lepsy(1:m), lepsz(1:m) )
      call feedepsintoam(l_xm, l_ym, l_zm, lepsx(1:l_xmymzm),  &
                         lepsy(1:l_xmymzm), lepsz(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      ! This is not accurate but I haven't figured out how to solve this problem.

      lfactor = l_epsout*(l_h*l_pbkappa)**2
      l_iv(1:l_xmymzm) = p_iv(1:l_xmymzm)

      i = 1
      m = mg_index(i)
      n = mg_index_ext(i)
      lxmym = mg_size(1,i)*mg_size(2,i)
      call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                     l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                     l_ad(m), l_bz(m), &
                     mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)

      do i = 2, mg_nlevel
         l = mg_index(i-1)
         m = mg_index(i)
         n = mg_index_ext(i)
         call restrict_eps_map(lepsx(l),lepsy(l),lepsz(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1), &
                               lepsx(m),lepsy(m),lepsz(m),mg_size(1,i),mg_size(2,i),mg_size(3,i))
         call restrict_iv(l_iv(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1),&
                          l_iv(m),mg_size(1,i),mg_size(2,i),mg_size(3,i) )
         lfactor = lfactor * 4
         lxmym = mg_size(1,i)*mg_size(2,i)

         call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                        l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                        l_ad(m), l_bz(m), &
                        mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)
      end do
      deallocate(lepsx,lepsy,lepsz)
   end select

contains

subroutine feedepsintoam(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ am1(1:xm,1:ym,1:zm)
    _REAL_ am2(1:xm,1:ym,1:zm)
    _REAL_ am3(1:xm,1:ym,1:zm)
    _REAL_ eps1(0:xm,1:ym,1:zm)
    _REAL_ eps2(1:xm,0:ym,1:zm)
    _REAL_ eps3(1:xm,1:ym,0:zm)
    am1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
    am2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
    am3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)
end subroutine feedepsintoam

subroutine feedepsintoad(xm, ym, zm, ad, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ ad(1:xm,1:ym,1:zm)
    _REAL_ eps1(0:xm,1:ym,1:zm)
    _REAL_ eps2(1:xm,0:ym,1:zm)
    _REAL_ eps3(1:xm,1:ym,0:zm)
    ad(1:xm,1:ym,1:zm) =                    eps1(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)
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
      l_tv(i) = xs(i) + l_am1(i     )*l_tv(i+1   ) &
                    + l_am2(i     )*l_tv(i+l_xm  ) &
                    + l_am3(i     )*l_tv(i+l_xmym)
   end do
   do i = 1, l_xmymzm
      l_zv(i) = xs(i) + l_ad (i     )*l_tv(i     ) &
                    + l_am1(i-1   )*l_zv(i-1   ) &
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

      l_tv(i) = l_pv(i) + l_am1(i     )*l_tv(i+1   ) &
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
      l_tv(i)  = xs(i) + l_am1(i)*l_tv(i+1   ) &
                     + l_am2(i)*l_tv(i+l_xm  ) &
                     + l_am3(i)*l_tv(i+l_xmym)

      phi(i) = l_tv(i)*l_rd(i)
   end do

end subroutine pb_iccg

!===========================================================================

subroutine pb_sor(phi,xs)

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer ii
   _REAL_ wsor1

!  calculate initial norm

   l_inorm = sum( abs( l_bv(1:l_xmymzm) ) )
   wsor1 = ONE - l_wsor
   uconvg = .true.
   l_itn = 0

   l_zv(1:l_xmymzm) = ONE/l_ad(1:l_xmymzm)

   do while ( uconvg )

      do ii = 1, l_xmymzm
         xs(ii) = wsor1*xs(ii) + l_wsor * (l_am1(ii-1     ) * xs(ii-1   ) + &
                                           l_am1(ii       ) * xs(ii+1   ) + &
                                           l_am2(ii-l_xm  ) * xs(ii-l_xm  ) + &
                                           l_am2(ii       ) * xs(ii+l_xm  ) + &
                                           l_am3(ii-l_xmym) * xs(ii-l_xmym) + &
                                           l_am3(ii       ) * xs(ii+l_xmym) + &
                                           l_bv(ii        )  ) * l_zv(ii)
      end do

      l_itn = l_itn + 1

      do ii = 1,l_xmymzm
         phi(ii) = l_am1(ii-1) * xs(ii-1) + l_am1(ii)* xs(ii+1) + &
                   l_am2(ii-l_xm) * xs(ii-l_xm) + l_am2(ii)*xs(ii+l_xm) + &
                   l_am3(ii-l_xmym) * xs(ii-l_xmym) + l_am3(ii)*xs(ii+l_xmym) + &
                   l_bv(ii) - l_ad(ii)* xs(ii)
      end do
      l_norm = sum(abs(phi(1:l_xmymzm)))

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then
         uconvg = .false.
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
         end if
      end if

   end do

   phi(1:l_xmymzm) = xs(1:l_xmymzm)

end subroutine pb_sor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_pcg
!
! CG core routine for linearized FDPB equation.  
!
! Authors:
! Ray Luo, Jun Wang
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine pb_pcg ( phi, xs )
! Passed variables
      _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)
!
! Local variables
!
      logical uconvg
      integer i, j, k, ii, jj
      _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

!
! begin code
!
!
! initialization
!
!   compute ||b||
!
      l_inorm = sum(ABS(l_bv(1:l_xmymzm)))
!     write(6, *)  'itn & norm ', 0, inorm
!
!   compute b - A * x(0) and save it in r(0)
!   p(0) = r(0)
!
! iteration 0:
!   compute <r(0),r(0)>
!
      l_itn = 0
      bdotb1 = ZERO
! WJ
      if ( l_bcopt == 10 ) then
         do j = 1, l_ym; do k = 1, l_zm
            ii = 1+(j-1)*l_xm+(k-1)*l_xmym
            jj = ii + l_xm - 1
            l_bv(ii)=l_bv(ii)+l_am4(ii)*xs(jj)
            l_bv(jj)=l_bv(jj)+l_am4(ii)*xs(ii)
         end do; end do
         do i = 1, l_xm; do k = 1, l_zm
            ii = i+(k-1)*l_xmym
            jj = ii + l_xmym - l_xm
            l_bv(ii)=l_bv(ii)+l_am5(ii)*xs(jj)
            l_bv(jj)=l_bv(jj)+l_am5(ii)*xs(ii)
         end do; end do
         do i = 1, l_xm; do j = 1, l_ym
            ii = i+(j-1)*l_xm
            jj = ii + l_xmymzm - l_xmym
            l_bv(ii)=l_bv(ii)+l_am6(ii)*xs(jj)
            l_bv(jj)=l_bv(jj)+l_am6(ii)*xs(ii)
         end do; end do
      end if
!
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
! WJ
         l_zv = ZERO
         if ( l_bcopt == 10 ) then
            do j = 1, l_ym; do k = 1, l_zm
               ii = 1+(j-1)*l_xm+(k-1)*l_xmym
               jj = ii + l_xm - 1
               l_zv(ii)=l_zv(ii)-l_am4(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am4(ii)*l_pv(ii)
            end do; end do
            do i = 1, l_xm; do k = 1, l_zm
               ii = i+(k-1)*l_xmym
               jj = ii + l_xmym - l_xm
               l_zv(ii)=l_zv(ii)-l_am5(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am5(ii)*l_pv(ii)
            end do; end do
            do i = 1, l_xm; do j = 1, l_ym
               ii = i+(j-1)*l_xm
               jj = ii + l_xmymzm - l_xmym
               l_zv(ii)=l_zv(ii)-l_am6(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am6(ii)*l_pv(ii)
            end do; end do
         end if
!
         do i = 1,l_xmymzm
            l_zv(i) =   l_ad(i)      *l_pv(i)      &
                    - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                    - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                    - l_am1(i-1   )*l_pv(i-1   ) &
                    - l_am1(i     )*l_pv(i+1   ) &
                    - l_am2(i     )*l_pv(i+l_xm  ) &
                    - l_am3(i     )*l_pv(i+l_xmym) &
! WJ
                    + l_zv(i)
!
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
            xs(i)       = xs(i)       + alpha*l_pv(i)
            l_bv(i)       = l_bv(i)       - alpha*l_zv(i)
            l_norm        = l_norm        +   ABS(l_bv(i))
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
            bdotb2      = bdotb2      + l_bv(i)*l_bv(i)
         end do
!        write(6, *)  'l_itn & l_norm ',l_itn, l_norm
!
!   check convergence
!
         if ( l_itn .ge. l_maxitn .or. l_norm .le. l_accept*l_inorm ) then

            uconvg = .false.
            if ( l_itn .ge. l_maxitn ) then
               write(6, *) 'PBMD WARNING: CG maxitn exceeded!'
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
!
!
      return
      end subroutine pb_pcg

!===========================================================================

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
   phi(1:l_xmymzm) = xs(1:l_xmymzm)

end subroutine pb_cg

!===========================================================================

subroutine pb_mg(phi,xs)

   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables
   logical uconvg
   integer j,lxly,lxlylz,p1,p2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

   l_itn = 0
   uconvg = .true.

   do while ( uconvg )
      mg_onorm = 9.9d99
      do j = 1, mg_nlevel-1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_before,l_accept,mg_onorm(j) )
         call restrict_bv(l_rv(mg_index(j)), mg_size(1,j), mg_size(2,j), mg_size(3,j), &
                          l_bv(mg_index(j+1)), mg_size(1,j+1),mg_size(2,j+1),mg_size(3,j+1) )
         l_xv(mg_x_idx(j+1):mg_x_idx(j+2)-1) = ZERO
      end do
      lxly = mg_size(1,j)*mg_size(2,j)
      lxlylz = lxly*mg_size(3,j)
      call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                 l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                 l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                 lxly,lxlylz,-1,l_accept,mg_onorm(j) )
      do j = mg_nlevel-1, 1, -1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         p1 = mg_x_idx(j+1)+mg_size(1,j+1)*mg_size(2,j+1)
         p2 = mg_x_idx(j)+lxly
         call interpolate(l_xv(p1), mg_size(1,j+1), mg_size(2,j+1), mg_size(3,j+1),&
                          l_xv(p2), mg_size(1,j) ,  mg_size(2,j), mg_size(3,j) , &
                          l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)), &
                          l_am3(mg_index_ext(j)),l_bz(mg_index(j)),l_epsout)
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)), l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_after,l_accept,mg_onorm(j) )
      end do
      l_itn = l_itn + 1
      l_norm = sum(abs(l_rv(1:l_xmymzm)))
!  
!  convergence
!  
      if ( l_itn .ge. l_maxitn .or. l_norm .le. l_inorm*l_accept ) then
         uconvg = .false.
!        write(6, *)  '   Multigrid itn & norm', litn, lnorm
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: Multigrid maxitn exceeded!'
         endif
      end if
   end do

   xs (1:l_xmymzm) = l_xv(l_xmym+1:l_xmym+l_xmymzm)
   phi(1:l_xmymzm) = xs (1:l_xmymzm)

end subroutine pb_mg

!===========================================================================

subroutine relax(xs,nx,ny,nz,lam1,lam2,lam3,lad,lbv,lrv,nxny,nxnynz,ncyc,accept,onorm)
 
   implicit none

   integer nx,ny,nz,nxny,nxnynz,ncyc
   _REAL_ xs(1-nxny:nxnynz+nxny), lam1(1-nxny:nxnynz), lam2(1-nxny:nxnynz), lam3(1-nxny:nxnynz)
   _REAL_ lad(1:nxnynz), lbv(1:nxnynz), lrv(1:nxnynz)
   _REAL_ accept, onorm

   logical luconvg
   integer ii,litn,itmax
   _REAL_ wsor, wsor1, linorm, lnorm
   integer itn_checknorm

   if (ncyc>0) then
      itn_checknorm = ncyc
      wsor = 1.0d0
      wsor1 = ONE - wsor
   else
      itn_checknorm = 10
      wsor = 1.9d0
      wsor1 = ONE - wsor
   end if

   linorm = sum(abs(lbv(1:nxnynz)))
!  write(6, *)  '      relax itn & norm ', litn, linorm!, onorm
   l_zv(1:nxnynz) = ONE/lad(1:nxnynz)

   litn = 0
!  lnorm = 0
   itmax = 1000
   luconvg = .true.

   do while ( luconvg )

!     start the sor iteration ...
 
      do ii = 1, nxnynz, 1
         xs(ii) = wsor1*xs(ii) + wsor * (lam1(ii-1   ) * xs(ii-1   ) + &
                                         lam1(ii     ) * xs(ii+1   ) + &
                                         lam2(ii-nx  ) * xs(ii-nx  ) + &
                                         lam2(ii     ) * xs(ii+nx  ) + &
                                         lam3(ii-nxny) * xs(ii-nxny) + &
                                         lam3(ii     ) * xs(ii+nxny) + &
                                          lbv(ii     )             ) * l_zv(ii)
      end do

      litn = litn + 1

      if ( mod(litn,itn_checknorm) == 0 ) then
         do ii = 1,nxnynz
            lrv(ii) = lam1(ii-1) * xs(ii-1) + lam1(ii)* xs(ii+1) + &
                     lam2(ii-nx) * xs(ii-nx) + lam2(ii)*xs(ii+nx) + &
                     lam3(ii-nxny) * xs(ii-nxny) + lam3(ii)*xs(ii+nxny) + &
                     lbv(ii) - lad(ii)* xs(ii)
         end do
         lnorm = sum(abs(lrv(1:nxnynz)))

!        write(6, *)  '      relax itn & norm ', litn, lnorm!, onorm
!  
!        check convergence
!  
         if ( litn .ge. itmax .or. ( ncyc .gt. 0 .and. (litn .ge. ncyc .and. lnorm < onorm ) ) &
              .or. lnorm .le. accept*linorm ) then

            luconvg = .false.
! WJ
            if ( ncyc .gt. 0 .and. litn .ge. ncyc .and. lnorm > onorm ) then
               write(6,*) "MG FAILED: ncyc, itn, norm, onorm", ncyc, litn, lnorm, onorm
               stop
            end if
!
            if ( ncyc .gt. 0 ) onorm = lnorm 
            if ( litn .ge. itmax ) then
               write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
            endif
         end if
      end if

   end do

end subroutine relax

!===========================================================================
      
subroutine set_am_ad( epsx,epsy,epsz,iv,lam1,lam2,lam3,lad,lbz, &
                      xn,yn,zn,lfactor,epsout )

    implicit none

   _REAL_  lfactor, epsout

   integer xn,yn,zn
   _REAL_  epsx(xn,yn,zn), epsy(xn,yn,zn), epsz(xn,yn,zn), iv(xn,yn,zn)
   _REAL_  lam1(xn,yn,zn),lam2(xn,yn,zn),lam3(xn,yn,zn)
   _REAL_  lad(xn,yn,zn),lbz(xn,yn,zn)

   integer  i,j,k,i1,j1,k1
   _REAL_ lam1t,lam2t,lam3t

   lam1(1:xn,1:yn,1:zn) = epsx(1:xn,1:yn,1:zn)
   lam2(1:xn,1:yn,1:zn) = epsy(1:xn,1:yn,1:zn)
   lam3(1:xn,1:yn,1:zn) = epsz(1:xn,1:yn,1:zn)

   do k = 1, zn
      k1 = k-1
      do j = 1, yn
         j1 = j-1
         do i = 1, xn
            i1 = i-1
            lam1t = epsout
            if ( i1 /= 0 ) lam1t = lam1(i1,j,k)
            lam2t = epsout
            if ( j1 /= 0 ) lam2t = lam2(i,j1,k)
            lam3t = epsout
            if ( k1 /= 0 ) lam3t = lam3(i,j,k1)
            lad(i,j,k) = lam1t + lam1(i,j,k) + lam2t + lam2(i,j,k) &
                       + lam3t + lam3(i,j,k) 
         end do
      end do
   end do

   lbz = lfactor*iv
   lad = lad + lbz

   do k = 1, zn; do j = 1, yn
      lam1(xn,j,k) = ZERO
   end do; end do
   do k = 1, zn; do i = 1, xn
      lam2(i,yn,k) = ZERO
   end do; end do
   do j = 1, yn; do i = 1, xn
      lam3(i,j,zn) = ZERO
   end do; end do

end subroutine set_am_ad

!===========================================================================

subroutine restrict_eps_map( epsx, epsy, epsz, xn, yn, zn, epsxr, epsyr, epszr, xnr, ynr, znr )

   implicit none
 
   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ epsx(xn,yn,zn),epsxr(xnr,ynr,znr)
   _REAL_ epsy(xn,yn,zn),epsyr(xnr,ynr,znr)
   _REAL_ epsz(xn,yn,zn),epszr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*),"Restriction of eps map Failed because of incorrect dimension"
      stop
   end if
   do k = 1, znr
      k2 = 2*k
      do j = 1, ynr
         j2 = 2*j
         do i = 1, xnr
            i2 = 2*i
            epsxr(i,j,k) = hmav(epsx(i2  ,j2  ,k2  ),epsx(i2+1,j2  ,k2  ))/4.d0 + &
                          (hmav(epsx(i2  ,j2-1,k2  ),epsx(i2+1,j2-1,k2  )) + &
                           hmav(epsx(i2  ,j2+1,k2  ),epsx(i2+1,j2+1,k2  )) + &
                           hmav(epsx(i2  ,j2  ,k2-1),epsx(i2+1,j2  ,k2-1)) + &
                           hmav(epsx(i2  ,j2  ,k2+1),epsx(i2+1,j2  ,k2+1)))/8.d0 + &
                          (hmav(epsx(i2  ,j2-1,k2-1),epsx(i2+1,j2-1,k2-1)) + &
                           hmav(epsx(i2  ,j2+1,k2-1),epsx(i2+1,j2+1,k2-1)) + &
                           hmav(epsx(i2  ,j2-1,k2+1),epsx(i2+1,j2-1,k2+1)) + &
                           hmav(epsx(i2  ,j2+1,k2+1),epsx(i2+1,j2+1,k2+1)))/16.d0 

            epsyr(i,j,k) = hmav(epsy(i2  ,j2  ,k2  ),epsy(i2  ,j2+1,k2  ))/4.d0 + &
                          (hmav(epsy(i2-1,j2  ,k2  ),epsy(i2-1,j2+1,k2  )) + &
                           hmav(epsy(i2+1,j2  ,k2  ),epsy(i2+1,j2+1,k2  )) + &
                           hmav(epsy(i2  ,j2  ,k2-1),epsy(i2  ,j2+1,k2-1)) + &
                           hmav(epsy(i2  ,j2  ,k2+1),epsy(i2  ,j2+1,k2+1)))/8.d0 + &
                          (hmav(epsy(i2-1,j2  ,k2-1),epsy(i2-1,j2+1,k2-1)) + &
                           hmav(epsy(i2+1,j2  ,k2-1),epsy(i2+1,j2+1,k2-1)) + &
                           hmav(epsy(i2-1,j2  ,k2+1),epsy(i2-1,j2+1,k2+1)) + &
                           hmav(epsy(i2+1,j2  ,k2+1),epsy(i2+1,j2+1,k2+1)))/16.d0 

            epszr(i,j,k) = hmav(epsz(i2  ,j2  ,k2  ),epsz(i2  ,j2  ,k2+1))/4.d0 + &
                          (hmav(epsz(i2  ,j2-1,k2  ),epsz(i2  ,j2-1,k2+1)) + &
                           hmav(epsz(i2  ,j2+1,k2  ),epsz(i2  ,j2+1,k2+1)) + &
                           hmav(epsz(i2-1,j2  ,k2  ),epsz(i2-1,j2  ,k2+1)) + &
                           hmav(epsz(i2+1,j2  ,k2  ),epsz(i2+1,j2  ,k2+1)))/8.d0 + &
                          (hmav(epsz(i2-1,j2-1,k2  ),epsz(i2-1,j2-1,k2+1)) + &
                           hmav(epsz(i2-1,j2+1,k2  ),epsz(i2-1,j2+1,k2+1)) + &
                           hmav(epsz(i2+1,j2-1,k2  ),epsz(i2+1,j2-1,k2+1)) + &
                           hmav(epsz(i2+1,j2+1,k2  ),epsz(i2+1,j2+1,k2+1)))/16.d0 

         end do
      end do
   end do

contains 

function hmav(a,b)

   _REAL_ hmav, a, b

   hmav = 2.d0*a*b/(a+b)

end function hmav

function hmav4(a,b,c,d)

   _REAL_ hmav4, a,b,c,d

   hmav4 = 4.d0/(1.d0/a+1.d0/b+1.d0/c+1.d0/d)

end function hmav4

end subroutine restrict_eps_map

!===========================================================================

subroutine restrict_iv(ivf, xn, yn, zn, ivr, xnr, ynr, znr)

   implicit none

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ ivf(xn,yn,zn),ivr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*),"Restriction of iv Failed because of incorrect dimension"
      stop
   end if

   do k = 1, znr
      k2 = k * 2
      do j = 1, ynr
         j2 = j * 2
         do i = 1, xnr
            i2 = i * 2
            ivr(i,j,k) =       ( ivf(i2-1,j2-1,k2-1) + TWO * ivf(i2  ,j2-1,k2-1) + ivf(i2+1,j2-1,k2-1) ) +  &
                        TWO  * ( ivf(i2-1,j2  ,k2-1) + TWO * ivf(i2  ,j2  ,k2-1) + ivf(i2+1,j2  ,k2-1) ) +  &
                               ( ivf(i2-1,j2+1,k2-1) + TWO * ivf(i2  ,j2+1,k2-1) + ivf(i2+1,j2+1,k2-1) ) +  &
                        TWO  * ( ivf(i2-1,j2-1,k2  ) + TWO * ivf(i2  ,j2-1,k2  ) + ivf(i2+1,j2-1,k2  ) ) +  &
                        FOUR * ( ivf(i2-1,j2  ,k2  ) + TWO * ivf(i2  ,j2  ,k2  ) + ivf(i2+1,j2  ,k2  ) ) +  &
                        TWO  * ( ivf(i2-1,j2+1,k2  ) + TWO * ivf(i2  ,j2+1,k2  ) + ivf(i2+1,j2+1,k2  ) ) +  &
                               ( ivf(i2-1,j2-1,k2+1) + TWO * ivf(i2  ,j2-1,k2+1) + ivf(i2+1,j2-1,k2+1) ) +  &
                        TWO  * ( ivf(i2-1,j2  ,k2+1) + TWO * ivf(i2  ,j2  ,k2+1) + ivf(i2+1,j2  ,k2+1) ) +  &
                               ( ivf(i2-1,j2+1,k2+1) + TWO * ivf(i2  ,j2+1,k2+1) + ivf(i2+1,j2+1,k2+1) )
            ivr(i,j,k) = ivr(i,j,k) / 64.d0
         end do
      end do
   end do

end subroutine restrict_iv

!===========================================================================

subroutine restrict_bv(bvf, xn, yn, zn, bvr, xnr, ynr, znr)
 
   implicit none 

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ bvf(xn,yn,zn),bvr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*),"Restriction of bv Failed because of incorrect dimension "
      stop
   end if

   do k = 1, znr
      k2 = k * 2
      do j = 1, ynr
         j2 = j * 2
         do i = 1, xnr
            i2 = i * 2
            bvr(i,j,k) =       ( bvf(i2-1,j2-1,k2-1) + TWO * bvf(i2  ,j2-1,k2-1) + bvf(i2+1,j2-1,k2-1) ) +  &
                        TWO  * ( bvf(i2-1,j2  ,k2-1) + TWO * bvf(i2  ,j2  ,k2-1) + bvf(i2+1,j2  ,k2-1) ) +  &
                               ( bvf(i2-1,j2+1,k2-1) + TWO * bvf(i2  ,j2+1,k2-1) + bvf(i2+1,j2+1,k2-1) ) +  &
                        TWO  * ( bvf(i2-1,j2-1,k2  ) + TWO * bvf(i2  ,j2-1,k2  ) + bvf(i2+1,j2-1,k2  ) ) +  &
                        FOUR * ( bvf(i2-1,j2  ,k2  ) + TWO * bvf(i2  ,j2  ,k2  ) + bvf(i2+1,j2  ,k2  ) ) +  &
                        TWO  * ( bvf(i2-1,j2+1,k2  ) + TWO * bvf(i2  ,j2+1,k2  ) + bvf(i2+1,j2+1,k2  ) ) +  &
                               ( bvf(i2-1,j2-1,k2+1) + TWO * bvf(i2  ,j2-1,k2+1) + bvf(i2+1,j2-1,k2+1) ) +  &
                        TWO  * ( bvf(i2-1,j2  ,k2+1) + TWO * bvf(i2  ,j2  ,k2+1) + bvf(i2+1,j2  ,k2+1) ) +  &
                               ( bvf(i2-1,j2+1,k2+1) + TWO * bvf(i2  ,j2+1,k2+1) + bvf(i2+1,j2+1,k2+1) )
            bvr(i,j,k) = bvr(i,j,k) / 16.d0
         end do
      end do
   end do

end subroutine restrict_bv

!===========================================================================

   subroutine interpolate(v, xn, yn ,zn, vi, xni, yni, zni, lam1, lam2, lam3, lbz, epsout )


   integer xn,yn,zn,xni,yni,zni
   _REAL_ v(1:xn*yn*zn),vi(1:xni*yni*zni), epsout
   _REAL_ lam1(1-xni*yni:xni*yni*zni),lam2(1-xni*yni:xni*yni*zni),lam3(1-xni*yni:xni*yni*zni)
   _REAL_ lbz(1:xni*yni*zni)

   integer i,j,k,xniyni,xniynizni,ii,ii2

   if ( xn*2+1 /= xni .or. yn*2+1 /= yni .or. zn*2+1 /= zni ) then
      write(6,*), "Interpolation Failed because of incorrect dimension"
      stop
   end if

   xniyni = xni*yni
   xniynizni = xni*yni*zni

   lam1(1-xniyni:0) = epsout
   lam2(1-xniyni:0) = epsout
   lam3(1-xniyni:0) = epsout
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = epsout
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = epsout
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = epsout
   end do; end do

   ii = 0
   do k = 1, zn; do j=1 ,yn; do i=1, xn
      ii = ii + 1
      ii2 = (k*2-1)*xniyni + (j*2-1)*xni + i*2
      vi(ii2) = vi(ii2) + v(ii)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,-1,lam2,xni,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,+1,lam2,xni,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,-xni,lam1,1,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,+xni,lam1,1,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,-xniyni,lam2,xni,lam1,1)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,+xniyni,lam2,xni,lam1,1)
   end do; end do; end do

   lam1(1-xniyni:0) = ZERO
   lam2(1-xniyni:0) = ZERO
   lam3(1-xniyni:0) = ZERO
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = ZERO
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = ZERO
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = ZERO
   end do; end do

contains

subroutine ipl_chain(vi,xnyn,xnynzn,l,v,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)

   integer xnyn,xnynzn,l,shift_1,shift_2,shift_3
   _REAL_ vi(1:xnynzn),v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l1
   _REAL_ v1

   v1 = ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1)
   l1 = l + shift_1 
   vi(l1) = vi(l1) + v1

   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,-shift_2,am_3,shift_3)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,-shift_3,am_2,shift_2)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,shift_3,am_2,shift_2)

end subroutine ipl_chain

subroutine ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)

   integer xnyn,xnynzn,l1,shift_1,shift_2,shift_3
   _REAL_ vi(1:xnynzn),v1,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l2
   _REAL_ v2

   v2 = ipl_comp2(v1,l1,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2)
   l2 = l1 + shift_2
   vi(l2) = vi(l2) + v2
   vi(l2-shift_3)= vi(l2-shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,-shift_3)
   vi(l2+shift_3)= vi(l2+shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,+shift_3)

end subroutine ipl_chain2

function ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1) 

   integer l,xnyn,xnynzn,shift_1
   _REAL_ ipl_comp1,v,am_1(1-xnyn:xnynzn),lbz(1:xnynzn)
   if ( shift_1 < 0 ) then
      ipl_comp1 = v * am_1(l+shift_1) / ( lbz(l+shift_1) + am_1(l+2*shift_1) + am_1(l+shift_1) )
   else
      ipl_comp1 = v * am_1(l) / ( lbz(l+shift_1) + am_1(l) + am_1(l+shift_1) )
   end if
   return

end function ipl_comp1

function ipl_comp2(v,l,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2) 

   integer l,xnyn,xnynzn,shift_1,shift_2
   _REAL_ ipl_comp2,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_2) + am_1(l+shift_2-abs(shift_1)) + lbz(l+shift_2)
   if ( shift_2 < 0 ) then
      ipl_comp2 = v * am_2(l+shift_2) / ( am_2(l+2*shift_2) + am_2(l+shift_2) + lad)
   else
      ipl_comp2 = v * am_2(l) / ( am_2(l) + am_2(l+shift_2) + lad )
   end if
   return

end function ipl_comp2
            
function ipl_comp3(v,l,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,shift_3) 

   integer l,xnyn,xnynzn,shift_1,shift_2,shift_3
   _REAL_ ipl_comp3,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_3) + am_1(l+shift_3-abs(shift_1)) + am_2(l+shift_3) + &
         am_2(l+shift_3-abs(shift_2)) + lbz(l+shift_3)
   if ( shift_3 < 0 ) then
      ipl_comp3 = v * am_3(l+shift_3) / ( am_3(l+2*shift_3) + am_3(l+shift_3) + lad)
   else
      ipl_comp3 = v * am_3(l) / ( am_3(l) + am_3(l+shift_3) + lad )
   end if
   return

end function ipl_comp3

end subroutine interpolate

!===========================================================================

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

subroutine pb_setupper2( l_am1, l_am2, l_am3, l_am4, l_am5, l_am6 )

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)
   _REAL_ l_am4(l_xm,l_ym,l_zm), l_am5(l_xm,l_ym,l_zm), l_am6(l_xm,l_ym,l_zm)
   integer i, j, k

   do j = 1, l_ym; do k = 1, l_zm
      if ( l_bcopt == 10 ) l_am4(1,j,k) = l_am1(l_xm,j,k)
      l_am1(l_xm,j,k) = ZERO
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      if ( l_bcopt == 10 ) l_am5(i,1,k) = l_am2(i,l_ym,k)
      l_am2(i,l_ym,k) = ZERO
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      if ( l_bcopt == 10 ) l_am6(i,j,1) = l_am3(i,j,l_zm)
      l_am3(i,j,l_zm) = ZERO
   end do; end do
end subroutine pb_setupper2
end module pb_lsolver
