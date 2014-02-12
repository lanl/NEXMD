#include "copyright.h"
#  define _REAL_ double precision
!#include "timer.h"

module pb_nlsolver

   implicit none

#  include "pb_constants.h"

   logical uconvg
   integer xm, ym, zm, xmym, xmymzm
   integer itn, maxitn, idamp
   integer itn_checknorm
   integer npbstep, npbgrid
   _REAL_ norm, inorm, savnorm, norm2, savnorm2
   _REAL_ factor, ktinv, factor1, maxbf
   _REAL_ wsor, lwsor, damp
   _REAL_ temp

!  MG
   integer :: mg_nlevel, ncyc_before, ncyc_after
   integer,allocatable :: mg_index(:), mg_index_ext(:), mg_x_idx(:), mg_size(:,:)
   _REAL_,allocatable :: mg_onorm(:), bz(:)

!  ICCG
   _REAL_,allocatable :: sav_am1(:), sav_am2(:), sav_am3(:), sav_ad(:)
   _REAL_,allocatable :: rd(:), zzv(:)

!  All
   _REAL_,allocatable :: am1(:), am2(:), am3(:), ad(:), ad1(:)
   _REAL_,allocatable :: bv(:), xv(:), rv(:), iv(:), xcorr(:) 
   _REAL_,allocatable :: pv(:), tv(:), zv(:)

contains
!===============================================================================
subroutine init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_npbstep,p_npbgrid, &
                      ivalence,h,pbkb,pbtemp,istrng,p_wsor,p_lwsor)

   integer nx, ny, nz, nxny, nxnynz, p_maxitn, p_npbstep, p_npbgrid
   _REAL_ ivalence, h, pbkb, pbtemp, istrng, p_wsor, p_lwsor

   uconvg = .true.   
 
   xm = nx
   ym = ny
   zm = nz
   xmym = nxny
   xmymzm = nxnynz
   mg_nlevel = 4
   maxitn = p_maxitn
   npbstep = p_npbstep
   npbgrid = p_npbgrid
   wsor = p_wsor
   lwsor = p_lwsor

   idamp = 0
   itn_checknorm = 10
   factor = TWO*istrng*h*h/ivalence
   ktinv = ONE/pbkb/pbtemp*ivalence
   factor1 = factor * ktinv
   damp = ONE
   inorm = ZERO
   norm = ZERO
   norm2 = ZERO
   savnorm = ZERO
   savnorm2 = ZERO
   maxbf = ZERO 
   temp = ZERO

end subroutine init_param
!=====================================================================================================================
subroutine allocate_array(solvopt)
   implicit none
   integer i,l,m,n
   integer solvopt

   select case ( solvopt )
   case (1)
      m = xmymzm
      l = m + xmym
      n = l + xmym
      allocate ( ad1(1:m), rv(1:m), iv(1:m),zzv(1:m), bv(1:l) )
      allocate ( am1(1:n), am2(1:n), am3(1:n), ad(1:l) )
      allocate ( sav_am1(1:n), sav_am2(1:n), sav_am3(1:n), sav_ad(1:l) )
      allocate ( rd(1:n), pv(1:n), tv(1:n), zv(1:n) )
      allocate ( xv(1:n), xcorr(1:n) )
   case (2)
      allocate ( mg_index(1:mg_nlevel+1), mg_index_ext(1:mg_nlevel+1),mg_x_idx(1:mg_nlevel+1),mg_size(1:3,1:mg_nlevel) )
      allocate ( mg_onorm(1:mg_nlevel) )

      mg_index_ext(1) = 1
      mg_index(1) = 1
      mg_x_idx(1) = 1
      mg_size(1,1) = xm
      mg_size(2,1) = ym
      mg_size(3,1) = zm
      m = xmymzm
      l = m + xmym
      n = l + xmym

      allocate ( zv(m) )

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

      allocate ( ad(1:m), ad1(1:m), bv(1:m), rv(1:m), iv(1:m), bz(1:m) )
      allocate ( am1(1:l), am2(1:l), am3(1:l) )
      allocate ( xv(1:n), xcorr(1:n) )
   case (3)
      m = xmymzm
      allocate ( ad(1:m), ad1(1:m), bv(1:m), zv(1:m), iv(1:m), tv(1:m), rv(1:m) )
      allocate ( am1(1-xmym:xmymzm), am2(1-xmym:xmymzm), am3(1-xmym:xmymzm) )
      allocate ( pv(1-xmym:m+xmym) )
   case (4)
      m = xmymzm
      allocate ( ad(1:m), ad1(1:m), bv(1:m), zv(1:m), iv(1:m) )
      allocate ( am1(1-xmym:xmymzm), am2(1-xmym:xmymzm), am3(1-xmym:xmymzm) )
   case (5)
      m = xmymzm
      allocate ( ad(1:m), ad1(1:m), bv(1:m), zv(1:m), iv(1:m) )
      allocate ( am1(1-xmym:xmymzm), am2(1-xmym:xmymzm), am3(1-xmym:xmymzm) )
   case (6)
      m = xmymzm
      allocate ( ad(1:m), ad1(1:m), bv(1:m), zv(1:m), iv(1:m), tv(1:m) )
      allocate ( am1(1-xmym:xmymzm), am2(1-xmym:xmymzm), am3(1-xmym:xmymzm) )
   case default
      m = xmymzm
      l = m + xmym
      n = l + xmym
      allocate ( ad(1:m), ad1(1:m), bv(1:m), rv(1:m), iv(1:m) )
      allocate ( tv(1:m), zv(1:m) )
      allocate ( am1(1-xmym:xmymzm), am2(1-xmym:xmymzm), am3(1-xmym:xmymzm) )
      allocate ( xv(1:n), xcorr(1:n), pv(1-xmym:m+xmym) )
   end select

end subroutine allocate_array
!=====================================================================================================================
subroutine deallocate_array( solvopt )

   integer solvopt 

   select case ( solvopt )
   case (1)
      deallocate( am1, am2, am3, xv, xcorr )
      deallocate( sav_am1, sav_am2, sav_am3, sav_ad )
      deallocate( ad, ad1, bv, rv, iv )
      deallocate( rd, pv, tv, zv, zzv )
   case (2)
      deallocate( am1, am2, am3, xv, xcorr )
      deallocate( ad, ad1, bv, rv, iv, zv, bz )
      deallocate( mg_index, mg_index_ext, mg_x_idx, mg_size, mg_onorm )
   case (3) 
      deallocate( am1, am2, am3 )
      deallocate( ad, ad1, bv, zv, iv, pv, tv, rv )
   case (4)
      deallocate( am1, am2, am3 )
      deallocate( ad, ad1, bv, zv, iv )
   case (5)
      deallocate( am1, am2, am3 )
      deallocate( ad, ad1, bv, zv, iv )
   case (6)
      deallocate( am1, am2, am3 )
      deallocate( ad, ad1, bv, zv, iv, tv )
   case default
      deallocate( am1, am2, am3, xv, xcorr )
      deallocate( ad, ad1, bv, rv, iv )
      deallocate( pv, tv, zv )
   end select

end subroutine deallocate_array
!==============================================================================
subroutine init_array(xs,epsx,epsy,epsz,p_iv,p_bv, &
                      solvopt,npbopt,h,epsout,eps0,pbtemp,pbkappa,nbnd,iepsav)

   integer npbopt, solvopt, nbnd
   integer iepsav(4,1:xmymzm)
   _REAL_ xs(1-xmym:xmymzm+xmym)
   _REAL_ epsx(1:xmymzm+ym*zm)
   _REAL_ epsy(1:xmymzm+xm*zm)
   _REAL_ epsz(1:xmymzm+xm*ym)
   _REAL_ p_iv(1:xmymzm), p_bv(1:xmymzm)
   _REAL_ h, epsout, eps0, pbtemp, pbkappa

   integer i,k,l,m,n,n1,n2,lxmym
   _REAL_ lfactor
   _REAL_,allocatable :: lepsx(:),lepsy(:),lepsz(:)

   select case ( solvopt )
   case (1)
      am1 = ZERO
      am2 = ZERO
      am3 = ZERO
      ad  = ZERO
      ad1 = ZERO
      rv  = ZERO
      bv  = ZERO
      pv  = ZERO
      tv  = ZERO
      zv  = ZERO
      rd  = ZERO
      zzv = ZERO
      xv  = ZERO
      xcorr  = ZERO
      iv(1:xmymzm) = p_iv(1:xmymzm)
      xv(1+xmym:xmymzm+xmym) = xs(1:xmymzm) 

      call feedepsintoam(epsx,epsy,epsz,am1(xmym+1),am2(xmym+1),am3(xmym+1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1:xmym) = ZERO
      am2(1:xmym) = ZERO
      am3(1:xmym) = ZERO
      call pb_setupper(am1(xmym+1),am2(xmym+1),am3(xmym+1))
      sav_am1 = am1
      sav_am2 = am2
      sav_am3 = am3
      sav_ad = ad

      if ( npbstep == 1 ) then
         bv(1+xmym:xmymzm+xmym) = p_bv(1:xmymzm)
      else
         call resid(xv,sav_am1,sav_am2,sav_am3, &
                    sav_ad,ad1,p_bv,bv(1+xmym),iv,npbopt, &
                    nbnd,iepsav)
      end if
   case (2)
      am1 = ZERO
      am2 = ZERO
      am3 = ZERO
      ad  = ZERO
      ad1 = ZERO
      iv  = ZERO
      bv  = ZERO
      rv  = ZERO
      zv  = ZERO
      xv  = ZERO
      xcorr  = ZERO
      iv(1:xmymzm) = p_iv(1:xmymzm)
      xv(1+xmym:xmymzm+xmym) = xs(1:xmymzm) 

      lxmym = mg_size(1,1)*mg_size(2,1)
      call feedepsintoam(epsx,epsy,epsz,am1(lxmym+1),am2(lxmym+1),am3(lxmym+1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      bz (1:xmymzm) = ad(1:xmymzm) - ad1(1:xmymzm)
      am1(1:lxmym) = ZERO
      am2(1:lxmym) = ZERO
      am3(1:lxmym) = ZERO
      call pb_setupper(am1(lxmym+1),am2(lxmym+1),am3(lxmym+1))

      if ( npbstep == 1 ) then
         bv(1:xmymzm) = p_bv(1:xmymzm)
      else  
         call resid(xv,am1,am2,am3,ad,ad1,p_bv,bv,iv,npbopt,nbnd,iepsav)
      end if

      m = 0
      do i = 1, mg_nlevel
         m = m + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) 
      end do
      allocate ( lepsx(1:m), lepsy(1:m), lepsz(1:m) )
      call feedepsintoam(epsx,epsy,epsz,lepsx(1),lepsy(1),lepsz(1))
      lfactor = factor1

      do i = 2, mg_nlevel
         l = mg_index(i-1)
         m = mg_index(i)
         n = mg_index_ext(i)
         call restrict_eps_map(lepsx(l),lepsy(l),lepsz(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1), &
                               lepsx(m),lepsy(m),lepsz(m),mg_size(1,i),mg_size(2,i),mg_size(3,i) )
         call restrict_iv(iv(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1),&
                          iv(m),mg_size(1,i),mg_size(2,i),mg_size(3,i) )
         lfactor = lfactor * 4
         lxmym = mg_size(1,i)*mg_size(2,i)

         call set_am_ad(lepsx(m),lepsy(m),lepsz(m),iv(m), &
                        am1(n+lxmym), am2(n+lxmym), am3(n+lxmym), &
                        ad(m), ad1(m), bz(m), &
                        mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,epsout,npbopt )

      end do

      deallocate( lepsx, lepsy, lepsz )
   case (3)
      pv  = ZERO
      zv  = ZERO
      tv  = ZERO
      rv  = ZERO
      iv (1:xmymzm) = p_iv(1:xmymzm)
      bv (1:xmymzm) = p_bv(1:xmymzm)

      call feedepsintoam(epsx,epsy,epsz,am1(1),am2(1),am3(1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1-xmym:0) = ZERO
      am2(1-xmym:0) = ZERO
      am3(1-xmym:0) = ZERO
      call pb_setupper(am1(1),am2(1),am3(1))
   case (4)
      zv  = ZERO
      iv (1:xmymzm) = p_iv(1:xmymzm)
      bv (1:xmymzm) = p_bv(1:xmymzm)
      call feedepsintoam(epsx,epsy,epsz,am1(1),am2(1),am3(1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1-xmym:0) = ZERO
      am2(1-xmym:0) = ZERO
      am3(1-xmym:0) = ZERO
      call pb_setupper(am1(1),am2(1),am3(1))
   case (5)
      zv  = ZERO
      iv (1:xmymzm) = p_iv(1:xmymzm)
      bv (1:xmymzm) = p_bv(1:xmymzm)
      call feedepsintoam(epsx,epsy,epsz,am1(1),am2(1),am3(1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1-xmym:0) = ZERO
      am2(1-xmym:0) = ZERO
      am3(1-xmym:0) = ZERO
      call pb_setupper(am1(1),am2(1),am3(1))
   case (6)
      zv  = ZERO
      tv  = ZERO
      iv (1:xmymzm) = p_iv(1:xmymzm)
      bv (1:xmymzm) = p_bv(1:xmymzm)
      call feedepsintoam(epsx,epsy,epsz,am1(1),am2(1),am3(1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1-xmym:0) = ZERO
      am2(1-xmym:0) = ZERO
      am3(1-xmym:0) = ZERO
      call pb_setupper(am1(1),am2(1),am3(1))
   case default
      am1 = ZERO
      am2 = ZERO
      am3 = ZERO
      ad  = ZERO
      ad1 = ZERO
      iv  = ZERO
      bv  = ZERO
      rv  = ZERO
      pv  = ZERO
      tv  = ZERO
      zv  = ZERO
      xv  = ZERO
      xcorr  = ZERO
      xv(1+xmym:xmymzm+xmym) = xs(1:xmymzm) 

      call feedepsintoam(epsx,epsy,epsz,am1(1),am2(1),am3(1))
      call pb_initad( epsx,epsy,epsz,ad(1),ad1(1),iv(1) )
      am1(1-xmym:0) = ZERO
      am2(1-xmym:0) = ZERO
      am3(1-xmym:0) = ZERO
      call pb_setupper(am1(1),am2(1),am3(1))
   end select

end subroutine init_array
!=====================================================================================================================
subroutine pb_setupper( am1, am2, am3 )

   _REAL_ am1(xm,ym,zm), am2(xm,ym,zm), am3(xm,ym,zm)

   integer i,j,k

   do j = 1, ym; do k = 1, zm
      am1(xm,j,k) = ZERO
   end do; end do
   do i = 1, xm; do k = 1, zm
      am2(i,ym,k) = ZERO
   end do; end do
   do i = 1, xm; do j = 1, ym
      am3(i,j,zm) = ZERO
   end do; end do
 
end subroutine pb_setupper
!=====================================================================================================================
subroutine feedepsintoam(eps1,eps2,eps3,lam1,lam2,lam3)
   implicit none
   _REAL_  eps1(0:xm,1:ym,1:zm), eps2(1:xm,0:ym,1:zm), eps3(1:xm,1:ym,0:zm)
   _REAL_  lam1(1:xm,1:ym,1:zm), lam2(1:xm,1:ym,1:zm), lam3(1:xm,1:ym,1:zm)

   lam1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
   lam2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
   lam3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)
end subroutine feedepsintoam
!==============================================================================
subroutine pb_initad(eps1,eps2,eps3,l_ad,lad1,l_iv)
   implicit none
   _REAL_  eps1(0:xm,1:ym,1:zm), eps2(1:xm,0:ym,1:zm), eps3(1:xm,1:ym,0:zm)
   _REAL_  l_ad(1:xm,1:ym,1:zm), lad1(1:xm,1:ym,1:zm), l_iv(1:xm,1:ym,1:zm)

   lad1(1:xm,1:ym,1:zm) =                      eps1(1:xm,  1:ym,  1:zm)
   lad1(1:xm,1:ym,1:zm) = lad1(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
   lad1(1:xm,1:ym,1:zm) = lad1(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
   lad1(1:xm,1:ym,1:zm) = lad1(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
   lad1(1:xm,1:ym,1:zm) = lad1(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
   lad1(1:xm,1:ym,1:zm) = lad1(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)
   l_ad = l_iv*factor1 + lad1
end subroutine pb_initad
!==============================================================================
subroutine set_am_ad( epsx,epsy,epsz,iv,lam1,lam2,lam3,lad,lad1,lbz, &
                      xn,yn,zn,lfactor,epsout,npbopt )

   _REAL_  lfactor, epsout

   integer xn,yn,zn,npbopt
   _REAL_  epsx(xn,yn,zn), epsy(xn,yn,zn), epsz(xn,yn,zn), iv(xn,yn,zn)
   _REAL_  lam1(xn,yn,zn),lam2(xn,yn,zn),lam3(xn,yn,zn)
   _REAL_  lad(xn,yn,zn),lad1(xn,yn,zn),lbz(xn,yn,zn)

   integer  i,j,k,i1,j1,k1
   _REAL_ lam1t,lam2t,lam3t

   lam1(1:xn,1:yn,1:zn) = epsx
   lam2(1:xn,1:yn,1:zn) = epsy
   lam3(1:xn,1:yn,1:zn) = epsz

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

   lad1 = lad
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

!=====================================================================================================================

subroutine restrict_eps_map( epsx, epsy, epsz, xn, yn, zn, epsxr, epsyr, epszr, xnr, ynr, znr )

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

!=====================================================================================================================

subroutine restrict_iv(ivf, xn, yn, zn, ivr, xnr, ynr, znr)

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
!=====================================================================================================================
subroutine restrict_bv(bvf, xn, yn, zn, bvr, xnr, ynr, znr)

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

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_nmg( phi, xs, p_bv, epsout, accept, npbopt, nbnd, iepsav )

   integer npbopt, nbnd
   _REAL_ accept, epsout
   _REAL_ phi(1:xmymzm), xs(1-xmym:xmymzm+xmym), p_bv(1:xmymzm) 
   integer iepsav(4,1:xmymzm)

   integer p1,p2,itmax,l
   integer ditn, maxditn
   _REAL_ vaccept

   ncyc_before = 5 
   ncyc_after = 5 

   itn = 0
   idamp = 1
   ditn = 0
   maxditn = 10
   itmax = 1000
   vaccept = accept

   uconvg = .true.
   inorm = sum(abs(p_bv(1:xmymzm)))
!  write (6,*) "NMG itn & norm", 0, inorm

   p1 = 1+xmym 
   p2 = xmymzm+xmym

   if ( npbopt == 0 ) then
      norm = inorm * accept
   else
      norm = inorm
   end if
   savnorm2 = ZERO

   do while ( uconvg )

      xcorr = ZERO

      if ( npbopt /= 0 .and. itn > 0 ) call update_ad()
     
!     if ( mod(npbstep+1, npbgrid) == 0 ) then
!        call resid(xv,am1,am2,am3,ad,ad1,p_bv,bv,iv,npbopt,nbnd,iepsav)
!     end if

      call pb_fmg( epsout, itmax, vaccept, norm )

      itn = itn + 1
      xv(p1:p2) = xv(p1:p2) + xcorr(p1:p2)

      call resid(xv,am1,am2,am3,ad,ad1,p_bv,bv,iv,npbopt,nbnd,iepsav)

      call innt(xv,xcorr,am1,am2,am3, &
                ad,ad1,p_bv,bv,rv,iv, &
                ditn,maxditn,npbopt,nbnd,iepsav)

      savnorm2 = norm2
!     mg_xcorr(p1:p2) = damp*mg_xcorr(p1:p2)
!     write (6,*) "NMG itn & norm", itn, norm
         
      if ( itn .ge. maxitn .or. norm .le. accept*inorm ) then
         uconvg = .false.
!        write (6,*) "NMG itn & norm", itn, norm
         if ( itn .ge. maxitn ) then
            write(6, *) 'PBMD WARNING: Multigrid maxitn exceeded!'
         endif
      end if

   end do

   xs(1:xmymzm) = xv(p1:p2)

end subroutine pb_nmg
!===========================================================================================================
subroutine pb_fmg ( epsout, itmax, vaccept, nnorm )

!  use pbtimer_module
   
   integer litn, itmax
   _REAL_  vaccept,epsout,nnorm

   logical luconvg
   integer j,p1,p2
   _REAL_ linorm, lnorm, raccept, errtol
!  integer lx, ly, lz, lxly, lxlylz
   integer lxly, lxlylz

! calculate initial norm

   linorm = sum(abs(bv(1:xmymzm)))
!  write (6,*) "   Multigrid itn & norm", 0, linorm
!  inorm = ONE

   litn = 0
   luconvg = .true.
   itmax = 1000
   raccept = 1.0D-3
   errtol = 1.0d0

!  if ( nnorm > 0.9d0 ) then
!     errtol = 0.9d0 * nnorm
!  else 
!     errtol = nnorm * nnorm
!  end if

   do while ( luconvg )
      mg_onorm = 9.9d99
      do j = 1, mg_nlevel-1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)

!        call pbtimer_start(PBTIME_MGRLX)
         call relax(xcorr(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    am1(mg_index_ext(j)),am2(mg_index_ext(j)),am3(mg_index_ext(j)), &
                    ad(mg_index(j)),bv(mg_index(j)),rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_before,raccept,mg_onorm(j) )
!        call pbtimer_stop(PBTIME_MGRLX)

!        call pbtimer_start(PBTIME_MGRES)
         call restrict_bv(rv(mg_index(j)), mg_size(1,j), mg_size(2,j), mg_size(3,j), &
                          bv(mg_index(j+1)), mg_size(1,j+1),mg_size(2,j+1),mg_size(3,j+1) )
!        call pbtimer_stop(PBTIME_MGRES)

         xcorr(mg_x_idx(j+1):mg_x_idx(j+2)-1) = ZERO
      end do

!     call pbtimer_start(PBTIME_MGRLX)
      lxly = mg_size(1,j)*mg_size(2,j)
      lxlylz = lxly*mg_size(3,j)
      call relax(xcorr(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                 am1(mg_index_ext(j)),am2(mg_index_ext(j)),am3(mg_index_ext(j)), &
                 ad(mg_index(j)),bv(mg_index(j)),rv(mg_index(j)), &
                 lxly,lxlylz,-1,raccept,mg_onorm(j) )
!     call pbtimer_stop(PBTIME_MGRLX)

      do j = mg_nlevel-1, 1, -1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         p1 = mg_x_idx(j+1)+mg_size(1,j+1)*mg_size(2,j+1)
         p2 = mg_x_idx(j)+lxly

!        call pbtimer_start(PBTIME_MGINT)
         call interpolate(xcorr(p1), mg_size(1,j+1), mg_size(2,j+1), mg_size(3,j+1),&
                          xcorr(p2), mg_size(1,j) ,  mg_size(2,j),   mg_size(3,j) , &
                          am1(mg_index_ext(j)),am2(mg_index_ext(j)),am3(mg_index_ext(j)),bz(mg_index(j)),epsout)
!        call pbtimer_stop(PBTIME_MGINT)

!        call pbtimer_start(PBTIME_MGRLX)
         call relax(xcorr(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    am1(mg_index_ext(j)),am2(mg_index_ext(j)),am3(mg_index_ext(j)), &
                    ad(mg_index(j)),bv(mg_index(j)), rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_after,raccept,mg_onorm(j) )
!        call pbtimer_stop(PBTIME_MGRLX)
      end do

      litn = litn + 1
      lnorm = sum(abs(rv(1:xmymzm)))
!     write(6, *)  '   Multigrid itn & norm', litn, lnorm

!  
!  convergence
!  
      if ( litn .ge. itmax .or. lnorm .le. errtol ) then 
         luconvg = .false.
!        write(6, *)  '   Multigrid itn & norm', litn, lnorm
         if ( litn .ge. itmax ) then
            write(6, *) 'PBMD WARNING: Multigrid maxitn exceeded!'
         endif
      end if

   end do

end subroutine pb_fmg
!=====================================================================================================================
subroutine update_ad

   integer j

   bz(1:xmymzm) = iv(1:xmymzm) * cosh( iv(1:xmymzm)*xv(1+xmym:xmymzm+xmym)*ktinv ) * factor1  
   ad(1:xmymzm) = ad1(1:xmymzm) + bz(1:xmymzm)

   do j = 2, mg_nlevel
      call restrict_ad(mg_size(1,j-1),mg_size(2,j-1),mg_size(3,j-1),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                       bz(mg_index(j-1)),bz(mg_index(j)),ad(mg_index(j)),ad1(mg_index(j)))
   end do

end subroutine update_ad
!=====================================================================================================================
subroutine restrict_ad(xf,yf,zf,xc,yc,zc,bzf,bzc,adc,ad1c)

   integer xf,yf,zf,xc,yc,zc
   _REAL_ bzf(xf,yf,zf), bzc(xc,yc,zc)
   _REAL_ adc(xc,yc,zc), ad1c(xc,yc,zc)

   integer i,j,k,i2,j2,k2

   do k = 1, zc
      k2 = 2*k
      do j = 1, yc
         j2 = 2*j
         do i = 1, xc
            i2 = 2*i
            bzc(i,j,k) =       ( bzf(i2-1,j2-1,k2-1) + &
                                 TWO * bzf(i2  ,j2-1,k2-1) + &
                                 bzf(i2+1,j2-1,k2-1) ) + &
                        TWO  * ( bzf(i2-1,j2  ,k2-1) + &
                                 TWO * bzf(i2  ,j2  ,k2-1) + &
                                 bzf(i2+1,j2  ,k2-1) ) + &
                               ( bzf(i2-1,j2+1,k2-1) + &
                                 TWO * bzf(i2  ,j2+1,k2-1) + &
                                 bzf(i2+1,j2+1,k2-1) ) + &
                        TWO  * ( bzf(i2-1,j2-1,k2  ) + &
                                 TWO * bzf(i2  ,j2-1,k2  ) + &
                                 bzf(i2+1,j2-1,k2  ) ) + &
                        FOUR * ( bzf(i2-1,j2  ,k2  ) + &
                                 TWO * bzf(i2  ,j2  ,k2  ) + &
                                 bzf(i2+1,j2  ,k2  ) ) + &
                        TWO  * ( bzf(i2-1,j2+1,k2  ) + &
                                 TWO * bzf(i2  ,j2+1,k2  ) + &
                                 bzf(i2+1,j2+1,k2  ) ) + &
                               ( bzf(i2-1,j2-1,k2+1) + &
                                 TWO * bzf(i2  ,j2-1,k2+1) + & 
                                 bzf(i2+1,j2-1,k2+1) ) + &
                        TWO  * ( bzf(i2-1,j2  ,k2+1) + &
                                 TWO * bzf(i2  ,j2  ,k2+1) + &
                                 bzf(i2+1,j2  ,k2+1) ) + &
                               ( bzf(i2-1,j2+1,k2+1) + &
                                 TWO * bzf(i2  ,j2+1,k2+1) + & 
                                 bzf(i2+1,j2+1,k2+1) ) 
            bzc(i,j,k) = bzc(i,j,k) / 16.d0
            adc(i,j,k) = ad1c(i,j,k) + bzc(i,j,k)
            !adc(i,j,k) = ad1c(i,j,k) + adf(i2+1,j2+1,k2+1) - ad1f(i2+1,j2+1,k2+1)
         end do
      end do
   end do 

end subroutine restrict_ad
!=====================================================================================================================

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

   lad = am_1(l+shift_3) + am_1(l+shift_3-abs(shift_1)) + am_2(l+shift_3) + am_2(l+shift_3-abs(shift_2)) + lbz(l+shift_3)
   if ( shift_3 < 0 ) then
      ipl_comp3 = v * am_3(l+shift_3) / ( am_3(l+2*shift_3) + am_3(l+shift_3) + lad)
   else
      ipl_comp3 = v * am_3(l) / ( am_3(l) + am_3(l+shift_3) + lad )
   end if
   return

end function ipl_comp3

end subroutine interpolate
!=====================================================================================================================

subroutine relax(xs,nx,ny,nz,lam1,lam2,lam3,lad,lbv,lrv,nxny,nxnynz,ncyc,accept,onorm)

   integer nx,ny,nz,nxny,nxnynz,ncyc
   _REAL_ xs(1-nxny:nxnynz+nxny), lam1(1-nxny:nxnynz), lam2(1-nxny:nxnynz), lam3(1-nxny:nxnynz)
   _REAL_ lad(1:nxnynz), lbv(1:nxnynz), lrv(1:nxnynz)
   _REAL_ accept, onorm

   logical luconvg
   integer ii,litn,itmax
   _REAL_ wsor, wsor1, linorm, lnorm

   if (ncyc>0) then
      itn_checknorm = ncyc
      wsor = 1.5d0
      wsor1 = ONE - wsor
   else
      itn_checknorm = 10
      wsor = 1.9d0
      wsor1 = ONE - wsor
   end if

   linorm = sum(abs(lbv(1:nxnynz)))
!  write(6, *)  '      relax itn & norm ', litn, linorm!, onorm
   zv(1:nxnynz) = ONE/lad(1:nxnynz)

   litn = 0
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
                                          lbv(ii     )             ) * zv(ii)
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
         if ( litn .ge. itmax .or. ( ncyc .gt. 0 .and. (litn .ge. ncyc .and. lnorm < onorm ) )  .or. lnorm .le. accept*linorm ) then

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_ncg ( phi, xs, accept, npbopt )

!  Passed variables

   integer npbopt
   _REAL_ phi(xmymzm), xs(1-xmym:xmymzm+xmym), accept

!  Local variables

   integer i
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2, bdotp
!  _REAL_ bdotc, cdotc
   
!
!  begin code
!
!
!  initialization
!
!  compute ||b||

   inorm = sum(ABS(bv(1:xmymzm)))
!  write(6, *)  'itn & inorm ', 0, inorm
!
!  compute b - A * x(0) and save it in r(0)
!  p(0) = r(0)
!
!  iteration 0:
!  compute <r(0),r(0)>
!
   if ( npbopt /= 0 ) ad(1:xmymzm) = ad1(1:xmymzm) 
   itn = 0
   do i = 1,xmymzm
       bv(i)  = bv(i)  + am3(i-xmym)*xs(i-xmym) &
                       + am2(i-xm  )*xs(i-xm  ) &
                       + am1(i-1   )*xs(i-1   ) &
                       + am1(i     )*xs(i+1   ) &
                       + am2(i     )*xs(i+xm  ) &
                       + am3(i     )*xs(i+xmym) &
                       - ad (i     )*xs(i     ) 
   end do
   bdotb1 = dot_product(bv(1:xmymzm),bv(1:xmymzm))

   if ( npbopt /= 0 ) then 
      rv(1:xmymzm) = sinh( iv(1:xmymzm) * xs(1:xmymzm) * ktinv )
      bv(1:xmymzm) = bv(1:xmymzm) - factor * rv(1:xmymzm) 
      bdotp = dot_product( bv(1:xmymzm), bv(1:xmymzm) )
   end if
   
   pv(1:xmymzm) = bv(1:xmymzm)

   norm = sum(abs(bv(1:xmymzm)))
!  write(6, *)  'itn & norm, bdotb1 ', itn, norm, bdotb1
!  write(6, *)  'bcg', bcg

!
!  the main loop of the CG solver
!
   uconvg = .true.
   do while ( uconvg )
!
!  iteration i:
!
!     compute Ap(i) = A * p(i)
!     compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
!
      do i = 1,xmymzm
         zv(i) = - am3(i-xmym)*pv(i-xmym) &
                 - am2(i-xm  )*pv(i-xm  ) &
                 - am1(i-1   )*pv(i-1   ) &
                 - am1(i     )*pv(i+1   ) &
                 - am2(i     )*pv(i+xm  ) &
                 - am3(i     )*pv(i+xmym) & 
                 + ad (i)     *pv(i)      
      end do
      pdotz = dot_product( pv(1:xmymzm), zv(1:xmymzm) )
!
!  iteration i+1:
!
      itn = itn + 1
!
!     Newton's root finding method to solve for alpha //CQ
!
      alpha  = bdotb1/pdotz
      if ( npbopt /= 0 ) call ntalpha(xs(1))
!
!     update x(i+1) = x(i) + alpha(i) p(i)
!            r(i+1) = r(i) - alpha(i) Ap(i)
!
!     revise bv (residual) for nonlinear PBE case //CQ
!
      xs(1:xmymzm) = xs(1:xmymzm) + alpha*pv(1:xmymzm)
!     cbv(1:xmymzm) = -alpha*zv(1:xmymzm)
       bv(1:xmymzm) = bv(1:xmymzm) - alpha*zv(1:xmymzm)
      if ( npbopt /= 0 ) then 
         tv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv
!        maxbf = maxval(tv(1:xmymzm))
         tv(1:xmymzm) = sinh( tv(1:xmymzm) )
!        cbv(1:xmymzm) = cbv(1:xmymzm) - factor * ( tv(1:xmymzm) - rv(1:xmymzm) )
         bv(1:xmymzm) = bv(1:xmymzm) - factor * ( tv(1:xmymzm) - rv(1:xmymzm) )
         rv(1:xmymzm) = tv(1:xmymzm)
      end if 
!     bv(1:xmymzm) = bv(1:xmymzm) + cbv(1:xmymzm)
      
!
!     compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
      norm = sum(abs(bv(1:xmymzm)))
      bdotb2 = dot_product(bv(1:xmymzm),bv(1:xmymzm))

!     if (npbopt /= 0) write(6, *)  'max boltzmann factor', maxbf
!     write(6, *)  'itn & norm',itn, norm
!
!     check convergence
!
      if ( itn .ge. maxitn .or. norm .le. accept*inorm ) then

         uconvg = .false.
         ad(1:xmymzm) = factor * tv(1:xmymzm)
         if ( itn .ge. maxitn ) then
            write(6, *) 'PBMD WARNING: CG maxitn exceeded!'
         endif

      else
!
!     compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two
!
!        if ( bcg == 1 ) then
            beta = bdotb2/bdotb1
!        else
!           bdotc = dot_product(bv(1:xmymzm),cbv(1:xmymzm))
!           cdotc = dot_product(cbv(1:xmymzm),cbv(1:xmymzm))
!           beta = bdotc/bdotb1
!           if ( beta < ZERO ) beta = ZERO 
!        end if 
         bdotb1 = bdotb2
!
!     update p(i+1) = r(i+1) + beta(i) p(i)
!
         pv(1:xmymzm) = bv(1:xmymzm) + beta*pv(1:xmymzm)
         if ( npbopt /= 0 ) bdotp = dot_product(pv(1:xmymzm),bv(1:xmymzm))
      end if

   end do !  end of the main CG loop
!  write(6, *)  'itn & norm ',itn, norm

contains 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Newton's root finding method to solve for alpha  //CQ
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ntalpha(xs)

   _REAL_ xs(1:xmymzm)
   
!
! Local variables
!
   integer alpcnt,alpimx
   _REAL_ alpcor,alptol,snhtmp,cshtmp
   
   alpimx = 5  ! alpimx = 5 -> 50   //CQ
   alptol = 0.1d0

   alpcor = 1.0d0

   alpcnt = 0
   do while ( (abs(alpcor) > alptol) .and. (alpcnt < alpimx) )

      alpcnt = alpcnt + 1

      tv(1:xmymzm) = iv(1:xmymzm) * ( xs(1:xmymzm) + alpha * pv(1:xmymzm) ) * ktinv
      snhtmp = factor * dot_product( pv(1:xmymzm), sinh(tv(1:xmymzm)) - rv(1:xmymzm) )
      cshtmp = factor * dot_product( pv(1:xmymzm)**2, iv(1:xmymzm)*cosh(tv(1:xmymzm)) ) * ktinv

      alpcor = - (alpha*pdotz-bdotp+snhtmp)/(pdotz+cshtmp)
      alpha = alpha + alpcor
!     write(6,*) 'alpcnt alpha ', alpcnt, alpha
      ! why do you need to do this correction below?
!     if ( alpcnt == alpimx ) then
!        alpha = alpha - alpcor/TWO
!     end if

   end do

end subroutine ntalpha

end subroutine pb_ncg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dsor ( phi, xs, wsor, accept, npbopt )

!  use pbtimer_module

!  Passed variables
    
   integer npbopt
   double precision phi(xmymzm), xs(1-xmym:xmymzm+xmym), accept
   _REAL_ wsor

!  local variables
   integer ii, ditn, maxditn
   _REAL_ wsor1

   wsor1 = ONE - wsor
   ditn = 0
   maxditn = 10
   idamp = 0
   savnorm2 = 9.9d99

!  if ( npbstep == 0 .and. npbopt /= 0 ) then
!     w = wsor
!     wsor = 1.5d0
!     wsor = ONE 
!  end if

!  calculate initial norm

   inorm = sum( abs( bv(1:xmymzm) ) )
!  write(6, *) 'itn & norm ', 0, inorm

!  the main loop

   uconvg = .true. 
   itn = 0
   if ( npbopt == 0 ) zv(1:xmymzm) = ONE/ad(1:xmymzm)
   if ( npbopt /= 0 .and. wsor > ONE ) idamp = 1

   do while ( uconvg )
      itn = itn + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     set up inverse diagonal elements, with salt if possible
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( npbopt /= 0 ) then
      !if (itn == 1) call pbtimer_start(PBTIME_SINH)

         !Nicholas without division by xs
         zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv 
!        maxbf = maxval(abs(zv(1:xmymzm)))
         ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
         zv(1:xmymzm) = ad1(1:xmymzm) * xs(1:xmymzm) + ad(1:xmymzm)
         zv(1:xmymzm) = xs(1:xmymzm)/zv(1:xmymzm)

      !if (itn == 1) call pbtimer_stop(PBTIME_SINH)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     start the sor iteration ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Nicholas with damping without division by xs
      if (npbopt == 0) then
         do ii = 1, xmymzm
            xs(ii) = wsor1*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
         end do
      else
         do ii = 1, xmymzm
            tv(ii) = -wsor*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
            xs(ii) = xs(ii) + tv(ii)
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     check convergence every ten steps ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      !Nicholas without division by xs
      if ( npbopt /= 0 ) then
         do ii = 1, xmymzm
            phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                      am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                      am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                          bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
         end do
         norm2 = dot_product(phi(1:xmymzm),phi(1:xmymzm))
         norm = sum(abs(phi(1:xmymzm)))

         ! with damping
         if ( idamp == 1 ) then
            damp = ONE 
            ditn = 0

            do while( norm2 > savnorm2 .and. ditn < maxditn )
               damp = damp / TWO
               xs(1:xmymzm) = xs(1:xmymzm) - damp * tv(1:xmymzm)

               zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv
               ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
               do ii = 1, xmymzm
                  phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                            am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                            am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                                bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
               end do

               norm = sum(abs(phi(1:xmymzm)))
               norm2 = dot_product(phi(1:xmymzm),phi(1:xmymzm))
               ditn = ditn + 1
            end do

            if ( ditn < maxditn ) then

               damp = damp / TWO
               xs(1:xmymzm) = xs(1:xmymzm) - damp * tv(1:xmymzm)

               zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv
               ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
               do ii = 1, xmymzm
                  phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                            am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                            am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                                bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
               end do

               savnorm = norm
               savnorm2 = norm2
               norm = sum(abs(phi(1:xmymzm)))
               norm2 = dot_product(phi(1:xmymzm),phi(1:xmymzm))
               ditn = ditn + 1

               do while ( norm2 < savnorm2 .and. ditn < maxditn ) 
                  damp = damp / TWO
                  xs(1:xmymzm) = xs(1:xmymzm) - damp * tv(1:xmymzm)

                  zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv
                  ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
                  do ii = 1, xmymzm
                     phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                               am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                               am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                                   bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
                  end do

                  savnorm = norm
                  savnorm2 = norm2
                  norm = sum(abs(phi(1:xmymzm)))
                  norm2 = dot_product(phi(1:xmymzm),phi(1:xmymzm))
                  ditn = ditn + 1
               end do 

               if ( norm2 >= savnorm2 ) then
                  xs(1:xmymzm) = xs(1:xmymzm) + damp * tv(1:xmymzm)
                  norm = savnorm
                  norm2 = savnorm2
                  ditn = ditn - 1
               end if

            end if

!           if ( ditn == 0 ) idamp = 0

!           write(6, *) '   ditn, norm2', ditn, norm2
         end if 
!        write(6, *) '   norm', norm

         savnorm2 = norm2
         
      end if

      if ( mod(itn, itn_checknorm) == 0 ) then

         if ( npbopt /= 0 ) then
            zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv
            ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
            end do
         else
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) -  ad(ii)* xs(ii     )
            end do
         end if
         norm = sum(abs(phi(1:xmymzm)))

!        write(6, *) 'itn, real norm', itn, norm

         if ( itn >= maxitn .or. norm <= accept*inorm ) then
            uconvg = .false.
!           write(6, *) 'itn & norm', itn, norm
            if ( itn .ge. maxitn ) write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
         end if

      end if

   end do

!  write(6, *) 'itn & norm ', itn, norm
!  if ( npbstep == 0 .and. npbopt /= 0 ) wsor = w

end subroutine pb_dsor
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_nticcg( phi, xs, p_bv, accept, npbopt, nbnd, iepsav )

!  use pb_mg

   integer npbopt, nbnd
   _REAL_ accept
   _REAL_ phi(1:xmymzm), xs(1-xmym:xmymzm+xmym)
   _REAL_ p_bv(1:xmymzm) 
   integer iepsav(4,1:xmymzm)

   integer itmax, ditn, maxditn, l
   _REAL_ vaccept, pnorm1, pnorm2

   itn = 0
   idamp = 1
   ditn = 0
   maxditn = 10
   vaccept = 0.1d0

   uconvg = .true.
   inorm = sum(abs(p_bv(1:xmymzm)))
   pnorm1 = inorm
   pnorm2 = 0

!  write (6,*) "NT itn & norm", 0, inorm

   itmax = 1000
   if ( npbopt == 0 ) then
      norm = accept * inorm
   else
      norm = inorm
   end if
   savnorm2 = ZERO

   do while ( uconvg )

      xcorr = ZERO

      am1 = sav_am1
      am2 = sav_am2
      am3 = sav_am3

      ! update ad
      if ( npbopt /= 0 .and. itn > 0 ) then
         ad(1:xmymzm) = ad1(1:xmymzm) + iv(1:xmymzm) * cosh( iv(1:xmymzm)*xv(1+xmym:xmymzm+xmym)*ktinv ) * factor1 
         sav_ad = ad
      end if

!     if ( mod(npbstep+1,npbgrid) == 0 ) then
!        call resid(xv,sav_am1,sav_am2,sav_am3, &
!                   sav_ad,ad1,p_bv,bv(1+xmym),iv,npbopt, &
!                   nbnd,iepsav)
!     end if

      call nt_iccg( xcorr, am1, am2, am3, ad, &
                    bv, rd, pv, tv, zv, &
                    vaccept, pnorm1, pnorm2)   

      itn = itn + 1
      xv = xv + xcorr

      ad = sav_ad
      call resid(xv,sav_am1,sav_am2,sav_am3, &
                 ad,ad1,p_bv,bv(1+xmym),iv,npbopt, &
                 nbnd,iepsav)
!     vaccept = abs(norm-pnorm2) / pnorm1
!     pnorm1 = norm

      if ( npbopt /= 0 ) then
         call innt(xv,xcorr,sav_am1,sav_am2,sav_am3, &
                   ad,ad1,p_bv,bv(1+xmym),rv,iv, &
                   ditn,maxditn,npbopt,nbnd,iepsav)
      end if

      savnorm2 = norm2
!     write (6,*) "NT itn & norm", itn, norm
         
      if ( itn .ge. maxitn .or. norm .le. accept*inorm ) then
         uconvg = .false.
!        write (6,*) "NT itn & norm", itn, norm
         if ( itn .ge. maxitn ) then
            write(6, *) 'PBMD WARNING: Multigrid maxitn exceeded!'
         endif
      end if

   end do

   xs(1:xmymzm) = xv(1+xmym:xmymzm+xmym)

end subroutine pb_nticcg
!=====================================================================================================================
subroutine innt(lxv,lxcorr,lam1,lam2,lam3,lad,lad1,ibv,lbv,lrv,liv, &
                ditn,maxditn,npbopt,nbnd,iepsav)

   integer ditn, maxditn, npbopt, nbnd
   integer iepsav(4,1:xmymzm)
   _REAL_  lxv(1-xmym:xmymzm+xmym), lxcorr(1-xmym:xmymzm+xmym)
   _REAL_  lam1(1-xmym:xmymzm), lam2(1-xmym:xmymzm), lam3(1-xmym:xmymzm)
   _REAL_  lad(1:xmymzm), lad1(1:xmymzm)
   _REAL_  ibv(1:xmymzm), lbv(1:xmymzm), lrv(1:xmymzm), liv(1:xmymzm)

   ! with damping
   if ( idamp == 1 ) then
      damp = ONE 
      ditn = 0

      do while( norm  > savnorm  .and. ditn < maxditn .and. itn > 1 )
         damp = damp / TWO
         lxv = lxv - damp * lxcorr

         call resid(lxv,lam1,lam2,lam3, &
                    lad,lad1,ibv,lbv,liv,npbopt, &
                    nbnd,iepsav)
         ditn = ditn + 1
      end do

      if ( ditn < maxditn ) then

         damp = damp / TWO
         lxv = lxv - damp * lxcorr

         savnorm = norm
         savnorm2 = norm2

         call resid(lxv,lam1,lam2,lam3, &
                    lad,lad1,ibv,lrv,liv,npbopt, &
                    nbnd,iepsav)

         ditn = ditn + 1

         do while ( norm  < savnorm  .and. ditn < maxditn ) 
            damp = damp / TWO
            lxv = lxv - damp * lxcorr
            lbv = lrv

            savnorm = norm
            savnorm2 = norm2

            call resid(lxv,lam1,lam2,lam3, &
                       lad,lad1,ibv,lrv,liv,npbopt, &
                       nbnd,iepsav)

            ditn = ditn + 1
         end do 

         if ( norm  >= savnorm  ) then
            lxv = lxv + damp * lxcorr
            norm = savnorm
            norm2 = savnorm2
            ditn = ditn - 1
         else
            lbv = lrv
         end if

      end if

!     write (6,*) "   ditn & norm2", ditn, norm2
      if ( ditn == 0 ) idamp = 0
   end if

end subroutine innt
!=====================================================================================================================
subroutine resid(lxs,lam1,lam2,lam3,lad,lad1,lbv,lrv,liv,npbopt,nbnd,iepsav)

!  use pbtimer_module
   
   integer nbnd,npbopt
   _REAL_ lxs(1-xmym:xmymzm+xmym), lam1(1-xmym:xmymzm), lam2(1-xmym:xmymzm), lam3(1-xmym:xmymzm)
   _REAL_ lad(1:xmymzm), lad1(1:xmymzm), lbv(1:xmymzm), lrv(1:xmymzm), liv(1:xmymzm)
   integer iepsav(4,1:xmymzm)

   integer ii,ip,i,j,k

   if ( npbopt == 0 ) then
      do ii = 1,xmymzm 
         lrv(ii) = lam1(ii-1) * lxs(ii-1) + lam1(ii)* lxs(ii+1) + &
                  lam2(ii-xm) * lxs(ii-xm) + lam2(ii)*lxs(ii+xm) + &
                  lam3(ii-xmym) * lxs(ii-xmym) + lam3(ii)*lxs(ii+xmym) + &
                  lbv(ii) - lad(ii)* lxs(ii)
      end do
   else
      do ii = 1,xmymzm 
         lrv(ii) = lam1(ii-1) * lxs(ii-1) + lam1(ii)* lxs(ii+1) + &
                   lam2(ii-xm) * lxs(ii-xm) + lam2(ii)*lxs(ii+xm) + &
                   lam3(ii-xmym) * lxs(ii-xmym) + lam3(ii)*lxs(ii+xmym) + &
                   lbv(ii) - lad1(ii)* lxs(ii)
      end do
      lad(1:xmymzm) = factor * sinh( liv(1:xmymzm) * lxs(1:xmymzm) * ktinv )
      lrv(1:xmymzm) = lrv(1:xmymzm) - lad(1:xmymzm)
   end if

   norm = sum(abs(lrv(1:xmymzm)))
!  write(6,*) 'resid norm',norm
   norm2 = dot_product(lrv(1:xmymzm),lrv(1:xmymzm))
!  write(6,*) 'nbnd',nbnd

end subroutine resid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ MICCG core iteration routine
subroutine nt_iccg ( xs, lam1, lam2, lam3, lad, &
                     lbv, lrd, lpv, ltv, lzv, &
                     vaccept, pnorm1, pnorm2 )   

   ! Passed variables
   
   _REAL_ xs(1-xmym:xmymzm+xmym) 
   _REAL_ lam1(1-xmym:xmymzm+xmym), lam2(1-xmym:xmymzm+xmym), lam3(1-xmym:xmymzm+xmym)
   _REAL_ lad(1:xmymzm+xmym), lbv(1-xmym:xmymzm)
   _REAL_ lrd(1-xmym:xmymzm+xmym), ltv(1-xmym:xmymzm+xmym)
   _REAL_ lzv(1-xmym:xmymzm+xmym), lpv(1-xmym:xmymzm+xmym)
   _REAL_ vaccept, pnorm1, pnorm2
    
   ! Local variables
    
   logical luconvg
   integer i,j,litn,mitn
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2
   _REAL_ fmiccg, linorm!, lnorm
   _REAL_ lnorm1,lnorm2!,la

!  if ( itn == 0 .and. npbstep > 0 .and. mod(npbstep+1,npbgrid) > 0 ) then
!     la = 0.01d0
!  else
!     la = 0.01d0
!     la = 1.0d0
!  end if

   fmiccg = -0.30d0
   mitn = 10000

!  write(6, *)  '   iccg itn & norm & vaccept', 0, norm, vaccept
 
   ! initialization
!  write(6,*) "fmiccg", fmiccg
   do i = 1, xmymzm
      lrd(i) = ONE/( lad(i) - &
         lam1(i-1   )*(       lam1(i-1   )+fmiccg*lam2(i-1   )+fmiccg*lam3(i-1   ))*lrd(i-1   ) - &
         lam2(i-xm  )*(fmiccg*lam1(i-xm  )+       lam2(i-xm  )+fmiccg*lam3(i-xm  ))*lrd(i-xm  ) - &
         lam3(i-xmym)*(fmiccg*lam1(i-xmym)+fmiccg*lam2(i-xmym)+       lam3(i-xmym))*lrd(i-xmym) )
   end do
    
   do i = 1, xmymzm
      lad(i) = lad(i)*lrd(i)
      lrd(i) = sqrt(lrd(i))
      lbv(i) = lbv(i)*lrd(i)
      lam1(i-1   ) = lam1(i-1   )*lrd(i)*lrd(i-1   )
      lam2(i-xm  ) = lam2(i-xm  )*lrd(i)*lrd(i-xm  )
      lam3(i-xmym) = lam3(i-xmym)*lrd(i)*lrd(i-xmym)
   end do
    
   linorm = ZERO
   do i = 1, xmymzm
      lad(i) = lad(i) - TWO

      lbv(i) = lbv(i) + lam1(i-1   )*lbv(i-1   ) &
                      + lam2(i-xm  )*lbv(i-xm  ) &
                      + lam3(i-xmym)*lbv(i-xmym)
      linorm = linorm + abs(lbv(i))
   end do
    
   do i = xmymzm, 1, -1
      ltv(i) = xs(i) + lam1(i     )*ltv(i+1   ) &
                     + lam2(i     )*ltv(i+xm  ) &
                     + lam3(i     )*ltv(i+xmym)
   end do
   do i = 1, xmymzm
      lzv(i) = xs(i) + lad (i     )*ltv(i     ) &
                     + lam1(i-1   )*lzv(i-1   ) &
                     + lam2(i-xm  )*lzv(i-xm  ) &
                     + lam3(i-xmym)*lzv(i-xmym)
   end do
   bdotb1 = ZERO
   do i = xmymzm, 1, -1
      lzv(i) = lzv(i) + ltv(i)
      lbv(i) = lbv(i) - lzv(i)
      
      ! iteration 0.
      
      bdotb1 = bdotb1 + lbv(i)*lbv(i)
      lpv(i) = lbv(i)
      
      ! first step of the matrix vector multiplication, see below
      
      ltv(i) = lpv(i) + lam1(i     )*ltv(i+1   ) &
                      + lam2(i     )*ltv(i+xm  ) &
                      + lam3(i     )*ltv(i+xmym)
   end do
    
   litn = 0
   luconvg = .true.
    
   ! the main loop of iccg solver
    
   do while ( luconvg )
       
      ! second and thilrd steps of the matrix vector multiplication
       
      pdotz = ZERO
      do i = 1, xmymzm+xmym
         lzv(i) = lpv(i) + lad (i     )*ltv(i     ) &
                         + lam1(i-1   )*lzv(i-1   ) &
                         + lam2(i-xm  )*lzv(i-xm  ) &
                         + lam3(i-xmym)*lzv(i-xmym)
         
         j = i - xmym
         lzv(j) = lzv(j) + ltv(j)
         
         pdotz = pdotz + lpv(j)*lzv(j)
      end do
      alpha = bdotb1/pdotz
       
!     lnorm  = ZERO
      lnorm1 = ZERO
      lnorm2 = ZERO
      bdotb2 = ZERO

      do i = 1, xmymzm
         xs(i)  = xs(i) + alpha*lpv(i)
         lbv(i) = lbv(i) - alpha*lzv(i)
!        lnorm  = lnorm  + abs(lbv(i))

         zzv(i) = lbv(i) - lam1(i-1   )*lbv(i-1   ) &
                         - lam2(i-xm  )*lbv(i-xm  ) &
                         - lam3(i-xmym)*lbv(i-xmym)
         zzv(i) = zzv(i)/lrd(i)
!        lnorm2: scaled norm
!        lnorm1: orgininal norm
         lnorm1 = lnorm1 + abs(zzv(i))
         lnorm2 = lnorm2 + abs(lbv(i))
 
         bdotb2 = bdotb2+ lbv(i)*lbv(i)
      end do
       
      litn = litn + 1
!     write(6, *)  '   iccg itn & norm',litn, lnorm1, lnorm2
      
      ! check convergence
       
      if ( litn >= mitn .or. lnorm1 < vaccept * norm ) then
!     if ( litn >= mitn .or. lnorm2 < norm ) then
          
         luconvg = .false.
         pnorm2 = lnorm1

!        write(6, *)  '   iccg itn & norm',litn, lnorm1
         if ( litn >= mitn ) then
            write(6, *) 'PB warning in pb_miccg(): CG maxitn exceeded!'
         end if
          
      else
          
         beta = bdotb2/bdotb1
         bdotb1 = bdotb2
          
         ! first step of the matrix vector multiplication
          
         do i = xmymzm, 1, -1
            lpv(i) = lbv(i) + beta*lpv(i)
             
            ltv(i) = lpv(i) + lam1(i)*ltv(i+1   ) &
                            + lam2(i)*ltv(i+xm  ) &
                            + lam3(i)*ltv(i+xmym)
         end do
      end if

   end do  !  while ( uconvg ), end of the main iccg loop
    
   ! back scaling of the solution

   do i = xmymzm, 1, -1
      ltv(i)  = xs(i) + lam1(i)*ltv(i+1   ) &
                      + lam2(i)*ltv(i+xm  ) &
                      + lam3(i)*ltv(i+xmym)
       
      xs(i) = ltv(i)*lrd(i)
   end do
    
end subroutine nt_iccg
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_nsor ( phi, xs, wsor, accept, npbopt )

!  use pbtimer_module

! Passed variables
    
   integer npbopt
   double precision phi(xmymzm), xs(1-xmym:xmymzm+xmym), accept
   _REAL_ wsor

! local variables
   integer ii
   _REAL_ wsor1

!  calculate initial norm

   inorm = sum( abs( bv(1:xmymzm) ) )
!  write(6, *) 'itn & norm ', 0, inorm

!  the main loop

   wsor1 = ONE - wsor
   uconvg = .true. 
   itn = 0

   if ( npbopt == 0 ) then 
      zv(1:xmymzm) = ONE/ad(1:xmymzm)
   else
      zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv ; 
!     maxbf = maxval(abs(zv(1:xmymzm)))
      ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
      zv(1:xmymzm) = ad1(1:xmymzm) * xs(1:xmymzm) + ad(1:xmymzm)
      zv(1:xmymzm) = xs(1:xmymzm)/zv(1:xmymzm)
   end if

   do while ( uconvg )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     start the sor iteration ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Nicholas without division by xs
      if (npbopt == 0) then
         do ii = 1, xmymzm
            xs(ii) = wsor1*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
         end do
      else
         do ii = 1, xmymzm
            xs(ii) = wsor1*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
         end do
      end if

      itn = itn + 1

!     set up inverse diagonal elements, with salt if possible

      if ( npbopt /= 0 ) then
      !if (itn == 1) call pbtimer_start(PBTIME_SINH)

         !Nicholas without division by xs
         zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv 
!        maxbf = maxval(abs(zv(1:xmymzm)))
         ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
         zv(1:xmymzm) = ad1(1:xmymzm) * xs(1:xmymzm) + ad(1:xmymzm)
         zv(1:xmymzm) = xs(1:xmymzm)/zv(1:xmymzm)

      !if (itn == 1) call pbtimer_stop(PBTIME_SINH)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     check convergence every ten steps ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      if ( mod(itn, itn_checknorm) == 0 ) then
         !Nicholas without division by xs
         if (npbopt == 0) then
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) -  ad(ii)* xs(ii     )
            end do
         else
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
            end do
         end if

         norm = sum(abs(phi(1:xmymzm)))

!        write(6, *) 'itn, norm', itn, norm
!        if (npbopt /= 0) write(6, *) 'max boltzmann factor', maxbf

         if ( itn >= maxitn .or. norm <= accept*inorm ) then
            uconvg = .false.
            if ( itn .ge. maxitn ) write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
         end if
      end if

   end do

!  write(6, *) 'itn & norm ', itn, norm

end subroutine pb_nsor
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_asor ( phi, xs, wsor, accept, npbopt )

!  use pbtimer_module

!  Passed variables
    
   integer npbopt
   double precision phi(xmymzm), xs(1-xmym:xmymzm+xmym), accept
   _REAL_ wsor

!  local variables
   integer ii
   _REAL_ wsor1

!  calculate initial norm

   inorm = sum( abs( bv(1:xmymzm) ) )
!  write(6, *) 'itn & norm ', 0, inorm

!  the main loop

   wsor1 = ONE - wsor
   uconvg = .true. 
   itn = 0
   if (npbopt == 0 ) then
      zv(1:xmymzm) = ONE/ad(1:xmymzm)
   else
      zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv 
!     maxbf = maxval(abs(zv(1:xmymzm)))
      ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
      zv(1:xmymzm) = ad1(1:xmymzm) * xs(1:xmymzm) + ad(1:xmymzm)
      zv(1:xmymzm) = xs(1:xmymzm)/zv(1:xmymzm)
   end if

   do while ( uconvg )
      itn = itn + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     start the sor iteration ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Nicholas without division by xs
      if (npbopt == 0) then
         do ii = 1, xmymzm
            xs(ii) = wsor1*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
         end do
      else
         do ii = 1, xmymzm
            xs(ii) = wsor1*xs(ii) + wsor * (am1(ii-1   ) * xs(ii-1   ) + &
                                            am1(ii     ) * xs(ii+1   ) + &
                                            am2(ii-xm  ) * xs(ii-xm  ) + &
                                            am2(ii     ) * xs(ii+xm  ) + &
                                            am3(ii-xmym) * xs(ii-xmym) + &
                                            am3(ii     ) * xs(ii+xmym) + &
                                             bv(ii     )             ) * zv(ii)
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     set up inverse diagonal elements, with salt if possible
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( npbopt /= 0 ) then
      !if (itn == 1) call pbtimer_start(PBTIME_SINH)

         !Nicholas without division by xs
         zv(1:xmymzm) = iv(1:xmymzm) * xs(1:xmymzm) * ktinv 
!        maxbf = maxval(abs(zv(1:xmymzm)))
         ad(1:xmymzm) = factor * sinh( zv(1:xmymzm) )
         zv(1:xmymzm) = ad1(1:xmymzm) * xs(1:xmymzm) + ad(1:xmymzm)
         zv(1:xmymzm) = xs(1:xmymzm)/zv(1:xmymzm)

      !if (itn == 1) call pbtimer_stop(PBTIME_SINH)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     check convergence every ten steps ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      if ( mod(itn, itn_checknorm) == 0 ) then

         !Nicholas without division by xs
         if (npbopt == 0) then
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) -  ad(ii)* xs(ii     )
            end do
         else
            do ii = 1, xmymzm
               phi(ii) = am1(ii-1   ) * xs(ii-1   ) + am1(ii)* xs(ii+1   ) + &
                         am2(ii-xm  ) * xs(ii-xm  ) + am2(ii)* xs(ii+xm  ) + &
                         am3(ii-xmym) * xs(ii-xmym) + am3(ii)* xs(ii+xmym) + &
                                             bv(ii) - ad1(ii)* xs(ii     ) - ad(ii)
            end do
         end if

         norm = sum(abs(phi(1:xmymzm)))
!        write(6, *) 'itn, wsor, norm', itn, wsor, norm
!        if (npbopt /= 0) write(6, *) '   max boltzmann factor', maxbf

         if ( npbopt /= 0 .and. itn > itn_checknorm .and. wsor > 1.0d0 ) then
            if ( norm > savnorm ) then
               wsor = wsor - 0.30d0
               wsor1 = ONE - wsor
            end if
         end if
         savnorm = norm  
   
         if ( itn >= maxitn .or. norm <= accept*inorm ) then
            uconvg = .false.
            if ( itn .ge. maxitn ) write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
         end if

      end if

   end do
!  write(6, *) 'itn & norm ', itn, norm

end subroutine pb_asor
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_ionene(nfocus,level,savgox,savgoy,savgoz, &
                     savxm,savym,savzm,savh,phi,ionene)

   integer nfocus, level
   integer savxm(nfocus), savym(nfocus), savzm(nfocus)
   _REAL_ savgox(nfocus), savgoy(nfocus), savgoz(nfocus), savh(nfocus)
   _REAL_ phi(xmymzm)
   _REAL_ ionene

   ad (1:xmymzm) = savh(level)*ad(1:xmymzm)
   ad1(1:xmymzm) = cosh(iv(1:xmymzm)*phi(1:xmymzm)*ktinv)
   factor = savh(level)*factor/ktinv
   if ( level == nfocus ) then 
      ionene = ionene - sum(ad(1:xmymzm)*phi(1:xmymzm) + &
                           (ad1(1:xmymzm)-iv(1:xmymzm))*factor)
   else
      call ion_ene(nfocus,level,savgox,savgoy,savgoz, &
                  savxm,savym,savzm,savh,phi,ad,ad1,iv,ionene)
   end if

end subroutine pb_ionene
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ion_ene(nfocus,level,savgox,savgoy,savgoz, &
                     savxm,savym,savzm,savh,phi,ad,ad1,iv,ionene)

   integer nfocus, level
   integer savxm(nfocus), savym(nfocus), savzm(nfocus)
   _REAL_ savgox(nfocus), savgoy(nfocus), savgoz(nfocus), savh(nfocus)
   _REAL_ phi(xm,ym,zm), ad(xm,ym,zm), ad1(xm,ym,zm), iv(xm,ym,zm)
   _REAL_ ionene

   integer i, j, k
   integer ix1, ix2, iy1, iy2, iz1, iz2
   _REAL_  aa1, aa2, bb1, bb2, cc1, cc2
   _REAL_  h1, h2, tmp, tmp1, tmp2, tmp3, tmp4
  
!  find in which coarse grid point the edge of fine grid lies

   h1 = savh(level); h2 = savh(level+1)

   tmp1 = (savgox(level+1) + (h2-h1)*HALF - savgox(level))/h1
!  print *, tmp1
   ix1  = floor(tmp1)
   aa1  = real(ix1+1) - tmp1

   tmp2 = (savgoy(level+1) + (h2-h1)*HALF - savgoy(level))/h1
!  print *, tmp2
   iy1  = floor(tmp2)
   bb1  = real(iy1+1) - tmp2

   tmp3 = (savgoz(level+1) + (h2-h1)*HALF - savgoz(level))/h1
!  print *, tmp3
   iz1  = floor(tmp3)
   cc1  = real(iz1+1) - tmp3

   tmp1 = tmp1 + real(savxm(level+1))*h2/h1
!  print *, tmp1
   ix2  = floor(tmp1)
   aa2  = tmp1 - real(ix2)

   tmp2 = tmp2 + real(savym(level+1))*h2/h1
!  print *, tmp2
   iy2  = floor(tmp2)
   bb2  = tmp2 - real(iy2)

   tmp3 = tmp3 + real(savzm(level+1))*h2/h1
!  print *, tmp3
   iz2  = floor(tmp3)
   cc2  = tmp3 - real(iz2)

!  print *, ix1, ix2, iy1, iy2, iz1, iz2
!  print *, aa1, aa2, bb1, bb2, cc1, cc2

!  do i = 1, xm; do j = 1, ym; do k = 1, zm
!     write(101,*) phi(i,j,k)
!  end do; end do; end do

!  the grid points inside

   do i = 1, savxm(level)
      if ( i >= ix1 .or. i <= ix2 ) cycle
      do j = 1, savym(level)
         if ( j >= iy1 .or. j <= iy2 ) cycle
         do k = 1, savzm(level)
            if ( k >= iz1 .or. k <= iz2 ) cycle
            ionene = ionene - ad(i,j,k)*phi(i,j,k) - &
                     (ad1(i,j,k)-iv(i,j,k))*factor                
         end do
      end do
   end do
         
!  the grid points on the six faces

   do j = iy1+1, iy2-1; do k = iz1+1, iz2-1
      i = ix1
      ionene = ionene - (ONE-aa1)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix2
      ionene = ionene - (ONE-aa2)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do; end do
   do i = ix1+1, ix2-1; do k = iz1+1, iz2-1
      j = iy1
      ionene = ionene - (ONE-bb1)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
      j = iy2
      ionene = ionene - (ONE-bb2)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do; end do
   do i = ix1+1, ix2-1; do j = iy1+1, iy2-1
      k = iz1
      ionene = ionene - (ONE-cc1)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
      k = iz2
      ionene = ionene - (ONE-cc2)*(ad(i,j,k)*phi(i,j,k) + &
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do; end do

!  the grid points on the twelve edges

   tmp1=ONE-bb1*cc1; tmp2=ONE-bb1*cc2; tmp3=ONE-bb2*cc1; tmp4=ONE-bb2*cc2
   do i = ix1+1, ix2-1
      j = iy1; k = iz1
      ionene = ionene - tmp1*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      j = iy1; k = iz2
      ionene = ionene - tmp2*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      j = iy2; k = iz1
      ionene = ionene - tmp3*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      j = iy2; k = iz2
      ionene = ionene - tmp4*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do
   tmp1=ONE-aa1*cc1; tmp2=ONE-aa1*cc2; tmp3=ONE-aa2*cc1; tmp4=ONE-aa2*cc2
   do j = iy1+1, iy2-1
      i = ix1; k = iz1
      ionene = ionene - tmp1*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix1; k = iz2
      ionene = ionene - tmp2*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix2; k = iz1
      ionene = ionene - tmp3*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix2; k = iz2
      ionene = ionene - tmp4*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do
   tmp1=ONE-aa1*bb1; tmp2=ONE-aa1*bb2; tmp3=ONE-aa2*bb1; tmp4=ONE-aa2*bb2
   do k = iz1+1, iz2-1
      i = ix1; j = iy1
      ionene = ionene - tmp1*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix1; j = iy2
      ionene = ionene - tmp2*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix2; j = iy1
      ionene = ionene - tmp3*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
      i = ix2; j = iy2
      ionene = ionene - tmp4*(ad(i,j,k)*phi(i,j,k) + & 
               (ad1(i,j,k)-iv(i,j,k))*factor)
   end do

!  if the fine grid spans only one coarse grid point in x, y, or z direction

   if ( ix1 == ix2 ) then
      aa1 = aa1 - HALF
      aa2 = aa2 - HALF
   end if
   if ( iy1 == iy2 ) then
      aa1 = aa1 - HALF
      aa2 = aa2 - HALF
   end if
   if ( iz1 == iz2 ) then
      aa1 = aa1 - HALF
      aa2 = aa2 - HALF
   end if

!  the grid points at eight corners

   i=ix1; j=iy1; k=iz1; tmp=ONE-aa1*bb1*cc1
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix1; j=iy1; k=iz2; tmp=ONE-aa1*bb1*cc2
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix1; j=iy2; k=iz1; tmp=ONE-aa1*bb2*cc1
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix1; j=iy2; k=iz2; tmp=ONE-aa1*bb2*cc2
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix2; j=iy1; k=iz1; tmp=ONE-aa2*bb1*cc1
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix2; j=iy1; k=iz2; tmp=ONE-aa2*bb1*cc2
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix2; j=iy2; k=iz1; tmp=ONE-aa2*bb2*cc1
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)
   i=ix2; j=iy2; k=iz2; tmp=ONE-aa2*bb2*cc2
   ionene = ionene - tmp*(ad(i,j,k)*phi(i,j,k) + &
            (ad1(i,j,k)-iv(i,j,k))*factor)

end subroutine ion_ene

end module pb_nlsolver
