! <compile=optimized>
#include "copyright.h"
#  define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ IIM algoritm driver
subroutine pb_iimdrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs&
                    )

   implicit none

#include "pb_constants.h"

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.

   ! passed variables

   integer npbstep, npbgrid, nstlim, atmfirst, atmlast,natom
   integer npbopt, solvopt, level, nfocus, bcopt
   _REAL_  h, savh(nfocus), gox, goy, goz, savgox(nfocus), savgoy(nfocus), savgoz(nfocus)
   integer xm, ym, zm, xmym, xmymzm, savxm(nfocus), savym(nfocus), savzm(nfocus)
   integer maxitn, itn
   _REAL_  fmiccg, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_  pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout, ionene
   _REAL_  gcrd(3,natom), acrg(natom)
   integer nbnd, iepsav(4,xmymzm)

   integer insas(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)
   _REAL_  lvlset(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)

   _REAL_  chgrd(xmymzm), saltgrd(xmymzm), phi(xmymzm)
   _REAL_  bv(xmymzm), cphi(xmymzm), xs(xmymzm+2*xmym)
   integer ii

   ! local varialbes

   _REAL_ rh

   rh = ONE/h

   ! for singularity-free PB equation we need to compute cphi() at the
   ! dielectric interface grid points

   if ( bcopt > 5 .and. bcopt < 10 ) then
      cphi(1:xmymzm) = ZERO
      call pb_dbcgrd( cphi(1), insas )
   end if

   ! for singular PB equation bv() is initialized as the grid charges

   if ( bcopt < 6 ) then
      bv = chgrd
   end if

   ! now we put effective charges at the space boundary grid points into bv()
   ! to take care of the boundary conditions

   call pb_bndcnd( bv(1), chgrd(1) )

   call IIM(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1),phi(1),xs(1),accept,h,atmfirst,atmlast,bcopt,solvopt)

contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbcgrd( cphi, insas )

   _REAL_ cphi(xm,ym,zm)
   integer insas(0:xm+1,0:ym+1,0:zm+1)

   integer i, j, k
   integer i0, j0, k0
   integer ip
   _REAL_ tmp

   tmp = ZERO

   do ip = 1, nbnd
      i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)

      i = i0 - 1
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      i = i0 + 1
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      j = j0 - 1
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp);
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      j = j0 + 1
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp)
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      k = k0 - 1
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      k = k0 + 1
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      if ( insas(i0,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k0) == ZERO ) then
            call get_coulpot(i0,j0,k0,tmp)
            cphi(i0,j0,k0) = tmp/epsin
         end if
      end if

   end do

end subroutine pb_dbcgrd
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulpot(i,j,k,pot)

   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz

   _REAL_ factor,qtmp,rinv,xtmp,ytmp,ztmp,dx,dy,dz
   _REAL_ a,a1,b,b1,c,c1

   factor = ONE/(FOURPI)/h

   pot = ZERO
   do iatm = atmfirst, atmlast
      xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
      qtmp = factor*acrg(iatm)
           
      dx = abs(i-xtmp); dy = abs(j-ytmp); dz = abs(k-ztmp)
      if (dx < 20.d0 .and. dy < 20.d0 .and. dz < 20.d0) then
         idx = floor(dx); idy = floor(dy); idz = floor(dz)
         a=dx-idx;b=dy-idy;c=dz-idz
         a1 = 1 - a; b1 = 1 - b; c1 = 1 - c
         rinv = a1*b1*c1*green(idx  ,idy  ,idz  ) &
               +a *b1*c1*green(idx+1,idy  ,idz  ) &
               +a1*b *c1*green(idx  ,idy+1,idz  ) &
               +a *b *c1*green(idx+1,idy+1,idz  ) &
               +a1*b1*c *green(idx  ,idy  ,idz+1) &
               +a *b1*c *green(idx+1,idy  ,idz+1) &
               +a1*b *c *green(idx  ,idy+1,idz+1) &
               +a *b *c *green(idx+1,idy+1,idz+1)
      else
         rinv = ONE/sqrt(dble(dx**2 + dy**2 + dz**2))
      end if
      pot = pot + qtmp*rinv
   end do  !  iatm = atmfirst, atmlast

end subroutine get_coulpot
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv, chgrd )
   
   ! Common variables
    
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
     
   ! Passed variables
    
   _REAL_ bv(xm,ym,zm), chgrd(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, iatm, ii
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp, factor
   _REAL_ x, y, z, dx, dy, dz, xtmp, ytmp, ztmp
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv

   ! part a: level = 1 cases ::::: 
   ! bcopt = 1
   ! zero potential in the singular PB.
   ! the boundary will be all solvent
    bv=0.0d0 !XP: This can make imin=6 work with ipb4
   if ( level == 1 .and. bcopt == 1 ) then
      write(6, *) "PB bomb in pb_bndcnd(): zero BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 2
   ! molecule dipolar debye-huckel potential in the singular PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 2 ) then
      write(6, *) "PB bomb in pb_iimdrv(): molecular dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 3
   ! sum of residue dipolar debye-huckel potentials in the singular PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 3 ) then
      write(6, *) "PB bomb in pb_iimdrv(): residue dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 4 .or. bcopt = 6
   ! sum of atom charge deby-huckel potentials in the singular (4) or singularity-free (6) PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. ( bcopt == 4 .or. bcopt == 6 ) ) then
      do iatm = atmfirst, atmlast
         xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
         qtmp = INV_FOURPI*acrg(iatm)
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            dx = i-xtmp; dy = j-ytmp; dz = ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
             
            dz = zm+1-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
          
         ! j=0 and ym+1 faces
          
         do k = 1, zm; do i = 1, xm
            dx = i-xtmp; dy  = ytmp; dz  = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
             
            dy = ym+1-ytmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
          
         ! i=0 and i=xm+1 faces
          
         do k = 1, zm; do j = 1, ym
            dx = xtmp; dy = j-ytmp; dz = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
             
            dx = xm+1-xtmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
      end do  !  iatm = atmfirst, atmlast
    
   ! bcopt = 5 .or. bcopt = 7
   ! sum of grid charge debye-huckel potentials in the singular (5) or singularity-free (7) PB.
   ! the boundary will be all solvent.

   else if ( level == 1 .and. (bcopt == 5 .or. bcopt == 7) ) then
       
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle

         qtmp = INV_FOURPI*qtmp 

         ! k=0 and k=zm+1 faces

         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do

      end do; end do; end do  !  itmp = 1, xm; jtmp = 1, ym; ktmp = 1, zm
       
   ! bcopt = 8
   ! sum of atom charge reaction field potentials in the singularity-free PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 8 ) then
       
      factor = INV_FOURPI*epsout*(ONE/epsout - ONE/epsin)
      do iatm = atmfirst, atmlast
         xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
         qtmp = factor*acrg(iatm)
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            dx = i-xtmp; dy = j-ytmp; dz = ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,1 ) = bv(i,j,1 ) + qtmp/r
             
            dz = zm+1-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,zm) = bv(i,j,zm) + qtmp/r
         end do; end do
          
         ! j=0 and ym+1 faces
          
         do k = 1, zm; do i = 1, xm
            dx = i-xtmp; dy  = ytmp; dz  = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,1 ,k) = bv(i,1 ,k) + qtmp/r
             
            dy = ym+1-ytmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,ym,k) = bv(i,ym,k) + qtmp/r
         end do; end do
          
         ! i=0 and i=xm+1 faces
          
         do k = 1, zm; do j = 1, ym
            dx = xtmp; dy = j-ytmp; dz = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(1 ,j,k) = bv(1 ,j,k) + qtmp/r
             
            dx = xm+1-xtmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(xm,j,k) = bv(xm,j,k) + qtmp/r
         end do; end do
      end do  !  iatm = atmfirst, atmlast
    
   ! bcopt = 9
   ! sum of grid charge reaction field potentials in the singularity-free PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 9 ) then
       
      factor = INV_FOURPI*epsout*(ONE/epsout - ONE/epsin)
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle

         qtmp = factor*qtmp
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + qtmp/r
            end if
         end do; end do

      end do; end do; end do  !  itmp = 1, xm; jtmp = 1, ym; ktmp = 1, zm

   ! part b: level > 1 case 
   ! electrostatic focusing
    
   else if ( level > 1 ) then
      xmtmp  = savxm(level-1) ; ymtmp  = savym(level-1) ; zmtmp  = savzm(level-1)
      htmp   = savh(level-1)
      goxtmp = savgox(level-1); goytmp = savgoy(level-1); goztmp = savgoz(level-1)
       
      ! k=0 and k=zm+1 faces
       
      do j = 1, ym; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy + h*j        ; z  = goz
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(i,j,1 ) = bv(i,j,1 ) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         z  = goz + h*(zm+1)
         zi = (z - goztmp)/htmp
         iz = int( zi )
         cc  = zi - dble( iz )
         cc1 = ONE - cc
         bv(i,j,zm) = bv(i,j,zm) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
   
      ! j=0 and j=ym+1 faces
   
      do k = 1, zm; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy              ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa = xi - dble( ix ); bb = yi - dble( iy ); cc = zi - dble( iz )
         aa1 = ONE - aa      ; bb1 = ONE - bb      ; cc1 = ONE - cc
         bv(i,1 ,k) = bv(i,1 ,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         y  = goy + h*(ym+1)
         yi = (y - goytmp)/htmp
         iy = int( yi )
         bb  = yi - dble( iy )
         bb1 = ONE - bb
         bv(i,ym,k) = bv(i,ym,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
        
      ! i=0 and i=xm+1 faces
        
      do k = 1, zm; do j = 1, ym
          
         x  = gox              ; y  = goy + h*j        ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(1 ,j,k) = bv(1 ,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         x  = gox + h * (xm+1)
         xi = (x - goxtmp)/htmp
         ix = int( xi )
         aa  = xi - dble( ix )
         aa1 = ONE - aa
         bv(xm,j,k) = bv(xm,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
       
   else
       
      ! unknown bcopt
       
      write(6, *) 'PB bomb in pb_iimdrv(): unknown BC option', bcopt
      call mexit(6, 1)
   end if 
 
 
end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ phi interpretation, xm, ym, zm are for the previous phi map
_REAL_ function phintp( xmtmp,ymtmp,zmtmp,ix,iy,iz,aa,bb,cc,aa1,bb1,cc1 )
    
   ! Passed variables
    
   integer, intent(in) :: xmtmp, ymtmp, zmtmp, ix, iy, iz
   _REAL_, intent(in) :: aa, bb, cc, aa1, bb1, cc1
    
   ! Local Variables
    
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! determine the position of the point w.r.t. the map
    
   bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1; bb1cc  = bb1*cc ; bb_cc  = bb *cc

   ! triliner interpolation
    
   phintp = aa1*bb1cc1*phi( ix   + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa *bb1cc1*phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa1*bb_cc1*phi( ix   + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa *bb_cc1*phi( ix+1 + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa1*bb1cc *phi( ix   + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa *bb1cc *phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa1*bb_cc *phi( ix   + xmtmp*( iy   + ymtmp*( iz   ) ) ) + &
            aa *bb_cc *phi( ix+1 + xmtmp*( iy   + ymtmp*( iz   ) ) )
                
end function phintp


end subroutine pb_iimdrv
