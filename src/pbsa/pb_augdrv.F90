! <compile=optimized>
#include "copyright.h"
#  define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ IIM algoritm driver
subroutine pb_augdrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs,&
                      pbverbose, &
                      augtoltype, augctf, augtol &
                    )

   implicit none

#include "pb_constants.h"

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.

   ! passed variables
   logical pbverbose
   integer npbstep, npbgrid, nstlim, atmfirst, atmlast
   integer npbopt, solvopt, level, nfocus, bcopt
   ! level, nfocus presumably is always 1
   integer natom
   integer maxitn, itn
   integer xm, ym, zm, xmym, xmymzm
   integer savxm(nfocus), savym(nfocus), savzm(nfocus)
   integer nbnd, iepsav(4,xmymzm)
   ! iepsav is allocated as iepsav(4,xmymzm_max)
   _REAL_  fmiccg, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_  pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout, ionene
   _REAL_  h, savh(nfocus), gox, goy, goz, savgox(nfocus), savgoy(nfocus), savgoz(nfocus)
   _REAL_  gcrd(3,natom), acrg(natom)
   !_REAL_  qgrdcrg(*)
   !   _REAL_,allocatable :: cirreg(:,:),ct(:,:),cepsav(:,:)
   !_REAL_,allocatable :: x(:),y(:),z(:),lvlst(:,:,:)
   !integer,allocatable :: index(:,:,:),index2(:,:,:)

   integer insas(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)
   _REAL_  lvlset(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)

   _REAL_  chgrd(xmymzm), saltgrd(xmymzm), phi(xmymzm)
   _REAL_  bv(xmymzm), cphi(xmymzm), xs(xmymzm+2*xmym)
   integer ii

   integer augtoltype !WMBS 0 for absolute 1 for relative
   _REAL_  augctf !WMBS cluster radius factor <= 0 -> h*h,
                  !     > 0 -> h*augctf
   _REAL_  augtol !WMBS tolerance for gmres routine

   ! local varialbes

   _REAL_ rh
  
   ! for cluster cutoff
 
! _REAL_ ctf,dist
! integer ctn,icn,i,j,k
  
   rh = ONE/h

   ! for singularity-free PB equation we need to compute cphi() at the
   ! dielectric interface grid points

   if ( bcopt > 5 .and. bcopt < 10 ) then
      cphi(1:xmymzm) = ZERO
      call pb_dbcgrd( cphi(1), insas )
   end if

   ! for singular PB equation bv() is initialized as the grid charges

   if ( bcopt < 6 .or. bcopt > 9 ) then
      bv = chgrd
   end if

   ! now we put effective charges at the space boundary grid points into bv()
   ! to take care of the boundary conditions

 
  !!! needs new nbnd and iepsav here before calling aug
  !  variables: integer ctn,icn XP: attempt to rule out bad irreg points
 
!  allocate (cirreg(1:nbnd,1:15),ct(1:nbnd,1:3),cepsav(1:4,1:xmymzm))
!  allocate (x(0:xm+1), y(0:ym+1), z(0:zm+1))
!  allocate (index(1:xm,1:ym,1:zm),index2(1:xm,1:ym,1:zm))
!  allocate (lvlst(0:xm+1,0:ym+1,0:zm+1))
!           
!  call  set_index(xm,ym,zm,h,gox,goy,goz,insas,lvlset,index)
!  
!  

!  do i = 0, xm+1
!      x(i) = gox + i*h
!  end do
!  do j = 0, ym+1
!     y(j) = goy + j*h
!  end do
!  do k = 0, zm+1
!     z(k) = goz + k*h
!  end do
!  
!
!  call pre_indexg(nbnd,x,y,z,index,lvlset,index2,cirreg)
!  do ii=1,nbnd
!  write(66,*) cirreg(ii,1),cirreg(ii,2),cirreg(ii,3)
!  end do 
!   stop
!  ctf=0.25d0;ctn=0

!  do ii=1,nbnd
!   write(6,*) "iter ii=",ii
!   if( ii==1) then
!    ct(1,1:3)=cirreg(ii,1:3)
!    cepsav(1:3,1)=iepsav(1:3,ii)
!     ctn=ctn+1
!   else 
!   do icn=1,ctn
!    dist=sqrt((ct(icn,1)-cirreg(ii,1))**2+(ct(icn,2)-cirreg(ii,2))**2+(ct(icn,3)-cirreg(ii,3))**2)
!    write(6,*) "ct",ct(icn,1:3)
!    write(6,*) "cirreg",cirreg(ii,1:3)
!    write(6,*) "icn",icn,"dist",dist
!    if(dist <= ctf) exit
!    if(icn==ctn) then
!     ctn=ctn+1
!    ct(ctn,1:3)=cirreg(ii,1:3)
!    cepsav(1:3,ctn)=iepsav(1:3,ii)
!    end if
!   end do
!    end if
!  write(6,*) "ctn=",ctn
!  end do
!  
!  do ii=1,ctn
!  write(99,*) cepsav(1:3,ii)
!  end do 
! 
! 
!    
!   write(6,*) 'nbnd=',nbnd,'ctn=',ctn
!     stop
   

   call pb_bndcnd( bv(1), chgrd(1) )
   call aug(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1),phi(1),accept,h,atmfirst,atmlast,bcopt,solvopt,&
            augtoltype,augctf,augtol)

!  deallocate (cirreg,ct,cepsav,index,index2)

!  deallocate (x, y, z)
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  subroutine pre_indexg(nbnd,x,y,z,index,lvlset,index2,cirreg)

!!XP: below attenpts to rule out bad irreg points,now set that job to aug.f
!   
!   integer nbnd,index(1:xm,1:ym,1:zm),index2(1:xm,1:ym,1:zm)
!   _REAL_  lvlset(0:xm+1,0:ym+1,0:zm+1)
!   _REAL_  cirreg(1:nbnd,1:15)
!   _REAL_  x(:),y(:),z(:)
!   !_REAL_  x(0,xm+1),y(0,ym+1),z(0,zm+1)
!   
!  lvlst(1:xm,1:ym,1:zm) = lvlset(1:xm,1:ym,1:zm)
! write(8001,*) nbnd,x
! write(8002,*) nbnd,y
! write(8003,*) nbnd,z
! write(8004,*) index
! write(8005,*) lvlst
! call  indexg(nbnd,x,y,z,index,lvlst,index2,cirreg)
!  do ii=1,nbnd
!  write(55,*) cirreg(ii,1),cirreg(ii,2),cirreg(ii,3)
!  end do 
!  
!   end subroutine pre_indexg

!ubroutine  set_index(xm,ym,zm,h,gox,goy,goz,insas,lvlset,index)
!  integer insas(0:xm+1,0:ym+1,0:zm+1)
!  _REAL_  lvlset(0:xm+1,0:ym+1,0:zm+1)
!  integer index(xm,ym,zm)
!  integer i,j,k,xm,ym,zm
!  _REAL_ h,gox,goy,goz
!  
!  do k = 1,zm 
!     do j = 1, ym
!        do i = 1, xm
!           if ( insas(i,j,k) > 0 ) then
!              index(i,j,k) = 1
!           else
!              index(i,j,k) = 5
!           end if
!        end do
!     end do
!  end do
!  do ii = 1, nbnd
!     i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)
!     if ( lvlset(i,j,k) > 0 ) then
!        index(i,j,k) = 4
!     else if ( lvlset(i,j,k) < 0 ) then
!        index(i,j,k) = 2
!     else
!        index(i,j,k) = 3
!     end if
!  end do
! 
!end subroutine set_index

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
   if (.not. ( (bcopt == 10) .and. (solvopt == 7 &
        .or. solvopt == 3))) then
           bv=0.0d0 !XP: This can make imin=6 work with ipb4 5 
        !WMBS: This will cause problems for fft or pcg solver
   end if

   if ( level == 1 .and. bcopt == 10 ) then
      if (.not. (solvopt == 7 .or. solvopt == 3)) then  
        write(6, *) "PB bomb in pb_bndcnd(): zero BC only allowed for fft or pcg"
        call mexit(6, 1)
      end if
    
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
!aug
subroutine aug(gox,goy,goz,l,m,n,phi,insas,nbnd,iepsav,epsin,epsout,bv,u,accept,h,&
               atmfirst,atmlast,bcopt,solvopt,augtoltype,augctf,augtol)
   ! the AUG method used to solve the PBE with different delectric constants of
   ! two regions separated by the interface   

   !    xs <- gox, ys <- goy, zs <- goz
   !    l  <- xm , m  <- ym , n  <- zm
   !    bi <- epsin, bo <- epsout
   !    u  <- phi,u0 <- xs

   implicit none

   ! passed variables
  
   !  go[x-z]                starting crds of the box r
   !  l,m,n                  dimension of the 3-D box
   !  phi                    level set function   
   !  insas                  index of grids to differentiate inner and outer
   !                         sider, the defination varies to the sasopt
   !  nbnd                   number of irregular grid points   
   !  iepsav                 array storing the grid crds of irregular ptns, it
   !                         originally also store the number of atom the grid
   !                         point belongs to, yet it does not apply the AUG
   !  epsin,epsout           dielectric constants of inside and outside 
   !  bv                     initial rhs of the linear eqns of the whole grids
   !  u                      potential of whole grids
   !  accept                 convergence criteria of the linear solver
   !  h                      grid spacing
   !  atmfirst,atmlast       the first and last atm index

   _REAL_  gox,goy,goz,h,accept,epsout,epsin,tol,dist
   integer l,m,n,nbnd,atmfirst,atmlast,bcopt,solvopt,iter
   integer imax,mm
   integer insas(0:l+1,0:m+1,0:n+1),iepsav(1:4,1:l*m*n)
   _REAL_  bv(l,m,n),u(l,m,n)
   _REAL_  phi(0:l+1,0:m+1,0:n+1)
   _REAL_  proj(1:1432,1:3)

   integer augtoltype ! 0 absolute, 1 relative
   _REAL_  augctf !cluster radius factor
   _REAL_  augtol !gmres tolerance
    
#include "pb_constants.h"

   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo,error
   integer  nirreg,nind3

   ! local variables
   ! x(),y(),z()             coordinates of grid points
   ! cirreg cirreg1          geometrical parameters at irregular points
   ! index()                 flag of all grid points with 1-5
   ! index2()                number index for irregular grid points
   ! index3()                number index for inner side irregular points
   ! wp,qp,wp1,qp1           [u]and[betaUn]for irregular and inner irregular points           
   ! unj(),unjf()            [Un] for irregular points
   ! bf(), fvec()            rhs of Ag=b,output of matrix-vector multiply
   ! w0p,wyp,wzp,w0p,        tangential derivatives of jump conditions
   ! wyp,wzp                
   ! wcoe,wxcoe,wycoe,wzcoe  coefficients obtained by SVD
   ! wxxcoe,wyycoe wzzcoe
   ! wxycoe,wxzcoe,wyzcoe
   ! c,c2                    correction terms
   ! sss1,sss2               dummy arrays
    

   _REAL_  ctf
   integer ctn,icn
   integer it,jt,kt,iit,it1
   integer i1,j1,k1,flag,i2,j2,k2
   _REAL_, parameter :: eps0 = 8.8542D-12 / (1.6022D-19)**2 /  &
                               (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer, parameter :: nq=27
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ rhs
   _REAL_,allocatable :: x(:),y(:),z(:),cirreg(:,:),cirreg1(:,:),ct(:,:)! geometrical parameters at irregular points
   integer,allocatable :: index(:,:,:),index2(:,:,:)!,index3(:,:,:)
   _REAL_,allocatable :: wp(:),qp(:),unj(:),unjf(:),bf(:),fvec(:) 
   _REAL_,allocatable :: q0p(:),qyp(:),qzp(:)
   _REAL_,allocatable :: w0p(:),wyp(:),wzp(:)
   _REAL_,allocatable :: wyyp(:),wzzp(:),wyzp(:)
   _REAL_,allocatable :: wcoe(:,:),wxcoe(:, :),wycoe(:, :)
   _REAL_,allocatable :: wzcoe(:,:),wxxcoe(:, :),wyycoe(:, :)
   _REAL_,allocatable :: wzzcoe(:, :),wxycoe(:, :)
   _REAL_,allocatable :: wxzcoe(:, :),wyzcoe(:, :)
   _REAL_,allocatable :: sss1(:),sss2(:)
   _REAL_,allocatable :: c(:,:,:,:), c2(:, :)
   _REAL_,allocatable :: f(:,:,:)
   !_REAL_,allocatable :: hg(:,:),v(:,:)
   character(10) str
   integer icall 

   character*11 potentialdxfilename
   character*9  potentialdxdataname
   integer      potentialdxfilenum 

   logical use_average_tolerance

   if (augtoltype == 1) then
      use_average_tolerance = .true.
   else
      use_average_tolerance = .false.
   endif

   potentialdxfilename= "pbsa_phi.dx"
   potentialdxdataname= "Potential"
   potentialdxfilenum= 90885

   ! interface variables 

   nirreg = nbnd
   xs = gox; ys = goy; zs = goz
   hx = h; hy = h; hz = h; hmax = h


   bi = epsin; bo = epsout

   ! allocating working arrays

   allocate (x(0:l+1), y(0:m+1), z(0:n+1))
   allocate (index(1:l,1:m,1:n), index2(1:l,1:m,1:n))
   allocate (cirreg(1:nbnd,1:15),cirreg1(1:nbnd,1:15),ct(1:nbnd,1:3))
   allocate (wp(nbnd),qp(nbnd))
   allocate (unj(nbnd),unjf(nbnd),bf(nbnd),fvec(nbnd))
   allocate (q0p(nbnd),qyp(nbnd),qzp(nbnd)) 
   allocate (w0p(nbnd),wyp(nbnd),wzp(nbnd))
   allocate (wyyp(nbnd),wzzp(nbnd),wyzp(nbnd))
   allocate (wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd, nq))
   allocate (wzcoe(nbnd,nq),wxxcoe(nbnd, nq),wyycoe(nbnd, nq))
   allocate (wzzcoe(nbnd, nq),wxycoe(nbnd, nq))
   allocate (wxzcoe(nbnd, nq),wyzcoe(nbnd, nq))
   allocate (sss1(nbnd),sss2(nbnd))
   allocate (c(l,m,n,7), c2(nbnd, 27))
   allocate (f(l,m,n))

   !setting up the box
   do i = 0, l+1
       x(i) = xs + i*hx
   end do
   do j = 0, m+1
      y(j) = ys + j*hy
   end do
   do k = 0, n+1
      z(k) = zs + k*hz
   end do

   ! setting up index
   ! index = 1 for interior regular grid points;
   ! index = 2 for interior irregular grid points;
   ! index = 3 for boundary grid points right on the interface;
   ! index = 4 for exterior irregular grid points;
   ! index = 5 for exterior regular grid points;

   index=0
   do k = 1, n
      do j = 1, m
         do i = 1, l
            if ( insas(i,j,k) > 0 ) then
               index(i,j,k) = 1
            else
               index(i,j,k) = 5
            end if
         end do
      end do
   end do
   
   do ii = 1, nbnd
      i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)

      if ( phi(i,j,k) > 0 ) then
         index(i,j,k) = 4
      else if ( phi(i,j,k) < 0 ) then
         index(i,j,k) = 2
      else
         index(i,j,k) = 3
      end if
   end do

   call indexg(l,m,n,h,hx,hy,hz,nbnd,x,y,z,index,phi,index2,cirreg) 

   !XP: index2 just to fill the blank of the arrays

   !XP: Below is to calculate the 1st derivatives of lvlset  for debuging
   !    call gradu(xm,ym,zm,-ONE,27,10,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,insas)

   ! XP: start to cluster  to get the new index2 and cirreg
   ! variable: ctn, ciireg1,index3

   if ( augctf > 0 .and. augctf <= 1) then
      ctf = h*augctf
   else
      ctf=h*h             ! cut off
   endif 
   ctn=0;index2=0         ! initialization


   do ii=1,nbnd
    if( ii==1) then       !1st one
     ct(1,1:3)=cirreg(ii,1:3)
     ctn=ctn+1
     cirreg1(ctn,1:15)=cirreg(ii,1:15)
     i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
     index2(i1,j1,k1)=ctn
    else 

    do icn=1,ctn         !not the 1st one, compare to the nearest current point

     dist=sqrt((ct(icn,1)-cirreg(ii,1))**2+(ct(icn,2)-cirreg(ii,2))**2+(ct(icn,3)-cirreg(ii,3))**2)

     if(dist <= ctf) then ! fall into a current one
     i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
     index2(i1,j1,k1)=icn
       exit
     end if

     if(icn==ctn) then    ! not fall into a current one
      ctn=ctn+1
     ct(ctn,1:3)=cirreg(ii,1:3)
     cirreg1(ctn,1:15)=cirreg(ii,1:15)
     i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
      index2(i1,j1,k1)=ctn
     end if
    end do
     end if
    end do

    nirreg=ctn ! ctn the after-cluster dimension

     !XP: print out all the projection points and see their distance with others.
    !!bpcont=0;bpct4=0;bpct5=0
    !!do i=1,nbnd
    !!   do j=i+1,nbnd
    !!    dist=sqrt((cirreg1(i,1)-cirreg1(j,1))**2+(cirreg1(i,2)-cirreg1(j,2))**2+(cirreg1(i,3)-cirreg1(j,3))**2)
    !!    if( dist < h*h) then
    !!     bpcont=bpcont+1
    !!     if(dist < h*h*h*h) bpct4=bpct4+1
    !!     if(dist < h*h*h*h*h) bpct5=bpct5+1
    !!   !write(6,*) 'close points:',iepsav(1:3,i),cirreg1(i,1:3),&
    !!    !'and',iepsav(1:3,i),cirreg1(i,1:3),'dist=',dist
    !!   end if 
    !!    end do 
    !!end do 

   ! setting up mm for gmres by experience
    if ( ctn < 1000 ) then
    mm= ceiling(real(ctn)/16.0d0)
    else 
    mm= ceiling(real(ctn)/60.0d0)
    end if
   !mm= ceiling(real(ctn)/8.0d0)
   
   !print dimension and clustering results
   !write(10,*) 'ctn=',ctn,'nbnd=',nbnd,'mm=',mm
   !write(10,*) 'Clustering Reduced:',100.0d0*(nbnd-ctn)/nbnd,'%'
   !write(6,*) 'ctn=',ctn,'nbnd=',nbnd,'mm=',mm
   !write(6,*) 'Clustering Reduced:',100.0d0*(nbnd-ctn)/nbnd,'%'
   !write(6,*) "close point pair number is",bpcont,"while nbnd = ",nbnd
   !write(6,*) "h^4 points:",bpct4,"h^5 points:",bpct5,"while nbnd = ",nbnd

   !write (60, *) nbnd
   !write (59, '("DOTS")')
   !write (60, '("DOTS")')
   !do i = 1, nbnd
   !  i1 = iepsav(1,i); j1 = iepsav(2,i) ; k1 = iepsav(3,i)
   !if(index(i1,j1,k1) == 2 .or. index(i1,j1,k1) == 3) then 
   !  write(str,'(i10)') i
   !   str=adjustl(str)
   !  write(59,'(a,3(f20.15),2x)') "H"//str, cirreg(i, 1:3)
   !else if(index(i1,j1,k1) == 4) then 
   ! write(str,'(i10)') i
   !   str=adjustl(str)
   !  write(60,'(a,3(f20.15),2x)') "H"//str, cirreg(i, 1:3)
   !end if
   !end do
   ! write(109,*) cirreg

   ! Below is the former attempts to use the inside irregular point to reduce
   ! the dimension. However, it turns out the double sides would be better.

   ! setting up index3, the address index of irregular grid points inside
   ! we don't use index3 in this two side approach,it is kept for later
   !nind3=0
   !index3=0
   !do k=1,n
   !  do j=1,m
   !    do i=1,l
   !       if(index2(i,j,k)>0) then
   !          if(phi(i,j,k)<0) then
   !             nind3=nind3+1
   !             index3(i,j,k)=nind3
   !          end if
   !       end if
   !    end do
   !  end do

   ! setting up the jump conditions at the projection points of the irreular (boundary) grid points
 
   wp=0.0d0;qp=0.0d0  ! Initialization

   call prodis(l,m,n,nirreg,bcopt, bi,bo,atmfirst, atmlast,nbnd, cirreg1, wp, qp) !wp[u] qp [beta_un]

   !  do i =1, ctn
   !  proj(i,1:3)= cirreg1(i,1:3)
   !   write(1310,*) cirreg(i,1:15)
   !   write(1311,*) cirreg1(i,1:3)
   !   write(1312,*) proj(i,1:3)
   !   write(1001,*) i,wp(i)
   !   write(1002,*) i,qp(i)
   !  end do

   ! Attempt to save the 2nd copy
   !do ii = 1, nbnd
   !   i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)
   !   if ( index3(i,j,k) == 0 ) cycle
   !   wp1(index3(i,j,k)) = wp(ii)        !wp1 qp1 index for 1 side
   !   qp1(index3(i,j,k)) = qp(ii)
   !   cirreg1(index3(i,j,k),1:15) = cirreg(ii,1:15)
   !end do

   !Below is used to test the boundary grid points
   !write(19,*) "nbnd = ",nbnd
   !flag=3
   ! do ii=1,nbnd
   !  i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)
   !  if( index(i,j,k) == 2 .or. index(i,j,k) == 3) then
   !     flag = 0
   !    do i1=i-1,i+1
   !       do j1=j-1,j+1
   !          do k1=k-1,k+1
   !          if(index(i1,j1,k1) == 4) then
   !              write(014,*) x(i1),y(j1),z(k1) !print out points
   !               flag = 1
   !            end if
   !          end do
   !        end do 
   !      end do 
   !        if(flag == 0) then
   !         write(017,*) x(i),y(j),z(k)
   !              write(18,*)  i,j,k
   !       end if
   ! end if
   !end do

   ! prepare to call matvec the first time to compute F2
   ! set g=0(unjf)  and also b(bf)=0 to get F2 as in IIM_augment

   do i=1, ctn ! 
      unj(i)=0.0d0
      unjf(i)=0.0d0
      bf(i)=0.0d0
   end do                             

   icall = 1

   !WMBS-Write projection points to .xyz format file for display in VMD 
   !Flag should be added to allow user control from mdin if desired
 ! if (.true.) then
 !    write(str,'(i10)') ctn
 !    write(100881,'(a)') str 
 !    write(100881,*) ' '
 !    do ii=1,ctn
 !       write(str,'(i10)') ii
 !       str=adjustl(str)
 !       write(100881,'(a,a6,f9.5,2f15.5)') "H",str,cirreg1(ii,1:3)
 !    end do
 !    flush(100881)
 ! end if 
 
   ! 1st time to call matvec to get the bf as in the book ch 6.1
   ! To solve (T-EA^-1 B)G=F2,  matvec(0)=(T-EA^-1 B)0-F2=-F2

   call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,& 
           wp,qp,& 
           cirreg1,& 
           unj,unjf,fvec,bf,epsin,epsout,index,index2,& 
           q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,& 
           wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,& 
           wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,& 
           accept,bv,icall,solvopt,norm,inorm)

   do i=1,  ctn
      bf(i)= -fvec(i) 
      unjf(i)=qp(i)   
   end do

   imax= 3*max(l,m,n)! ;stop                             !set values for imax
 
   if (augtol > 0 .and. augtol < 1) then
      tol = augtol
   else
      tol = 1.0d-5
   endif

   if (use_average_tolerance) then
      !write(10,*) 'Using average residual tolerance:tol=avg tol=abs tol/ctn'
      !write(10,*) 'abs tol = ',tol*ctn,'avg tol = ',tol
      tol = tol*ctn
   end if

   !inferface for gmres
   !write(6001,*) mm,l,n,nbnd,nind3,imax,h
   !write(6001,*)   phi
   !write(6001,*)   x,y,z
   !write(6001,*)   index2,index3,iepsav
   !write(6001,*)   wp,qp,wp1,qp1
   !write(6001,*)   cirreg,cirreg1
   !write(6001,*)   unj,u,f,unjf,v,hg,bf
   !write(6001,*)   tol,iter,error,bout,bin
   !write(6001,*)   q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp
   !write(6001,*)   wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe
   !write(6001,*)   wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2
   !write(6001,*)   accept,bv

   ! Check the projection pnts and projection pnts from the exmol are on the same
   ! interface
   !do ii = 1, nbnd
   !   i = iepsav(1,ii); j = iepsav(2,ii); k = iepsav(3,ii)
   !   write(131,*) gox+i*h,goy+j*h,goz+k*h
   ! ! write(134,*) gox+i*h,goy+i*h,goz+i*h,u(i,j,k)
   !end do
   !  stop
  
   unj=0.0d0; unjf=0.0d0;u=0.0d0;f=0.0d0  !Initialization..

   call gmresx(mm,l,m,n,nbnd,ctn,imax,h,phi,x,y,z,xs,ys,zs,&
            index,index2,&
            wp,qp,& 
            cirreg1,&
            unj,u,f,unjf,bf,&
            tol,iter,error,epsout,epsin,&
            q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp ,&
            wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
            wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&                         
            accept,bv,solvopt,pbverbose,norm,inorm)                                  

   ! converting back to the PBSA unit for potential
      
   do k=1, n
      do j=1, m
         do i=1, l
           u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
         end do
      end do
   end do  

   !write(6,*) "-debug: Writting potential map dx file";flush(6)
  !if ( solvopt == 7 ) then
  !!WMBS - writting potential data to pbsa_phi.dx for visualization / testing
  !     call gen_dx_file(xm,ym,zm,h,gox,goy,goz,u(1:xm,1:ym,1:zm),&
  !            potentialdxfilename,potentialdxfilenum,&
  !            potentialdxdataname)
  !end if
   
   !write(6,*) "-debug: Deallocating arrays";flush(6)
   deallocate (x, y, z)
   deallocate (index, index2)
   deallocate (cirreg,cirreg1)
   deallocate (wp,qp)
   deallocate (unj,unjf,bf,fvec)
   deallocate (q0p,qyp,qzp) 
   deallocate (w0p,wyp,wzp)
   deallocate (wyyp,wzzp,wyzp)
   deallocate (wcoe,wxcoe,wycoe)
   deallocate (wzcoe,wxxcoe,wyycoe)
   deallocate (wzzcoe,wxycoe)
   deallocate (wxzcoe,wyzcoe)
   deallocate (sss1,sss2)
   deallocate (c, c2)
   deallocate (f)

end subroutine  aug
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


end subroutine pb_augdrv
