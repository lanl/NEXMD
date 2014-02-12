! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ partition of atoms into internal and external portions according to a sphere
subroutine pb_atmpart( verbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                       inatm,outwat,oution,ipres,outflag,xctr,yctr,zctr,rdiel,sepbuf,x,ifcap )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   logical verbose, pbprint
   integer natom, ibgwat, ienwat, ibgion, ienion, inatm, outwat, oution, ifcap
   integer ipres(*), outflag(natom)
   _REAL_ xctr, yctr, zctr, rdiel, sepbuf
   _REAL_ x(3,natom)

   ! Local variables

   integer ires, num
   _REAL_ sepr, atmr, xtmp, ytmp, ztmp
 
   sepr = (rdiel-sepbuf)**2
   outflag = 0
 
   ! internal portion always include protein atoms
   ! external portion only include water atoms
 
   inatm = ipres(ibgwat)-1
   outwat = 0
   oution = 0

   if( ifcap == 2 .and. ibgion /= 0 ) then
      ! Consider monoatomic ions separately
      inatm = inatm - (ienion - ibgion + 1)
      do ires = ibgion, ienion
         num = ipres(ires)
         xtmp = x(1,num); ytmp = x(2,num); ztmp = x(3,num)
         atmr = (xtmp-xctr)**2 + (ytmp-yctr)**2 + (ztmp-zctr)**2
         if ( atmr > sepr ) then
            oution = oution + 1
            outflag(num) = 1
         end if
         ! Always count ions as "in" even if they are "out"
         inatm = inatm + 1
      end do
   end if

   ! Consider water molecules
   do ires = ibgwat, ienwat
      num = ipres(ires)
      xtmp = x(1,num); ytmp = x(2,num); ztmp = x(3,num)
      atmr = (xtmp-xctr)**2 + (ytmp-yctr)**2 + (ztmp-zctr)**2
      if ( atmr > sepr ) then
         outwat = outwat + 3
         outflag(num) = 1
         outflag(num+1) = 1
         outflag(num+2) = 1
      else
         inatm = inatm + 3
      end if
   end do

   if ( verbose .and. pbprint ) then
      write(6,'(a,2i6,a,f6.3)') ' Atoms are partitioned into two regions', inatm-oution, outwat+oution, ' with a buffer of', sepbuf
   end if

end subroutine pb_atmpart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ partition of atoms into internal and external portions according to a shell
subroutine pb_atmpart2( verbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                        inatm,outwat,oution,ipres,outflag,distance,x )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   logical verbose, pbprint
   integer natom, ibgwat, ienwat, ibgion, ienion, inatm, outwat, oution
   integer ipres(*), outflag(natom)
   _REAL_ distance
   _REAL_ x(3,natom)

   ! Local variables

   integer ires, jatm, num, proatm
   logical isin 
   _REAL_ dist2, dist2tmp, xtmp, ytmp, ztmp
 
   dist2 = distance**2
   outflag = 0
 
   ! Internal portion always include protein atoms
   if ( ibgion /= 0 ) then
      proatm = (ipres(ibgwat) - 1) - (ienion - ibgion + 1)
   else
      proatm = (ipres(ibgwat) - 1)
   end if
   inatm = proatm

   ! External portion include water atoms and ions
   outwat = 0
   oution = 0

   ! Consider monoatomic ions
   if ( ibgion /= 0 ) then
      do ires = ibgion, ienion
         isin = .false.
         num = ipres(ires)
         xtmp = x(1,num); ytmp = x(2,num); ztmp = x(3,num)

         do jatm = 1, proatm
            dist2tmp = (x(1,jatm) - xtmp)**2 + (x(2,jatm) - ytmp)**2 + (x(3,jatm) - ztmp)**2
            if ( dist2tmp < dist2 ) then
               isin = .true.
               exit
            end if
         end do ! jatm
         if(.not. isin) then
            oution = oution + 1
            outflag(num) = 1
         end if

         ! Always count ions as "in" even if they are "out"
         inatm = inatm + 1
      end do ! ires
   end if

   ! Consider water molecules
   do ires = ibgwat, ienwat
      isin = .false.
      num = ipres(ires)
      xtmp = x(1,num); ytmp = x(2,num); ztmp = x(3,num)

      do jatm = 1, proatm
         dist2tmp = (x(1,jatm) - xtmp)**2 + (x(2,jatm) - ytmp)**2 + (x(3,jatm) - ztmp)**2
         if ( dist2tmp < dist2 ) then
            isin = .true.
            exit
         end if
      end do ! jatm
      if(.not. isin) then
         outwat = outwat + 3
         outflag(num) = 1
         outflag(num+1) = 1
         outflag(num+2) = 1
      else
         inatm = inatm + 3
      end if

   end do ! ires

   if ( verbose .and. pbprint ) then
      write(6,'(a,2i6,a,f6.3)') ' Atoms are partitioned into two regions', inatm-oution, outwat+oution
   end if

end subroutine pb_atmpart2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ partition of atoms into internal and external ones according to ligand mask
subroutine pb_atmpart3(verbose,pbprint,natom,buffer,xmin,xmax,ymin,ymax,zmin,zmax,liveflag,realflag,outflag,x)

   use poisson_boltzmann, only : savh, nfocus, solvopt, nbuffer
   use solvent_accessibility, only : radi, dprob
   implicit none

   ! PASSED VARIABLES

   logical verbose, pbprint
   integer natom
   integer liveflag(natom)
   integer realflag(natom)
   integer outflag(natom)
   _REAL_ buffer,xmin,xmax,ymin,ymax,zmin,zmax
   _REAL_ x(3,natom)

#  include "pb_constants.h"

   ! LOCAL VARIABLES

   integer iatm, xtmp, ytmp, ztmp
   _REAL_  lxmin, lymin, lzmin
   _REAL_  lxmax, lymax, lzmax
   _REAL_  htmp, xbox, ybox, zbox
   _REAL_  xlength, ylength, zlength
   _REAL_  pos1(3), pos2(3)
   _REAL_  myradi

   ! find the bounding box large enough surrounding the ligand atoms with a
   ! prescribed buffer

   ! liveflag is initialized in pb_init if ligandmask specified.
   if ( xmin == ZERO .or. xmax == ZERO .or. &
        ymin == ZERO .or. ymax == ZERO .or. &
        zmin == ZERO .or. zmax == ZERO      ) then
      xmin=9999;  ymin=9999;  zmin=9999
      xmax=-9999; ymax=-9999; zmax=-9999

      do iatm = 1, natom
         if ( liveflag(iatm) == 0 ) cycle
         pos1 = x(:,iatm)
         pos2 = pos1 - radi(iatm)
         pos1 = pos1 + radi(iatm)
         if (pos1(1) > xmax) xmax = pos1(1)
         if (pos1(2) > ymax) ymax = pos1(2)
         if (pos1(3) > zmax) zmax = pos1(3)
         if (pos2(1) < xmin) xmin = pos2(1)
         if (pos2(2) < ymin) ymin = pos2(2)
         if (pos2(3) < zmin) zmin = pos2(3)
      end do
   end if

   ! After the bounding box is determined, the list of atoms inside
   ! of it can be stored in liveflag.

   liveflag = 1
   do iatm = 1, natom
      pos1 = x(:,iatm)
      if ( pos1(1) > xmax .or. pos1(1) < xmin .or. &
           pos1(2) > ymax .or. pos1(2) < ymin .or. &
           pos1(3) > zmax .or. pos1(3) < zmin ) then
         liveflag(iatm) = 0
      endif
   enddo

   if ( verbose .and. pbprint ) then
      write(6, '(1x,a,i10)') 'Bounding box identified for live atoms: '
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
      write(6, '(a)') '------- VMD goodie --------'
      write(6, '(a)') '#liveflag boundary Box'
      write(6, '(a)') 'draw materials off'
      write(6, '(a)') 'draw color white'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymin, zmin,'" "', &
        xmax, ymin, zmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymin, zmin,'" "', &
        xmin, ymax, zmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymin, zmin,'" "', &
        xmin, ymin, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmax, ymin, zmin,'" "', &
        xmax, ymax, zmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmax, ymin, zmin,'" "', &
        xmax, ymin, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymax, zmin,'" "', &
        xmax, ymax, zmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymax, zmin,'" "', &
        xmin, ymax, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymin, zmax,'" "', &
        xmax, ymin, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmin, ymin, zmax,'" "', &
        xmin, ymax, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmax, ymax, zmax,'" "', &
        xmax, ymax, zmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmax, ymax, zmax,'" "', &
        xmin, ymax, zmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        xmax, ymax, zmax,'" "', &
        xmax, ymin, zmax,'"'
      write(6, '(a)') '--- End VMD goodie --------'
   end if

   ! align box center and box dimension on grid

   xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO
   htmp = savh(nfocus)
   xbox = nint(xbox/htmp)*htmp; ybox = nint(ybox/htmp)*htmp; zbox = nint(zbox/htmp)*htmp

   ! set real atom box, to be exactly consistent with the FDPB focus box

   xlength = xmax-xmin; ylength = ymax-ymin; zlength = zmax-zmin

   xtmp = nint(xlength/htmp) + nbuffer
   ytmp = nint(ylength/htmp) + nbuffer
   ztmp = nint(zlength/htmp) + nbuffer

   if ( solvopt == 2 ) then
      xtmp = 16*ceiling( dble(xtmp)/16.0d0 ) - 1
      ytmp = 16*ceiling( dble(ytmp)/16.0d0 ) - 1
      ztmp = 16*ceiling( dble(ztmp)/16.0d0 ) - 1
   else
      xtmp = 2*nint( dble(xtmp)*HALF ) + 1
      ytmp = 2*nint( dble(ytmp)*HALF ) + 1
      ztmp = 2*nint( dble(ztmp)*HALF ) + 1
   end if

   xlength = htmp*(xtmp - 1)
   ylength = htmp*(ytmp - 1)
   zlength = htmp*(ztmp - 1)

   lxmin = - (xlength) * HALF + xbox
   lymin = - (ylength) * HALF + ybox
   lzmin = - (zlength) * HALF + zbox

   lxmax = lxmin + xlength 
   lymax = lymin + ylength
   lzmax = lzmin + zlength

   ! flag any atoms inside the focus box

!  liveflag = 0 ! 1 if inside the block  w/o pad
   realflag = 0 ! 1 if inside the block with pad
   outflag  = 1 ! 1 if outside the geometry for surface
   myradi = MAXVAL(radi(1:natom))+dprob*2
   do iatm = 1, natom
      pos1 = x(:,iatm)
      if ( pos1(1) < lxmax .and. pos1(1) > lxmin .and. &
           pos1(2) < lymax .and. pos1(2) > lymin .and. &
           pos1(3) < lzmax .and. pos1(3) > lzmin ) then
         realflag(iatm) = 1
      end if
      if ( pos1(1) <= lxmax+myradi .and. pos1(1) >= lxmin-myradi .and. &
           pos1(2) <= lymax+myradi .and. pos1(2) >= lymin-myradi .and. &
           pos1(3) <= lzmax+myradi .and. pos1(3) >= lzmin-myradi ) then
         outflag(iatm) = 0
      end if
   end do
   if ( verbose .and. pbprint ) then
      write(6, '(1x,a,i10)') 'Corrected bounding box identified for real atoms: '
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', lxmin, lxmax, lxmax-lxmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', lymin, lymax, lymax-lymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', lzmin, lzmax, lzmax-lzmin
      write(6, '(1x,a,i10)') 'Corrected bounding box identified for outside atoms: '
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', lxmin-myradi, lxmax+myradi, lxmax-lxmin+myradi+myradi
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', lymin-myradi, lymax+myradi, lymax-lymin+myradi+myradi
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', lzmin-myradi, lzmax+myradi, lzmax-lzmin+myradi+myradi
      write(6, '(a)') '------- VMD goodie --------'
      write(6, '(a)') '#realflag boundary Box'
      write(6, '(a)') 'draw materials off'
      write(6, '(a)') 'draw color red'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymin, lzmin,'" "', &
        lxmax, lymin, lzmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymin, lzmin,'" "', &
        lxmin, lymax, lzmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymin, lzmin,'" "', &
        lxmin, lymin, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax, lymin, lzmin,'" "', &
        lxmax, lymax, lzmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax, lymin, lzmin,'" "', &
        lxmax, lymin, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymax, lzmin,'" "', &
        lxmax, lymax, lzmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymax, lzmin,'" "', &
        lxmin, lymax, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymin, lzmax,'" "', &
        lxmax, lymin, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin, lymin, lzmax,'" "', &
        lxmin, lymax, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax, lymax, lzmax,'" "', &
        lxmax, lymax, lzmin,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax, lymax, lzmax,'" "', &
        lxmin, lymax, lzmax,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax, lymax, lzmax,'" "', &
        lxmax, lymin, lzmax,'"'

      write(6, '(a)') '#outflag boundary Box'
      write(6, '(a)') 'draw color yellow'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymin-myradi, lzmin-myradi,'" "', &
        lxmax+myradi, lymin-myradi, lzmin-myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymin-myradi, lzmin-myradi,'" "', &
        lxmin-myradi, lymax+myradi, lzmin-myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymin-myradi, lzmin-myradi,'" "', &
        lxmin-myradi, lymin-myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax+myradi, lymin-myradi, lzmin-myradi,'" "', &
        lxmax+myradi, lymax+myradi, lzmin-myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax+myradi, lymin-myradi, lzmin-myradi,'" "', &
        lxmax+myradi, lymin-myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymax+myradi, lzmin-myradi,'" "', &
        lxmax+myradi, lymax+myradi, lzmin-myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymax+myradi, lzmin-myradi,'" "', &
        lxmin-myradi, lymax+myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymin-myradi, lzmax+myradi,'" "', &
        lxmax+myradi, lymin-myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmin-myradi, lymin-myradi, lzmax+myradi,'" "', &
        lxmin-myradi, lymax+myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax+myradi, lymax+myradi, lzmax+myradi,'" "', &
        lxmax+myradi, lymax+myradi, lzmin-myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax+myradi, lymax+myradi, lzmax+myradi,'" "', &
        lxmin-myradi, lymax+myradi, lzmax+myradi,'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "', &
        lxmax+myradi, lymax+myradi, lzmax+myradi,'" "', &
        lxmax+myradi, lymin-myradi, lzmax+myradi,'"'
      write(6, '(a)') '--- End VMD goodie --------'
   end if

end subroutine pb_atmpart3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Grid Partition subroutine by Mengjuei Hsieh
!+ For slave node to claim ownership from the focus sub-grids.
!+ This has been generalized that it doesn't need to know about the MPI
!+ If we want to add print out for multiple thread (MPI), we need to pass
!+ the logical variable "master" for the print out control.
subroutine pb_atmpart_mb(verbose, pbprint, natom, istart, iend, h,   &
              blknx,  blkny,  blknz,  blkxlo, blkylo, blkzlo,        &
              blkxup, blkyup, blkzup, blkgox, blkgoy, blkgoz,        &
              x, lfocuswpad, lfocus0pad, lfocus4s3f,                 &
                 nfirstwpad, nfirst0pad, nfirst4s3f, master)
   use solvent_accessibility, only : radi, dprob
   implicit none

   ! Passed variables

   logical verbose, pbprint, master
   integer natom
   integer istart, iend
   integer blkxlo(*), blkylo(*), blkzlo(*)
   integer blkxup(*), blkyup(*), blkzup(*)
   integer blknx(*),  blkny(*),  blknz(*)
   integer lfocuswpad(*), lfocus0pad(*), lfocus4s3f(*)
   integer nfirstwpad(*), nfirst0pad(*), nfirst4s3f(*)
   _REAL_ h
   _REAL_ blkgox(*), blkgoy(*), blkgoz(*)
   _REAL_ x(3,natom)

   ! Local variables

   integer i, j, k, l, m, n
   !logical blockhasatom
   logical myinside(natom)
   _REAL_ dimxlower, dimylower, dimzlower
   _REAL_ dimxupper, dimyupper, dimzupper
   _REAL_ blkxlower, blkylower, blkzlower
   _REAL_ blkxupper, blkyupper, blkzupper
   _REAL_ myradi

   !check starts.
   k = 1
   l = 1
   m = 1
   n = 1
   myinside = .false.
   myradi = MAXVAL(radi(1:natom))+dprob*2
   blkloop: do i = istart, iend
      dimxlower = blkgox(i) + h
      blkxlower = blkgox(i) + blkxlo(i)*h
      dimxupper = blkgox(i) +  blknx(i)*h
      blkxupper = blkgox(i) + blkxup(i)*h
      dimylower = blkgoy(i) + h
      blkylower = blkgoy(i) + blkylo(i)*h
      dimyupper = blkgoy(i) +  blkny(i)*h
      blkyupper = blkgoy(i) + blkyup(i)*h
      dimzlower = blkgoz(i) + h
      blkzlower = blkgoz(i) + blkzlo(i)*h
      dimzupper = blkgoz(i) +  blknz(i)*h
      blkzupper = blkgoz(i) + blkzup(i)*h

      atmloop: do j = 1, natom
         !since the outer edge does not contain atoms, it should be safe.
         if ( x(1,j) <  dimxlower - myradi - h ) cycle atmloop
         if ( x(1,j) >= dimxupper + myradi + h ) cycle atmloop
         if ( x(2,j) <  dimylower - myradi - h ) cycle atmloop
         if ( x(2,j) >= dimyupper + myradi + h ) cycle atmloop
         if ( x(3,j) <  dimzlower - myradi - h ) cycle atmloop
         if ( x(3,j) >= dimzupper + myradi + h ) cycle atmloop
!print *, dimxlower - myradi - h, dimxupper + myradi + h, dimylower - myradi - h
!print *, dimyupper + myradi + h, dimzlower - myradi - h, dimzupper + myradi + h;stop
         lfocus4s3f(n) = j
         n = n + 1
         if ( x(1,j) <  dimxlower ) cycle atmloop
         if ( x(1,j) >= dimxupper ) cycle atmloop
         if ( x(2,j) <  dimylower ) cycle atmloop
         if ( x(2,j) >= dimyupper ) cycle atmloop
         if ( x(3,j) <  dimzlower ) cycle atmloop
         if ( x(3,j) >= dimzupper ) cycle atmloop
         lfocuswpad(k) = j
         k = k + 1
         if ( x(1,j) <  blkxlower ) cycle atmloop
         if ( x(1,j) >= blkxupper ) cycle atmloop
         if ( x(2,j) <  blkylower ) cycle atmloop
         if ( x(2,j) >= blkyupper ) cycle atmloop
         if ( x(3,j) <  blkzlower ) cycle atmloop
         if ( x(3,j) >= blkzupper ) cycle atmloop
         if (myinside(j)) then
            write(6,*) "atom:",j," is assigned at least twice in the block ",i
            cycle atmloop
         endif
         lfocus0pad(l) = j
         l = l + 1
         myinside(j) = .true.
      end do atmloop
      nfirst4s3f(m+1) = nfirst4s3f(1)+n-1
      nfirstwpad(m+1) = nfirstwpad(1)+k-1
      nfirst0pad(m+1) = nfirst0pad(1)+l-1

      m = m + 1
   end do blkloop
   !print *,"mjhsieh: k,l,m",k,l,m
   do j = 1, natom
      if (.not.myinside(j)) then
         write(6,*) "PBerror: At least one atom is outside ", &
           "the focus fine grids."
         call mexit(6,1)
      endif
   end do

   ! TODO: add a permanent check to make sure the atom numbers sum up to natom
   ! This has to be after the MPI_REDUCE
   ! write(6,*) m-1,nfirst0pad(m)-1,natom
end subroutine pb_atmpart_mb
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Atom-based nblist for PBMD/PBDOCK
!+ This is also used by multi-block focusing
subroutine pb_atmlist( verbose,pbprint,maxnba,natom,ntypes,iac,ico,natex, &
                        nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,&
                        cn2,cn1pb,cn2pb,cn3pb,cg,acrd)
   use poisson_boltzmann, only : liveflag, outflag
   implicit none
    
#  include "pb_constants.h"
    
   ! Passed variables
    
   integer maxnba,natom,ntypes               ! should be readonly
   integer natex(*),nshrt(0:natom),nex(*),iex(64,natom),iac(*),ico(*) !readonly
   integer iar1pb(4,0:natom)                 ! will be updated
   integer iprshrt(*)                        ! will be updated
   _REAL_ acrd(3,*), cn1(*), cn2(*), cg(*)   ! should be readonly
   _REAL_ cn1pb(*), cn2pb(*), cn3pb(*)       ! will be updated
   _REAL_ cutnb, cutsa, cutfd                ! should be readonly
    
   ! Local variables
    
   logical verbose, pbprint, i_is_excluded, j_is_excluded
   integer iaci, ic
   integer iatm, jatm
   integer i, j, jp, k, kp
   integer cntr, eclose, pclose, sclose, nclose
   integer tmpex(MAXNEI), tmppb(MAXNEI), tmpsa(MAXNEI), tmpnb(MAXNEI)
   _REAL_ xi, yi, zi, dx, dy, dz, d2
   _REAL_ cgi
    
   ! set up zeroth atom for limits
    
   iar1pb(1, 0) = 0
   iar1pb(2, 0) = 0
   iar1pb(3, 0) = 0
   iar1pb(4, 0) = 0
   iar1pb(1, natom) = 0
   iar1pb(2, natom) = 0
   iar1pb(3, natom) = 0
   iar1pb(4, natom) = 0
   nex(1:natom) = 0
   iprshrt(1:maxnba) = 0
   cn1pb(1:maxnba) = 0
   cn2pb(1:maxnba) = 0
   cn3pb(1:maxnba) = 0
    
   ! this is the global index of atom-based pair
    
   cntr = 0
    
   ! we shall assume that only ligand atoms are moving for now,
   ! and they are continuously located in the coordinate array
   ! mjhsieh: the original code was looping ligand atoms over
   !          all atoms in the docking box, is that correct?
    
   do i = 1, natom
      if ( outflag(i) == 1 ) then
         iar1pb(1, i) = cntr
         iar1pb(2, i) = cntr
         iar1pb(3, i) = cntr
         iar1pb(4, i) = cntr
         cycle
      end if
      iaci = ntypes*(iac(i)-1)
      cgi  = -cg(i)

      ! the inner loop is over all atoms in the docking box 
       
      ! part a: save nonboneded pairs into tmppb, tmpsa, and tmpnb
       
      xi = acrd(1, i)
      yi = acrd(2, i)
      zi = acrd(3, i)
      eclose = 0
      pclose = 0
      sclose = 0
      nclose = 0
      do j = 1, natom

         if ( outflag(j) == 0 .and. j <= i )  cycle

         ! save excluded pairs into tmpex

         j_is_excluded = .false.
         do kp = nshrt(i-1) + 1, nshrt(i)
            k = natex(kp)
            if (k == 0) cycle
            if (k == j) then
               eclose = eclose + 1
               tmpex(eclose) = j
               j_is_excluded = .true.
               exit
            end if
         end do ! kp = nshrt(i-1) + 1, nshrt(i)
         if (j_is_excluded) cycle

         i_is_excluded = .false.
         do kp = nshrt(j-1) + 1, nshrt(j)
            k = natex(kp)
            if (k == 0) cycle
            if (k == i) then
               eclose = eclose + 1
               tmpex(eclose) = j
               i_is_excluded = .true.
               exit
            end if
         end do ! kp = nshrt(j-1) + 1, nshrt(j)
         if (i_is_excluded) cycle

         ! save nbonded pairs into the following ...

         dx = xi - acrd(1, j)
         dy = yi - acrd(2, j)
         dz = zi - acrd(3, j)
         d2 = dx**2 + dy**2 + dz**2
         if (d2 <= cutfd) then
            pclose = pclose + 1
            tmppb(pclose) = j
         else if (d2 <= cutsa) then
            sclose = sclose + 1
            tmpsa(sclose) = j
         else if (d2 <= cutnb) then
            nclose = nclose + 1
            tmpnb(nclose) = j
         end if
      end do ! j = 1, natom
!write(2012,*) '====', i,eclose,cntr
!write(2012,*) '====', i,pclose,cntr, cutfd
!write(2012,*) '====', i,sclose,cntr, cutsa
!write(2012,*) '====', i,nclose,cntr, cutnb
       
      ! part b: pack them into the new atom-based nblist
      ! since there is no separation of h-atom and other atoms, we need to
      ! set up CN1, CN2 properly in every situation, though we may not use it
      ! in docking. Once confirmed in the later phase of the project,
      ! these extra arrays can be removed. 
       
      if (eclose > MAXNEI .or. pclose > MAXNEI .or. sclose > MAXNEI .or. &
          nclose > MAXNEI) then
         write(6, *) 'PB bomb in pb_atmlist(): MAXNEI too short'
         call mexit(6, 1)
      end if
      if (eclose + pclose + sclose + nclose + cntr > maxnba) then
         write(6, '(4x,a,i8)') 'PB Bomb in pb_atmlist(): maxnba too short'
         call mexit(6, 1)
      end if

      do j = 1, eclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpex(j)
         iatm = i
         jatm = tmpex(j)
         nex(iatm) = nex(iatm) + 1
         iex(nex(iatm),iatm) = jatm
         nex(jatm) = nex(jatm) + 1
         iex(nex(jatm),jatm) = iatm
      end do ! j = 1, eclose
      iar1pb(1, i) = cntr ! for exclusion and FDCOUL calculation
      do j = 1, pclose
         cntr = cntr + 1
         iprshrt(cntr) = tmppb(j)
         ic = ico(iaci+iac(tmppb(j)))
         if ( ic > 0 ) then ! normal nonbonded pairs
            cn1pb(cntr) = cn1(ic)
            cn2pb(cntr) = cn2(ic)
         else               ! h-bonded pairs
            cn1pb(cntr) = ZERO
            cn2pb(cntr) = ZERO
         end if
         cn3pb(cntr) = cgi*cg(tmppb(j))
      end do ! j = 1, pclose
      iar1pb(2, i) = cntr ! for FDCOUL calculation
      do j = 1, sclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpsa(j)
         ic = ico(iaci+iac(tmpsa(j)))
         if ( ic > 0 ) then ! normal nonbonded pairs
            cn1pb(cntr) = cn1(ic)
            cn2pb(cntr) = cn2(ic)
         else               ! h-bonded pairs
            cn1pb(cntr) = ZERO
            cn2pb(cntr) = ZERO
         end if
         cn3pb(cntr) = cgi*cg(tmpsa(j))
      end do ! j = 1, sclose
      iar1pb(3, i) = cntr ! for SA calculation
      do j = 1, nclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpnb(j)
         ic = ico(iaci+iac(tmpnb(j)))
         if ( ic > 0 ) then ! normal nonbonded pairs
            cn1pb(cntr) = cn1(ic)
            cn2pb(cntr) = cn2(ic)
         else               ! h-bonded pairs
            cn1pb(cntr) = ZERO
            cn2pb(cntr) = ZERO
         end if
         cn3pb(cntr) = cgi*cg(tmpnb(j))
      end do ! j = 1, nclose
      iar1pb(4, i) = cntr ! for VDW calculation
   end do ! i = 1, natom 
!do i = 1, natom
!   write(2013,*) i,outflag(i),iar1pb(3,i)-iar1pb(4,i-1)
!end do
   if ( verbose .and. pbprint ) write(6,'(2x,a,i9)') 'NB-update: atom-based nb list', cntr
   
end subroutine pb_atmlist 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up FDPB grid for a de novo call
subroutine pb_setgrd( ipb, verbose, pbprint, initial, pbgrid, ifcap, atmfirst, atmlast, xcap, ycap, zcap, cutcap )
    
   use poisson_boltzmann
   use solvent_accessibility, only : radi, dprob

   implicit none

   ! Passed variables
    
   logical verbose, pbprint, initial, pbgrid
   integer ipb,ifcap, atmfirst, atmlast
   _REAL_ xcap, ycap, zcap, cutcap
    
   ! Local variables
    
   integer totsavxmymzm
   integer xm_max, ym_max, zm_max
   integer xmymzm_max,xmymzm_ext
   integer l, alloc_err(32)
!  integer block

   ! get center and dimension information of the current molecule
    
   !call setgrd( verbose, pbprint, initial, ifcap, atmfirst, atmlast, xcap, ycap, zcap, cutcap )
   call setgrd( verbose, pbprint, pbgrid, ifcap, atmfirst, atmlast, xcap, ycap, zcap, cutcap )
   
   ! Mengjuei> xmymzm_ext = (xm + 2)*(ym + 2)*(zm + 2)
   xmymzm_max = maxval(savxmymzm(1:nfocus))
   xm_max     = maxval(savxm(1:nfocus))
   ym_max     = maxval(savym(1:nfocus))
   zm_max     = maxval(savzm(1:nfocus))
   if ( multiblock ) then
      xmymzm_max = max(xmymzm_max,maxval(blknxnynz))
      xm_max = max(xm_max,maxval(blknx))
      ym_max = max(ym_max,maxval(blkny))
      zm_max = max(zm_max,maxval(blknz))
   end if
   xmymzm_ext = (xm_max + 2) * (ym_max + 2) * (zm_max + 2)
    
   ! if allocate working arrays for fdpb
    
   alloc_err(1:32) = 0
   if ( .not. initial ) then
      deallocate(    phi, stat = alloc_err(1 ) )
      deallocate(  chgrd, stat = alloc_err(2 ) )
      deallocate(   epsx, stat = alloc_err(4 ) )
      deallocate(   epsy, stat = alloc_err(5 ) )
      deallocate(   epsz, stat = alloc_err(6 ) )
      if(ipb /= 4 .and. ipb /= 5) deallocate(saltgrd, stat = alloc_err(7 ) )
      deallocate(     bv, stat = alloc_err(3 ) )

      deallocate(  insas, stat = alloc_err(8 ) )
      deallocate( atmsas, stat = alloc_err(9 ) )
      deallocate( lvlset, stat = alloc_err(10) )
      deallocate(     zv, stat = alloc_err(11) )

      deallocate(   cphi, stat = alloc_err(13) )
      if(ipb /= 4 .and. ipb /= 5) deallocate( fedgex, stat = alloc_err(14) )
      if(ipb /= 4 .and. ipb /= 5) deallocate( fedgey, stat = alloc_err(15) )
      if(ipb /= 4 .and. ipb /= 5) deallocate( fedgez, stat = alloc_err(16) )
      deallocate( iepsav, stat = alloc_err(17) )
      if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavx, stat = alloc_err(18) )
      if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavy, stat = alloc_err(19) )
      if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavz, stat = alloc_err(20) )

      deallocate(     xs, stat = alloc_err(30) )
      if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
           alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
           alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
           alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
           alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
           alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30) /= 0 ) then
        write(6, *) 'PB bomb in pb_setgrd(): Deallocation aborted', alloc_err(1:30)
         call mexit(6, 1)
      end if
   end if

   ! physical constant maps for numerical solutions

   allocate(    phi( xmymzm_max ), stat = alloc_err(1 ) )
   allocate(  chgrd( xmymzm_max ), stat = alloc_err(2 ) )
   allocate(   epsx( xmymzm_max+ym_max*zm_max ), stat = alloc_err(3 ) )
   allocate(   epsy( xmymzm_max+xm_max*zm_max ), stat = alloc_err(4 ) )
   allocate(   epsz( xmymzm_max+xm_max*ym_max ), stat = alloc_err(5 ) )
   if(ipb /= 4 .and. ipb /= 5) allocate(saltgrd( xmymzm_max ), stat = alloc_err(6 ) )
   allocate(     bv( xmymzm_max ), stat = alloc_err(7 ) )

   ! geometry propery maps and auxiliary arrays for mapping dielectric and stern interfaces

   allocate(  insas( xmymzm_ext ), stat = alloc_err(8 ) )
   allocate( atmsas( xmymzm_ext ), stat = alloc_err(9 ) )
   allocate( lvlset( xmymzm_ext ), stat = alloc_err(10) )
   allocate(     zv( xmymzm_ext ), stat = alloc_err(11) )

   ! physical property maps for forces

   allocate(   cphi   (1:xmymzm_max), stat = alloc_err(13) )
   if(ipb /= 4 .and. ipb /= 5) allocate( fedgex   (1:xmymzm_max), stat = alloc_err(14) )
   if(ipb /= 4 .and. ipb /= 5) allocate( fedgey   (1:xmymzm_max), stat = alloc_err(15) )
   if(ipb /= 4 .and. ipb /= 5) allocate( fedgez   (1:xmymzm_max), stat = alloc_err(16) )
   allocate( iepsav (4,1:xmymzm_max), stat = alloc_err(17) )
   if(ipb /= 4 .and. ipb /= 5) allocate( iepsavx(4,1:xmymzm_max), stat = alloc_err(18) )
   if(ipb /= 4 .and. ipb /= 5) allocate( iepsavy(4,1:xmymzm_max), stat = alloc_err(19) )
   if(ipb /= 4 .and. ipb /= 5) allocate( iepsavz(4,1:xmymzm_max), stat = alloc_err(20) )

   ! the saved phi array for pbmd/pbdock

   totsavxmymzm = 0
   do l = 1, nfocus
      totsavxmymzm = totsavxmymzm + savxmymzm(l)+2*SAVXMYM(l)
   end do
   if ( multiblock ) then
      totsavxmymzm = max(totsavxmymzm, (maxval(blknxnynz) + &
         savxmymzm(1) + 2*savxmym(1) + 2*maxval(blknxny)))
   end if
   allocate( xs(1:totsavxmymzm), stat = alloc_err(30) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30) /= 0 ) then
      write(6, *) 'PB bomb in pb_setgrd(): Allocation aborted', alloc_err(1:30)
      call mexit(6, 1)
   end if

   ! initialize saved phi map
    
   xs = ZERO
    
   ! save fine grid limits for checking of atom-out-of-grid situation
   
   gxmin = savgox(nfocus) + 3*savh(nfocus);
   gxmax = savgox(nfocus) + (savxm(nfocus)-2)*savh(nfocus)
   gymin = savgoy(nfocus) + 3*savh(nfocus);
   gymax = savgoy(nfocus) + (savym(nfocus)-2)*savh(nfocus)
   gzmin = savgoz(nfocus) + 3*savh(nfocus);
   gzmax = savgoz(nfocus) + (savzm(nfocus)-2)*savh(nfocus)

   ! Auxiliary setup for multiblock

!  if ( multiblock ) then
!     block = 0
!     do l = 2, nfocus
!        block = block + levelblock(l)
!     end do
!  end if

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set up FDPB grid for a de novo call
subroutine setgrd( verbose,pbprint,initial,ifcap,atmfirst,atmlast,xcap,ycap,zcap,cutcap )
   implicit none
    
#  include "pb_constants.h"
#  include "extra.h"
    
   ! Passed variables
    
   logical verbose, pbprint, initial
   integer ifcap, atmfirst, atmlast 
   _REAL_ xcap, ycap, zcap, cutcap
    
   ! Local variables
    
   integer l, iatm
   logical isfillratiosane, newbox
   _REAL_ xlength, ylength, zlength
   _REAL_ xbox, ybox, zbox, htmp

   ! Local multi-block variables
   integer tmpxm,  tmpym,  tmpzm
   integer tmpgox, tmpgoy, tmpgoz
   integer block,  blkpad
   integer iblk,   jblk,   kblk
   integer alloc_err(19)
    
   ! set bounding box center for selected atoms
    
   if ( verbose .and. pbprint ) then
      write(6, *)
      write(6, *)
      write(6, *) '======== Setting up Grid Parameters ========'
      write(6, *) 'Using bounding box for grid setup'
   end if

   ! :::: set the second/fine grid ::::

   l = nfocus

   ! find the bounding box if it has not been initialized from
   ! pb_read() or from pb_atmpart3() for focusing on the ligand.

   newbox = .false.
   if ( xmin == ZERO .and. xmax == ZERO .and. &
        ymin == ZERO .and. ymax == ZERO .and. &
        zmin == ZERO .and. zmax == ZERO ) newbox = .true.

   if ( .not. ligand ) then
      if ( ifcap == 0 .or. ifcap == 5 ) then
         xmin =  9999.0; ymin =  9999.0; zmin =  9999.0
         xmax = -9999.0; ymax = -9999.0; zmax = -9999.0
         do iatm = atmfirst, atmlast
            if ( ifcap == 5 .and. outflag(iatm) == 1 ) cycle
            if ( acrd(1,iatm)-radi(iatm) .lt. xmin ) xmin = acrd(1,iatm)-radi(iatm)
            if ( acrd(1,iatm)+radi(iatm) .gt. xmax ) xmax = acrd(1,iatm)+radi(iatm)
            if ( acrd(2,iatm)-radi(iatm) .lt. ymin ) ymin = acrd(2,iatm)-radi(iatm)
            if ( acrd(2,iatm)+radi(iatm) .gt. ymax ) ymax = acrd(2,iatm)+radi(iatm)
            if ( acrd(3,iatm)-radi(iatm) .lt. zmin ) zmin = acrd(3,iatm)-radi(iatm)
            if ( acrd(3,iatm)+radi(iatm) .gt. zmax ) zmax = acrd(3,iatm)+radi(iatm)
         end do
!        xmin = xmin - dprob; xmax = xmax + dprob
!        ymin = ymin - dprob; ymax = ymax + dprob
!        zmin = zmin - dprob; zmax = zmax + dprob
      else
         xmin = xcap - cutcap; xmax = xcap + cutcap
         ymin = ycap - cutcap; ymax = ycap + cutcap
         zmin = zcap - cutcap; zmax = zcap + cutcap
      end if
   end if

   ! this is the box center
   ! round it to the nearest h unit for easy restarting

   xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO
   htmp = savh(nfocus)
   xbox = nint(xbox/htmp)*htmp; ybox = nint(ybox/htmp)*htmp; zbox = nint(zbox/htmp)*htmp

   if ( verbose .and. pbprint .and. nfocus > 1 ) then
      write(6, '(1x,a,i10)') 'Bounding Box at level:  ', l
      write(6, '(1x,a,3f10.3)') 'Bounding Box Center:  ', xbox, ybox, zbox
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
   end if

   ! this is for updating the bounding box information 

   if ( initial ) then
      cxbox(l) = xbox; cybox(l) = ybox; czbox(l) = zbox
      if ( verbose .and. pbprint .and. nfocus > 1 ) write(6, '(a,1x,i5,1x,3f10.3)') &
         '   beginning box center at level ', l, cxbox(l), cybox(l), czbox(l)
   else
      savxbox(l) = cxbox(l); savybox(l) = cybox(l); savzbox(l) = czbox(l)
      if ( verbose .and. pbprint .and. nfocus > 1 ) write(6, '(a,1x,i5,1x,3f10.3)') &
         '   previous box center at level ', l, savxbox(l), savybox(l), savzbox(l)
   end if
    
   ! set the grid dimension
    
   xlength = xmax-xmin; ylength = ymax-ymin; zlength = zmax-zmin
   if ( outphi .and. phiform == 0 ) then
      xlength = max(xlength, ylength, zlength)
      ylength = xlength; zlength = xlength
   endif
   savxm(l) = nint( xlength/savh(l) ) + nbuffer
   savym(l) = nint( ylength/savh(l) ) + nbuffer
   savzm(l) = nint( zlength/savh(l) ) + nbuffer
    
   !    adjust for the multigrid solver (four levels) or other solver options ...
    
   if ( solvopt == 2 ) then 
      savxm(l) = 16*ceiling( dble(savxm(l))/16.0d0 ) - 1
      savym(l) = 16*ceiling( dble(savym(l))/16.0d0 ) - 1
      savzm(l) = 16*ceiling( dble(savzm(l))/16.0d0 ) - 1   
   else
      savxm(l) = 2*nint( dble(savxm(l))*HALF ) + 1
      savym(l) = 2*nint( dble(savym(l))*HALF ) + 1
      savzm(l) = 2*nint( dble(savzm(l))*HALF ) + 1
   end if
    
   !    additional stored data and printing
    
   savxmym(l) = savxm(l)*savym(l)
   savxmymzm(l) = savxmym(l)*savzm(l)
   if ( verbose .and. pbprint .and. nfocus > 1 ) write(6, '(a,i5,1x,3i5)') &
      ' Grid dimension at level ', l, savxm(l), savym(l), savzm(l)
    
   ! set grid origin
    
   if ( initial ) then
      savgox(l) = - dble(savxm(l)+1) * savh(l) * HALF + xbox
      savgoy(l) = - dble(savym(l)+1) * savh(l) * HALF + ybox
      savgoz(l) = - dble(savzm(l)+1) * savh(l) * HALF + zbox
   else
      cxbox(l) = savxbox(l) + nint( (xbox - savxbox(l))/savh(l) )*savh(l)
      cybox(l) = savybox(l) + nint( (ybox - savybox(l))/savh(l) )*savh(l)
      czbox(l) = savzbox(l) + nint( (zbox - savzbox(l))/savh(l) )*savh(l)
      if ( verbose .and. pbprint .and. nfocus > 1 ) write(6, '(a,i5,1x,3f10.3)') &
         ' Box center corrected at level ', l, cxbox(l), cybox(l), czbox(l)
      savgox(l) = - dble(savxm(l)+1) * savh(l) * HALF + cxbox(l)
      savgoy(l) = - dble(savym(l)+1) * savh(l) * HALF + cybox(l)
      savgoz(l) = - dble(savzm(l)+1) * savh(l) * HALF + czbox(l)
   end if
   if ( verbose .and. pbprint .and. nfocus > 1 ) write(6, '(a,i5,1x,3f10.3)') &
      ' Grid origin corrected at level ', l, savgox(l), savgoy(l), savgoz(l)

   ! Multi-block section, mjhsieh
   if ( .not. multiblock ) then
      ngrdblkx=savxm(l)
      ngrdblky=savym(l)
      ngrdblkz=savzm(l)
   else if ( l /= 2 ) then
      write(6,*) 'setgrd(): nfocus should be 2 with multiblock focusing'
      call mexit(6,1)
   else
      levelblock(1) = 1
      h = savh(l)
      blkpad = nint(buffer/h)
      tmpxm=savxm(l)
      tmpym=savym(l)
      tmpzm=savzm(l)
      if (xmblk == 0) xmblk = ceiling(REAL(tmpxm-1)/(ngrdblkx-1))
      if (ymblk == 0) ymblk = ceiling(REAL(tmpym-1)/(ngrdblky-1))
      if (zmblk == 0) zmblk = ceiling(REAL(tmpzm-1)/(ngrdblkz-1))
      block = xmblk*ymblk*zmblk
      allocate(     blkxo(block), stat = alloc_err( 1) )
      allocate(     blkyo(block), stat = alloc_err( 2) )
      allocate(     blkzo(block), stat = alloc_err( 3) )
      allocate(    blkxlo(block), stat = alloc_err( 4) )
      allocate(    blkylo(block), stat = alloc_err( 5) )
      allocate(    blkzlo(block), stat = alloc_err( 6) )
      allocate(    blkxup(block), stat = alloc_err( 7) )
      allocate(    blkyup(block), stat = alloc_err( 8) )
      allocate(    blkzup(block), stat = alloc_err( 9) )
      allocate(     blknx(block), stat = alloc_err(10) )
      allocate(     blkny(block), stat = alloc_err(11) )
      allocate(     blknz(block), stat = alloc_err(12) )
      allocate(   blknxny(block), stat = alloc_err(13) )
      allocate(   blknynz(block), stat = alloc_err(14) )
      allocate(   blknxnz(block), stat = alloc_err(15) )
      allocate( blknxnynz(block), stat = alloc_err(16) )
      allocate(    blkgox(block), stat = alloc_err(17) )
      allocate(    blkgoy(block), stat = alloc_err(18) )
      allocate(    blkgoz(block), stat = alloc_err(19) )
      if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
           alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
           alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
           alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)               &
           /= 0 ) then
         write(6, *) 'PB bomb in setgrd(): Allocation aborted', alloc_err(1:19)
         call mexit(6, 1)
      end if
      if ( mod(ngrdblkx-1,fscale) /= 0 .or. &
           mod(ngrdblkz-1,fscale) /= 0 .or. &
           mod(ngrdblkz-1,fscale) /= 0      ) then
         write(6,*) 'PB bomb in setgrd(): unfavorable ngrdbl[x-z] setting.'
         call mexit(6, 1)
      end if
      tmpxm = xmblk*(ngrdblkx-1)+1
      tmpym = ymblk*(ngrdblky-1)+1
      tmpzm = zmblk*(ngrdblkz-1)+1
      tmpgox = xbox-(nint(REAL(tmpxm-1)*HALF)-1)*h
      tmpgoy = ybox-(nint(REAL(tmpym-1)*HALF)-1)*h
      tmpgoz = zbox-(nint(REAL(tmpzm-1)*HALF)-1)*h
      savxm(l)  = tmpxm
      savym(l)  = tmpym
      savzm(l)  = tmpzm
      savgox(l) = tmpgox
      savgoy(l) = tmpgoy
      savgoz(l) = tmpgoz
      levelblock(l) = xmblk * ymblk * zmblk
      if ( master ) then
         write(6,'(a)') "multiblock range:"
         write(6,'(3f10.3)') tmpgox +h, tmpgox +tmpxm*h
         write(6,'(3f10.3)') tmpgoy +h, tmpgoy +tmpym*h
         write(6,'(3f10.3)') tmpgoz +h, tmpgoz +tmpzm*h
      end if
!     gxmin = tmpgox+h-tmpxm*h
!     gymin = tmpgoy+h-tmpym*h
!     gzmin = tmpgoz+h-tmpzm*h
!     gxmax = tmpgox  +tmpxm*h
!     gymax = tmpgoy  +tmpym*h
!     gzmax = tmpgoz  +tmpzm*h

      block = 0
      tmpxm = ngrdblkx + 2*blkpad
      tmpym = ngrdblky + 2*blkpad
      tmpzm = ngrdblkz + 2*blkpad
      do kblk = 1, zmblk; do jblk = 1, ymblk; do iblk = 1, xmblk
         block = block + 1
         ! set up the default dimension, origin, and grid limits within
         ! which the atoms' forces will be computed. That is to say
         ! atoms within the padding regions are discarded because they
         ! are too close to the block boundary, i.e. accuracy too low.
         gox = tmpgox + (iblk - 1)*(ngrdblkx - 1)*h - buffer
         goy = tmpgoy + (jblk - 1)*(ngrdblky - 1)*h - buffer
         goz = tmpgoz + (kblk - 1)*(ngrdblkz - 1)*h - buffer
         blkxo(block) = (iblk - 1)*(ngrdblkx - 1)   - blkpad
         blkyo(block) = (jblk - 1)*(ngrdblky - 1)   - blkpad
         blkzo(block) = (kblk - 1)*(ngrdblkz - 1)   - blkpad
         ilower = 1 + blkpad
         jlower = 1 + blkpad
         klower = 1 + blkpad
         
         iupper = ilower + (ngrdblkx - 1)
         jupper = jlower + (ngrdblky - 1)
         kupper = klower + (ngrdblkz - 1)
         blknx(block) = tmpxm
         blkny(block) = tmpym
         blknz(block) = tmpzm
         blknxny(block) = tmpxm*tmpym
         blknynz(block) = tmpym*tmpzm
         blknxnz(block) = tmpxm*tmpzm
         blknxnynz(block) = tmpxm*tmpym*tmpzm
         blkgox(block) = gox
         blkgoy(block) = goy
         blkgoz(block) = goz
         blkxlo(block) = ilower; blkxup(block) = iupper
         blkylo(block) = jlower; blkyup(block) = jupper
         blkzlo(block) = klower; blkzup(block) = kupper
         !if ( verbose .and. pbprint .and. master) then
         if ( master .and. nfocus > 1 ) then
            write(6, '(a,2i5,1x,3i5)') ' Grid dimension at level/block ',&
               l, block,  tmpxm, tmpym, tmpzm
            write(6, '(a,2i5,1x,3f10.3)') ' Grid orig at level/block', &
               l, block, gox, goy, goz
            if ( ilower == 1 ) then
               write(6,'(a,2f10.3)') ' Block inner limits ', gox+h*ilower, gox+h*iupper
            else
               write(6,'(a,2f10.3)') ' Block inner limits ', gox+h*ilower+h, gox+h*iupper
            end if
            if ( klower == 1 ) then
               write(6,'(a,2f10.3)') ' Block inner limits ', goy+h*jlower, goy+h*jupper
            else
               write(6,'(a,2f10.3)') ' Block inner limits ', goy+h*jlower+h, goy+h*jupper
            end if
            if ( klower == 1 ) then
               write(6,'(a,2f10.3)') ' Block inner limits ', goz+h*klower, goz+h*kupper
            else
               write(6,'(a,2f10.3)') ' Block inner limits ', goz+h*klower+h, goz+h*kupper
            end if
         end if
      end do; end do; end do
   end if
   ! end of multi-block

   ! if no focusing is used, we are done ...
 
   !if ( nfocus > 1 ) then

   ! :::::: set the first/coarse grid ::::::
    
   ! find the bounding box of all atoms
    
   if ( ifcap == 0 .or. ifcap == 5 ) then
      xmin =  9999.0; ymin =  9999.0; zmin =  9999.0
      xmax = -9999.0; ymax = -9999.0; zmax = -9999.0
      do iatm = atmfirst, atmlast
         if ( ifcap == 5 .and. outflag(iatm) == 1 ) cycle
         if ( acrd(1,iatm)-radi(iatm) .lt. xmin ) xmin = acrd(1,iatm)-radi(iatm)
         if ( acrd(1,iatm)+radi(iatm) .gt. xmax ) xmax = acrd(1,iatm)+radi(iatm)
         if ( acrd(2,iatm)-radi(iatm) .lt. ymin ) ymin = acrd(2,iatm)-radi(iatm)
         if ( acrd(2,iatm)+radi(iatm) .gt. ymax ) ymax = acrd(2,iatm)+radi(iatm)
         if ( acrd(3,iatm)-radi(iatm) .lt. zmin ) zmin = acrd(3,iatm)-radi(iatm)
         if ( acrd(3,iatm)+radi(iatm) .gt. zmax ) zmax = acrd(3,iatm)+radi(iatm)
      end do
      xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO

      ! rounding to nearest h unit for easy restarting

      htmp = savh(nfocus)
      xbox = nint(xbox/htmp)*htmp; ybox = nint(ybox/htmp)*htmp; zbox = nint(zbox/htmp)*htmp
   else
      xmin = xcap - cutcap; xmax = xcap + cutcap
      ymin = ycap - cutcap; ymax = ycap + cutcap
      zmin = zcap - cutcap; zmax = zcap + cutcap
      xbox = xcap; ybox = ycap; zbox = zcap
   end if
    
   if ( verbose .and. pbprint ) then
      write(6, '(1x,a,i10)') 'Bounding Box at level:  ', 1
      write(6, '(1x,a,3f10.3)') 'Bounding Box Center:  ', xbox, ybox, zbox
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
   end if

   ! this is for updating of box information 

   if ( initial ) then
      cxbox(1) = xbox; cybox(1) = ybox; czbox(1) = zbox
      if ( verbose .and. pbprint ) write(6, '(a,1x,i5,1x,3f10.3)') &
         '   beginning box center at level ', 1, cxbox(1), cybox(1), czbox(1)
   else
      savxbox(1) = cxbox(1); savybox(1) = cybox(1); savzbox(1) = czbox(1)
      if ( verbose .and. pbprint ) write(6, '(a,1x,i5,1x,3f10.3)') &
         '   previous box center at level ', 1, savxbox(1), savybox(1), savzbox(1)
   end if
    
   ! set the grid dimension
    
   isfillratiosane = .false.
   do while ( .not. isfillratiosane ) 
    
   !    at least fillratio as large as the solute bounding box to search the
   !    large enough grid

   xlength = (xmax-xmin)*fillratio; ylength = (ymax-ymin)*fillratio; zlength = (zmax-zmin)*fillratio
   if ( outphi .and. phiform == 0 ) then
      xlength = max(xlength, ylength, zlength)
      ylength = xlength; zlength = xlength
   endif
   savxm(1) = nint(xlength/savh(1))
   savym(1) = nint(ylength/savh(1))
   savzm(1) = nint(zlength/savh(1))

   !    adjust for the multigrid solver (four levels) or other solver options ...
    
   if ( solvopt == 2 ) then
      savxm(1) = 16*ceiling( dble(savxm(1))/16.0d0 ) - 1
      savym(1) = 16*ceiling( dble(savym(1))/16.0d0 ) - 1
      savzm(1) = 16*ceiling( dble(savzm(1))/16.0d0 ) - 1   
   else
      savxm(1) = 2*nint( dble(savxm(1))*HALF ) + 1
      savym(1) = 2*nint( dble(savym(1))*HALF ) + 1
      savzm(1) = 2*nint( dble(savzm(1))*HALF ) + 1
   end if

   !    additional stored data and printing
    
   savxmym(1) = savxm(1)*savym(1)
   savxmymzm(1) = savxmym(1)*savzm(1)
    
   ! set grid origin
    
   if ( initial ) then
      savgox(1) = - dble(savxm(1)+1) * savh(1) * HALF + xbox
      savgoy(1) = - dble(savym(1)+1) * savh(1) * HALF + ybox
      savgoz(1) = - dble(savzm(1)+1) * savh(1) * HALF + zbox
   else
      cxbox(1) = savxbox(1) + nint( (xbox - savxbox(1))/savh(1) )*savh(1)
      cybox(1) = savybox(1) + nint( (ybox - savybox(1))/savh(1) )*savh(1)
      czbox(1) = savzbox(1) + nint( (zbox - savzbox(1))/savh(1) )*savh(1)
      if ( verbose .and. pbprint ) write(6, '(a,i5,1x,3f10.3)') &
         ' Box center corrected at level ', 1, cxbox(1), cybox(1), czbox(1)
      savgox(1) = - dble(savxm(1)+1) * savh(1) * HALF + cxbox(1)
      savgoy(1) = - dble(savym(1)+1) * savh(1) * HALF + cybox(1)
      savgoz(1) = - dble(savzm(1)+1) * savh(1) * HALF + czbox(1)
   end if
    
   ! sanity checks for electrostatic focussing
    
   isfillratiosane = .true.
   do l = 2, nfocus
      if ( savgox(l)+(savxm(l)+1)*savh(l) > savgox(l-1)+savxm(l-1)*savh(l-1) .or. &
           savgoy(l)+(savym(l)+1)*savh(l) > savgoy(l-1)+savym(l-1)*savh(l-1) .or. &
           savgoz(l)+(savzm(l)+1)*savh(l) > savgoz(l-1)+savzm(l-1)*savh(l-1) ) then
         write(6, '(a,i5)') 'PB Warn in setgrd(): focusing grid too large', l
         isfillratiosane = .false.
         fillratio = fillratio + 1
         write(6, '(a,f6.3)') ' fillratio is automatically increased to', fillratio
      end if
   end do

   if ( multiblock .and. isfillratiosane ) then ! multiblock section
      ! NOTE: tmp[x-z]m already include the paddings
      if ( minval(blkgox)         <= savgox(1)                  .or. &
           minval(blkgoy)         <= savgoy(1)                  .or. &
           minval(blkgoz)         <= savgoz(1)                  .or. &
           minval(blkgox)+tmpxm*h >= savgox(1)+savxm(1)*savh(1) .or. &
           minval(blkgoy)+tmpym*h >= savgoy(1)+savym(1)*savh(1) .or. &
           minval(blkgoz)+tmpzm*h >= savgoz(1)+savzm(1)*savh(1) ) then
         write(6, '(a)'  ) 'setgrd(): fine grid larger than coarse grid'
         write(6, '(a,f6.3)') &
              'We need a fillratio larger than', max(&
              (h*(tmpxm-1))/(xmax-xmin), &
              (h*(tmpym-1))/(ymax-ymin), &
              (h*(tmpzm-1))/(zmax-zmin))
         isfillratiosane = .false.
         fillratio = fillratio + 1
         write(6, '(a,f6.3)') 'Automatically increased fillratio to', fillratio
      end if
   end if ! end of multiblock section

   end do ! while not isfillratiosane
   if ( verbose .and. pbprint ) write(6, '(a,i5,1x,3i5)') &
      ' Grid dimension at level ', 1, savxm(1),savym(1),savzm(1)
   if ( verbose .and. pbprint ) write(6, '(a,i5,1x,3f10.3)') &
      ' Grid origin corrected at level ', 1, savgox(1),savgoy(1),savgoz(1)
 
!  end if ! if ( nfocus > 1 )

   ! :::::: for all grids ::::::

   ! if requested offseting grid, do it here
    
   if ( offx + offy + offz /= ZERO ) then
      do l = 1, nfocus
         savgox(l) = savgox(l) + offx
         savgoy(l) = savgoy(l) + offy
         savgoz(l) = savgoz(l) + offz
         if ( verbose .and. pbprint ) write(6, '(a,i5,1x,3f10.3)') &
         ' Grid origin at level after offset', l, savgox(l), savgoy(l), savgoz(2)
      end do
   end if
   if ( verbose .and. pbprint .and. .not. multiblock ) then
      write(6, '(a)') '------- VMD goodie --------'
      write(6, '(a)') '#First Level Geometric Box'
      write(6, '(a)') 'draw materials off'
      write(6, '(a)') 'draw color green'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1),savgoz(1),'" "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1),savgoz(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1),savgoz(1),'" "',&
           savgox(1),savgoy(1)+savym(1)*savh(1),savgoz(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1),savgoz(1),'" "',&
           savgox(1),savgoy(1),savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1),&
           savgoz(1),'" "',savgox(1)+savxm(1)*savh(1),&
           savgoy(1)+savym(1)*savh(1),savgoz(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1),&
           savgoz(1),'" "',savgox(1)+savxm(1)*savh(1),&
           savgoy(1),savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1)+savym(1)*savh(1),savgoz(1),&
           '" "',savgox(1)+savxm(1)*savh(1),savgoy(1)+savym(1)*savh(1),savgoz(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1)+savym(1)*savh(1),savgoz(1),&
           '" "',savgox(1),savgoy(1)+savym(1)*savh(1),savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1),savgoz(1)+savzm(1)*savh(1),&
           '" "',savgox(1)+savxm(1)*savh(1),savgoy(1),savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1),savgoy(1),savgoz(1)+savzm(1)*savh(1),&
           '" "',savgox(1),savgoy(1)+savym(1)*savh(1),savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1)+savym(1)*savh(1),&
           savgoz(1)+savzm(1)*savh(1),'" "',savgox(1)+savxm(1)*savh(1),&
           savgoy(1)+savym(1)*savh(1),savgoz(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1)+savym(1)*savh(1),&
           savgoz(1)+savzm(1)*savh(1),'" "',savgox(1),savgoy(1)+savym(1)*savh(1),&
           savgoz(1)+savzm(1)*savh(1),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(1)+savxm(1)*savh(1),savgoy(1)+savym(1)*savh(1),&
           savgoz(1)+savzm(1)*savh(1),'" "',savgox(1)+savxm(1)*savh(1),savgoy(1),&
           savgoz(1)+savzm(1)*savh(1),'"'

      write(6, '(a)') '#2nd Level Geometric Box'
      write(6, '(a)') 'draw color blue'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2),savgoz(2),'" "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2),savgoz(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2),savgoz(2),'" "',&
           savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2),savgoz(2),'" "',&
           savgox(2),savgoy(2),savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2),&
           savgoz(2),'" "',savgox(2)+savxm(2)*savh(2)+savh(2),&
           savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2),&
           savgoz(2),'" "',savgox(2)+savxm(2)*savh(2)+savh(2),&
           savgoy(2),savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),&
           '" "',savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),&
           '" "',savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2),savgoz(2)+savzm(2)*savh(2)+savh(2),&
           '" "',savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2),savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2),savgoy(2),savgoz(2)+savzm(2)*savh(2)+savh(2),&
           '" "',savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2)+savym(2)*savh(2)+savh(2),&
           savgoz(2)+savzm(2)*savh(2)+savh(2),'" "',savgox(2)+savxm(2)*savh(2)+savh(2),&
           savgoy(2)+savym(2)*savh(2)+savh(2),savgoz(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2)+savym(2)*savh(2)+savh(2),&
           savgoz(2)+savzm(2)*savh(2)+savh(2),'" "',savgox(2),savgoy(2)+savym(2)*savh(2)+savh(2),&
           savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a,3f9.3,a,3f9.3,a)') 'draw line "',&
           savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2)+savym(2)*savh(2)+savh(2),&
           savgoz(2)+savzm(2)*savh(2)+savh(2),'" "',savgox(2)+savxm(2)*savh(2)+savh(2),savgoy(2),&
           savgoz(2)+savzm(2)*savh(2)+savh(2),'"'
      write(6, '(a)') '--- End VMD goodie --------'
   endif
   
end subroutine setgrd

end subroutine pb_setgrd
