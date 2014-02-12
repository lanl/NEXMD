! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine plane here]
subroutine plane (i1, i2, i3, x, rn, a, cent)
   
   ! Subroutine determine PLANE:
   
   !  --- given three atoms i1,i2 and i3, and coordinates x
   !      returns rn(i) [i=1,2,3] components of the normalized
   !          vector normal to the plane containing the three
   !          points;  and a(i,j) [i=1,2,3, j=1..9] which
   !          are the derivatives of rn(i) with respect to
   !          the nine cartesian coordinates of points i1,i2,i3.
   !          Also returns cent[1..3], the center of the three atoms.
   
   implicit none
   integer :: i1p, i1, i2p, i2, i3p, i3, i, j
   _REAL_ :: x1, x, y1, z1, x2, y2, z2, x3, y3, z3
   _REAL_ :: tempa, tempb, tempc, ax, ay, az, anorm, rn, a, anx, cent
   dimension x(*)
   dimension a(3,9),rn(9),anx(9),cent(3)
   
   i1p = 3*i1-3
   i2p = 3*i2-3
   i3p = 3*i3-3
   
   x1 = x(i1p+1)
   y1 = x(i1p+2)
   z1 = x(i1p+3)
   x2 = x(i2p+1)
   y2 = x(i2p+2)
   z2 = x(i2p+3)
   x3 = x(i3p+1)
   y3 = x(i3p+2)
   z3 = x(i3p+3)
   
   !       ----- coefficients of the equation for the plane of atoms 1-3
   
   tempa = y1*z2 - y2*z1
   tempb = y3*z1 - y1*z3
   tempc = y2*z3 - y3*z2
   ax = tempa + tempb + tempc
   ay = -(x1*z2 - x2*z1 + x3*z1 - x1*z3 + x2*z3 - x3*z2)
   az = x1*y2 - x2*y1 + x3*y1 - x1*y3 + x2*y3 - x3*y2
   anorm = 1.d0/sqrt(ax*ax + ay*ay + az*az)
   
   !       ----- normalize to standard form for plane equation (i.e. such
   !       ----- that length of the vector "a" is unity
   
   rn(1) = ax*anorm
   rn(2) = ay*anorm
   rn(3) = az*anorm
   
   ! ----- first derivatives of ax,ay,ax w/resp. to coords. of atoms 1-3
   ! ----- first index holds ax, ay or az;
   ! ----- second index holds x1,y1...y3,z3
   ! ----- (note: these are the derivatives of the "unnormalized" a vector)
   
   a(1,1) = 0.0d0
   a(1,2) = z2 - z3
   a(1,3) = y3 - y2
   a(1,4) = 0.0d0
   a(1,5) = z3 - z1
   a(1,6) = y1 - y3
   a(1,7) = 0.0d0
   a(1,8) = z1 - z2
   a(1,9) = y2 - y1
   a(2,1) = z3 - z2
   a(2,2) = 0.0d0
   a(2,3) = x2 - x3
   a(2,4) = z1 - z3
   a(2,5) = 0.0d0
   a(2,6) = x3 - x1
   a(2,7) = z2 - z1
   a(2,8) = 0.0d0
   a(2,9) = x1 - x2
   a(3,1) = y2 - y3
   a(3,2) = x3 - x2
   a(3,3) = 0.0d0
   a(3,4) = y3 - y1
   a(3,5) = x1 - x3
   a(3,6) = 0.0d0
   a(3,7) = y1 - y2
   a(3,8) = x2 - x1
   a(3,9) = 0.0d0
   
   ! ----- first derivatives of normalization const. w/resp. to atoms 1-3
   
   do i=1,9
      anx(i) = (ax*a(1,i) + ay*a(2,i) + az*a(3,i))*anorm
   end do
   
   !       ---- finally, derivates of normalized vector:
   
   do i=1,3
      do j=1,9
         a(i,j) = (a(i,j) - rn(i)*anx(j))*anorm
      end do
   end do
   
   cent(1) = (x1+x2+x3)/3.d0
   cent(2) = (y1+y2+y3)/3.d0
   cent(3) = (z1+z2+z3)/3.d0
   
   return
end subroutine plane 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cshf here]
subroutine cshf(natom,x,f)
   
   !  Subroutine Chemical SHiFt
   
   !  -- computes constraint energies and derivatives based on
   !     chemical shift calculations outlined in Osapay and Case,
   !     J. Am. Chem. Soc. 113:9436-9444 (1991).
   
   use file_io_dat

   implicit none
   integer:: i, i1, i2, i3, iatr, ifind, ihelp, iin, iloc, ilochelp, &
        inc, innit, ino, ip, ipear, iprot, iring, j, jres, k, kp1, llist, &
        m, maxllist, n, natom, natr, nc, nca, ncount, ncurlist, nhn, &
        nmrcut, nmrcut2, nnit, nnnxloc, no, nprot, nres, nrescur, &
        nrescursave, nresprot, nring, nsum, nter, numxloc
   _REAL_ :: anorm, anx1, anx2, anx3, anx4, anx5, anx6, anx7, anx8, &
        anx9, ax, ay, az, b, bx1, bx2, bx3, bx4, bx5, bx6, bx7, bx8, bx9, &
        cent, charge, const, constca, constcagp, constnh, constside, d, &
        d1, d1cx, d1cy, d1cz, d2, d2cx, d2cy, d2cz, ddelta, deldx1, &
        deldx10, deldx11, deldx12, deldx2, deldx3, deldx4, deldx5, &
        deldx6, deldx7, deldx8, deldx9, delta, dist1, dist2, disthh, &
        dr1m3, dr2m3, drn, drph, ds12, ds12r, dshifr, dshift, eshfi, f, &
        fac, gansum, gesum, helst, hflygare, obs, obsp, ph1, ph2, ph3, &
        pp1, pp2, pp3, prob, r1, r1sq, r2, r2sq, r5ph, relc, rf, rms, rn, &
        rph2, rphi, s12p, sh_p, shav, shavp, shcut, shhm, shi, shift, &
        shp, shrang, signr, str, vp1, vp2, vp3, wt, x, x0, x1, x2, x3, &
        xrfac6, xrfact, y0, y1, y2, y3, z, z0, z1, z2, z3

#  include "nmr.h"
#  include "md.h"
#  include "def_time.h"
   integer cter
   character(len=8) namr(mring)
   logical first,nh_calc,details
   parameter (constca=-0.754d0, constside=-0.041d0, constnh=-0.55d0)
   parameter (constcagp=-0.51d0)
   parameter (hflygare=13.09d0, helst=-5.99d0)
   parameter (nh_calc=.false.)
   parameter (nmrcut=100.0d0,nmrcut2=nmrcut*nmrcut)
   
   dimension d(3,matom),signr(mring)
   dimension x(3,*),f(3,*),dshifr(9)
   dimension rn(3,mring),cent(3,mring),drn(3,9,mring),dshift(9)
   dimension natr(mring),iatr(16,mring),str(mring), &
         iprot(mshf),obs(mshf),wt(mshf),shhm(mshf),shrang(mshf), &
         shp(mshf)
   dimension pp1(mxr),pp2(mxr),pp3(mxr)
   dimension ds12r(9)
   dimension r1(3),r2(3),ds12(9),dr1m3(3),dr2m3(3)
   dimension nc(mxr),no(mxr),nca(mxr),nhn(mxr),nnit(mxr), &
         charge(matom),llist(mxr,5*mxr),numxloc(mshf),const(mshf), &
         maxllist(mxr),nresprot(mshf)
   save nring,natr,iatr,str,namr,nprot,iprot,obs,wt,shrang, &
         nc,no,nca,nhn,nnit,charge,nter,cter,ipear,numxloc, &
         llist,const,maxllist,nresprot,first,shcut,details
   dimension obsp(mshf),shavp(mshf)
   namelist/shf/ nring,natr,iatr,str,nprot,obs,wt,iprot,shcut,namr, &
         shrang,nter,cter,details
   data first /.true./
   data shcut /0.3/
   data details/.false./
   
   !=====================================================================
   
   !   First time through, set up a number of arrays, and read input data:
   
   !=====================================================================
   
   if (first) then
      
      ! (If restraint input has been redirected, open the appropriate file)
      
      iin = 5
      if (iredir(5) /= 0) then
         call amopen(36,redir(5)(1:iredir(5)),'O','F','R')
         iin = 36
         write(6,10) redir(5)(1:iredir(5))
         10 format(' Chemical shifts will be read from file: ',a)
      end if
      
      !  --- zero out the intial shrang array:
      
      do i=1,mshf
         shrang(i) = 0.0d0
      end do
      cter = 0
      nter = 1
      call nmlsrc('shf',iin,ifind)
      if (ifind == 0) goto 30
      read(iin,nml=shf,err=30,end=30)
      goto 40
      30 write(6,*) 'namelist reports error reading &shf'
      call mexit(6,1)
      40 continue
      if (cter < nter) then
         write(6,*) 'nter,cter error: ',nter,cter
         call mexit(6,1)
      end if
      if(nring > mring .or. nprot > mshf .or. nring < 1 &
            .or. nprot < 1) then
         write(6,70)
         70 format('Errors in nring or nprot')
         call mexit(6,1)
      end if
      
      !  --- multiply input str values by 4.06/0.7442943 = 5.4548
      !      which converts conventional "ring current intensities"
      !      into the value by which the geometric factor must be
      !      mutliplied to obtain shifts in ppm.
      
      do iring=1,nring
         str(iring) = str(iring)*5.4548d0
      end do
      
      !   set up arrays containing C,O,N,CA,H atom numbers
      !       for peptide contribution calculations
      
      do i=1,mxr
         nc(i)=0
         no(i)=0
         nca(i)=0
         nnit(i)=0
         nhn(i)=0
         maxllist(i)=0
         do j=1,5*mxr
            llist(i,j)=0
         end do
      end do
      !       write(6,*)'nter cter ',nter,cter
      do i=1,natom
         if (resat(i)(1:4) == 'C   ')then
            read(resat(i)(10:13),'(i4)')nres
            nc(nres)=i
            charge(i)=0.55d0
         else if (resat(i)(1:4) == 'O   ')then
            read(resat(i)(10:13),'(i4)')nres
            no(nres)=i
            charge(i)=-0.55d0
         else if (resat(i)(1:4) == 'CA  ')then
            read(resat(i)(10:13),'(i4)')nres
            nca(nres)=i
            charge(i)=0.10d0
         else if (resat(i)(1:4) == 'H   ' &
               .or. resat(i)(1:4) == 'HN  ') then
            read(resat(i)(10:13),'(i4)')nres
            nhn(nres)=i
            charge(i)=0.25d0
         else if (resat(i)(1:4) == 'N   ')then
            read(resat(i)(10:13),'(i4)')nres
            nnit(nres)=i
            charge(i)=-0.35d0
            if (resat(i)(6:8) == 'PRO')charge(i)=-0.20d0
         else
            charge(i)=0.0d0
         end if
      end do
      nrescursave=0
      ncount=0
      do i=1,nprot
         ihelp=iprot(i)
         read(resat(ihelp)(10:13),'(i4)') nrescur
         nresprot(i) = nrescur
         numxloc(i)=0
         ihelp=ihelp-1
         120 if (resat(ihelp)(1:1) == 'C' .or. &
               resat(ihelp)(1:1) == 'N' )then
            numxloc(i)=ihelp
         else
            ihelp=ihelp-1
            if (ihelp == 0)write(6,*) &
                  'Error in finding heavy atom connected to H ',ihelp
            goto 120
         end if
         
         !   Check distance of C-H or N-H
         
         disthh=0.0d0
         ihelp=iprot(i)
         ilochelp=numxloc(i)
         do iloc=1,3
            disthh=disthh+(x(iloc,ihelp)-x(iloc,ilochelp))**2
         end do
         disthh=sqrt(disthh)
         if (disthh < 0.7 .or. disthh > 1.5)write(6,*) &
               'Error:  X-H distance is ',disthh,' protnum: ',iprot(i)
#ifdef DNA_SHIFT
         const(i) = 0.0d0
#else
         const(i)=constside
         if (resat(ilochelp)(1:4) == 'CA  ') then
            if (resat(ilochelp)(6:8) == 'GLY'  .or. &
                  resat(ilochelp)(6:8) == 'PRO') then
               const(i)=constcagp
            else
               const(i)=constca
            end if
         end if
         !          if (nh_calc)const(i)=constNH
#endif
         
         !   update the llist if proton of new residue
         
         if (nrescur == nrescursave) cycle
         jres=nter
         if (jres /= nrescur) then
            ncount=ncount+1
            llist(nrescur,ncount)=nca(jres)
            ncount=ncount+1
            llist(nrescur,ncount)=nc(jres)
            ncount=ncount+1
            llist(nrescur,ncount)=no(jres)
         end if
         do jres=nter+1,cter-1
            !             if(.not.nh_calc.and.jres.eq.nrescur) go to 35
            !             if(nh_calc.and.jres.eq.nrescur-1) go to 35
            if(jres == nrescur) cycle
            ncount=ncount+1
            llist(nrescur,ncount)=nnit(jres)
            if (nhn(jres) /= 0) then
               ncount=ncount+1
               llist(nrescur,ncount)=nhn(jres)
            end if
            ncount=ncount+1
            llist(nrescur,ncount)=nca(jres)
            ncount=ncount+1
            llist(nrescur,ncount)=nc(jres)
            ncount=ncount+1
            llist(nrescur,ncount)=no(jres)
         end do
         jres=cter
         if (jres /= nrescur) then
            ncount=ncount+1
            llist(nrescur,ncount)=nnit(jres)
            if (nhn(jres) /= 0) then
               ncount=ncount+1
               llist(nrescur,ncount)=nhn(jres)
            end if
         end if
         nrescursave=nrescur
         maxllist(nrescur)=ncount
         ncount=0
      end do
      do jres=nter+1,cter
         nnit(jres-1) = nnit(jres)
      end do
      
      !     at this point, the following is valid:
      !     nc(i),no(i),nnit(i)  =  peptide group (i)  = i-th C=O and i+1-th N
      !     e.g. 1st peptide group means  C=O of 1st and N of 2nd residue
      !     loop from nter to cter-1 will go through all peptide groups
      
      
      !   Check peptide group geometry:
      
      do j=nter,cter-1
         inc=nc(j)
         ino=no(j)
         innit=nnit(j)
         dist1=0.0d0
         dist2=0.0d0
         do i=1,3
            dist1 = dist1 + (x(i,inc)-x(i,ino))**2
         end do
         do i=1,3
            dist2 = dist2 + (x(i,inc)-x(i,innit))**2
         end do
         dist1=sqrt(dist1)
         dist2=sqrt(dist2)
         if (dist1 > 2.0 .or. dist2 > 2.0)then
            write(6,*)'Error in peptide group bonds, C is ',inc
            write (6,*)'dist(C=O) = ',dist1,' dist(C-N) = ',dist2
         end if
      end do
      first = .false.
   end if
   
   !  End of "first" loop
   
   !=====================================================================
   
   !   --Initialize variables, etc.:
   
   eshf = 0.0d0
   nsum = 0
   shi = 0.0d0
   ipear = 0
   do n=1,3
      do m=1,natom
         d(n,m) = 0.0d0
      end do
   end do
   if (iprint /= 0) rewind (82)
   if (iprint /= 0) write(82,220)
   220 format(13x,'Proton       Calc.     Observed   Penalty'/ &
         10x,'---------------------------------------------')
   
   !=====================================================================
   
   ! --      loop over all protons to be calculated:
   
   !=====================================================================
   
   do i=1,nprot
      ip = iprot(i)
      if( details ) write(6,'(i5,a13)') i,resat(ip)(1:13)
      
      !-------------------------------------------------------------------
      
      !=====================================================================
      ! --    Ring currents: loop over all rings:
      !=====================================================================
      
      call timer_start(TIME_RINGCURR)
      shhm(i) = 0.0d0
      do j=1,nring
         sh_p = 0.0d0
         
         !  --   get normal to this ring
         
         if (i == 1) then
            if(namr(j)(1:3) == 'HEM'.and.namr(j)(8:8) == 'M') then
               call plane(iatr(4,j),iatr(8,j),iatr(12,j),x,rn(1,j), &
                     drn(1,1,j),cent(1,j))
            else
               call plane(iatr(1,j),iatr(3,j),iatr(5,j),x,rn(1,j), &
                     drn(1,1,j),cent(1,j))
            end if
            
            !  --   check on signr of normal vector
            
            signr(j) = 1.0d0
            d1cx = x(1,iatr(1,j)) - cent(1,j)
            d1cy = x(2,iatr(1,j)) - cent(2,j)
            d1cz = x(3,iatr(1,j)) - cent(3,j)
            d2cx = x(1,iatr(3,j)) - cent(1,j)
            d2cy = x(2,iatr(3,j)) - cent(2,j)
            d2cz = x(3,iatr(3,j)) - cent(3,j)
            vp1 = d1cy*d2cz - d1cz*d2cy
            vp2 = d1cz*d2cx - d1cx*d2cz
            vp3 = d1cx*d2cy - d1cy*d2cx
            if ((vp3*rn(3,j)+vp2*rn(2,j)+vp1*rn(1,j)) &
                  > 0.0d0) signr(j) = -1.0d0
         end if
         
         !  --     skip atoms in this ring
         
         !            (needs special code to include porphyrin rings
         !             in calculation of proximal histidine)
         
         if (resat(ip)(6:8) /= 'HEM') then
#ifdef DNA_SHIFT
            if(namr(j)(4:6) == resat(ip)(11:13) &
                  .and. resat(ip)(3:3) /= '''') cycle
#else
            if(namr(j)(4:6) == resat(ip)(11:13)) cycle
#endif
         else
            if(namr(j)(1:3) == 'HEM'  .and. &
                  namr(j)(8:8) == 'H') cycle
         end if
         
         !  --     skip rings whose center is more than NMRCUT away
         
         if (namr(j)(1:3) /= 'HEM' .or. namr(j)(8:8) /= 'M') then
            relc = (x(1,ip)-cent(1,j))**2 + (x(2,ip)-cent(2,j))**2 &
                  +(x(3,ip)-cent(3,j))**2
            if (relc > nmrcut2) cycle
         end if
         
         !=====================================================================
         
         !  --   loop over pairs of bonded atoms in the ring, and compute
         !         contribution to the shift and its derivatives:
         
         !=====================================================================
         
         do k=1,natr(j)
            kp1 = k + 1
            if (kp1 > natr(j)) kp1 = 1
            
            !          -- given proton ip, and ring atoms k and kp1, with
            !                coordinates x and ring normal vector rn;
            !                and derivatives drn of the ring normal with respect
            !                to the cartesian coordinates of atoms k and kp1;
            !          -- compute Haigh-Maillion shift contribution, and its
            !                derivatives with respect to the cartesian coordinates
            !                of ip, k and kp1
            
            do n=1,3
               r1(n) = x(n,iatr(k,j)) - x(n,ip)
               r2(n) = x(n,iatr(kp1,j)) - x(n,ip)
            end do
            
            !          -- compute triple scalar product:
            
            s12p = r1(1)*(r2(2)*rn(3,j) - r2(3)*rn(2,j)) &
                  + r1(2)*(r2(3)*rn(1,j) - r2(1)*rn(3,j)) &
                  + r1(3)*(r2(1)*rn(2,j) - r2(2)*rn(1,j))
            
            !         -- get the derivatives of s12p with respect to the cartesian
            !              coordinates of ip
            
            ds12(1) = -(rn(3,j)*(r2(2)-r1(2)) + rn(2,j)*(r1(3)-r2(3)))
            ds12(2) = -(rn(1,j)*(r2(3)-r1(3)) + rn(3,j)*(r1(1)-r2(1)))
            ds12(3) = -(rn(2,j)*(r2(1)-r1(1)) + rn(1,j)*(r1(2)-r2(2)))
            
            ds12(4) = rn(3,j)*r2(2) - rn(2,j)*r2(3)
            ds12(5) = rn(1,j)*r2(3) - rn(3,j)*r2(1)
            ds12(6) = rn(2,j)*r2(1) - rn(1,j)*r2(2)
            ds12(7) =-rn(3,j)*r1(2) + rn(2,j)*r1(3)
            ds12(8) =-rn(1,j)*r1(3) + rn(3,j)*r1(1)
            ds12(9) =-rn(2,j)*r1(1) + rn(1,j)*r1(2)
            do n=1,9
               ds12r(n) = r1(1)*(r2(2)*drn(3,n,j) - r2(3)*drn(2,n,j)) &
                     + r1(2)*(r2(3)*drn(1,n,j) - r2(1)*drn(3,n,j)) &
                     + r1(3)*(r2(1)*drn(2,n,j) - r2(2)*drn(1,n,j))
            end do
            
            !          -- get radial factors
            
            r1sq = r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3)
            d1 = sqrt(r1sq)
            r2sq = r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3)
            d2 = sqrt(r2sq)
            temp = 0.5d0*(1.d0/(r1sq*d1) + 1.0d0/(r2sq*d2))
            shift = s12p*temp
            
            do n=1,9
               dshift(n) = temp*ds12(n)
               dshifr(n) = temp*ds12r(n)
            end do
            
            do n=1,3
               dr1m3(n) = 3.0d0*r1(n)/(r1sq*r1sq*d1)
               dr2m3(n) = 3.0d0*r2(n)/(r2sq*r2sq*d2)
               dshift(n) = dshift(n) + s12p*(dr1m3(n) + dr2m3(n))*0.5d0
            end do
            do n=4,6
               dshift(n) = dshift(n) -s12p*dr1m3(n-3)*0.5d0
               dshift(n+3) = dshift(n+3) -s12p*dr2m3(n-3)*0.5d0
            end do
            
            shhm(i) = shhm(i)+ signr(j)*str(j)*shift
            sh_p    = sh_p   + signr(j)*str(j)*shift
            fac = signr(j)*str(j)
            do n=1,9
               dshift(n) = fac*dshift(n)
               dshifr(n) = fac*dshifr(n)
            end do
            if(namr(j)(1:3) == 'HEM'.and.namr(j)(8:8) == 'M') then
               do n=1,3
                  d(n,ip) = d(n,ip) - dshift(n)
                  d(n,iatr(k,j)) = d(n,iatr(k,j)) - dshift(n+3)
                  d(n,iatr(kp1,j)) = d(n,iatr(kp1,j)) - dshift(n+6)
                  d(n,iatr(4,j)) = d(n,iatr(4,j)) - dshifr(n)
                  d(n,iatr(8,j)) = d(n,iatr(8,j)) - dshifr(n+3)
                  d(n,iatr(12,j)) = d(n,iatr(12,j)) - dshifr(n+6)
               end do
            else
               do n=1,3
                  d(n,ip) = d(n,ip) - dshift(n)
                  d(n,iatr(k,j)) = d(n,iatr(k,j)) - dshift(n+3)
                  d(n,iatr(kp1,j)) = d(n,iatr(kp1,j)) - dshift(n+6)
                  d(n,iatr(1,j)) = d(n,iatr(1,j)) - dshifr(n)
                  d(n,iatr(3,j)) = d(n,iatr(3,j)) - dshifr(n+3)
                  d(n,iatr(5,j)) = d(n,iatr(5,j)) - dshifr(n+6)
               end do
            end if
         end do
         if( details ) write(6,'(20x,a8,f8.3)') namr(j), sh_p
      end do
      call timer_stop_start(TIME_RINGCURR,TIME_ELECNOE)
      if( details ) write(6,'(20x,a8,f8.3)') 'total   ', shhm(i)
      
      !   ---end of loop over rings; end of ring current section
      
      
      !=====================================================================
      
      !   ---Electrostatics contribution:
      
      !=====================================================================
      
#ifdef DNA_SHIFT
      gesum = 0.0d0
#else
      nnnxloc=numxloc(i)
      nrescur = nresprot(i)
      ncurlist=maxllist(nrescur)
      call elstat(ip,nrescur,x,nnnxloc,charge,gesum,d, &
            mxr,helst,llist,ncurlist)
      if( details ) write(6,'(20x,a8,f8.3)') 'elec    ', gesum
#endif
      call timer_stop_start(TIME_ELECNOE,TIME_ANISO)
      
      !=====================================================================
      
      !   ---Peptide group anistropy calculation:
      
      !=====================================================================
      
      gansum=0.0d0
#ifndef DNA_SHIFT
      if (i == 1) then
         do jres=nter,cter-1
            !         if (nh_calc) then
            !           if (jres.eq.nrescur) goto 30
            !         end if
            inc=nc(jres)
            ino=no(jres)
            innit=nnit(jres)
            pp1(jres)=0.5586d0*x(1,innit)+ &
                  0.5955d0*x(1,ino)-0.1542d0*x(1,inc)
            pp2(jres)=0.5586d0*x(2,innit)+ &
                  0.5955d0*x(2,ino)-0.1542d0*x(2,inc)
            pp3(jres)=0.5586d0*x(3,innit)+ &
                  0.5955d0*x(3,ino)-0.1542d0*x(3,inc)
         end do
      end if
      
      !   --- get distance from this proton to peptide group of jres:
      
      !forcevector
      do jres=nter,cter-1
         ph1=x(1,ip)-pp1(jres)
         ph2=x(2,ip)-pp2(jres)
         ph3=x(3,ip)-pp3(jres)
         rph2 = ph1*ph1 + ph2*ph2 + ph3*ph3
         
         
         !   --- skip peptide groups more than NMRCUT away:
         
         if (rph2 > nmrcut2) cycle
         
         rphi = 1.0d0/sqrt(rph2)
         r5ph = rphi/(rph2*rph2)
         
         !   calculate distance from the O,C,N plane
         
         i1=nc(jres)
         i2=no(jres)
         i3=nnit(jres)
         x0 = x(1,ip)
         y0 = x(2,ip)
         z0 = x(3,ip)
         x1 = x(1,i1)
         y1 = x(2,i1)
         z1 = x(3,i1)
         x2 = x(1,i2)
         y2 = x(2,i2)
         z2 = x(3,i2)
         x3 = x(1,i3)
         y3 = x(2,i3)
         z3 = x(3,i3)
         
         !       ----- coefficients of the equation for the plane of atoms 1-3
         
         bx1 = y2*z3 - y3*z2
         bx2 = x3*z2 - x2*z3
         bx3 = x2*y3 - x3*y2
         bx4 = y3*z1 - y1*z3
         bx5 = x1*z3 - x3*z1
         bx6 = x3*y1 - x1*y3
         bx7 = y1*z2 - y2*z1
         bx8 = x2*z1 - x1*z2
         bx9 = x1*y2 - x2*y1

         ax = bx1 + bx4 + bx7
         ay = bx2 + bx5 + bx8
         az = bx3 + bx6 + bx9
         b = x1*bx1 + x2*bx4 + x3*bx7
         anorm = sqrt(ax*ax + ay*ay + az*az)
         
         !       ----- normalize to standard form for plane equation (i.e. such
         !       ----- that length of the vector "a" is unity
         
         ax = ax/anorm
         ay = ay/anorm
         az = az/anorm
         b  =  b/anorm
         
         !       ----- delta is the desired rxn. coordinate
         
         delta = b - ax*x0 - ay*y0 - az*z0
         
         !       ----- first derivatives of ax,ay,ax w/resp. to coords. of atoms 1-3
         !       ----- first index holds ax, ay or az;
         !       ----- second index holds x1,y1...y3,z3
         !       ----- (note: these are the derivatives of the "unnormalized" a vector)
         anx1=ay*(z3 - z2) + az*(y2 - y3)
         anx2=ax*(z2 - z3) + az*(x3 - x2)
         anx3=ax*(y3 - y2) + ay*(x2 - x3)
         anx4=ay*(z1 - z3) + az*(y3 - y1)
         anx5=ax*(z3 - z1) + az*(x1 - x3)
         anx6=ax*(y1 - y3) + ay*(x3 - x1)
         anx7=ay*(z2 - z1) + az*(y1 - y2)
         anx8=ax*(z1 - z2) + az*(x2 - x1)
         anx9=ax*(y2 - y1) + ay*(x1 - x2)
         
         !       ----- first derivatives of delta w/resp. to cartesians
         
         deldx1 = - ax
         deldx2 = - ay
         deldx3 = - az
         deldx4 = - delta*anx1 + bx1 &
               - y0 * (z3-z2) - z0 * (y2-y3)
         deldx5 = - delta*anx2 + bx2 &
               - x0 * (z2 - z3) - z0 * (x3 - x2)
         deldx6 = - delta*anx3 + bx3 &
               - x0 * (y3 - y2) - y0 * (x2 - x3)
         deldx7 = - delta*anx4 + bx4 &
               - y0 * (z1 - z3) - z0 * (y3 - y1)
         deldx8 = - delta*anx5 + bx5 &
               - x0 * (z3 - z1) - z0 * (x1 - x3)
         deldx9 = - delta*anx6 + bx6 &
               - x0 * (y1 - y3) - y0 * (x3 - x1)
         deldx10 = - delta*anx7 + bx7 &
               - y0 * (z2 - z1) - z0 * (y1 - y2)
         deldx11 = - delta*anx8 + bx8 &
               - x0 * (z1 - z2) - z0 * (x2 - x1)
         deldx12 = - delta*anx9 + bx9 &
               - x0 * (y2 - y1) - y0 * (x1 - x2)
         
         gansum=gansum+rphi/(3.0d0*rph2) - delta*delta*r5ph
         
         !   Calculate derivatives
         
         !      d gan/d x = d gan/d rph * d rph/d x + d gan/d delta * d delta/d x
         !                       1.            2.          3.                4.
         
         drph = hflygare*(( -1.0d0/(rph2*rph2)) + &
               5.0d0*(delta*delta)/(rph2*rph2*rph2))
         ddelta = -2.0d0*hflygare*delta*r5ph
         d(1,ip)=d(1,ip) -drph*ph1*rphi - ddelta*deldx1
         d(1,i1)=d(1,i1) +drph*ph1*(-0.1542d0)*rphi - ddelta*deldx4/anorm
         d(1,i2)=d(1,i2) +drph*ph1*( 0.5955d0)*rphi - ddelta*deldx7/anorm
         d(1,i3)=d(1,i3) +drph*ph1*( 0.5586d0)*rphi -ddelta*deldx10/anorm
         d(2,ip)=d(2,ip) -drph*ph2*rphi - ddelta*deldx2
         d(2,i1)=d(2,i1) +drph*ph2*(-0.1542d0)*rphi - ddelta*deldx5/anorm
         d(2,i2)=d(2,i2) +drph*ph2*( 0.5955d0)*rphi - ddelta*deldx8/anorm
         d(2,i3)=d(2,i3) +drph*ph2*( 0.5586d0)*rphi -ddelta*deldx11/anorm
         d(3,ip)=d(3,ip) -drph*ph3*rphi - ddelta*deldx3
         d(3,i1)=d(3,i1) +drph*ph3*(-0.1542d0)*rphi - ddelta*deldx6/anorm
         d(3,i2)=d(3,i2) +drph*ph3*( 0.5955d0)*rphi - ddelta*deldx9/anorm
         d(3,i3)=d(3,i3) +drph*ph3*( 0.5586d0)*rphi -ddelta*deldx12/anorm
      end do
      
      !   end of loop over peptide groups contributing to proton "i"
      
      gansum=gansum*hflygare
      if( details ) write(6,'(20x,a8,f8.3)') 'peptide ', gansum
#endif
      call timer_stop(TIME_ANISO)
      
      !---------------------------------------------------------------------
      
      !   --- get anisotropy, electrostatic and constant terms into shp:
      
      shp(i) = gesum + gansum + const(i)
      if( details ) write(6,'(20x,a8,f8.3)') 'constant', const(i)
      
      !=====================================================================
      
      ! --  now get the constraint energies and update the forces
      
      !=====================================================================
      
      call timer_start(TIME_SHFDER)
      shi = shi + shhm(i) + shp(i)
      nsum = nsum + 1
      if (wt(i) > 0.0) then
         shav = shi/nsum
         if( details ) write(6,'(20x,a8,f8.3)') 'TOTALav ', shav
         if (ipnlty == 1 .or. ipnlty == 3) then
            if (shav > obs(i)+shrang(i)) then
               eshfi = wshift*wt(i)*(shav-(obs(i)+shrang(i)))
            else if (shav < obs(i)-shrang(i)) then
               eshfi = wshift*wt(i)*((obs(i)-shrang(i))-shav)
            else
               eshfi = 0.0d0
            end if
         else if (ipnlty == 2) then
            if (shav > obs(i)+shrang(i)) then
               eshfi = wshift*wt(i)*(shav-(obs(i)+shrang(i)))**2
            else if (shav < obs(i)-shrang(i)) then
               eshfi = wshift*wt(i)*(shav-(obs(i)-shrang(i)))**2
            else
               eshfi = 0.0
            end if
         else
            write(6,*) 'bad values for ipnlty:', ipnlty
            call mexit(6,1)
         end if
         eshf = eshf + eshfi
         if (ipnlty == 1 .or. ipnlty == 3) then
            if (shav > obs(i)+shrang(i)) then
               fac= wshift*wt(i)/nsum
            else if (shav < obs(i)-shrang(i)) then
               fac= -wshift*wt(i)/nsum
            else
               fac = 0.0d0
            end if
         else if (ipnlty == 2) then
            if (shav > obs(i)+shrang(i)) then
               fac= 2.0d0*wshift*wt(i)*(shav-(obs(i)+shrang(i)))/nsum
            else if (shav < obs(i)-shrang(i)) then
               fac= 2.0d0*wshift*wt(i)*(shav-(obs(i)-shrang(i)))/nsum
            else
               fac = 0.0d0
            end if
         end if
         if (iprint /= 0 .and. abs(shav-obs(i)) > shcut) &
               write(82,350) ip,resat(ip),shav,obs(i),eshfi
         350 format(5x,i4,1x,a13,3f10.2)
         ipear = ipear + 1
         if (ipear > mshf) then
            write(6,*) 'ipear: ',ipear
            call mexit(6,1)
         end if
         shavp(ipear) = shav
         obsp(ipear) = obs(i)
         shi = 0.0d0
         nsum = 0
         
         ! --  update the force array
         
         do n=1,3
            do m=1,natom
               f(n,m) = f(n,m) + fac*d(n,m)
               d(n,m) = 0.0d0
            end do
         end do
      end if

      call timer_stop(TIME_SHFDER)
      
      !=====================================================================
      
      !  --- End grand loop over all protons whose shift is to estimated:
      
   end do
   
   !=====================================================================
   
   if (iprint /= 0) then
      write(82,44) eshf
      44 format(40x,'Total shift    constraint:',f8.2)
      write(82,390)
      390 format(/21x,'#  Pearson r   rms error'/18x,35('-'))
      call pearsn(obsp,shavp,ipear,rf,prob,z,rms,xrfact,xrfac6,iuse)
      write(82,'(a18,i5,2f10.5)') 'Shift correlation:',iuse,rf,rms
   end if
   
   
   return
end subroutine cshf 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine elstat here]
subroutine elstat(ip,nrescur,x,nnnxloc,charge,gesum,d, &
      mxr,helst,llist,ncurlist)
   
   !   Subroutine calculating electrostatics contribution
   !      to chemical shifts
   
   !   - given proton ip, coordinates in x,
   !   - return shift contribution from electrostatics and
   !            its derivatives with respect to the coordinates
   
   implicit none
#  include "md.h"

   integer:: i, ip, jatom, llist, maxcurlist, mxr, ncurlist, &
        nnnxloc, nrescur
   _REAL_ :: ch1, ch2, ch3, charge, d, derconst, dot, dud1, dud2, &
        dud3, dvd1, dvd2, dwd1, dwd3, gesum, helst, hx1, hx2, hx3, rch, &
        rch2, rhx, rhx2, upar, vpar, wpar, x

   dimension x(3,*),charge(*),d(3,*)
   dimension llist(mxr,5*mxr)
   parameter (maxcurlist=1000)
   dimension rhx2(maxcurlist),rch2(maxcurlist),dot(maxcurlist), &
         ch1(maxcurlist),hx1(maxcurlist), &
         ch2(maxcurlist),hx2(maxcurlist), &
         ch3(maxcurlist),hx3(maxcurlist), &
         rch(maxcurlist),rhx(maxcurlist)
   
   if (ncurlist > maxcurlist) then
      write(6,*) 'Need to re-dimension elstat: ',ncurlist,maxcurlist
      call mexit(6,1)
   end if
   gesum = 0.0d0
   do i=1,ncurlist
      jatom=llist(nrescur,i)
      
      !  --- estimate the contribution of electrostics to the proton
      !         chemical shift
      
      !    jatom is the atom whose charge is creating the field (X)
      !    ip    is the proton whose shift is under consideration (H)
      !    nnnxloc is the heavy atom bonded to ip  (C)
      
      !    charge is the array of charges
      !    x is the array of coordinates
      
      !    ge is returned as the current contribution to the shift,
      !        but has not yet been multiplied by helst
      !    d is the derivative array which is updated by this subroutine
      
      !    helst is the constant mutliplying the electrostatic contibution
      
      ch1(i)=x(1,ip)-x(1,nnnxloc)
      hx1(i)=x(1,jatom)-x(1,ip)
      ch2(i)=x(2,ip)-x(2,nnnxloc)
      hx2(i)=x(2,jatom)-x(2,ip)
      ch3(i)=x(3,ip)-x(3,nnnxloc)
      hx3(i)=x(3,jatom)-x(3,ip)
      rhx2(i) = hx1(i)*hx1(i) + hx2(i)*hx2(i) + hx3(i)*hx3(i)
      rch2(i) = ch1(i)*ch1(i) + ch2(i)*ch2(i) + ch3(i)*ch3(i)
      dot(i) =  ch1(i)*hx1(i) + ch2(i)*hx2(i) + ch3(i)*hx3(i)
      rch(i) = sqrt(rch2(i))
      rhx(i) = sqrt(rhx2(i))
      gesum = gesum + (charge(jatom)*dot(i))/(rhx2(i)*rch(i)*rhx(i))
      !      write(6,'(4i5,2f10.3)')i,jatom,ip,nnnxloc,charge(jatom),dot(i)
   end do
   
   !  Calculate derivatives
   
   !   d ge/d x = (helst*charge/v**2 * w**4)*(vw du/dx - 3uv dw/dx -uw dv/dx)
   !                    derconst           upar       wpar       vpar
   
   
   !forcevector
   do i=1,ncurlist
      jatom=llist(nrescur,i)
      derconst = helst*charge(jatom)/(rch2(i)*rhx2(i)*rhx2(i))
      upar = rch(i)*rhx(i)
      wpar = -3.0d0*dot(i)*rch(i)
      vpar = -dot(i)*rhx(i)
      
      !      1: ip (h)   2: nnnxloc (c)  3: jatom (x)
      
      dud1 = hx1(i) - ch1(i)
      dud2 = -hx1(i)
      dud3 = ch1(i)
      dvd1 = ch1(i)/rch(i)
      dvd2 = -dvd1
      dwd1 = -hx1(i)/rhx(i)
      dwd3 = -dwd1
      d(1,ip)=d(1,ip) &
            -derconst* (upar*dud1 + wpar*dwd1 + vpar*dvd1)
      d(1,nnnxloc)=d(1,nnnxloc) &
            -derconst* (upar*dud2 + vpar*dvd2)
      d(1,jatom)=d(1,jatom) &
            -derconst* (upar*dud3 + wpar*dwd3)
      
      dud1 = hx2(i) - ch2(i)
      dud2 = -hx2(i)
      dud3 = ch2(i)
      dvd1 = ch2(i)/rch(i)
      dvd2 = -dvd1
      dwd1 = -hx2(i)/rhx(i)
      dwd3 = -dwd1
      d(2,ip)=d(2,ip) &
            -derconst* (upar*dud1 + wpar*dwd1 + vpar*dvd1)
      d(2,nnnxloc)=d(2,nnnxloc) &
            -derconst* (upar*dud2 + vpar*dvd2)
      d(2,jatom)=d(2,jatom) &
            -derconst* (upar*dud3 + wpar*dwd3)
      
      dud1 = hx3(i) - ch3(i)
      dud2 = -hx3(i)
      dud3 = ch3(i)
      dvd1 = ch3(i)/rch(i)
      dvd2 = -dvd1
      dwd1 = -hx3(i)/rhx(i)
      dwd3 = -dwd1
      d(3,ip)=d(3,ip) &
            -derconst* (upar*dud1 + wpar*dwd1 + vpar*dvd1)
      d(3,nnnxloc)=d(3,nnnxloc) &
            -derconst* (upar*dud2 + vpar*dvd2)
      d(3,jatom)=d(3,jatom) &
            -derconst* (upar*dud3 + wpar*dwd3)
      
   end do
   
   gesum=gesum*helst
   return
end subroutine elstat 
