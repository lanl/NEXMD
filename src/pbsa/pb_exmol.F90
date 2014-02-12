! <compile=optimized>
#include "copyright.h"
#  define _REAL_ double precision
#include "pb_def.h"
#include "timer.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular ion map assignment
subroutine pb_ionmap( pbverbose,ifcap,natom,iprob,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              outflag,gcrd,radi,atmsas,insas,zv,saltgrd)

   implicit none

   ! passed variables
 
   logical pbverbose
   integer ifcap
   integer natom
   _REAL_ iprob, h, gox, goy, goz
   integer xm, ym, zm, xmymzm
   integer outflag(*)
   _REAL_ gcrd(3,*), radi(*)
   integer atmsas(*), insas(*)
   _REAL_ zv(*), saltgrd(xmymzm)

   ! local variables

   integer iatm, xmymzm_ext
   _REAL_ rh, range0, range1, xi, yi, zi

   rh = 1.0d0/h

   ! local array setup
   xmymzm_ext = xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8

   ! resetting atmsas and insas

   insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0

   ! mark grid points within Stern layer as -3
    
   zv(1:xmymzm_ext) = 9999.0d0
   do iatm = 1, natom
      if (ifcap == 5 .and. outflag(iatm) == 1) cycle
      range0 = radi(iatm)
      if ( range0 == 0.0d0 ) cycle
      range1 = (range0+iprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exstsph( -3, insas, atmsas, zv )
   end do

   ! set up the ion exclusion map

   saltgrd(1:xmymzm) = 1.0d0
   call ionmap( insas, saltgrd )


contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Exlusion of ions from protein interior.
subroutine ionmap ( insas,saltgrd )

   ! Passed variables

   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ saltgrd(xm,ym,zm)

   ! Local variables

   integer i, j, k
   ! for InsightII display
   !_REAL_ g(3)

   ! for InsightII display
   !open (unit=55, file='ions.dot')
   !write (55, '("DOTS")')
   do k = 1, zm; do j = 1, ym; do i = 1, xm
      if ( insas(i,j,k) /= -4 ) then
         saltgrd(i,j,k) = 0.0d0
         ! for InsightII display
         !g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
         !write (55,'(4(f8.3,2x))') g(1:3), 300.
      end if
   end do; end do; end do
   ! for InsightII display
   !close(55)

end subroutine ionmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic stern spheres
subroutine exstsph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer dielsph
   integer insph(0:xm+1,0:ym+1,0:zm+1)
   integer inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer i, j, k
   integer lowi, lowj, lowk
   integer highi, highj, highk
   _REAL_ range2, range3!, d, d2, r

   ! mjhsieh: iatm is defined in the parent routine, which
   !          should be put into the interface instead.
   !r = range0

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               !d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
               !if ( d > r ) then
               !   d = d - r
               !   if ( insph(i,j,k) == dielsph ) then
               !      if ( d < dst(i,j,k) ) then
               !         inatm(i,j,k) = iatm; dst(i,j,k) = d
               !      end if
               !      cycle
               !   end if
               !   insph(i,j,k) = dielsph;
               !   inatm(i,j,k) = iatm; dst(i,j,k) = d
               !else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = 0.0d0
               !end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exstsph


end subroutine pb_ionmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment
subroutine pb_exmol_ses( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
              smoothopt,dprob,epsin,epsout,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              level,nfocus,&
              narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
              outflag,gcrd,acrd,radi,radip3,&
              marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
              atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
              iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)

   use poisson_boltzmann, only: ligand, multiblock, membraneopt, mthick,&
                              mctrdz, outlvlset, outmlvlset, poretype,&
                              poreradius

   use pbtimer_module
   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap, ipb, natom, smoothopt
   _REAL_ dprob, epsin, epsout
   _REAL_ h, gox, goy, goz
   integer xm, ym, zm, xmymzm, level, nfocus
   integer narcdot, maxarc
   integer nbnd, nbndx, nbndy, nbndz
   integer outflag(*)
   _REAL_ gcrd(3,*), acrd(3,*), radi(*), radip3(*)
   integer marc(*), m2narc(maxarc,*), fstarc(*), arcatm(2,*), dotarc(*)
   _REAL_ arccrd(3,*), savarc(3,*)
   integer atmsas(*), insas(*)
   _REAL_ lvlset(*), zv(*)
   _REAL_ epsx(xmymzm+ym*zm), epsy(xmymzm+xm*zm), epsz(xmymzm+xm*ym)
   integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer savbcopt(nfocus), saopt, sasopt

   ! Local variables

   integer ip, iatm, buf, nwarn, xmymzm_ext
   integer i, j, k
   integer newown
   integer ierr
   integer, allocatable :: atmx(:,:,:)
   integer, allocatable :: atmy(:,:,:)
   integer, allocatable :: atmz(:,:,:)
   _REAL_ xi, yi, zi
   _REAL_ range1, rh
   _REAL_, allocatable :: clsx(:,:,:)
   _REAL_, allocatable :: clsy(:,:,:)
   _REAL_, allocatable :: clsz(:,:,:)

        !membrane variables
        !_REAL_ mthick

        !levelset debugging variables
!		character*14 		lvlsetfilename
!		character*9 	lvlsetdataname
!		integer			lvlsetf

   allocate(atmx(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   allocate(atmy(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   allocate(atmz(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   allocate(clsx(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   allocate(clsy(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   allocate(clsz(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
   call pbtimer_start(PBTIME_PBEXMOL)
   call pbtimer_start(PBTIME_PBEXMOL_SETUP)
   ! local array setup

    !levelset debugging
    !lvlsetfilename="pbsa_lvlset.dx"
    !lvlsetdataname="level set"
    !lvlsetfn=3334

   xmymzm_ext = xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8
   if(ipb /= 4 .and. ipb /= 5) then
   epsx(1:xmymzm+ym*zm) = epsout; epsy(1:xmymzm+xm*zm) = epsout; epsz(1:xmymzm+xm*ym) = epsout
   end if
   clsx = 1.d0; clsy = 1.d0; clsz = 1.d0 
   atmx = 0; atmy = 0; atmz = 0 
   newown = 1

   call pbtimer_stop(PBTIME_PBEXMOL_SETUP)
   ! in uniform dielectric systems, nothing to do here
   ! this is most likely for some reference state calculation or development work

   if ( epsin == epsout ) return

   ! in heterogeneous dielectrics systems ...

   rh = 1.0d0/h

   if ( sasopt ==  1 ) then

      ! part a: reset atmsas and insas

      call pbtimer_start(PBTIME_PBEXMOL_PARTA)
      insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0
      call pbtimer_stop(PBTIME_PBEXMOL_PARTA)

      ! part b: mark grid points just a bit larger than SAS by 2 grid points as -2

      ! WJ:
      ! 1. if use the level set-based SES, we first need to take care of the possibly very small
      ! solvent probe since we need level set function values on not just
      ! boundary grid points but its neighobrs as well. Adding 2 grids beyound
      ! dprob makes it safer for the later search of neighbors
      ! 2. this is also done for classical SES for comparison only ...

      call pbtimer_start(PBTIME_PBEXMOL_PARTB)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+4.0d0; xi = gcrd(1,iatm);yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTB)
!write(600+mytaskid,*)insas(1:xmymzm_ext),atmsas(1:xmymzm_ext);call mexit(0,0)

      ! part c, mark grid points within SAS as 1, outside remains to be -2 from part b

      call pbtimer_start(PBTIME_PBEXMOL_PARTC)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( 2, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTC)

   ! Mengjuei: I've tested multiblock with first level ses enabled. The
   !           accuracy was improved by not much. Not worth it. (01262011)

   !XP: add density surface using sasopt == 2

   else if ( sasopt == 2 ) then
    !XP: resetting lvlst and insas
       lvlset(1:xmymzm_ext) = 0.0d0;insas(1:xmymzm_ext) = 0 
    call set_lvlst(insas,lvlset,membraneopt,mthick,mctrdz,outlvlset,outmlvlset,&
                   poretype,poreradius)


   else  if ( sasopt == 0 .and. level == nfocus ) then

      ! SES at the fine level

      ! part a: reset atmsas and insas

      call pbtimer_start(PBTIME_PBEXMOL_PARTA)
      insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0
      call pbtimer_stop(PBTIME_PBEXMOL_PARTA)
      !iepsav=0; iepsavx=0; iepsavy=0; iepsavz=0

      ! part b: mark grid points just a bit larger than SAS by 2 grid points as -2

      ! WJ:
      ! 1. if use the level set-based SES, we first need to take care of the possibly very small
      ! solvent probe since we need level set function values on not just
      ! boundary grid points but its neighobrs as well. Adding 2 grids beyound
      ! dprob makes it safer for the later search of neighbors
      ! 2. this is also done for classical SES for comparison only ...

      call pbtimer_start(PBTIME_PBEXMOL_PARTB)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+4.0d0; xi = gcrd(1,iatm);yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTB)

      ! part c, mark grid points within SAS as 1, outside remains to be -2 from part b

      call pbtimer_start(PBTIME_PBEXMOL_PARTC)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( 1, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTC)
    
      ! part d, mark grid points within VDW as 2, outside is 1 from part c
    
      call pbtimer_start(PBTIME_PBEXMOL_PARTD)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTD)
    
      ! part e, mark grid points in the contact region accessible to the solvent probe as -2
    
      call pbtimer_start(PBTIME_PBEXMOL_PARTE)
      call contact( insas, atmsas )
      call pbtimer_stop(PBTIME_PBEXMOL_PARTE)
    
      ! part f, mark grid points in the reentry region accessible to the solvent probe as -1
         
      call pbtimer_start(PBTIME_PBEXMOL_PARTF)
      buf = 2; range1 = dprob*rh + buf
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = narcdot, 1, -1
         xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
         call exresph( -1, insas, atmsas, zv(1) )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTF)
!   call extractinsas2(-1, insas)
!   stop

   ! modified VDW to approximate SES at the coarse level
 
   else
       
      ! part a, reset grid points everywhere as -2
 
      insas(1:xmymzm_ext) = -2

      ! part b, mark volume within VDW as 2
      ! 1. since outside is -2, there are only contact fractional edges
      ! 2. note that zv is set to zero because we don't want to store any
      ! atom/distance info at this level and atmsas won't be used for this level
      ! 3. we are using the modified VDW radii that have been
      ! augmented by radinc in sa_driver()

      zv(1:xmymzm_ext) = 0.0d0
      do iatm = 1, natom
         range1 = radip3(iatm)
         if ( range1 == 0.0d0 ) cycle
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv(1) )
      end do
 
   end if
 
   ! save boundary edges for eneopt=2 energy and forces
   ! Qin: always call epsbnd because bcopt == 6 requires this 
   ! WJ: epsbnd has been revised. should be free of warning in normal calls

   call pbtimer_start(PBTIME_PBEPSBND)
   call epsbnd( atmsas, insas )
   call pbtimer_stop(PBTIME_PBEPSBND)

   ! WJ: if requested, set up the level set function as signed distance to SES

   if ( ( ipb == 2 .or. ipb == 4 .or. ipb == 5 ) .and. level == nfocus .and. &
       sasopt /= 2 ) then
      lvlset(1:xmymzm_ext) = 9999.0d0
      call assignlvlset( atmsas, insas, lvlset )
   end if


   ! set up epsx, epsy and epsz maps
   ! boundary edges for eneopt=1 and frcopt=1&3 are saved inside.

   call pbtimer_start(PBTIME_PBEPSMAP)
   if( ipb /= 4 .and. ipb /= 5) call epsmap( ipb, insas, atmsas, lvlset, epsx, epsy, epsz )
   call pbtimer_stop(PBTIME_PBEPSMAP)
   call pbtimer_stop(PBTIME_PBEXMOL)
!write(600+mytaskid,*)epsx(1:xmymzm+ym*zm),epsy(1:xmymzm+xm*zm),epsz(1:xmymzm+xm*ym);call mexit(0,0)

   call pbtimer_start(PBTIME_PBCALSA)
   if ( level == nfocus .and. (.not. ligand .or. .not. multiblock)  ) then
      if ( abs(saopt) == 1 ) &
         call calc_sa1(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                       iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt)
      if ( abs(saopt) == 2 ) &
         call calc_sa2(acrd,xm,ym,zm,xmymzm,nbnd,iepsav,iepsavx,iepsavy,&
                       iepsavz,gox,goy,goz,h,sasopt,smoothopt, &
                       epsin,epsout,insas,epsx,epsy,epsz)
   end if

    !Membrane dielectric grid edge array setup. Only useable with gaussian
    !level set surface for now.
    if ( ( ipb == 2 .or. ipb == 4 .or. ipb == 5 ) .and. level == nfocus ) then
      if (membraneopt >0) then
      !   write(6,*) 'writting protein level set'; flush(6)
      !   call printlvlset(lvlset,9000)
      !  write(6,*) 'Building membrane level set'; flush(6)
         !call membrane_lvlset(lvlset,atmsas,insas)
        !Adding membrane region for insas, lvlset, and atmsas
        !Setting void checking to .true. for now. Will later be wired
        !to the appropriate variable from mdin
        !setting membrane thickness to 5 angstroms, and lvlset layer to iprob
        ! write(6,*) 'exmol dprob:',dprob
      !   call setup_membrane(lvlset,insas,&
      !          atmsas,xm,ym,zm,gox,goy,goz,h,&
      !          dprob,.true.,epsx,epsy,epsz,ipb)
      !   write(6,*) 'setting membrane epsbnd';flush(6)
      !   call epsbnd(atmsas,insas)
      !  write(6,*) 'writting new level set';flush(6)        
      !  call printlvlset(lvlset,9001)
         call set_membrane_eps(epsx,epsy,epsz,epsin,epsout,epsin,&
                        lvlset,insas,gox,goy,goz,xm,ym,zm,h,mthick,&
                        mctrdz,ipb)
      end if
    end if        

   call pbtimer_stop(PBTIME_PBCALSA)


   deallocate(atmx, stat = ierr )
   deallocate(atmy, stat = ierr )
   deallocate(atmz, stat = ierr )
   deallocate(clsx, stat = ierr )
   deallocate(clsy, stat = ierr )
   deallocate(clsz, stat = ierr )
!if ( level == nfocus ) then
!if ( ligand ) then
!   call printinsas(insas,8881,0,xm+1,0,ym+1,0,zm+1)
!   call printepsx(  epsx,9991,0,xm  ,1,ym  ,1,zm  )
!   call printepsy(  epsy,9992,1,xm  ,0,ym  ,1,zm  )
!   call printepsz(  epsz,9993,1,xm  ,1,ym  ,0,zm  )
!else
!   call prntinsas_(insas,8881,-13.5d0,11.5d0,-44.5d0,-20.5d0,10.5d0,35.5d0)
!   call prntepsx_(  epsx,9991,-13.5d0,11.0d0,-44.0d0,-21.0d0,11.0d0,35.0d0)
!   call prntepsy_(  epsy,9992,-13.0d0,11.0d0,-44.5d0,-21.0d0,11.0d0,35.0d0)
!   call prntepsz_(  epsz,9993,-13.0d0,11.0d0,-44.0d0,-21.0d0,10.5d0,35.0d0)
!   call prntinsas_(insas,8881,-16.5d0, 4.5d0,-18.5d0,  2.5d0,-18.0d0, 4.0d0)
!   call prntepsx_(  epsx,9991,-16.5d0, 4.0d0,-18.0d0,  2.0d0,-17.5d0, 3.5d0)
!   call prntepsy_(  epsy,9992,-16.0d0, 4.0d0,-18.5d0,  2.0d0,-17.5d0, 3.5d0)
!   call prntepsz_(  epsz,9993,-16.0d0, 4.0d0,-18.0d0,  2.0d0,-18.0d0, 3.5d0)
!end if
!end if
contains
subroutine printlvlset(lvlset,fn)
    !Passed Variables
    integer fn
    _REAL_ lvlset(0:xm+1,0:ym+1,0:zm+1)

    !Local Variables
    integer xi,yi,zi

    write(fn,*) "#xi,    yi,    zi,    lvlset"
    do zi=1,zm; do yi=1,ym; do xi=1,xm
        write(fn,*) xi,yi,zi,lvlset(xi,yi,zi)
    end do; end do; end do
       
end subroutine printlvlset

!the Gaussian surface
subroutine set_lvlst(ins,lvl,membraneopt,mthick,mctrdz,outlvlset,outmlvlset,&
                     poretype,porerad)

   integer ins(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ lvl(0:xm+1,0:ym+1,0:zm+1)

   integer i,j,k
   _REAL_  zi,yi,xi

   logical outlvlset, outmlvlset
   
   !membrane related vars-------------------------------------------
   !see membrane.f for membrane related routines
   _REAL_   dval
   _REAL_   lvlsum
   _REAL_   meps 
   _REAL_   mthick,mprobe,porerad
   _REAL_   poredata(0:zm+1,3) !will need to be passed later
   _REAL_   mctrdz
   integer  poretype
   integer  membraneopt   
   _REAL_, allocatable ::  mlvl(:,:,:)
   integer  mlvlAllocErr
    !debuging - for lvlset file writting
    character*9     lvlsetdataname, mlvldataname
    character*14    lvlsetfilename
    character*15    mlvlfilename
    integer         lvlsetfn,       mlvlfn
    
    if (membraneopt /= 0) then
        allocate(mlvl(0:xm+1,0:ym+1,0:zm+1),STAT=mlvlAllocErr)

        if (mlvlAllocErr /= 0) then
                write(6,*) 'PB ERROR in set_lvlst(): failed to allocate mlvl';flush(6)
                write(6,*) 'allocation error type ',mlvlAllocErr
                call mexit(6,1)
        end if
    end if
    lvlsetdataname="Level Set"; mlvldataname="Level Set"
    lvlsetfilename="pbsa_lvlset.dx"; mlvlfilename="pbsa_mlvlset.dx"
    lvlsetfn=3334
    mlvlfn=3335

    !!!write(6,*) '-debug: writting levelset'; flush(6)
   !write(6,*) '-debug: poretype = ',poretype; flush(6)

   !Hard coding membrane variables here
   !will need to make these passed in later
   mprobe = h*2.0d0
   !poretype=0 !only use 0 or 1 for now, 0= no pore, 1= cylinder
   meps=epsin
   porerad = 6.0d0 !for cylindrical pore

   if ( membraneopt /= 0 .and. zm*h < mthick ) then
      write(6,*) "PB BOMB in set_lvl(): membrane thickness wider than grid"
      write(6,*) "mthick = ",sngl(mthick),"box z dim=",sngl(zm*h)
      call mexit(6,1)  
   end if

   if (poretype ==1) then !will need to add more when other poretypes 
                          !are implemented
   ! write(6,*) '-debug: generating cylinder exclusion region data'; flush(6)
      call gen_cylinder_data(xm,ym,zm,gox,goy,goz,h,mthick,porerad,&
                    natom,poredata) 
    !!!write(6,*) '-debug: finished generating exclusion region data'; flush(6)
      !above subroutine located in membrane.f
    !do k=0,zm+1
    !write(6,*) 'k,poredata(k,1:3)',k,poredata(k,1:3);flush(6)
    !end do
   end if
   !-----------------------------------------------------------------

   do k = 0, zm+1
      zi = real(k)
        !!!write(6,*) '-debug: set_lvl(): k = ',k;flush(6)
      do j = 0, ym+1 
         yi = real(j)
         do i = 0, xm+1
            xi = real(i)
            call density_calc(xi,yi,zi,lvl(i,j,k))
            !membrane related code ------------------------------------
            if (membraneopt /= 0) then
               call membrane_density_calc(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                        mthick,mprobe,dprob,poretype,&
                        poredata,dval,mctrdz) !see membrane.f for subroutine

                mlvl(i,j,k) = 1.0d0-dval
            end if
            !---------------------------------------------------------
!           write (90,*),i,j,k
!           write (90,*),xi,yi,zi,lvl(i,j,k)
            !we will set ins before adding membrane levelset. Otherwise
            !we would have trouble locating which regions belong to solute
            !and which ones belong to membrane... not real important at this
            !point since dielectrics are the same, but will be very important
            !if they are not.
            if (membraneopt == 0) then !no membrane region to worry about
               if ( lvl(i,j,k) < 0 ) then
                  ins(i,j,k) = 1
               else
                  ins(i,j,k) = -1
               end if
            else !membrane region present
               lvlsum = 1.0d0-(dval-1.0d0*(lvl(i,j,k)-1.0d0))
               if (lvlsum < 0) then !inside either membrane or solute
                  if (meps == epsin) then !can use solute ins for both
                  !there will be no true membrane-solute interface to worry
                  !about
                     ins(i,j,k) = 1
                  else !we will need to mark grids in membrane as seperate from
                  !grids in solute. Will use ins=0 for grids in membrane
                     if (lvl(i,j,k) < 0) then !we are in the solute region
                        ins(i,j,k) = 1
                     else !we are in membrane region
                        ins(i,j,k) = 0
                     end if
                  end if
               else  !in solvent region
                  ins(i,j,k) = -1
               end if 
            end if
               !WMBS- Need to tac on dval (membrane levelset) to u
               !we do this after setting ins since it would otherwise cause
               !needless complications. The process will still be
               !a little convoluted since it isn't done directly in
               !density calc. We want to add membrane lvlset value (dval)
               !just like summing atom lvlsets.
               !But, since density_calc sets u = 1-u after summing u for atoms,
               !we will need to undo this, add dval, and then redo it.
            if (membraneopt /= 0) then
               lvl(i,j,k) = -1.0d0*(lvl(i,j,k)-1.0d0) + dval
               lvl(i,j,k) = 1.0d0-lvl(i,j,k)
            end if
         end do
      end do
   end do

   !!write(6,*) '-debug: writting level set files'; flush(6)
   if (outlvlset) then
      call gen_dx_file(xm,ym,zm,h,gox,goy,goz,&
            lvl(1:xm,1:ym,1:zm),lvlsetfilename,lvlsetfn,lvlsetdataname)
   end if
   if (outmlvlset .and. membraneopt>0) then
    !!!write(6,*) '-debug: done with set_lvlst()'; flush(6)
   !stop   

    !WMBS - Write dx format output file for level sets. Should be linked
    !to flags in the input file later to allow user control
  !   if (.true. .and. membraneopt>0) then !write membrane levelset
        call gen_dx_file(xm,ym,zm,h,gox,goy,goz,&
                   mlvl(1:xm,1:ym,1:zm),mlvlfilename,mlvlfn,mlvldataname)
  !   end if
  !   if (.true.) then !write total level set
  !      call gen_dx_file(xm,ym,zm,h,gox,goy,goz,&
  !                 lvl(1:xm,1:ym,1:zm),lvlsetfilename,lvlsetfn,lvlsetdataname)
   end if 
  ! !stop

end subroutine set_lvlst


subroutine density_calc(i,j,k,u)

   _REAL_ i,j,k
   _REAL_ u

   _REAL_ spcoef(4,4) , dash(5)
   integer l,m
   _REAL_ ia,ja,ka,dist

   data spcoef /1.000000,0.1000000,6.4999998E-02,2.9999999E-02, &
                -4.527143,-1.745714,0.2900000,-0.2542857,       &
                0.0000000,11.12571,-2.982857,0.8057143,         &
                14.83429,-18.81143,5.051429,-1.074286 /
   data dash /0.d0,0.25d0,0.5d0,0.75d0,1.d0 /

   save spcoef,dash

   u = 0.d0
   !do l = 1, nacrg ! natom
   do l = 1, natom
      ia = gcrd(1,l); ja = gcrd(2,l); ka = gcrd(3,l)
!     write(90,*),ia,ja,ka
      dist = sqrt((ia-i)**2+(ja-j)**2+(ka-k)**2) - radip3(l)/h
      dist = dist * h / 2.d0 / dprob
      
      if ( dist > 1.d0 ) cycle
      if ( dist < 0.d0 ) then
         u = u + 1.00000 - 4.527143 * dist
!        write(90,*) 0
      else
         do m = 1, 4
            if ( dist > dash(m) .and. dist <= dash(m+1) ) then
!              write(90,*) m,dash(m)
               u = u                             +&
                   spcoef(m,1)                   +&
                   spcoef(m,2)*(dist-dash(m))    +&
                   spcoef(m,3)*(dist-dash(m))**2 +&
                   spcoef(m,4)*(dist-dash(m))**3
            endif
         end do
      end if
!     write(90,*) dist,u
   end do 
   u = 1.0-u
end subroutine density_calc

subroutine printinsas(mymap,fn,lowi,highi,lowj,highj,lowk,highk)
   implicit none
   integer mymap(0:xm+1,0:ym+1,0:zm+1)
   integer i,j,k,fn,lowi,highi,lowj,highj,lowk,highk
   do i=lowi,highi
      do j=lowj,highj
         do k=lowk,highk
            write(fn,*) mymap(i,j,k)
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine printinsas
subroutine prntinsas_(mymap,fn,minx,maxx,miny,maxy,minz,maxz)
   implicit none
   integer mymap(0:xm+1,0:ym+1,0:zm+1)
   integer i,j,k,fn
   _REAL_ minx,maxx,miny,maxy,minz,maxz
   do i=0,xm
      if ( i*h+gox < minx .or. i*h+gox > maxx ) cycle
      do j=1,ym
         if ( j*h+goy < miny .or. j*h+goy > maxy ) cycle
         do k=1,zm
            if ( k*h+goz < minz .or. k*h+goz > maxz ) cycle
            write(fn,*) mymap(i,j,k)
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine prntinsas_
subroutine printepsx(mymap,fn,lowi,highi,lowj,highj,lowk,highk)
   implicit none
   _REAL_ mymap(0:xm,1:ym,1:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn,lowi,highi,lowj,highj,lowk,highk
   do i=lowi,highi
      do j=lowj,highj
         do k=lowk,highk
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine printepsx
subroutine printepsy(mymap,fn,lowi,highi,lowj,highj,lowk,highk)
   implicit none
   _REAL_ mymap(1:xm,0:ym,1:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn,lowi,highi,lowj,highj,lowk,highk
   do i=lowi,highi
      do j=lowj,highj
         do k=lowk,highk
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine printepsy
subroutine printepsz(mymap,fn,lowi,highi,lowj,highj,lowk,highk)
   implicit none
   _REAL_ mymap(1:xm,1:ym,0:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn,lowi,highi,lowj,highj,lowk,highk
   do i=lowi,highi
      do j=lowj,highj
         do k=lowk,highk
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine printepsz
subroutine prntepsx_(mymap,fn,minx,maxx,miny,maxy,minz,maxz)
   implicit none
   _REAL_ mymap(0:xm,1:ym,1:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn
   _REAL_ minx,maxx,miny,maxy,minz,maxz
   do i=0,xm
      if ( i*h+gox < minx .or. i*h+gox > maxx ) cycle
      do j=1,ym
         if ( j*h+goy < miny .or. j*h+goy > maxy ) cycle
         do k=1,zm
            if ( k*h+goz < minz .or. k*h+goz > maxz ) cycle
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine prntepsx_
subroutine prntepsy_(mymap,fn,minx,maxx,miny,maxy,minz,maxz)
   implicit none
   _REAL_ mymap(1:xm,0:ym,1:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn
   _REAL_ minx,maxx,miny,maxy,minz,maxz
   do i=0,xm
      if ( i*h+gox < minx .or. i*h+gox > maxx ) cycle
      do j=1,ym
         if ( j*h+goy < miny .or. j*h+goy > maxy ) cycle
         do k=1,zm
            if ( k*h+goz < minz .or. k*h+goz > maxz ) cycle
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine prntepsy_
subroutine prntepsz_(mymap,fn,minx,maxx,miny,maxy,minz,maxz)
   implicit none
   _REAL_ mymap(1:xm,1:ym,0:zm)
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer i,j,k,fn
   _REAL_ minx,maxx,miny,maxy,minz,maxz
   do i=0,xm
      if ( i*h+gox < minx .or. i*h+gox > maxx ) cycle
      do j=1,ym
         if ( j*h+goy < miny .or. j*h+goy > maxy ) cycle
         do k=1,zm
            if ( k*h+goz < minz .or. k*h+goz > maxz ) cycle
            write(fn,*) mymap(i,j,k)/eps0
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine prntepsz_
subroutine extractinsas(imesh, array)
   implicit none
   integer imesh
   integer array(0:xm+1,0:ym+1,0:zm+1)
   integer i, j, k
   do i = 0, xm+1; do j = 0, ym+1; do k = 0, zm+1
      if (array(i,j,k) == imesh) then
         !  H        0.00000        2.49029        0.00000
         write(2000,'(2x,a,3f15.5)') "H",gox+i*h,&
                                         goy+j*h,&
                                         goz+k*h
      endif
   end do; end do; end do
end subroutine extractinsas

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exsasph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d, d2, r

   r = rh*radi(iatm)

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
               if ( d > r ) then
                  d = d - r
                  if ( insph(i,j,k) == dielsph ) then
                     if ( d < dst(i,j,k) ) then
                        inatm(i,j,k) = iatm; dst(i,j,k) = d
                     end if
                     cycle
                  end if
                  insph(i,j,k) = dielsph;
                  inatm(i,j,k) = iatm; dst(i,j,k) = d
               else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = 0.0d0
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
          
end subroutine exsasph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
    
   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3
   _REAL_ projection(1:3),tmp(1:3),xl(1:3),dist,d2
   integer iatml

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               if ( insph(i,j,k) == dielsph ) then
                  if ( level /= nfocus ) cycle
                  iatml = inatm(i,j,k)
                  xl = gcrd(1:3,iatml)
                  tmp(1) = i - xl(1)
                  tmp(2) = j - xl(2)
                  tmp(3) = k - xl(3)
                  dist = sqrt(sum(tmp*tmp))
                  if ( dist == 0.d0 ) then
                     inatm(i,j,k) = iatm
                     cycle
                  end if
                  tmp = tmp / dist
                  projection = xl+tmp*radi(iatml)*rh
                  tmp(1) = projection(1)-xi
                  tmp(2) = projection(2)-yi 
                  tmp(3) = projection(3)-zi
                  d2 = sum(tmp*tmp)
                  if ( d2 < range1*range1 ) then
                     inatm(i,j,k) = iatm
                  end if
               else
                  insph(i,j,k) = dielsph
                  inatm(i,j,k) = iatm
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 >= 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
   
 
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute density at grid points within mol sas surf
subroutine exdensph(ip,insph,packing,dens,atmctr )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer nn
   integer insph(xm,ym,zm)
   _REAL_ packing(xm,ym,zm), dens(xm,ym,zm)
   _REAL_ atmctr(xm,ym,zm)

   ! Local variables
    
   integer  i, j, k, ii, jj, kk
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3,step
   _REAL_ d, r, d2, d2g, point, density_tmp
   integer  l , flag, ip
   
   step = 0.01d0
   d2g = 1.0d0/((2.0d0*dprob)*rh )

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 > 0.0d0 ) then
             
            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(1,ceiling(xi - range3)); highi = min(xm,floor(xi + range3))
            do i = lowi, highi

               ! no need to spline if it is outside sas

               if ( insph(i,j,k) == 1) then
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
                  d = d - r; point = d*d2g ! in the unit of (2*dprob/h)

                   ! get density by bicubic interpolation
                   
                   if( point > 1.0d0 ) cycle
                   if( point < -0.4 ) then
                      density_tmp = 3.352d0
                   else   
                      density_tmp = density( point, packing(i,j,k) )
                   end if
                   
                   dens(i,j,k) = dens(i,j,k) + density_tmp
                   !RLflag = i+xm*(j-1)+xmym*(k-1)
                   !RLif( density > density_value(1,flag) .and. density > density_value(2,flag)) then
                   !RL          density_list(2,flag)  = density_list(1,flag)   
                   !RL          density_list(1,flag)  = ip   
                   !RL          density_value(2,flag) = density_value(1,flag)   
                   !RL          density_value(1,flag) = density
                   !RLend if   
                   !RL          
                   !RLif( density < density_value(1,flag) .and. density > density_value(2,flag)) then
                   !RL          density_list(2,flag)  = ip   
                   !RL          density_value(2,flag) = density
                   !RLend if 
   !####################################################################################
   !               do ii =-40, 100
   !                 point = step*ii
   !                 if (abs(point) < 0.0000001) point = 0.0000001
   !                 call density_func(point,5.0d0,density, allc)
   !                 write(189,*)step*ii, density
   !               end do
   !               print *, 'OK-density' 
   !               stop 
                   
               end if ! ( insph(i,j,k) == 1 )

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exdensph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark contact grid points between vdw and sas surfaces
subroutine contact( insas,atmsas )
    
   integer insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   integer i, j, k, buffer, iatm, ii, jj, iarc, inside
   _REAL_ xg(3), xi(3), xj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji

   _REAL_, parameter :: small = 0.01d0
    
   buffer = 1
   !do k = buffer, zm+1-buffer; do j = buffer, ym+1-buffer; do i = buffer, xm+1-buffer
   !do k = 1, zm; do j = 1, ym; do i = 1, xm
   do k = 0, zm+1; do j = 0, ym+1; do i = 0, xm+1
       
      if ( insas(i,j,k) < 0 ) cycle
       
      xg(1) = gox + i*h; xg(2) = goy + j*h; xg(3) = goz + k*h
       
      ! this is the atom that marked this grid within sasrf, so it will be the
      ! grid's contact atom if it is marked so.
       
      iatm = atmsas(i,j,k)
      xi(1) = acrd(1,iatm); xi(2) = acrd(2,iatm); xi(3) = acrd(3,iatm)
       
      ! go through all arcs that this atom generates in circle()
      if ( fstarc(iatm) == 0 ) then
         print *, 'debug: in contact() 0 fstarc bomb!'
         print *, 'atom info:', iatm, xi(1:3), radi(iatm), radi(iatm) + dprob
         print *, 'grid info:', xg(1:3), i, j, k, insas(i,j,k) 
         stop
      end if
      inside = -2
!if ( marc(iatm) == 0 ) print *,"debug: marc == 0", i, j, k, iatm
      do ip = 1, marc(iatm)
         iarc = m2narc(ip,iatm)
!if ( iarc == 0 ) print *,"debug: iarc == 0", i, j, k, iatm, ip          

         ! generated by outer loop, i.e. the atom is iatm in circle()
          
         if ( iarc >= fstarc(iatm) ) then
            jj = arcatm(1,iarc)
            xj(1) = acrd(1,jj)
            xj(2) = acrd(2,jj)
            xj(3) = acrd(3,jj)
            cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
          
         ! generated by inner loop, i.e. the atom is jatm in circle()
          
         else
            xj = xi
            ii = arcatm(2,iarc)
            xj(1) = acrd(1,ii)
            xj(2) = acrd(2,ii)
            xj(3) = acrd(3,ii)
            cosaji = savarc(1,iarc); cosaij = savarc(2,iarc)
         end if
         rxij = savarc(3,iarc)
         xij = rxij*(xj - xi)
          
         dx = xg - xi; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         dx = xg - xj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         ! if gij < aij .and. gji < aji, this is a reentry grid
          
         if ( cosgij <= cosaij .or. cosgji <= cosaji ) cycle
         if ( insas(i,j,k) == 1 ) then
            inside = 1
         else if ( insas(i,j,k) == 2 ) then
            inside = 1
         end if
         exit
      end do
       
      if ( inside == -2 .and. insas(i,j,k) == 2 ) inside = 2
      insas(i,j,k) = inside
       
   end do; end do; end do
   
 
end subroutine contact
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within reentry surf
subroutine exresph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a renentry sphere (dielectric constant dielsph)
   ! of index iatm as dielsph. Modified from UHBD (Comp. Phys. Comm. 91:57-95,
   ! 1995) routines excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)

   ! Local variables

   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2, front, aa

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi

               if ( insph(i,j,k) == 2 ) cycle
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == -2 ) then
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) == 1 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) == 1 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
               else if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) - 1.d-9 ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
               else if ( insph(i,j,k) == 1 ) then
!print *,"debug1:",insph(26,17,7),insph(26,17,8),d2,range1,buf
                  if ( d2 < (range1 - buf)**2 ) then
                     insph(i,j,k) = dielsph;
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
                  else if ( d2 < dst(i,j,k) - 1.d-9 ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
               end if
!              insph(i,j,k) = dielsph;
!              inatm(i,j,k) = iatm; dst(i,j,k) = d2

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
   
 
end subroutine exresph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd ( atmsas,insas )
    
   use poisson_boltzmann, only: savgox, savgoy, savgoz
   implicit none
    
   integer atmsas(0:xm+1,0:ym+1,0:zm+1), insas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   logical boundary
   integer buffer, i, j, k, clstmp
   integer iarc
! WJ: used for contact surface test
!  integer iatml
!  _REAL_ xl(1:3),projection(1:3),tmp(1:3),dist
!
    
   nwarn = 0
   nbnd = 0
   buffer = 1
   !     !!!write(6,*) '-debug: pb_exmol.f: epsbnd: entering main loop';flush(6)
   !do k = buffer, zm+1-buffer; do j = buffer, ym+1-buffer; do i = buffer, xm+1-buffer
   do k = 1, zm; do j = 1, ym; do i = 1, xm
       
      ! set up condition for a boundary grid point
!if( dot_product((/savgox(2)+i*.5,savgoy(2)+j*.5,savgoz(2)+k*.5/)-(/10.5,0.,-12./),&
!                (/savgox(2)+i*.5,savgoy(2)+j*.5,savgoz(2)+k*.5/)-(/10.5,0.,-12./))&
!     < .001 ) then
!   if ( abs(insas(i,j,k)) == 2 ) then
!      write(6,*) "mjhsieh:2 ", insas(i,j,k), atmsas(i,j,k), &
!         savgox(2)+i*.5,savgoy(2)+j*.5,savgoz(2)+k*.5
!      write(6,*) "mjhsieh:2 ", i, j, k
!   else if ( abs(insas(i,j,k)) == 1 ) then
!      write(6,*) "mjhsieh:1 ", insas(i,j,k), atmsas(i,j,k), &
!         savgox(2)+i*.5,savgoy(2)+j*.5,savgoz(2)+k*.5
!      write(6,*) "mjhsieh:2 ", i, j, k
!      iarc = dotarc(atmsas(i,j,k)) ! this leads to the circle/arc
!      write(6,*) 'the two atoms forming the arc', iarc, arcatm(1:2,iarc)
!   else
!      write(6,*) "bomb.."; stop
!   end if
!end if
     ! write(6,*) 'i,j,k',i,j,k;flush(6)
      if (membraneopt==1 .and. sasopt/=2) then
         WRITE(6,*) 'PB Bomb in epsbnd(): Implicit membrane only allowed with sasopt=2'
         call mexit(6,1)
      end if
      boundary = .false.
      if ( (insas(i,j,k) ==  1 .or. insas(i,j,k) ==  2) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) == -1 .or. insas(i,j,k) == -2) .and.&
           (insas(i-1,j,k) ==  1 .or. insas(i-1,j,k) ==  2 .or. insas(i+1,j,k) ==  1 .or.&
            insas(i+1,j,k) ==  2 .or. insas(i,j-1,k) ==  1 .or. insas(i,j-1,k) ==  2 .or.&
            insas(i,j+1,k) ==  1 .or. insas(i,j+1,k) ==  2 .or. insas(i,j,k-1) ==  1 .or.&
            insas(i,j,k-1) ==  2 .or. insas(i,j,k+1) ==  1 .or. insas(i,j,k+1) ==  2 .or.&
            insas(i-1,j,k) ==  0 .or. insas(i+1,j,k) ==  0 .or. insas(i,j-1,k) ==  0 .or.&
            insas(i,j+1,k) ==  0 .or. insas(i,j,k-1) ==  0 .or. insas(i,j,k+1) ==  0) ) then
            boundary = .true.
      else if ( (insas(i,j,k) ==  0) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2 &
           ) ) then 
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
        !!!!write(6,*) '-debug: pb_exmol.f: epsbnd: found boundary node:',i,j,k;flush(6) 
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k
! WJ test
!     if ( insas(i,j,k) == 2 ) then 
!        iatml = atmsas(i,j,k)
!        xl = gcrd(1:3,iatml)
!        tmp(1) = i - xl(1)
!        tmp(2) = j - xl(2)
!        tmp(3) = k - xl(3)
!        dist = sqrt(sum(tmp*tmp))
!        tmp = tmp / dist
!        projection = xl+tmp*radi(iatml)*rh
!        do iatm = 1, natom
!           if ( iatm == iatml ) cycle
!           xl = gcrd(1:3,iatm)
!           tmp = projection - xl
!           dist = sqrt(sum(tmp*tmp))
!           if ( dist < radi(iatm)*rh ) then
!              write(6,*) 'Warning: Contact error'
!              write(6,*) i,j,k
!              write(6,*) insas(i-1,j,k), insas(i+1,j,k)
!              write(6,*) insas(i,j-1,k), insas(i,j+1,k)
!              write(6,*) insas(i,j,k-1), insas(i,j,k+1)
!              write(6,*) projection
!              write(6,'(i3,4f10.6)') iatml,gcrd(1:3,iatml),radi(iatml)*rh
!              write(6,'(i3,4f10.6)') iatm,gcrd(1:3,iatm),radi(iatm)*rh
!           end if
!        end do
!     end if
! test end
      if ( ifcap /= 0 .and. ifcap /= 5 ) then
         iepsav(4,nbnd) = 0
         cycle
      end if
 
      ! for a grid point in contact region +/- 2 or in a  solvent probe, simply use the atom/probe that
      ! marks it
    if( sasopt/=2) then  !in lvlst iepsav( 4 ) is no longer needed.
      clstmp = 0
      if ( abs(insas(i,j,k)) == 2 ) then
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
      elseif ( insas(i,j,k) == -1 ) then
         clstmp = -atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
 
      ! for a buried reentry grid point, find the atom that marked its neighoring exposed reentry
      ! grid points. Note that this may not be possible when grid spacing is large
 
      else if ( insas(i,j,k) == 1 ) then
         clstmp = -atmsas(i,j,k)
         if ( clstmp == 0 ) then
            write(6,*) 'PB Info: Close contact cannot be found. fndcls() is called'
            clstmp = fndcls( i, j, k, insas, atmsas )
            nwarn = nwarn + 1
         end if

      end if
      iepsav(4,nbnd) = clstmp
    end if
   end do; end do; end do
   !!!write(6,*) '-debug: pb_exmol.f: epsbnd: done with main loop';flush(6)
   if ( nwarn > 0 ) then
      if ( pbverbose .and. level == nfocus ) write(6, '(a,i4)') &
      'PB Warning in epsbnd(): No neighbor found for boundary grids total:', nwarn
   end if
 
 
end subroutine epsbnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls( i,j,k,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls, i, j, k
   integer  insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - dprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls = clsatm(ip)
   end do

 
end function fndcls
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine assignlvlset( atmsas,insas,u )

   implicit none

   ! passed variables

   integer atmsas(0:xm+1,0:ym+1,0:zm+1), insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer i,j,k,l,buffer
   _REAL_ dist,d2,d

   rh = 1.0d0 / h
   dist = dprob*rh
   buffer = 2
   do l = 1, nbnd
      do k = iepsav(3,l) - buffer, iepsav(3,l) + buffer
         do j = iepsav(2,l) - buffer, iepsav(2,l) + buffer
            do i = iepsav(1,l) - buffer, iepsav(1,l) + buffer
               if ( atmsas(i,j,k) == 0 ) then
                  write(6,'(a,3i5)') 'PB Bomb in assignlvlset(): no atmsas', i,j,k
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l)+1,iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l)-1,iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)+1,iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)-1,iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)+1)
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)-1)
!                 call mexit(6,1)
               end if
               if ( abs(insas(i,j,k)) == 2 ) then 
                  iatm = atmsas(i,j,k)
                  xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d = sqrt(d2)
!                 u(i,j,k) = dist - ( (radi(iatm)+dprob)*rh - d )
                  u(i,j,k) = - radi(iatm)*rh + d 
               else if ( abs(insas(i,j,k)) == 1 ) then
                  iatm = atmsas(i,j,k)
                  xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d =sqrt(d2)
                  u(i,j,k) = dist - d
               else
                  write(6,'(a,4i5)') 'PB Bomb in assignlvlset(): illegal insas flag', i,j,k, insas(i,j,k)
                  call mexit(6,1)
               end if
            end do
         end do
      end do
   end do


end subroutine assignlvlset
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap( ipb,insas,atmsas,u,epsx,epsy,epsz )

   implicit none

   ! passed variables

   integer ipb
   integer insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsx(0:xm,1:ym,1:zm)
   _REAL_ epsy(1:xm,0:ym,1:zm)
   _REAL_ epsz(1:xm,1:ym,0:zm)

   ! local variables

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   integer x_flag, y_flag, z_flag
   _REAL_ epsint, epsint0

   ! debugging variables
   character*12   epszfilename
   character*4    epszdataname
   integer        epszfn


   epszfilename="pbsa_epsz.dx"
   epszdataname="epsz"
   epszfn = 3333
   ! set default value for simple harmonic average

   epsint0 = 2.0d0*epsin*epsout/(epsin+epsout)
   epsint = epsint0

   ! initialize boundary edge counters

   nbndx = 0; nbndy = 0; nbndz = 0
   x_flag = 0; y_flag = 0; z_flag = 0 ! local copies for the level set version

   ! the fraction caused by atom or probe

   do k = 0, zm; do j = 0, ym; do i = 0, xm
      a = insas(i,j,k)
      b = insas(i+1,j,k)
      a1 = atmsas(i,j,k)
      b1 = atmsas(i+1,j,k)
      if ( j == 0 .or. k == 0 ) then
         continue
      else if ( sign(a,b) == a ) then
         if ( a > 0 ) then
            epsx(i,j,k) = epsin
         end if
      else
!        if ( smoothopt > 0 .and. level == nfocus ) then
         if ( level == nfocus ) then
            if ((ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. b == 2)) .or. ipb == 4 .or. ipb ==5)&
             .and. sasopt /= 2) & 
            call epsfracx  (i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout  )
            if (((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. b == 1)) .or. ipb == 4 .or. ipb == 5 ).and. nbndx > x_flag) &
             .or. (ipb == 2  .and. sasopt ==2)) then
            call epsfracx_r(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout,u)
            x_flag = x_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
            if ( smoothopt == 0 ) epsint = epsint0
         end if
         epsx(i,j,k) = epsint
      end if
      c = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( i == 0 .or. k == 0 ) then
         continue
      else if ( sign(a,c) == a ) then
         if ( a > 0 ) then
            epsy(i,j,k) = epsin
         end if
      else
!        if ( smoothopt > 0 .and. level == nfocus ) then
         if ( level == nfocus ) then
            if ((ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. c == 2)) .or. ipb == 4  .or. ipb == 5)&
            .and. sasopt /= 2) &
            call epsfracy  (i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout  )
            if (((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. c == 1)) .or. ipb == 4 .or. ipb == 5) .and. nbndy > y_flag)&
              .or. (ipb ==2 .and. sasopt ==2)) then
            call epsfracy_r(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout,u)
            y_flag = y_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
            if ( smoothopt == 0 ) epsint = epsint0
         end if
         epsy(i,j,k) = epsint
      end if
      d = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( i == 0 .or. j == 0 ) then
         continue
      else if ( sign(a,d) == a ) then
         if ( a > 0 ) then
            epsz(i,j,k) = epsin
         end if
      else
!        if ( smoothopt > 0 .and. level == nfocus ) then
         if ( level == nfocus ) then
            if ((ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. d == 2)) .or. ipb == 4  .or. ipb == 5)&
               .and. sasopt /= 2) &
            call epsfracz  (i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout  )
            if (((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. d == 1)) .or. ipb == 4  .or. ipb == 5) .and. nbndz > z_flag)&
              .or. (ipb == 2 .and. sasopt ==2 )) then
            call epsfracz_r(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout,u)
            z_flag = z_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
            if ( smoothopt == 0 ) epsint = epsint0
         end if
         epsz(i,j,k) = epsint
      end if

        !write(6,*) 'i,j,k, epsz',i,j,k,epsz(i,j,k) 

      ! checking the sanity of epsx, epsy, and epsz

      if ( j == 0 .or. k == 0 ) then
         continue
      else if ( epsx(i,j,k) < epsin .or. epsx(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsx out of range', i,j,k
         write(6,*) 'epsin,epsout,epsx: ',epsin,epsout,epsx(i,j,k)
         call mexit(6,1)
      end if
      if ( i == 0 .or. k == 0 ) then
         continue
      else if ( epsy(i,j,k) < epsin .or. epsy(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsy out of range', i,j,k
         write(6,*) 'epsin,epsout,epsy:',epsin,epsout,epsy(i,j,k)
         call mexit(6,1)
      end if
      if ( i == 0 .or. j == 0 ) then
         continue
      else if ( epsz(i,j,k) < epsin .or. epsz(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsz out of range', i,j,k
         write(6,*) 'epsin,epsout,epsz:',epsin,epsout,epsz(i,j,k)
         call mexit(6,1)
      end if

   end do; end do; end do

   ! legacy checking output 
   !do k = 1, zm
   !   write(20, *) 'plane', k
   !do j = 1, ym
   !   write(20, '(100f6.1)') epsx(1:xm,j,k)/eps0
   !end do
   !end do
   !do k = 1, zm
   !   write(21, *) 'plane', k
   !do i = 1, xm
   !   write(21, '(100f6.1)') epsy(i,1:ym,k)/eps0
   !end do
   !end do
   !do j = 1, ym
   !   write(22, *) 'plane', j
   !do i = 1, xm
   !   write(22, '(100f6.1)') epsz(i,j,1:zm)/eps0
   !end do
   !end do

end subroutine epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)
   _REAL_ eps0
   eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = range3 + xi - REAL(i)
         fedgex(nbndx) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = range3 - xi + REAL(i+1)
         fedgex(nbndx) = 1 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + xi - REAL(i)
         xg(1) = gox + h*i + h*aa
         xg(2) = goy + h*j
         xg(3) = goz + h*k
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndx = nbndx + 1
            fedgex(nbndx) =  aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
            iepsavx(4,nbndx) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         nbndx = nbndx + 1
         aa = xi-range3-REAL(i)
         fedgex(nbndx) = aa
         if ( newown == 1 ) then
            fedgex(nbndx) = clsx(i,j,k)
            iatm = atmx(i,j,k)
            aa = clsx(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm x3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - xi + REAL(i+1)
         xg(1) = gox + h*i + h*(1-aa)
         xg(2) = goy + h*j
         xg(3) = goz + h*k
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndx = nbndx + 1
            fedgex(nbndx) = 1-aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
            iepsavx(4,nbndx) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         nbndx = nbndx + 1
         aa = REAL(i+1)-xi-range3
         fedgex(nbndx) = 1-aa
         if ( newown == 1 ) then
            fedgex(nbndx) = 1.d0-clsx(i,j,k)
            iatm = atmx(i,j,k)
            aa = clsx(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm x4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = xi+range3-REAL(i)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i + h*aa
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(yi-j)**2
         if ( front >= 0.0d0 ) then
            range3 = sqrt(front)
            aa = xi-range3-REAL(i)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) =  aa
               if ( newown == 1 ) then
                  fedgex(nbndx) = clsx(i,j,k)
                  iatm = atmx(i,j,k)
                  aa = clsx(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm x5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(i+1)-xi+range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i + h*(1-aa)
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) = 1 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(yi-j)**2
         if ( front >= 0.0d0 ) then
            range3 = sqrt(front)
            aa = REAL(i+1)-xi-range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndx = nbndx + 1
               fedgex(nbndx) = 1 - aa
               if ( newown == 1 ) then
                  fedgex(nbndx) = 1.d0-clsx(i,j,k)
                  iatm = atmx(i,j,k)
                  aa = clsx(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm x6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(yi-j)**2
      nbndx = nbndx + 1
      range3 = sqrt(front)
      aa = REAL(i+1)-xi-range3
      fedgex(nbndx) = 1-aa

      if ( newown == 1 ) then
         fedgex(nbndx) = 1.d0-clsx(i,j,k)
         iatm = atmx(i,j,k)
         aa = clsx(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm x1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
      iepsavx(4,nbndx) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(yi-j)**2
      nbndx = nbndx + 1
      range3 = sqrt(front)
      aa = xi-range3-REAL(i)
      fedgex(nbndx) = aa

      if ( newown == 1 ) then
         fedgex(nbndx) = clsx(i,j,k)
         iatm = atmx(i,j,k)
         aa = clsx(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm x2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
      iepsavx(4,nbndx) = -iatm
   end if

end subroutine epsfracx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)
   _REAL_ eps0
   eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = range3 + yi - REAL(j)
         fedgey(nbndy) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = range3 - yi + REAL(j+1)
         fedgey(nbndy) = 1 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + yi - REAL(j)
         xg(1) = gox + h*i
         xg(2) = goy + h*j + h*aa
         xg(3) = goz + h*k
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndy = nbndy + 1
            fedgey(nbndy) =  aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
            iepsavy(4,nbndy) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         nbndy = nbndy + 1
         aa = yi-range3-REAL(j)
         fedgey(nbndy) = aa
         if ( newown == 1 ) then
            fedgey(nbndy) = clsy(i,j,k)
            iatm = atmy(i,j,k)
            aa = clsy(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm y3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - yi + REAL(j+1)
         xg(1) = gox + h*i
         xg(2) = goy + h*j + h*(1-aa)
         xg(3) = goz + h*k
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndy = nbndy + 1
            fedgey(nbndy) = 1-aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
            iepsavy(4,nbndy) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         nbndy = nbndy + 1
         aa = REAL(j+1)-yi-range3
         fedgey(nbndy) = 1-aa
         if ( newown == 1 ) then
            fedgey(nbndy) = 1.d0-clsy(i,j,k)
            iatm = atmy(i,j,k)
            aa = clsy(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm y4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(xi-i)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = yi+range3-REAL(j)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j + h*aa
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(xi-i)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = yi-range3-REAL(j)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) =  aa
               if ( newown == 1 ) then
                  fedgey(nbndy) = clsy(i,j,k)
                  iatm = atmy(i,j,k)
                  aa = clsy(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm y5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(xi-i)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(j+1)-yi+range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j + h*(1-aa)
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = 1 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(xi-i)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = REAL(j+1)-yi-range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndy = nbndy + 1
               fedgey(nbndy) = 1 - aa
               if ( newown == 1 ) then
                  fedgey(nbndy) = 1.d0-clsy(i,j,k)
                  iatm = atmy(i,j,k)
                  aa = clsy(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm y6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(xi-i)**2
      nbndy = nbndy + 1
      range3 = sqrt(front)
      aa = REAL(j+1)-yi-range3
      fedgey(nbndy) = 1-aa

      if ( newown == 1 ) then
         fedgey(nbndy) = 1.d0-clsy(i,j,k)
         iatm = atmy(i,j,k)
         aa = clsy(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm y1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
      iepsavy(4,nbndy) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(xi-i)**2
      nbndy = nbndy + 1
      range3 = sqrt(front)
      aa = yi-range3-REAL(j)
      fedgey(nbndy) = aa

      if ( newown == 1 ) then
         fedgey(nbndy) = clsy(i,j,k)
         iatm = atmy(i,j,k)
         aa = clsy(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm y2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
      iepsavy(4,nbndy) = -iatm
   end if

end subroutine epsfracy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)
   _REAL_ eps0
   eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = range3 + zi - REAL(k)
         fedgez(nbndz) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = range3 - zi + REAL(k+1)
         fedgez(nbndz) = 1 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + zi - REAL(k)
         xg(1) = gox + h*i
         xg(2) = goy + h*j
         xg(3) = goz + h*k + h*aa
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndz = nbndz + 1
            fedgez(nbndz) =  aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
            iepsavz(4,nbndz) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
         nbndz = nbndz + 1
         aa = zi-range3-REAL(k)
         fedgez(nbndz) = aa
         if ( newown == 1 ) then
            fedgez(nbndz) = clsz(i,j,k)
            iatm = atmz(i,j,k)
            aa = clsz(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm z3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - zi + REAL(k+1)
         xg(1) = gox + h*i
         xg(2) = goy + h*j
         xg(3) = goz + h*k + h*(1-aa)
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndz = nbndz + 1
            fedgez(nbndz) = 1-aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
            iepsavz(4,nbndz) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
         nbndz = nbndz + 1
         aa = REAL(k+1)-zi-range3
         fedgez(nbndz) = 1-aa
         if ( newown == 1 ) then
            fedgez(nbndz) = 1.d0-clsz(i,j,k)
            iatm = atmz(i,j,k)
            aa = clsz(i,j,k)
            if ( iatm == 0 ) then
               write(6,*) "wrong iatm z4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(xi-i)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = zi+range3-REAL(k)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j
            xg(3) = goz + h*k + h*aa
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(xi-i)**2-(yi-j)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = zi-range3-REAL(k)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) =  aa
               if ( newown == 1 ) then
                  fedgez(nbndz) = clsz(i,j,k)
                  iatm = atmz(i,j,k)
                  aa = clsz(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm z5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(xi-i)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(k+1)-zi+range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j
            xg(3) = goz + h*k + h*(1-aa)
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = 1 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(xi-i)**2-(yi-j)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = REAL(k+1)-zi-range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndz = nbndz + 1
               fedgez(nbndz) = 1 - aa
               if ( newown == 1 ) then
                  fedgez(nbndz) = 1.d0-clsz(i,j,k)
                  iatm = atmz(i,j,k)
                  aa = clsz(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,*) "wrong iatm z6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(xi-i)**2-(yi-j)**2
      nbndz = nbndz + 1
      range3 = sqrt(front)
      aa = REAL(k+1)-zi-range3
      fedgez(nbndz) = 1-aa

      if ( newown == 1 ) then
         fedgez(nbndz) = 1.d0-clsz(i,j,k)
         iatm = atmz(i,j,k)
         aa = clsz(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm z1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
      iepsavz(4,nbndz) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(xi-i)**2-(yi-j)**2
      nbndz = nbndz + 1
      range3 = sqrt(front)
      aa = zi-range3-REAL(k)
      fedgez(nbndz) = aa

      if ( newown == 1 ) then
         fedgez(nbndz) = clsz(i,j,k)
         iatm = atmz(i,j,k)
         aa = clsz(i,j,k)
         if ( iatm == 0 ) then
            write(6,*) "wrong iatm z2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
      iepsavz(4,nbndz) = -iatm
   end if

end subroutine epsfracz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flag_value(xg,iatm,flag)

   use solvent_accessibility, only: ntri, triarc

   implicit none

   ! passed variables

   _REAL_ xg(3)
   integer iatm, flag

   ! local variables

   integer ii, jj, iarc
   _REAL_ xii(3), xjj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji, arcpos(3)
     
   arcpos(1:3) = arccrd(1:3,iatm)

   if ( iatm > narcdot - ntri ) then
      do ip = 1, 3
         iarc = triarc(ip,iatm-narcdot+ntri)
         jj = arcatm(1,iarc)
         xjj(1) = acrd(1,jj)
         xjj(2) = acrd(2,jj)
         xjj(3) = acrd(3,jj)
         ii = arcatm(2,iarc)
         xii(1) = acrd(1,ii)
         xii(2) = acrd(2,ii)
         xii(3) = acrd(3,ii)
         cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
         rxij = savarc(3,iarc)
         xij = rxij*(xjj - xii)

         dx = xg - xii; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

         dx = xg - xjj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

         if ( cosgij > cosaij .and. cosgji > cosaji ) then
            flag = 1
            exit
         else
            flag = 0
         end if
      end do
   else
      iarc = dotarc(iatm)

      jj = arcatm(1,iarc)
      xjj(1) = acrd(1,jj)
      xjj(2) = acrd(2,jj)
      xjj(3) = acrd(3,jj)
      ii = arcatm(2,iarc)
      xii(1) = acrd(1,ii)
      xii(2) = acrd(2,ii)
      xii(3) = acrd(3,ii)
      cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
      rxij = savarc(3,iarc)
      xij = rxij*(xjj - xii)

      dx = xg - xii; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
      cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

      dx = xg - xjj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
      cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

      if ( cosgij > cosaij .and. cosgji > cosaji ) then
         flag = 1
      else
         flag = 0
      end if
   end if

end subroutine flag_value
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for x-edges with the level set function
subroutine epsfracx_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0

   if ( a > 0 ) then
      x1 = dble(i-1)
      x2 = dble(i  )
      x3 = dble(i+1)
      f1 = u(i-1,j,k)
      f2 = u(i  ,j,k)
      f3 = u(i+1,j,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i+1,j,k) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(i  )
      x2 = dble(i+1)
      x3 = dble(i+2)
      f1 = u(i  ,j,k)
      f2 = u(i+1,j,k)
      f3 = u(i+2,j,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(i)
   else
      aa = dble(i+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if(sasopt==2) nbndx=nbndx+1
   fedgex(nbndx) = t - dble(i)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)


end subroutine epsfracx_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for y-edges with the level set function
subroutine epsfracy_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(j-1)
      x2 = dble(j  )
      x3 = dble(j+1)
      f1 = u(i,j-1,k)
      f2 = u(i,j  ,k)
      f3 = u(i,j+1,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j+1,k) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(j  )
      x2 = dble(j+1)
      x3 = dble(j+2)
      f1 = u(i,j  ,k)
      f2 = u(i,j+1,k)
      f3 = u(i,j+2,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(j)
   else
      aa = dble(j+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if(sasopt==2) nbndy=nbndy+1
   fedgey(nbndy) = t - dble(j)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)

end subroutine epsfracy_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for z-edges with the level set function
subroutine epsfracz_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(k-1)
      x2 = dble(k)
      x3 = dble(k+1)
      f1 = u(i,j,k-1)
      f2 = u(i,j,k)
      f3 = u(i,j,k+1)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j,k+1) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(k)
      x2 = dble(k+1)
      x3 = dble(k+2)
      f1 = u(i,j,k)
      f2 = u(i,j,k+1)
      f3 = u(i,j,k+2)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(k)
   else
      aa = dble(k+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if(sasopt == 2) nbndz=nbndz+1
   fedgez(nbndz) =  t - dble(k)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)


end subroutine epsfracz_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ root returns the approximated root between x0 and x1 if f0*f1 <=0 using
!+ quadratic interpolation
subroutine root(x0,x1,x2,f0,f1,f2,t0)

   implicit none

   ! passed variables

   _REAL_ x0,x1,x2,f0,f1,f2,t0

   ! local variables

   _REAL_ b,c,a0,b0,c0,t,r1,r2

   b = (f0-f1)/(x0-x1)
   c = f2 - f1 - b*(x2-x1)
   c = c/( (x2-x0)*(x2-x1))

   a0 = c
   b0 = b - c*(x0+x1)
   c0 = f1 -b*x1 + c*x0*x1

   if ( a0 == 0 ) then
      t0 = -c0/b0
      return
   end if

   t = b0*b0 - 4.0d0*a0*c0

   ! If t <=0, must be double root t is close to zero

   if ( t <= 0.0d0 ) then
      t0 = -b0/(2.0d0*a0)
      return
   end if

   t = sqrt(t)
   if ( b0 >= 0.0d0 ) then
      r1 = (-b0-t)/(2.0d0*a0)
   else
      r1 = (-b0+t)/(2.0d0*a0)
   end if

   r2 = -b0/a0-r1

   if ( x0 <= r1 + 1.0d-7 .and. r1 <= x1+1.0d-7 ) then
      t0 = r1
   else
      t0 = r2
   end if

   if ( x0 > t0 ) t0 = x0
   if ( x1 < t0 ) t0 = x1


end subroutine root
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute density function given distance and packing index
function density( distance,packing )

   implicit none

   ! common block variables
    
   _REAL_ allc(4,4,-2:4,0:14)
   common /density_coefficient/ allc

   ! passed variables

   _REAL_ density
   _REAL_ distance, packing

   ! local variables

   integer i, j, k, l
   _REAL_ u, v

   packing = 1.0d0 ! no need to use packing for now ...
 
   i = floor( distance/0.20d0 )
   j = floor( packing )

   if ( distance >= 0.0d0 ) then
      u = mod(distance/0.2d0, 1.0d0)
   else
      u = 1.0d0 - mod(abs(distance/0.2d0), 1.0d0)
   end if

   if ( packing >= 0.0d0 ) then
      v = mod(packing, 1.d0)
   else
      v = 1.0d0 - mod(abs(packing), 1.0d0)
   end if

   density = 0.0d0
   do l = 1, 4
      do k = 1, 4
         density = density + allc(k,l,i,j)*( u**(k-1) )*( v**(l-1) )
      end do
   end do


end function density


end subroutine pb_exmol_ses

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set up coefficients for bicubic interpolation
subroutine density_init( )
    
   implicit none

   ! common block variables
    
   _REAL_ allc(4,4,-2:4,0:14)
   common /density_coefficient/ allc

   ! local variables

   integer i, j, k, l
   _REAL_ coef(4,4)
   _REAL_ f(4), f1(4), f12(4), f2(4)
   _REAL_ surf(-3:6,-1:16)
   _REAL_ surf1d(-2:5,-1:16)
   _REAL_ surf2d(-3:6,0:15)
   _REAL_ surf12d(-2:5,0:15)
   _REAL_ delta_dist, delta_pack

   !RLwrite(6,*) 'PB Info: Setting up spline coefficients for the density function'
   !RLopen(999,file='coef.dat')
   !RL
   !RL set up surf()
   !RL
   !RLdo j = -1, 16
   !RL   read(999,*) surf(-3:6,j)
   !RLend do
   !RL
   !RLclose(999)
 
   delta_dist = 0.20d0
   delta_pack = 1.0d0
 
   ! set up surf1d()
    
   do j = -1, 16
      do i = -2, 5
         surf1d(i,j) = (surf(i+1,j)-surf(i-1,j))/(delta_dist*2.0d0)
      end do
   end do
    
   ! set up surf2d()
    
   do j = 0, 15
      do i = -3, 6
         surf2d(i,j) = (surf(i,j+1)-surf(i,j-1))/(delta_pack*2.0d0)
      end do
   end do
    
   ! set up surf12d()
    
   do j = 0, 15
      do i = -2, 5
         surf12d(i,j) = (surf(i+1,j+1)-surf(i+1,j-1)-surf(i-1,j+1)+surf(i-1,j-1))/(delta_dist*delta_pack*4.0d0)
      end do
   end do
   
   ! now set up coefficients for each interpolation rectangle of x and y ...

   ! note the following order for surf(x, y): x changes first and y changes second, i.e. the fortran way
   ! both the orders for i,j and k,l are changed

   do j = 0, 14
      do i = -2, 4
          
         ! set up f
          
         f(1) = surf(i,j)
         f(2) = surf(i+1,j)
         f(3) = surf(i+1,j+1)
         f(4) = surf(i,j+1)
           
         ! set up f first derivitive of i
          
         f1(1) = surf1d(i,j)
         f1(2) = surf1d(i+1,j)
         f1(3) = surf1d(i+1,j+1)
         f1(4) = surf1d(i,j+1)
          
         ! set up f first derivitive of j
          
         f2(1) = surf2d(i,j)
         f2(2) = surf2d(i+1,j)
         f2(3) = surf2d(i+1,j+1)
         f2(4) = surf2d(i,j+1)
          
         ! set up f cross derivitive of i and j
         
         f12(1) = surf12d(i,j)
         f12(2) = surf12d(i+1,j)
         f12(3) = surf12d(i+1,j+1)
         f12(4) = surf12d(i,j+1)
         
         ! ready to call the set up routine ...
          
         call bicubic_coef(f,f1,f2,f12,delta_dist,delta_pack,coef)
          
         ! store the coefficients in the common block for later ...
          
         do l = 1,4
            do k = 1,4
               allc(l,k,i,j) = coef(l,k)
            end do
         end do
          
      end do
   end do


contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute coefficients for the bicubic interpolation
subroutine bicubic_coef(f,f1,f2,f12,d1,d2,coef)
   !
   ! the algorithm is documented in Numerical Recipes, 2nd edition, 1992.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   _REAL_ f(4), f1(4), f12(4), f2(4)
   _REAL_ d1, d2
   _REAL_ coef(4,4)

   ! local variables

   integer i, j, k
   _REAL_ d1d2, coef_tmp, ctmp(16), tmp(16)
   _REAL_ weight(16,16)

   ! weighting factors from Numerical Recipes

   !data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*0, &
   !9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4, &
   !1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0, &
   !-6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2, &
   !10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4, &
   !-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0, &
   !2,-2,2*0,-1,1/
   data weight/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,&
                0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,&
                0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,&
                0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,&
                0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,&
                0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,&
                0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,&
                0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,&
                0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,&
                0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1,&
                0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1 /

   ! setting up the working arrays

   d1d2 = d1*d2

   do i = 1, 4
      tmp(i) = f(i)
      tmp(i+4) = f1(i)*d1
      tmp(i+8) = f2(i)*d2
      tmp(i+12) = f12(i)*d1d2
   end do

   ! computing the coefficients

   do i = 1, 16
      ctmp(i) = 0.0d0
      do j = 1, 16
         ctmp(i) = ctmp(i) + weight(i,j)*tmp(j)
      end do
   end do

   ! saving it into the 2-d arragy for later

   k = 0
   do i = 1, 4
      do j = 1, 4
         k = k + 1
         coef(i,j) = ctmp(k)
      end do
   end do


end subroutine bicubic_coef


end subroutine density_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment.
subroutine pb_exmol_cap( pbverbose,ifcap )

   use poisson_boltzmann    
   use solvent_accessibility

   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap
 
   ! Local variables

   logical ses 
   integer ip, iatm, nwarn, xmymzm_ext
   _REAL_ xi, yi, zi
   _REAL_ range1, rh

   ! local array setup
   xmymzm_ext = xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8
 
   epsx(1:xmymzm+ym*zm) = epsout; epsy(1:xmymzm+xm*zm) = epsout; epsz(1:xmymzm+xm*ym) = epsout

   ! mark volume within vdw srf as 2, outside -2, so there are only contact-type
   ! boundary grid points

   rh = 1/h
   insas(1:xmymzm_ext) = -2
   zv(1:xmymzm_ext) = 9999.0d0
   iatm = -1
   range1 = radi(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   iatm = 1
   call exvwsph_cap( 2, insas, atmsas, zv(1) )

   ! finally, save boundary edges for db energy and forces105 
   call epsbnd_cap( atmsas, insas )

   ! use the insas grid to setup epsx, epsy and epsz maps

   call epsmap_cap( insas, atmsas, epsx, epsy, epsz )


contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph_cap( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none
    
   ! Passed variables
    
   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2
    
   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
               end if
               insph(i,j,k) = dielsph;
               inatm(i,j,k) = iatm; dst(i,j,k) = d2
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exvwsph_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd_cap ( atmsas,insas )
    
   implicit none
    
   integer atmsas(0:xm+1,0:ym+1,0:zm+1)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   logical boundary
   integer buffer, i, j, k, clstmp
    
   nwarn = 0
   nbnd = 0
   buffer = 1
   do k = buffer, zm+1-buffer; do j = buffer, ym+1-buffer; do i = buffer, xm+1-buffer
       
      ! set up condition for a boundary grid point
       
      boundary = .false.
      if ( (insas(i,j,k) ==  1 .or. insas(i,j,k) ==  2) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) == -1 .or. insas(i,j,k) == -2) .and.&
           (insas(i-1,j,k) ==  1 .or. insas(i-1,j,k) ==  2 .or. insas(i+1,j,k) ==  1 .or.&
            insas(i+1,j,k) ==  2 .or. insas(i,j-1,k) ==  1 .or. insas(i,j-1,k) ==  2 .or.&
            insas(i,j+1,k) ==  1 .or. insas(i,j+1,k) ==  2 .or. insas(i,j,k-1) ==  1 .or.&
            insas(i,j,k-1) ==  2 .or. insas(i,j,k+1) ==  1 .or. insas(i,j,k+1) ==  2) ) then
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
 
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k
      if ( ifcap /= 0 ) then
         iepsav(4,nbnd) = -1
         cycle
      end if
 
      ! for a grid point in contact region +/- 2 or in a  solvent probe, simply use the atom/probe that
      ! marks it
 
      clstmp = 0
      if ( abs(insas(i,j,k)) == 2 ) then
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
      elseif ( insas(i,j,k) == -1 ) then
         clstmp = -atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
 
      ! for a buried reentry grid point, find the atom that marked its neighoring exposed reentry
      ! grid points. Note that this may not be possible when grid spacing is large
 
      else if ( insas(i,j,k) == 1 ) then
         clstmp = -fndcls_cap( i, j, k, insas, atmsas )
         if ( clstmp == 0 ) then
            nwarn = nwarn + 1
         end if
      end if
 
      iepsav(4,nbnd) = clstmp
   end do; end do; end do
   if ( nwarn > 0 ) then
      if ( pbverbose .and. level == nfocus ) write(6, '(a,i4)') &
      'PB Warning in epsbnd(): No neighbor found for boundary grids total:', nwarn
   end if
 
end subroutine epsbnd_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls_cap( i,j,k,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls_cap, i, j, k
   integer  insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls_cap = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - dprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls_cap = clsatm(ip)
   end do
 
end function fndcls_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap_cap( insas,atmsas,epsx,epsy,epsz )

   implicit none
   integer insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   _REAL_ epsint

   epsint = 2.0d0*epsin*epsout/(epsin+epsout)
 
   do k = 0, zm; do j = 0, ym; do i = 0, xm
      a = insas(i,j,k)
      b = insas(i+1,j,k)
      a1 = atmsas(i,j,k)
      b1 = atmsas(i+1,j,k)
      if ( j == 0 .or. k == 0 ) then
         ! do nothing
      else if ( sign(a,b) == a ) then
         if ( a > 0 ) then
            epsx(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracx_cap(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout)
         epsx(i,j,k) = epsint
      end if
      c = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( i == 0 .or. k == 0 ) then
         ! do nothing
      else if ( sign(a,c) == a ) then
         if ( a > 0 ) then
            epsy(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracy_cap(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout)
         epsy(i,j,k) = epsint
      end if
      d = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( i == 0 .or. j == 0 ) then
         ! do nothing
      else if ( sign(a,d) == a ) then
         if ( a > 0 ) then
            epsz(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracz_cap(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout)
         epsz(i,j,k) = epsint
      end if
   end do; end do; end do
 
!   do k = 1, zm
!      write(20, *) 'plane', k
!   do j = 1, ym
!      write(20, '(100f6.1)') epsx(1:xm,j,k)/eps0
!   end do
!   end do
!   do k = 1, zm
!      write(21, *) 'plane', k
!   do i = 1, xm
!      write(21, '(100f6.1)') epsy(i,1:ym,k)/eps0
!   end do
!   end do
!   do j = 1, ym
!      write(22, *) 'plane', j
!   do i = 1, xm
!      write(22, '(100f6.1)') epsz(i,j,1:zm)/eps0
!   end do
!   end do

end subroutine epsmap_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this x-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracx_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - xi + dble(i+1)
      else
         aa = range3 + xi - dble(i)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracx_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this y-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracy_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - yi + dble(j+1)
      else
         aa = range3 + yi - dble(j)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracy_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this z-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracz_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - zi + dble(k+1)
      else
         aa = range3 + zi - dble(k)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracz_cap

end subroutine pb_exmol_cap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calc_sa1(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt)

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   implicit none

#  include "pb_constants.h"
#  include "md.h"

   ! Passed variables
   _REAL_ acrd(3,*)
   integer xm,ym,zm,xmymzm,nbndx,nbndy,nbndz
   !integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   integer iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer sasopt

   ! Local variables

   integer i, j, k, iatm, ip
   integer cnt_dot, rnt_dot, dim_dot, tri_dot
   _REAL_ x(3), crd(3)
   _REAL_ rn(1:3), rsphere, dr, r1, r2, r3, h2, hh
   _REAL_ ds1, total_s1, ds2, total_s2
   _REAL_ ds, total_s , cnt_s, rnt_s, dim_s, tri_s
   _REAL_ dss, tss, rx, ry, rz, ess, e0, e1

   ! mjhsieh: warning eliminator
   rsphere = -1d0
   h2 = h*h
   hh = HALF*h

   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO

   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
!     crd(1) = gox + h*i + fedgex(ip)*h; crd(2) = goy + h*j; crd(3) = goz + h*k
      crd(1) = gox + h*i + hh; crd(2) = goy + h*j; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         write(6,*) 'PBMD FATAL ERROR: cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(1))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     rx = ONE/rn(1)
!     dss = atan(rx*(rn(2)+hh)*(rn(3)+hh)/sqrt(rn(1)**2+(rn(2)+hh)**2+(rn(3)+hh)**2)) &
!         - atan(rx*(rn(2)+hh)*(rn(3)-hh)/sqrt(rn(1)**2+(rn(2)+hh)**2+(rn(3)-hh)**2)) &
!         - atan(rx*(rn(2)-hh)*(rn(3)+hh)/sqrt(rn(1)**2+(rn(2)-hh)**2+(rn(3)+hh)**2)) &
!         + atan(rx*(rn(2)-hh)*(rn(3)-hh)/sqrt(rn(1)**2+(rn(2)-hh)**2+(rn(3)-hh)**2)) 
!     dss = dss*rx*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndx

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
!     crd(1) = gox + h*i; crd(2) = goy + h*j + fedgey(ip)*h; crd(3) = goz + h*k
      crd(1) = gox + h*i; crd(2) = goy + h*j + hh; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         write(6,*) 'PBMD FATAL ERROR: cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(2))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3 
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(2)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     ry = ONE/rn(2)
!     dss = atan(ry*(rn(1)+hh)*(rn(3)+hh)/sqrt(rn(2)**2+(rn(1)+hh)**2+(rn(3)+hh)**2)) &
!         - atan(ry*(rn(1)+hh)*(rn(3)-hh)/sqrt(rn(2)**2+(rn(1)+hh)**2+(rn(3)-hh)**2)) &
!         - atan(ry*(rn(1)-hh)*(rn(3)+hh)/sqrt(rn(2)**2+(rn(1)-hh)**2+(rn(3)+hh)**2)) &
!         + atan(ry*(rn(1)-hh)*(rn(3)-hh)/sqrt(rn(2)**2+(rn(1)-hh)**2+(rn(3)-hh)**2)) 
!     dss = dss*ry*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndy

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
!     crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+fedgez(ip)*h
      crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+hh

      if ( iatm == 0 ) then
         write(6,*) 'PBMD FATAL ERROR: cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(3))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3 
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(3)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     rz = ONE/rn(3)
!     dss = atan(rz*(rn(1)+hh)*(rn(2)+hh)/sqrt(rn(3)**2+(rn(1)+hh)**2+(rn(2)+hh)**2)) &
!         - atan(rz*(rn(1)+hh)*(rn(2)-hh)/sqrt(rn(3)**2+(rn(1)+hh)**2+(rn(2)-hh)**2)) &
!         - atan(rz*(rn(1)-hh)*(rn(2)+hh)/sqrt(rn(3)**2+(rn(1)-hh)**2+(rn(2)+hh)**2)) &
!         + atan(rz*(rn(1)-hh)*(rn(2)-hh)/sqrt(rn(3)**2+(rn(1)-hh)**2+(rn(2)-hh)**2)) 
!     dss = dss*rz*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndz

   total_s = total_s*h2
!  total_s1 = total_s1*h2
!  total_s2 = total_s2*h2
   cnt_s = cnt_s*h2
   rnt_s = rnt_s*h2
   tri_s = tri_s*h2
   dim_s = dim_s*h2
!  e1 = e1/(nbndx+nbndy+nbndz)
!  if ( saopt < 0 .or. imin == 6 ) then
      write(6,'(a,f12.4)') ' Total molecular surface',total_s
!     write(6,'(a,f20.10)') 'Total molecular surface',total_s
!     write(6,'(a,f20.10)') 'Total contact surface',cnt_s
!     write(6,'(a,f20.10)') 'Total reentry surface',rnt_s
!     write(6,'(a,f20.10)') 'Total dimer surface',dim_s
!     write(6,'(a,f20.10)') 'Total trimer surface',tri_s
!     write(6,'(a,i10)') 'contact boundary point',cnt_dot
!     write(6,'(a,i10)') 'reentrant boundary point',rnt_dot
!     write(6,'(a,i10)') 'dimer boundary point',dim_dot
!     write(6,'(a,i10)') 'trimer boundary point',tri_dot

!     write(6,'(a,f20.10)') 'zero order',total_s1
!     write(6,'(a,f20.10)') 'second order',total_s2
!     write(6,'(a,f20.10)') 'all order',tss
!     write(6,'(a,f20.10)') 'mode oo error',e0
!     write(6,'(a,f20.10)') 'mode 1 error',e1
!     print *,total_s,cnt_s,rnt_s
!  end if
!  print *, total_s

end subroutine calc_sa1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calc_sa2(acrd,xm,ym,zm,xmymzm,nbnd,iepsav,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,sasopt,smoothopt,&
                    epsin,epsout,insas,epsx,epsy,epsz)

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   implicit none

#  include "pb_constants.h"
#  include "md.h"

   ! Passed variables

   _REAL_ acrd(3,*)
   integer xm,ym,zm,xmymzm,nbnd
   integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   integer sasopt,smoothopt
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsin, epsout
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)

   ! Local variables

   integer ip, i, j, k, iatm
   _REAL_ g(3), x(3), dx(3), hh(3), repsp(3), repsm(3)
   _REAL_ dist, dx2, ds, total_s, sgn, half_h, epsth, repsin, repsout
   _REAL_ cnt_s, rnt_s, term, dim_s, tri_s

   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   repsin = ONE/epsin; repsout = ONE/epsout
   epsth = TWO/(repsin+repsout)
!  half_h = HALF*h

   do ip = 1, nbnd

      ! collecting boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the surface grid point on to the molecular surface, crd() is the
      ! new coord, and x() is the atom/probe coord, fx/y/z0 is the grid version
      ! of crd()

      if      ( iatm == 0 ) then
         write(6,*) 'PBMD FATAL ERROR: can not find owner of boundary grid points'
         call mexit(6, 1)
      else if ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            dist = radi(iatm)+dprob
         else
            dist = radi(iatm)
         end if
         sgn = ONE
      else
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         dist = dprob
         sgn = -ONE
      end if

      dx = g - x; dx2 = dx(1)**2 + dx(2)**2 + dx(3)**2
!     hh = half_h

      dist = sqrt(dx2)
!     ds = (dx(1)+hh(1))/epsx(i,j,k)-(dx(1)-hh(1))/epsx(i-1,j,k) + &
!          (dx(2)+hh(2))/epsy(i,j,k)-(dx(2)-hh(2))/epsy(i,j-1,k) + &
!          (dx(3)+hh(3))/epsz(i,j,k)-(dx(3)-hh(3))/epsz(i,j,k-1)  
!     ds = dx(1)/epsx(i,j,k)-dx(1)/epsx(i-1,j,k) + &
!          dx(2)/epsy(i,j,k)-dx(2)/epsy(i,j-1,k) + &
!          dx(3)/epsz(i,j,k)-dx(3)/epsz(i,j,k-1) 
      if ( smoothopt == 1 .or. smoothopt == 2 ) then
         repsp = repsin
         repsm = repsin
         if ( epsx(i,j,k) > epsth ) repsp(1) = repsout
         if ( epsy(i,j,k) > epsth ) repsp(2) = repsout
         if ( epsz(i,j,k) > epsth ) repsp(3) = repsout
         if ( epsx(i-1,j,k) > epsth ) repsm(1) = repsout
         if ( epsy(i,j-1,k) > epsth ) repsm(2) = repsout
         if ( epsz(i,j,k-1) > epsth ) repsm(3) = repsout
         ds = dx(1)*repsp(1)-dx(1)*repsm(1) + &
              dx(2)*repsp(2)-dx(2)*repsm(2) + &
              dx(3)*repsp(3)-dx(3)*repsm(3) 
      else
         ds = dx(1)/epsx(i,j,k)-dx(1)/epsx(i-1,j,k) + &
              dx(2)/epsy(i,j,k)-dx(2)/epsy(i,j-1,k) + &
              dx(3)/epsz(i,j,k)-dx(3)/epsz(i,j,k-1) 
      end if
      term = sgn*ds/dist
      total_s = total_s + term
      if ( iatm > 0  ) then
         cnt_s = cnt_s + term
      else
         rnt_s = rnt_s + term
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + term
         else
            dim_s = dim_s + term
         end if
      end if
      
   end do

   term = h*h/(repsout-repsin)
   total_s = total_s*term
   cnt_s = cnt_s*term
   rnt_s = rnt_s*term
   dim_s = dim_s*term
   tri_s = tri_s*term
!  if ( saopt < 0 .or. imin == 6 ) then
!     print *, nbnd, nbndx+nbndy+nbndz
      write(6,'(1x,a,f12.4)') 'Total molecular surface',total_s
!     write(6,'(1x,a,f12.4)') 'Total contact surface',cnt_s
!     write(6,'(1x,a,f12.4)') 'Total reentry surface',rnt_s
!     write(6,'(1x,a,f12.4)') 'Total dimer surface',dim_s
!     write(6,'(1x,a,f12.4)') 'Total trimer surface',tri_s
!     print *,total_s,cnt_s,rnt_s
!  end if
!  print *, total_s

end subroutine calc_sa2
