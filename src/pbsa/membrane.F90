#define _REAL_ double precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Membrane Levelset Setup Functions !!!!!!!!!!!!!!!!!!!!!!!!!!

!membrane_density_calc
!calculates the levelset value contribution due to the membrane
!   location specified by i,j,k
!   poretype specifies an exclusion region
!   handled by multiplying the membrane density value
!   by an appropriate factor if in or near the excluded zone
!   poretype meanings:
!       value       meaning
!       -----       ------
!       0           no pore
!       1           cylindrical pore
!       2-4         reserved for future implementations
!   This file should be made into a module later. For now
!   poredata will have to be calculated before hand, stored
!   externally, and passed in.
!   see exclusion_calculation for more info.

!This routine is being modularized to (hopefully) make adding
!new membrane definitions easier. 
!Membrane density at a given point i,j,k will be computed by a call to
!routine integrate_point_density.
!This subroutine will calculate the membrane density for the point at i,j,k
! by integrating the membrane density in a ball of radius
! mprobe+2*dprobe centered on i,j,k.
! This is a bit klunky for now since it will simply do this point wise
! by itterating over all grid nodes within the ball.
! Each point will be checked using routine is_in_membrane to see if it
! should contribute to the density at i,j,k. For each node that does
! contribute, routine pointwise_density_calc(xi,yi,zi,xj,yj,zj,u) is called
! to compute the density contribution.

subroutine membrane_density_calc(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                        mthick,mprobe,dprobe,poretype,&
                        poredata,dval,mctrdz)
    implicit none
    !passed variables
    integer     i,j,k       !node location for calculation
    integer     xm,ym,zm    !grid array index dimensions
    integer     poretype
    _REAL_      gox,goy,goz !box origin
    _REAL_      h           !grid spacing
    _REAL_      mprobe      !membrane probe radius (just a dummy for now)
    _REAL_      dprobe      !water probe radius
    _REAL_      mthick      !membrane thickness
    _REAL_      poredata(*) !dimensions will depend on poretype
                            !will be reshaped later
    _REAL_      dval        !output
    _REAL_      mctrdz


    !local variables
    _REAL_      u
    _REAL_      mdist       !distance to membrane z-center at k
    _REAL_      mzctr       !z grid coordinate of membrane center 
    !_REAL_      eval        !distance to surface of exclusion region
    _REAL_      dcutoff     !cutoff for the magnitude of the membrane
                            !density near the membrane center.

    mzctr = goz + (1+zm)/2.0d0*h + mctrdz
    dval = 0;

    !we need to set cutoff to the minimum between the membrane probe
    !radius and the membrane thicknes. This will ensure the magnitude
    !of the membrane levelset will
    !1) Not have a cusp at the center (ensured by smoothing function)
    !2) Will not be larger than mprobe
    dcutoff = MIN(mprobe,mthick)

    If (mprobe > mthick) then
        write(6,*) 'PB Warning: membrane thickness is less than',&
                    ' membrane probe radius'
    end if

    If (mthick > zm*h) then
        write(6,*) 'PB Bomb in membrane_density_calc(): membrane is thicker than',&
                   'computation grid'
        call mexit(6,1)
    end if

    !mdist is signed distance to slab region boundary

    !write(6,*) '-debug: computing membrane slab levelset'
    mdist = abs(k*h+goz-mzctr)-mthick/2.0d0

    dval = mdist

    IF (poretype == 0) then
    !Slab like membrane, only need to rescale distance to match
    !membrane profile. Distance profile set such that interior
    !density will taper off to a value of mprobe for large membranes
    !this should remove any cusps from the center of the membrane.

        call distance_to_density(dval,u,dcutoff,dprobe)
        dval = u
            
    else IF (poretype == 1) then
    !Membrane has a cylinderical pore region. We use
    !routine get_pore_surf_dist to find the distance from i,j,k to
    !the surface of this region. If this distance is 

       !write(6,*) '-debug: adding cylinder exclusion'; flush(6)
       
       call add_cylinder_exclusion(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                                   poredata,dcutoff,dval)

       !write(6,*) '-debug: computing density from distance';flush(6)
       call distance_to_density(dval,u,dcutoff,dprobe)
       dval = u
        
    else
        write(6,*) 'PB Bomb in membrane_density_calc(): unknown poretype'
        call mexit(6,1)
    end if

    !dval = 1-dval
    !write(6,*) '-debug: done'    

end subroutine membrane_density_calc

!This routine modifies the slab-like membrane levelset to include
!a cylinderical exclusion region. The poredata array is to contain the
!radius of this region at z-constant plane of the grid. This is to
!potentially allow more complicated surfaces, such as slanted solids
!of revolution in later versions but modification of the cylinder
!surface distance calculation will be needed if they are used.
!For now, only orthoganal uniform radius cylinders will work correctly
!see gen_cylinder_data for generation of poredata array
subroutine add_cylinder_exclusion(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                poredata,mprobe,ddist)

    implicit none
    !Passed Data
    integer xm,ym,zm                !grid nodes per dim
    integer i,j,k                   !grid location indices for calculation
    _REAL_  gox,goy,goz,h           !grid origin and spacing
    _REAL_  poredata(0:zm+1,3)      !data for cylinderical pore
    _REAL_  mprobe                  !membrane probe radius
    _REAL_  ddist                   !output

    !Local Variables
    _REAL_  edist                   !distance to cylinder surface
    _REAL_  mdist                   !scratch variable
    _REAL_  x,y,z                   !grid cooridinates from indicies
    integer zi

    x = gox+i*h; y = goy + j*h; z = goz + k*h
    mdist = ddist / 2.0 / mprobe

    !mprobe is intended to allow definition of a membrane probe radius
    !that is different from that of the solvent...
    !mprobe is not actually doing anything here as of yet but is
    !included for potential use in future implementations.

    !write(6,*) '-debug: poredata = ';flush(6)
!        do zi=0,zm+1
!           write(6,*) zi,poredata(zi,1:3)
!        end do

    !write(6,*) '-debug: computing distance to cylinder surface';flush(6)
    call cylinder_surf_distance(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                                poredata,mprobe,edist)

    !write(6,*) '-debug: adjusting level set'; flush(6)
    edist = edist

    If ( mdist < 0 ) THEN !Interior to membrane bounds
        IF ( edist > 0 ) THEN 
            !Inside the membrane but not in the exclusion region. Here,
            !'distance' is calculated as geometric mean of membrane and
            !cylinder distance when cylinder surface is closer than
            !membrane surface. If membrane is closest cylinder is
            !defacto ignored
            ddist = -Sqrt(abs(mdist)*min(abs(mdist),abs(edist)))
        ELSE
            !We are inside the exclusion region so we use the distance
            !to the cylinder surface
            ddist = -edist
        END IF
    ELSE !Exterior to membrane bounds
        !we assume cylinderical surface ends at membrane bounds 
        !and does not have faces, so we will need to compute the distance
        !to the edge if we are directly above the cylinder 
        IF (edist <= 0) THEN
        !We are directly above the cylinder so we will need to compute the
        !distance to the edge of the cylinder. This is the root square
        !sum of the distance to the extended cylinder edgen and the distance
        !to the membrane surface... we assume that the pore data for regions
        !above/below the membrane will contain the radius of the cylinder at
        !the corresponding membrane edge... if not there may be problems
        !so modify at your own risk
            ddist = sqrt(mdist**2.0d0 + edist**2.0d0)
        ELSE
        !We are closer to membrane surface here so we use its distance directly
            ddist = mdist
        END IF
    END IF

    ddist = ddist*2.0d0*mprobe

end subroutine add_cylinder_exclusion

!Calculates the distance for node i,j,k to the surface of the cylinderical
!exclusion region. poredata(k,1) and poredata(k,2) are x/y location of center
!of the cylinderical region at k. poredata(k,3) is the radius of the region
!at k. Although only orthoganal cylinders are supported presently, this
!should allow future implementations to construct more complex regions
!such as slanted cylinders. The distance computation will need to be modified,
!here and in the above add_cylinder_exclusion routines to do so properly.
subroutine cylinder_surf_distance(i,j,k,xm,ym,zm,gox,goy,goz,h,&
                                    poredata,mprobe,dist)

    implicit none
    !Passed Variables
    integer xm,ym,zm
    integer i,j,k
    _REAL_  gox,goy,goz,h
    _REAL_  poredata(0:zm+1,3)
    _REAL_  mprobe
    _REAL_  dist                !output

    !Internal Variables
    _REAL_  x,y,z
    _REAL_  dx,dy,dr
    integer zi

    !write(6,*) '-debug: writting poredata';flush(6)
!    do zi=0,zm+1
!        write(6,*) zi,poredata(zi,1:3);flush(6)
!    end do

    !write(6,*) '-debug: computing i,j,k->x,y,z';flush(6)

    x = i*h+gox;y=j*h+goy;z=k*h+goz

    !write(6,*) '-debug: computing distance to pore center';flush(6)

    dx = Abs(x-poredata(k,1))
    dy = ABS(y-poredata(k,2))
    dr = sqrt(dx**2.0d0 + dy**2.0d0)

    !write(6,*) '-debug: getting distance to cylinder edge';flush(6)

    dist = dr - poredata(k,3)
    dist = dist / 2.0d0 / mprobe

    !write(6,*) '-debug: done';flush(6)
    
end subroutine cylinder_surf_distance


!Rescales the given distance value to density value for the membrane interior
!Units of distance here should be in # of solvent probe radii, not angstroms
!rprobe should be membrane probe radius in terms of # of solvent probe radii
!the distance will be transformed from a linear function to a cubic polynomial
!that smoothly transitions from 0 to plateaue at a minimum of -rprobe 
!over a range of rprobe/2 units. The initial slope will match that of the spline
!function for the external density.
!!!Alternative option may be to hijack the exterior spline function and apply
!!!it symmetrically about zero instead... this may make matching derivatives
!!!tricky when the membrane probe radius differs from the solvent probe radius.
subroutine interior_density(dist,u,rprobe)

    implicit none
    
    !Passed Variables
    _REAL_ dist
    _REAL_ u
    _REAL_ rprobe

    !Local Variables
    _REAL_ sdata(4)
    _REAL_ temp

    data sdata /0,4.527143,-6.108572,2.108572/

    save sdata    

    u = dist

    temp = -1.0d0*rprobe

    IF (u < temp/2.0d0) THEN
        u = temp
    Else
        u = sdata(1)                            + &
            sdata(2)*u                          + &
            sdata(3)*(u**2.0d0)/(temp)          + &
            sdata(4)*(u**3.0d0)/(temp**2.0d0)
    end if

    !u = 1.0 - u

end subroutine interior_density

!Similar to the density calc in pb_exmol. The input value
!ddist should be grid unit distance divided by 2*dprobe 
!will return output value u using the same spline based fitting
!as in density_calc, however, unlike pb_exmol, u WILL NOT
!be set to u=1.00-u at the end. This is to allow it to be
!more easily added to existing u value.
subroutine distance_to_density(dist,u,mprobe,dprobe)

    implicit none

    !Passed Variables
    _REAL_  dist    !psuedo signed distance
    _REAL_  mprobe  !membrane probe radius
    _REAL_  dprobe  !solvent probe radius
    _REAL_  u       !output

    !Internal Variables
    _REAL_  mspcoef(4,4), mdash(5)

    integer m       !counter

   
   !Duplicate of spline data for density_calc in pb_exmol 
   data mspcoef /1.000000,0.1000000,6.4999998E-02,2.9999999E-02, &
                -4.527143,-1.745714,0.2900000,-0.2542857,       &
                0.0000000,11.12571,-2.982857,0.8057143,         &
                14.83429,-18.81143,5.051429,-1.074286 /
   data mdash /0.d0,0.25d0,0.5d0,0.75d0,1.d0 / 

   save mspcoef,mdash

    u=0.0d0

    dist = dist / 2.0 / dprobe !rescale from angstroms to # of probe radii

    if (dist <= 1.0d0) then
        if ( dist <= 0.0d0 ) then
        !mprobe/dprobe fed to probe radius to ensure that the minimum
        !attained is the membrane prode radius (in terms # of solvent probe
        !radii)
            !write(6,*) '-debug: setting interior density';flush(6)
            call interior_density(dist,u,1.0d0)
            u =  1.0d0-u
            !u = 1.0d0
        else
        !    call interior_density(-dist,u,1.0d0)
        !    u = -u

            !write(6,*) '-debug: setting exterior density';flush(6)
            do m=1,4
                if ( dist > mdash(m) .and. dist <= mdash(m+1) ) then
                    u = u                               +&
                        mspcoef(m,1)                    +&
                        mspcoef(m,2)*(dist-mdash(m))    +&
                        mspcoef(m,3)*(dist-mdash(m))**2 +&
                        mspcoef(m,4)*(dist-mdash(m))**3
                end if
            end do
        end if

    end if

end subroutine distance_to_density

!gen_cylinder_data
!Generates the poredata needed to describe a right cylinder exclusion
!region, centered on centeroid of the solute inside the membrane region.
!be sure that the poredata array has been allocated as a 1D array
!-> use allocate( poredata(3*(zm+2)) )
subroutine gen_cylinder_data(xm,ym,zm,gox,goy,goz,h,mthick,rad,&
                    natm,poredata)

    use poisson_boltzmann, only : acrd, mctrdz

    implicit none
    !passed variables
    integer     xm,ym,zm
    integer     natm
    _REAL_      gox,goy,goz,h
    _REAL_      mthick,rad
    _REAL_      poredata(0:zm+1,3) !output

    !local variables
    integer     zi,iatm,acount
    _REAL_      centeroid(2)
    _REAL_      mzctr

    
    mzctr = goz+(1+zm)/2.0d0*h+mctrdz   

    acount = 0
    write(6,*) "building cylindrical pore exclusion region data"
    write(6,*) "radius = ",rad
    !find centeroid of membrane bound solute
    do iatm = 1,natm
        if ( abs(mzctr-acrd(3,iatm)) <= mthick/2.0d0 ) then
            centeroid(1) = centeroid(1) + acrd(1,iatm)
            centeroid(2) = centeroid(2) + acrd(2,iatm)
            acount = acount + 1
        end if  
    end do
    centeroid = centeroid / acount
    write(6,*) "center = ",centeroid; flush(6)
    !build cylinder data
    !!write(6,*) '-debug: writting poredata';flush(6)
    do zi = 0,zm+1
        poredata(zi,1) = centeroid(1)
        poredata(zi,2) = centeroid(2)
        poredata(zi,3) = rad
       ! write(6,*) zi,poredata(zi,1:3);flush(6)
       ! write(6,*) 'zi,poredata(zi,1-3)',zi,poredata(zi,1:3); flush(6)
    end do


end subroutine gen_cylinder_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! MEMBRANE EPS ARRAY SETUP ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_membrane_eps(epsx,epsy,epsz,epsin,epsout,epsmemb,&
               lvlset,insas,gox,goy,goz,xm,ym,zm,h,mthick,mctrdz,&
               ipb)
    implicit none
    !passed variables
    integer ipb
    integer xm,ym,zm
    _REAL_ gox,goy,goz,h,mthick,mctrdz
    _REAL_ epsin,epsout,epsmemb
        !the eps grids apparently need to store on a single
        !additional outer edge, so take care with indices
        !if you acces their unshaped 1D versions
    _REAL_ epsx(0:xm,1:ym,1:zm)
    _REAL_ epsy(1:xm,0:ym,1:zm)
    _REAL_ epsz(1:xm,1:ym,0:zm)
    _REAL_ lvlset(0:xm+1,0:ym+1,0:zm+1)
    integer insas(0:xm+1,0:ym+1,0:zm+1)
    
    !local variables
    integer xi,yi,zi
    _REAL_ epstemp,epsuv0,epsum0,epsvm0
    _REAL_ mgzctr,mgzmax,mgzmin
    integer mizmax,mizmin
    integer p0,px,py,pz

    !since grids outside the membrane bounds are essentially unchanged,
    !we can save some time by only itterating over relevant z-indices
    !to do so, we will need to again locate the membrane bounds.
    mgzctr = goz + (1+zm)/2.0d0*h + mctrdz
    mgzmax = mgzctr + mthick/2.0d0
    mgzmin = mgzctr - mthick/2.0d0
    mizmax = Ceiling(1+(mgzmax-goz)/h)
    mizmin = Floor(1+(mgzmin-goz)/h)
    
    if (mizmax > zm .or. mizmin < 0) then
      !sanity check on the membrane bounds. Shouldn't be a 
      !problem, but here just in case the earlier membrane
      !setup failed to catch it.
      write(6,*) 'PB Bomb in set_membrane_eps(): Membrane boundaries exceed grid boundaries'
      call mexit(6,1)
    end if

    !!!write(6,*) '-debug: writting computed membrane topology data:'
    write(6,*) CHAR(09),'gox,goy,goz,h: ',sngl(gox),sngl(goy),sngl(goz),sngl(h)
    write(6,*) CHAR(09),'mgzctr: ',sngl(mgzctr)
    write(6,*) CHAR(09),'mgzmin: ',sngl(mgzmin)
    write(6,*) CHAR(09),'mgzmax: ',sngl(mgzmax)
    write(6,*) CHAR(09),'mizmin: ',mizmin
    write(6,*) CHAR(09),'mizmax: ',mizmax
    flush(6)

    !Since we now have 3 seperate interfacial types to consider
    ! *solute-solvent* (same as before)
    ! *solute-membrane* 
    ! *solvent-membrane*
    !we will need to hold more than just 1 reduced mass eps value
    !so instead of just epsint0 we now have
    ! *epsuv0 : solute-solvent - same as epsint0 in pb_exmol
    ! *epsum0 : solute-membrane
    ! *epsvm0 : solvent-membrane
    
    !epsuv0 is probably unnecessary since the only thing we need to
    !change is grids in or adjacent to the membrane. It is included
    !for completeness
    epsuv0 = 2.0d0*epsin*epsout/(epsin+epsout)
    epsum0 = 2.0d0*epsmemb*epsin/(epsin+epsmemb)
    epsvm0 = 2.0d0*epsmemb*epsout/(epsmemb+epsout)
    !in pb_exmol, epstint (epstemp here) would be initizlized now...
    !unfortunately we cant do that yet since the value we want will
    !depend on the locations of the grid nodes we choose.
    

    !this is the main work loop
    !eps arrays are set as if smoothopt=1 and sasopt=2 are being used
    !may add in alternatives later. Since I dont like to reinvent the wheel,
    !we make calls to epsmembfracx_r, epsmembfracy_r and epsmembfracz_r
    !which just modified copies of the equivalents in pb_exmol.
    !Direct calls were infeasible since the level set does not behave
    !properly across the membrane-solute interface.
    !!!write(6,*) '-debug: starting outer eps array z loop';flush(6)        

    do zi = mizmin,mizmax !only need to loop over z values near membrane
        !!!write(6,*) '-debug: setting eps arrays for zi = ',zi,' zg = ',goz+zi*h;flush(6)
    do xi = 0,xm; do yi = 0,ym
        p0 = insas(xi,yi,zi)
        px = insas(xi+1,yi,zi)
        py = insas(xi,yi+1,zi)
        pz = insas(xi,yi,zi+1)
        if ( ((p0*px*py*pz) ==0) .and. ((abs(p0)+abs(px)+abs(py)+abs(pz)) /= 0 ) ) then
                !write(6,*) CHAR(09),'-debug: p0,px,py,pz: ',p0,px,py,pz;flush(6)
        end if
        !since most of eps arrays have been set already, we only
        !need to focus on grids within or adjacent to the membrane.
        !The insas array is 0 in the membrane, so we can ignore any pair
        !of adjacent nodes whose insas product is nonzero.
        !We can also exploit insas to tell which dielectric boundary has 
        !been crossed by looking at the sum.

        !set epsx
        if ( yi>0 .and. zi>0 ) then !epsx undefined for yi or zi<=0
        if (p0*px == 0) then
           if (p0+px == 0) then !membrane regular point
                epsx(xi,yi,zi) = epsmemb
            else if (p0+px > 0) then !membrane-solute boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsx at',xi,yi,zi;flush(6) 
                epstemp = epsum0
                !write(6,*) CHAR(09),CHAR(09),'-debug: interpolating epsx value';flush(6)
                call epsmembfracx_r( xi, yi, zi, p0, px, 0, 0, 1/h, epstemp,&
                                     epsin, epsout, lvlset, insas,xm,ym,zm,ipb )
                !write(6,*) CHAR(09),CHAR(09),'-debug: setting epsx';flush(6)
                epsx(xi,yi,zi) = epstemp
            else if (p0+px < 0) then !membrane-solvent boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsx at',xi,yi,zi;flush(6)  
                epstemp = epsvm0
                 call epsmembfracx_r( xi, yi, zi, p0, px, 0, 0, 1/h, epstemp,&
                                     epsin, epsout, lvlset, insas,xm,ym,zm,ipb )
                 epsx(xi,yi,zi) = epstemp
            else !this should be mathematically impossible
            !but just in case it does occur,
            ! we cowardly refuse to go any further
                  write(6,*) 'PB Bomb in set_membrane_eps(): insas array memory is corrputed'
                  call mexit(6,1)
            end if
        end if
        end if
        !set epsy
        if ( xi > 0 .and. zi>0 ) then
            if (p0*py == 0) then
               if (p0+py == 0) then !membrane regular point
                  epsy(xi,yi,zi) = epsmemb
               else if (p0+py > 0) then !membrane-solute boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsy at',xi,yi,zi;flush(6)
                  epstemp = epsum0
                  call epsmembfracy_r(xi,yi,zi,p0,py,0,0,1/h,epstemp,&
                                       epsin,epsout,lvlset,insas,xm,ym,zm,ipb)
                  epsy(xi,yi,zi) = epstemp
               else if (p0+py < 0) then !membrane-solvent boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsy at',xi,yi,zi;flush(6)
                  epstemp = epsvm0 
                  call epsmembfracy_r(xi,yi,zi,p0,py,0,0,1/h,epstemp,&
                                       epsin,epsout,lvlset,insas,xm,ym,zm,ipb)
                  epsy(xi,yi,zi) = epstemp
               else
                  ! This shouldn't be possibe, but just in case it does
                  ! we will bravely run away
                  write(6,*) 'PB Bomb in set_membrane_eps(): insas array memory is corrputed'
                  call mexit(6,1)
               end if
            end if
        end if
        !set epsz
        if ( xi > 0 .and. yi>0) then
            if (p0*pz == 0) then
                if (p0+pz == 0) then !membrane regular point
                    epsz(xi,yi,zi) = epsmemb
                else if (p0+pz > 0) then !membrane-solute boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsz at',xi,yi,zi;flush(6)
                    epstemp = epsmemb
                    call epsmembfracz_r(xi,yi,zi,p0,pz,0,0,1/h,epstemp,&
                                        epsin,epsout,lvlset,insas,xm,ym,zm,ipb)
                    epsz(xi,yi,zi) = epstemp
                else if (p0+pz<0) then !membrane-solvent boundary
            !write(6,*) CHAR(09),'-debug: setting membrane epsz at',xi,yi,zi;flush(6)
                    epstemp=epsmemb
                    call epsmembfracz_r(xi,yi,zi,p0,pz,0,0,1/h,epstemp,&
                                        epsin,epsout,lvlset,insas,xm,ym,zm,ipb)
                    epsz(xi,yi,zi) = epstemp
                else !this shouldnt be possible. But just in case
                !We bail with a shout of memory corruption
                  write(6,*) 'PB Bomb in set_membrane_eps(): insas array memory is corrputed'
                  call mexit(6,1)
                end if
            end if
        end if
        !Eps array sanity checks. We assume that solvent has the highest eps value,
        !and that either solute or membrane have the lowest eps value.
        !If any of the eps array's have a value that falls outside of this range
        !we will complain angerly and bail out.
        if ( ((p0*px*py*pz) ==0) .and. ((abs(p0)+abs(px)+abs(py)+abs(pz)) /= 0 ) ) then
 
            if ((yi/=0 .and. zi/=0).and.(p0*px == 0)) then
            if (epsx(xi,yi,zi) > epsout .or. &
                (epsx(xi,yi,zi) < epsin .and. epsx(xi,yi,zi) < epsmemb)) then
                write(6,*) 'PB Bomb in set_membrane_eps: epsx array value out of bounds'
                write(6,*) 'value: ',sngl(epsx(xi,yi,zi)),' anomalous at ',xi,yi,zi
                call mexit(6,1)
            end if
            end if
            if ((xi/=0 .and. zi/=0).and.(p0*py==0)) then
            if (epsy(xi,yi,zi) > epsout .or. &
                (epsy(xi,yi,zi) < epsin .and. epsy(xi,yi,zi) < epsmemb)) then
                write(6,*) 'PB Bomb in set_membrane_eps: epsy array value out of bounds'
                write(6,*) 'value: ',sngl(epsy(xi,yi,zi)),' anomalous at ',xi,yi,zi
                call mexit(6,1)
            end if
            end if
            if ((xi/=0 .and. yi/=0).and.(p0*pz==0)) then
            if (epsz(xi,yi,zi) > epsout .or. &
                (epsz(xi,yi,zi) < epsin .and. epsz(xi,yi,zi) < epsmemb)) then
                write(6,*) 'PB Bomb in set_membrane_eps: epsz array value out of bounds'
                write(6,*) 'value: ',sngl(epsz(xi,yi,zi)),' anomalous at ',xi,yi,zi
                call mexit(6,1)
            end if
            end if
        end if
    end do; end do
    end do

        !!!write(6,*) '-debug: done setting membrane eps';flush(6)
        !stop
end subroutine set_membrane_eps

!epsmembfracx_r
!Computes needed value for epsx array at the given node (i,j,k) using levelset
!and insas information.
!Modified version of equivalent in pb_exmol. 
!To use it here, we must consider:
!1) insas = 0 in membrane
!2) At the membrane-solute interface levelset goes to zero, but does
!   not cross. We will need to send "membroot" a fake value if want to hijack
!   it for our purposes
!There may still be a problem if we encounter a case where the three nodes
!fed to membroot each occupy a different insas region (i.e. membrane, solute, 
!and solvent). We should still be able to get a solution from membroot, but it will
!likely be inaccurate. Luckily, this should be a very rare occurance unless
!we encounter molecules with very unusual topologies.
subroutine epsmembfracx_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout,&
                         u, insas, xm, ym, zm,ipb )
!There seems to be several variables fed in that never get used. They are
!kept for now in the interest of preserving the original signature as much as
!possible, but can probably be removed eventually.
    use poisson_boltzmann, only: nbndx,nbndy,nbndz,fedgex,fedgey,fedgez

   implicit none

   ! passed variables
   integer ipb
   integer  i, j, k, a, b, a1, b1, xm, ym, zm
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)
   integer insas(0:xm+1,0:ym+1,0:zm+1) 

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t
   integer a0,ap,app,am

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0

   !just in case, we ensure that our insas region values
   !are all either -1,0,or 1
   !write(6,*) CHAR(09),CHAR(09),'-debug: epsmembfracx_r(): initializing'     
   if (abs(a) > 1) then
      a0 = sign(1,a)
   else
      a0 = a
   end if
   if (abs(b) > 1) then
      ap = sign(1,b)
   else
      ap = b
   end if   

   !We must be careful with what we feed to the call to membroot as it
   !requires a sign change in the levelset (u) to function.
   !Unfortunately, the levelset will not change sign across the
   !membrane-solute interface, so there may be a problem if that is
   !the case. We will therefore need to force a sign change in u.
   !To deal with this, we will use a convention: Since, in the
   !original insas / lvlset scheme, insas had the opposite sign
   !of lvlset, and lvlset was positive when insas was negative,
   !we will enforce that if insas at a > insas at b then lvlset
   !at a will be negative and lvlset at b will be positive
   
   !lastly, since we are feeding 3 nodes to membroot for quadratic
   !membroot finding, we need to also consider the third nodes levelset
   !value. As convention we will always compare its insas value with
   !the insas value of the central node to determine whether its
   !levelset value will be made positive or negative.

   if ( a0 > ap ) then
    !write(6,*) CHAR(09),CHAR(09),'-debug: a0>ap, setting am';flush(6)
   !here we have the central node with a larger insas value than
   !its neighbor to the right. We will thus force the central node
   !to have a negative levelset value and the node to the right
   !will be given a positive levelset value. The third point,
   !which is the left side neighbor of the central node, will have
   !its insas value compared to the central node. If it is greater
   !than or equal to the insas value at the central node, the left
   !node will be given a negative levelset value (same as central
   !node) otherwise it will be given a positive levelset value
      x1 = dble(i-1)
      x2 = dble(i  )
      x3 = dble(i+1)
      if (abs(insas(i-1,j,k)) > 1) then
         am = sign(1,insas(i-1,j,k))
      else
         am = insas(i-1,j,k)
      end if
    !write(6,*) CHAR(09),CHAR(09),'-debug: a0,ap,am ',a0,ap,am;flush(6)

      if (am >= a0) then
        f1=-abs(u(i-1,j,k))
      else
        f1= abs(u(i-1,j,k))
      end if
        f2=-abs(u(i  ,j,k))
        f3= abs(u(i+1,j,k))
    !write(6,*) CHAR(09),CHAR(09),'-debug: f1,f2,f3',sngl(f1),sngl(f2),sngl(f3);flush(6)
    !write(6,*) CHAR(09),CHAR(09),'-debug: finding edge-interface intersect';flush(6) 
      call membroot(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i+1,j,k) is inside, f2 < 0 and f1 > 0

   if ( ap > a0 ) then
   !This case has the central node with a lower insas value than the
   !node to the right. We will thus force it to have a positive levelset
   !value and give the next node a negative levelset value.
   !Keeping to the above convention, we attempt to place our third
   !interpolation value in the same region as node with the higher insas
   !value. In this case, two nodes to the right of the central node.
   !We again compare the third nodes insas value to the central node
   !to determine the sign of its levelset value to be fed to membroot.
    !write(6,*) CHAR(09),CHAR(09),'-debug: ap>a0, setting app';flush(6)
      x1 = dble(i  )
      x2 = dble(i+1)
      x3 = dble(i+2)
      if (abs(insas(i+2,j,k)) > 1) then
         app=sign(1,insas(i+2,j,k))
      else
         app=insas(i+2,j,k)
      end if
    !write(6,*) CHAR(09),CHAR(09),'-debug: a0,ap,app ',a0,ap,app;flush(6)
      f1 = abs(u(i  ,j,k))
      f2 =-abs(u(i+1,j,k))
      if (app <= a0) then
         f3 = abs(u(i+2,j,k))
      else
         f3 =-abs(u(i+2,j,k))
      end if
    !write(6,*) CHAR(09),CHAR(09),'-debug: f1,f2,f3',sngl(f1),sngl(f2),sngl(f3);flush(6)
    !write(6,*) CHAR(09),CHAR(09),'-debug: finding edge-interface intersect';flush(6)
      call membroot(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(i)
   else
      aa = dble(i+1) - t
   end if
   !This seems to be using variables from pb_exmol.
   !As there is no module available for pb_exmol, we will have keep it
   !this way... It should probably be fixed in the future in the interest
   !of a more object oriented programming style.
   !if(sasopt==2) 
   !write(6,*) CHAR(09),CHAR(09),'-debug: incrmenting nbndx';flush(6)
        nbndx=nbndx+1
   if (ipb /=4 .and. ipb/=5) then
    !write(6,*) CHAR(09),CHAR(09),'-debug: updating fedgx @',nbndx;flush(6)  
      fedgex(nbndx) = t - dble(i)
      !write(6,*) CHAR(09),CHAR(09),'-debug: returning epsx interpolated value';flush(6) 
   end if
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   !write(6,*) CHAR(09),CHAR(09),'-debug: done with epsmembfracx_r()';flush(6)
 

end subroutine epsmembfracx_r

!Just like epsmembfracx_r but along y
!see epsmembfracx_r for more information on method / modifications
subroutine epsmembfracy_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout,&
                       u, insas, xm, ym, zm,ipb )

     use poisson_boltzmann, only: nbndx,nbndy,nbndz,fedgex,fedgey,fedgez

   implicit none

   ! passed variables
   integer ipb
   integer  i, j, k, a, b, a1, b1, xm, ym, zm
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)
   integer insas(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t
   integer am,a0,ap,app

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   !force insas values to be -1,0, or 1 by sign for comparisons
   if (abs(a) > 1) then
      a0 = sign(1,a)
   else
      a0 = a
   end if
   if (abs(b) > 1) then
      ap = sign(1,b)
   else
      ap = b
   end if


   if ( a0 > ap ) then
      x1 = dble(j-1)
      x2 = dble(j  )
      x3 = dble(j+1)
      if (abs(insas(i,j-1,k)) > 1) then
         am = sign(1,insas(i,j-1,k))
      else
         am = insas(i,j-1,k)
      end if
      
      if (am >= a0) then
         f1 =-abs(u(i,j-1,k))
      else
         f1 = abs(u(i,j-1,k))
      end if
      f2 =-abs(u(i,j  ,k))
      f3 = abs(u(i,j+1,k))
      call membroot(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j+1,k) is inside, f2 < 0 and f1 > 0

   if ( ap > a0 ) then
      x1 = dble(j  )
      x2 = dble(j+1)
      x3 = dble(j+2)
      if (insas(i,j+2,k) <= ap) then
         app = sign(1,insas(i,j+2,k))
      else
         app = insas(i,j+2,k)
      end if      
      f1 = abs(u(i,j  ,k))
      f2 =-abs(u(i,j+1,k))
      if (app <= a0) then
         f3 = abs(u(i,j+2,k))
      else
         f3 =-abs(u(i,j+2,k))
      end if
      call membroot(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(j)
   else
      aa = dble(j+1) - t
   end if
   !if(sasopt==2) 
        nbndy=nbndy+1
   if (ipb /=4 .and. ipb /=5) then
      fedgey(nbndy) = t - dble(j)
   end if
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)

end subroutine epsmembfracy_r

!like epsmembfracx_r and epsmembfracy_r, but along z axis. In the
!interest of space I will ommit comments for this one
subroutine epsmembfracz_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout,&
                           u, insas, xm, ym, zm, ipb )

     use poisson_boltzmann, only: nbndx,nbndy,nbndz,fedgex,fedgey,fedgez

     

   implicit none

   ! passed variables
   integer ipb
   integer  i, j, k, a, b, a1, b1, xm, ym, zm
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)
   integer insas(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t
   integer am,a0,ap,app

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 
   if (abs(a)>1) then
      a0 = sign(1,a)
   else
      a0 = a
   end if
   if (abs(b)>1) then
      ap = sign(1,b)
   else
      ap = b
   end if

   if ( a0 > ap ) then
      x1 = dble(k-1)
      x2 = dble(k)
      x3 = dble(k+1)
      if (abs(insas(i,j,k-1))>1) then
         am = sign(1,insas(i,j,k-1))
      else
         am = insas(i,j,k-1)
      end if

      if (am >= a0) then
         f1 =-abs(u(i,j,k-1))
      else
         f1 = abs(u(i,j,k-1))
      end if
      f2 =-abs(u(i,j,k))
      f3 = abs(u(i,j,k+1))
      call membroot(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j,k+1) is inside, f2 < 0 and f1 > 0

   if ( ap > a0 ) then
      x1 = dble(k)
      x2 = dble(k+1)
      x3 = dble(k+2)
      if (abs(insas(i,j,k+2))>1) then
         app = sign(1,insas(i,j,k+2))
      else
         app = insas(i,j,k+2)
      end if
      f1 =-abs(u(i,j,k))
      f2 = abs(u(i,j,k+1))
      if ( app <= a0 ) then
         f3 =-abs(u(i,j,k+2))
      else
         f3 = abs(u(i,j,k+2))
      end if
      call membroot(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(k)
   else
      aa = dble(k+1) - t
   end if
   !if(sasopt == 2) 
        nbndz=nbndz+1
   if (ipb /= 4 .and. ipb /= 5) then
   fedgez(nbndz) =  t - dble(k)
   end if
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)


end subroutine epsmembfracz_r

subroutine membroot(x0,x1,x2,f0,f1,f2,t0)

   implicit none

   ! passed variables

   _REAL_ x0,x1,x2,f0,f1,f2,t0

   ! local variables

   _REAL_ b,c,a0,b0,c0,t,r1,r2

   !write(6,*) CHAR(09),CHAR(09),'-debug: starting membroot()';flush(6)    

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

   ! If t <=0, must be double membroot t is close to zero

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

   !write(6,*) CHAR(09),CHAR(09),'-debug: finished membroot()';flush(6)

end subroutine membroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! OLD EXCLUSION FUNCTIONs, REMOVE IF NOT REIMPLEMENTED LATER !!!!!!!!!!!!!

!This is just a wrapper for the calculation right now.
!its main job is to determine which type of pore region is being used
!and call the appropriate routine. This is necessary to allow flexibility
!in shaping the poredata array.
subroutine exclusion_calculation(i,j,k,gox,goy,goz,h,xm,ym,zm,&
                                 mprobe,poretype,poredata,eval)
    implicit none
    !passed variables
    integer     i,j,k
    integer     xm,ym,zm
    integer     poretype
    _REAL_      gox,goy,goz
    _REAL_      mprobe
    _REAL_      poredata(*)
    _REAL_      eval        !output
    _REAL_      h

    if (poretype == 1) then !cylinderical region
        call cylinder_exclusion(i,j,k,gox,goy,goz,h,xm,ym,zm,&
                mprobe,poredata,eval)
    else if (poretype == 2) then    !automated void checking
        write(6,*) 'PB Warning: Automated void checking not yet implemented'
        eval = 1.0d0
    else if (poretype == 3) then    !spherical exclusion region
        write(6,*) 'PB Warning: Spherical exlusion regions not yet implemented'
        eval = 1.0d0 
    else if (poretype == 4) then    !load from dx file
        write(6,*) 'PB Warning: External void region files not yet implemented'
        eval = 1.0d0
    else
        write(6,*) 'PB Bomb in exclusion_calculatin(): Unknown exclusion type'
        call mexit(6,1)
    end if

end subroutine exclusion_calculation

!cylinder_exclusion
!Calculates exclusion region factor for a cylinder - like region
!Inteneded to allow definition of orthoganal or skewed cylinderical
!regions, but could be used for any contiguos region defined 
!as an orthoganal or skewed solid of revolution.
!(R0(z),r(z)) -R0(z) is xy coords of region center at z
!               r(z) is radius of region at z
!The exclusion factor is then calculated using the arctan function
!such that if the grid coordinates at j,k fall within
! r(z)-mprobe/2 the exclusion value will yield ~0
!if j,k falls outside r(z)+mprobe/2 the value will yield ~1
!otherwise, there will be a fast change from ~.031 to ~.968
!For most purposes, a small mprobe should suffice since a sharp
!drop off will not make a big difference if the cylinder edges are
!all well within the solute interior. If, for some reason, this is
! or cant be the case, a lareger mprobe will give smoother transitions
!to avoid introducing additional instabilities due to ill conditioned
!derivative calculations.
subroutine cylinder_exclusion(i,j,k,gox,goy,goz,h,xm,ym,zm,mprobe,&
                                poredata,eval)

    implicit none
#include "pb_constants.h"
    !passed variables
    integer     i,j,k
    integer     xm,ym,zm
    _REAL_      gox,goy,goz,h
    _REAL_      mprobe
    _REAL_      poredata(0:zm+1,3)
                    !dimensions 1,2 store x/y coords of exclusion
                    !region centeroid at given z index
                    !dimension 3 holds radius at currend z index
    _REAL_      eval    !output

    !local variables
    _REAL_      xcrd,ycrd,dist

    xcrd = i*h+gox; ycrd = j*h+goy
    dist = sqrt( (xcrd-poredata(k,1))**2 + (ycrd-poredata(k,2))**2 )
    eval = .5 + TWO*INV_TWOPI*atan(20.0d0/mprobe*(dist-poredata(k,3)&
                -mprobe))

    
        !if ( (i==Floor((1+xm)/2.0d0)) .and. (j==Floor((1+ym)/2.0d0)) ) then
        !    write(6,*) 'k,pdata(1-3),eval:',k,&
        !                poredata(k,1:3),eval; flush(6)
        !end if
        !if ( (i==Floor((1+xm)/2.0d0)) .and. (k==Floor((1+zm)/2.0d0))) then
        !    if (j==0) then
        !        write(12345,*) 'cylinder reginon center y-line data'
        !    end if
        !    write(12345,*) 'j,pdata(1-3),eval',j,poredata(k,1:3),eval;flush(12345)
        !end if 

end subroutine cylinder_exclusion
