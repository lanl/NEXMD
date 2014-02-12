!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     +++ start of pb_fft driver +++                           !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!pb_fftsolv -   uses fast fourier transforms to solve poisson's equations for coulombic
!               potential. Coulombic potential is then used to generate surface charge
!               grid (bv). This will be used with bc option 9 and ene option 2 combination
subroutine pb_fftsolv(cggrid,gridorder,nx,ny,nz,cphi,gridh,& 
                      flags)
!phi,bv,epsx,epsy,epsz,iepsav,insas,nbnd,epsin,epsout,&
   use poisson_boltzmann, only : acrd, gox, goy, goz, fillratio, grdcrg, qgrdcrg, ngrdcrg&
                                ,bv,epsx,epsy,epsz,iepsav,insas,nbnd
   implicit none
#  include "md.h"
#  include "pb_constants.h"

    ! Common variables
    _REAL_ green(0:20, 0:20, 0:20)
    common /blk_green/ green

    ! passed variables
    _REAL_ cggrid(nx,ny,nz) !charge grid
    integer nx,ny,nz    !number gridnodes in x,y,z
    integer gridorder   !unused for now, but here as a place holder for future use, set to 1
    !integer nbnd        !number of boundar grid nodes     
    !integer insas(0:nx+1,0:ny+1,0:nz+1) !solvent accesisble surface mask 1 if in sas, else 0
    _REAL_  cphi(nx,ny,nz)  !coulombic potential grid
    !_REAL_  phi(nx,ny,nz)   !potential grid, may mean total or reaction field depending on
                            !location of corresponding nodes and choice of solver option combo
    !_REAL_  bv(nx,ny,nz)    !boundary value grid, used to store boundary conditions. I.e. charges
                            !stored on boundaries, exact usage varies with solver options
    !_REAL_  epsx(:,:,:),epsy(:,:,:),epsz(:,:,:) !dielectric constant at grid edges in x,y,z orientation
    !_REAL_  iepsav(4,*) !boundary grid node indices (iepsav(1-3,i)) and types (iepsav(4,i))
    _REAL_  epsin, epsout !dielectric constants of solute and solvent
    _REAL_  gridh       !grid spacing 
    _REAL_  fft_kappa   !ionic strength, kappa, for helmholtz formulation of linear poisson boltzmann equation
                        !currently unimplmented. See pb_fft.f
    _REAL_  gridhx, gridhy, gridhz

    ! internal variables
    integer, allocatable :: grdcrg_l(:,:)
    _REAL_, allocatable :: cgmag(:)
    _REAL_, allocatable :: cggridtemp(:,:,:)
    integer tempn
    integer ngrdcrg_l
    integer ist,counter
    _REAL_, allocatable :: fd_phi_grid(:,:,:) !store fd potential grid. only used for debugging / numerical testing
    integer xi, yi, zi
    character(LEN=34) :: phi_grid_format !high precision output format for writting phi.dat file
    character(LEN=34) :: phi_node_format !high precision output format for writting phi.dat file
    logical fdcheck_flag !run fd potential and output energies 
    logical genfile_phigrid_flag !Will write phi grid for fft and, if fdcheck_flag is on, fd as well 
    logical genfile_cgrid_flag, grid_based_fd, genfile_grdcrg_flag, genfile_grdphi_flag !output control flags
    logical internal_grdcrg, gen_surface_charge, genfile_surface_charge, write_progress
    integer flags       !Control flag code. See below for description and usage   
    _REAL_  tempcrg 
    !flag meanings:
    !   bit id      flag name               description
    !   (to turn
    !   on add #)
    !   --------    ---------               -----------
    !   1 (+1)      fdcheck_flag        -   controls whether or not to perform fd potential caluculations
    !                                        this will be stored in a seperate grid from fft calculations
    !                                        and will not be passed to the calling superroutine. It is primarily
    !                                        for the purpose of checking numerical accuracy of the routine against
    !                                        the known fd implementation. If enabled, potential grid output files
    !                                        will contain a second entry for potentials at each node. The first
    !                                        entry will be fft potential, the second will be fd potential
    !   2 (+2)      genfile_phigrid_flag-   controls output of the entire potential grid to the output 
    !                                        file: phi_grid.dat. This file will list each node's x, y and z
    !                                        indices, followed by the potential at that node. This file can quickly
    !                                        become very large so the flag should be turned off when not needed,
    !                                        particularly for very large grids as this will slow the program considerably
    !                                        and the file output can be in excess of 100 MB for grids over 200x200x200.
    !                                        In most cases, the grdphi file will be sufficient.
    !   3 (+4)      genfile_cgrid_fla   if (write_progress) then
    !   4 (+8)      grid_based_fd       -   holdover from older code. leave set to true
    !   5 (+16)     genfile_grdcrg_flag -   controls output of locations charged grid nodes and their
    !                                        corresponding charges to the outupt file: charged_nodes.dat
    !   6 (+32)     genfile_grdphi_flag -   generates a file containing locations of charged nodes and the
    !                                        potential stored in phi at the corresponding location.
    !   7 (+64)     internal_grdcrg     -   if set to true, will use internal arrays for grdcrg, qgrdcrg
    !                                        and ngrdcrg rather than modifying the existing arrays from the
    !                                        poisson boltzmann module. Should be left off if these arrays have already
    !                                        been generated elsewhere
    !   8 (+128)    gen_surface_charge  -   Generate surface charge density grid. Turn on to generate surface charge
    !                                        density grid based on phi generated in pb_fft.
    !   9 (+256)  genfile_surface_charge-   Write the surface charge density to a file
    !  10 (+512)    write_progress      -   Turns on / off write(6,*)ing of progress to mdout

    !flags 8 and 9 (surface density related) are currently non-functional. Will be implemented when / if the surface
    !charge density generation routine is completed. Also, if flag 8 is off, flag 9 will automatically be turned off
    !as well (this should be obvious, but is in place as a precaution).

        write(6,*) '-debug: pb_fftsolv: initializing';Flush(6)

        allocate(cggridtemp(nx,ny,nz))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                       BEGIN FLAG INITIALIZATION BLOCK                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write(6,*) 'Setting flags'
    if (modulo((flags/1),2) == 1) then
        fdcheck_flag = .true.
    else
        fdcheck_flag = .false.
    endif

    if (modulo((flags/2),2) == 1) then
        genfile_phigrid_flag = .true.
    else
        genfile_phigrid_flag = .false.
    endif

    if (modulo((flags/4),2) == 1) then
        genfile_cgrid_flag = .true.
    else
        genfile_cgrid_flag = .false.
    endif

    if (modulo((flags/8),2) == 1) then
        grid_based_fd = .true.
    else
        grid_based_fd = .false.
    endif

    if (modulo((flags/16),2) == 1) then
        genfile_grdcrg_flag = .true.
    else
        genfile_grdcrg_flag = .false.
    endif

    if (modulo((flags/32),2) == 1) then
        genfile_grdphi_flag = .true.
    else
        genfile_grdphi_flag = .false.
    endif

    if (modulo((flags/64),2) == 1) then
        internal_grdcrg = .true.
    else
        internal_grdcrg = .false.
    endif

    if (modulo((flags/128),2) == 1) then
        gen_surface_charge = .true.
        if (modulo((flags/256),2) == 1) then 
            genfile_surface_charge = .true.
        else
            genfile_surface_charge = .false.
        endif
    else
        gen_surface_charge = .false.
        genfile_surface_charge = .false. !only enable writting if generation
                                         !is turned on as well
    endif

    if (modulo((flags/512),2) == 1) then
        write_progress = .true.
    else
        write_progress = .false.
    endif

    if (write_progress) then
        Write(6,*) 'FD Green function checking Flag set to: ',fdcheck_flag
        Write(6,*) 'Phi grid writting flag set to: ',genfile_phigrid_flag
        Write(6,*) 'Charge grid writting Flag set to: ',genfile_cgrid_flag
        Write(6,*) 'Use grid based fd Flag set to: ',grid_based_fd
        Write(6,*) 'Charged node arrays writting Flag set to: ',genfile_grdcrg_flag
        Write(6,*) 'Charged node potential writting Flag set to: ',genfile_grdphi_flag
        Write(6,*) 'Charged node array overwritting disablement Flag set to: ',internal_grdcrg
        Write(6,*) 'Surface charge density generation Flag set to: ',gen_surface_charge
        Write(6,*) 'Surface charge density writting Flag set to: ',genfile_surface_charge
        flush(6)
    end if

    if (fdcheck_flag) then
       allocate(fd_phi_grid(nx,ny,nz))
    endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                        END FLAG INITIALIZATION BLOCK                        !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                   BEGIN GRID & ARRAY INITIALIZATION BLOCK                   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !allocate and initialize internal charge grid data arrays
   !grdcrg, qgrdcrg, and ngrdcrg are read directly from poisson boltzmann module
   !internal versions, grdcrg_l, cgmag, and ngrdcrg_l are allocated.
   !care must be taken as these were originally assumed to be those
   !used by the poisson boltzmann module. If internal_grdcrg flag is set to off
   !then we proceed under the assumption that the they will have shapes identicle
   !to their external counterparts. However, since this method is now used outside
   !of pb_force, this may not be the case:
   !  In the case where internal_grdcrg flag is turned on, the internal arrays will
   !  need to be generated based on the total number of charged nodes in the grid
   !  that has been fed in. If boundary charges have been added, this will likely
   !  be different than the original charge grid which was generated directly from
   !  the charged atoms in the solute.
   !In either case, the build_charge_arrays will be used to construct the internal
   !arrays. If the internal_grdcrg flag is off, the grdcrg and qgrdcrg arrays
   !will then be overwritten with the internally generated versions.
   if (write_progress) then
    write(6,*) 'entering array allocation block'
    flush(6)
   end if
   if (internal_grdcrg) then
    tempn = 0
    do xi=1,nx; do yi = 1,ny; do zi=1,nz !count number of charged grid nodes
        if (cggrid(xi,yi,zi) /= 0.0d0) then
            tempn = tempn + 1
        end if
    end do; end do; end do 
    if (write_progress) then
        write(6,*) 'there are ',tempn,' charged gridnodes'
        write(6,*) 'allocating memory for local charged grid node data arrays'
        flush(6)
    end if
      allocate ( grdcrg_l(3,tempn), stat=ist)
      if ( ist /= 0 ) then 
         write(6,*),'allocate grdcrg_l(3,tempn) failed',tempn,ist
         stop
        else
        if (write_progress) then
         write(6,*),'allocation of grdcrg_l succeeded'
         flush(6)
        endif
      end if
     allocate( cgmag(tempn), stat=ist )
      if ( ist /= 0 ) then 
         write(6,*),'allocate cgmag(tempn) failed',tempn,ist
         stop
        else
        if (write_progress) then
         write(6,*),'allocation of cgmag succeeded'
         flush(6)
        endif
      end if
     if (write_progress) then
        write(6,*) 'building charged grid node arrays from charge grid'
        flush(6)
     end if
    call build_charge_arrays(cggrid,grdcrg_l,cgmag,ngrdcrg_l,nx,ny,nz)
   else
     tempn = SIZE(grdcrg,1)
     if (write_progress) then
        write(6,*) 'grdcrg holds ',tempn,' node index sets'
        write(6,*) 'allocating memory for local charge grid node data arrays'
        flush(6)
     end if
      allocate ( grdcrg_l(3,tempn), stat=ist)
      if ( ist /= 0 ) then 
         write(6,*),'allocate grdcrg_l(3,tempn) failed',tempn,ist
         stop
      else
        if (write_progress) then
         write(6,*),'allocation of grdcrg_l succeeded'
         flush(6)
        endif
      end if
     allocate( cgmag(tempn), stat=ist )
      if ( ist /= 0 ) then 
         write(6,*),'allocate cgmag(tempn) failed',tempn,ist
         stop
        else
        if (write_progress) then
         write(6,*),'allocation of cgmag succeeded'
         flush(6)
        endif
      end if
     if (write_progress) then
        write(6,*) 'building charged grid node arrays from charge grid'
        flush(6)
     end if
     call build_charge_arrays(cggrid,grdcrg_l,cgmag,ngrdcrg_l,nx,ny,nz)
     if (write_progress) then
        write(6,*) 'writting local arrays to global copies'
        flush(6)
     end if
     grdcrg = grdcrg_l
     ngrdcrg = ngrdcrg_l
     qgrdcrg = cgmag
   endif

   if (write_progress) then
    write(6,*) 'initializing phi grid format'
    flush(6)
   end if

   if (fdcheck_flag) then

    phi_grid_format="(I7,I7,I7,ES24.12E3,ES24.12E3)"
    phi_node_format="(I7,I7,I7,I7,ES24.12E3,ES24.12E3)"
   else
    phi_grid_format="(I7,I7,I7,ES24.12E3)"
    phi_node_format="(I7,I7,I7,I7,ES24.12E3)"
   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !             BEGIN CHARGE GRID DATA ANALYSIS / WRITTING BLOCKS               !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tempcrg=SUM(cggrid)
    if (write_progress) then
        write(6,*) 'Analyzing grid data'
        write(6,*) 'Net charge of cggrid: ', tempcrg
        write(6,*) 'fill ratio set to: ', fillratio
        write(6,*) 'nodes per dim. x,y,z ',nx,ny,nz
        write(6,*) 'max net charge magnitued threshold is:',nx*ny*nz*1.0d-15 
        flush(6)
    end if
   if (abs(tempcrg) >= (nx*ny*nz*1.0d-15)) then
        if (write_progress) then
                write(6,*) 'net charge is greater than threshold'
                write(6,*) 'neutralizing grid'
                flush(6)
        end if
        cggridtemp = cggrid-(tempcrg/(nx*ny*nz))
       ! do xi=1,nx; do yi=1,ny; do zi=1,nz
       !         cggridtemp(xi,yi,zi)=cggrid(xi,yi,zi)-tempcrg/(nx*ny*nz)
       ! end do; end do; end do;
        if (write_progress) then
            write(6,*) 'added ',tempcrg/(nx*ny*nz),'to each grid node'
            write(6,*) 'net charge is now: ',SUM(cggridtemp)
        end if 
   else
        cggridtemp = cggrid
   end if
   gridhx = gridh
   gridhy = gridh
   gridhz = gridh
    if (write_progress) then
        write(6,*) 'grid cell volume= gridhx*gridhy*gridhz = ',sngl(gridhx*gridhy*gridhz)
        write(6,*) 'gridspacing (x,y,z)', sngl(gridhx),sngl(gridhy),sngl(gridhz)
        flush(6)
    end if
   if ( genfile_cgrid_flag ) then !write charge grid to file if debugging
    if (write_progress) then 
        write(6,*) 'Writting Charge Grid Data File'; flush(6) 
    end if
    call genfile_cgrid(cggridtemp,nx,ny,nz) 
    write(6,*) 'done writting'; flush(6)
   end if
   if ( genfile_grdcrg_flag ) then !write charged node data arrays to file
        if (write_progress) then
            write(6,*) 'Writting charged node data file'; flush(6)
            flush(6)
        end if
        call genfile_grdcrg_l(grdcrg_l,cgmag,ngrdcrg_l)
   end if 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                  END GRID DATA ANALYSIS / WRITTING BLOCKS                   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                   BEGIN FFT PHI GRID CALCULATION BLOCK                      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !generate coulombic potential

    if ( write_progress ) then
        write(6,*) 'calling pb_fft'; flush(6)
    end if
    write(6,*) '-debug: pb_fftsolv: calling pb_fft';flush(6)
    fft_kappa = 0.0d0 !turned off helmholtz solve capabilities till needed
    epsin = 1 !arbitrarily set to vaccuum case for testing purposes - fd assumes eps=1
#ifdef FFTW
    call pb_fft(cggridtemp,nx,ny,nz,cphi,gridh,epsin,fft_kappa)
#else
    write(6,*) "PBSA FFT solver is disabled, probably because gcc <= 4.2"
#endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                    END FFT PHI GRID CALCULATION BLOCK                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                   BEGIN FD PHI GRID CALCULATION BLOCK                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !+ FD Calculation based on FD Green's function                               +!
   !+ Computes FD Phi grid at each grid node by looping over all charged nodes  +!
   !+ and applying the FD Green's Function if charge is less than 20 grid units +!
   !+ away. If charge is 20 or more grid units away, the standard 1/r potential +!
   !+ is applied. This method was previously implemented here directly, but has +!
   !+ been plut in its own subroutine now for the sake of modularity/readability+!
   if ( fdcheck_flag ) then !fd calculation, used for debugging
    if (write_progress) then
        write(6,*) 'generating fd potential map'
        call flush(6)
    end if
    call gen_fd_phi_grid(cggridtemp,fd_phi_grid,grdcrg_l,cgmag,ngrdcrg_l,gridh,&
                         nx,ny,nz,tempcrg,grid_based_fd)
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                     BEGIN FD PHI GRID CALCULATION BLOCK                     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                        BEGIN PHI GRID WRITTING BLOCKS                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( genfile_phigrid_flag ) then
        if (write_progress) then
            write(6,*) 'writting potential map data files'
            call flush(6)
        end if
        if ( fdcheck_flag ) then
            call genfile_joint_phigrid(cphi,fd_phi_grid,nx,ny,nz,phi_grid_format)
        else
            call genfile_fft_phigrid(cphi,nx,ny,nz,phi_grid_format)
        end if
   end if
   if ( genfile_grdphi_flag ) then
        if (write_progress) then
            write(6,*), 'writting charged node potential data files'
            call flush(6)
        end if
        if ( fdcheck_flag ) then
            call genfile_joint_grdphi(cphi,fd_phi_grid,grdcrg_l,ngrdcrg_l,phi_node_format)
        else
            call genfile_fft_grdphi(cphi,grdcrg_l,ngrdcrg_l,phi_node_format)
        endif
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                         END PHI GRID WRITTING BLOCKS                        !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                BEGIN SURFACE CHARGE DENSITY GENERATION BLOCKS               !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !interpolate coulombic potential to generate dielectric boundary surface charges
        !last term is for grid interpolation order. Used as a place holder for future usage
    if (gen_surface_charge) then
        call gen_dbcrg( grdcrg,qgrdcrg,insas,epsx,epsy,epsz,bv,cphi,nx,ny,nz,nbnd,1 )
    end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                 END SURFACE CHARGE DENSITY GENERATION BLOCKS               !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (write_progress) then
        write(6,*), 'deallocating arrays'
        call flush(6)
    end if 
    deallocate(grdcrg_l)
    deallocate(cgmag)
    if (write_progress) then
        write(6,*), 'exiting pb_fftsolv'
        call flush(6)
    end if
contains
!++++++++++++++++++++++++++Charge array building subroutine+++++++++++++++++++
!+ Builds up arrays to store information on charged grid nodes
!+ grdcrg_l holds indices of each node with nonzero charge
!+ cgmag (alias of grdcrg_l in pb_force) contains magnitudes of charges at the
!+       corresponding node
!+ ngrdcrg_l is the total number of nodes with nonzero charge
!+ These arrays are not essential to building phi grid with fft method, but
!+ will be employed in energy calculations and, if fdcheck_flag is enabled
!+ they will be used by the fd phi grid building method
!+ This method is essentialy an internal copy of block of code in pb_fddrv.f
!+ used to do this task.
subroutine build_charge_arrays(cggrid,grdcrg_l,cgmag,ngrdcrg_l,nx,ny,nz)
    ! Passed Variables:
    _REAL_ cggrid(nx,ny,nz)
    integer grdcrg_l(3,*)
    _REAL_ cgmag(*)
    integer ngrdcrg_l, nx, ny, nz
    
    !Internal Variables:
    integer xi, yi, zi !counters

    ngrdcrg_l = 0 !initializ number of charged grid nodes to zero
    do xi=1,nx; do yi=1,ny; do zi=1,nz
        if ( cggrid(xi,yi,zi) .ne. 0.0 ) then
        !loop over all grid nodes and check for nonzero charge
            !if charge at grid node is nonzero, add values to arrays accordingly
            ngrdcrg_l = ngrdcrg_l + 1
            grdcrg_l(1,ngrdcrg_l) = xi
            grdcrg_l(2,ngrdcrg_l) = yi
            grdcrg_l(3,ngrdcrg_l) = zi
            cgmag(ngrdcrg_l) = cggrid(xi,yi,zi)
        end if
    end do; end do; end do

end subroutine build_charge_arrays

!+++++++++++++++++Charge grid file writting subroutine++++++++++++++++++++++++
!+ Writes charge grid to the data file to 'pb_cg.dat'
subroutine genfile_cgrid(cggrid,nx,ny,nz)

    !Passed Variables
    _REAL_ cggrid(nx,ny,nz)
    integer nx, ny, nz
      open(64,file='pb_cg.dat',form="formatted")
      write(64,*) '# charge grid data'
      write(64,*) '# col 1 - 3 x,y,z locations, col 4 charge value'
      write(64,*) '# written in ascending order of z, y, x'
      do zi = 1,nz
         do yi = 1,ny
            do xi = 1,nx
               write(64,*) xi,yi,zi,cggrid(xi,yi,zi)
            end do
         end do 
      end do
      close(64)

end subroutine genfile_cgrid

!++++++++++++++FD Based Phi Grid Calculation Subroutine+++++++++++++++++++++++++
!+ Computes Phi grid using the fd green's funtion
subroutine gen_fd_phi_grid(cggrid,fd_phi_grid,grdcrg_l,cgmag,ngrdcrg_l,gridh,&
                           nx,ny,nz,tempcrg,grid_based_fd)
    !Passed Variable
    _REAL_ fd_phi_grid(nx,ny,nz)
    integer grdcrg_l(3,*)
    _REAL_ cgmag(*)
    _REAL_ gridh
    _REAL_ cggrid(nx,ny,nz), tempcrg
    integer ngrdcrg_l, nx, ny, nz
    logical grid_based_fd

    !Internal Variables
    integer xi, yi, zi, xj, yj, zj, qi !counters
    _REAL_ dr, qq, qpot
    integer qxi, qyi, qzi, dxi, dyi, dzi

    if ( (abs(tempcrg) <= nx*ny*nz*1.0d-32) .AND. grid_based_fd ) then 
        write(6,*) 'net grid charge magnitude within threshold and grid based fd on.',&
                   ' Using grd arrays'
    do xi=1,nx; do yi=1,ny; do zi=1,nz
        fd_phi_grid(xi,yi,zi) = 0.0 !initialize to zero
        qpot = 0.0d0
        call flush(6)
        do qi=1,ngrdcrg_l !phi accumulation. Loop over all charges in grdcrg_l
            !grab node indices qxi, qyi, qzi, and node charge qq
            qxi=grdcrg_l(1,qi)
            qyi=grdcrg_l(2,qi)
            qzi=grdcrg_l(3,qi)
            qq=cgmag(qi)
            call flush(6)
            !compute distance from charged node to phi grid node
            dxi = xi - qxi
            dyi = yi - qyi
            dzi = zi - qzi
            dr = Sqrt(dxi**2.0d0 + dyi**2.0d0 + dzi**2.0d0)
            if (dr < 20) then
                qpot = qq*green(abs(dxi),abs(dyi), abs(dzi))
            else
                qpot = qq/(dr)
            end if
            fd_phi_grid(xi,yi,zi) = fd_phi_grid(xi,yi,zi) + qpot
        end do
        fd_phi_grid(xi,yi,zi) = fd_phi_grid(xi,yi,zi)*INV_FOURPI/gridh
    end do; end do; end do
    call flush(6)
    else
        write(6,*) 'net grid charge maginitude too high or grid based fd off.',&
                   ' Using full grid fd'
        do xi=1,nx; do yi=1,ny; do zi=1,nz
           fd_phi_grid(xi,yi,zi) = 0.0
           qpot = 0.0d0
           do xj=1,nx; do yj=1,ny; do zj=1,nz
                dxi = xi-xj
                dyi = yi-yj
                dzi = zi-zj
                qq=cggrid(xj,yj,zj)
                dr = Sqrt(dxi**2.0d0 + dyi**2.0d0 + dzi**2.0d0)
                if (dr < 20) then
                    qpot = qq*green(abs(dxi),abs(dyi), abs(dzi))
                else
                    qpot = qq/(dr)
                end if
            fd_phi_grid(xi,yi,zi) = fd_phi_grid(xi,yi,zi) + qpot
           end do; end do; end do
           fd_phi_grid(xi,yi,zi) = fd_phi_grid(xi,yi,zi)*INV_FOURPI/gridh
        end do; end do; end do
    end if
        write(6,*) 'done making fd_phi_grid'; flush(6)
end subroutine gen_fd_phi_grid

!++++++++++++++++++++Phi Grid Data File Writting Subroutines++++++++++++++++++++!
!++++++++++++++++++++++++++joint FD & FFT Phi Grid++++++++++++++++++++++++++++++
subroutine genfile_joint_phigrid(phi_grid,fd_phi_grid,nx,ny,nz,phi_grid_format)
    !Passed Variables
    _REAL_ phi_grid(nx,ny,nz)
    _REAL_ fd_phi_grid(nx,ny,nz)
    integer nx, ny, nz
    character(LEN=33) phi_grid_format

    !Internal Variables
    integer xi, yi, zi
    
    open(101,file="phi_grid.dat",form="formatted")
    write(101,*) '# X index, Y index, Z index, fft phi, fd phi'
    do xi=1,nx; do yi=1,ny; do zi=1,nz
        write(101,phi_grid_format) xi, yi, zi, phi_grid(xi,yi,zi), fd_phi_grid(xi,yi,zi)
    end do; end do; end do
    close(101)

end subroutine genfile_joint_phigrid

!++++++++++++++++++++++++fft Phi Grid++++++++++++++++++++++++++++++++++++++++++
subroutine genfile_fft_phigrid(phi_grid,nx,ny,nz,phi_grid_format)
    !Passed Variables
    _REAL_ phi_grid(nx,ny,nz)
    integer nx, ny, nz
    character(LEN=33) phi_grid_format

    !Internal Variables
    integer xi, yi, zi
    
    open(101,file="phi_grid.dat",form="formatted")
    write(101,*) '# X index, Y index, Z index, fft phi'
    do xi=1,nx; do yi=1,ny; do zi=1,nz
        write(101,phi_grid_format) xi, yi, zi, phi_grid(xi,yi,zi)
    end do; end do; end do
    close(101)


end subroutine genfile_fft_phigrid

!+++++++++++++++++++Charged Node Data File Writting Subroutines+++++++++++++++++!
!++++Write Data File for  Joint FD and FFT Potentials at charged nodes++++++++++
subroutine genfile_joint_grdphi(phi_grid,fd_phi_grid,grdcrg_l,ngrdcrg_l,phi_node_format)
    !Passed Variables
    _REAL_ phi_grid(nx,ny,nz)
    _REAL_ fd_phi_grid(nx,ny,nz)
    integer grdcrg_l(3,ngrdcrg_l)
    integer ngrdcrg_l
    character(LEN=33) phi_node_format

    !Internal Variables
    integer icrg

    open(101,file="charge_node_phi.dat",form="formatted")
    write(101,*) "# grdcrg_l id, x index, yindex, fft phi, fd phi"
    do icrg=1,ngrdcrg_l
        write(101,phi_node_format) icrg,grdcrg_l(1,icrg),grdcrg_l(2,icrg),grdcrg_l(3,icrg),&
            phi_grid(grdcrg_l(1,icrg),grdcrg_l(2,icrg),grdcrg_l(3,icrg)),&
            fd_phi_grid(grdcrg_l(1,icrg),grdcrg_l(2,icrg),grdcrg_l(3,icrg))
    end do
    close(101)
end subroutine genfile_joint_grdphi

!++++Write Data File for FFT Potentials at charged nodes++++++++++
subroutine genfile_fft_grdphi(phi_grid,grdcrg_l,ngrdcrg_l,phi_node_format)
    !Passed Variables
    _REAL_ phi_grid(nx,ny,nz)
    integer grdcrg_l(3,ngrdcrg_l)
    integer ngrdcrg_l
    character(LEN=33) phi_node_format

    !Internal Variables
    integer icrg

    open(101,file="charge_node_phi.dat",form="formatted")
    write(101,*) "# grdcrg_l id, x index, yindex, fft phi"
    do icrg=1,ngrdcrg_l
        write(101,phi_node_format) icrg,grdcrg_l(1,icrg),grdcrg_l(2,icrg),grdcrg_l(3,icrg),&
            phi_grid(grdcrg_l(1,icrg),grdcrg_l(2,icrg),grdcrg_l(3,icrg))
    end do
    close(101)
end subroutine genfile_fft_grdphi

!+++++++++++++++++Write Data File for Charged Grid Nodes++++++++++++++++++++++++
subroutine genfile_grdcrg_l(grdcrg_l,cgmag,ngrdcrg_l)
    !Passed Variables
    integer grdcrg_l(3,ngrdcrg_l)
    _REAL_ cgmag(ngrdcrg_l)
    integer ngrdcrg_l

    !Internal Variables
    integer icrg

    open(101,file="charged_nodes.dat",form="formatted")
    write(101,*) "# grdcrg_l id, x index, y index, z index, charge"
    do icrg=1,ngrdcrg_l
        write(101,*) icrg, grdcrg_l(1,icrg), grdcrg_l(2,icrg), grdcrg_l(3,icrg),cgmag(icrg)
    end do
    close(101)
end subroutine genfile_grdcrg_l

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ adds potential for atom to eel
subroutine add_atom_potential(eel,enAt,cggrid,atom,atomcoords,cgmag,phi_grid,pb_bounds,gridhx,gridhy,gridhz,nx,ny,nz)

   ! Passed variables:
   integer nx,ny,nz, atom
   _REAL_ eel,gridhx,gridhy,gridhz
   _REAL_ atomcoords(3,*),phi_grid(nx,ny,nz),pb_bounds(3,2),cgmag(*)
   _REAL_ cggrid(nx,ny,nz)
   _REAL_ enAt
   _REAL_ self_energy

   ! Local variables
   integer xindex,yindex,zindex
   _REAL_ hvol, dconst

   ! dconst is used a coefficient to subtract the atoms self energy     +!
   ! from the potential (inclusion of self energy is a necessary result +!
   ! of implementing discretized fft methods to solve poisson's eqn.    +!
   hvol = (gridhx * gridhy * gridhz);
   dconst = (hvol ** (1.0d0 / 3.0d0)) * 3.1488d-1;
   
   xindex = floor((atomcoords(1,atom)-pb_bounds(1,1))/gridhx)+ 1;
   yindex = floor((atomcoords(2,atom)-pb_bounds(2,1))/gridhy)+ 1;
   zindex = floor((atomcoords(3,atom)-pb_bounds(3,1))/gridhz)+ 1;

   enAt = (18.2223d0 ** 2)*cgmag(atom)*(phi_grid(xindex,yindex,zindex))! - cgmag(atom) / (dconst*hvol));
   eel  = eel + enAt
   !subtract out discretized self energies if not comparing to fd total energies
   self_energy = (18.2223d0 ** 2)*cgmag(atom)*(cggrid(xindex,yindex,zindex)) / (dconst)
   if ( .not. fdcheck_flag) then
      eel = eel - self_energy
   endif
   
   !write potentials of each atom to output file. Used to debug with fdcheck_flag
   !if ( fdcheck_flag ) then
   write(6,*) 'adding potential for atom ', atom
   write(6,*) '   grid node: ', xindex, yindex, zindex
   write(6,*) '   location: ', sngl(atomcoords(1,atom)),sngl(atomcoords(2,atom)),sngl(atomcoords(3,atom))
   write(6,*) '   phi: ', phi_grid(xindex,yindex,zindex);
   write(6,*) '   node self potential: ', self_energy
   !write(6,*) '   corrected phi: ', (phi_grid(xindex,yindex,zindex) - cgmag(atom) / (dconst*hvol))
   write(6,*) '   atom charge: ', cgmag(atom)
   write(6,*) '   node charge: ', cggrid(xindex,yindex,zindex)
   write(6,*) '   uncorrected energy: ', phi_grid(xindex,yindex,zindex)*cgmag(atom)&
                                         *AMBER_ELECTROSTATIC2
   write(6,*) '   corrected energy:   ', phi_grid(xindex,yindex,zindex)*cgmag(atom)&
                                         *AMBER_ELECTROSTATIC2 - self_energy
   !endif
end subroutine add_atom_potential

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculates force on given atom... remains to be implemented
subroutine calc_force(pbfrc,iatm,grdcrg_l,phi_grid,pb_bounds,gridhx,gridhy,gridhz,nx,ny,nz);

   _REAL_ pbfrc(3,*)
   integer iatm, nx, ny, nz
   _REAL_ grdcrg_l(3)
   _REAL_ phi_grid(nx,ny,nz)
   _REAL_ gridhx, gridhy, gridhz
   _REAL_ pb_bounds(3,2)

end subroutine calc_force

!                      +++end of pb_fft driver+++                             !

!+ Dielectric boundary charge assignment routine. Adds charges to the boundary
!+ (bv) according to singularity free charge boundary formulation.
!+ See paper: Qin Cai et. al. "On Removal of Charge Singularity in Poisson-Boltzmann Equation"
!+ In brief:
!+          *itterate over all boundary charges, and for each boundary charge bdi:
!+             sum over all adjacent nodes, adji, lying inside solute region,
!+                 adding the inverse of the product of the coulombic potential at adji 
!+                     times the sum of all adjacent grid edge dielectric values for adji
!+             if bdi is not in ion accesisble region, add product of coulombic potential at bdi
!+                times sum of all connected grid edge dielectric values for bi                               
subroutine gen_dbcrg( grdcrg,qgrdcrg,insas,epsx,epsy,epsz,bv,cphi,nx,ny,nz,nbnd,gridorder )

    !passed variables
    integer insas(0:nx+1,0:ny+1,0:nz+1) !marks whether grid nodes are in solvent region
    integer nx, ny, nz 
    integer checksum
    integer nbnd                        !number of grid nodes on dielectric boundary
    integer gridorder                   !grid interpolation order, just a place holder for now
    integer grdcrg(3,*)
    _REAL_  qgrdcrg(*)
    _REAL_  epsx(0:nx,1:ny,1:nz) 
    _REAL_  epsy(0:nx,1:ny,1:nz)
    _REAL_  epsz(0:nx,1:ny,1:nz)
    _REAL_  bv(nx,ny,nz) !output
    _REAL_  cphi(nx,ny,nz) !input, coulombic potential
    _REAL_  gridh

    !internal variables
    integer xi, yi, zi, xi0, yi0, zi0, bndi, dxi, dyi, dzi
    _REAL_  qtemp, qtemp0, epstemp
    
    do bndi = 1,bndi !itterate over each node, bndi, in dielectric boundary
        xi0 = iepsav(1,bndi); yi0 = iepsav(2,bndi); zi0 = iepsav(3,bndi)
        qtemp = 0
        qtemp0 = 0
        do dxi = -1,1; do dyi = -1,1; do dzi = -1,1 !itterate over bndi and adjacent nodes
            checksum = abs(dxi) + abs(dyi) + abs(dzi) !use checksum to exclude corners
            if ( checksum <= 1 ) then
                xi = xi0 + dxi
                yi = yi0 + dyi
                zi = zi0 + dzi
                if ( insas(xi,yi,zi) == 0 ) then !only include effects of nodes inside solute region
                    epstemp = epsx(xi,yi,zi) + epsx(xi-1,yi,zi) + epsy(xi,yi,zi)&
                            + epsy(xi,yi-1,zi) + epsz(xi,yi,zi) + epsz(xi,yi,zi-1)
                    if ( checksum == 0 ) then !check if this is central node, bndi, or an adjacent node
                        qtemp = qtemp + cphi(xi,yi,zi)*epstemp !add central node product
                    else 
                        qtemp = qtemp - cphi(xi,yi,zi)*epstemp !subtract adjacent node products
                    end if
                end if
            end if
        bv(xi0,yi0,zi0) = bv(xi0,yi0,zi0) + qtemp
        end do; end do; end do
    end do    

end subroutine gen_dbcrg

end subroutine pb_fftsolv

