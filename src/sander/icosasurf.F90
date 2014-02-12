! <compile=optimized>
#include "copyright.h"
module icosasurf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! --- Module for icosasurf algorithm
!     Icosasurf calculates approximations of atomic SAS's
!     The implementation follows an idea of M. Rarey
!
!     Holger Gohlke
!     31.12.2003
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "assert.fh"
#include "dprec.fh"

! --- Icosahedron used for approximation
integer, private, dimension(20,3), parameter :: &
   ictri = reshape( (/0, 4, 1,  0, 9, 4,  9, 5, 4,  4, 5, 8,  4, 8, 1, &
                      8,10, 1,  8, 3,10,  5, 3, 8,  5, 2, 3,  2, 7, 3, &
                      7,10, 3,  7, 6,10,  7,11, 6, 11, 0, 6,  0, 1, 6, &
                      6, 1,10,  9, 0,11,  9,11, 2,  9, 2, 5,  7, 2,11/), &
                     shape(ictri), order = (/2, 1/) )

_REAL_, private, dimension(12,3), parameter :: &
   icosasym = reshape( (/-0.525731112119133606, 0.0, 0.850650808352039932, &
                          0.525731112119133606, 0.0, 0.850650808352039932, &
                         -0.525731112119133606, 0.0,-0.850650808352039932, &
                          0.525731112119133606, 0.0,-0.850650808352039932, &
      
                          0.0, 0.850650808352039932, 0.525731112119133606, &
                          0.0, 0.850650808352039932,-0.525731112119133606, &
                          0.0,-0.850650808352039932, 0.525731112119133606, &
                          0.0,-0.850650808352039932,-0.525731112119133606, &
      
                          0.850650808352039932, 0.525731112119133606, 0.0, &
                         -0.850650808352039932, 0.525731112119133606, 0.0, &
                          0.850650808352039932,-0.525731112119133606, 0.0, &
                         -0.850650808352039932,-0.525731112119133606, 0.0/), &
                       shape(icosasym), order = (/2, 1/) )

_REAL_, private, dimension(12,3) :: icosa

! --- Very low level approximation 
!     SAS_MAXLEVEL          2     SAS_MINLEVEL          0
!     MAX_SURFACE_PARTS   320     MIN_SURFACE_PARTS    20  
! --- Low level approximation 
!     SAS_MAXLEVEL          3     SAS_MINLEVEL          1
!     MAX_SURFACE_PARTS  1280     MIN_SURFACE_PARTS    80 
! --- High level approximation 
!     SAS_MAXLEVEL          4     SAS_MINLEVEL          3
!     MAX_SURFACE_PARTS  5120     MIN_SURFACE_PARTS  1280
integer, private :: ismin
integer, private :: ismax
integer, private :: ipmin
integer, private :: ipmax

! --- Water radius
_REAL_, private :: srad

! --- Max nof neighboring atoms including an icosahedra point
integer, private, parameter :: iemax = 200

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initializes icosa variables, checks consistency
subroutine icosa_init(isasmin, isasmax, sasrad)

      implicit none
      integer i
      integer isasmin, isasmax
      _REAL_ sasrad, angle, axis, rotmat, point
      dimension axis(3), rotmat(3,3), point(3)

      angle = 11.5
      axis(1) = 1.0
      axis(2) = 1.0
      axis(3) = 1.0

      ! --- Rotate icosahedra so that symmetry wrt coordinate axes is broken

      call gen_rot_mat(angle, axis, rotmat)
      do i=1,12
        point(1) = icosasym(i,1)
        point(2) = icosasym(i,2)
        point(3) = icosasym(i,3)
        call rot_point(point, rotmat)
        icosa(i,1) = point(1) 
        icosa(i,2) = point(2) 
        icosa(i,3) = point(3) 
      end do       

      if(isasmin.lt.0 .or. isasmax.lt.isasmin .or. &
         sasrad.lt.0.0d0) then
        write(6,*) 'Wrong input parameters for icosa algorithm'
        write(6,*) 'isasmin: ', isasmin
        write(6,*) 'isasmax: ', isasmax
        write(6,*) 'sasrad: ', sasrad
        call mexit(6,1)
      end if

      ismin = isasmin
      ismax = isasmax
      srad = sasrad

      ipmin = 20
      do i=1,ismin
        ipmin = ipmin * 4
      end do
      ipmax = 20
      do i=1,ismax
        ipmax = ipmax * 4
      end do

      return
end subroutine icosa_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Extracts actual atom environment from ineighbor
subroutine icosa_atom_environment(ineighborpt, ineighbor, iatenvcnt, iatenv)

      implicit none
      integer ineighborpt, ineighbor
      integer iatenvcnt, iatenv
      dimension ineighbor(*), iatenv(*)

      iatenvcnt = 0

      do while(ineighbor(ineighborpt).gt.0)
        iatenvcnt = iatenvcnt + 1
        if(iatenvcnt .gt. iemax) then
          write(6,*) 'iatenvcnt .gt. iemax in icosasurf'
          call mexit(6,1)
        end if
        iatenv(iatenvcnt) = ineighbor(ineighborpt)
        ineighborpt = ineighborpt + 1
      end do
      ineighborpt = ineighborpt + 1

      return
end subroutine icosa_atom_environment

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determines if point (x,y,z) outside any neighboring atom
logical function icosa_point_exclusion(x, vdwrad, &
                                       iatenvcnt, iatenv, &
                                       iatinclcnt, iatincl, &
                                       xp, yp, zp)

      implicit none
      integer iatenvcnt, iatenv
      integer iatinclcnt, iatincl   ! Contains pointers to atoms that include (x,y,z) on output
      integer i, j
      _REAL_ x, vdwrad, r2, vdw
      _REAL_ xp, yp, zp
      _REAL_ xj, yj, zj
      logical binside
      dimension x(*), vdwrad(*), iatenv(*), iatincl(*)

      binside = .false.
      iatinclcnt = 0

      do i=1, iatenvcnt
        j = iatenv(i)
        xj = x(3*j-2)
        yj = x(3*j-1)
        zj = x(3*j  )
        r2 = (xp-xj)**2 + (yp-yj)**2 + (zp-zj)**2
        vdw = vdwrad(j) + srad
        if(r2.lt.(vdw*vdw)) then
          binside = .true.
          iatinclcnt = iatinclcnt + 1
          iatincl(iatinclcnt) = j
        end if
      end do

      icosa_point_exclusion = .not. binside

      return
end function icosa_point_exclusion

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Recursively calcs contribution of surface patch
recursive &
function icosa_patch_approx(i, x, vdwrad, &
                            v1, v2, v3, &
                            exc1, exc2, exc3, &
                            iatinclcnt1, iatincl1, &
                            iatinclcnt2, iatincl2, &
                            iatinclcnt3, iatincl3, &
                            iatenvcnt, iatenv, &
                            nof_parts, &
                            idecomp) &
       result(ipapprox)

      use decomp, only: decsasa_icosa
      use constants, only : FOURPI
      implicit none
      integer i, idecomp
      integer iatinclcnt1, iatincl1
      integer iatinclcnt2, iatincl2
      integer iatinclcnt3, iatincl3
      integer iatinclcnt12, iatincl12
      integer iatinclcnt23, iatincl23
      integer iatinclcnt31, iatincl31
      integer iatenvcnt, iatenv
      integer nof_parts, nof_exc
      _REAL_ x, vdwrad
      _REAL_ v1, v2, v3
      _REAL_ v12, v23, v31
      _REAL_ sas, total_sas
      _REAL_ norm_x, norm_y, norm_z, normsum, normi
      _REAL_ xi, yi, zi, r
      _REAL_ ipapprox
      logical exc1, exc2, exc3
      logical exc12, exc23, exc31
      dimension x(*), vdwrad(*)
      dimension v1(*), v2(*), v3(*)
      dimension iatincl1(*), iatincl2(*), iatincl3(*)
      dimension iatincl12(iemax), iatincl23(iemax), &
                iatincl31(iemax)
      dimension iatenv(*)
      dimension v12(3), v23(3), v31(3)  

      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      r = vdwrad(i)
      total_sas= FOURPI*r*r / dble(nof_parts)

      if(nof_parts .ge. ipmax) then

        ! ---   finish recursion

        !write(6,*) 'EOR ', i, iatinclcnt1, iatinclcnt2, iatinclcnt3

        nof_exc = 0
        if(exc1) then
          nof_exc = nof_exc + 1
        else if(idecomp.gt.0) then
          call decsasa_icosa(i,iatinclcnt1,iatincl1, &
                             -total_sas/3.0d0, &
                             idecomp)
        endif
        if(exc2) then
          nof_exc = nof_exc + 1
        else if(idecomp.gt.0) then
          call decsasa_icosa(i,iatinclcnt2,iatincl2, &
                             -total_sas/3.0d0,&
                             idecomp)
        endif
        if(exc3) then
          nof_exc = nof_exc + 1
        else if(idecomp.gt.0) then
          call decsasa_icosa(i,iatinclcnt3,iatincl3, &
                             -total_sas/3.0d0, &
                             idecomp)
        endif

        sas = 0.0d0
        if(nof_exc .eq. 1) then
          sas = total_sas * 0.3333333333d0  ! before: * 0.25d0
        else if(nof_exc .eq. 2) then
          sas = total_sas * 0.6666666666d0  ! before: * 0.75d0
        else if(nof_exc .eq. 3) then
          sas = total_sas
        end if
      else

        ! ---   next level of approximation

        norm_x = v1(1) + v2(1)
        norm_y = v1(2) + v2(2)
        norm_z = v1(3) + v2(3)
        normsum = norm_x**2 + norm_y**2 + norm_z**2
        normi  = 1.0d0 / sqrt(normsum)
        v12(1) = norm_x * normi * r
        v12(2) = norm_y * normi * r
        v12(3) = norm_z * normi * r

        norm_x = v2(1) + v3(1)
        norm_y = v2(2) + v3(2)
        norm_z = v2(3) + v3(3)
        normsum = norm_x**2 + norm_y**2 + norm_z**2
        normi  = 1.0d0 / sqrt(normsum)
        v23(1) = norm_x * normi * r
        v23(2) = norm_y * normi * r
        v23(3) = norm_z * normi * r

        norm_x = v3(1) + v1(1)
        norm_y = v3(2) + v1(2)
        norm_z = v3(3) + v1(3)
        normsum = norm_x**2 + norm_y**2 + norm_z**2
        normi  = 1.0d0 / sqrt(normsum)
        v31(1) = norm_x * normi * r
        v31(2) = norm_y * normi * r
        v31(3) = norm_z * normi * r

        exc12 = icosa_point_exclusion(x, vdwrad, &
                                      iatenvcnt, iatenv, &
                                      iatinclcnt12, iatincl12, &
                                      xi + v12(1), yi + v12(2), zi + v12(3))
        exc23 = icosa_point_exclusion(x, vdwrad, &
                                      iatenvcnt, iatenv, &
                                      iatinclcnt23, iatincl23, &
                                      xi + v23(1), yi + v23(2), zi + v23(3))
        exc31 = icosa_point_exclusion(x, vdwrad, &
                                      iatenvcnt, iatenv, &
                                      iatinclcnt31, iatincl31, &
                                      xi + v31(1), yi + v31(2), zi + v31(3))

        sas = 0.0d0
        ! ---   triangle 1
        if( exc1 .and. exc12 .and. exc31 ) then
          sas = sas + total_sas / 4.0d0
        else if( (nof_parts .lt. ipmin) .or. &
                 exc1 .or. exc12 .or. exc31 ) then
          sas = sas + icosa_patch_approx(i, x, vdwrad, v1, v12, v31, &
                                         exc1, exc12, exc31, &
                                         iatinclcnt1, iatincl1, &
                                         iatinclcnt12, iatincl12, &
                                         iatinclcnt31, iatincl31, &
                                         iatenvcnt, iatenv, &
                                         4 * nof_parts, &
                                         idecomp) 
        else if(idecomp.gt.0) then
          !write(6,*) 'TR1 ', i, iatinclcnt1, iatinclcnt12, iatinclcnt31
          call decsasa_icosa(i,iatinclcnt1,iatincl1, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt12,iatincl12, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt31,iatincl31, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
        end if
        ! ---   triangle 2
        if( exc2 .and. exc23 .and. exc12 ) then 
          sas = sas + total_sas / 4.0d0
        else if( (nof_parts .lt. ipmin) .or. &
                 exc2 .or. exc23 .or. exc12 ) then
          sas = sas + icosa_patch_approx(i, x, vdwrad, v2, v23, v12, &
                                         exc2, exc23, exc12, &
                                         iatinclcnt2, iatincl2, &
                                         iatinclcnt23, iatincl23, &
                                         iatinclcnt12, iatincl12, &
                                         iatenvcnt, iatenv, &
                                         4 * nof_parts, &
                                         idecomp)
        else if(idecomp.gt.0) then
          !write(6,*) 'TR2 ', i, iatinclcnt2, iatinclcnt23, iatinclcnt12
          call decsasa_icosa(i,iatinclcnt2,iatincl2, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt23,iatincl23, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt12,iatincl12, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
        end if
        ! ---   triangle 3
        if( exc3 .and. exc31 .and. exc23 ) then
          sas = sas + total_sas / 4.0d0
        else if( (nof_parts .lt. ipmin) .or. &
                 exc3 .or. exc31 .or. exc23 ) then
          sas = sas + icosa_patch_approx(i, x, vdwrad, v3, v31, v23, &
                                         exc3, exc31, exc23, &
                                         iatinclcnt3, iatincl3, &
                                         iatinclcnt31, iatincl31, &
                                         iatinclcnt23, iatincl23, &
                                         iatenvcnt, iatenv, &
                                         4 * nof_parts, &
                                         idecomp)
        else if(idecomp.gt.0) then
          !write(6,*) 'TR3 ', i, iatinclcnt3, iatinclcnt31, iatinclcnt23
          call decsasa_icosa(i,iatinclcnt3,iatincl3, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt31,iatincl31, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt23,iatincl23, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
        end if
        ! ---   triangle c
        if( exc12 .and. exc23 .and. exc31 ) then
          sas = sas + total_sas / 4.0d0
        else if( (nof_parts .lt. ipmin) .or. &
                 exc12 .or. exc23 .or. exc31 ) then
          sas = sas + icosa_patch_approx(i, x, vdwrad, v12, v23, v31, &
                                         exc12, exc23, exc31, &
                                         iatinclcnt12, iatincl12, &
                                         iatinclcnt23, iatincl23, &
                                         iatinclcnt31, iatincl31, &
                                         iatenvcnt, iatenv, &
                                         4 * nof_parts, &
                                         idecomp)
        else if(idecomp.gt.0) then
          !write(6,*) 'TRC ', i, iatinclcnt12, iatinclcnt23, iatinclcnt31
          call decsasa_icosa(i,iatinclcnt12,iatincl12, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt23,iatincl23, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
          call decsasa_icosa(i,iatinclcnt31,iatincl31, &
                             -total_sas/4.0d0/3.0d0, &
                             idecomp)
        end if

      end if

      ipapprox = sas
      
      return
end function icosa_patch_approx
    
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Approximates sas of a sphere
_REAL_ function icosa_sphere_approx(i, x, vdwrad, &
                                    ineighborpt, ineighbor, &
                                    idecomp)
      use constants, only : FOURPI
      use decomp, only : decsasa
      implicit none
      integer i, j, idecomp
      integer ineighborpt, ineighbor
      integer iatenvcnt, iatenv
      integer iatinclcnt, iatincl
      logical exc
      _REAL_ x, vdwrad, r
      _REAL_ total_sas, sas, sastmp
      _REAL_ v
      _REAL_ xi, yi, zi
      dimension x(*), vdwrad(*), ineighbor(*), exc(12), v(12*3)
      dimension iatenv(iemax)
      dimension iatinclcnt(12), iatincl(12*iemax)

      r = vdwrad(i)
      xi = x(3*i-2)   
      yi = x(3*i-1)   
      zi = x(3*i  )   
      total_sas = FOURPI*r*r / 20.0d0
      sas = 0.0d0

      call icosa_atom_environment(ineighborpt, ineighbor, &
                                  iatenvcnt, iatenv)

      if(iatenvcnt .gt. 0) then
        do j= 1, 12
          v(3*(j-1)+1) = icosa(j,1) * r
          v(3*(j-1)+2) = icosa(j,2) * r
          v(3*(j-1)+3) = icosa(j,3) * r
          exc(j) = icosa_point_exclusion(x, vdwrad, &
                                         iatenvcnt, iatenv, &
                                         iatinclcnt(j), iatincl((j-1)*iemax+1), &
                                         xi + v(3*(j-1)+1), yi + v(3*(j-1)+2), &
                                         zi + v(3*(j-1)+3)) 
        end do

        ! ---   Decomposition used here:
        !       First, the complete unburied SASA of an atom is given as self term.
        !       Subsequently, for each buried triangle, the amount of
        !       burial is added as negative SASA contribution.
        !       Hence, summing over self, indirect and direct term yields
        !       the correct total SASA of this atom. Furthermore, this allows
        !       to estimate to what extent atoms of neighboring residues contribute 
        !       to the burial.

        if(idecomp.eq.1 .or. idecomp.eq.2) then
          call decsasa(1,i,0,0,FOURPI*r*r)
        else if(idecomp.eq.3 .or. idecomp.eq.4) then
          call decsasa(-1,i,0,0,FOURPI*r*r)
        endif

        ! ---   Caveat: ictri values run from 0 .. 11

        do j= 1,20
          if( exc(ictri(j,1)+1) .and. exc(ictri(j,2)+1) .and. &
              exc(ictri(j,3)+1) ) then
            sas = sas + total_sas
          else
            sastmp = icosa_patch_approx (i, x, vdwrad, &
                           v(3*ictri(j,1)+1), v(3*ictri(j,2)+1), &
                           v(3*ictri(j,3)+1), &
                           exc(ictri(j,1)+1), exc(ictri(j,2)+1), &
                           exc(ictri(j,3)+1), &
                           iatinclcnt(ictri(j,1)+1), &
                           iatincl(ictri(j,1)*iemax+1), &
                           iatinclcnt(ictri(j,2)+1), &
                           iatincl(ictri(j,2)*iemax+1), &
                           iatinclcnt(ictri(j,3)+1), &
                           iatincl(ictri(j,3)*iemax+1), &
                           iatenvcnt, iatenv, 20, &
                           idecomp)
            sas = sas + sastmp
          end if
        end do
      end if

      icosa_sphere_approx = sas

      return
end function icosa_sphere_approx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Generates rotation matrix for rotation by angle around axis.
subroutine gen_rot_mat(angle, axis, rotmat)

      ! --- The angle must be given in degree, the axis is normalized here.
      !     Rotmat is a (3,3) matrix.
      use constants, only : DEG_TO_RAD
      implicit none
      _REAL_ axis, angle, rotmat
      _REAL_ norm
      _REAL_ c, s, t
      dimension axis(*), rotmat(*)


      angle = angle * DEG_TO_RAD
      c = cos(angle)
      s = sin(angle)
      t = 1 - c
      norm = 1.0d0 / sqrt(axis(1)**2 + axis(2)**2 + axis(3)**2)
      axis(1) = axis(1) * norm
      axis(2) = axis(2) * norm
      axis(3) = axis(3) * norm

      rotmat(3*1-2) = t * axis(1)**2 + c
      rotmat(3*1-1) = t * axis(1) * axis(2) - s * axis(3)
      rotmat(3*1  ) = t * axis(1) * axis(3) + s * axis(2)

      rotmat(3*2-2) = t * axis(1) * axis(2) + s * axis(3)
      rotmat(3*2-1) = t * axis(2)**2 + c
      rotmat(3*2  ) = t * axis(2) * axis(3) - s * axis(1)

      rotmat(3*3-2) = t * axis(1) * axis(3) - s * axis(2)
      rotmat(3*3-1) = t * axis(2) * axis(3) + s * axis(1)
      rotmat(3*3  ) = t * axis(3)**2 + c

      return
end  subroutine gen_rot_mat

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Rotates point according to rotmat.
subroutine rot_point(point, rotmat)

      implicit none
      integer i, j
      _REAL_ point, rotmat
      _REAL_ ptmp
      dimension point(*), rotmat(*), ptmp(3)

      do i=1,3
        ptmp(i) = 0.0d0
        do j=1,3
          ptmp(i) = ptmp(i) + rotmat(3*(i-1)+j) * point(j)
        end do
      end do

      do i=1,3
        point(i) = ptmp(i)
      end do

      return
end subroutine rot_point

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module icosasurf
