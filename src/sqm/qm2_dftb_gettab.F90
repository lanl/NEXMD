! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"


! ========================================
! Subroutine Gettab
! -----------------
!
! Reads the Slater-Kosher tables from the
! files indicated. In Amber, the SK files
! must be stored in:
!         $(AMBERHOME)/dat/slko
! ========================================
subroutine gettab(ntype)

  use qm2_dftb_module, only: MAXINT, sktab, espin, mcharge, spltab

  implicit none

!! Passed in
  integer, intent(in) :: ntype

!! Locals
  character(1024) :: skfile
  character(128) :: chdummy
  _REAL_  :: qzeroh(3),uhubbh(3)
  _REAL_  :: at_mass, polyn_coeff(8), polyn_cutoff, bondmap_cutoff
  _REAL_  :: dp_dip_int, ps_dip_int, extra_info(7), real_junk(20)
  integer :: ppp, i, j, k, l

  do i = 1,ntype
     mcharge%qzero(i) = 0.0d0
     mcharge%uhubb(i) = 0.0d0
     do j = 1,ntype
        skfile = sktab%skfiles(i,j)

        ! Opens Slater-Koster file
        open (3,file=skfile,status='unknown')
        rewind 3

        !------------------
        ! Read file header
        !------------------

        ! --> sr: step width
        ! --> dimens: number of distances
        read (3,*) sktab%sr(i,j),sktab%dimens(i,j)

        ! For homonuclear bonds only:
        if (i == j) then

           read (3,*) (sktab%skself(l,i),l = 1,3),& ! --> skself: d, p and s (self) energies.     
                      espin(i), &                   ! --> espin:  atomic spin-polarization energy 
                      (uhubbh(4-l), l=1,3), &       ! --> uhubbh: d, p, s hubbard parameters      
                      (qzeroh(4-l), l=1,3)          ! --> qzeroh: d, p, s atomic charges          

           mcharge%uhubb(i) = uhubbh(1)

           ! Atomic charge (qzero) is the sum of
           ! the s, p and d charges [qzeroh(1--3)].
           do k = 1,3
              mcharge%qzero(i) = mcharge%qzero(i) + qzeroh(k)
           ENDDO

           ! Atomic Spin Polarization Energy:
           do l=1,3
              espin(i) = espin(i) + sktab%skself(l,i)*qzeroh(4-l)
           ENDDO

        endif  !! if(i == j)
        ! End reading header

        ! -----------------
        ! Intermediary part
        ! -----------------

        ! Every one of the S-K files contains, after the mentioned data, a list with
        ! 20 lines supposed to contain:
        ! --> (1) atomic mass
        ! --> (8) 8 polynomial coefficients
        ! --> (1) cutoff radius of polynomial
        ! --> (1) cutoff radius of bondmap
        ! --> (1) d, p dipole integral
        ! --> (1) p, s dipole integral
        ! --> (7) 7 "unused" positions.
        ! (A total of 20 values)
        !
        ! --> plus more 19 lines of the type:
        !     20*1.0  or
        !     20*0.0
        !
        ! For some reason, those 20 lines are being read into H and S, so that 
        ! when the first 20 distances in H ans S are *NOT* what they were supposed to 
        ! be from the S-K tables!!
        !
        ! 
     


        ! ------------------------
        ! Reads H and S parameters
        ! ------------------------

        do k = 1,sktab%dimens(i,j)
           read (3,*) (sktab%skhtab(l,k,i,j),l = 1,10), & ! 10 H and
                      (sktab%skstab(l,k,i,j),l = 1,10)    ! 20 S values for each distance.
        ENDDO



        ! ---------------------------------------
        ! Reads repulsive potential 'Spline' data
        ! ---------------------------------------

        ! The repulsive potential is written in the form of
        ! a number of cubic splines, that each span a small
        ! interval of the distance. At each interval, the
        ! corresponding spline can be writen as:
        ! Y(t) = a + b.t + c.t^2 + d.t^3
        !
        ! The LAST interval is actually fit to a 5th order spline:
        ! Y(t) = a + b.t + c.t^2 + d.t^3 + e.t^4 + f.t^5

        ! First, locate it's beggining. It is
        ! marked by the 'Spline' keyword.
23016   if( .TRUE. )then
           read(3,'(A)') chdummy
           if(chdummy == 'Spline')then
              goto 23017
           endif
           goto 23016
        endif
23017   continue


        ! The first line after 'Spline':
        ! --> numint: number of intervals
        ! --> cutoff distance
        read(3,*) spltab%numint(i,j),spltab%cutoff(i,j)
        if(spltab%numint(i,j) > MAXINT)then
           write(6,*) 'Too many intervalls!'
           goto 99
        endif

        ! --> efkt: Exponential function parameters (small distances)
        read(3,*) (spltab%efkt(ppp,i,j),ppp=1,3)

        ! For each interval, read the spline data:
        ! --> xr(1,ppp,i,j)        : start distance of interval
        ! --> xr(2,ppp,i,j)        : end   distance of interval
        ! --> coeff(1--4, ppp,i,j) : Coefficients for cubic spline
        do  ppp=1,spltab%numint(i,j)
           if(ppp < spltab%numint(i,j))then
              read (3,*) spltab%xr(1,ppp,i,j),   spltab%xr(2,ppp,i,j),    &
                         spltab%coeff(1,ppp,i,j),spltab%coeff(2,ppp,i,j), &
                         spltab%coeff(3,ppp,i,j),spltab%coeff(4,ppp,i,j)
           else
              ! The last line has 2 extra numbers because the last
              ! spline is of 5th order.
              read (3,*) spltab%xr(1,ppp,i,j),   spltab%xr(2,ppp,i,j),    &
                         spltab%coeff(1,ppp,i,j),spltab%coeff(2,ppp,i,j), &
                         spltab%coeff(3,ppp,i,j),spltab%coeff(4,ppp,i,j), &
                         spltab%coeff(5,ppp,i,j),spltab%coeff(6,ppp,i,j)
           endif
        ENDDO ! do ppp=1,spltab%numint(i,j)

        ! Checks that the end of the last interval coincides with the 
        ! cutoff. Otherwise, the file must be corrupted.
        if(spltab%xr(2,spltab%numint(i,j),i,j) /= spltab%cutoff(i,j))then
           write(6,*) 'Error in Datafile'
           goto 99
        endif

        ! OK, we are done here. Time for doughnuts.
        close (3)
     ENDDO ! do j = 1,ntype
  ENDDO ! do i = 1,ntype

99 continue
end subroutine gettab
