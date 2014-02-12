! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#include "pb_def.h"
#include "timer.h"

module solvent_accessibility

   implicit none

   integer               :: radiopt       ! radi option
   integer               :: maxsph        ! max no of surface dots per atom
   integer               :: maxarc        ! max no of arcs per atom
   integer               :: maxarcdot     ! max no of arc dots per atom
   integer               :: nsrfdot       ! no of sa surface dots
   integer               :: narcdot       ! no of sa arc dots
   integer               :: narc          ! no of sa arcs
   integer               :: triopt 
   _REAL_              :: radiscale     ! scale factor of radii
   _REAL_              :: sprob         ! solvent probe for dispersion surface
   _REAL_              :: vprob         ! solvent probe for cavity volume
   _REAL_              :: dprob         ! solvent probe for dielectric surface
   _REAL_              :: iprob         ! ion probe for stern surface
   _REAL_              :: arcres        ! arc dot resolution, with respect to the grid spacing
   _REAL_              :: prtsas        ! total sa surface area for printout
   _REAL_              :: prtsav        ! total volume within sa surface area for printout

   integer , allocatable ::  nzratm(:)    ! atomic index of atoms with nonzero radii
   _REAL_, allocatable ::   mdsig(:)    ! atomic radii (np vdw radii)
   _REAL_, allocatable ::    rmin(:)    ! atomic radii (np vdw radii)
   _REAL_, allocatable ::    radi(:)    ! atomic radii (pb cavity radii)
   _REAL_, allocatable ::   radip(:)    ! atomic radii + dprob/sprob
   _REAL_, allocatable ::  radip2(:)    ! squared radip
   _REAL_, allocatable ::    scrd(:,:)  ! coordinates of maxsph evenly distributed dots on the unit sphere

   integer , allocatable ::  fstsdot(:)   ! first surface dot of an atom
   integer , allocatable ::  lstsdot(:)   ! last surface dot of an atom
   integer , allocatable ::  fstadot(:)   ! first arc dot of an atom
   integer , allocatable ::  lstadot(:)   ! last arc dot of an atom
   integer , allocatable ::   fstarc(:)   ! first arc of an atom
   integer , allocatable ::   lstarc(:)   ! last arc of an atom
   integer , allocatable ::     marc(:)   ! no. of arcs of an atom
   integer , allocatable ::   m2narc(:,:) ! marc to narc index for an atom
   integer , allocatable ::   arcatm(:,:) ! two atoms that generate an arc
   integer , allocatable ::   dotarc(:)   ! arc index of an arc dot
   _REAL_, allocatable ::   savarc(:,:) ! saved arc geometries
   _REAL_, allocatable ::  savactr(:,:) ! saved arc center coordinates
   _REAL_, allocatable ::   srfcrd(:,:) ! surface dot coordinates
   _REAL_, allocatable ::   arccrd(:,:) ! arc dot coordinates

   integer               ::   ntri
   integer               ::   maxtri
   integer , allocatable ::   triatm(:,:) ! 
   integer , allocatable ::  triatm1(:,:) ! 
   integer , allocatable ::   triarc(:,:) ! 
   integer , allocatable ::  triarc1(:,:) ! 
   _REAL_  , allocatable ::   tricrd(:,:) ! 
   _REAL_  , allocatable ::  tricrd1(:,:) ! 
   logical , allocatable :: knocktri(:)   ! 

   logical , allocatable :: knockout(:)   ! auxiliary flags for checking surface/arcs
   integer , allocatable ::   spharc(:)   ! auxiliary arc index
   integer , allocatable ::  spharc1(:)   ! auxiliary arc index
   _REAL_, allocatable ::   sphcrd(:,:) ! auxiliary surface/arc dot coordinates
   _REAL_, allocatable ::  sphcrd1(:,:) ! auxiliary surface/arc dot coordinates

   ! PBMD dynamics radii setting

   _REAL_              :: radinc
   _REAL_              :: expthresh
   _REAL_, allocatable ::  radip3(:)
   _REAL_, allocatable ::    nmax(:)
   _REAL_, allocatable ::    nexp(:)
   _REAL_, allocatable :: sumnmax(:)
   _REAL_, allocatable :: sumnexp(:)
   _REAL_, allocatable ::  avnmax(:)
   _REAL_, allocatable ::  avnexp(:)

   ! Ligand Mask / Traditional Focusing / Multiple-block Focusing
   logical done_sa_driver
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Cleanup routine of solvent accessible surface and arcs
subroutine sa_free ( dosas,ndosas,npsa )

   integer dosas, ndosas, alloc_err(16)
   logical npsa

   alloc_err = 0
   if ( dosas .eq. ndosas) then
   deallocate( fstsdot, stat = alloc_err(1 ) )
   deallocate( lstsdot, stat = alloc_err(2 ) )
   deallocate(  srfcrd, stat = alloc_err(3 ) )
   if ( .not. npsa ) then
      deallocate( fstadot, stat = alloc_err(4 ) )
      deallocate( lstadot, stat = alloc_err(5 ) )
      deallocate(  fstarc, stat = alloc_err(6 ) ) 
      deallocate(  lstarc, stat = alloc_err(7 ) )
      deallocate(    marc, stat = alloc_err(8 ) )
      deallocate(  m2narc, stat = alloc_err(9 ) )
      deallocate(  arcatm, stat = alloc_err(10) )
      deallocate(  savarc, stat = alloc_err(11) )
      deallocate( savactr, stat = alloc_err(12) )
      deallocate(  dotarc, stat = alloc_err(13) )
      deallocate(  arccrd, stat = alloc_err(14) )
      deallocate(  triatm, stat = alloc_err(15) )
      deallocate(  triarc, stat = alloc_err(16) )
   end if

   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err( 5)+&
        alloc_err(6)+alloc_err(7)+alloc_err(8)+alloc_err(9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+&
        alloc_err(15)+alloc_err(16) /= 0 ) then
      write(6, *) 'SA Bomb in sa_free(): DeAllocation aborted', alloc_err(1:16)
      call mexit(6, 1)
   end if
   end if

end subroutine sa_free
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Cleanup routine of solvent accessible surface and arcs
subroutine sa_free_mb ( dosas,ndosas )

   integer dosas, ndosas, alloc_err(21)

   alloc_err(1:21) = 0
   if ( dosas .eq. ndosas) then
   !deallocate( fstsdot, stat = alloc_err(1 ) )
   !deallocate( lstsdot, stat = alloc_err(2 ) )
   !deallocate(  srfcrd, stat = alloc_err(3 ) )
   deallocate( fstadot, stat = alloc_err(4 ) )
   deallocate( lstadot, stat = alloc_err(5 ) )
   deallocate(  fstarc, stat = alloc_err(6 ) ) 
   deallocate(  lstarc, stat = alloc_err(7 ) )
   deallocate(    marc, stat = alloc_err(8 ) )
   deallocate(  m2narc, stat = alloc_err(9 ) )
   deallocate(  arcatm, stat = alloc_err(10) )
   deallocate(  savarc, stat = alloc_err(11) )
   deallocate( savactr, stat = alloc_err(12) )
   deallocate(  dotarc, stat = alloc_err(13) )
   deallocate(  arccrd, stat = alloc_err(14) )
   deallocate(  triatm, stat = alloc_err(15) )
!  deallocate( triatm1, stat = alloc_err(16) )
   deallocate(  triarc, stat = alloc_err(17) )
!  deallocate( triarc1, stat = alloc_err(18) )
!  deallocate(  tricrd, stat = alloc_err(19) )
!  deallocate( tricrd1, stat = alloc_err(20) )
!  deallocate(knocktri, stat = alloc_err(21) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err( 5)+&
        alloc_err(6)+alloc_err(7)+alloc_err(8)+alloc_err(9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+&
        alloc_err(15)+alloc_err(16)+alloc_err(17)+alloc_err(18)+&
        alloc_err(19)+alloc_err(20)+alloc_err(21) /= 0 ) then
      write(6, *) 'SA Bomb in sa_free_mb(): DeAllocation aborted', alloc_err(1:21)
      call mexit(6, 1)
   end if
   end if

end subroutine sa_free_mb
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ NP dispersion/cavity radi setup based on nonelectrostatic solvation free energies in TIP3P
subroutine sa_init( verbose, pbprint, natom, atmlast, ifcap, prob, r, rp, rp2, outflag )

   implicit none

#  include "pb_constants.h"

   ! Passed variables

   logical verbose, pbprint
   integer natom, atmlast, ifcap
   integer outflag(*)
   !character (len=4) :: isymbl(natom)
   _REAL_ prob, r(natom), rp(natom), rp2(natom)

   ! Local variables
    
   integer iatm
    
   ! for InsightII surface display only

   !do iatm = 1, natom
   !   else if ( isymbl(iatm)(1:1) == 'C' .or. isymbl(iatm)(1:1) == 'c' ) then
   !      r(iatm) = 1.55d0
   !   else if ( isymbl(iatm)(1:1) == 'H' .or. isymbl(iatm)(1:1) == 'h' ) then
   !      r(iatm) = 1.10d0
   !   else if ( isymbl(iatm)(1:1) == 'N' .or. isymbl(iatm)(1:1) == 'n' ) then
   !      r(iatm) = 1.40d0
   !   else if ( isymbl(iatm)(1:1) == 'O' .or. isymbl(iatm)(1:1) == 'o' ) then
   !      r(iatm) = 1.35d0
   !   else if ( isymbl(iatm)(1:1) == 'P' .or. isymbl(iatm)(1:1) == 'p' ) then
   !      r(iatm) = 1.88d0
   !   else if ( isymbl(iatm)(1:1) == 'S' .or. isymbl(iatm)(1:1) == 's' ) then
   !      r(iatm) = 1.81d0
   !   else
   !      write(6, *) 'SA Bomb in sa_init(): No radius assigned for atom', iatm, isymbl(iatm)
   !      call mexit(6, 1)
   !   end if
   !end do
   
!  if ( verbose .and. pbprint ) write(6, *) ' SA surface: setting up working radii'
   if ( verbose .and. pbprint ) then 
      write(6,*)
      write(6,*) '======== Setting up Solvent Accessibility Data ========'
      write(6,*) 'Setting up working radii'
   end if
   rp  = ZERO
   rp2 = ZERO
   do iatm = 1, atmlast
      if(ifcap == 5 .and. outflag(iatm) == 1) cycle
      if ( r(iatm) /= ZERO ) then
         rp(iatm) = r(iatm) + prob
         rp2(iatm) = rp(iatm)**2
      endif
   end do

   ! Ligand Mask / Traditional Focusing / Multiple-block Focusing
   done_sa_driver = .false.

end subroutine sa_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ VDW radi setup for Pymol surface potential using Pymol vdw radii
subroutine phi_aaradi( natom, isymbl, radi )

   implicit none

   ! Passed variables

   integer natom
   character (len=4) :: isymbl(natom)
   _REAL_ radi(natom)

   ! Local varialbes

   integer iatm

   do iatm = 1, natom
      if      (isymbl(iatm)(1:2) == 'Br' .or. isymbl(iatm)(1:2) == 'br') then
         radi(iatm) = 1.85d0
      else if ( isymbl(iatm)(1:2) == 'C0' ) then
         radi(iatm) = 1.80d0
      else if (isymbl(iatm)(1:2) == 'Cl' .or. isymbl(iatm)(1:2) == 'cl') then
         radi(iatm) = 1.75d0
      else if ( isymbl(iatm)(1:2) == 'CU' ) then
         radi(iatm) = 1.40d0
      else if ( isymbl(iatm)(1:2) == 'Cs' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'C' .or. isymbl(iatm)(1:1) == 'c' ) then
         radi(iatm) = 1.70d0
      else if ( isymbl(iatm)(1:2) == 'FE' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'F' .or. isymbl(iatm)(1:1) == 'f' ) then
         radi(iatm) = 1.47d0
      else if ( isymbl(iatm)(1:1) == 'H' .or. isymbl(iatm)(1:1) == 'h' ) then
         radi(iatm) = 1.20d0
      else if ( isymbl(iatm)(1:2) == 'IB' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:2) == 'IM' ) then
         radi(iatm) = 1.75d0
      else if ( isymbl(iatm)(1:2) == 'IP' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:1) == 'I' .or. isymbl(iatm)(1:1) == 'i' ) then
         radi(iatm) = 1.98d0
      else if ( isymbl(iatm)(1:2) == 'K ' ) then
         radi(iatm) = 2.75d0
      else if ( isymbl(iatm)(1:2) == 'Li' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'LP' ) then
         radi(iatm) = 1.20d0
      else if ( isymbl(iatm)(1:2) == 'Mg' ) then
         radi(iatm) = 1.73d0
      else if ( isymbl(iatm)(1:2) == 'Mn' ) then
         radi(iatm) = 1.73d0
      else if ( isymbl(iatm)(1:2) == 'Na' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:1) == 'N' .or. isymbl(iatm)(1:1) == 'n' ) then
         radi(iatm) = 1.55d0
      else if ( isymbl(iatm)(1:1) == 'O' .or. isymbl(iatm)(1:1) == 'o' ) then
         radi(iatm) = 1.52d0
      else if ( isymbl(iatm)(1:1) == 'P' .or. isymbl(iatm)(1:1) == 'p' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'Rb' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'S' .or. isymbl(iatm)(1:1) == 's' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'Zn' ) then
         radi(iatm) = 1.39d0
      else
         write(6, *) 'SA Bomb in phi_aaradi(): No radius assigned for atom', iatm, isymbl(iatm)
         call mexit(6, 1)
      end if
   end do

   radi = radi -0.10d0

end subroutine phi_aaradi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PB cavity radi setup based on electrostatic solvation free energies in TIP3P
subroutine pb_aaradi( natom, nbonh, ibh, jbh, radi, acrg, ucrgh, ucrga, resid, igraph, isymbl, rin )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Chuck Tan, Luo Research Group, UC-Irvine
   !
   ! The assigned cavity radii were opitimzed from reproducing TIP3P electrostatic
   ! solvation free energies of all Amber standard residues/fragments (42) and single
   ! ions in the database files.
   !
   ! However, small organic molecules generated by Antechamber use parm file radii!
   ! 
   ! Dielectric assignment use a solvent probe of 1.6 A and solvent excluded surface.
   ! fdpb used a grid spacing of 0.5 A with the harmonic smoothing method to
   ! define the dielectric boundary. fdpb energy was averaged over 27 different
   ! grid positions. See Tan, Yang and Luo, In Prep. 2006 for more information.
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "pb_constants.h"

   ! Passed variables

   integer natom, nbonh
   integer ibh(*), jbh(*)
   character (len=4) :: resid(*), igraph(*), isymbl(*)
   _REAL_ radi(*), acrg(*), ucrgh(*), ucrga(*), rin(*)

   ! Local variables
    
   integer nhn(natom)
   integer iatm, jatm, idum
    
   ! accumulate no. of H attached to N

   nhn(1:natom) = 0

   do idum = 1, nbonh

      iatm = ibh(idum)/3 + 1
      jatm = jbh(idum)/3 + 1

      if ( isymbl(iatm)(1:1) == 'N' )  then
         nhn(iatm) = nhn(iatm) + 1
      endif

      if ( isymbl(jatm)(1:1) == 'N' )  then
         nhn(jatm) = nhn(jatm) + 1
      endif

   enddo

   ! ::::: radii are zero by default :::::
   ! this takes care of All zero H's
    
   radi(1:natom) = ZERO
    
   ! ::::: radii for heavy atoms without H :::::

   do iatm = 1, natom

      if     ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! ----- radii of CT without H -----
      if     ( isymbl(iatm)(1:2) == 'CT' ) then
         radi(iatm) = 2.00d0
      ! ----- radii of CX without H -----
      elseif ( isymbl(iatm)(1:2) == 'CX' ) then
         radi(iatm) = 2.00d0
      ! ----- radii of aromatic C without H -----
      elseif ( isymbl(iatm)(1:2) == 'CA' .or. &
               isymbl(iatm)(1:2) == 'CB' .or. &
               isymbl(iatm)(1:2) == 'CC' .or. &
               isymbl(iatm)(1:2) == 'CK' .or. &
               isymbl(iatm)(1:2) == 'CM' .or. &
               isymbl(iatm)(1:2) == 'CN' .or. &
               isymbl(iatm)(1:2) == 'CQ' .or. &
               isymbl(iatm)(1:2) == 'CR' .or. &
               isymbl(iatm)(1:2) == 'CV' .or. &
               isymbl(iatm)(1:2) == 'CW' .or. &
               isymbl(iatm)(1:2) == 'C*' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 1.85d0
         else
            radi(iatm) = 1.75d0
         endif
      ! ----- radii of C  -----
      elseif ( isymbl(iatm)(1:2) == 'C ' ) then
         if ( ucrga(iatm) < -0.5d0 ) then
            radi(iatm) = 1.65d0
         else
            radi(iatm) = 2.00d0
         endif

      ! ----- radii of NB, NC, and N* without H -----
      ! two kinds of charge environments
      elseif ( isymbl(iatm)(1:2) == 'NB' .or. &
               isymbl(iatm)(1:2) == 'NC' .or. &
               isymbl(iatm)(1:2) == 'N*' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 1.80d0
         else
            radi(iatm) = 1.70d0
         endif
      ! ----- radii of N without H -----
      elseif ( isymbl(iatm)(1:2) == 'N ' ) then
         radi(iatm) = 1.70d0

      ! ----- radii of O, O2, OS without H -----
      elseif ( isymbl(iatm)(1:2) == 'O ' ) then
         radi(iatm) = 1.57d0
      elseif ( isymbl(iatm)(1:2) == 'O2' ) then
         radi(iatm) = 1.25d0
      elseif ( isymbl(iatm)(1:2) == 'OS' ) then
         radi(iatm) = 1.30d0

      ! ----- radii of S atoms without H -----
      elseif ( isymbl(iatm)(1:1) == 'S' ) then
         if ( ucrga(iatm) < -0.5d0 ) then
            radi(iatm) = 1.45d0
         else
            radi(iatm) = 2.00d0
         endif

      ! ----- radii of P atoms without H -----
      elseif ( isymbl(iatm)(1:1) == 'P' ) then
         radi(iatm) = 1.95d0

      ! ----- radii of Ions -----
      elseif ( isymbl(iatm)(1:2) == 'Li' ) then
         radi(iatm) = 1.481d0
      elseif ( isymbl(iatm)(1:2) == 'Na' ) then
         radi(iatm) = 1.875d0
      elseif ( isymbl(iatm)(1:2) == 'IP' ) then
         radi(iatm) = 1.875d0
      elseif ( isymbl(iatm)(1:2) == 'K ' ) then
         radi(iatm) = 2.288d0
      elseif ( isymbl(iatm)(1:2) == 'Rb' ) then
         radi(iatm) = 2.448d0
      elseif ( isymbl(iatm)(1:2) == 'Cs' ) then
         radi(iatm) = 2.712d0
      elseif ( isymbl(iatm)(1:2) == 'F ' ) then
         radi(iatm) = 1.223d0
      elseif ( isymbl(iatm)(1:2) == 'Cl' ) then
         radi(iatm) = 1.516d0
      elseif ( isymbl(iatm)(1:2) == 'IM' ) then
         radi(iatm) = 1.815d0
      elseif ( isymbl(iatm)(1:2) == 'Br' ) then
         radi(iatm) = 1.745d0
      elseif ( isymbl(iatm)(1:2) == 'I ' ) then
         radi(iatm) = 1.870d0
      elseif ( isymbl(iatm)(1:2) == 'MG' ) then
         radi(iatm) = 1.515d0
      elseif ( isymbl(iatm)(1:2) == 'C0' ) then
         radi(iatm) = 2.126d0
      elseif ( isymbl(iatm)(1:2) == 'Zn' ) then
         radi(iatm) = 1.469d0

      endif

   end do

   ! ::::: radii for heavy atoms with H :::::

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! ----- radii of CT/HC_n, CT/H1_n, and CT/HP_n -----

      if     ( isymbl(iatm)(1:2) == 'CT' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) then
            if ( ucrga(iatm) > 0.9d0 ) then      ! LYS
               radi(iatm) = 2.50d0
            elseif ( ucrga(iatm) < -0.9d0 ) then ! CYM
               radi(iatm) = 2.10d0
            else
               radi(iatm) = 2.25d0
            endif
         endif

      ! ----- radii of CX/HC_n, CX/H1_n, and CX/HP_n -----

      elseif ( isymbl(iatm)(1:2) == 'CX' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) then
            if ( ucrga(iatm) > 0.9d0 ) then      ! LYS
               radi(iatm) = 2.50d0
            elseif ( ucrga(iatm) < -0.9d0 ) then ! CYM
               radi(iatm) = 2.10d0
            else
               radi(iatm) = 2.25d0
            endif
         endif

     ! ----- radii of aromatic C/H -----

      elseif ( isymbl(iatm)(1:2) == 'CA' .or. &
               isymbl(iatm)(1:2) == 'CB' .or. &
               isymbl(iatm)(1:2) == 'CC' .or. &
               isymbl(iatm)(1:2) == 'CK' .or. &
               isymbl(iatm)(1:2) == 'CM' .or. &
               isymbl(iatm)(1:2) == 'CN' .or. &
               isymbl(iatm)(1:2) == 'CQ' .or. &
               isymbl(iatm)(1:2) == 'CR' .or. &
               isymbl(iatm)(1:2) == 'CV' .or. &
               isymbl(iatm)(1:2) == 'CW' .or. &
               isymbl(iatm)(1:2) == 'C*' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) radi(iatm) = 2.15d0

      ! ----- radii of NA always with H -----

      elseif ( isymbl(iatm)(1:2) == 'NA' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 2.15d0
         else
            radi(iatm) = 2.10d0
         endif

      ! ----- radii of N, N2, and N3 with H -----

      elseif ( nhn(iatm) /= 0 .and. &
             ( isymbl(iatm)(1:2) == 'N ' .or. &
               isymbl(iatm)(1:2) == 'N2' .or. &
               isymbl(iatm)(1:2) == 'N3' ) ) then

         if ( ucrga(iatm) < 0.00d0 ) then
            radi(iatm) = 1.80d0 ! LYN only
         elseif ( ucrga(iatm) > 0.80d0 ) then
            radi(iatm) = 2.50d0 ! LYS and ARG only
         else
            if     ( nhn(iatm) == 2 .or. nhn(iatm) == 3 ) then
               radi(iatm) = 2.10d0
            elseif ( nhn(iatm) == 1 ) then
               radi(iatm) = 1.95d0
            endif
         endif

     ! ----- radii of OH always with H -----
     ! OH on rings, ASH/GLH behave like OH on rings

      elseif ( isymbl(iatm)(1:2) == 'OH' ) then
         radi(iatm) = 1.925d0

      ! ----- radii of OW always with H -----
      ! in TIP3P water, TIP4P/TIP5P not checked

      elseif ( isymbl(iatm)(1:2) == 'OW' ) then
         radi(iatm) = 1.85d0

      ! ----- radii of SH with H -----

      elseif ( isymbl(iatm)(1:2) == 'SH' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) radi(iatm) = 2.30d0

      endif

   enddo

   ! exceptions:

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! Charged aromatic C/N: HIP only

      if ( resid(iatm)(1:3) == 'HIP' ) then
         if ( isymbl(iatm)(1:2) == 'CW' .or. &
              isymbl(iatm)(1:2) == 'CR' .or. &
              isymbl(iatm)(1:2) == 'NA' ) radi(iatm) = 2.65d0
      endif

      ! this is not exception
      ! OH not on rings, i.e. OH on linear alphatic CT

      if ( isymbl(iatm)(1:2) == 'OH' ) then
         if ( resid(iatm)(1:3) == 'SER' .or. &
              resid(iatm)(1:3) == 'THR' .or. &
              igraph(iatm)(1:2) == 'O5' ) radi(iatm) = 1.79d0
      endif

      ! this is not exception
      ! OS not on rings, they behave like O2

      if ( isymbl(iatm)(1:2) == 'OS' ) then
         if ( igraph(iatm)(1:2) == 'O4' ) radi(iatm) = 1.25d0 
      endif

   enddo

   ! for Antechamber ligands with atom types of lower cases
   ! use radii from the parm file

   do iatm = 1, natom
      if ( isymbl(iatm)(1:1) == 'h' .or. &
           isymbl(iatm)(1:1) == 'c' .or. &
           isymbl(iatm)(1:1) == 'n' .or. &
           isymbl(iatm)(1:1) == 'o' .or. &
           isymbl(iatm)(1:1) == 'p' .or. &
           isymbl(iatm)(1:1) == 'f' .or. &
           isymbl(iatm)(1:1) == 'b' .or. &
           isymbl(iatm)(1:1) == 'i' .or. &
           isymbl(iatm)(1:1) == 's' ) radi(iatm) = rin(iatm)
   end do

   ! safe-guarding, very important!

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      if ( radi(iatm) == ZERO ) then
         write(6, *) 'PB Bomb in pb_aaradi(): No radius assigned for atom', iatm, igraph(iatm), isymbl(iatm)
         call mexit(6, 1)
      endif

   enddo

end subroutine pb_aaradi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize solvent exposed dots of a unit sphere
subroutine sa_sphere(maxsph, scrd)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Generates a set of about maxsph roughly evenly spaced points on the surface of a
   ! unit sphere, for use in generating probe-accessible surface. Number of points
   ! generated is maxsph, approximately equal to requested number. Modified from UHBD
   ! (Comp. Phys. Comm. 91:57-95, 1995) routine sphere() by Michael Gilson.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "pb_constants.h"

   ! Passed variables

   integer maxsph
   _REAL_ scrd(3, *)

   ! Local variables

   integer ntheta,npsi,npsimax,i,nt,np
   _REAL_ thtstp,theta,stheta,ctheta,psistp,psi,cpsi,spsi

   ! begin code

   ntheta = int(sqrt(PI*maxsph/FOUR))
   npsimax = int(TWO*ntheta)
   thtstp = PI/ntheta

   i = 1
   do nt = 1, ntheta
      theta = thtstp*nt
      stheta = sin(theta)
      ctheta = cos(theta)
      npsi = nint(stheta*npsimax)
      if (npsi == 0) cycle
      psistp = TWOPI/npsi
      do np = 1, npsi
         psi = np*psistp
         cpsi = cos(psi)
         spsi = sin(psi)
         scrd(1,i) = cpsi * stheta
         scrd(2,i) = spsi * stheta
         scrd(3,i) = ctheta
         i = i + 1
      enddo
   enddo
   maxsph = i - 1

end subroutine sa_sphere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of solvent accessible surface and arcs
subroutine sa_driver( verbose,pbprint,ipb,inp,natom,atmlast,dosas,ndosas,&
              npbstep,nsaslag,ligand,outflag,acrd,iar1pb,iprshrt,nex,iex,npsa )
 
   use pbtimer_module
   implicit none
#  include "pb_constants.h"

   ! Passed variables
    
   logical verbose, pbprint, ligand
   integer ipb, inp, natom, atmlast, dosas, ndosas, npbstep, nsaslag
   integer iar1pb(4,0:natom), iprshrt(*), nex(natom), iex(64,natom)
   integer outflag(natom)
   _REAL_ acrd(3,1:natom)
   logical npsa
    
   ! Local variables
    
   integer iatm
   _REAL_ wf0, wf1, exposure, increase, smallsas
 
   prtsas = ZERO
   prtsav = ZERO

   ! compute sa every "ndosas" FDPB calculations

   if ( (.not. ligand) .and. done_sa_driver ) then
      write(6,*) "sa_driver() is not allowed to be called twice without ligand/mb."
      call mexit(6,1)
   end if
   if ( dosas .eq. ndosas ) then
      dosas = 1
!  DO THIS FOR THE FIRST TIME SA_DRIVER GOT CALLED OR WHEN NO LIGAND/MULTIBLOCK SPECIFIED
      if ( .not. done_sa_driver ) then
          call pbtimer_start(PBTIME_PBSASRF)
          call sa_srf( verbose, pbprint, natom, 1, atmlast )
          call pbtimer_stop(PBTIME_PBSASRF)
      end if
!  DO THIS AFTER THE FIRST TIME SA_DRIVER GOT CALLED OR WHEN NO LIGAND/MULTIBLOCK SPECIFIED
      if ( .not. ligand .or. done_sa_driver ) then

         call pbtimer_start(PBTIME_PBSAARC)
         if ( .not. npsa .and. (ipb == 1 .or. ipb == 2 .or. ipb == 4 .or. ipb == 5) ) call sa_arc( verbose,pbprint,natom,1,atmlast )
         call pbtimer_stop(PBTIME_PBSAARC)
         if ( inp == 2 ) call sa_vol( verbose, pbprint, atmlast )
      endif
   end if
! mjhsieh: 1st level doesn't need sa_arc
    
   ! compute time-averaged atomic exposures over the last nsaslag steps
    
!  if ( npbstep <= nsaslag ) then
!     do ip = 1, nsatm
!        sumnmax(ip) = sumnmax(ip) + nmax(ip)
!        avnmax(ip) = sumnmax(ip)/dble(npbstep)
!        sumnexp(ip) = sumnexp(ip) + nexp(ip)
!        avnexp(ip) = sumnexp(ip)/dble(npbstep)
!     end do
!  else
!     wf0 = ONE/dble(nsaslag)
!     wf1 = ONE - wf0
!     do ip = 1, nsatm
!        avnmax(ip) = wf1*avnmax(ip) + wf0*nmax(ip)
!        avnexp(ip) = wf1*avnexp(ip) + wf0*nexp(ip)
!     end do
!  end if

!  GET RADIP3 FOR THE FIRST TIME SA_DRIVER GOT CALLED OR WHEN NO LIGAND/MULTIBLOCK SPECIFIED
   if ( done_sa_driver ) then
!write(600+mytaskid,*)radi(1:natom),radip2(1:natom),radip3(1:natom),avnmax(1:natom),avnexp(1:natom);call mexit(0,0)
      continue
   else
      ! CQ, the average values are fixed to easy reproduction of different runs
      do iatm = 1, natom
         avnmax(iatm) = nmax(iatm)
         avnexp(iatm) = nexp(iatm)
      end do
    
      ! use the exposure information to determine the effective dynamic cavity radii
      ! for modified van der Waals surface for MD
       
      smallsas = 0.05d0*maxsph
      if ( verbose .and. pbprint ) then
         write(6, '(a)')
         write(6, '(a)') 'Atomic solvent accessible surface area:'
      end if
      
      ! In the future radip3 should be read-in for docking
      
      radip3 = ZERO
      do iatm = 1, natom
         if ( radip(iatm) == ZERO ) cycle
      
         if ( avnmax(iatm) <= smallsas ) then
            exposure = ONE ! buried by 1-4 pairs only, no need to increase radius
         else if ( avnmax(iatm) <= 2*smallsas ) then
            exposure = exp(-(avnmax(iatm)-smallsas)**2/(2*smallsas))
            exposure = exposure + avnexp(iatm)/avnmax(iatm) ! transition zone
         else
            exposure = avnexp(iatm)/avnmax(iatm) ! normal exposed atoms
         end if
         if ( exposure > expthresh ) then
            increase = ZERO
         else 
!           increase = ONE - exposure/expthresh
            increase = HALF*( cos( PI*(exposure/expthresh) ) + ONE )
         end if
      
         radip3(iatm) = radi(iatm) + increase*radinc
          
         ! compute SASA for printout
          
         prtsas = prtsas + FOURPI*radip2(iatm)*nexp(iatm)/dble(maxsph)
         if ( verbose .and. pbprint ) write(6, '(i12,f19.13)') &
            iatm, FOURPI*radip2(iatm)*nexp(iatm)/dble(maxsph)
         !           1   126.721289834435
         !1234567890123456789012345678901
         !           1     126.721289834
      end do
      ! This is for future dry-run functionality
      !if (dumpsas) then
      !endif
   end if

!  print *, prtsas
   ! Ligand Mask / Traditional Focusing / Multiple-block Focusing
   done_sa_driver = .true.
    
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Solvent accessible surface calculation
subroutine sa_srf( verbose,pbprint,natom,atmfirst,atmlast )
    
   implicit none
    
   ! Passed variables
    
   logical verbose, pbprint
   integer natom, atmfirst, atmlast
    
   ! Local variables
    
   integer alloc_err(4)
   integer iatm, jatm, jfirst, jlast, ilast, isph, jsph, jp
   integer nsrf, nsurf
   integer fstsph(natom), lstsph(natom)
   _REAL_ sx, sy, sz, xi, yi, zi, dx, dy, dz, dxij, dyij, dzij, d2
   _REAL_, parameter :: small = 0.0001d0
    
   allocate(knockout(  maxsph*natom), stat = alloc_err(1) )
   allocate(  sphcrd(3,maxsph*natom), stat = alloc_err(2) )
   if ( alloc_err( 1)+alloc_err( 2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Allocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if
    
   ! generate solvent accessible points using excluded atoms only
   ! then save these points per atoms for later.
   ! this should be the maximum sas area of this molecule

   nsrf = 0
   do iatm = atmfirst, atmlast
      if ( radip(iatm) == ZERO ) cycle
      fstsph(iatm) = nsrf + 1
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
      jfirst = 1; jlast  = nex(iatm)
        
      ! working on iatm's points
        
      do isph = 1, maxsph
         sx = radip(iatm)*scrd(1,isph) + xi
         sy = radip(iatm)*scrd(2,isph) + yi
         sz = radip(iatm)*scrd(3,isph) + zi
            
         ! loop over exclusion pairs:
            
         do jp = jfirst, jlast
            jatm = iex(jp,iatm)
            if ( radip(jatm) == ZERO )  cycle
            dx = acrd(1,jatm) - sx
            dy = acrd(2,jatm) - sy
            dz = acrd(3,jatm) - sz
            d2 = dx*dx + dy*dy + dz*dz + small
            if (d2 < radip2(jatm)) goto 10
         enddo

         ! if sphere point is not "knocked out" by all its neighbors
         ! update global srf point counter

         nsrf = nsrf + 1

         ! save its coord, wrt the ctr of iatm

         sphcrd(1,nsrf) = sx - xi
         sphcrd(2,nsrf) = sy - yi
         sphcrd(3,nsrf) = sz - zi

10    continue
      end do

      ! remember the last srf point of the atom in global srf list

      lstsph(iatm) = nsrf
   end do  ! iatm = 1, atmlast
 
   knockout(1:nsrf) = .false.

   allocate( fstsdot(         natom), stat = alloc_err(1) )
   allocate( lstsdot(         natom), stat = alloc_err(2) )
   allocate(  srfcrd(3,       nsrf ), stat = alloc_err(3) )
   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Allocation aborted', alloc_err(1:3)
      call mexit(6, 1)
   end if

   ! loop over nblist with cutsa2 (currently cutsa)

   do iatm = atmfirst, atmlast
      if ( radip(iatm) == ZERO ) cycle
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
      jfirst = iar1pb(1,iatm) + 1; jlast  = iar1pb(3,iatm)
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if ( radip(jatm) == ZERO ) cycle

         dxij = xi - acrd(1,jatm)
         dyij = yi - acrd(2,jatm)
         dzij = zi - acrd(3,jatm)
         d2 = dxij**2 + dyij**2 + dzij**2
         if (d2 > (radip(iatm) + radip(jatm))**2) cycle
          
         ! working on iatm's points
          
         do isph = fstsph(iatm), lstsph(iatm)
            if ( knockout(isph) ) cycle
            dx = sphcrd(1,isph) + dxij
            dy = sphcrd(2,isph) + dyij
            dz = sphcrd(3,isph) + dzij
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < radip2(jatm) ) knockout(isph) = .true.
         enddo
          
         ! working on jatm's points
          
         do jsph = fstsph(jatm), lstsph(jatm)
            if ( knockout(jsph) ) cycle
            dx = sphcrd(1,jsph) - dxij
            dy = sphcrd(2,jsph) - dyij
            dz = sphcrd(3,jsph) - dzij
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < radip2(iatm) ) knockout(jsph) = .true.
         enddo
       
      enddo  ! jatm = iprshrt(ip), jp = jfirst, jlast
 
   enddo  ! iatm = atmfirst, atmlast

   ! determine the exposed atoms and buried atoms
   ! save surface dots for eps calculations

   lstsdot = 0
   nsrfdot = 0 ! sa dots for the whole molecule
   nmax = ZERO; nexp = ZERO ! required to initialized for later reading
   do iatm = atmfirst, atmlast
      fstsdot(iatm) = nsrfdot+1
      if ( radip(iatm) == ZERO ) cycle

      nsurf = 0 ! sa dots for this atom
      do isph = fstsph(iatm), lstsph(iatm)
         if ( knockout(isph) ) cycle
         nsurf = nsurf + 1
         nsrfdot = nsrfdot + 1
         srfcrd(1,nsrfdot) = sphcrd(1,isph)! + acrd(1,iatm) for display
         srfcrd(2,nsrfdot) = sphcrd(2,isph)! + acrd(2,iatm) for display
         srfcrd(3,nsrfdot) = sphcrd(3,isph)! + acrd(3,iatm) for display
      end do
      lstsdot(iatm) = nsrfdot

      ! get maximum no. surface dots of an atom
      ! the actual no. surface dots of an atom

      if ( lstsph(iatm) <= fstsph(iatm) ) then
         nmax(iatm) = ZERO
         nexp(iatm) = dble(nsurf)
      else
         nmax(iatm) = dble(lstsph(iatm) - fstsph(iatm) + 1)
         nexp(iatm) = dble(nsurf)
      end if
   enddo

   deallocate(knockout, stat = alloc_err(1) )
   deallocate(  sphcrd, stat = alloc_err(2) )
   if ( alloc_err( 1)+alloc_err( 2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Deallocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if
   
   if ( verbose .and. pbprint ) then
!     write(6,'(a,i6)') 'Number of SA srf points exposed', nsrfdot
      write(6,*) 'Number of SA srf points exposed', nsrfdot
   end if

   ! for InsightII/Sybyl display
   !print *, 'Number of SA srf points exposed', nsrfdot
   !open (unit=55, file='sasrf.pdb')
   !do isph = 1, nsrfdot
   !   write (55,'("ATOM  ",x,i4,2x,"O   ALA",3x,i3,4x,3(f8.3))') isph,1,srfcrd(1:3,isph)
   !   write (55, '("TER")')
   !end do
   !close(55)
   !open (unit=55, file='sasrf.dot')
   !write (55, '("DOTS")')
   !do isph = 1, nsrfdot
   !   write (55,'(4(f8.3,2x))') srfcrd(1,isph), srfcrd(2,isph), srfcrd(3,isph), 300.
   !end do
   !close(55)

end subroutine sa_srf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Solvent accessible arc calculation
subroutine sa_arc( verbose,pbprint,natom,atmfirst,atmlast )
    
   implicit none

   logical verbose, pbprint
   integer natom, atmfirst, atmlast
    
   integer alloc_err(8)
   integer ksrf, lsrf, msrf, isph
   integer fstsph(natom), lstsph(natom), fstsph1(natom), lstsph1(natom)
   character(10) str
   integer iatm, mtri
   integer fsttri(natom),lsttri(natom),fsttri1(natom),lsttri1(natom)

   call pbtimer_start(PBTIME_PBSAARC_SETUP)
   if ( triopt > 0 ) then
      allocate(triatm (3,maxtri*natom))
      allocate(triatm1(3,maxtri*natom))
      allocate(triarc (3,maxtri*natom))
      allocate(triarc1(3,maxtri*natom))
      allocate(tricrd (3,maxtri*natom))
      allocate(tricrd1(3,maxtri*natom))
      allocate(knocktri(maxtri*natom))
   else
      allocate(triatm  (3,1))
      allocate(triatm1 (3,1))
      allocate(triarc  (3,1))
      allocate(triarc1 (3,1))
      allocate(tricrd  (3,1))
      allocate(tricrd1 (3,1))
      allocate(knocktri(1))
   end if

   allocate( fstadot(            natom), stat = alloc_err(1) )
   allocate( lstadot(            natom), stat = alloc_err(2) )
   allocate(  fstarc(            natom), stat = alloc_err(3) ) 
   allocate(  lstarc(            natom), stat = alloc_err(4) )
   allocate(    marc(            natom), stat = alloc_err(5) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc()1: Allocates aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
    
   ! safe initialization to silence potential compiler error messages 

   triatm = 0; triatm1 = 0; triarc = 0; triarc1 = 0
   fsttri = 0; lsttri = 0; fsttri1 = 0; lsttri1 =0 
   ntri = 0; mtri = 0
   fstarc = 0; lstarc = 0; marc = 0; fstsph = 0; lstsph = 0 
   fstsph1 = 0; lstsph1 = 0; ksrf = 0; msrf = 0; lsrf = 0
   tricrd = ZERO; tricrd1 = ZERO
    
   ! get all possible arcs and arc dots

   allocate(  m2narc(  maxarc   ,natom), stat = alloc_err(1) )
   allocate(  arcatm(2,maxarc   *natom), stat = alloc_err(2) )
   allocate(  savarc(3,maxarc   *natom), stat = alloc_err(3) )
   allocate( savactr(3,maxarc   *natom), stat = alloc_err(4) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc()2: Allocates aborted', alloc_err(1:4)
      call mexit(6, 1)
   end if
    
   allocate(knockout(  maxarcdot*natom), stat = alloc_err(1) )
   allocate(  sphcrd(3,maxarcdot*natom), stat = alloc_err(2) )
   allocate( sphcrd1(3,maxarcdot*natom), stat = alloc_err(3) )
   allocate(  spharc(  maxarcdot*natom), stat = alloc_err(4) )
   allocate( spharc1(  maxarcdot*natom), stat = alloc_err(5) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc()3: Allocates aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
   call pbtimer_stop(PBTIME_PBSAARC_SETUP)

   !spharc=0
   !spharc1=0
   !m2narc=0
    
   call pbtimer_start(PBTIME_PBCIRCLE)
   call circle( natom,atmfirst,atmlast,ksrf,arcres,fstsph,lstsph, & 
                ntri,mtri,triatm,tricrd,fsttri,lsttri,triarc, &
                triatm1,tricrd1,fsttri1,lsttri1,triarc1 )
   call pbtimer_stop(PBTIME_PBCIRCLE)

   ! knock out points overlapped by other atoms
   ! step 1: loop over exclusion list
    
   call pbtimer_start(PBTIME_PBEXCLUDE)
   call exclud( atmfirst,atmlast,4,1,ksrf,lsrf,fstsph,lstsph,fstsph1,lstsph1,spharc1,spharc,sphcrd,sphcrd1, & 
                mtri,ntri,triatm1,tricrd1,fsttri1,lsttri1,triarc1, &
                triatm,tricrd,fsttri,lsttri,triarc )
   call pbtimer_stop(PBTIME_PBEXCLUDE)
    
   ! step 2: loop over short nblist with cutsa1 (currently cutfd) of 5 A
    
   call pbtimer_start(PBTIME_PBEXCLUDE)
   call exclud( atmfirst,atmlast,1,2,lsrf,msrf,fstsph1,lstsph1,fstsph,lstsph,spharc,spharc1,sphcrd1,sphcrd, & 
                ntri,mtri,triatm,tricrd,fsttri,lsttri,triarc, &
                triatm1,tricrd1,fsttri1,lsttri1,triarc1 )
   call pbtimer_stop(PBTIME_PBEXCLUDE)
    
   ! step 3: loop over long nblist with cutsa2 (currently cutsa) 9 A
   ! note that narcdot, fstadot, lstadot, dotarc, arccrd are returned and used elsewhere
    
   call pbtimer_start(PBTIME_PBSAARC_SETUP)
   allocate(  dotarc(   msrf), stat = alloc_err(1) )
   allocate(  arccrd(3, msrf+mtri), stat = alloc_err(2) )
   if ( alloc_err(1)+alloc_err(2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Allocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if
   call pbtimer_stop(PBTIME_PBSAARC_SETUP)
   !dotarc=0
   !arccrd=0d0
    
   call pbtimer_start(PBTIME_PBEXCLUDE)
   call exclud( atmfirst,atmlast,2,3,msrf,narcdot,fstsph,lstsph,fstadot,lstadot,spharc1,dotarc,sphcrd,arccrd, & 
                mtri,ntri,triatm1,tricrd1,fsttri1,lsttri1,triarc1, &
                triatm,tricrd,fsttri,lsttri,triarc )
   call pbtimer_stop(PBTIME_PBEXCLUDE)
    
   call pbtimer_start(PBTIME_PBSAARC_SETUP)
   if ( triopt > 0 ) then
      do isph = 1, ntri
         narcdot = narcdot + 1
         arccrd(1,narcdot) = tricrd(1,isph)
         arccrd(2,narcdot) = tricrd(2,isph)
         arccrd(3,narcdot) = tricrd(3,isph)
!        write(6,*) 'trimer',triatm(1:3,isph)
      end do
   end if

   deallocate(knockout, stat = alloc_err(1) )
   deallocate(  sphcrd, stat = alloc_err(2) )
   deallocate( sphcrd1, stat = alloc_err(3) )
   deallocate(  spharc, stat = alloc_err(4) )
   deallocate( spharc1, stat = alloc_err(5) )
   deallocate(triatm1)
   deallocate(triarc1)
   deallocate(tricrd )
   deallocate(tricrd1)
   deallocate(knocktri)
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Deallocation aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
   call pbtimer_stop(PBTIME_PBSAARC_SETUP)
    
   if ( verbose .and. pbprint ) then
      write(6,*) 'Number of SA arcs generated', narc
!     write(6,*) 'Number of trimer points', ntri
      write(6,'(A,I12,A,E11.6)') &
             ' Number of SA arc points exposed', narcdot, &
             '  with resolution (A) = ', arcres
   end if
    
   ! for InsightII/Sybyl display
   !print *, 'Number of SA arcs generated', narc
   !print *, 'Number of SA arc points exposed', narcdot, ' with resolution (A)', arcres
   !open (unit=55, file='saarc.pdb')
   !do isph = 1, narcdot
   !   write (55,'("ATOM  ",x,i4,2x,"O   ALA",3x,i3,4x,3(f8.3))') isph,1,arccrd(1:3,isph)
   !   write (55, '("TER")')
   !enddo
   !close(55)

!  dumping arcs and trimer points (for visualization)
!  open(unit=58,file='saarc.dot')
!  write (58, *) narcdot-ntri
!  write (58, '("DOTS")')
!  do isph = 1, narcdot-ntri
!     write(str,'(i10)') isph
!     str = adjustl(str)
!     write(58,'(a,3(f20.15),2x)') "H"//str, arccrd(1:3,isph)
!  enddo
!  close(58)
!  write (59, *) ntri
!  write (59, '("DOTS")')
!  do isph = 1, ntri
!     write(str,'(i10)') isph
!     str = adjustl(str)
!     write(59,'(a,3(f20.15),2x)') "H"//str, tricrd(1:3,isph)
!  enddo

   !open(unit=55,file='saarc.dot')
   !write (55, '("DOTS")')
   !do isph = 1, narcdot
   !   write(55,'(4(f8.3,2x))') arccrd(1:3,isph), 300.
   !enddo
   !close(55)

end subroutine sa_arc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Generate maximum solvent-accessible circles
subroutine circle( natom,atmfirst,atmlast,msrf,cstep,fstsph,lstsph, &
                   ntri,mtri,triatm,tricrd,fsttri,lsttri,triarc, &
                   triatm1,tricrd1,fsttri1,lsttri1,triarc1 )
    
   implicit none
    
   ! Passed variables
    
   integer natom, atmfirst, atmlast, msrf, fstsph(*), lstsph(*)
   _REAL_ cstep
   integer ntri, mtri, triatm(3,*), triatm1(3,*), triarc(3,*), triarc1(3,*)
   integer fsttri(natom), lsttri(natom), fsttri1(natom), lsttri1(natom)
   _REAL_  tricrd(3,*), tricrd1(3,*)
    
   ! Local variables
    
   integer  jp, iatm, ifirst, ilast, jatm, jfirst, jlast, isph, nspha
   _REAL_ ri, xi, yi, zi, radius, aphi, rj
   _REAL_ x, y, dxij, dyij, dzij, d2, d, rd, d1, d3, side2(3)
   _REAL_ sintheta, costheta, cosaij, cosaji
   _REAL_ arcctr(3), frc_ctr
   _REAL_ psi, cosphi, sinphi
   integer katm, kp, atms(3), arc
   _REAL_ xj, yj,zj , xk, yk, zk, rk
   _REAL_ dxjk, dyjk, dzjk, dxik, dyik, dzik
   _REAL_ crn1(3), crn2(3)
   logical isexist, iscycle
   integer kfirst, klast, pp, atm
   integer n1, n2, n3
   _REAL_ xl, yl, zl, rl2, dx, dy, dz
   integer lp, lfirst, llast, latm
   logical knock

   _REAL_, parameter :: small = 0.000001d0

   n1 = 0; n2 = 0; n3 = 0
   narc = 0; msrf = 0
   do iatm = atmfirst, atmlast
       
      ri = radip(iatm)
      if ( ri == ZERO ) then
         fstsph(iatm) = lstsph(iatm) + 1
         cycle
      end if
      fstarc(iatm) = narc + 1
      fstsph(iatm) = msrf + 1
      fsttri1(iatm) = mtri + 1
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
       
      jfirst = iar1pb(4,iatm-1) + 1; jlast  = iar1pb(3,iatm)
      do jp = jfirst, jlast
          
         jatm = iprshrt(jp)
         rj = radip(jatm)
         if ( rj == ZERO ) cycle
          
         xj = acrd(1,jatm); yj = acrd(2,jatm); zj = acrd(3,jatm)
         dxij = xj - xi
         dyij = yj - yi
         dzij = zj - zi
         d2 = dxij**2 + dyij**2 + dzij**2
         if ( d2 >= (ri + rj)**2 .or. d2 <= (ri - rj)**2 ) cycle
          
         ! setting up indexes ...
          
         narc = narc + 1
         marc(iatm) = marc(iatm) + 1; marc(jatm) = marc(jatm) + 1
         m2narc(marc(iatm),iatm) = narc
         m2narc(marc(jatm),jatm) = narc
          
         ! getting the center of the arc circle
          
         d = sqrt(d2); rd = ONE/d
         cosaij = (d2 + radip2(iatm) - radip2(jatm))/(TWO*d*ri)
         cosaji = (d2 + radip2(jatm) - radip2(iatm))/(TWO*d*rj)
         arcatm(1,narc) = jatm; arcatm(2,narc) = iatm
         savarc(1,narc) = cosaij; savarc(2,narc) = cosaji; savarc(3,narc) = rd

         radius = ri*sqrt(ONE - cosaij**2)
         frc_ctr = ri*cosaij*rd
         arcctr(1) = xi + frc_ctr*dxij
         arcctr(2) = yi + frc_ctr*dyij
         arcctr(3) = zi + frc_ctr*dzij
         savactr(1:3,narc) = arcctr(1:3)
          
         ! determining the orientation of the arc w.r.t. the z-axis
          
         if ( abs(dyij) < small .and. abs(dxij) < small ) then
            ! the vector is aligned with the z-axis
            cosphi = ONE
            sinphi = ZERO

         else if ( abs(dyij) < small ) then
            ! the vector is aligned with the x-axis
            cosphi = ZERO
            if ( dxij > ZERO ) then
               sinphi = ONE
            else
               sinphi = - ONE
            end if

            ! otherwise ...
         else if ( dyij > ZERO ) then
            aphi = atan(dxij/dyij)
            cosphi = cos(aphi)
            sinphi = sin(aphi)
         else
            aphi = PI + atan(dxij/dyij)
            cosphi = cos(aphi)
            sinphi = sin(aphi)
         end if
         costheta = dzij/d
         sintheta = sqrt(ONE - costheta**2)
          
         ! generating dots on the arc
          
         nspha = int(TWOPI*radius/cstep/EIGHT)*8
         psi = ZERO
         do isph = 1, nspha
            msrf = msrf + 1
            x = radius*cos(psi)
            y = radius*sin(psi)
            sphcrd(1,msrf) = arcctr(1) + cosphi*(x) + sinphi*(   costheta*y)
            sphcrd(2,msrf) = arcctr(2) - sinphi*(x) + cosphi*(   costheta*y)
            sphcrd(3,msrf) = arcctr(3) +                                        ( - sintheta*y)
            spharc1(msrf) = narc
            psi = psi + TWOPI/nspha
         end do

         if ( triopt > 0 ) then
            ! looking for all groups of three atoms that can possibly form trimer
            ! points

            do kp = jp+1, jlast

               katm = iprshrt(kp)
               rk = radip(katm)
               if ( rk == ZERO ) cycle
             
               iscycle = .true.
               if ( jatm < katm ) then
                  kfirst = iar1pb(4,jatm-1) + 1
                  klast  = iar1pb(3,jatm)
                  pp = kfirst 
                  do while ( pp <= klast )
                     atm = iprshrt(pp)
                     if ( katm == atm ) then
                        iscycle = .false.
                        exit
                     end if
                     pp = pp + 1
                  end do
               else 
                  kfirst = iar1pb(4,katm-1) + 1
                  klast  = iar1pb(3,katm)
                  pp = kfirst 
                  do while ( pp <= klast )
                     atm = iprshrt(pp)
                     if ( jatm == atm ) then
                        iscycle = .false.
                        exit
                     end if
                     pp = pp + 1
                  end do
               end if
!              iscycle = .false.
               if ( iscycle ) then
                  n1 = n1 + 1
               else
                  n2 = n2 + 1
               end if
               if ( iscycle ) cycle

               xk = acrd(1,katm); yk = acrd(2,katm); zk = acrd(3,katm)
               dxjk = xk - xj
               dyjk = yk - yj
               dzjk = zk - zj
               d1 = dxjk**2 + dyjk**2 + dzjk**2
               if ( d1 >= (rj + rk)**2 .or. d1 <= (rj - rk)**2 ) cycle

               dxik = xk - xi
               dyik = yk - yi
               dzik = zk - zi
               d3 = dxik**2 + dyik**2 + dzik**2
               if ( d3 >= (ri + rk)**2 .or. d3 <= (ri - rk)**2 ) cycle

               n3 = n3 + 1
               atms(1) = iatm
               atms(2) = jatm
               atms(3) = katm
               side2(1) = d2
               side2(2) = d3
               side2(3) = d1
               ! ruling out groups of three atoms 
               ! whose three arcs can't cross each other
               call get_arccrn(natom,acrd,atms,isexist,crn1,crn2)
               if ( isexist ) then
                  lfirst = iar1pb(4,iatm-1) + 1; llast  = iar1pb(3,iatm)
                  knock = .false.
                  do lp = lfirst, llast
                     latm = iprshrt(lp)
                     if ( latm == jatm .or. latm == katm ) cycle
                     rl2 = radip2(latm)
                     if ( rl2 == ZERO ) cycle
                     xl = acrd(1,latm); yl = acrd(2,latm); zl = acrd(3,latm)
                     dx = crn1(1) - xl
                     dy = crn1(2) - yl
                     dz = crn1(3) - zl
                     d2 = dx*dx + dy*dy + dz*dz + 1.d-9
                     if ( d2 < rl2 ) then
                        knock = .true.
                        exit
                     end if
                  end do
                  if ( .not. knock ) then
                     mtri = mtri + 1
                     triatm1(1:3,mtri) = atms
                     tricrd1(1:3,mtri) = crn1
                     triarc1(1,mtri) = narc
                  end if
                  knock = .false.
                  do lp = lfirst, llast
                     latm = iprshrt(lp)
                     if ( latm == jatm .or. latm == katm ) cycle
                     rl2 = radip2(latm)
                     if ( rl2 == ZERO ) cycle
                     xl = acrd(1,latm); yl = acrd(2,latm); zl = acrd(3,latm)
                     dx = crn2(1) - xl
                     dy = crn2(2) - yl
                     dz = crn2(3) - zl
                     d2 = dx*dx + dy*dy + dz*dz + 1.d-9
                     if ( d2 < rl2 ) then
                        knock = .true.
                        exit
                     end if
                  end do
                  if ( .not. knock ) then
                     mtri = mtri + 1
                     triatm1(1:3,mtri) = atms
                     tricrd1(1:3,mtri) = crn2
                     triarc1(1,mtri) = narc
                  end if
               end if
            end do 
         end if
          
      end do  ! jatm = iprshrt(ip), jp =  jfirst, jlast
      lstarc(iatm) = narc
      lstsph(iatm) = msrf
      lsttri1(iatm) = mtri
       
      if ( narc > maxarc*natom ) then
         write(6,*) 'SA Bomb in circle(): Stored surface arcs over limit', iatm, narc
         call mexit(6,1)
      end if
      if ( marc(iatm) > maxarc ) then
         write(6,*) 'SA Bomb in circle(): Stored surface arcs over limit', iatm, marc(iatm)
         call mexit(6,1)
      end if
      if ( msrf > maxarcdot*natom ) then
         write(6,*) 'SA Bomb in circle(): Stored surface points over limit', msrf
         call mexit(6,1)
      end if
        
   end do  ! iatm = atmfirst, atmlast
   fstarc(atmlast) = narc + 1
   lstarc(atmlast) = narc + 1
   fstsph(atmlast) = lstsph(atmlast) + 1

!  write(6,*) 'n1n2n3',n1,n2,n3
!  write(6,*) 'msrf',msrf
!  write(6,*) 'mtri',mtri
   if ( triopt > 0 ) then
      ! searching for the other two arcs a trimer point belongs to
      do jp = 1, mtri
         katm = triatm1(3,jp)
         do kp = 1, marc(katm)
            arc = m2narc(kp,katm)
            if ( triatm1(1,jp) == arcatm(1,arc) .or. &
                 triatm1(1,jp) == arcatm(2,arc)) then
               triarc1(2,jp) = arc
               cycle
            end if
            if ( triatm1(2,jp) == arcatm(1,arc) .or. &
                 triatm1(2,jp) == arcatm(2,arc) ) then
               triarc1(3,jp) = arc
               cycle
            end if
         end do
      end do
   end if

end subroutine circle
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_arccrn(natom,crd,idx,isexist,crn1,crn2)

   implicit none

#  include "pb_constants.h"

   ! Passed
   integer natom, idx(3)
   _REAL_ crd(3,natom) 
   _REAL_ prob, r(natom)
   _REAL_ crn1(3), crn2(3)
   logical isexist

   ! Local
   integer m, n, mp, np, j
   _REAL_ a1(3), a2(3), a3(3)
   _REAL_ p1(3), p2(3), q1(3), q2(3), p(3), d(3)
   _REAL_ u(3,3), v(3,3), w(3), b(3), t(3)
   _REAL_ TOL, wmax, thresh, h, d1, d2
   _REAL_ m1,m2,m3,t1,t2,t3,det,tt(3)
 
   isexist = .true.
   crn1 = ZERO; crn2 = ZERO
   a1 = crd(1:3,idx(1)); a2 = crd(1:3,idx(2)); a3 = crd(1:3,idx(3))
   p1 = a2 - a1; p2 = a3 - a1
   q1 = a2 + a1; q2 = a3 + a1
   d(1) = radip(idx(1))
   d(2) = radip(idx(2))
   d(3) = radip(idx(3))
   p(1) = p1(2)*p2(3) - p2(2)*p1(3)
   p(2) = p2(1)*p1(3) - p1(1)*p2(3)
   p(3) = p1(1)*p2(2) - p2(1)*p1(2)
   b(1) = (dot_product(p1,q1) + (d(1)+d(2))*(d(1)-d(2)))*HALF
   b(2) = (dot_product(p2,q2) + (d(1)+d(3))*(d(1)-d(3)))*HALF
   b(3) = dot_product(p,a1)
   u(1,1:3) = p1
   u(2,1:3) = p2
   u(3,1:3) = p
  
   TOL = 1.d-5
   m1 = u(2,1)*u(3,2)-u(2,2)*u(3,1)
   m2 = u(1,1)*u(3,2)-u(1,2)*u(3,1)
!  m3 = u(1,1)*u(2,2)-u(1,2)*u(2,1)
   m3 = p(3)
   det = u(1,3)*m1-u(2,3)*m2+u(3,3)*m3
   if ( abs(det) < TOL ) then
      isexist = .false.
      return
   end if 
   t3 = b(1)*m1-b(2)*m2+b(3)*m3
   t(3) = t3/det
   if ( abs(m1) > TOL ) then
      t(2) = (b(3)-u(3,3)*t(3))*u(2,1) - (b(2)-u(2,3)*t(3))*u(3,1)
      t(2) = t(2)/m1
   elseif ( abs(m2) > TOL ) then
      t(2) = (b(3)-u(3,3)*t(3))*u(1,1) - (b(1)-u(1,3)*t(3))*u(3,1)
      t(2) = t(2)/m2
   elseif ( abs(m3) > TOL ) then
      t(2) = (b(2)-u(2,3)*t(3))*u(1,1) - (b(1)-u(1,3)*t(3))*u(2,1)
      t(2) = t(2)/m3
   else
      isexist = .false.
      return
   end if
   if ( abs(u(1,1)) > TOL ) then
      t(1) = b(1) - u(1,3)*t(3) - u(1,2)*t(2)
      t(1) = t(1)/u(1,1)
   elseif ( abs(u(2,1)) > TOL ) then
      t(1) = b(2) - u(2,3)*t(3) - u(2,2)*t(2)
      t(1) = t(1)/u(2,1)
   elseif ( abs(u(3,1)) > TOL ) then
      t(1) = b(1) - u(3,3)*t(3) - u(3,2)*t(2)
      t(1) = t(1)/u(3,1)
   else
      isexist = .false.
      return
   end if
      
!  tt = t
!  mp = 3; np = 3
!  m = 3; n = 3
!  TOL = 1.d-5
!  call svdcmp(u,m,n,mp,np,w,v)
!  wmax = ZERO
!  do j = 1, n
!    if ( w(j) > wmax ) wmax = w(j)
!  end do
!  thresh = TOL * wmax
!  do j = 1, n
!    if ( w(j) < thresh ) w(j) = ZERO
!  end do
!  call svbksb(u,w,v,m,n,mp,np,b,t)

!  if ( dot_product(tt-t,tt-t) > 1.d-14 ) then 
!     write(6,*) "inconsistent",tt(1:3),t(1:3)
!     write(6,*) "u1",u(1,1:3)
!     write(6,*) "u2",u(2,1:3)
!     write(6,*) "u3",u(3,1:3)
!     write(6,*) "b",b(1:3)
!     write(6,*) "det",det,t3
!     stop
!  end if

!  t = tt
!  p = p/sqrt(p(1)**2+p(2)**2+p(3)**2)
   d1 = d(1)**2
   d2 = dot_product(t-a1,t-a1)
!  t2 = dot_product(tt-a1,tt-a1)
!  if ( d1 <= d2 .and. d1 > t2 ) then
!     write(6,*) "inconsistent",d1,d2,t2
!  end if
   if ( d1 <= d2 ) then
      isexist = .false.
      return
   end if 
   h = sqrt((d1-d2)/dot_product(p,p))
   crn1 = t + h*p
   crn2 = t - h*p

end subroutine get_arccrn
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Knock out dots on the circles that are overlapped by close pairs
!  mjhsieh: spharc1 was dotarc outside.
subroutine exclud( atmfirst,atmlast,istart,istop,lsrf,msrf,fstsph,lstsph,fstsph1,lstsph1,spharc,spharc1,sphcrd,sphcrd1, & 
                   ntri,mtri,triatm,tricrd,fsttri,lsttri,triarc, &
                   triatm1,tricrd1,fsttri1,lsttri1,triarc1 )
    
   implicit none
    
   integer atmfirst, atmlast, istart, istop, lsrf, msrf
   integer fstsph(*), lstsph(*), fstsph1(*), lstsph1(*), spharc(*), spharc1(*)
   _REAL_ sphcrd(3,*), sphcrd1(3,*)
    
   integer ntri, mtri, triatm(3,*), triatm1(3,*), triarc(3,*), triarc1(3,*)
   integer fsttri(natom),lsttri(natom),fsttri1(natom),lsttri1(natom)
   _REAL_  tricrd(3,*), tricrd1(3,*)

   integer jp, iatm, jatm, ifirst, ilast, jfirst, jlast, isph, jsph
   _REAL_ xi, yi, zi, xj, yj, zj
   _REAL_ ri2, rj2, dx, dy, dz, d2
!integer i, myoutflag(natom),katm

   _REAL_, parameter :: small = 1.d-9

   knockout(1:lsrf) = .false.
   knocktri(1:ntri) = .false.
   
!if ( ligand ) then
!  myoutflag=outflag
!else
!myoutflag = 0
!do iatm = 1 , natom
!   read(5555,'(i4)') myoutflag(i)
!   write(9999,*)radip2(iatm)
!enddo
!rewind
!endif
   do iatm = atmfirst, atmlast
       
      ri2 = radip2(iatm)
      if ( ri2 == ZERO ) cycle 
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
       
      if ( istart == 4 ) then
         jfirst = iar1pb(istart,iatm-1) + 1
      else
         jfirst = iar1pb(istart,iatm) + 1
      end if
      jlast  = iar1pb(istop,iatm)
      do jp = jfirst, jlast
          
         jatm = iprshrt(jp); rj2 = radip2(jatm)
         if ( rj2 == ZERO ) cycle
         xj = acrd(1,jatm); yj = acrd(2,jatm); zj = acrd(3,jatm)
           
         ! working on iatm's points
          
         do isph = fstsph(iatm), lstsph(iatm)
            if ( knockout(isph) ) cycle
            dx = sphcrd(1,isph) - xj
            dy = sphcrd(2,isph) - yj
            dz = sphcrd(3,isph) - zj
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < rj2 ) knockout(isph) = .true.
         end do

!if (ligand)then
!write(3333,*) iatm,jatm,outflag(iatm)
!else
!write(3333,*) iatm,jatm,myoutflag(iatm)
!endif
          
         ! working on jatm's points
           
         do jsph = fstsph(jatm), lstsph(jatm)
            if ( knockout(jsph) ) cycle
            dx = sphcrd(1,jsph) - xi
            dy = sphcrd(2,jsph) - yi
            dz = sphcrd(3,jsph) - zi
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < ri2 ) knockout(jsph) = .true.
         end do

         if ( triopt > 0 ) then
            do jsph = fsttri(jatm), lsttri(jatm)
               if ( knocktri(jsph) ) cycle
               dx = tricrd(1,jsph) - xi
               dy = tricrd(2,jsph) - yi
               dz = tricrd(3,jsph) - zi
               d2 = dx*dx + dy*dy + dz*dz + small
               if ( d2 < ri2 ) knocktri(jsph) = .true.
            end do
         end if 
       
      end do  ! jatm = iprshrt(ip), jp = jfirst, jlast
    
   end do  ! iatm = atmfirst, atmlast
    
   ! condense exposed points
    
   msrf = 0
   do iatm = atmfirst, atmlast
      if ( radip2(iatm) == ZERO ) then
         fstsph1(iatm) = lstsph1(iatm) + 1
         cycle
      end if
      fstsph1(iatm) = msrf + 1
      do isph = fstsph(iatm), lstsph(iatm)
         if ( knockout(isph) ) cycle
         msrf = msrf + 1
         sphcrd1(1,msrf) = sphcrd(1,isph)
         sphcrd1(2,msrf) = sphcrd(2,isph)
         sphcrd1(3,msrf) = sphcrd(3,isph)
         spharc1(msrf) = spharc(isph)
!jatm = arcatm(1,spharc(isph)) 
!katm = arcatm(2,spharc(isph)) 
!if ( myoutflag(jatm) == 0 .or. myoutflag(katm) == 0 ) &
!write(3333,*) min(jatm,katm), max(jatm,katm), iatm, sngl(sphcrd1(1:3, msrf))
      end do
      lstsph1(iatm) = msrf
!write(2020,*) msrf, fstsph1(iatm), lstsph1(iatm), iatm, outflag(iatm)
   end do
   fstsph1(atmlast) = lstsph1(atmlast) + 1
    
   if ( triopt > 0 ) then
      mtri = 0
      do iatm = atmfirst, atmlast
         if ( radip2(iatm) == ZERO ) then
            fsttri1(iatm) = lsttri1(iatm) + 1
            cycle
         end if
         fsttri1(iatm) = mtri + 1
         do isph = fsttri(iatm), lsttri(iatm)
            if ( knocktri(isph) ) cycle
            mtri = mtri + 1
            tricrd1(1:3,mtri) = tricrd(1:3,isph)
            triatm1(1:3,mtri) = triatm(1:3,isph)
            triarc1(1:3,mtri) = triarc(1:3,isph)
         end do
         lsttri1(iatm) = mtri
      end do
      fsttri1(atmlast) = lsttri1(atmlast) + 1
   end if

end subroutine exclud
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Volume within the solvent accessible surface
subroutine sa_vol( verbose,pbprint,atmlast )
    
   implicit none
    
   ! Passed variables
    
   logical verbose, pbprint
   integer atmlast
    
   ! Local variables
    
   integer alloc_err
   integer iatm, volnum, nbuffer, xm, ym, zm, xmymzm!, i, j, k
   integer, allocatable :: insas(:) 
    
   _REAL_ xmin, xmax, ymin, ymax, zmin, zmax, xbox, ybox, zbox
   _REAL_ htmp, rh, gox, goy, goz, range1, xi, yi, zi
   _REAL_ gcrd(3, atmlast)
    
   ! to silence valgrind errors
   volnum = 0

   htmp = 0.5d0; rh = ONE/htmp
   nbuffer = 2*(int(TWO*vprob*rh)+1)+1
    
   ! set bounding box center for all atoms
    
   xmin = 9999.0d0; ymin = 9999.0d0; zmin = 9999.0d0
   xmax = -9999.0d0; ymax = -9999.0d0; zmax = -9999.0d0
   do iatm = 1, atmlast
      if ( acrd(1,iatm)-radip(iatm) .lt. xmin ) xmin = acrd(1,iatm)-radip(iatm)
      if ( acrd(1,iatm)+radip(iatm) .gt. xmax ) xmax = acrd(1,iatm)+radip(iatm)
      if ( acrd(2,iatm)-radip(iatm) .lt. ymin ) ymin = acrd(2,iatm)-radip(iatm)
      if ( acrd(2,iatm)+radip(iatm) .gt. ymax ) ymax = acrd(2,iatm)+radip(iatm)
      if ( acrd(3,iatm)-radip(iatm) .lt. zmin ) zmin = acrd(3,iatm)-radip(iatm)
      if ( acrd(3,iatm)+radip(iatm) .gt. zmax ) zmax = acrd(3,iatm)+radip(iatm)
   enddo
   xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO
   xbox = nint(xbox*rh)*htmp; ybox = nint(ybox*rh)*htmp; zbox = nint(zbox*rh)*htmp
   if ( verbose .and. pbprint ) then
      write(6, '(1x,a,3f10.3)') ' SAV: Bounding Box Center:  ', xbox, ybox, zbox
      write(6, '(1x,a,3f10.3)') ' SAV: Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') ' SAV: Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') ' SAV: Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
   end if
    
   xm = nint( (xmax - xmin)*rh ) + nbuffer; xm = 2*nint( dble(xm)*HALF ) + 1
   ym = nint( (ymax - ymin)*rh ) + nbuffer; ym = 2*nint( dble(ym)*HALF ) + 1
   zm = nint( (zmax - zmin)*rh ) + nbuffer; zm = 2*nint( dble(zm)*HALF ) + 1
   xmymzm = xm*ym*zm
   if ( verbose .and. pbprint ) write(6, '(a,1x,3i5)') ' SAV: Grid dimension ', xm, ym, zm
   gox = - dble(xm+1)*htmp*HALF + xbox
   goy = - dble(ym+1)*htmp*HALF + ybox
   goz = - dble(zm+1)*htmp*HALF + zbox
   if ( verbose .and. pbprint ) write(6, '(a,1x,3f10.3)') ' SAV: Grid origin ', gox, goy, goz
    
   do iatm = 1, atmlast
      gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
      gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
      gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
   end do
    
   allocate( insas(xmymzm), stat = alloc_err )
    
   insas(1:xmymzm) = -1
   do iatm = 1, atmlast
      range1 = radip(iatm)
      if ( range1 == ZERO ) cycle
      range1 = (range1-sprob+vprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exsasph( xm, ym, zm, range1, xi, yi, zi, insas )
   end do
    
   call calsav( volnum, xm, ym, zm, insas)
    
   deallocate( insas, stat = alloc_err )
    
   prtsav = volnum*(htmp)**3
    
end subroutine sa_vol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Count molecular volume within solvent accessible surface
subroutine calsav( volnum, xm, ym, zm, insas )

   integer volnum, xm, ym, zm
   integer insas(xm, ym, zm)

   integer i, j, k
    
!   open (unit=57, file='sav.dot')
!   write (57, '("DOTS")')
   volnum = 0
   do i = 1, xm
   do j = 1, ym
   do k = 1, zm
      if (insas(i, j, k) == 1 ) then
         volnum = volnum + 1
!         xtmp = i * htmp + gox
!         ytmp = j * htmp + goy
!         ztmp = k * htmp + goz
!         write (57,'(4(f8.3,2x))') xtmp, ytmp, ztmp, 300.
      end if
   end do
   end do
   end do
!  close(57)

end subroutine calsav   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark atomic volume within solvent accessible surface
subroutine exsasph( xm,ym,zm,range1,xi,yi,zi,insph )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom of index iatm as 1. Modified from UHBD
   ! (Comp. Phys. Comm. 91:57-95, 1995) routines excrx() and exsrfx() by
   ! Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   implicit none

   ! Passed variables

   integer  xm, ym, zm
   integer  insph(xm,ym,zm)
   _REAL_ range1, xi, yi, zi
 
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3

   lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = max(1,ceiling(xi - range3)); highi = min(xm,floor(xi + range3))
            do i = lowi, highi
               insph(i,j,k) = 1 
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exsasph

end subroutine sa_driver

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
subroutine saslave_init(natom)
   implicit none
   include "mpif.h"
#  include "parallel.h"
#  include "extra.h"
   integer natom
   !LOCAL
   integer ierr
   !MPI initialization VERY BAD IMPLEMENTATION
   if ( .not. master ) then
      allocate(    radi(    natom),STAT=ierr)
      allocate(   radip(    natom),STAT=ierr)
      allocate(  radip2(    natom),STAT=ierr)
      allocate(  radip3(    natom),STAT=ierr)
      allocate(  nzratm(    natom),STAT=ierr)
      allocate(    nmax(    natom),STAT=ierr)
      allocate(    nexp(    natom),STAT=ierr)
      allocate( sumnmax(    natom),STAT=ierr)
      allocate( sumnexp(    natom),STAT=ierr)
      allocate(  avnmax(    natom),STAT=ierr)
      allocate(  avnexp(    natom),STAT=ierr)
      allocate(   mdsig(    natom),STAT=ierr)!dummy for multiblock so far
      allocate(    rmin(    natom),STAT=ierr)!dummy for multiblock so far
   end if
   call MPI_BCAST(   maxsph,       1,         MPI_INTEGER,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(   maxtri,       1,         MPI_INTEGER,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(   maxarc,       1,         MPI_INTEGER,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(maxarcdot,       1,         MPI_INTEGER,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(   arcres,       1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(    dprob,       1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(   radinc,       1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(    iprob,       1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(expthresh,       1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BCAST(  radi(1),   natom,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
   if ( .not. master ) then
      allocate(scrd(3,maxsph),STAT=ierr);REQUIRE(ierr==0)
   end if
   call MPI_BCAST(scrd(1,1),3*maxsph,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr);REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
!print *,radinc
end subroutine saslave_init
#endif /*def MPI*/
#endif /*ndef SANDER or LIBPBSA*/

end module solvent_accessibility
