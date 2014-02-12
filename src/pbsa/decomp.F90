
! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
#define ASSERT(e) if(.not.(e)) call croak(__FILE__,__LINE__)

module decomp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! --- Module for energy decomposition
!
!     Holger Gohlke
!     5.1.2004
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer, public, dimension(:), allocatable  :: jgroup, index, irespw

integer, private, parameter                 :: ndectype = 28
!  number of "enum" variables following

integer, private, parameter                 :: &

            ! DECOMP contributions for sidechains
            sideintsel =  1, &  ! DECOMP -ene int sel
            sideintind =  2, &  ! DECOMP  ene int ind
            sidevdwsel =  3, &  ! DECOMP -ene vdw sel
            sidevdwind =  4, &  ! DECOMP  ene vdw ind
            sidevdwdir =  5, &  ! DECOMP  ene vdw dir
            sideeelsel =  6, &  ! DECOMP -ene eel sel
            sideeelind =  7, &  ! DECOMP  ene eel ind
            sideeeldir =  8, &  ! DECOMP  ene eel dir
            sidepolsel =  9, &  ! DECOMP -ene pol sel
            sidepolind = 10, &  ! DECOMP  ene pol ind
            sidepoldir = 11, &  ! DECOMP  ene pol dir
            sidesassel = 12, &  ! DECOMP -ene sas sel
            sidesasind = 13, &  ! DECOMP  ene sas ind
            sidesasdir = 14, &  ! DECOMP  ene sas dir

            ! DECOMP contributions for backbone
            backintsel = 15, &  ! DECOMP -ene int sel
            backintind = 16, &  ! DECOMP  ene int ind
            backvdwsel = 17, &  ! DECOMP -ene vdw sel
            backvdwind = 18, &  ! DECOMP  ene vdw ind
            backvdwdir = 19, &  ! DECOMP  ene vdw dir
            backeelsel = 20, &  ! DECOMP -ene eel sel
            backeelind = 21, &  ! DECOMP  ene eel ind
            backeeldir = 22, &  ! DECOMP  ene eel dir
            backpolsel = 23, &  ! DECOMP -ene pol sel
            backpolind = 24, &  ! DECOMP  ene pol ind
            backpoldir = 25, &  ! DECOMP  ene pol dir
            backsassel = 26, &  ! DECOMP -ene sas sel
            backsasind = 27, &  ! DECOMP  ene sas ind
            backsasdir = 28     ! DECOMP  ene sas dir
                                               
integer, private, dimension(ndectype)       :: ndecind 
!  index into the dec array

integer, private                            :: ndecno  
!  length of the dec array minus 1

_REAL_,  private, dimension(:), allocatable :: dec     
!  stores decomposition values

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for integer arrays of module decomp
subroutine allocate_int_decomp( natom, nres )

   implicit none
   integer, intent(in) :: natom
   integer, intent(in) :: nres
   integer ier

   allocate( jgroup(natom), stat = ier )
   REQUIRE( ier == 0 )
   jgroup(:) = 0

   allocate( index(nres), stat = ier )
   REQUIRE( ier == 0 )
   index(:) = 0

   allocate( irespw(nres), stat = ier )
   REQUIRE( ier == 0 )
   irespw(:) = 0

   return

end subroutine allocate_int_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for real array of module decomp
subroutine allocate_real_decomp( nin )

   implicit none
   integer, intent(in) :: nin
   integer i, j, ncnt, n, ier

   ndecno = nin

   n = nin + 1  !  +1 for "rest" storage
   allocate( dec(ndectype * n), stat = ier )
   REQUIRE( ier == 0 )

   ! Init ndecind and dec arrays
   ncnt = 1
   do i=1,ndectype
      ndecind(i) = ncnt
      do j=1,n
         dec(ncnt) = 0.0d0
         ncnt = ncnt + 1
      end do
   end do

   return

end subroutine allocate_real_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for integer arrays of module decomp
subroutine deallocate_int_decomp( )

   implicit none
   integer ier

   if ( allocated( jgroup ) ) then
      deallocate( jgroup, stat = ier )
      REQUIRE( ier == 0 )
   else
      ASSERT( .false. )  ! cannot deallocate un-allocated array
   end if

   if ( allocated( index ) ) then
      deallocate( index, stat = ier )
      REQUIRE( ier == 0 )
   else
      ASSERT( .false. )  ! cannot deallocate un-allocated array
   end if

   if ( allocated( irespw ) ) then
      deallocate( irespw, stat = ier )
      REQUIRE( ier == 0 )
   else
      ASSERT( .false. )  ! cannot deallocate un-allocated array
   end if
   return

end subroutine deallocate_int_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for real array of module decomp
subroutine deallocate_real_decomp( )

   implicit none
   integer ier

   if ( allocated( dec ) ) then
      deallocate( dec, stat = ier )
      REQUIRE( ier == 0 )
   else
      ASSERT( .false. )  ! cannot deallocate un-allocated array
   end if

   return

end subroutine deallocate_real_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine decpair here]
subroutine decpair(nty,nat1,nat2,fval)
   
   ! Accumulates contributions fval into ENE_XXX_SEL/IND/DIR,
   !   where XXX stands for EGB, VDW, EEL or
   !   INT in the case of bond energies
   
   ! Holger Gohlke
   !   12.11.2001
   
   implicit none
   
   _REAL_  fval,hfval
   integer nty, ntype
   integer nat1, nat2
   integer nres1, nres2
   integer nind1, nind2
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir
   integer ntmp
   logical isside1, isside2
   logical isprot1, isprot2
   logical pairwise

#include "memory.h"
   ! mjhsieh: warnings eliminator
   nssel = -1; nsdir = -1; nsind = -1
   nbsel = -1; nbind = -1; nbdir = -1
   
   hfval = 0.5d0 * fval

   if(nty < 0) then
      pairwise = .true.
      ntype = -nty
   else
      pairwise = .false.
      ntype = nty
   end if

   if (ntype == 1) then
      !       --- GB decomposition
      nssel = ndecind(sidepolsel)
      nsind = ndecind(sidepolind)
      nsdir = ndecind(sidepoldir)
      nbsel = ndecind(backpolsel)
      nbind = ndecind(backpolind)
      nbdir = ndecind(backpoldir)
   else if (ntype == 2) then
      !       --- EEL decomposition
      !           (+ 1-4 EEL if idecomp == 2 or 4;
      !              this follows decomp. for mm_pbsa)
      nssel = ndecind(sideeelsel)
      nsind = ndecind(sideeelind)
      nsdir = ndecind(sideeeldir)
      nbsel = ndecind(backeelsel)
      nbind = ndecind(backeelind)
      nbdir = ndecind(backeeldir)
   else if (ntype == 3) then
      !       --- VDW decomposition
      !           (+ 1-4 VDW if idecomp == 2 or 4;
      !              this follows decomp for mm_pbsa)
      nssel = ndecind(sidevdwsel)
      nsind = ndecind(sidevdwind)
      nsdir = ndecind(sidevdwdir)
      nbsel = ndecind(backvdwsel)
      nbind = ndecind(backvdwind)
      nbdir = ndecind(backvdwdir)
   else if (ntype == 4) then
      !       --- BOND, UREY-BRADLEY decomposition
      !           (+ 1-4 EEL + 1-4 VDW if idecomp == 1 or 3;
      !              this follows decomp as done within sander)
      nssel = ndecind(sideintsel)
      nsind = ndecind(sideintind)
      nsdir = nsind                  ! "Hack" that assures that cases such
      nbsel = ndecind(backintsel)    !   as nres1 = 0, nres2 != 0 will work
      nbind = ndecind(backintind)
      nbdir = nbind                  ! dito
   else
      write(6,*) 'Wrong input for ntype: ',ntype
      call mexit(6,1)
   end if  ! (ntype == 1)
   
   ! --- Determine if side or back and if prot or lig
   
   isside1 = .true.
   isprot1 = .true.
   nres1 = jgroup(nat1)
   if(nres1 < 0) then
      isprot1 = .false.
      nres1 = -nres1
   end if
   if(nres1 > nres) then
      isside1 = .false.
      nres1 = nres1 - nres
   end if

   isside2 = .true.
   isprot2 = .true.
   nres2 = jgroup(nat2)
   if(nres2 < 0) then
      isprot2 = .false.
      nres2 = -nres2
   end if
   if(nres2 > nres) then
      isside2 = .false.
      nres2 = nres2 - nres
   end if
   
   ! --- Echo result
   
   !      write(6,*) nat1,nres1,isside1,isprot1, &
   !                 nat2,nres2,isside2,isprot2
   
   ! --- Decompose
   
   if (nres1 == nres2) then
      !       --- self-energy of residue
      if (isside1) then
         nind1 = nssel
      else
         nind1 = nbsel
      end if
      if (isside2) then
         nind2 = nssel
      else
         nind2 = nbsel
      end if
   else
      if ((     isprot1 .and.      isprot2) .or. &
            (.not.isprot1 .and. .not.isprot2)) then
         !         --- indirect energy from within own molecule
         if (isside1) then
            nind1 = nsind
         else
            nind1 = nbind
         end if
         if (isside2) then
            nind2 = nsind
         else
            nind2 = nbind
         end if
      else
         !         --- direct energy between molecules
         if (isside1) then
            nind1 = nsdir
         else
            nind1 = nbdir
         end if
         if (isside2) then
            nind2 = nsdir
         else
            nind2 = nbdir
         end if
      end if
   end if  ! (nres1 == nres2)
   
   ! --- Accumulate results
   
   if(nind1 == 0 .or. nind2 == 0) then
      write(6,*) 'NIND1 or NIND2 wrong in DECPAIR.'
      write(6,*) 'NTYPE = ',ntype
      write(6,*) nat1,nres1,isside1,isprot1, &
            nat2,nres2,isside2,isprot2
      call mexit(6,1)
   end if

   if(nres1 > 0 .and. nres2 > 0) then
      !       --- Interactions "wanted" to be considered
      if(pairwise) then
         !         --- Pairwise decomp
         ntmp = (index(nres1) - 1) * npdec + index(nres2)
         dec(nind1 + ntmp - 1) = dec(nind1 + ntmp - 1) + hfval
         ntmp = (index(nres2) - 1) * npdec + index(nres1)
         dec(nind2 + ntmp - 1) = dec(nind2 + ntmp - 1) + hfval
      else
         !         --- Residue-based decomp
         dec(nind1 + nres1 - 1) = dec(nind1 + nres1 - 1) + hfval
         dec(nind2 + nres2 - 1) = dec(nind2 + nres2 - 1) + hfval
      end if
   else
      !       --- Collect as "rest"
      dec(nind1 + ndecno) = dec(nind1 + ndecno) + hfval
      dec(nind2 + ndecno) = dec(nind2 + ndecno) + hfval
   end if
   
   return
end subroutine decpair 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine decangle here]
subroutine decangle(nat1,nat2,nat3,fval)
   
   ! Accumulates contributions fval into ENE_INT_SEL/IND,
   !   for ANGLE energies
   
   ! Holger Gohlke
   !   12.11.2001
   
   implicit none
   
   _REAL_  fval,hfval
   integer nat1, nat2, nat3
   integer nres1, nres2, nres3
   integer nind1, nind2, nind3
   logical isside1, isside2, isside3
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir

#include "memory.h"

   hfval = 0.3333333333333333334d0 * fval
   nssel = ndecind(sideintsel)
   nsind = ndecind(sideintind)
   nsdir = 0
   nbsel = ndecind(backintsel)
   nbind = ndecind(backintind)
   nbdir = 0
   
   ! --- Determine if side or back and if prot or lig
   
   isside1 = .true.
   nres1 = iabs(jgroup(nat1))
   if(nres1 > nres) then
      isside1 = .false.
      nres1 = nres1 - nres
   end if

   isside2 = .true.
   nres2 = iabs(jgroup(nat2))
   if(nres2 > nres) then
      isside2 = .false.
      nres2 = nres2 - nres
   end if

   isside3 = .true.
   nres3 = iabs(jgroup(nat3))
   if(nres3 > nres) then
      isside3 = .false.
      nres3 = nres3 - nres
   end if
   
   ! --- Echo result
   
   !      write(6,*) nat1,nres1,isside1,
   !     +           nat2,nres2,isside2,
   !     +           nat3,nres3,isside3
   
   ! --- Decompose
   
   if ((nres1 == nres2) .and. (nres1 == nres3)) then
      !       --- self-energy of residue
      if (isside1) then
         nind1 = nssel
      else
         nind1 = nbsel
      end if
      if (isside2) then
         nind2 = nssel
      else
         nind2 = nbsel
      end if
      if (isside3) then
         nind3 = nssel
      else
         nind3 = nbsel
      end if
   else
      !       --- indirect energy from within own molecule
      if (isside1) then
         nind1 = nsind
      else
         nind1 = nbind
      end if
      if (isside2) then
         nind2 = nsind
      else
         nind2 = nbind
      end if
      if (isside3) then
         nind3 = nsind
      else
         nind3 = nbind
      end if
   end if  ! ((nres1 == nres2) .and. (nres1 == nres3))
   
   ! --- Accumulate results
   
   if(nind1 == 0 .or. nind2 == 0 .or. nind3 == 0) then
      write(6,*) 'NIND1 or NIND2 or NIND3 wrong in DECPAIR.'
      write(6,*) nat1,nres1,isside1, &
            nat2,nres2,isside2, &
            nat3,nres3,isside3
      call mexit(6,1)
   end if

   if(nres1 > 0 .and. nres2 > 0 .and. nres3 > 0) then
      !       --- Interactions "wanted" to be considered
      dec(nind1 + nres1 - 1) = dec(nind1 + nres1 - 1) + hfval
      dec(nind2 + nres2 - 1) = dec(nind2 + nres2 - 1) + hfval
      dec(nind3 + nres3 - 1) = dec(nind3 + nres3 - 1) + hfval
   else
      !       --- Collect as "rest"
      dec(nind1 + ndecno) = dec(nind1 + ndecno) + hfval
      dec(nind2 + ndecno) = dec(nind2 + ndecno) + hfval
      dec(nind3 + ndecno) = dec(nind3 + ndecno) + hfval
   end if
   
   return
end subroutine decangle 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine decphi here]
subroutine decphi(nat1,nat2,nat3,nat4,fval)
   
   ! Accumulates contributions fval into ENE_INT_SEL/IND,
   !   for PHI energies
   
   ! Holger Gohlke
   !   12.11.2001
   
   implicit none
   
   _REAL_  fval,hfval
   integer nat1, nat2, nat3, nat4
   integer nres1, nres2, nres3, nres4
   integer nind1, nind2, nind3, nind4
   logical isside1, isside2, isside3, isside4
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir

#include "memory.h"

   hfval = 0.25d0 * fval
   nssel = ndecind(sideintsel)
   nsind = ndecind(sideintind)
   nsdir = 0
   nbsel = ndecind(backintsel)
   nbind = ndecind(backintind)
   nbdir = 0
   
   ! --- Determine if side or back and if prot or lig
   
   isside1 = .true.
   nres1 = iabs(jgroup(nat1))
   if(nres1 > nres) then
      isside1 = .false.
      nres1 = nres1 - nres
   end if

   isside2 = .true.
   nres2 = iabs(jgroup(nat2))
   if(nres2 > nres) then
      isside2 = .false.
      nres2 = nres2 - nres
   end if

   isside3 = .true.
   nres3 = iabs(jgroup(nat3))
   if(nres3 > nres) then
      isside3 = .false.
      nres3 = nres3 - nres
   end if

   isside4 = .true.
   nres4 = iabs(jgroup(nat4))
   if(nres4 > nres) then
      isside4 = .false.
      nres4 = nres4 - nres
   end if
   
   ! --- Echo result
   
   !      write(6,*) nat1,nres1,isside1,
   !     +           nat2,nres2,isside2,
   !     +           nat3,nres3,isside3,
   !     +           nat4,nres4,isside4
   
   ! --- Decompose
   
   if ((nres1 == nres2) .and. &
         (nres1 == nres3) .and. &
         (nres1 == nres4)) then
      !       --- self-energy of residue
      if (isside1) then
         nind1 = nssel
      else
         nind1 = nbsel
      end if
      if (isside2) then
         nind2 = nssel
      else
         nind2 = nbsel
      end if
      if (isside3) then
         nind3 = nssel
      else
         nind3 = nbsel
      end if
      if (isside4) then
         nind4 = nssel
      else
         nind4 = nbsel
      end if
   else
      !       --- indirect energy from within own molecule
      if (isside1) then
         nind1 = nsind
      else
         nind1 = nbind
      end if
      if (isside2) then
         nind2 = nsind
      else
         nind2 = nbind
      end if
      if (isside3) then
         nind3 = nsind
      else
         nind3 = nbind
      end if
      if (isside4) then
         nind4 = nsind
      else
         nind4 = nbind
      end if
   end if  !  ((nres1 == nres2) .and.
   
   ! --- Accumulate results
   
   if(nind1 == 0 .or. nind2 == 0 .or. &
         nind3 == 0 .or. nind4 == 0) then
      write(6,*) &
            'NIND1 or NIND2 or NIND3 or NIND4 wrong in DECPAIR.'
      write(6,*) nat1,nres1,isside1, &
            nat2,nres2,isside2, &
            nat3,nres3,isside3, &
            nat4,nres4,isside4
      call mexit(6,1)
   end if

   if(nres1 > 0 .and. nres2 > 0 .and. &
         nres3 > 0 .and. nres4 > 0) then
      !       --- Interactions "wanted" to be considered
      dec(nind1 + nres1 - 1) = dec(nind1 + nres1 - 1) + hfval
      dec(nind2 + nres2 - 1) = dec(nind2 + nres2 - 1) + hfval
      dec(nind3 + nres3 - 1) = dec(nind3 + nres3 - 1) + hfval
      dec(nind4 + nres4 - 1) = dec(nind4 + nres4 - 1) + hfval
   else
      !       --- Collect as "rest"
      dec(nind1 + ndecno) = dec(nind1 + ndecno) + hfval
      dec(nind2 + ndecno) = dec(nind2 + ndecno) + hfval
      dec(nind3 + ndecno) = dec(nind3 + ndecno) + hfval
      dec(nind4 + ndecno) = dec(nind4 + ndecno) + hfval
   end if
   
   return
end subroutine decphi 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine decsasa here]
subroutine decsasa(nty,nat1,nat2,nat3,fval)
   
   ! Accumulates contributions fval into ENE_SAS_SEL/IND/DIR
   
   ! Holger Gohlke
   !   12.11.2001
   
   implicit none
   
   _REAL_  fval,hfval,addfval
   integer nty, ntype
   integer nat1, nat2, nat3
   integer nres1, nres2, nres3
   integer nind1, nind2, nind3
   integer ntmp
   logical isside1, isside2, isside3
   logical isprot1, isprot2, isprot3
   logical pairwise
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir

#include "memory.h"
   ! mjhsieh: warning eliminator
   ntype = -1; nres3 = -1; isprot3 = .false.
   addfval = -1; nind3 = -1; nind1 = -1
   nres2 = -1; nind2 = -1; pairwise = .false.

   hfval = 0.5d0 * fval
   nssel = ndecind(sidesassel)
   nsind = ndecind(sidesasind)
   nsdir = ndecind(sidesasdir)
   nbsel = ndecind(backsassel)
   nbind = ndecind(backsasind)
   nbdir = ndecind(backsasdir)
   
   ! --- Check NTYPE
   
   if(nty == 0 .or. abs(nty) >= 4) then
      write(6,*) 'Wrong input for NTYPE: ',nty
      call mexit(6,1)
   else if(nty < 0) then
      pairwise = .true.
      ntype = -nty
   else
      pairwise = .false.
      ntype = nty
   end if
   
   ! --- Determine if side or back and if prot or lig
   
   isside1 = .true.
   isprot1 = .true.
   nres1 = jgroup(nat1)
   if(nres1 < 0) then
      isprot1 = .false.
      nres1 = -nres1
   end if
   if(nres1 > nres) then
      isside1 = .false.
      nres1 = nres1 - nres
   end if

   if(ntype >= 2) then
      isside2 = .true.
      isprot2 = .true.
      nres2 = jgroup(nat2)
      if(nres2 < 0) then
         isprot2 = .false.
         nres2 = -nres2
      end if
      if(nres2 > nres) then
         isside2 = .false.
         nres2 = nres2 - nres
      end if

      if(ntype == 3) then
         isside3 = .true.
         isprot3 = .true.
         nres3 = jgroup(nat3)
         if(nres3 < 0) then
            isprot3 = .false.
            nres3 = -nres3
         end if
         if(nres3 > nres) then
            isside3 = .false.
            nres3 = nres3 - nres
         end if
      end if
   end if
   
   ! --- Echo result
   
   !      if(ntype.eq.1) then
   !        write(6,*) nat1,nres1,isside1, isprot1, fval
   !      else if(ntype.eq.2) then
   !        write(6,*) nat1,nres1,isside1, isprot1, &
   !                   nat2,nres2,isside2, isprot2, &
   !                   fval
   !      else if(ntype.eq.3) then
   !        write(6,*) nat1,nres1,isside1, isprot1, &
   !                   nat2,nres2,isside2, isprot2, &
   !                   nat3,nres3,isside3, isprot3, &
   !                   fval
   !      end if
   
   ! --- Decompose
   !     here: nat1 solely determines if sidechain or backbone

   if(ntype == 1) then
      !       --- self-energy
      addfval = fval
      if (isside1) then
         nind1 = nssel
      else
         nind1 = nbsel
      end if
   else if(ntype >= 2) then
      if(ntype == 2) then
         addfval = fval
      else if(ntype == 3) then
         addfval = hfval
      end if

      if(nres2 == nres1) then
         !         --- self-energy
         if (isside1) then
            nind2 = nssel
         else
            nind2 = nbsel
         end if
      else if ((     isprot1 .and.      isprot2) .or. &
               (.not.isprot1 .and. .not.isprot2)) then
         !         --- indirect energy
         if (isside1) then
            nind2 = nsind
         else
            nind2 = nbind
         end if
      else
         !         --- direct energy
         if (isside1) then
            nind2 = nsdir
         else
            nind2 = nbdir
         end if
      end if

      if(ntype == 3) then
         if(nres3 == nres1) then
            !           --- self-energy
            if (isside1) then
               nind3 = nssel
            else
               nind3 = nbsel
            end if
         else if ((     isprot1 .and.      isprot3) .or. &
                  (.not.isprot1 .and. .not.isprot3)) then
            !           --- indirect energy
            if (isside1) then
               nind3 = nsind
            else
               nind3 = nbind
            end if
         else
            !           --- direct energy
            if (isside1) then
               nind3 = nsdir
            else
               nind3 = nbdir
            end if
         end if
      end if
      
   end if
   
   ! --- Accumulate results
   
   if((ntype == 1 .and. nres1 > 0) .or. &
      (ntype == 2 .and. nres1 > 0 .and. nres2 > 0) .or. &
      (ntype == 3 .and. &
        nres1 > 0 .and. nres2 > 0 .and. nres3 > 0)) then
      !       --- Interactions "wanted" to be considered
      !           here: all contributions are collected for nres1(!)
      if(ntype == 1) then
         if(pairwise) then
            ntmp = (index(nres1) - 1) * npdec + index(nres1)
            dec(nind1 + ntmp - 1) = dec(nind1 + ntmp - 1) + addfval
         else
            dec(nind1 + nres1 - 1) = dec(nind1 + nres1 - 1) + addfval
         end if
      else if(ntype >= 2) then
         if(pairwise) then
            ntmp = (index(nres1) - 1) * npdec + index(nres2)
            dec(nind2 + ntmp - 1) = dec(nind2 + ntmp - 1) + addfval
         else
            dec(nind2 + nres1 - 1) = dec(nind2 + nres1 - 1) + addfval
         end if
         if(ntype == 3) then
            if(pairwise) then
               ntmp = (index(nres1) - 1) * npdec + index(nres3)
               dec(nind3 + ntmp - 1) = dec(nind3 + ntmp - 1) + addfval
            else
               dec(nind3 + nres1 - 1) = dec(nind3 + nres1 - 1) + addfval
            end if
         end if
      end if
   else
      !       --- Collect as "rest"
      if(ntype == 1) then
         dec(nind1 + ndecno) = dec(nind1 + ndecno) + addfval
      else if(ntype >= 2) then
         dec(nind2 + ndecno) = dec(nind2 + ndecno) + addfval
         if(ntype == 3) then
            dec(nind3 + ndecno) = dec(nind3 + ndecno) + addfval
         end if
      end if
   end if  ! ((ntype == 1 .and. nres1 > 0) .or.
   
   return
end subroutine decsasa 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine checkdec here]
subroutine checkdec(idecomp)
   
   ! Sums contents in ENE_XXX_SEL/IND/DIR
   
   ! Holger Gohlke
   !   13.11.2001
   
   implicit none
   
   integer, intent(in) :: idecomp

   _REAL_ &
         esintsel, esintind, esvdwsel, esvdwind, esvdwdir, &
         eseelsel, eseelind, eseeldir, espolsel, espolind, &
         espoldir, essassel, essasind, essasdir
   _REAL_ &
         ebintsel, ebintind, ebvdwsel, ebvdwind, ebvdwdir, &
         ebeelsel, ebeelind, ebeeldir, ebpolsel, ebpolind, &
         ebpoldir, ebsassel, ebsasind, ebsasdir

   integer i, ndec
   integer nsintsel, nsintind, nsvdwsel, nsvdwind, nsvdwdir, &
         nseelsel, nseelind, nseeldir, nspolsel, nspolind, &
         nspoldir, nssassel, nssasind, nssasdir
   integer nbintsel, nbintind, nbvdwsel, nbvdwind, nbvdwdir, &
         nbeelsel, nbeelind, nbeeldir, nbpolsel, nbpolind, &
         nbpoldir, nbsassel, nbsasind, nbsasdir

#include "memory.h"
   
   nsintsel = ndecind(sideintsel)
   nsintind = ndecind(sideintind)
   nsvdwsel = ndecind(sidevdwsel)
   nsvdwind = ndecind(sidevdwind)
   nsvdwdir = ndecind(sidevdwdir)
   nseelsel = ndecind(sideeelsel)
   nseelind = ndecind(sideeelind)
   nseeldir = ndecind(sideeeldir)
   nspolsel = ndecind(sidepolsel)
   nspolind = ndecind(sidepolind)
   nspoldir = ndecind(sidepoldir)
   nssassel = ndecind(sidesassel)
   nssasind = ndecind(sidesasind)
   nssasdir = ndecind(sidesasdir)
                                 
   nbintsel = ndecind(backintsel)
   nbintind = ndecind(backintind)
   nbvdwsel = ndecind(backvdwsel)
   nbvdwind = ndecind(backvdwind)
   nbvdwdir = ndecind(backvdwdir)
   nbeelsel = ndecind(backeelsel)
   nbeelind = ndecind(backeelind)
   nbeeldir = ndecind(backeeldir)
   nbpolsel = ndecind(backpolsel)
   nbpolind = ndecind(backpolind)
   nbpoldir = ndecind(backpoldir)
   nbsassel = ndecind(backsassel)
   nbsasind = ndecind(backsasind)
   nbsasdir = ndecind(backsasdir)
   
   esintsel = 0.0d0
   esintind = 0.0d0
   esvdwsel = 0.0d0
   esvdwind = 0.0d0
   esvdwdir = 0.0d0
   eseelsel = 0.0d0
   eseelind = 0.0d0
   eseeldir = 0.0d0
   espolsel = 0.0d0
   espolind = 0.0d0
   espoldir = 0.0d0
   essassel = 0.0d0
   essasind = 0.0d0
   essasdir = 0.0d0

   ebintsel = 0.0d0
   ebintind = 0.0d0
   ebvdwsel = 0.0d0
   ebvdwind = 0.0d0
   ebvdwdir = 0.0d0
   ebeelsel = 0.0d0
   ebeelind = 0.0d0
   ebeeldir = 0.0d0
   ebpolsel = 0.0d0
   ebpolind = 0.0d0
   ebpoldir = 0.0d0
   ebsassel = 0.0d0
   ebsasind = 0.0d0
   ebsasdir = 0.0d0

   if(idecomp < 3) then
      ndec = nres
   else if(idecomp > 2) then
      ndec = npdec*npdec
   endif

   do i=0,ndec-1  ! here: not till nres -> keeps "rest" apart
      esintsel = esintsel + dec(nsintsel + i)
      esintind = esintind + dec(nsintind + i)
      esvdwsel = esvdwsel + dec(nsvdwsel + i)
      esvdwind = esvdwind + dec(nsvdwind + i)
      esvdwdir = esvdwdir + dec(nsvdwdir + i)
      eseelsel = eseelsel + dec(nseelsel + i)
      eseelind = eseelind + dec(nseelind + i)
      eseeldir = eseeldir + dec(nseeldir + i)
      espolsel = espolsel + dec(nspolsel + i)
      espolind = espolind + dec(nspolind + i)
      espoldir = espoldir + dec(nspoldir + i)
      essassel = essassel + dec(nssassel + i)
      essasind = essasind + dec(nssasind + i)
      essasdir = essasdir + dec(nssasdir + i)

      ebintsel = ebintsel + dec(nbintsel + i)
      ebintind = ebintind + dec(nbintind + i)
      ebvdwsel = ebvdwsel + dec(nbvdwsel + i)
      ebvdwind = ebvdwind + dec(nbvdwind + i)
      ebvdwdir = ebvdwdir + dec(nbvdwdir + i)
      ebeelsel = ebeelsel + dec(nbeelsel + i)
      ebeelind = ebeelind + dec(nbeelind + i)
      ebeeldir = ebeeldir + dec(nbeeldir + i)
      ebpolsel = ebpolsel + dec(nbpolsel + i)
      ebpolind = ebpolind + dec(nbpolind + i)
      ebpoldir = ebpoldir + dec(nbpoldir + i)
      ebsassel = ebsassel + dec(nbsassel + i)
      ebsasind = ebsasind + dec(nbsasind + i)
      ebsasdir = ebsasdir + dec(nbsasdir + i)
   end do
   
   ! --- Output total energies (w/ rest)
   
   write(6,300)
   write(6,350) esintsel + esintind + ebintsel + ebintind + &
                dec(nsintsel+ndecno) + dec(nsintind+ndecno) + &
                dec(nbintsel+ndecno) + dec(nbintind+ndecno)
   write(6,360) esvdwsel + esvdwind + esvdwdir + &
                ebvdwsel + ebvdwind + ebvdwdir + &
                dec(nsvdwsel+ndecno) + dec(nsvdwind+ndecno) + &
                dec(nsvdwdir+ndecno) + &
                dec(nbvdwsel+ndecno) + dec(nbvdwind+ndecno) + &
                dec(nbvdwdir+ndecno), &
                eseelsel + eseelind + eseeldir + &
                ebeelsel + ebeelind + ebeeldir + &
                dec(nseelsel+ndecno) + dec(nseelind+ndecno) + &
                dec(nseeldir+ndecno) + &
                dec(nbeelsel+ndecno) + dec(nbeelind+ndecno) + &
                dec(nbeeldir+ndecno)
   write(6,370) espolsel + espolind + espoldir + &
                ebpolsel + ebpolind + ebpoldir + &
                dec(nspolsel+ndecno) + dec(nspolind+ndecno) + &
                dec(nspoldir+ndecno) + &
                dec(nbpolsel+ndecno) + dec(nbpolind+ndecno) + &
                dec(nbpoldir+ndecno), &
                essassel + essasind + essasdir + &
                ebsassel + ebsasind + ebsasdir + &
                dec(nssassel+ndecno) + dec(nssasind+ndecno) + &
                dec(nssasdir+ndecno) + &
                dec(nbsassel+ndecno) + dec(nbsasind+ndecno) + &
                dec(nbsasdir+ndecno)
   
   ! --- Output self energies (w/o rest)
   
   write(6,302)
   write(6,350) esintsel + ebintsel
   write(6,360) esvdwsel + ebvdwsel, &
                eseelsel + ebeelsel
   write(6,370) espolsel + ebpolsel, &
                essassel + ebsassel
   
   ! --- Output indirect energies (w/o rest)
   
   write(6,304)
   write(6,350) esintind + ebintind
   write(6,360) esvdwind + ebvdwind, &
                eseelind + ebeelind
   write(6,370) espolind + ebpolind, &
                essasind + ebsasind
   
   ! --- Output direct energies (w/o rest)
   
   write(6,306)
   write(6,350) 0.0d0    ! No direct interactions for internal energies
   write(6,360) esvdwdir + ebvdwdir, &
                eseeldir + ebeeldir
   write(6,370) espoldir + ebpoldir, &
                essasdir + ebsasdir
   
   ! --- Output rest energies
   
   write(6,308)
   write(6,350) dec(nsintsel+ndecno) + dec(nsintind+ndecno) + &
                dec(nbintsel+ndecno) + dec(nbintind+ndecno)
   write(6,360) dec(nsvdwsel+ndecno) + dec(nsvdwind+ndecno) + &
                dec(nsvdwdir+ndecno) + &
                dec(nbvdwsel+ndecno) + dec(nbvdwind+ndecno) + &
                dec(nbvdwdir+ndecno), &
                dec(nseelsel+ndecno) + dec(nseelind+ndecno) + &
                dec(nseeldir+ndecno) + &
                dec(nbeelsel+ndecno) + dec(nbeelind+ndecno) + &
                dec(nbeeldir+ndecno)
   write(6,370) dec(nspolsel+ndecno) + dec(nspolind+ndecno) + &
                dec(nspoldir+ndecno) + &
                dec(nbpolsel+ndecno) + dec(nbpolind+ndecno) + &
                dec(nbpoldir+ndecno), &
                dec(nssassel+ndecno) + dec(nssasind+ndecno) + &
                dec(nssasdir+ndecno) + &
                dec(nbsassel+ndecno) + dec(nbsasind+ndecno) + &
                dec(nbsasdir+ndecno)
   
   write(6,'(/)')
   
   300 format(/ /20x,'CHECK DECOMP - TOTAL ENERGIES (w/ REST)',/)
   302 format(/ /20x,'CHECK DECOMP - SELF ENERGIES (w/o REST)',/)
   304 format(/ /20x,'CHECK DECOMP - INDIRECT ENERGIES (w/o REST)',/)
   306 format(/ /20x,'CHECK DECOMP - DIRECT ENERGIES (w/o REST)',/)
   308 format(/ /20x,'CHECK DECOMP - REST ENERGIES',/)
   350 format(1x,'INTERNAL= ',f13.4)
   360 format(1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4)
   370 format(1x,'EGB     = ',f13.4,2x,'ESURF   = ',f13.4)
   
   return
end subroutine checkdec 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine printdec here]
subroutine printdec(ix)
   
   ! Outputs contents in ENE_XXX_SEL/IND/DIR with respect to residues
   
   ! Holger Gohlke
   !   13.11.2001
   
   implicit none
   
   _REAL_  energy
   dimension energy(5)

   integer i, j, ipresst, igst, istart
   integer ix(*)
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir
   dimension nssel(5), nsind(5), nsdir(5)
   dimension nbsel(5), nbind(5), nbdir(5)
   
#  include "memory.h"
   
   ipresst  = i02
   igst     = ibellygp

   nssel(1) = ndecind(sideintsel)
   nsind(1) = ndecind(sideintind)
   nsdir(1) = 0                       
   nssel(2) = ndecind(sidevdwsel)
   nsind(2) = ndecind(sidevdwind)
   nsdir(2) = ndecind(sidevdwdir)
   nssel(3) = ndecind(sideeelsel)
   nsind(3) = ndecind(sideeelind)
   nsdir(3) = ndecind(sideeeldir)
   nssel(4) = ndecind(sidepolsel)
   nsind(4) = ndecind(sidepolind)
   nsdir(4) = ndecind(sidepoldir)
   nssel(5) = ndecind(sidesassel)
   nsind(5) = ndecind(sidesasind)
   nsdir(5) = ndecind(sidesasdir)
                                      
   nbsel(1) = ndecind(backintsel)
   nbind(1) = ndecind(backintind)
   nbdir(1) = 0                       
   nbsel(2) = ndecind(backvdwsel)
   nbind(2) = ndecind(backvdwind)
   nbdir(2) = ndecind(backvdwdir)
   nbsel(3) = ndecind(backeelsel)
   nbind(3) = ndecind(backeelind)
   nbdir(3) = ndecind(backeeldir)
   nbsel(4) = ndecind(backpolsel)
   nbind(4) = ndecind(backpolind)
   nbdir(4) = ndecind(backpoldir)
   nbsel(5) = ndecind(backsassel)
   nbind(5) = ndecind(backsasind)
   nbdir(5) = ndecind(backsasdir)
   
   ! --- Output total energies
   
   write(6,300)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0) then
         do j=1,5
            energy(j) = 0.0d0
            if(nssel(j) > 0 .and. nbsel(j) > 0) &
                  energy(j) = energy(j) + dec(nssel(j)+i) + dec(nbsel(j)+i)
            if(nsind(j) > 0 .and. nbind(j) > 0) &
                  energy(j) = energy(j) + dec(nsind(j)+i) + dec(nbind(j)+i)
            if(nsdir(j) > 0 .and. nbdir(j) > 0) &
                  energy(j) = energy(j) + dec(nsdir(j)+i) + dec(nbdir(j)+i)
         end do
         write(6,360) i+1,(energy(j),j=1,5)
      end if
   end do
   
   ! --- Output sidechain energies
   
   write(6,302)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0) then
         do j=1,5
            energy(j) = 0.0d0
            if(nssel(j) > 0) &
                  energy(j) = energy(j) + dec(nssel(j)+i)
            if(nsind(j) > 0) &
                  energy(j) = energy(j) + dec(nsind(j)+i)
            if(nsdir(j) > 0) &
                  energy(j) = energy(j) + dec(nsdir(j)+i)
         end do
         write(6,362) i+1,(energy(j),j=1,5)
      end if
   end do
   
   ! --- Output backbone energies
   
   write(6,304)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0) then
         do j=1,5
            energy(j) = 0.0d0
            if(nbsel(j) > 0) &
                  energy(j) = energy(j) + dec(nbsel(j)+i)
            if(nbind(j) > 0) &
                  energy(j) = energy(j) + dec(nbind(j)+i)
            if(nbdir(j) > 0) &
                  energy(j) = energy(j) + dec(nbdir(j)+i)
         end do
         write(6,364) i+1,(energy(j),j=1,5)
      end if
   end do
   
   write(6,'(/)')
   
   300 format(/ /20x,'PRINT DECOMP - TOTAL ENERGIES',/)
   302 format(/ /20x,'PRINT DECOMP - SIDECHAIN ENERGIES',/)
   304 format(/ /20x,'PRINT DECOMP - BACKBONE ENERGIES',/)
   350 format(3x,1x,'resid |internal |vdw      |eel      |pol      |sas')
   355 format(6('=========='))
   360 format('TDC',1x,i6,5(1x,f9.3))
   362 format('SDC',1x,i6,5(1x,f9.3))
   364 format('BDC',1x,i6,5(1x,f9.3))
   
   return
end subroutine printdec 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine printpdec here]
subroutine printpdec(ix)
   
   ! Outputs contents in ENE_XXX_SEL/IND/DIR with
   !   respect to pairwise decomp
   
   ! Holger Gohlke
   !   29.11.2001
   
   implicit none
   
   _REAL_  energy
   dimension energy(5)

   integer i, j, k, ipresst, igst, istart
   integer ix(*)
   integer nssel, nsind, nsdir
   integer nbsel, nbind, nbdir
   integer ntmp
   dimension nssel(5), nsind(5), nsdir(5)
   dimension nbsel(5), nbind(5), nbdir(5)
   
#  include "memory.h"
   
   ipresst  = i02
   igst     = ibellygp

   nssel(1) = ndecind(sideintsel)
   nsind(1) = ndecind(sideintind)
   nsdir(1) = 0                       
   nssel(2) = ndecind(sidevdwsel)
   nsind(2) = ndecind(sidevdwind)
   nsdir(2) = ndecind(sidevdwdir)
   nssel(3) = ndecind(sideeelsel)
   nsind(3) = ndecind(sideeelind)
   nsdir(3) = ndecind(sideeeldir)
   nssel(4) = ndecind(sidepolsel)
   nsind(4) = ndecind(sidepolind)
   nsdir(4) = ndecind(sidepoldir)
   nssel(5) = ndecind(sidesassel)
   nsind(5) = ndecind(sidesasind)
   nsdir(5) = ndecind(sidesasdir)
                                      
   nbsel(1) = ndecind(backintsel)
   nbind(1) = ndecind(backintind)
   nbdir(1) = 0                       
   nbsel(2) = ndecind(backvdwsel)
   nbind(2) = ndecind(backvdwind)
   nbdir(2) = ndecind(backvdwdir)
   nbsel(3) = ndecind(backeelsel)
   nbind(3) = ndecind(backeelind)
   nbdir(3) = ndecind(backeeldir)
   nbsel(4) = ndecind(backpolsel)
   nbind(4) = ndecind(backpolind)
   nbdir(4) = ndecind(backpoldir)
   nbsel(5) = ndecind(backsassel)
   nbind(5) = ndecind(backsasind)
   nbdir(5) = ndecind(backsasdir)

   ! --- Output total energies
   
print *, dec(186)
   write(6,300)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0 .and. &
            index(i + 1) > 0) then
         do k=0,nres-1
            istart = ix(ipresst + k)
            if(ix(igst + istart - 1) > 0 .and. index(k + 1) > 0) then
               ntmp = (index(i + 1) - 1) * npdec + index(k + 1) - 1
               do j=1,5
                  energy(j) = 0.0d0
                  if(nssel(j) > 0 .and. nbsel(j) > 0) &
                        energy(j) = energy(j) + &
                        dec(nssel(j)+ntmp) + dec(nbsel(j)+ntmp)
                  if(nsind(j) > 0 .and. nbind(j) > 0) &
                        energy(j) = energy(j) + &
                        dec(nsind(j)+ntmp) + dec(nbind(j)+ntmp)
                  if(nsdir(j) > 0 .and. nbdir(j) > 0) &
                        energy(j) = energy(j) + &
                        dec(nsdir(j)+ntmp) + dec(nbdir(j)+ntmp)
               end do
               write(6,360) i+1,k+1,(energy(j),j=1,5)
            end if
         end do
      end if
   end do
   
   ! --- Output sidechain energies
   
   write(6,302)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0 .and. index(i + 1) > 0) then
         do k=0,nres-1
            istart = ix(ipresst + k)
            if(ix(igst + istart - 1) > 0 .and. index(k + 1) > 0) then
               ntmp = (index(i + 1) - 1) * npdec + index(k + 1) - 1
               do j=1,5
                  energy(j) = 0.0d0
                  if(nssel(j) > 0) &
                        energy(j) = energy(j) + dec(nssel(j)+ntmp)
                  if(nsind(j) > 0) &
                        energy(j) = energy(j) + dec(nsind(j)+ntmp)
                  if(nsdir(j) > 0) &
                        energy(j) = energy(j) + dec(nsdir(j)+ntmp)
               end do
               write(6,362) i+1,k+1,(energy(j),j=1,5)
            end if
         end do
      end if
   end do
   
   ! --- Output backbone energies
   
   write(6,304)
   write(6,350)
   write(6,355)
   do i=0,nres-1
      istart = ix(ipresst + i)
      if(ix(igst + istart - 1) > 0 .and. index(i + 1) > 0) then
         do k=0,nres-1
            istart = ix(ipresst + k)
            if(ix(igst + istart - 1) > 0 .and. index(k + 1) > 0) then
               ntmp = (index(i + 1) - 1) * npdec + index(k + 1) - 1
               do j=1,5
                  energy(j) = 0.0d0
                  if(nbsel(j) > 0) &
                        energy(j) = energy(j) + dec(nbsel(j)+ntmp)
                  if(nbind(j) > 0) &
                        energy(j) = energy(j) + dec(nbind(j)+ntmp)
                  if(nbdir(j) > 0) &
                        energy(j) = energy(j) + dec(nbdir(j)+ntmp)
               end do
               write(6,364) i+1,k+1,(energy(j),j=1,5)
            end if
         end do
      end if
   end do
   
   write(6,'(/)')
   
   300 format(/ /20x,'PRINT PAIR DECOMP - TOTAL ENERGIES',/)
   302 format(/ /20x,'PRINT PAIR DECOMP - SIDECHAIN ENERGIES',/)
   304 format(/ /20x,'PRINT PAIR DECOMP - BACKBONE ENERGIES',/)
   350 format(3x,1x,'resid1 ->resid2 ', &
         '|internal |vdw      |eel      |pol      |sas')
   355 format(7('=========='))
   360 format('TDC',1x,i7,'->',i7,5(1x,f9.3))
   362 format('SDC',1x,i7,'->',i7,5(1x,f9.3))
   364 format('BDC',1x,i7,'->',i7,5(1x,f9.3))
   
   return
end subroutine printpdec 

end module decomp
