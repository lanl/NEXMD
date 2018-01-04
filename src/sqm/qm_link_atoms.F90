! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!Written by Ross Walker & Mike Crowley (TSRI, 2005)

subroutine identify_link_atoms(qmmm_nml,qmmm_struct, nbona,ib,jb)
!Link atom methods in Amber:
!
! 3) In this method the link atom sees the MM field. Here to avoid overcounting the 
!    MM link pair atoms are zeroed in the main charge array. Thus MM-MM interactions both
!    direct and reciprocal are between all MM atoms that are NOT part of a link pair.
!    QM-MM interactions are between ALL QM atoms and all non-link pair MM atoms. In the reciprocal
!    calculation QM-QM interactions are between ALL QM atoms. QM-MM interactions are between 
!    ALL QM atoms and all non-link pair MM atoms.
!
! Current link atom restrictions:
!
! Link atoms are not placed between bonds to hydrogen. So cutting across a C-H bond will not give
! you a link atom across that bond. This is not presently tested for. It is just hoped that the
! user will not try such a non-sensicle QM-MM boundary.
!
! Link atoms are restricted to one per MM link pair atom. This is tested for during the detection
! of link atoms and an error is generated if this requirement is violated. This would seem to be
! a sensible policy otherwise you could have two link atoms too close together. E.g. consider a 
! ring:
!
!  C---C              C---C
! /     \            /     \
!C       C--- ----> C       C*--- Where c* is the MM atom and C are all QM. This would give
! \     /            \     /      two link atoms very close together and is thus not allowed.
!  C---C              C---C
!
!
! Since there are no resp charges available for link atoms qmgb=1 does not work with link atoms.
! if link atoms are detected below and qmgb=1 then the program will quit.

  use qmmm_struct_module, only : new, qmmm_struct_type
  use qmmm_nml_module   , only : qmmm_nml_type

  implicit none

!passed in
  type(qmmm_nml_type), intent(inout) :: qmmm_nml
  type(qmmm_struct_type), intent(inout) :: qmmm_struct
  integer, intent(in) :: nbona,ib(nbona),jb(nbona)
!locals
  integer mm, i, j, kk, imm
  integer ier
  logical ii,jj

!     This subroutine identifies link atoms which are placed along bonds
!     connecting quantum mechanical and molecular mechanical atoms.
!     From this routine we need the number of link atoms and the pairs
!     of atoms they are between.

  qmmm_struct%nlink = 0
  do i=1,nbona !Loop over bonds without H and see if one atom is QM and one MM
      ii = .false.
      jj = .false.
      do j=1,qmmm_struct%nquant
         ii=ii .or. ((ib(i)/3 + 1).eq.qmmm_struct%iqmatoms(j))
         jj=jj .or. ((jb(i)/3 + 1).eq.qmmm_struct%iqmatoms(j))
      end do
      if (ii .neqv. jj) then
         qmmm_struct%nlink = qmmm_struct%nlink+1
      end if
  end do

  !We now know how many link atoms there are and have stored it in qmmm_struct%nlink
  ier=0
  if (qmmm_struct%nlink > 0 ) then
    allocate ( qmmm_struct%link_pairs(2,qmmm_struct%nlink), stat=ier )
    !3,x position is filled on every call to qm_mm when we fill the pair list.
    REQUIRE(ier == 0)
  
    !now we redo the loop to fill our link_pairs array with the MM-QM bonds that were broken for link atoms
    mm = 1
    do i=1,nbona
       ii = .false.
       jj = .false.
       do j=1,qmmm_struct%nquant
         if ((ib(i)/3 + 1).eq.qmmm_struct%iqmatoms(j)) then
            kk=j
            ii = .true.
         end if
         if ((jb(i)/3 + 1).eq.qmmm_struct%iqmatoms(j)) then
            kk=j
            jj = .true.
         end if 
       enddo
       if(ii .neqv. jj) then ! we are at a MM-QM boundary bond
           if(ii) then
             qmmm_struct%link_pairs(1,mm)=(jb(i)/3 + 1)
           else
             qmmm_struct%link_pairs(1,mm)=(ib(i)/3 + 1)
           end if
           qmmm_struct%link_pairs(2,mm)=kk
           mm = mm + 1
       end if
    enddo

    !EGB assumes that the link_pairs list is sorted on the first index which is
    !the MM link pair id.
    call qmsort_link_pairs(qmmm_struct%link_pairs, qmmm_struct%nlink) 
  end if

  ! We now know how many link atoms there are so re-allocate our qmmm_struct%iqm_atomic_numbers to include them
  qmmm_struct%nquant_nlink = qmmm_struct%nquant + qmmm_struct%nlink !for speed to save subroutines doing the add.

  qmmm_struct%mm_link_mask = .false. !Note, sets entire natom long array to false

  !Check that our link atoms satisfy the requirements.
  !First test - link atoms must be in numerical order.
  !Second test - MM link pairs must be unique.
  do i = 1, qmmm_struct%nlink-1
    if (qmmm_struct%link_pairs(1,i) >= qmmm_struct%link_pairs(1,i+1)) then
      write(6,*) 'QMMM: ERROR with link atom ',i
      write(6,*) 'QMMM: MM id is: ',qmmm_struct%link_pairs(1,i)
      if (i<qmmm_struct%nlink) write(6,*) 'QMMM: id of next MM atom is: ',qmmm_struct%link_pairs(1,i+1)
      call sander_bomb('identify_link_atoms','Illegal Link Atom Specification.', &
     'either a MM atom has been substituted by more than one link atom or the list of &
     &MM link pairs is not in numerical order. Cannot continue.')
    end if
  end do

  if (qmmm_struct%nlink /= 0 ) then
     ! Raellocate iqmatoms and iqm_atomic_numbers arrays to make them nquant_nlink long
     call new(qmmm_struct, qmmm_nml%qmmm_int, qmmm_nml%qmmm_switch)
     ! Now fill the link atoms atomic numbers
     do i = 1, qmmm_struct%nlink
        qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant+i) = qmmm_nml%lnk_atomic_no
     end do
     ! Fill the link atom section of iqmatoms with the mm link pair numbers
     ! At the same time we will set the mm_link_mask to true for mm link pair atoms
     do i = 1, qmmm_struct%nlink
        imm = qmmm_struct%link_pairs(1,i)
        qmmm_struct%iqmatoms(qmmm_struct%nquant+i) = imm
        qmmm_struct%mm_link_mask(imm) = .true. !This sets the MM atom flag to true to show it is part of an MM-QM link pair.
     end do 
  end if

!LINK ATOM RESTRICTION - QMGB=1 does not work with link atoms.
  if (qmmm_struct%nlink>0 .and. qmmm_nml%qmgb == 1) &
     call sander_bomb('identify_link_atoms','Link atoms are not currently available for qm_gb=1', &
                        'No resp charges for link atoms - abort')

  return
end subroutine identify_link_atoms

!------------------------------------------------------------------------

subroutine position_link_atoms(qmmm_nml,qmmm_struct, unimaged_coords)
!This routine positions link atoms a distance of lnk_dis along the bond
!vector connecting the two atoms.

!If lnk_dis < 0.0d0 then the link atom is placed on top of the MM link pair atom.

!At this point qm_coords from the qmmm_structure should contain the real QM atom
!coordinates. We will place the link atom coordinates at the end of this array.
!Hence it must have been allocated as 3*nquant+nlink long.

  use qmmm_nml_module   , only : qmmm_nml_type
   use qmmm_struct_module, only : qmmm_struct_type
   implicit none
!Passed in
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    type(qmmm_struct_type), intent(inout) :: qmmm_struct
   _REAL_, intent(in) :: unimaged_coords(3,*)

!Local
   _REAL_ :: xij, yij, zij, rij2, onerij, qmx, qmy, qmz
   integer :: i, offset, im, iq, it

   offset = qmmm_struct%nquant
   if (qmmm_nml%lnk_dis < 0.0d0) then
     !The link atom will simply be given the
     !coordinates of the MM_link_pair atom.
     do i = 1, qmmm_struct%nlink
        !Get atom identities
        im = qmmm_struct%link_pairs(1,i) !Position of mm atom in main amber coords array
        it = qmmm_struct%link_pairs(2,i)
        iq = qmmm_struct%iqmatoms(it)
        qmx = unimaged_coords(1,iq)
        qmy = unimaged_coords(2,iq)
        qmz = unimaged_coords(3,iq)
        xij = unimaged_coords(1,im) - qmx
        yij = unimaged_coords(2,im) - qmy
        zij = unimaged_coords(3,im) - qmz
        offset = offset+1
        qmmm_struct%qm_coords(1,offset) = xij + qmmm_struct%qm_coords(1,it)
        qmmm_struct%qm_coords(2,offset) = yij + qmmm_struct%qm_coords(2,it)
        qmmm_struct%qm_coords(3,offset) = zij + qmmm_struct%qm_coords(3,it)
     end do 
   else
     !The link atom will get placed at lnk_dis along the bond vector.
     do i=1,qmmm_struct%nlink
        !Get atom identities
        im = qmmm_struct%link_pairs(1,i) !Position of mm atom in main amber coords array
        it = qmmm_struct%link_pairs(2,i)
        iq = qmmm_struct%iqmatoms(it)
        qmx = unimaged_coords(1,iq)
        qmy = unimaged_coords(2,iq)
        qmz = unimaged_coords(3,iq)
        xij = unimaged_coords(1,im) - qmx
        yij = unimaged_coords(2,im) - qmy
        zij = unimaged_coords(3,im) - qmz
        rij2 = xij*xij + yij*yij + zij*zij
        onerij = 1.0d0/sqrt(rij2) * qmmm_nml%lnk_dis
        offset = offset+1
        qmmm_struct%qm_coords(1,offset) = xij * onerij + qmmm_struct%qm_coords(1,it)
        qmmm_struct%qm_coords(2,offset) = yij * onerij + qmmm_struct%qm_coords(2,it)
        qmmm_struct%qm_coords(3,offset) = zij * onerij + qmmm_struct%qm_coords(3,it)
     end do 
   end if

   return
end subroutine position_link_atoms

!------------------------------------------------------------------------

subroutine distribute_lnk_f(forcemod,flink,mmcoord,unimaged_qmcoord,link_distance)
!Written by Ross Walker (TSRI, 2005)
!This subroutine uses the chain rule to distribute the force
!that was on the link atom flink(x,y,z) between the MM and QM
!pair that the link atom sat between.

!The re-calculated force is returned in forcemod. The new QM atom force
!should then be:

!FQM(x,y,z)=FQM(x,y,z)+Flink(x,y,z)-FORCEMOD

!On the MM atom it should be:

!FMM(x,y,z)=FMM(x,y,z)+FORCEMOD

!Note in the case of lnk_dis < 0.0d0 then the link atom distance
!is no longer constant and so the differential is wrong since
!d(L-QM) is no longer a constant. However, the constraint is
!that it always has to sit on top of the MM atom thus the
!force modification is simply that the force on the link
!atom should just be placed entirely on the MMlink pair atom.

     implicit none

!Passed in
     _REAL_, intent(out) :: forcemod(3)
     _REAL_, intent(in) :: flink(3) !Force on link atom as calculated in QM calc
     _REAL_, intent(in) :: mmcoord(3) !Coordinates of MM atom in link pair (Unimaged from x array)
     _REAL_, intent(in) :: unimaged_qmcoord(3) !Unimaged coordinate of qm atom (from x array)
     _REAL_, intent(in) :: link_distance !The link atom distance from the QM atom as it was placed 
                                         !along the QM-MM bond vector (in angstroms).

!Local
     _REAL_ :: R2, oneR, dotprod, vec(3)
     _REAL_ :: lnk_dis_oneR

!If the link_distance is < 0.0d0 then the link atom was placed on top of the MM link pair atom.
     if (link_distance < 0.0d0) then
       forcemod(1:3) = flink(1:3)
     else
       vec(1:3) = mmcoord(1:3) - unimaged_qmcoord(1:3)

       R2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)

       oneR = 1.0d0 / sqrt(R2)
  
       lnk_dis_oneR = link_distance*oneR

       vec(1:3) = vec(1:3)*oneR

       dotprod = flink(1)*vec(1)+flink(2)*vec(2)+flink(3)*vec(3)

       forcemod(1) = lnk_dis_oneR*(flink(1)-(dotprod*vec(1))) 
       forcemod(2) = lnk_dis_oneR*(flink(2)-(dotprod*vec(2))) 
       forcemod(3) = lnk_dis_oneR*(flink(3)-(dotprod*vec(3))) 
     end if

     return
end subroutine distribute_lnk_f

!------------------------------------------------------------------------

subroutine print_link_atom_info(qmmm_nml,qmmm_struct, qmcoords, atom_type )

!Writes out link atom info to mdout file

  use qmmm_nml_module   , only : qmmm_nml_type
  use qmmm_struct_module, only : qmmm_struct_type
  use constants, only : INV_AMBER_ELECTROSTATIC
  implicit none

!Passed in
  type(qmmm_nml_type), intent(in) :: qmmm_nml
  type(qmmm_struct_type), intent(in) :: qmmm_struct
  _REAL_ qmcoords(3,qmmm_struct%nquant_nlink)
  character(len=4), intent(in):: atom_type(*)

  integer i,j

  write(6,'(/"QMMM: Link Atom Information")')
  write(6,'("QMMM: ------------------------------------------------------------------------")')
  !ONLY PRINT XYZ COORDS IF THERE ARE LINK ATOMS
  if (qmmm_struct%nlink>0) then
    write(6,'("QMMM:  nlink = ",i5,"                   Link Coords              Resp Charges")') qmmm_struct%nlink
    write(6,'("QMMM:    MM(typ)  -  QM(typ)      X         Y         Z         MM        QM")')
    !Charges of MM link pairs has been zeroed in main array so we have to use the saved charges
    do j = 1, qmmm_struct%nlink
      write(6,'("QMMM:",i6," ",A4," ",i6," ",A4,5F10.3)') &
                qmmm_struct%link_pairs(1,j),atom_type(qmmm_struct%link_pairs(1,j)), &
                 qmmm_struct%iqmatoms(qmmm_struct%link_pairs(2,j)), &
                 atom_type(qmmm_struct%iqmatoms(qmmm_struct%link_pairs(2,j))), &
                ( qmcoords(i,j+qmmm_struct%nquant),i=1,3 ), &
                qmmm_struct%mm_link_pair_resp_charges(j)*INV_AMBER_ELECTROSTATIC, &
                qmmm_struct%qm_resp_charges(qmmm_struct%link_pairs(2,j))*INV_AMBER_ELECTROSTATIC
    end do
  else
    write(6,'("QMMM: nlink = ",i5)') qmmm_struct%nlink
  end if
  write(6,'("QMMM: ------------------------------------------------------------------------")')

  return

end subroutine print_link_atom_info

subroutine adj_mm_link_pair_crd(qmmm_nml,qmmm_struct, unimaged_coords)

  !This routine will save the coordinates of mm link pair atoms from the main amber 
  !coordinate array (unimaged) and then replace the coordinates of these atoms with
  !the corresponding link atom coordinate - this is used for things like qm_ewald ktable 
  !and GB.

  !If lnk_dis < 0.0d0 then the link atom is placed on top of the MM link pair atom which effectively
  !means there is no change.

  use qmmm_nml_module   , only : qmmm_nml_type
  use qmmm_struct_module, only : qmmm_struct_type
 implicit none

!passed in
  type(qmmm_nml_type), intent(inout) :: qmmm_nml
  type(qmmm_struct_type), intent(inout) :: qmmm_struct
  _REAL_, intent(inout) :: unimaged_coords(3,*) !natom long - amber's unimaged coordinate array (x)

!Local
  integer :: ier=0
  integer :: i, im, iq
  _REAL_ :: qmx, qmy, qmz, xij, yij, zij, rij2, onerij

   if (qmmm_nml%lnk_dis < 0.0d0) return !link atom sits on top of MMlink pair atom.

!Allocate mm_link_pair_saved_coords if first call
   if (qmmm_struct%adj_mm_link_pair_crd_first_call) then
     qmmm_struct%adj_mm_link_pair_crd_first_call = .false.
     qmmm_struct%mmcoords_contains_lnk_coords = .false.
     allocate(qmmm_struct%mm_link_pair_saved_coords(3,qmmm_struct%nlink), stat=ier)
     REQUIRE(ier == 0)
   end if

   if (qmmm_struct%mmcoords_contains_lnk_coords) then
     call sander_bomb('adj_mm_link_pair_crd','Call to set the mm link pair atom coordinates to be link atom coords', &
                      'when the main amber coordinate array already contains the link atom coordinates.')
   end if

   !link atom sees MM region then we need to save the coordinates of the MM link
   !pairs from amber's MAIN coordinate array (unimaged) and then work out what the link atom
   !coordinates should be in this reference frame and then insert them in.

   do i = 1, qmmm_struct%nlink
     !Save the MM atom's coordinates first
     im = qmmm_struct%link_pairs(1,i) !Position of mm atom in main amber coords array
     qmmm_struct%mm_link_pair_saved_coords(1:3,i) = unimaged_coords(1:3,im)

     iq = qmmm_struct%iqmatoms(qmmm_struct%link_pairs(2,i)) !Position of qm atom in main amber coords array
     qmx = unimaged_coords(1,iq)
     qmy = unimaged_coords(2,iq)
     qmz = unimaged_coords(3,iq)
     xij = unimaged_coords(1,im) - qmx
     yij = unimaged_coords(2,im) - qmy
     zij = unimaged_coords(3,im) - qmz
     rij2 = xij*xij + yij*yij + zij*zij
     onerij = 1.0d0/sqrt(rij2) * qmmm_nml%lnk_dis
     unimaged_coords(1,im) = xij * onerij + qmx
     unimaged_coords(2,im) = yij * onerij + qmy
     unimaged_coords(3,im) = zij * onerij + qmz
   end do

   qmmm_struct%mmcoords_contains_lnk_coords = .true.

   return

end subroutine adj_mm_link_pair_crd

subroutine rst_mm_link_pair_crd(qmmm_nml,qmmm_struct, unimaged_coords)

  !This routine replaces the mm coordinates from the mm_link_pair_saved array back
  !into the main amber array.

  use qmmm_nml_module   , only : qmmm_nml_type
  use qmmm_struct_module, only : qmmm_struct_type
  implicit none

!Passed in
  type(qmmm_nml_type), intent(inout) :: qmmm_nml
  type(qmmm_struct_type), intent(inout) :: qmmm_struct
  _REAL_, intent(out) :: unimaged_coords(3,*) !natom long - amber's unimaged coord array

!Local
  integer :: i, im

!Note if lnk_dis < 0.0d0 then the link atom sits on top of the MMlink atom so there is nothing to do.
  if (qmmm_nml%lnk_dis < 0.0d0) return !link atom sits on top of MMlink pair atom.

  if (.not. qmmm_struct%mmcoords_contains_lnk_coords) then
    call sander_bomb('rst_mm_link_pair_crd','Call to restore the mm link pair atom coordinates with the mm coords', &
                     'when the main amber coordinate array does not contain link atom coordinates.')
  end if

  do i = 1, qmmm_struct%nlink
    im = qmmm_struct%link_pairs(1,i) !Position of mm atom in main amber coords array
    unimaged_coords(1:3,im) = qmmm_struct%mm_link_pair_saved_coords(1:3,i)

  end do

  qmmm_struct%mmcoords_contains_lnk_coords = .false.

  return

end subroutine rst_mm_link_pair_crd

!------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+Sorts the array of link atoms pairs into numerical order
subroutine qmsort_link_pairs( link_pairs, nlink)

   implicit none
   integer, intent(in) :: nlink
   integer, intent(inout) :: link_pairs(2,nlink)

! Local
   integer i,j,lcurrent, qmcurrent
!
!   sort array in ascending order based on first index
!
   do i = 1, nlink
      lcurrent = link_pairs(1,i)
      qmcurrent = link_pairs(2,i)
      do j = i+1,nlink
         if (lcurrent.gt.link_pairs(1,j)) then
            link_pairs(1,i) = link_pairs(1,j)
            link_pairs(2,i) = link_pairs(2,j)
            link_pairs(1,j) = lcurrent
            link_pairs(2,j) = qmcurrent
            lcurrent = link_pairs(1,i)
            qmcurrent = link_pairs(2,i)
         endif
      end do
   end do
end subroutine qmsort_link_pairs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

