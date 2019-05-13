! <compile=optimized>

#include "dprec.fh"
subroutine slkmatrices_a_to_bohrs(i,j,qm_coords,ham,over,matsize)

   ! Passed in:
   !   matsize: Dimension of the (square) arrays 'ham' and 'over'
   !   i,j: Atoms between which the interaction is being calculated
   !   qm_coords: the coordinates matrix

   ! Output:
   !   ham:  The Hamiltonian matrix
   !   over: The Overlap matrix

   use constants, only: A_TO_BOHRS
   use qm2_dftb_module, only: lmax, izp_str

   implicit none

   ! Passed in:
   integer, intent(in)  :: i,j           ! Atoms between which the interaction is being calculated
   integer, intent(in)  :: matsize           ! Dimension of the (square) arrays 'ham' and 'over'
   _REAL_ , intent(in)  :: qm_coords(3,*)      ! Coordinates matrix in angstroms
   _REAL_ , intent(out) :: ham(matsize,matsize)  ! The Hamiltonian matrix
   _REAL_ , intent(out) :: over(matsize,matsize) ! The Overlap matrix


   ! Locals:
   integer :: izpi,izpj                  ! IZP pf each atom
   _REAL_  :: dif(3)                     ! Distance vector between the 2 atoms
   external skspar,skhpar                ! Subroutines to calculate the terms of:
                                         !    -> Overlap matrix     (skspar)
                                         !    -> Hamiltonian matrix (skhpar)


   ! Distance vector between the 2 atoms
   dif(1:3)=(qm_coords(1:3,j)-qm_coords(1:3,i))*A_TO_BOHRS

   izpi= izp_str%izp(i)
   izpj= izp_str%izp(j)

   ! Zero all elements of the matrices

   ham  = 0.0d0
   over = 0.0d0

   ! Get hamiltonian (use skhpar)
   call slkode(dif,izpi,izpj,ham,matsize,skhpar)

   ! Get overlap (use skspar)
   call slkode(dif,izpi,izpj,over,matsize,skspar)

end subroutine slkmatrices_a_to_bohrs

subroutine slkmatrices(i,j,qm_coords,ham,over,matsize)

   ! Passed in:
   !   matsize: Dimension of the (square) arrays 'ham' and 'over'
   !   i,j: Atoms between which the interaction is being calculated
   !   qm_coords: the coordinates matrix

   ! Output:
   !   ham:  The Hamiltonian matrix
   !   over: The Overlap matrix

   use qm2_dftb_module, only: lmax, izp_str
   use constants, only: A_TO_BOHRS

   implicit none

   ! Passed in:
   integer, intent(in)  :: i,j           ! Atoms between which the interaction is being calculated
   integer, intent(in)  :: matsize           ! Dimension of the (square) arrays 'ham' and 'over'
   _REAL_ , intent(in)  :: qm_coords(3,*)      ! Coordinates matrix in angstroms
   _REAL_ , intent(out) :: ham(matsize,matsize)  ! The Hamiltonian matrix
   _REAL_ , intent(out) :: over(matsize,matsize) ! The Overlap matrix


   ! Locals:
   integer :: izpi,izpj                  ! IZP pf each atom
   _REAL_  :: dif(3)                     ! Distance vector between the 2 atoms
   external skspar,skhpar                ! Subroutines to calculate the terms of:
                                         !    -> Overlap matrix     (skspar)
                                         !    -> Hamiltonian matrix (skhpar)


   ! Distance vector between the 2 atoms
   dif(1:3)=(qm_coords(1:3,j)-qm_coords(1:3,i)) * A_TO_BOHRS

   izpi= izp_str%izp(i)
   izpj= izp_str%izp(j)

   ! Zero all elements of the matrices

   ham  = 0.0d0
   over = 0.0d0

   ! Get hamiltonian (use skhpar)
   call slkode(dif,izpi,izpj,ham,matsize,skhpar)

   ! Get overlap (use skspar)
   call slkode(dif,izpi,izpj,over,matsize,skspar)

end subroutine slkmatrices

subroutine slkode(dist,i,j,em,LDIM,iovpar)

  use qm2_dftb_module, only: lmax, dummy

  implicit none

  !Passed in:
  _REAL_ , intent(in)  :: dist(3)       ! Vector distance between the 2 atoms
  integer, intent(in)  :: i,j
  integer, intent(in)  :: LDIM          ! Dimension of the EM and dummy matrices
  _REAL_ , intent(out) :: em(LDIM,LDIM) ! Container for the H or S matrix (depending on iovpar)
  external iovpar                       ! Subroutine to use:
                                        !     skspar for overlaps
                                        !     skhpar for hamiltonian

  !Locals
  _REAL_  :: x(6),x2(6),r2,r2i,ri
  integer :: l,k, maxmax, minmax



  r2=0.0
  do l=1,3
     x(l)=dist(l)
     x2(l)=x(l)*x(l)
     r2=r2+x2(l)
  end do

  if(r2 >= 1.0e-8)then
     r2i = 1.0 / r2
     ri = sqrt(r2i)
     do l=1,3
        x(l)   = x(l) * ri
        x(l+3) = x(l)
        x2(l)  = x2(l) * r2i
        x2(l+3)= x2(l)
     end do

     ! Number of orbitals
     maxmax = max(lmax(i),lmax(j))
     minmax = min(lmax(i),lmax(j))

     ! Get s-s interactions
     call skss(x,x2,i,j,r2,iovpar,em,ldim)
     if (maxmax <= 1) return

     ! Get s-p and p-p interactions
     if(minmax >= 2)then
        call skpp(x,x2,i,j,r2,iovpar,em(2,2),ldim)
        call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),ldim)
        if(i /= j)then
           call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),ldim)
        endif
     else
        if(lmax(j) >= 2)then
           call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),ldim)
        else
           call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),ldim)
        endif
     endif
     if(maxmax <= 2) return

     ! Get s-d, p-d and d-d
     if(minmax == 3)then
        call skdd(x,x2,i,j,r2,iovpar,em(5,5),ldim)
        call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
        call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),ldim)
        if(i /= j)then
           call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
           call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),ldim)
        endif
     else
        if(lmax(i) == 1)then
           call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
        else
           if(lmax(i) == 2)then
              call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
              call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),ldim)
           else
              if(lmax(j) == 1)then
                 call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
              else
                 call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
                 call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),ldim)
              endif !  if(lmax(j) == 1)then
           endif ! if(lmax(i) == 2)then
        endif ! if(lmax(i) == 1)then
     endif ! if(minmax == 3)then

  else ! if(r2 >= 1.0e-8)then

     ! The interactions are on the same atom.

     em = 0.0d0

     if(i /= j) return

     call selfs(i,j,r2,iovpar,em,ldim)
     if(lmax(i) <= 1) return

     call selfp(i,j,r2,iovpar,em(2,2),ldim)
     if(lmax(i) <= 2) return

     call selfd(i,j,r2,iovpar,em(5,5),ldim)

  endif ! if(r2 >= 1.0e-8)then

end subroutine slkode
