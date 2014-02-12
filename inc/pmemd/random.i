!*******************************************************************************
!
! Type: random_state
!
! Description:
!
! Variables necessary for Marsaglia random number stream. We create a 
! random_state type so that we can have independent random number streams that
! do not take from each others' sequence. This allows you to desynchronize
! random number sequences in one case (i.e. for Langevin dynamics to improve
! scaling), but keep another sequence synchronized if necessary (i.e. for REMD).
!
!*******************************************************************************

type :: random_state

  sequence

  ! Real variables in Marsaglia algorithm

  double precision, dimension(97) :: u
  double precision                :: c
  double precision                :: cd
  double precision                :: cm

  ! Pointers into u() in Marsaglia algorithm

  integer                         :: i97
  integer                         :: j97

  ! Determines if this random generator has been set or not

  logical                         :: set = .false.

end type random_state
