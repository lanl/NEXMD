#include "copyright.i"

!*******************************************************************************
!
! Module:  random_mod
!
! Description: 
!
!*******************************************************************************

module random_mod

use file_io_dat_mod

  implicit none

! Real variables in Marsaglia algorithm:

  double precision, private, save       :: u(97)
  double precision, private, save       :: c
  double precision, private, save       :: cd
  double precision, private, save       :: cm

! Pointers into u() in Marsaglia algorithm:

  integer, private, save                :: i97
  integer, private, save                :: j97

! Set is .true. if amrset has been called:

  logical, private, save                :: set = .false.

contains

!*******************************************************************************
!
! Subroutine:  amrset
!
! Description: 
!
! Initialization routine for Marsaglias random number generator
! as implemented in Amber 3.0 Rev A by George Seibel.  See doc in amrand.
!
! Testing:  Call amrset with iseed = 54185253.  This should result
! in is1 = 1802 and is2 = 9373.  Call amrand 20000 times, then six
! more times, printing the six random numbers * 2**24 (ie, 4096*4096)
! They should be: (6f12.1)
! 6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
!              
! INPUT:        iseed
!
! OUTPUT:       to module private variables
!
!*******************************************************************************

subroutine amrset(iseed)

  implicit none

! Formal arguments:

  integer           iseed   ! integer seed greater than zero.

! Local variables:

! Two internal seeds used in Marsaglia algorithm:

  integer           is1
  integer           is2

! Max value of first seed (is1), 31328:

  integer           is1max

! Max value of second seed (is2), 30081

  integer           is2max

  integer           i, ii, j, jj, k, l, m

  double precision  s, t

  data is1max, is2max /31328, 30081/

! Construct two internal seeds from single unbound Amber seed: 
!
! is1 and is2 are quotient and remainder of iseed/IS2MAX.  We add
! one to keep zero and one results from both mapping to one.
! max and min functions keep is1 and is2 in required bounds.

  is1 = max((iseed / is2max) + 1, 1)
  is1 = min(is1, is1max)

  is2 = max(1, mod(iseed, is2max) + 1)
  is2 = min(is2, is2max)

  i = mod(is1/177, 177) + 2
  j = mod(is1    , 177) + 2
  k = mod(is2/169, 178) + 1
  l = mod(is2    , 169)

  do ii = 1, 97
    s = 0.0d0
    t = 0.5d0
    do jj = 1, 24
      m = mod(mod(i*j, 179)*k, 179)
      i = j
      j = k
      k = m
      l = mod(53*l + 1, 169)
      if (mod(l*m, 64) .ge. 32) s = s + t
      t = 0.5d0 * t
    end do
    u(ii) = s
  end do

  c  = 362436.0d0   / 16777216.0d0
  cd = 7654321.0d0  / 16777216.0d0
  cm = 16777213.0d0 / 16777216.0d0

  i97 = 97
  j97 = 33

  set = .true.

  return

end subroutine amrset

#ifdef  NEED_AMRAND
!*******************************************************************************
!
! Subroutine:  amrand
!
! Description: 
!
! Portable Random number generator by George Marsaglia
! Amber 3.0 Rev A implementation by George Seibel
!
! This random number generator originally appeared in *Toward a Universal
! Random Number Generator* by George Marsaglia and Arif Zaman.  Florida
! State University Report: FSU-SCRI-87-50 (1987)
!
! It was later modified by F. James and published in *A Review of Pseudo-
! random Number Generators*
!
! This is claimed to be the best known random number generator available.
! It passes ALL of the tests for random number generators and has a
! period of 2^144, is completely portable (gives bit identical results on
! all machines with at least 24-bit mantissas in the floating point
! representation).
!
! The algorithm is a combination of a Fibonacci sequence (with lags of 97
! and 33, and operation "subtraction plus one, modulo one") and an
! "arithmetic sequence" (using subtraction).
!
! INPUT:        from module private variables
!
! OUTPUT:       y:  A random number between 0.0 and 1.0
!
!*******************************************************************************

subroutine amrand(y)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision  y

! Local variables:

  double precision  uni

  if (.not. set) then
    write(mdout, '(a)') 'amrand not initialized'
    call mexit(6, 1)
  end if

  uni = u(i97) - u(j97)
  if (uni .lt. 0.0d0) uni = uni + 1.0d0
  u(i97) = uni
  i97 = i97 - 1
  if (i97 .eq. 0) i97 = 97
  j97 = j97 - 1
  if (j97 .eq. 0) j97 = 97
  c = c - cd
  if (c .lt. 0.0d0) c = c + cm
  uni = uni - c
  if (uni .lt. 0.0d0) uni = uni + 1.0d0
  y = uni

  return

end subroutine amrand
#endif /* NEED AMRAND */

!*******************************************************************************
!
! Subroutine:   gauss
!
! Description:  Generate a pseudo-random Gaussian sequence, with mean am and
!               std. dev. sd.  This is a version of amrand() that adds the
!               constraint of a gaussian distribution, with mean "AM" and
!               standard deviation "SD".  Output is to variable "V". It also
!               requires amrset to have been called first, and "uses up" the
!               same sequence that amrand() does.

!*******************************************************************************

subroutine gauss(am, sd, v)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)                  :: am
  double precision, intent(in)                  :: sd
  double precision, optional, intent(out)       :: v

! Local variables:

  double precision              :: tmp1, tmp2
  double precision              :: uni
  double precision              :: zeta1, zeta2

  if (.not. set) then
    write(mdout, '(a)') 'amrand not initialized!'
    call mexit(6, 1)
  end if

  ! Use the method of Box and Muller

  ! for some applications, one could use both "v" and "veven" in random
  ! sequence; but this won't work for most things we need (e.g. Langevin
  ! dynamics,) since the two adjacent variables are themselves highly
  ! correlated.  Hence we will just use the first ("v") variable.

  ! get two random numbers, even on (-1,1):

  do

    uni = u(i97) - u(j97)
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    u(i97) = uni
    i97 = i97 - 1
    if (i97 .eq. 0) i97 = 97
    j97 = j97 - 1
    if (j97 .eq. 0) j97 = 97
    c = c - cd
    if (c .lt. 0.0d0) c = c + cm
    uni = uni - c
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    zeta1 = uni + uni - 1.d0

    uni = u(i97) - u(j97)
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    u(i97) = uni
    i97 = i97 - 1
    if (i97 .eq. 0) i97 = 97
    j97 = j97 - 1
    if (j97 .eq. 0) j97 = 97
    c = c - cd
    if (c .lt. 0.0d0) c = c + cm
    uni = uni - c
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    zeta2 = uni + uni - 1.d0

    tmp1 = zeta1 * zeta1 + zeta2 * zeta2

    if (tmp1 .lt. 1.d0 .and. tmp1 .ne. 0.d0) then

      if (present(v)) then
        tmp2 = sd * sqrt(-2.d0 * log(tmp1)/tmp1)
        v = zeta1 * tmp2 + am
      end if

      return

    end if

  end do

end subroutine gauss

end module random_mod
