#include "copyright.h"
#include "dprec.h"

!-----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize common block for default random stream
block data raset

#include "random.h"
   type (rand_gen_state) :: def_gen
   !     Initializes random number initialization flag to false.
   common /raset1/ def_gen
   
   data def_gen%set /.false./
end   
!-----------------------------------------------------------------------



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper to initialize default random generator
subroutine amrset(iseed)
#include "random.h"
   integer, intent (in) :: iseed  
   type (rand_gen_state) :: def_gen
   common /raset1/ def_gen
   call amrset_gen(def_gen, iseed)
end subroutine amrset
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initializes the random number process
subroutine amrset_gen(st, iseed)
   implicit none
   !     Initialization routine for Marsaglias random number generator
   !     as implemented in Amber 3.0 Rev A by George Seibel
   !     See doc in amrand.
   !     Testing:  Call amrset with iseed = 54185253.  This should result
   !     in is1 = 1802 and is2 = 9373.  Call amrand 20000 times, then six
   !     more times, printing the six random numbers * 2**24 (ie, 4096*4096)
   !     They should be: (6f12.1)
   !     6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
   
   !     INPUT:
#include "random.h"
   type (rand_gen_state), intent (out) :: st

   integer, intent (in) :: iseed
   !        ... integer seed greater than zero
   
   !     INTERNAL:
   
   integer is1, is2
   !        ... the two internal seeds used in Marsaglia algorithm
   integer, parameter :: is1max = 31328, is2max = 30081
   !        ... max values of internal seeds
   integer i,j,k,l,m
   !        ... used in generation of st%u()
   double precision s,t
   !        ... used in generation of st%u()
   integer ii, jj
   !        ... loop indices
   

   !     --- construct two internal seeds from single unbound Amber seed ---
   
   !         is1 and is2 are quotient and remainder of iseed/IS2MAX.  We add
   !         one to keep zero and one results from both mapping to one.
   !         max and min functions keep is1 and is2 in required bounds.
   
   is1 = max((iseed / is2max)+1, 1)
   is1 = min(is1, is1max)
   
   is2 = max(1, mod(iseed, is2max)+1)
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
         l = mod(53*l+1, 169)
         if (mod(l*m, 64) >= 32) s = s + t
         t = 0.5d0 * t
      end do
      st%u(ii) = s
   end do
   
   st%c  = 362436.0d0   / 16777216.0d0
   st%cd = 7654321.0d0  / 16777216.0d0
   st%cm = 16777213.0d0 / 16777216.0d0
   
   st%i97 = 97
   st%j97 = 33
   
   st%set = .true.
   return
end subroutine amrset_gen
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper to retreive uniform variate from default generator
subroutine amrand(y)
#include "random.h"
   _REAL_, intent (out) :: y
   type (rand_gen_state) :: def_gen
   common /raset1/ def_gen
   call amrand_gen(def_gen, y)
end subroutine amrand
!-----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Generate uniform variate y based on state variables st
subroutine amrand_gen(st, y)
   implicit none
   !     Portable Random number generator by George Marsaglia
   !     Amber 3.0 Rev A implementation by George Seibel
   
#include "random.h"
   type (rand_gen_state), intent (inout) :: st
   
   !     OUTPUT:
   
   !     y:  A random number between 0.0 and 1.0
   
   _REAL_, intent (out) :: y
   
   !     INTERNAL:
   
   double precision uni
   !        ... working var. for random number
   
   !     This random number generator originally appeared in *Toward a Universal
   !     Random Number Generator* by George Marsaglia and Arif Zaman.  Florida
   !     State University Report: FSU-SCRI-87-50 (1987)
   
   !     It was later modified by F. James and published in *A Review of Pseudo-
   !     random Number Generators*
   
   !     This is claimed to be the best known random number generator available.
   !     It passes ALL of the tests for random number generators and has a
   !     period of 2^144, is completely portable (gives bit identical results on
   !     all machines with at least 24-bit mantissas in the floating point
   !     representation).
   
   !     The algorithm is a combination of a Fibonacci sequence (with lags of 97
   !     and 33, and operation "subtraction plus one, modulo one") and an
   !     "arithmetic sequence" (using subtraction).
   
   if ( .not. st%set ) then
      write(6,'(a)') 'amrand not initd'
      call mexit(6, 1)
   end if
   
   uni = st%u(st%i97) - st%u(st%j97)
   if ( uni < 0.0d0 ) uni = uni + 1.0d0
   st%u(st%i97) = uni
   st%i97 = st%i97 - 1
   if (st%i97 == 0) st%i97 = 97
   st%j97 = st%j97 - 1
   if (st%j97 == 0) st%j97 = 97
   st%c = st%c - st%cd
   if ( st%c < 0.0d0 ) st%c = st%c + st%cm
   uni = uni - st%c
   if ( uni < 0.0d0 ) uni = uni + 1.0d0
#ifdef DPREC
   y = uni
#else
   y = real(uni)
#endif
   
   return
end subroutine amrand_gen 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Default generator wrapper for gauss_gen
subroutine gauss(am,sd,v)
   implicit none
#include "random.h"
   _REAL_, intent (in) :: am, sd
   _REAL_, intent (out) :: v
   type (rand_gen_state) :: def_gen
   common /raset1/ def_gen  
   call gauss_gen(def_gen, am, sd, v)
end subroutine gauss
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Generates a pseudo-random Gaussian sequence, with mean am and std. dev. sd
subroutine gauss_gen(st,am,sd,v)
   implicit none
   
   !     This is a version of amrand() that adds the constraint of
   !     a gaussian distribution, with mean "AM" and standard deviation "SD".
   !     Output is to variable "V". It also requires amrset to have
   !     been called first, and "uses up" the same sequence that
   !     amrand() does.

#include "random.h"
   type (rand_gen_state), intent (inout) :: st
   
   _REAL_, parameter :: zero = 0.0d0, six = 6.0d0, twopi = 6.28318531d0
   _REAL_, intent (in) :: am, sd
   _REAL_, intent (out) :: v
   double precision :: a, uni, zeta(2), tmp1, tmp2 
   double precision, save :: veven
   integer :: i
   logical :: odd = .true.
   
   if ( .not. st%set ) then
      write(6,'(a)') 'amrand not initd'
      call mexit(6, 1)
   end if
   !  Use the method of Box and Muller

   !  for some applications, one could use both "v" and "veven" in random
   !  sequence; but this won't work for most things we need (e.g. Langevin
   !  dynamics,) since the two adjacent variables are themselves highly
   !  correlated.  Hence we will just use the first ("v") variable.  The
   !  code commented out here is maintained in case someone wants to
   !  re-examine this issue later.

!  if( .not. odd ) then
!#ifdef DPREC
!    v = veven
!#else
!    v = real(veven)
!#endif
!    odd = .true.
!    return
!  end if

   ! get two random numbers, even on (-1,1):

10 do i=1,2
      uni = st%u(st%i97) - st%u(st%j97)
      if ( uni < 0.0d0 ) uni = uni + 1.0d0
      st%u(st%i97) = uni
      st%i97 = st%i97 - 1
      if (st%i97 == 0) st%i97 = 97
      st%j97 = st%j97 - 1
      if (st%j97 == 0) st%j97 = 97
      st%c = st%c - st%cd
      if ( st%c < 0.0d0 ) st%c = st%c + st%cm
      uni = uni - st%c
      if ( uni < 0.0d0 ) uni = uni + 1.0d0
      zeta(i) = uni + uni - 1.d0
    end do

   tmp1 = zeta(1)*zeta(1) + zeta(2)*zeta(2)
   if( tmp1 >= 1.d0 .or. tmp1 == 0.d0 ) goto 10
   tmp2 = sd*sqrt( -2.d0*log(tmp1)/tmp1 )
#ifdef DPREC
   v =        zeta(1)*tmp2 + am
#else
   v =        real(zeta(1)*tmp2 + am)
#endif
!  veven =    zeta(2)*tmp2 + am
!  odd = .false.

   return
end subroutine gauss_gen 
!-----------------------------------------------------------------------

