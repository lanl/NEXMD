! State variables required by random number generator
type :: rand_gen_state
   sequence
   !        real variables in Marsaglia algorithm
   double precision :: u(97), c, cd, cm
   !        pointers into u() in Marsaglia algorithm
   integer :: i97, j97
   !        set is true if amrset has been called; rand_dummy to make the
   !        size of the type a multiple of its largest element (important
   !        for Intel compilers, at least).
   logical :: set, rand_dummy
end type rand_gen_state
