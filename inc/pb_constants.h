!-------------------------------------------------------
! Physical Constants

! The conversion factor from electron/A to kcal/mol/A.
_REAL_, parameter :: AMBER_ELECTROSTATIC = 18.2223d0
_REAL_, parameter :: AMBER_ELECTROSTATIC2 =  AMBER_ELECTROSTATIC*AMBER_ELECTROSTATIC
_REAL_, parameter :: INV_AMBER_ELECTROSTATIC = 1.0d0/AMBER_ELECTROSTATIC
_REAL_, parameter :: INV_AMBER_ELECTROSTATIC2 = INV_AMBER_ELECTROSTATIC*INV_AMBER_ELECTROSTATIC

! The Bohr radius of 1998, physics.nist.gov/constants.
_REAL_, parameter :: BOHR_RADIUS = 0.5291772083d0
_REAL_, parameter :: INV_BOHR_RADIUS = 1.0d0/BOHR_RADIUS


!-------------------------------------------------------
! Numeric Constants

_REAL_, parameter ::      PI = 3.14159265358979323846d0 
_REAL_, parameter ::   TWOPI = 2.0d0*PI
_REAL_, parameter ::  FOURPI = 4.0d0*PI
_REAL_, parameter :: EIGHTPI = 8.0d0*PI
_REAL_, parameter :: INV_TWOPI = 1.0d0/TWOPI
_REAL_, parameter :: INV_FOURPI = 1.0d0/FOURPI
_REAL_, parameter :: INV_EIGHTPI = 1.0d0/EIGHTPI

_REAL_, parameter :: SQRT2 = 1.41421356237309504880d0
_REAL_, parameter :: SQRT3 = 1.73205080756887729352d0

!-------------------------------------------------------
! Generic Floating Point Constants

_REAL_, parameter :: TEN_TO_MINUS1  = 1.0d-1
_REAL_, parameter :: TEN_TO_MINUS2  = 1.0d-2
_REAL_, parameter :: TEN_TO_MINUS3  = 1.0d-3
_REAL_, parameter :: TEN_TO_MINUS4  = 1.0d-4
_REAL_, parameter :: TEN_TO_MINUS5  = 1.0d-5
_REAL_, parameter :: TEN_TO_MINUS6  = 1.0d-6
_REAL_, parameter :: TEN_TO_MINUS10 = 1.0d-10
_REAL_, parameter :: TEN_TO_PLUS1   = 1.0d+1
_REAL_, parameter :: TEN_TO_PLUS2   = 1.0d+2
_REAL_, parameter :: TEN_TO_PLUS3   = 1.0d+3
_REAL_, parameter :: TEN_TO_PLUS10  = 1.0d+10

_REAL_, parameter :: ZERO     = 0.0d0
_REAL_, parameter :: ONE      = 1.0d0
_REAL_, parameter :: TWO      = 2.0d0
_REAL_, parameter :: THREE    = 3.0d0
_REAL_, parameter :: FOUR     = 4.0d0
_REAL_, parameter :: FIVE     = 5.0d0
_REAL_, parameter :: SIX      = 6.0d0
_REAL_, parameter :: SEVEN    = 7.0d0
_REAL_, parameter :: EIGHT    = 8.0d0
_REAL_, parameter :: NINE     = 9.0d0
_REAL_, parameter :: TEN      = 10.0d0
_REAL_, parameter :: ELEVEN   = 11.0d0
_REAL_, parameter :: TWELVE   = 12.0d0

_REAL_, parameter :: HALF     = ONE/TWO
_REAL_, parameter :: THIRD    = ONE/THREE
_REAL_, parameter :: FOURTH   = ONE/FOUR
_REAL_, parameter :: FIFTH    = ONE/FIVE
_REAL_, parameter :: SIXTH    = ONE/SIX
_REAL_, parameter :: SEVENTH  = ONE/SEVEN
_REAL_, parameter :: EIGHTH   = ONE/EIGHT
_REAL_, parameter :: NINTH    = ONE/NINE
_REAL_, parameter :: TENTH    = ONE/TEN
_REAL_, parameter :: ELEVENTH = ONE/ELEVEN
_REAL_, parameter :: TWELFTH  = ONE/TWELVE
