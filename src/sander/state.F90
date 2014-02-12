! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"


module state
  ! Mark James Williamson (MJW) Feb 2010
  !
  ! A universal state term for storing energy and other values in SANDER.
  ! Inspired partly from some of Bob Duke's examples in the PMEMD code.

  ! Motivations
  ! -----------
  ! The legacy nature of ene() and ener() arrays in the old code 
  ! limited the ability to add in new charmm related energy terms.
  ! These legacy arrays were initially used to account for the potential
  ! energy terms within the system in question, but over time seemed to have
  ! been hijacked for storing and transporting (via mpi) other variables.
  ! In addition to this, these arrays had fixed sizes with the element's
  ! position determining its meaning. Layered on top of this were arbitrary 
  ! regions in the middle of the array that were used as scratch space.
  ! Finally, some subroutines would be passed part of this array whilst
  ! other would access it via another equivalenced array.
  !
  ! All in all, a totally confusing beast.
  
  ! Overview
  ! --------
  ! This module provides a state derived type for containing various 
  ! properties of a system. The state_rec type in turn is composed of 
  ! other derived type such as potential_energy and kinetic_energy.
  ! The embedding of sub derived types enables certain operations
  ! to be carried out on parts of the system's state, such as a 
  ! mpi_reduce on just the potential energy terms of the system.
  !
  ! Below is a partial ASCII example of derived types interplay:
  !
  !   state_rec
  !      |
  !       \-kinetic_energy_rec
  !      |         |
  !      |          \total           == stateA%kin%tot
  !      |         |
  !      |          \Temp_solute     
  !      |      
  !       \-potential_energy_rec
  !      |         |
  !      |          \total
  !      |         |
  !      |          \bond            == stateA%pot%bond
  !       \total 
  !

  !
  ! Specific interfaces for the derived types have been coded to 
  ! allow array like operations on the various types such as
  ! stateA + stateB. These operations are generally used for 
  ! generating statistics. To save on handcoding, most of the
  ! interfaces make use of the transfer intrinsic to project
  ! the type in question onto a real array, carry out the maths
  ! easily using the fact that many intrinsics are available
  ! for real arrays and then finally projecting back onto the
  ! derived type after the operation has been carried out.

  !=============
  !| TODO List |
  !=============
  !
  ! - zero_neg_values_state() 
  ! - zero_pot_energy()       could both be improved
  !
  ! - Wide spread duplication of routines for calculating energy statistics,
  !   perhaps they should be consolidated here?
  !
  ! - Bookeeping using the various *_rec_len offsets for the various module
  !   interfaces is a bit hacky and needs peer review.
  !
  ! - Replace the need for the local ene() array in force()
  !
  ! - Write a function like fdist() in parallel.f, that is just for 
  !   energy terms.
  !
  ! - Consolidate the passing of ener%vir with ener in force()

  !===========================
  !| TODO before release List |
  !===========================
  ! 
  ! - Newbalance in fdist() is broken and USE_MPI_IN_PLACE may also be broken
  !
  ! - xref.f::cns_xref_run() is broken, I cannot understand what is being done
  !   in this routine to map it over to the new state_rec.


  implicit none
  private 

!--------------- Public variables and subroutines -----------
  ! Variables and derived types
  public :: potential_energy_rec     ! Potential energy type
  public :: potential_energy_rec_len ! Its length (used in mpi_* calls)
  public :: kinetic_energy_rec       ! Kinetic energy type
  public :: kinetic_energy_rec_len   ! Its length
  public :: sgld_rec       ! SGLD properties
  public :: sgld_rec_len   ! SGLD length
  public :: state_rec      ! State type; composed of the two types above 
                           ! and some other floats
  public :: state_rec_len  ! Its length

  public :: null_state_rec ! Used to zero all values in an
                           ! instance of a state type

  public :: null_potential_energy_rec
  public :: null_sgld_rec

  ! Subroutines / functions

  public :: zero_neg_values_state
  public :: zero_pot_energy

  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: sqrt
!------------- End Public variables and subroutines ---------


  integer, parameter :: kinetic_energy_rec_len = 7

  type kinetic_energy_rec
    sequence

                               ! What is it           What it was historically
                               !=============         ========================
                               
    _REAL_  :: tot             ! Total Kinetic          ener(2)
    _REAL_  :: solt            ! Temp_solute            ener(3)*onefac(1)
    _REAL_  :: solv            ! Temp_solvent           ener(4)*onefac(2)
    _REAL_  :: pres_scale_solt ! Press_SCALE_solute     ener(5)*onefac(3)
    _REAL_  :: pres_scale_solv ! Press_SCALE_solvent    ener(6)
    _REAL_  :: sgft            ! Guiding factor of
                               ! SGLD simulation        ener(48), ene(?)
    _REAL_  :: tempsg          ! Guiding temperature
                               ! of SGLD simulation     ener(49), ene(?)
  end type kinetic_energy_rec

  type(kinetic_energy_rec), parameter      :: null_kinetic_energy_rec = &
                                 kinetic_energy_rec(       &
                                 0.d0,0.d0,0.d0,0.d0,0.d0, &
                                 0.d0,0.d0)







  integer, parameter  :: potential_energy_rec_len = 25

  type potential_energy_rec
    sequence
                               ! What is it            What it was historically
                               ! ============          ========================

    _REAL_  :: tot             ! potential energy       ener(23), ene(1)
    _REAL_  :: vdw             ! van der Waals          ener(24), ene(2)
    _REAL_  :: elec            ! electrostatic energy   ener(25), ene(3)
    _REAL_  :: gb              ! GB (when igb.gt.0)     ener(26), ene(4)
    _REAL_  :: bond            ! Bond                   ener(27), ene(5)
    _REAL_  :: angle           ! Angle                  ener(28), ene(6)
    _REAL_  :: dihedral        ! Torsion                ener(29), ene(7)
    _REAL_  :: vdw_14          ! 1-4 non bonded         ener(30), ene(8)
    _REAL_  :: elec_14         ! 1-4 electrostatic      ener(31), ene(9)
    _REAL_  :: constraint      ! Constraint             ener(32), ene(10)

    _REAL_  :: polar           ! polarization           ener(33), ene(none)

    _REAL_  :: hbond           ! 10-12 (hb)
                               ! This was actually      ener(26), ene(4)
                               ! when igb == 0

    _REAL_  :: surf            ! surface-area dependent
                               ! energy, or cavity energy  
                               !                        ener(37), ene(23)

    _REAL_  :: scf             ! QMMM SCF Energy        ener(47), ene(25)
    _REAL_  :: disp            ! implicit solvation
                               ! dispersion energy      ener(50), ene(26)

    _REAL_  :: dvdl            ! charging free energy
                               ! result                 ener(39), ene(21)


! The next 3 are totally new
    _REAL_  :: angle_ub        ! CHARMM Urey-Bradley
    _REAL_  :: imp             ! CHARMM Improper
    _REAL_  :: cmap            ! CHARMM CMAP

    _REAL_  :: emap            ! EMAP constraint energy

    _REAL_  :: les             !  TODO
    _REAL_  :: noe             !  TODO
    _REAL_  :: pb              ! Poisson Boltzmann
    _REAL_  :: rism            ! 3D-RISM total solvation free energy
    _REAL_  :: ct              ! charge transfer

  end type potential_energy_rec

  type(potential_energy_rec), parameter      :: null_potential_energy_rec = &
                                 potential_energy_rec(&
                                 0.d0,0.d0,0.d0,0.d0,0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,&  ! 10
                                 0.d0,0.d0,0.d0,0.d0,0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,&  ! 20
                                 0.d0,0.d0,0.d0,0.d0,0.d0)

  integer, parameter :: sgld_rec_len = 14

  type sgld_rec
    sequence

                               ! What is it           What it was historically
                               !=============         ========================
                               
    _REAL_  :: sgft            ! momentum guiding factor
    _REAL_  :: sgff            ! force guiding factor
    _REAL_  :: tempsg          ! guiding temperature
    _REAL_  :: sgscal          ! energy conservation factor
    _REAL_  :: templf          ! low frequency temperature
    _REAL_  :: temphf          ! high frequency temperature
    _REAL_  :: treflf          ! reference low frequency temperature
    _REAL_  :: trefhf          ! reference high frequency temperature
    _REAL_  :: frclf           ! low frequency force factor
    _REAL_  :: frchf           ! high frequency force factor
    _REAL_  :: epotlf          ! low frequency potential energy
    _REAL_  :: epothf          ! high frequency potential energy
    _REAL_  :: virsg           ! guiding force virial
    _REAL_  :: sgwt            ! weighting power

  end type sgld_rec

  type(sgld_rec), parameter      :: null_sgld_rec = &
                                 sgld_rec(       &
                                 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                                 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

  ! The main type
  integer, parameter  :: state_rec_len = 1 + &
                                         kinetic_energy_rec_len   + &
                                         3 + 1 + 4 + 4 + 4        + &
                                         potential_energy_rec_len + &
                                         sgld_rec_len + &
                                         12
                         ! I wish there was a better way to calculate this :(


  type state_rec
    sequence
                               ! What is it            What it was historically
                               ! ============          ========================


    _REAL_  :: tot             ! Total energy           ener(1)

    type(kinetic_energy_rec) &
            :: kin             ! Kinetic energy         ener(2:4) ?

    _REAL_  :: box(3)          ! Box_{x,y,z}            ener(7:9)
    _REAL_  :: volume          ! Volume of box          ener(10)
    _REAL_  :: pres(1:4)       ! pres_{x,y,z,total}     ener(11:14)
    _REAL_  :: vir(1:4)        ! vir_{x,y,z,total}      ener(19:22)
    _REAL_  :: cmt(1:4)        ! cmt_{x,y,z,total}      ener(15:18)

    type(potential_energy_rec) &
            :: pot             ! Total potential        ener(23:), ene(1:30)

    !More bits?
    _REAL_  :: aveper          ! ave perm moment        ener(34), ene(none)
    _REAL_  :: aveind          ! ave ind  moment        ener(35), ene(none)
    _REAL_  :: avetot          ! ave total moment       ener(36), ene(none)
    _REAL_  :: exray           ! xray related           ener(38), ene(?)
    _REAL_  :: eptot           ! potential enery
                               ! for a subset of atoms
                               ! also used by xray as "r"  ener(40), ene(24) 
                               ! I think this is REMD specific?
    _REAL_  :: rfree           ! xray related           ener(41), ene(?)
    _REAL_  :: density         !                        ener(42), ene(?)
    _REAL_  :: virvsene        ! Ewald error estimate?  ener(43), ene(?)
    _REAL_  :: diprms          !                        ener(44)
    _REAL_  :: dipiter         !                        ener(45)
    _REAL_  :: dipole_temp     !                        ener(46)
    _REAL_  :: surface_ten     ! Surface tension when running with csurften>0

    type(sgld_rec) &
            :: sgld             ! SGLD properties

end type state_rec

  type(state_rec), parameter      :: null_state_rec = state_rec(&
                                 0.d0,null_kinetic_energy_rec,0.d0,0.d0,0.d0, &  ! 5
                                 0.d0,0.d0,null_potential_energy_rec,0.d0,0.d0,& ! 10
                                 0.d0,0.d0,0.d0,0.d0,0.d0, &                     ! 15
                                 0.d0,0.d0,0.d0,0.d0,0.d0,null_sgld_rec)         ! 21

  interface assignment(=)
    module procedure assign_state_with_float
  end interface
  interface operator(+)
    module procedure add_state
  end interface
  interface operator(*)
    module procedure mul_state
    module procedure mul_state_by_float
  end interface
  interface operator(-)
    module procedure sub_state
  end interface
  interface operator(/)
    module procedure div_state
    module procedure div_state_by_float
    module procedure div_state_by_int
  end interface
  interface sqrt
    module procedure sqrt_state
  end interface





  contains

  ! This is how it should be done offically....
  !
  !function add_state(a, b)
  !  type(state_rec), intent(in)    :: a,b
  !  type(state_rec)                :: add_state

  !  add_state%tot      = a%tot  + b%tot
  !  add_state%kin      = a%kin  + b%kin
  !  add_state%box      = a%box  + b%box
  !  ...etc.....

  !end function add_state

  ! But this is how to do it quicker (and legally) by
  !  making use of the transfer intrinsic:

  function add_state(stateA, stateB)
    type(state_rec)                :: add_state
    type(state_rec), intent(in)    :: stateA, stateB

    _REAL_              :: arrayA(1:state_rec_len)
    _REAL_              :: arrayB(1:state_rec_len)

  ! Use the TRANSFER intrinsic over EQUIVALENCE.
  ! http://www.liv.ac.uk/HPC/HTMLF90Course/HTMLF90CourseNotesnode262.html
  ! lengthData = size(transfer(r, data))
  !  Seems to work :

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA) 
    arrayB = transfer(stateB, arrayB) 

    ! Now do some array math
    arrayA = arrayA + arrayB

    ! Now project back for the output
    add_state = transfer(arrayA, stateA)
  end function add_state


!-------------------------------------------------
  function mul_state(stateA, stateB)
    type(state_rec)                :: mul_state
    type(state_rec), intent(in)    :: stateA, stateB

    _REAL_              :: arrayA(1:state_rec_len)
    _REAL_              :: arrayB(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)
    arrayB = transfer(stateB, arrayB)

    ! Now do some array math
    arrayA = arrayA * arrayB

    ! Now project back for the output
    mul_state = transfer(arrayA, stateA)
  end function mul_state

!-------------------------------------------------
  function mul_state_by_float(stateA, myFloat)
    type(state_rec)                :: mul_state_by_float
    type(state_rec), intent(in)    :: stateA
    _REAL_, intent(in)             :: myFloat

    _REAL_              :: arrayA(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    ! Now do some array math
    arrayA = arrayA * myFloat

    ! Now project back for the output
    mul_state_by_float = transfer(arrayA, stateA)
  end function mul_state_by_float


!-------------------------------------------------
  function sub_state(stateA, stateB)
    type(state_rec)                :: sub_state
    type(state_rec), intent(in)    :: stateA, stateB

    _REAL_              :: arrayA(1:state_rec_len)
    _REAL_              :: arrayB(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)
    arrayB = transfer(stateB, arrayB)

    ! Now do some array math
    arrayA = arrayA - arrayB

    ! Now project back for the output
    sub_state = transfer(arrayA, stateA)
  end function sub_state

!-------------------------------------------------
  function div_state(stateA, stateB)
    type(state_rec)                :: div_state
    type(state_rec), intent(in)    :: stateA, stateB

    _REAL_              :: arrayA(1:state_rec_len)
    _REAL_              :: arrayB(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)
    arrayB = transfer(stateB, arrayB)

    ! Now do some array math
    ! arrayA = arrayA / arrayB

    ! RCW: Use vector functions
    call vdinv(state_rec_len, arrayB, arrayB)
    arrayA = arrayA * arrayB

    ! Now project back for the output
    div_state = transfer(arrayA, stateA)
  end function div_state

!-------------------------------------------------
  function div_state_by_float(stateA, myFloat)
    type(state_rec)                :: div_state_by_float
    type(state_rec), intent(in)    :: stateA
    _REAL_, intent(in)             :: myFloat

    _REAL_              :: arrayA(1:state_rec_len)
    _REAL_              :: one_float

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    one_float = 1.0d0 / myFloat

    ! Now do some array math
    arrayA = arrayA * one_float

    ! Now project back for the output
    div_state_by_float = transfer(arrayA, stateA)
  end function div_state_by_float

!-------------------------------------------------
  function div_state_by_int(stateA, myInt)
    type(state_rec)                :: div_state_by_int
    type(state_rec), intent(in)    :: stateA
    integer, intent(in)            :: myInt

    _REAL_              :: arrayA(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    ! Now do some array math
    arrayA = arrayA / myInt

    ! Now project back for the output
    div_state_by_int = transfer(arrayA, stateA)
  end function div_state_by_int

!----------------------------------------------
  subroutine assign_state_with_float(stateA, myFloat)
    type(state_rec), intent(inout) :: stateA
    _REAL_, intent(in)             :: myFloat
    _REAL_                         :: arrayA(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    ! Now do some array math
    arrayA = myFloat

    ! Now project back for the output
    stateA = transfer(arrayA, stateA)
  end subroutine assign_state_with_float

!--------------------------------------------------------
  subroutine zero_pot_energy(Mypot_energy)
    type(potential_energy_rec), intent(inout) :: Mypot_energy
    _REAL_                         :: arrayA(1:potential_energy_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(Mypot_energy, arrayA)

    ! Now do some array math
    arrayA = 0.d0

    ! Now project back for the output
    Mypot_energy = transfer(arrayA, Mypot_energy)
  end subroutine zero_pot_energy


!--------------------------------------------------------
  subroutine zero_neg_values_state(stateA)
    type(state_rec), intent(inout) :: stateA
    _REAL_                         :: arrayA(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    ! Now do some array math
    where (arrayA < 0.d0)
      arrayA = 0.d0
    endwhere

    ! Now project back for the output
    stateA = transfer(arrayA, stateA)
  end subroutine zero_neg_values_state


!-----------------------------------------------
  function sqrt_state(stateA)
    type(state_rec)                :: sqrt_state
    type(state_rec), intent(in)    :: stateA

    _REAL_              :: arrayA(1:state_rec_len)

    ! Project the state_recs onto arrays
    arrayA = transfer(stateA, arrayA)

    ! Now do some array math
    ! arrayA = sqrt(arrayA)

    ! RCW use vector functions
    call vdsqrt(state_rec_len, arrayA, arrayA)

    ! Now project back for the output
    sqrt_state = transfer(arrayA, stateA)
  end function sqrt_state


end module state
