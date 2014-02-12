#include "dprec.fh"
module qm2_pm6_hof_module
! ----------------------------------------------------------------------
! PURPOSE: PM6 specific corrections to the Heat of Formation (HOF)
!          and corresponding gradient
! 
! The correction is calculated in kcal/mol
! It is not documented in any publication. I took the information from
! Jimmy Stewart's code and double checked results by comparison against
! HOFs obtained with MOPAC2009 and Gaussian09
!
! For the moment the covalent radii are stored in this subroutine
! They should be moved into a qm2_parameters module when the code
! will be cleaned up
! 
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : April 2010
! ----------------------------------------------------------------------
  
  implicit none

  private
#if 0
  public :: pm6_correction
#endif
  public :: corInfoType, print, hofCorrection, hofCorrectionGradient
  public :: cct, nsp2
  public :: strlen

  interface print
     module procedure printCorInfoType
  end interface

#if 0
  interface pm6_correction
     module procedure hofCorrection
     module procedure hofCorrectionGradient
  end interface
#endif

  ! Data type collecting information on PM6 HOF corrections
  type corInfoType
     logical :: inUse  ! correction in use ?
     integer :: natom  ! correction for how many atoms ?
     _REAL_  :: energy ! correction to HOF
  end type corInfoType

  type(corInfoType) :: cct, nsp2
  integer, parameter :: strlen = 80

contains

  ! -------------------------------
  ! Calculate PM6 correction to HOF
  ! -------------------------------
  _REAL_ function hofCorrection()

    use constants, only : zero
    use qmmm_module, only : qmmm_struct
    implicit none

    integer :: numBonds(qmmm_struct%nquant_nlink)
    integer :: bondedAtoms(qmmm_struct%nquant_nlink,qmmm_struct%nquant_nlink)
    integer :: natom
    logical, parameter :: debug = .false.
    integer :: iat, jat
    
    hofCorrection = zero
    natom = qmmm_struct%nquant_nlink

    ! Determine which atoms are bonded to each other
    call setupDentate(natom, qmmm_struct%qm_coords, numBonds, bondedAtoms)

    if (debug ) then
       ! print info on bonded atoms
       write (6,'(a)') ' ** BONDS **' 
       do iat = 1, natom
          write (6, '(/a,i3)') ' atom number :', iat
          write (6, '(a,i3)') ' bonds       :', numBonds(iat)
          write (6, '(a,10i3)') ' partners    :', (bondedAtoms(jat,iat), jat=1,numBonds(iat))
       end do
    end if

    ! HOF CC triple bond correction
    call ccTripleBond(natom, qmmm_struct%qm_coords, numBonds, bondedAtoms, cct)
    hofCorrection = hofCorrection + cct%energy

    ! HOF MM correction for nitrogen atoms with three ligands
    call nsp2Correction(natom, qmmm_struct%qm_coords, numBonds, bondedAtoms, nsp2)
    hofCorrection = hofCorrection + nsp2%energy

  end function hofCorrection

  ! ----------------------------------------------
  ! Determine which atoms are bonded to each other
  ! ----------------------------------------------
  subroutine setupDentate(natom, coord, numBonds, bondedAtoms)

    use qmmm_module, only : qmmm_struct
    implicit none
    integer, intent(in)  :: natom
    _REAL_,  intent(in)  :: coord(3,natom)
    integer, intent(out) :: numBonds(natom)
    integer, intent(out) :: bondedAtoms(natom,natom)

    integer :: iat, jat, iatnum, jatnum
    _REAL_ :: distance, threshold
    _REAL_ :: vec(3)
    !  Atoms are assumed attached if they are within
    !  1.1  times the sum of their covalent radii.
    _REAL_, parameter :: factor = 1.1d0
    _REAL_ :: covalentAtomRadius(107)
    data covalentAtomRadius / &
    !      H       He
         & 0.37d0, 0.32d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 1.34d0, 0.90d0, 0.82d0, 0.77d0, 0.75d0, 0.73d0, 0.71d0, 0.69d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 1.54d0, 1.30d0, 1.18d0, 1.11d0, 1.06d0, 1.02d0, 0.99d0, 0.97d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 1.96d0, 1.74d0, 1.44d0, 1.36d0, 1.25d0, 1.27d0, 1.39d0, 1.25d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 1.26d0, 1.21d0, 1.38d0, 1.31d0, 1.26d0, 1.22d0, 1.19d0, 1.16d0, &
    !      Br      Kr
         & 1.14d0, 1.10d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 2.11d0, 1.92d0, 1.62d0, 1.48d0, 1.37d0, 1.45d0, 1.56d0, 1.26d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 1.35d0, 1.31d0, 1.53d0, 1.48d0, 1.44d0, 1.41d0, 1.38d0, 1.35d0, &
    !      I       Xe
         & 1.33d0, 1.30d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 2.25d0, 1.98d0, 1.69d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 1.60d0, 1.50d0, 1.38d0, 1.46d0, 1.59d0, 1.28d0, 1.37d0, 1.28d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 1.44d0, 1.49d0, 1.48d0, 1.47d0, 1.46d0, 999.d0, 999.d0, 1.45d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lr      Rf      Db      Sg       Tv
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0 /

    numBonds = 0
    bondedAtoms = 0

    do iat = 1, natom
       do jat = 1, iat - 1
          vec(1:3) = coord(1:3,iat) - coord(1:3,jat)
          distance = vec(1)**2 + vec(2)**2 + vec(3)**2
          iatnum = qmmm_struct%iqm_atomic_numbers(iat)
          jatnum = qmmm_struct%iqm_atomic_numbers(jat)
          threshold = ( factor * (covalentAtomRadius(iatnum)+covalentAtomRadius(jatnum)) )**2
          if (distance < threshold) then
             numBonds(iat) = numBonds(iat) + 1
             numBonds(jat) = numBonds(jat) + 1
             bondedAtoms(numBonds(iat), iat) = jat
             bondedAtoms(numBonds(jat), jat) = iat
          end if
       end do
    end do
    
  end subroutine setupDentate

  ! ----------------------------------------------------------------
  ! Evaluate energy contribution from acetylenic bonds - this is
  ! a correction to account for the extra stabilization of yne bonds
  ! (written according to Jimmy Stewart's code)
  ! ----------------------------------------------------------------
  subroutine ccTripleBond(natom, coord, numBonds, bondedAtoms, cct)

    use qmmm_module, only : qmmm_struct
    use constants, only : zero
    implicit none

    integer, intent(in) :: natom
    _REAL_,  intent(in) :: coord(3,natom)
    integer, intent(in) :: numBonds(natom)
    integer, intent(in) :: bondedAtoms(natom,natom)
    type(corInfoType), intent(out) :: cct

    integer :: iat, jat, iatnum, jatnum
    _REAL_ :: rab, vec(3)
    logical :: ccBond
    _REAL_, parameter :: threshold = 1.65d0
    _REAL_, parameter :: factor = 6.0d0  !  (The value "6" was determined empirically by Jimmy Stewart)

    call init(cct)
    do iat = 1, natom

       iatnum = qmmm_struct%iqm_atomic_numbers(iat)
       if ( iatnum == 6 .and. numBonds(iat) == 2) then

          ! Possible CC triple bond
          ccBond = .false.
          rab = 10.d0
          jat = bondedAtoms(1, iat)
          jatnum = qmmm_struct%iqm_atomic_numbers(jat)
          if (jatnum == 6) then
             ccBond = .true.
             vec(1:3) = coord(1:3,iat) - coord(1:3,jat)
             rab = min(rab, vec(1)**2 + vec(2)**2 + vec(3)**2)
          end if
          jat = bondedAtoms(2, iat)
          jatnum = qmmm_struct%iqm_atomic_numbers(jat)
          if(jatnum == 6) then
             ccBond = .true.
             vec(1:3) = coord(1:3,iat) - coord(1:3,jat)
             rab = min(rab, vec(1)**2 + vec(2)**2 + vec(3)**2)
          end if
          if (.not. ccBond) cycle
          if (rab > threshold) cycle

          cct%natom = cct%natom + 1
       end if

    end do

    if (cct%natom > 0) then
       cct%inUse = .true.
       cct%energy = dble(cct%natom) * factor
    end if

  end subroutine  ccTripleBond

  ! ----------------------------------------------------
  ! Molecular mechanics correction to all nitrogen atoms
  ! that have exactly three ligands
  ! ----------------------------------------------------
  subroutine nsp2Correction(natom, coord, numBonds, bondedAtoms, nsp2)

    use qmmm_module, only : qmmm_struct
    implicit none

    integer, intent(in) :: natom
    _REAL_,  intent(in) :: coord(3,natom)
    integer, intent(in) :: numBonds(natom)
    integer, intent(in) :: bondedAtoms(natom,natom)
    type(corInfoType), intent(out) :: nsp2

    integer :: iat, iatnum, jat(3), jatnum(3), numHydrogen

    call init(nsp2)

    do iat = 1, natom
       iatnum = qmmm_struct%iqm_atomic_numbers(iat)
       if (iatnum == 7 .and. numBonds(iat) == 3) then
          jat(1:3) = bondedAtoms(1:3,iat)
          jatnum(1) = qmmm_struct%iqm_atomic_numbers(jat(1))
          jatnum(2) = qmmm_struct%iqm_atomic_numbers(jat(2))
          jatnum(3) = qmmm_struct%iqm_atomic_numbers(jat(3))
          numHydrogen = 0
          if ( jatnum(1) == 1) numHydrogen = numHydrogen + 1
          if ( jatnum(2) == 1) numHydrogen = numHydrogen + 1
          if ( jatnum(3) == 1) numHydrogen = numHydrogen + 1
          if ( numHydrogen < 2) then
             nsp2%natom = nsp2%natom + 1
             nsp2%energy = nsp2%energy &
                  + nsp2AtomCorrection(coord, iat, jat(1), jat(2), jat(3))
          end if
       end if
    end do

    if (nsp2%natom > 0) then
       nsp2%inUse = .true.
    end if

  end subroutine nsp2Correction

  ! --------------------------------------------------------
  ! Penalty for triple coordinated nitrogen being non-planar
  ! Routine from Jimmy Stewart
  ! --------------------------------------------------------
  _REAL_ function nsp2AtomCorrection(coord,n,i,j,k)

    use constants, only : half, one, four, ten
    implicit none
    _REAL_,  intent(in) :: coord(3,*)
    integer, intent(in) :: n, i, j, k

    _REAL_ :: a, b, c, ab, ac, bc, tot, cosa, cosb, cosc
    _REAL_ :: vec(3)

    ! Evaluate the penalty for non-planarity
    ! - done by working out the three angles about the central atom
    ! (here "n") subtended by the lines to atoms "i", "j", and "k".
    vec(1:3) = coord(1:3,n) - coord(1:3,i)
    a = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    vec(1:3) = coord(1:3,n) - coord(1:3,j)
    b = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    vec(1:3) = coord(1:3,n) - coord(1:3,k)
    c = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    vec(1:3) = coord(1:3,j) - coord(1:3,i)
    ab = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    vec(1:3) = coord(1:3,k) - coord(1:3,i)
    ac = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    vec(1:3) = coord(1:3,j) - coord(1:3,k)
    bc = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

    cosa = acos((b**2 +c**2 -bc**2)/(2.d0*b*c))
    cosb = acos((a**2 +c**2 -ac**2)/(2.d0*a*c))
    cosc = acos((b**2 +a**2 -ab**2)/(2.d0*b*a))

    ! tot = difference between the sum of the three angles and
    !       360 degrees, expressed as radians.
    tot = four*asin(one) - (cosa + cosb + cosc)
    nsp2AtomCorrection = -half*exp(-ten*tot)

  end function nsp2AtomCorrection

  ! ---------------------------------
  ! Gradient of PM6 correction to HOF
  ! ---------------------------------
  subroutine hofCorrectionGradient(natom, dxyz)

    use qmmm_module, only : qmmm_struct
    implicit none

    integer, intent(in)    :: natom
    _REAL_,  intent(inout) :: dxyz(3,natom)

    integer :: numBonds(natom)
    integer :: bondedAtoms(natom,natom)

    ! Determine which atoms are bonded to each other
    call setupDentate(natom, qmmm_struct%qm_coords, numBonds, bondedAtoms)

    ! gradient due to 
    call nsp2CorrectionGrad(natom, qmmm_struct%qm_coords, numBonds, bondedAtoms, dxyz)

  end subroutine hofCorrectionGradient

  ! -------------------------------------------------
  ! Gradient of molecular mechanics correction to all 
  ! nitrogen atoms that have exactly three ligands
  ! -------------------------------------------------
  subroutine nsp2CorrectionGrad(natom, coord, numBonds, bondedAtoms, dxyz)

    use qmmm_module, only : qmmm_struct
    implicit none

    integer, intent(in)   :: natom
    _REAL_,  intent(in)   :: coord(3,natom)
    integer, intent(in)   :: numBonds(natom)
    integer, intent(in)   :: bondedAtoms(natom,natom)
    _REAL_, intent(inout) :: dxyz(3,natom)

    integer :: iat, iatnum, jat(3), jatnum(3), numHydrogen

    do iat = 1, natom
       iatnum = qmmm_struct%iqm_atomic_numbers(iat)
       if (iatnum == 7 .and. numBonds(iat) == 3) then
          jat(1:3) = bondedAtoms(1:3,iat)
          jatnum(1) = qmmm_struct%iqm_atomic_numbers(jat(1))
          jatnum(2) = qmmm_struct%iqm_atomic_numbers(jat(2))
          jatnum(3) = qmmm_struct%iqm_atomic_numbers(jat(3))
          numHydrogen = 0
          if ( jatnum(1) == 1) numHydrogen = numHydrogen + 1
          if ( jatnum(2) == 1) numHydrogen = numHydrogen + 1
          if ( jatnum(3) == 1) numHydrogen = numHydrogen + 1
          if ( numHydrogen < 2) then
             ! call nsp2AtomCorrectionGrad(coord, iat, jat(1), jat(2), jat(3), dxyz)
             call nsp2AtomCorrectionGradNum(coord, iat, jat(1), jat(2), jat(3), dxyz)
          end if
       end if
    end do

  end subroutine nsp2CorrectionGrad

  ! --------------------------------------------------------------
  ! Analytical gradient of penalty function for triple coordinated 
  ! nitrogen being non-planar
  ! CURRENTLY DEFUNCT, DO NOT USE!!!
  ! NEEDS DEBUGGING
  ! --------------------------------------------------------------
  subroutine nsp2AtomCorrectionGrad(coord, n, i, j, k, dxyz)

    use constants, only : zero, one, two, four, five, ten
    implicit none
    _REAL_,  intent(in)    :: coord(3,*)
    integer, intent(in)    :: n, i, j, k
    _REAL_,  intent(inout) :: dxyz(3,*)

    _REAL_ :: ea, eb, ec, ab, ac, bc, tot, cosa, cosb, cosc
    _REAL_ :: vea(3), veb(3), vec(3), vab(3), vac(3), vbc(3)
    _REAL_ :: tmp, tmp1, tmp2, ua, ub, uc
    _REAL_ :: grad(3,4) ! n/e=1, i/a=2, j/b=3, k/c=4

    grad(3,4) = zero

    ! PROBABLY THERE IS A SIGN ERROR HERE OR LATER
    vea(1:3) = coord(1:3,n) - coord(1:3,i)
    veb(1:3) = coord(1:3,n) - coord(1:3,j)
    vec(1:3) = coord(1:3,n) - coord(1:3,k)
    vab(1:3) = coord(1:3,j) - coord(1:3,i)
    vac(1:3) = coord(1:3,k) - coord(1:3,i)
    vbc(1:3) = coord(1:3,k) - coord(1:3,j)

    ea = sqrt(vea(1)**2 + vea(2)**2 + vea(3)**2)
    eb = sqrt(veb(1)**2 + veb(2)**2 + veb(3)**2)
    ec = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    ab = sqrt(vab(1)**2 + vab(2)**2 + vab(3)**2)
    ac = sqrt(vac(1)**2 + vac(2)**2 + vac(3)**2)
    bc = sqrt(vbc(1)**2 + vbc(2)**2 + vbc(3)**2)

    ua = (eb**2 +ec**2 -bc**2)/(two*eb*ec)
    ub = (ea**2 +ec**2 -ac**2)/(two*ea*ec)
    uc = (eb**2 +ea**2 -ab**2)/(two*eb*ea)

    cosa = acos(ua)
    cosb = acos(ub)
    cosc = acos(uc)

    ! part 1: exponential prefactor, common to all
    tmp = four*asin(one) - (cosa + cosb + cosc)
    tmp1 = five*exp(-ten*tmp)

    ! part 2a: angle alpha derivative, that is, cosa
    tmp2 = -one / sqrt(one-ua**2)
    tmp = tmp1 * tmp2 * one / (eb*ec)
    grad(1:3,3) = grad(1:3,3) + tmp*(-vbc(1:3) -veb(1:3))
    grad(1:3,4) = grad(1:3,4) + tmp*( vbc(1:3) -vec(1:3))
    tmp = tmp1 * tmp2 * ua / (eb**2)
    grad(1:3,3) = grad(1:3,3) + tmp*veb(1:3)
    tmp = tmp1 * tmp2 * ua / (ec**2)
    grad(1:3,4) = grad(1:3,4) + tmp*vec(1:3)

    ! part 2b: angle beta derivative, that is, cosb
    tmp2 = -one / sqrt(one-ub**2)
    tmp = tmp1 * tmp2 * one / (ea*ec)
    grad(1:3,2) = grad(1:3,2) + tmp*(-vac(1:3) -vea(1:3))
    grad(1:3,4) = grad(1:3,4) + tmp*( vac(1:3) -vec(1:3))
    tmp = tmp1 * tmp2 * ub / (ea**2)
    grad(1:3,2) = grad(1:3,2) + tmp*vea(1:3)
    tmp = tmp1 * tmp2 * ub / (ec**2)
    grad(1:3,4) = grad(1:3,4) + tmp*vec(1:3)

    ! part 2c: angle gamma derivative, that is, cosc
    tmp2 = -one / sqrt(one-uc**2)
    tmp = tmp1 * tmp2 * one / (ea*eb)
    grad(1:3,2) = grad(1:3,2) + tmp*(-vab(1:3) -vea(1:3))
    grad(1:3,3) = grad(1:3,3) + tmp*( vab(1:3) -veb(1:3))
    tmp = tmp1 * tmp2 * uc / (ea**2)
    grad(1:3,2) = grad(1:3,2) + tmp*vea(1:3)
    tmp = tmp1 * tmp2 * uc / (eb**2)
    grad(1:3,3) = grad(1:3,3) + tmp*veb(1:3)

    ! get n/e gradient from ijk/abc gradient
    grad(1:3,1) = - (grad(1:3,2) + grad(1:3,3) + grad(1:3,4))

    write(6,'(a)') ' *** Nsp2 MM GRADIENT CORRECTION ***'
    write(6,'(i3,3f20.4)') n, grad(1:3,1)
    write(6,'(i3,3f20.4)') i, grad(1:3,2)
    write(6,'(i3,3f20.4)') j, grad(1:3,3)
    write(6,'(i3,3f20.4)') k, grad(1:3,4)

    ! add gradient to gradient vector dxyz
    dxyz(1:3,n) = dxyz(1:3,n) + grad(1:3,1)
    dxyz(1:3,i) = dxyz(1:3,i) + grad(1:3,2)
    dxyz(1:3,j) = dxyz(1:3,j) + grad(1:3,3)
    dxyz(1:3,k) = dxyz(1:3,k) + grad(1:3,4)

  end subroutine nsp2AtomCorrectionGrad

  ! --------------------------------------------------------------
  ! Numerical gradient of penalty function for triple coordinated 
  ! nitrogen being non-planar
  ! using 3 point central differences with a step size of 1.d-5 A
  ! -------------------------------------------------------------
  subroutine nsp2AtomCorrectionGradNum(coord, n, i, j, k, dxyz)

    use constants, only : zero, half, one, two, four, ten
    implicit none
    _REAL_,  intent(in)    :: coord(3,*)
    integer, intent(in)    :: n, i, j, k
    _REAL_,  intent(inout) :: dxyz(3,*)

    integer :: iat, ixyz, istep
    _REAL_  :: fac
    _REAL_  :: xyz(3,4), grad(3,4) ! n/e=1, i/a=2, j/b=3, k/c=4
    _REAL_  :: backup
    _REAL_  :: vec(3)
    _REAL_  :: a, b, c, ab, ac, bc, cosa, cosb, cosc, tot
    _REAL_  :: ener(2)
    _REAL_, parameter :: stepsize = 1.0d-5

    xyz(1:3,1) = coord(1:3,n)
    xyz(1:3,2) = coord(1:3,i)
    xyz(1:3,3) = coord(1:3,j)
    xyz(1:3,4) = coord(1:3,k)

    grad(:,:) = zero
    ener(:)   = zero

    ! gradient for atoms i, j, k
    do iat = 2, 4
       do ixyz = 1, 3
          do istep = 1, 2

             fac = one
             if (istep == 2) then
                fac = -one
             end if

             backup = xyz(ixyz,iat)
             xyz(ixyz,iat) = xyz(ixyz,iat) + fac*stepsize
          
             vec(1:3) = xyz(1:3,1) - xyz(1:3,2)
             a = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             vec(1:3) = xyz(1:3,1) - xyz(1:3,3)
             b = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             vec(1:3) = xyz(1:3,1) - xyz(1:3,4)
             c = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             vec(1:3) = xyz(1:3,3) - xyz(1:3,2)
             ab = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             vec(1:3) = xyz(1:3,4) - xyz(1:3,2)
             ac = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             vec(1:3) = xyz(1:3,3) - xyz(1:3,4)
             bc = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
             
             cosa = acos((b**2 +c**2 -bc**2)/(two*b*c))
             cosb = acos((a**2 +c**2 -ac**2)/(two*a*c))
             cosc = acos((b**2 +a**2 -ab**2)/(two*b*a))
             
             ! tot = difference between the sum of the three angles and
             !       360 degrees, expressed as radians.
             tot = four*asin(one) - (cosa + cosb + cosc)
             ener(istep) = -half*exp(-ten*tot)

             !write(6,*) 'atom=',iat, '  coord=',ixyz, '  step=',istep, '  fac=',fac, 'ener=',ener(istep)

             ! restore geometry
             xyz(ixyz,iat) = backup
             
          end do
          
          grad(ixyz,iat) = (ener(1) - ener(2)) / (two*stepsize)

       end do
    end do

    ! get n/e gradient from ijk/abc gradient
    grad(1:3,1) = - (grad(1:3,2) + grad(1:3,3) + grad(1:3,4))

    ! write(6,'(a)') ' *** Nsp2 MM GRADIENT CORRECTION ***'
    ! write(6,'(i3,3f20.16)') n, grad(1:3,1)
    ! write(6,'(i3,3f20.16)') i, grad(1:3,2)
    ! write(6,'(i3,3f20.16)') j, grad(1:3,3)
    ! write(6,'(i3,3f20.16)') k, grad(1:3,4)

    ! add gradient to gradient vector dxyz
    dxyz(1:3,n) = dxyz(1:3,n) + grad(1:3,1)
    dxyz(1:3,i) = dxyz(1:3,i) + grad(1:3,2)
    dxyz(1:3,j) = dxyz(1:3,j) + grad(1:3,3)
    dxyz(1:3,k) = dxyz(1:3,k) + grad(1:3,4)

  end subroutine nsp2AtomCorrectionGradNum

  ! ----------------------
  ! initialize corInfoType
  ! ----------------------
  subroutine init(self)
    use constants, only : zero
    implicit none
    type(corInfoType), intent(out) :: self
    self%inUse  = .false.
    self%natom  = 0
    self%energy = zero
  end subroutine init

  ! -----------------
  ! print CorInfoType
  ! -----------------
  subroutine printCorInfoType(self, string)
    use constants, only : KCAL_TO_EV
    implicit none
    type(corInfoType), intent(in) :: self
    character(len=strlen), intent(in) :: string
    if (self%inUse) then
       write(6,'(/3a,i3,a)') ' PM6: ', trim(string), ' in use for', &
            self%natom, ' atoms'
       write(6,'(a,f20.8,a,f18.8,a)') ' Correction          =', &
            self%energy, ' kcal/mol  (', self%energy*KCAL_TO_EV, ' eV)'
    end if
  end subroutine printCorInfoType

end module qm2_pm6_hof_module
