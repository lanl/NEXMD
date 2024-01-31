#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++
!This module contains constants and indexes related
!to elements, orbitals, and two center integrals
!Created by Taisung Lee, 01/27/2011
!+++++++++++++++++++++++++++++++++++++++

module ElementOrbitalIndex

    use UtilitiesModule, only : Upcase

    implicit none

    integer, parameter :: NumberElements = 86

    character(2), dimension(1:NumberElements), parameter :: ElementSymbol = (/ &
        'H ', 'He', &
        'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
        'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
        'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
        'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
        'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
        'Cs', 'Ba', 'La', &
        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
        'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'/)
!
!  The following constants are for semiempirical methods only.  Should be separated
!  into another module later
!  N=3, L=2 --  Up to d-orbital for now
!
    integer, parameter :: MaxPrincipleQuantumNumber = 6
    integer, parameter :: MaxAngularQuantumNumber = 2
    integer, parameter :: MaxValenceOrbitals = (MaxAngularQuantumNumber + 1)**2
    integer, parameter :: MaxValenceDimension = MaxValenceOrbitals*(MaxValenceOrbitals + 1)/2

    integer, parameter :: MaxGaussianExpansion = 6 ! STO-6G expansion

! Principal Quantum Numbers for AO.
    integer, parameter, dimension(1:NumberElements) :: SPPrincipalQuantumNumber = (/ &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, &
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, &
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, &
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
        6, 6, 6, 6, 6, 6/)
    integer, parameter, dimension(1:NumberElements) :: DPrincipalQuantumNumber = (/ &
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
        4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
        6, 6, 6, 6, 6, 6/)

! Indexes used in arranging 2-electron integrals
    integer, parameter :: Index1_2Electron = 243

    integer, parameter, dimension(1:Index1_2Electron) :: IntIJ = (/ &
        1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, &
        4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, &
        9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, &
        13, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, &
        16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, &
        20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, &
        22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, &
        26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, &
        29, 29, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, &
        34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, &
        36, 37, 37, 37, 37, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 40, 40, 40, 41, 42, &
        42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 45, &
        45, 45, 45/)

    integer, parameter, dimension(1:Index1_2Electron) :: IntKL = (/ &
        15, 21, 28, 36, 45, 12, 19, 23, 39, 11, 15, 21, 22, 26, 28, 36, 45, 13, 24, 32, &
        38, 34, 37, 43, 11, 15, 21, 22, 26, 28, 36, 45, 17, 25, 31, 16, 20, 27, 44, 29, &
        33, 35, 42, 15, 21, 22, 28, 36, 45, 3, 6, 11, 21, 26, 36, 2, 12, 19, 23, 39, &
        4, 13, 24, 32, 38, 14, 17, 31, 1, 3, 6, 10, 15, 21, 22, 28, 36, 45, 8, 16, &
        20, 27, 44, 7, 14, 17, 25, 31, 18, 30, 40, 2, 12, 19, 23, 39, 8, 16, 20, 27, &
        44, 1, 3, 6, 10, 11, 15, 21, 22, 26, 28, 36, 45, 3, 6, 10, 15, 21, 22, 28, &
        36, 45, 2, 12, 19, 23, 39, 4, 13, 24, 32, 38, 7, 17, 25, 31, 3, 6, 11, 21, &
        26, 36, 8, 16, 20, 27, 44, 1, 3, 6, 10, 15, 21, 22, 28, 36, 45, 9, 29, 33, &
        35, 42, 18, 30, 40, 7, 14, 17, 25, 31, 4, 13, 24, 32, 38, 9, 29, 33, 35, 42, &
        5, 34, 37, 43, 9, 29, 33, 35, 42, 1, 3, 6, 10, 11, 15, 21, 22, 26, 28, 36, &
        45, 5, 34, 37, 43, 4, 13, 24, 32, 38, 2, 12, 19, 23, 39, 18, 30, 40, 41, 9, &
        29, 33, 35, 42, 5, 34, 37, 43, 8, 16, 20, 27, 44, 1, 3, 6, 10, 15, 21, 22, &
        28, 36, 45/)

    integer, parameter, dimension(1:Index1_2Electron) :: IntRep = (/ &
        1, 1, 1, 1, 1, 3, 3, 8, 3, 9, 6, 6, 12, 14, 13, 7, 6, 15, 8, 3, &
        3, 11, 9, 14, 17, 6, 7, 12, 18, 13, 6, 6, 3, 2, 3, 9, 11, 10, 11, 9, &
        16, 10, 11, 7, 6, 4, 5, 6, 7, 9, 17, 19, 32, 22, 40, 3, 33, 34, 27, 46, &
        15, 33, 28, 41, 47, 35, 35, 42, 1, 6, 6, 7, 29, 38, 22, 31, 38, 51, 9, 19, &
        32, 21, 32, 3, 35, 33, 24, 34, 35, 35, 35, 3, 34, 33, 26, 34, 11, 32, 44, 37, &
        49, 1, 6, 7, 6, 32, 38, 29, 21, 39, 30, 38, 38, 12, 12, 4, 22, 21, 19, 20, &
        21, 22, 8, 27, 26, 25, 27, 8, 28, 25, 26, 27, 2, 24, 23, 24, 14, 18, 22, 39, &
        48, 45, 10, 21, 37, 36, 37, 1, 13, 13, 5, 31, 30, 20, 29, 30, 31, 9, 19, 40, &
        21, 32, 35, 35, 35, 3, 42, 34, 24, 33, 3, 41, 26, 33, 34, 16, 40, 44, 43, 50, &
        11, 44, 32, 39, 10, 21, 43, 36, 37, 1, 7, 6, 6, 40, 38, 38, 21, 45, 30, 29, &
        38, 9, 32, 19, 22, 3, 47, 27, 34, 33, 3, 46, 34, 27, 33, 35, 35, 35, 52, 11, &
        32, 50, 37, 44, 14, 39, 22, 48, 11, 32, 49, 37, 44, 1, 6, 6, 7, 51, 38, 22, &
        31, 38, 29/)

    integer, parameter, dimension(1:Index1_2Electron) ::  IntRf1 = (/ &
        19, 19, 19, 19, 19, 3, 3, 8, 3, 3, 33, 33, 8, 27, 25, 35, 33, 15, 8, 3, &
        3, 34, 3, 27, 15, 33, 35, 8, 28, 25, 33, 33, 3, 2, 3, 3, 34, 24, 35, 3, &
        41, 26, 35, 35, 33, 2, 23, 33, 35, 3, 15, 1, 32, 22, 40, 3, 6, 11, 14, 0, &
        15, 6, 18, 16, 0, 7, 11, 16, 19, 33, 33, 35, 29, 44, 22, 48, 44, 52, 3, 1, &
        32, 21, 32, 3, 11, 6, 10, 11, 7, 11, 11, 3, 11, 6, 10, 11, 34, 32, 38, 37, &
        50, 19, 33, 35, 33, 32, 44, 29, 21, 37, 36, 44, 44, 8, 8, 2, 22, 21, 1, 20, &
        21, 22, 8, 14, 10, 13, 14, 8, 18, 13, 10, 14, 2, 10, 5, 10, 27, 28, 22, 37, &
        31, 43, 24, 21, 37, 30, 39, 19, 25, 25, 23, 48, 36, 20, 29, 36, 48, 3, 1, 40, &
        21, 32, 11, 7, 11, 3, 16, 11, 10, 6, 3, 16, 10, 6, 11, 41, 40, 38, 45, 49, &
        34, 38, 32, 37, 26, 21, 45, 30, 37, 19, 35, 33, 33, 40, 44, 44, 21, 43, 36, 29, &
        44, 3, 32, 1, 22, 3, 0, 14, 11, 6, 3, 0, 11, 14, 6, 11, 11, 7, 51, 35, &
        32, 49, 37, 38, 27, 37, 22, 31, 35, 32, 50, 39, 38, 19, 33, 33, 35, 52, 44, 22, &
        48, 44, 29/)

    integer, parameter, dimension(1:Index1_2Electron) ::  IntRf2 = (/ &
        19, 19, 19, 19, 19, 9, 9, 12, 9, 3, 33, 33, 8, 27, 25, 35, 33, 17, 12, 9, &
        9, 35, 3, 27, 15, 33, 35, 8, 28, 25, 33, 33, 9, 4, 9, 3, 35, 26, 34, 3, &
        42, 24, 34, 35, 33, 2, 23, 33, 35, 3, 15, 19, 32, 22, 40, 9, 33, 35, 27, 47, &
        17, 33, 28, 42, 46, 35, 34, 41, 19, 33, 33, 35, 29, 44, 22, 48, 44, 52, 3, 19, &
        32, 21, 32, 9, 34, 33, 26, 35, 35, 34, 34, 9, 35, 33, 24, 35, 35, 32, 44, 39, &
        0, 19, 33, 35, 33, 32, 44, 29, 21, 37, 36, 44, 44, 8, 8, 2, 22, 21, 19, 20, &
        21, 22, 12, 27, 24, 25, 27, 12, 28, 25, 24, 27, 4, 26, 23, 26, 27, 28, 22, 37, &
        48, 43, 26, 21, 39, 36, 37, 19, 25, 25, 23, 48, 36, 20, 29, 36, 48, 3, 19, 40, &
        21, 32, 34, 35, 34, 9, 41, 35, 26, 33, 9, 42, 24, 33, 35, 42, 40, 44, 43, 0, &
        35, 44, 32, 37, 24, 21, 43, 36, 39, 19, 35, 33, 33, 40, 44, 44, 21, 43, 36, 29, &
        44, 3, 32, 19, 22, 9, 46, 27, 35, 33, 9, 47, 35, 27, 33, 34, 34, 35, 52, 34, &
        32, 0, 39, 44, 27, 37, 22, 48, 34, 32, 0, 37, 44, 19, 33, 33, 35, 52, 44, 22, &
        48, 44, 29/)

contains

    function GetAtomicNumber(inputSymbol) result(atomicNumber)

        implicit none

        character(2), intent(in) :: inputSymbol
        integer :: atomicNumber

        integer :: i

        atomicNumber = -1
        do i = 1, NumberElements
            if (Upcase(ElementSymbol(i)) == Upcase(inputSymbol)) then
                atomicNumber = i
                exit
            end if
        end do

    end function GetAtomicNumber

    subroutine Is2CenterIntegralZero(r2center)

        implicit none

        logical, intent(out) :: r2center(45, 45)

        logical, save :: initialized = .false.
        logical, save :: r2saved(45, 45)

        if (initialized) then
            r2center = r2saved
            return
        else
            r2center = .false.

            r2center(1, 1) = .TRUE.
            r2center(2, 1) = .TRUE.
            r2center(3, 1) = .TRUE.
            r2center(6, 1) = .TRUE.
            r2center(10, 1) = .TRUE.
            r2center(11, 1) = .TRUE.
            r2center(12, 1) = .TRUE.
            r2center(15, 1) = .TRUE.
            r2center(18, 1) = .TRUE.
            r2center(21, 1) = .TRUE.
            r2center(25, 1) = .TRUE.
            r2center(28, 1) = .TRUE.
            r2center(36, 1) = .TRUE.
            r2center(45, 1) = .TRUE.
            r2center(1, 2) = .TRUE.
            r2center(2, 2) = .TRUE.
            r2center(3, 2) = .TRUE.
            r2center(6, 2) = .TRUE.
            r2center(10, 2) = .TRUE.
            r2center(11, 2) = .TRUE.
            r2center(12, 2) = .TRUE.
            r2center(15, 2) = .TRUE.
            r2center(18, 2) = .TRUE.
            r2center(21, 2) = .TRUE.
            r2center(25, 2) = .TRUE.
            r2center(28, 2) = .TRUE.
            r2center(36, 2) = .TRUE.
            r2center(45, 2) = .TRUE.
            r2center(1, 3) = .TRUE.
            r2center(2, 3) = .TRUE.
            r2center(3, 3) = .TRUE.
            r2center(6, 3) = .TRUE.
            r2center(10, 3) = .TRUE.
            r2center(11, 3) = .TRUE.
            r2center(12, 3) = .TRUE.
            r2center(15, 3) = .TRUE.
            r2center(18, 3) = .TRUE.
            r2center(21, 3) = .TRUE.
            r2center(25, 3) = .TRUE.
            r2center(28, 3) = .TRUE.
            r2center(36, 3) = .TRUE.
            r2center(45, 3) = .TRUE.
            r2center(4, 4) = .TRUE.
            r2center(5, 4) = .TRUE.
            r2center(13, 4) = .TRUE.
            r2center(16, 4) = .TRUE.
            r2center(17, 4) = .TRUE.
            r2center(20, 4) = .TRUE.
            r2center(31, 4) = .TRUE.
            r2center(34, 4) = .TRUE.
            r2center(40, 4) = .TRUE.
            r2center(43, 4) = .TRUE.
            r2center(4, 5) = .TRUE.
            r2center(5, 5) = .TRUE.
            r2center(13, 5) = .TRUE.
            r2center(16, 5) = .TRUE.
            r2center(17, 5) = .TRUE.
            r2center(20, 5) = .TRUE.
            r2center(31, 5) = .TRUE.
            r2center(34, 5) = .TRUE.
            r2center(40, 5) = .TRUE.
            r2center(43, 5) = .TRUE.
            r2center(1, 6) = .TRUE.
            r2center(2, 6) = .TRUE.
            r2center(3, 6) = .TRUE.
            r2center(6, 6) = .TRUE.
            r2center(10, 6) = .TRUE.
            r2center(11, 6) = .TRUE.
            r2center(12, 6) = .TRUE.
            r2center(15, 6) = .TRUE.
            r2center(18, 6) = .TRUE.
            r2center(21, 6) = .TRUE.
            r2center(25, 6) = .TRUE.
            r2center(28, 6) = .TRUE.
            r2center(29, 6) = .TRUE.
            r2center(33, 6) = .TRUE.
            r2center(36, 6) = .TRUE.
            r2center(45, 6) = .TRUE.
            r2center(7, 7) = .TRUE.
            r2center(8, 7) = .TRUE.
            r2center(14, 7) = .TRUE.
            r2center(22, 7) = .TRUE.
            r2center(23, 7) = .TRUE.
            r2center(26, 7) = .TRUE.
            r2center(32, 7) = .TRUE.
            r2center(35, 7) = .TRUE.
            r2center(39, 7) = .TRUE.
            r2center(42, 7) = .TRUE.
            r2center(7, 8) = .TRUE.
            r2center(8, 8) = .TRUE.
            r2center(14, 8) = .TRUE.
            r2center(22, 8) = .TRUE.
            r2center(23, 8) = .TRUE.
            r2center(26, 8) = .TRUE.
            r2center(32, 8) = .TRUE.
            r2center(35, 8) = .TRUE.
            r2center(39, 8) = .TRUE.
            r2center(42, 8) = .TRUE.
            r2center(9, 9) = .TRUE.
            r2center(27, 9) = .TRUE.
            r2center(37, 9) = .TRUE.
            r2center(41, 9) = .TRUE.
            r2center(1, 10) = .TRUE.
            r2center(2, 10) = .TRUE.
            r2center(3, 10) = .TRUE.
            r2center(6, 10) = .TRUE.
            r2center(10, 10) = .TRUE.
            r2center(11, 10) = .TRUE.
            r2center(12, 10) = .TRUE.
            r2center(15, 10) = .TRUE.
            r2center(18, 10) = .TRUE.
            r2center(21, 10) = .TRUE.
            r2center(25, 10) = .TRUE.
            r2center(28, 10) = .TRUE.
            r2center(29, 10) = .TRUE.
            r2center(33, 10) = .TRUE.
            r2center(36, 10) = .TRUE.
            r2center(45, 10) = .TRUE.
            r2center(1, 11) = .TRUE.
            r2center(2, 11) = .TRUE.
            r2center(3, 11) = .TRUE.
            r2center(6, 11) = .TRUE.
            r2center(10, 11) = .TRUE.
            r2center(11, 11) = .TRUE.
            r2center(12, 11) = .TRUE.
            r2center(15, 11) = .TRUE.
            r2center(18, 11) = .TRUE.
            r2center(21, 11) = .TRUE.
            r2center(25, 11) = .TRUE.
            r2center(28, 11) = .TRUE.
            r2center(36, 11) = .TRUE.
            r2center(45, 11) = .TRUE.
            r2center(1, 12) = .TRUE.
            r2center(2, 12) = .TRUE.
            r2center(3, 12) = .TRUE.
            r2center(6, 12) = .TRUE.
            r2center(10, 12) = .TRUE.
            r2center(11, 12) = .TRUE.
            r2center(12, 12) = .TRUE.
            r2center(15, 12) = .TRUE.
            r2center(18, 12) = .TRUE.
            r2center(21, 12) = .TRUE.
            r2center(25, 12) = .TRUE.
            r2center(28, 12) = .TRUE.
            r2center(36, 12) = .TRUE.
            r2center(45, 12) = .TRUE.
            r2center(4, 13) = .TRUE.
            r2center(5, 13) = .TRUE.
            r2center(13, 13) = .TRUE.
            r2center(16, 13) = .TRUE.
            r2center(17, 13) = .TRUE.
            r2center(20, 13) = .TRUE.
            r2center(31, 13) = .TRUE.
            r2center(34, 13) = .TRUE.
            r2center(40, 13) = .TRUE.
            r2center(43, 13) = .TRUE.
            r2center(7, 14) = .TRUE.
            r2center(8, 14) = .TRUE.
            r2center(14, 14) = .TRUE.
            r2center(22, 14) = .TRUE.
            r2center(23, 14) = .TRUE.
            r2center(26, 14) = .TRUE.
            r2center(32, 14) = .TRUE.
            r2center(35, 14) = .TRUE.
            r2center(39, 14) = .TRUE.
            r2center(42, 14) = .TRUE.
            r2center(1, 15) = .TRUE.
            r2center(2, 15) = .TRUE.
            r2center(3, 15) = .TRUE.
            r2center(6, 15) = .TRUE.
            r2center(10, 15) = .TRUE.
            r2center(11, 15) = .TRUE.
            r2center(12, 15) = .TRUE.
            r2center(15, 15) = .TRUE.
            r2center(18, 15) = .TRUE.
            r2center(21, 15) = .TRUE.
            r2center(25, 15) = .TRUE.
            r2center(28, 15) = .TRUE.
            r2center(36, 15) = .TRUE.
            r2center(45, 15) = .TRUE.
            r2center(4, 16) = .TRUE.
            r2center(5, 16) = .TRUE.
            r2center(13, 16) = .TRUE.
            r2center(16, 16) = .TRUE.
            r2center(17, 16) = .TRUE.
            r2center(20, 16) = .TRUE.
            r2center(31, 16) = .TRUE.
            r2center(34, 16) = .TRUE.
            r2center(40, 16) = .TRUE.
            r2center(43, 16) = .TRUE.
            r2center(4, 17) = .TRUE.
            r2center(5, 17) = .TRUE.
            r2center(13, 17) = .TRUE.
            r2center(16, 17) = .TRUE.
            r2center(17, 17) = .TRUE.
            r2center(20, 17) = .TRUE.
            r2center(31, 17) = .TRUE.
            r2center(34, 17) = .TRUE.
            r2center(40, 17) = .TRUE.
            r2center(43, 17) = .TRUE.
            r2center(1, 18) = .TRUE.
            r2center(2, 18) = .TRUE.
            r2center(3, 18) = .TRUE.
            r2center(6, 18) = .TRUE.
            r2center(10, 18) = .TRUE.
            r2center(11, 18) = .TRUE.
            r2center(12, 18) = .TRUE.
            r2center(15, 18) = .TRUE.
            r2center(18, 18) = .TRUE.
            r2center(21, 18) = .TRUE.
            r2center(25, 18) = .TRUE.
            r2center(28, 18) = .TRUE.
            r2center(29, 18) = .TRUE.
            r2center(33, 18) = .TRUE.
            r2center(36, 18) = .TRUE.
            r2center(45, 18) = .TRUE.
            r2center(4, 20) = .TRUE.
            r2center(5, 20) = .TRUE.
            r2center(13, 20) = .TRUE.
            r2center(16, 20) = .TRUE.
            r2center(17, 20) = .TRUE.
            r2center(20, 20) = .TRUE.
            r2center(31, 20) = .TRUE.
            r2center(34, 20) = .TRUE.
            r2center(40, 20) = .TRUE.
            r2center(43, 20) = .TRUE.
            r2center(1, 21) = .TRUE.
            r2center(2, 21) = .TRUE.
            r2center(3, 21) = .TRUE.
            r2center(6, 21) = .TRUE.
            r2center(10, 21) = .TRUE.
            r2center(11, 21) = .TRUE.
            r2center(12, 21) = .TRUE.
            r2center(15, 21) = .TRUE.
            r2center(18, 21) = .TRUE.
            r2center(21, 21) = .TRUE.
            r2center(25, 21) = .TRUE.
            r2center(28, 21) = .TRUE.
            r2center(29, 21) = .TRUE.
            r2center(33, 21) = .TRUE.
            r2center(36, 21) = .TRUE.
            r2center(45, 21) = .TRUE.
            r2center(7, 22) = .TRUE.
            r2center(8, 22) = .TRUE.
            r2center(14, 22) = .TRUE.
            r2center(22, 22) = .TRUE.
            r2center(23, 22) = .TRUE.
            r2center(26, 22) = .TRUE.
            r2center(32, 22) = .TRUE.
            r2center(35, 22) = .TRUE.
            r2center(39, 22) = .TRUE.
            r2center(42, 22) = .TRUE.
            r2center(7, 23) = .TRUE.
            r2center(8, 23) = .TRUE.
            r2center(14, 23) = .TRUE.
            r2center(22, 23) = .TRUE.
            r2center(23, 23) = .TRUE.
            r2center(26, 23) = .TRUE.
            r2center(32, 23) = .TRUE.
            r2center(35, 23) = .TRUE.
            r2center(39, 23) = .TRUE.
            r2center(42, 23) = .TRUE.
            r2center(1, 25) = .TRUE.
            r2center(2, 25) = .TRUE.
            r2center(3, 25) = .TRUE.
            r2center(6, 25) = .TRUE.
            r2center(10, 25) = .TRUE.
            r2center(11, 25) = .TRUE.
            r2center(12, 25) = .TRUE.
            r2center(15, 25) = .TRUE.
            r2center(18, 25) = .TRUE.
            r2center(21, 25) = .TRUE.
            r2center(25, 25) = .TRUE.
            r2center(28, 25) = .TRUE.
            r2center(29, 25) = .TRUE.
            r2center(33, 25) = .TRUE.
            r2center(36, 25) = .TRUE.
            r2center(45, 25) = .TRUE.
            r2center(7, 26) = .TRUE.
            r2center(8, 26) = .TRUE.
            r2center(14, 26) = .TRUE.
            r2center(22, 26) = .TRUE.
            r2center(23, 26) = .TRUE.
            r2center(26, 26) = .TRUE.
            r2center(32, 26) = .TRUE.
            r2center(35, 26) = .TRUE.
            r2center(39, 26) = .TRUE.
            r2center(42, 26) = .TRUE.
            r2center(9, 27) = .TRUE.
            r2center(27, 27) = .TRUE.
            r2center(37, 27) = .TRUE.
            r2center(41, 27) = .TRUE.
            r2center(1, 28) = .TRUE.
            r2center(2, 28) = .TRUE.
            r2center(3, 28) = .TRUE.
            r2center(6, 28) = .TRUE.
            r2center(10, 28) = .TRUE.
            r2center(11, 28) = .TRUE.
            r2center(12, 28) = .TRUE.
            r2center(15, 28) = .TRUE.
            r2center(18, 28) = .TRUE.
            r2center(21, 28) = .TRUE.
            r2center(25, 28) = .TRUE.
            r2center(28, 28) = .TRUE.
            r2center(29, 28) = .TRUE.
            r2center(33, 28) = .TRUE.
            r2center(36, 28) = .TRUE.
            r2center(45, 28) = .TRUE.
            r2center(6, 29) = .TRUE.
            r2center(10, 29) = .TRUE.
            r2center(18, 29) = .TRUE.
            r2center(21, 29) = .TRUE.
            r2center(25, 29) = .TRUE.
            r2center(28, 29) = .TRUE.
            r2center(29, 29) = .TRUE.
            r2center(33, 29) = .TRUE.
            r2center(4, 31) = .TRUE.
            r2center(5, 31) = .TRUE.
            r2center(13, 31) = .TRUE.
            r2center(16, 31) = .TRUE.
            r2center(17, 31) = .TRUE.
            r2center(20, 31) = .TRUE.
            r2center(31, 31) = .TRUE.
            r2center(34, 31) = .TRUE.
            r2center(40, 31) = .TRUE.
            r2center(43, 31) = .TRUE.
            r2center(7, 32) = .TRUE.
            r2center(8, 32) = .TRUE.
            r2center(14, 32) = .TRUE.
            r2center(22, 32) = .TRUE.
            r2center(23, 32) = .TRUE.
            r2center(26, 32) = .TRUE.
            r2center(32, 32) = .TRUE.
            r2center(35, 32) = .TRUE.
            r2center(39, 32) = .TRUE.
            r2center(42, 32) = .TRUE.
            r2center(6, 33) = .TRUE.
            r2center(10, 33) = .TRUE.
            r2center(18, 33) = .TRUE.
            r2center(21, 33) = .TRUE.
            r2center(25, 33) = .TRUE.
            r2center(28, 33) = .TRUE.
            r2center(29, 33) = .TRUE.
            r2center(33, 33) = .TRUE.
            r2center(4, 34) = .TRUE.
            r2center(5, 34) = .TRUE.
            r2center(13, 34) = .TRUE.
            r2center(16, 34) = .TRUE.
            r2center(17, 34) = .TRUE.
            r2center(20, 34) = .TRUE.
            r2center(31, 34) = .TRUE.
            r2center(34, 34) = .TRUE.
            r2center(40, 34) = .TRUE.
            r2center(43, 34) = .TRUE.
            r2center(7, 35) = .TRUE.
            r2center(8, 35) = .TRUE.
            r2center(14, 35) = .TRUE.
            r2center(22, 35) = .TRUE.
            r2center(23, 35) = .TRUE.
            r2center(26, 35) = .TRUE.
            r2center(32, 35) = .TRUE.
            r2center(35, 35) = .TRUE.
            r2center(39, 35) = .TRUE.
            r2center(42, 35) = .TRUE.
            r2center(1, 36) = .TRUE.
            r2center(2, 36) = .TRUE.
            r2center(3, 36) = .TRUE.
            r2center(6, 36) = .TRUE.
            r2center(10, 36) = .TRUE.
            r2center(11, 36) = .TRUE.
            r2center(12, 36) = .TRUE.
            r2center(15, 36) = .TRUE.
            r2center(18, 36) = .TRUE.
            r2center(21, 36) = .TRUE.
            r2center(25, 36) = .TRUE.
            r2center(28, 36) = .TRUE.
            r2center(36, 36) = .TRUE.
            r2center(45, 36) = .TRUE.
            r2center(9, 37) = .TRUE.
            r2center(27, 37) = .TRUE.
            r2center(37, 37) = .TRUE.
            r2center(41, 37) = .TRUE.
            r2center(7, 39) = .TRUE.
            r2center(8, 39) = .TRUE.
            r2center(14, 39) = .TRUE.
            r2center(22, 39) = .TRUE.
            r2center(23, 39) = .TRUE.
            r2center(26, 39) = .TRUE.
            r2center(32, 39) = .TRUE.
            r2center(35, 39) = .TRUE.
            r2center(39, 39) = .TRUE.
            r2center(42, 39) = .TRUE.
            r2center(4, 40) = .TRUE.
            r2center(5, 40) = .TRUE.
            r2center(13, 40) = .TRUE.
            r2center(16, 40) = .TRUE.
            r2center(17, 40) = .TRUE.
            r2center(20, 40) = .TRUE.
            r2center(31, 40) = .TRUE.
            r2center(34, 40) = .TRUE.
            r2center(40, 40) = .TRUE.
            r2center(43, 40) = .TRUE.
            r2center(9, 41) = .TRUE.
            r2center(27, 41) = .TRUE.
            r2center(37, 41) = .TRUE.
            r2center(41, 41) = .TRUE.
            r2center(7, 42) = .TRUE.
            r2center(8, 42) = .TRUE.
            r2center(14, 42) = .TRUE.
            r2center(22, 42) = .TRUE.
            r2center(23, 42) = .TRUE.
            r2center(26, 42) = .TRUE.
            r2center(32, 42) = .TRUE.
            r2center(35, 42) = .TRUE.
            r2center(39, 42) = .TRUE.
            r2center(42, 42) = .TRUE.
            r2center(4, 43) = .TRUE.
            r2center(5, 43) = .TRUE.
            r2center(13, 43) = .TRUE.
            r2center(16, 43) = .TRUE.
            r2center(17, 43) = .TRUE.
            r2center(20, 43) = .TRUE.
            r2center(31, 43) = .TRUE.
            r2center(34, 43) = .TRUE.
            r2center(40, 43) = .TRUE.
            r2center(43, 43) = .TRUE.
            r2center(1, 45) = .TRUE.
            r2center(2, 45) = .TRUE.
            r2center(3, 45) = .TRUE.
            r2center(6, 45) = .TRUE.
            r2center(10, 45) = .TRUE.
            r2center(11, 45) = .TRUE.
            r2center(12, 45) = .TRUE.
            r2center(15, 45) = .TRUE.
            r2center(18, 45) = .TRUE.
            r2center(21, 45) = .TRUE.
            r2center(25, 45) = .TRUE.
            r2center(28, 45) = .TRUE.
            r2center(36, 45) = .TRUE.
            r2center(45, 45) = .TRUE.

            r2saved = r2center
            initialized = .true.
        end if

    end subroutine Is2CenterIntegralZero

end module ElementOrbitalIndex
