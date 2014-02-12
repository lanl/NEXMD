#include "dprec.fh"
#include "assert.fh"

module naesmd_constants
    implicit none
! evkcal converts eV in Kcal/mol
    _REAL_, parameter :: evkcal=23.060529d0
    _REAL_, parameter :: kcalev=1.0d0/evkcal
! kcalau converts Kcal/mol in a.u.
    _REAL_, parameter :: kcalau=1.5936013d-3
! feVmdqt converts a.u. in eV
    _REAL_, parameter :: feVmdqt=27.2116d0
! convl transform atomic units in amstrong
    _REAL_, parameter :: convl=0.529177249d0
! definition of phi (didn't find any references to it - KGB)
!   _REAL_, parameter :: phi=dacos(-1.0d0)
! boltzman   Boltzmann constant (kB) in  au ! g*Ang**2/ps**2/K/mole
    _REAL_, parameter :: boltzman=0.83143435d0*1.59360d-3/4.184d2
!  1/convtf transform the time from femtoseconds to a.u.
    _REAL_, parameter :: convtf=2.41888D-2
! 1/convt transform from picoseconds to atomic units
    _REAL_, parameter :: convt=2.41888d-5
! convm transfrom masses to atomic units
    _REAL_, parameter :: convm=1822.8885d0


    CHARACTER*2, parameter :: ELEMNT(107) = (/' H','He', &
     'Li','Be',' B',' C',' N',' O',' F','Ne', &
     'Na','Mg','Al','Si',' P',' S','Cl','Ar', &
     ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu', &
     'Zn','Ga','Ge','As','Se','Br','Kr', &
     'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag', &
     'Cd','In','Sn','Sb','Te',' I','Xe', &
     'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
     'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt', &
     'Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
     'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','XX', &
     'Fm','Md','Cb','++',' +','--',' -','Tv'/)
end module 
