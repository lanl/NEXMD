!+ Specification and control of Amber's working precision

#ifndef DPREC_H
#define DPREC_H

! Description:
! Preprocessor directives that characterize the floating-point
! working precision as single or double precision.
! The current scheme guarantees internal consistency at the expense
! of flexibility.  A need for flexibility has yet to appear.
! The preprocessor guard DPREC_H prevents multiple, and thus
! inconsistent, definitions.
! The default working precision is double precision.
! User control of the working precision at build time should be
! exercised via the preprocessor name _REAL_.
! To build a single precision Amber use
!     make -e AMBERBUILDFLAGS=' -D_REAL_ '
! The preprocessor names that characterize the precision are

!   _REAL_     precision type specifier.
!              Use  _REAL_ foo  as the precision independent
!              notation for  double precision foo  and  real foo.

!   AMBER_MPI_REAL
!              MPI precision type specifier.
!              Use AMBER_MPI_REAL as the precision independent
!              notation for MPI_DOUBLE_PRECISION and MPI_REAL.

!   D_OR_S()   precision prefix for the BLAS and LAPACK Library routines.
!              Use, e.g.,  D_OR_S()axpy(...)  as the precision independent
!              notation for daxpy(...) and saxpy(...).

!   DPREC      defined when the working precision is double;
!              undefined when the working precision is single.

!   VD_OR_VS() precision prefix for the Intel Vector Math Library routines.
!              Use, e.g.,  VD_OR_VS()exp(...)  as the precision independent
!              notation for vdexp(...) and vsexp(...).

!   WIDE_REAL  defined when a platform's single precision is wider than the
!              IEEE 754 single precision (32-bit).  This is the case on Cray
!              machines.  Note that on these machines _REAL_ is defined to
!              double precision, and the compiler flag -dp is used to disable
!              double precision.  The result is consistent 64-bit working
!              precision.

! History:
! $Id: dprec.h,v 9.0 2006/04/03 23:35:47 case Exp $

! Code Description:
!   Languages:          C Preprocessor and Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code":
!     http://nsipp.gsfc.nasa.gov/infra/eurorules.html

! References:
!   IEEE 754: Standard for Binary Floating-Point Arithmetic.
!     http://grouper.ieee.org/groups/754/

!   The unused Fortran 90 module numerics.f90.


#ifndef _REAL_
#  define _REAL_ double precision
#  define AMBER_MPI_REAL MPI_DOUBLE_PRECISION
#  define DPREC 1
#else
#  undef  _REAL_
#  define _REAL_ real
#  undef  AMBER_MPI_REAL
#  define AMBER_MPI_REAL MPI_REAL
#  undef  DPREC
#endif

#ifdef CRAY_PVP
#  undef  AMBER_MPI_REAL
#  define AMBER_MPI_REAL MPI_REAL8
#endif

#if defined WIDE_REAL || ! defined DPREC
#  define D_OR_S() s
#  define VD_OR_VS() s
#else
#  define D_OR_S() d
#  define VD_OR_VS() d
#endif

! DPREC_H
#endif

