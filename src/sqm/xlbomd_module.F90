#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

! ----------------------------------------------------------------------
! Written by Josiah A. Bjorgaard at LANL March 2014
! Holds parameters and allocates variables for XL-BOMD Calculations
! Phi:auxiliarly variable which contributes to the initial scf guess
! ----------------------------------------------------------------------
module xlbomd_module

	implicit none
	
	!Following variables are for pushing a stack containing all phi
	type pointers_to_phi
		_REAL_, pointer	::	guess(:)
	end type pointers_to_phi
	
	type xlbomd_structure
		!Variables initialized from qmmm namelist
		integer			::	xlbomd_flag,K,Kpassable
		_REAL_			::	dt2w2,xlalpha
		
		!Variables not from qmmm namelist
		integer 		::	mat_size,num_calls
		_REAL_, pointer	::	coef(:),phi(:,:)
		!_REAL_, allocatable	::	P(:)
		type(pointers_to_phi), dimension(:),allocatable	::	phi_point
	end type xlbomd_structure



end module xlbomd_module
