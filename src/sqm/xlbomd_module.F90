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

contains
	
	!Initialize and Allocate Variables
	subroutine init_xlbomd(xlbomd_struct,mat_size)

		implicit none
		
		type(xlbomd_structure), intent(inout) :: xlbomd_struct
		integer		::	n,mat_size
		!Allocate and assign xlbomd_struct%coefficients
		allocate(xlbomd_struct%coef(xlbomd_struct%K+1))
		allocate(xlbomd_struct%phi(mat_size,xlbomd_struct%K+1))
		allocate(xlbomd_struct%phi_point(xlbomd_struct%K+2)) !xlbomd_struct%K+2 is here for rotating the pointers
		xlbomd_struct%phi=0.d0 !Clear xlbomd_struct%phi
		xlbomd_struct%Kpassable=xlbomd_struct%K
		select case (xlbomd_struct%K)	
			case(3)
				xlbomd_struct%dt2w2=1.69 
				xlbomd_struct%xlalpha=150d-3
				xlbomd_struct%coef=(/-2,3,0,-1/)
			case(4)
				xlbomd_struct%dt2w2=1.75
				xlbomd_struct%xlalpha=57d-3
				xlbomd_struct%coef=(/-3,6,-2,-2,1/)
			case(5)
				xlbomd_struct%dt2w2=1.82
				xlbomd_struct%xlalpha=18d-3
				xlbomd_struct%coef=(/-6 ,14 ,-8 ,-3 ,4,-1/) 
			case(6)
				xlbomd_struct%dt2w2=1.84
				xlbomd_struct%xlalpha=5.5d-3
				xlbomd_struct%coef=(/-14,36,-27,-2,12,-6,1/)
			case(7)
				xlbomd_struct%dt2w2=1.86
				xlbomd_struct%xlalpha=1.6d-3
				xlbomd_struct%coef=(/-36,99,-88,11,32,-25,8,-1/)
			case(8)
				xlbomd_struct%dt2w2=1.88
				xlbomd_struct%xlalpha=0.44d-3
				xlbomd_struct%coef=(/-99,286,-286,78,78,-90,42,-10,1/)
			case(9)
				xlbomd_struct%dt2w2=1.89
				xlbomd_struct%xlalpha=0.12d-3
				xlbomd_struct%coef=(/-286,858,-936,364,168,-300,184 ,-63 ,12,-1/)
		end select
		!Initial pointer assignments
		do n=1,xlbomd_struct%K+1
			xlbomd_struct%phi_point(n)%guess=>xlbomd_struct%phi(1:mat_size,n)
		enddo
	end subroutine init_xlbomd
	
	!Predict density matrix from previous guesses. Could also be used for wavefunctions. 
	subroutine predictdens_xlbomd(xlbomd_struct,num_calls,P)

		implicit none

		integer		::	n,num_calls
		_REAL_		::	P(:)
                type(xlbomd_structure),intent(inout) :: xlbomd_struct

		if (num_calls<xlbomd_struct%K+2) then !Initial burnin
                        !Add the solved density matrix as the initial guess for 'burn in'
                        write(6,*)'Full convergence for XL-BOMD Burn in:',num_calls
                        if (num_calls>0) then
                        	xlbomd_struct%phi_point(xlbomd_struct%K+2-num_calls)%guess=P
			endif
                endif   
		if (num_calls>xlbomd_struct%K+1) then !this is after the initial 'burn in' and on the first step
			!Calculate new initial guess
			write(6,*)'XL-BOMD Predicted Initial Guess'
			P=2.0*xlbomd_struct%phi_point(1)%guess-xlbomd_struct%phi_point(2)%guess+xlbomd_struct%dt2w2*(P-xlbomd_struct%phi_point(1)%guess)
			do n=1,xlbomd_struct%K+1
				P=P+xlbomd_struct%xlalpha*xlbomd_struct%coef(n)*xlbomd_struct%phi_point(n)%guess
			enddo
			
			xlbomd_struct%phi_point(xlbomd_struct%K+2)%guess=>xlbomd_struct%phi_point(xlbomd_struct%K+1)%guess !extra pointer for rotating
			!Increment the stack of initial guesses
			do n=xlbomd_struct%K,1,-1 !Increment stack
				xlbomd_struct%phi_point(n+1)%guess=>xlbomd_struct%phi_point(n)%guess
			enddo
                        xlbomd_struct%phi_point(1)%guess=>xlbomd_struct%phi_point(xlbomd_struct%K+2)%guess!points to last member of the previous stack for replacement
			
			xlbomd_struct%phi_point(1)%guess=P !Push new guess onto stack
			
		endif	
	end subroutine predictdens_xlbomd


end module xlbomd_module
