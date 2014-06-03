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
	
	!Variables initialized from qmmm namelist
	integer			::	xlbomd_flag,K,Kpassable
	_REAL_			::	dt2w2,xlalpha
	
	!Variables not from qmmm namelist
	integer 		::	mat_size,num_calls
	_REAL_, allocatable,target	::	coef(:),phi(:,:)
	!_REAL_, allocatable	::	P(:)

	!Following variables are for pushing a stack containing all phi
	type pointers_to_phi
		_REAL_, pointer	::	guess(:)
	end type pointers_to_phi
	type(pointers_to_phi), dimension(:),allocatable	::	phi_point

contains
	
	!Initialize and Allocate Variables
	subroutine init_xlbomd(mat_size)

		implicit none

		integer		::	n,mat_size
		!Allocate and assign coefficients
		allocate(coef(K+1))
		allocate(phi(mat_size,K+1))
		allocate(phi_point(K+2)) !K+2 is here for rotating the pointers
		phi=0.d0 !Clear phi
		Kpassable=K
		select case (K)	
			case(3)
				dt2w2=1.69 
				xlalpha=150d-3
				coef=(/-2,3,0,-1/)
			case(4)
				dt2w2=1.75
				xlalpha=57d-3
				coef=(/-3,6,-2,-2,1/)
			case(5)
				dt2w2=1.82
				xlalpha=18d-3
				coef=(/-6 ,14 ,-8 ,-3 ,4,-1/) 
			case(6)
				dt2w2=1.84
				xlalpha=5.5d-3
				coef=(/-14,36,-27,-2,12,-6,1/)
			case(7)
				dt2w2=1.86
				xlalpha=1.6d-3
				coef=(/-36,99,-88,11,32,-25,8,-1/)
			case(8)
				dt2w2=1.88
				xlalpha=0.44d-3
				coef=(/-99,286,-286,78,78,-90,42,-10,1/)
			case(9)
				dt2w2=1.89
				xlalpha=0.12d-3
				coef=(/-286,858,-936,364,168,-300,184 ,-63 ,12,-1/)
		end select
		!Initial pointer assignments
		do n=1,K+1
			phi_point(n)%guess=>phi(1:mat_size,n)
		enddo
	end subroutine init_xlbomd
	
	!Predict density matrix from previous guesses. Could also be used for wavefunctions. 
	subroutine predictdens_xlbomd(num_calls,P)

		implicit none

		integer		::	n,num_calls
		_REAL_		::	P(:)

		if (num_calls<K+2) then !Initial burnin
                        !Add the solved density matrix as the initial guess for 'burn in'
                        write(6,*)'Full convergence for XL-BOMD Burn in:',num_calls
                        if (num_calls>0) then
                        	phi_point(K+2-num_calls)%guess=P
			endif
                        !write(6,*)'Density matrix from burnin', num_calls
                        !write(6,*)P
			!write(6,*)'***********************************************************************'
                endif   
		if (num_calls>K+1) then !this is after the initial 'burn in' and on the first step
			!Test
			!write(6,*)'Input density matrix D'
			!write(6,*)P
			!write(6,*)'----------------------------------------------------------------------'
			!do n=1,K+1
			!	write(6,*)'Initial Density Matrix',n
			!	write(6,*)phi_point(n)%guess
			!	write(6,*)'------------------------------------------------------------------'
			!enddo
			!Calculate new initial guess
			write(6,*)'XL-BOMD Predicted Initial Guess'
			P=2.0*phi_point(1)%guess-phi_point(2)%guess+dt2w2*(P-phi_point(1)%guess)
			do n=1,K+1
				P=P+xlalpha*coef(n)*phi_point(n)%guess
			enddo
			
			phi_point(K+2)%guess=>phi_point(K+1)%guess !extra pointer for rotating
			!Increment the stack of initial guesses
			do n=K,1,-1 !Increment stack
				phi_point(n+1)%guess=>phi_point(n)%guess
			enddo
                        phi_point(1)%guess=>phi_point(K+2)%guess!points to last member of the previous stack for replacement
			
			!Test 0
			!write(6,*)'Predicted Density Matrix'
			!write(6,*)P
			!write(6,*)'--------------------------------------------------------------------'
			!Test 1
                        !do n=1,K+1
                        !        write(6,*)'Final Density Matrix 1',n
                        !        write(6,*)phi_point(n)%guess
                        !        write(6,*)'--------------------------------------------------------------------'
                        !enddo

			phi_point(1)%guess=P !Push new guess onto stack
			
			!Test 2
			!do n=1,K+1
			!	write(6,*)'Final Density Matrix 2',n
			!	write(6,*)phi_point(n)%guess
			!	write(6,*)'--------------------------------------------------------------------'
			!enddo
			!write(6,*)'P'
			!write(6,*)P
			!write(6,*)'---------------------------------------------------------------------'
		endif	
	end subroutine predictdens_xlbomd

end module xlbomd_module
