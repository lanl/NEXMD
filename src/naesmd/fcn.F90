#include "dprec.fh"
#include "assert.fh"

!
!********************************************************************      
!********************************************************************      
!   
subroutine fcn(x,yg,ygprime,naesmd_struct)
	use naesmd_module
	implicit none
      _REAL_, intent(in) :: x                                  !indep!
      _REAL_, dimension(:), intent(in) :: yg                   !dep!
      _REAL_, dimension(:), intent(in) :: ygprime                   !dep!
       type(naesmd_structure), intent(inout) :: naesmd_struct 

	    call interpolate(size(yg,1),x,naesmd_struct)
	    call vqcalc(size(yg,1),yg,ygprime,naesmd_struct)
       return
end subroutine fcn
 
!********************************************************************
!
!  Subroutine to interpolate the values of the terms <phi_k | d phi_j/ d t >
!  and the adiabatic energies
!
!********************************************************************
!
subroutine interpolate(n,x,naesmd_struct)
    use naesmd_module
    use md_module
    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct 
    integer n,k,j
    _REAL_ x
 
    do k=1,n/2
        do j=1,n/2
            naesmd_struct%cadiab(k,j)=naesmd_struct%cadiabmiddleold(k,j) &
                +naesmd_struct%bcoeffcadiab(k,j)*(x-naesmd_struct%tini0)
        end do
    end do

    do k=1,n/2
        naesmd_struct%vmdqt(k)=naesmd_struct%vmdqtmiddleold(k)+naesmd_struct%bcoeffvmdqt(k)*(x-naesmd_struct%tini0)
    end do
        
    return
end subroutine
!
!********************************************************************
!
!  Subroutine to calculate de d(yg)/dt and dqq/dt
!  That is, the derivatives of the modulus and phase of the
!  electronic coefficients
!
!********************************************************************
!

subroutine vqcalc(n,yg,yprime,naesmd_struct)
    use naesmd_module
    use md_module
    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct 
    integer n,k,j
    _REAL_ yg(n),yprime(n)

    do k=1,n/2
        yprime(k)=0.d0
        yprime(k+n/2)=-1.0d0*naesmd_struct%vmdqt(k)

        do j=1,n/2
            naesmd_struct%vnqcorrhop(k,j)=-1.0d0*yg(j)* &
                dcos(yg(j+n/2)-yg(k+n/2))*naesmd_struct%cadiab(k,j)

            yprime(k)=yprime(k)+naesmd_struct%vnqcorrhop(k,j)
            naesmd_struct%vnqcorrhop(k,j)=naesmd_struct%vnqcorrhop(k,j)*2.0d0*yg(k)
        end do

        do j=1,n/2
            yprime(k+n/2)=yprime(k+n/2)-yg(j)/(1.0d-6+yg(k))* &
                dsin(yg(j+n/2)-yg(k+n/2))*naesmd_struct%cadiab(k,j)
        end do
    end do
    return
end subroutine
