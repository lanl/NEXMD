#include "dprec.fh"
#include "assert.fh"

!
!********************************************************************      
!********************************************************************      
!   
subroutine fcn(x,yg_new,ygprime_new,naesmd_struct)
	use naesmd_module
	implicit none
      _REAL_, intent(in) :: x                                  !indep!
      _REAL_, dimension(:), intent(in) :: yg_new                   !dep!
      _REAL_, dimension(:), intent(in) :: ygprime_new                   !dep!
       type(naesmd_structure), intent(inout) :: naesmd_struct 

	    call interpolate(size(yg_new,1),x,naesmd_struct)
	    call vqcalc(size(yg_new,1),yg_new,ygprime_new,naesmd_struct)
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
    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct 
    integer n,k,j
    _REAL_ x
 
    do k=1,n/3
        do j=1,n/3
            naesmd_struct%cadiab(k,j)=naesmd_struct%cadiabmiddleold(k,j) &
                +naesmd_struct%bcoeffcadiab(k,j)*(x-naesmd_struct%tini0)
        end do
    end do

    do k=1,n/3
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

subroutine vqcalc(n,yg_new,yprime_new,naesmd_struct)
    use naesmd_module
    use md_module
    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct 
    integer n,k,j
    _REAL_ yg_new(n),yprime_new(n)

    do k=1,n/3
        yprime_new(k)=0.d0
        yprime_new(k+n/3)=0.d0
        yprime_new(k+2*n/3)=-1.0d0*naesmd_struct%vmdqt(k)

        do j=1,n/3
            yprime_new(k)=yprime_new(k)-1.0d0*(yg_new(j) &
               *dcos(yg_new(j+2*n/3)-yg_new(k+2*n/3))-yg_new(j+n/3) &
               *dsin(yg_new(j+2*n/3)-yg_new(k+2*n/3)))*naesmd_struct%cadiab(k,j)
        end do

        do j=1,n/3
            yprime_new(k+n/3)=yprime_new(k+n/3)-1.0d0*(yg_new(j+n/3) &
                *dcos(yg_new(j+2*n/3)-yg_new(k+2*n/3))+yg_new(j) &
                *dsin(yg_new(j+2*n/3)-yg_new(k+2*n/3)))*naesmd_struct%cadiab(k,j)

            naesmd_struct%vnqcorrhop(k,j)=2.0d0*yg_new(k) &
                *(-1.0d0)*(yg_new(j)*dcos(yg_new(j+2*n/3)-yg_new(k+2*n/3))-yg_new(j+n/3) &
                *dsin(yg_new(j+2*n/3)-yg_new(k+2*n/3)))*naesmd_struct%cadiab(k,j) &
                +2.0d0*yg_new(k+n/3)*(-1.0d0)*(yg_new(j+n/3)*dcos(yg_new(j+2*n/3)-yg_new(k+2*n/3))+yg_new(j) &
                *dsin(yg_new(j+2*n/3)-yg_new(k+2*n/3)))*naesmd_struct%cadiab(k,j)
        end do
    end do
    return
end subroutine
