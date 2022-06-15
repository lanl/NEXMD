#include "dprec.fh"
#include "assert.fh"

!*******************************************************************************
!*This file contains auxiliary subroutines for calculating nuclear normal modes*
!*******************************************************************************


!Subroutine for the calculation of the forces
subroutine hessian_calc(sim,deltaX)
    use communism !for numerical derivatives
    use naesmd_constants
    use constants, only: TWOPI 

    implicit none

    type(simulation_t),target :: sim !communism module
    type(simulation_t),pointer :: simpoint

!Pased
    _REAL_, intent(in) :: deltaX
!Internal
    integer i,j,k
    _REAL_, allocatable, dimension(:,:) :: fxPlus, fxMinus !Forces for calculation the Hessian
    _REAL_, allocatable, dimension(:,:) :: fyPlus, fyMinus
    _REAL_, allocatable, dimension(:,:) :: fzPlus, fzMinus
    _REAL_, allocatable, dimension(:,:) :: hessian
    _REAL_, allocatable, dimension(:,:) :: a !Normal modes
    _REAL_, allocatable, dimension(:) :: work, iwork !For the diagonlization of the Hessian
    integer n, lda, lwork, liwork, info !For the diagonlization of the Hessian
    _REAL_, allocatable, dimension(:) :: w, freq2 !Frequencies
    
    simpoint => sim

!Allocating internal variables
    n = 3*simpoint%naesmd%natom
    lda = 3*simpoint%naesmd%natom
    lwork = 2*(3*simpoint%naesmd%natom)**2+6*(3*simpoint%naesmd%natom)+1
    liwork = 5*(3*simpoint%naesmd%natom)+3
    allocate(fxPlus(simpoint%naesmd%natom,3*simpoint%naesmd%natom), fxMinus(simpoint%naesmd%natom,3*simpoint%naesmd%natom))
    allocate(fyPlus(simpoint%naesmd%natom,3*simpoint%naesmd%natom), fyMinus(simpoint%naesmd%natom,3*simpoint%naesmd%natom))
    allocate(fzPlus(simpoint%naesmd%natom,3*simpoint%naesmd%natom), fzMinus(simpoint%naesmd%natom,3*simpoint%naesmd%natom))
    allocate(hessian(3*simpoint%naesmd%natom,3*simpoint%naesmd%natom))
    allocate(w(3*simpoint%naesmd%natom))
    allocate(freq2(3*simpoint%naesmd%natom))
    allocate(a(lda,n))
    allocate(work(lwork),iwork(liwork))

    write(6,*) "Calculating nuclear normal modes, deltaX = ", deltaX

!Calculating fPlus
    do i=1,simpoint%naesmd%natom
        !fxPlus
        simpoint%naesmd%rx(i) = simpoint%naesmd%rx(i) + deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest 
        call deriv(simpoint,simpoint%naesmd%ihop)
        fxPlus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%rx(i) = simpoint%naesmd%rx(i) - deltaX/convl

        !fyPlus
        simpoint%naesmd%ry(i) = simpoint%naesmd%ry(i) + deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call deriv(simpoint,simpoint%naesmd%ihop)
        fyPlus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%ry(i) = simpoint%naesmd%ry(i) - deltaX/convl

        !fzPlus
        simpoint%naesmd%rz(i) = simpoint%naesmd%rz(i) + deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call deriv(simpoint,simpoint%naesmd%ihop)
        fzPlus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%rz(i) = simpoint%naesmd%rz(i) - deltaX/convl

        !fxMinus
        simpoint%naesmd%rx(i) = simpoint%naesmd%rx(i) - deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call deriv(simpoint,simpoint%naesmd%ihop)
        fxMinus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%rx(i) = simpoint%naesmd%rx(i) + deltaX/convl

        !fyMinus
        simpoint%naesmd%ry(i) = simpoint%naesmd%ry(i) - deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call deriv(simpoint,simpoint%naesmd%ihop)
        fyMinus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%ry(i) = simpoint%naesmd%ry(i) + deltaX/convl

        !fzMinus
        simpoint%naesmd%rz(i) = simpoint%naesmd%rz(i) - deltaX/convl
        call naesmd2qmmm_r(simpoint)
        if(simpoint%excn>0) then
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%excn+1)
        else
            call do_sqm_davidson_update(simpoint,0,simpoint%naesmd%cmdqt,simpoint%naesmd%vmdqt,simpoint%naesmd%vgs, &
                statelimit=simpoint%qmmm%state_of_interest)
        end if
        sim%naesmd%ihop=sim%qmmm%state_of_interest
        call deriv(simpoint,simpoint%naesmd%ihop)
        fzMinus(i,:) = simpoint%deriv_forces*convl/feVmdqt !Converted to a.u.
        simpoint%naesmd%rz(i) = simpoint%naesmd%rz(i) + deltaX/convl
    enddo
    !Calculating the hessian
    do i=1,simpoint%naesmd%natom
        hessian(3*i-2,:) = -(fxPlus(i,:) - fxMinus(i,:))/(2.0d0*deltaX/convl)
        hessian(3*i-1,:) = -(fyPlus(i,:) - fyMinus(i,:))/(2.0d0*deltaX/convl)
        hessian(3*i,:) = -(fzPlus(i,:) - fzMinus(i,:))/(2.0d0*deltaX/convl)
    enddo
    !Symmetrizing the hessian
    do i=1,3*simpoint%naesmd%natom
        do j=i,3*simpoint%naesmd%natom
            hessian(i,j)=0.5d0*(hessian(i,j)+hessian(j,i))
        enddo
    enddo 
    !Weighting by mass
    do i=1,simpoint%naesmd%natom
        do j=i,simpoint%naesmd%natom
            hessian(3*i-2,3*j-2)=hessian(3*i-2,3*j-2)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i-2,3*j-1)=hessian(3*i-2,3*j-1)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i-2,3*j)=hessian(3*i-2,3*j)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i-1,3*j-2)=hessian(3*i-1,3*j-2)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i-1,3*j-1)=hessian(3*i-1,3*j-1)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i-1,3*j)=hessian(3*i-1,3*j)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i,3*j-2)=hessian(3*i,3*j-2)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i,3*j-1)=hessian(3*i,3*j-1)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
            hessian(3*i,3*j)=hessian(3*i,3*j)/dsqrt(simpoint%naesmd%massmdqt(i))/dsqrt(simpoint%naesmd%massmdqt(j))
        enddo
    enddo
    !Finishing symmetrizig of the Hessian
    do i=1,3*simpoint%naesmd%natom
        do j=i,3*simpoint%naesmd%natom
            hessian(j,i)=hessian(i,j)
        enddo
    enddo
    do i=1,3*simpoint%naesmd%natom
        do j=1,3*simpoint%naesmd%natom
            a(i,j)=hessian(i,j)
        enddo
    enddo

    call dsyevd('V','U',n,a,lda,w,work,lwork,iwork,liwork,info)

!    print *, 'dsyevd for the Hessian diagonalization finished, info: ', info
       

!Writing data
    open(34, file='nma_modes.out') !Normal modes eigenvalues
    do i=1,3*simpoint%naesmd%natom
        write(34,198) (a(i,j),j=1,3*simpoint%naesmd%natom)
    enddo
    close(34)
    open(35, file='nma_freq.out') !Normal modes frequencies
    do i=1,3*simpoint%naesmd%natom
        freq2(i)=dsqrt(dabs(w(i)))
        if (w(i).lt.0.d0) then
           freq2(i)=-freq2(i)
        endif 
        write(35,*) i,freq2(i)
    enddo
    close(35)
    open(10,file='hessian.out') !Hessian
    do i=1,3*simpoint%naesmd%natom
        write(10,198) (hessian(i,j),j=1,3*simpoint%naesmd%natom)
    enddo
    close(10)

198   format(1000(e18.10,1x))

end subroutine hessian_calc
