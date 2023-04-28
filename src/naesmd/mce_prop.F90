#include "dprec.fh"
#include "assert.fh"

!module mce_prop
!    use naesmd_constants
!    use communism
!    use cadiab_module
!    use additional_subroutines

!contains

subroutine update_nuclear(nuclear,sims,Nsim,w,wr)
    use naesmd_constants
    use communism
    implicit none
    type(mce) :: nuclear
    type(sim_pointer_array), dimension(Nsim), intent(in) :: sims
    integer, intent(in) :: Nsim
    integer, intent(in) :: w  !for writing mce data
    integer, intent(in) :: wr !for writing restart (last files)

    integer :: n,m,i,j, k, info, nexc
    _REAL_, dimension(3*sims(1)%sim%naesmd%natom) :: x1, x2 !vectors with positions
    _REAL_, dimension(3*sims(1)%sim%naesmd%natom) :: p1, p2 !vectors with momenta
    _REAL_, dimension(3*sims(1)%sim%naesmd%natom) :: mass !masses for each dimension
    double complex, dimension(Nsim,Nsim) :: S !to store complete overlap
    double complex, dimension(Nsim,Nsim) :: sE_all !to store the complete electronic overlap
    double complex, dimension(Nsim,Nsim) :: sek !kinetic matrix elements
    double complex, dimension(Nsim,Nsim) :: Sinv !To invert the complete overlap matrix
    double complex, dimension(Nsim) :: tmp4 !Temporary array for matrix multiplication
    double complex, dimension(sims(1)%sim%excN) :: a1, a2, tmp1, tmp3 !electronic amplitudes 
    double complex ek,alpha,beta !to store excpectation value of the nuclear kinetic operator T
    double complex S_dot !to store the time derivative of the nuclear overlap
    double complex S_dotE !to store the time derivative of the electronic overlaps
    double complex, dimension(3*sims(1)%sim%naesmd%natom) :: Sdx, Sdp
    _REAL_, dimension(3*sims(1)%sim%naesmd%natom,sims(1)%sim%excN,sims(1)%sim%excN) :: grad12
    double complex, dimension(3*sims(1)%sim%naesmd%natom) :: gradtr, gradd
    double complex, dimension(sims(1)%sim%excN,sims(1)%sim%excN) :: hen2,tmp2
    double complex, dimension(sims(1)%sim%excN,sims(1)%sim%excN) :: he12
    double complex, dimension(Nsim,Nsim) :: H
    double complex, dimension(Nsim,Nsim,sims(1)%sim%excN) :: AmpE_traj !To calculate populations
    _REAL_, dimension(sims(1)%sim%excN) :: pop
    logical :: cloned

    AmpE_traj = 0.0d0
    sE_all = 0.0d0
    sek =0.0d0
    ek = 0.0d0
    hen2 = 0.0d0
    he12 = 0.0d0
    S_dotE = 0.0d0
    S_dot = 0.0d0
    H = 0.0d0
    AmpE_traj = 0.0d0
    alpha = cmplx(1.0d0,0.0d0)
    beta = cmplx(0.0d0,0.0d0)
    tmp1 = 0.0d0
    tmp2 = 0.0d0
    tmp3 = 0.0d0
    tmp4 = 0.0d0
    nexc = sims(1)%sim%excN



    cloned = .false.
    do i=1,Nsim
        if(nuclear%cloned(i)) cloned =.true.
    enddo
    do n=1,Nsim
        do i=1,sims(n)%sim%naesmd%natom
            x1(3*i-2)=sims(n)%sim%naesmd%rx(i)
            x1(3*i-1)=sims(n)%sim%naesmd%ry(i)
            x1(3*i)=sims(n)%sim%naesmd%rz(i)
            mass(3*i-2)=sims(n)%sim%naesmd%massmdqt(i)
            mass(3*i-1)=sims(n)%sim%naesmd%massmdqt(i)
            mass(3*i)=sims(n)%sim%naesmd%massmdqt(i)
            p1(3*i-2)=mass(3*i-2)*sims(n)%sim%naesmd%vx(i)
            p1(3*i-1)=mass(3*i-1)*sims(n)%sim%naesmd%vy(i)
            p1(3*i)=mass(3*i)*sims(n)%sim%naesmd%vz(i)
        enddo
        a1=(sims(n)%sim%naesmd%yg_new(1:sims(n)%sim%excN)+ &
           (0.0d0,1.0d0)*sims(n)%sim%naesmd%yg_new(sims(n)%sim%excN+1:2*sims(n)%sim%excN))* &
           (cos(sims(n)%sim%naesmd%yg_new(2*sims(n)%sim%excN+1:3*sims(n)%sim%excN))+ &
           (0.0d0,1.0d0)*(sin(sims(n)%sim%naesmd%yg_new(2*sims(n)%sim%excN+1:3*sims(n)%sim%excN))))
        do m=1,Nsim
            a2=(sims(m)%sim%naesmd%yg_new(1:sims(m)%sim%excN)+ &
               (0.0d0,1.0d0)*sims(m)%sim%naesmd%yg_new(sims(m)%sim%excN+1:2*sims(m)%sim%excN))* &
               (cos(sims(m)%sim%naesmd%yg_new(2*sims(m)%sim%excN+1:3*sims(m)%sim%excN))+ &
               (0.0d0,1.0d0)*(sin(sims(m)%sim%naesmd%yg_new(2*sims(m)%sim%excN+1:3*sims(m)%sim%excN))))
            do i=1,sims(n)%sim%naesmd%natom
                x2(3*i-2)=sims(m)%sim%naesmd%rx(i)
                x2(3*i-1)=sims(m)%sim%naesmd%ry(i)
                x2(3*i)=sims(m)%sim%naesmd%rz(i)
                p2(3*i-2)=mass(3*i-2)*sims(m)%sim%naesmd%vx(i)
                p2(3*i-1)=mass(3*i-1)*sims(m)%sim%naesmd%vy(i)
                p2(3*i)=mass(3*i)*sims(m)%sim%naesmd%vz(i)
            enddo
            call overlap_nuclear(x1,p1,sims(n)%sim%naesmd%pha, &           !coordinates, momenta and gamma for n
                                 x2,p2,sims(m)%sim%naesmd%pha, &           !coordinates, momenta and gamma for m
                                 sims(m)%sim%deriv_forces*convl/feVmdqt, & !forces for m
                                 3*sims(1)%sim%naesmd%natom,  &            !number of dimensions
                                 S(n,m), S_dot, ek, Sdx, Sdp,  &           !output
                                 mass, sims(1)%sim%naesmd%w)               !masses and widths
            do i=1, nexc
                do j=1, nexc
                    tmp2(i,j) = cmplx(nuclear%sE(n,m,i,j),0.0d0)
                enddo
            enddo
            call zgemm('N', 'N',  nexc, 1, nexc, alpha, tmp2, nexc, a2, nexc, beta, tmp1, nexc)

            sE_all(n,m)=dot_product(a1,tmp1) !Electronic overlap (including electronic amplitudes)
            sek(n,m)=ek*sE_all(n,m) !MCE kinetic energy matrix
            do i=1,sims(n)%sim%excN
                do j=1,sims(m)%sim%excN
                    grad12(:,i,j)=(sims(n)%sim%deriv_forces_state(i,:)+sims(m)%sim%deriv_forces_state(j,:))*convl/feVmdqt
                enddo
            enddo
            do i=1,3*sims(1)%sim%naesmd%natom
                do j=1, nexc
                    do k=1, nexc
                        tmp2(j,k) = cmplx(nuclear%sE(n,m,j,k)*grad12(i,j,k),0.0d0)
                    enddo
                enddo

                call zgemm('N', 'N', nexc, 1, nexc, alpha, tmp2, nexc, a2, nexc, beta, tmp1, nexc)

                gradtr(i)=dot_product(a1,tmp1)
            enddo

            do i=1,sims(n)%sim%excN
                do j=1,sims(n)%sim%excN
                    grad12(:,i,j)=(sims(n)%sim%deriv_forces_state(i,:)-sims(m)%sim%deriv_forces_state(j,:))*convl/feVmdqt
                enddo
            enddo
            do i=1,3*sims(1)%sim%naesmd%natom
                do j=1, nexc
                    do k=1, nexc
                        tmp2(j,k) = cmplx(nuclear%sE(n,m,j,k)*grad12(i,j,k),0.0d0)
                    enddo
                enddo
                call zgemm('N', 'N', nexc, 1, nexc, alpha, tmp2, nexc, a2, nexc, beta, tmp1, nexc)

                gradd(i)=dot_product(a1,tmp1)
            enddo
            do i=1,sims(n)%sim%excN
                do j=1,sims(m)%sim%excN
                    hen2(i,j)=sims(n)%sim%naesmd%vmdqt(i)+sims(m)%sim%naesmd%vmdqt(j)
                enddo
            enddo
            he12=nuclear%sE(n,m,:,:)*hen2
            call zgemm('N', 'N', nexc, 1, nexc, alpha ,he12, nexc, a2, nexc, beta, tmp1, nexc)
            H(n,m)=ek*sE_all(n,m)+0.5d0*dot_product(a1,tmp1)*S(n,m)+ &
                   0.5d0*(0.0d0,0.25d0)*sum((p2-p1)*gradtr/sims(n)%sim%naesmd%w)*S(n,m)- &
                   0.25d0*sum((x1-x2)*gradd)*S(n,m)
            do j=1, nexc
                do k=1, nexc
                        tmp2(j,k) = cmplx(nuclear%sE(n,m,j,k),0.0d0)
                enddo
                tmp3(j) = a2(j)*cmplx(sims(m)%sim%naesmd%vmdqt(j),0.0d0)
            enddo

            call zgemm('N', 'N', nexc, 1, nexc, alpha, tmp2, nexc, tmp3, nexc, beta, tmp1, nexc)

            S_dotE=-(0.d0,1.d0)*dot_product(a1,tmp1)
            S_dot = S_dot*sE_all(n,m)+S_dotE*S(n,m)

            H(n,m)=H(n,m)-(0.d0,1.d0)*S_dot
            do i=1, nexc
                do j=1, nexc
                    tmp2(i,j) = cmplx(nuclear%sE(n,m,i,j),0.0d0)
                enddo
            enddo

            call zgemm('N', 'N', nexc, 1, nexc, alpha, tmp2, nexc, a2, nexc, beta, tmp1, nexc)

            AmpE_traj(n,m,:)=conjg(a1)*tmp1
        enddo!m
    enddo!n

!debugging
!            write(nuclear%outfile_0,*) sims(1)%sim%naesmd%tfemto, (sims(i)%sim%aimc%nclones, n=1, Nsim)
!debugging
    
    Sinv(1:Nsim,1:Nsim)=S(1:Nsim,1:Nsim)*sE_all(1:Nsim,1:Nsim)!To invert the overlap matrix

    call zpotrf('U',Nsim,Sinv,Nsim,info)
    if(info.ne.0) then
        print*,'Error inverting overlap matrix for MCE propagation: 1', info
        stop
    endif
    call zpotri('U',Nsim,Sinv,Nsim,info)
    if(info.ne.0) then
        print*,'Error inverting overlap matrix for MCE propagation: 2', info
        stop
    endif
    do i=2,Nsim
        do j=1,i
            Sinv(i,j)=conjg(Sinv(j,i))
        enddo
    enddo

    do i=1,Nsim
        do j=1,Nsim
            nuclear%Heff_old(i,j)=nuclear%Heff(i,j)
        enddo
    enddo
    call zgemm('N', 'N',Nsim, Nsim, Nsim, alpha, Sinv, Nsim, H(1:Nsim,1:Nsim), Nsim, beta, nuclear%Heff(1:Nsim,1:Nsim), Nsim)

    if(cloned) then
        nuclear%Heff_old(1:Nsim,1:Nsim)=nuclear%Heff(1:Nsim,1:Nsim)
    endif
    nuclear%cloned=.false.

    call propagate_c(nuclear,sims,Nsim,cloned)

    do i=1,sims(1)%sim%excN
        call zgemm('N', 'N', Nsim, 1, Nsim, alpha, S(1:Nsim,1:Nsim)*AmpE_traj(1:Nsim,1:Nsim,i), Nsim, nuclear%D(1:Nsim), &
            Nsim, beta, tmp4, Nsim)

        pop(i)=real(dot_product(nuclear%D(1:Nsim),tmp4))
    enddo
    nuclear%pop=pop
    nuclear%D=nuclear%D/sqrt(sum(pop)) !Renormalization (sum(pop) should be 1)

!Writing MCE output
    if(w.eq.1) then
        write(nuclear%outfile_1,999) Nsim, sims(1)%sim%naesmd%tfemto, (pop(i),i=1,sims(1)%sim%excN), sum(pop)
        write(nuclear%outfile_2,999) Nsim, sims(1)%sim%naesmd%tfemto, nuclear%D(1:Nsim)
        call flush(nuclear%outfile_1)
        call flush(nuclear%outfile_2)
        do n=1,Nsim
            do m=n+1,Nsim
                do i=1,sims(1)%sim%excN
                    do j=1,sims(1)%sim%excN
                        write(nuclear%outfile_3,998) Nsim, sims(1)%sim%naesmd%tfemto, n, m, i, j, nuclear%sE(n,m,i,j)
                    enddo
                enddo
            enddo
        write(sims(n)%sim%outfile_29,997) sims(n)%sim%naesmd%tfemto, sims(n)%sim%naesmd%pha 
        enddo
        call flush(nuclear%outfile_3)
        call flush(sims(n)%sim%outfile_29)
        if(wr.eq.1) then !Writing last files
            open(nuclear%outfile_4,file='Heff_last.out')
                write(nuclear%outfile_4,*) sims(1)%sim%naesmd%tfemto
                do i = 1, Nsim
                    write(nuclear%outfile_4,997) nuclear%Heff(i,1:Nsim)
                enddo
            close(nuclear%outfile_4)
            open(nuclear%outfile_5,file='Pha_last.out')
                write(nuclear%outfile_5,999) Nsim, sims(1)%sim%naesmd%tfemto,(sims(i)%sim%naesmd%pha,i=1,Nsim)
            close(nuclear%outfile_5)
        endif
    endif

999     FORMAT(I3,1000(1X,F18.10))
997     FORMAT(1000(1X,F18.10))
998     FORMAT(I3,1(1X,F18.10),4I3,1000(1X,F18.10))

end subroutine update_nuclear


subroutine propagate_c(nuclear,sims,Nsim,cloned)
    use communism
    implicit none
    type(mce) :: nuclear
    type(sim_pointer_array), dimension(Nsim) :: sims
    integer, intent(in) :: Nsim
    logical, intent(in) :: cloned

    _REAL_ :: tini,tfin
    integer :: i
    double complex, dimension(Nsim) :: D_fin

    tini=(sims(1)%sim%naesmd%tfemtoquantum-sims(1)%sim%naesmd%dtquantum*convtf)/convtf
    do i=1,10
        tfin = tini + sims(1)%sim%naesmd%dtquantum/10
!        if(.not.cloned) then !Do not propagate in the step after cloning
            call integrate_1step_c(tini,tfin,nuclear%D,D_fin,nuclear,Nsim)
            nuclear%D = D_fin
!        endif
        tini = tfin
    enddo

end subroutine propagate_c


subroutine integrate_1step_c(tini,tfin,c,c_fin,nuclear,Nsim)
    use communism
    implicit none
    _REAL_, intent(in) :: tini, tfin
    double complex, intent(in), dimension(Nsim) :: c
    double complex, intent(out), dimension(Nsim) :: c_fin
    type(mce) :: nuclear
    integer, intent(in) :: Nsim

    _REAL_ dt
    _REAL_, dimension(4) :: pfac1
    _REAL_, dimension(4) :: pfac2
    integer irk4, i, j
    double complex, dimension(Nsim) :: krk4, krk4_aux, tmp1,tmp2
    _REAL_ :: tini_aux
    double complex :: alpha, beta
    double complex, dimension(Nsim,Nsim) :: Heff_t

    dt = tfin - tini
    pfac1(1) = dt*0.d0
    pfac1(2) = dt*5.d-1
    pfac1(3) = dt*5.d-1
    pfac1(4) = dt*1.d0
    pfac2(1) = dt/6.d0
    pfac2(2) = dt/3.d0
    pfac2(3) = dt/3.d0
    pfac2(4) = dt/6.d0
    alpha = cmplx(1.0d0,0.0d0)
    beta = cmplx(0.0d0,0.0d0)

    c_fin = c
    krk4 = 0.0d0
    tini_aux = tini
    do irk4 = 1, 4
!        krk4 = fp(t0 + pfac1(irk4), f0 + krk4 * pfac1(irk4))
        tini_aux = tini + pfac1(irk4)
        call interpolate_Heff(tini,tini_aux,tfin,nuclear%Heff_old(1:Nsim,1:Nsim),Heff_t,nuclear%Heff(1:Nsim,1:Nsim),Nsim)
        krk4_aux = krk4
        do j = 1, Nsim
                tmp2(j) = c(j)+krk4(j)*pfac1(irk4)
        enddo
        call zgemm('N', 'N', Nsim, 1,  Nsim, alpha, Heff_t, Nsim, tmp2, Nsim, beta, tmp1, Nsim)

        krk4 = (0.0d0,-1.0d0)*tmp1
!        ff = ff + pfac2(irk4) * krk4
        c_fin = c_fin + pfac2(irk4) * krk4
    enddo !irk4
    
end subroutine integrate_1step_c


subroutine interpolate_Heff(tini,t,tfin,Heff_old,Heff_t,Heff,Nsim)
    implicit none
    _REAL_, intent(in) :: tini,t,tfin
    double complex, dimension(Nsim,Nsim), intent(in) :: Heff, Heff_old
    double complex, dimension(Nsim,Nsim), intent(out) :: Heff_t
    integer, intent(in) :: Nsim
    integer :: i,j

    Heff_t=(Heff*(t-tini)-Heff_old*(t-tfin))/(tfin-tini)

end subroutine interpolate_Heff


subroutine write_mce_data_0(sims, nuclear, Nsim)
    use communism
    implicit none   
    type(mce), intent(in) :: nuclear
    type(sim_pointer_array), dimension(Nsim), intent(in) :: sims
    integer, intent(in) :: Nsim
    _REAL_, dimension(sims(1)%sim%excN) :: pop

!Writing MCE output at t = 0 fs
    pop = sims(1)%sim%naesmd%yg(1:sims(1)%sim%excN)**2
    write(nuclear%outfile_1,999) Nsim, sims(1)%sim%naesmd%tfemto, pop, sum(pop)
    write(nuclear%outfile_2,999) Nsim, sims(1)%sim%naesmd%tfemto, nuclear%D(1:Nsim)

999     FORMAT(I3,1000(1X,F18.10))

end subroutine write_mce_data_0


subroutine overlap_nuclear(x1, p1, Ph1, x2, p2, Ph2, f2, ndim, s, s_dot, ek, Sdx, Sdp, mass, w)
implicit none

interface
  function overlap_CG(x1, p1, w1, x2, p2, w2)
    double complex :: overlap_CG
    _REAL_ :: x1, p1, w1, x2, p2, w2
  end function overlap_CG

  function overlap_dx_CG(x1, p1, w1, x2, p2, w2, p)
    double complex :: overlap_dx_CG
    _REAL_ :: x1, p1, w1, x2, p2, w2
    logical, optional :: p
  end function overlap_dx_CG

  function overlap_dp_CG(x1, p1, w1, x2, p2, w2, p)
    double complex :: overlap_dp_CG
    _REAL_ :: x1, p1, w1, x2, p2, w2
    logical, optional :: p
  end function overlap_dp_CG

  function overlap_d2x_CG(x1, p1, w1, x2, p2, w2, p)
    double complex :: overlap_d2x_CG
    _REAL_ :: x1, p1, w1, x2, p2, w2
    logical, optional :: p
  end function overlap_d2x_CG
end interface

_REAL_, dimension (ndim) :: x1, p1, x2, p2, f2
double complex s, s0(ndim) ,s1, ek, ek0, s_dot, Sdx(ndim), Sd2x(ndim), Sdp(ndim)
integer i, ndim
double precision, intent(in) :: mass(ndim), w(ndim)
_REAL_ :: MassToAU, PhaseDot, Ph1, Ph2
MassToAu=1822.887

do i=1,ndim
  S0(i) =  overlap_CG(x1(i), p1(i), w(i), x2(i), p2(i), w(i))
  Sdx(i) = overlap_dx_CG(  x1(i), p1(i), w(i), x2(i), p2(i), w(i),.true.)
  Sd2x(i)= overlap_d2x_CG( x1(i), p1(i), w(i), x2(i), p2(i), w(i),.true.)
  Sdp(i) = overlap_dp_CG(  x1(i), p1(i), w(i), x2(i), p2(i), w(i),.true.)
end do
S=product(S0)!*exp(cmplx(0.,-Ph1+Ph2))

ek=-0.5*S*sum(Sd2x/mass)

PhaseDot=0.0d0
do i=1,ndim
  PhaseDot=PhaseDot+p2(i)**2/mass(i)/2.0d0
enddo

S_dot = (dot_product(p2/mass, Sdx) + dot_product( f2, Sdp)) * S! + (0.d0,1.d0)* PhaseDot * S
return
end subroutine overlap_nuclear

function overlap_CG( x_i, p_i, alpha_i, &
                     x_j, p_j, alpha_j  ) result ( S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the overlap between two complex gaussian
_REAL_, intent(in) ::  x_i, p_i, alpha_i, &
                                  x_j, p_j, alpha_j

double complex :: S_ij

_REAL_ :: delta_x, delta_p, &  ! difference is position & momentum
                     x_cent,           &  ! the centroid position
                     prefactor,        &
                     real_part,        &
                     imag_part

delta_x = x_i-x_j
delta_p = p_i-p_j

real_part = ( alpha_i*alpha_j*delta_x**2+ 0.25d0*delta_p**2 ) / &
                         ( alpha_i + alpha_j )

   prefactor = sqrt( 2.0 * sqrt( alpha_i*alpha_j ) / &
                              ( alpha_i + alpha_j ) )

   x_cent  = ( alpha_i*x_i + alpha_j*x_j ) / ( alpha_i + alpha_j )

   imag_part = ( p_i*x_i - p_j*x_j ) - x_cent * delta_p

   S_ij = prefactor * exp( -real_part + (0.,1.)*imag_part )

return
end function overlap_CG



function overlap_d2x(x1,p1,x2,p2,w) result (s2_ij)
!
implicit none
double complex preex, sup, s2_ij
_REAL_ :: x1,p1,x2,p2,w, q, p, p0
q=x1-x2
p=p1-p2
p0=(p1+p2)*0.5
preex=-0.5*w+(0.5*w*q-(0.,1.)*p0)**2
sup=-0.25*(q**2*w+p**2/w)+(0.,1.)*p0*q
s2_ij=preex*exp(sup)
return
end function overlap_d2x

function overlap_d2x_CG( x_i, p_i, alpha_i, &
                         x_j, p_j, alpha_j, &
                         prefactor         ) result ( d2x_S_ij )
implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expectation of the d^2/dx^2  between two complex gaussians
_REAL_, intent(in) :: x_i, p_i, alpha_i, &
                                    x_j, p_j, alpha_j
logical,optional                 :: prefactor   ! pre calculated overlap
double complex            :: d2x_S_ij, overlap_CG

double complex :: S_ij
_REAL_ :: P_ij, delta_x

delta_x = x_i - x_j
P_ij    = alpha_i*p_j + alpha_j*p_i


d2x_S_ij = -(+4.d0*(0.,1)* alpha_i*alpha_j * delta_x * P_ij    &
             +2.d0*        alpha_i*alpha_j *(alpha_i+alpha_j)  &
             -4.d0*delta_x**2 * alpha_i**2 * alpha_j**2        &
              + P_ij**2 ) / ( alpha_i + alpha_j )**2

if( present(prefactor) )then
   if( .not. prefactor )then
      d2x_S_ij = d2x_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j  )
   endif
else
   d2x_S_ij = d2x_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                     x_j, p_j, alpha_j  )
endif

return
end function overlap_d2x_CG

function overlap_dx_CG( x_i, p_i, alpha_i, &
                        x_j, p_j, alpha_j, &
                        prefactor        ) result (dx_S_ij )
implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expection of d/dx between two complex gaussians 
_REAL_, intent(in) ::  x_i, p_i, alpha_i,&
                                     x_j, p_j, alpha_j
logical,optional                 ::  prefactor
double complex            ::  dx_S_ij, overlap_CG

_REAL_ :: P_ij, delta_x

delta_x = x_i - x_j
P_ij    = alpha_i*p_j + alpha_j*p_i

! CAUTION SIGN FLIPPED
dx_S_ij =-( (0.,1.)*P_ij - 2.d0 * alpha_j*alpha_i*delta_x ) &
                 / ( alpha_i + alpha_j )

if( present(prefactor) )then
   if( .not. prefactor)then
      dx_S_ij = dx_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                      x_j, p_j, alpha_j  )
   endif
else
   dx_S_ij = dx_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                   x_j, p_j, alpha_j  )
endif

return
end function overlap_dx_CG


function overlap_dp_CG( x_i, p_i, alpha_i, &
                        x_j, p_j, alpha_j, &
                        prefactor        ) result (dp_S_ij )
implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expection of d/dp between two complex gaussians 
! right acting
_REAL_, intent(in) ::  x_i, p_i, alpha_i,&
                                     x_j, p_j, alpha_j
logical,optional                 ::  prefactor ! pre calculated overlap
double complex            ::  dp_S_ij, overlap_CG

_REAL_ :: delta_x, delta_p

delta_x = x_i - x_j
delta_p = p_i - p_j

dp_S_ij =   ( delta_p + 2.d0 * (0.,1.) * alpha_i*delta_x ) &
                / ( 2.d0* (alpha_i + alpha_j ) )

if( present(prefactor) )then
   if( .not. prefactor )then
      dp_S_ij = dp_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                      x_j, p_j, alpha_j  )
   endif
else
   dp_S_ij = dp_S_ij * overlap_CG( x_i, p_i, alpha_i, &
                                   x_j, p_j, alpha_j  )
endif

return
end function overlap_dp_CG

!end module
