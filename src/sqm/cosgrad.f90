
!Calculate COSMO Gradient
subroutine cosgrad(cosmo_c_struct, rho1,rho2)

use cosmo_C, only :: cosmo_C_structure
use qmmm_module
implicit none
type(cosmo_C_structure), intent(inout) :: cosmo_c_struct
integer         ::      n,i,j,iatom,jatom,INFO
_REAL_          ::      GradA(3,cosmo_c_struct%nps,cosmo_c_struct%nps),GradB(3,cosmo_c_struct%lm61,cosmo_c_struct%nps),GradM(3,cosmo_c_struct%lm61,cosmo_c_struct%lm61),&
                        AinvBtGradA(3,cosmo_c_struct%lm61,cosmo_c_struct%nps)
_REAL_          ::      A(cosmo_c_struct%nps,cosmo_c_struct%nps),Ainv(cosmo_c_struct%nps,cosmo_c_struct%nps),AinvB(cosmo_c_struct%nps,cosmo_c_struct%lm61),IPIV(cosmo_c_struct%nps)
_REAL_          ::      GradE,ALPHA,BETA
_REAL_,intent(in)::     rho1(qm2_struct%norbs**2),rho2(qm2_struct%norbs**2)

!Calculate Grad(A)
do i=1,cosmo_c_struct%nps
    iatom=iatsp(i)
    do j=1,i
        jatom=iatsp(j)
        if (iatom.ne.jatom)
        do n=1,3
        GradA(n,i,j)=abs(cosmo_c_struct%amat((i-1)*i/2+j))**(-3)*(cosmo_c_struct%cosurf(n,i)-cosmo_c_struct%cosurf(n,j)
        GradA(n,j,i)=-GradA(n,i,j)
        A(i,j)=cosmo_c_struct%amat((i-1)*i/2)
        A(j,i)=A(i,j)
    enddo
enddo

!Calculate Grad(B)
do i= 1,cosmo_c_struct%nps
     iatom=iatsp(i)
     do j=1,cosmo_c_struct%lm61
        jatom=ATOM THAT J IS ON
        if (iatom.ne.jatom)
        do n=1,3
        GradB(n,i,j)=abs(cosmo_c_struct%bmat((j,i)))**(-3)*(cosmo_c_struct%cosurf(n,i)-qmmm_struct%qm_coords(n,jatom))
     enddo
enddo

!Calculate A^(-1)
Ainv(1:cosmo_c_struct%nps,1:cosmo_c_struct%nps)=0.0;
     do i=1,cosmo_c_struct%nps
           Ainv(i,i)=1.0; ! Unitary Diagonal matix 
     end do
CALL DGESV( cosmo_c_struct%nps, cosmo_c_struct%nps, A, cosmo_c_struct%nps, IPIV, Ainv, cosmo_c_struct%nps, INFO);

!Calculate Gradient
! GradB*AinvB + AinvB^t*GradB + AinvB^t*GradA*AinvB
BETA=0.0; ALPHA=1.0;
CALL DGEMM('N','N',cosmo_c_struct%nps, cosmo_c_struct%nps, cosmo_c_struct%lm61, ALPHA, Ainv, cosmo_c_struct%nps, B, cosmo_c_struct%lm61,BETA, AinvB, cosmo_c_struct%lm61 );

CALL DGEMM('T','N',cosmo_c_struct%lm61, cosmo_c_struct%nps, cosmo_c_struct%nps, ALPHA, AinvB, cosmo_c_struct%lm61, GradA(n,:,:),cosmo_c_struct%nps,BETA,AinvBtGradA(n,:,:), cosmo_c_struct%lm61 );
CALL DGEMM('T','N',cosmo_c_struct%lm61, cosmo_c_struct%nps, cosmo_c_struct%lm61, ALPHA, GradB, cosmo_c_struct%lm61, AinvB, cosmo_c_struct%lm61, BETA,GradM, cosmo_c_struct%lm61 ); 
!Transpose GradB? Equiv of first and second terms?
BETA=1.0; 
CALL DGEMM('T','N',cosmo_c_struct%lm61, cosmo_c_struct%nps, cosmo_c_struct%lm61, ALPHA, AinvB, cosmo_c_struct%lm61, GradB, cosmo_c_struct%lm61, BETA, GradM, cosmo_c_struct%lm61 );
CALL DGEMM('N','N',cosmo_c_struct%lm61, cosmo_c_struct%nps, cosmo_c_struct%lm61, ALPHA, AinvBtGradA, cosmo_c_struct%lm61, AinvB, cosmo_c_struct%lm61,BETA,GradM, cosmo_c_struct%lm61 );

!Calculate With Density Matrices
do i=1,cosmo_c_struct%lm61
     do j=1,cosmo_c_struct%lm61
          GradE=GradE+rho1(i)*rho2(j)*GradM*(i,j)
     enddo
enddo
GradE=GradE*cosmo_c_struct%a0*cosmo_c_struct%ev !Gradient in eV/A or eV/Bohr??

return
end subroutine
