
!Calculate COSMO Gradient
subroutine cosgrad(rho1,rho2)

use cosmo_C
use qmmm_module
implicit none

integer         ::      n,i,j,iatom,jatom,INFO
_REAL_          ::      GradA(3,nps,nps),GradB(3,lm61,nps),GradM(3,lm61,lm61),&
                        AinvBtGradA(3,lm61,nps)
_REAL_          ::      A(nps,nps),Ainv(nps,nps),AinvB(nps,lm61),IPIV(nps)
_REAL_          ::      GradE,ALPHA,BETA
_REAL_,intent(in)::     rho1(qm2_struct%norbs**2),rho2(qm2_struct%norbs**2)

!Calculate Grad(A)
do i=1,nps
    iatom=iatsp(i)
    do j=1,i
        jatom=iatsp(j)
        if (iatom.ne.jatom)
        do n=1,3
        GradA(n,i,j)=abs(amat((i-1)*i/2+j))**(-3)*(cosurf(n,i)-cosurf(n,j)
        GradA(n,j,i)=-GradA(n,i,j)
        A(i,j)=amat((i-1)*i/2)
        A(j,i)=A(i,j)
    enddo
enddo

!Calculate Grad(B)
do i= 1,nps
     iatom=iatsp(i)
     do j=1,lm61
        jatom=ATOM THAT J IS ON
        if (iatom.ne.jatom)
        do n=1,3
        GradB(n,i,j)=abs(bmat((j,i)))**(-3)*(cosurf(n,i)-qmmm_struct%qm_coords(n,jatom))
     enddo
enddo

!Calculate A^(-1)
Ainv(1:nps,1:nps)=0.0;
     do i=1,nps
           Ainv(i,i)=1.0; ! Unitary Diagonal matix 
     end do
CALL DGESV( nps, nps, A, nps, IPIV, Ainv, nps, INFO);

!Calculate Gradient
! GradB*AinvB + AinvB^t*GradB + AinvB^t*GradA*AinvB
BETA=0.0; ALPHA=1.0;
CALL DGEMM('N','N',nps, nps, lm61, ALPHA, Ainv, nps, B, lm61,BETA, AinvB, lm61 );

CALL DGEMM('T','N',lm61, nps, nps, ALPHA, AinvB, lm61, GradA(n,:,:),nps,BETA,AinvBtGradA(n,:,:), lm61 );
CALL DGEMM('T','N',lm61, nps, lm61, ALPHA, GradB, lm61, AinvB, lm61, BETA,GradM, lm61 ); 
!Transpose GradB? Equiv of first and second terms?
BETA=1.0; 
CALL DGEMM('T','N',lm61, nps, lm61, ALPHA, AinvB, lm61, GradB, lm61, BETA, GradM, lm61 );
CALL DGEMM('N','N',lm61, nps, lm61, ALPHA, AinvBtGradA, lm61, AinvB, lm61,BETA,GradM, lm61 );

!Calculate With Density Matrices
do i=1,lm61
     do j=1,lm61
          GradE=GradE+rho1(i)*rho2(j)*GradM*(i,j)
     enddo
enddo
GradE=GradE*a0*ev !Gradient in eV/A or eV/Bohr??

return
end subroutine
