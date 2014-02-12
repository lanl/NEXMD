! <compile=optimized>

#include "dprec.fh"

subroutine dftb_matrix_diag(NA,N,A,EW,IER)
   use qm2_dftb_module, only: scr_space
   implicit none
!! Passed in
   integer, intent(in)    :: NA        ! linear dimension of A
   integer, intent(in)    :: N         ! Dimension of the problem
   integer, intent(out)   :: IER       ! Error code
   _REAL_ , intent(inout) :: A(NA,N)   ! (I):Matrix to be diagonalized
                                       ! (O): Eigenvectors
   _REAL_ , intent(out)   :: EW(N)     ! Eigenvalues

   call DSYEV('V','L',N,A,NA,EW,scr_space,3*NA,IER)
   return
end subroutine dftb_matrix_diag

!=============================================
SUBROUTINE EWEVGE (NA,NB,N,A,B,EW,IEV,IER)
   use qm2_dftb_module, only: scr_space

   ! Solves a general eigenvalue problem

   implicit none 

!! Passed in
   integer, intent(in)    :: NA       ! Linear dimension of A
   integer, intent(in)    :: NB       ! Linear dimension of B
   integer, intent(in)    :: N        ! Dimension of Problem
   integer, intent(in)    :: IEV      ! Wether to calculate eigenvalues (1==yes)
   integer, intent(out)   :: IER      ! Error code
   _REAL_ , intent(inout) :: A(NA,N)  ! (I):Matrix A, (O):Eigenvectors
   _REAL_ , intent(inout) :: B(NB,N)  ! Matrix B
   _REAL_ , intent(out)   :: EW(N)    ! Eigenvalues

   CALL DSYGV(IEV,'V','L',N,A,NA,B,NB,EW,scr_space,3*NA,IER)
   
   RETURN
end SUBROUTINE EWEVGE

    ! ** To use DSYGVD (divide-and-conquer) **
    !DIMENSION  A(NA,N),B(NB,N),EW(N),H(N),AUX(3*NA)
    !integer :: IWORK(5*N+3)
    !integer :: LWORK, LIWORK

    !LWORK = 1 + 6*N + 2*N*N
    !CALL DSYGVD(IEV,'V','L',N,A,NA,B,NB,EW,AUX,LWORK,IWORK,5*N+3,IER)
