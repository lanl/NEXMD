! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
! *********************************************************************         
! The following soubroutines calculate the one electron and two electron
! contribution to the fock matrix (with d-orbital implementation.
! They are NOT optimized yet.
!
!  By Taisung Lee (Rutgers, 2011)                                                                              
!                                                                               
! *********************************************************************

module qm2_fock_d

    public qm2_fock1_d, qm2_fock2_d, W2Fock_atompair
   
    private InitializeWPosition

    
contains 

  
    subroutine qm2_fock2_d(qm2_struct, qmmm_struct, F, PTOT, W)
        !***********************************************************************
        !
        ! FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
        ! MATRIX
        ! ON INPUT  PTOT = TOTAL DENSITY MATRIX.
        !           W    = TWO-ELECTRON INTEGRAL MATRIX.
        !
        !  ON OUTPUT F   = PARTIAL FOCK MATRIX
        !
        !  Why isn't this used for excited states? !JAB !FIXME
        !***********************************************************************
        use qmmm_module, only : qm2_params, qmmm_mpi, qm2_structure
        use qmmm_struct_module, only : qmmm_struct_type
	implicit none

         type(qmmm_struct_type), intent(inout) :: qmmm_struct
         type(qm2_structure), intent(inout) :: qm2_struct
        _REAL_, intent(inout) :: F(:)
        _REAL_, intent(in) :: ptot(:)
        _REAL_, intent(in) :: W(:)

        !Local

        integer:: m,i,j, ij, ji, k, l, kl, lk, kk, ii, ia, ib, jk, kj, jj, ja, jb
        integer:: starting, size_ij

        if (.not.qmmm_struct%w_position_initialized) call InitializeWPosition(qmmm_struct)
 
        do ii=1,qmmm_struct%nquant_nlink
            IA=qm2_params%orb_loc(1,ii)
            IB=qm2_params%orb_loc(2,ii)
            M=0
            do J=IA,IB
                do K=IA,IB
                    M=M+1
                    JK=MIN(J,K)
                    KJ=K+J-JK
                    JK=JK+qm2_params%pascal_tri1(KJ)
                    qm2_struct%fock2_PTOT2(M,ii)=PTOT(JK)
                end do
            end do
        end do
          
        do ii=1,qmmm_struct%nquant_nlink
            do jj=1, qmmm_struct%nquant_nlink
        
                k=qm2_params%orb_loc(1,ii)
                l=qm2_params%orb_loc(1,jj)
                if (ii.ne.jj) then
            
                    i=qm2_params%natomic_orbs(ii)
                    j=qm2_params%natomic_orbs(jj)
                    starting=qmmm_struct%w_position(ii,jj)
                    size_ij=( i*(i+1)*j*(j+1) ) /4
                    !write(*,*) "starting ", starting,size_ij,i,j,k,l
                    call W2Fock_atompair(W(starting:starting+size_ij-1), F, ptot, &
                        i, j, k, l, qmmm_struct%W2Fock_atompair_initialized,qmmm_struct%w_index)
                end if
            end do
        end do

        return
    end subroutine qm2_fock2_d



    subroutine qm2_fock1_d(qmmm_struct, F, PTOT)

        use constants          , only : fourth
        use ElementOrbitalIndex, only : MaxValenceDimension, &
            Index1_2Electron, IntRep,  &
            IntRf1, IntRf2, IntIJ, IntKL
        use MNDOChargeSeparation, only: GetOneCenter2Electron
        use qmmm_module         , only : qmmm_mpi, qm2_params
        use qmmm_struct_module, only : qmmm_struct_type

        implicit none

        type(qmmm_struct_type), intent(inout) :: qmmm_struct

        _REAL_, intent(inout) :: F(:)
        _REAL_, intent(in) :: PTOT(:)
 
        ! local
 
        _REAL_::F_local(MaxValenceDimension), P_local(MaxValenceDimension)
    
        integer::i,j,k,i1,i2,ij,kl,counter
        integer::qmType
  
  
        if (.not.qmmm_struct%w_position_initialized) call InitializeWPosition(qmmm_struct)

        ! first calculate the SP contributions
        call qm2_fock1(qmmm_struct, F,PTOT)
    
        do i=1,qmmm_struct%nquant_nlink
            qmType=qmmm_struct%qm_atom_type(i)
            k=qm2_params%natomic_orbs(i)
        
            if (k.ge.9) then !only do those atoms w/ d-orbitals
        
                ! the integrals
                if ((.not.qmmm_struct%qm2_fock1_d_initialized).or. (qmType.ne.qmmm_struct%qmType_saved)) then
                    do j=1,Index1_2Electron
                        i1 = IntRf1(j)
                        i2 = IntRf2(j)
                        qmmm_struct%W(j) = GetOneCenter2Electron(qmType, IntRep(j))
                        if(i1>0) qmmm_struct%W(j) = qmmm_struct%W(j)-fourth*GetOneCenter2Electron(qmType,i1)
                        if(i2>0) qmmm_struct%W(j) = qmmm_struct%W(j)-fourth*GetOneCenter2Electron(qmType,i2)
                    end do
                    qmmm_struct%qmType_saved=qmType
                    qmmm_struct%qm2_fock1_d_initialized=.true.
                end if
            
                ! Copy the density matrix to local
            
                i1=qm2_params%orb_loc(1,i)
                i2=qm2_params%orb_loc(2,i)
                counter=0
                do j=i1, i2
                    do k=qm2_params%pascal_tri1(j)+i1, qm2_params%pascal_tri2(j)-1
                        counter=counter+1
                        P_local(counter)=PTOT(k)*2.d0  ! off-diag terms need to be
                                ! doubled to account for the half matrix summation
                    end do
                    counter=counter+1
                    P_local(counter)=PTOT(qm2_params%pascal_tri2(j))
                end do
            
                ! the coulombic contribution
                F_local=0.0D0
                do j=1,Index1_2Electron
                    ij=IntIJ(j)
                    kl=IntKL(j)
                    F_local(ij)=F_local(ij)+P_local(kl)*qmmm_struct%W(j)
                end do
            
                ! the exchange contribution
                ! no exchange for the RHF case
                ! the UHF case--to be done later
            
                ! add local contribution back to the Fock matrix
                i1=qm2_params%orb_loc(1,i)
                i2=qm2_params%orb_loc(2,i)
                counter=0
                do j=i1, i2
                    do k=qm2_params%pascal_tri1(j)+i1, qm2_params%pascal_tri2(j)
                        counter=counter+1
                        F(k)=F(k)+F_local(counter)
                    end do
                end do
           
            end if
        end do

        return
    end subroutine qm2_fock1_d

    subroutine W2Fock_atompair(W, F, D, norbs_a, norbs_b,  &
        na_starting, nb_starting, initialized, w_index)

        use constants, only : half

        implicit none

        integer, intent(in)::norbs_a, norbs_b
        integer, intent(in)::na_starting
        integer, intent(in)::nb_starting

        _REAL_, intent(in)::W(*), D(*)
        _REAL_, intent(inout)::F(*)
        logical, intent(inout)::initialized

        integer, intent(inout):: w_index(:,:,:,:,:,:)
        !the w_index needs "lots" of memory but significantly improves
        !the size and readability of the code
        
        !local
        integer, parameter::orbital_length(3)=(/ 1, 4, 9 /)
        integer, parameter::pair_length(3)=(/1, 10, 45 /)
        integer::starting(9)

        integer::i,j,k,l,a1,a2,b1,b2, ii,jj
        integer::na, nb, nn, counter_F, counter_ij, counter_kl, location
        integer::index_Fa(norbs_a,norbs_a), index_Fb(norbs_b,norbs_b)
        integer::index_Fab(norbs_a,norbs_b)
    
        _REAL_::temp1, temp2, tt1(81), tt2(81), wlocal(9,9,9,9)
        logical::toCalc(81)

        !write(*,*) "norbs_a = ",norbs_a,"norbs_b = ", norbs_b,"na_starting = ", na_starting,"nb_starting = ", nb_starting
        if (.not. initialized) then
    
            w_index=-1
            do na=1,3
                do nb=1,3
                    counter_ij=0
                    do i=1, orbital_length(na)
                        do j=1, i
                            counter_ij=counter_ij+1
                            counter_kl=0
                            do k=1, orbital_length(nb)
                                do l=1, k
                                    counter_kl=counter_kl+1
                                    location=(counter_ij-1)*pair_length(nb)+counter_kl
                                    w_index(i,j,k,l, na, nb)=location
                                    w_index(i,j,l,k, na, nb)=location
                                    w_index(j,i,k,l, na, nb)=location
                                    w_index(j,i,l,k, na, nb)=location
                                end do !l
                            end do !k
                        end do !j
                    end do !i
       
                end do ! nb
            end do ! na

            initialized=.true.
        endif
   
        do i=1, norbs_a
            ii=na_starting+i-1
            jj=ii*(ii-1)/2+na_starting-1
            do j=1, i
                location=jj+j
                !write(*,*)"Fa i = ",i,"j = " ,j,"location =",location;
                index_Fa(i, j)=location
                index_Fa(j, i)=location
            end do !nb
        end do
        do k=1, norbs_b
            ii=nb_starting+k-1
            jj=ii*(ii-1)/2+nb_starting-1
            do l=1, k
                location=jj+l
                !write(*,*)"Fb k = ",k,"l = " ,l,"location =",location;
                index_Fb(k, l)=location
                index_Fb(l, k)=location
            end do !nb
        end do
    
        if (na_starting.ge.nb_starting) then
            k=0
            a1=na_starting
            b1=nb_starting
            a2=norbs_a
            b2=norbs_b
        else
            k=1
            a1=nb_starting
            b1=na_starting
            a2=norbs_b
            b2=norbs_a
        end if
        do i=1, a2
            ii=a1+i-1
            jj=ii*(ii-1)/2
            do j=1, b2
                location=jj+b1+j-1
                if (k==0) then
                    !write(*,*)"Fab k=0  i = ",i,"j = " ,j,"location =",location;
                    index_Fab(i, j)=location
                else
                    !write(*,*)"Fab k=1  j = ",j,"i = " ,i,"location =",location;
                    index_Fab(j, i)=location
                end if
            end do !nb
        end do
 
        ! real stuff

        na=sqrt(norbs_a+1.0D-5)
        nb=sqrt(norbs_b+1.0D-5)
          
        ! for the coloumbic part, calculate all na/nb pairs but be carful about the index for W
        ii=1
        do k=1, norbs_b
            do l=1, k-1
                tt1(ii)=D(index_Fb(l,k))*2
                ii=ii+1
            end do
            tt1(ii)=D(index_Fb(k,k))
            ii=ii+1
        end do
                 
        do i=1, norbs_a
            do j=1, i
       
                temp1=0.0D0
                if (na_starting.gt.nb_starting) then
                    ii=1
                    do k=1, norbs_b
                        do l=1, k
                            temp1=temp1+tt1(ii)* &
                                W(w_index(i,j,k,l,na,nb))
                            !write(*,*)"w_index_(",i,",",j,",",k,",",l,") = ", w_index(i,j,k,l,na,nb)
                            !write(*,*) "W_(",w_index(i,j,k,l,na,nb),")=", W(w_index(i,j,k,l,na,nb))
                            ii=ii+1
                        end do
                    end do
                else
                    ii=1
                    do k=1, norbs_b
                        do l=1, k
                            temp1=temp1+tt1(ii)* &
                                W(w_index(k,l,i,j,nb,na));
                            !write(*,*)"_w_index(",k,",",l,",",i,",",j,") = ", w_index(k,l,i,j,na,nb);
                            !write(*,*) "_W(",w_index(k,l,i,j,na,nb),")=", W(w_index(k,l,i,j,na,nb));
                            ii=ii+1;
                        end do
                    end do
                end if  !na_starting.lt.nb_starting
                !write(*,*)"F(",i,",",j,") = ", F(index_Fa(i,j))
                !write(*,*)"D(",i,",",j,") =	", D(index_Fa(i,j))
                !write(*,*)"G(",i,",",j,") = ", temp1;
                F(index_Fa(i,j))=F(index_Fa(i,j))+temp1
            end do
        end do

        ! for the exchange part, only calculates when na_starting.gt.nb_starting
        ! to avoid double counting
        if (na_starting.gt.nb_starting) then
            !write(*,*)"na_starting = ", na_starting,"nb_starting = ", nb_starting
            ii=1
            do l=1, norbs_b
                do j=1, norbs_a
                    tt1(ii)=D(index_Fab(j,l))
                    if (abs(tt1(ii))>1.0d-6) then
                        toCalc(ii)=.true.
                    else
                        toCalc(ii)=.false.
                    end if
                    ii=ii+1
                end do
            end do
               
            do i=1, norbs_a
                do k=1, norbs_b
                    temp1=0.0D0
                    ii=1
                    do l=1, norbs_b
                        do j=1, norbs_a
                            !if (toCalc(ii)) then
                            temp1=temp1+tt1(ii)* &
                                W(w_index(i,j,k,l,na,nb))
                            !write(*,*)"Fab w_index(",i,",",j,",",k,",",l,") = ", w_index(i,j,k,l,na,nb);
                            !write(*,*) "W(",w_index(i,j,k,l,na,nb),")=", W(w_index(i,j,k,l,na,nb));
                            ii=ii+1
                            !end if
                        end do
                    end do
                    F(index_Fab(i,k))=F(index_Fab(i,k))-temp1*half

                end do
            end do
        end if

    end subroutine W2Fock_atompair

    subroutine InitializeWPosition(qmmm_struct)

        use qmmm_struct_module, only : qmmm_struct_type
	use qmmm_module, only : qm2_params
    
        type(qmmm_struct_type), intent(inout) :: qmmm_struct

        integer::i,j,k,ii,jj,kk,n
    
        if (qmmm_struct%w_position_initialized) return
    
        n=qmmm_struct%nquant_nlink
        allocate(qmmm_struct%w_position(n,n))
        
        kk=1
        qmmm_struct%w_position=-1
        do ii=1,n
            do jj=1, ii-1
                i=qm2_params%natomic_orbs(ii)
                j=qm2_params%natomic_orbs(jj)
                
                qmmm_struct%w_position(ii,jj)=kk
                qmmm_struct%w_position(jj,ii)=kk

                kk=kk+(i*(i+1)/2)*(j*(j+1)/2)
            end do
        end do
        qmmm_struct%w_position_initialized=.true.
    
    end subroutine InitializeWPosition


end module qm2_fock_d
