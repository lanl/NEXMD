#include "dprec.fh"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!This file contains subroutines to include COSMO and Onsager type continuum
!models in ground and excited state calculations. Routines specific to COSMO
!cavity construction, cholesky factorization, and COSMO gradients are found in
!cosmo.f90
!
! Josiah A. Bjorgaard, Vasyl Kuzmenko, Kirill Velizhanin 
! 2013-2014 Los Alamos National Laboratory
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!THIS SUBROUTINE REPLACES Lxi IN THE DAVIDSON ROUTINE WHEN A SOLVENT MODEL IS
!USED, IT INCLUDES ALL SOLVENT MODELS DESIGNATED BY solvent_model and
!potential_type
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine Lxi_testing(u1,v1,solvent_model)
    use qm2_davidson_module
    use qmmm_module,only:qm2_struct;
    use cosmo_C, only: v_solvent_difdens,v_solvent_xi,potential_type, EF;

    implicit none

    _REAL_ u1(qm2ds%Nrpa),v1(qm2ds%Nrpa)
    _REAL_ f,f1,f2,fs1,ddot
    integer i,j,p,h
    integer one
    integer solvent_model
    _REAL_  tmp(qm2ds%nb,qm2ds%nb),tmp2(qm2ds%nb,qm2ds%nb);

    parameter (one=1)

    fs1=0;
    qm2ds%xi=0.d0
    call mo2site(u1,qm2ds%xi,qm2ds%eta) !Change basis of guess vector of Davidson from M.O to A.O
    qm2ds%eta=0.0;
    call Vxi(qm2ds%xi,qm2ds%eta);    !Calculate Vacuum Electron Correlation
    !!SELECT SOLVENT MODEL AND POTENTIAL TYPE
    if ((solvent_model.eq.1)) then !1:Linear Response
        tmp=0.d0;
        if (potential_type.eq.3) then !COSMO Potential
            call VxiM(qm2ds%xi,tmp);
        elseif (potential_type.eq.2) then !Onsager Potential
            call rcnfld(tmp,qm2ds%xi,qm2ds%nb)
        elseif (potential_type.eq.1) then !testing
            do i=1,qm2ds%nb; tmp(i,i)=qm2ds%eta(qm2ds%nb*(i-1)+i); enddo !double diag vac correlation
            endif
            tmp=tmp
            call VxiM_end(qm2ds%eta,tmp);        !Add selected potential to vacuum correlation
        elseif (solvent_model.eq.99) then !For Z-vector equation if different
            tmp=0.d0;
            if (potential_type.eq.3) then !COSMO Potential
                call VxiM(qm2ds%xi,tmp);
            elseif (potential_type.eq.2) then !Onsager Potential
                call rcnfld(tmp,qm2ds%xi,qm2ds%nb)
            elseif (potential_type.eq.1) then !testing
                do i=1,qm2ds%nb; tmp(i,i)=qm2ds%eta(qm2ds%nb*(i-1)+i); enddo !double diag vac correlation
                endif
                tmp=2*tmp
                call VxiM_end(qm2ds%eta,tmp);            !Add selected potential to vacuum correlation
            elseif (solvent_model.eq.98) then !For Z-vector equation with SS model
                tmp=0.d0;
                if (potential_type.eq.3) then !COSMO Potential
                    call VxiM(qm2ds%xi,tmp);
                elseif (potential_type.eq.2) then !Onsager Potential
                    call rcnfld(tmp,qm2ds%xi,qm2ds%nb)
                elseif (potential_type.eq.1) then !testing
                    do i=1,qm2ds%nb; tmp(i,i)=qm2ds%eta(qm2ds%nb*(i-1)+i); enddo !double diag vac correlation
                    endif
                    call VxiM_end(qm2ds%eta,tmp);                    !Add selected potential to vacuumcorrelation
                    !Commutator is performed here for State Specific Solvent Routines
                    tmp=0.d0;
                    call commutator(qm2ds%xi,v_solvent_difdens,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                elseif(solvent_model.eq.2) then ! 2: State Specific [V_s(T+Z),xi]
                    tmp=0.d0;
                    !Commutator is performed here for State Specific Solvent Routines
                    call commutator(qm2ds%xi,v_solvent_difdens,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                elseif(solvent_model.eq.3) then !3: State Specific [V_s(xi),xi]
                    call commutator(v_solvent_xi,qm2ds%xi,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                elseif(solvent_model.eq.5) then !5: Variational State Specific term
                    call commutator(qm2ds%xi,v_solvent_difdens,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                    write(6,*)'Adding variational term in Liouville operator'
                elseif(solvent_model.eq.6) then!6: Solve nonlinear Liouville equation testing
                    call commutator(qm2ds%eta,qm2ds%xi,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                elseif(solvent_model.eq.10) then!10: NO GS Solvent test
                    !call addnuc(tmp(1,1));
                    tmp=0.d0; tmp2=0.d0;
                    call rcnfld_fock(tmp,qm2_struct%den_matrix,qm2ds%Nb)
                    call unpacking(qm2ds%Nb,tmp,tmp2,'s'); tmp=0.d0;
                    call commutator(tmp2,qm2ds%xi,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                elseif(solvent_model.eq.7) then!7: combined LR and VE
                endif
 
                !add constant electric field potential [V_E,xi]
                if(EF.eq.2) then!Constant Electric Field in ES only
                    tmp=0.d0; tmp2=0.d0
                    call efield_fock(tmp,qm2ds%Nb)
                    call unpacking(qm2ds%Nb,tmp,tmp2,'s'); tmp=0.d0
                    call commutator(qm2ds%xi,tmp2,qm2ds%Nb,tmp,.false.)
                    call VxiM_end(qm2ds%eta,tmp)
                endif

                call site2mo(qm2ds%xi,qm2ds%eta,v1);                !Change basis of xi again to M.O.

                i=0
                do p=1,qm2ds%Np
                    do h=qm2ds%Np+1,qm2ds%Nb
                        i=i+1
                        f=qm2ds%ehf(h)-qm2ds%ehf(p);
                        !JAB Test
                        v1(i)=v1(i)+f*u1(i)
                        v1(i+qm2ds%Ncis)=-(v1(i+qm2ds%Ncis)+f*u1(i+qm2ds%Ncis))
                    end do
                enddo
                do i=1,qm2ds%Nrpa
                    qm2ds%temp1(i)=qm2ds%v0(i,1)
                    qm2ds%temp2(i)=u1(i)
                end do
 
                if(qm2ds%Mj.gt.0) then
                    write(6,*) 'In the weird block' !!JAB
                    fs1=qm2ds%fs+qm2ds%e0(qm2ds%Mj)
                    do j=1,qm2ds%Mj
                        f1=fs1*(ddot(qm2ds%Ncis,qm2ds%v0(1,j),one,u1(1),one) &
                            -ddot(qm2ds%Ncis,qm2ds%v0(qm2ds%Ncis+1,j),one,u1(qm2ds%Ncis+1),one))
                        f2=fs1*(ddot(qm2ds%Ncis,qm2ds%v0(qm2ds%Ncis+1,j),one,u1(1),one) &
                            -ddot(qm2ds%Ncis,qm2ds%v0(1,j),one,u1(qm2ds%Ncis+1),one))

                        call daxpy(qm2ds%Ncis,f1,qm2ds%v0(1,j),one,v1(1),one)
                        call daxpy(qm2ds%Ncis,f2,qm2ds%v0(1+qm2ds%Ncis,j),one,v1(1),one)
                        call daxpy(qm2ds%Ncis,f1,qm2ds%v0(1+qm2ds%Ncis,j), &
                            one,v1(1+qm2ds%Ncis),one)
                        call daxpy(qm2ds%Ncis,f2,qm2ds%v0(1,j),one,v1(1+qm2ds%Ncis),one)

                    end do
                endif
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !THIS SUBROUTINE CALCULATES THE SOLVENT POTENTIAL OPERATOR V_S FOR COSMO USING CHOLESKY FACTORIZATION
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine VxiM(xi,Vximmat)
                use qm2_davidson_module
                use cosmo_C,only:fepsi,nps,lm61,a0,ev,amat,bmat,ipiden,gden,nsetf
                implicit none
                _REAL_, intent(in)::xi(qm2ds%nb,qm2ds%nb)
                _REAL_, intent(inout)::Vximmat(qm2ds%nb,qm2ds%nb)
                _REAL_ :: density(lm61),charges(nps),phi(nps)
                _REAL_ :: p(qm2ds%nb*(qm2ds%nb+1)/2),fcon
                integer:: i,im

                fcon=fepsi*a0*ev !scaling factor
                call packing(qm2ds%nb,xi,p,'s') !Note that the factor of two for diagonal elements which happens in this subroutine is corrected for by gden
                ! FIRST CALCULATE QDENEL FROM DENSITY MATRIX
                ! gden is -2 for diagonal and -1 for nondiagonal
                do i=1,lm61
                    density(i)=gden(i)*p(ipiden(i))
                end do
                !  NOW CALCULATE PHIEL FROM BMAT*QDENEL
                call  DGEMV ( 'T', lm61, nps, 1.d0, bmat, lm61, density, 1, 0.d0, phi, 1 ) !LAPACK subroutine for above
                ! NOW CALCULATE CHARGES
                call coscl2(amat,nsetf,charges,phi,nps) !Use cholesky factorization routine
 
                ! NOW ADD BMAT*QSCEL TO FOCK MATRIX
                call  DGEMV ( 'N', lm61, nps, fcon, bmat, lm61, charges, 1, 0.d0, density, 1 ) !LAPACK subroutine for above
                p=0.d0
                do i=1,lm61;
                    im=ipiden(i); p(im)=density(i);                !in A.U.
                enddo
                call unpacking(qm2ds%nb,p,Vximmat,'s')
            end subroutine VxiM

            !^^^^^^^^^^^^^^^^^^^^^^^
            !!ADD tmp to Vximmat...
            !^^^^^^^^^^^^^^^^^^^^^^^
            subroutine VxiM_end(Vximmat,tmp)
                use qmmm_module,only:qm2_struct;
                implicit none
                _REAL_,intent(inout)::tmp(qm2_struct%norbs,qm2_struct%norbs);
                _REAL_,intent(inout)::Vximmat(qm2_struct%norbs,qm2_struct%norbs);
                      !Modifing Vxi=Vxi+Mxi;
                Vximmat(1:qm2_struct%norbs,1:qm2_struct%norbs)= &
                    Vximmat(1:qm2_struct%norbs,1:qm2_struct%norbs)+ &
                    tmp(1:qm2_struct%norbs,1:qm2_struct%norbs);
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !!ADD tmp to Vximmat if tmp is the vector of a diagonal matrix
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine VxiM_end45(Vximmat,tmp)
                use qmmm_module,only:qm2_struct;
                implicit none
                _REAL_,intent(inout)::tmp(qm2_struct%norbs);
                _REAL_,intent(inout)::Vximmat(qm2_struct%norbs,qm2_struct%norbs);
                integer i;
                do i=1,qm2_struct%norbs
                    Vximmat(i,i)=Vximmat(i,i)+2*tmp(i);
                end do
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !SUBROUTINE FOR DIRECT CALCULATION OF A COMMUTATOR
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine commutator(A,B,matrix_size,res,flag);
                integer matrix_size;
                logical flag;
                _REAL_, dimension(matrix_size,matrix_size) :: A,B,Res;
                _REAL_, dimension(:,:), allocatable :: tmp
                _REAL_ ALPHA;
                _REAL_ BETA;

                ALPHA=1.0D+0;
                BETA=0.0D+0;
                allocate(tmp(matrix_size,matrix_size));
                if (flag) then !Transpose A if flag
                    CALL DGEMM('T','N', matrix_size , matrix_size , matrix_size , ALPHA , A,  &
                        matrix_size , B , matrix_size , BETA , res , matrix_size);
                    CALL DGEMM('N','T', matrix_size , matrix_size , matrix_size , ALPHA , B,  &
                        matrix_size , A , matrix_size , BETA , tmp , matrix_size);
                else !Otherwise no transpose
                    CALL DGEMM('N','N', matrix_size , matrix_size , matrix_size , ALPHA , A,  &
                        matrix_size , B , matrix_size , BETA , res , matrix_size);
                    CALL DGEMM('N','N', matrix_size , matrix_size , matrix_size , ALPHA , B, &
                        matrix_size , A , matrix_size , BETA , tmp , matrix_size);
                end if;
                res=res-tmp;
                deallocate(tmp);
                return
            end


            !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            ! Onsager-type solvent model (single point and gradient)functions
            ! and Electric field potentials are below
            !
            ! Josiah A. Bjorgaard, Vasyl Kuzmenko, Kirill Velizhanin
            ! 2013-2014 Los Alamos National Laboratory
            !
            !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE DIPOLE CAVITY SCREENING WITHOUT NUCLEAR TERMS FOR
            !EXCITED STATE CALCULATIONS
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfld(f,p,n)
                use qmmm_module,only:qmmm_struct;
                use cosmo_C, only: fepsi,onsager_radius;
                use constants, only : one, BOHRS_TO_A, AU_TO_EV
                implicit none;
                integer n,i,j;
                _REAL_ f(n,n);
                _REAL_ p(n,n);
                _REAL_ tmp(n,n);
                _REAL_ dipx(n,n),dipy(n,n),dipz(n,n);
                _REAL_ dip(3,n*(n+1)/2);
                _REAL_ scaled(3),elec_dip(3);
                tmp=0.d0
                dip=0.d0
                dipx=0.d0
                dipy=0.d0
                dipz=0.d0
                !!GET ELEC DIP MATRIX
                call get_dipole_matrix(qmmm_struct%qm_coords, dip)
                call unpacking(n,dip(1,:),dipx,'s')
                call unpacking(n,dip(2,:),dipy,'s')
                call unpacking(n,dip(3,:),dipz,'s')

                !!CALC TOTAL AND ELECTRIC DIP
                elec_dip(1:3)=0.0
                do i=1,n
                    do j=1,n
                        elec_dip(1) = elec_dip(1)+P(i,j)*dipx(i,j)
                        elec_dip(2) = elec_dip(2)+P(i,j)*dipy(i,j)
                        elec_dip(3) = elec_dip(3)+P(i,j)*dipz(i,j)
                    end do
                end do

                !!CALCULATE SCALING FACTOR AND INDUCED DIPOLE
                scaled=(fepsi/onsager_radius**3)*elec_dip*BOHRS_TO_A*AU_TO_EV;                ! for units of eV
        
                !!CALCULATE REACTION FIELD POTENTIAL
                f=f-2.d0*(scaled(1)*dipx+scaled(2)*dipy+scaled(3)*dipz);                 !in eV
                return
            end

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE SOLVENT MODEL GRADIENT FOR TWO TRIANGULAR MATRICES i.e. tr(F(rho)Z)
            !WITH OR WITHOUT NUCLEAR PART OF FOCK MATRIX
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfldgrad2(dxyz,p1,p2,n,calc_nuc)
                use qmmm_module,only: qm2_params,qmmm_struct;
                use cosmo_C, only: fepsi,onsager_radius,numat;
                use constants, only : one, BOHRS_TO_A, AU_TO_EV,EV_TO_KCAL
                implicit none;

                integer::  i,i1,i2,i3,j,k;
                integer::  n
                _REAL_ ::  dxyz(3,numat),origin(3,numat)
                _REAL_ ::  p1(n*(n+1)/2),p2(n*(n+1)/2)
                _REAL_ ::  dip(3,n*(n+1)/2);
                _REAL_ ::  scaled1(3),scaled2(3),scaled_nuc(3)&
                    ,elec_dip1(3),elec_dip2(3),nuc_dip(3),q_elec1,q_elec2
                logical::  calc_nuc

                dxyz=0.d0; elec_dip1=0.d0; elec_dip2=0.d0
                nuc_dip=0.d0;dip=0.d0;q_elec1=0.d0;q_elec2=0.d0

                call centercoords(origin)
                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                nuc_dip(1)=sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=sum(dip(2,:))
                nuc_dip(3)=sum(dip(3,:))
                dip=0.d0

                !ELECTRIC DIPOLE
                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)
                i=1
                do j=1,n
                    do k=1,j
                        if (j.eq.k) then
                            elec_dip1 = elec_dip1 + P1(i)*dip(:,i) !diagonal
                            elec_dip2 = elec_dip2 + P2(i)*dip(:,i)
                        else
                            elec_dip1 = elec_dip1 + P1(i)*dip(:,i)*2.d0 !off-diagonal
                            elec_dip2 = elec_dip2 + P2(i)*dip(:,i)*2.d0
                        endif
                        i=i+1
                    enddo
                end do

                !! CALCULATE SCALING FACTOR AND INDUCED DIPOLES
                scaled1=(fepsi/onsager_radius**3)*elec_dip1*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL
                scaled2=(fepsi/onsager_radius**3)*elec_dip2*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL
                scaled_nuc=(fepsi/onsager_radius**3)*nuc_dip*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL

                !!CALCULATE DERIVATIVES
                i3=1
                do k=1,numat
                    i1=qm2_params%orb_loc(1,k)
                    i2=qm2_params%orb_loc(2,k)
                    do i=i1,i2
                        do j=1,3
                            dxyz(j,k)=dxyz(j,k)+2*(p1(i3)*scaled2(j)+scaled1(j)*p2(i3)) !E-E parts
                            if (calc_nuc) dxyz(j,k)=dxyz(j,k)+2*scaled_nuc(j)*p2(i3) !N-E part
                        enddo
                        q_elec1=q_elec1+p1(i3)
                        q_elec2=q_elec2+p2(i3)
                        i3=i3+i+1 !diagonal indices
                    enddo
                    if(calc_nuc) dxyz(:,k)=dxyz(:,k)-2*qm2_params%core_chg(k)*scaled2 !N-E part
                enddo
                do j=1,3 !Translation
                    dxyz(j,:)=dxyz(j,:)-2*(q_elec1*scaled2(j)+scaled1(j)*q_elec2)/numat !E-E
                    if (calc_nuc) dxyz(j,:)=dxyz(j,:)+2*sum(qm2_params%core_chg)*scaled2(j)/numat&
                        -2*scaled_nuc(j)*q_elec2/numat !N-E
                enddo
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE POTENTIAL GRADIENT FOR GS, i.e. tr(F(rho)rho) INCLUDING
            !NUCLEAR-ELECTRON PART AND NUCLEAR-NUCLEAR GRADIENT
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfldgrad(dxyz,p,n)
                use qmmm_module,only: qm2_params,qmmm_struct;
                use cosmo_C, only: fepsi,onsager_radius,numat;
                use constants, only : one, BOHRS_TO_A, AU_TO_EV,EV_TO_KCAL

                implicit none;

                integer::  i,i1,i2,i3,j,k;
                integer::  n
                _REAL_ ::  dxyz(3,numat),origin(3,numat)
                _REAL_ ::  p(n*(n+1)/2)
                _REAL_ ::  dip(3,n*(n+1)/2);
                _REAL_ ::  scaled(3),elec_dip(3),nuc_dip(3),q_elec

                dxyz=0.d0; elec_dip=0.d0
                nuc_dip=0.d0;dip=0.d0;q_elec=0.d0

                call centercoords(origin)
                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                nuc_dip(1)=sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=sum(dip(2,:))
                nuc_dip(3)=sum(dip(3,:))
                dip=0.d0

                !ELECTRIC DIPOLE
                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)
                i=1
                do j=1,n
                    do k=1,j
                        if (j.eq.k) elec_dip = elec_dip + P(i)*dip(:,i) !diagonal
                        if (j.ne.k) elec_dip = elec_dip + P(i)*dip(:,i)*2.d0 !off-diagonal
                        i=i+1
                    enddo
                end do

                !! CALCULATE SCALING FACTOR
                scaled=(fepsi/onsager_radius**3)*(elec_dip+nuc_dip)*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL;
                !!CALCULATE DERIVATIVES
                i3=1
                do k=1,numat
                    i1=qm2_params%orb_loc(1,k)
                    i2=qm2_params%orb_loc(2,k)
                    do i=i1,i2
                        do j=1,3
                            dxyz(j,k)=dxyz(j,k)+2*p(i3)*scaled(j) !E-E,E-N parts
                        enddo
                        q_elec=q_elec+p(i3)
                        i3=i3+i+1 !diagonal indices
                    enddo
                    dxyz(:,k)=dxyz(:,k)-2*qm2_params%core_chg(k)*scaled !N-N,N-E parts
                enddo
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !GRADIENT FOR ONSAGER POTENTIAL FOR FULLL MATRICES, i.e. tr(V(xi)xi^(+))
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfldgrad_full(dxyz,p,n)
                use qmmm_module,only: qm2_params,qmmm_struct;
                use cosmo_C, only: fepsi,onsager_radius,numat;
                use constants, only : one, BOHRS_TO_A, AU_TO_EV,EV_TO_KCAL

                implicit none;

                integer::  i,i1,i2,j,k;
                integer::  n
                _REAL_ ::  dxyz(3,numat),origin(3,numat)
                _REAL_ ::  p(n**2)
                _REAL_ ::  dip(3,n*(n+1)/2),dip2(3,n**2);
                _REAL_ ::  scaled(3),elec_dip(3),nuc_dip(3),q_elec

                dxyz=0.d0; elec_dip=0.d0
                nuc_dip=0.d0;dip=0.d0; dip2=0.d0; q_elec=0.d0
                call centercoords(origin)

                !!GET ELEC DIP MATRIX
                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)
                call unpacking(n,dip(1,:),dip2(1,:),'s')
                call unpacking(n,dip(2,:),dip2(2,:),'s')
                call unpacking(n,dip(3,:),dip2(3,:),'s')

                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)

                do i=1,3
                    elec_dip(i) = sum(P*dip2(i,:)) !diagonal
                enddo

                !! CALCULATE SCALING FACTOR
                scaled=(fepsi/onsager_radius**3)*elec_dip*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL;
                !!CALCULATE DERIVATIVES
                do k=1,numat
                    i1=qm2_params%orb_loc(1,k)
                    i2=qm2_params%orb_loc(2,k)
                    do i=i1,i2;
                        do j=1,3
                            dxyz(j,k)=dxyz(j,k)+4*2*P((i-1)*n+i)*scaled(j)
                        enddo
                        q_elec=q_elec+P((i-1)*n+i)
                    enddo
                enddo
                do j=1,3 !Translation
                    dxyz(j,:)=dxyz(j,:)-2*q_elec/numat*scaled(j) !E-E
                enddo
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE POTENTIAL FOR GROUND STATE CALCULATION INCLUDING NUCLEAR-ELECTRONIC PART
            !THIS POTENTIAL IS ADDED TO THE FOCK OPERATOR
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfld_fock(f,p,n)
                use qmmm_module,only:qmmm_struct;
                use cosmo_C, only: numat,fepsi,onsagE,onsager_radius;
                use constants, only : BOHRS_TO_A, AU_TO_EV,CODATA08_AU_TO_DEBYE
                implicit none;
                _REAL_ f0((n+1)*n/2),f((n+1)*n/2),p((n+1)*n/2)
                _REAL_ dip(3,n*(n+1)/2);
                _REAL_ nuc_dip(3),scaled(3),elec_dip(3),Etest;
                integer n,i,j,k
                _REAL_ origin(3,numat)
                elec_dip=0.d0; dip=0.d0; Etest=0.d0

                call centercoords(origin) !origin to center of charge distribution

                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                nuc_dip(1)=sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=sum(dip(2,:))
                nuc_dip(3)=sum(dip(3,:))
                dip=0.d0
       
        
                !ELECTRIC DIPOLE
                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)
                i=1
                do j=1,n
                    do k=1,j
                        if (j.eq.k) elec_dip = elec_dip + P(i)*dip(:,i) !diagonal
                        if (j.ne.k) elec_dip = elec_dip + P(i)*dip(:,i)*2.d0 !off-diagonal
                        i=i+1
                    enddo
                end do
        
                !!CALCULATE SCALING FACTOR * INDUCED DIPOLE
                scaled=(fepsi/onsager_radius**3)*elec_dip*BOHRS_TO_A*AU_TO_EV

                !!CALCULATE REACTION FIELD POTENTIAL OPERATOR AND ADD TO FOCK MATRIX
                do i=1,(n+1)*n/2
                    f0(i)=2.d0*(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i));
                enddo
                f=f-f0
                onsagE=-sum((elec_dip+nuc_dip)**2)*fepsi/onsager_radius**3*BOHRS_TO_A*AU_TO_EV
                return
            end

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE POTENTIAL FOR GROUND STATE NUCLEAR-NUCLEAR PART
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfldnuc(enuclr)
                use qmmm_module,only:qm2_struct,qmmm_struct;
                use cosmo_C, only: numat,fepsi,onsager_radius
                use constants, only : one, BOHRS_TO_A, AU_TO_EV
                implicit none;
                _REAL_ dip(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)
                _REAL_ nuc_dip(3),scaled(3)
                _REAL_ enuclr,origin(3,numat)
                integer n
                n=qm2_struct%norbs
                dip=0.d0
        
                call centercoords(origin)
                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                nuc_dip(1)=sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=sum(dip(2,:))
                nuc_dip(3)=sum(dip(3,:))

                !!CALCULATE SCALING FACTOR * INDUCED DIPOLE
                scaled=(fepsi/onsager_radius**3)*nuc_dip*BOHRS_TO_A*AU_TO_EV

                !CALCULATE NUCLEAR-NUCLEAR + NUCLEAR-ELECTRONIC ENERGIES
                enuclr=enuclr-sum(scaled*nuc_dip)
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !ONSAGER-TYPE POTENTIAL FOR NUCLEAR-ELECTRONIC PART
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine rcnfldhcr(h)
                use qmmm_module,only:qm2_struct,qmmm_struct;
                use cosmo_C, only: numat,fepsi,onsager_radius
                use constants, only : one, BOHRS_TO_A, AU_TO_EV
                implicit none;
                _REAL_ dip(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)
                _REAL_ nuc_dip(3),scaled(3),elec_dip(3),Etest
                _REAL_ h(qm2_struct%norbs*(qm2_struct%norbs+1)/2)
                _REAL_ h0(qm2_struct%norbs*(qm2_struct%norbs+1)/2)
                _REAL_ origin(3,numat)
                integer n,i,j,k
                n=qm2_struct%norbs
                dip=0.d0
                Etest=0.d0

                call centercoords(origin)
                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                nuc_dip(1)=sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=sum(dip(2,:))
                nuc_dip(3)=sum(dip(3,:))
                dip=0.d0

                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip)
                i=1
                do j=1,n
                    do k=1,j
                        if (j.eq.k) elec_dip = elec_dip + qm2_struct%den_matrix(i)*dip(:,i) !diagonal
                        if (j.ne.k) elec_dip = elec_dip + qm2_struct%den_matrix(i)*dip(:,i)*2.d0 !off-diagonal
                        i=i+1
                    enddo
                end do

                !!CALCULATE SCALING FACTOR * INDUCED DIPOLE
                scaled=(fepsi/onsager_radius**3)*nuc_dip*BOHRS_TO_A*AU_TO_EV

                !CALCULATE NUCLEAR-ELECTRONIC and ELECTRONIC-NUCLEAR ENERGIES
                do i=1,(n+1)*n/2
                    h0(i)=-2*(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i));
                enddo
                h=h+h0
            end subroutine

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !SUBROUTINE FOR CONSTANT ELECTRIC FIELD SCREENING
            !OF ELECTRON-ELECTRON INTERACTION
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine efield_fock(f,n)
                use qmmm_module,only:qm2_struct,qmmm_struct;
                use cosmo_C, only: Ex,Ey,Ez
                use constants, only : one, BOHRS_TO_A, AU_TO_EV
                use qm2_davidson_module
                implicit none;
                integer n,k;
                _REAL_ f(n*(n+1)/2);
                _REAL_ tmp(n*(n+1)/2);
                _REAL_ dip(3,n*(n+1)/2);
                _REAL_ origin(3,qmmm_struct%nquant_nlink)
                _REAL_ ddot
                _REAL_ GSDM(n,n),mu_gr(3)
                tmp=0.d0
                dip=0.d0
                origin=0.d0;
                call centercoords(origin)
                !!GET ELEC DIP MATRIX
                call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip) !dipole in au * A
                call unpacking(qm2ds%Nb,qm2_struct%den_matrix,GSDM,'s')
                do k=1,3  ! loop over x,y,z
                    call unpacking(qm2ds%Nb,dip(k,:),qm2ds%eta_scratch,'s')
                    mu_gr(k)=ddot(qm2ds%Nb**2,GSDM,1,qm2ds%eta_scratch,1)/BOHRS_TO_A
                enddo
                !!CALCULATE ELECTRIC FIELD POTENTIAL OPERATOR
                tmp=Ex*dip(1,:)+Ey*dip(2,:)+Ez*dip(3,:);                !now in hartree
                f=f+tmp/BOHRS_TO_A
                return
            end

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !SUBROUTINE FOR CONSTANT ELECTRIC FIELD SCREENING
            !OF NUCLEAR-NUCLEAR INTERACTION
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine efield_nuc(enuclr)
                use qmmm_module,only:qm2_struct,qmmm_struct;
                use cosmo_C, only: Ex,Ey,Ez
                use constants, only : BOHRS_TO_A, AU_TO_EV
                implicit none;
                _REAL_ dip(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)
                _REAL_ nuc_dip(3)
                _REAL_ enuclr,origin(3,qmmm_struct%nquant_nlink)
                origin=0.d0; dip=0.d0
                call centercoords(origin)
                !!GET NUC DIP MAT and CALC NUC DIP
                call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
                write(6,*)'Nuclear Dipole:',sum(dip(1,:))/BOHRS_TO_A, &
                    sum(dip(2,:))/BOHRS_TO_A, &
                    sum(dip(3,:))/BOHRS_TO_A
                nuc_dip(1)=Ex*sum(dip(1,:))!NUCLEAR DIPOLE IN a.u. * angstrom
                nuc_dip(2)=Ey*sum(dip(2,:)) !Now in hartree
                nuc_dip(3)=Ez*sum(dip(3,:))
                !CALCULATE NUCLEAR-NUCLEAR ENERGY
                enuclr=enuclr+sum(nuc_dip)/BOHRS_TO_A!/2.d0
                return
            end


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !CALCULATE COORDINATE MOMENT
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            subroutine centercoords(origin)
                use qmmm_module,only: qmmm_struct
                use cosmo_C,only:numat
                implicit none;
                integer :: j
        
                _REAL_  :: origin(3,numat)
                origin=0.d0
                do j=1,3
                    origin(j,:)=origin(j,:)+sum(qmmm_struct%qm_coords(j,:))/numat
                enddo
            end subroutine

