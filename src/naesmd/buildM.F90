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

   _REAL_ u1(qm2ds%Nrpa),v1(qm2ds%Nrpa),v2(qm2ds%Nrpa),xi2(qm2ds%Nb*qm2ds%Nb)
   _REAL_ f,f1,f2,fs1,ddot,fmin
   integer i,j,p,h,k
   logical check_symmetry;
   logical flag_;
   integer one,junk
   integer solvent_model
   _REAL_  tmp(qm2ds%nb,qm2ds%nb),tmp2(qm2ds%nb,qm2ds%nb);
   _REAL_ WORK(qm2ds%Nb*3),dotprod1,dotprod2,testdip1,testdip2
   INTEGER LWORK,INFO

   parameter (one=1)

   fs1=0;
   call mo2site(u1,qm2ds%xi,qm2ds%eta) !Change basis of guess vector of Davidson from M.O to A.O
   qm2ds%eta=0.0;
   call Vxi(qm2ds%xi,qm2ds%eta); !Calculate Vacuum Electron Correlation
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
	tmp=2*tmp
        call VxiM_end(qm2ds%eta,tmp); !Add selected potential to vacuum correlation
   elseif (solvent_model.eq.99) then !For Z-vector equation if different
        tmp=0.d0;
        if (potential_type.eq.3) then !COSMO Potential
        call VxiM(qm2ds%xi,tmp);
        elseif (potential_type.eq.2) then !Onsager Potential
        call rcnfld(tmp,qm2ds%xi,qm2ds%nb)
        elseif (potential_type.eq.1) then !testing
        do i=1,qm2ds%nb; tmp(i,i)=qm2ds%eta(qm2ds%nb*(i-1)+i); enddo !double diag vac correlation
        endif
        tmp=2*tmp !linear response
        call VxiM_end(qm2ds%eta,tmp); !Add selected potential to vacuum correlation
   elseif(solvent_model.eq.2) then ! 2: State Specific [V_s(T+Z),xi]
        !Commutator is performed here for State Specific Solvent Routines

	!TESTING BLOCK
        !v_solvent_difdens=0.d0
        !i=0
        !do i=1,qm2ds%nb
	!	write(6,*)v_solvent_difdens(:,i)
            	!v_solvent_difdens(i,i)=2
            	!v_solvent_difdens(i*2,i*2)=-2
        !enddo; !i=0
	tmp=0.d0
        call commutator(v_solvent_difdens,qm2ds%xi,qm2ds%Nb,tmp,.false.)
	!tmp=-tmp/2 !test
        call VxiM_end(qm2ds%eta,tmp)
   elseif(solvent_model.eq.3) then !3: State Specific [V_s(xi),xi]
        call commutator(v_solvent_xi,qm2ds%xi,qm2ds%Nb,tmp,.false.)
        call VxiM_end(qm2ds%eta,tmp)
   elseif(solvent_model.eq.6) then!6: Solve nonlinear Liouville equation testing
        call commutator(qm2ds%eta,qm2ds%xi,qm2ds%Nb,tmp,.false.)
        call VxiM_end(qm2ds%eta,tmp)
   elseif(solvent_model.eq.10) then!10: NO GS Solvent test
        !call addnuc(tmp(1,1)); 
        tmp=0.d0; tmp2=0.d0;
        !call rcnfldhcr(tmp)
        !write(6,*)'hcr',tmp
        call rcnfld_fock(tmp,qm2_struct%den_matrix,qm2ds%Nb)
        !write(6,*)'solvfock=',tmp
        !stop
        call unpacking(qm2ds%Nb,tmp,tmp2,'s'); tmp=0.d0;
        call commutator(tmp2,qm2ds%xi,qm2ds%Nb,tmp,.false.)
        call VxiM_end(qm2ds%eta,tmp)
   endif
   
   !add constant electric field potential [V_EF,xi]
   if(EF.eq.2) then!Constant Electric Field in ES onl
        !write(6,*)'Using Constant EF in Excited State only'
        tmp=0.d0; tmp2=0.d0
        call efield_fock(tmp,qm2ds%Nb)
	call unpacking(qm2ds%Nb,tmp,tmp2,'s'); tmp=0.d0
        call commutator(qm2ds%xi,tmp2,qm2ds%Nb,tmp,.false.)
	!write(6,*)'Efield Operator:'
	!do p=1,qm2ds%Nb
	!	write(6,*)tmp2(:,p)
	!	do h=1,qm2ds%Nb
	!		tmp2(p,h)=(tmp2(p,p)-tmp2(h,h))*qm2ds%xi((p-1)*qm2ds%Nb+h)
	!	enddo
	!enddo
        call VxiM_end(qm2ds%eta,tmp)
   endif
   call site2mo(qm2ds%xi,qm2ds%eta,v1); !Change basis of xi again to M.O.

   i=0
   do p=1,qm2ds%Np
    do h=qm2ds%Np+1,qm2ds%Nb
       i=i+1
       f=qm2ds%ehf(h)-qm2ds%ehf(p);
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
	use cosmo_C,only:fepsi,nps,lm61,numat,a0,ev,amat,bmat,ipiden,gden,nsetf
	implicit none
	_REAL_, intent(in)::xi(qm2ds%nb,qm2ds%nb)
	_REAL_, intent(inout)::Vximmat(qm2ds%nb,qm2ds%nb)
	_REAL_ :: density(lm61),charges(nps),phi(nps)
  	_REAL_ :: p(qm2ds%nb*(qm2ds%nb+1)/2),fcon,INFO
  	integer:: i,im,j

  	fcon=fepsi*a0*ev !scaling factor
  	call packing(qm2ds%nb,xi,p,'s') !Note that the factor of two for diagonal elements which happens in this subroutine is corrected for by gden
  	! FIRST CALCULATE QDENEL FROM DENSITY MATRIX
  	! gden is -2 for diagonal and -1 for nondiagonal
  	do i=1,lm61
      		density(i)=gden(i)*p(ipiden(i))
  	end do
  	!  NOW CALCULATE PHIEL FROM BMAT*QDENEL
	if (1==0) then
  		phi=0.d0; p=0.d0;
  		do i=1,nps ! running over SAS tiles
    			do j=1,lm61
      				phi(i)=phi(i)+bmat(j,i)*density(j)
    			end do
  		end do
	else
  		call  DGEMV ( 'T', lm61, nps, 1.d0, bmat, lm61, density, 1, 0.d0, phi, 1 ) !LAPACK subroutine for above
	endif
  	! NOW CALCULATE CHARGES
  	call coscl2(amat,nsetf,charges,phi,nps) !Use cholesky factorization routine
 
 	! NOW ADD BMAT*QSCEL TO FOCK MATRIX
	if (1==0) then
  		do i=1,lm61
    			im=ipiden(i)
    			do j=1,nps
      				p(im)=p(im)+fcon*bmat(i,j)*charges(j)
    			enddo
  		enddo
	else
  		call  DGEMV ( 'N', lm61, nps, fcon, bmat, lm61, charges, 1, 0.d0, density, 1 ) !LAPACK subroutine for above
  		p=0.d0
  		do i=1,lm61;
			im=ipiden(i); p(im)=density(i); !in A.U.
		 enddo
	endif
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
  integer matrix_size,i,j,k;
	logical flag;
	real(8), dimension(matrix_size,matrix_size) :: A,B,Res;
	real(8), dimension(:,:), allocatable :: tmp,tmp2
	real(8) ALPHA;
	real(8) BETA;

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
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: fepsi,onsagE,onsager_radius;
        use constants, only : one, BOHRS_TO_A, AU_TO_EV
        implicit none;
        integer n;
        _REAL_ f(n,n);
        _REAL_ p(n,n);
        _REAL_ tmp(n,n);
        _REAL_ dipx(n,n),dipy(n,n),dipz(n,n);
        _REAL_ dip(3,n*(n+1)/2);
        _REAL_ scaled(3),elec_dip(3);
        integer i,j,k;
        logical check_symmetry
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
        scaled=(fepsi/onsager_radius**3)*elec_dip*BOHRS_TO_A*AU_TO_EV; ! for units of eV

        !!CALCULATE REACTION FIELD POTENTIAL
        f=f-2.d0*(scaled(1)*dipx+scaled(2)*dipy+scaled(3)*dipz); !in eV
        !f=f-(scaled(1)*dipx+scaled(2)*dipy+scaled(3)*dipz); !in eV test

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

        !call unpacking(n,p,p_full,'s') !Surely this can be done without unpacking
        !!CALCULATE DERIVATIVES
        i3=1
        do k=1,numat
                i1=qm2_params%orb_loc(1,k)
                i2=qm2_params%orb_loc(2,k)
        !        p_full(i1,i1)=p_full(i1,i1)-qm2_params%core_chg(k)
        !        do i=i1,i2;
        !               do j=1,3 
        !               dxyz(j,k)=dxyz(j,k)+2*p_full(i,i)*scaled(j)
        !               enddo
        !        enddo
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
              !dxyz(:,k)=dxyz(:,k)-2*qm2_params%core_chg(k)*scaled_nuc !N-N part
        enddo
           do j=1,3 !Translation
              dxyz(j,:)=dxyz(j,:)-2*(q_elec1*scaled2(j)+scaled1(j)*q_elec2)/numat !E-E
              if (calc_nuc) dxyz(j,:)=dxyz(j,:)+2*sum(qm2_params%core_chg)*scaled2(j)/numat&
                           -2*scaled_nuc(j)*q_elec2/numat !N-E
              !dxyz(j,:)=dxyz(j,:)-2*sum(qm2_params%core_chg)*scaled_nuc(j)/numat !N-N
           enddo
        !write(6,*)'nuc_dip',nuc_dip
        !write(6,*)'elec_dip 2',elec_dip2/BOHRS_TO_A
        !write(6,*)'sum_dip',nuc_dip+elec_dip1+elec_dip2/BOHRS_TO_A
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
        !scaled=(fepsi/onsager_radius**3)*elec_dip*BOHRS_TO_A*AU_TO_EV*EV_TO_KCAL; ! Testing

        !call unpacking(n,p,p_full,'s') !Surely this can be done without unpacking
        !!CALCULATE DERIVATIVES
        i3=1
        do k=1,numat 
                i1=qm2_params%orb_loc(1,k)
                i2=qm2_params%orb_loc(2,k)
        !        p_full(i1,i1)=p_full(i1,i1)-qm2_params%core_chg(k)
        !        do i=i1,i2;
        !               do j=1,3 
        !               dxyz(j,k)=dxyz(j,k)+2*p_full(i,i)*scaled(j)
        !               enddo
        !        enddo
        do i=i1,i2
           do j=1,3
              dxyz(j,k)=dxyz(j,k)+2*p(i3)*scaled(j) !E-E,E-N parts
           enddo
           q_elec=q_elec+p(i3)
           i3=i3+i+1 !diagonal indices
        enddo
        dxyz(:,k)=dxyz(:,k)-2*qm2_params%core_chg(k)*scaled !N-N,N-E parts
        enddo
           do j=1,3 !Translation
              !dxyz(j,:)=dxyz(j,:)-2*q_elec/numat*scaled(j) !E-E
              !dxyz(j,:)=dxyz(j,:)-(sum(qm2_params%core_chg)+q_elec)*scaled(j)/numat !E-N + N-E
              !dxyz(j,:)=dxyz(j,:)+2*sum(qm2_params%core_chg)/numat*scaled(j) !N-N
           enddo
        !write(6,*)'nuc_dip',nuc_dip
        !write(6,*)'elec_dip',elec_dip
        !write(6,*)'sum_dip',nuc_dip+elec_dip
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!GRADIENT FOR ONSAGER POTENTIAL FOR FULLL MATRICES, i.e. tr(V(xi)xi^(+))
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine rcnfldgrad_full(dxyz,p,n)
        use qmmm_module,only: qm2_params,qmmm_struct;
        use cosmo_C, only: fepsi,onsager_radius,numat;
        use constants, only : one, BOHRS_TO_A, AU_TO_EV,EV_TO_KCAL

        implicit none;

        integer::  i,i1,i2,i3,j,k;
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
        !write(6,*)'Dipole:',elec_dip/BOHRS_TO_A*sqrt(2.d0)
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
        use qmmm_module,only:qm2_struct,qmmm_struct;
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
        !f=qm2_struct%hmatrix !test
        do i=1,(n+1)*n/2
           f0(i)=2.d0*(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i));
           !f0(i)=(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i)); !test

           !Etest=Etest+f0(i)*P(i)
        enddo
           f=f-f0
        !onsagE=sum((elec_dip+nuc_dip)**2)*fepsi/onsager_radius**3*BOHRS_TO_A*AU_TO_EV/2
        onsagE=-sum((elec_dip+nuc_dip)**2)*fepsi/onsager_radius**3*BOHRS_TO_A*AU_TO_EV
        !WRITE RESUTLS FOR TESTING
        !write(6,*)'fs',scaled(1)*dipx+scaled(2)*dipy+scaled(3)*dipz
        !write(6,*)'N-N,E-E,N-E,total',(fepsi/onsager_radius**3)*BOHRS_TO_A*AU_TO_EV*sum(nuc_dip**2),&
        !                        (fepsi/onsager_radius**3)*BOHRS_TO_A*AU_TO_EV*sum(elec_dip**2),&
        !                        (fepsi/onsager_radius**3)*BOHRS_TO_A*AU_TO_EV*sum(elec_dip*nuc_dip),&
        !                        onsagE
        !write(6,*)'totale',sum(scaled*elec_dip)+(fepsi/onsager_radius**3)*sum(nuc_dip**2)*BOHRS_TO_A*AU_TO_EV
        !write(6,*)'totale',sum((elec_dip+nuc_dip)**2)*(fepsi/onsager_radius**3)*BOHRS_TO_A*AU_TO_EV
        !write(6,*)'E-E=',Etest,'E-E2=',sum(scaled*elec_dip)
        !write(6,*)'E-N+N-N',sum(scaled*nuc_dip),OnsagE
        !write(6,*)'nuc:',nuc_dip
        !write(6,*)'elec:',elec_dip
        !write(6,*)'sum:',(nuc_dip+elec_dip)
        !write(6,*)'den_matrix=',qm2_struct%den_matrix
return
end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!ONSAGER-TYPE POTENTIAL FOR GROUND STATE NUCLEAR-NUCLEAR PART
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine rcnfldnuc(enuclr)
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: numat,fepsi,onsager_radius,tri_2d,onsagE
        use constants, only : one, BOHRS_TO_A, AU_TO_EV
        implicit none;
        _REAL_ dip(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)
        _REAL_ nuc_dip(3),scaled(3),elec_dip(3)
        _REAL_ enuclr,origin(3,numat)
        integer n,i
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
        !enuclr=enuclr-sum(scaled*nuc_dip)/2.d0 !test

       
        !WRITE RESULTS AND TEST
        !write(6,*)'N-N=',sum(scaled*nuc_dip)
        !enuclr=-(fepsi/onsager_radius**3)*sum(elec_dip**2)*BOHRS_TO_A*AU_TO_EV
        !write(6,*)'addnuc: E-N + N-N',enuclr
        !write(6,*)'nuc_dip:',nuc_dip
        !write(6,*)'elec_dip:',elec_dip
        !write(6,*)'sum_dip:',elec_dip+nuc_dip
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!ONSAGER-TYPE POTENTIAL FOR NUCLEAR-ELECTRONIC PART
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine rcnfldhcr(h)
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: numat,fepsi,onsager_radius,tri_2d,onsagE
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
        !h=0.d0 !test
        do i=1,(n+1)*n/2
           h0(i)=-2*(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i));
           !h0(i)=-(scaled(1)*dip(1,i)+scaled(2)*dip(2,i)+scaled(3)*dip(3,i)); !test

        enddo
         h=h+h0
        !TEST
        !write(6,*)elec_dip
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!SUBROUTINE FOR CONSTANT ELECTRIC FIELD SCREENING 
!OF ELECTRON-ELECTRON INTERACTION
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine efield_fock(f,n)
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: fepsi,Ex,Ey,Ez
        use constants, only : one, BOHRS_TO_A, AU_TO_EV
        implicit none;
        integer n,j,k;
        _REAL_ f(n*(n+1)/2);
        _REAL_ tmp(n*(n+1)/2);
        _REAL_ dip(3,n*(n+1)/2);
        _REAL_ efv(3)
	_REAL_ origin(3,qmmm_struct%nquant_nlink)
        tmp=0.d0
        dip=0.d0
	origin=0.d0;
        call centercoords(origin)
        !!GET ELEC DIP MATRIX
        call get_dipole_matrix(qmmm_struct%qm_coords-origin, dip) !dipole in au * A
        !!CALCULATE ELECTRIC FIELD POTENTIAL OPERATOR
        tmp=Ex*dip(1,:)+Ey*dip(2,:)+Ez*dip(3,:); !now in hartree
	!write(6,*)'Efield Operator in Subroutine'
	!do j=1,n
	!	write(6,*)tmp(:,j)
	!enddo
        f=f+tmp/BOHRS_TO_A
return
end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!SUBROUTINE FOR CONSTANT ELECTRIC FIELD SCREENING 
!OF NUCLEAR-NUCLEAR INTERACTION
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine efield_nuc(enuclr)
        use qmmm_module,only:qm2_struct,qmmm_struct;
        use cosmo_C, only: fepsi,Ex,Ey,Ez
        use constants, only : BOHRS_TO_A, AU_TO_EV
        implicit none;
        _REAL_ dip(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)
        _REAL_ nuc_dip(3)
        _REAL_ enuclr,origin(3,qmmm_struct%nquant_nlink)
	origin=0.d0; dip=0.d0
        call centercoords(origin)
        !!GET NUC DIP MAT and CALC NUC DIP
        call get_nuc_dip(qmmm_struct%qm_coords-origin, dip);
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
        use qmmm_module,only: qmmm_struct,qm2_struct,qm2_params
        use cosmo_C,only:numat
        implicit none;
        integer :: i,i1,i2,i3,j,k
        
        _REAL_  :: origin(3,numat)
           origin=0.d0
           do j=1,3
              origin(j,:)=origin(j,:)+sum(qmmm_struct%qm_coords(j,:))/numat
           enddo
           !write(6,*)'origin',origin
end subroutine


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!The following functions aren't called by the program but are included here 
!testing purposes.
!
! Josiah A. Bjorgaard, Vasyl Kuzmenko, Kirill Velizhanin 
! 2013-2014 Los Alamos National Laboratory
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!THIS SUBROUTINE CALCULATES THE SOLVENT POTENTIAL OPERATOR V_s FOR COSMO DIRECTLY USING MATRIX INVERSION
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine VxiM_full(xi,Vximmat)
  use constants, only : half,AU_TO_EV,BOHRS_TO_A;
  use qm2_davidson_module;
  use qmmm_module,only:qm2_struct;
  use cosmo_C,only:ipiden,lm61,tri_2D,mmat;
  implicit none
  _REAL_,intent(inout)::xi(qm2_struct%norbs,qm2_struct%norbs);
  _REAL_,intent(inout)::Vximmat(qm2_struct%norbs,qm2_struct%norbs);
  _REAL_ e0;
  integer i,j,m,n,k,l,one;
  _REAL_::mmat_1,mmat_2;


        do i=1,lm61
                do j=1,i
                        k=tri_2D(1,ipiden(j)); l=tri_2D(2,ipiden(j));
                        m=tri_2D(1,ipiden(i)); n=tri_2D(2,ipiden(i));

      mmat_1=mmat(i,j);

      if (m.eq.n) mmat_1=mmat_1/2.d0;
      if (k.eq.l) mmat_1=mmat_1/2.d0;
      if (i.eq.j) mmat_1=mmat_1/2.d0;
      
      mmat_2=mmat_1*1.d0;
      !mmat_2=mmat_1*2.d0     
        
                        Vximmat(m,n)=Vximmat(m,n)+mmat_2*xi(l,k);
                        Vximmat(m,n)=Vximmat(m,n)+mmat_2*xi(k,l);
                        Vximmat(n,m)=Vximmat(n,m)+mmat_2*xi(l,k);
                        Vximmat(n,m)=Vximmat(n,m)+mmat_2*xi(k,l);

                        Vximmat(l,k)=Vximmat(l,k)+mmat_2*xi(m,n);
                        Vximmat(l,k)=Vximmat(l,k)+mmat_2*xi(n,m);
                        Vximmat(k,l)=Vximmat(k,l)+mmat_2*xi(m,n);
                        Vximmat(k,l)=Vximmat(k,l)+mmat_2*xi(n,m);

      !exchange       

                        !Vximmat(m,l)=Vximmat(m,l)-mmat_1*xi(n,k);
                        !Vximmat(m,k)=Vximmat(m,k)-mmat_1*xi(n,l);
                        !Vximmat(n,l)=Vximmat(n,l)-mmat_1*xi(m,k);
                        !Vximmat(n,k)=Vximmat(n,k)-mmat_1*xi(m,l);

                        !Vximmat(k,n)=Vximmat(k,n)-mmat_1*xi(l,m);
                        !Vximmat(k,m)=Vximmat(k,m)-mmat_1*xi(l,n);
                        !Vximmat(l,n)=Vximmat(l,n)-mmat_1*xi(k,m);
                        !Vximmat(l,m)=Vximmat(l,m)-mmat_1*xi(k,n);

                end do
        end do
write(6,*)'Vximmat',Vximmat
stop
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!CALCULATE M=B^t A^(-1) B in units of eV and scaled by dielectric factor by brute force matrix inversion
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine buildM_testing()
	use cosmo_C,only:nps,lm61,amat,bmat,nsetf,mmat,fepsi,qdenet,a0,ev
	integer i;
	integer j;
	real(8) scalarM;
	_REAL_, dimension (:,:), allocatable::Atmp;
	_REAL_, dimension (:,:), allocatable::A;
	_REAL_, dimension (:,:), allocatable::B;
	_REAL_, dimension (:,:), allocatable::tmp;
	integer INFO;
	integer IPIV(nps);
	real(8) ALPHA;
	real(8) BETA;
	integer LWORK;
	!integer M;
        !integer LDB; 
        !integer NRHS;
        !integer KD;
        !integer LDAB;
        

        write(6,*) 'BuildM_testing called' !!JAB Testing
	!write(6,*)'amat:',amat
        !write(6,*)'bmat:',bmat
	LWORK =nps;
	ALPHA = 1.0;
	BETA = 0.0;

	allocate(A(nps,nps));
	do i=1,nps
		do j = 1, i
			A(i,j) = amat((i-1)*i/2+j);
			A(j,i) = A(i,j);
		end do
	end do
	allocate(Atmp(nps,nps));
	Atmp(1:nps,1:nps)=0.0;
	do i=1,nps
               Atmp(i,i)=1.0; ! Unitary Diagonal matix 
	end do
	CALL DGESV( nps, nps, A, nps, IPIV, Atmp, nps, INFO);
	deallocate(A)
	ALPHA=1.0;
	BETA=0.0;
        allocate(B(lm61,nps));
	B(1:lm61,1:nps)=bmat(1:lm61,1:nps);
	allocate(tmp(nps,lm61));
	CALL DGEMM('N','T',nps, lm61, nps, ALPHA, Atmp, nps, B, lm61, BETA, tmp, nps );
    	deallocate(Atmp);
	scalarM=-fepsi*a0*ev;
	B(1:lm61,1:nps)= scalarM * bmat(1:lm61,1:nps);
        CALL DGEMM('N','N',lm61, lm61, nps, ALPHA, B, lm61, tmp, nps, BETA, mmat, lm61 );
	deallocate(B,tmp);
end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!return indexes of full F/P/KSI/E <==what is this?
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine init_tri_2D(tri_2D,T2DSize)
        use qmmm_module,only:qm2_struct;
        implicit none;
        integer  T2DSize;
        integer, dimension(4,T2DSize), intent(inout) :: tri_2D;
        integer, parameter :: im = 1;
        integer, parameter :: in = 2;
        integer, parameter :: imu = 3;
        integer, parameter :: inu = 4;
        integer m, n, i;
        integer Nb;
        Nb = qm2_struct%norbs;
	i=0;
	do m=1,Nb
                do n=1,m
                        i=i+1;
                        tri_2D(im,i) = m;
                        tri_2D(in,i) = n;
                end do
        end do
	i =0;
	do m=1,Nb
                do n=m,Nb
                        i=i+1;
                        tri_2D(imu,i) = m;
                        tri_2D(inu,i) = n;
                end do
        end do
end subroutine init_tri_2D

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!return orbital indexes for NDDO style charges
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine index_lm61(i_in,row_out,column_out)
	use constants, only : half;
	use cosmo_C,only:ipiden,lm61;
	integer, intent(in)::i_in;
	integer, intent(out)::row_out;
	integer, intent(out)::column_out;
	integer alpha;
	real(8) a;
	alpha = ipiden(i_in);
	!solving alpha=a(a+1)/2
	a=-half*(1-sqrt(1.+8.*alpha));
	row_out=CEILING(a);
	!if ((k_out*(k_out+1)/2>alpha) .and.(k_out/(k_out+1)/2.<))
	column_out=row_out;

	!	search for column
	!	i-row,j-column
	!	k -triangular index
	!	while (k=(i-1)*i/2+j > our 1D array index ) then j=j-1;	

	do while (((row_out-1)*row_out/2+column_out) >ipiden(i_in)) 
		column_out=column_out-1;
	end do
end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Check matrix symmetry
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
logical function check_symmetry(M,matrix_size)
        implicit none;
        integer matrix_size;
        _REAL_ M(matrix_size,matrix_size)
        _REAL_ sum,tmp_scalar,norm;
        integer i;
        integer j;
        check_symmetry=.true.;
        sum=0.d0;
        do i=1,matrix_size
          do j=1,i-1
                tmp_scalar=M(i,j)-M(j,i);
               	sum=sum+tmp_scalar**2.;
          end do
       	end do
	norm= sqrt(sum)
        write(6,*)"Norm:",norm
        if (norm> 1.0d-10) check_symmetry=.false.;
        end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Calculate COSMO Gradient Brute Force Test (Broken)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine cosgrad(GradE,rho1,rho2)

use cosmo_C, only:nps,lm61,iatsp,amat,bmat,cosurf,idenat,ipiden,a0,ev,numat,tri_2D
use qmmm_module
implicit none

integer         ::      n,i,j,iatom,jatom,im1,im2,jm1,jm2,iin,ifi,jin,jfi,INFO
integer         ::      AOtolm61(qm2_struct%norbs**2,qm2_struct%norbs**2)
_REAL_          ::      GradA(3,nps,nps),GradB(3,lm61,nps),GradM(3,lm61,lm61),&
                        AinvBtGradA(3,lm61,nps)
_REAL_          ::      A(nps,nps),Ainv(nps,nps),AinvB(nps,lm61),IPIV(nps)
_REAL_          ::      GradE(3,qmmm_struct%nquant_nlink),ALPHA,BETA
_REAL_,intent(in)::     rho1(qm2_struct%norbs,qm2_struct%norbs),&
                        rho2(qm2_struct%norbs,qm2_struct%norbs)


!convert lm61 indices to AO indices
do i=1,lm61
       AOtolm61(tri_2D(1,ipiden(i)),tri_2D(2,ipiden(i)))=i;
       AOtolm61(tri_2D(2,ipiden(i)),tri_2D(1,ipiden(i)))=i;
enddo

!Calculate Grad(A)
do i=1,nps
    iatom=iatsp(i)
    do j=1,i
        jatom=iatsp(j)
        if (iatom.ne.jatom) then
        write(6,*)'AMAT:',amat((i-1)*i/2+j),1/sqrt(sum((cosurf(:,i)-cosurf(:,j))**2))
        do n=1,3
        GradA(n,i,j)=abs(amat((i-1)*i/2+j))**(3)*(cosurf(n,i)-cosurf(n,j))
        GradA(n,j,i)=-GradA(n,i,j)
        A(i,j)=amat((i-1)*i/2)
        A(j,i)=A(i,j)
        enddo
        endif
    enddo
enddo

!Calculate Grad(B)
do i= 1,nps
     iatom=iatsp(i)
     do jatom=1,numat
        if (iatom.ne.jatom) then
          jin=qm2_params%orb_loc(1,jatom)
          jfi=qm2_params%orb_loc(2,jatom)
             do jm1=jin,jfi
             do jm2=jin,jm1
             !write(6,*)i,AOtolm61(jm1,jm2),size(bmat)
   write(6,*)'BMAT:',bmat(AOtolm61(jm1,jm2),i),1/sqrt(sum((cosurf(:,i)-qmmm_struct%qm_coords(:,jatom))**2)),i,AOtolm61(jm1,jm2)
                do n=1,3
                GradB(n,i,AOtolm61(jm1,jm2))=abs(bmat(i,AOtolm61(jm1,jm2)))**(3)*&
                        (cosurf(n,i)-qmmm_struct%qm_coords(n,jatom))
                enddo
             enddo
             enddo
        endif
     enddo
enddo
stop
!Calculate A^(-1)
Ainv(1:nps,1:nps)=0.0;
     do i=1,nps
           Ainv(i,i)=1.0; ! Unitary Diagonal matix 
     end do
CALL DGESV( nps, nps, A, nps, IPIV, Ainv, nps, INFO);

!Calculate Grad(M)
! GradB*AinvB + AinvB^t*GradB + AinvB^t*GradA*AinvB
BETA=0.0; ALPHA=1.0;
CALL DGEMM('N','N',nps, nps, lm61, ALPHA, Ainv, nps, bmat, lm61,BETA, AinvB, lm61);

do n=1,3
BETA=0.0
CALL DGEMM('T','N',lm61, nps, nps, ALPHA, AinvB, lm61,GradA(n,:,:),nps,BETA,AinvBtGradA(n,:,:), lm61 );
write(6,*)AinvBtGradA
stop
CALL DGEMM('T','N',lm61, nps, lm61, ALPHA, GradB(n,:,:), lm61, AinvB, lm61, BETA,GradM(n,:,:),lm61 );
!Transpose GradB? Equiv of first and second terms?
BETA=1.0;
CALL DGEMM('T','N',lm61, nps, lm61, ALPHA, AinvB, lm61, GradB(n,:,:), lm61, BETA,GradM(n,:,:), lm61 );
CALL DGEMM('N','N',lm61, nps, lm61, ALPHA, AinvBtGradA(n,:,:), lm61, AinvB,lm61,BETA,GradM(n,:,:), lm61 );
enddo
!Note: GradM should be symmetric

!Calculate  <rho1|GradM|rho2>=Tr(V_s(rho2)*rho1^t) for each atom
do i=1,numat
          iin=qm2_params%orb_loc(1,i)
          ifi=qm2_params%orb_loc(2,i)
     do j=1,numat
          if(i.ne.j) then
          jin=qm2_params%orb_loc(1,j)
          jfi=qm2_params%orb_loc(2,j)
          do im1=iin,ifi
          do im2=iin,ifi
             do jm1=jin,jfi
             do jm2=jin,jfi
                do n=1,3
                GradE(n,i)=GradE(n,i)+rho1(im1,im2)*rho2(jm1,jm2)*&
                                GradM(n,AOtolm61(im1,im2),AOtolm61(jm1,jm2))
                enddo
             enddo
             enddo
          enddo
          enddo
          endif
     enddo
enddo
GradE=GradE*a0*ev !Gradient in eV/A or eV/Bohr??

return
end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!A GENERAL SUBROUTINE FOR CALCULATING THE SOLVENT ENERGY AS Tr(V_s(rho_1)*rho_2)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
subroutine calc_Esolv(Esolv,sol_dens,dens,docore)
        use qmmm_module,only:qm2_struct
        use cosmo_C
        implicit none
        logical        ::      docore
        integer        ::      n,m,p
        _REAL_         ::      Esolv,Enuclr
        _REAL_         ::      sol_dens(qm2_struct%norbs**2)
        _REAL_         ::      sol_dens_p(qm2_struct%norbs*(qm2_struct%norbs+1)/2)
        _REAL_         ::      dens(qm2_struct%norbs**2)
        _REAL_         ::      V_s(qm2_struct%norbs**2)
        _REAL_         ::      V_s_p(qm2_struct%norbs*(qm2_struct%norbs+1)/2)

        Esolv=0.d0; Enuclr=0.d0
        V_s=0.d0; V_s_p=0.d0

!calculate potential with or without core charges
if (.not.docore) then !Calculate potential without core charges as in ES

        call VxiM(sol_dens,V_s) 


else !Calculate potential with core charges as in GS
        call packing(qm2_struct%norbs,sol_dens,sol_dens_p, 's')
        call addnuc(Enuclr)
        call addhcr(V_s_p); V_s_p=V_s_p*2.d0 !add 1 e- interactions
        call addfck(V_s_p,sol_dens_p) !add 2 e- interactions
        call unpacking(qm2_struct%norbs,V_s_p,V_s,'s')
endif

!calculate tr(VxiM*dens)
do n=1,qm2_struct%norbs
       do m=1,qm2_struct%norbs
              p=(n-1)*qm2_struct%norbs+m
              Esolv=Esolv+V_s(p)*dens(p)
              !write(6,*)V_s(p)*dens(p),V_s(p),dens(p)
       enddo
enddo
Esolv=Esolv/2.d0
write(6,*)'Solvent Energies (Electronic/Nuclear/Total)',Esolv,Enuclr,Esolv+Enuclr

!Note that addfck and VxiM give identical results if solute charges due to cores
!are eliminated in addfck
end subroutine


