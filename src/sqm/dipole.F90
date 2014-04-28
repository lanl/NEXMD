! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

!*************************************************
! Produces the dipole operator matrix for semiemirical methods
!*****************************************************
subroutine get_dipole_matrix(coord, dpmat)
 	use qmmm_module, only : qm2_params, qm2_struct, qmmm_struct, qmmm_nml
 	use constants, only : BOHRS_TO_A
	! use findmask
      	implicit none

	! Calculates the dipole matrices
      	_REAL_, intent(inout) :: coord(*)
	! CML TEST 7/17/12
      	_REAL_, intent(inout) :: dpmat(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2) 
      	integer :: loop_count,orb_beg,orb_end,orb_size,i,nquant

      	! AWG: Currently not supported with d orbitals for NDDO methods
      	if ( .not. qmmm_nml%qmtheory%DFTB ) then
         	do i = 1,qmmm_struct%nquant_nlink
            		if ( qm2_params%natomic_orbs(i) > 5 ) return
         	end do
      	end if

	dpmat = 0.d0

      	! Increasing the number of quantum atoms and determining orbital numbers
      	do nquant=1,qmmm_struct%nquant_nlink
         	! Asking for the Method, if DFTB is used it does not compute the changes in charge due to orbitals.             
         	if (.not. qmmm_nml%qmtheory%DFTB) then
            		orb_beg=qm2_params%orb_loc(1,nquant)
            		orb_end=qm2_params%orb_loc(2,nquant)

            		! Determination of the values difference between the orbitals
            		orb_size=orb_end-orb_beg

			! Fill diagonal elements of dipole matrices with coordinate * electron charge (-1)
			do loop_count = orb_beg, orb_end
				do i = 1,3
	           		  dpmat(i,qm2_params%pascal_tri2(loop_count)) = -coord((nquant-1)*3+i)
				end do
			end do

			! Fill the off-diagonal elements for s-p interactions on the same atom
            		if (orb_size>0) then 
               			i=0
               			do loop_count=orb_beg+1,orb_beg+3 ! only for p-orbitals
                  			i=i+1 ! X, Y, Z
				 	dpmat(i,qm2_params%pascal_tri1(loop_count)+orb_beg) = &
				  		-qm2_params%multip_2c_elec_params(1,nquant)*BOHRS_TO_A  !//CODATA08_A_TO_BOHRS ! Convert to electron-Angstroms
               			end do
            		end if 
         	end if
      	end do	
      	
	
      	return
endsubroutine get_dipole_matrix

!******************************************************************
! Calculates the transition dipole
!******************************************************************

subroutine trans_dipole(mu, alpha)
	use qm2_davidson_module
	use qmmm_module, only: qmmm_struct
	use constants, only : BOHRS_TO_A, SQRT2
	implicit none
	
	_REAL_, intent(out) :: mu(3, qm2ds%Mx)		! Dipole moment
	_REAL_, intent(out) :: alpha(3)		! Polarizability

	_REAL_ :: dip(3,qm2ds%Nb*(qm2ds%Nb+1)/2) ! Dipole matrix elements in 3D
	_REAL_ :: ddot ! This is a function
	integer :: j, one, Lt

	Lt = qm2ds%Nb*(qm2ds%Nb+1)/2
	one = 1
	mu = 0.d0
	alpha = 0.d0

	call get_dipole_matrix(qmmm_struct%qm_coords, dip)

	do j = 1,qm2ds%Mx	
		call mo2site(qm2ds%v0(1,j), qm2ds%xi_scratch, qm2ds%eta_scratch)
		call unpacking(qm2ds%Nb,dip(1,:),qm2ds%eta_scratch,'s')
		mu(1,j) = ddot(qm2ds%Nb**2,qm2ds%xi_scratch,one,qm2ds%eta_scratch,one)*SQRT2
		call unpacking(qm2ds%Nb,dip(2,:),qm2ds%eta_scratch,'s')
		mu(2,j) = ddot(qm2ds%Nb**2,qm2ds%xi_scratch,one,qm2ds%eta_scratch,one)*SQRT2
		call unpacking(qm2ds%Nb,dip(3,:),qm2ds%eta_scratch,'s')
		mu(3,j) = ddot(qm2ds%Nb**2,qm2ds%xi_scratch,one,qm2ds%eta_scratch,one)*SQRT2
		alpha(1) = alpha(1)+2*mu(1,j)**2/qm2ds%e0(j)
		alpha(2) = alpha(2)+2*mu(2,j)**2/qm2ds%e0(j)
		alpha(3) = alpha(3)+2*mu(3,j)**2/qm2ds%e0(j)
	end do

	mu(1:3,1:qm2ds%Mx)=mu(1:3,1:qm2ds%Mx)/BOHRS_TO_A ! convert to BOHRS

	return

end subroutine trans_dipole

!*************************************************************************
!Calcualtes the ground and excited state dipoles
!*************************************************************************
subroutine qm2_calc_molecular_dipole_in_excited_state()
	use qm2_davidson_module
	use qmmm_module, only: qmmm_struct, qm2_struct, qm2_params, qmmm_nml
	use constants, only : light_speed, charge_on_elec
	implicit none
	_REAL_ :: mu(3,qm2ds%Mx),mu_relaxed(3,qm2ds%Mx),mu_unrelaxed(3,qm2ds%Mx)          ! Dipole moment
	_REAL_ :: nuc_dipole(3), mu_gr(3),summc
        _REAL_ :: ddot ! This is a function
	_REAL_, dimension(:),allocatable :: ex_mchg;
	_REAL_, dimension(:,:), allocatable :: GSDM,ESDM, T, TZ, dip, tmp
	_REAL_, parameter :: convert_to_debye = charge_on_elec*(1.0d-10)*(light_speed/(1.0d-21)); !This conversion is from MOPAC manual
	_REAL_, parameter :: convert_debye_to_AU = 2.541746
	integer :: Lt,i,j,k,l,nquant,state
	logical flag;

	allocate(dip(3,qm2ds%Nb*(qm2ds%Nb+1)/2));
        allocate(GSDM(qm2ds%Nb,qm2ds%Nb),tmp(qm2ds%Nb,qm2ds%Nb));

        if(qmmm_nml%printcharges) then
        allocate(ex_mchg(qmmm_struct%nquant_nlink));
        end if

!nuclear part
	call get_nuc_dip(qmmm_struct%qm_coords, dip) !return nuclear dipole operator in angstroms
	do k=1,3 ! loop over x,y,z
		call unpacking(qm2ds%Nb,dip(k,:),tmp,'s')
		summc=0.d0
		do i=1,qm2ds%Nb
			summc=summc+tmp(i,i)
		end do
		nuc_dipole(k)=summc
	end do
!end nuclear part

!ground state
	call get_dipole_matrix(qmmm_struct%qm_coords, dip);
	call unpacking(qm2ds%Nb,qm2_struct%den_matrix,GSDM,'s')
        do k=1,3  ! loop over x,y,z
		call unpacking(qm2ds%Nb,dip(k,:),qm2ds%eta_scratch,'s')
                mu_gr(k) = ddot(qm2ds%Nb**2,GSDM,1,qm2ds%eta_scratch,1) + nuc_dipole(k)
        enddo
	mu_gr=mu_gr*convert_to_debye/convert_debye_to_AU
        write(6,*)
        write(6,*) 'Ground State Molecular Dipole Moment (A.U.)'
        write(6,"(25x,'dx',14x,'dy',14x,'dz',10x,'ftotal')")
        write(6,"(19x,4g15.7)") mu_gr(1),mu_gr(2),mu_gr(3),sqrt(mu_gr(1)**2 + mu_gr(2)**2 + mu_gr(3)**2)

!excited states
if (qm2ds%Mx>0) then
	allocate(TZ(qm2ds%Nb,qm2ds%Nb));
        allocate(T(qm2ds%Nb,qm2ds%Nb));
        allocate(ESDM(qm2ds%Nb,qm2ds%Nb));

	do state=1,qm2ds%Mx
	        call calc_rhotz(state,qm2ds%rhoTZ,.true.); 
                call calc_rhotz(state,qm2ds%rhoT,.false.);
		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoTZ,TZ,qm2ds%tz_scratch)		
                call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%rhoT,T,tmp)

		ESDM=GSDM+TZ !Construct excited state density matrix

		!The mulliken charges are here but this is a bad place for them
		if(qmmm_nml%printcharges) then		
	        	do i=1,qmmm_struct%nquant_nlink
				call qm2_calc_mulliken(i,ex_mchg(i),ESDM)
	        	end do
			write (6,*)
			write (6,'("QMMM: Mulliken Charges")')
			call qm2_print_charges(state,1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                                        ex_mchg,qmmm_struct%iqm_atomic_numbers)
		end if

                !tr(Mu Rho) calculate dipole
		do k=1,3  ! loop over x,y,z
	 	   	call unpacking(qm2ds%Nb,dip(k,:),qm2ds%eta_scratch,'s')
		   	mu(k,state)=ddot(qm2ds%Nb**2,ESDM,1,qm2ds%eta_scratch,1)
                   	mu(k,state)=mu(k,state)+nuc_dipole(k);
                   	mu_unrelaxed(k,state)=ddot(qm2ds%Nb**2,T,1,qm2ds%eta_scratch,1) 
                   	mu_relaxed(k,state)=ddot(qm2ds%Nb**2,TZ,1,qm2ds%eta_scratch,1)
                enddo

	end do

	deallocate(T,TZ,ESDM,GSDM,dip,tmp);

	if(qmmm_nml%printcharges) then
		deallocate(ex_mchg);
	end if

	mu(1:3,:)=mu(1:3,:) * convert_to_debye;
	if (qmmm_nml%printdipole > 0) then
		write(6,*)
        	write(6,*) 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)'
	        write(6,"(8x,'Omega',12x,'dx',14x,'dy',14x,'dz',10x,'ftotal')")
		do i=1,qm2ds%Mx
	          write(6,"(i4,5g15.7)") i,qm2ds%e0(i), &
        	        mu(1,i), &
                	mu(2,i), &
	                mu(3,i), &
			sqrt(mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
	        end do
		mu(1:3,:)=mu(1:3,:)/convert_debye_to_AU
	        write(6,*)
        	write(6,*) 'Frequencies (eV) and Total Molecular Dipole Moments (AU)'
	        write(6,"(8x,'Omega',12x,'dx',14x,'dy',14x,'dz',10x,'ftotal')")
        	do i=1,qm2ds%Mx
	          write(6,"(i4,5g15.7)") i,qm2ds%e0(i), &
        	        mu(1,i), &
                	mu(2,i), &
	                mu(3,i), &
        	        sqrt(mu(1,i)**2 + mu(2,i)**2 + mu(3,i)**2)
	        end do
                
                !PRINTING RELAXED AND UNRELAXED DIPOLES
                if (qmmm_nml%printdipole > 1) then
                mu_relaxed(1:3,:)=mu_relaxed(1:3,:) * convert_to_debye/convert_debye_to_AU;
                mu_unrelaxed(1:3,:)=mu_unrelaxed(1:3,:) * convert_to_debye/convert_debye_to_AU;

                write(6,*)
                write(6,*) 'Frequencies (eV) Unrelaxed Difference Dipole Moments (AU)'
                write(6,"(8x,'Omega',12x,'dx',14x,'dy',14x,'dz',10x,'ftotal')")
                do i=1,qm2ds%Mx
                  write(6,"(i4,5g15.7)") i,qm2ds%e0(i), &
                        mu_unrelaxed(1,i), &
                        mu_unrelaxed(2,i), &
                        mu_unrelaxed(3,i), &
                        sqrt(mu_unrelaxed(1,i)**2 + mu_unrelaxed(2,i)**2 + mu_unrelaxed(3,i)**2)
                 end do
                write(6,*)
                write(6,*) 'Frequencies (eV) Relaxed Difference Dipole Moments (AU)'
                write(6,"(8x,'Omega',12x,'dx',14x,'dy',14x,'dz',10x,'ftotal')")
                do i=1,qm2ds%Mx
                  write(6,"(i4,5g15.7)") i,qm2ds%e0(i), &
                        mu_relaxed(1,i), &
                        mu_relaxed(2,i), &
                        mu_relaxed(3,i), &
                        sqrt(mu_relaxed(1,i)**2 + mu_relaxed(2,i)**2 + mu_relaxed(3,i)**2)
                 end do
                end if
	end if
end if
end 

!subroutine md_trans_dipole(mu, alpha)
!	use qm2_davidson_module
!	use qmmm_module, only: qmmm_struct
!
!	implicit none
!
!	_REAL_, intent(out) :: mu(3, qm2ds%Mx)		! Dipole moment
!	_REAL_, intent(out) :: alpha(3, qm2ds%Mx)		! Polarizability
!
!	_REAL_ :: drho(qm2ds%Nb*(qm2ds%Nb+1)/2)
!	_REAL_ :: dip(3, qm2ds%Nb*(qm2ds%Nb+1)/2) ! Dipole matrix elements in 3D
!	_REAL_ :: v0_ao(qm2ds%Nb,qm2ds%Nb)	! Lanczos eigenvector in AO
!	_REAL_ :: fn	! Normalization constant
!
!	integer :: i, j, k, m, n, Lt
!
!	fn = 1/sqrt(2.0)
!	Lt=qm2ds%Nb*(qm2ds%Nb+1)/2
!	
!	mu = 0.d0
!	alpha = 0.d0
!
!!	do i=1,Lt
!!		rrwork(i)=eta(i)
!!	end do
!
!!	call response11 (Nb,Nb_M,Mx,M2,e0,dip,rrwork(1), &
!!		eta,xi,rrwork(Lt_M+1),rrwork(2*Lt_M+1),mu,polar)
!
!!	do i=1,Lt
!!		xi(i)=rrwork(Lt_M+i)
!!	end do
!
!	call wrb_arr(dip(1), Lt, 'dip0x.b', 'r', 's')
!	call wrb_arr(dip(Lt+1), Lt, 'dip0y.b', 'r', 's')
!	call wrb_arr(dip(2*Lt+1), Lt, 'dip0z.b', 'r', 's')
!
!	do i = 1,qm2ds%Mx
!		call mo2site(qm2ds%v0(:,i), v0_ao, qm2ds%eta_scratch)
!		do j = 1,3
!			drho = 0.d0
!			do m = 1,qm2ds%Nb
!				do n = 1,m-1
!					mu(j,i) = mu(j,i)+2*(v0_ao(n,m)+v0_ao(m,n))*dip(j,m*(m-1)/2+n)*fn
!				end do
!				mu(j,i) = mu(j,i)+2*v0_ao(m,m)*dip(j,m*(m-1)/2+m)*fn
!			end do
!			do m = 1,qm2ds%Nb
!				do n = 1,m
!					k=m*(m-1)/2+n
!					drho(k) = drho(k)+mu(j,i)*(v0_ao(n,m)+v0_ao(m,n))*fn/qm2ds%e0(i)
!				end do
!			end do
!	
!		end do
!	end do
!
!		
!
!	return
!
!end subroutine md_trans_dipole


subroutine get_nuc_dip(coord, dpmat)

	 use qmmm_module, only : qm2_params, qm2_struct, qmmm_struct, qmmm_nml
	 use constants, only : CODATA08_A_TO_BOHRS, bohr_radius
! use findmask
      implicit none

! Calculates the dipole matrices

      _REAL_ :: coord(*)
        ! CML TEST 7/17/12
      _REAL_ :: dpmat(3,qm2_struct%norbs*(qm2_struct%norbs+1)/2)

      integer :: loop_count,orb_beg,orb_end,orb_size,i,nquant
      ! AWG: Currently not supported with d orbitals for NDDO methods
      if ( .not. qmmm_nml%qmtheory%DFTB ) then
         do i = 1,qmmm_struct%nquant_nlink
            if ( qm2_params%natomic_orbs(i) > 5 ) return
         end do
      end if
	dpmat=0.0
      do nquant=1,qmmm_struct%nquant_nlink
         ! Asking for the Method, if DFTB is used it does not compute the changes in charge due to orbitals.
         if (.not. qmmm_nml%qmtheory%DFTB) then
            orb_beg=qm2_params%orb_loc(1,nquant)
            orb_end=qm2_params%orb_loc(2,nquant)

            ! Determination of the values difference between the orbitals
            orb_size=orb_end-orb_beg

                        ! Fill diagonal elements of dipole matrices with coordinate * Nuclear  charge (Z)
                        do loop_count = orb_beg, orb_beg
                                do i = 1,3
                                   dpmat(i,qm2_params%pascal_tri2(loop_count)) = &
					qm2_params%core_chg(nquant)*coord((nquant-1)*3+i)
                                end do
                        end do

         end if
      end do
end




subroutine write_down(m,m2,sizem)
	integer sizem,i,j;
	_REAL_ m(sizem,sizem),m2(sizem,sizem)
	_REAL_ summ
	summ=0.
	open(unit=243,file="dip");
	do i=1,sizem
	  do j=1,sizem
		if (i==j) summ=summ+m(i,j)*m2(i,j)
		if (i==j) write(243,*)i,j,m(i,j)*m2(i,j),summ
	  end do
	end do	

	close(243)
end

subroutine write_down_correct(m,sizem)
	use constants, only : BOHRS_TO_A, SQRT2
        integer	sizem,i,j;
        _REAL_ m(sizem,sizem)
        open(unit=244,file="dip");
        do i=1,sizem
          do j=1,sizem
                 write(244,*)i,j, 0.5*(m(i,j)-m(j,i))


         end do
       	end do
        close(244)
end

subroutine gaddt(sizem,m1,m2)
	integer sizem
	_REAL_ m1(sizem,sizem),m2(sizem,sizem)
	m1(1:sizem,1:sizem)=m1(1:sizem,1:sizem)+m2(1:sizem,1:sizem)
end 
