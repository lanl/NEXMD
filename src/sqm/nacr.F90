#include "dprec.fh"

subroutine nacR_analytic(xyz_in, ihop, icheck)
	use qmmm_module, only : qmmm_struct, qm2_struct
	use qm2_davidson_module
	use constants, only : KCAL_TO_EV

	implicit none

	! NEW OR MODIFIED VARIABLES
	integer :: M4_M, M2_M, Np, Nh, Nb, Mx_M, Mx
	integer :: ideriv ! This will come from a module later
	_REAL_, intent(in) :: xyz_in(3*qmmm_struct%nquant_nlink)
    integer, intent(in):: ihop,icheck

	! OLD VARIABLES (MAY BE DEPRECATED)
	_REAL_ fPi,fbar,f,t,d,d1,ff,ff0,ff1,ff11,ddot
    integer i,ii,j,k,im,one,istate,mdflag,ip,ih,Nm,N3

	parameter (fPi = 3.1415926535898d0)
	parameter (fbar = 0.05d0)  ! maximum dE in numerical derivative, eV.
    integer Na,npot
!	_REAL_ xx(qmmm_struct%nquant_nlink),yy(qmmm_struct%nquant_nlink),zz(qmmm_struct%nquant_nlink)
!	_REAL_ xx1(qmmm_struct%nquant_nlink),yy1(qmmm_struct%nquant_nlink),zz1(qmmm_struct%nquant_nlink)
    _REAL_ xyz(3,qmmm_struct%nquant_nlink),dxyz1(3,qmmm_struct%nquant_nlink)

	Na = qmmm_struct%nquant_nlink
    N3 = 3*qmmm_struct%nquant_nlink
    one=1
    ff0=0.0
    ff1=1.0
    ff11=-1.0
	
	Na = qmmm_struct%nquant_nlink
	Np = qm2ds%Np
	Nh = qm2ds%Nh
	Nb = qm2ds%Nb
	M4_M = qm2ds%Np*qm2ds%Nh
	M2_M = M4_M*2
	Mx = qm2ds%Mx 
	Mx_M = Mx
	npot = Mx

! This will need to be put outside of this routine because cmdqt is important for NAD calcs
!	do j=1,npot
!		f= ddot(nbasis,cmdqt(1,j),1,v0(1,kx(j)),1)
!		if (f.lt.0.0d0) then
!			do i=1,nbasis
!				cmdqt(i,j)=-v0(i,kx(j))
!			end do
!		else
!			do i=1,nbasis
!				cmdqt(i,j)=v0(i,kx(j))
!			enddo
!		end if
!	end do
! But, it's here now to test this function
! Note, Davidson must be run prior to this for these assignments
	if (qmmm_struct%qm_mm_first_call) then
		write(6,*)  'sqm_energy() must be run once before executing this procedure!'
		if (qm2ds%Mx == 0) write(6,*)  'excN must be > 0 to run this procedure!'
		call mexit(6,1)
	end if

!	do j=1,qm2ds%Mx
!		f= ddot(M4_M,qm2ds%cmdqt(1,j),1,qm2ds%v0(1,qm2ds%kx(j)),1)
!		if (f.lt.0.0d0) then
!			do i=1,M4_M
!				qm2ds%cmdqt(i,j) = -qm2ds%v0(i,qm2ds%kx(j))
!			end do
!		else
!			do i=1,M4_M
!				qm2ds%cmdqt(i,j) = qm2ds%v0(i,qm2ds%kx(j))
!			enddo
!		end if
!	end do

    !KGB commented 2 lines out
	!ihop = 2  ! initial state
	!icheck = 1	! final state
	
! Form transition density martix between state ihop and icheck
!  xi_ih_ic = (I-2rho)(xi_ih xi_ic + xi_ic xi_ih)

	call getmodef(M2_M,Mx_M,Np,Nh,ihop,qm2ds%cmdqt,qm2ds%nacr_scratch)
	call getmodef(M2_M,Mx_M,Np,Nh,icheck,qm2ds%cmdqt,qm2ds%eta_scratch)
	
	call dgemm('N','T',Nb,Nb,Nb,ff1,qm2ds%nacr_scratch,Nb,qm2ds%eta_scratch,Nb,ff0,qm2ds%eta,Nb)
	call dgemm('T','N',Nb,Nb,Nb,ff1,qm2ds%eta_scratch,Nb,qm2ds%nacr_scratch,Nb,ff1,qm2ds%eta,Nb)
	call Iminus2rho(Nb,Np,qm2ds%eta,qm2ds%xi)
	call mo2sitef (Nb,qm2ds%vhf,qm2ds%xi,qm2ds%eta,qm2ds%xi_scratch)
	
! Above eta contains transition density martix between state ihop and icheck in AO

	call getmodef(M2_M,Mx_M,Np,Nh,ihop,qm2ds%cmdqt,qm2ds%nacr_scratch)
	call mo2sitef (Nb,qm2ds%vhf,qm2ds%nacr_scratch,qm2ds%xi,qm2ds%eta_scratch)
	call getmodef(M2_M,Mx_M,Np,Nh,icheck,qm2ds%cmdqt,qm2ds%nacr_scratch)
	call mo2sitef (Nb,qm2ds%vhf,qm2ds%nacr_scratch,qm2ds%xi_scratch_2,qm2ds%eta_scratch)

! Above xi and xi_scratch_2 contain transition density martices for states ihop and icheck in AO
	ideriv = 3
	if (ideriv.eq.3) then  ! Fast GS MOPAC derivatives
	    do i = 1,qmmm_struct%nquant_nlink
	        do j = 1,3
       		    xyz(j,i) = xyz_in((i-1)*3+j)
				dxyz1(j,i) = 0.d0
	        end do
	    end do
! Term Tr(F^x rho_ij) (only symmetric part contributes)
		call packing(Nb,qm2ds%eta,qm2ds%nacr_scratch,'s')
		call DCART1(dxyz1,qm2_struct%den_matrix,qm2ds%nacr_scratch,xyz)
! FIXME Why den_matrix is twice of that in old version?

! Convert from kcal/A to eV/A
		do j = 3,N3,3
			qm2ds%dij(j-2)= dxyz1(1,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(icheck))-qm2ds%e0(qm2ds%kx(ihop)))
			qm2ds%dij(j-1)= dxyz1(2,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(icheck))-qm2ds%e0(qm2ds%kx(ihop)))
			qm2ds%dij(j)  = dxyz1(3,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(icheck))-qm2ds%e0(qm2ds%kx(ihop)))
		end do
	end if

	!write(6,*)  'dij'
	!do i = 1,Na*3
    !write(6,*)  qm2ds%dij(i)
	!end do

	flush(6)
	!STOP

	return

end subroutine nacR_analytic
