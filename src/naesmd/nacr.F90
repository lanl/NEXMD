#include "dprec.fh"

subroutine nacR_analytic(qm2_params,qmmm_nml,qmmm_mpi,qm2_rij_eqns,qm2_struct,qm2ds,qmmm_struct, xyz_in, ihop, icheck)
    use qmmm_module, only : qm2_structure, qmmm_mpi_structure, qm2_rij_eqns_structure
    use qm2_davidson_module
    use constants, only : KCAL_TO_EV
    use qmmm_struct_module, only : qmmm_struct_type
    use qmmm_nml_module   , only : qmmm_nml_type
    use qm2_params_module,  only : qm2_params_type

    implicit none
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi


    ! NEW OR MODIFIED VARIABLES
    type(qm2_davidson_structure_type), intent(inout) :: qm2ds
    integer :: M4_M, M2_M, Np, Nh, Nb, Mx_M, Mx
    integer :: ideriv ! This will come from a module later
    _REAL_, intent(in) :: xyz_in(3*qmmm_struct%nquant_nlink)
    integer, intent(in):: ihop,icheck

    ! OLD VARIABLES (MAY BE DEPRECATED)
    _REAL_ fPi,fbar,f,t,d,d1,ff,ff0,ff1,ff11,ddot
    integer i,ii,j,k,im,one,istate,mdflag,ip,ih,Nm,N3

    parameter (fPi = 3.1415926535898d0)
    parameter (fbar = 0.05d0)  ! maximum dE in numerical derivative, eV.
    integer Na
    _REAL_ xyz(3,qmmm_struct%nquant_nlink),dxyz1(3,qmmm_struct%nquant_nlink)
    if (qmmm_nml%verbosity.gt.4) write(6,*)'nacR_analytic called'

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


    if (qmmm_struct%qm_mm_first_call) then
        write(6,*)  'sqm_energy() must be run once before executing this procedure!'
        if (qm2ds%Mx == 0) write(6,*)  'excN must be > 0 to run this procedure!'
        call mexit(6,1)
    end if

    qm2ds%nacr_scratch=0
    call getmodef(M2_M,Mx_M,Np,Nh,ihop,qm2ds%cmdqt,qm2ds%nacr_scratch)
    call getmodef(M2_M,Mx_M,Np,Nh,icheck,qm2ds%cmdqt,qm2ds%eta_scratch)
    call dgemm('N','T',Nb,Nb,Nb,ff1,qm2ds%nacr_scratch,Nb,qm2ds%eta_scratch,Nb,ff0,qm2ds%eta,Nb)
    call dgemm('T','N',Nb,Nb,Nb,ff1,qm2ds%eta_scratch,Nb,qm2ds%nacr_scratch,Nb,ff1,qm2ds%eta,Nb)
    call Iminus2rho(Nb,Np,qm2ds%eta,qm2ds%xi)
    call mo2sitef (Nb,qm2ds%vhf,qm2ds%xi,qm2ds%eta,qm2ds%xi_scratch)
    ! Above eta contains transition density matrix between state ihop and icheck in AO

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
        call DCART1(qmmm_nml, qm2_params, qmmm_mpi, qm2_rij_eqns, qm2_struct,qm2ds,qmmm_struct, &
		 dxyz1,qm2_struct%den_matrix,qm2ds%nacr_scratch,xyz)

        ! Convert from kcal/A to eV/A
        do j = 3,N3,3
            qm2ds%dij(j-2)= dxyz1(1,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(ihop))-qm2ds%e0(qm2ds%kx(icheck)))
            qm2ds%dij(j-1)= dxyz1(2,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(ihop))-qm2ds%e0(qm2ds%kx(icheck)))
            qm2ds%dij(j)  = dxyz1(3,j/3)*KCAL_TO_EV/(qm2ds%e0(qm2ds%kx(ihop))-qm2ds%e0(qm2ds%kx(icheck)))
        end do
    end if

	    write(6,*) "no"

    return

end subroutine nacR_analytic
