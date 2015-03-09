#include "dprec.fh"

subroutine calc_rhotz(state, rhoTZ,calc_Z)
        use qm2_davidson_module
        use qmmm_module
        use cosmo_C,only:ceps,potential_type,solvent_model,v_solvent_difdens,v_solvent_xi
        implicit none
	logical calc_Z; !This flag doesn't seem to work
        integer, intent(inout) :: state
        _REAL_, intent(out) :: rhoTZ(qm2ds%Nb**2)
        _REAL_ tmp(qm2ds%Nb,qm2ds%Nb)
        _REAL_ f,t,d,d1,ff,ff0,ff1,ff11
        integer i,ii,j,k,im,one,istate,mdflag,ip,ih,solvent_model_1

	if((solvent_model.gt.0).and.(solvent_model.ne.4)) then
		solvent_model_1=99 !Use linear response solvent model on lhs of Z-vector equation
                !should actually be able to just use the solvent effects on the
                !molecular orbital energies here. See paper 2 for discussion.
	else
		solvent_model_1=0
	endif
        tmp=0.d0
        one=1
    !KGB mdflag -> 
        qm2ds%mdflag=0
        ff0=0.0
        ff1=1.0
        ff11=-1.0
! Verify Davidson has been run, if not, kill the program...perhaps calling
! sqm_energy() might be wise...but it requires so many input arguments that aren't here...CML 7/17/12
        if (qmmm_struct%qm_mm_first_call) then
                write(6,*)  'sqm_energy() must be run once before executing this procedure!'
                if (qm2ds%Mx == 0) write(6,*)  'excN must be > 0 to run this procedure!'
                call mexit(6,1)
        end if

! Calculate unrelaxed density of a given excited state
! T=[[xi^+ rho] xi]= (I-2rho)(xi^+ xi+xi xi^+)
        call getmodef(2*qm2ds%Np*qm2ds%Nh,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
                state,qm2ds%v0,qm2ds%tz_scratch)
        call dgemm('N','T',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,ff1,qm2ds%tz_scratch, &
                qm2ds%Nb,qm2ds%tz_scratch,qm2ds%Nb,ff0,qm2ds%eta_tz,qm2ds%Nb)
        call dgemm('T','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,ff1,qm2ds%tz_scratch, &
                qm2ds%Nb,qm2ds%tz_scratch,qm2ds%Nb,ff1,qm2ds%eta_tz,qm2ds%Nb)
        call Iminus2rho(qm2ds%Nb,qm2ds%Np,qm2ds%eta_tz,rhoTZ) 

if (calc_Z) then 
! Start putting together left hand side of LZ-equation
! Form [V(rhoT), rho]=(I-2rho)V(rhoT) 
        call mo2sitef (qm2ds%Nb,qm2ds%vhf,rhoTZ, &
                qm2ds%eta_tz,qm2ds%tz_scratch(qm2ds%Nb**2+1))
        qm2ds%tz_scratch(:)=0.d0;
        call Vxi (qm2ds%eta_tz,qm2ds%tz_scratch(1))
!**************SOLVENT BLOCK !!JAB
        if(solvent_model.gt.0) then
        	tmp=0.d0
        if (potential_type.eq.3) then !COSMO Potential
        	call VxiM(qm2ds%eta_tz,tmp); 
        else if (potential_type.eq.2) then !Onsager Potential
        	call rcnfld(tmp,qm2ds%eta_tz,qm2ds%nb); 
        end if
		tmp=2.d0*tmp !linear response
        	call VxiM_end(qm2ds%tz_scratch(1),tmp); !Add selected potential to vacuum correlation
        end if
!!************END SOLVENT BLOCK

!!************GAS PHASE BLOCK
        call site2mof (qm2ds%Nb,qm2ds%vhf,qm2ds%tz_scratch(1),qm2ds%xi_tz, &
                     qm2ds%tz_scratch(qm2ds%Nb**2+1))
        call Iminus2rho(qm2ds%Nb,qm2ds%Np,qm2ds%xi_tz,qm2ds%rhoLZ) 
       
        call project(qm2ds%Nb,qm2ds%Np,qm2ds%Nh,qm2ds%rhoLZ)

! Form [[xi^+, rho], V(xi)], rho] +cc= (I-2rho)[(I-2rho)xi^+, V(xi)]+cc
        call getmodef(2*qm2ds%Np*qm2ds%Nh,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
                        state,qm2ds%v0,qm2ds%xi_tz) !This may be redundant with the call above FIXME
        call mo2sitef (qm2ds%Nb,qm2ds%vhf,qm2ds%xi_tz,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
                qm2ds%eta_tz(:)=0.d0;
        call Vxi(qm2ds%tz_scratch(1),qm2ds%eta_tz)

!**************END GAS PHASE BLOCK

!**************SOLVENT BLOCK to add V_s(xi)
if((solvent_model.eq.1)) then !Linear Response solvent
        tmp=0.d0;
        if (potential_type.eq.3) then !COSMO Potential
        	call VxiM(qm2ds%tz_scratch(1),tmp); 
        elseif (potential_type.eq.2) then !Onsager Potential
        	call rcnfld(tmp,qm2ds%tz_scratch(1),qm2ds%nb);
        end if
        call VxiM_end(qm2ds%eta_tz,tmp); !Add selected potential to vacuum correlation
endif
!!************END SOLVENT BLOCK

!!************GAS PHASE BLOCK
        call site2mof (qm2ds%Nb,qm2ds%vhf,qm2ds%eta_tz,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
        call transp1(qm2ds%Nb,qm2ds%xi_tz)
        call Iminus2rho(qm2ds%Nb,qm2ds%Np,qm2ds%xi_tz,qm2ds%eta_tz)
        call dgemm('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,ff1,qm2ds%eta_tz, &
                        qm2ds%Nb,qm2ds%tz_scratch(1),qm2ds%Nb,ff0,qm2ds%xi_tz,qm2ds%Nb)
        call dgemm('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,ff11,qm2ds%tz_scratch(1), &
                        qm2ds%Nb,qm2ds%eta_tz,qm2ds%Nb,ff1,qm2ds%xi_tz,qm2ds%Nb)   
        call symmetr(qm2ds%Nb,qm2ds%xi_tz) 
        !call Iminus2rho(qm2ds%Nb,qm2ds%Np,qm2ds%xi_tz,qm2ds%eta_tz) !Now these are done below after a solvent block
        !call project(qm2ds%Nb,qm2ds%Np,qm2ds%Nh,qm2ds%eta_tz)
!************END GAS PHASE BLOCK

!!***********BEGIN SOLVENT BLOCK
if((solvent_model.eq.2).or.(solvent_model.eq.4)) then !VE and SS solvent
! Add [[[xi^+_k, V_S(T)], xi_k],rho] + cc by calculating commutators FIXME it
! currently does T_k and xi_n so it will only work for the state of interest
        call getmodef(2*qm2ds%Np*qm2ds%Nh,qm2ds%Mx,qm2ds%Np,qm2ds%Nh, &
                        state,qm2ds%v0,qm2ds%eta_tz)
        call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%eta_tz,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
        call commutator(qm2ds%tz_scratch(1),v_solvent_difdens,qm2ds%Nb,tmp,.true.)!inner commutator
        call commutator(tmp,qm2ds%tz_scratch(1),qm2ds%Nb,qm2ds%eta_tz,.false.) !second commutator with transpose
        call site2mof(qm2ds%Nb,qm2ds%vhf,qm2ds%eta_tz,qm2ds%tz_scratch(1), &
                        qm2ds%tz_scratch(qm2ds%Nb**2+1))
        qm2ds%xi_tz=qm2ds%xi_tz-1.0*qm2ds%tz_scratch(1:qm2ds%Nb**2)
endif
!!***********END SOLVENT BLOCK 

!************BEGIN VACUUM BLOCK
        call Iminus2rho(qm2ds%Nb,qm2ds%Np,qm2ds%xi_tz,qm2ds%eta_tz)
        call project(qm2ds%Nb,qm2ds%Np,qm2ds%Nh,qm2ds%eta_tz)

! Finally put together left hand side of LZ-equation
! 2 is coming from the complex comjugate for eta
        do i=1,qm2ds%Nb**2
                qm2ds%rhoLZ(i)= qm2ds%rhoLZ(i) - 2*qm2ds%eta_tz(i)
        enddo

! Now solve for Z equation LZ=rhoLZ
! Start with 0-order
        call continter(qm2ds%Nb,qm2ds%Np,qm2ds%Nh,qm2ds%Np*qm2ds%Nh, &
                qm2ds%tz_scratch(3*2*qm2ds%Np*qm2ds%Nh+1),qm2ds%rhoLZ)
        call dcopy(2*qm2ds%Np*qm2ds%Nh,qm2ds%tz_scratch(3*2*qm2ds%Np*qm2ds%Nh+1), &
                one,qm2ds%tz_scratch(1),one)

                f=1.0
! Start loop here
		!'RhoTZ min iter set to at least 100 for testing... FIXME'    
                do im=1,max(qm2ds%icount_M,100)    ! JAB switched to at least 100
                        i = 0
                        do ip = 1,qm2ds%Np
                                do ih = qm2ds%Np+1,qm2ds%Nb
                                        i = i + 1
                                        ff = 1/(qm2ds%ehf(ih) - qm2ds%ehf(ip))
                                        qm2ds%eta_tz(i) = ff*qm2ds%tz_scratch(3*2*qm2ds%Np*qm2ds%Nh+i)
                                        qm2ds%eta_tz(i+qm2ds%Np*qm2ds%Nh) = &
                                                -ff*qm2ds%tz_scratch(7*qm2ds%Np*qm2ds%Nh+i) ! 3*M2+M4+i
                                end do
                        end do
                        if (im.eq.1) then
                                call dcopy(2*qm2ds%Np*qm2ds%Nh,qm2ds%eta_tz,one, &
                                        qm2ds%tz_scratch(2*qm2ds%Np*qm2ds%Nh+1),one)
                        else
                                call summing(2*qm2ds%Np*qm2ds%Nh,qm2ds%eta_tz, &
                                        qm2ds%tz_scratch(2*qm2ds%Np*qm2ds%Nh+1))
                        end if

                        call Lxi_testing(qm2ds%tz_scratch(2*qm2ds%Np*qm2ds%Nh+1), &
                                qm2ds%tz_scratch(4*qm2ds%Np*qm2ds%Nh+1),solvent_model_1)
! Check for convergency
                        f=0.0
                        do j=1,2*qm2ds%Np*qm2ds%Nh
                                qm2ds%tz_scratch(6*qm2ds%Np*qm2ds%Nh+j)= &
                                        qm2ds%tz_scratch(j)-qm2ds%tz_scratch(4*qm2ds%Np*qm2ds%Nh+j)
                                f=f+abs(qm2ds%tz_scratch(6*qm2ds%Np*qm2ds%Nh+j)**2)
                        end do
! If converged, exit    
                        if(f.lt.qm2ds%ftol0**2) then ! converged
                                call expinter(qm2ds%Nb,qm2ds%Np,qm2ds%Nh,qm2ds%Np*qm2ds%Nh, &
                                        qm2ds%tz_scratch(2*qm2ds%Np*qm2ds%Nh+1),qm2ds%eta_tz)
                                do j=1,qm2ds%Nb**2
                                        rhoTZ(j)=rhoTZ(j)+qm2ds%eta_tz(j)
                                end do
                                exit
                        end if 
! Quit if convergence is not achieved
                        if(im.eq.qm2ds%icount_M*2) then !JAB *2 more for solvent
                                write(6,*)  'Eq. for Z: LZ=rhoLZ did not converge'
                                write(6,*)  'Achieved convergence= ', f
                                call mexit(6,1)
                        end if 
                end do
        end if
        return
end subroutine calc_rhotz
