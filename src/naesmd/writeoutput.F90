#include "dprec.fh"
#include "assert.fh"
module writeoutput_module
    use naesmd_constants
    !***********************************************
    ! write output information
    ! at the initial step
    !***********************************************
    implicit none

contains

    subroutine writeoutputini(sim,ibo,yg,yg_new,lprint)
        use naesmd_constants
        use naesmd_module
        use communism
        use constants
        implicit none

        type(simulation_t), target :: sim
        type(simulation_t), pointer :: simpoint
        integer  i,j,k,l,n,ibo,kki,lprint
        _REAL_ ntot
        _REAL_ xcm,ycm,zcm
        _REAL_ NACR(sim%naesmd%natom*3)
        character*1000 card, tmpcard
        
        character*10 file_status
        character*10 file_access

        _REAL_ yg(:)
        _REAL_ yg_new(:)
        _REAL_ energy !!JAB VE model
        integer unitcard
        character*1000 line

        simpoint => sim

        if(sim%cosmo%solvent_model.eq.2) then
            call calc_excsolven(sim%cosmo,sim%dav,sim%qmmm,energy) !JAB Test
            sim%naesmd%vmdqt(sim%naesmd%ihop)=sim%naesmd%vmdqt(sim%naesmd%ihop)-0.5*energy/feVmdqt !JAB Test
        endif

        sim%naesmd%kin=0.0d0
        do i =1, sim%naesmd%natom
            sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(i)*(sim%naesmd%vx(i)**2+sim%naesmd%vy(i)**2+sim%naesmd%vz(i)**2)/2
        end do

        ntot=0
        do j=1,sim%excN
            ntot=ntot+yg(j)**2
        end do

        sim%naesmd%kinini=sim%naesmd%kin
        if(sim%naesmd%state.eq.'fund') then
            sim%naesmd%vini=sim%naesmd%vgs
            sim%naesmd%etotini=sim%naesmd%kin+sim%naesmd%vgs
            write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%vgs*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt-sim%naesmd%etotini*feVmdqt
        end if

        if(sim%naesmd%state.eq.'exct') then
            if(sim%naesmd%dynam_type.eq.'tully') then
              sim%naesmd%vini=sim%naesmd%vmdqt(sim%naesmd%ihop)
              sim%naesmd%etotini=sim%naesmd%kin+sim%naesmd%vmdqt(sim%naesmd%ihop)
              write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt, &
                sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt,&
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt- &
                sim%naesmd%etotini*feVmdqt
            elseif(sim%naesmd%dynam_type.eq.'aimc'.or.sim%naesmd%dynam_type.eq.'mf') then
              sim%naesmd%vmfini=0.0d0
              do j=1,sim%excN
                sim%naesmd%vmfini=sim%naesmd%vmfini+yg(j)**2*sim%naesmd%vmdqt(j)
              enddo
              sim%naesmd%etotini=sim%naesmd%kin+sim%naesmd%vmfini
              write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vmfini*feVmdqt, &
                sim%naesmd%vmfini*feVmdqt-sim%naesmd%vmfini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmfini*feVmdqt,&
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmfini*feVmdqt- &
                sim%naesmd%etotini*feVmdqt
            endif
        end if

        if(sim%naesmd%dynam_type.eq.'aimc') then
            write(sim%outfile_29,889) sim%naesmd%tfemto, 0.0d0 
        endif
        
        if(lprint.ge.1) then
            if(sim%excN.gt.0) then
                write(sim%outfile_4,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                    (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                call flush(sim%outfile_4)
            endif
        endif

        if(sim%naesmd%state.eq.'exct') then
            if(lprint.ge.1) then
                if(ibo.ne.1) then
                    write(sim%outfile_8,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot
                    do j=1,sim%excN
                      do k=j+1,sim%excN
                        call nacR_analytic_wrap(simpoint,j,k,NACR)
                        NACR=sim%naesmd%sgn(j,k)*NACR
                        sim%naesmd%cadiabnew(j,k)=0.0d0
                        do l=1,sim%naesmd%natom
                          sim%naesmd%cadiabnew(j,k)=sim%naesmd%cadiabnew(j,k)+sim%naesmd%vx(l)*NACR(3*l-2)*convl+sim%naesmd%vy(l)*NACR(3*l-1)*convl+sim%naesmd%vz(l)*NACR(3*l)*convl
                        enddo
                        sim%naesmd%cadiabnew(k,j)=-1.0d0*sim%naesmd%cadiabnew(j,k)
                      enddo
                    enddo
                    write(sim%outfile_7,888) sim%naesmd%tfemto,((sim%naesmd%cadiabnew(j,k),k=1,sim%excN),j=1,sim%excN)
                end if
            elseif(sim%naesmd%dynam_type.eq.'aimc') then
               write(sim%outfile_8,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot 
            end if

            if(lprint.ge.3.and.ibo.ne.1) then
                !write(sim%outfile_16,889) sim%naesmd%tfemto,(dsin(yg(j+sim%excN)),j=1,sim%excN)
                write(sim%outfile_16,889) sim%naesmd%tfemto,(dmod(yg(j+sim%excN),TWOPI),j=1,sim%excN)
                call flush(sim%outfile_16)
            end if

            call flush(sim%outfile_5)
        endif

        !State dependent forces needed for Multiconfigurational post-processing
        if(sim%naesmd%dynam_type.eq.'aimc'.or.sim%naesmd%dynam_type.eq.'mf') then
            write (sim%outfile_25,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do n=1,sim%excN
            do k=1,sim%naesmd%natom
                write(sim%outfile_25,999) sim%naesmd%atomtype(k),sim%deriv_forces_state(n,1+3*(k-1)) &
                    ,sim%deriv_forces_state(n,2+3*(k-1)),sim%deriv_forces_state(n,3+3*(k-1))
            end do
            end do
            if(sim%naesmd%dynam_type.eq.'aimc') write(sim%outfile_26,*) '  time = ',sim%naesmd%tfemto, ' weight = ', sim%aimc%Weight 
        endif
        call flush(sim%outfile_25)
        if (sim%naesmd%dynam_type.eq.'aimc') call flush(sim%outfile_26)
        !Old Force writing code
        if((lprint.ge.3).or.(sim%naesmd%dynam_type.eq.'aimc')) then
            write (sim%outfile_12,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(sim%outfile_12,999) sim%naesmd%atomtype(k),sim%deriv_forces(1+3*(k-1)) &
                    ,sim%deriv_forces(2+3*(k-1)),sim%deriv_forces(3+3*(k-1))
            end do
            call flush(sim%outfile_12)
        endif
        if(lprint.ge.3) then
            ! Check the position of the center of mass
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0
        
            do j=1,sim%naesmd%natom
                xcm=xcm+sim%naesmd%rx(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                ycm=ycm+sim%naesmd%ry(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                zcm=zcm+sim%naesmd%rz(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            end do
        
            write(sim%outfile_17,889) sim%naesmd%tfemto,xcm-sim%naesmd%xcmini,ycm-sim%naesmd%ycmini,zcm-sim%naesmd%zcmini
            call flush(sim%outfile_17)
            call flush(sim%outfile_18)
        end if

        if((sim%naesmd%state.eq.'exct').and.(lprint.ge.1).and. &
           (sim%naesmd%dynam_type.ne.'mf').and.(sim%naesmd%dynam_type.ne.'aimc')) then
            if(sim%naesmd%printTDM.eq.0) then
                write(sim%outfile_5,889) sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
            else
                write(sim%outfile_5,889) sim%naesmd%tfemto,(sim%dav%v2(j,sim%naesmd%ihop),j=1,sim%dav%Nb**2)
            endif
            call flush(sim%outfile_5)
        else if((sim%naesmd%state.eq.'exct').and.(lprint.ge.1).and. &
                ((sim%naesmd%dynam_type.eq.'mf').or.(sim%naesmd%dynam_type.eq.'aimc'))) then
            do k=1,sim%excN
                sim%naesmd%ihop=k
                if(sim%naesmd%printTDM.eq.0) then
                    write(sim%outfile_5,999) k,sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
                else
                    write(sim%outfile_5,999) k,sim%naesmd%tfemto,(sim%dav%v2(j,sim%naesmd%ihop),j=1,sim%dav%Nb**2)
                endif
                call flush(sim%outfile_5)
            enddo
            ! in order to print the initial transition density of all states
            !Unopened file
        end if

        
        !Determine file status based on job sim%naesmd%state
        if(sim%naesmd%tfemto.eq.0.d0) then
            file_access='sequential'
            file_status='unknown'
        else
            file_access='append'
            file_status='old'
        end if

        
        if(sim%naesmd%tfemto.eq.0.d0) then
        
            write (sim%outfile_20,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            write(sim%outfile_20,557) '$VELOC'
            
            do k=1,sim%naesmd%natom
                write(sim%outfile_20,223) sim%naesmd%vx(k)*convl/convt, &
                    sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt
            end do
            
            write(sim%outfile_20,556) '$ENDVELOC'
            
            if(ibo==0) then
                write (sim%outfile_21,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                    '  time = ',sim%naesmd%tfemto
                write(sim%outfile_21,557) '$COEFF'
                
                do k=1,sim%excN
                    write(sim%outfile_21,223) yg(k)*cos(yg(k+sim%excN)),yg(k)*sin(yg(k+sim%excN)),yg_new(k),yg_new(k+sim%excN),dmod(yg_new(k+2*sim%excN),TWOPI)
                enddo
                
                write(sim%outfile_21,556) '$ENDCOEFF'
            end if
            
            write (sim%outfile_22,*) sim%naesmd%natom
            write (sim%outfile_22,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(sim%outfile_22,302) ELEMNT(sim%naesmd%atomtype(k)),sim%naesmd%rx(k)*convl, &
                    sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
            end do
        end if


        if(sim%naesmd%iview.eq.1) then

            ! to be used in case we want to print the transition densities of all the states at t=0
            do kki=1,sim%excN
            if(sim%Nsim.eq.1) then 
                card='view' // sim%naesmd%ktbig(sim%naesmd%icontini) // '-' //  sim%naesmd%ktbig(kki) // '.DATA'
            else 
                    write (tmpcard, "(I4.4)") sim%id 
                    card='view' // sim%naesmd%ktbig(sim%naesmd%icontini) // '-' //  &
                          & sim%naesmd%ktbig(kki) // '_' // trim(tmpcard) // '.DATA'
            end if 
            unitcard=90*10000+sim%id
                !************************************************************************************
                !       card='view' // sim%naesmd%ktbig(sim%naesmd%icontini) // '-' //  sim%naesmd%ktbig(sim%naesmd%ihop) // '.DATA'
                !************************************************************************************
                OPEN(unitcard,FILE=card)
                write(unitcard,440) ' Number of atoms:'
                write(unitcard,99) sim%naesmd%natom
                write(unitcard,441) ' Number of orbitals:'
                write(unitcard,99) sim%naesmd%nbasis
                write(unitcard,445) ' Number of occupied orbitals:'
                write(unitcard,222) ' 1'
                write(unitcard,442) ' Number of eigenvectors printed:'
                write(unitcard,222) ' 1'
                write(unitcard,441) ' Atomic coordinates:'

                do k=1,sim%naesmd%natom
                    write(unitcard,999) sim%naesmd%atomtype(k), &
                        sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
                end do

                write(unitcard,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
                do k=1,sim%dav%Nb
                    ! to be used in case we want to print the transition densities of all the states at t=0
                    write(unitcard,*) sim%dav%v2(sim%dav%Nb*(k-1)+k,kki)
                end do
                write(unitcard,*) ""
                close(unitcard)
            ! to be used in case we want to print the transition densities of all the states at t=0
            end do
        end if

        call flush(sim%outfile_1)

        if(sim%naesmd%state.eq.'exct'.and.lprint.ge.2) then
            write(sim%outfile_11,688) sim%naesmd%tfemto,(sim%naesmd%iorden(j),j=1,sim%excN)
            call flush(sim%outfile_11)
        end if


222     FORMAT(A2,3(1x,F12.6))
223     FORMAT(6(1x,F16.10))
688     FORMAT(F18.10,10000(1X,I4))
440     format(a17)
441     format(a20)
442     format(a32)
443     format(a41)
445     format(a29)
449     format(a28,F18.10,a9,F18.10)
889     FORMAT(100000(1X,F18.10))
999     FORMAT(I3,1X,100000(1X,F18.10))
99      FORMAT(I5)
555     format(a90)
556     format(a9)
557     format(a6)
600     format(a7,a4)
302     format(A3,3f16.10)
888     FORMAT(30000(1X,F18.10)) 

        return
    end subroutine

    !***********************************************
    ! write output information
    ! in the classical loop
    !***********************************************

    subroutine writeoutput(sim,ibo,yg,yg_new,lprint,cross,Nsim)
        use naesmd_constants
        use naesmd_module
        use langevin_temperature
        use communism
        use constants
        implicit none
 
        type(simulation_t),pointer::sim
        INTEGER j,k,n,ibo,kki,lprint,Nsim
        _REAL_ ntot
        _REAL_ xcm,ycm,zcm
        _REAL_ vmf
        _REAL_ randre
        integer randint
        integer cross(sim%excN)
        character*1000 card, tmpcard
        _REAL_ yg(2*sim%excN)
        _REAL_ yg_new(3*sim%excN)
        _REAL_ energy !JAB Test
        integer unitcard
        if(sim%cosmo%solvent_model.eq.2) then
            call calc_excsolven(sim%cosmo,sim%dav,sim%qmmm,energy) !JAB Test
            sim%naesmd%vmdqt(sim%naesmd%ihop)=sim%naesmd%vmdqt(sim%naesmd%ihop)-0.5*energy/feVmdqt !JAB Test
        endif

        if(sim%naesmd%state.eq.'fund') then
            write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%vgs*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt,&
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt-sim%naesmd%etotini*feVmdqt
        end if

        if(sim%naesmd%state.eq.'exct') then
          if(sim%naesmd%dynam_type.eq.'tully') then
            if(ibo.eq.1) then
                write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,&
                    sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                    sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt, &
                    sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                    sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt,sim%naesmd%kin*feVmdqt &
                    +sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt &
                    -sim%naesmd%etotini*feVmdqt
            else
                write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,&
                    sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                    sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt, &
                    sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                    sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt,sim%naesmd%kin*feVmdqt &
                    +sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt &
                    -sim%naesmd%etotini*feVmdqt
            end if
          elseif(sim%naesmd%dynam_type.eq.'aimc'.or.sim%naesmd%dynam_type.eq.'mf') then
            vmf=0.0d0
            do j=1,sim%excN
              vmf=vmf+yg(j)**2*sim%naesmd%vmdqtnew(j)
            enddo
            write(sim%outfile_1,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,&
                    sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                    vmf*feVmdqt, &
                    vmf*feVmdqt-sim%naesmd%vmfini*feVmdqt, &
                    sim%naesmd%kin*feVmdqt+vmf*feVmdqt,sim%naesmd%kin*feVmdqt &
                    +vmf*feVmdqt &
                    -sim%naesmd%etotini*feVmdqt 
          endif
        end if

        !ntot is the variable to check the norm conservation
        if(sim%naesmd%state.eq.'exct') then
            ntot=0
            do j=1,sim%excN
                ntot=ntot+yg(j)**2
            end do

            if(lprint.ge.1) then
                if(ibo.eq.1) then
                    write(sim%outfile_4,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                        (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                else
                    write(sim%outfile_4,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                        (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                    write(sim%outfile_8,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot
                    write(sim%outfile_7,888) sim%naesmd%tfemto,((sim%naesmd%cadiabnew(j,k),k=1,sim%excN),j=1,sim%excN)

                    call flush(sim%outfile_8)
                    call flush(sim%outfile_7)
                end if
                call flush(sim%outfile_4)
            elseif(sim%naesmd%dynam_type.eq.'aimc') then
                write(sim%outfile_8,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot
            end if

            if(lprint.ge.3.and.ibo.ne.1) then
                write(sim%outfile_16,889) sim%naesmd%tfemto,(dmod(yg(j+sim%excN),TWOPI),j=1,sim%excN)
                call flush(sim%outfile_16)
            end if
        else
            if(sim%excN.gt.0.and.lprint.ge.1) then
                write(sim%outfile_4,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                        (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
            endif
        end if

        write(sim%outfile_2,889) sim%naesmd%tfemto,sim%naesmd%tempi,sim%naesmd%tempf
        call flush(sim%outfile_2)
        !
        !
        !State dependent forces needed for Multiconfigurational post-processing
        if(sim%naesmd%dynam_type.eq.'aimc'.or.sim%naesmd%dynam_type.eq.'mf') then
            write (sim%outfile_25,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do n=1,sim%excN
            do k=1,sim%naesmd%natom
                write(sim%outfile_25,999) sim%naesmd%atomtype(k),sim%deriv_forces_state(n,1+3*(k-1)) &
                    ,sim%deriv_forces_state(n,2+3*(k-1)),sim%deriv_forces_state(n,3+3*(k-1))
            end do
            end do
            if(sim%naesmd%dynam_type.eq.'aimc') write(sim%outfile_26,*) '  time = ',sim%naesmd%tfemto, ' weight = ', sim%aimc%Weight 
        endif
        call flush(sim%outfile_25)
        if(sim%naesmd%dynam_type.eq.'aimc') call flush(sim%outfile_26)
        if((lprint.ge.3).or.(sim%naesmd%dynam_type.eq.'aimc')) then
            write (sim%outfile_12,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(sim%outfile_12,999) sim%naesmd%atomtype(k),sim%deriv_forces(1+3*(k-1)) &
                    ,sim%deriv_forces(2+3*(k-1)),sim%deriv_forces(3+3*(k-1))
            end do
            call flush(sim%outfile_12)
        endif
        if(lprint.ge.3) then
            ! Check the position of the center of mass
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0

            do j=1,sim%naesmd%natom
                xcm=xcm+sim%naesmd%rx(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                ycm=ycm+sim%naesmd%ry(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                zcm=zcm+sim%naesmd%rz(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            end do

            write(sim%outfile_17,889) sim%naesmd%tfemto,xcm-sim%naesmd%xcmini,ycm-sim%naesmd%ycmini,zcm-sim%naesmd%zcmini
            call flush(sim%outfile_17)
            call flush(sim%outfile_18)
            call flush(sim%outfile_19)
        end if

        if(sim%naesmd%state.eq.'exct'.and.lprint.ge.2) then
            write(sim%outfile_11,688) sim%naesmd%tfemto,(sim%naesmd%iorden(j),j=1,sim%excN),cross
            call flush(sim%outfile_11)
        end if

        if((sim%naesmd%state.eq.'exct').and.(lprint.ge.1).and. &
           (sim%naesmd%dynam_type.ne.'mf').and.(sim%naesmd%dynam_type.ne.'aimc')) then
            if(sim%naesmd%printTDM.eq.0) then
                write(sim%outfile_5,889) sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
            else
                write(sim%outfile_5,889) sim%naesmd%tfemto,(sim%dav%v2(j,sim%naesmd%ihop),j=1,sim%dav%Nb**2)
            endif
            call flush(sim%outfile_5)
        else if((sim%naesmd%state.eq.'exct').and.(lprint.ge.1).and. &
                ((sim%naesmd%dynam_type.eq.'mf').or.(sim%naesmd%dynam_type.eq.'aimc'))) then
            do k=1,sim%excN
                sim%naesmd%ihop=k
                if(sim%naesmd%printTDM.eq.0) then
                    write(sim%outfile_5,999) k,sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
                else
                    write(sim%outfile_5,999) k,sim%naesmd%tfemto,(sim%dav%v2(j,sim%naesmd%ihop),j=1,sim%dav%Nb**2)
                endif
                call flush(sim%outfile_5)
            enddo
        end if

        !if(sim%naesmd%icont.ne.sim%naesmd%nstepcoord) then
        !    sim%naesmd%icont=sim%naesmd%icont+1
        !else
        !    sim%naesmd%icont=1
        !    sim%naesmd%icontpdb=sim%naesmd%icontpdb+1
            
            write (sim%outfile_20,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            write(sim%outfile_20,557) '$VELOC'
    
            do k=1,sim%naesmd%natom
                write(sim%outfile_20,223) sim%naesmd%vx(k)*convl/convt, &
                    sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt
            end do
    
            write(sim%outfile_20,556) '$ENDVELOC'
            
            if(ibo==0) then
                write (sim%outfile_21,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                    '  time = ',sim%naesmd%tfemto
                write(sim%outfile_21,557) '$COEFF'
                
                do k=1,sim%excN
                    write(sim%outfile_21,223) yg(k)*cos(yg(k+sim%excN)),yg(k)*sin(yg(k+sim%excN)),yg_new(k),yg_new(k+sim%excN),dmod(yg_new(k+2*sim%excN),TWOPI)
                enddo
                
                write(sim%outfile_21,556) '$ENDCOEFF'
            end if

            write (sim%outfile_22,*) sim%naesmd%natom
            write (sim%outfile_22,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(sim%outfile_22,302) ELEMNT(sim%naesmd%atomtype(k)),sim%naesmd%rx(k)*convl, &
                    sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
            end do

            if(sim%naesmd%iview.eq.1) then
                do kki=1,sim%excN
                    if(sim%Nsim.eq.1) then
                        card='view' // sim%naesmd%ktbig(sim%naesmd%icontpdb) // '-' //  sim%naesmd%ktbig(kki) // '.DATA'
                    else 
                        write (tmpcard, "(I4.4)") sim%id 
                        card='view' // sim%naesmd%ktbig(sim%naesmd%icontpdb) // '-' //  &
                             & sim%naesmd%ktbig(kki) // '_' // trim(tmpcard) // '.DATA'
                    end if 
                    unitcard=90*10000+sim%id
                    OPEN(unitcard,FILE=card)
                    write(unitcard,440) ' Number of atoms:'
                    write(unitcard,99) sim%naesmd%natom
                    write(unitcard,441) ' Number of orbitals:'
                    write(unitcard,99) sim%naesmd%nbasis
                    write(unitcard,445) ' Number of occupied orbitals:'
                    write(unitcard,222) ' 1'
                    write(unitcard,442) ' Number of eigenvectors printed:'
                    write(unitcard,222) ' 1'
                    write(unitcard,441) ' Atomic coordinates:'

                    do k=1,sim%naesmd%natom
                        write(unitcard,999) sim%naesmd%atomtype(k), &
                            sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
                    end do

                    write(unitcard,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
                    do k=1,sim%dav%Nb
                        write(unitcard,*) sim%dav%v2(sim%dav%Nb*(k-1)+k,kki)
                    end do
                    write(unitcard,*) ""
                    close(unitcard)
                end do
            end if
        call flush(sim%outfile_1)
222     FORMAT(A2,3(1x,F12.6))
223     FORMAT(6(1x,F16.10))
440     format(a17)
441     format(a20)
442     format(a32)
443     format(a41)
445     format(a29)
889     FORMAT(100000(1X,F18.10))
999     FORMAT(I3,1X,100000(1X,F18.10))
99      FORMAT(I5)
555     format(a90)
556     format(a9)
557     format(a6)
302     format(A3,3f16.10)
449     format(a28,F18.10,a9,F18.10)
450     format(a8,I10)
688     FORMAT(F18.10,10000(1X,I4))
888     FORMAT(30000(1X,F18.10))
600     format(a7,a4)

        return
    end subroutine

!Suroutine for writing the restart file    
    subroutine write_restart(sim,Nsim,drops,Nsim_max)
        use naesmd_module
        use communism
        use constants
        
        implicit none

!Pased in
        type(simulation_t), target :: sim
        type(simulation_t), pointer :: simpoint
        integer, intent(in) :: Nsim
        integer, intent(in), dimension(:) :: drops(Nsim_max)
        integer, intent(in) :: Nsim_max
!Internal 
        integer unitcard
        character*150 card
        character*150 card_tmp
        character*150 line
        integer i,inputSteps,excN,outData
        character*150 dynam_type

        simpoint => sim

        card_tmp = ''
        card = ''

        unitcard = 10*1000 + simpoint%id
        if(Nsim.eq.1) then
            card_tmp='restart.tmp'
            card='restart.out'
        else
            write (card_tmp, "(A8,I4.4,A4)") "restart_", simpoint%id + sum(drops(:simpoint%id + 1)), ".tmp"
            write (card, "(A8,I4.4,A4)") "restart_", simpoint%id + sum(drops(:simpoint%id + 1)), ".out"
        endif

        OPEN(unitcard,FILE=card_tmp(1:16),ACTION='write',STATUS='replace')
        do i=1,size(simpoint%input_line)
            if(simpoint%input_line(i)(1:18).eq.'   exc_state_init=') then
                write(unitcard,"(A18,I2,A48)") '   exc_state_init=',simpoint%naesmd%ihop,', ! Initial excited state (0 - ground state) [0]'
            elseif(simpoint%input_line(i)(1:26).eq.'   n_exc_states_propagate=') then
                write(unitcard,"(A26,I3,A32)") '   n_exc_states_propagate=', simpoint%excN, ', ! Number of excited states [0]'
            elseif(simpoint%input_line(i)(1:13).eq.'   time_init=') then
                write(unitcard,"(A13,F9.2,A26)") '  time_init=',simpoint%naesmd%tfemto,', ! Initial time, fs [0.0]'
            elseif(simpoint%input_line(i)(1:17).eq.'   n_class_steps=') then
                write(unitcard,"(A17,I7,A33)") '  n_class_steps=',simpoint%naesmd%nstep0-nint(simpoint%naesmd%tfemto/(simpoint%naesmd%dtmdqt*convtf)),', ! Number of classical steps [1]'
            elseif(simpoint%input_line(i)(1:12).eq.'   nclones0=') then
                write(unitcard,"(A12,I2,A76)") '  nclones0=',simpoint%aimc%nclones,", ! Clones count for 'aimc' (must be declared here for restarting 'aimc') [0]"
            else
                write(unitcard,'(A)') simpoint%input_line(i)(1:len_trim(simpoint%input_line(i)))
            endif
        enddo
        do i=1,simpoint%naesmd%natom
            write(unitcard,"(I4,3F16.10)")  simpoint%naesmd%atomtype(i),simpoint%naesmd%rx(i)*convl,simpoint%naesmd%ry(i)*convl,sim%naesmd%rz(i)*convl
        enddo
        write(unitcard,'(A)') '&endcoord\n'
        write(unitcard,'(A)') '&veloc'
        do i=1,simpoint%naesmd%natom
            write(unitcard,"(3F16.10)") simpoint%naesmd%vx(i)*convl/convt,simpoint%naesmd%vy(i)*convl/convt,simpoint%naesmd%vz(i)*convl/convt
        enddo
        write(unitcard,'(A)') '&endveloc'
        write(unitcard,'(A)') '&coeff'
        do i=1,simpoint%excN
            write(unitcard,"(2F16.10)") simpoint%naesmd%yg(i), dmod(simpoint%naesmd%yg(i+simpoint%excN),TWOPI)
        enddo
        write(unitcard,'(A)') '&endcoeff'
        call flush(unitcard)
        close(unitcard)
        call system('mv '//card_tmp(1:16)//' '//card(1:16))

    endsubroutine write_restart

end module

