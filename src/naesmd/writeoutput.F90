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

    subroutine writeoutputini(sim,ibo,yg,lprint)
        use naesmd_constants
        use naesmd_module
        use communism
        implicit none

        type(simulation_t) sim
        integer  i,j,k,ibo,kki,lprint
        _REAL_ ntot
        _REAL_ xcm,ycm,zcm
        character*1000 card
        
        character*10 file_status
        character*10 file_access

        _REAL_ yg(:)
        _REAL_ energy !!JAB VE model

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
            write(98,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%vgs*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt-sim%naesmd%etotini*feVmdqt
        end if

        if(sim%naesmd%state.eq.'exct') then
            sim%naesmd%vini=sim%naesmd%vmdqt(sim%naesmd%ihop)
            sim%naesmd%etotini=sim%naesmd%kin+sim%naesmd%vmdqt(sim%naesmd%ihop)
            write(98,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt, &
                sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt,&
                sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt- &
                sim%naesmd%etotini*feVmdqt
        end if

        if(sim%naesmd%state.eq.'exct') then
            if(lprint.ge.1) then
                write(96,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                    (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                if(ibo.ne.1) then
                    write(95,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot
                end if
            end if

            if(lprint.ge.3.and.ibo.ne.1) then
                write(94,889) sim%naesmd%tfemto,(dsin(yg(j+sim%excN)),j=1,sim%excN)
                call flush(94)
            end if

            call flush(96)
            call flush(95)
        endif
        ! use for vibrations**************************************
        !      write(99,779) sim%naesmd%tfemto,(sim%naesmd%atomtype(k),
        !     $sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl,k=1,sim%naesmd%natom)
        !      write(80,779) sim%naesmd%tfemto,(sim%naesmd%atomtype(k),
        !     $sim%naesmd%vx(k)*convl/convt,sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt,k=1,sim%naesmd%natom)
        !***********************************************************

        !Old Force writing code
        if(lprint.ge.3) then
            write (85,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(85,999) sim%naesmd%atomtype(k),sim%deriv_forces(1+3*(k-1)) &
                    ,sim%deriv_forces(2+3*(k-1)),sim%deriv_forces(3+3*(k-1))
            end do
        !    write(85,889) sim%naesmd%tfemto,(sim%deriv_forces(1+3*(j-1)),j=1,sim%naesmd%natom)
        !    write(84,889) sim%naesmd%tfemto,(sim%deriv_forces(2+3*(j-1)),j=1,sim%naesmd%natom)
        !    write(83,889) sim%naesmd%tfemto,(sim%deriv_forces(3+3*(j-1)),j=1,sim%naesmd%natom)
            ! Check the position of the center of mass
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0
        
            do j=1,sim%naesmd%natom
                xcm=xcm+sim%naesmd%rx(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                ycm=ycm+sim%naesmd%ry(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                zcm=zcm+sim%naesmd%rz(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            end do
        
            write(91,889) sim%naesmd%tfemto,xcm-sim%naesmd%xcmini,ycm-sim%naesmd%ycmini,zcm-sim%naesmd%zcmini
            call flush(85)
            call flush(84)
            call flush(83)
            call flush(91)
            call flush(99)
        end if

        if(sim%naesmd%state.eq.'exct'.and.lprint.ge.1) then
            write(89,889) sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
            call flush(89)
            ! in order to print the initial transition density of all states
            do k=1,sim%excN
                write(77,889) sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,k),j=1,sim%dav%Nb)
                call flush(77)
            end do
        end if

        ! write the atomic positions to trajectory file
!        card='coordinates.out'
!        OPEN(201,FILE=card,action='write',STATUS='replace')
!        
!        write (201,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
!            '  time = ',sim%naesmd%tfemto
!        write(201,557) '$COORD'
!        do i=1,sim%naesmd%txtinicoord
!            write(201,555) sim%naesmd%txtinput(i)
!        end do
!
!        do i=1,sim%naesmd%natom
!            write(201,999) sim%naesmd%atomtype(i),sim%naesmd%rx(i)*convl &
!                ,sim%naesmd%ry(i)*convl,sim%naesmd%rz(i)*convl
!        end do
!
!        write(201,556) '$ENDCOORD'
        
        !Determine file status based on job sim%naesmd%state
        if(sim%naesmd%tfemto.eq.0.d0) then
            file_access='sequential'
            file_status='unknown'
        else
            file_access='append'
            file_status='old'
        end if
        
        OPEN(202,FILE='velocity.out',action='write',STATUS=file_status,ACCESS=file_access)
        if(ibo==0) OPEN(203,FILE='coefficient.out',action='write',STATUS=file_status,ACCESS=file_access)
        OPEN(9,file='coords.xyz',action='write',STATUS=file_status,ACCESS=file_access)
        
        if(sim%naesmd%tfemto.eq.0.d0) then
        
            write (202,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            write(202,557) '$VELOC'
            
            do k=1,sim%naesmd%natom
                write(202,223) sim%naesmd%vx(k)*convl/convt, &
                    sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt
            end do
            
            write(202,556) '$ENDVELOC'
            
            if(ibo==0) then
                write (203,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                    '  time = ',sim%naesmd%tfemto
                write(203,557) '$COEFF'
                
                do k=1,sim%excN
                    write(203,223) yg(k)**2,dsin(yg(k+sim%excN))
                enddo
                
                write(203,556) '$ENDCOEFF'
            end if
            
            write (9,*) sim%naesmd%natom
            write (9,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(9,302) ELEMNT(sim%naesmd%atomtype(k)),sim%naesmd%rx(k)*convl, &
                    sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
            end do
        end if


        if(sim%naesmd%iview.eq.1) then

            ! to be used in case we want to print the transition densities of all the states at t=0
            do kki=1,sim%excN
                card='view' // sim%naesmd%ktbig(sim%naesmd%icontini) // '-' //  sim%naesmd%ktbig(kki) // '.DATA'
                !************************************************************************************
                !       card='view' // sim%naesmd%ktbig(sim%naesmd%icontini) // '-' //  sim%naesmd%ktbig(sim%naesmd%ihop) // '.DATA'
                !************************************************************************************
                OPEN(90,FILE=card)
                write(90,440) ' Number of atoms:'
                write(90,99) sim%naesmd%natom
                write(90,441) ' Number of orbitals:'
                write(90,99) sim%naesmd%nao
                write(90,445) ' Number of occupied orbitals:'
                write(90,222) ' 1'
                write(90,442) ' Number of eigenvectors printed:'
                write(90,222) ' 1'
                write(90,441) ' Atomic coordinates:'

                do k=1,sim%naesmd%natom
                    write(90,999) sim%naesmd%atomtype(k), &
                        sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
                end do

                write(90,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
                do k=1,sim%dav%Nb
                    ! to be used in case we want to print the transition densities of all the states at t=0
                    write(90,*) sim%dav%v2(sim%dav%Nb*(k-1)+k,kki)
                end do
                close(90)
            ! to be used in case we want to print the transition densities of all the states at t=0
            end do
        end if

        call flush(98)

        !111   FORMAT(A2,14x,3(1x,F16.12,1x,I1))
222     FORMAT(A2,3(1x,F12.6))
223     FORMAT(3(1x,F16.10))
        !333   format(a5,I3,1X,a80)
440     format(a17)
441     format(a20)
442     format(a32)
443     format(a41)
        !444   format(a80)
445     format(a29)
449     format(a28,F18.10,a9,F18.10)
889     FORMAT(10000(1X,F18.10))
        !887   FORMAT(F18.10,10000(1X,I1))
999     FORMAT(I3,1X,1000(1X,F18.10))
        !998   FORMAT(I4)
99      FORMAT(I5)
        !777   FORMAT(A7,I4,2X,A2,2X,A5,I4,4X,3(1X,F7.3))
        !778   format(1000000(f10.6,1x))
        !779   format(F18.10,10000(1X,I3,3(1X,F18.10)))
555     format(a90)
556     format(a9)
557     format(a6)
        !558   format(a7)
        !88    format(a1)
600     format(a7,a4)
302     format(A3,3f16.10)

        return
    end subroutine

    !***********************************************
    ! write output information
    ! in the classical loop
    !***********************************************

    subroutine writeoutput(sim,ibo,yg,lprint,cross)
        use naesmd_constants
        use naesmd_module
        use langevin_temperature
        use communism
        implicit none
 
        type(simulation_t),pointer::sim
        INTEGER j,k,ibo,kki,lprint
        _REAL_ ntot
        _REAL_ xcm,ycm,zcm
        _REAL_ randre
        integer randint
        integer cross(sim%excN)
        character*1000 card
        _REAL_ yg(2*sim%excN)
        _REAL_ energy !JAB Test
        logical first
        data first /.true./
        save first

        if(sim%cosmo%solvent_model.eq.2) then
            call calc_excsolven(sim%cosmo,sim%dav,sim%qmmm,energy) !JAB Test
            sim%naesmd%vmdqt(sim%naesmd%ihop)=sim%naesmd%vmdqt(sim%naesmd%ihop)-0.5*energy/feVmdqt !JAB Test
        endif

        if(sim%naesmd%state.eq.'fund') then
            write(98,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                sim%naesmd%vgs*feVmdqt, &
                sim%naesmd%vgs*feVmdqt-sim%naesmd%vini*feVmdqt, &
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt,&
                sim%naesmd%kin*feVmdqt+sim%naesmd%vgs*feVmdqt-sim%naesmd%etotini*feVmdqt
        end if

        if(sim%naesmd%state.eq.'exct') then
            if(ibo.eq.1) then
                write(98,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                    sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt, &
                    sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                    sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt,sim%naesmd%kin*feVmdqt &
                    +sim%naesmd%vmdqt(sim%naesmd%ihop)*feVmdqt &
                    -sim%naesmd%etotini*feVmdqt
            else
                write(98,889) sim%naesmd%tfemto,sim%naesmd%kin*feVmdqt,sim%naesmd%kin*feVmdqt-sim%naesmd%kinini*feVmdqt, &
                    sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt, &
                    sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt-sim%naesmd%vini*feVmdqt, &
                    sim%naesmd%kin*feVmdqt+sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt,sim%naesmd%kin*feVmdqt &
                    +sim%naesmd%vmdqtnew(sim%naesmd%ihop)*feVmdqt &
                    -sim%naesmd%etotini*feVmdqt

            end if
        end if

        !ntot is the variable to check the norm conservation
        if(sim%naesmd%state.eq.'exct') then
            ntot=0
            do j=1,sim%excN
                ntot=ntot+yg(j)**2
            end do

            if(lprint.ge.1) then
                if(ibo.eq.1) then
                    write(96,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                        (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                else
                    write(96,889) sim%naesmd%tfemto,sim%naesmd%vgs*feVmdqt, &
                        (sim%naesmd%vmdqt(j)*feVmdqt,j=1,sim%excN)
                    write(95,999) sim%naesmd%ihop,sim%naesmd%tfemto,(yg(j)**2,j=1,sim%excN),ntot
                    write(93,888) sim%naesmd%tfemto,((sim%naesmd%cadiabnew(j,k),k=1,sim%excN),j=1,sim%excN)

                    call flush(95)
                    call flush(93)
                end if

                call flush(96)
            end if

            if(lprint.ge.3.and.ibo.ne.1) then
                write(94,889) sim%naesmd%tfemto,(dsin(yg(j+sim%excN)),j=1,sim%excN)
                call flush(94)
            end if
        end if

        write(92,889) sim%naesmd%tfemto,sim%naesmd%tempi,sim%naesmd%tempf
        call flush(92)
        !
        ! use for vibrations*****************************
        !      write(99,779) sim%naesmd%tfemto,(sim%naesmd%atomtype(k),
        !     $sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl,k=1,sim%naesmd%natom)
        !      write(80,779) sim%naesmd%tfemto,(sim%naesmd%atomtype(k),
        !     $sim%naesmd%vx(k)*convl/convt,sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt,k=1,sim%naesmd%natom)
        !******************************************************
        !
        if(lprint.ge.3) then
            write (85,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(85,999) sim%naesmd%atomtype(k),sim%deriv_forces(1+3*(k-1)) &
                    ,sim%deriv_forces(2+3*(k-1)),sim%deriv_forces(3+3*(k-1))
            end do
!            write(85,889) sim%naesmd%tfemto,(sim%deriv_forces(1+3*(j-1)),j=1,sim%naesmd%natom)
!            write(84,889) sim%naesmd%tfemto,(sim%deriv_forces(2+3*(j-1)),j=1,sim%naesmd%natom)
!            write(83,889) sim%naesmd%tfemto,(sim%deriv_forces(3+3*(j-1)),j=1,sim%naesmd%natom)

            write(125,889) sim%naesmd%tfemto,(sim%naesmd%ax(j),j=1,sim%naesmd%natom)
            write(126,889) sim%naesmd%tfemto,(sim%naesmd%ax(j),j=1,sim%naesmd%natom)
            ! Check the position of the center of mass
            xcm=0.0d0
            ycm=0.0d0
            zcm=0.0d0

            do j=1,sim%naesmd%natom
                xcm=xcm+sim%naesmd%rx(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                ycm=ycm+sim%naesmd%ry(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
                zcm=zcm+sim%naesmd%rz(j)*sim%naesmd%massmdqt(j)/sim%naesmd%masstot
            end do

            write(91,889) sim%naesmd%tfemto,xcm-sim%naesmd%xcmini,ycm-sim%naesmd%ycmini,zcm-sim%naesmd%zcmini
            call flush(85)
            call flush(84)
            call flush(83)
            call flush(91)
            call flush(99)
            call flush(80)
        end if

        if(sim%naesmd%state.eq.'exct'.and.lprint.ge.2) then
            write(100,688) sim%naesmd%tfemto,(sim%naesmd%iorden(j),j=1,sim%excN),cross
            write(120,688) sim%naesmd%tfemto,(cross(j),j=1,sim%excN)
            call flush(120)
            call flush(100)
        end if

        if(sim%naesmd%state.eq.'exct'.and.lprint.ge.1) then
            write(89,889) sim%naesmd%tfemto,(sim%dav%v2(sim%dav%Nb*(j-1)+j,sim%naesmd%ihop),j=1,sim%dav%Nb)
            call flush(89)
        end if

        if(sim%naesmd%icont.ne.sim%naesmd%nstepcoord) then
            sim%naesmd%icont=sim%naesmd%icont+1
        else
            sim%naesmd%icont=1
            sim%naesmd%icontpdb=sim%naesmd%icontpdb+1
            
!            write (201,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
!                '  time = ',sim%naesmd%tfemto
!            write(201,557) '$COORD'
!            do k=1,sim%naesmd%txtinicoord
!                write(201,555) sim%naesmd%txtinput(k)
!            end do
!    
!            do k=1,sim%naesmd%natom
!                write(201,999) sim%naesmd%atomtype(k),sim%naesmd%rx(k)*convl &
!                    ,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
!            end do
!    
!            write(201,556) '$ENDCOORD'
            
            
            write (202,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            write(202,557) '$VELOC'
    
            do k=1,sim%naesmd%natom
                write(202,223) sim%naesmd%vx(k)*convl/convt, &
                    sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt
            end do
    
            write(202,556) '$ENDVELOC'
            
            if(ibo==0) then
                write (203,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                    '  time = ',sim%naesmd%tfemto
                write(203,557) '$COEFF'
                
                do k=1,sim%excN
                    write(203,223) yg(k)**2,dsin(yg(k+sim%excN))
                enddo
                
                write(203,556) '$ENDCOEFF'
            end if
    
            card='restart.out'
            OPEN(10,FILE=card,ACTION='write',STATUS='replace')
            write (10,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            randint=int(rranf1(sim%naesmd%iseedmdqt)*1d8)
            write (10,*) 'State = ', sim%naesmd%ihop
            write (10,*) 'Seed  = ', randint
            write(10,557) '$COORD'
            do k=1,sim%naesmd%txtinicoord
                write(10,555) sim%naesmd%txtinput(k)
            end do
    
            do k=1,sim%naesmd%natom
                write(10,999) sim%naesmd%atomtype(k),sim%naesmd%rx(k)*convl &
                    ,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
            end do
    
            write(10,556) '$ENDCOORD'
            write(10,557) '$VELOC'
    
            do k=1,sim%naesmd%natom
                write(10,223) sim%naesmd%vx(k)*convl/convt, &
                    sim%naesmd%vy(k)*convl/convt,sim%naesmd%vz(k)*convl/convt
            end do
    
            write(10,556) '$ENDVELOC'
            write(10,557) '$COEFF'
    
            do k=1,sim%excN
                write(10,223) yg(k)**2,dsin(yg(k+sim%excN))
            enddo
    
            write(10,556) '$ENDCOEFF'
            close(10)
            

            write (9,*) sim%naesmd%natom
            write (9,449) 'FINAL HEAT OF FORMATION =   ', (sim%naesmd%kin+sim%naesmd%vgs)*feVmdqt, &
                '  time = ',sim%naesmd%tfemto
            do k=1,sim%naesmd%natom
                write(9,302) ELEMNT(sim%naesmd%atomtype(k)),sim%naesmd%rx(k)*convl, &
                    sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
            end do

            if(sim%naesmd%iview.eq.1) then
                do kki=1,sim%excN
                    card='view' // sim%naesmd%ktbig(sim%naesmd%icontpdb) // '-' //  sim%naesmd%ktbig(kki) // '.DATA'
                    OPEN(90,FILE=card)
                    write(90,440) ' Number of atoms:'
                    write(90,99) sim%naesmd%natom
                    write(90,441) ' Number of orbitals:'
                    write(90,99) sim%naesmd%nao
                    write(90,445) ' Number of occupied orbitals:'
                    write(90,222) ' 1'
                    write(90,442) ' Number of eigenvectors printed:'
                    write(90,222) ' 1'
                    write(90,441) ' Atomic coordinates:'

                    do k=1,sim%naesmd%natom
                        write(90,999) sim%naesmd%atomtype(k), &
                            sim%naesmd%rx(k)*convl,sim%naesmd%ry(k)*convl,sim%naesmd%rz(k)*convl
                    end do

                    write(90,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
                    do k=1,sim%dav%Nb
                        write(90,*) sim%dav%v2(sim%dav%Nb*(k-1)+k,kki)
                    end do
                    close(90)
                end do
            end if
        endif
        call flush(98)
        !111   FORMAT(A2,14x,3(1x,F16.12,1x,I1))
222     FORMAT(A2,3(1x,F12.6))
223     FORMAT(3(1x,F16.10))
        !333   format(a5,I3,1X,a80)
440     format(a17)
441     format(a20)
442     format(a32)
443     format(a41)
        !444   format(a80)
445     format(a29)
889     FORMAT(10000(1X,F18.10))
        !887   FORMAT(F18.10,10000(1X,I1))
999     FORMAT(I3,1X,1000(1X,F18.10))
        !998   FORMAT(I4)
99      FORMAT(I5)
        !777   FORMAT(A7,I4,2X,A2,2X,A5,I4,4X,3(1X,F7.3))
        !778   format(1000000(f10.6,1x))
        !779   format(F18.10,10000(1X,I3,3(1X,F18.10)))
555     format(a90)
556     format(a9)
557     format(a6)
        !558   format(a7)
        !88    format(a1)
302     format(A3,3f16.10)
        !333   format(a5,I3,1X,a80)
        !440   format(a17)
        !441   format(a20)
        !442   format(a32)
        !443   format(a41)
        !444   format(a80)
        !445   format(a29)
449     format(a28,F18.10,a9,F18.10)
450     format(a8,I10)
688     FORMAT(F18.10,10000(1X,I4))
        !889   FORMAT(20000(1X,F18.10))
888     FORMAT(30000(1X,F18.10))
        !887   FORMAT(F18.10,10000(1X,I1))
        !999   FORMAT(I3,1X,1000(1X,F18.10))
        !998   FORMAT(I4)
        !99    FORMAT(I5)
        !777   FORMAT(A7,I4,2X,A2,2X,A5,I4,4X,3(1X,F7.3))
        !778   format(1000000(f10.6,1x))
        !779   format(F18.10,10000(1X,I3,3(1X,F18.10)))
        !555   format(a90)
        !556   format(a9)
        !557   format(a6)
        !558   format(a7)
        !88    format(a1)
600     format(a7,a4)

        return
    end subroutine
end module

