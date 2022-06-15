#include "dprec.fh"
#include "assert.fh"

module quantum_prop_add
    use naesmd_constants
    use communism
    use cadiab_module
    use additional_subroutines

contains

    subroutine checkcrossing(sim,cross,lprint)
        implicit none

        type(simulation_t),pointer::sim
        integer k,j,i,lprint
        integer cross(sim%excN)
        
        !TODO ascpr and z are for use with apc subroutine min-cost algorithm written
        !in F77. Should be updated to have an allocatable option.
        integer ascpr(260,260),z
        integer iordenapc(260)

        !--------------------------------------------------------------------
        !
        !  Following the nonavoiding crossing of states
        !
        !  sim%naesmd%scpr is the overlap between new and old transition densities
        !  if there is no crossing, the matrix is supposed to be close
        !  to identity matrix
        !
        !--------------------------------------------------------------------
        if (sim%excN>260) then
            write(6,*)'Error:quantum propagator does not yet handle more than 260 states'
        endif

        sim%naesmd%scpr(1:sim%excN,1:sim%excN)=0.d0

        ! kav: FIXME
        ! Here and below, in the old code the size Ncis is used as a size of
        ! sim%naesmd%cmdqt vectors. This is not entirely correct, since the size is
        ! Nrpa=2*Ncis (check it!!!). However, it can be very accurate becuase
        ! the contribution of negative energies to a positive excitation
        ! energy can be very small (would be great to check it)
        do i=1,sim%excN
            do j=1,sim%excN
                do k=1,sim%dav%Ncis
                    sim%naesmd%scpr(i,j)=sim%naesmd%scpr(i,j)+sim%naesmd%cmdqtold(k,i)*sim%naesmd%cmdqtnew(k,j)
                end do
            end do
        end do

        do i=1,sim%excN
            do j=1,sim%excN
                ascpr(i,j)=int(sim%naesmd%scpr(i,j)**2*1.d5)
            end do
        end do
        ! window**************************************
        do i=1,sim%excN
            do j=1,sim%excN
                if((j.lt.(i-2)).or.(j.gt.(i+2))) then
                    ascpr(i,j)=-1*ifix(sngl(1.d5))
                end if
            end do
        end do

        !************************************

        do i=1,sim%excN
            do j=1,sim%excN
                ascpr(i,j)=-1*ascpr(i,j)
            enddo
        end do
        
        !Translate to old APC min-cost algorithm routine in F77
        !Requires static matrix sizes
        iordenapc(1:sim%excN)=sim%naesmd%iorden
        call apc(sim%naesmd%apc_common,sim%excN,ascpr,iordenapc,z)
        sim%naesmd%iorden=iordenapc(1:sim%excN)

        do i=1,sim%excN
            sim%naesmd%iordenhop(i)=sim%naesmd%iorden(i)
            if(sim%naesmd%iorden(i).ne.i) then
                if(i.lt.sim%naesmd%iorden(i).or.i.eq.sim%naesmd%ihop) then
                    if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).lt.0.9d0) then
                        cross(i)=1
                    else
                        cross(i)=2
                    end if
                    if(lprint.ge.2) then
                    endif
                else
                    cross(i)=0
                end if
            else
                cross(i)=0
            end if
        end do


        ! Store the overlap matrix to print it if hops or cross

        do i=1,sim%excN
            do j=1,sim%excN
                sim%naesmd%scprreal(i,j)=sim%naesmd%scpr(i,j)
            end do
        end do

        !***********************************


        return
    end subroutine

    subroutine checkcrossingmiddle(sim,cross)
        implicit none

        type(simulation_t),pointer::sim
        integer k,j,i
        integer cross(sim%excN)

        !***************************************************
        ! following the nonavoiding crossing of states
        !***************************************************
        do i=1,sim%excN
            do j=1,sim%excN
                sim%naesmd%scpr(i,j)=0.0d0
                do k=1,sim%dav%Ncis
                    sim%naesmd%scpr(i,j)=sim%naesmd%scpr(i,j)+sim%naesmd%cmdqtmiddleold(k,i)*sim%naesmd%cmdqtmiddle(k,j)
                end do
            end do
        end do
        do i=1,sim%excN
            if((i.lt.sim%naesmd%iorden(i)).or.(i.eq.sim%naesmd%ihop)) then
                    if(sim%naesmd%ihop.gt.0) then
                         if(i.ne.sim%naesmd%iorden(sim%naesmd%ihop)) then
                             if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).ge.0.9d0) then
                                 cross(i)=2
                             end if
                         end if
                    else
                         if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).ge.0.9d0) then
                             cross(i)=2
                         end if
                    end if
            end if
        end do
        return
    end subroutine

    subroutine vmdqtmiddlecalc(sim, iimdqt,Na)
        implicit none

        type(simulation_t), pointer :: sim
        integer i,j,iimdqt
        integer Na
        _REAL_ xx(Na),yy(Na),zz(Na)

        if(iimdqt.eq.1) then
            do i=1,sim%excN
                sim%naesmd%vmdqtmiddleold(i)=sim%naesmd%vmdqtold(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    sim%naesmd%cmdqtmiddleold(i,j)=sim%naesmd%cmdqtold(i,j)
                end do
            end do

        else
            do i=1,sim%excN
                sim%naesmd%vmdqtmiddleold(i)=sim%naesmd%vmdqtmiddle(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    sim%naesmd%cmdqtmiddleold(i,j)=sim%naesmd%cmdqtmiddle(i,j)
                end do
            end do
        end if

        if(iimdqt.eq.sim%naesmd%nquantumstep) then

            do i=1,sim%excN
                sim%naesmd%vmdqtmiddle(i)=sim%naesmd%vmdqtnew(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    sim%naesmd%cmdqtmiddle(i,j)=sim%naesmd%cmdqtnew(i,j)
                end do
            end do

        else

            if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
                do j=1,sim%naesmd%natom
                    xx(j)=sim%naesmd%rxold(j) &
                        +sim%naesmd%vxold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%axold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2

                    yy(j)=sim%naesmd%ryold(j) &
                        +sim%naesmd%vyold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%ayold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2

                    zz(j)=sim%naesmd%rzold(j) &
                        +sim%naesmd%vzold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%azold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2
                end do
            end if

            if(sim%naesmd%ensemble.eq.'langev') then
                do j=1,sim%naesmd%natom

                    xx(j)=sim%naesmd%rxold(j) &
                        +sim%naesmd%vxold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%axold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        +sim%naesmd%prand(1,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)

                    yy(j)=sim%naesmd%ryold(j) &
                        +sim%naesmd%vyold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%ayold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        +sim%naesmd%prand(2,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)

                    zz(j)=sim%naesmd%rzold(j) &
                        +sim%naesmd%vzold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        +sim%naesmd%azold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        +sim%naesmd%prand(3,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)
                end do
            end if


            call do_sqm_davidson_update(sim,0,cmdqt=sim%naesmd%cmdqtmiddle, &
                vmdqt=sim%naesmd%vmdqtmiddle,vgs=sim%naesmd%vgs,rx=xx,ry=yy,rz=zz)
        end if
        write(6,*)'sim%naesmd%iordenhop:',sim%naesmd%iordenhop
        write(6,*)'sizevmdqtmiddle',shape(sim%naesmd%vmdqtmiddle)
! option 1: check trivial unavoided crossingfor all the states
! Note: Currently broken due to sim%naesmd%iordenhop initialized as 0
        do j=1,sim%excN
            if(sim%naesmd%lowvalue(j).gt.dabs(sim%naesmd%vmdqtmiddle(j)- &
                sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))) then

                sim%naesmd%lowvalue(j)=dabs(sim%naesmd%vmdqtmiddle(j)- &
                    sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))
                sim%naesmd%lowvaluestep(j)=iimdqt
            end if
        end do
! option 2: check trivial unavoided crossingfor for the current sim%naesmd%state
        return
    end subroutine

    subroutine vmdqtlowvalue(sim,win,iimdqt,Na)
        implicit none

        type(simulation_t), pointer :: sim
        integer j,win,iimdqt
        integer Na
        _REAL_ xx(Na),yy(Na),zz(Na)

        if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
            do j=1,sim%naesmd%natom
                xx(j)=sim%naesmd%rxold(j) &
                    +sim%naesmd%vxold(j)*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%axold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2

                yy(j)=sim%naesmd%ryold(j) &
                    +sim%naesmd%vyold(j)*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%ayold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2

                zz(j)=sim%naesmd%rzold(j) &
                    +sim%naesmd%vzold(j)*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%azold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2
            end do
        end if

        if(sim%naesmd%ensemble.eq.'langev') then
            do j=1,sim%naesmd%natom
                xx(j)=sim%naesmd%rxold(j) &
                    +sim%naesmd%vxold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%axold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2 &
                    +sim%naesmd%prand(1,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win)

                yy(j)=sim%naesmd%ryold(j) &
                    +sim%naesmd%vyold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%ayold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2 &
                    +sim%naesmd%prand(2,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win)

                zz(j)=sim%naesmd%rzold(j) &
                    +sim%naesmd%vzold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win) &
                    +sim%naesmd%azold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt+win))**2 &
                    +sim%naesmd%prand(3,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt+win)
            end do
        end if
         


        call do_sqm_davidson_update(sim,0,cmdqt=sim%naesmd%cmdqtmiddle, &
            vmdqt=sim%naesmd%vmdqtmiddle,vgs=sim%naesmd%vgs,rx=xx,ry=yy,rz=zz)

        if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
            do j=1,sim%naesmd%natom
                xx(j)=sim%naesmd%rxold(j) &
                    +sim%naesmd%vxold(j)*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%axold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2

                yy(j)=sim%naesmd%ryold(j) &
                    +sim%naesmd%vyold(j)*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%ayold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2

                zz(j)=sim%naesmd%rzold(j) &
                    +sim%naesmd%vzold(j)*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%azold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2
            end do
        end if

        if(sim%naesmd%ensemble.eq.'langev') then
            do j=1,sim%naesmd%natom
                xx(j)=sim%naesmd%rxold(j) &
                    +sim%naesmd%vxold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%axold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2 &
                    +sim%naesmd%prand(1,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win)

                yy(j)=sim%naesmd%ryold(j) &
                    +sim%naesmd%vyold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%ayold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2 &
                    +sim%naesmd%prand(2,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win)

                zz(j)=sim%naesmd%rzold(j) &
                    +sim%naesmd%vzold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win) &
                    +sim%naesmd%azold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                    *(sim%naesmd%dtquantum*dfloat(iimdqt-win))**2 &
                    +sim%naesmd%prand(3,j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt-win)
            end do
        end if

        call do_sqm_davidson_update(sim,0,sim%naesmd%cmdqtmiddleold,rx=xx,ry=yy,rz=zz)

        return
    end subroutine

    subroutine do_trivial_wrap(sim)
        type(simulation_t), pointer :: sim
        integer win,i,j
        character(*), parameter :: fmt1 = "(F10.4,1X,F10.4,1X,I3,1X,I3,1X,I3,4(1X,F18.10))"

        if(sim%dotrivial.eq.1.and.conthop.eq.0) then
            call checkcrossing(sim,sim%naesmd%cross,sim%lprint)
            else
               write(6,*)'WARNING:TRIVIAL UNAVOIDED CROSSING DETECTION IS OFF'
               do i=1,sim%naesmd%npot
                  sim%naesmd%cross(i)=0
               enddo
            endif

            sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
            sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)
            !*******************************************
            ! loop to detect the crossing point
            !******************************************
            sim%naesmd%crosstot=0
            do i=1,sim%excN
                if(sim%naesmd%cross(i).eq.1) sim%naesmd%crosstot=1
            end do
            if(sim%naesmd%crosstot.eq.1) then
                if (sim%lprint.gt.1) then
                    write(6,*)'there is crossing'
                    write(6,*)'cross(:)=',sim%naesmd%cross(1:sim%excN)
                endif
                if(sim%qmmm%state_of_interest.gt.0) then
                sim%naesmd%cadiabhop=sim%naesmd%cadiabnew(sim%qmmm%state_of_interest, &
                                                          sim%naesmd%iorden(sim%qmmm%state_of_interest))
                endif
                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal*sim%naesmd%nstepcross
                sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)
                do j=1,sim%excN
                    sim%naesmd%lowvalue(j)=1000.0d0
                end do
                do iimdqt=1,sim%naesmd%nquantumstep
                    sim%naesmd%tfemtoquantum=sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                        +iimdqt*sim%naesmd%dtquantum*convtf
                    call vmdqtmiddlecalc(sim,iimdqt,sim%Na)
                    call checkcrossingmiddle(sim,sim%naesmd%cross)
                    if(sim%lprint.ge.2) then
                        do j=1,sim%excN
                            if(sim%naesmd%cross(j).eq.1) then
                                write(sim%outfile_15,fmt1) sim%naesmd%tfemto, &
                                    sim%naesmd%tfemtoquantum, &
                                    sim%naesmd%cross(j),j,sim%naesmd%iordenhop(j), &
                                    sim%naesmd%scpr(j,sim%naesmd%iordenhop(j)), &
                                    sim%naesmd%scprreal(j,sim%naesmd%iordenhop(j)), &
                                    sim%naesmd%vmdqtmiddle(j)- &
                                    sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j))
                                call flush(sim%outfile_15)
                            end if
                        end do
                    end if
                end do
                do win=1,3
                    do j=1,sim%excN
                        if(sim%naesmd%cross(j).eq.1) then
                            if(j.lt.sim%naesmd%iordenhop(j).or.j.eq.sim%naesmd%ihop) then
                                    call vmdqtlowvalue(sim,win,sim%naesmd%lowvaluestep(j),sim%Na)
                                    call checkcrossingmiddle(sim,sim%naesmd%cross)
                                    if(sim%lprint.ge.2) then
                                        write(sim%outfile_15,fmt1) sim%naesmd%tfemto, &
                                             sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf + &
                                             sim%naesmd%lowvaluestep(j)*sim%naesmd%dtquantum*convtf, &
                                             sim%naesmd%cross(j),j,sim%naesmd%iordenhop(j),&
                                             sim%naesmd%scpr(j,sim%naesmd%iordenhop(j)), &
                                             sim%naesmd%scprreal(j,sim%naesmd%iordenhop(j)),&
                                             sim%naesmd%lowvalue(j)
                                        call flush(sim%outfile_15)
                                    end if
                            end if
                        end if
                    end do
                end do
                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)        
            end if
            !
            !  remove the couplings if cross=2
            !
            if(sim%qmmm%state_of_interest.gt.0) then
            do i=1,sim%excN
                if(i.eq.sim%qmmm%state_of_interest) then

                    if(sim%naesmd%cross(i).eq.2) then
                        sim%naesmd%nquantumstep=sim%naesmd%nquantumreal

                        !if(sim%naesmd%conthop.gt.0) then
                        !    if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop=0
                        !end if

                        !if(sim%naesmd%conthop2.gt.0) then
                        !    if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop2=0
                        !end if

                        if(sim%naesmd%conthop.eq.0) then
                            sim%naesmd%cadiabold(i,sim%naesmd%iorden(i))=0.0d0
                            sim%naesmd%cadiabold(sim%naesmd%iorden(i),i)=0.0d0
                            sim%naesmd%cadiabnew(i,sim%naesmd%iorden(i))=0.0d0
                            sim%naesmd%cadiabnew(sim%naesmd%iorden(i),i)=0.0d0
                        end if
                    end if
                else
                    if(i.ne.sim%naesmd%iorden(sim%qmmm%state_of_interest)) then
                        if(i.lt.sim%naesmd%iorden(i)) then
                            if(sim%naesmd%cross(i).eq.2) then
                                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                                sim%naesmd%cadiabold(i,sim%naesmd%iorden(i))=0.0d0
                                sim%naesmd%cadiabold(sim%naesmd%iorden(i),i)=0.0d0
                                sim%naesmd%cadiabnew(i,sim%naesmd%iorden(i))=0.0d0
                                sim%naesmd%cadiabnew(sim%naesmd%iorden(i),i)=0.0d0
                            end if
                        end if
                    end if
                end if
            end do
            else
            do i=1,sim%excN
                if(i.lt.sim%naesmd%iorden(i)) then
                    if(sim%naesmd%cross(i).eq.2) then
                        sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                        sim%naesmd%cadiabold(i,sim%naesmd%iorden(i))=0.0d0
                        sim%naesmd%cadiabold(sim%naesmd%iorden(i),i)=0.0d0
                        sim%naesmd%cadiabnew(i,sim%naesmd%iorden(i))=0.0d0
                        sim%naesmd%cadiabnew(sim%naesmd%iorden(i),i)=0.0d0
                    end if
                end if
            end do
            endif
    end subroutine

    subroutine quantum_propagation(sim,imdqt,nuclear,Nsim)
        use communism, only: MCE
        implicit none

        type(simulation_t),pointer::sim
        type(MCE) :: nuclear
        integer, intent(in) :: Nsim
        integer,intent(in) :: imdqt
        external :: fcn ! function call for interpolation (has to be passed)
        character(*), parameter :: fmt1 = "(10000(1X,F18.10))"
        integer :: iimdqt, j, k
        _REAL_  :: tini,tfin
        _REAL_, allocatable, dimension(:) :: yg_new_fin
        _REAL_ :: pha_fin

        allocate(yg_new_fin(3*sim%excN))

            do iimdqt=1,sim%naesmd%nquantumstep
                write(6,*)'Begin quantum step #',iimdqt
                sim%naesmd%tfemtoquantum=sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                    +iimdqt*sim%naesmd%dtquantum*convtf
                ! Definition of initial and final time for the quantum propagator

                tini=(sim%naesmd%tfemtoquantum-sim%naesmd%dtquantum*convtf)/convtf
                tfin=tini+sim%naesmd%dtquantum
                sim%naesmd%tini0=tini!used for interpolation
                !--------------------------------------------------------------------
                ! Calculation of the coordinates at the middle of the step.
                ! This intermediate structure is obtained using Newton equation
                ! with constant velocities and accelerations.
                ! and calculation of the energies(=sim%naesmd%vmdqtmiddle) and nacT(=sim%naesmd%cadiabmiddle) values
                ! at intermediate times of the classical step
                !--------------------------------------------------------------------
                call cadiabmiddlecalc(sim,iimdqt,sim%Na,sim%naesmd%cross)

                !--------------------------------------------------------------------
                !
                !  Calculation of parameters to fit the values of sim%naesmd%cadiab and
                !  sim%naesmd%vmdqt during propagation
                !
                !--------------------------------------------------------------------
                call fitcoef(sim)
                !--------------------------------------------------------------------
                ! Runge-Kutta-Verner propagator
                !--------------------------------------------------------------------
                do k = 1, 10 !loop to improve RK
                    tfin = tini + sim%naesmd%dtquantum/10
                    call integrate_1step_a(tini,tfin,sim%naesmd%yg_new,yg_new_fin,sim)
                    sim%naesmd%yg_new = yg_new_fin
                    call integrate_1step_gamma(tini,tfin,sim%naesmd%pha,pha_fin,sim)
                    sim%naesmd%pha = pha_fin
                    tini = tfin
                enddo
                ! Check the norm
                !******************************************************
                ! values for hop probability
                do k=1,sim%excN
                    do j=1,sim%excN
                        sim%naesmd%vnqcorrhop(k,j)=2.0d0*sim%naesmd%yg_new(k)*(-1.0d0)*(sim%naesmd%yg_new(j) &
                            *dcos(sim%naesmd%yg_new(j+2*sim%excN)-sim%naesmd%yg_new(k+2*sim%excN))-sim%naesmd%yg_new(j+sim%excN) &
                            *dsin(sim%naesmd%yg_new(j+2*sim%excN)-sim%naesmd%yg_new(k+2*sim%excN)))*sim%naesmd%cadiabmiddle(k,j) &
                            +2.0d0*sim%naesmd%yg_new(k+sim%excN)*(-1.0d0)*(sim%naesmd%yg_new(j+sim%excN) &
                            *dcos(sim%naesmd%yg_new(j+2*sim%excN)-sim%naesmd%yg_new(k+2*sim%excN))+sim%naesmd%yg_new(j) &
                            *dsin(sim%naesmd%yg_new(j+2*sim%excN)-sim%naesmd%yg_new(k+2*sim%excN)))*sim%naesmd%cadiabmiddle(k,j)
                    end do
                end do

                do k=1,sim%excN
                    do j=1,sim%excN
                        sim%naesmd%vnqcorrhoptot(k,j)=sim%naesmd%vnqcorrhoptot(k,j) &
                            +sim%naesmd%vnqcorrhop(k,j)*sim%naesmd%dtquantum
                    end do
                end do
                write(6,*)'End quantum step #',iimdqt
            end do
    end subroutine

    subroutine check_ntotcoher(sim)
        type(simulation_t),pointer::sim
        integer :: j
        _REAL_ :: ntotcoher
            if(sim%naesmd%icontw.eq.sim%naesmd%nstepw) then
                if(sim%naesmd%state.eq.'exct') then
                    if(sim%lprint.ge.1) then
                        ntotcoher=0.0d0

                        do j=1,sim%excN
                            ntotcoher=ntotcoher+sim%naesmd%yg_new(j)**2+sim%naesmd%yg_new(j+sim%excN)**2
                        end do
!BTN: removed file coeff-n-before.out. grep this line to undo
                        call flush(105)
                    end if
                end if
            end if
        
    end subroutine

!
!********************************************************************      
!
!  Subroutine for integrating the TDS for electronic amplitudes "a"
!  by means of a 4th order Runge-Kutta propagator
!
!********************************************************************      
!   
subroutine integrate_1step_a(tini,tfin,yg_new,yg_new_fin,sim)
    implicit none

    _REAL_, intent(in) :: tini
    _REAL_, intent(in) :: tfin
    _REAL_, intent(in), dimension(:) :: yg_new
    _REAL_, intent(out), dimension(:) :: yg_new_fin
    type(simulation_t), pointer :: sim

    _REAL_ dt
    _REAL_, dimension(4) :: pfac1
    _REAL_, dimension(4) :: pfac2
    integer irk4, j, k
    _REAL_, allocatable, dimension(:) :: krk4, krk4_aux
    _REAL_ tini_aux

    allocate(krk4(3*sim%excN),krk4_aux(3*sim%excN))

    dt = tfin - tini
    pfac1(1) = dt*0.d0
    pfac1(2) = dt*5.d-1
    pfac1(3) = dt*5.d-1
    pfac1(4) = dt*1.d0
    pfac2(1) = dt/6.d0
    pfac2(2) = dt/3.d0
    pfac2(3) = dt/3.d0
    pfac2(4) = dt/6.d0

    yg_new_fin = yg_new
    krk4 = 0.0d0
    tini_aux = tini
    do irk4 = 1, 4
        tini_aux = tini + pfac1(irk4)
        call interpolate(size(yg_new,1),tini_aux,sim%naesmd)
        krk4_aux = krk4
        do k=1,sim%excN
            krk4(k) = 0.d0
            krk4(k+sim%excN) = 0.d0
            krk4(k+2*sim%excN)=-sim%naesmd%vmdqt(k)
            do j=1,sim%excN
                krk4(k)=krk4(k)-((yg_new(j)+krk4_aux(j)*pfac1(irk4)) &
                    *dcos(yg_new(j+2*sim%excN)+krk4_aux(j+2*sim%excN)*pfac1(irk4)-yg_new(k+2*sim%excN)-krk4_aux(k+2*sim%excN)*pfac1(irk4))-(yg_new(j+sim%excN)+krk4_aux(j+sim%excN)*pfac1(irk4)) &
                    *dsin(yg_new(j+2*sim%excN)+krk4_aux(j+2*sim%excN)*pfac1(irk4)-yg_new(k+2*sim%excN)-krk4_aux(k+2*sim%excN)*pfac1(irk4)))*sim%naesmd%cadiab(k,j)

                krk4(k+sim%excN) = krk4(k+sim%excN)-((yg_new(j+sim%excN)+krk4_aux(j+sim%excN)*pfac1(irk4)) &
                    *dcos(yg_new(j+2*sim%excN)+krk4_aux(j+2*sim%excN)*pfac1(irk4)-yg_new(k+2*sim%excN)-krk4_aux(k+2*sim%excN)*pfac1(irk4))+(yg_new(j)+krk4_aux(j)*pfac1(irk4)) &
                    *dsin(yg_new(j+2*sim%excN)+krk4_aux(j+2*sim%excN)*pfac1(irk4)-yg_new(k+2*sim%excN)-krk4_aux(k+2*sim%excN)*pfac1(irk4)))*sim%naesmd%cadiab(k,j)
            end do
        end do
        yg_new_fin(:) = yg_new_fin(:) + pfac2(irk4) * krk4(:)
    enddo !irk4

end subroutine integrate_1step_a

!
!********************************************************************
!
!  Subroutine to propagate the ODE corresponding to the nuclear
!  phases (gamma) by means of a 4th order RK algorithm
!
!********************************************************************
!
subroutine integrate_1step_gamma(tini,tfin,pha,pha_fin,sim)
    implicit none
    _REAL_, intent(in) :: tini, tfin
    _REAL_, intent(in) :: pha
    _REAL_, intent(out) :: pha_fin
    type(simulation_t) :: sim

    _REAL_ dt
    _REAL_, dimension(4) :: pfac1
    _REAL_, dimension(4) :: pfac2
    integer irk4, i
    _REAL_ :: krk4, krk4_aux
    _REAL_ :: tini_aux
    _REAL_ :: vxt(sim%naesmd%natom)
    _REAL_ :: vyt(sim%naesmd%natom)
    _REAL_ :: vzt(sim%naesmd%natom)

    dt = tfin - tini
    pfac1(1) = dt*0.d0
    pfac1(2) = dt*5.d-1
    pfac1(3) = dt*5.d-1
    pfac1(4) = dt*1.d0
    pfac2(1) = dt/6.d0
    pfac2(2) = dt/3.d0
    pfac2(3) = dt/3.d0
    pfac2(4) = dt/6.d0

    pha_fin = pha
    krk4 = 0.0d0
    tini_aux = tini
    do irk4 = 1, 4
        !krk4 = fp(t0 + pfac1(irk4), f0 + krk4 * pfac1(irk4))
        tini_aux = tini + pfac1(irk4)
        !interpalation for r and p
        call interpolate_veloc(tini,tini_aux,tfin,sim%naesmd%vx,vxt,sim%naesmd%vxold)
        call interpolate_veloc(tini,tini_aux,tfin,sim%naesmd%vy,vyt,sim%naesmd%vyold)
        call interpolate_veloc(tini,tini_aux,tfin,sim%naesmd%vz,vzt,sim%naesmd%vzold)
        krk4_aux = krk4
        krk4 = 0.0d0
        do i = 1, sim%naesmd%natom
            krk4 = krk4 + sim%naesmd%massmdqt(i)*vxt(i)*vxt(i)/2.0d0
            krk4 = krk4 + sim%naesmd%massmdqt(i)*vyt(i)*vyt(i)/2.0d0
            krk4 = krk4 + sim%naesmd%massmdqt(i)*vzt(i)*vzt(i)/2.0d0
        enddo
        pha_fin = pha_fin + pfac2(irk4) * krk4
    enddo !irk4

end subroutine integrate_1step_gamma

!
!**************************************************************
!
!  This subroutine interpolates linearly the velocity 
!
!**************************************************************
!

subroutine interpolate_veloc(tini,t,tfin,v0,v,vf)
    implicit none
    _REAL_, intent(in) :: tini,t,tfin
    _REAL_, dimension(:), intent(in) :: v0, vf
    _REAL_, dimension(:), intent(out) :: v
    
    v = (vf*(t-tini)-v0*(t-tfin))/(tfin-tini)

end subroutine interpolate_veloc


!********************************************************************
!
!  Subroutine to interpolate the values of the terms <phi_k | d phi_j/ d t >
!  and the adiabatic energies
!
!********************************************************************
!
subroutine interpolate(n,x,naesmd_struct)
    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct
    integer n,k,j
    _REAL_ x

    do k=1,n/3
        do j=1,n/3
            naesmd_struct%cadiab(k,j)=naesmd_struct%cadiabmiddleold(k,j) &
                +naesmd_struct%bcoeffcadiab(k,j)*(x-naesmd_struct%tini0)
        end do
        naesmd_struct%cadiab(k,k)=0.0d0
    end do

    do k=1,n/3
        naesmd_struct%vmdqt(k)=naesmd_struct%vmdqtmiddleold(k)+naesmd_struct%bcoeffvmdqt(k)*(x-naesmd_struct%tini0)
    end do

    return
end subroutine

!
!*********************************************************************
!
!  Subroutine for integrating the ODEs for electronic overlpas "sE"
!
!*********************************************************************
!

subroutine propagate_sE(sims,nuclear,Nsim)
    use communism
    use AIMC_type_module
    use naesmd_constants
    implicit none
    type(sim_pointer_array), dimension(:), intent(in) :: sims
    type(MCE) :: nuclear
    integer, intent(in) :: Nsim

    integer i,j,k,n,m
    _REAL_, allocatable, dimension(:,:) :: sE_fin
    _REAL_ :: tini,tfin

    allocate(sE_fin(sims(1)%sim%excN,sims(1)%sim%excN))

    do n=1,Nsim
    do m=1,Nsim
        tini=(sims(1)%sim%naesmd%tfemtoquantum-sims(1)%sim%naesmd%dtquantum*convtf)/convtf
        do k=1,10
            tfin = tini + sims(1)%sim%naesmd%dtquantum/10
            call integrate_1step_sE(tini,tfin,nuclear%sE(n,m,:,:),sE_fin,sims(n)%sim,sims(m)%sim)
            nuclear%sE(n,m,:,:) = sE_fin
            tini = tfin
        enddo
    enddo
    enddo

end subroutine propagate_sE

!
!********************************************************************
!
!   Subroutine for integrating one step of the RK4 method for the
!   electronic overlaps "sE"
!
!********************************************************************
!

subroutine integrate_1step_sE(tini,tfin,sE,sE_fin,sim1,sim2)
    implicit none
    _REAL_, intent(in) :: tini
    _REAL_, intent(in) :: tfin
    _REAL_, intent(in), dimension(:,:) :: sE !Electronic overlaps between two trajectroies
    _REAL_, intent(out), dimension(:,:) :: sE_fin
    type(simulation_t), pointer :: sim1
    type(simulation_t), pointer :: sim2

    _REAL_ dt
    _REAL_, dimension(4) :: pfac1
    _REAL_, dimension(4) :: pfac2
    integer irk4, i,j,m,n,k
    _REAL_, allocatable, dimension(:,:) :: krk4, krk4_aux
    _REAL_ tini_aux

    allocate(krk4(sim1%excN,sim2%excN), krk4_aux(sim1%excN,sim2%excN))

    dt = tfin - tini
    pfac1(1) = dt*0.d0
    pfac1(2) = dt*5.d-1
    pfac1(3) = dt*5.d-1
    pfac1(4) = dt*1.d0
    pfac2(1) = dt/6.d0
    pfac2(2) = dt/3.d0
    pfac2(3) = dt/3.d0
    pfac2(4) = dt/6.d0

    sE_fin = sE
    krk4 = 0.0d0
    tini_aux = tini
    do irk4 = 1, 4
!        krk4 = fp(t0 + pfac1(irk4), f0 + krk4 * pfac1(irk4))
        tini_aux = tini + pfac1(irk4)
        call interpolate(3*sim1%excN,tini_aux,sim1%naesmd)
        call interpolate(3*sim2%excN,tini_aux,sim2%naesmd)
        krk4_aux = krk4
        do i=1,sim1%excN
            do j=1,sim2%excN
                krk4(i,j) = 0.0d0
                do k=1,sim1%excN
                    krk4(i,j)=krk4(i,j)+sim1%naesmd%cadiab(k,i)*(sE(k,j)+krk4_aux(k,j)*pfac1(irk4))+&
                                        sim2%naesmd%cadiab(k,j)*(sE(i,k)+krk4_aux(i,k)*pfac1(irk4))
                enddo
            enddo
        enddo
!        ff = ff + pfac2(irk4) * krk4
        sE_fin = sE_fin + pfac2(irk4) * krk4
    enddo!irk4

end subroutine integrate_1step_sE

end module
