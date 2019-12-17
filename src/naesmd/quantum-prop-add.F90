#include "dprec.fh"
#include "assert.fh"

module quantum_prop_add
    use naesmd_constants
    use communism
    use cadiab_module
    use additional_subroutines
    use rksuite_90

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
        integer tmp(sim%excN)

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
            if(sim%naesmd%iorden(i).ne.i) then
                if(i.lt.sim%naesmd%iorden(i).or.i.eq.sim%naesmd%ihop) then
                    if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).lt.0.9d0) then
                        cross(i)=1
                    else
                        cross(i)=2
                        tmp = sim%naesmd%dbtorden
                        if(i.lt.sim%naesmd%iorden(i)) then
                            sim%naesmd%dbtorden(i) = tmp(sim%naesmd%iorden(i))
                            sim%naesmd%dbtorden(i+1) = tmp(sim%naesmd%iorden(i+1))
                        endif
                    end if
                    sim%naesmd%iordenhop(i)=sim%naesmd%iorden(i)
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
        integer tmp(sim%excN)

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
            if(i.lt.sim%naesmd%iorden(i).or.i.eq.sim%naesmd%ihop) then
                if(i.ne.sim%naesmd%iorden(sim%naesmd%ihop)) then
                    if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).ge.0.9d0) then
                        cross(i)=2
                        tmp = sim%naesmd%dbtorden
                        if(i.lt.sim%naesmd%iorden(i)) then
                            sim%naesmd%dbtorden(i) = tmp(sim%naesmd%iorden(i))
                            sim%naesmd%dbtorden(i+1) = tmp(sim%naesmd%iorden(i+1))
                        endif
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
        integer temp(1:sim%excN)

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


            call do_sqm_davidson_update(sim,cmdqt=sim%naesmd%cmdqtmiddle, &
                vmdqt=sim%naesmd%vmdqtmiddle,vgs=sim%naesmd%vgs,rx=xx,ry=yy,rz=zz)
        end if
        write(6,*)'sim%naesmd%iordenhop:',sim%naesmd%iordenhop
        write(6,*)'sizevmdqtmiddle',shape(sim%naesmd%vmdqtmiddle)
! option 2: check trivial unavoided crossingfor for the current sim%naesmd%state

        temp = sim%naesmd%iordenhop
        do j=1,sim%excN
            if(sim%naesmd%iordenhop(j) .eq. 0) then
                sim%naesmd%iordenhop(j) = j
            endif
        enddo
        do j=1,sim%excN
            if(sim%naesmd%lowvalue(j).gt.dabs(sim%naesmd%vmdqtmiddle(j)- &
                sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))) then

                sim%naesmd%lowvalue(j)=dabs(sim%naesmd%vmdqtmiddle(j)- &
                    sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))
                sim%naesmd%lowvaluestep(j)=iimdqt
            end if
        end do
        sim%naesmd%iordenhop = temp
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
         


        call do_sqm_davidson_update(sim,cmdqt=sim%naesmd%cmdqtmiddle, &
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

        call do_sqm_davidson_update(sim,sim%naesmd%cmdqtmiddleold,rx=xx,ry=yy,rz=zz)

        return
    end subroutine

    subroutine do_trivial_wrap(sim)
        type(simulation_t), pointer :: sim
        integer win,i,j
        character(*), parameter :: fmt1 = "(F7.4,1X,F7.4,1X,I3,1X,I3,1X,I3,4(1X,F18.10))"

        if(sim%dotrivial.eq.1) then
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
                sim%naesmd%cadiabhop=sim%naesmd%cadiabnew(sim%qmmm%state_of_interest, &
                                                          sim%naesmd%iorden(sim%qmmm%state_of_interest))
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
                                if(j.ne.sim%naesmd%iordenhop(j)) then
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
                        end if
                    end do
                end do
                sim%naesmd%nquantumstep=sim%naesmd%nquantumreal
                sim%naesmd%dtquantum=sim%naesmd%dtmdqt/dfloat(sim%naesmd%nquantumstep)        
            end if
            !
            !  remove the couplings if cross=2
            !
            do i=1,sim%excN
                if(i.eq.sim%qmmm%state_of_interest) then

                    if(sim%naesmd%cross(i).eq.2) then
                        sim%naesmd%nquantumstep=sim%naesmd%nquantumreal

                        if(sim%naesmd%conthop.gt.0) then
                            if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop=0
                        end if

                        if(sim%naesmd%conthop2.gt.0) then
                            if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop2=0
                        end if

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
    end subroutine

    subroutine quantum_propagation(sim,imdqt)
        implicit none

        type(simulation_t),pointer::sim
        integer,intent(in) :: imdqt
        external :: fcn ! function call for interpolation (has to be passed)
        character(*), parameter :: fmt1 = "(10000(1X,F18.10))"
        integer :: iimdqt, j, k
        _REAL_  :: tini

            do iimdqt=1,sim%naesmd%nquantumstep
                write(6,*)'Begin quantum step #',iimdqt
                sim%naesmd%tfemtoquantum=sim%naesmd%tfemto-sim%naesmd%dtmdqt*convtf &
                    +iimdqt*sim%naesmd%dtquantum*convtf
                ! Definition of initial and final time for the quantum propagator
                if(imdqt.eq.1.and.iimdqt.eq.1) then
                    tini=0.0d0
                else
                    tini=sim%rk_comm%tend
                end if

                sim%rk_comm%tend=sim%naesmd%tfemtoquantum/convtf
                sim%naesmd%tini0=tini
                !--------------------------------------------------------------------
                ! Calculation of the coordinates at the middle of the step.
                ! This intermediate structure is obtained using Newton equation
                ! with constant velocities and accelerations.
                ! and calculation of the energies(=sim%naesmd%vmdqtmiddle) and nacT(=sim%naesmd%cadiabmiddle) values
                ! at intermediate times of the classical step
                !--------------------------------------------------------------------
                call cadiabmiddlecalc(sim,iimdqt,sim%Na,sim%naesmd%cross)

                if(sim%lprint.ge.3) then
                    if(iimdqt.ne.sim%naesmd%nquantumstep) then
                        write(sim%outfile_4,fmt1) sim%naesmd%tfemtoquantum, &
                            sim%naesmd%vgs*feVmdqt,(sim%naesmd%vmdqtmiddle(j)*feVmdqt,j=1,sim%excN)

                        write(sim%outfile_7,fmt1) sim%naesmd%tfemtoquantum, &
                            ((sim%naesmd%cadiabmiddle(j,k),k=1,sim%excN),j=1,sim%excN)
                        call flush(96)
                        call flush(93)
                    end if
                end if
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
	 	sim%rk_comm%tend=min(sim%rk_comm%tend,sim%rk_comm%tmax)	
		call range_integrate(sim%rk_comm,fcn,sim%rk_comm%tend,sim%rk_comm%tgot,&
                                     sim%naesmd%yg,sim%naesmd%ygprime,&
                                     sim%naesmd,sim%rk_comm%rk_flag)
		sim%rk_comm%tend=sim%rk_comm%tgot
                ! Check the norm
                call checknorm(sim%rk_comm,sim%excN,sim%rk_comm%tend,sim%rk_comm%tmax,&
                               sim%rk_comm%rk_tol,sim%rk_comm%thresholds,sim%naesmd%yg)
                !******************************************************
                ! values for hop probability
                do k=1,sim%excN
                    do j=1,sim%excN
                        sim%naesmd%vnqcorrhop(k,j)=-1.0d0*sim%naesmd%yg(j) &
                            *dcos(sim%naesmd%yg(j+sim%excN)-sim%naesmd%yg(k+sim%excN))*sim%naesmd%cadiabmiddle(k,j)
                        sim%naesmd%vnqcorrhop(k,j)=sim%naesmd%vnqcorrhop(k,j)*2.0d0*sim%naesmd%yg(k)
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
                            ntotcoher=ntotcoher+sim%naesmd%yg(j)**2
                        end do
!BTN: removed file coeff-n-before.out. grep this line to undo
                        call flush(105)
                    end if
                end if
            end if
        
    end subroutine


end module

