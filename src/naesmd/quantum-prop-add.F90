#include "dprec.fh"
#include "assert.fh"

module quantum_prop_add
    use naesmd_constants
    use communism
contains

    subroutine checkcrossing(sim,cross,lprint)
        use md_module
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
        !write(6,*)'sim%naesmd%iorden:',sim%naesmd%iorden
        
        !Translate to old APC min-cost algorithm routine in F77
        !Requires static matrix sizes
        iordenapc(1:sim%excN)=sim%naesmd%iorden
        call apc(sim%excN,ascpr,sim%naesmd%iorden,z)
        sim%naesmd%iorden=iordenapc(1:sim%excN)

        do i=1,sim%excN
            !write(6,*)'i,sim%naesmd%iorden(i),sim%naesmd%ihop,sim%naesmd%scpr(i,sim%naesmd%iorden(i))',i,sim%naesmd%iorden(i),sim%naesmd%ihop,sim%naesmd%scpr(i,sim%naesmd%iorden(i))
            if(sim%naesmd%iorden(i).ne.i) then
                if(i.lt.sim%naesmd%iorden(i).or.i.eq.sim%naesmd%ihop) then
                    if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).lt.0.9d0) then
                        cross(i)=1
                    else
                        cross(i)=2
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
        use md_module
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
            if(i.lt.sim%naesmd%iorden(i).or.i.eq.sim%naesmd%ihop) then
                if(i.ne.sim%naesmd%iorden(sim%naesmd%ihop)) then
                    if(dabs(sim%naesmd%scpr(i,sim%naesmd%iorden(i))).ge.0.9d0) then
                        cross(i)=2
                    end if
                end if
            end if
        end do
        return
    end subroutine

    subroutine vmdqtmiddlecalc(sim, iimdqt,Na)
        use md_module
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


            call do_sqm_davidson_update(sim,cmdqt=sim%naesmd%cmdqtmiddle, &
                vmdqt=sim%naesmd%vmdqtmiddle,vgs=sim%naesmd%vgs,rx=xx,ry=yy,rz=zz)
        end if
        write(6,*)'sim%naesmd%iordenhop:',sim%naesmd%iordenhop
        write(6,*)'sizevmdqtmiddle',shape(sim%naesmd%vmdqtmiddle)
! option 1: check trivial unavoided crossingfor all the states
! Note: Currently broken due to sim%naesmd%iordenhop initialized as 0
!        do j=1,sim%excN
!            if(sim%naesmd%lowvalue(j).gt.dabs(sim%naesmd%vmdqtmiddle(j)- &
!                sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))) then
!
!                sim%naesmd%lowvalue(j)=dabs(sim%naesmd%vmdqtmiddle(j)- &
!                    sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(j)))
!                sim%naesmd%lowvaluestep(j)=iimdqt
!            end if
!        end do
! option 2: check trivial unavoided crossingfor for the current sim%naesmd%state
        if(sim%naesmd%lowvalue(sim%naesmd%ihop).gt.dabs(sim%naesmd%vmdqtmiddle(sim%naesmd%ihop)- &
                sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(sim%naesmd%ihop)))) then
                sim%naesmd%lowvalue(sim%naesmd%ihop)=dabs(sim%naesmd%vmdqtmiddle(sim%naesmd%ihop)- &
                        sim%naesmd%vmdqtmiddle(sim%naesmd%iordenhop(sim%naesmd%ihop)))
                sim%naesmd%lowvaluestep(sim%naesmd%ihop)=iimdqt
        end if
        return
    end subroutine

    subroutine vmdqtlowvalue(sim,win,iimdqt,Na)
        use md_module
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
end module
