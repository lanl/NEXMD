#include "dprec.fh"
#include "assert.fh"

module quantum_prop_add
    use naesmd_constants
    use communism
contains

    subroutine checkcrossing(sim,cross,lprint)
        use naesmd_module
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
        !  scpr is the overlap between new and old transition densities
        !  if there is no crossing, the matrix is supposed to be close
        !  to identity matrix
        !
        !--------------------------------------------------------------------
        if (sim%excN>260) then
            write(6,*)'Error:quantum propagator does not yet handle more than 260 states'
        endif

        scpr(1:sim%excN,1:sim%excN)=0.d0

        ! kav: FIXME
        ! Here and below, in the old code the size Ncis is used as a size of
        ! cmdqt vectors. This is not entirely correct, since the size is
        ! Nrpa=2*Ncis (check it!!!). However, it can be very accurate becuase
        ! the contribution of negative energies to a positive excitation
        ! energy can be very small (would be great to check it)
        do i=1,sim%excN
            do j=1,sim%excN
                do k=1,sim%dav%Ncis
                    scpr(i,j)=scpr(i,j)+cmdqtold(k,i)*cmdqtnew(k,j)
                end do
            end do
        end do

        do i=1,sim%excN
            do j=1,sim%excN
                ascpr(i,j)=int(scpr(i,j)**2*1.d5)
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
        !write(6,*)'iorden:',iorden
        
        !Translate to old APC min-cost algorithm routine in F77
        !Requires static matrix sizes
        iordenapc(1:sim%excN)=iorden
        call apc(sim%excN,ascpr,iorden,z)
        iorden=iordenapc(1:sim%excN)

        do i=1,sim%excN
            !write(6,*)'i,iorden(i),ihop,scpr(i,iorden(i))',i,iorden(i),ihop,scpr(i,iorden(i))
            if(iorden(i).ne.i) then
                if(i.lt.iorden(i).or.i.eq.ihop) then
                    if(dabs(scpr(i,iorden(i))).lt.0.9d0) then
                        cross(i)=1
                    else
                        cross(i)=2
                    end if
                    iordenhop(i)=iorden(i)
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
                scprreal(i,j)=scpr(i,j)
            end do
        end do

        !***********************************


        return
    end subroutine

    subroutine checkcrossingmiddle(sim,cross)
        use naesmd_module
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
                scpr(i,j)=0.0d0
                do k=1,sim%dav%Ncis
                    scpr(i,j)=scpr(i,j)+cmdqtmiddleold(k,i)*cmdqtmiddle(k,j)
                end do
            end do
        end do
        do i=1,sim%excN
            if(i.lt.iorden(i).or.i.eq.ihop) then
                if(i.ne.iorden(ihop)) then
                    if(dabs(scpr(i,iorden(i))).ge.0.9d0) then
                        cross(i)=2
                    end if
                end if
            end if
        end do
        return
    end subroutine

    subroutine vmdqtmiddlecalc(sim, iimdqt,Na)
        use naesmd_module
        use md_module
        implicit none

        type(simulation_t), pointer :: sim
        integer i,j,iimdqt
        integer Na
        _REAL_ xx(Na),yy(Na),zz(Na)

        if(iimdqt.eq.1) then
            do i=1,sim%excN
                vmdqtmiddleold(i)=vmdqtold(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    cmdqtmiddleold(i,j)=cmdqtold(i,j)
                end do
            end do

        else
            do i=1,sim%excN
                vmdqtmiddleold(i)=vmdqtmiddle(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    cmdqtmiddleold(i,j)=cmdqtmiddle(i,j)
                end do
            end do
        end if

        if(iimdqt.eq.nquantumstep) then

            do i=1,sim%excN
                vmdqtmiddle(i)=vmdqtnew(i)
            end do

            do i=1,sim%dav%Ncis
                do j=1,sim%excN
                    cmdqtmiddle(i,j)=cmdqtnew(i,j)
                end do
            end do

        else

            if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
                do j=1,natom
                    xx(j)=rxold(j) &
                        +vxold(j)*dtquantum*dfloat(iimdqt) &
                        +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

                    yy(j)=ryold(j) &
                        +vyold(j)*dtquantum*dfloat(iimdqt) &
                        +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

                    zz(j)=rzold(j) &
                        +vzold(j)*dtquantum*dfloat(iimdqt) &
                        +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2
                end do
            end if

            if(ensemble.eq.'langev') then
                do j=1,natom

                    xx(j)=rxold(j) &
                        +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
                        +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
                        *(dtquantum*dfloat(iimdqt))**2 &
                        +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt)

                    yy(j)=ryold(j) &
                        +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
                        +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
                        *(dtquantum*dfloat(iimdqt))**2 &
                        +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt)

                    zz(j)=rzold(j) &
                        +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
                        +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
                        *(dtquantum*dfloat(iimdqt))**2 &
                        +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt)
                end do
            end if


            call do_sqm_davidson_update(sim,cmdqt=cmdqtmiddle, &
                vmdqt=vmdqtmiddle,vgs=vgs,rx=xx,ry=yy,rz=zz)
        end if
        write(6,*)'iordenhop:',iordenhop
        write(6,*)'sizevmdqtmiddle',shape(vmdqtmiddle)
! option 1: check trivial unavoided crossingfor all the states
! Note: Currently broken due to iordenhop initialized as 0
!        do j=1,sim%excN
!            if(lowvalue(j).gt.dabs(vmdqtmiddle(j)- &
!                vmdqtmiddle(iordenhop(j)))) then
!
!                lowvalue(j)=dabs(vmdqtmiddle(j)- &
!                    vmdqtmiddle(iordenhop(j)))
!                lowvaluestep(j)=iimdqt
!            end if
!        end do
! option 2: check trivial unavoided crossingfor for the current state
        if(lowvalue(ihop).gt.dabs(vmdqtmiddle(ihop)- &
                vmdqtmiddle(iordenhop(ihop)))) then
                lowvalue(ihop)=dabs(vmdqtmiddle(ihop)- &
                        vmdqtmiddle(iordenhop(ihop)))
                lowvaluestep(ihop)=iimdqt
        end if
        return
    end subroutine

    subroutine vmdqtlowvalue(sim,win,iimdqt,Na)
        use naesmd_module
        use md_module
        implicit none

        type(simulation_t), pointer :: sim
        integer j,win,iimdqt
        integer Na
        _REAL_ xx(Na),yy(Na),zz(Na)

        if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
            do j=1,natom
                xx(j)=rxold(j) &
                    +vxold(j)*dtquantum*dfloat(iimdqt+win) &
                    +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2

                yy(j)=ryold(j) &
                    +vyold(j)*dtquantum*dfloat(iimdqt+win) &
                    +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2

                zz(j)=rzold(j) &
                    +vzold(j)*dtquantum*dfloat(iimdqt+win) &
                    +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2
            end do
        end if

        if(ensemble.eq.'langev') then
            do j=1,natom
                xx(j)=rxold(j) &
                    +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
                    +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt+win))**2 &
                    +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)

                yy(j)=ryold(j) &
                    +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
                    +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt+win))**2 &
                    +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)

                zz(j)=rzold(j) &
                    +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
                    +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt+win))**2 &
                    +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)
            end do
        end if
         


        call do_sqm_davidson_update(sim,cmdqt=cmdqtmiddle, &
            vmdqt=vmdqtmiddle,vgs=vgs,rx=xx,ry=yy,rz=zz)

        if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
            do j=1,natom
                xx(j)=rxold(j) &
                    +vxold(j)*dtquantum*dfloat(iimdqt-win) &
                    +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2

                yy(j)=ryold(j) &
                    +vyold(j)*dtquantum*dfloat(iimdqt-win) &
                    +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2

                zz(j)=rzold(j) &
                    +vzold(j)*dtquantum*dfloat(iimdqt-win) &
                    +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2
            end do
        end if

        if(ensemble.eq.'langev') then
            do j=1,natom
                xx(j)=rxold(j) &
                    +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
                    +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt-win))**2 &
                    +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)

                yy(j)=ryold(j) &
                    +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
                    +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt-win))**2 &
                    +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)

                zz(j)=rzold(j) &
                    +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
                    +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
                    *(dtquantum*dfloat(iimdqt-win))**2 &
                    +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)
            end do
        end if

        call do_sqm_davidson_update(sim,cmdqtmiddleold,rx=xx,ry=yy,rz=zz)

        return
    end subroutine
end module
