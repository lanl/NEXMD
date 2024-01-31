#include "dprec.fh"
#include "assert.fh"
!*******************************************************
! The analytic NAC for t.
! are calculated inside of cadiaboldcalc,cadiabmiddlecalc, and cadiabnewcalc
! Input should be xxp,yyp,zzp 'plus' xyz at t+dt
! xxm,yym,zzm 'minus' - xyz at t-dt, dt should correspond ~10^-4 10^-5 A shift
! Vectors and frequencies from the above ceo, at xx,yy,zz geometry
! The routine will output sim%naesmd%cadiab_analt - an array of NA couplings for T
! from sim%naesmd%state sim%naesmd%ihop to all other states sim%naesmd%npot
! xxp,yyp, and zzp are xyz at t + sim%naesmd%dtnact
! xxm,yym, and zzm are xyz at t - sim%naesmd%dtnact

module cadiab_module
    use naesmd_constants
    use communism
    use nacT_analytic_module
contains
    subroutine cadiaboldcalc(sim, imdqt)
        use qm2_davidson_module
        implicit none

        type(simulation_t), pointer :: sim
        integer k, j, i, imdqt, l
        type(xstep_t), target::xs
        type(xstep_t), pointer::xstep
        type(realp_t) :: rold(3), vold(3)
        type(realp_t) :: deltaRp(3), deltaRm(3), deltaRpold(3), deltaRmold(3)
        type(simulation_t), pointer::simpoint
        _REAL_ NACR(sim%naesmd%natom*3)

        simpoint => sim

        allocate (xs%Rp(3, sim%Na))
        allocate (xs%Rm(3, sim%Na))
        allocate (xs%R(3, sim%Na))
        rold(1)%p => sim%naesmd%rxold
        rold(2)%p => sim%naesmd%ryold
        rold(3)%p => sim%naesmd%rzold
        vold(1)%p => sim%naesmd%vxold
        vold(2)%p => sim%naesmd%vyold
        vold(3)%p => sim%naesmd%vzold
        deltaRp(1)%p => sim%naesmd%deltaxxpnew
        deltaRp(2)%p => sim%naesmd%deltayypnew
        deltaRp(3)%p => sim%naesmd%deltazzpnew
        deltaRm(1)%p => sim%naesmd%deltaxxmnew
        deltaRm(2)%p => sim%naesmd%deltayymnew
        deltaRm(3)%p => sim%naesmd%deltazzmnew
        deltaRpold(1)%p => sim%naesmd%deltaxxpold
        deltaRpold(2)%p => sim%naesmd%deltayypold
        deltaRpold(3)%p => sim%naesmd%deltazzpold
        deltaRmold(1)%p => sim%naesmd%deltaxxmold
        deltaRmold(2)%p => sim%naesmd%deltayymold
        deltaRmold(3)%p => sim%naesmd%deltazzmold

        xstep => xs

        if (imdqt .eq. 1) then
            call do_sqm_davidson_update(sim, 0, cmdqt=sim%naesmd%cmdqtnew, &
                vmdqt=sim%naesmd%vmdqtnew, vgs=sim%naesmd%vgs, r=rold)!, sim%naesmd%cmdqt=sim%naesmd%cmdqtnew) !JAKB, not necessary bc updated already in
            call new_xstep_dtnact_r3(sim, rold, xstep)
            do k = 1, 3
                deltaRp(k)%p(1:sim%naesmd%natom) = xstep%Rp(k, 1:sim%naesmd%natom) &
                    - xstep%R(k, 1:sim%naesmd%natom)
                deltaRm(k)%p(:sim%naesmd%natom) = xstep%Rm(k, 1:sim%naesmd%natom) &
                    - xstep%R(k, 1:sim%naesmd%natom)
            end do
            if (sim%naesmd%dynam_type .eq. 'tsh') then
                call nacT_analytic(sim, sim%naesmd%cadiabnew, xstep)
            elseif (sim%naesmd%dynam_type .eq. 'aimc' .or. sim%naesmd%dynam_type .eq. 'mf') then
                do j = 1, sim%excN
                    sim%naesmd%cadiabnew(j, j) = 0.0d0
                    do k = j + 1, sim%excN
                        call nacR_analytic_wrap(sim, j, k, NACR)
                        NACR = sim%naesmd%sgn(j, k)*NACR
                        sim%naesmd%cadiabnew(j, k) = 0.0d0
                        do l = 1, sim%naesmd%natom
                            sim%naesmd%cadiabnew(j, k) = sim%naesmd%cadiabnew(j, k) + sim%naesmd%vxold(l)*NACR(3*l - 2)*convl + &
                                sim%naesmd%vyold(l)*NACR(3*l - 1)*convl + sim%naesmd%vzold(l)*NACR(3*l)*convl
                        end do
                        sim%naesmd%cadiabnew(k, j) = -1.0d0*sim%naesmd%cadiabnew(j, k)
                    end do
                end do
            end if
        end if

        do k = 1, 3
            deltaRpold(k)%p(1:sim%naesmd%natom) = deltaRp(k)%p(:sim%naesmd%natom)
            deltaRmold(k)%p(1:sim%naesmd%natom) = deltaRm(k)%p(:sim%naesmd%natom)
        end do

        do j = 1, sim%excN
            sim%naesmd%vmdqtold(j) = sim%naesmd%vmdqtnew(j)
        end do

        do i = 1, sim%dav%Ncis
            do j = 1, sim%excN
                sim%naesmd%cmdqtold(i, j) = sim%naesmd%cmdqtnew(i, j)
            end do
        end do

        do j = 1, sim%excN
            do k = 1, sim%excN
                sim%naesmd%cadiabold(j, k) = sim%naesmd%cadiabnew(j, k)
            end do
        end do

        deallocate (xs%Rp)
        deallocate (xs%Rm)
        deallocate (xs%R)

        return
    end subroutine cadiaboldcalc
    !
    !********************************************************************
    !
    !  Adiabatic middle calculation
    !
    !********************************************************************
    !
    subroutine cadiabmiddlecalc(sim, iimdqt, Na, cross)
        use qm2_davidson_module
        implicit none
        type(simulation_t), pointer::sim
        integer j, i, iimdqt, Na, k, l
        _REAL_ xx(Na), yy(Na), zz(Na)
        _REAL_ xxp(Na), yyp(Na), zzp(Na)
        _REAL_ xxm(Na), yym(Na), zzm(Na)
        integer cross(sim%excN)
        type(xstep_t), target::xs
        type(xstep_t), pointer::xstep
        type(simulation_t), pointer::simpoint
        _REAL_ NACR(sim%naesmd%natom*3)

        simpoint => sim

        allocate (xs%Rp(3, sim%Na))
        allocate (xs%Rm(3, sim%Na))
        allocate (xs%R(3, sim%Na))

        xstep => xs

        if (iimdqt .eq. 1) then
            do i = 1, sim%excN
                do j = 1, sim%excN
                    sim%naesmd%cadiabmiddleold(i, j) = sim%naesmd%cadiabold(i, j)
                end do
                sim%naesmd%vmdqtmiddleold(i) = sim%naesmd%vmdqtold(i)
            end do
            do i = 1, sim%dav%Ncis
                do j = 1, sim%excN
                    sim%naesmd%cmdqtmiddleold(i, j) = sim%naesmd%cmdqtold(i, j)
                end do
            end do
        else
            do i = 1, sim%excN
                do j = 1, sim%excN
                    sim%naesmd%cadiabmiddleold(i, j) = sim%naesmd%cadiabmiddle(i, j)
                end do
                sim%naesmd%vmdqtmiddleold(i) = sim%naesmd%vmdqtmiddle(i)
            end do
            do i = 1, sim%dav%Ncis
                do j = 1, sim%excN
                    sim%naesmd%cmdqtmiddleold(i, j) = sim%naesmd%cmdqtmiddle(i, j)
                end do
            end do
        end if
        if (iimdqt .eq. sim%naesmd%nquantumstep) then
            do i = 1, sim%excN
                sim%naesmd%vmdqtmiddle(i) = sim%naesmd%vmdqtnew(i)
            end do
            do i = 1, sim%excN
                do j = 1, sim%excN
                    sim%naesmd%cadiabmiddle(i, j) = sim%naesmd%cadiabnew(i, j)
                end do
            end do
            do i = 1, sim%dav%Ncis
                do j = 1, sim%excN
                    sim%naesmd%cmdqtmiddle(i, j) = sim%naesmd%cmdqtnew(i, j)
                end do
            end do
        else
            if (sim%naesmd%ensemble .eq. 'energy' .or. sim%naesmd%ensemble .eq. 'temper') then
                do j = 1, sim%naesmd%natom

                    xx(j) = sim%naesmd%rxold(j) + sim%naesmd%vxold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%axold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2

                    yy(j) = sim%naesmd%ryold(j) + sim%naesmd%vyold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%ayold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2

                    zz(j) = sim%naesmd%rzold(j) + sim%naesmd%vzold(j)*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%azold(j)*0.5d0*(sim%naesmd%dtquantum*dfloat(iimdqt))**2
                end do
            else if (sim%naesmd%ensemble .eq. 'langev') then
                do j = 1, sim%naesmd%natom
                    xx(j) = sim%naesmd%rxold(j) &
                        + sim%naesmd%vxold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%axold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        + sim%naesmd%prand(1, j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)

                    yy(j) = sim%naesmd%ryold(j) &
                        + sim%naesmd%vyold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%ayold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        + sim%naesmd%prand(2, j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)

                    zz(j) = sim%naesmd%rzold(j) &
                        + sim%naesmd%vzold(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt) &
                        + sim%naesmd%azold(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt) &
                        *(sim%naesmd%dtquantum*dfloat(iimdqt))**2 &
                        + sim%naesmd%prand(3, j)/sim%naesmd%dtmdqt*sim%naesmd%dtquantum*dfloat(iimdqt)
                end do
            end if
            call do_sqm_davidson_update(sim, 0, cmdqt=sim%naesmd%cmdqtmiddle, &
                vmdqt=sim%naesmd%vmdqtmiddle, vgs=sim%naesmd%vgs, rx=xx, ry=yy, rz=zz)

            ! xxp,yyp, and zzp are xyz at t + sim%naesmd%dtnact
            ! xxm,xym, and zzm are xyz at t - sim%naesmd%dtnact

            do j = 1, sim%naesmd%natom
                xxp(j) = xx(j) + sim%naesmd%deltaxxpold(j) &
                    + (sim%naesmd%deltaxxpnew(j) - sim%naesmd%deltaxxpold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)

                yyp(j) = yy(j) + sim%naesmd%deltayypold(j) &
                    + (sim%naesmd%deltayypnew(j) - sim%naesmd%deltayypold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)

                zzp(j) = zz(j) + sim%naesmd%deltazzpold(j) &
                    + (sim%naesmd%deltazzpnew(j) - sim%naesmd%deltazzpold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)

                xxm(j) = xx(j) + sim%naesmd%deltaxxmold(j) &
                    + (sim%naesmd%deltaxxmnew(j) - sim%naesmd%deltaxxmold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)

                yym(j) = yy(j) + sim%naesmd%deltayymold(j) &
                    + (sim%naesmd%deltayymnew(j) - sim%naesmd%deltayymold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)

                zzm(j) = zz(j) + sim%naesmd%deltazzmold(j) &
                    + (sim%naesmd%deltazzmnew(j) - sim%naesmd%deltazzmold(j))/sim%naesmd%dtmdqt &
                    *sim%naesmd%dtquantum*dfloat(iimdqt)
            end do

            call new_xstep(sim, xx, yy, zz, xxp, yyp, zzp, xxm, yym, zzm, xstep)
            if (sim%naesmd%dynam_type .eq. 'tsh') then
                call nacT_analytic(sim, sim%naesmd%cadiab, xstep)
            elseif (sim%naesmd%dynam_type .eq. 'aimc' .or. sim%naesmd%dynam_type .eq. 'mf') then
                do j = 1, sim%excN
                    do k = j + 1, sim%excN
                        call nacR_analytic_wrap(sim, j, k, NACR)
                        NACR = sim%naesmd%sgn(j, k)*NACR
                        sim%naesmd%cadiab(j, k) = 0.0d0
                        do l = 1, sim%naesmd%natom
                            sim%naesmd%cadiab(j, k) = sim%naesmd%cadiab(j, k) + sim%naesmd%vxold(l)*NACR(3*l - 2)*convl + &
                                sim%naesmd%vyold(l)*NACR(3*l - 1)*convl + sim%naesmd%vzold(l)*NACR(3*l)*convl
                        end do
                        sim%naesmd%cadiab(k, j) = -1.0d0*sim%naesmd%cadiab(j, k)
                    end do
                end do
            end if
            do i = 1, sim%excN
                do j = 1, sim%excN
                    sim%naesmd%cadiabmiddle(i, j) = sim%naesmd%cadiab(i, j)
                end do
            end do

        end if

        do i = 1, sim%excN
            if (sim%naesmd%dynam_type .eq. 'aimc' .or. sim%naesmd%dynam_type .eq. 'mf') then
                if (cross(i) .eq. 2) then
                    sim%naesmd%cadiabmiddle(i, sim%naesmd%iorden(i)) = 0.0d0
                    sim%naesmd%cadiabmiddle(sim%naesmd%iorden(i), i) = 0.0d0
                end if
            elseif (i .eq. sim%naesmd%ihop) then
                if (cross(i) .eq. 2) then
                    !if(sim%naesmd%conthop.gt.0) then
                    !    if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop=0
                    !end if

                    !if(sim%naesmd%conthop2.gt.0) then
                    !    if(sim%naesmd%iordenhop(i).ne.sim%naesmd%ihopprev) sim%naesmd%conthop2=0
                    !end if

                    if (sim%naesmd%conthop .eq. 0) then
                        sim%naesmd%cadiabmiddle(i, sim%naesmd%iorden(i)) = 0.d0
                        sim%naesmd%cadiabmiddle(sim%naesmd%iorden(i), i) = 0.d0
                    end if
                end if
            else
                if (i .ne. sim%naesmd%iorden(sim%naesmd%ihop)) then
                    if (i .lt. sim%naesmd%iorden(i)) then
                        if (cross(i) .eq. 2) then
                            sim%naesmd%cadiabmiddle(i, sim%naesmd%iorden(i)) = 0.d0
                            sim%naesmd%cadiabmiddle(sim%naesmd%iorden(i), i) = 0.d0
                        end if
                    end if
                end if
            end if
        end do

        if (iimdqt .eq. sim%naesmd%nquantumstep) then
            do i = 1, sim%excN
                do j = 1, sim%excN
                    sim%naesmd%cadiabnew(i, j) = sim%naesmd%cadiabmiddle(i, j)
                end do
            end do
        end if

        deallocate (xs%Rp)
        deallocate (xs%Rm)
        deallocate (xs%R)

        return
    end subroutine cadiabmiddlecalc
    !
    !********************************************************************
    !

    subroutine cadiabnewcalc(sim)
        use qm2_davidson_module
        implicit none

        type(simulation_t), pointer :: sim
        integer k, j, i, l
        type(xstep_t), target :: xs
        type(xstep_t), pointer::xstep
        type(realp_t) :: r(3)
        type(realp_t) :: deltaRp(3), deltaRm(3)
        type(simulation_t), pointer::simpoint
        _REAL_ NACR(sim%naesmd%natom*3)

        simpoint => sim

        r(1)%p => sim%naesmd%rx
        r(2)%p => sim%naesmd%ry
        r(3)%p => sim%naesmd%rz
        deltaRp(1)%p => sim%naesmd%deltaxxpnew
        deltaRp(2)%p => sim%naesmd%deltayypnew
        deltaRp(3)%p => sim%naesmd%deltazzpnew
        deltaRm(1)%p => sim%naesmd%deltaxxmnew
        deltaRm(2)%p => sim%naesmd%deltayymnew
        deltaRm(3)%p => sim%naesmd%deltazzmnew
        allocate (xs%Rp(3, sim%Na))
        allocate (xs%Rm(3, sim%Na))
        allocate (xs%R(3, sim%Na))

        xstep => xs

        call do_sqm_davidson_update(sim, sim%naesmd%printTdipole, cmdqt=sim%naesmd%cmdqtnew, &
            vmdqt=sim%naesmd%vmdqtnew, vgs=sim%naesmd%vgs)

        !  xxp,yyp, and zzp are xyz at t + sim%naesmd%dtnact
        !  xxm,yym, and zzm are xyz at t - sim%naesmd%dtnact

        call new_xstep_dtnact_r3(sim, r, xstep)

        do k = 1, 3
            deltaRp(k)%p(:sim%naesmd%natom) = xstep%Rp(k, :sim%naesmd%natom) &
                - xstep%R(k, :sim%naesmd%natom)
            deltaRm(k)%p(:sim%naesmd%natom) = xstep%Rm(k, :sim%naesmd%natom) &
                - xstep%R(k, :sim%naesmd%natom)
        end do

        if (sim%naesmd%dynam_type .eq. 'tsh') then
            call nacT_analytic(sim, sim%naesmd%cadiab, xstep)
        elseif (sim%naesmd%dynam_type .eq. 'aimc' .or. sim%naesmd%dynam_type .eq. 'mf') then
            do j = 1, sim%excN
                do k = j + 1, sim%excN
                    call nacR_analytic_wrap(sim, j, k, NACR)
                    NACR = sim%naesmd%sgn(j, k)*NACR
                    sim%naesmd%cadiab(j, k) = 0.0d0
                    do l = 1, sim%naesmd%natom
                        sim%naesmd%cadiab(j, k) = sim%naesmd%cadiab(j, k) + sim%naesmd%vxold(l)*NACR(3*l - 2)*convl + &
                            sim%naesmd%vyold(l)*NACR(3*l - 1)*convl + sim%naesmd%vzold(l)*NACR(3*l)*convl
                    end do
                    sim%naesmd%cadiab(k, j) = -1.0d0*sim%naesmd%cadiab(j, k)
                end do
            end do
        end if
        do i = 1, sim%excN
            do j = 1, sim%excN
                sim%naesmd%cadiabnew(i, j) = sim%naesmd%cadiab(i, j)
            end do
        end do

        deallocate (xs%Rp)
        deallocate (xs%Rm)
        deallocate (xs%R)

        return
    end subroutine cadiabnewcalc

end module cadiab_module
