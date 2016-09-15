#include "dprec.fh"
#include "assert.fh"

module cadiab_module
    use naesmd_constants
    use communism
    use nacT_analytic_module
contains


    subroutine cadiaboldcalc(sim,imdqt,Na,Nm)
        use naesmd_module
        use md_module
        use qm2_davidson_module
        implicit none
 
        type(simulation_t), pointer :: sim

        integer k,j,i,ii,iii,imdqt
        integer mdflag,Na,Nm
        _REAL_ E0
        type(xstep_t)::xstep
      
        if(imdqt.eq.1) then
            call do_sqm_davidson_update(sim,cmdqt=cmdqtnew, &
                vmdqt=vmdqtnew,vgs=vgs,r=sim%naesmd%r%vold)!, cmdqt=cmdqtnew) !JAKB, not necessary bc updated already in
            xstep=new_xstep_dtnact_r3(sim,sim%naesmd%r%vold)
            do k=1,3
                sim%naesmd%deltaRp%v(k)%p(1:natom)=xstep%Rp(k,1:natom) &
                    -xstep%R(k,1:natom)
                sim%naesmd%deltaRm%v(k)%p(:natom)=xstep%Rm(k,1:natom) &
                    -xstep%R(k,1:natom)
            end do
            call nacT_analytic(sim,cadiabnew,xstep)
        end if

        do k=1,3
            sim%naesmd%deltaRp%vold(k)%p(1:natom)=sim%naesmd%deltaRp%v(k)%p(:natom)
            sim%naesmd%deltaRm%vold(k)%p(1:natom)=sim%naesmd%deltaRm%v(k)%p(:natom)
        end do

        do j=1,qm2ds%Mx
            vmdqtold(j)=vmdqtnew(j)
        enddo

        do i=1,sim%dav%Ncis
            do j=1,qm2ds%Mx
                cmdqtold(i,j)=cmdqtnew(i,j)
            enddo
        enddo

        do j=1,qm2ds%Mx
            do k=1,qm2ds%Mx
                cadiabold(j,k)=cadiabnew(j,k)
            end do
        end do

        return
    end subroutine
    !
    !********************************************************************
    !
    !  Adiabatic middle calculation
    !
    !********************************************************************
    !
    subroutine cadiabmiddlecalc(sim,iimdqt,Na,Nm,cross)
        use qm2_davidson_module
        use naesmd_module
        use md_module
        implicit none
        type(simulation_t),pointer::sim
        integer k,j,i,ii,iii,iimdqt
        integer mdflag,Na,Nm
        real(8) E0
        real*8 Omega(qm2ds%Mx),fosc(qm2ds%Mx)
        real*8 xx(Na),yy(Na),zz(Na)
        real*8 xxp(Na),yyp(Na),zzp(Na)
        real*8 xxm(Na),yym(Na),zzm(Na)
        integer cross(sim%excN)
        type(xstep_t)::xstep

        if(iimdqt.eq.1) then
            do i=1,qm2ds%Mx
                do j=1,qm2ds%Mx
                    cadiabmiddleold(i,j)=cadiabold(i,j)
                end do
                vmdqtmiddleold(i)=vmdqtold(i)
            end do
            do i=1,sim%dav%Ncis
                do j=1,qm2ds%Mx
                    cmdqtmiddleold(i,j)=cmdqtold(i,j)
                end do
            end do
        else
            do i=1,qm2ds%Mx
                do j=1,qm2ds%Mx
                    cadiabmiddleold(i,j)=cadiabmiddle(i,j)
                end do
                vmdqtmiddleold(i)=vmdqtmiddle(i)
            end do
            do i=1,sim%dav%Ncis
                do j=1,qm2ds%Mx
                    cmdqtmiddleold(i,j)=cmdqtmiddle(i,j)
                end do
            end do
        end if
        if(iimdqt.eq.nquantumstep) then
            do i=1,qm2ds%Mx
                vmdqtmiddle(i)=vmdqtnew(i)
            end do
            do i=1,qm2ds%Mx
                do j=1,qm2ds%Mx
                    cadiabmiddle(i,j)=cadiabnew(i,j)
                end do
            end do
            do i=1,sim%dav%Ncis
                do j=1,qm2ds%Mx
                    cmdqtmiddle(i,j)=cmdqtnew(i,j)
                end do
            end do
        else
            if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
                do j=1,natom

                    xx(j)=rxold(j) + vxold(j)*dtquantum*dfloat(iimdqt) &
                        +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

                    yy(j)=ryold(j) + vyold(j)*dtquantum*dfloat(iimdqt) &
                        +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

                    zz(j)=rzold(j) + vzold(j)*dtquantum*dfloat(iimdqt) &
                        +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2
                end do
            else if (ensemble.eq.'langev') then
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
            endif
            call do_sqm_davidson_update(sim,cmdqt=cmdqtmiddle, &
                vmdqt=vmdqtmiddle,vgs=vgs,rx=xx,ry=yy,rz=zz)

            ! xxp,yyp, and zzp are xyz at t + dtnact
            ! xxm,xym, and zzm are xyz at t - dtnact

            do j=1,natom
                xxp(j)=xx(j)+deltaxxpold(j) &
                    +(deltaxxpnew(j)-deltaxxpold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)

                yyp(j)=yy(j)+deltayypold(j) &
                    +(deltayypnew(j)-deltayypold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)

                zzp(j)=zz(j)+deltazzpold(j) &
                    +(deltazzpnew(j)-deltazzpold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)

                xxm(j)=xx(j)+deltaxxmold(j) &
                    +(deltaxxmnew(j)-deltaxxmold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)

                yym(j)=yy(j)+deltayymold(j) &
                    +(deltayymnew(j)-deltayymold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)

                zzm(j)=zz(j)+deltazzmold(j) &
                    +(deltazzmnew(j)-deltazzmold(j))/dtmdqt &
                    *dtquantum*dfloat(iimdqt)
            end do

            xstep=new_xstep(sim,xx,yy,zz,xxp,yyp,zzp,xxm,yym,zzm)
            call nacT_analytic(sim,cadiab,xstep)

            do i=1,qm2ds%Mx
                do j=1,qm2ds%Mx
                    cadiabmiddle(i,j)=cadiab(i,j)
                end do
            end do

        end if

        do i=1,qm2ds%Mx
            if(i.eq.ihop) then
                if(cross(i).eq.2) then
                    if(conthop.gt.0) then
                        if(iordenhop(i).ne.ihopprev) conthop=0
                    end if

                    if(conthop2.gt.0) then
                        if(iordenhop(i).ne.ihopprev) conthop2=0
                    end if

                    if(conthop.eq.0) then
                        cadiabmiddle(i,iorden(i))=0.d0
                        cadiabmiddle(iorden(i),i)=0.d0
                    end if
                end if
            else
                if(i.ne.iorden(ihop)) then
                    if(i.lt.iorden(i)) then
                        if(cross(i).eq.2) then
                            cadiabmiddle(i,iorden(i))=0.d0
                            cadiabmiddle(iorden(i),i)=0.d0
                        end if
                    end if
                end if
            end if
        end do

        if(iimdqt.eq.nquantumstep) then
            do i=1,qm2ds%Mx
                do j=1,qm2ds%Mx
                    cadiabnew(i,j)=cadiabmiddle(i,j)
                end do
            end do
        end if

889     FORMAT(10000(1X,F18.10))

        return
    end subroutine
    !
    !********************************************************************
    !

    subroutine cadiabnewcalc(sim,Na,Nm)
        use naesmd_module
        use md_module
        use qm2_davidson_module
        implicit none

        type(simulation_t), pointer :: sim
        integer k,j,i,ii,iii,iimdqt
        integer mdflag,Na,Nm
        double precision E0

        type(xstep_t) :: xstep
        call do_sqm_davidson_update(sim,cmdqt=cmdqtnew, &
            vmdqt=vmdqtnew,vgs=vgs)

        !  xxp,yyp, and zzp are xyz at t + dtnact
        !  xxm,yym, and zzm are xyz at t - dtnact

        xstep=new_xstep_dtnact_r3(sim, sim%naesmd%r%v)

        do k=1,3
            sim%naesmd%deltaRp%v(k)%p(:natom)=xstep%Rp(k,:natom) &
                -xstep%R( k,:natom)
            sim%naesmd%deltaRm%v(k)%p(:natom)=xstep%Rm(k,:natom) &
                -xstep%R( k,:natom)
        end do

        call nacT_analytic(sim,cadiab,xstep)

        do i=1,qm2ds%Mx
            do j=1,qm2ds%Mx
                cadiabnew(i,j)=cadiab(i,j)
            end do
        end do

889     format(1000(1x,f18.10))

        return
    end subroutine

end module
