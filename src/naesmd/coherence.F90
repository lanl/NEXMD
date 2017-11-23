#include "dprec.fh"
#include "assert.fh"

! Subroutine to loose coherence

subroutine coherence(sim,rkcomm, Na, tend, tmax,rk_tolerance, &
    thresholds,yg,constcoherE0,constcoherC,cohertype)
!subroutine coherence(sim,Na,ido,neq,tini,tend,toldivprk, &
!    param,yg,constcoherE0,constcoherC,cohertype,idocontrol)
    use rksuite_90, only: rk_comm_real_1d, setup 
    use qm2_davidson_module
    use communism
    use md_module
    implicit none

    type(simulation_t),pointer  ::  sim
    type(rk_comm_real_1d) :: rkcomm 
   
    integer cohertype
    !integer cohertype,idocontrol
    integer Na,kk
    _REAL_ xx(Na),yy(Na),zz(Na)
    integer k,i,j
    !integer ido,neq
    !_REAL_ tini,tend,toldivprk,param(50)
    _REAL_ tend,tmax,rk_tolerance
    _REAL_ constcoherE0,constcoherC
    _REAL_ yg(sim%excN*2), thresholds(sim%excN*2) 
    _REAL_ taocoher(sim%excN)
    _REAL_ norm,norm1,norm2,normdij
    _REAL_ vect1(Na*3) &
        ,vect2(Na*3) &
        ,vecs(Na*3)
    _REAL_ kinec(sim%excN)
    _REAL_ dij(Na*3)

    !external fcn

    do j=1,sim%excN
        kinec(j)=0.0d0
    end do

    if(cohertype.eq.0) then

        sim%naesmd%kin=0.0d0
        do i=1,sim%naesmd%natom
            sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(i)*(sim%naesmd%vx(i)**2+sim%naesmd%vy(i)**2+sim%naesmd%vz(i)**2)/2
        end do

        do j=1,sim%excN
            if(j.ne.sim%naesmd%ihop) then
                taocoher(j)=1.0d0/dabs((sim%naesmd%vmdqtnew(j)-sim%naesmd%vmdqtnew(sim%naesmd%ihop))) &
                    *(constcoherC+constcoherE0/sim%naesmd%kin)
            end if
        end do

        do j=1,sim%excN
            if(j.ne.sim%naesmd%ihop) then
                if(yg(j).ne.0.0d0) then
                    yg(j)=yg(j)*dexp(-sim%naesmd%dtmdqt/taocoher(j))
                end if
            end if
        end do

    end if

    if(cohertype.eq.1) then

        do j=1,sim%naesmd%natom
            xx(j)=sim%naesmd%rx(j)*convl
            yy(j)=sim%naesmd%ry(j)*convl
            zz(j)=sim%naesmd%rz(j)*convl
        end do

        sim%dav%mdflag=2
        call do_sqm_davidson_update(sim,cmdqt=sim%naesmd%cmdqt, &
            vmdqt=sim%naesmd%vmdqt,vgs=sim%naesmd%vgs,rx=xx,ry=yy,rz=zz)

        do j=1,sim%excN
            if(j.ne.sim%naesmd%ihop) then
                if(yg(j).ne.0.d0) then
                    call nacR_analytic_wrap(sim,sim%naesmd%ihop,j,dij)
                    ! calculate the magnitude of nacR ina a.u.
                    do k=1,sim%naesmd%natom*3
                        dij(k)=dij(k)*convl
                    end do

                    normdij=0.0d0
                    do k=1,sim%naesmd%natom*3
                        normdij=normdij+dij(k)**2
                    end do
                    if(normdij.ne.0.0d0) then
                        normdij=dsqrt(normdij)
                        ! inner product of normalized nacR and P
                        kk=1
                        norm=0.d0

                        do k=1,sim%naesmd%natom
                            norm=norm+sim%naesmd%massmdqt(k)*sim%naesmd%vx(k)*dij(kk)/normdij &
                                +sim%naesmd%massmdqt(k)*sim%naesmd%vy(k)*dij(kk+1)/normdij &
                                +sim%naesmd%massmdqt(k)*sim%naesmd%vz(k)*dij(kk+2)/normdij

                            kk=kk+3
                        end do
                        ! vector projection of P on nacR
                        do k=1,sim%naesmd%natom*3
                            dij(k)=dij(k)*norm
                        end do
                        ! sum up of vector projection of P on nacR and vector P
                        ! checking which sign leads to a summation
                        kk=1
                        do k=1,sim%naesmd%natom
                            vect1(kk)=dij(kk)+sim%naesmd%massmdqt(k)*sim%naesmd%vx(k)
                            vect1(kk+1)=dij(kk+1)+sim%naesmd%massmdqt(k)*sim%naesmd%vy(k)
                            vect1(kk+2)=dij(kk+2)+sim%naesmd%massmdqt(k)*sim%naesmd%vz(k)
                            kk=kk+3
                        end do

                        kk=1
                        do k=1,sim%naesmd%natom
                            vect2(kk)=dij(kk)-sim%naesmd%massmdqt(k)*sim%naesmd%vx(k)
                            vect2(kk+1)=dij(kk+1)-sim%naesmd%massmdqt(k)*sim%naesmd%vy(k)
                            vect2(kk+2)=dij(kk+2)-sim%naesmd%massmdqt(k)*sim%naesmd%vz(k)
                            kk=kk+3
                        end do

                        norm1=0.0d0
                        norm2=0.0d0
                        do kk=1,sim%naesmd%natom*3
                            norm1=norm1+vect1(kk)**2
                            norm2=norm2+vect2(kk)**2
                        end do

                        ! final decoherence direction versor s
                        if(norm1.gt.norm2) then
                            do kk=1,sim%naesmd%natom*3
                                vecs(kk)=vect1(kk)/dsqrt(norm1)
                            end do
                        else
                            do kk=1,sim%naesmd%natom*3
                                vecs(kk)=vect2(kk)/dsqrt(norm2)
                            end do
                        end if
                        ! projection of P on s
                        kk=1
                        norm=0.0d0
                        do k=1,sim%naesmd%natom
                            norm=norm+sim%naesmd%massmdqt(k)*sim%naesmd%vx(k)*vecs(kk) &
                                +sim%naesmd%massmdqt(k)*sim%naesmd%vy(k)*vecs(kk+1) &
                                +sim%naesmd%massmdqt(k)*sim%naesmd%vz(k)*vecs(kk+2)
                            kk=kk+3
                        end do
                        ! vector projection of P on s
                        do k=1,sim%naesmd%natom*3
                            vecs(k)=vecs(k)*norm
                        end do
                        ! Kinetic energy from this projection
                        sim%naesmd%kin=0.0d0
                        kk=1
                        do k =1,sim%naesmd%natom
                            sim%naesmd%kin=sim%naesmd%kin+vecs(kk)**2/(2.0d0*sim%naesmd%massmdqt(k)) &
                                +vecs(kk+1)**2/(2.0d0*sim%naesmd%massmdqt(k)) &
                                +vecs(kk+2)**2/(2.0d0*sim%naesmd%massmdqt(k))
                            kk=kk+3
                        end do
                        kinec(j)=sim%naesmd%kin*feVmdqt

                    else

                        sim%naesmd%kin=0.0d0
                        do k =1,sim%naesmd%natom
                            sim%naesmd%kin=sim%naesmd%kin+ &
			sim%naesmd%massmdqt(k)*(sim%naesmd%vx(k)**2+sim%naesmd%vy(k)**2+sim%naesmd%vz(k)**2)/2
                        end do
                        kinec(j)=sim%naesmd%kin*feVmdqt
                    end if

                    ! caclulate decoherence time ofr sim%naesmd%state j

                    taocoher(j)=1.0d0/dabs((sim%naesmd%vmdqtnew(j)-sim%naesmd%vmdqtnew(sim%naesmd%ihop))) &
                        *(constcoherC+constcoherE0/sim%naesmd%kin)

                    ! apply decoherence for sim%naesmd%state j

                    yg(j)=yg(j)*dexp(-sim%naesmd%dtmdqt/taocoher(j))

                end if
            end if
        end do

        write(111,888) sim%naesmd%tfemto,(kinec(k),k=1,sim%excN)
        call flush(111)
    end if

    norm=0.d0
    do k=1,sim%excN
        if(k.ne.sim%naesmd%ihop) then
            norm=norm+yg(k)*yg(k)
        end if
    end do

    yg(sim%naesmd%ihop)=yg(sim%naesmd%ihop)*dsqrt((1-norm)/(yg(sim%naesmd%ihop)*yg(sim%naesmd%ihop)))

    call    setup(rkcomm, tend, yg, tmax, rk_tolerance, thresholds, &
	    'M','R')
    !if(idocontrol.eq.0) then
    !    ido=3
    !    call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
    !    ido=1
    !    idocontrol=1
    !end if
    ! end modified by Seba


888 FORMAT(10000(1X,F18.10))
    return
end
!
