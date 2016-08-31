! Subroutine to loose coherence

! modified by Seba
!        SUBROUTINE coherence(ido,neq,tini,tend,toldivprk, &
!      param,yg,constcoherE0,constcoherC)
   subroutine coherence(sim,Na,Nm,mdflag,d,E0,Omega,fosc, &
      ido,neq,tini,tend,toldivprk, &
      param,yg,constcoherE0,constcoherC,cohertype,idocontrol)
! end modified by Seba
   use qm2_davidson_module
   use communism
   
   implicit none

   type(simulation_t),pointer::sim

! modified by Seba
   integer cohertype,idocontrol
! end modified by Seba

! modified by Seba
!        integer Na
   integer mdflag,Na,Nm,kk
   double precision E0,d
   include 'md.par'
   include 'parH.par'
   real*8 Omega(sim%excN),fosc(sim%excN)
   real*8 xx(Na),yy(Na),zz(Na)
! end modified by Seba
   integer k,i,j,icheck,itest,ini,ihopavant
   integer ido,neq
   include 'sizes'
! modified by Seba
   include 'md.cmn'
! end modified by Seba
   double precision tini,tend,toldivprk,param(50)
   double precision iseedhop,eavant, eapres 
   double precision constcoherE0,constcoherC 
   include 'common'
   double precision yg(sim%excN),ytemp,ytemp2 
   double precision taocoher(sim%excN) 
! modified by Seba
!        double precision norm 
   double precision norm,norm1,norm2,normdij
   double precision vect1(Na*3) &
                   ,vect2(Na*3) &
                   ,vecs(Na*3)
   double precision dij(Na*3),kinec(sim%excN)
! end modified by Seba

   external fcn

! modified by Seba
! We need to modify it in order to call the proper ceo and nacR_analitic subroutines
!
!        kin=0.0d0
!
!        do i =1, natom
!            kin = kin + massmdqt(I)*(VX(I)**2+VY(I)**2+VZ(I)**2)
!        enddo
!        kin = 0.5d0 * kin
!
!        do j=1,sim%excN
!           if(j.ne.ihop) then
!             taocoher(j)=1.0d0/dabs((vmdqtnew(j)-vmdqtnew(ihop))) &
!      *(constcoherC+constcoherE0/kin)
!           endif
!        enddo
!
!        do j=1,sim%excN
!           if(j.ne.ihop) then
!              if(yg(j).ne.0.0d0) then
!                  yg(j)=yg(j)*dexp(-dtmdqt/taocoher(j))
!              endif
!           endif
!        enddo
!
!        norm=0.0d0
!        do k = 1,sim%excN
!           if(k.ne.ihop) then
!              norm=norm+yg(k)*yg(k)
!           endif
!        enddo
!
!        yg(ihop)=yg(ihop)*dsqrt((1-norm)/(yg(ihop)*yg(ihop)))
!
!        ido=3
!        call divprk(ido,neq,fcn,tini,tend, &
!      toldivprk,param,yg)
!        ido=1

   do j=1,sim%excN
      kinec(j)=0.0d0
   end do

   if(cohertype.eq.0) then

      kin=0.0d0
      do i=1,natom
         kin=kin+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)/2
      end do

      do j=1,sim%excN
         if(j.ne.ihop) then
            taocoher(j)=1.0d0/dabs((vmdqtnew(j)-vmdqtnew(ihop))) &
               *(constcoherC+constcoherE0/kin)
         end if
      end do

      do j=1,sim%excN
         if(j.ne.ihop) then
            if(yg(j).ne.0.0d0) then
               yg(j)=yg(j)*dexp(-dtmdqt/taocoher(j))
            end if
         end if
      end do

   end if

   if(cohertype.eq.1) then

      do j=1,natom
         xx(j)=rx(j)*convl
         yy(j)=ry(j)*convl
         zz(j)=rz(j)*convl
      end do

      !mdflag=2
      !call ceo(Na,xx,yy,zz,atoms,sim%excN,E0,Omega,fosc,mdflag)
      sim%dav%mdflag=2
      call do_sqm_davidson_update(sim,cmdqt=cmdqt, &
         vmdqt=vmdqt,vgs=vgs,rx=xx,ry=yy,rz=zz)

      do j=1,sim%excN
         if(j.ne.ihop) then
            if(yg(j).ne.0.d0) then
! calculate the non-adiabatic coupling vector(nacR)
               !call nacR_analytic(Na,Nm,d,xx,yy,zz,sim%excN,Omega, &
               !   mdflag,ihop,j,cmdqt,dij)

               call nacR_analytic_wrap(sim,ihop,j,dij)
! calculate the magnitude of nacR ina a.u. 
               do k=1,natom*3
                  dij(k)=dij(k)*convl
               end do

               normdij=0.0d0
               do k=1,natom*3
                  normdij=normdij+dij(k)**2
               end do
!                    write(107,*) tfemto,j,normdij
!                    call flush(107)
               if(normdij.ne.0.0d0) then
                  normdij=dsqrt(normdij)
! inner product of normalized nacR and P
                  kk=1
                  norm=0.d0

                  do k=1,natom
                     norm=norm+massmdqt(k)*vx(k)*dij(kk)/normdij &
                        +massmdqt(k)*vy(k)*dij(kk+1)/normdij &
                        +massmdqt(k)*vz(k)*dij(kk+2)/normdij

                     kk=kk+3
                  end do
!                    write(108,*) tfemto,j,norm
!                    call flush(108)
! vector projection of P on nacR
                  do k=1,natom*3
                     dij(k)=dij(k)*norm
                  end do
! sum up of vector projection of P on nacR and vector P
! checking which sign leads to a summation
                  kk=1
                  do k=1,natom
                     vect1(kk)=dij(kk)+massmdqt(k)*vx(k)
                     vect1(kk+1)=dij(kk+1)+massmdqt(k)*vy(k)
                     vect1(kk+2)=dij(kk+2)+massmdqt(k)*vz(k)
                     kk=kk+3
                  end do

                  kk=1
                  do k=1,natom
                     vect2(kk)=dij(kk)-massmdqt(k)*vx(k)
                     vect2(kk+1)=dij(kk+1)-massmdqt(k)*vy(k)
                     vect2(kk+2)=dij(kk+2)-massmdqt(k)*vz(k)
                     kk=kk+3
                  end do

                  norm1=0.0d0
                  norm2=0.0d0
                  do kk=1,natom*3
                     norm1=norm1+vect1(kk)**2
                     norm2=norm2+vect2(kk)**2
                  end do

!                    write(109,*) tfemto,j,norm1,norm2
!                    call flush(109)
! final decoherence direction versor s
                  if(norm1.gt.norm2) then
                     do kk=1,natom*3
                        vecs(kk)=vect1(kk)/dsqrt(norm1)
                     end do
                  else
                     do kk=1,natom*3
                        vecs(kk)=vect2(kk)/dsqrt(norm2)
                     end do
                  end if
! projection of P on s
                  kk=1
                  norm=0.0d0
                  do k=1,natom
                     norm=norm+massmdqt(k)*vx(k)*vecs(kk) &
                        +massmdqt(k)*vy(k)*vecs(kk+1) &
                        +massmdqt(k)*vz(k)*vecs(kk+2)
                     kk=kk+3
                  end do
!                    write(110,*) tfemto,j,norm
!                    call flush(110)
! vector projection of P on s
                  do k=1,natom*3
                     vecs(k)=vecs(k)*norm
                  end do
! Kinetic energy from this projection
                  kin=0.0d0
                  kk=1
                  do k =1,natom
                     kin=kin+vecs(kk)**2/(2.0d0*massmdqt(k)) &
                        +vecs(kk+1)**2/(2.0d0*massmdqt(k)) &
                        +vecs(kk+2)**2/(2.0d0*massmdqt(k))
                     kk=kk+3
                  end do
                  kinec(j)=kin*feVmdqt

               else

                  kin=0.0d0
                  do k =1,natom
                     kin=kin+massmdqt(k)*(vx(k)**2+vy(k)**2+vz(k)**2)/2
                  end do
                  kinec(j)=kin*feVmdqt
               end if

!                    write(110,*) tfemto,j,kin
!                    call flush(110)

! caclulate decoherence time ofr state j

               taocoher(j)=1.0d0/dabs((vmdqtnew(j)-vmdqtnew(ihop))) &
                  *(constcoherC+constcoherE0/kin)

! apply decoherence for state j

               yg(j)=yg(j)*dexp(-dtmdqt/taocoher(j))

            end if
         end if
      end do

      write(111,888) tfemto,(kinec(k),k=1,sim%excN)
      call flush(111)
   end if

   norm=0.d0
   do k=1,sim%excN
      if(k.ne.ihop) then
         norm=norm+yg(k)*yg(k)
      end if
   end do

   yg(ihop)=yg(ihop)*dsqrt((1-norm)/(yg(ihop)*yg(ihop)))

   if(idocontrol.eq.0) then
      ido=3
      call divprk(ido,neq,fcn,tini,tend,toldivprk,param,yg)
      ido=1
   end if
! end modified by Seba


889   FORMAT(I3,10000(1X,F18.10))
888   FORMAT(10000(1X,F18.10))
887   FORMAT(F18.10,1X,I2,1X,I3,1X,I3,10000(1X,F18.10))

   return
   end
!
