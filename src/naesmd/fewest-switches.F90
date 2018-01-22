#include "dprec.fh"
#include "assert.fh"

! Subroutine to calculate the probability to hop and the velocity adjustment after an effective hop
module fewest_switches
    use naesmd_constants
    use random
    use langevin_temperature
    use communism
    implicit none
contains
    subroutine evalhop(sim, rkcomm, lprint, tend,tmax,rk_tolerance,thresholds, &
        Na,yg,cross)
	use rksuite_90, only: rk_comm_real_1d, setup 
        implicit none
        type(simulation_t), pointer :: sim
        type(rk_comm_real_1d) :: rkcomm 
        integer Na,lprint,excNtemp,ipn,indx(Na)
        integer k,i,j,jj,icheck,itest,ini,ihopavant
        _REAL_ tmax,tend,rk_tolerance,pn,kavant
        _REAL_ g(sim%excN),gacum(sim%excN)
        _REAL_ iseedhop,eavant, eapres
        _REAL_ xx(Na),yy(Na),zz(Na)
        _REAL_ yg(sim%excN*2),ytemp,ytemp2, thresholds(sim%excN*2)
        _REAL_ t_start,t_finish
        integer cross(sim%excN),crosstemp,ininonhop
        if(sim%naesmd%conthop.eq.3) sim%naesmd%conthop=0
        if(sim%naesmd%conthop2.eq.3) sim%naesmd%conthop2=0
        if(sim%naesmd%conthop.gt.0) sim%naesmd%conthop=sim%naesmd%conthop+1
        if(sim%naesmd%conthop2.gt.0) sim%naesmd%conthop2=sim%naesmd%conthop2+1
        if(sim%naesmd%conthop.gt.0) then
            if(cross(sim%naesmd%ihop).eq.2) then
                if(sim%naesmd%iordenhop(sim%naesmd%ihop).ne.sim%naesmd%ihopprev) sim%naesmd%conthop=0
            end if
        end if
        iseedhop=rranf1(sim%naesmd%iseedmdqt)
        ihopavant=sim%naesmd%ihop
        eavant=sim%naesmd%vmdqtnew(sim%naesmd%ihop)
        do j=1,sim%naesmd%natom
            eavant=eavant+0.5d0*sim%naesmd%massmdqt(j)*(sim%naesmd%vx(j)**2+sim%naesmd%vy(j)**2+sim%naesmd%vz(j)**2)
        end do
        do j=1,sim%excN
            g(j)=0.d0
            gacum(j)=0.d0
        end do
        ! g(j) is the probability to hop from the
        ! current sim%naesmd%state to the sim%naesmd%state j
        do j=1,sim%excN
            if(j.ne.sim%naesmd%ihop) then
                g(j)=sim%naesmd%vnqcorrhoptot(j,sim%naesmd%ihop)/(sim%naesmd%nqold**2)
                if(g(j).lt.0.0d0) g(j)=0.0d0
            end if
        end do
        icheck=0
        do j=1,sim%excN
            gacum(j)=0.0d0
            if(j.ne.sim%naesmd%ihop) then
                do k=1,j
                    if(k.ne.sim%naesmd%ihop) gacum(j)=gacum(j)+g(k)
                end do
            end if
        end do
        itest=0
        do j=1,sim%excN
            if(j.ne.sim%naesmd%ihop) then
                if(iseedhop.le.gacum(j).and.itest.eq.0) then
                    icheck=j
                    itest=1
                end if
            end if
        end do
        crosstemp=cross(sim%naesmd%ihop)
        if(icheck.ne.0.and.cross(sim%naesmd%ihop).ne.2) then
            ini=0
            ! adjustment of velocities
            if(sim%naesmd%conthop2.gt.0) then
                if(sim%naesmd%ihopprev.ne.icheck) sim%naesmd%conthop2=0
            end if
            if(sim%naesmd%conthop2.eq.0) then
                call veladjustment(sim, lprint,Na,icheck,ini)
                if(lprint.ge.2) then
                    write(sim%outfile_10,*) sim%naesmd%tfemto,icheck,ini
                    call flush(sim%outfile_10)
                end if
            else
                ini=1
            end if
            if(ini.eq.0) then
                sim%naesmd%ihop=icheck
                call naesmd2qmmm_r(sim)
                call cpu_time(t_start)
                call deriv(sim,sim%naesmd%ihop)
                call cpu_time(t_finish)
                sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
                call do_sqm_davidson_update(sim,vmdqt=sim%naesmd%vmdqtnew,vgs=sim%naesmd%vgs)

!                if(sim%naesmd%decorhop.eq.1) then
!                    do j=1,sim%excN
!                        yg(j)=0.0d0
!                        yg(j+sim%excN)=rranf1(sim%naesmd%iseedmdqt)
!                    end do
!                    yg(sim%naesmd%ihop)=1.d0
!                end if

!Added for patch JAKB
               if(sim%naesmd%decorhop.ne.0) then
                  do j=1,sim%excN
                     yg(j)=0.0d0
!                     call random(sim%naesmd%iseedmdqt,iseedhop)
!                     yg(j+sim%excN)=iseedhop
                     yg(j+sim%excN)=rranf1(sim%naesmd%iseedmdqt)
                  enddo
                  yg(sim%naesmd%ihop)=1.0d0
	          call    setup(rkcomm, tend, yg, tmax, rk_tolerance, thresholds, 'M','R')
! after kirill
! check the reduction of sim%excN
                 if(sim%naesmd%iredpot.eq.1) then
                   if((sim%naesmd%ihop+sim%naesmd%nstates).lt.sim%excN) then
                      sim%excN=sim%naesmd%ihop+sim%naesmd%nstates
                      write(sim%outfile_23,887) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN
                      call flush(sim%outfile_23)
                   endif
                 endif
                 if(sim%naesmd%iredpot.eq.2) then
                   excNtemp=sim%excN
                   do i=1,sim%excN
                      if((sim%naesmd%vmdqtnew(sim%naesmd%ihop)+sim%naesmd%deltared).lt.sim%naesmd%vmdqtnew(sim%excN+1-i)) then
                         excNtemp=sim%excN+1-i
                      endif
                   enddo
                   sim%excN=excNtemp
                   write(sim%outfile_23,887) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN
                   call flush(sim%outfile_23)
                 endif
! after Josiah
                 if(sim%naesmd%iredpot.eq.3) then
                   jj=1
                   pn=0.0d0
                   do i=1,sim%naesmd%natom
                      if(sim%naesmd%atomtype(i).eq.1) then
                        sim%naesmd%cicoeffao3(i)=sim%naesmd%cicoeffao2(jj,sim%naesmd%ihop)**2
                        jj=jj+1
                      endif
                      if(sim%naesmd%atomtype(i).eq.6) then
                        sim%naesmd%cicoeffao3(i)=sim%naesmd%cicoeffao2(jj,sim%naesmd%ihop)**2 &
			    + sim%naesmd%cicoeffao2(jj+1,sim%naesmd%ihop)**2&
                            +sim%naesmd%cicoeffao2(jj+2,sim%naesmd%ihop)**2 + sim%naesmd%cicoeffao2(jj+3,sim%naesmd%ihop)**2
                          jj=jj+4
                      endif
                      pn=pn+ sim%naesmd%cicoeffao3(i)**2
                   enddo
                   pn=1.0d0/pn
                   ipn=int(pn)+1
                   call indexx(sim%naesmd%natom,sim%naesmd%cicoeffao3,indx)
                   sim%naesmd%deltared=0.0d0
                   do jj=1,sim%naesmd%natom
                      write(sim%outfile_24,*) ipn,sim%naesmd%atomtype(indx(jj)),sim%naesmd%cicoeffao3(indx(jj))
                   enddo
                   call flush(sim%outfile_24)
                   do jj=1,sim%naesmd%natom
                      if(indx(jj).ge.(sim%naesmd%natom-ipn)) then
                        sim%naesmd%deltared = sim%naesmd%deltared + &
			0.5d0*sim%naesmd%massmdqt(jj)*(sim%naesmd%vx(jj)**2 &
			+sim%naesmd%vy(jj)**2+sim%naesmd%vz(jj)**2)
                      endif
                   enddo
                   excNtemp=sim%excN
                   do i=1,sim%excN
                     if((sim%naesmd%vmdqtnew(sim%naesmd%ihop)+sim%naesmd%deltared).lt.sim%naesmd%vmdqtnew(sim%excN+1-i)) then
                         excNtemp=sim%excN+1-i
                     endif
                   enddo
                   sim%excN=excNtemp
                   write(sim%outfile_23,884) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN, &
				ipn,sim%naesmd%deltared*feVmdqt,kavant*feVmdqt
                   call flush(sim%outfile_23)
                 endif
! end after josiah
! after kirill
               endif
               sim%naesmd%conthop=1
! sim%naesmd%ihopprev keep the value of the previous sim%naesmd%state from where we hop
! in order to allow new hops for the next 2 steps only in case that
! we do not want to hop again to the same sim%naesmd%state
!               sim%naesmd%ihopprev=ihopavant
!######################################### 
            endif
!######################################### 
! added after Kirill
!######################################### 
           if (ini.eq.2) then
!######################################### 
               do j=1,sim%excN
                  yg(j)=0.0d0
                  yg(j+sim%excN)=rranf1(sim%naesmd%iseedmdqt)
               enddo
               yg(ihopavant)=1.0d0
		    call    setup(rkcomm, tend, yg, tmax, rk_tolerance, thresholds, &
			    'M','R')
               sim%naesmd%conthop=1
! check the reduction of sim%excN
               if(sim%naesmd%iredpot.eq.1) then
                   if((ihopavant+sim%naesmd%nstates).lt.sim%excN) then
                      sim%excN=ihopavant+sim%naesmd%nstates
                      write(sim%outfile_23,887) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN
                      call flush(sim%outfile_23)
                   endif
               endif
               if(sim%naesmd%iredpot.eq.2) then
                  excNtemp=sim%excN
                  do i=1,sim%excN
                     if((sim%naesmd%vmdqtnew(ihopavant)+sim%naesmd%deltared).lt.sim%naesmd%vmdqtnew(sim%excN+1-i)) then
                         excNtemp=sim%excN+1-i
                     endif
                  enddo
                  sim%excN=excNtemp
                  write(sim%outfile_23,887) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN
                  call flush(sim%outfile_23)
               endif
! after Josiah
               if(sim%naesmd%iredpot.eq.3) then
                  jj=1
                  pn=0.0d0
                  do i=1,sim%naesmd%natom
                     if(sim%naesmd%atomtype(i).eq.1) then
                        sim%naesmd%cicoeffao3(i)=sim%naesmd%cicoeffao2(jj,ihopavant)**2
                        jj=jj+1
                     endif
                     if(sim%naesmd%atomtype(i).eq.6) then
                        sim%naesmd%cicoeffao3(i)=sim%naesmd%cicoeffao2(jj,ihopavant)**2 &
      + sim%naesmd%cicoeffao2(jj+1,ihopavant)**2+sim%naesmd%cicoeffao2(jj+2,ihopavant)**2 &
      + sim%naesmd%cicoeffao2(jj+3,ihopavant)**2
                        jj=jj+4
                     endif
                     pn=pn+ sim%naesmd%cicoeffao3(i)**2
                  enddo
                  pn=1.0d0/pn
                  ipn=int(pn)+1
                  call indexx(sim%naesmd%natom,sim%naesmd%cicoeffao3,indx)
                  sim%naesmd%deltared=0.0d0
              do jj=1,sim%naesmd%natom
                write(sim%outfile_24,*) ipn,sim%naesmd%atomtype(indx(jj)),sim%naesmd%cicoeffao3(indx(jj))
             enddo
              call flush(sim%outfile_24)
                  do jj=1,sim%naesmd%natom
                     if(indx(jj).ge.(sim%naesmd%natom-ipn)) then
                        sim%naesmd%deltared = sim%naesmd%deltared + 0.5d0*sim%naesmd%massmdqt(jj)*&
                                (sim%naesmd%vx(jj)**2+sim%naesmd%vy(jj)**2+sim%naesmd%vz(jj)**2)
                     endif
                  enddo
                  excNtemp=sim%excN
                  do i=1,sim%excN
                     if((sim%naesmd%vmdqtnew(ihopavant)+sim%naesmd%deltared).lt.sim%naesmd%vmdqtnew(sim%excN+1-i)) then
                         excNtemp=sim%excN+1-i
                     endif
                  enddo
                  sim%excN=excNtemp
                  write(sim%outfile_23,884) sim%naesmd%tfemto,ihopavant,sim%naesmd%ihop,sim%excN,ipn,&
                                sim%naesmd%deltared*feVmdqt,kavant*feVmdqt
                  call flush(sim%outfile_23)
               endif
!###############################################3       
            endif
!###############################################3       
! end after josiah
! end added after Kirill
! end added for patch JAKB

!                sim%naesmd%conthop=1
!            end if

        endif

        if(crosstemp.eq.2.and.sim%naesmd%conthop.eq.0) then
            sim%naesmd%conthop2=1
            ! sim%naesmd%ihopprev keep the value of the previous sim%naesmd%state from where we hop
            ! in order to allow new hops for the next 2 steps only in case that
            ! we do not want to hop again to the same sim%naesmd%state
            sim%naesmd%ihopprev=sim%naesmd%ihop
            ini=0
            sim%naesmd%ihop=sim%naesmd%iordenhop(sim%naesmd%ihop)
            icheck=sim%naesmd%ihop
            ! after the hop, we reinitialize the variables
            call naesmd2qmmm_r(sim)
            call cpu_time(t_start)
            call deriv(sim,sim%naesmd%ihop)
            call cpu_time(t_finish)
            sim%time_deriv_took=sim%time_deriv_took+t_finish-t_start
            do j=1,sim%naesmd%natom
                xx(j)=sim%naesmd%rx(j)
                yy(j)=sim%naesmd%ry(j)
                zz(j)=sim%naesmd%rz(j)
            end do
            call do_sqm_davidson_update(sim,vmdqt=sim%naesmd%vmdqtnew,vgs=sim%naesmd%vgs)
            ytemp=yg(ihopavant)
            ytemp2=yg(ihopavant+sim%excN)
            yg(ihopavant)=yg(sim%naesmd%ihop)
            yg(ihopavant+sim%excN)=yg(sim%naesmd%ihop+sim%excN)
            yg(sim%naesmd%ihop)=ytemp
            yg(sim%naesmd%ihop+sim%excN)=ytemp2
	    call    setup(rkcomm, tend, yg, tmax, rk_tolerance, thresholds, &
	            'M','R')
        end if
        ! Evaluation of other crossings that do not involve the sim%naesmd%ihop sim%naesmd%state
        !*************************************************************
        ininonhop=1
        do i=1,sim%excN
            if(i.lt.sim%naesmd%iorden(i).and.i.ne.ihopavant.and.i.ne.sim%naesmd%iorden(ihopavant)) then
                if(cross(i).eq.2) then
                    ininonhop=0
                    icheck=i
                    ! after the hop, we reinicialize the variables
                    ytemp=yg(i)
                    ytemp2=yg(i+sim%excN)
                    yg(i)=yg(sim%naesmd%iorden(i))
                    yg(i+sim%excN)=yg(sim%naesmd%iorden(i)+sim%excN)
                    yg(sim%naesmd%iorden(i))=ytemp
                    yg(sim%naesmd%iorden(i)+sim%excN)=ytemp2
		    call    setup(rkcomm, tend, yg, tmax, rk_tolerance, thresholds, &
	            	'M','R')
                end if
            end if
        enddo
        ! Check the conservation of the total energy in the hop
        ! and recalculate the sim%naesmd%kin to be printed in writeoutput.f
        if(icheck.ne.0) then
            if(ini.eq.0.or.ininonhop.eq.0) then
                eapres=sim%naesmd%vmdqtnew(sim%naesmd%ihop)
                sim%naesmd%kin=0.0d0
                do j=1,sim%naesmd%natom
                    sim%naesmd%kin=sim%naesmd%kin+sim%naesmd%massmdqt(j)*(sim%naesmd%vx(j)**2+&
			sim%naesmd%vy(j)**2+sim%naesmd%vz(j)**2)/2
                end do
                eapres=eapres+sim%naesmd%kin
                do j=1,sim%excN
                    if(j.ne.sim%naesmd%ihop) then
                        if(j.lt.sim%naesmd%iorden(j)) then
                            if(cross(j).ne.0) then
                                write(sim%outfile_3,887) sim%naesmd%tfemto,cross(j),j,sim%naesmd%iorden(j),eavant,eapres
                                call flush(sim%outfile_3)
                            end if
                        end if
                    else
                        if(sim%naesmd%ihop.ne.ihopavant) then
                            write(sim%outfile_3,887) sim%naesmd%tfemto,cross(sim%naesmd%ihop), &
                                                     & ihopavant,sim%naesmd%ihop,eavant,eapres
                            call flush(sim%outfile_3)
                        end if
                    end if
                end do
            end if
        end if
        !***********************************
        ! end analyze the hopping
        !**********************************
887     FORMAT(F18.10,1X,I2,1X,I3,1X,I3,10000(1X,F18.10))
884   FORMAT(F18.10,1X,I2,1X,I3,1X,I3,1X,I3,10000(1X,F18.10))
        return
    end subroutine
    ! At the point of hop, in general, the value of the potential energy in the new
    ! surface is different to the one in the older.
    ! In order to conserve the energy, we adjust the velocities
    subroutine veladjustment(sim, lprint,Na,icheck,ini)
        implicit none
        type(simulation_t),pointer::sim
        integer Na,lprint
        integer i,j,icheck,ini,ihoptemp
        _REAL_ dij(Na*3),vicheck
        _REAL_ alpha,racine,ctehop1,dctehop1
        _REAL_ vtemp(sim%excN),vgstemp
        !********************************************************
        ! adjustment of velocities
        !********************************************************
        ! Added by ST: calculate here NAC <psi| d psi/dR> in one step:
        !   Current energy calculation
        call do_sqm_davidson_update(sim,vmdqt=vtemp,vgs=vgstemp)
        !    Feed here energies and wavefunctions and geometry, get back dij
        ! if necessary here the signs of the CI coefficient matrix can be checked right here
        ! analytical calculation of nacR
        call nacR_analytic_wrap(sim, sim%naesmd%ihop, icheck, dij)
        ! end of the calculation of the non-adiabatic coupling vector dij
        !*********************************
        if(lprint.ge.1) then
            j=1
            do i=1,sim%naesmd%natom
                write(sim%outfile_6,*) i,dij(j),dij(j+1),dij(j+2)
                j=j+3
            end do
            call flush(sim%outfile_6)
        end if
        ! calculation of the current energy
        ! and the velocities adjustment
        ihoptemp=icheck
        vicheck=vtemp(ihoptemp)
        alpha=sim%naesmd%vmdqtnew(sim%naesmd%ihop)-vicheck
        racine = 0.0d0

        j=1
        do i=1,sim%naesmd%natom
            racine=racine+sim%naesmd%vx(i)*dij(j)+sim%naesmd%vy(i)*dij(j+1) &
                +sim%naesmd%vz(i)*dij(j+2)
            j=j+3
        end do
        racine=racine**2
        j=1
        do i=1,sim%naesmd%natom
            racine=racine+2.0d0*alpha/sim%naesmd%massmdqt(i) &
                *(dij(j)**2+dij(j+1)**2+dij(j+2)**2)
            j=j+3
        end do
        if(racine.le.0.0d0) then
            ini=1
! change made after Kirill
            if(sim%naesmd%decorhop.eq.2) ini=2
! end change made after Kirill
            goto 4321
        end if
        ctehop1=0.0d0
        j=1
        do i=1,sim%naesmd%natom
            ctehop1=ctehop1+sim%naesmd%vx(i)*dij(j)+sim%naesmd%vy(i)*dij(j+1)+sim%naesmd%vz(i)*dij(j+2)
            j=j+3
        end do
        ctehop1=ctehop1+dsqrt(racine)
        dctehop1=0.d0
        j=1
        do i=1,sim%naesmd%natom
            dctehop1=dctehop1+1.0d0/sim%naesmd%massmdqt(i) &
                *(dij(j)**2+dij(j+1)**2+dij(j+2)**2)
            j=j+3
        end do
        ctehop1=ctehop1/dctehop1

        ! option to adjust the velocities in the direction of
        ! the nonadiabatic coupling vector
        j=1
        do i=1,sim%naesmd%natom
            sim%naesmd%vx(i)=sim%naesmd%vx(i)-ctehop1*dij(j)/sim%naesmd%massmdqt(i)
            sim%naesmd%vy(i)=sim%naesmd%vy(i)-ctehop1*dij(j+1)/sim%naesmd%massmdqt(i)
            sim%naesmd%vz(i)=sim%naesmd%vz(i)-ctehop1*dij(j+2)/sim%naesmd%massmdqt(i)
            j=j+3
        end do
4321 continue
     !********************************************************
     ! end of adjustment of velocities
     !********************************************************
     return
 end subroutine

subroutine decoherence_E0_and_C(sim)
        type(simulation_t),pointer::sim
            
            if(sim%constcoherE0.ne.0.d0.and.sim%constcoherC.ne.0.d0) then
                if(sim%naesmd%conthop.ne.1.and.sim%naesmd%conthop2.ne.1) then
                    !BTN: I have not had time to investigate this fully, but the simulation object gets garbled when passed to coherence. 
                    !My guess is that one of the variables passed (though it does not appear to be sim) is not defined identically correctly inside coherence.
                    !For now, just exit.
                    write(6,*) 'Your choices for decoher_e0 and decoher_c &
                        are not currently available'
                    stop
                    call coherence(sim,sim%rk_comm, sim%Na, sim%rk_comm%tend,sim%rk_comm%tmax, &
                                   sim%rk_comm%rk_tol, &
                         sim%rk_comm%thresholds,sim%naesmd%yg,sim%constcoherE0,sim%constcoherC,sim%cohertype)
                end if
            end if
end subroutine

subroutine indexx(n,arr,indx)
! ************************************************************************

      integer n,indx(n),M,NSTACK
      double precision arr(n)
      parameter (M=7, NSTACK=50)
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
! is in ascending order for j=1,2,...N. The input quantities n and arr are not changed.

      integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      double precision a

      do j=1,n
        indx(j)=j
      enddo
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do j=l+1,ir
           indxt=indx(j)
           a=arr(indxt)
           do i=j-1,l,-1
               if(arr(indx(i)).le.a)goto 2
               indx(i+1)=indx(i)
           enddo
           i=l-1
2          indx(i+1)=indxt
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
         itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
           i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
           j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
        istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1

      return
      end
 end module
