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
   use md_module
   use communism
   use qm2_davidson_module,only:qm2ds
   use Cosmo_C, only : solvent_model
   implicit none      

   type(simulation_t) sim
   integer  i,j,k,kk,slen,readstring,ibo,kki,lprint
   double precision ntot
   double precision xcm,ycm,zcm 
   character*1000 card
   character*1000 cardmopac

   double precision yg(2*sim%excN) 
   double precision energy !!JAB VE model

   if(solvent_model.eq.2) then
        call calc_excsolven(energy) !JAB Test
        vmdqt(ihop)=vmdqt(ihop)-0.5*energy/feVmdqt !JAB Test
   endif

   kin=0.0d0
   do i =1, natom
      kin=kin+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)/2
   end do

   ntot=0
   do j=1,sim%excN
      ntot=ntot+yg(j)**2
   end do

   kinini=kin
   if(state.eq.'fund') then
      vini=vgs
      etotini=kin+vgs
      write(98,889) tfemto,kin*feVmdqt,kin*feVmdqt-kinini*feVmdqt, &
         vgs*feVmdqt, &
         vgs*feVmdqt-vini*feVmdqt, &
         kin*feVmdqt+vgs*feVmdqt,kin*feVmdqt+vgs*feVmdqt-etotini*feVmdqt
   end if

   if(state.eq.'exct') then
      vini=vmdqt(ihop)
      etotini=kin+vmdqt(ihop)
      write(98,889) tfemto,kin*feVmdqt,kin*feVmdqt-kinini*feVmdqt, &
         vmdqt(ihop)*feVmdqt, &
         vmdqt(ihop)*feVmdqt-vini*feVmdqt, &
         kin*feVmdqt+vmdqt(ihop)*feVmdqt,kin*feVmdqt+vmdqt(ihop)*feVmdqt- &
         etotini*feVmdqt
   end if

   if(state.eq.'exct') then
      if(lprint.ge.1) then
         write(96,889) tfemto,vgs*feVmdqt, &
            (vmdqt(j)*feVmdqt,j=1,sim%excN)
         if(ibo.ne.1) then
            write(95,999) ihop,tfemto,(yg(j)**2,j=1,sim%excN),ntot
         end if
      end if

      if(lprint.ge.3.and.ibo.ne.1) then
         write(94,889) tfemto,(dsin(yg(j+sim%excN)),j=1,sim%excN)
         call flush(94)
      end if

      call flush(96)
      call flush(95)
   endif
! use for vibrations**************************************
!      write(99,779) tfemto,(atomtype(k),
!     $rx(k)*convl,ry(k)*convl,rz(k)*convl,k=1,natom)
!      write(80,779) tfemto,(atomtype(k),
!     $vx(k)*convl/convt,vy(k)*convl/convt,vz(k)*convl/convt,k=1,natom)
!***********************************************************

      if(lprint.ge.3) then
         write(85,889) tfemto,(sim%deriv_forces(1+3*(j-1)),j=1,natom)
         write(84,889) tfemto,(sim%deriv_forces(2+3*(j-1)),j=1,natom)
         write(83,889) tfemto,(sim%deriv_forces(3+3*(j-1)),j=1,natom)
        ! Check the position of the center of mass
         xcm=0.0d0
         ycm=0.0d0
         zcm=0.0d0

         do j=1,natom
            xcm=xcm+rx(j)*massmdqt(j)/masstot
            ycm=ycm+ry(j)*massmdqt(j)/masstot
            zcm=zcm+rz(j)*massmdqt(j)/masstot
         end do

         write(91,889) tfemto,xcm-xcmini,ycm-ycmini,zcm-zcmini
         call flush(85)
         call flush(84)
         call flush(83)
         call flush(91)
         call flush(99)
      end if

      if(state.eq.'exct'.and.lprint.ge.1) then
         write(89,889) tfemto,(qm2ds%v2(qm2ds%Nb*(j-1)+j,ihop),j=1,qm2ds%Nb)
         call flush(89)
! in order to print the initial transition density of all states
         do k=1,sim%excN
            write(77,889) tfemto,(qm2ds%v2(qm2ds%Nb*(j-1)+j,k),j=1,qm2ds%Nb)
            call flush(77)
         end do
      end if

! write the atomic positions in angstoms
      card='coordinates' // ktbig(icontini) // '.dat'
      OPEN(10,FILE=card)
      write(10,557) '$COORD'
      do i=1,txtinicoord
           write(10,555) txtinput(i)
      end do

      do i=1,natom
         write(10,999) atomtype(i),rx(i)*convl &
            ,ry(i)*convl,rz(i)*convl
      end do

      write(10,556) '$ENDCOORD'
      write(10,557) '$VELOC'

      do k=1,natom
         write(10,223) vx(k)*convl/convt, &
            vy(k)*convl/convt,vz(k)*convl/convt
      end do

      write(10,556) '$ENDVELOC'
      write(10,557) '$COEFF'

      do k=1,sim%excN
         write(10,223) yg(k)**2,dsin(yg(k+sim%excN))
      enddo

      write(10,556) '$ENDCOEFF'

      !do i=txtendcoord+7,jend
      !   write(10,555) txtinput(i)
      !end do
      !close(10)

      if(iview.eq.1) then

! to be used in case we want to print the transition densities of all the states at t=0
      do kki=1,sim%excN
         card='view' // ktbig(icontini) // '-' //  ktbig(kki) // '.DATA'
!************************************************************************************
!       card='view' // ktbig(icontini) // '-' //  ktbig(ihop) // '.DATA'
!************************************************************************************
         OPEN(90,FILE=card)
         write(90,440) ' Number of atoms:'
         write(90,99) natom
         write(90,441) ' Number of orbitals:'
         write(90,99) nao
         write(90,445) ' Number of occupied orbitals:'
         write(90,222) ' 1'
         write(90,442) ' Number of eigenvectors printed:'
         write(90,222) ' 1'
         write(90,441) ' Atomic coordinates:'

         do k=1,natom
            write(90,999) atomtype(k), &
               rx(k)*convl,ry(k)*convl,rz(k)*convl
         end do

         write(90,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
         do k=1,qm2ds%Nb
! to be used in case we want to print the transition densities of all the states at t=0
            write(90,*) qm2ds%v2(qm2ds%Nb*(k-1)+k,kki)
         end do
         close(90)
! to be used in case we want to print the transition densities of all the states at t=0
      end do
   end if

   call flush(98)

111   FORMAT(A2,14x,3(1x,F16.12,1x,I1))
222   FORMAT(A2,3(1x,F12.6))
223   FORMAT(3(1x,F16.10))
333   format(a5,I3,1X,a80)
440   format(a17)
441   format(a20)
442   format(a32)
443   format(a41)
444   format(a80)
445   format(a29)
889   FORMAT(10000(1X,F18.10))
887   FORMAT(F18.10,10000(1X,I1))
999   FORMAT(I3,1X,1000(1X,F18.10))
998   FORMAT(I4)
99    FORMAT(I5)
777   FORMAT(A7,I4,2X,A2,2X,A5,I4,4X,3(1X,F7.3))
778   format(1000000(f10.6,1x))
779   format(F18.10,10000(1X,I3,3(1X,F18.10)))
555   format(a90)
556   format(a9)
557   format(a6)
558   format(a7)
88    format(a1)

   return
   end subroutine

!***********************************************
! write output information 
! in the classical loop
!***********************************************

   subroutine writeoutput(sim,i,ibo,yg,lprint,cross)
   use naesmd_constants
   use naesmd_module
   use md_module
   use communism
   use qm2_davidson_module,only:qm2ds
   use Cosmo_C,only: solvent_model
   implicit none      
 
   type(simulation_t),pointer::sim
   INTEGER l,i,j,jj,k,kk,slen,readstring,Nb,ibo,kki,lprint
   double precision ntot
   double precision xcm,ycm,zcm 
   integer cross(sim%excN)
   character*1000 card
   character*1000 cardmopac
   double precision yg(2*sim%excN) 
   double precision poblacring1,poblacring4
   double precision energy !JAB Test
   integer nring,indx(sim%excN)
   logical first
   data first /.true./
   save first

   if(solvent_model.eq.2) then
        call calc_excsolven(energy) !JAB Test
        vmdqt(ihop)=vmdqt(ihop)-0.5*energy/feVmdqt !JAB Test
   endif

   if(state.eq.'fund') then
      write(98,889) tfemto,kin*feVmdqt,kin*feVmdqt-kinini*feVmdqt, &
         vgs*feVmdqt, &
         vgs*feVmdqt-vini*feVmdqt, &
         kin*feVmdqt+vgs*feVmdqt,kin*feVmdqt+vgs*feVmdqt-etotini*feVmdqt
   end if

   if(state.eq.'exct') then
      if(ibo.eq.1) then
         write(98,889) tfemto,kin*feVmdqt,kin*feVmdqt-kinini*feVmdqt, &
            vmdqt(ihop)*feVmdqt, &
            vmdqt(ihop)*feVmdqt-vini*feVmdqt, &
            kin*feVmdqt+vmdqt(ihop)*feVmdqt,kin*feVmdqt &
            +vmdqt(ihop)*feVmdqt &
            -etotini*feVmdqt
      else
         write(98,889) tfemto,kin*feVmdqt,kin*feVmdqt-kinini*feVmdqt, &
            vmdqtnew(ihop)*feVmdqt, &
            vmdqtnew(ihop)*feVmdqt-vini*feVmdqt, &
            kin*feVmdqt+vmdqtnew(ihop)*feVmdqt,kin*feVmdqt &
            +vmdqtnew(ihop)*feVmdqt &
            -etotini*feVmdqt

      end if
   end if

   !ntot is the variable to check the norm conservation
   if(state.eq.'exct') then
      ntot=0
      do j=1,sim%excN
         ntot=ntot+yg(j)**2
      end do

      if(lprint.ge.1) then
         if(ibo.eq.1) then
            write(6,*)tfemto,shape(tfemto)
            write(6,*)vgs,shape(vgs)
            write(6,*)feVmdqt,shape(feVmdqt)
            write(6,*)vmdqt,shape(vmdqt)
            write(6,*)sim%excN,shape(sim%excN)
            write(6,*) tfemto,vgs*feVmdqt, &
               (vmdqt(j)*feVmdqt,j=1,sim%excN)
            write(6,*)'done0'
            write(96,889) tfemto,vgs*feVmdqt, &
               (vmdqt(j)*feVmdqt,j=1,sim%excN)
            write(6,*)'done'
         else    
            write(96,889) tfemto,vgs*feVmdqt, &
               (vmdqt(j)*feVmdqt,j=1,sim%excN)
            write(95,999) ihop,tfemto,(yg(j)**2,j=1,sim%excN),ntot
            write(93,888) tfemto,((cadiabnew(j,k),k=1,sim%excN),j=1,sim%excN)

            call flush(95)
            call flush(93)
         end if

         call flush(96)
      end if

      if(lprint.ge.3.and.ibo.ne.1) then
         write(94,889) tfemto,(dsin(yg(j+sim%excN)),j=1,sim%excN)
         call flush(94)
      end if
   end if

   write(92,889) tfemto,tempi,tempf
   call flush(92)
!
! use for vibrations*****************************
!      write(99,779) tfemto,(atomtype(k),
!     $rx(k)*convl,ry(k)*convl,rz(k)*convl,k=1,natom)
!      write(80,779) tfemto,(atomtype(k),
!     $vx(k)*convl/convt,vy(k)*convl/convt,vz(k)*convl/convt,k=1,natom)
!******************************************************
!
   if(lprint.ge.3) then
      !write(85,889) tfemto,(fxmdqt(j)/convl*feVmdqt,j=1,natom)
      !write(84,889) tfemto,(fymdqt(j)/convl*feVmdqt,j=1,natom)
      !write(83,889) tfemto,(fzmdqt(j)/convl*feVmdqt,j=1,natom)
      write(85,889) tfemto,(sim%deriv_forces(1+3*(j-1)),j=1,natom)
      write(84,889) tfemto,(sim%deriv_forces(2+3*(j-1)),j=1,natom)
      write(83,889) tfemto,(sim%deriv_forces(3+3*(j-1)),j=1,natom)

      write(125,889) tfemto,(sim%naesmd%a%x(j),j=1,natom)
      write(126,889) tfemto,(ax(j),j=1,natom)
      ! Check the position of the center of mass
      xcm=0.0d0
      ycm=0.0d0
      zcm=0.0d0

      do j=1,natom
         xcm=xcm+rx(j)*massmdqt(j)/masstot
         ycm=ycm+ry(j)*massmdqt(j)/masstot
         zcm=zcm+rz(j)*massmdqt(j)/masstot
      end do

      write(91,889) tfemto,xcm-xcmini,ycm-ycmini,zcm-zcmini
      call flush(85)
      call flush(84)
      call flush(83)
      call flush(91)
      call flush(99)
      call flush(80)
   end if

   if(state.eq.'exct'.and.lprint.ge.2) then
      write(100,688) tfemto,(iorden(j),j=1,sim%excN),cross
      write(120,688) tfemto,(cross(j),j=1,sim%excN)
      call flush(120)
      call flush(100)
   end if

   if(state.eq.'exct'.and.lprint.ge.1) then
      write(89,889) tfemto,(qm2ds%v2(qm2ds%Nb*(j-1)+j,ihop),j=1,qm2ds%Nb)
      call flush(89)
   end if

   if(icont.ne.nstepcoord) then
      icont=icont+1
   else
      icont=1
      icontpdb=icontpdb+1

      card='coordinates' // ktbig(icontpdb) // '.dat'
      open(10,FILE=card)
      write(10,557) '$COORD'
      do k=1,txtinicoord
         write(10,555) txtinput(k)
      end do

      do k=1,natom
         write(10,999) atomtype(k),rx(k)*convl, &
         ry(k)*convl,rz(k)*convl
      end do
      write(10,556) '$ENDCOORD'
      write(10,557) '$VELOC'

      do k=1,natom
         write(10,223) vx(k)*convl/convt, &
            vy(k)*convl/convt,vz(k)*convl/convt
      end do

      write(10,556) '$ENDVELOC'
      write(10,557) '$COEFF'

      do k=1,sim%excN
         write(10,223) yg(k)**2,dsin(yg(k+sim%excN))
      end do

      write(10,556) '$ENDCOEFF'

      !do k=txtendcoord+8,txtendcoord+9
      !   write(10,555) txtinput(k)
      !end do

      !write(10,*)int(nstep-(icontpdb-icontini)*nstepcoord*nstepw) 
      !do k=txtendcoord+11,txtendcoord+21
      !   write(10,555) txtinput(k)
      !end do

      !write(10,*) iseedmdqt
      !do k=txtendcoord+23,jend
      !   write(10,555) txtinput(k)
      !end do
      !close(10)

      open (9,file='coords.xyz',access='append')
      write (9,*) natom
      write (9,449) 'FINAL HEAT OF FORMATION =   ', (kin+vgs)*feVmdqt, &
         '  time = ',tfemto
      do k=1,natom
         write(9,302) ELEMNT(atomtype(k)),rx(k)*convl, &
            ry(k)*convl,rz(k)*convl
      end do
      close (9)

      if(iview.eq.1) then
         do kki=1,sim%excN
         card='view' // ktbig(icontpdb) // '-' //  ktbig(kki) // '.DATA'
         OPEN(90,FILE=card)
         write(90,440) ' Number of atoms:'
         write(90,99) natom
         write(90,441) ' Number of orbitals:'
         write(90,99) nao
         write(90,445) ' Number of occupied orbitals:'
         write(90,222) ' 1'
         write(90,442) ' Number of eigenvectors printed:'
         write(90,222) ' 1'
         write(90,441) ' Atomic coordinates:'

         do k=1,natom
            write(90,999) atomtype(k), &
               rx(k)*convl,ry(k)*convl,rz(k)*convl
         end do

         write(90,443) ' Eigenvector:   1  with Eigenvalue:   0.0'
         do k=1,qm2ds%Nb
            write(90,*) qm2ds%v2(qm2ds%Nb*(k-1)+k,kki)
         end do
         close(90)
      end do
   end if
   endif
   call flush(98)

111   FORMAT(A2,14x,3(1x,F16.12,1x,I1))
222   FORMAT(A2,3(1x,F12.6))
223   FORMAT(3(1x,F16.10))
302   format(A3,3f12.6)
333   format(a5,I3,1X,a80)
440   format(a17)
441   format(a20)
442   format(a32)
443   format(a41)
444   format(a80)
445   format(a29)
449   format(a28,F18.10,a9,F18.10)
688   FORMAT(F18.10,10000(1X,I4))
889   FORMAT(20000(1X,F18.10))
888   FORMAT(30000(1X,F18.10))
887   FORMAT(F18.10,10000(1X,I1))
999   FORMAT(I3,1X,1000(1X,F18.10))
998   FORMAT(I4)
99    FORMAT(I5)
777   FORMAT(A7,I4,2X,A2,2X,A5,I4,4X,3(1X,F7.3))
778   format(1000000(f10.6,1x))
779   format(F18.10,10000(1X,I3,3(1X,F18.10)))
555   format(a90)
556   format(a9)
557   format(a6)
558   format(a7)
88    format(a1)

   return
   end subroutine
end module

