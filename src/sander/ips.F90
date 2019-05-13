! <compile=optimized>

module nbips

#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"
!
!      variables for Isotropic Periodic Sum (IPS) calculation
!
!      IPS         IPS options: 1--for both ele and l-j using 3D IPS
!                               2--for ele only using 3D IPS
!                               3--for l-j only using 3D IPS
!                               4--for both ele and l-j using 3D IPS/DFFT
!                               5--for ele only using 3D IPS/DFFT
!                               6--for l-j only using 3D IPS/DFFT
!      NNBIPS      Number of nonbonded atom pairs
!      NNBIPST     Provious Number of nonbonded atom pairs
!
integer,save:: IPS,NNBIPST,NNBIPS,IPSSIZ
!


!    IPS parameters
!    AIPSE(0:3)    Electrostatic IPS coefficients
!    BIPSE(3)    Electrostatic IPS derivative coefficients
!    AIPSVA(0:3)   Lennard-Jones repulsion IPS coefficients
!    BIPSVA(3)   Lennard-Jones repulsion IPS derivative coefficients
!    AIPSVC(0:3)   Lennard-Jones dispersion IPS coefficients
!    BIPSVC(3)   Lennard-Jones dispersion IPS derivative coefficients
!    RIPS      Radius of IPS local region 
!    RAIPS     Radius of IPS extensive local region 
!    PIPS*0    Self IPS pair energies 
!    PIPS*C    IPS boundary energies 
!    EIPSS*C   IPS boundary energies 
!    VIRIPS*C  IPS boundary virials 
!    VBOXIPS   IPS local region volume
!

_REAL_,save::  AIPSE(0:3),BIPSE(3), &
       AIPSVA(0:3),BIPSVA(3), &
       AIPSVC(0:3),BIPSVC(3), &
       RIPS,RIPS2,RIPSR,RIPS2R,RIPS6R,RIPS12R, &
       RAIPS,GRIDIPS,VBOXIPS,DVBIPS, &
       PIPSE0,PIPSVA0,PIPSVC0,PIPSEC,PIPSVAC,PIPSVCC,  &
       EIPSSNB,EIPSSEL,VIRIPS,EIPSANB,EIPSAEL,VIRAIPS

_REAL_,save::  CGSUM,CIJSUM,AIJSUM,VIREXIPS(3,3)


INTEGER,save::  MIPSX,MIPSY,MIPSZ,MIPSO

!     TEIPS   ! Perform IPS for electrostatic interaction
!     TVIPS   ! Perform IPS for Lennard-Jones interaction
!
LOGICAL,save:: TEIPS,TVIPS,TEAIPS,TVAIPS,TIPSMEM,TIPSUPD,TIPSFIN

_REAL_,allocatable,dimension(:),save :: wnb

integer,save :: ierr_allocate,alloc_err


_REAL_,allocatable,dimension(:),save :: elearray,vdwarray
_REAL_,allocatable,dimension(:),save :: elexx,elexy, &
       elexz,eleyy,eleyz,elezz,vdwxx,vdwxy,vdwxz,vdwyy,vdwyz,vdwzz

_REAL_,allocatable,dimension(:),save ::prefac1,prefac2,prefac3
!
integer,save :: forder
!
INTEGER,save :: NFFT1,NFFT2,NFFT3,nfftdim1,nfftdim2,nfftdim3, &
        nfftable,nffwork,SIZ_Q,SIZFFTAB,SIZFFWRK
_REAL_,allocatable,dimension(:),save :: fftable,ffwork

!===========================================================================
contains
!===========================================================================


!Performance updates by Ross Walker (TSRI, 2005)


      SUBROUTINE IPSSYS(NATOM,NTYPE,NTB,CG,CUT,CN1,CN2,IAC,ICO,X)
!-----------------------------------------------------------------------
!    Calculate system energy and force in 3D periodic systems
!
      use constants

      implicit none

#include "extra.h"

      INTEGER NTYPE,I3,J3
      INTEGER NATOM,NTB
      _REAL_ CG(*),CN1(*), CN2(*)
      _REAL_ X(*)
      INTEGER IAC(*),ICO(*)
      _REAL_  CUT
      INTEGER I,J
      INTEGER ITI,ITJ,ITMAX,IACI,IC,NITI,NITJ
      _REAL_  CISUMI
      _REAL_  XI,YI,ZI,XIJ,YIJ,ZIJ,R2
      _REAL_  AIJ,CIJ,ANIJ
      _REAL_  CGIJSUM

        INTEGER,allocatable,dimension(:) :: NSUMIT
        _REAL_,allocatable,dimension(:) :: CISUM
#ifdef MPI
   teips=.false.
   tvips=.false.
   teaips=.false.
   tvaips=.false.
   if((ips-4)*(ips-6) == 0 )tvaips =.true.
   if ( (ips-4)*(ips-5) == 0 )teaips =.true.
   if( tvaips.OR.( (ips -1)*(ips-3) == 0 ))tvips =.true.
   if( teaips.OR.((ips -1)*(ips-2) == 0 ))teips =.true.
#endif

!  IPS Radius:

      RIPS2=cut
      RIPS2R=one/cut
      RIPS6R=RIPS2R*RIPS2R*RIPS2R
      RIPS12R=RIPS6R*RIPS6R
      RIPS=SQRT(RIPS2)
      RIPSR=one/RIPS

!  AIPS parameters
      IPSSIZ=0

!  Ele IPS parameters:
          AIPSE(0)=-35.0/16.0
          AIPSE(1)=35.0/16.0
          AIPSE(2)=-21.0/16.0D0
          AIPSE(3)=5.0D0/16.0D0
!          PIPSEC=ZERO
        PIPSEC=ONE+AIPSE(0)+AIPSE(1)+AIPSE(2)+AIPSE(3)
        PIPSE0=AIPSE(0)-PIPSEC
        BIPSE(1)=TWO*AIPSE(1)
        BIPSE(2)=FOUR*AIPSE(2)
        BIPSE(3)=SIX*AIPSE(3)

!  Dispersion IPS parameters:
          AIPSVC(0)=7.0D0/16.0D0
          AIPSVC(1)=9.0D0/14.0D0
          AIPSVC(2)=-3.0D0/28.0D0
          AIPSVC(3)=6.0D0/7.0D0
        PIPSVCC=ONE+AIPSVC(0)+AIPSVC(1)+AIPSVC(2)+AIPSVC(3)
        PIPSVC0=AIPSVC(0)-PIPSVCC
        BIPSVC(1)=TWO*AIPSVC(1)
        BIPSVC(2)=FOUR*AIPSVC(2)
        BIPSVC(3)=SIX*AIPSVC(3)
!  Repulsion IPS parameters:
        AIPSVA(0)=5.0D0/787.0D0
        AIPSVA(1)=9.0D0/26.0D0
        AIPSVA(2)=-3.0D0/13.0D0
        AIPSVA(3)=27.0D0/26.0D0
        PIPSVAC=ONE+AIPSVA(0)+AIPSVA(1)+AIPSVA(2)+AIPSVA(3)
        PIPSVA0=AIPSVA(0)-PIPSVAC
        BIPSVA(1)=FOUR*AIPSVA(1)
        BIPSVA(2)=EIGHT*AIPSVA(2)
        BIPSVA(3)=TWO*SIX*AIPSVA(3)

!
      EIPSSNB=0.0D0
      EIPSSEL=0.0D0

!=======================================================================
!   Main loop begin
!=======================================================================
      if(allocated(NSUMIT))deallocate(NSUMIT,CISUM)
      allocate(NSUMIT(NTYPE),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate NSUMIT"
      allocate(CISUM(NTYPE),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate CISUM"

      DO ITI=1,NTYPE
        NSUMIT(ITI)=0
        CISUM(ITI)=0
      ENDDO
      CGSUM=0.0D0
      ITMAX=0
      DO I=1,NATOM
        CGSUM=CGSUM+CG(I)
        ITI=IAC(I)
        IF(ITI>NTYPE)STOP "problem in IAC!"
        NSUMIT(ITI)=NSUMIT(ITI)+1
      ENDDO
      CIJSUM=0.0D0
      AIJSUM=0.0D0
      NNBIPST=0

!  system energy is calculated based on all pairs:
      DO ITI=1,NTYPE
        iaci = ntype * (ITI - 1)
        ic = ico(iaci+ITI)
        AIJ=CN1(IC)
        CIJ=CN2(IC)
        NITI=NSUMIT(ITI)
        ANIJ=NITI*NITI*half
        CISUMI=CISUM(ITI)+CIJ*NITI
        CIJSUM=CIJSUM+CIJ*ANIJ
        AIJSUM=AIJSUM+AIJ*ANIJ
        NNBIPST=NNBIPST+NITI*NITI
        DO ITJ=ITI+1,NTYPE
          NITJ=NSUMIT(ITJ)
          ic = ico(iaci+ITJ)
          AIJ=CN1(IC)
          CIJ=CN2(IC)
          ANIJ=NITI*NITJ
          CIJSUM=CIJSUM+CIJ*ANIJ
          AIJSUM=AIJSUM+AIJ*ANIJ
          CISUMI=CISUMI+CIJ*NITJ
          CISUM(ITJ)=CISUM(ITJ)+CIJ*NITI
          NNBIPST=NNBIPST+2*NITI*NITJ
        ENDDO
        CISUM(ITI)=CISUMI
      ENDDO
      CGIJSUM=CGSUM*CGSUM*half
      if( cgijsum < 1.D-15 ) cgijsum = 0.d0
      IF(TVIPS.and.master)THEN
        EIPSSNB=AIJSUM*PIPSVAC*RIPS12R
        IF(.NOT.TVAIPS)EIPSSNB=EIPSSNB-CIJSUM*PIPSVCC*RIPS6R
      ELSE
        EIPSSNB=0.0D0
      ENDIF
      IF(TEIPS.and.master)THEN
        EIPSSEL=0.0D0
        IF(.NOT.TEAIPS)EIPSSEL=CGIJSUM*PIPSEC*RIPSR
      ELSE
        EIPSSEL=0.0D0
      ENDIF
      IF(TVAIPS.OR.TEAIPS)THEN
      if(allocated(wnb))deallocate(wnb)
      allocate(wnb(natom),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate wnb"
        DO I=1,NATOM
          ITI=IAC(I)
!   pair interaction part
          WNB(I)=CISUM(ITI)/SQRT(CIJSUM*TWO)
        ENDDO
      ENDIF
      deallocate(NSUMIT,CISUM)
 
!=======================================================================
!   Main loop end
!=======================================================================

! Calculate volume virials:

      VIRIPS=-(EIPSSNB+EIPSSEL)
      VBOXIPS=FOURPI*RIPS2*RIPS*third
      tipsfin=.false.
      IF(NTB.EQ.0)THEN
        tipsfin=.true.
        ! For non-periodic system, system energy is calculated based on cutoff
        NNBIPS=0
        DO I=1,NATOM
          ! Atom i long-range reference and self-interaction:
          I3=3*I-2
          XI=X(I3)
          YI=X(I3+1)
          ZI=X(I3+2)
          NNBIPS=NNBIPS+1
          DO J=I+1,NATOM
            J3=3*J-2
            XIJ=XI-X(J3)
            YIJ=YI-X(J3+1)
            ZIJ=ZI-X(J3+2)
            R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            !  RIPS2 is infinite except for non periodic conditions
            IF(R2<RIPS2)NNBIPS=NNBIPS+2
          ENDDO
        ENDDO
      ENDIF

      if( master ) then
         WRITE(6,'(" ----------------------------------")')
         IF(TEAIPS.OR.TVAIPS)THEN
            WRITE(6,'(" Using 3D-IPS/DFFT algorithm")')
            IF(TIPSFIN)THEN
                WRITE(6,'("   An infinite IPS Radius is used")')
            ELSE IF(RAIPS>0)THEN
                WRITE(6,'("   IPS Radius: ",F6.2," A")')RAIPS
            ELSE
                WRITE(6,'("   IPS Radius is set to the box size")')
            ENDIF
            IF(MIPSX>0)THEN
                WRITE(6,'("   FFT dimensions:",3I6)')MIPSX,MIPSY,MIPSZ
            ELSE
                WRITE(6,'("   FFT dimensions are set based on grid size:",F6.2," A")')GRIDIPS
            ENDIF
            WRITE(6,'("   FFT b-spline order is:",I4)')MIPSO
         ELSE
            WRITE(6,'(" Using 3D-IPS algorithm")')
            WRITE(6,'("   IPS Radius: ",F6.2," A")')RIPS
         ENDIF
         IF(TEIPS)THEN
            WRITE(6,'("   Using IPS for electrostatic energy")')
         ENDIF
         IF(TVIPS)THEN
            WRITE(6,'("   Using IPS for L-J energy")')
         ENDIF
         WRITE(6,'(" ----------------------------------")')
      else 
        EIPSSNB=0.0D0
        EIPSSEL=0.0D0
        VIRIPS=0.0D0
        NNBIPS=0
      end if
      RETURN 

end subroutine ipssys

SUBROUTINE IPSUPDATE(NTB)

!-----------------------------------------------------------------------
!     Update parameters once IPS radius or the box size changed
!-----------------------------------------------------------------------

      use nblist, only : volume
      implicit none
#ifdef MPI
#  include "parallel.h"
   include "mpif.h"
   integer ierr
   INTEGER NNBTMP
#endif
      INTEGER NTB
      _REAL_ FIPS,CHANGE 

      IF(NTB>0)THEN
        FIPS=VBOXIPS/volume
      ELSE
#ifdef MPI
        call mpi_allreduce(NNBIPS,NNBTMP,1,MPI_INTEGER, &
         mpi_sum,commsander,ierr)
        NNBIPS=NNBTMP 
#endif
        FIPS=NNBIPS*1.0D0/NNBIPST
      ENDIF
      CHANGE=ABS(FIPS-1.0D0)

      ! Update system energies and forces:
      IF(CHANGE>DVBIPS)THEN
         VIRIPS=VIRIPS*FIPS
         EIPSSNB=EIPSSNB*FIPS
         EIPSSEL=EIPSSEL*FIPS
         TIPSUPD=NTB>0
         VBOXIPS=volume
         NNBIPST=NNBIPS
      ELSE
         TIPSUPD=.FALSE.
      ENDIF
      RETURN
end subroutine ipsupdate

SUBROUTINE EEXIPS(ENB,EEL,IFRSTA,ILASTA,NTB,NTYPES,IAC,ICO,numex,natex,CG, &
                        CN1,CN2,DX,X)

!-----------------------------------------------------------------------
!   3D IPS interaction between excluded atom pairs
!   This routine must be called first to update IPS parameters when needed
!       and to initalize electrostatic and vdw energies
!
!   by Xiongwu Wu  - 9/10/2004
!
!-----------------------------------------------------------------------

      use constants, only : zero,half,ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE
!      use crg_reloc, only : ifcr, cr_add_dcdr_factor

      implicit none
      _REAL_ ENB, EEL
      INTEGER IFRSTA,ILASTA, NTB,NTYPES,IAC(*),ICO(*),numex(*),natex(*)
      _REAL_ CG(*),CN1(*), CN2(*)
      _REAL_ DX(*),X(*)
      ! local:
      INTEGER I,J,K,IC,IACI,ITI,ITJ
      INTEGER I3,I31,I32,J3,J31,J32
      INTEGER  ILAST, IFIRST
      _REAL_ DXI, DYI, DZI,DIJ,DIJX,DIJY,DIJZ
      _REAL_ ENBIJ,EELIJ
      _REAL_  XI,YI,ZI,XIJ,YIJ,ZIJ
      _REAL_  CGI,CGIJ,AIJ,CIJ
      _REAL_  R2,R6,R12,U1,U2,U4
      _REAL_  PE,DEU,PVA,DVAU,PVC,DVCU

      ! check to see if volume or atom pair changed 
      CALL IPSUPDATE(NTB)

      ! Setup constants for use in inner loops
      ENB=EIPSSNB
      EEL=EIPSSEL
      DEU=zero
      DVAU=zero
      DVCU=zero
      CGIJ=zero
      AIJ=zero
      CIJ=zero
      VIREXIPS = zero
      VIREXIPS(1,1) = virips
      VIREXIPS(2,2) = virips
      VIREXIPS(3,3) = virips

!=======================================================================
!   Main loop begin
!=======================================================================

      ILAST = 0
      DO I=1,IFRSTA-1
         ILAST=ILAST+numex(I)
      ENDDO
      NNBIPS=0
      DO I=IFRSTA,ILASTA !1,natom
         NNBIPS=NNBIPS+1
         IFIRST = ILAST+1
         ILAST=ILAST+numex(I)
         IF(TEIPS)THEN
            CGI=CG(I)
            CGIJ=CGI*CGI
            EELIJ=HALF*CGIJ*PIPSE0*RIPSR
            EEL=EEL+EELIJ
         ENDIF
         IF(TVIPS)THEN
            ITI=IAC(I)
            iaci = ntypes * (iti - 1)
            ic = ico(iaci+ITI)
            AIJ=CN1(IC)
            CIJ=CN2(IC)
            ! Atom i long-range reference and self-interaction
            ENBIJ=HALF*(AIJ*PIPSVA0*RIPS6R-CIJ*PIPSVC0)*RIPS6R
            ENB=ENB+ENBIJ
         ENDIF
         I32=3*I
         I31=I32-1
         I3=I32-2
         DXI=DX(I3)
         DYI=DX(I31)
         DZI=DX(I32)
         XI=X(I3)
         YI=X(I31)
         ZI=X(I32)
         DO K=IFIRST,ILAST !1,No interactions for this atom.
            J = natex(K) !negative if atom J should be excluded from interaction
                         ! with I.  This deals with 1-2,1-3 and 1-4 only. 
                        ! Does not deal with cutoff.
           IF(J<=0) cycle
           J32=3*J
           J31=J32-1
           J3=J32-2
           XIJ=X(J3)-XI
           YIJ=X(J31)-YI
           ZIJ=X(J32)-ZI
           R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
           NNBIPS=NNBIPS+2
           R6=R2*R2*R2
           R12=R6*R6
           U2=R2*RIPS2R
           IF(TEIPS)THEN
            U1=SQRT(U2)
            CGIJ=CGI*CG(J)/RIPS
            PE=AIPSE(0)+U2*(AIPSE(1)+U2*(AIPSE(2)+U2*(AIPSE(3))))
            DEU=U2*(BIPSE(1)+U2*(BIPSE(2)+U2*(BIPSE(3))))
            EELIJ=CGIJ*(PE-PIPSEC)
            EEL=EEL+EELIJ
            DIJ=-CGIJ*DEU/R2
           ELSE
            DIJ=ZERO
           ENDIF
           IF(TVIPS)THEN
              ITJ=IAC(J)
              ic = ico(iaci+ITJ)
              AIJ=CN1(IC)*RIPS12R
              CIJ=CN2(IC)*RIPS6R
              U4=U2*U2

!
            PVC=CIJ*(AIPSVC(0)+U2*(AIPSVC(1)+U2*(AIPSVC(2)+U2*(AIPSVC(3))))-PIPSVCC)
            DVCU=CIJ*(U2*(BIPSVC(1)+U2*(BIPSVC(2)+U2*BIPSVC(3))))
!
            PVA=AIJ*(AIPSVA(0)+U4*(AIPSVA(1)+U4*(AIPSVA(2) +U4*(AIPSVA(3))))-PIPSVAC)
            DVAU=AIJ*(U4*(BIPSVA(1)+U4*(BIPSVA(2)+U4*BIPSVA(3))))
            ENBIJ=PVA-PVC
            ENB=ENB+ENBIJ
            DIJ=DIJ-(DVAU-DVCU)/R2
           ENDIF

           DIJX=DIJ*XIJ
           DIJY=DIJ*YIJ
           DIJZ=DIJ*ZIJ
           DXI=DXI-DIJX
           DYI=DYI-DIJY
           DZI=DZI-DIJZ
           DX(J3)=DX(J3)+DIJX
           DX(J31)=DX(J31)+DIJY
           DX(J32)=DX(J32)+DIJZ
           VIREXIPS(1,1) = VIREXIPS(1,1) - DIJX*XIJ
           VIREXIPS(1,2) = VIREXIPS(1,2) - DIJX*YIJ
           VIREXIPS(1,3) = VIREXIPS(1,3) - DIJX*ZIJ
           VIREXIPS(2,2) = VIREXIPS(2,2) - DIJY*YIJ
           VIREXIPS(2,3) = VIREXIPS(2,3) - DIJY*ZIJ
           VIREXIPS(3,3) = VIREXIPS(3,3) - DIJZ*ZIJ
         end do
         DX(I3) =  DXI
         DX(I31) =  DYI
         DX(I32) =  DZI
      ENDDO

!=======================================================================
!   Main loop end
!=======================================================================
      VIREXIPS(2,1)=VIREXIPS(1,2)
      VIREXIPS(3,1)=VIREXIPS(1,3)
      VIREXIPS(3,2)=VIREXIPS(2,3)
      RETURN
end subroutine eexips


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine aipspbc here]
subroutine aipspbc( eerw,eerq,natom,crd,charge,frcx,frc,rec_vir)
   use stack
   use ew_bspline,only:get_grid_weights
   use nblist,only:recip
   
   implicit none
#  include "flocntrl.h"
      
#  include "def_time.h"
#  include "box.h"

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
    include "mpif.h"
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif


   ! INPUT
   !       natom:  number of atoms
   !       crd   atomic coords
   !       charge  atomic charges

#ifdef MPI
   integer ierr
#endif
   
   integer natom,num_ks_trial
   _REAL_ crd(3,natom),charge(natom)
   integer nmine,nderiv
   
   ! OUTPUT
   !       eerq:  ips remaining electrostatic  energy
   !       eerw:  ips remaining Lennard-Jones  energy
   !       frc forces incremented by the remaining sum
   
   _REAL_ eerq,eerw,eeaips,evaips
   _REAL_ frcx(3),frc(3,natom)
   _REAL_ rec_vir(3,3),aips_vir(3,3)
   
   character(kind=1,len=7) :: routine="do_aips"
   ! HEAP STORAGE:  These arrays need to be preserved throughout simulation
   
   integer l_d2th1,l_d2th2,l_d2th3
   
   integer l_fr1,l_fr2,l_fr3
   integer l_th1,l_th2,l_th3,l_dth1,l_dth2,l_dth3
   integer l_fftw,l_q,l_w
   integer imy_cg
   integer num_ks
   integer m1,m2
   save num_ks
   data num_ks/0/
   integer l_tmpy,l_alpha,l_beta
   logical not_done

   nderiv = 1
   
   call aipscell(natom,crd)
!
!     get some integer array dimensions
   call get_fftdims(nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,SIZFFTAB,SIZFFWRK)
# ifdef MPI
#   ifdef LES
   num_ks = natom
#   else
   if(num_ks == 0) &
         num_ks = min((((nxyslab(0)+forder-1)*natom*4)/ &
         (3*nfft3)),natom)
#   endif
# else
   num_ks = natom
# endif
      IF(TIPSMEM.OR.TIPSUPD)THEN
!    Calculate IPS FUNCTIONS
        CALL AIPS_SETUP(NATOM)
      ENDIF

   call timer_start(TIME_AIPS_GRID)

   num_ks_trial = 0
   not_done = .true.
   do while(not_done)
      call get_stack(l_fftw,sizffwrk,routine)
      call get_stack(l_q,siz_q,routine)
      call get_stack(l_w,siz_q,routine)
      call get_stack(l_fr1,num_ks,routine)
      call get_stack(l_fr2,num_ks,routine)
      call get_stack(l_fr3,num_ks,routine)
      call get_stack(l_th1,num_ks*forder,routine)
      call get_stack(l_th2,num_ks*forder,routine)
      call get_stack(l_th3,num_ks*forder,routine)
      if ( nderiv == 1 )then
         call get_stack(l_dth1,num_ks*forder,routine)
         call get_stack(l_dth2,num_ks*forder,routine)
         call get_stack(l_dth3,num_ks*forder,routine)
         call get_stack(l_d2th1,forder,routine)
         call get_stack(l_d2th2,forder,routine)
         call get_stack(l_d2th3,forder,routine)
      end if
      if(.not. rstack_ok)then
         deallocate(r_stack)
         allocate(r_stack(1:lastrst),stat=alloc_ier)
         call reassign_rstack(routine)
      endif
      REQUIRE(rstack_ok)

      call get_istack(imy_cg,num_ks,routine)
      if(.not. istack_ok)then
         deallocate(i_stack)
         allocate(i_stack(1:lastist),stat=alloc_ier)
         call reassign_istack(routine)
      endif
      REQUIRE(istack_ok)

      call get_grid_weights( &
            natom,crd,recip,nfft1,nfft2,nfft3, &
            r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3),forder, &
            r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
            r_stack(l_dth1),r_stack(l_dth2),r_stack(l_dth3), &
            r_stack(l_d2th1),r_stack(l_d2th2),r_stack(l_d2th3), &
            i_stack(imy_cg),nmine,nderiv,num_ks)
      if(nmine > num_ks)then
         if( num_ks_trial >= 2 ) call mexit( 6,1 )
         if ( nderiv == 1 )then
            call free_stack(l_d2th3,routine)
            call free_stack(l_d2th2,routine)
            call free_stack(l_d2th1,routine)
            call free_stack(l_dth3,routine)
            call free_stack(l_dth2,routine)
            call free_stack(l_dth1,routine)
         end if
         call free_stack(l_th3,routine)
         call free_stack(l_th2,routine)
         call free_stack(l_th1,routine)
         call free_stack(l_fr3,routine)
         call free_stack(l_fr2,routine)
         call free_stack(l_fr1,routine)
         call free_stack(l_w,routine)
         call free_stack(l_q,routine)
         call free_stack(l_fftw,routine)
         call free_istack(imy_cg,routine)
         num_ks= nmine*4/3
         num_ks_trial = num_ks_trial + 1
      else
         not_done = .false.
      end if
   enddo
   
   !........Fill Charge Grid
   !          charges are approximated on an even grid
#ifdef MPI
   call fill_ips_grid(natom,charge,wnb, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         r_stack(l_q),r_stack(l_w),i_stack(imy_cg),nmine)
#else
   call fill_ips_grid(natom,charge,wnb, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         r_stack(l_q),r_stack(l_w),nmine)
#endif

   call timer_stop_start(TIME_AIPS_GRID,TIME_AIPS_FFT)
   call get_stack(l_tmpy,2*nfftdim1,routine)
   call get_stack(l_alpha,nfft1,routine)
   call get_stack(l_beta,nfft1,routine)
   if(.not. rstack_ok)then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
   endif
   REQUIRE(rstack_ok)

      IF(teaips)THEN
        call fft_backrc( &
         r_stack(l_q),fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ELSE
        r_stack(l_q:l_q+SIZ_Q-1)=0.0d0
      ENDIF
      IF(tvaips)THEN
           call fft_backrc( &
         r_stack(l_w),fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ELSE
        r_stack(l_w:l_w+SIZ_Q-1)=0.0d0
      ENDIF


   call timer_stop_start(TIME_AIPS_FFT,TIME_AIPS_SUM)
   
   !           -------------SCALAR SUM------------------
   
      call aips_sumrc( &
            r_stack(l_q), r_stack(l_w),elearray,vdwarray,elexx,elexy, &
       elexz,eleyy,eleyz,elezz,vdwxx,vdwxy,vdwxz,vdwyy,vdwyz,vdwzz,&
       nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eeaips,evaips,aips_vir)
#ifndef noVIRIAL
   do m2 = 1,3
      do m1 = 1,3
         rec_vir(m1,m2) = rec_vir(m1,m2)+aips_vir(m1,m2)
      end do
   end do
#endif

   call timer_stop_start(TIME_AIPS_SUM,TIME_AIPS_FFT)
   
   !           -----------FFT FORWARD--------------------
   
      IF(teaips)THEN
        eerq=eerq+eeaips
   call fft_forwardrc( &
         r_stack(l_q),fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ENDIF
      IF(tvaips)THEN
        eerw=eerw+evaips
   call fft_forwardrc( &
         r_stack(l_w),fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ENDIF
   call free_stack(l_beta,routine)
   call free_stack(l_alpha,routine)
   call free_stack(l_tmpy,routine)

   call timer_stop_start(TIME_AIPS_FFT,TIME_AIPS_FRC)
   
   !           -----------GRAD SUM--------------------
#ifdef MPI
   call aips_grad_sumrc( &
         natom,charge,wnb,recip, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_dth1),r_stack(l_dth2),r_stack(l_dth3), &
         frcx,frc,r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         r_stack(l_q), r_stack(l_w),&
         i_stack(imy_cg),nmine)
#else
   call aips_grad_sumrc( &
         natom,charge,wnb,recip, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_dth1),r_stack(l_dth2),r_stack(l_dth3), &
         frcx,frc,r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         r_stack(l_q), r_stack(l_w),&
         nmine)
#endif
         
   call timer_stop(TIME_AIPS_FRC)
   if ( nderiv == 1 )then
      call free_stack(l_d2th3,routine)
      call free_stack(l_d2th2,routine)
      call free_stack(l_d2th1,routine)
      call free_stack(l_dth3,routine)
      call free_stack(l_dth2,routine)
      call free_stack(l_dth1,routine)
   end if
   call free_stack(l_th3,routine)
   call free_stack(l_th2,routine)
   call free_stack(l_th1,routine)
   call free_stack(l_fr3,routine)
   call free_stack(l_fr2,routine)
   call free_stack(l_fr1,routine)
   call free_stack(l_w,routine)
   call free_stack(l_q,routine)
   call free_stack(l_fftw,routine)

   call free_istack(imy_cg,routine)
#ifdef MPI
   call mpi_barrier(recip_comm,ierr)
#endif
   return
end subroutine aipspbc


!****************************************************************
!                        AIPS_SETUP
!****************************************************************
      SUBROUTINE AIPS_SETUP(NATOM)
!-----------------------------------------------------------------------
!     This routine allocates space and defines variables for
!     the FFT calculation
!-----------------------------------------------------------------------
!
      use ew_bspline,only:load_prefacs
      use stack
#ifdef MPI
  use fft,only:fft_init,column_fft_flag
#endif
     implicit none
     INTEGER NATOM
!
     character(kind=1,len=10) :: routine="aips_setup"
      _REAL_ dummy
      INTEGER kbot,ktop
      integer l_fftw,l_tmpy,l_alpha,l_beta,opt_infl
      integer siztheta,sizstack

      INTEGER MAXORDER,MAXNFFT
      PARAMETER (MAXORDER=25, MAXNFFT=2000)
#  include "def_time.h"
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
      INTEGER kbot0
#endif

#ifndef MPI
   siztheta = natom*forder
   siz_q = 2*nfftdim1*nfftdim2*nfftdim3
#else
   siztheta = natom*forder*(nxyslab(0)+6)/nfft3
   siz_q = max( ntxyslab*(nxyslab(0)), ntxzslab*nxzslab(0))
#endif
   sizstack = siz_q+6*siztheta+sizffwrk+3*natom
!
      if(allocated(prefac1)) deallocate(prefac1,prefac2,prefac3)
      allocate(prefac1(nfft1),prefac2(nfft2),prefac3(nfft3))

     if(allocated(fftable))deallocate(fftable,ffwork)
      allocate(fftable(sizfftab),ffwork(sizffwrk))

   opt_infl=1
   call load_prefacs(prefac1,prefac2,prefac3, &
         nfft1,nfft2,nfft3,forder,opt_infl)
   call fft_setup(dummy,fftable,ffwork, &
         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
         nfftable,nffwork)
#ifdef MPI
   if(column_fft_flag)then
      call fft_init(nfft1,nfft2,nfft3)
   endif
#endif
!     IPS energy functions
      IF(IPSSIZ>0)THEN
          deallocate(elearray,vdwarray,stat=alloc_err)
          if(alloc_err /= 0 ) write(0,*)"unable to deallocate elearray"
          deallocate(elexx,elexy,elexz,eleyy,eleyz,elezz,stat=alloc_err)
          if(alloc_err /= 0 ) write(0,*)"unable to deallocate elexx"
          deallocate(vdwxx,vdwxy,vdwxz,vdwyy,vdwyz,vdwzz,stat=alloc_err)
          if(alloc_err /= 0 ) write(0,*)"unable to deallocate vdwxx"
      ENDIF
      allocate(elearray(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate elearray"
      allocate(elexx(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate elexx"
      allocate(elexy(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate elexy"
      allocate(elexz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate elexz"
      allocate(eleyy(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate eleyy"
      allocate(eleyz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate eleyz"
      allocate(elezz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate elezz"
      allocate(vdwarray(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwarray"
      allocate(vdwxx(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwxx"
      allocate(vdwxy(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwxy"
      allocate(vdwxz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwxz"
      allocate(vdwyy(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwyy"
      allocate(vdwyz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwyz"
      allocate(vdwzz(siz_q),stat=alloc_err)
      if(alloc_err /= 0 ) write(0,*)"unable to allocate vdwzz"
!  end of memory allocation
      IPSSIZ=SIZ_Q
! Update energy function grid. stored at the image part of Q and W
        elearray(1:SIZ_Q)=0.0d0
        elexx(1:SIZ_Q)=0.0d0
        elexy(1:SIZ_Q)=0.0d0
        elexz(1:SIZ_Q)=0.0d0
        eleyy(1:SIZ_Q)=0.0d0
        eleyz(1:SIZ_Q)=0.0d0
        elezz(1:SIZ_Q)=0.0d0
        vdwarray(1:SIZ_Q)=0.0d0
        vdwxx(1:SIZ_Q)=0.0d0
        vdwxy(1:SIZ_Q)=0.0d0
        vdwxz(1:SIZ_Q)=0.0d0
        vdwyy(1:SIZ_Q)=0.0d0
        vdwyz(1:SIZ_Q)=0.0d0
        vdwzz(1:SIZ_Q)=0.0d0
#ifdef MPI
   kbot0 = mxystart(mytaskid)
   kbot = kbot0 + 1
   ktop = kbot0 + mxyslabs
#else
      kbot = 1
      ktop = nfft3
#endif
  call timer_start(TIME_AIPS_FUNC)
        CALL PBC_IPS_ENG(kbot, ktop)
   call timer_stop_start(TIME_AIPS_FUNC,TIME_AIPS_FFT)
   call get_stack(l_fftw,sizffwrk,routine)
   call get_stack(l_tmpy,2*nfftdim1,routine)
   call get_stack(l_alpha,nfft1,routine)
   call get_stack(l_beta,nfft1,routine)
      if(.not. rstack_ok)then
         deallocate(r_stack)
         allocate(r_stack(1:lastrst),stat=alloc_ier)
         call reassign_rstack(routine)
      endif
      REQUIRE(rstack_ok)
      IF(TEAIPS)THEN
           call fft_backrc( &
         elearray,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         elexx,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         elexy,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         elexz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         eleyy,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         eleyz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         elezz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ENDIF
      IF(TVAIPS)THEN
           call fft_backrc( &
         vdwarray,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwxx,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwxy,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwxz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwyy,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwyz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
           call fft_backrc( &
         vdwzz,fftable,r_stack(l_fftw), &
         nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork &
         ,r_stack(l_tmpy), &
         r_stack(l_alpha),r_stack(l_beta) &
         )
      ENDIF
  call timer_stop(TIME_AIPS_FFT)
   call free_stack(l_beta,routine)
   call free_stack(l_alpha,routine)
   call free_stack(l_tmpy,routine)
   call free_stack(l_fftw,routine)
      RETURN
      END subroutine aips_setup


      SUBROUTINE PBC_IPS_ENG(KBOT,KTOP)    
! Calculate grid IPS potentials

      use constants
      use nblist, only: ucell,volume
      implicit none
#  include "extra.h"
      INTEGER KBOT,KTOP,MRCX,MRCY,MRCZ

      INTEGER IBINX,IBINY,IBINZ,IBINX0,IBINY0,IBINZ0
      INTEGER IBINX1,IBINY1,IBINZ1
      INTEGER I000
      _REAL_ RIPSC,RIPSC2,RIPSCR,RIPSC2R,RIPSC6R,FIPS
      _REAL_ WRK11,WRK21,WRK31,WRK12,WRK22,WRK32
      _REAL_ WRK13,WRK23,WRK33
      _REAL_ XBIN,YBIN,ZBIN
      _REAL_ EELIJ,ENBIJ
      _REAL_  XI,YI,ZI,XIJ,YIJ,ZIJ,PC(3)
      _REAL_  R2,R2R
      _REAL_ U1,U2,U4,U8,U6R,U12R
      _REAL_ PE,DPE,DEIJ,PEC,PDEC
      _REAL_ PVC,DPVC,DVIJ,PVCC,PDVCC
      _REAL_ CIPSE(0:3),CIPSDE(0:3),DIPSE(3),DIPSDE(3)
      _REAL_ CIPSVC(0:3),CIPSDVC(0:3),DIPSVC(3),DIPSDVC(3)
! RCFFT addition
      integer nfftdimrc
      nfftdimrc=2*nfftdim1
! Convert to fractional coordinates
      WRK11=ucell(1,1)
      WRK21=ucell(2,1)
      WRK31=ucell(3,1)
      WRK12=ucell(1,2)
      WRK22=ucell(2,2)
      WRK32=ucell(3,2)
      WRK13=ucell(1,3)
      WRK23=ucell(2,3)
      WRK33=ucell(3,3)
      XBIN=1.0d0/NFFT1
      YBIN=1.0d0/NFFT2
      ZBIN=1.0d0/NFFT3
!  Calculate anisotropic radius
      IF(RAIPS>ZERO)THEN
        RIPSC=RAIPS
      ELSE
!  Define RIPSC the longest box side
        RIPSC=ABS(wrk11)
        IF(RIPSC<ABS(wrk22))RIPSC=ABS(wrk22)
        IF(RIPSC<ABS(wrk33))RIPSC=ABS(wrk33)
!   Make RIPSC constant within certain range of box sizes
        RIPSC=ANINT(RIPSC)
        IF(RIPSC<RIPS+20)RIPSC=RIPS+20
      ENDIF
      RIPSCR=ONE/RIPSC
      RIPSC2=RIPSC*RIPSC
      RIPSC2R=ONE/RIPSC2
      RIPSC6R=RIPSC2R*RIPSC2R*RIPSC2R
      IF(TIPSFIN)THEN
        MRCX=NFFT1/2-1
        MRCY=NFFT2/2-1
        MRCZ=NFFT3/2-1
!  electrostatic parameters
        CIPSE(0)=ZERO
        CIPSE(1)=ZERO
        CIPSE(2)=ZERO
        CIPSE(3)=ZERO
        PEC=ZERO
!   r6 parameters
        CIPSVC(0)=ZERO
        CIPSVC(1)=ZERO
        CIPSVC(2)=ZERO
        CIPSVC(3)=ZERO
        PVCC=ZERO
      ELSE
! Find the lattice range that enclose the cutoff sphere
!    x-y plan is: pxy=cross(vx,vy); 
!        plan distance is: rz=dot(vz,pxy);
!        plan numbers: 2*ripsc/rz + 1
        CALL CROSS(ucell(1,1),ucell(1,2),PC)
        PC=PC/SQRT(PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3))
        CALL DOT(PC,ucell(1,3),ZIJ)
        MRCZ=INT(RIPSC*NFFT3/ABS(ZIJ)+ONE)
        CALL CROSS(ucell(1,2),ucell(1,3),PC)
        PC=PC/SQRT(PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3))
        CALL DOT(PC,ucell(1,1),XIJ)
        MRCX=INT(RIPSC*NFFT1/ABS(XIJ)+ONE)
        CALL CROSS(ucell(1,3),ucell(1,1),PC)
        PC=PC/SQRT(PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3))
        CALL DOT(PC,ucell(1,2),YIJ) 
        MRCY=INT(RIPSC*NFFT2/ABS(YIJ)+ONE)
!  electrostatic parameters
        CIPSE(0)=-35.0/16.0
        CIPSE(1)=35.0/16.0
        CIPSE(2)=-21.0/16.0
        CIPSE(3)=5.0/16.0
        PEC=ZERO
        PEC=ONE+CIPSE(0)+CIPSE(1)+CIPSE(2)+CIPSE(3) 
!   r6 parameters
        CIPSVC(0)=7.0/16.0
        CIPSVC(1)=9.0/14.0
        CIPSVC(2)=-3.0/28.0
        CIPSVC(3)=6.0/7.0
        PVCC=ONE+CIPSVC(0)+CIPSVC(1)+CIPSVC(2)+CIPSVC(3) 
      ENDIF
! Precalculate IPS constants
      U2=RIPS2/RIPSC2
      U4=U2*U2
      U8=U4*U4
      DIPSE(1)=TWO*CIPSE(1)
      DIPSE(2)=FOUR*CIPSE(2)
      DIPSE(3)=SIX*CIPSE(3)
      CIPSDE(0)=CIPSE(0)*RIPSCR-AIPSE(0)*RIPSR
      CIPSDE(1)=CIPSE(1)*U2*RIPSCR-AIPSE(1)*RIPSR
      CIPSDE(2)=CIPSE(2)*U4*RIPSCR-AIPSE(2)*RIPSR
      CIPSDE(3)=CIPSE(3)*U4*U2*RIPSCR-AIPSE(3)*RIPSR
      PDEC=PIPSEC*RIPSR-PEC*RIPSCR
      DIPSDE(1)=TWO*CIPSDE(1)
      DIPSDE(2)=FOUR*CIPSDE(2)
      DIPSDE(3)=SIX*CIPSDE(3)
      DIPSVC(1)=TWO*CIPSVC(1)
      DIPSVC(2)=FOUR*CIPSVC(2)
      DIPSVC(3)=SIX*CIPSVC(3)
      CIPSDVC(0)=CIPSVC(0)*RIPSC6R-AIPSVC(0)*RIPS6R
      CIPSDVC(1)=CIPSVC(1)*U2*RIPSC6R-AIPSVC(1)*RIPS6R
      CIPSDVC(2)=CIPSVC(2)*U4*RIPSC6R-AIPSVC(2)*RIPS6R
      CIPSDVC(3)=CIPSVC(3)*U4*U2*RIPSC6R-AIPSVC(3)*RIPS6R
      DIPSDVC(1)=TWO*CIPSDVC(1)
      DIPSDVC(2)=FOUR*CIPSDVC(2)
      DIPSDVC(3)=SIX*CIPSDVC(3)
      PDVCC=PIPSVCC*RIPS6R-PVCC*RIPSC6R
!  Cutoff grids
      DO 300 IBINZ1=-MRCZ,MRCZ
        ZI=-IBINZ1*ZBIN
        IBINZ=MOD(IBINZ1+100*NFFT3,NFFT3)+1
        IF(IBINZ<KBOT .OR. IBINZ > KTOP)GOTO 300    
        IBINZ0=(IBINZ-KBOT)*NFFTDIMRC*NFFTDIM2
        DO 400 IBINY1=-MRCY,MRCY
          YI=-IBINY1*YBIN
          IBINY=MOD(IBINY1+100*NFFT2,NFFT2)
          IBINY0=IBINZ0+IBINY*NFFTDIMRC
          DO 500 IBINX1=-MRCX,MRCX
            XI=-IBINX1*XBIN
            IBINX=MOD(IBINX1+100*NFFT1,NFFT1)
            IBINX0=IBINY0+IBINX
            I000=IBINX0+1
            XIJ=XI*WRK11+YI*WRK12+ZI*WRK13
            YIJ=XI*WRK21+YI*WRK22+ZI*WRK23
            ZIJ=XI*WRK31+YI*WRK32+ZI*WRK33
            R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            EELIJ=ZERO
            ENBIJ=ZERO
            IF(R2>RIPSC2)GOTO 500
            IF(R2<RIPS2)THEN
              U2=R2*RIPS2R
              R2R=R2/(R2*R2+TEN_TO_MINUS10)
!  Electrostatic IPS
! 
              PE=CIPSDE(0)+U2*(CIPSDE(1)+U2*(CIPSDE(2)+U2*(CIPSDE(3))))+PDEC
              DPE=U2*(DIPSDE(1)+U2*(DIPSDE(2)+U2*(DIPSDE(3))))
              EELIJ=EELIJ+PE
              DEIJ=DPE
!  Lennard-Jones IPS
              U4=U2*U2
              U6R=ONE/U4/U2
              U12R=U6R*U6R
!  L-J r6 term
!
        PVC=CIPSDVC(0)+U2*(CIPSDVC(1)+U2*(CIPSDVC(2)+U2*(CIPSDVC(3))))+PDVCC
              DPVC=U2*(DIPSDVC(1)+U2*(DIPSDVC(2)+U2*(DIPSDVC(3))))
!  L-J r12 term --neglected
!
              ENBIJ=ENBIJ-PVC
              DVIJ=-DPVC
          ELSE 
            R2R=ONE/R2
            U2=R2*RIPSC2R
!  Electrostatic IPS
!
            U1=SQRT(U2)
            PE=ONE/U1+CIPSE(0)+U2*(CIPSE(1)+U2*(CIPSE(2)+U2*(CIPSE(3))))-PEC
            DPE=-ONE/U1+U2*(DIPSE(1)+U2*(DIPSE(2)+U2*(DIPSE(3))))
            EELIJ=EELIJ+PE*RIPSCR
            DEIJ=DPE*RIPSCR
!  Lennard-Jones IPS
            U4=U2*U2
            U6R=ONE/U4/U2
            U12R=U6R*U6R
!  L-J r6 term
!
        PVC=U6R+CIPSVC(0)+U2*(CIPSVC(1)+U2*(CIPSVC(2)+U2*(CIPSVC(3))))-PVCC
        DPVC=-SIX*U6R+U2*(DIPSVC(1)+U2*(DIPSVC(2)+U2*(DIPSVC(3))))
!  L-J r12 term --neglected
!
            ENBIJ=ENBIJ-PVC*RIPSC6R
            DVIJ=-DPVC*RIPSC6R
          ENDIF
!  Add boundary tensor
            DVIJ=(DVIJ)*R2R
            DEIJ=(DEIJ)*R2R
            VDWXX(I000)=VDWXX(I000)+XIJ*XIJ*DVIJ
            VDWXY(I000)=VDWXY(I000)+XIJ*YIJ*DVIJ
            VDWYY(I000)=VDWYY(I000)+YIJ*YIJ*DVIJ
            VDWXZ(I000)=VDWXZ(I000)+XIJ*ZIJ*DVIJ
            VDWYZ(I000)=VDWYZ(I000)+YIJ*ZIJ*DVIJ
            VDWZZ(I000)=VDWZZ(I000)+ZIJ*ZIJ*DVIJ
            ELEXX(I000)=ELEXX(I000)+XIJ*XIJ*DEIJ
            ELEXY(I000)=ELEXY(I000)+XIJ*YIJ*DEIJ
            ELEYY(I000)=ELEYY(I000)+YIJ*YIJ*DEIJ
            ELEXZ(I000)=ELEXZ(I000)+XIJ*ZIJ*DEIJ
            ELEYZ(I000)=ELEYZ(I000)+YIJ*ZIJ*DEIJ
            ELEZZ(I000)=ELEZZ(I000)+ZIJ*ZIJ*DEIJ
            VDWARRAY(I000)=VDWARRAY(I000)+ENBIJ
            ELEARRAY(I000)=ELEARRAY(I000)+EELIJ
500       CONTINUE 
400     CONTINUE 
300   CONTINUE
!
        FIPS=FOUR*PI*RIPSC2*RIPSC/THREE/VOLUME
        IF(TVAIPS.and.master.AND. .NOT.TIPSFIN)THEN 
!  Only the 6th-term uses ripsc as cutoff
          EIPSANB=-FIPS*CIJSUM*PVCC*RIPSC6R
        ELSE
          EIPSANB=ZERO
        ENDIF
        IF(TEAIPS.and.master.AND. .NOT.TIPSFIN)THEN 
          EIPSAEL=FIPS*CGSUM*CGSUM*PEC/RIPSC/TWO
        ELSE
          EIPSAEL=ZERO
        ENDIF
        VIRAIPS=-(EIPSANB+EIPSAEL)
      RETURN
      END SUBROUTINE PBC_IPS_ENG

!-------------------------------------------------------------------
!                                -----
!     --- FILL_CHARGE_GRID -- RC-------
!                                -----
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_charge_grid here]
#ifdef MPI
subroutine fill_ips_grid(natom,charge, vdws,&
      theta1,theta2,theta3,fr1,fr2,fr3, &
      q,w,my_cg,nmine)
#else
subroutine fill_ips_grid(natom,charge, vdws,&
      theta1,theta2,theta3,fr1,fr2,fr3, &
      q,w,nmine)
#endif


   !---------------------------------------------------------------------
   ! INPUT:
   !      natom:  number of atoms
   !      charge: the array of atomic charges
   !      theta1,theta2,theta3: the spline coeff arrays
   !      fr1,fr2,fr3 the scaled and shifted fractional coords
   !      nfft1,nfft2,nfft3: the charge grid dimensions
   !      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
   !      order: the order of spline interpolation
   ! OUTPUT:
   !      Q the charge grid
   !---------------------------------------------------------------------

  use ew_bspline,only:kbot,ktop
  
   implicit none
   integer natom
   _REAL_ fr1(natom),fr2(natom),fr3(natom)
   _REAL_ theta1(forder,natom),theta2(forder,natom), &
         theta3(forder,natom),charge(natom),vdws(natom)
   _REAL_ q(*),w(*)
   integer nmine

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
   integer kbot0
   integer my_cg(*)
#endif

   integer n,ith1,ith2,ith3,i0,j0,k0,k,i00,j00,iqk,iqj
   _REAL_ prod,prodw
   integer im,kq,ntot
   
   !........Zero the Charge grids
   
   
#ifdef MPI
   kbot0 = mxystart(mytaskid)
   kbot = kbot0 + 1
   ktop = kbot0 + mxyslabs
   ntot = 2*nfftdim1*nfftdim2*mxyslabs
#else
   ntot = 2*nfftdim1*nfftdim2*nfftdim3
#endif
   q(1:ntot)=0.0d0
   w(1:ntot)=0.0d0
   
   do im = 1,nmine
#ifdef MPI
      n = my_cg(im)
#else
      n = im
#endif
      k0 = int(fr3(im)) - forder
      j00 = int(fr2(im)) - forder
      i00 = int(fr1(im)) - forder
      do ith3 = 1,forder
         k0 = k0 + 1
         if (k0 >= 0) then
            k = k0 + 1
         else
            k = k0 + 1 + nfft3
         end if

#ifdef MPI
         if ( k >= kbot .and. k <= ktop ) then
            kq = k - kbot0
#else
            kq = k
#endif
            iqk = (kq-1)*2*nfftdim1*nfftdim2
            j0 = j00
            do ith2 = 1,forder
               j0 = j0 + 1

               iqj = iqk + j0*2*nfftdim1
               if( j0 < 0 ) iqj = iqj + 2*nfft2*nfftdim1

               prod = theta2(ith2,im)*theta3(ith3,im)*charge(n)
               prodw = theta2(ith2,im)*theta3(ith3,im)*vdws(n)
               i0 = i00 + 1
               do ith1 = 1,forder
                  i0 = i0 + 1
                  if (i0 >= 1) then
                     q(i0+iqj) = q(i0+iqj) + theta1(ith1,im)*prod
                     w(i0+iqj) = w(i0+iqj) + theta1(ith1,im)*prodw
                  else
                     q(i0+nfft1+iqj) = q(i0+nfft1+iqj) + theta1(ith1,im)*prod
                     w(i0+nfft1+iqj) = w(i0+nfft1+iqj) + theta1(ith1,im)*prodw
                  end if
               end do
            end do
#ifdef MPI
         end if
#endif
      end do  !  ith3 = 1,order
   end do  !  im = 1,nmine
   return
end subroutine fill_ips_grid 

!-------------------------------------------------------------------
!  AIPS_GRAD_SUM RC    **** REAL not complex ****
!-----------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grad_sumrc here]
#ifdef MPI
subroutine aips_grad_sumrc( &
      natom,charge,vdws,recip,theta1,theta2,theta3,  &
      dtheta1,dtheta2,dtheta3,frcx,frc,fr1,fr2,fr3, &
      q,w,my_cg,nmine)
#else
subroutine aips_grad_sumrc( &
      natom,charge,vdws,recip,theta1,theta2,theta3,  &
      dtheta1,dtheta2,dtheta3,frcx,frc,fr1,fr2,fr3, &
      q,w,nmine)
#endif
  use constants, only:zero,one
  use ew_bspline,only:kbot,ktop
  use decomp, only: decpr, decpair
  use file_io_dat
  
   implicit none
   integer natom
   integer kq
#include "box.h"
#include "md.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
   include "mpif.h"
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
   _REAL_ q(*),w(*)
#else
   _REAL_ q(*),w(*)
   
   
#endif

   _REAL_ recip(3,3),recip11,recip22,recip33
   _REAL_ fr1(natom),fr2(natom),fr3(natom)
   _REAL_ frcx(3),frc(3,natom)
   _REAL_ theta1(forder,natom),theta2(forder,natom), &
         theta3(forder,natom),charge(natom),vdws(natom)
   _REAL_ dtheta1(forder,natom),dtheta2(forder,natom), &
         dtheta3(forder,natom)
#ifdef MPI
   integer my_cg(*)
#endif
   integer nmine
   integer iqk,iqj,j00,i00

   integer n,ith1,ith2,ith3,i0,j0,k0,k,im
   _REAL_ f1,f2,f3,term,chargen,vdwn
   _REAL_ dfx,dfy,dfz,dnfft1,dnfft2,dnfft3
   _REAL_ cfact,f1fac,f2fac,f3fac
   _REAL_ dectmp

   CFACT=one/(NFFT1*NFFT2*NFFT3)

   dnfft1 = nfft1
   dnfft2 = nfft2
   dnfft3 = nfft3
   
   recip11 = recip(1,1)*dnfft1
   recip22 = recip(2,2)*dnfft2
   recip33 = recip(3,3)*dnfft3

      do im = 1,nmine
#ifdef MPI
         n = my_cg(im)
#else
         n = im
#endif
         f1 = zero
         f2 = zero
         f3 = zero
         i00 = int(fr1(im)) - forder
         j00 = int(fr2(im)) - forder
         k0 = int(fr3(im)) - forder
         chargen = cfact*charge(n)
         vdwn = cfact*vdws(n)
         do ith3 = 1,forder
            k0 = k0 + 1
            if (k0 >= 0) then
               k = k0 + 1
            else
               k = k0 + 1 + nfft3
            end if
            
#ifdef MPI
            if ( k >= kbot .and. k <= ktop )then
               kq = k - mxystart(mytaskid)
#else
               kq = k
#endif
               iqk = (kq-1)*2*nfftdim1*nfftdim2
               j0 = j00
               do ith2 = 1,forder
                  j0 = j0 + 1
                  
                  iqj = iqk + j0*2*nfftdim1
                  if( j0 < 0 ) iqj = iqj + 2*nfft2*nfftdim1
                  
                  i0 = i00 + 1
                  f1fac =  theta2(ith2,im) *  theta3(ith3,im) 
                  f2fac = dtheta2(ith2,im) *  theta3(ith3,im) 
                  f3fac =  theta2(ith2,im) * dtheta3(ith3,im) 
                  do ith1 = 1,forder
                     i0 = i0 + 1
                     if (i0 >= 1) then
                        term = q(i0 + iqj)*chargen+w(i0+iqj)*vdwn
                     else
                        term = q(i0 + nfft1 + iqj)*chargen+w(i0 + nfft1 + iqj)*vdwn
                     end if
                    
                     !               ---force is negative of grad
                     
                     f1 = f1 - term * dtheta1(ith1,im) * f1fac
                     f2 = f2 - term *  theta1(ith1,im) * f2fac
                     f3 = f3 - term *  theta1(ith1,im) * f3fac
                     ! -- ti decomp
                     if(decpr .and. idecomp > 0) then
                        dectmp = term*theta1(ith1,im)*theta2(ith2,im)*theta3(ith3,im)*0.5d0
                        call decpair(2,n,n,dectmp/(nstlim/ntpr))
                     endif
                  end do
               end do
#ifdef MPI
            end if  ! ( k >= kbot .and. k <= ktop )
#endif
         end do  !  ith3 = 1,order
         
         if ( ifbox == 1 ) then  ! orthogonal unit cell
            dfx = recip11 * f1
            dfy = recip22 * f2
            dfz = recip33 * f3
         else
            f1 = dnfft1*f1
            f2 = dnfft2*f2
            f3 = dnfft3*f3
            dfx=recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
            dfy=recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
            dfz=recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
         end if
         frc(1,n) = frc(1,n) + dfx
         frc(2,n) = frc(2,n) + dfy
         frc(3,n) = frc(3,n) + dfz
         frcx(1)=frcx(1)+dfx
         frcx(2)=frcx(2)+dfy
         frcx(3)=frcx(3)+dfz
         
      end do  !  im = 1,nmine
   return
end subroutine aips_grad_sumrc 
!=================================================================
!     --- AIPS_SUM RC---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine scalar_sumrc here]
subroutine aips_sumrc( &
      q, w,elen,vdwn,qxx,qxy,qxz,qyy,qyz,qzz,wxx,wxy,wxz,wyy,wyz,wzz,&
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eerq,eerw,rec_vir)
   implicit none
#  include "extra.h"
#  include "box.h"
#  include "md.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
   integer ier,i
#endif
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   _REAL_ eerq,eerw,rec_vir(3,3)
   integer k2q
   _REAL_ q(2,nfft3,nfftdim1,nfft2),w(2,nfft3,nfftdim1,nfft2)
   _REAL_ elen(2,nfft3,nfftdim1,nfft2),vdwn(2,nfft3,nfftdim1,nfft2)
   _REAL_ qxx(2,nfft3,nfftdim1,nfft2),wxx(2,nfft3,nfftdim1,nfft2)
   _REAL_ qxy(2,nfft3,nfftdim1,nfft2),wxy(2,nfft3,nfftdim1,nfft2)
   _REAL_ qxz(2,nfft3,nfftdim1,nfft2),wxz(2,nfft3,nfftdim1,nfft2)
   _REAL_ qyy(2,nfft3,nfftdim1,nfft2),wyy(2,nfft3,nfftdim1,nfft2)
   _REAL_ qyz(2,nfft3,nfftdim1,nfft2),wyz(2,nfft3,nfftdim1,nfft2)
   _REAL_ qzz(2,nfft3,nfftdim1,nfft2),wzz(2,nfft3,nfftdim1,nfft2)
   _REAL_ enq,enw
   integer k1,k2,k3,m1,m2,m3,nff
   integer nf1,nf2,nf3
   _REAL_ struc2q,struc2w
   _REAL_ tmp1,tmp2
   _REAL_ cfact,cfact1,cfact2,cfact3,cfact4
   _REAL_ qi,qj,wi,wj,pqi,pqj,pwi,pwj

   CFACT=0.5d0/(NFFT1*NFFT2*NFFT3)
   nff = nfft1*nfft2
   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1
   enq = 0.0d0
   enw = 0.0d0
#ifndef noVIRIAL
   do m2 = 1,3
      do m1 = 1,3
         rec_vir(m1,m2) = 0.d0
      end do
   end do
#endif

   !======================================================================
   !        BIG LOOP
   !======================================================================

#ifdef MPI
   do k2q = 1, mxzslabs
      if(master)then
         k2=k2q
      else
         k2 = k2q + mxzstart(mytaskid)
      end if
#else
   do k2q = 1, nfft2
      k2=k2q
#endif
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      
      CFACT1=prefac2(K2)
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         CFACT2=CFACT1*prefac3(K3)
         do k1 = 1, nf1+1
            m1 = mod(nfft1-k1 + 1,nfft1)+1
            CFACT3=CFACT2*prefac1(K1)
            if ( k1 == m1 ) then
                CFACT4=CFACT3
            else
                CFACT4=2.0d0*CFACT3
            endif
            QI=Q(1,k3,k1,k2q) 
            QJ=Q(2,k3,k1,k2q)
            PQI=ELEn(1,k3,k1,k2q) 
            PQJ=ELEn(2,k3,k1,k2q)
            Q(1,k3,k1,k2q) = CFACT3*(QI*PQI-QJ*PQJ)
            Q(2,k3,k1,k2q) = CFACT3*(QI*PQJ+QJ*PQI)
            WI=W(1,k3,k1,k2q)
            WJ=W(2,k3,k1,k2q)
            PWI=VDWn(1,k3,k1,k2q) 
            PWJ=VDWn(2,k3,k1,k2q)
            W(1,k3,k1,k2q) = CFACT3*(WI*PWI-WJ*PWJ)
            W(2,k3,k1,k2q) =CFACT3*(WI*PWJ+WJ*PWI)
            struc2q = CFACT4*(QI*QI + QJ*QJ)
            tmp1 = PQI *struc2q
            enq = enq + tmp1
            struc2w = CFACT4*(WI*WI + WJ*WJ)
            tmp2 = PWI *struc2w
            enw = enw + tmp2
#ifndef noVIRIAL
            rec_vir(1,1) = rec_vir(1,1) + (qxx(1,k3,k1,k2q)*struc2q+wxx(1,k3,k1,k2q)*struc2w)
            rec_vir(1,2) = rec_vir(1,2) + (qxy(1,k3,k1,k2q)*struc2q+wxy(1,k3,k1,k2q)*struc2w)
            rec_vir(1,3) = rec_vir(1,3) + (qxz(1,k3,k1,k2q)*struc2q+wxz(1,k3,k1,k2q)*struc2w)
            rec_vir(2,2) = rec_vir(2,2) + (qyy(1,k3,k1,k2q)*struc2q+wyy(1,k3,k1,k2q)*struc2w)
            rec_vir(2,3) = rec_vir(2,3) + (qyz(1,k3,k1,k2q)*struc2q+wyz(1,k3,k1,k2q)*struc2w)
            rec_vir(3,3) = rec_vir(3,3) + (qzz(1,k3,k1,k2q)*struc2q+wzz(1,k3,k1,k2q)*struc2w)
#endif
            

         end do
      end do
   end do

   eerq = eipsael+cfact * enq
   eerw = eipsanb+cfact * enw

#ifndef noVIRIAL
   do m2 = 1,3
      do m1 = 1,m2
         rec_vir(m1,m2) = cfact*rec_vir(m1,m2)
      end do
      rec_vir(m2,m2) = rec_vir(m2,m2)+viraips
   end do
   rec_vir(2,1) = rec_vir(1,2)
   rec_vir(3,1) = rec_vir(1,3)
   rec_vir(3,2) = rec_vir(2,3)
#endif

   return
end subroutine aips_sumrc 

      subroutine setndim(nin,nout)
!-----------------------------------------------------------------------
!  This routine set an even integer with factors of only 2,3, or 5
!
      implicit none
      INTEGER nin,nout,n
      nout=((nin+1)/2)*2
10    n=nout     
20    IF(N.EQ.2*(N/2))THEN
        N=N/2
        GOTO 20
      ENDIF
30    IF(N.EQ.3*(N/3))THEN
        N=N/3
        GOTO 30
      ENDIF
50    IF(N.EQ.5*(N/5))THEN
        N=N/5
        GOTO 50
      ENDIF
      if(n>1)then
        nout=nout+2
        goto 10
      endif 
      return
      end subroutine setndim

      SUBROUTINE AIPSCELL(NATOM,X)
!-----------------------------------------------------------------------
!  This routine define large box to enclude a finite system
!
      use nblist,only:a,b,c,ucell,recip,volume
      implicit none
      INTEGER NATOM
      _REAL_ X(*)
      _REAL_ extents(3,2)
   _REAL_ u23(3),u31(3),u12(3)
   _REAL_ onevolume
      INTEGER i,j,n,nx,ny,nz
   TIPSMEM=.FALSE.
   FORDER=MIPSO
   if(FORDER<4)FORDER=4
   if(tipsfin)then
     do i=1,3
       extents(i,1)=x(i)
       extents(i,2)=x(i)
     end do
     do i=3,3*natom-1,3
       do j=1,3
         extents(j,1)=min(extents(j,1),x(i+j))
         extents(j,2)=max(extents(j,2),x(i+j))
       end do
     end do
     n=2*INT((extents(1,2)-extents(1,1))/GRIDIPS+2+MIPSO)
     call setndim(n,nx)
     a=nx*GRIDIPS
     n=2*INT((extents(2,2)-extents(2,1))/GRIDIPS+2+MIPSO)
     call setndim(n,ny)
     b=ny*GRIDIPS
     n=2*INT((extents(3,2)-extents(3,1))/GRIDIPS+2+MIPSO)
     call setndim(n,nz)
     c=nz*GRIDIPS
     ucell(1,1) = a
     ucell(2,1) = 0.d0
     ucell(3,1) = 0.d0
     ucell(1,2) = 0.d0
     ucell(2,2) = b
     ucell(3,2) = 0.d0
     ucell(1,3) = 0.d0
     ucell(2,3) = 0.d0
     ucell(3,3) = c
   !  now get reciprocal vectors

     call cross(ucell(1,2),ucell(1,3),u23)
     call cross(ucell(1,3),ucell(1,1),u31)
     call cross(ucell(1,1),ucell(1,2),u12)
     call dot(ucell(1,1),u23,volume)
     onevolume=1.0d0/volume
     do j = 1,3
       recip(j,1) = u23(j)*onevolume
       recip(j,2) = u31(j)*onevolume
       recip(j,3) = u12(j)*onevolume
     end do
   else
     if(mipsx>0)then
        call setndim(mipsx,nx)
     else
        i=anint(a/gridips)
        call setndim(i,nx)
     endif
     if(mipsy>0)then
        call setndim(mipsy,ny)
     else
        i=anint(b/gridips)
        call setndim(i,ny)
     endif
     if(mipsz>0)then
        call setndim(mipsz,nz)
     else
        i=anint(c/gridips)
        call setndim(i,nz)
     endif
   endif
   if(nfft1/=nx)then
     nfft1=nx
     tipsmem=.true.
   endif
   if(nfft2/=ny)then
     nfft2=ny
     tipsmem=.true.
   endif
   if(nfft3/=nz)then
     nfft3=nz
     tipsmem=.true.
   endif
   RETURN
   END SUBROUTINE AIPSCELL


end module nbips
