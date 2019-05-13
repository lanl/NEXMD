! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine to compute energies and forces
subroutine force(xx,ix,ih,ipairs,x,f,ener,vir, &
      fs, rborn, reff, onereff, qsetup,do_list_update,nstep)

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : ncsu_on_force => on_force
#  if !defined(LES) && defined(MPI)
   use remd, only : rem, mdloop
#  endif
#  ifdef MPI
   use ncsu_sander_proxy, only : ncsu_remember_initremd => remember_initremd
#  endif
#endif
   
   use file_io_dat
   use genborn
   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : npopt, np_force
   use pbtimer_module, only : pbtimer_init, pbtimer_summary
#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm, rism_force
#endif
#ifdef APBS
   use apbs
#endif /* APBS */
   use trace
   use stack
   use pimd_vars, only: ipimd,nbead,bnd_vir,Epot_spring,Epot_deriv,real_mass,equal_part,nrg_all,nebrms,itimass
   use neb_vars, only: ineb, neb_force
   use full_pimd_vars, only: totener
   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, &
                           qmmm_mpi, qmewald, qmmm_scratch
   use constants, only : zero
   use relax_mat
   use ew_recip
   use parms, only:cn1,cn2,asol,bsol,pk,rk,tk,numbnd,numang,nptra,nphb,nimprp,&
                   cn3,cn4,cn5!mjhsieh: for another vdwmodel
#ifdef PUPIL_SUPPORT

   ! Using the ucell variable in nblist
   use nblist, only:nonbond_list,a,b,c,alpha,beta,gamma,ucell
#else
   ! Only used in qmmm call
   use nblist, only:nonbond_list,a,b,c,alpha,beta,gamma
#endif /*PUPIL_SUPPORT*/

#ifdef DSSP
   use dssp, only: fdssp, edssp, idssp
#endif


   use amoeba_interface,only: AM_VAL_eval,AM_NonBond_eval
   use amoeba_mdin, only : iamoeba,am_nbead

   use amd_mod
   use nbips,only: ips,eexips
   use emap,only:temap,emapforce

#if defined(MPI)
   use evb_data, only: nrg_frc
   use softcore, only: sc_ener
#  ifdef LES
      use miller, only: dlnQ_dl
      use remd, only : rem, mdloop
#  endif
#endif /* MPI */

#ifdef PUPIL_SUPPORT
   use pupildata
#endif /*PUPIL_SUPPORT*/

   use linear_response, only: ilrt, ee_linear_response, energy_m0, energy_w0, &
        energy_vdw0, cn1_lrt, cn2_lrt, crg_m0, crg_w0, do_lrt, f_scratch, &
        lrt_solute_sasa

   use cns_xref
#ifdef _XRAY
   use xray_interface_module, only: xray_get_derivative, xray_active
#endif

!CHARMM Force Field Support
   use charmm_mod, only : charmm_active, charmm_calc_impropers, &
                          charmm_calc_cmap, charmm_calc_urey_bradley,&
                          charmm_dump_gold, &
                          do_charmm_dump_gold
   use ff11_mod, only : cmap_active, calc_cmap
   use decomp, only: &
#ifdef MPI
      collect_dec2,  &
#endif
      init_dec
   use state

   use crg_reloc, only: ifcr, cr_reassign_charge, cr_calc_force

   implicit none
  
   integer, intent(in) :: nstep

#ifdef PUPIL_SUPPORT
   character(kind=1,len=5) :: routine="force"
#endif
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ fs(*),rborn(*),reff(*),dvdl
   _REAL_, intent(out) :: onereff(*)
#include "def_time.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "extra_pts.h"
#include "parallel.h"

#ifdef MPI
#include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
   integer gb_rad_mpistart
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer ierr
#ifdef MPI
   _REAL_ etmp
#endif 
   logical belly
#include "md.h"
#include "pb_md.h"
#include "box.h"
#include "nmr.h"
#include "memory.h"
#include "extra.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "flocntrl.h"
#include "les.h"
   integer istart,iend
   _REAL_ evdwex, eelex
   _REAL_ escf
   _REAL_ enemap

   logical, intent(inout) :: qsetup
   logical, intent(out) :: do_list_update

   _REAL_  enmr(6),devdis(4),devang(4),devtor(4),devpln(4),devplpt(4),devgendis(4),entr,ecap
   _REAL_  x(*),f(*),vir(4)
   type(state_rec)  ener

   !Local
   _REAL_                     :: ene(30)    !Used locally ONLY
   type(potential_energy_rec) :: pot        !Used locally ONLY

#ifdef LES
   _REAL_  :: nrg_bead(nbead)
#endif

   integer ibead,bead_per_node
   _REAL_  virtmp(4),am_Ebnd,am_Eang,am_Edih,am_Evdw,am_Eelt,am_Enb14,am_Ee14,am_Epolar, &
           am_Espring,am_Ederiv

   integer i,m,i3,j,j3
   _REAL_  virvsene,eelt,epol,esurf,edisp,enpol
#ifdef RISMSANDER
   _REAL_ erism
#endif /*RISMSANDER*/

   _REAL_ ect ! charge transfer

   _REAL_  epolar,aveper,aveind,avetot,emtot,dipiter,dipole_temp
   integer l_r2x,l_rjx,l_tmp1,l_tmp2,l_tmp3,l_tmp4,l_tmp5
   integer l_tmp6,l_tmp7,l_tmp8,l_jj,l_skipv, l_kvls,l_jvls,l_psi
   integer l_da,l_sumd
   integer newbalance
   save newbalance
   
!AMD variables
   _REAL_ amd_totdih, temp_amd_totdih
#if defined(MPI)
   integer :: status( MPI_STATUS_SIZE )
   _REAL_ :: vel0_nrg_sum
#endif /* MPI */

   ect = 0.0

   call trace_enter( 'force' )

   call timer_start(TIME_FORCE)
   if ( idecomp /= 0 .and. icfe == 0) call init_dec
   ene(:) = ZERO 
   call zero_pot_energy(pot)

   belly = ibelly > 0

#ifdef MPI
   if (mpi_orig) then

      !     Check to see if we are done yet in mpi_orig case (tec3).
      !     This is done by monitoring the status of an integer notdone.
      !     If notdone .eq. 1 then we keep going.  notdone is set to zero
      !     when we no longer want to call force().  This perhaps is not the
      !     most efficient means to implement this check...

      call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
      if (notdone /= 1) return

      !       Send copies of xyz coords, setbox common block, vir array
      !       and NTNB value to all nodes from master with a broadcast.

      if (numtasks > 1) then

         call mpi_bcast(box,BC_BOXR,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(ntb,BC_BOXI,mpi_integer,0,commsander,ierr)
         call mpi_bcast(vir,3,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(xx(lcrd),3*natom,MPI_DOUBLE_PRECISION, &
               0,commsander,ierr)
         call mpi_bcast(ntnb,1,mpi_integer,0,commsander,ierr)
         if (iabs(ntb) >= 2) then
            call mpi_bcast(xx(l45),3*natom,MPI_DOUBLE_PRECISION, &
                  0,commsander,ierr)
         end if
      end if
   end if

   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = natom
#endif

   !-----------------------------------
   !   QMMM Variable QM Solvent Scheme
   !-----------------------------------
   if (qmmm_nml%ifqnt .and. qmmm_nml%vsolv > 0) then
     call timer_start(TIME_QMMM)
     call timer_start(TIME_QMMMVARIABLESOLVCALC)
     ! For the moment this is NOT parallel - all threads need to call this.
     call qmmm_vsolv_update(nstep, natom, nres, ix(i70), xx(l15), x, xx(lmass), &
                             nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                             ix(iibh),ix(ijbh),ix(iicbh), &
                             ix(iiba),ix(ijba),ix(iicba), &
                             ix(i24),ix(i26),ix(i28),ix(i30),ix(i32),ix(i34),ix(i36),ix(i38), &
                             ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),ix(i50),ix(i52),ix(i54), &
                             ix(i56),ix(i58), ix(ibellygp))
     call timer_stop(TIME_QMMMVARIABLESOLVCALC)
     call timer_stop(TIME_QMMM)
   end if 
   !-----------------------------------
   ! END QMMM Variable QM Water Scheme
   !-----------------------------------

   if(iamoeba.eq.1) then
      REQUIRE(am_nbead.eq.ncopy)
   end if

   !     ----- ZERO OUT THE ENERGIES AND FORCES -----
   enoe = 0.d0
   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   dipiter=0.d0
   dvdl=0.d0
   dipole_temp=0.d0
   enmr(1:6) = 0.d0
   vir(1:4) = 0.d0
   virvsene = 0.d0
   f(1:3*natom+iscale) = 0.d0
#ifdef LES
   if( ipimd>0) nrg_all(1:nbead)=0.d0
#endif

#ifndef PUPIL_SUPPORT
   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_NONBON)
      call timer_start(TIME_LIST)

      call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
               ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
               ix(ibellygp),belly,newbalance,cn1, &
               xx(lvel),xx(lvel2),ntp,xx(l45), qsetup, &
               do_list_update)
      call timer_stop(TIME_LIST)
      call timer_stop(TIME_NONBON)
   end if
   ! charge reassign here !
   if ( ifcr /= 0 ) then
      call cr_reassign_charge( x, f, pot%ct, xx(l15), natom )
   end if
#endif
   
#ifndef DISABLE_NCSU
#  ifdef MPI
   call ncsu_remember_initremd(rem.gt.0.and.mdloop.eq.0)
#  endif
   call ncsu_on_force(x, f, vir)
#endif

   ! ----------------------------------------------------------------
   ! Do weight changes, if requested
   ! ----------------------------------------------------------------

   if (nmropt > 0) &
         call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
         devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb,xx(lnmr01), &
         ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
         cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
         nimprp,nphb,natom,natom,ntypes,nres, &
         rad,wel,radhb,welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'WEIT')
   ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints

   epolar = 0.d0

   ! -----------------------------------------------------------------
   ! EGB if igb>0 and /=6 then we need to calculate the GB radii for this
   ! structure
   ! -----------------------------------------------------------------

   ! If we are doing qm_gb=2 then we need to calculate the GB radii
   ! before calling qm_mm
   if( igb > 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. &
      ( irespa < 2 .or. mod(irespa,nrespai) == 0) ) then
#ifdef MPI
      gb_rad_mpistart = mytaskid+1
#endif
      call timer_start(TIME_EGB)
      call timer_start(TIME_GBRAD1)
      !If qmmm and then this will calculate us the radii for the
      !link atoms and not the mm link pair atoms. The initial radii used are those
      !of the mm link pair's atom type though.
      !The MMlink pair atoms must be the coordinates of the link atoms here.
      if (qmmm_nml%ifqnt) call adj_mm_link_pair_crd(qmmm_nml,qmmm_struct,x)

      call egb_calc_radii(igb,natom,x,fs,reff, &
                     onereff,fsmax,rgbmax, rborn, offset, &
                     rbornstat,xx(l188),xx(l189),         &
                     xx(l186),xx(l187), gbneckscale, ncopy, rdt, &
                     xx(l2402),xx(l2403),xx(l2404),  & ! Hai Nguyen: gbvalpha,gbvbeta,gbvgamma 
                     gbalpha, gbbeta, gbgamma & !Hai Nguyen: keep for sander.LES
#ifdef MPI
                     ,gb_rad_mpistart &
#endif
                       )
      if (qmmm_nml%ifqnt) call rst_mm_link_pair_crd(qmmm_nml,qmmm_struct, x)

      call timer_stop(TIME_GBRAD1)
      call timer_stop(TIME_EGB)

   end if
 
   ! QM/MM Contributions are now calculated before the NON-Bond info.
   ! ----------------------------------------------------------------
   ! Calculate the qm/mm contributions
   ! ----------------------------------------------------------------

   if(qmmm_nml%ifqnt) then

      ! If we are doing periodic boundaries with QM/MM PME then we need to
      ! do the PME calculation twice. First here to get the potential at
      ! each QM atom due to the PME and then again after all the QM is done
      ! to get the MM-MM potential and all of the gradients.

      if(qmmm_nml%qm_pme) then
         ! Ewald force will put the potential into the qm_ewald%mmpot array.
         call timer_start(TIME_EWALD)
         call ewald_force(x,natom,ix(i04),ix(i06), &
               xx(l15),cn1,cn2,asol,bsol,eelt,epolar, &
!               f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),.true. &
! Modified by WJM, YD, mjhsieh
               f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),xx(ldf), &
               xx(lpol2),.true., cn3, cn4, cn5 )
         call timer_stop(TIME_EWALD)
      endif

      call timer_start(TIME_QMMM)

#ifndef LES
        !========================================================
        !                      REGULAR QMMM
        !========================================================

         
      call qm_mm(x, natom,qmmm_struct%scaled_mm_charges, &
                 f,escf,periodic,reff,onereff, &
                intdiel,extdiel,Arad, cut,qm2_struct%scf_mchg,ntypes, &
                 ih(m06),xx(lmass), ix(i04),nstep)
        pot%scf = escf
#endif

        !========================================================
        !                  END REGULAR QMMM
        !========================================================
      call timer_stop(TIME_QMMM)
   end if 

   !---------------------------------------------------------------
   !END qm/mm contributions
   !---------------------------------------------------------------

#ifdef PUPIL_SUPPORT

   !*****************************************************
   !     Getting the Quantum forces with PUPIL package
   !*****************************************************

!  Reconstruct the simulation cell if there is any change
   do iPup=1,3    !vector loop
     do jPup=1,3  !Component loop
       qcell((iPup-1)*3+jPup) = ucell(jPup,iPup)
     enddo
   enddo
!  minimum point of the box ..... we assume (0,0,0) ????
   qcell(10) = 0.0d0
   qcell(11) = 0.0d0
   qcell(12) = 0.0d0

!  temporary vector to wrap the real coordinates to pass through
!  PUPIL interface.
   call get_stack(l_puptmp,3*natom,routine)
   if(.not. rstack_ok)then
     deallocate(r_stack)
     allocate(r_stack(1:lastrst),stat=alloc_ier)
     call reassign_rstack(routine)
   endif
   REQUIRE(rstack_ok)
   do iPup=1,3*natom
     r_stack(l_puptmp + iPup - 1) = x(iPup)
   end do

   if(ntb > 0) then
     call wrap_molecules(nspm,ix(i70),r_stack(l_puptmp))
     if(ifbox == 2) call wrap_to(nspm,ix(i70),r_stack(l_puptmp),box)
   end if


!  Preparing the coordinates, velocity and classic forces
!  to get quantum force
   do iPup=1,natom
     bs1 = (iPup-1)*9
     bs2 = (iPup-1)*3
     do jPup=1,3
       qcdata(bs1+jPup) = r_stack(l_puptmp + bs2 + jPup - 1)
       qcdata(bs1+3+jPup) = realStack(lvel+bs2+jPup-1)
       qcdata(bs1+6+jPup) = f(bs2+jPup)
     enddo
   enddo

!  Deallocating temporary stack
   call free_stack(l_puptmp,routine)

!  We are going to use the qmmm_nml and qmmm_struct variables to skip quantum atoms
!  in the force calculation
   qmmm_nml%ifqnt = .true.
   if(pupStep .EQ. 0) then
!!!  To keep initial values from the MD step 1
     ierr = 0
     allocate ( pupnb14(numnb14*3),stat=ierr)
     REQUIRE(ierr == 0)

     allocate ( pupbonh(nbonh*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupbona(nbona*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheth(ntheth*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheta(ntheta*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphih(nphih*5),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphia(nphia*5),stat=ierr)
     REQUIRE(ierr == 0)
     pupnbonh   = nbonh
     pupnbona   = nbona
     pupntheth  = ntheth
     pupntheta  = ntheta
     pupnphih   = nphih
     pupnphia   = nphia
     pupnumnb14 = numnb14
     call copy_14nb(ix(inb_14),pupnb14,numnb14)
     do iPup = 1,nbonh
       bs1 = iPup-1
       bs2 = bs1*3
       pupbonh(bs2+1)  = ix(iibh +bs1)
       pupbonh(bs2+2)  = ix(ijbh +bs1)
       pupbonh(bs2+3)  = ix(iicbh+bs1)
     enddo
     do iPup = 1,nbona
       bs1 = iPup-1
       bs2 = bs1*3
       pupbona(bs2+1)  = ix(iiba +bs1)
       pupbona(bs2+2)  = ix(ijba +bs1)
       pupbona(bs2+3)  = ix(iicba+bs1)
     enddo
     do iPup = 1,ntheth
       bs1 = iPup-1
       bs2 = bs1*4
       puptheth(bs2+1) = ix(i24  +bs1)
       puptheth(bs2+2) = ix(i26  +bs1)
       puptheth(bs2+3) = ix(i28  +bs1)
       puptheth(bs2+4) = ix(i30  +bs1)
     enddo
     do iPup = 1,ntheta
       bs1 = iPup-1
       bs2 = bs1*4
       puptheta(bs2+1) = ix(i32  +bs1)
       puptheta(bs2+2) = ix(i34  +bs1)
       puptheta(bs2+3) = ix(i36  +bs1)
       puptheta(bs2+4) = ix(i38  +bs1)
     enddo
     do iPup = 1,nphih
       bs1 = iPup-1
       bs2 = bs1*5
       pupphih(bs2+1)  = ix(i40  +bs1)
       pupphih(bs2+2)  = ix(i42  +bs1)
       pupphih(bs2+3)  = ix(i44  +bs1)
       pupphih(bs2+4)  = ix(i46  +bs1)
       pupphih(bs2+5)  = ix(i48  +bs1)
     enddo
     do iPup = 1,nphia
       bs1 = iPup-1
       bs2 = bs1*5
       pupphia(bs2+1)  = ix(i50  +bs1)
       pupphia(bs2+2)  = ix(i52  +bs1)
       pupphia(bs2+3)  = ix(i54  +bs1)
       pupphia(bs2+4)  = ix(i56  +bs1)
       pupphia(bs2+5)  = ix(i58  +bs1)
     enddo
   endif

!  Getting the quantum forces for a specific quantum domain
   pupStep  = pupStep + 1
   puperror = 0
   pupLevelData = 3
   call getquantumforces(natom,pupLevelData,pupStep,puperror,qcdata,qcell)
   if (puperror.ne.0) then
     write (6,*) 'Fatal error: Could not obtain quantum forces!'
     call mexit(6,1)
   endif
 ! Quantum energy treatment....
   pot%scf = qmEnergy

!  Deleting interactions CL-QZ if a new list of quantum atoms is given
   if (pupQZchange .ne. 0) then

!    ********** Rebuilding the nonbonding 14 list ***************
!    deleting all connectivity between the QM atoms

!      reinitializing internal nb 14 list structures from the beginning
       numnb14= pupnumnb14
       nbonh  = pupnbonh
       nbona  = pupnbona
       ntheth = pupntheth
       ntheta = pupntheta
       nphih  = pupnphih
       nphia  = pupnphia
       call copy_14nb(pupnb14,ix(inb_14),pupnumnb14)
       do iPup = 1,nbonh
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iibh +bs1) = pupbonh(bs2+1)
         ix(ijbh +bs1) = pupbonh(bs2+2)
         ix(iicbh+bs1) = pupbonh(bs2+3)
       enddo
       do iPup = 1,nbona
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iiba +bs1) = pupbona(bs2+1)
         ix(ijba +bs1) = pupbona(bs2+2)
         ix(iicba+bs1) = pupbona(bs2+3)
       enddo
       do iPup = 1,ntheth
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i24  +bs1) = puptheth(bs2+1)
         ix(i26  +bs1) = puptheth(bs2+2)
         ix(i28  +bs1) = puptheth(bs2+3)
         ix(i30  +bs1) = puptheth(bs2+4)
       enddo
       do iPup = 1,ntheta
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i32  +bs1) = puptheta(bs2+1)
         ix(i34  +bs1) = puptheta(bs2+2)
         ix(i36  +bs1) = puptheta(bs2+3)
         ix(i38  +bs1) = puptheta(bs2+4)
       enddo
       do iPup = 1,nphih
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i40  +bs1) = pupphih(bs2+1)
         ix(i42  +bs1) = pupphih(bs2+2)
         ix(i44  +bs1) = pupphih(bs2+3)
         ix(i46  +bs1) = pupphih(bs2+4)
         ix(i48  +bs1) = pupphih(bs2+5)
       enddo
       do iPup = 1,nphia
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i50  +bs1) = pupphia(bs2+1)
         ix(i52  +bs1) = pupphia(bs2+2)
         ix(i54  +bs1) = pupphia(bs2+3)
         ix(i56  +bs1) = pupphia(bs2+4)
         ix(i58  +bs1) = pupphia(bs2+5)
       enddo

       call deleting_qm_atoms()
       qsetup=.true.

!      Setting as current quantum zone
       pupQZchange = 0

   endif
  ! For PUPIL, rebuild the neighbour list 
  ! and zero the charges on QM atoms at every step
  if ( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then
         
    ! (for GB: do all nonbondeds together below)
    call timer_start(TIME_NONBON)
    call timer_start(TIME_LIST)
    !do_list_update=.true.         
    call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
                      ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
                      ix(ibellygp),belly,newbalance,cn1, &
                      xx(lvel),xx(lvel2),ntp,xx(l45), qsetup, &
                      do_list_update)
    !call qm_zero_charges(x(L15))
    call timer_stop(TIME_LIST)
    call timer_stop(TIME_NONBON)
  end if
  ! charge reassign here !
  if ( ifcr /= 0 ) then
     call cr_reassign_charge( x, f, pot%ct, xx(l15), natom )
  end if

#endif /*PUPIL_SUPPORT*/


  if(iamoeba.eq.1) then
      vir(1:4)=0.0
   end if
   ! ----------------------------------------------------------------
   ! Calculate the non-bonded contributions
   ! ----------------------------------------------------------------
   call timer_start(TIME_NONBON)


   call timer_start(TIME_EEXIPS)
   if( ips > 0 ) then
      call eexips(evdwex,eelex,istart,iend, ntb,ntypes, &
           ix(i04),ix(i06),ix(i08),ix(i10),xx(l15),cn1,cn2,f,x)
   endif
   call timer_stop(TIME_EEXIPS)

   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_EWALD)

      if ( iamoeba == 1 )then
         call AM_NonBond_eval(natom,x,f,vir,xx,ipairs, &
                               evdw,eelt,epolar,&
                               enb14,ee14,diprms,dipiter)
      else
! M-WJ
!         if ( induced == 1 )then
         if ( induced > 0 ) then
!
              call handle_induced(x,natom,ix(i04),ix(i06), &
                 xx(l15),cn1,cn2,asol,bsol, &
!                eelt,epolar,f,xx,ix,ipairs,xx(lpol), &
! Modified by WJM, YD, mjhsieh
                 eelt,epolar,f,xx,ix,ipairs,xx(lpol),xx(ldf),xx(lpol2), &
                 xx(l45),virvsene,ix(i02),ibgwat,nres, &
                 aveper,aveind,avetot,emtot,diprms,dipiter,dipole_temp,dt, &
                 ntb,cn3,cn4,cn5)
         else
            if ( ilrt /= 0 ) then
               ! Modifications for computing interaction energy
               ! according to the Linear Response Theory, LIE module
               if ( do_lrt ) then
                  ! call with molecule charges set to zero
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       crg_m0,cn1,cn2,asol,bsol,energy_m0,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),xx(ldf), &
                       xx(lpol2), .false. , cn3, cn4, cn5 )
                  ! call with water charges set to zero
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       crg_w0,cn1,cn2,asol,bsol,energy_w0,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),xx(ldf), &
                       xx(lpol2), .false. , cn3, cn4, cn5 )
                  ! call with full charges but no vdw interaction between solute and solvent
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       xx(l15),cn1_lrt,cn2_lrt,asol,bsol,eelt,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),xx(ldf), &
                       xx(lpol2), .false. , cn3, cn4, cn5 )
                  energy_vdw0 = evdw
                  call lrt_solute_sasa(x,natom, xx(l165))
               end if
               ! call normal_ewald force this will overwrite everything 
               ! computed above except energy_m0 and energy_w0
               call ewald_force(x,natom,ix(i04),ix(i06), &
                    xx(l15),cn1,cn2,asol,bsol,eelt,epolar, &
                    f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),xx(ldf), &
                    xx(lpol2), .false. , cn3, cn4, cn5 )
               energy_vdw0 = evdw - energy_vdw0
               ! count call to ltr, maybe calculate Eee and print it
               call ee_linear_response(eelt, master)
            else ! just call ewald_force normally
             call ewald_force(x,natom,ix(i04),ix(i06), &
                xx(l15),cn1,cn2,asol,bsol,eelt,epolar, &
                f,xx,ix,ipairs,xx(l45),virvsene,xx(lpol),xx(ldf), &
                 xx(lpol2), .false. , cn3, cn4, cn5 )
  
            end if 
         end if 
      end if

      call timer_stop(TIME_EWALD)

#ifdef MPI
      if(mytaskid == 0)then
#endif
         pot%vdw      = evdw
         pot%elec     = eelt
         pot%hbond    = ehb  !whereis ehb?
#ifdef MPI
      else
         ! energies have already been reduced to the master
         ! node in ewald_force, so here we zero out elements
         ! on non-master nodes:
         pot%vdw      = 0.d0
         pot%elec     = 0.d0
         pot%hbond    = 0.d0

      end if
#endif

      if( ips > 0 )then
         pot%vdw   = pot%vdw   + evdwex
         pot%elec  = pot%elec  + eelex
      endif

   end if  

   call timer_stop(TIME_NONBON)

   ! ----------------------------------------------------------------
   ! Calculate the other contributions
   ! ----------------------------------------------------------------

   !     -- when igb==10, all nonbonds are done in routine pb_force, and
   !                      all nonpolar interactions are done in np_force:
   !
   !     -- HG put this part here such that "outflag" is known from a call
   !        of pb_force; outflag is needed in the "bond" routine in the case
   !        of ifcap == 2,5 (i.e., ivcap == 1,5)

#ifdef MPI
   if(mytaskid == 0)then
#endif
      if( igb == 10 .or. ipb /= 0 ) then
         call timer_start(TIME_PBFORCE)
         call pbtimer_init
         call pb_force(natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),ix(i10), &
                 cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
         if ( pbgrid ) pbgrid = .false.
         if ( pbinit ) pbinit = .false.
         pot%vdw  = evdw
         pot%elec = eelt
         pot%pb   = epol
         call timer_stop(TIME_PBFORCE)

         call timer_start(TIME_NPFORCE)
         esurf = 0.0d0; edisp = 0.0d0
         if ( ifcap == 0  .and. npopt /= 0 ) &
            call np_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),cn1,cn2,x,f,esurf,edisp);
         if ( pbprint ) pbprint = .false.
         pot%surf = esurf
         pot%disp = edisp
         call pbtimer_summary
         call timer_stop(TIME_NPFORCE)

      end if 

#ifdef MPI
   end if
#endif

!  +---------------------------------------------------------------+
!  |  Bonds with H                                                 |
!  +---------------------------------------------------------------+

   call timer_start(TIME_BOND)

   ! initialize bond virial
   if(ipimd>0) bnd_vir = zero

#ifdef MPI /* SOFT CORE */
   ! zero only once, sc bond energy is sum of H and non-H terms
   sc_ener(1) = 0.0d0
#endif

   if( ntf < 2 ) then

      ebdev = 0.d0
      call bond(nbonh,ix(iibh),ix(ijbh),ix(iicbh),x,xx,ix,f,ene(6))
      pot%bond = pot%bond + ene(6)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesb
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Bonds without H                                              |
!  +---------------------------------------------------------------+

   if( ntf < 3 ) then

      call bond(nbona+nbper,ix(iiba),ix(ijba),ix(iicba),x,xx,ix,f,ene(7))
      pot%bond = pot%bond + ene(7)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesb
      endif
#  endif
#endif
      if (nbonh+nbona > 0) ebdev = sqrt( ebdev/(nbonh+nbona) )
   end if

!  +---------------------------------------------------------------+
!  |  Angles with H                                                |
!  +---------------------------------------------------------------+

   if( ntf < 4 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(2) = 0.0d0
#endif

      eadev = 0.d0

      call angl(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),x,xx,ix,f,ene(8))
      pot%angle = pot%angle + ene(8)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesa
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Angles without H                                             |
!  +---------------------------------------------------------------+

   if( ntf < 5 ) then

      call angl(ntheta+ngper,ix(i32),ix(i34),ix(i36),ix(i38),x,xx,ix,f, &
           ene(9)) 
      pot%angle = pot%angle + ene(9)

#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesa
      endif
#  endif
#endif
      if (ntheth+ntheta > 0) eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
   end if

!  +--------------------------------------------------------------------+
!  | AMD calculate dihedral energy first to estimate dihedral weight,   |
!  | then use it in the regular ephi function. Added by Romelia Salomon |
!  | Fix me later, AMD dihedral weight does NOT support AMOEBA          |
!  +--------------------------------------------------------------------+

   if(iamd .gt. 1)then
!    |  Dihedrals with H 
     if( ntf < 6 ) then
       call ephi_ene_amd(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
            x,amd_dih_noH)
     endif
!    |  Dihedrals without H 
     if( ntf < 7 ) then
       call ephi_ene_amd(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
            x,amd_dih_H)
     endif
#ifdef MPI 
     temp_amd_totdih = amd_dih_noH + amd_dih_H
#ifdef USE_MPI_IN_PLACE
     call mpi_allreduce(MPI_IN_PLACE,temp_amd_totdih,1,MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
     amd_totdih = temp_amd_totdih
#else
     call mpi_allreduce(temp_amd_totdih, amd_totdih, &
                  1, MPI_DOUBLE_PRECISION, &
                  mpi_sum, commsander, ierr)
#endif
#else
     amd_totdih = amd_dih_noH + amd_dih_H
#endif
     call calculate_amd_dih_weights(amd_totdih,temp0)
   endif

!  +---------------------------------------------------------------+
!  |  Dihedrals with H                                             |
!  +---------------------------------------------------------------+
   ! initialize 14 nb energy virial
   if(ipimd>0) e14vir = zero

   if( ntf < 6 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(3) = 0.0d0
#endif

      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
           xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(10),ene(11),ene(12),xx(l190))

      pot%dihedral = pot%dihedral + ene(10) ! Combine contributions from dihedrals with H
      pot%vdw_14   = pot%vdw_14   + ene(11) ! Combine 1-4 vdw contributions from dihedrals with H
      pot%elec_14  = pot%elec_14  + ene(12) ! Combine 1-4 elec contributions from dihedrals with H


#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesd
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Dihedrals without H                                          |
!  +---------------------------------------------------------------+

   if( ntf < 7 ) then

      call ephi(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
           xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(13),ene(14),ene(15),xx(l190))

      pot%dihedral = pot%dihedral + ene(13) ! Combine contributions from dihedrals without H
      pot%vdw_14   = pot%vdw_14   + ene(14) ! Combine 1-4 vdw contributions from dihedrals without H
      pot%elec_14  = pot%elec_14  + ene(15) ! Combine 1-4 elec contributions from dihedrals without H

#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesd
      endif
#  endif
#endif

!  +---------------------------------------------------------------+
!  |  CHARMM IMPROPERS IF CHARMM FORCEFIELD IS IN USE              |
!  +---------------------------------------------------------------+
!  Note: CHARMM does not distinguish between impropers with and 
!        without hydrogen hence it is not possible to strictly
!        conform to the ntf options of sander. Here CHARMM impropers
!        are calculated as long as ntf < 7 - so CHARMM impropers are
!        essentially considered to be in the same set as dihedrals
!        NOT involving hydrogen.
     if (charmm_active) then
       call charmm_calc_urey_bradley(x,pot%angle_ub,f)
       call charmm_calc_impropers(x,pot%imp,f)
       call charmm_calc_cmap(x,pot%cmap,f)
     end if

! --- END CHARMM IMPROPERS IF REQUIRED ---
     if (cmap_active) then
        call calc_cmap(x,pot%cmap,f)
     end if

   end if !ntf < 7

   if(iamoeba==1) then
      call AM_VAL_eval(x,f,vir,ene(6),ene(8),ene(10))
                              !ebond, eangle,etors
      pot%bond     = pot%bond     + ene(6)
      pot%angle    = pot%angle    + ene(8)
      pot%dihedral = pot%dihedral + ene(10)

   end if

   call timer_stop(TIME_BOND)

   ! --- calculate the EMAP constraint energy ---

   if(temap) then   ! ntr=1 (positional restraints)
       call emapforce(natom,enemap,xx(lmass),x,f )
       pot%emap = enemap
   end if

   ! --- calculate the position constraint energy ---

   if(natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
       call xconst(natc,entr,ix(icnstrgp),x,f,xx(lcrdr),xx(l60),natom )
       pot%constraint = entr
   end if

   if ( itgtmd==1 .and. (nattgtfit > 0 .or. nattgtrms > 0) ) then

      ! Calculate rmsd for targeted md (or minimization) if requested.
      ! All nodes do rms fit, could just be master then broadcast.
      ! All nodes need all coordinates for this.

      call rmsfit(xx(lcrdr),x,xx(lmass),ix(itgtfitgp),  &
                  ix(itgtrmsgp),rmsdvalue,nattgtrms,nattgtfit,rmsok)

      if (.not.rmsok) then
         if (master) write (6,*) 'Fatal error: Error calculating RMSD!'
         call mexit(6, 1)
      end if

      call xtgtmd(entr,ix(itgtrmsgp),x,f,xx(lcrdr),xx(lmass),tgtrmsd,tgtmdfrc,rmsdvalue,nattgtrms)
      pot%constraint = entr
   else if(itgtmd == 2) then
      call mtmdcall(entr,xx(lmtmd01),ix(imtmd02),x,f,ih(m04),ih(m02),ix(i02),&
                    ih(m06),xx(lmass),natom,nres,'CALC')
      pot%constraint = entr
   end if

   if(ifcap == 1 .or. ifcap == 2) then
      call capwat(natom,x,f,ecap)
      pot%constraint = pot%constraint + ecap
   else if(ifcap == 3) then
      write(6,*) 'No energy expression for spherical boundary known yet'
      call mexit(6,1)
   else if(ifcap == 4) then
      write(6,*) 'No energy expression for orthorhombic boundary known yet'
      call mexit(6,1)
   end if
   ! No energy expression for ifcap == 5 given because only
   !    one step of minimization is allowed with this.

   !  (this seems very weird: we have already done an allreduce on molvir
   !  in ewald_force(); this just collects it on processor 0 (with zeroes
   !  on all slave nodes), then later does an allreduce...)

   if( mytaskid == 0 .and. iamoeba == 0 ) then
      vir(1) = vir(1)+0.5d0*molvir(1,1)
      vir(2) = vir(2)+0.5d0*molvir(2,2)
      vir(3) = vir(3)+0.5d0*molvir(3,3)
   end if

   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then
      ener%virvsene    = virvsene
      ener%diprms      = diprms
      ener%dipiter     = dipiter
      ener%dipole_temp = dipole_temp
   end if

   !     ---- get the noesy volume penalty energy: ------

   pot%noe = 0.d0
   if( iredir(4) /= 0 ) then
      call timer_start(TIME_NOE)
      call noecalc(x,f,xx,ix)
      call timer_stop(TIME_NOE)
   end if
   ! Do we need a pot%noe here?  mjw TODO

   !     -- when igb!=0 and igb!=10, all nonbonds are done in routine egb:

   esurf = 0.d0
   if( igb /= 0 .and. igb /= 10 .and. ipb == 0 ) then
      call timer_start(TIME_EGB)
      call egb( x,f,rborn,fs,reff,onereff,xx(l15),ix(i04),ix(i06), &
            ix(i08),ix(i10),xx(l190), &
            cut,ntypes,natom,natbel,epol,eelt,evdw, &
            esurf,dvdl,xx(l165),ix(i82),xx(l170),xx(l175),xx(l180), &
            xx(l185), xx(l186),xx(l187),xx(l188),xx(l189),ncopy, &
            xx(l2402),xx(l2403),xx(l2404)  )   ! last three are GB8 arrays

      pot%vdw  = evdw
      pot%elec = eelt
      pot%gb   = epol
      pot%surf = esurf
      pot%dvdl = dvdl
      call timer_stop(TIME_EGB)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
        pot%les = pot%les + elesp
      endif
#  endif
#endif

   end if

#ifdef RISMSANDER
   
   if(rismprm%irism == 1) then
      call timer_start(TIME_RISM)
      call rism_force(x,f,erism,irespa)
      pot%rism = erism
      call timer_stop(TIME_RISM)
   endif
#endif

#ifdef APBS
! APBS forces
      if( mdin_apbs ) then
         if (igb /= 6) then
            write(6, '(a)') '&apbs keyword requires igb=6.'
            call mexit(6,1)
         end if
         call timer_start(TIME_PBFORCE)
! in: coords, radii, charges
! out: updated forces (via apbs_params) and solvation energy (pol + apolar)
         if (sp_apbs) then
            call apbs_spenergy(natom, x, f, eelt, enpol)
         else
            call apbs_force(natom, x, f, pot%vdw, eelt, enpol)
         end if
         pot%pb   = eelt 
         pot%surf = enpol
         call timer_stop(TIME_PBFORCE)

      end if  ! ( mdin_apbs )
#endif /* APBS */

   if( master ) then
      !  These parts of the NMR energies are not parallelized, so only
      !  are done on the master node:
      eshf = 0.d0
      epcshf = 0.d0
      ealign = 0.d0
      ecsa = 0.d0
      if (iredir(5) /= 0) call cshf(natom,x,f)
      if (iredir(7) /= 0) call pcshift(natom,x,f)
      if (iredir(9) /= 0) call csa1(natom,x,f)
      if (iredir(8) /= 0) call align1(natom,x,f,xx(lmass))
   end if

   ! additional force due to charge relocation
   if ( ifcr /= 0 ) then
      call cr_calc_force( f )
   end if

#ifdef MPI

   call timer_barrier( commsander )
   call timer_start(TIME_COLLFRC)

   !     add force, ene, vir, copies from all nodes
   !            also add up newbalance for nonperiodic.

   ! Remember to work on the local instance of the
   ! potential energy array, i.e. pot and NOT the global one,
   ! i.e. ener%pot

   call fdist(f,xx(lfrctmp),pot,vir,newbalance)


   call timer_stop(TIME_COLLFRC)

#endif

   ! ---- at this point, the parallel part of the force calculation is
   !      finished, and the forces have been distributed to their needed
   !      locations.  All forces below here are computed redundantly on
   !      all processors, and added into the force vector.  Hence, below
   !      is the place to put any component of the force calculation that
   !      has not (yet) been parallelized.

   ! Calculate the NMR restraint energy contributions, if requested.
   ! (Even though this is not parallelized, it needs to be run on all
   ! threads, since this code is needed for weight changes as well as
   ! for NMR restraint energy analysis.  The whole thing could stand a
   ! major re-write....)

   if (nmropt > 0) &
      call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
         devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb, &
         xx(lnmr01),ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
         cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
         nimprp,nphb,natom,natom,ntypes,nres, &
         rad,wel,radhb,welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'CALC')
#ifdef MPI
   call mpi_reduce(enoe,pot%noe,1,MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
   enoe = pot%noe ! so all processors now have the full enoe value
#else
   pot%noe = enoe
#endif

#ifdef DSSP
   if( idssp > 0 ) then
      call fdssp( natom,x,f,edssp )
   else
      edssp = 0.d0
   end if
#endif

   !     ----- CALCULATE TOTAL ENERGY AND GROUP THE COMPONENTS -----

#ifndef LES
   if( igb == 0 .and. ipb == 0 ) then
      pot%vdw_14   = pot%vdw_14   + enb14
      pot%elec_14  = pot%elec_14  + ee14

   endif

#endif

   pot%constraint = pot%constraint+ eshf+epcshf+ pot%noe + &
                    sum(enmr(1:6))+ealign+ecsa



#ifdef DSSP
   pot%constraint = pot%constraint + edssp
#endif

   pot%polar = epolar

   pot%tot=      pot%vdw        + &
                 pot%elec       + &
                 pot%gb         + &
                 pot%pb         + &
                 pot%bond       + &
                 pot%angle      + &
                 pot%dihedral   + &
                 pot%vdw_14     + &
                 pot%elec_14    + &
                 pot%hbond      + &
                 pot%constraint + &
                 pot%emap       + &
                 pot%rism       + &
                 pot%ct

   pot%tot = pot%tot + pot%polar + pot%surf + pot%scf + pot%disp

   !Charmm related
   pot%tot = pot%tot + pot%angle_ub + pot%imp + pot%cmap 

   !The handover
   ener%pot = pot

!  +---------------------------------------------------------------+
!  | AMD calculate total potential energy weight, then apply it to |
!  | all the force elements f=f*fwgt. Added by Romelia Salomon     |
!  +---------------------------------------------------------------+

   if(iamd .gt. 0)then
     call calculate_amd_total_weights(natom,pot%tot,amd_totdih,f,temp0)
   end if


   ener%aveper = aveper
   ener%aveind = aveind
   ener%avetot = avetot



   
   ! This is now historical; MJW Feb 2010
   !
   !    Here is a summary of how the ene array is used.  For parallel runs,
   !    these values get summed then rebroadcast to all nodes (via
   !    mpi_allreduce).

   !    ene(1):      total energy
   !    ene(2):      van der Waals
   !    ene(3):      electrostatic energy
   !    ene(4):      10-12 (hb) energy, or GB energy when igb.gt.0
   !    ene(5):      bond energy
   !    ene(6):      angle energy
   !    ene(7):      torsion angle energy
   !    ene(8):      1-4 nonbonds
   !    ene(9):      1-4 electrostatics
   !    ene(10):     constraint energy
   !    ene(11-19):  used as scratch, but not needed further below
   !    ene(20):     position constraint energy + cap energy
   !    ene(21):     charging free energy result
   !    ene(22):     noe volume penalty
   !    ene(23):     surface-area dependent energy, or cavity energy
   !    ene(24):     potential energy for a subset of atoms
   !    ene(25):     SCF Energy when doing QMMM
   !    ene(26):     implicit solvation dispersion energy


#ifdef PUPIL_SUPPORT
   !*****************************************************
   !     Closing the qmmm structure consideration
   !*****************************************************
!  Adding the quantum forces from last QM calculation
   do iPup=1,pupqatoms
     bs1 = (abs(pupqlist(iPup))-1)*3
     do jPup=1,3
       bs2    = bs1    + jPup
       f(bs2) = f(bs2) + qfpup(bs2)
     enddo
   enddo


!  Disconnecting qmmmm interactions
   qmmm_nml%ifqnt = .false.

#endif

   ! ----ADD X-RAY TARGET FUNCTION AND GRADIENT
   call cns_xref_run(natom,ih(m04), x,f,ener)

#ifdef _XRAY
   ! ---- BUILT-IN X-RAY TARGET FUNCTION AND GRADIENT
   if (xray_active) call xray_get_derivative(x,f)
#endif

   !     if freezemol has been set, zero out all of the forces for the
   !     real atoms; (no longer necessary to set ibelly).
   if( ifreeze > 0 ) then
      do i=1,3*natom
         f(i) = 0.d0
      end do
   end if

   !     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----
   if (belly) call bellyf(natom,ix(ibellygp),f)

!  +---------------------------------------------------------------+
!  |  Interface to EVB                                             |
!  +---------------------------------------------------------------+

#if defined(MPI)
#ifdef LES
   if( nbead > 0 ) then
      call mpi_allreduce ( nrg_all, nrg_bead, nbead, MPI_DOUBLE_PRECISION &
                         , MPI_SUM, commsander, ierr )
      nrg_all(:) = nrg_bead(:)
   end if
#endif
   if( ievb /= 0 ) call evb_ntrfc ( x, f, ener, xx(lmass), ix, ipairs, vel0_nrg_sum )
#endif /* MPI */

   if( ipimd>0 ) then
#ifdef LES
      call pimd_part_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
#ifdef MPI
      if( ievb /= 0 ) then
         nrg_frc(3)= vel0_nrg_sum
         nrg_frc(2)= equal_part + Epot_deriv
         nrg_frc(1)= nrg_frc(3) + nrg_frc(2)
         dlnQ_dl = dvdl
      endif
#endif
#else
      ener = ener/nbead 
      f(1:natom*3) = f(1:natom*3)/nbead
      vir(1:3) = vir(1:3) /nbead
      atvir = atvir/nbead
      e14vir = e14vir/nbead
      bnd_vir = bnd_vir/nbead
      call pimd_full_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
# ifdef MPI
      if (master) call mpi_reduce(ener,totener,state_rec_len,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,commmaster,ierr)
# endif
#endif
      ! Pass dvdl = dV/dl for TI w.r.t. mass.
      if (itimass > 0) ener%pot%dvdl = dvdl
   end if
#ifdef MPI

     !CARLOS: ONLY MPI SUPPORTS NEB (MULTISANDER)
     if(ineb>0) then
        if(sanderrank.eq.0) then  !only masters do NEB
           call full_neb_forces( xx(lmass), x, f, ener%pot%tot, ix(itgtfitgp),ix(itgtrmsgp))
        endif

        ! now master will broadcast the neb forces for the rmsgp atoms
        ! all nodes will add this to the current force total
        ! if we wanted all nodes to calculate the neb forces, they
        ! would all need access to the neighbor bead coordinates, which is
        ! probably more expensive than having the master send out the modified
        ! forces

        ! master broadcasts the force update for the rmsgp atoms only

        call mpi_bcast(neb_force,nattgtrms*3,MPI_DOUBLE_PRECISION,0,commsander,ierr)

        do i=1,nattgtrms
           j3 = 3*(i - 1) !pointer into the packed neb force array
           j=ix(itgtrmsgp+i-1) !actual atom # for atom in rms group
           i3 = 3*(j - 1) !pointer into real force array

           f(i3+1)=f(i3+1)+neb_force(j3+1)
           f(i3+2)=f(i3+2)+neb_force(j3+2)
           f(i3+3)=f(i3+3)+neb_force(j3+3)
        enddo

! CARLOS: what is this doing? ener(27) is neb energy
! looks like it reduces entire energy array
! needs MUCH better documentation on details (such as 28)

        if(sanderrank.eq.0) then
           call mpi_reduce(ener,totener,state_rec_len,MPI_DOUBLE_PRECISION,MPI_SUM,0,commmaster,ierr)
        end if

   end if !ineb>0
#endif /*MPI*/


   if (charmm_active) then
     if ( do_charmm_dump_gold == 1 ) then
       call charmm_dump_gold(f,natom,ener)
     endif 
   end if


#ifdef MPI

   if(icfe == 0) then
      if(idecomp == 1 .or. idecomp == 2) then
         call collect_dec2(nres)
      end if
          
      if(idecomp >= 3) then
         call collect_dec2(npdec*npdec)
      end if
   end if

#endif

   call timer_stop(TIME_FORCE)
   call trace_exit( 'force' )
   return


end subroutine force


