! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!----------------------------------------------------------------------
!                   EWALD_FORCE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine ewald_force(crd,numatoms,iac,ico,charge, &
      cn1,cn2,asol,bsol,eelt,epol,frc,x,ix,ipairs, &
!      xr,virvsene,pol,qm_pot_only &
! Modified by WJM, YD
      xr,virvsene,pol,dampfactor,pol2,qm_pot_only &
! by mjhsieh
      , cn3, cn4, cn5 &
!! 
      )
   use trace
   use ew_recip
   use ew_dipole_recip
   use stack
   use decomp, only: decpr
   use nblist, only: recip, cutoffnb,nbfilter, volume, &
                     maxnblst,nucgrd1,nucgrd2,nucgrd3, &
                     adjust_imagcrds, map_coords, &
                     nvdwcls
   use constants, only : INV_AMBER_ELECTROSTATIC2, INV_AMBER_ELECTROSTATIC, &
                         zero, AMBER_ELECTROSTATIC
   use qmmm_module, only : qmmm_struct, qm2_struct, qmmm_nml, qmewald, &
                           qmmm_scratch
   use parms, only : one_scee, one_scnb
   use charmm_mod, only : charmm_active, charmm_cn114, charmm_cn214
   use file_io_dat
#ifdef LES
   use pimd_vars, only : ipimd, nrg_all,nbead
   use part_pimd_vars, only: frcx_copy,frcx_temp
#endif
   use nbips, only :aipspbc, ips,teips,tvips,teaips,tvaips,virexips,wnb

   implicit none
   character(kind=1,len=11) :: routine="ewald_force"
#  include "extra.h"
#  include "ew_cntrl.h"
#  include "def_time.h"
#  include "ew_time.h"
#  include "ew_pme_recip.h"
#  include "ew_mpole.h"
#  include "ew_erfc_spline.h"
#  include "memory.h"
#  include "md.h"
#  include "extra_pts.h"
#  include "box.h"
#  include "les.h"
#  include "ew_frc.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer numatoms,iac(*),ico(*)
   _REAL_ crd(3,*),charge(*),cn1(*),cn2(*),asol(*), &
         bsol(*),eelt,epol,frc(3,*),xr(3,*),virvsene
!   _REAL_ pol(*)
! Modified by WJM, YD
   _REAL_ pol(*), dampfactor(*), pol2(*)
! by mjhsieh
   _REAL_ cn3(*), cn4(*), cn5(*)
!!
   logical, intent(in) :: qm_pot_only ! Flag for whether we do a regular 
                                      ! full pme (gradients as well) or just 
                                      ! get the potential for qm scf procedure.
#ifdef MPI
   integer l_temp,ndel
   _REAL_ rl_temp(BC_EW_COMM3)
! JMS: xdist doesn't consider mxvar from nmr.h, so we won't bother with it here
   _REAL_ x_temp(numatoms*3)
#endif

   integer nstart,ntop
   _REAL_ x(*)
   integer ix(*),ipairs(*)
   integer mm_no,qm_no, lnk_no
   _REAL_ forcemod(3)
   _REAL_ boltz2
#ifdef LES
   integer icopy
#endif

   !     ARGUMENTS:
   !     CRD (INPUT) is the array holding atomic coords.
   !     NUMATOMS (INPUT) is the number of atoms
   !     IAC, ICO and NTYPES (INPUT) are used to look up vdw and hbond
   !       interactions.  ICO is an NTYPES by NTYPES array giving lookup
   !       into VDW and HBOND coefficients. IAC(i) is the atom type of
   !       atom i. The coefficients for 6-12 interactions for atoms i
   !       and j are indexed by ICO(iac(i),iac(j)).  In practice ICO
   !       is unrolled to a 1 dimensional array.  They are needed here
   !       to split nonbond list into vdw,hbonds (see pack_nb_list())
   !     CHARGE (INPUT) is the array of atomic partial charges.
   !     CN1 and CN2 (INPUT) are the VDW coefficients
   !     ASOL and BSOL (INPUT) are the HBOND coefficients
   !     EELT,EVDW,EHB (OUTPUT) are the energies computed and
   !     FRC (OUTPUT) is the force array
   !     X and IX are real and integer arrays which comprise the total
   !       "dynamic" memory  for amber; coords,forces,bond lists etc are
   !       accessed as offsets in them.
   
   integer indtop, i,j,commsander_mytaskid,commsander_numtasks, qm_temp_count

   call trace_enter( 'ewald_force' )
#ifdef MPI
   commsander_mytaskid = mytaskid
   commsander_numtasks = numtasks
#else
   commsander_mytaskid = 0
   commsander_numtasks = 1
#endif
 call timer_start(TIME_EWFSTRT)
   
   ! Clear all energies and virials plus our scratch force array
   
   eelt = 0.d0
   eer = 0.d0
   eed = 0.d0
   eea = 0.d0
   ees = 0.d0
   epol = 0.d0
   epold = 0.d0
   epola = 0.d0
   epols = 0.d0
   evdw = 0.d0
   ehb = 0.d0
   eedvir = 0.d0
   evdwr = 0.d0
   !     if( numextra.gt.0 ) then
   ee14 = 0.d0
   enb14 = 0.d0
   epol14 = 0.d0
   !     endif
   eeles=0.d0
   rec_vir(1:3,1:3) = 0.d0
   rec_vird(1:3,1:3) = 0.d0
   dir_vir(1:3,1:3) = 0.d0
   adj_vir(1:3,1:3) = 0.d0
   self_vir(1:3,1:3) = 0.d0
   les_vir(1:3,1:3) = 0.d0
   e14vir(1:3,1:3) = 0.d0
   framevir(1:3,1:3) = 0.d0
   virvsene = 0.d0
   if(ips>0)then
     adj_vir=virexips
   endif
! M-WJ
!   if ( mpoltype == 1 )then
   if ( mpoltype > 0 ) then
!
      call zero_array(frc,3*numatoms) ! note that ips and ipol are inconsistent
      call zero_array(x(lfield),3*numatoms)
   end if

#ifdef LES
   if(ipimd>0)   frcx_copy(1:3,1:nbead) = 0.d0
#endif

   frcx(1:3) = 0.d0
   
   call timer_stop(TIME_EWFSTRT)
   if( numextra > 0 ) call local_to_global(crd,x,ix)
#ifdef MPI
! M-WJ
!   if ( induced == 1 ) then
   if ( induced > 0 ) then
!
      call timer_start(TIME_DISTDIP)
      call xdist(x(linddip),x_temp,numatoms)
      call timer_stop(TIME_DISTDIP)
   end if
#endif

   
   if( use_pme /= 0 .and. mod(irespa,nrespa) == 0) then
      
      !--------------------------------------------------------
      ! SELF ENERGY
      !--------------------------------------------------------
      
      call timer_start(TIME_SELF)
      if ( master ) call self(charge,numatoms,ees,ew_coeff,volume,self_vir,frc)
      call timer_stop(TIME_SELF)

      !--------------------------------------------------------
      ! RECIPROCAL ENERGY
      !--------------------------------------------------------


      call timer_start(TIME_REC)
      if (.not. qmmm_nml%ifqnt) then
#ifdef MPI
        if(i_do_recip)then
           call mpi_comm_size(recip_comm,numtasks,ierr)
           call mpi_comm_rank(recip_comm,mytaskid,ierr)
           master = mytaskid.eq.0
#endif

           if ( ew_type == 0 )then
              call do_pme_recip(mpoltype,numatoms,crd,charge,frc,     &
                    X(linddip),x(lfield),x(lprefac1),x(lprefac2),    &
                    X(lprefac3),x(lfftable),qm_pot_only )
           else if ( ew_type == 1 )then
              call map_coords(crd,numatoms,recip)
              if ( mpoltype == 0 )then
                 call recip_reg(numatoms,charge,eer,rec_vir, &
                       mlimit,volume,recip,frc, &
                       ew_coeff,maxexp,commsander_mytaskid,commsander_numtasks)
! M-WJ
!              else if ( mpoltype == 1 )then
              else if ( mpoltype > 0 )then
!
                 call recip_reg_dipole(numatoms,charge,eer,rec_vir, &
                       mlimit,volume,recip,frc, &
                       ew_coeff,maxexp,x(linddip),x(lfield), &
                       commsander_mytaskid,commsander_numtasks,forctime, &
                       iscsum,stack_time)
              end if
           else
              ASSERT( .false. )  ! ew_type input validation occurs in mdread.f
           end if
        
           ! Get the dipole specific part of virial, which depends on recip 
           ! efield.  Note this logic relies on the recip being done first; 
           ! i.e. total efield is the recip efield at this point
         
! M-WJ
!           if ( mpoltype == 1 )then
           if ( mpoltype > 0 ) then
!
              call dip_field_corr(x(lfield),x(linddip),rec_vir, &
                    numatoms, commsander_mytaskid,commsander_numtasks)
            
              ! Get the energy and field contributions from the self-images 
              ! of dipoles
            
              call self_dipole(numatoms,epols,ew_coeff,x(lfield), &
                    x(linddip),commsander_mytaskid,commsander_numtasks,frc)
           end if
         
           ! Scaleup forces, field for respa
         
           if ( nrespa > 1 ) &
                 call respa_scale(numatoms,frc,x(lfield),nrespa)
#ifdef MPI
           numtasks = commsander_numtasks
           mytaskid = commsander_mytaskid
           master = mytaskid.eq.0
        end if         !loop for ido_recip
#endif
      else !if .not ifqnt
        !--- QUANTUM PME ---
        if (qmmm_nml%qm_pme .and. .not. qm_pot_only) then
          !copy the mulliken charges out of scf_mchg and into main charge array.
          !all threads do this - 1,nquant_nlink
          do qm_temp_count = 1, qmmm_struct%nquant_nlink
             charge(qmmm_struct%iqmatoms(qm_temp_count)) = &
                        qm2_struct%scf_mchg(qm_temp_count)*AMBER_ELECTROSTATIC
          end do
        end if
#ifdef MPI
        if(i_do_recip)then
           call mpi_comm_size(recip_comm,numtasks,ierr)
           call mpi_comm_rank(recip_comm,mytaskid,ierr)
#endif
           if (qm_pot_only) then
             call adj_mm_link_pair_crd(qmmm_nml,qmmm_struct, crd) !Replace the MM atom's coordinates 
                                          !      with the link atom.

             call do_pme_recip(mpoltype,numatoms,crd,charge, &
                   qmmm_scratch%qm_real_scratch,  &
                   X(linddip),x(lfield),x(lprefac1),x(lprefac2),    &
                   X(lprefac3),x(lfftable),qm_pot_only )

             call rst_mm_link_pair_crd(qmmm_nml,qmmm_stuct, crd)
             !-- store the reciprocal energy so that we can subtract off 
             !   the QM-QM reciprocal energy
             !   when we next come through.
             qmewald%mm_recip_e=eer
             call timer_stop(TIME_REC)
  
             ! If qm_pot_only is true then do not do any more stuff here.
  
             return
           end if

           qmmm_scratch%qm_real_scratch(1:3*natom) = zero
           call adj_mm_link_pair_crd(qmmm_nml,qmmm_struct, crd) !Replace the MM atom's coordinates 
                                          !      with the link atom.
           call do_pme_recip(mpoltype,numatoms,crd,charge, &
                 qmmm_scratch%qm_real_scratch,  &
                 X(linddip),x(lfield),x(lprefac1),x(lprefac2),    &
                 X(lprefac3),x(lfftable),qm_pot_only )


           call rst_mm_link_pair_crd(qmmm_nml,qmmm_stuct, crd)
           !We need to redistribute the force on link atoms
           do i=1,qmmm_struct%nlink
             mm_no = qmmm_struct%link_pairs(1,i)  !location of atom in x array
             lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom 
                                                  !      bound to link atom
             qm_no = qmmm_struct%iqmatoms(lnk_no)
             !Note this routine uses the flink in the form -flink.
             call distribute_lnk_f(forcemod, &
                         qmmm_scratch%qm_real_scratch(3*mm_no-2), &
                         crd(1,mm_no), crd(1,qm_no), qmmm_nml%lnk_dis)

             frc(1:3,mm_no) = frc(1:3,mm_no) + forcemod(1:3)
             frc(1:3,qm_no) = frc(1:3,qm_no) &
                   + qmmm_scratch%qm_real_scratch(3*mm_no-2:3*mm_no) &
                   - forcemod(1:3)

             !Zero out the link atom force in the scratch array.
             qmmm_scratch%qm_real_scratch(3*mm_no-2) = zero
             qmmm_scratch%qm_real_scratch(3*mm_no-1) = zero
             qmmm_scratch%qm_real_scratch(3*mm_no) = zero
           end do
 
           do i = 1, natom
             frc(1,i)=frc(1,i)+qmmm_scratch%qm_real_scratch(3*i-2)
             frc(2,i)=frc(2,i)+qmmm_scratch%qm_real_scratch(3*i-1)
             frc(3,i)=frc(3,i)+qmmm_scratch%qm_real_scratch(3*i)
           end do

           if (qmmm_nml%qm_pme) then
             ! Zero the mulliken charges in main charge array. 
             ! Probably belongs after do_pme_recip call above
             !     but here it avoids the need for an extra if statement.
             ! All threads do this.
 
             do qm_temp_count = 1, qmmm_struct%nquant_nlink
                charge(qmmm_struct%iqmatoms(qm_temp_count)) = zero
             end do
 
             ! Build a charge array that is natom long but contains only QM 
             ! scf charges.

             qmmm_scratch%qm_pme_scratch(1:natom) = zero
             qmmm_scratch%qm_real_scratch(1:3*natom) = zero
             do qm_temp_count = 1, qmmm_struct%nquant_nlink
                qmmm_scratch%qm_pme_scratch(qmmm_struct%iqmatoms(qm_temp_count)) &
                     = qm2_struct%scf_mchg(qm_temp_count)*AMBER_ELECTROSTATIC
             end do
 
             ! Calculate QM-QM PME forces to remove from 
             !   the full PME we did above.
             call adj_mm_link_pair_crd(qmmm_nml,qmmm_struct, crd)
             call do_pme_recip(mpoltype,numatoms,crd, &
                  qmmm_scratch%qm_pme_scratch, &
                  qmmm_scratch%qm_real_scratch, &
                  X(linddip),x(lfield),x(lprefac1),x(lprefac2),    &
                  X(lprefac3),x(lfftable),qm_pot_only )
 
             ! Now subtract out the QM-QM forces from the pme.
  
             ! Redistribute the force on link atoms
             call rst_mm_link_pair_crd(qmmm_nml,qmmm_struct, crd)
             do i=1,qmmm_struct%nlink
                mm_no = qmmm_struct%link_pairs(1,i)  !location in x array
                lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom 
                                                     !      bound to link atom
                qm_no = qmmm_struct%iqmatoms(lnk_no)
                ! Note this routine uses the flink in the form -flink.
                call distribute_lnk_f(forcemod, &
                     qmmm_scratch%qm_real_scratch(3*mm_no-2),crd(1,mm_no), &
                     crd(1,qm_no),qmmm_nml%lnk_dis)
                
               frc(1:3,mm_no) = frc(1:3,mm_no) - forcemod(1:3)
               frc(1:3,qm_no) = frc(1:3,qm_no) &
                    - qmmm_scratch%qm_real_scratch(3*mm_no-2:3*mm_no) &
                    + forcemod(1:3)
             end do
 
             do qm_temp_count = 1, qmmm_struct%nquant
               i = qmmm_struct%iqmatoms(qm_temp_count)
               frc(1,i)=frc(1,i)-qmmm_scratch%qm_real_scratch(3*i-2)
               frc(2,i)=frc(2,i)-qmmm_scratch%qm_real_scratch(3*i-1)
               frc(3,i)=frc(3,i)-qmmm_scratch%qm_real_scratch(3*i)
             end do

             ! Put back the energy from recip without the QM-QMrecip 
             ! contributions; these are already in ESCF.

             eer = qmewald%mm_recip_e
           end if !qmmm_nml%qm_pme
#ifdef MPI
           numtasks = commsander_numtasks
           mytaskid = commsander_mytaskid 
        end if         !loop for ido_recip
#endif
      end if !if .not. ifqnt
      call timer_stop(TIME_REC)
      
      !-------------------------------------------------------
      ! LONG RANGE DISPERSION CONTRIBUTION
      !-------------------------------------------------------

      call timer_start(TIME_SELF)
      if ( vdwmeth == 1)then
         if ( master ) then
               call vdw_correction(ico,ntypes,nvdwcls, &
               volume,evdwr,rec_vird,cn2,cutoffnb)
               ! -- ti decomp
               if(decpr .and. idecomp > 0) call vdwdec_correction(natom,ntypes,iac,ico,nvdwcls,cn2,volume,cutoffnb)
         end if
      end if
      call timer_stop(TIME_SELF)
      
   end if      ! (use_pme and respa check)
   if( teaips .or. tvaips ) then
      !--------------------------------------------------------
      ! Remaining IPS ENERGY
      !--------------------------------------------------------

      call timer_stop_start(TIME_EWALD,TIME_AIPS)
#ifdef MPI
        if(i_do_recip)then
           call mpi_comm_size(recip_comm,numtasks,ierr)
           call mpi_comm_rank(recip_comm,mytaskid,ierr)
           master = mytaskid.eq.0
#endif

            call aipspbc(evdwr,eer,numatoms,crd,charge,frcx,frc,rec_vir)

           ! Scaleup forces, field for respa
         
#ifdef MPI
           numtasks = commsander_numtasks
           mytaskid = commsander_mytaskid
           master = mytaskid.eq.0
        end if         !loop for ido_recip
#endif
      call timer_stop_start(TIME_AIPS,TIME_EWALD)
      
   end if      ! (do aips )
   
   if( ips>0 .and. vdwmeth==1 .and. master ) then

      !-------------------------------------------------------
      ! LONG RANGE DISPERSION CONTRIBUTION
      !-------------------------------------------------------
      
      call timer_start(TIME_SELF)
      call vdw_correction(ico,ntypes,nvdwcls, &
               volume,evdwr,rec_vird,cn2,cutoffnb)
      call timer_stop(TIME_SELF)

   end if

   !--------------------------------------------------------
   ! DIRECT PART OF EWALD PLUS VDW, HBOND
   !--------------------------------------------------------
   
   call timer_start(TIME_DIR)
#ifdef MPI
   if(i_do_direct)then
      call mpi_comm_size(direct_comm,numtasks,ierr)
      call mpi_comm_rank(direct_comm,mytaskid,ierr)
#endif

      !         Get fractional coordinates, and adjust the imaged cartesians
      !         for pre-imaged short range interactions

      call map_coords(crd,numatoms,recip)
      call adjust_imagcrds(crd,natom )

      indtop = nucgrd1*nucgrd2*nucgrd3
      call get_nb_energy(iac,ico,ntypes,charge, &
            cn1,cn2,asol,bsol,frc,numatoms, &
            indtop, &
            ipairs, &
            ew_coeff,eedtbdns,x(leed_cub),x(leed_lin), &
            maxnblst,eed,evdw,ehb,dir_vir,eedvir, &
            nbfilter,ee_type,eedmeth,dxdr, &
! M-WJ
!            pol,dipdamp,&
! Modified by WJM
            pol,dampfactor,pol2, &
!
!
! mjhsieh
            cn3, cn4, cn5, &
!
            epold,x(linddip),x(lfield),mpoltype)


#ifdef MPI
       numtasks = commsander_numtasks
       mytaskid = commsander_mytaskid
   end if
#endif
   call timer_stop(TIME_DIR)

   !--------------------------------------------------------
   ! ADJUST ENERGIES,FORCES FOR MASKED OUT PAIRS
   !--------------------------------------------------------
   if( use_pme /= 0 ) then
      call timer_start(TIME_ADJ)
      if ( mpoltype == 0 )then

         call nb_adjust(charge,eea,crd, &
               ix(imask1),ix(imask2),numadjst, &
               ew_coeff,eedtbdns, &
               x(leed_cub),x(leed_lin),frc,numatoms,adj_vir,cutoffnb, &
               ee_type,eedmeth)

! M-WJ
!      else if ( mpoltype == 1 )then
      else if ( mpoltype > 0 )then
!
         call nb_adjust_dipole(charge,eea,crd, &
               ix(imask1),ix(imask2),numadjst,ew_coeff,eedtbdns, &
               x(leed_cub),x(leed_lin),frc,numatoms,adj_vir,cutoffnb, &
               ee_type,eedmeth,x(linddip),x(lfield),epola)
      end if

      call timer_stop(TIME_ADJ)
   end if

#ifdef LES
   
   ! ADJUST ENERGIES,FORCES FOR LES INTRA-COPY PAIR SCALING FACTOR
   
   if (nlesadj > 0.and.ipimd==0) then
      call timer_start(TIME_LESADJ)
      call nb_adjust_les(charge,eeles,crd, frc, &
            numatoms,les_vir,use_pme)
      call timer_stop(TIME_LESADJ)
   end if
  
#endif /* LES */
   
#ifndef LES
   !--------------------------------------------------------
   ! do the 1-4 here & remove them from ephi
   !--------------------------------------------------------

   if ( mpoltype == 0 )then
      if (charmm_active) then
        call get_14_cg(charge,crd,frc,iac,ico,ntypes, &
              charmm_cn114,charmm_cn214, ee14,enb14,one_scee,one_scnb, &
              e14vir,ix,commsander_mytaskid,commsander_numtasks,eedmeth)
      ! mjhsieh: note that ff11cmap doesn't need to touch this
      else
        call get_14_cg(charge,crd,frc,iac,ico,ntypes, &
              cn1,cn2, ee14,enb14,one_scee,one_scnb, &
              e14vir,ix,commsander_mytaskid,commsander_numtasks,eedmeth)
      end if
! M-WJ
!   else if ( mpoltype == 1 ) then
   else if ( mpoltype > 0 ) then
!
      call get_14_dipole(charge,crd,frc, &
            x(linddip),x(lfield), &
            iac,ico,ntypes, &
! M-WJ
!            cn1,cn2, ee14,enb14,epol14,one_scee,one_scnb,e14vir,ix, &
!            cn1,cn2, pol,dipdamp, ee14,enb14,epol14,one_scee,one_scnb,e14vir,ix, &
! Modified by WJM, YD
            cn1,cn2, pol,dampfactor, pol2, &
            ee14,enb14,epol14,one_scee,one_scnb,e14vir,ix, &
!
            commsander_mytaskid,commsander_numtasks)
   end if

   ! Now transfer force and torque from extra points:
   if (numextra > 0) call orient_frc(crd,frc,framevir,ix)
#endif

   !--------------------------------------------------------
   ! ADJUST FORCES FOR net force (comes from recip work)
   !--------------------------------------------------------

   call timer_start(TIME_EWFRCADJ)
#ifdef MPI
   ndel=(numatoms+numtasks-1)/numtasks
   nstart = ndel*mytaskid+1
   ntop = nstart+ndel-1
   if(ntop > numatoms)ntop=numatoms
   call trace_mpi('mpi_allreduce',3,'MPI_DOUBLE_PRECISION',mpi_sum)

#ifdef LES
   if(ipimd>0) then
      call mpi_allreduce(frcx_copy,frcx_temp,3*ncopy,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
      frcx_copy = frcx_temp
   end if
#endif


# ifdef USE_MPI_IN_PLACE
   call mpi_allreduce(MPI_IN_PLACE,frcx,3,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
#else
   call mpi_allreduce(frcx,rl_temp,3,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
   call putm_back(frcx,rl_temp,3)
#endif

#else

   !  ---- Put force back in frc array from temp. force array ----

   nstart=1
   ntop=numatoms
#endif

#ifdef LES
   if(ipimd>0) frcx_copy=frcx_copy/numatoms
#endif

   frcx(1:3)=frcx(1:3)/dble(numatoms - numextra)
   if ( verbose > 0 )then
      if ( master )write(6,33)frcx(1),frcx(2),frcx(3)
      33 format(1x,'NET FORCE PER ATOM: ',3(1x,e12.4))
   end if
   !     ---do not remove net force if netfrc = 0; e.g. in minimization
   
   if ( netfrc > 0 )then
      do i = nstart,ntop

#ifdef LES
         if(ipimd>0) then
            if(cnum(i).eq.0) then
               do icopy=1,ncopy
                  frc(1:3,i) = frc(1:3,i)-frcx_copy(1:3,icopy)
               enddo
            else
               frc(1:3,i)=frc(1:3,i)-frcx_copy(1:3,cnum(i))
            endif
         else
            frc(1:3,i)=frc(1:3,i)-frcx(1:3)
         end if
#else
         frc(1,i) = frc(1,i) - frcx(1)
         frc(2,i) = frc(2,i) - frcx(2)
         frc(3,i) = frc(3,i) - frcx(3)
#endif
      end do
#ifndef LES
      ! Re-zero any extra points forces that got whacked above...
      if (numextra > 0) call zero_extra_pnts_vec(frc, ix)
#endif
   end if
   call timer_stop_start(TIME_EWFRCADJ,TIME_EWVIRIAL)   

   !--------------------------------------------------------
   !            NB  VIRIAL
   !--------------------------------------------------------

   !only get the virial if this is a periodic simulation.
   if (ntb>0) then
     call get_nonbond_virial(atvir,molvir,subvir, &
           rec_vir,dir_vir,adj_vir,rec_vird,self_vir, &
           les_vir, e14vir,framevir, &
           eelt,eedvir,virvsene,frc,crd,xr,numatoms)
  end if

#ifdef MPI
   
   !     Accumulate field; similar to accum force
   
! M-WJ
!   if ( induced == 1 ) call fsum(x(lfield),x(lfrctmp),r_stack(1),r_stack(1))
   if ( induced > 0 ) call fsum(x(lfield),x(lfrctmp),r_stack(1),r_stack(1))
!
#endif
   
   !     ---Now get dipole self energy if induced
   
! M-WJ
!   if ( induced == 1 )then
   if ( induced > 0 )then
!
      call get_dipinfo(numatoms,pol,x(lfield),x(linddip), &
            x(ldipvel),dipself,dipkine,diprms,dipndf,indmeth)
   end if
#ifdef MPI
   
   !     Accumulate the virials and energies
   
   call trace_mpi('mpi_allreduce', &
         BC_EW_COMM3,'MPI_DOUBLE_PRECISION',mpi_sum)
#  ifdef USE_MPI_IN_PLACE
   call mpi_allreduce(MPI_IN_PLACE,eer,BC_EW_COMM3, &
         MPI_DOUBLE_PRECISION, mpi_sum,commsander,ierr)
#  else
   call mpi_allreduce(eer,rl_temp,BC_EW_COMM3, &
         MPI_DOUBLE_PRECISION, mpi_sum,commsander,ierr)
   call putm_back(eer,rl_temp,BC_EW_COMM3)
#  endif

#  ifdef LES
#   ifdef USE_MPI_IN_PLACE
   call mpi_allreduce(MPI_IN_PLACE,eeles,10,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
#   else
   call mpi_allreduce(eeles,rl_temp,10,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
   call putm_back(eeles,rl_temp,10)
#   endif
#  endif
#endif /* MPI */
   
   !     finish up dipole stuff. Note that we just need to scale by
   !     degrees of freedom, which had to be summed over MPI procs.
   !     hence not done in get_dipinfo
   
! M-WJ
!   if ( induced == 1 )then
   if ( induced > 0 )then
!
      if ( indmeth == 3 )then
         boltz2 = 8.31441d-3 * 0.5d0
         boltz2 = boltz2/4.184d0
         dipkine = 0.5d0*dipmass*dipkine*INV_AMBER_ELECTROSTATIC2
         diptemp = dipkine/(boltz2*3.d0*dipndf)
      end if
! M-WJ
!     diprms = sqrt(diprms/dipndf)
      if (dipndf /= 0) diprms = sqrt(diprms/dipndf)
!
      
      !       ---change to Debye units
      
      diprms = diprms*4.8d0*INV_AMBER_ELECTROSTATIC
   end if

#ifdef LES
   eelt = ees + eer + eed + eea + eeles
#else
   eelt = ees + eer + eed + eea
#endif
   epol = epold + epola + epols + dipself + dipkine
   epol = epol + epol14
   
   evdw = evdw + evdwr
   
#ifdef LES

!  +---------------------------------------------------------------+
!  |  NOTE: evdwr (which is only calculated for the master) has    |
!  |  been mpi_allreduced above but this gets overcounted when     |
!  |  nrg_all is later mpi_allreduced before evb_ntrfc.  Only      |
!  |  accumulate evdwr on the master PE.                           | 
!  +---------------------------------------------------------------+

!   if(ipimd>0) nrg_all(:) = nrg_all(:) + evdwr/ncopy
   if( ( ipimd > 0 ) .and. master ) &
      nrg_all(:) = nrg_all(:) + evdwr/ncopy

#endif
   !  [need to do this here, after reduction of the virials, rather than
   !  in get_nonbond_vir().  (should figure out why....)]

   virvsene = rec_vir(1,1) + rec_vir(2,2) + rec_vir(3,3) + &
         eedvir + adj_vir(1,1)+adj_vir(2,2)+adj_vir(3,3) + &
         les_vir(1,1)+les_vir(2,2)+les_vir(3,3)
   
   !  Compare the virvsene number to the electrostatic energy:
   
   if (eelt /= 0.0d0) then
      virvsene = abs(virvsene+eelt)/abs(eelt)
   else
      virvsene = 0.0d0
   end if
   
   !     --- Dump info if verbose
   
   call force_info(ees,eer,eed,eea, &
         epol,epold,epola,epols,dipself,dipkine,induced, &
         eeles,les_vir, &
         evdw,ehb,eelt,virvsene,molvir,eedvir, &
         rec_vir,dir_vir,adj_vir,self_vir, &
         rec_vird,atvir,subvir,verbose)

   call timer_stop(TIME_EWVIRIAL)  
   call trace_exit( 'ewald_force' )

   return
end subroutine ewald_force 

!--------------------------------------------------------------------
subroutine assign_copy_charge(numatoms,charge,chg_single,id)
   implicit none
   integer numatoms,id
#  include "memory.h"
#  include "les.h"
   _REAL_  charge(*), chg_single(*)
   integer i

   do i=1,numatoms
      if( cnum(i).eq.0 ) then
         chg_single(i) = charge(i)
      else if( cnum(i).eq.id ) then
         chg_single(i) = charge(i)*ncopy
      else 
         chg_single(i) = 0.0
      end if
    end do
end subroutine

!--------------------------------------------------------------------
!         DO_PME_RECIP
!--------------------------------------------------------------------
subroutine do_pme_recip(mpoltype,numatoms,crd,charge,frc,dipole,   &
      field,prefac1,prefac2,prefac3,fftable,qm_pot_only)
   use ew_recip
   use ew_recip_spatial
   use ew_dipole_recip
   use nblist, only: recip, volume
   use fft,only:column_fft_flag
#ifdef LES
   use pimd_vars, only: ipimd, nrg_all
   use part_pimd_vars, only: ftmp, pimd_mmchg, frcx_copy
#endif
   implicit none
#  include "memory.h"

   integer mpoltype, numatoms,ierr
   _REAL_ crd(3,numatoms),charge(numatoms),frc(3,numatoms)
   _REAL_ dipole(3,numatoms),field(3,numatoms)
   integer nmine,nderiv,icopy
   logical,intent(in) :: qm_pot_only

#  include "ew_frc.h"
#  include "md.h"
#  include "ew_pme_recip.h"

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
#endif
   ! OUTPUT
   !       eer:  ewald reciprocal or k-space  energy
   !       frc forces incremented by k-space sum
   !       virial:  virial due to k-space sum (valid for atomic scaling;
   !                rigid molecule virial needs a correction term not
   !                computed here
   
   ! HEAP STORAGE:  These arrays need to be preserved throughout simulation
   _REAL_ prefac1(*),prefac2(*),prefac3(*),fftable(*)

#ifdef LES
   integer iatom, aid
   _REAL_ eer_sum
#endif

#ifdef LES
   if( ipimd > 0 ) then
      eer_sum = 0.0
      ftmp = 0.0
      pimd_mmchg = 0.0
   end if
#endif

   if ( mpoltype == 0 )then
#ifdef LES      
   if( ipimd > 0 ) then
      do icopy = 1, ncopy
         frcx = frcx_copy(:,icopy)
         call assign_copy_charge(numatoms, charge, pimd_mmchg, icopy)
         call do_pmesh_kspace(numatoms,crd,pimd_mmchg,ftmp,&
            prefac1,prefac2,prefac3,fftable,qm_pot_only)
         frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa
         eer_sum = eer_sum + eer
         nrg_all(icopy) = nrg_all(icopy) + eer/ncopy
         frcx_copy(:,icopy) = frcx
      end do
   else
      call do_pmesh_kspace(numatoms,crd,charge,frc, &
           prefac1,prefac2,prefac3,fftable,qm_pot_only)
           frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa
   endif

#else /* LES */

# ifdef MPI
      if(column_fft_flag)then
         call spatial_do_pmesh_kspace(numatoms,crd,charge,frc, &
              prefac1,prefac2,prefac3,fftable,qm_pot_only)
      else
# endif         
         call do_pmesh_kspace(numatoms,crd,charge,frc, &
             prefac1,prefac2,prefac3,fftable,qm_pot_only)

# ifdef MPI
      endif
# endif
      frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa

#endif /* LES */

! M-WJ
!   else if ( mpoltype == 1 )then
   else if ( mpoltype > 0 )then
!
      call do_pmesh_dipole_kspace( &
            numatoms,crd,charge,recip,volume,ew_coeff, &
            eer,frc,rec_vir,dipole,field,prefac1,prefac2,prefac3, &
            fftable,frcx(1),frcx(2),frcx(3))
      frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa
   end if /* mpoltype */

#ifdef LES
   if( ipimd > 0 ) then
      eer = eer_sum/ncopy
      frc = frc + ftmp/ncopy
      frcx_copy = frcx_copy/ncopy
   end if
#endif

end subroutine do_pme_recip

#ifdef MPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine putm_back here]
subroutine putm_back(x,y,n)
   implicit none
   integer, intent(in) :: n
   _REAL_, intent(out) :: x(n)
   _REAL_, intent(in) :: y(n)

!   integer :: i
   x(1:n) = y(1:n)
!   do i=1,n
!      x(i)=y(i)
!   end do
   return
end subroutine putm_back 
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit verbose information about the PME calculation.
subroutine force_info(ees,eer,eed,eea, &
      epol,epold,epola,epols,dipself,dipkine,induced, &
      eeles,les_vir, &
      evdw,ehb,eelt,virvsene,molvir,eedvir, &
      rec_vir,dir_vir,adj_vir,self_vir, &
      rec_vird,atvir,subvir,verbose)
   implicit none
   _REAL_ ees,eer,eed,eea,evdw,ehb,eelt, &
         epol,epold,epola,epols,dipself,dipkine, &
         eedvir,virvsene,molvir(3,3), &
         rec_vir(3,3),dir_vir(3,3),adj_vir(3,3),self_vir(3,3), &
         rec_vird(3,3),atvir(3,3),subvir(3,3)
   integer verbose,induced
#  include "extra.h"
   _REAL_ eeles,les_vir(3,3)
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
   
   if (.not. master) return

   if ( verbose >= 1 ) then
      write(6, '(4(/,5x,a,f24.12))' ) &
            'Evdw                   = ', evdw, &
            'Ehbond                 = ', ehb, &
            'Ecoulomb               = ', eelt, &
            'E polarization         = ', epol
      write(6, '(2(/,5x,a,f24.12))') &
            'Iso virial             = ', &
            molvir(1,1)+molvir(2,2)+molvir(3,3), &
            'Eevir vs. Ecoulomb     = ', virvsene
   end if
   if ( verbose >= 2 )then
      write(6, '(4(/,5x,a,f24.12),/)') &
            'E electrostatic (self) = ', ees, &
            '                (rec)  = ', eer, &
            '                (dir)  = ', eed, &
            '                (adj)  = ', eea &
#ifdef LES
            ,'             (LESadj)  = ', eeles &
#endif
            ;

! M-WJ
!      if ( induced == 1 )then
      if ( induced > 0 )then
!
         write(6, '(5x,a,a)') &
               'E polarization  (rec)  : ', &
               'Included in Electrostatic (rec)'
         write(6, '(5(/,5x,a,f24.12),/)') &
               'E polarization  (self) = ', epols, &
               '                (dir)  = ', epold, &
               '                (adj)  = ', epola, &
               'E dip induction (pot)  = ', dipself, &
               '            (kinetic)  = ', dipkine
      end if
      write(6,30)molvir(1,1),molvir(1,2),molvir(1,3)
      write(6,30)molvir(2,1),molvir(2,2),molvir(2,3)
      write(6,30)molvir(3,1),molvir(3,2),molvir(3,3)
      30 format(5x,'MOLECULAR VIRIAL: ',3(1x,e14.8))
   end if
   if ( verbose == 3 )then
      write(6,*)'--------------------------------------------'
      write(6,31)rec_vir(1,1),rec_vir(1,2),rec_vir(1,3)
      write(6,31)rec_vir(2,1),rec_vir(2,2),rec_vir(2,3)
      write(6,31)rec_vir(3,1),rec_vir(3,2),rec_vir(3,3)
      write(6,*)'..................'
      31 format(5x,'Reciprocal VIRIAL: ',3(1x,e14.8))
      write(6,32)dir_vir(1,1),dir_vir(1,2),dir_vir(1,3)
      write(6,32)dir_vir(2,1),dir_vir(2,2),dir_vir(2,3)
      write(6,32)dir_vir(3,1),dir_vir(3,2),dir_vir(3,3)
      write(6,*)'..................'
      32 format(5x,'Direct VIRIAL: ',3(1x,e14.8))
      write(6,38)eedvir
      write(6,*)'..................'
      38 format(5x,'Dir Sum EE vir trace: ',e14.8)
      write(6,33)adj_vir(1,1),adj_vir(1,2),adj_vir(1,3)
      write(6,33)adj_vir(2,1),adj_vir(2,2),adj_vir(2,3)
      write(6,33)adj_vir(3,1),adj_vir(3,2),adj_vir(3,3)
      write(6,*)'..................'
      33 format(5x,'Adjust VIRIAL: ',3(1x,e14.8))
      write(6,34)rec_vird(1,1),rec_vird(1,2),rec_vird(1,3)
      write(6,34)rec_vird(2,1),rec_vird(2,2),rec_vird(2,3)
      write(6,34)rec_vird(3,1),rec_vird(3,2),rec_vird(3,3)
      write(6,*)'..................'
      34 format(5x,'Recip Disp. VIRIAL: ',3(1x,e14.8))
      write(6,35)self_vir(1,1),self_vir(1,2),self_vir(1,3)
      write(6,35)self_vir(2,1),self_vir(2,2),self_vir(2,3)
      write(6,35)self_vir(3,1),self_vir(3,2),self_vir(3,3)
      write(6,*)'..................'
      35 format(5x,'Self VIRIAL: ',3(1x,e14.8))
#ifdef LES
      write(6,36)les_vir(1,1),les_vir(1,2),les_vir(1,3)
      write(6,36)les_vir(2,1),les_vir(2,2),les_vir(2,3)
      write(6,36)les_vir(3,1),les_vir(3,2),les_vir(3,3)
      write(6,*)'..................'
      36 format(5x,'LES VIRIAL: ',3(1x,e14.8))
#endif
      write(6,37)atvir(1,1),atvir(1,2),atvir(1,3)
      write(6,37)atvir(2,1),atvir(2,2),atvir(2,3)
      write(6,37)atvir(3,1),atvir(3,2),atvir(3,3)
      37 format(5x,'Atomic VIRIAL: ',3(1x,e14.8))
      write(6,*)'--------------------------------------------'
      write(6,39)subvir(1,1),subvir(1,2),subvir(1,3)
      write(6,39)subvir(2,1),subvir(2,2),subvir(2,3)
      write(6,39)subvir(3,1),subvir(3,2),subvir(3,3)
      39 format(5x,'Sub    VIRIAL: ',3(1x,e14.8))
      write(6,*)'--------------------------------------------'
   end if  ! ( verbose == 3 )
   return
end subroutine force_info 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_nonbond_virial here]
subroutine get_nonbond_virial(atvir,molvir,subvir, &
      rec_vir,dir_vir,adj_vir,rec_vird,self_vir, &
      les_vir, e14vir,framevir, &
      eelt,eedvir,virvsene,frc,crd,xr,numatoms)
   implicit none
   _REAL_ les_vir(3,3)
   _REAL_ e14vir(3,3),framevir(3,3)
   _REAL_ atvir(3,3),molvir(3,3),subvir(3,3)
   _REAL_ rec_vir(3,3),dir_vir(3,3), &
         adj_vir(3,3),rec_vird(3,3)
   _REAL_ self_vir(3,3)
   _REAL_ eedvir
   _REAL_ eelt,virvsene
   _REAL_ frc(3,*),crd(3,*),xr(3,*)
   integer numatoms

   integer i,j,n

#ifndef noVIRIAL
   !     The atomic virial needs to be corrected by terms from
   !     bonds,angles etc.  Rather than getting the force due to
   !     those terms alone without the nonbonds and computing the
   !     force*crd tensor we will compute the full force*crd
   !     tensor and then correct it with the atvir tensor which is the
   !     force*crd tensor computed  with force given by the nonbond
   !     terms.  Always compute atomic virial. Only compute molecular
   !     virial if xr is valid, allowing a correction term to be
   !     computed...

   !     rec_vir()  == reciprocal virial, computed in scalar_sumrc()
   !     rec_vird() == long-range vdw (dispersion) correction, computed in
   !                   vdw_correction()
   !     dir_vir()  == direct-space nonbond virial, computed in short_ene()
   !     adj_vir()  == corrections for 1-2, 1-3, etc, computed in nb_adjust()
   !     self_vir() == correction computed in self()
   !     e14vir()   == 1-4 term (correction?), computed in get_14_cg()
   !     framevir() == extra-point transfers, computed in orient_force()
   
   do j = 1,3
      do i = 1,3
         atvir(i,j) = rec_vir(i,j) + dir_vir(i,j) + adj_vir(i,j) + &
               rec_vird(i,j) + self_vir(i,j) + e14vir(i,j) + framevir(i,j)
#ifdef LES
         atvir(i,j) = atvir(i,j) + les_vir(i,j)
#endif
         ! subvir(i,j) = 0.d0
         molvir(i,j) = atvir(i,j)
      end do
   end do
   do n = 1,numatoms
      do j = 1,3
         do i = 1,3
            ! subvir(i,j) = subvir(i,j) - frc(i,n)*crd(j,n)
            molvir(i,j) = molvir(i,j) + frc(i,n)*xr(j,n)
         end do
      end do
   end do

#endif
   return
end subroutine get_nonbond_virial 


!     --- NB_ADJUST ---

!     ...the part of ewald due to gaussian counterion about an atom
!     you are bonded to or otherwise for which you do not compute the
!     nonbond pair force. NECESSARY since you ARE computing this pair
!     in the reciprocal sum...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust here]
subroutine nb_adjust(charge,eea,crd, &
      mask1,mask2, &
      numadjst,ewaldcof,eedtbdns, &
      eed_cub,eed_lin,frc,numatoms, &
      adj_vir,cutoffnb,ee_type,eedmeth)


   use constants, only : INVSQRTPI, third, half
   use decomp, only: decpr, decpair
   use crg_reloc, only: ifcr, cr_add_dcdr_factor
   use file_io_dat
#ifdef LES
   use pimd_vars, only: ipimd, nrg_all, nbead, nbead_inv
#endif
   implicit none

   _REAL_ cutoffnb
   _REAL_ charge(*),eea,crd(3,*)
   integer numadjst,numatoms
   integer mask1(*),mask2(*)
   integer ee_type,eedmeth
   _REAL_ ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
   _REAL_ adj_vir(3,3),frc(3,numatoms)

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer numlo,numhi
#ifdef MPI
   integer numdel,numleft
#endif

   integer i,k,nx,ind
   _REAL_ r2,r,x,df,cgi,cgk,erfc,derfc
   _REAL_ dfx,dfy,dfz
   _REAL_ del,dx,xx,delr2inv
   _REAL_ delx,dely,delz
   _REAL_ vxx,vxy,vxz,vyy,vyz,vzz
   _REAL_ d0,d1
   _REAL_ ecur
#include "box.h"
#include "md.h"

#ifdef LES

   _REAL_ xbig
   _REAL_ xsmall,rsmall,esmall,ezero,efit, &
         fsmall,fzero,ffit,prefac
   
   ! for LES, we do NOT use the scale factor for nb_adjust.
   ! this is because we are correcting for the interactions calculated
   ! in the reciprocal sum that should not have been included.
   ! so, we want to subtract exactly what we calculated- the non-scaled value.
   ! we do not scale the direct sum part since we CANNOT scale
   ! the reciprocal component.
   ! for part PIMD, we do use the scale factor for nb_adjust.
   ! since we have calculate the reciprocal energy correctly
   ! by call ew_recip P times
#include "les.h"

#endif


#  include "flocntrl.h"
   ! FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return

   !
   !     initialize virial
   !
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   !
   !     return if no terms to adjust!
   !
   if ( numadjst .eq. 0) return

   
   !     Get numlo, numhi for MPI and 1processor version
   
#ifdef MPI
   
   !     numdel is the amount of numadjst that each PE gets
   
   numdel = numadjst / numtasks
   if ( numdel == 0 )numdel = 1
   
   !--MFC-- numleft is the amout of numadjst left after numadjst is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra numadjst.
   
   numleft = mod(numadjst,numtasks)
   
   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible bu number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else    /*   MPI   */
   numlo = 1
   numhi = numadjst
#endif
#ifdef LES
   
   !     Carlos/Tom D. modifications for LES
   !     this is to account for possible short and long distances
   !     between atoms that are excluded
   
   prefac = 2.d0 * INVSQRTPI
   xsmall = 0.01d0
   rsmall = xsmall / ewaldcof
   call erfcfun(xsmall,erfc)
   derfc = prefac*exp(-xsmall**2)
   esmall = (erfc - 1.d0)/rsmall
   ezero = -ewaldcof * prefac
   
   fsmall = (erfc - 1.d0)/rsmall + ewaldcof*derfc
   fsmall = fsmall/(rsmall*rsmall)

   fzero = 0.d0

#endif /* LES */
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   eea = 0.d0
   del = 1.d0 / eedtbdns
   do nx = numlo,numhi
      i = mask1(nx)
      k = mask2(nx)


      cgi = charge(i)
      cgk = charge(k)
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)
      r2 = delx*delx+dely*dely+delz*delz
      r=sqrt(r2)
#ifdef TEST_CLOSEATOMS
      if(r < .8) &
            write(6,'("<NB_adjust> Atoms close:",2i8,e10.3)')i,k,r
#endif
      x = ewaldcof*r
#ifdef LES
      xbig=cutoffnb*ewaldcof
      
      !       special cases below for LES
      !       xbig is for excluded atoms that are far apart in space
      !       xsmall is for copies that may be overlapping
      
      if ( x > xbig )then
         
         !         put 1/r2 here so it won't cause NaN in xsmall routine
         
         r2=1.d0/r2
         eea = eea - cgi*cgk/r
         if ( ifcr /= 0 ) then
            call cr_add_dcdr_factor( i, -cgk/r )
            call cr_add_dcdr_factor( k, -cgi/r )
         end if
         df = -cgi*cgk*r2/r
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
      else if (x <= xsmall)then
         efit = (x/xsmall)*esmall + ((xsmall - x)/xsmall)*ezero
         ffit = (x/xsmall)*fsmall + ((xsmall - x)/xsmall)*fzero
         eea = eea + cgi*cgk*efit
         if ( ifcr /= 0 ) then
            call cr_add_dcdr_factor( i, cgk*efit )
            call cr_add_dcdr_factor( k, cgi*efit )
         end if
         
         !         *r2 accounted for in fsmall calc above
         
         !         df = cgi*cgk*(ffit)*r2
         df = cgi*cgk*(ffit)
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
      else

#endif /* LES */
         
         ! similar code to that in short_ene; however the only valid option
         !         here is the erfc switch, and dxdr = ewaldcof.
         
         if ( eedmeth == 1 )then
            
            !         ---cubic spline on erfc derfc
            
            ind = eedtbdns*x + 1
            dx = x - (ind-1)*del
            erfc = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                  dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
            derfc = (eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                  dx*eed_cub(4,ind)*half))
         else if ( eedmeth == 2 )then
            
            !         ---linear lookup on erfc, deriv
            
            xx = eedtbdns*x + 1
            ind = xx
            dx = xx - ind
            erfc = (1.d0 - dx)*eed_lin(1,ind) + &
                  dx*eed_lin(1,ind+1)
            derfc = (1.d0 - dx)*eed_lin(2,ind) + &
                  dx*eed_lin(2,ind+1)
         else if ( eedmeth == 3 )then
            call get_ee_func(x,erfc,derfc,ee_type)
         end if
         delr2inv = 1.d0 / r2
         d0 = (erfc - 1.d0)/r
         d1 = (d0 - ewaldcof*derfc)*delr2inv
#ifdef LES
         if( ipimd>0 ) then
            lfac = lesfac((lestyp(i)-1)*nlesty+lestyp(k))
            ecur = lfac*cgi*cgk*d0
            eea = eea + ecur
            if ( ifcr /= 0 ) then
               call cr_add_dcdr_factor( i, lfac*cgk*d0 )
               call cr_add_dcdr_factor( k, lfac*cgi*d0 )
            end if
            df = lfac * cgi * cgk * d1
             
            if(cnum(i).eq.0.and.cnum(k).eq.0) then
               nrg_all(1:nbead)=nrg_all(1:nbead) + ecur*nbead_inv
            else 
               if(cnum(i).ne.0) then
                  nrg_all(cnum(i)) = nrg_all(cnum(i)) + ecur
               else
                  nrg_all(cnum(k)) = nrg_all(cnum(k)) + ecur
               end if
            end if

         else
            eea = eea + cgi*cgk*d0
            if ( ifcr /= 0 ) then
               call cr_add_dcdr_factor( i, cgk*d0 )
               call cr_add_dcdr_factor( k, cgi*d0 )
            end if
            df = cgi*cgk*d1
         end if       
#else
         eea = eea + cgi*cgk*d0
         if ( ifcr /= 0 ) then
            call cr_add_dcdr_factor( i, cgk * d0 )   
            call cr_add_dcdr_factor( k, cgi * d0 )   
         end if
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,k,cgi*cgk*d0/(nstlim/ntpr))
         df = cgi*cgk*d1
#endif
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
#ifdef LES
         
         !       endif for (x.gt.xbig) above
         
      end if  ! ( x > xbig )
#endif
      vxx = vxx - dfx*delx
      vxy = vxy - dfx*dely
      vxz = vxz - dfx*delz
      vyy = vyy - dfy*dely
      vyz = vyz - dfy*delz
      vzz = vzz - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
   end do  !  nx = numlo,numhi
   adj_vir(1,1) = vxx
   adj_vir(1,2) = vxy
   adj_vir(2,1) = vxy
   adj_vir(1,3) = vxz
   adj_vir(3,1) = vxz
   adj_vir(2,2) = vyy
   adj_vir(2,3) = vyz
   adj_vir(3,2) = vyz
   adj_vir(3,3) = vzz
   return
end subroutine nb_adjust 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust_dipole here]
subroutine nb_adjust_dipole(charge,eea,crd, &
      mask1,mask2,numadjst,ewaldcof,eedtbdns, &
      eed_cub,eed_lin,frc,numatoms, &
      adj_vir,cutoffnb,ee_type,eedmeth, &
      dipole,field,epola)
   use constants, only : third, half
   implicit none

   _REAL_ cutoffnb
   _REAL_ charge(*),eea,crd(3,*)
   integer mask1(*),mask2(*),numadjst,numatoms
   integer ee_type,eedmeth
   _REAL_ ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
   _REAL_ adj_vir(3,3),frc(3,numatoms), &
         dipole(3,*),field(3,*),epola

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer numlo,numhi
#ifdef MPI
   integer numdel,numleft
#endif

   integer i,j,k,nx,ind
   _REAL_ r2,r,x,cgi,cgk,erfc,derfc
   _REAL_ dfx,dfy,dfz
   _REAL_ del,dx,xx,delr2inv
   _REAL_ delx,dely,delz
   _REAL_ d0,d1,d2,d3,fac,dotir,dotkr,dotik,dotp
   _REAL_ term,term0,term1,termi,termk,cgp
   _REAL_ dphii_dx,dphii_dy,dphii_dz, &
         dphik_dx,dphik_dy,dphik_dz

#  include "flocntrl.h"
   !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return
   
   !     Get numlo, numhi for MPI and 1processor version
   
#ifdef MPI
   
   !     numdel is the amount of numadjst that each PE gets
   
   numdel = numadjst / numtasks
   if ( numdel == 0 )numdel = 1

   !--MFC-- numleft is the amout of numadjst left after numadjst is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra numadjst.
   
   numleft = mod(numadjst,numtasks)
   
   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible bu number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else
   numlo = 1
   numhi = numadjst
#endif
   do j = 1,3
      do i = 1,3
         adj_vir(i,j) = 0.d0
      end do
   end do
   eea = 0.d0
   epola = 0.d0
   del = 1.d0 / eedtbdns
   fac = 2.d0*ewaldcof*ewaldcof
   do nx = numlo,numhi
      i = mask1(nx)
      k = mask2(nx)
      cgi = charge(i)
      cgk = charge(k)
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)
      r2 = delx*delx+dely*dely+delz*delz
      r=sqrt(r2)
      x = ewaldcof*r
      
      !       similar code to that in short_ene; however the only valid option
      !       here is the erfc switch, and dxdr = ewaldcof.
      
      if ( eedmeth == 1 )then
         
         !         ---cubic spline on erfc derfc
         
         ind = eedtbdns*x + 1
         dx = x - (ind-1)*del
         erfc = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
               dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
         derfc = (eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
               dx*eed_cub(4,ind)*half))
      else if ( eedmeth == 2 )then
         
         !         ---linear lookup on erfc, deriv
         
         xx = eedtbdns*x + 1
         ind = xx
         dx = xx - ind
         erfc = (1.d0 - dx)*eed_lin(1,ind) + &
               dx*eed_lin(1,ind+1)
         derfc = (1.d0 - dx)*eed_lin(2,ind) + &
               dx*eed_lin(2,ind+1)
      else if ( eedmeth == 3 )then
         call get_ee_func(x,erfc,derfc,ee_type)
      end if
      delr2inv = 1.d0 / r2
      d0 = (erfc - 1.d0)/r
      d1 = (d0 - ewaldcof*derfc)*delr2inv
      d2 = (3.d0*d1 - fac*ewaldcof*derfc)*delr2inv
      d3 = (5.d0*d2 - fac*fac*ewaldcof*derfc)*delr2inv
      dotkr = dipole(1,k)*delx+dipole(2,k)*dely+dipole(3,k)*delz
      dotir = dipole(1,i)*delx+dipole(2,i)*dely+dipole(3,i)*delz
      dotp = dotkr*dotir
      cgp = cgi*cgk
      dotik = dipole(1,i)*dipole(1,k)+dipole(2,i)*dipole(2,k)+ &
            dipole(3,i)*dipole(3,k)
      term = cgk*dotir-cgi*dotkr+dotik
      term0 = cgp*d0 + term*d1 - dotp*d2
      term1 = cgp*d1 + term*d2 - dotp*d3
      termi = cgi*d1+dotir*d2
      termk = cgk*d1 - dotkr*d2
      eea = eea + cgp*d0
      epola = epola + term*d1 - dotp*d2
      dfx = term1*delx + termi*dipole(1,k) - termk*dipole(1,i)
      dfy = term1*dely + termi*dipole(2,k) - termk*dipole(2,i)
      dfz = term1*delz + termi*dipole(3,k) - termk*dipole(3,i)
      adj_vir(1,1) = adj_vir(1,1) - dfx*delx
      adj_vir(1,2) = adj_vir(1,2) - dfx*dely
      adj_vir(1,3) = adj_vir(1,3) - dfx*delz
      adj_vir(2,1) = adj_vir(2,1) - dfy*delx
      adj_vir(2,2) = adj_vir(2,2) - dfy*dely
      adj_vir(2,3) = adj_vir(2,3) - dfy*delz
      adj_vir(3,1) = adj_vir(3,1) - dfz*delx
      adj_vir(3,2) = adj_vir(3,2) - dfz*dely
      adj_vir(3,3) = adj_vir(3,3) - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
      dphii_dx = termk*delx + d1*dipole(1,k)
      dphii_dy = termk*dely + d1*dipole(2,k)
      dphii_dz = termk*delz + d1*dipole(3,k)
      dphik_dx = -termi*delx + d1*dipole(1,i)
      dphik_dy = -termi*dely + d1*dipole(2,i)
      dphik_dz = -termi*delz + d1*dipole(3,i)
      field(1,i) = field(1,i) - dphii_dx
      field(2,i) = field(2,i) - dphii_dy
      field(3,i) = field(3,i) - dphii_dz
      field(1,k) = field(1,k) - dphik_dx
      field(2,k) = field(2,k) - dphik_dy
      field(3,k) = field(3,k) - dphik_dz
   end do  !  nx = numlo,numhi
   return
end subroutine nb_adjust_dipole 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine self here]
subroutine self(cg,numatoms,ene,ewaldcof,volume,self_vir,frc)
   use constants, only : PI, INVSQRTPI
   use decomp, only : decpair
   use crg_reloc, only : ifcr, cr_add_dcdr_factor
   use linear_response, only: ilrt
   use file_io_dat
#ifdef LES
   use pimd_vars, only: ipimd,nrg_all,nbead,nbead_inv
#endif
   implicit none
#  include "memory.h"
   integer numatoms
   _REAL_ cg(*),ene,ewaldcof,ee_plasma,volume,self_vir(3,3)
   _REAL_ sumq,sumq2
   integer i,j
   _REAL_ factor, d0, d1
   integer issetup
   _REAL_ frc(3,numatoms)
#  include "flocntrl.h"
#  include "md.h"

#ifdef LES
#  include "les.h"
#endif

   save issetup
   save sumq,sumq2
   data sumq,sumq2/0.d0,0.d0/
   data issetup/0/
   
   if ( do_self == 0 .and. ilrt == 0)return
   !     ---only compute sumq and sumq2 at beginning. They don't change
   ! if LIE module is in use, this is called repeatedly with different charge arrays, so recompute

   d0 = -ewaldcof*INVSQRTPI
   
   if ( issetup == 0 )then
      sumq = 0.d0
      sumq2 = 0.d0
      factor = sqrt(pi/(ewaldcof*ewaldcof*volume))*sqrt(0.5d0)
      do i = 1,numatoms
         sumq = sumq + cg(i)
#ifdef LES
         if(ipimd>0) then
            lfac=lesfac((lestyp(i)-1)*nlesty+lestyp(i))
            sumq2 = sumq2 + lfac*cg(i)*cg(i)
            if ( ifcr /= 0 ) then
               call cr_add_dcdr_factor( i, 2.0*cg(i)*lfac*d0 )
            end if
         else
            sumq2 = sumq2 + cg(i)*cg(i)
            if ( ifcr /= 0 ) then
               call cr_add_dcdr_factor( i, 2.0*cg(i)*d0 )
            end if
         endif
#else
         sumq2 = sumq2 + cg(i)*cg(i)
         if ( ifcr /= 0 ) then
            call cr_add_dcdr_factor( i, 2.0*cg(i)*d0 )
         end if
         ! -- ti decomp
         if(idecomp == 1 .or. idecomp == 2) then
            factor = pi/(ewaldcof*ewaldcof*volume)
            call decpair(2,i,i,cg(i)*cg(i)*d0)
         endif
#endif
      end do
      ! -- ti decomp
      if((idecomp == 1 .or. idecomp == 2) .and. -0.5d0*factor*sumq*sumq /= 0.0d0) then
         do i = 1,numatoms
            call decpair(2,i,i,-0.5d0*factor*sumq*sumq/numatoms) 
         end do
      end if
      if ( ilrt == 0 .and. ifcr == 0 ) then
         issetup = 1
      else
         ! recompute sum every time because total charge can change
      end if
   end if
   ene = sumq2*d0
   factor = pi / (ewaldcof*ewaldcof*volume)
   ee_plasma = -0.5d0*factor*sumq*sumq
   ! force related to ee_plasma
   if ( ifcr /= 0 ) then
      d1 = -factor * sumq
      do i = 1,numatoms
         call cr_add_dcdr_factor( i, d1 )
      end do
   end if
   ene = ene + ee_plasma
   do j = 1,3
      do i = 1,3
         self_vir(i,j) = 0.d0
      end do
      self_vir(j,j) = -ee_plasma
   end do
#ifdef LES
   if(ipimd>0) then
      do i=1,numatoms
         lfac=lesfac((lestyp(i)-1)*nlesty+lestyp(i))
 
         if(cnum(i).eq.0) then
            nrg_all(1:nbead)=nrg_all(1:nbead) &
                 + lfac*cg(i)*cg(i)*d0*nbead_inv
         else 
            nrg_all(cnum(i)) = nrg_all(cnum(i)) &
                 + lfac*cg(i)*cg(i)*d0
         endif
      end do
   end if
#endif

   return
end subroutine self 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine self_dipole here]
subroutine self_dipole(numatoms,epol,ewaldcof,efield,dipole, &
      iproc,nproc)
   use constants, only : SQRTPI
   implicit none
   integer numatoms,iproc,nproc
   _REAL_ efield(3,*),dipole(3,*),epol,ewaldcof

   _REAL_ fac_dip
   integer i
#  include "flocntrl.h"

   if ( do_self == 0 )return
   epol = 0.d0
   fac_dip = 4.d0*ewaldcof**3/(3.d0*SQRTPI)
   do i = 1,numatoms
      if ( iproc == mod(i,nproc) )then
         epol = epol - fac_dip* &
               (dipole(1,i)**2+dipole(2,i)**2+dipole(3,i)**2)
         efield(1,i) = efield(1,i) + fac_dip*dipole(1,i)
         efield(2,i) = efield(2,i) + fac_dip*dipole(2,i)
         efield(3,i) = efield(3,i) + fac_dip*dipole(3,i)
      end if
   end do
   epol = 0.5d0*epol
   return
end subroutine self_dipole 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dip_field_corr here]
subroutine dip_field_corr(efield,dipole,rec_vir,numatoms, &
      iproc,nproc)
   implicit none
   _REAL_ efield(3,*),dipole(3,*),rec_vir(3,3)
   integer numatoms,iproc,nproc
   
   ! The virial needs to be adjusted when dipoles are present
   ! the structure factors contribute a term, which is equal to
   ! the tensor product of the dipole and the reciprocal sum electric field
   ! this routine computes this term
   
   integer i,j,n
   do n = 1,numatoms
      !      if ( iproc .eq. mod(n,nproc) )then
      do j = 1,3
         do i = 1,3
            rec_vir(i,j) = rec_vir(i,j) + efield(i,n)*dipole(j,n)
         end do
      end do
      !      endif
   end do

   return
end subroutine dip_field_corr 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine respa_scale here]
subroutine respa_scale(numatoms,frc,field,nrespa)
   implicit none
   integer numatoms,nrespa
   _REAL_ frc(3,*),field(3,*)
   integer n
   do n = 1, numatoms
      frc(1,n) = nrespa*frc(1,n)
      frc(2,n) = nrespa*frc(2,n)
      frc(3,n) = nrespa*frc(3,n)
      !       field(1,n) = nrespa*field(1,n)
      !       field(2,n) = nrespa*field(2,n)
      !       field(3,n) = nrespa*field(3,n)
   end do
   return
end subroutine respa_scale 
!-------------------------------------------------------------------

#ifdef LES

!-------------------------------------------------------------------
!     --- NB_ADJUST_LES ---

!     Author: Carlos Simmerling

!     The electrostatic energies have not been corrected for the
!     intra-copy scaling factors. Here we go through all pairs
!     that need correcting (calculated at program start) and
!     correct the energy. The calculated energy was cgi*cgk/r,
!     and the proper LES energy is cgi*cgj/r*lfac where lfac
!     corresponds to the number of copies of that pair interaction.
!     The PME electrostatic energy was therefore underestimated by
!     the standard code.

!     This routine is based on the nb_adjust routine.



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust_les here]
subroutine nb_adjust_les(charge,ene,crd, &
      frc,numatoms,adj_vir,use_pme)

   implicit none
   _REAL_ charge(*),ene,crd(3,*)
   integer numatoms,use_pme
   _REAL_ frc(3,numatoms)
   _REAL_ adj_vir(3,3)
   integer numlo,numhi,numdel


#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer numleft, izero
   !-------------------------------------------------------------------
   _REAL_ tmpreal

   integer i,k,nx,ind
   _REAL_ r2,r,x,df,cgi,cgk,erfc,derfc
   _REAL_ dfx,dfy,dfz
   _REAL_ del,dx
   _REAL_ delx,dely,delz
   _REAL_ vxx,vxy,vxz,vyy,vyz,vzz
   
#  include "les.h"
#  include "box.h"
#  include "flocntrl.h"

   !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return
   
   ! if we are not doing PME, then return because the only calculation
   ! was for the direct sum and therefore does not need correcting here
   ! since it was done at the time the direct sum was calculated
   
   if (use_pme == 0) return

#ifdef MPI
   
   ! numdel is the amount of nlesadj that each PE gets
   
   numdel = nlesadj / numtasks
   if ( numdel == 0 )numdel = 1
   
   !--MFC-- numleft is the amout of nlesadj left after nlesadj is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra nlesadj.
   
   numleft = mod(nlesadj,numtasks)

   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible bu number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else
   numlo = 1
   numhi = nlesadj
#endif
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   eeles = 0.d0
   do nx = numlo,numhi
      i = ileslst(nx)
      k = jleslst(nx)
      cgi = charge(i)
      cgk = charge(k)
      lfac=lesfac((lestyp(i)-1)*nlesty+lestyp(k))
      
      !           should check here to make sure lfac is not 1.0! if it
      !           _is_ 1.0, then there is no need for the correction, and
      !           the list is incorrect.
      
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)

      !       r is really 1/r
      
      r = 1.d0/sqrt(delx*delx+dely*dely+delz*delz)
      r2 = r*r
      
      !       (lfac-1) adds the correct term and subtracts the incorrect
      !       one already added
      
      !       routine made more efficient with a
      !       temp variable for cgi*cgk*(lfac - 1.d0)*r
      
      tmpreal=cgi*cgk*(lfac - 1.d0)*r
      !if ( ifcr /= 0 ) then 
      !   call cr_add_dcdr_factor( i, cgk*(lfac-1.0d0)*r )
      !   call cr_add_dcdr_factor( k, cgi*(lfac-1.0d0)*r )
      !end if
      
      eeles = eeles + tmpreal
      
      !       have to correct the force in the same way, adding
      !       correct force and subtracting incorrect one
      
      df = tmpreal*r2
      dfx = delx*df
      dfy = dely*df
      dfz = delz*df
      vxx = vxx - dfx*delx
      vxy = vxy - dfx*dely
      vxz = vxz - dfx*delz
      vyy = vyy - dfy*dely
      vyz = vyz - dfy*delz
      vzz = vzz - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
   end do  !  nx = numlo,numhi
   ene = eeles
   adj_vir(1,1) = vxx
   adj_vir(1,2) = vxy
   adj_vir(2,1) = vxy
   adj_vir(1,3) = vxz
   adj_vir(3,1) = vxz
   adj_vir(2,2) = vyy
   adj_vir(2,3) = vyz
   adj_vir(3,2) = vyz
   adj_vir(3,3) = vzz
   return
end subroutine nb_adjust_les 
!-------------------------------------------------------------------
#endif /*LES*/
