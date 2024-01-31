#include "dprec.fh"

!-------------------------------------------------------------------------------
module amoeba_interface
    implicit none
    private
    public AMOEBA_readparm, &
        AMOEBA_deallocate, AM_VAL_eval, AM_NonBond_eval
contains
!-------------------------------------------------------------------------------
    subroutine AM_NONBOND_readparm(nf, natom, mass)
        use amoeba_multipoles, only : AM_MPOLE_readparm
        use amoeba_adjust, only : AM_ADJUST_readparm
        use amoeba_vdw, only : AM_VDW_read_parm
        use amoeba_induced, only : AM_INDUCED_readparm
        use amoeba_recip, only : AM_RECIP_allocate
        use amoeba_runmd, only : AM_RUNMD_init
        integer, intent(in) :: nf, natom

        integer :: mpole_valid, adjust_valid, vdw_valid, polar_valid
        _REAL_ mass(*)
        mpole_valid = AM_MPOLE_readparm(nf, natom)
        adjust_valid = AM_ADJUST_readparm(nf)
        vdw_valid = AM_VDW_read_parm(nf)
        polar_valid = AM_INDUCED_readparm(nf, natom, mass)
        ! alocate to initialize recip
        call AM_RECIP_allocate(natom)
        call AM_RUNMD_init(natom)
    end subroutine AM_NONBOND_readparm
!-------------------------------------------------------------------------------
    subroutine AM_VAL_readparm(nf, ntf, ntc)
        use amoeba_bonds, only : AM_BONDS_readparm
        use amoeba_ureyb, only : AM_UREYB_readparm
        use amoeba_reg_angles, only : AM_REG_ANGLES_readparm
        use amoeba_trig_angles, only : AM_TRIG_ANGLES_readparm
        use amoeba_opbend_angles, only : AM_OPBEND_ANGLES_readparm
        use amoeba_torsions, only : AM_TORSIONS_readparm
        use amoeba_stretch_torsions, only : AM_STRETCH_TORSIONS_readparm
        use amoeba_pitorsions, only : AM_PITORSIONS_readparm
        use amoeba_stretch_bend, only : AM_STRETCH_BEND_readparm
        use amoeba_torsion_torsion, only : AM_TOR_TOR_readparm
        integer, intent(in) :: nf
        integer, intent(inout) :: ntf, ntc !set these to 8,1 if amoeba-friendly prmtop
        !so regular amber valence & shake not run

        integer :: bonds_valid, ureyb_valid, angles_valid, trig_angles_valid, &
            opbend_valid, torsions_valid, pitorsions_valid, strechbend_valid, &
            torsion_torsion_valid, stretch_torsion_valid
        bonds_valid = AM_BONDS_readparm(nf)
        ureyb_valid = AM_UREYB_readparm(nf)
        angles_valid = AM_REG_ANGLES_readparm(nf)
        trig_angles_valid = AM_TRIG_ANGLES_readparm(nf)
        opbend_valid = AM_OPBEND_ANGLES_readparm(nf)
        torsions_valid = AM_TORSIONS_readparm(nf)
        stretch_torsion_valid = AM_STRETCH_TORSIONS_readparm(nf)
        pitorsions_valid = AM_PITORSIONS_readparm(nf)
        strechbend_valid = AM_STRETCH_BEND_readparm(nf)
        torsion_torsion_valid = AM_TOR_TOR_readparm(nf)
        if ((bonds_valid == 1) .or. (ureyb_valid == 1) .or. &
            (angles_valid == 1) .or. (trig_angles_valid == 1) .or. &
            (opbend_valid == 1) .or. (torsions_valid == 1) .or. &
            (pitorsions_valid == 1) .or. (strechbend_valid == 1) .or. &
            (torsion_torsion_valid == 1)) then
            ntf = 8 !make sure no amber valence calculations performed
            ntc = 1 !make sure SHAKE not performed
        end if
    end subroutine AM_VAL_readparm
!-------------------------------------------------------------------------------
    subroutine AMOEBA_readparm(nf, ntf, ntc, numatoms, mass)
        use amoeba_mdin, only : iamoeba
        integer, intent(in) :: nf, numatoms
        integer, intent(inout) :: ntf, ntc !set these to 8,1 if amoeba-friendly prmtop
        !so regular amber valence & shake not run
        _REAL_ mass(*)
        call AMOEBA_check_parm_legal(nf)
        if (iamoeba /= 1) return
        call AM_VAL_readparm(nf, ntf, ntc)
        call AM_NONBOND_readparm(nf, numatoms, mass)
    end subroutine AMOEBA_readparm
!-------------------------------------------------------------------------------
    subroutine AM_VAL_deallocate()
        use amoeba_bonds, only : AM_BONDS_deallocate
        use amoeba_ureyb, only : AM_UREYB_deallocate
        use amoeba_reg_angles, only : AM_REG_ANGLES_deallocate
        use amoeba_trig_angles, only : AM_TRIG_ANGLES_deallocate
        use amoeba_opbend_angles, only : AM_OPBEND_ANGLES_deallocate
        use amoeba_torsions, only : AM_TORSIONS_deallocate
        use amoeba_pitorsions, only : AM_PITORSIONS_deallocate
        use amoeba_stretch_bend, only : AM_STRETCH_BEND_deallocate
        use amoeba_torsion_torsion, only : AM_TOR_TOR_deallocate
        implicit none
        call AM_BONDS_deallocate()
        call AM_UREYB_deallocate()
        call AM_REG_ANGLES_deallocate()
        call AM_TRIG_ANGLES_deallocate()
        call AM_OPBEND_ANGLES_deallocate()
        call AM_TORSIONS_deallocate()
        call AM_PITORSIONS_deallocate()
        call AM_STRETCH_BEND_deallocate()
        call AM_TOR_TOR_deallocate()
    end subroutine AM_VAL_deallocate
!-------------------------------------------------------------------------------
    subroutine AM_NONBOND_deallocate()
        use amoeba_multipoles, only : AM_MPOLE_deallocate
        call AM_MPOLE_deallocate()
    end subroutine AM_NONBOND_deallocate
!-------------------------------------------------------------------------------
    subroutine AMOEBA_deallocate()
        call AM_VAL_deallocate()
        call AM_NONBOND_deallocate()
    end subroutine AMOEBA_deallocate
!-------------------------------------------------------------------------------
    subroutine AM_VAL_eval(crd, frc, sander_vir, ebond, eangle, etors)
        use amoeba_bonds, only : AM_BONDS_eval
        use amoeba_ureyb, only : AM_UREYB_eval
        use amoeba_reg_angles, only : AM_REG_ANGLES_eval
        use amoeba_trig_angles, only : AM_TRIG_ANGLES_eval
        use amoeba_opbend_angles, only : AM_OPBEND_ANGLES_eval
        use amoeba_torsions, only : AM_TORSIONS_eval
        use amoeba_stretch_torsions, only : AM_STRETCH_TORSIONS_eval
        use amoeba_pitorsions, only : AM_PITORSIONS_eval
        use amoeba_stretch_bend, only : AM_STRETCH_BEND_eval
        use amoeba_torsion_torsion, only : AM_TOR_TOR_eval
        use amoeba_mdin, only : iamoeba, do_amoeba_valence, amoeba_verbose
        implicit none
# include "extra.h"
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), sander_vir(4)
        _REAL_, intent(inout) :: ebond, eangle, etors

        _REAL_ :: ene(10), vir_tensor(3, 3)

        ene = 0.d0
        vir_tensor(:, :) = 0.d0
        if (iamoeba == 1 .and. do_amoeba_valence == 1) then
            call AM_BONDS_eval(crd, frc, ene(1), vir_tensor)
            call AM_UREYB_eval(crd, frc, ene(2), vir_tensor)
            call AM_REG_ANGLES_eval(crd, frc, ene(3), vir_tensor)
            call AM_TRIG_ANGLES_eval(crd, frc, ene(4), vir_tensor)
            call AM_OPBEND_ANGLES_eval(crd, frc, ene(5), vir_tensor)
            call AM_TORSIONS_eval(crd, frc, ene(6), vir_tensor)
            call AM_PITORSIONS_eval(crd, frc, ene(7), vir_tensor)
            call AM_STRETCH_BEND_eval(crd, frc, ene(8), vir_tensor)
            call AM_TOR_TOR_eval(crd, frc, ene(9), vir_tensor)
            call AM_STRETCH_TORSIONS_eval(crd, frc, ene(10), vir_tensor)
        end if
        if (master .and. amoeba_verbose > 0) then
            write (6, '(a,3(1x,e16.8))') &
                'valence energies: bond,ureyb,angle ', ene(1), ene(2), ene(3)
            write (6, '(a,3(1x,f14.4))') &
                'valence energies: trangle,opbend,tor ', ene(4), ene(5), ene(6)
            write (6, '(a,3(1x,f14.4))') &
                'valence energies: pitors,strbend,tortor ', ene(7), ene(8), ene(9)
            write (6, '(a,1x,f14.4)') &
                'valence energies: strtor ', ene(10)
            write (6, '(a,3(1x,f14.4))') &
                'valence vir = ', vir_tensor(1, 1), vir_tensor(1, 2), vir_tensor(1, 3)
            write (6, '(a,3(1x,f14.4))') &
                'valence vir = ', vir_tensor(2, 1), vir_tensor(2, 2), vir_tensor(2, 3)
            write (6, '(a,3(1x,f14.4))') &
                'valence vir = ', vir_tensor(3, 1), vir_tensor(3, 2), vir_tensor(3, 3)
        end if
        ebond = ene(1)!bonds energy
        eangle = ene(2) + ene(3) + ene(4) + ene(5) + ene(8)!angle energy
        etors = ene(6) + ene(7) + ene(9) + ene(10)!torsion energy
        sander_vir(1) = sander_vir(1) + vir_tensor(1, 1)
        sander_vir(2) = sander_vir(2) + vir_tensor(2, 2)
        sander_vir(3) = sander_vir(3) + vir_tensor(3, 3)
        sander_vir(4) = sander_vir(4) + &
            vir_tensor(1, 1) + vir_tensor(2, 2) + vir_tensor(3, 3)
    end subroutine AM_VAL_eval
!-------------------------------------------------------------------------------
    subroutine AM_NonBond_eval(numatoms, crd, frc, sander_vir, x, ipairs, &
        evdw, eelt, epolar, evdw_14, eelt_14, &
        diprms, dipiter)
        use amoeba_mdin, only : iamoeba, do_amoeba_nonbond, amoeba_verbose
        use amoeba_multipoles, only : AM_MPOLE_local_to_global, torque_field, &
            global_multipole, AM_MPOLE_torque_to_force
        use amoeba_induced, only : AM_INDUCED_eval
        use nblist, only : recip, adjust_imagcrds, map_coords
        integer, intent(in) :: numatoms
        _REAL_, intent(in) ::  crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), sander_vir(4)
        _REAL_, intent(in) :: x(*)
        integer, intent(in) :: ipairs(*)
        _REAL_, intent(out) :: evdw, eelt, epolar, evdw_14, eelt_14, dipiter, diprms

        integer :: j, k
        _REAL_ :: vir_tensor(3, 3)
        evdw = 0.d0
        eelt = 0.d0
        epolar = 0.d0
        evdw_14 = 0.d0
        eelt_14 = 0.d0
        do k = 1, 3
            do j = 1, 3
                vir_tensor(j, k) = 0.d0
            end do
        end do
        if (iamoeba == 1 .and. do_amoeba_nonbond == 1) then

            call zero_array(torque_field, 10*numatoms)
            call zero_array(global_multipole, 10*numatoms)

            ! update the imaged crds
            call map_coords(crd, numatoms, recip)
            call adjust_imagcrds(crd, numatoms)

            call AM_MPOLE_local_to_global(crd)

            call AM_INDUCED_eval(numatoms, crd, x, ipairs, diprms, dipiter)

            call AM_NonBond_ene_frc(numatoms, crd, x, ipairs, &
                eelt, epolar, evdw, evdw_14, frc, vir_tensor)

            !add the torque contributions
            call AM_MPOLE_torque_to_force(numatoms, crd, frc, vir_tensor)

            if (amoeba_verbose > 0) then
                write (6, '(a,3(1x,g16.8))') &
                    'nonbond vir = ', vir_tensor(1, 1), vir_tensor(1, 2), vir_tensor(1, 3)
                write (6, '(a,3(1x,g16.8))') &
                    'nonbond vir = ', vir_tensor(2, 1), vir_tensor(2, 2), vir_tensor(2, 3)
                write (6, '(a,3(1x,g16.8))') &
                    'nonbond vir = ', vir_tensor(3, 1), vir_tensor(3, 2), vir_tensor(3, 3)
            end if

        end if
        sander_vir(1) = sander_vir(1) + vir_tensor(1, 1)
        sander_vir(2) = sander_vir(2) + vir_tensor(2, 2)
        sander_vir(3) = sander_vir(3) + vir_tensor(3, 3)
        sander_vir(4) = sander_vir(4) + &
            vir_tensor(1, 1) + vir_tensor(2, 2) + vir_tensor(3, 3)

    end subroutine AM_NonBond_eval
!-------------------------------------------------------------------------------
end module amoeba_interface
!-------------------------------------------------------------------------------
subroutine AMOEBA_check_parm_legal(nf)
    use amoeba_mdin, only : iamoeba
    integer, intent(in) :: nf

    integer :: iok, ionerr, indicator
    character(len=80) :: fmt, fmtin, dtype

    fmtin = '(I5)'
    dtype = 'AMOEBA_FORCEFIELD'
    ionerr = 1 ! not fatal if missing
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    if (iok == 0) then !this data type found in prmtop
        read (nf, fmt) indicator
        if (iamoeba /= 1) then
            write (6, *) 'AMOEBA style prmtop, but amoeba NOT set to one'
            call mexit(6, 1)
        end if
        return ! amoeba prmtop, amoeba set to 1
    else
        if (iamoeba == 1) then
            write (6, *) 'NOT an AMOEBA style prmtop, but amoeba IS set to one'
            call mexit(6, 1)
        end if
        return ! NON amoeba prmtop, amoeba NOT set to 1
    end if
end subroutine AMOEBA_check_parm_legal
!-------------------------------------------------------------------------------
subroutine AMOEBA_check_newstyle_inpcrd(inpcrd, newstyle)
    character(len=*), intent(in) :: inpcrd
    logical, intent(out) :: newstyle

    integer :: nf, iok, ionerr
    character(len=80) :: fmt, fmtin, dtype
    nf = 30
    call amopen(nf, inpcrd, 'O', 'F', 'R')
    call nxtsec_crd_reset()
    fmtin = '(a)'
    dtype = 'TITLE'
    ionerr = 1 ! not fatal if missing
    call nxtsec_crd(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    if (iok == 0) then !this data type found in prmtop
        newstyle = .true.
    else
        newstyle = .false.
    end if
    call nxtsec_crd_reset()
    close (nf)
end subroutine AMOEBA_check_newstyle_inpcrd
!-------------------------------------------------------------------------------
subroutine AM_NonBond_perm_fields(numatoms, is_polarizable, crd, x, ipairs, &
    dip_field_d, dip_field_p, cart_dipole_field)
    use amoeba_recip, only : AM_RECIP_perm_field
    use amoeba_direct, only : AM_DIRECT_permfield
    use amoeba_adjust, only : AM_ADJUST_permfield
    use amoeba_self, only : AM_SELF_permfield
    integer, intent(in) :: numatoms
    logical, intent(in) :: is_polarizable(*)
    _REAL_, intent(in) ::  crd(3, *)
    _REAL_, intent(in) :: x(*)
    integer, intent(in) :: ipairs(*)
    _REAL_, intent(out) :: dip_field_d(3, *), dip_field_p(3, *), cart_dipole_field(3, *)
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
    logical master
#endif

    call zero_array(dip_field_d, 3*numatoms)
    call zero_array(dip_field_p, 3*numatoms)
    call zero_array(cart_dipole_field, 3*numatoms)

#ifdef MPI
    call mpi_comm_size(recip_comm, numtasks, ierr)
    call mpi_comm_rank(recip_comm, mytaskid, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_RECIP_perm_field(numatoms, crd, cart_dipole_field, x)

#ifdef MPI
    call mpi_comm_rank(commsander, mytaskid, ierr)
    call mpi_comm_size(commsander, numtasks, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_DIRECT_permfield(ipairs, x, cart_dipole_field)
    call AM_ADJUST_permfield(crd, x, dip_field_d, dip_field_p)
    call AM_SELF_permfield(numatoms, dip_field_d, dip_field_p)
    call AM_INDUCED_add_cart_to_dip(numatoms, is_polarizable, &
        cart_dipole_field, &
        dip_field_d, dip_field_p)
end subroutine AM_NonBond_perm_fields
!-------------------------------------------------------------------------------
subroutine AM_NonBond_dip_dip_fields(numatoms, x, &
    ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
    use amoeba_recip, only : AM_RECIP_dipole_field
    use amoeba_direct, only : AM_DIRECT_dip_dip_field
    use amoeba_adjust, only : AM_ADJUST_dip_dip_fields
    use amoeba_self, only : AM_SELF_dipole_field
    integer, intent(in) :: numatoms
    _REAL_, intent(in) :: x(*)
    _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
    _REAL_, intent(out) :: dip_field_d(3, *), dip_field_p(3, *)
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
    logical master
#endif

    call zero_array(dip_field_d, 3*numatoms)
    call zero_array(dip_field_p, 3*numatoms)

#ifdef MPI
    call mpi_comm_size(recip_comm, numtasks, ierr)
    call mpi_comm_rank(recip_comm, mytaskid, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_RECIP_dipole_field(numatoms, x, &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)

#ifdef MPI
    call mpi_comm_rank(commsander, mytaskid, ierr)
    call mpi_comm_size(commsander, numtasks, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_DIRECT_dip_dip_field( &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
    call AM_ADJUST_dip_dip_fields(ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
    call AM_SELF_dipole_field(numatoms, ind_dip_d, ind_dip_p, dip_field_d, &
        dip_field_p)
end subroutine AM_NonBond_dip_dip_fields
!-------------------------------------------------------------------------------
subroutine AM_NonBond_ene_frc(numatoms, crd, x, ipairs, &
    ene_perm, ene_ind, ene_vdw, &
    ene_vdw_14, frc, vir_tensor)
    use amoeba_recip, only : AM_RECIP_ene_frc
    use amoeba_direct, only : AM_DIRECT_ene_frc
    use amoeba_adjust, only : AM_ADJUST_ene_frc
    use amoeba_self, only : AM_SELF_ene_torque
    use amoeba_vdw, only : AM_VDW_longrange_ene
    use amoeba_induced, only : ind_dip_d, ind_dip_p
    use amoeba_mdin, only : amoeba_verbose
    integer, intent(in) :: numatoms
    _REAL_, intent(in) :: crd(3, *), x(*)
    integer, intent(in) :: ipairs(*)
    _REAL_, intent(out) :: ene_perm, ene_ind, ene_vdw, ene_vdw_14
    _REAL_, intent(inout) :: frc(3, *), vir_tensor(3, 3)

#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
    logical master
#endif

    _REAL_ :: e_rec_perm, e_rec_ind, e_dir_perm, e_dir_ind, e_adj_perm, &
        e_adj_ind, e_self_perm, e_self_ind, e_dir_vdw, e_adj_vdw, e_rec_vdw

#ifdef MPI
    call mpi_comm_size(recip_comm, numtasks, ierr)
    call mpi_comm_rank(recip_comm, mytaskid, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_RECIP_ene_frc(numatoms, crd, x, ind_dip_d, ind_dip_p, &
        e_rec_perm, e_rec_ind, frc, vir_tensor)

#ifdef MPI
    call mpi_comm_rank(commsander, mytaskid, ierr)
    call mpi_comm_size(commsander, numtasks, ierr)
    master = mytaskid .eq. 0
#endif

    call AM_NonBond_remove_net_force(numatoms, frc)
    call AM_DIRECT_ene_frc(ipairs, crd, x, ind_dip_d, ind_dip_p, &
        e_dir_perm, e_dir_ind, e_dir_vdw, frc, vir_tensor)
    call AM_ADJUST_ene_frc(crd, x, ind_dip_d, ind_dip_p, &
        e_adj_perm, e_adj_ind, e_adj_vdw, frc, vir_tensor)
    call AM_SELF_ene_torque(numatoms, ind_dip_d, ind_dip_p, &
        e_self_perm, e_self_ind)
    call AM_VDW_longrange_ene(e_rec_vdw, vir_tensor)

    if (amoeba_verbose > 0) then
        write (6, '(a,/,4(1x,f14.4))') &
            'e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm = ', &
            e_rec_perm, e_dir_perm, e_adj_perm, e_self_perm
        write (6, '(a,/,4(1x,f14.4))') &
            'e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind = ', &
            e_rec_ind, e_dir_ind, e_adj_ind, e_self_ind
        write (6, '(a,/,3(1x,f14.4))') &
            'e_dir_vdw,e_adj_vdw,e_rec_vdw = ', &
            e_dir_vdw, e_adj_vdw, e_rec_vdw
    end if
    ene_perm = e_rec_perm + e_dir_perm + e_adj_perm + e_self_perm
    ene_ind = e_rec_ind + e_dir_ind + e_adj_ind + e_self_ind
    ene_vdw = e_dir_vdw + e_rec_vdw
    ene_vdw_14 = e_adj_vdw
end subroutine AM_NonBond_ene_frc
!-------------------------------------------------------------------------------
subroutine AM_NonBond_remove_net_force(numatoms, frc)
    integer, intent(in) :: numatoms
    _REAL_, intent(inout) :: frc(3, *)

    _REAL_ :: frcx, frcy, frcz
    integer n
    frcx = 0.d0
    frcy = 0.d0
    frcz = 0.d0
    do n = 1, numatoms
        frcx = frcx + frc(1, n)
        frcy = frcy + frc(2, n)
        frcz = frcz + frc(3, n)
    end do
    frcx = frcx/numatoms
    frcy = frcy/numatoms
    frcz = frcz/numatoms
    do n = 1, numatoms
        frc(1, n) = frc(1, n) - frcx
        frc(2, n) = frc(2, n) - frcy
        frc(3, n) = frc(3, n) - frcz
    end do

end subroutine AM_NonBond_remove_net_force
!-------------------------------------------------------------
subroutine AM_NONBOND_set_user_bit(do_recip, do_adjust, do_direct, do_self, &
    do_vdw, do_induce)
    use amoeba_recip, only : AM_RECIP_set_user_bit
    use amoeba_adjust, only : AM_ADJUST_set_user_bit
    use amoeba_direct, only : AM_DIRECT_set_user_bit
    use amoeba_self, only : AM_SELF_set_user_bit
    use amoeba_vdw, only : AM_VDW_set_user_bit
    use amoeba_induced, only : AM_INDUCED_set_user_bit
    implicit none
    integer, intent(in) :: do_recip, do_adjust, do_direct, do_self, &
        do_vdw, do_induce
    call AM_RECIP_set_user_bit(do_recip)
    call AM_ADJUST_set_user_bit(do_adjust)
    call AM_DIRECT_set_user_bit(do_direct)
    call AM_SELF_set_user_bit(do_self)
    call AM_VDW_set_user_bit(do_vdw)
    call AM_INDUCED_set_user_bit(do_induce)
end subroutine AM_NONBOND_set_user_bit
!-------------------------------------------------------------------------------
subroutine AM_VAL_set_user_bit(do_bond, do_ureyb, do_reg_angle, do_trig_angle, &
    do_opbend_angle, do_torsions, do_str_torsions, do_pitorsions, &
    do_stretch_bend, do_torsion_torsion)
    use amoeba_bonds, only : AM_BONDS_set_user_bit
    use amoeba_ureyb, only : AM_UREYB_set_user_bit
    use amoeba_reg_angles, only : AM_REG_ANGLES_set_user_bit
    use amoeba_trig_angles, only : AM_TRIG_ANGLES_set_user_bit
    use amoeba_opbend_angles, only : AM_OPBEND_ANGLES_set_user_bit
    use amoeba_torsions, only : AM_TORSIONS_set_user_bit
    use amoeba_stretch_torsions, only : AM_STRETCH_TORSIONS_suser_bit
    use amoeba_pitorsions, only : AM_PITORSIONS_set_user_bit
    use amoeba_stretch_bend, only : AM_STRETCH_BEND_set_user_bit
    use amoeba_torsion_torsion, only : AM_TOR_TOR_set_user_bit
    implicit none
    integer, intent(in) :: do_bond, do_ureyb, do_reg_angle, do_trig_angle, &
        do_opbend_angle, do_torsions, do_pitorsions, &
        do_str_torsions, do_stretch_bend, do_torsion_torsion
    call AM_BONDS_set_user_bit(do_bond)
    call AM_UREYB_set_user_bit(do_ureyb)
    call AM_REG_ANGLES_set_user_bit(do_reg_angle)
    call AM_TRIG_ANGLES_set_user_bit(do_trig_angle)
    call AM_OPBEND_ANGLES_set_user_bit(do_opbend_angle)
    call AM_TORSIONS_set_user_bit(do_torsions)
    call AM_STRETCH_TORSIONS_suser_bit(do_str_torsions)
    call AM_PITORSIONS_set_user_bit(do_pitorsions)
    call AM_STRETCH_BEND_set_user_bit(do_stretch_bend)
    call AM_TOR_TOR_set_user_bit(do_torsion_torsion)
end subroutine AM_VAL_set_user_bit
!-------------------------------------------------------------------------------
subroutine AMOEBA_get_numlist(header, nf, num_list)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(out) :: num_list

    integer :: iok, ionerr
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    fmtin = '(10I8)'
    dtype = header//'NUM_LIST'
    ionerr = 1 ! not fatal if missing
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    if (iok == 0) then !this data type found in prmtop
        read (nf, fmt) num_list
    else !either old style prmtop or data not found
        num_list = 0 ! upon return this will invalidate valence_term
    end if
end subroutine AMOEBA_get_numlist
!------------------------------------------------------------------------
subroutine AMOEBA_read_list_data(header, nf, dim1, num_list, list)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf, dim1, num_list
    integer, intent(out) :: list(dim1, num_list)

    integer :: iok, ionerr, j, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(10I8)'
    dtype = header//'LIST'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) ((list(j, k), j=1, dim1), k=1, num_list)
end subroutine AMOEBA_read_list_data
!----------------------------------------------------------
subroutine AMOEBA_read_real_list_data(header, nf, dim1, num_list, list)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf, dim1, num_list
    _REAL_, intent(out) :: list(dim1, num_list)

    integer :: iok, ionerr, j, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'LIST'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) ((list(j, k), j=1, dim1), k=1, num_list)
end subroutine AMOEBA_read_real_list_data
!--------------------------------------------------------------------
subroutine AMOEBA_read_real_scalar(flag, nf, scalar_value)
    implicit none
    character(len=*), intent(in) :: flag
    integer, intent(in) :: nf
    _REAL_, intent(out) :: scalar_value

    integer :: iok, ionerr
    character(len=80) :: fmt
    character(len=80) :: fmtin

    ionerr = 0 !fatal if missing
    fmtin = '(E16.8)'
    call nxtsec(nf, 6, ionerr, fmtin, flag, fmt, iok)
    read (nf, fmt) scalar_value
end subroutine AMOEBA_read_real_scalar
!--------------------------------------------------------------------
subroutine AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
    implicit none
    integer, intent(in) :: num_list
    integer, intent(out) :: startlist, endlist, siztask

    integer piece
    integer numtasks, mytaskid

    numtasks = 1
    mytaskid = 0
    if (numtasks > 1) then
        piece = num_list/numtasks
        startlist = mytaskid*piece + 1
        endlist = mytaskid*piece + piece
        if (mytaskid == (numtasks - 1)) endlist = num_list
        siztask = endlist - startlist + 1
    else
        startlist = 1
        endlist = num_list
        siztask = endlist - startlist + 1
    end if
end subroutine AMOEBA_get_startlist_endlist
!-------------------------------------------------------------
subroutine array_add(a, b, num)
    _REAL_, intent(inout) :: a(*)
    _REAL_, intent(in) :: b(*)
    integer, intent(in) :: num

    integer :: n
    do n = 1, num
        a(n) = a(n) + b(n)
    end do
end subroutine array_add
!-------------------------------------------------------------
subroutine dump_dipoles(dip, nsites, nf)
    implicit none
    _REAL_ dip(3, *)
    integer nsites, nf
    integer j, n

    do n = 1, nsites
        write (nf, '(3f20.12)') (dip(j, n), j=1, 3)
    end do
    call amflsh(nf)
end subroutine dump_dipoles
