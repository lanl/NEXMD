#include "copyright.i"

!*******************************************************************************
!
! Module:  prmtop_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module prmtop_dat_mod

  use gbl_datatypes_mod
  use file_io_dat_mod

  implicit none

! Global parameter definitions:

! The following parameter MAY need to be increased for a very large system:

  integer, parameter    :: max_dihed_dups = 10000

! Global data.  This stuff should be broadcast:

  ! Starting with nttyp, the integer values are derived rather than read.
  ! nttyp is the number of 6-12 vdw parameters.

  integer, parameter    :: prmtop_int_cnt = 26

  integer               natom, ntypes, nbonh, ntheth, nphih, next, nres, &
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &
                        ifbox, ifcap, nspm, numextra, ncopy, nttyp, &
                        bonda_idx, anglea_idx, diheda_idx, &
                        gbl_bond_allocsize, gbl_angle_allocsize, &
                        gbl_dihed_allocsize

  common / prmtop_int / natom, ntypes, nbonh, ntheth, nphih, next, nres, &
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &
                        ifbox, ifcap, nspm, numextra, ncopy, nttyp, &
                        bonda_idx, anglea_idx, diheda_idx, &
                        gbl_bond_allocsize, gbl_angle_allocsize, &
                        gbl_dihed_allocsize

  save  :: / prmtop_int /

! Cap data not currently supported...

! integer, parameter    :: prmtop_dbl_cnt = 4

! double precision      cutcap, xcap, ycap, zcap

! common / prmtop_dbl / cutcap, xcap, ycap, zcap

! NOTE - Much of this data applies to amber forcefields used with pme or
! generalized Born simulation methods, and not the Amoeba forcefield.  Under
! Amoeba, we will read what we need of the previously defined prmtop content,
! and then get Amoeba-specific information for valence interactions (bond,
! angle, dihedral) and vdw interactions.  The charges are not even used...

  ! atm_qterm = the atom partial charge array for amber pme, and the
  !             sq_polinv array for Amoeba.  We will allocate/assign values to
  !             this later for Amoeba; we put this data here because it is
  !             roughly analogous to charge and has to be used in calcs at
  !             roughly equivalent points; we basically load up img%qterm from
  !             this array without having to check if we are in amber pme or
  !             Amoeba...

  ! atm_mass = the atom mass array.
  ! atm_iac = TBS
  ! typ_ico = TBS
  ! atm_nsp = the atom submolecule index array.
  ! atm_igraph = the atom name array.
  ! gbl_res_atms = The residue atoms index array.
  ! gbl_labres = the residue labels array.

  ! Atom and residue physical consts, names, indices into parameters,
  ! parameters, etc.

  double precision,     allocatable, save       :: atm_qterm(:)   ! Amoeba
  double precision,     allocatable, save       :: atm_mass(:)    ! Amoeba
  integer,              allocatable, save       :: atm_iac(:)     ! Amoeba
  integer,              allocatable, save       :: typ_ico(:)     ! Amoeba
  double precision,     allocatable, save       :: gbl_cn1(:)
  double precision,     allocatable, save       :: gbl_cn2(:)
  double precision,     allocatable, save       :: gbl_rk(:)
  double precision,     allocatable, save       :: gbl_req(:)
  double precision,     allocatable, save       :: gbl_tk(:)
  double precision,     allocatable, save       :: gbl_teq(:)
  double precision,     allocatable, save       :: gbl_asol(:)
  double precision,     allocatable, save       :: gbl_bsol(:)
  double precision,     allocatable, save       :: gbl_pk(:)
  double precision,     allocatable, save       :: gbl_pn(:)
  double precision,     allocatable, save       :: gbl_phase(:)

  ! Derived from atom paramaters:

  double precision,     allocatable, save       :: gbl_gamc(:)
  double precision,     allocatable, save       :: gbl_gams(:)
  double precision,     allocatable, save       :: gbl_fmn(:) 
  integer,              allocatable, save       :: gbl_ipn(:)

  character(4),         allocatable, save       :: atm_igraph(:)   ! Amoeba
  integer,              allocatable, save       :: gbl_res_atms(:) ! Amoeba
  character(4),         allocatable, save       :: gbl_labres(:)   ! Amoeba
  integer,              allocatable, save       :: atm_nsp(:)      ! Amoeba(?)

  ! For Generalized Born:

  double precision,     allocatable, save       :: atm_gb_fs(:)
  double precision,     allocatable, save       :: atm_gb_radii(:)

  ! This is the original atom bond, angle, and dihedral information.
  ! It is used at initialization and when atoms are reassigned to make
  ! the arrays that are actually used in the bonds, angles, and dihedrals
  ! modules.

  type(bond_rec),       allocatable, save       :: gbl_bond(:)
  type(angle_rec),      allocatable, save       :: gbl_angle(:)
  type(dihed_rec),      allocatable, save       :: gbl_dihed(:)

  ! NOTE! These are deallocated after setup, unless we are using
  !       Generalized Born:

  ! atm_numex = TBS
  ! gbl_natex = TBS
  ! atm_isymbl = the atom symbol array.
  ! atm_itree = TBS

  integer,              allocatable, save       :: atm_numex(:)    ! Amoeba
  integer,              allocatable, save       :: gbl_natex(:)    ! Amoeba
  character(4),         allocatable, save       :: atm_isymbl(:)   ! Amoeba
  character(4),         allocatable, save       :: atm_itree(:)    ! Amoeba

  ! Main title from prmtop.  No need to broadcast, since only master does i/o.
  
  character(80), save   :: prmtop_ititl

! Old prmtop variables that must be printed out but that are otherwise unused:

  integer, private      :: mbona
  integer, private      :: mphia
  integer, private      :: mtheta
  integer, private      :: natyp
  integer, private      :: nhparm
  integer, private      :: nmxrs
  integer, private      :: nparm

! Copies of original values of prmtop variables, the original values of which
! need to be printed out, but which are modified prior to printout:

  integer, private      :: orig_nphih
  integer, private      :: orig_nphia

! PRMTOP type definition:

  ! Enumerations for prmtop_type

  integer, parameter                    :: prmtop_type_undefined = 0
  integer, parameter                    :: prmtop_type_nonamoeba = 1
  integer, parameter                    :: prmtop_type_amoeba = 2

  integer, save                         :: prmtop_type = prmtop_type_undefined

! Hide internal routines:

#ifdef AMOEBA
  private        alloc_amoeba_prmtop_mem
#endif /* AMOEBA */
  private        alloc_amber_prmtop_mem, &
                 duplicate_dihedrals, calc_dihedral_parms

contains

!*******************************************************************************
!
! Subroutine:  init_prmtop_dat
!
! Description:  Read the molecular topology and parameters.  What is read
!               depends on whether we are expecting an amber or an amoeba
!               forcefield.
!
!*******************************************************************************

subroutine init_prmtop_dat(num_ints, num_reals, inpcrd_natom)

  use parallel_dat_mod
  use pmemd_lib_mod
  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use nextprmtop_section_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     
  integer, intent(in)           :: inpcrd_natom

! Local variables:

  integer               :: alloc_failed
  double precision      :: dielc_sqrt
  integer               :: i
  integer               :: ic

  ! Stuff used in support of the amber 7 prmtop format:

  integer               :: errcode
  character(4)          :: isymbl
  character(80)         :: fmt
  character(80)         :: fmtin
  character(80)         :: afmt, ifmt, rfmt
  character(80)         :: type

! Other unused prmtop variables:

  integer               :: ifpert, mbper, mdper, mgper, nbper, ndper, ngper

! BEGIN DBG
! integer               :: nxt_atm
! END DBG       

! Initialize stuff used in support of the amber 7 prmtop format:

  afmt = '(20A4)'
  ifmt = '(12I6)'
  rfmt = '(5E16.8)'

! Formatted input:

  call amopen(prmtop, prmtop_name, 'O', 'F', 'R')
  call nxtsec_reset()   ! insurance; the nextprmtop stuff caches the fileptr.

  fmtin = '(a80)'
  type = 'TITLE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  ! NOTE the hack below (copied from sander) that prevents the hollerith format
  ! returned from the prmtop from screwing things up (the read format is
  ! fmtin, not the returned fmt).

  read(prmtop, fmtin) prmtop_ititl

  fmtin = ifmt
  type = 'POINTERS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  numextra = 0
  ncopy = 0

  read(prmtop, fmt) natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, &
               nhparm, nparm, next, nres, nbona, ntheta, nphia, &
               numbnd, numang, nptra, natyp, nphb, ifpert, nbper, ngper, &
               ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap, numextra, &
               ncopy

  if (natom .ne. inpcrd_natom) then
    write(mdout, '(a,a)') error_hdr, &
      'natom mismatch in inpcrd/restrt and prmtop files!'
    call mexit(6, 1)
  end if

  if (ifpert .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files with perturbation information!'
    call mexit(6, 1)
  end if

  if (nparm .eq. 1) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files created by ADDLES with nparm=1!'
    write(mdout, '(a,a)') extra_line_hdr, &
      'Please use sander compiled with -DLES.'
    call mexit(6, 1)
  end if

  if (numextra .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files with extra points data!'
    write(mdout, '(a)') use_sander
    call mexit(6, 1)
  end if

  if (mbona .ne. nbona .or. mtheta .ne. ntheta .or. mphia .ne. nphia) then
    write(mdout,'(a,a)') error_hdr, &
      'PMEMD no longer allows constraints in prmtop!'
    write(mdout,'(a,a)') extra_line_hdr, &
      'The prmtop must have mbona=nbona, mtheta=ntheta, mphia=nphia.'
    call mexit(6,1)
  end if

  nttyp = ntypes * (ntypes + 1) / 2

! Here we set nscm and ndfmin. We need natom to do this...

  if (nscm .gt. 0 .and. ntb .eq. 0) then
    ndfmin = 6   ! both translation and rotation com motion removed
    if (natom == 1) ndfmin = 3
    if (natom == 2) ndfmin = 5
  else if (nscm .gt. 0) then
    ndfmin = 3    ! just translation com will be removed
  else
    ndfmin = 0
  end if
  if (ibelly .gt. 0) then   ! No COM Motion Removal, ever.
    nscm = big_int
    ndfmin = 0
  end if
  
  if (nscm .le. 0) nscm = big_int

  ! For Langevin Dynamics...

  if (gamma_ln .gt. 0.0d0) ndfmin = 0 ! No COM motion removal for LD simulation

  ! Enter amber ff or amoeba ff-specific code...

#ifdef AMOEBA
  if (iamoeba .eq. 0) then
    call init_amber_prmtop_dat
  else
    call init_amoeba_prmtop_dat
  end if
#else
  call init_amber_prmtop_dat
#endif /* AMOEBA */


  ! Reset cap data if necessary:

  if (ivcap .eq. 2) ifcap = 0

  if (ifcap .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files with a solvent cap (ifcap != 0)!'
    call mexit(6, 1)
  end if

  ! Check for single H residues, which are archaic, a pain, and no longer
  ! supported.

  do i = 1, nres
    if (gbl_res_atms(i + 1) - gbl_res_atms(i) .eq. 1) then
      if (atm_mass(gbl_res_atms(i)) .lt. 2.d0) then
        write(mdout, '(a,a)') error_hdr, &
          'PMEMD does not support single H residues!'
        write(mdout, '(a,a,i6,a,i6)') extra_line_hdr, &
          'Atom number ', gbl_res_atms(i), ', Residue number ', i
        call mexit(6, 1)
      end if
    end if
  end do

!BEGIN DBG
! nxt_atm = 1
! do i = 1, nspm
!   if (atm_nsp(i) .gt. 3) then
!     write(0,2000)'DBG: molecule, atm_cnt, atoms: ', &
!       i, atm_nsp(i), nxt_atm, '-', nxt_atm + atm_nsp(i) - 1
!   end if
!   nxt_atm = nxt_atm + atm_nsp(i)
! end do

! do i = 1, nres
!   if (gbl_res_atms(i+1) - gbl_res_atms(i) .gt. 3) then
!     write(0,2000)'DBG: res, atm_cnt, atoms: ', &
!       i, gbl_res_atms(i+1) - gbl_res_atms(i), &
!       gbl_res_atms(i),'-',gbl_res_atms(i+1) - 1
!   end if
! end do

!2000 format(a, 3i7, a1, i7)
!END DBG

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  init_amber_prmtop_dat
!
! Description:  Amber forcefield-specific prmtop processing.
!
!*******************************************************************************

subroutine init_amber_prmtop_dat

  implicit none

! Local variables:

  double precision      :: dumd                 ! scratch or trash variable
  integer               :: idum                 ! scratch or trash variable

  ! BEWARE! The bond, angle, and dihedral arrays get repacked, and the amount
  ! of space used changes.  For this reason, we store the sizes allocated
  ! by the master, and non-master processes just use those sizes.  Also, the
  ! allocated sizes are used to broadcast these arrays.  It is a little
  ! wasteful, but the safest approach for now.  The values for nbonh, nbona,
  ! ntheth, ntheta, nphih, and nphia are all the values read from the prmtop
  ! file; they may later be adjusted in eliminating contstrained bonds or
  ! handling dihedral duplicates.

  bonda_idx = nbonh + 1                   ! Bond array non-hydrogens idx.
  anglea_idx = ntheth + 1                 ! Angle array non-hydrogens idx.
  diheda_idx = nphih + max_dihed_dups + 1 ! Dihedral array non-hydrogens idx.

  gbl_bond_allocsize = nbonh + nbona
  gbl_angle_allocsize = ntheth + ntheta
  gbl_dihed_allocsize = nphih + nphia + 2 * max_dihed_dups

  call alloc_amber_prmtop_mem(num_ints, num_reals)

! Read the various prmtop arrays, making mods as required:

  fmtin = afmt
  type = 'ATOM_NAME'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_igraph(i), i = 1, natom)

  fmtin = rfmt
  type = 'CHARGE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_qterm(i), i = 1, natom)

! Scale the charges if dielc .gt. 1.0d0:

  if (using_pme_potential .and. dielc .ne. 1.0d0) then
    dielc_sqrt = sqrt(dielc)
    do i = 1, natom
      atm_qterm(i) = atm_qterm(i) / dielc_sqrt
    end do
  end if

  fmtin = rfmt
  type = 'MASS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_mass(i), i = 1, natom)

  fmtin = ifmt
  type = 'ATOM_TYPE_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_iac(i), i = 1, natom)
  
  fmtin = ifmt
  type = 'NUMBER_EXCLUDED_ATOMS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_numex(i), i = 1, natom)

  fmtin = ifmt
  type = 'NONBONDED_PARM_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (typ_ico(i), i = 1, ntypes * ntypes)

  fmtin = afmt
  type = 'RESIDUE_LABEL'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_labres(i), i = 1, nres)

  fmtin = ifmt
  type = 'RESIDUE_POINTER'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_res_atms(i), i = 1, nres)

  gbl_res_atms(nres + 1) = natom + 1

! Read the parameters:

  ! Throw away the solty array.  It is currently not used.

  fmtin = rfmt
  type = 'BOND_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_rk(i),    i = 1, numbnd)

  fmtin = rfmt
  type = 'BOND_EQUIL_VALUE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_req(i),   i = 1, numbnd)

  fmtin = rfmt
  type = 'ANGLE_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_tk(i),    i = 1, numang)

  fmtin = rfmt
  type = 'ANGLE_EQUIL_VALUE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_teq(i),   i = 1, numang)
  
  fmtin = rfmt
  type = 'DIHEDRAL_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_pk(i),    i = 1, nptra)

  fmtin = rfmt
  type = 'DIHEDRAL_PERIODICITY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_pn(i),    i = 1, nptra)

  fmtin = rfmt
  type = 'DIHEDRAL_PHASE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_phase(i), i = 1, nptra)

  fmtin = rfmt
  type = 'SOLTY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (dumd, i = 1, natyp)           ! solty

  fmtin = rfmt
  type = 'LENNARD_JONES_ACOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_cn1(i),   i = 1, nttyp)
  
  fmtin = rfmt
  type = 'LENNARD_JONES_BCOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_cn2(i),   i = 1, nttyp)

! Read the bonding information:

  fmtin = ifmt
  type = 'BONDS_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bond(i)%atm_i, &
                gbl_bond(i)%atm_j, &
                gbl_bond(i)%parm_idx, &
                i = 1, nbonh)

  do i = 1, nbonh
    gbl_bond(i)%atm_i = gbl_bond(i)%atm_i / 3 + 1
    gbl_bond(i)%atm_j = gbl_bond(i)%atm_j / 3 + 1
  end do

  fmtin = ifmt
  type = 'BONDS_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bond(bonda_idx + i)%atm_i, &
                gbl_bond(bonda_idx + i)%atm_j, &
                gbl_bond(bonda_idx + i)%parm_idx, &
                i = 0, nbona - 1)

  do i = 0, nbona - 1
    gbl_bond(bonda_idx + i)%atm_i = gbl_bond(bonda_idx + i)%atm_i / 3 + 1
    gbl_bond(bonda_idx + i)%atm_j = gbl_bond(bonda_idx + i)%atm_j / 3 + 1
  end do

! Read the angle information:

  fmtin = ifmt
  type = 'ANGLES_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_angle(i)%atm_i, &
                gbl_angle(i)%atm_j, &
                gbl_angle(i)%atm_k, &
                gbl_angle(i)%parm_idx, &
                i = 1, ntheth)

  do i = 1, ntheth
    gbl_angle(i)%atm_i = gbl_angle(i)%atm_i / 3 + 1
    gbl_angle(i)%atm_j = gbl_angle(i)%atm_j / 3 + 1
    gbl_angle(i)%atm_k = gbl_angle(i)%atm_k / 3 + 1
  end do

  fmtin = ifmt
  type = 'ANGLES_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_angle(anglea_idx + i)%atm_i, &
                gbl_angle(anglea_idx + i)%atm_j, &
                gbl_angle(anglea_idx + i)%atm_k, &
                gbl_angle(anglea_idx + i)%parm_idx, &
                i = 0, ntheta - 1)

  do i = 0, ntheta - 1
    gbl_angle(anglea_idx + i)%atm_i = gbl_angle(anglea_idx + i)%atm_i / 3 + 1
    gbl_angle(anglea_idx + i)%atm_j = gbl_angle(anglea_idx + i)%atm_j / 3 + 1
    gbl_angle(anglea_idx + i)%atm_k = gbl_angle(anglea_idx + i)%atm_k / 3 + 1
  end do

! Read the dihedral information:

  fmtin = ifmt
  type = 'DIHEDRALS_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_dihed(i)%atm_i, &
                gbl_dihed(i)%atm_j, &
                gbl_dihed(i)%atm_k, &
                gbl_dihed(i)%atm_l, &
                gbl_dihed(i)%parm_idx, &
                i = 1, nphih)

  do i = 1, nphih
    gbl_dihed(i)%atm_i = gbl_dihed(i)%atm_i / 3 + 1
    gbl_dihed(i)%atm_j = gbl_dihed(i)%atm_j / 3 + 1
    if (gbl_dihed(i)%atm_k .gt. 0) then
      gbl_dihed(i)%atm_k = gbl_dihed(i)%atm_k / 3 + 1
    else
      gbl_dihed(i)%atm_k = gbl_dihed(i)%atm_k / 3 - 1
    end if
    if (gbl_dihed(i)%atm_l .gt. 0) then
      gbl_dihed(i)%atm_l = gbl_dihed(i)%atm_l / 3 + 1
    else
      gbl_dihed(i)%atm_l = gbl_dihed(i)%atm_l / 3 - 1
    end if
  end do

  fmtin = ifmt
  type = 'DIHEDRALS_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_dihed(diheda_idx + i)%atm_i, &
                gbl_dihed(diheda_idx + i)%atm_j, &
                gbl_dihed(diheda_idx + i)%atm_k, &
                gbl_dihed(diheda_idx + i)%atm_l, &
                gbl_dihed(diheda_idx + i)%parm_idx, &
                i = 0, nphia - 1)

  do i = 0, nphia - 1
    gbl_dihed(diheda_idx + i)%atm_i = gbl_dihed(diheda_idx + i)%atm_i / 3 + 1
    gbl_dihed(diheda_idx + i)%atm_j = gbl_dihed(diheda_idx + i)%atm_j / 3 + 1
    if (gbl_dihed(diheda_idx + i)%atm_k .gt. 0) then
      gbl_dihed(diheda_idx + i)%atm_k = gbl_dihed(diheda_idx + i)%atm_k / 3 + 1
    else
      gbl_dihed(diheda_idx + i)%atm_k = gbl_dihed(diheda_idx + i)%atm_k / 3 - 1
    end if
    if (gbl_dihed(diheda_idx + i)%atm_l .gt. 0) then
      gbl_dihed(diheda_idx + i)%atm_l = gbl_dihed(diheda_idx + i)%atm_l / 3 + 1
    else
      gbl_dihed(diheda_idx + i)%atm_l = gbl_dihed(diheda_idx + i)%atm_l / 3 - 1
    end if
  end do

  fmtin = ifmt
  type = 'EXCLUDED_ATOMS_LIST'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_natex(i), i = 1, next)

! Read the h-bond parameters:

  fmtin = rfmt
  type = 'HBOND_ACOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_asol(i), i = 1, nphb)

  fmtin = rfmt
  type = 'HBOND_BCOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bsol(i), i = 1, nphb)

  fmtin = rfmt
  type = 'HBCUT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (dumd, i = 1, nphb)    ! Throw away hbcut. It is not used.

! Read isymbl, itree, join and irotat arrays.

  fmtin = afmt
  type = 'AMBER_ATOM_TYPE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_isymbl(i), i = 1, natom)

  fmtin = afmt
  type = 'TREE_CHAIN_CLASSIFICATION'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_itree(i), i = 1, natom)

  fmtin = ifmt
  type = 'JOIN_ARRAY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (idum, i = 1, natom) ! Throw away the join array. Not used.

  fmtin = ifmt
  type = 'IROTAT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (idum, i = 1, natom)! Throw away the irotat array. Not used.
  
! Read the boundary condition stuff if needed:

  if (ifbox .gt. 0) then

    fmtin = ifmt
    type = 'SOLVENT_POINTERS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) idum, nspm, idum

  else

    nspm = 1

  end if

  allocate(atm_nsp(nspm), stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_nsp)

  ! We read the prmtop box and angle data here but ignore it, per standard
  ! amber behaviour.  The ifbox value also influences atm_nsp() setup though.

  if (ifbox .gt. 0) then

    fmtin = ifmt
    type = 'ATOMS_PER_MOLECULE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (atm_nsp(i), i = 1, nspm)

    fmtin = rfmt
    type = 'BOX_DIMENSIONS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) dumd, dumd, dumd, dumd

  else

    atm_nsp(1) = natom

  end if

! Load the cap information if needed:

  if (ifcap .gt. 0) then

    fmtin = '(I6)'
    type = 'CAP_INFO'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) idum   ! natcap

    fmtin = '(4E16.8)'
    type = 'CAP_INFO2'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    ! At least for now, cap data is thrown away.  As far as I know, all
    ! cap simulations do not have PBC.

    read(prmtop, fmt) dumd, dumd, dumd, dumd ! cutcap, xcap, ycap, zcap

! else

!   natcap = 0
!   cutcap = 0.0d0
!   xcap = 0.0d0
!   ycap = 0.0d0
!   zcap = 0.0d0

  end if

! Generalized Born, if appropriate:

  if (using_gb_potential) then

    if (errcode .eq. -1) then
      write(mdout, '(a,a)') error_hdr, &
        'GB calculations now require a new-style prmtop file'
      call mexit(6, 1)
    end if

    fmtin = rfmt
    type = 'RADII'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (atm_gb_radii(i), i = 1, natom)

    fmtin = rfmt
    type = 'SCREEN'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (atm_gb_fs(i), i = 1, natom)

  end if

! Polarizabilities would need to be read here if they are supported.

! Fix up the typ-ico array for more efficient processing of vdw and 10_12
! interactions.  Basically, if the ico array entry is 0, then vdw processing
! is unnecessary for amber pme.

  do i = 1, ntypes * ntypes
    ic = typ_ico(i)
    if (ic .gt. 0) then
      if (gbl_cn1(ic) .eq. 0.d0 .and. gbl_cn2(ic) .eq. 0.d0) typ_ico(i) = 0
    else if (ic .lt. 0) then
      ic = -ic
      if (gbl_asol(ic) .eq. 0.d0 .and. gbl_bsol(ic) .eq. 0.d0) typ_ico(i) = 0
    end if
  end do

! Duplicate dihedral pointers for vector ephi.  nphih and nphia may be increased
! in these calls.

  orig_nphih = nphih        ! Saved for reporting later
  orig_nphia = nphia        ! Saved for reporting later

  call duplicate_dihedrals(nphih, gbl_dihed, gbl_pn)
  call duplicate_dihedrals(nphia, gbl_dihed(diheda_idx), gbl_pn)

! Pre-calculate and save some parameters for vector ephi:

  call calc_dihedral_parms(nptra, gbl_pk, gbl_pn, gbl_phase, gbl_gamc, &
                           gbl_gams, gbl_ipn, gbl_fmn)

! Print prmtop parameters header ('resource use'):

  write(mdout, 1020)

  ! This is after-the-fact, but done here for output consistency with sander 8:

  if (ntb .gt. 0) then
    write(mdout, '(a, /)')' getting new box info from bottom of inpcrd'
  end if

! Print prmtop parameters:

  write(mdout, 1030) natom, ntypes, nbonh, mbona, ntheth, mtheta, &
                     orig_nphih, mphia, nhparm, nparm, next, &
                     nres, nbona, ntheta, orig_nphia, numbnd, numang, &
                     nptra, natyp, nphb, ifbox, nmxrs, ifcap, numextra, ncopy

  if (using_gb_potential) then

    fmtin = afmt
    type = 'RADIUS_SET'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    ! Write implicit solvent radius and screening info to mdout:

    if (errcode .eq. 0) then
      read(prmtop, fmt) type    ! note re-use of type to hold result, a hack...
      write(mdout, '(a,a)') &
        ' Implicit solvent radii are ', type
    end if

    ! If igb .eq. 7 use special S_x screening params; here we overwrite
    ! the tinker values read from the prmtop.

    if (igb .eq. 7) then

      write(mdout,'(a)') &
        ' Replacing prmtop screening parameters with GBn (igb=7) values'

      do i = 1, natom

        isymbl = atm_isymbl(i)

        if (isymbl(1:1) .eq. 'C' .or. isymbl(1:1) .eq. 'c') then
          atm_gb_fs(i) = 4.84353823306d-1
        else if (isymbl(1:1) .eq. 'H' .or. isymbl(1:1) .eq. 'h') then
          atm_gb_fs(i) = 1.09085413633d0
        else if (isymbl(1:1) .eq. 'N' .or. isymbl(1:1) .eq. 'n') then
          atm_gb_fs(i) = 7.00147318409d-1
        else if (isymbl(1:1) .eq. 'O' .or. isymbl(1:1) .eq. 'o') then
          atm_gb_fs(i) = 1.06557401132d0
        else if (isymbl(1:1) .eq. 'S' .or. isymbl(1:1) .eq. 's') then
          atm_gb_fs(i) = 6.02256336067d-1
        else
          atm_gb_fs(i) = 5.d-1
        end if

      end do

    end if

    ! Put fs(i) * (rborn(i) - offset) into the "fs" array:

    gb_fs_max = 0.d0

    do i = 1, natom
      atm_gb_fs(i) = atm_gb_fs(i) * (atm_gb_radii(i) - offset)
      gb_fs_max = max(gb_fs_max, atm_gb_fs(i))
    end do

  end if

  close(unit = prmtop)

  return

 1020 format(80('-')/'   1.  RESOURCE   USE: ',/80('-')/)

 1030 format(t2, &
        'NATOM  = ', i7, ' NTYPES = ', i7, ' NBONH = ', i7, ' MBONA  = ', i7, &
      /' NTHETH = ', i7, ' MTHETA = ', i7, ' NPHIH = ', i7, ' MPHIA  = ', i7, &
      /' NHPARM = ', i7, ' NPARM  = ', i7, ' NNB   = ', i7, ' NRES   = ', i7, &
      /' NBONA  = ', i7, ' NTHETA = ', i7, ' NPHIA = ', i7, ' NUMBND = ', i7, &
      /' NUMANG = ', i7, ' NPTRA  = ', i7, ' NATYP = ', i7, ' NPHB   = ', i7, &
      /' IFBOX  = ', i7, ' NMXRS  = ', i7, ' IFCAP = ', i7, ' NEXTRA = ', i7, &
      /' NCOPY  = ', i7/)

end subroutine init_amber_prmtop_dat

#ifdef AMOEBA
!*******************************************************************************
!
! Subroutine:  init_amoeba_prmtop_dat
!
! Description:  Amoeba forcefield-specific prmtop processing.
!
!*******************************************************************************

subroutine init_amoeba_prmtop_dat

  implicit none

! Local variables:

  double precision      :: dumd                 ! scratch or trash variable
  integer               :: idum                 ! scratch or trash variable

  if (ntypes .ne. 1) then
    write(mdout, '(a,a)') error_hdr, &
      'Amoeba code expects prmtop to specify only one Amber atom type (ntypes)!'
    call mexit(6, 1)
  end if

  ! All the amber bond, angle, and dihedral info is not used by Amoeba; it has
  ! its own valence terms...

  bonda_idx = 0
  anglea_idx = 0
  diheda_idx = 0

  gbl_bond_allocsize = 0
  gbl_angle_allocsize = 0
  gbl_dihed_allocsize = 0

  call alloc_amoeba_prmtop_mem(num_ints, num_reals)

! Read the various prmtop arrays, making mods as required:

  fmtin = afmt
  type = 'ATOM_NAME'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_igraph(i), i = 1, natom)

  ! The atm_qterm array is allocated in code in this module and will be
  ! filled in later.  We may want to refactor this at some point...

  fmtin = rfmt
  type = 'MASS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_mass(i), i = 1, natom)

  ! There is only one amber atom type in amoeba, but we are dependent
  ! on atm_iac for pairlist processing at the moment...

  fmtin = ifmt
  type = 'ATOM_TYPE_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_iac(i), i = 1, natom)
  
  fmtin = ifmt
  type = 'NUMBER_EXCLUDED_ATOMS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_numex(i), i = 1, natom)

  ! There is only one amber atom type in amoeba, but we are dependent
  ! on typ_ico for pairlist processing at the moment...

  fmtin = ifmt
  type = 'NONBONDED_PARM_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (typ_ico(i), i = 1, ntypes * ntypes)

  fmtin = afmt
  type = 'RESIDUE_LABEL'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_labres(i), i = 1, nres)

  fmtin = ifmt
  type = 'RESIDUE_POINTER'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_res_atms(i), i = 1, nres)

  gbl_res_atms(nres + 1) = natom + 1

  ! Force the Amber bond, angle, dihedral param counts to 0...

  numbnd = 0
  numang = 0
  nptra = 0

  ! We don't use the amber LJ coefficients in Amoeba...

  ! Force the Amber bond information counts to 0...

  nbonh = 0
  nbona = 0
  mbona = 0

  ! Force the Amber angle information counts to 0...

  ntheth = 0
  ntheta = 0
  mtheta = 0

  ! Force the Amber dihedral information counts to 0...

  nphih = 0
  nphia = 0
  mphia = 0

  fmtin = ifmt
  type = 'EXCLUDED_ATOMS_LIST'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_natex(i), i = 1, next)

! We don't use the amber HBOND coefficients in Amoeba...

  nphb = 0

! Read isymbl, itree, join and irotat arrays.

  fmtin = afmt
  type = 'AMBER_ATOM_TYPE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_isymbl(i), i = 1, natom)

  fmtin = afmt
  type = 'TREE_CHAIN_CLASSIFICATION'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_itree(i), i = 1, natom)

! Read the boundary condition stuff if needed:

  if (ifbox .gt. 0) then
    fmtin = ifmt
    type = 'SOLVENT_POINTERS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop, fmt) idum, nspm, idum
  else
    nspm = 1
  end if

  allocate(atm_nsp(nspm), stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_nsp)

  ! We read the prmtop box and angle data here if appropriate.  I think amber
  ! 6 and 7 both ignore it, but there are scenarios (with binary restart files
  ! I believe) where there is no other place for the info to come from, so
  ! PMEMD will look at it as a last resort.

  ! The ifbox value only influences atm_nsp() setup in pmemd...

  if (ifbox .gt. 0) then
    fmtin = ifmt
    type = 'ATOMS_PER_MOLECULE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop, fmt) (atm_nsp(i), i = 1, nspm)
  else
    atm_nsp(1) = natom
  end if

! We skip the cap info - not used for explicit solvent simulations...

! For amoeba we force the typ_ico value (only 1 amber atom type, theoretically,
! since the typing is handled in amoeba) to 0 also to short-circuit some
! unnnecessary code.  This is a hack...

  do i = 1, ntypes * ntypes
    typ_ico(i) = 0
  end do


  ! Here we dink with amber dihedrals if dealing with an amber ff. No need
  ! under Amoeba.  Report 0's for the dihedral cnts for now.

  orig_nphih = 0        ! Saved for reporting later
  orig_nphia = 0        ! Saved for reporting later

! Print prmtop parameters header ('resource use'):

  write(mdout, 1020)

  ! This is after-the-fact, but done here for output consistency with sander 8:

  if (ntb .gt. 0) then
    write(mdout, '(a, /)')' getting new box info from bottom of inpcrd'
  end if

! Print prmtop parameters.  This needs to change for Amoeba...

  write(mdout, 1030) natom, ntypes, nbonh, mbona, ntheth, mtheta, &
                     orig_nphih, mphia, nhparm, nparm, next, &
                     nres, nbona, ntheta, orig_nphia, numbnd, numang, &
                     nptra, natyp, nphb, ifbox, nmxrs, ifcap, numextra, ncopy


  close(unit = prmtop)

  return

 1020 format(80('-')/'   1.  RESOURCE   USE: ',/80('-')/)

 1030 format(t2, &
        'NATOM  = ', i7, ' NTYPES = ', i7, ' NBONH = ', i7, ' MBONA  = ', i7, &
      /' NTHETH = ', i7, ' MTHETA = ', i7, ' NPHIH = ', i7, ' MPHIA  = ', i7, &
      /' NHPARM = ', i7, ' NPARM  = ', i7, ' NNB   = ', i7, ' NRES   = ', i7, &
      /' NBONA  = ', i7, ' NTHETA = ', i7, ' NPHIA = ', i7, ' NUMBND = ', i7, &
      /' NUMANG = ', i7, ' NPTRA  = ', i7, ' NATYP = ', i7, ' NPHB   = ', i7, &
      /' IFBOX  = ', i7, ' NMXRS  = ', i7, ' IFCAP = ', i7, ' NEXTRA = ', i7, &
      /' NCOPY  = ', i7/)

end subroutine init_amoeba_prmtop_dat
#endif /* AMOEBA */

end subroutine init_prmtop_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_amber_prmtop_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_amber_prmtop_dat

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

  integer               :: alloc_failed
  integer               :: bytes_per_unit
  integer               :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(natom, prmtop_int_cnt, mpi_integer, 0, mpi_comm_world, &
                 err_code_mpi)

! There is only cap data, which we currently don't use, so we comment this
! out:
!
! call mpi_bcast(cutcap, prmtop_dbl_cnt, mpi_double_precision, 0, &
!                mpi_comm_world, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_amber_prmtop_mem(num_ints, num_reals)
  end if

  call mpi_bcast(atm_qterm, natom, mpi_double_precision, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_mass, natom, mpi_double_precision, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_iac, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(typ_ico, ntypes * ntypes, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_res_atms, nres + 1, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_igraph, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_labres, nres, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  if (nphb > 0) then
    call mpi_bcast(gbl_asol, nphb, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_bsol, nphb, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  if (nttyp > 0) then
    call mpi_bcast(gbl_cn1, nttyp, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_cn2, nttyp, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  if (numbnd > 0) then
    call mpi_bcast(gbl_rk, numbnd, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_req, numbnd, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  if (numang > 0) then
    call mpi_bcast(gbl_tk, numang, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_teq, numang, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  if (nptra > 0) then
    call mpi_bcast(gbl_pk, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_pn, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_phase, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_gamc, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_gams, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_fmn, nptra, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_ipn, nptra, mpi_integer, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  ! Generalized Born:

  if (using_gb_potential) then
    call mpi_bcast(atm_gb_fs, natom, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(atm_gb_radii, natom, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
  end if

  ! Polarizabilities would need to be broadcast if they are supported.

  call get_bytesize(gbl_bond(1), gbl_bond(2), bytes_per_unit)

  call mpi_bcast(gbl_bond, size(gbl_bond) * bytes_per_unit, mpi_byte, 0, &
                 mpi_comm_world, err_code_mpi)

  call get_bytesize(gbl_angle(1), gbl_angle(2), bytes_per_unit)

  call mpi_bcast(gbl_angle, size(gbl_angle) * bytes_per_unit, mpi_byte, 0, &
                 mpi_comm_world, err_code_mpi)

  call get_bytesize(gbl_dihed(1), gbl_dihed(2), bytes_per_unit)

  call mpi_bcast(gbl_dihed, size(gbl_dihed) * bytes_per_unit, mpi_byte, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_numex, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_natex, next, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  ! No need to broadcast atm_isymbl or atm_itree

  return

end subroutine bcast_amber_prmtop_dat
#endif

#ifdef MPI
#ifdef AMOEBA
!*******************************************************************************
!
! Subroutine:  bcast_amoeba_prmtop_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_amoeba_prmtop_dat

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

  integer               :: alloc_failed
  integer               :: bytes_per_unit
  integer               :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(natom, prmtop_int_cnt, mpi_integer, 0, mpi_comm_world, &
                 err_code_mpi)

! There is only cap data, which we currently don't use, so we comment this
! out:
!
! call mpi_bcast(cutcap, prmtop_dbl_cnt, mpi_double_precision, 0, &
!                mpi_comm_world, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_amoeba_prmtop_mem(num_ints, num_reals)
  end if

  call mpi_bcast(atm_mass, natom, mpi_double_precision, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_iac, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(typ_ico, ntypes * ntypes, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_res_atms, nres + 1, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_igraph, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_labres, nres, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(atm_numex, natom, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(gbl_natex, next, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  ! No need to broadcast atm_isymbl or atm_itree

  return

end subroutine bcast_amoeba_prmtop_dat
#endif /* AMOEBA */
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  alloc_amber_prmtop_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_amber_prmtop_mem(num_ints, num_reals)

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_qterm(natom), &
           atm_mass(natom), &
           atm_iac(natom), &
           typ_ico(ntypes * ntypes), &
           gbl_cn1(nttyp), &
           gbl_cn2(nttyp), &
           gbl_rk(numbnd), &
           gbl_req(numbnd), &
           gbl_tk(numang), &
           gbl_teq(numang), &
           gbl_res_atms(nres + 1), &
           atm_igraph(natom), &
           gbl_labres(nres), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_qterm) + size(atm_mass) + &
                          size(gbl_cn1) + size(gbl_cn2) + size(gbl_rk) + &
                          size(gbl_req) + size(gbl_tk) + size(gbl_teq)

  num_ints = num_ints + size(atm_iac) + size(typ_ico) + &
                       size(gbl_res_atms) + size(atm_igraph) + size(gbl_labres)

  if (nphb .ne. 0) then
    allocate(gbl_asol(nphb), &
             gbl_bsol(nphb), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_asol) + size(gbl_bsol)

  end if

  ! Allocate space for atom bond, angle and dihedral information:

  allocate(gbl_bond(gbl_bond_allocsize), &
           gbl_angle(gbl_angle_allocsize), &
           gbl_dihed(gbl_dihed_allocsize), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_bond) * bond_rec_ints + &
             size(gbl_angle) * angle_rec_ints + &
             size(gbl_dihed) * dihed_rec_ints

  if (nptra > 0) then
    allocate(gbl_pk(nptra), &
             gbl_pn(nptra), &
             gbl_phase(nptra), &
             gbl_gamc(nptra), &
             gbl_gams(nptra), &
             gbl_fmn(nptra), &
             gbl_ipn(nptra), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_pk) + &
                            size(gbl_pn) + &
                            size(gbl_phase) + &
                            size(gbl_gamc) + &
                            size(gbl_gams) + &
                            size(gbl_fmn)

    num_ints = num_ints + size(gbl_ipn)

    ! We know these get initialized properly, so no need to zero.

  end if

  ! Allocations for Generalized Born:

  if (using_gb_potential) then

    allocate(atm_gb_fs(natom), &
             atm_gb_radii(natom), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(atm_gb_fs) + size(atm_gb_fs)

  end if

  ! Polarizabilities would need to have space alloc'd if they are supported.

  ! Allocation of stuff that will later be deallocated:

  allocate(atm_numex(natom), &
           gbl_natex(next), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_numex) + size(gbl_natex)

  ! Allocation of stuff that will later be deallocated that is only used in
  ! the master:

  if (master) then

    allocate(atm_isymbl(natom), &
             atm_itree(natom), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_isymbl) + size(atm_itree)

  end if

  return

end subroutine alloc_amber_prmtop_mem

#ifdef AMOEBA
!*******************************************************************************
!
! Subroutine:  alloc_amoeba_prmtop_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_amoeba_prmtop_mem(num_ints, num_reals)

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_qterm(natom), &
           atm_mass(natom), &
           atm_iac(natom), &
           typ_ico(ntypes * ntypes), &
           gbl_res_atms(nres + 1), &
           atm_igraph(natom), &
           gbl_labres(nres), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_qterm) + size(atm_mass)

  atm_qterm(:) = 0.d0   ! Safety measure...

  num_ints = num_ints + size(atm_iac) + size(typ_ico) + &
                       size(gbl_res_atms) + size(atm_igraph) + size(gbl_labres)

  ! Allocation of stuff that will later be deallocated:

  allocate(atm_numex(natom), &
           gbl_natex(next), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_numex) + size(gbl_natex)

  ! Allocation of stuff that will later be deallocated that is only used in
  ! the master:

  if (master) then

    allocate(atm_isymbl(natom), &
             atm_itree(natom), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_isymbl) + size(atm_itree)

  end if

  return

end subroutine alloc_amoeba_prmtop_mem
#endif /* AMOEBA */

#ifdef AMOEBA
!*******************************************************************************!
! Subroutine:  get_prmtop_type
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_prmtop_type

  use file_io_dat_mod
  use file_io_mod
  use nextprmtop_section_mod

  implicit none

! Local variables:

  integer               :: iok
  integer               :: ionerr
  character(len = 80)   :: fmt
  character(len = 80)   :: fmtin
  character(len = 80)   :: dtype

  call amopen(prmtop, prmtop_name, 'O', 'F', 'R')
  call nxtsec_reset()   ! insurance; the nextprmtop stuff caches the fileptr.

  dtype = 'AMOEBA_FORCEFIELD'
  fmtin = '(I5)'
  ionerr = 1    ! not fatal if missing

  call nxtsec(prmtop, mdout, ionerr, fmtin, dtype, fmt, iok)

  if (iok .eq. 0) then          !this data type found in prmtop
    prmtop_type = prmtop_type_amoeba
  else
    prmtop_type = prmtop_type_nonamoeba
  end if

  close(unit = prmtop)

  return

end subroutine get_prmtop_type
#endif /* AMOEBA */

!*******************************************************************************
!
! Subroutine:  duplicate_dihedrals
!
! Description:  Duplicates pointers to multi-term dihedrals for vector ephi.
!               H-atom diheds are duplicated at the end of the h-atom lists,
!               but heavy atom diheds must be inserted between the end of the
!               heavy-atom diheds.  In order to use this code, extra space MUST
!               be allocated (2 * max_dihed_dups).
!
! Author:  George Seibel
!              
!*******************************************************************************

subroutine duplicate_dihedrals(nphi, dihed, pn)

  use gbl_constants_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: nphi         ! num real dihedrals
  type(dihed_rec)       :: dihed(*)     ! the dihedral record
  double precision      :: pn(*)        ! periodicity; negative if multi-term,
                                        ! read until + encountered:
! Local variables:

  integer               :: i
  integer               :: ic           ! working parameter pointer:
  integer               :: ndup         ! num of diheds duplicated:

  ndup = 0

  do i = 1, nphi

    ic = dihed(i)%parm_idx

    do while (pn(ic) .lt. 0)

      ndup = ndup + 1

      if (ndup .gt. max_dihed_dups) then

        write(mdout,'(a,a,i5,a)') error_hdr, &
          'max_dihed_dups = ',max_dihed_dups,' exceeded.'

        write(mdout,'(a,a,i5)') extra_line_hdr, &
          'Set max_dihed_dups =', ndup
        call mexit(6, 1)

      else 

     ! Duplicate pointers of multi-term dihedrals:

        dihed(nphi + ndup)%atm_i  = dihed(i)%atm_i
        dihed(nphi + ndup)%atm_j  = dihed(i)%atm_j
        dihed(nphi + ndup)%atm_k  = dihed(i)%atm_k
        dihed(nphi + ndup)%atm_l  = dihed(i)%atm_l

     ! But use the NEXT parameter pointer:

        dihed(nphi + ndup)%parm_idx = ic + 1

      end if

      ! Check for a third or higher term:

      ic = ic + 1

    end do

  end do

  nphi = nphi + ndup

  write(mdout,'(a,i5,a,/)') '| Duplicated',ndup,' dihedrals'

  return

end subroutine duplicate_dihedrals

!*******************************************************************************
!
! Subroutine:  calc_dihedral_parms
!
! Description: Routine to calc additional parameters for the vectorised ephi.
!
!*******************************************************************************

subroutine calc_dihedral_parms(numphi, pk, pn, phase, gamc, gams, ipn, fmn)

  use gbl_constants_mod

  implicit none

! Formal arguments:

  integer           numphi
  double precision  pk(*)
  double precision  pn(*)
  double precision  phase(*)
  double precision  gamc(*)
  double precision  gams(*)
  integer           ipn(*)
  double precision  fmn(*)

! Local variables:

  double precision  dum
  double precision  dumc
  double precision  dums
  double precision  four
  integer           i
  double precision  one
  double precision  pim
  double precision  tenm3
  double precision  tenm10
  double precision  zero

  data zero, one, tenm3, tenm10/0.0d+00, 1.0d+00, 1.0d-03, 1.0d-06/
  data four /4.0d+00/

  pim = four*atan(one)

  do i = 1, numphi
    dum = phase(i)
    if (abs(dum - PI) .le. tenm3) dum = sign(pim, dum)
    dumc = cos(dum)
    dums = sin(dum)
    if (abs(dumc) .le. tenm10) dumc = zero
    if (abs(dums) .le. tenm10) dums = zero
    gamc(i) = dumc*pk(i)
    gams(i) = dums*pk(i)
    fmn(i) = one
    if (pn(i) .le. zero) fmn(i) = zero
    pn(i) = abs(pn(i))
    ipn(i) = int(pn(i) + tenm3)
  end do

  return

end subroutine calc_dihedral_parms

end module prmtop_dat_mod
