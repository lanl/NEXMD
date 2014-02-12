module cosmo_C
  !USE vast_kind_param, ONLY:  double
  implicit none
  integer :: solvent_model, potential_type, EF
  logical :: iseps, noeps, useps, lpka
  logical::coserr=.false.
  integer :: nspa, nps, nps2, nden, lenabc, nppa = 1082, &
  amat_dim, isude_dim, nipc, ioldcv
  integer, dimension (2) :: n0
   integer lm61,numat,mpack
   real(8),parameter::a0=0.5291772083d0 ! Bohr in angstroms
   real(8),parameter::ev=27.2113834d0 ! Hartree in eV

  integer, dimension(:), allocatable :: &
  & iatsp,   & !
  & nar_csm, & !
  & nsetf,i, & !  Watch out for the "i"
  & ipiden,  & !
  & idenat,  & !
  & nset       !
  integer, dimension(:,:), allocatable :: &
  & isude,   & !
  & nn
  real(8) ::&
    fepsi,       & ! Dielectric factor =  (e-1)/(e+0.5), e = dielectric constant
                   !
    rds, disex2, &
                   !
    ediel,       & ! Dielectric energy, in eV
                   !
    solv_energy, & ! Solvation energy, in eV
                   !
    area,        & ! Surface area, in square Angstroms
                   !
    cosvol,      & ! Molecular volume, in cubic Angstroms
                   !
    cif1, cif2,  & !
    cdiagi,      &
   fnsq, rsolv,  &
   cosmo_scf_ftol, &		! Tolerance in SCF COSMO procedure
   index_of_refraction, & 	! Optical index of Refraction	
   ceps, &			! Dielectric permitivity	
   onsagE, &
   Ex, &
   Ey, &
   Ez, &
   onsager_radius	

  double precision, dimension(3,3,1082) :: tm
  double precision, dimension(4,1082) :: dirsm, dirvec
  integer, dimension(:,:), allocatable :: tri_2D;
  double precision, dimension(:), allocatable :: &
  & amat,    & !
  & cmat,    & !
  & gden,    & !
  & qscat,   & !
  & arat,    & !
  & srad,    & !
  & abcmat,  & !
  & qden,    & !
  & bh,      & ! 
  & cdiag,   & !
  & xi_k
  double precision, dimension(:,:), target, allocatable :: v_solvent_difdens;
  double precision, dimension(:,:), pointer :: v_solvent_xi;	
  double precision, dimension(:,:), allocatable :: &
  & bmat,    & !
  & phinet,  & !
  & qscnet,  & !
  & qdenet,  & !
  & cosurf,  & !
  & sude,    & !
  & xsp,     & !
  & cxy,     & !
  & mmat
end module cosmo_C
