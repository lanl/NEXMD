
#include "dprec.fh"

!------------------------------------------------------------------------------
subroutine transformation_matrix(nbead)
!------------------------------------------------------------------------------

    use constants, only : pi, kB, hbar

    use cmd_vars, only : cart_to_nmode, nmode_to_cart, fcart_to_fnmode, &
        lambda_nmode, omega_nmode, mass_nmode, adiab_param

    !..................................................

    implicit none

#include "md.h"
#include "extra.h"

    integer, intent(in) :: nbead

    _REAL_ :: kT, beta, omega_P

    !..................................................

    integer :: n, n1, n2, ibead

    _REAL_ :: tmp

    !..................................................

    kT = kB*temp0

    beta = 1.d0/kT

    omega_P = sqrt(dble(nbead))/hbar/beta

    !---------------------------------------

    !! Transformation matrix from normal modes to cartesian coordinates.

    nmode_to_cart(1:nbead, 1) = 1.d0

    do n1 = 1, nbead/2
        nmode_to_cart(2*n1, nbead) = -1.0d0
        nmode_to_cart(2*n1 - 1, nbead) = 1.0d0
    end do

    do n1 = 1, nbead
        tmp = 2.0d0*pi*dble(n1 - 1)/dble(nbead)
        do n2 = 1, nbead/2 - 1
            nmode_to_cart(n1, 2*n2) = sqrt(2.d0)*cos(tmp*dble(n2))
            nmode_to_cart(n1, 2*n2 + 1) = -sqrt(2.d0)*sin(tmp*dble(n2))
        end do
    end do

    !---------------------------------------

    !! Transformation matrix from cartesian coordinates to normal modes.

    do n2 = 1, nbead
        do n1 = 1, nbead
            cart_to_nmode(n2, n1) = nmode_to_cart(n1, n2)/dble(nbead)
            fcart_to_fnmode(n2, n1) = nmode_to_cart(n1, n2)
        end do
    end do

    !---------------------------------------

    !! Normal mode frequencies. n=1 --> path centroid.

    lambda_nmode(1) = 0.0d0

    do n = 1, nbead/2 - 1
        tmp = 2.0d0*pi*dble(n)/dble(nbead)
        tmp = 2.0d0*dble(nbead)*(1.0d0 - cos(tmp))
        lambda_nmode(2*n) = tmp
        lambda_nmode(2*n + 1) = tmp
    end do

    lambda_nmode(nbead) = 4.0d0*dble(nbead)

    !---------------------------------------

    !! Characteristic oscillation frequency.

    omega_nmode(1) = omega_P

    do n = 2, nbead
        if (adiab_param < 1.d0) then
            omega_nmode(n) = omega_P/adiab_param
        else
            omega_nmode(n) = omega_P
        end if
    end do

end subroutine transformation_matrix

!------------------------------------------------------------------------------
subroutine trans_pos_cart_to_nmode(x)

! Transform cartesian positions into normal mode positions.
!------------------------------------------------------------------------------

    use pimd_vars, only : nbead, tempbuf

    use cmd_vars, only : cart_to_nmode

    !..................................................
    implicit none

#include "memory.h"
#include "les.h"

    integer :: first, last

    _REAL_  :: x(3, natom)

    !..................................................

    integer :: iatom, ibead

    !..................................................

#ifdef LES
    call trans_part(x, tempbuf, cart_to_nmode)
#else
    call trans_full(x, tempbuf, cart_to_nmode)
#endif
    x(1:3, 1:natom) = tempbuf(1:3, 1:natom)
    return

end subroutine trans_pos_cart_to_nmode

!------------------------------------------------------------------------------
subroutine trans_pos_nmode_to_cart(x, cartpos)

! Transform normal mode positions into cartesian positions.
!------------------------------------------------------------------------------

    use pimd_vars, only : nbead

    use cmd_vars, only : nmode_to_cart
    !..................................................

    implicit none

#include "memory.h"
#include "les.h"

    integer :: first, last

    _REAL_ :: x(3, natom), cartpos(3, natom)
    !..................................................

    integer :: iatom, ibead

    !..................................................

#ifdef LES
    call trans_part(x, cartpos, nmode_to_cart)
#else

# ifdef MPI
    !RCW: Can cartpos be used here as a temporary 3 x natom scratch array?
    call xdist(x, cartpos, natom)
# endif
    call trans_full(x, cartpos, nmode_to_cart)
#endif

end subroutine trans_pos_nmode_to_cart

!------------------------------------------------------------------------------
subroutine trans_frc_cart_to_nmode(f)

! Transform cartesian forces into normal mode forces
!------------------------------------------------------------------------------

    use pimd_vars, only : nbead, tempbuf

    use cmd_vars, only : fcart_to_fnmode

    !..................................................

    implicit none

#include "memory.h"
#include "les.h"

    integer :: first, last

    _REAL_  :: f(3, natom)

    !..................................................

    integer :: iatom, ibead

    !..................................................

#ifdef LES
    call trans_part(f, tempbuf, fcart_to_fnmode)
#else

# ifdef MPI
    !RCW: Can tempbuf be used here as a temporary 3 x natom scratch array?
    call xdist(f, tempbuf, natom)
# endif

    call trans_full(f, tempbuf, fcart_to_fnmode)

#endif
    f(1:3, 1:natom) = tempbuf(1:3, 1:natom)
    return

end subroutine trans_frc_cart_to_nmode

!------------------------------------------------------------------------------
subroutine trans_vel_nmode_to_cart(vel_nmode, cartvel)

! Transform normal mode positions into cartesian positions.
!------------------------------------------------------------------------------

    use pimd_vars, only : nbead

    use cmd_vars, only : nmode_to_cart

    !..................................................

    implicit none
#include "les.h"
#include "memory.h"

    integer :: iatom, first, last

    _REAL_, intent(in) :: vel_nmode(3, natom), cartvel(3, natom)

    !..................................................
    integer :: iatomCL, ibead, idim
    !..................................................

#ifdef LES
    call trans_part(vel_nmode, cartvel, nmode_to_cart)
#else

# ifdef MPI
    !RCW: Can cartvel be used here as a temporary 3 x natom scratch array?
    call xdist(vel_nmode, cartvel, natom)
# endif
    call trans_full(vel_nmode, cartvel, nmode_to_cart)
    return
#endif
end subroutine trans_vel_nmode_to_cart

!------------------------------------------------------------------------------
subroutine trans_vel_from_cart_to_nmode(v)

! Transform normal mode positions into cartesian positions.
!------------------------------------------------------------------------------

    use pimd_vars, only : nbead, tempbuf

    use cmd_vars, only : cart_to_nmode

    !..................................................

    implicit none
#include "memory.h"
#include "extra.h"
#include "les.h"

    integer :: first, last

    _REAL_  :: v(3, natom)

    !..................................................

    integer :: iatom, ibead, idim

    !..................................................

#ifdef LES
    call trans_part(v, tempbuf, cart_to_nmode)
#else
    call trans_full(v, tempbuf, cart_to_nmode)
#endif
    v(1:3, 1:natom) = tempbuf(1:3, 1:natom)

end subroutine trans_vel_from_cart_to_nmode

!------------------------------------------------------------------------------
subroutine full_scale_vel_centroid(v, mass, istart, iend)
!
! scale centroid velocities and fix total momentum equal to zero
!------------------------------------------------------------------------------

    use constants, only : kB
    use full_pimd_vars, only : mybeadid
    !..................................................

    implicit none

#include "md.h"
#include "les.h"
#include "memory.h"
#ifdef MPI
#  include "parallel.h"
    include 'mpif.h'
#endif
    _REAL_ :: v(3, natom), mass(natom)

    integer iatom, istart, iend, ierr, ibead

    !..................................................

    _REAL_ :: scale_factor

    _REAL_ :: E_kin2, temp_ini

    _REAL_ :: p(3)
#ifndef USE_MPI_IN_PLACE
    _REAL_ :: ptmp(3)
#endif

    !..................................................

    if (mybeadid .ne. 1) return

    !! Set total momentum equal to zero for path centroid.
    p(1:3) = 0.0
    do iatom = istart, iend
        p(1) = p(1) + mass(iatom)*v(1, iatom)
        p(2) = p(2) + mass(iatom)*v(2, iatom)
        p(3) = p(3) + mass(iatom)*v(3, iatom)
    end do

    p(1) = p(1)/dble(natom)
    p(2) = p(2)/dble(natom)
    p(3) = p(3)/dble(natom)

#ifdef MPI
# ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, p, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
# else
    call mpi_allreduce(p, ptmp, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
    p = ptmp
# endif
#endif

    E_kin2 = 0.d0

    do iatom = istart, iend
        v(1, iatom) = v(1, iatom) - p(1)/mass(iatom)
        v(2, iatom) = v(2, iatom) - p(2)/mass(iatom)
        v(3, iatom) = v(3, iatom) - p(3)/mass(iatom)

        E_kin2 = E_kin2 + mass(iatom) &
            *(v(1, iatom)**2 + v(2, iatom)**2 + v(3, iatom)**2)
    end do

#ifdef MPI
# ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, E_kin2, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
# else
    call mpi_allreduce(E_kin2, ptmp, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
    E_kin2 = ptmp(1)
# endif
#endif

    temp_ini = E_kin2/3.d0/kB/dble(natom)

    scale_factor = sqrt(temp0/temp_ini)

    !! Scale velocity according to target temperature.
    do iatom = istart, iend
        v(1, iatom) = v(1, iatom)*scale_factor
        v(2, iatom) = v(2, iatom)*scale_factor
        v(3, iatom) = v(3, iatom)*scale_factor
    end do

end subroutine full_scale_vel_centroid

#ifdef LES
!------------------------------------------------------------------------------
subroutine part_scale_vel_centroid(v, mass, istart, iend)
!
! scale centroid velocities and fix total momentum equal to zero
!------------------------------------------------------------------------------

    use constants, only : kB

    use pimd_vars, only : natomCL

    !..................................................

    implicit none

#include "md.h"
#include "les.h"
#include "memory.h"
#ifdef MPI
#  include "parallel.h"
    include 'mpif.h'
#endif
    _REAL_ :: v(3, natom), mass(natom)

    integer iatom, istart, iend, ierr, ibead

    !..................................................

    _REAL_ :: scale_factor

    _REAL_ :: E_kin2, temp_ini

    _REAL_ :: p(3)
#ifndef USE_MPI_IN_PLACE
    _REAL_ :: ptmp(3)
#endif

    !..................................................

    !! Set total momentum equal to zero for path centroid.
    p(1:3) = 0.0
    do iatom = istart, iend
        if (cnum(iatom) .eq. 0 .or. cnum(iatom) .eq. 1) then
            p(1) = p(1) + mass(iatom)*v(1, iatom)
            p(2) = p(2) + mass(iatom)*v(2, iatom)
            p(3) = p(3) + mass(iatom)*v(3, iatom)
        end if
    end do

    p(1) = p(1)/dble(natomCL)
    p(2) = p(2)/dble(natomCL)
    p(3) = p(3)/dble(natomCL)

#ifdef MPI
# ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, p, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
# else
    call mpi_allreduce(p, ptmp, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
    p = ptmp
# endif
#endif

    E_kin2 = 0.d0

    do iatom = istart, iend
        if (cnum(iatom) .eq. 0 .or. cnum(iatom) .eq. 1) then
            v(1, iatom) = v(1, iatom) - p(1)/mass(iatom)
            v(2, iatom) = v(2, iatom) - p(2)/mass(iatom)
            v(3, iatom) = v(3, iatom) - p(3)/mass(iatom)

            E_kin2 = E_kin2 + mass(iatom) &
                *(v(1, iatom)**2 + v(2, iatom)**2 + v(3, iatom)**2)
        end if
    end do

#ifdef MPI
# ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, E_kin2, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
# else
    call mpi_allreduce(E_kin2, ptmp, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
    E_kin2 = ptmp(1)
# endif
#endif

    temp_ini = E_kin2/3.d0/kB/dble(natomCL)

    scale_factor = sqrt(temp0/temp_ini)

    !! Scale velocity according to target temperature.
    do iatom = istart, iend
        if (cnum(iatom) .eq. 0 .or. cnum(iatom) .eq. 1) then
            v(1, iatom) = v(1, iatom)*scale_factor
            v(2, iatom) = v(2, iatom)*scale_factor
            v(3, iatom) = v(3, iatom)*scale_factor
        end if
    end do

end subroutine part_scale_vel_centroid

!------------------------------------------------------------------------------
subroutine trans_part(a, b, a_to_b)
!------------------------------------------------------------------------------
    use pimd_vars, only : nbead
    implicit none

#include "les.h"
#include "memory.h"

    _REAL_ a(3, *), b(3, *), a_to_b(nbead, nbead)
    integer iatom, ibead, first, last

    do iatom = 1, natom
        if (cnum(iatom) .eq. 0) then
            b(1:3, iatom) = a(1:3, iatom)
        else
            ibead = cnum(iatom)
            first = iatom - ibead + 1
            last = iatom - ibead + nbead
            b(1, iatom) = &
                sum(a_to_b(ibead, 1:nbead)*a(1, first:last))

            b(2, iatom) = &
                sum(a_to_b(ibead, 1:nbead)*a(2, first:last))

            b(3, iatom) = &
                sum(a_to_b(ibead, 1:nbead)*a(3, first:last))
        end if
    end do
end subroutine trans_part

#endif

!------------------------------------------------------------------------------
subroutine trans_full(a, b, a_to_b)

! Transform cartesian positions into normal mode positions.
    use pimd_vars, only : nbead
    use full_pimd_vars, only : mybeadid, xall
    !..................................................
    implicit none

#include "les.h"
#include "memory.h"

    integer :: iatm3, ibead, ierr
    integer :: jatm3, jbead
    integer :: myfirst, mylast

#ifdef MPI
    include 'mpif.h'
# include "parallel.h"
#endif

    _REAL_, intent(in) :: a(3*natom), a_to_b(nbead, nbead)
    _REAL_  :: b(3, natom)
    !..................................................
#ifdef MPI
    if (sanderrank .eq. 0) then
        call mpi_allgather(a, 3*natom, MPI_DOUBLE_PRECISION, &
            xall, 3*natom, MPI_DOUBLE_PRECISION, &
            commmaster, ierr)
    end if

    call mpi_bcast(xall, 3*natom*nbead, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif

    b(1:3, 1:natom) = 0.d0

    do jbead = 1, nbead
        b(1:3, 1:natom) = b(1:3, 1:natom) + a_to_b(mybeadid, jbead)*xall(1:3, 1:natom, jbead)
    end do

end subroutine trans_full
