#include "dprec.fh"

!------------------------------------------------------------------------------
module nose_hoover_vars

! internal variables used for Nose'-Hoover thermostat
!------------------------------------------------------------------------------

    use nose_hoover_module, only : Thermostat_type

    implicit none

    save

    logical :: use_nose_hoover = .false.

    integer :: file_nhc = 1001

    integer :: nchain, nthermo

    type(Thermostat_type), allocatable :: thermo(:, :)

! type( Thermostat_type ), allocatable :: thermo_nmode( :, :, : )

    _REAL_ :: Econserved = 0.d0

    _REAL_ :: tau  !! characteristic time scale of the system

end module nose_hoover_vars
