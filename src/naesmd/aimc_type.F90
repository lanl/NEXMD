#include "dprec.fh"
#include "assert.fh"

module AIMC_type_module
implicit none

private

public :: AIMC_type

type AIMC_type
        _REAL_, dimension(:), allocatable  :: FM ! AIMC FM (PCCP18)
        _REAL_, dimension(:), allocatable  :: FE ! AIMC "second" force term (PCCP18)
        _REAL_, dimension(:), allocatable  :: Fmax ! AIMC force on max populated surface (PCCP18)
        _REAL_ :: Weight=1
        _REAL_ :: delta_clone_1
        _REAL_ :: delta_clone_2
        _REAL_ :: delta_clone_3
        _REAL_ :: force_pop_min
        integer :: nclones,nclones_max,imax
        logical :: new_clone=.false.
end type AIMC_type

end module
