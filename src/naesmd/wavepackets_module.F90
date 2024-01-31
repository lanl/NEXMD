#include "dprec.fh"
#include "assert.fh"

module wavepackets_module
    use naesmd_constants
    use communism
    implicit none

    private

    public :: MCE_structure
    complex*16 :: WP_coeff

    type wavepacket_structure
        _REAL_, allocatable :: width_fixed_diagonal(:)
        _REAL_, allocatable :: WP_Positions(:)
        _REAL_, allocatable :: WP_Momemntum(:)
    end type wavepacket_structure

    type MCE_structure
    end type MCE_structure

end module wavepackets_module
