#include "ncsu-config.h"

module ncsu_constants

implicit none

private

#ifdef NCSU_REAL_IS_DOUBLE
NCSU_REAL, public, parameter :: ZERO  = 0.d0
NCSU_REAL, public, parameter :: ONE   = 1.d0
NCSU_REAL, public, parameter :: TWO   = 2.d0
NCSU_REAL, public, parameter :: THREE = 3.d0
NCSU_REAL, public, parameter :: FOUR  = 4.d0
#else
NCSU_REAL, public, parameter :: ZERO  = 0.0
NCSU_REAL, public, parameter :: ONE   = 1.0
NCSU_REAL, public, parameter :: TWO   = 2.0
NCSU_REAL, public, parameter :: THREE = 3.0
NCSU_REAL, public, parameter :: FOUR  = 4.0
#endif /* NCSU_REAL_IS_DOUBLE */

integer, public, parameter :: STRING_LENGTH = 256

integer, public, parameter :: ERR_UNIT = SANDER_STDERR_UNIT
integer, public, parameter :: OUT_UNIT = SANDER_STDOUT_UNIT

integer, public, parameter :: ABMD_MONITOR_UNIT = SANDER_LAST_UNIT + 1
integer, public, parameter :: SMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 2
integer, public, parameter :: PMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 3

#ifdef MPI
integer, public, parameter :: REM_MDIN_UNIT = SANDER_LAST_UNIT + 4
integer, public, parameter :: PMD_REMLOG_UNIT = SANDER_LAST_UNIT + 5
integer, public, parameter :: ABMD_REMLOG_UNIT = SANDER_LAST_UNIT + 6
integer, public, parameter :: BBMD_MONITOR_UNIT = SANDER_LAST_UNIT + 7
integer, public, parameter :: BBMD_LOG_UNIT = SANDER_LAST_UNIT + 8
#endif /* MPI */

integer, public, parameter :: EVEC_UNIT1 = SANDER_LAST_UNIT + 9
integer, public, parameter :: CRD_UNIT1  = SANDER_LAST_UNIT + 10
integer, public, parameter :: REF_UNIT1 = SANDER_LAST_UNIT + 11
integer, public, parameter :: IDX_UNIT1 = SANDER_LAST_UNIT + 12

end module ncsu_constants
