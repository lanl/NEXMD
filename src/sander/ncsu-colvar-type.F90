! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar_type

implicit none

private

integer, public, parameter :: COLVAR_ANGLE           = 1
integer, public, parameter :: COLVAR_TORSION         = 2
integer, public, parameter :: COLVAR_DISTANCE        = 3
integer, public, parameter :: COLVAR_MULTI_RMSD      = 4
integer, public, parameter :: COLVAR_R_OF_GYRATION   = 5
integer, public, parameter :: COLVAR_HANDEDNESS      = 6
integer, public, parameter :: COLVAR_N_OF_BONDS      = 7
integer, public, parameter :: COLVAR_N_OF_STRUCTURES = 8
integer, public, parameter :: COLVAR_LCOD            = 9
integer, public, parameter :: COLVAR_COS_OF_DIHEDRAL = 10
integer, public, parameter :: COLVAR_COM_ANGLE       = 11
integer, public, parameter :: COLVAR_COM_TORSION     = 12
integer, public, parameter :: COLVAR_COM_DISTANCE    = 13
integer, public, parameter :: COLVAR_PCA             = 14



type, public :: colvar_t

   integer :: type = -1

   integer,   pointer :: i(:) => null()
   NCSU_REAL, pointer :: r(:) => null()
   
   integer :: tag ! (see ncsu-cv-priv.*)

   ! avgcrd : average crd of the trajectory 
   ! r      : reference crd 
   ! evec   : eigenvector from PCA 
   NCSU_REAL, pointer :: avgcrd(:) => null() 
   NCSU_REAL, pointer :: evec(:) => null() 
     
   ! state(:) stores the sate of reference part of ref.crd
   
   integer,  pointer :: state_ref(:) => null()
   integer,  pointer :: state_pca(:) => null()
   integer,  pointer :: ipca_to_i(:) => null()
 
end type colvar_t

end module ncsu_colvar_type
