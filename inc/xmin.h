!+ Specification of Amber's Interface to LMOD

!  Use defines instead of parameters because this file may
!  be included by both C and Fortran.

! Xmin return flags
#define DONE            0
#define CALCENRG        1
#define CALCGRAD        2
#define CALCBOTH        3
#define CALCENRG_NEWNBL 4
#define CALCGRAD_NEWNBL 5
#define CALCBOTH_NEWNBL 6
#define CALCENRG_OLDNBL 7
#define CALCGRAD_OLDNBL 8
#define CALCBOTH_OLDNBL 9

! Lmod return flags
#define UNIMPLEMENTED   -1
#define MINIMIZE        5
#define RELAX           6

! NonBond List Update control options
#define NBL_UPDATE_ALWAYS 569
#define NBL_UPDATE_NEVER  691
#define NBL_UPDATE_NORMAL 769

