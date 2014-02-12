!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module file_io_dat contains all of the information regarding file IO objects
module file_io_dat

! For now it'll just be files.h. This allows us to include only those properties
! that we want in various subroutines. My main reason for doing this for now is
! so that we can pass the MAX_FN_LEN parameter all over the place without having
! to include all of files.h

#include "files.h"

end module file_io_dat
