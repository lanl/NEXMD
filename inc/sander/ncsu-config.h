#ifndef NCSU_CONFIG_H
#define NCSU_CONFIG_H

#if defined(LES) || defined(QMMM)
#  define DISABLE_NCSU yes
#endif

#if defined(MPI) && defined(BINTRAJ) && !defined(DISABLE_NCSU)
#  define NCSU_ENABLE_BBMD sure,whynot
#endif

#ifndef _REAL_
#  include "dprec.fh"
#  define NCSU_REAL _REAL_
#  ifdef DPREC
#    define NCSU_REAL_IS_DOUBLE indeed
#  endif
#endif /* _REAL_ */

#ifdef NCSU_REAL_IS_DOUBLE
#  define NCSU_TO_REAL(x) dble(x)
#else
#  define NCSU_TO_REAL(x) real(x)
#endif /* NCSU_REAL_IS_DOUBLE */

#define SANDER_STDERR_UNIT 6
#define SANDER_STDOUT_UNIT 6

/* EVB uses 75th (see files.h) */
#define SANDER_LAST_UNIT 77

#ifndef BINTRAJ
#  define NCSU_NO_NETCDF entirely
#endif /* BINTRAJ */

#endif /* NCSU_CONFIG_H */
