#ifndef NCSU_UTILS_H
#define NCSU_UTILS_H

/*
 *  'afailed' is provided by 'ncsu_utils' module
 */

#ifndef NCSU_DISABLE_ASSERT
#  define ncsu_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define ncsu_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NCSU_PURE_EXCEPT_ASSERT
#  define NCSU_USE_AFAILED use ncsu_utils, only : afailed
#else
#  define ncsu_assert(s) 
#  define ncsu_assert_not_reached() 
#  define NCSU_PURE_EXCEPT_ASSERT pure
#  define NCSU_USE_AFAILED
#endif /* NCSU_DISABLE_ASSERT */

#define NCSU_OUT_OF_MEMORY call out_of_memory(__FILE__, __LINE__)

#ifdef MPI
#  define NCSU_MASTER_ONLY_BEGIN if (sanderrank.eq.0) then
#  define NCSU_MASTER_ONLY_END end if
#else
#  define NCSU_MASTER_ONLY_BEGIN
#  define NCSU_MASTER_ONLY_END
#endif /* MPI */

#define NCSU_ERROR   ' ** NCSU-Error ** : '
#define NCSU_WARNING ' ** NCSU-Warning ** : '
#define NCSU_INFO    ' NCSU : '

#endif /* NCSU_UTILS_H */
