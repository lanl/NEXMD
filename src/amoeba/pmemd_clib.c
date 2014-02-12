#include <stdlib.h>
#include <sys/time.h>
#ifndef NO_RLIMIT_STACK_CTRL
#include <sys/resource.h>
#endif

/*
 * CLINK_CAPS is used for unicos (cray)
 * NO_C_UNDERSCORE is used for ibm spx
 * DBL_C_UNDERSCORE is used for pathscale and g95
 */

#ifdef CLINK_CAPS
#define PROC_GET_BYTESIZE       GET_BYTESIZE
#define PROC_GET_WALL_TIME      GET_WALL_TIME
#define PROC_UNLIMIT_STACK      UNLIMIT_STACK
#else
#ifdef NO_C_UNDERSCORE
#define PROC_GET_BYTESIZE       get_bytesize
#define PROC_GET_WALL_TIME      get_wall_time
#define PROC_UNLIMIT_STACK      unlimit_stack
#else
#ifdef DBL_C_UNDERSCORE
#define PROC_GET_BYTESIZE       get_bytesize__
#define PROC_GET_WALL_TIME      get_wall_time__
#define PROC_UNLIMIT_STACK      unlimit_stack__
#else
#define PROC_GET_BYTESIZE       get_bytesize_
#define PROC_GET_WALL_TIME      get_wall_time_
#define PROC_UNLIMIT_STACK      unlimit_stack_
#endif
#endif
#endif

void
PROC_GET_BYTESIZE(void * start, void * end, int * bytes)
{
  *bytes = ((char *)end - (char *)start);
}

void
PROC_GET_WALL_TIME(int * sec, int * usec)
{
  struct timeval        wall_time;

  gettimeofday(&wall_time, NULL);

  *sec = wall_time.tv_sec;
  *usec = wall_time.tv_usec;

  return;
}

/*
 * This routine returns -1 to indicated that the stack has been reset to the
 * maximum.  For machines that define NO_RLIMIT_STACK_CTRL because they don't
 * have the rlimit interface, we are basically lying...  If the routine returns
 * a positive value, that is the actual limit up to 1,000,000,000.  Any limit
 * over 1,000,000,000 is considered "unlimited" for purposes of this code (this
 * actually helps keep 32/64 bit issues from unnecessarily setting off alarms
 * on systems like AIX.
 */

void
PROC_UNLIMIT_STACK(int * new_limit)
{
#ifdef NO_RLIMIT_STACK_CTRL
  *new_limit = -1;
#else
  struct rlimit         rlim;

  getrlimit(RLIMIT_STACK, &rlim);

  if (rlim.rlim_cur == RLIM_INFINITY)
  {
    /* the stack is already unlimited */
    *new_limit = -1;
  }
  else
  {
    if (rlim.rlim_max == RLIM_INFINITY)
    {
      /* the stack can be set to unlimited by anyone */
      rlim.rlim_cur = rlim.rlim_max;
      setrlimit(RLIMIT_STACK, &rlim);
      *new_limit = -1;
    }
    else
    {
      /* there is a hard limit in effect; bump up to it */

      if (rlim.rlim_max >= 1000000000) /* we consider 1 billion "unlimited" */
        *new_limit = -1;
      else
        *new_limit = (int)rlim.rlim_max;

      if (rlim.rlim_cur < rlim.rlim_max)
      {
        rlim.rlim_cur = rlim.rlim_max;
        setrlimit(RLIMIT_STACK, &rlim);
      }
    }
  }
#endif

  return;
}
