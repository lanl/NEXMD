#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#ifdef WINDOWS
/* tw struct timeval is in there on windows */
#include <Winsock2.h>
#endif

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

#ifdef WINDOWS
/* tw 
 * replacement for gettimeofday based on  
 * http://www.cpp-programming.net/c-tidbits/gettimeofday-function-for-windows/
 */
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
 
struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	unsigned __int64 tmpres = 0;
	static int tzflag;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		/*converting file time to unix epoch*/
		tmpres /= 10;  /*convert into microseconds*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS;
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}
/* tw makes deprecated error but unused anyway   
	if (NULL != tz)
	{
		if (!tzflag)
		{
			_tzset();
			tzflag++;
		}
		tz->tz_minuteswest = _get_timezone / 60;
		tz->tz_dsttime = _get_daylight;
	}
*/
	return 0;
}
/* end time.h from tw */
#endif 

#ifdef WINDOWS
extern "C" void
#else
void
#endif
PROC_GET_BYTESIZE(void * start, void * end, int * bytes)
{
  *bytes = ((char *)end - (char *)start);
}

#ifdef WINDOWS
extern "C" void
#else
void
#endif
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


#ifdef WINDOWS
extern "C" void
#else
void
#endif
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
