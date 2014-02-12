#ifdef CLINK_CAPS
#  define wallclock_ WALLCLOCK
#  define nwallclock_ NWALLCLOCK
#  define microsec_ MICROSEC
#endif
#ifdef CLINK_PLAIN
#  define wallclock_ wallclock
#  define nwallclock_ nwallclock
#  define microsec_ microsec
#endif

#ifndef F90_TIMER

#include <sys/time.h>
#include <unistd.h>

static int ncalls=0;

void wallclock_( double *wallc ){

#  ifdef NO_DETAILED_TIMINGS
	*wallc = 0.0;
#  else
	struct timeval tv;
	struct timezone tz;
	long rsec;

	gettimeofday( &tv, &tz );

	rsec = tv.tv_sec - 959800000;
	*wallc = (double)rsec + (double)(tv.tv_usec)/1000000.;
	ncalls++;
#  endif

	return;

}

void microsec_( int *ig ){

	struct timeval tv;
	struct timezone tz;

	gettimeofday( &tv, &tz );
	*ig = (int)tv.tv_usec;

	return;

}

void nwallclock_( int *n ){  
		/* provides the number of times wallclock was called */

	*n = ncalls;
	return;

}
#endif
