#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#ifndef CYGWIN
#  ifdef USE_AMBER_C9XCOMPLEX
#    include "complex.h"
#  else
#    include <complex.h>
#  endif
#endif
#include "sff.h"
#include "memutil.h"
#include "timer.h"

#ifdef SCALAPACK
#include "mpi.h"
#endif

/* Here are offsets into the ScaLAPACK descriptor array. */

#define DTYPE_ (0)
#define CTXT_ (1)
#define M_ (2)
#define N_ (3)
#define MB_ (4)
#define NB_ (5)
#define RSRC_ (6)
#define CSRC_ (7)
#define LLD_ (8)
#define DLEN_ (9)

REAL_T seconds(void);

INT_T get_blocksize(void);

INT_T get_mytaskid(void);

INT_T get_numtasks(void);

INT_T dspev_(char *, char *, int *, REAL_T *,
             REAL_T *, REAL_T *, int *, REAL_T *, int *);

INT_T dsyev_(char *, char *, INT_T *, REAL_T *, INT_T *,
             REAL_T *, REAL_T *, INT_T *, INT_T *);

INT_T dsyevd_(char *jobz, char *uplo, INT_T *, REAL_T *,
              INT_T *, REAL_T *, REAL_T *, INT_T *,
              INT_T *, INT_T *, INT_T *);

INT_T dpotrf_(char *, INT_T *, REAL_T *, INT_T *, INT_T *);

INT_T dpotrs_(char *, INT_T *, INT_T *, REAL_T *, INT_T *, REAL_T *,
              INT_T *, INT_T *);

void dcopy_(INT_T *, REAL_T *, INT_T *, REAL_T *, INT_T *);

INT_T dsaupd_(INT_T *, char *, INT_T *, char *, INT_T *,
              REAL_T *, REAL_T *, INT_T *, REAL_T *, INT_T *, INT_T *,
              INT_T *, REAL_T *, REAL_T *, INT_T *, INT_T *);

INT_T dseupd_(INT_T *, char *, INT_T *,
              REAL_T *, REAL_T *, INT_T *, REAL_T *,
              char *, INT_T *, char *, INT_T *, REAL_T *,
              REAL_T *, INT_T *, REAL_T *, INT_T *, INT_T *, INT_T *,
              REAL_T *, REAL_T *, INT_T *, INT_T *, short *, short *,
              short *);

INT_T dgetrs_(char *, INT_T *, INT_T *, REAL_T *, INT_T *, INT_T *,
              REAL_T *, INT_T *, INT_T *);

void blacs_pinfo_(INT_T *, INT_T *);

void sl_init_(INT_T *, INT_T *, INT_T *);

INT_T sl_gridreshape_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *,
                      INT_T *);

void blacs_gridinfo_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void blacs_barrier_(INT_T *, char *);

void blacs_gridexit_(INT_T *);

INT_T numroc_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void descinit_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *,
               INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void pdsyev_(char *, char *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *);

void pdgemr2d_(INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *, INT_T *);

INT_T myroc(int, int, int, int);

REAL_T *ptr2d(REAL_T *, int *, int, int);


FILE *safe_fopen(char *, char *);

static int c__3 = 3;

INT_T dgeev_(char *, char *, INT_T *, REAL_T *, INT_T *,
             REAL_T *, REAL_T *, REAL_T *, INT_T *,
             REAL_T *, INT_T *, REAL_T *, INT_T *, INT_T *);

REAL_T atof(const char *);

INT_T dgemv_(char *, INT_T *, INT_T *, REAL_T *, REAL_T *, INT_T *,
             REAL_T *, INT_T *, REAL_T *, REAL_T *, INT_T *);

static INT_T typerun;

FILE	*safe_fopen( char *filename, char *mode )
{
	FILE *f;
	f = fopen( filename, mode );
	if( f == NULL ){
		fprintf( stderr, "unable to open file %s\n", filename );
		exit(1);
	}
	return( f );
}

#ifdef CYGWIN
/*********************************************************************** 
 * Internal complex number math routines, for use with CYGWIN.
 * Daniel R. Roe - 2010-03-02
 **********************************************************************/
#define ONEOVERROOTTWO 0.707106781
typedef struct _complextype {
  double r;
  double i;
} complextype;
/*
 * complex_sqrt()
 * Calculate the square root of a complex number.
 *
 * Given a complex number c, there exists a complex number z such that:
 *   c = z^2
 * Every complex number except 0 will have two roots, positive and negative.
 * If 
 *   c = a + bI (a and b are both real)
 * and 
 *   z = p + qI
 * then
 *   a + bI = (p + qI)^2
 * The right hand side becomes: 
 *   p^2 + pqI + pqI +qI^2 --> p^2 - q^2 + 2pqI
 * Equating the real and imaginary parts to a + bI yields:
 *   a = p^2 - q^2 [eq 1]
 *   b = 2pq       [eq 2]
 * Solving the second eq. for q yields:
 *   q = b / 2p    [eq 3]
 * Which when substituted back into the first eq yields:
 *   a = p^2 - (b / 2p)^2, a = p^2 - b^2 / 4p^2
 * Multiple both sides by 4p^2:
 *   4ap^2 = 4p^4 - b^2
 * Which can be converted to quadratic form for p^2:
 *   4p^4 - 4ap^2 - b^2 = 0
 * p^2 can be solved by using the quadratic formula:
 *   p^2 = (a + (a^2 + b^2)^0.5) / 2
 *   p = 2^-0.5 (a + (a^2 + b^2)^0.5)^0.5 [eq 4]
 * The fourth equation can then be substituted back into the third:
 *   q = b / 2 (2^-0.5 (a + (a^2 + b^2)^0.5)^0.5)
 * which can then be rearranged to:
 *   q = (b / abs(b)) (2^-0.5) (( (a^2 + b^2)^0.5 ) - a)^0.5
 */
complextype complex_sqrt(complextype C) {
  double a2, b2, sq_a2b2, a_sq_a2b2, p, q, bsign;
  complextype answer;

  //fprintf(stderr,"Computing complex sqrt for %lf + %lf I)\n",C.r,C.i);

  a2=C.r * C.r;
  b2=C.i * C.i;
  if (b2>0)
    bsign=1.0;
  else
    bsign=-1.0;

  sq_a2b2 = sqrt(a2 + b2);

  a_sq_a2b2 = C.r + sq_a2b2;
     
  p = ONEOVERROOTTWO * sqrt(a_sq_a2b2);

  a_sq_a2b2 = sq_a2b2 - C.r;

  q = ONEOVERROOTTWO * bsign * sqrt(a_sq_a2b2);

  //fprintf(stderr,"Real: %lf   Imag: %lf\n",p,q);

  answer.r=p;
  answer.i=q;

  return answer;
}

/*
 * complex_div()
 * Perform division of two complex numbers
 */
complextype complex_div(complextype N, complextype D) {
  double a,b,c,d;
  double c2d2;
  complextype answer;

  //fprintf(stderr, "Computing complex div for %lf + %lf I / %lf + %lf I\n",
  //        N.r,N.i,D.r,D.i);

  a=N.r;
  b=N.i;
  c=D.r;
  d=D.i;
  answer.r=0.0;
  answer.i=0.0;

  if (c==0.0 && d==0.0) {
    fprintf(stderr,"Error: nmode.c: complex_div: Denominator is 0.0, answer is undefined.\n");
    return answer;
  }

  c2d2 = (c * c) + (d * d);
  
  answer.r = ((a * c) + (b * d)) / c2d2;

  answer.i = ((b * c) - (a * d)) / c2d2;

  //fprintf(stderr,"Real: %lf   Imag: %lf\n",answer.r, answer.i);

  return answer;
}

/*
 * complex_add()
 * Perform addition of two complex numbers
 */
complextype complex_add(complextype C1, complextype C2) {
  complextype answer;

  answer.r = C1.r + C2.r;
  answer.i = C1.i + C2.i;

  return answer;
}

/*
 * complex_mult()
 * Perform multiplication of two complex numbers
 */
complextype complex_mult(complextype C1, complextype C2) {
  complextype answer;

  answer.r = (C1.r * C2.r) - (C1.i * C2.i);
  answer.i = (C1.i * C2.r) + (C1.r * C2.i);

  return answer;
}
#endif //CYGWIN

/***********************************************************************
                            SYMNUM()
************************************************************************/

/*     ----- routine to give the symmetry number. only for linear */
/*           molecules. for others symmetry number is unity ----- */

static int symnum(int natoms, REAL_T * amass, REAL_T * sn, int *linear)
{


   /* Parameter adjustments */
   --amass;

   /* Function Body */
   *sn = 1.;
   if (natoms <= 2) {
      *linear = 1;
      if (amass[1] == amass[2]) {
         *sn = 2.;
      }
   }
   return 0;
}                               /* symnum */

/* ----------------------------------------------------------------------- */

static int mofi(int natoms, REAL_T * c, REAL_T * amass, REAL_T * pmom)
{
   /* Initialized data */

   static char uplo[1] = "U";
   static char jobz[1] = "V";
   static REAL_T zero = 0.;

   /* System generated locals */

   /* Local variables */
   static REAL_T t[9], x, y, z, e2[30], wt;
   static int i, iat;
   static REAL_T com[3];
   static int ier, iaind;
   static REAL_T totwt, eigvec[9];


/*     compute the principal moments of inertia. */
/*     units are amu-bohr**2 */

   /* Parameter adjustments */
   --pmom;
   --amass;
   --c;

   /* Function Body */


/*     compute the position of the center of mass and translate */
/*     it to the origin. */

   com[0] = zero;
   com[1] = zero;
   com[2] = zero;

   totwt = zero;
   for (iat = 1; iat <= natoms; ++iat) {
      iaind = (iat - 1) * 3;
      wt = amass[iat];
      totwt += wt;
      com[0] += wt * c[iaind + 1];
      com[1] += wt * c[iaind + 2];
      com[2] += wt * c[iaind + 3];
   }

   com[0] /= totwt;
   com[1] /= totwt;
   com[2] /= totwt;

/*     compute the principal moments. */

   for (i = 1; i <= 9; ++i) {
      t[i - 1] = zero;
   }

   for (iat = 1; iat <= natoms; ++iat) {
      wt = amass[iat];
      x = c[1 + 3 * (iat - 1)] - com[1 - 1];
      y = c[2 + 3 * (iat - 1)] - com[2 - 1];
      z = c[3 + 3 * (iat - 1)] - com[3 - 1];
      t[0] += wt * (y * y + z * z);
      t[2] += wt * (x * x + z * z);
      t[5] += wt * (x * x + y * y);
      t[1] -= wt * x * y;
      t[3] -= wt * x * z;
      t[4] -= wt * y * z;
   }
   ier = 0;
   dspev_(jobz, uplo, &c__3, t, &pmom[1], eigvec, &c__3, e2, &ier);
   return 0;
}                               /* mofi */

/***********************************************************************
                            THERMO()
************************************************************************/

static int thermo(int natoms, int n, int ilevel, REAL_T * c,
                  REAL_T * amass, REAL_T * freq, REAL_T t, REAL_T patm,
                  REAL_T energy)
{
   REAL_T *vtemp = NULL, *evibn = NULL, *cvibn = NULL, *svibn = NULL;

   /* Initialized data */

   static REAL_T zero = 0.;
   static REAL_T akilo = 1e3;
   static REAL_T pstd = 101325.;
   static REAL_T pt2 = .2;
   static REAL_T half = .5;
   static REAL_T one = 1.;
   static REAL_T onept5 = 1.5;
   static REAL_T two = 2.;
   static REAL_T twopt5 = 2.5;
   static REAL_T four = 4.;
   static REAL_T eight = 8.;

   /* System generated locals */
   REAL_T d__1;

   /* Local variables */
   REAL_T e;
   int i;
   REAL_T p, s, pi, cv, sn, rt, em1;
   int iff;
   REAL_T arg, gas;
   int iat;
   REAL_T con, dum, dum1, dum2, argd, cvib, evib;
   int ndof;
   REAL_T avog, ezkc, pipi, ezpe, tokg, svib, crot, pmom[10],
       erot, etot, ctot, srot, stot, tovt, jpcal, tocal, ccont, ctran,
       econt, etran, scont, stran, tomet, rtemp, boltz, etovt, rtemp1,
       rtemp2, rtemp3, planck;
   int linear;
   REAL_T tokcal, hartre, weight;
   int nimag;
   int mytaskid;

   /* Allocate temporary vectors. */

   vtemp = vector(1, n);
   evibn = vector(1, n);
   cvibn = vector(1, n);
   svibn = vector(1, n);

   /* Get the MPI task id. */

   mytaskid = get_mytaskid();

/*     given the structure of a molecule and its normal mode vibrational */
/*     frequencies this routine uses standard statistical mechanical */
/*     formulas for an ideal gas (in the canonical ensemble, see, */
/*     for example, d. a. mcquarrie, "statistical thermodynamics", */
/*     harper & row, new york, 1973, chapters 5, 6, and 8) to compute */
/*     the entropy, heat capacity, and internal energy. */

/*     the si system of units is used internally.  conversion to units */
/*     more familiar to most chemists is made for output. */


/*     amass:   atomic weights, in amu. */
/*     pmom:    principal moments of inertia, in amu-bohr**2 and */
/*              in ascending order. */
/*     freq:    vibrational frequencies, in cm**-1 and in ascending */
/*              order */
/*     c    :   coordinates in Angstroms */
/*     vtemp:   vibrational temperatures, in kelvin. */
/*     evibn:   contribution to e from the vibration n. */
/*     cvibn:   contribution to cv from the vibration n. */
/*     svibn:   contribution to s from the vibration n. */
/*     t:       temperature */
/*     patm:    pressure, in atmospheres */


   /* Parameter adjustments */
   --freq;
   --amass;
   --c;


/*     tokg:    kilograms per amu. */
/*     boltz:   boltzman constant, in joules per kelvin. */
/*     planck:  planck constant, in joule-seconds. */
/*     avog:    avogadro constant, in mol**(-1). */
/*     jpcal:   joules per calorie. */
/*     tomet:   metres per Angstrom. */
/*     hartre:  joules per hartree. */

   tokg = 1.660531e-27;
   boltz = 1.380622e-23;
   planck = 6.626196e-34;
   avog = 6.022169e23;
   jpcal = 4.18674;
   tomet = 1e-10;
   hartre = 4.35981e-18;

/*     compute the gas constant, pi, pi**2, and e. */
/*     compute the conversion factors cal per joule and kcal per joule. */

   gas = avog * boltz;
   pi = four * atan(one);
   pipi = pi * pi;
   e = exp(one);
   tocal = one / jpcal;
   tokcal = tocal / akilo;

/*     print the temperature and pressure. */

   p = pstd * patm;
   if (mytaskid == 0) {
      fprintf(nabout, "\n                - Thermochemistry -\n\n");
      fprintf(nabout, "Temperature:  %8.3f\n   Pressure:  %8.3f\n", t,
              patm);
   }
   rt = gas * t;

/*     compute and print the molecular mass in amu, then convert to */
/*     kilograms. */

   weight = zero;
   for (iat = 1; iat <= natoms; ++iat) {
      weight += amass[iat];
   }
   if (mytaskid == 0) {
      fprintf(nabout, "       Mass:  %8.3f\n", weight);
   }
   weight *= tokg;

/*     compute contributions due to translation: */
/*        etran-- internal energy */
/*        ctran-- constant v heat capacity */
/*        stran-- entropy */

   dum1 = boltz * t;
   d__1 = two * pi;
   dum2 = pow(d__1, onept5);
   arg = pow(dum1, onept5) / planck;
   arg = arg / p * (dum1 / planck);
   arg = arg * dum2 * (weight / planck);
   arg = arg * sqrt(weight) * pow(e, twopt5);
   stran = gas * log(arg);
   etran = onept5 * rt;
   ctran = onept5 * gas;

/*     Compute contributions due to electronic motion: */
/*        It is assumed that the first electronic excitation energy */
/*        is much greater than kt and that the ground state has a */
/*        degeneracy of one.  Under these conditions the electronic */
/*        partition function can be considered to be unity.  The */
/*        ground electronic state is taken to be the zero of */
/*        electronic energy. */

/*     for monatomics print and return. */

   if (natoms <= 1) {
      s = stran * tocal;
      e = etran * tokcal;
      cv = ctran * tocal;
      /* need print statement here */
      free_vector(vtemp, 1, n);
      free_vector(evibn, 1, n);
      free_vector(cvibn, 1, n);
      free_vector(svibn, 1, n);
      return 0;
   }

/*     compute contributions due to rotation. */

/*     Compute the principal moments of inertia, get the rotational */
/*     symmetry number, see if the molecule is linear, and compute */
/*     the rotational temperatures.  Note the imbedded conversion */
/*     of the moments to SI units. */

   mofi(natoms, &c[1], &amass[1], pmom);
   if (mytaskid == 0) {
      fprintf(nabout, "Principal moments of inertia in amu-A**2:\n");
      fprintf(nabout, "     %12.2f%12.2f%12.2f\n", pmom[0], pmom[1],
              pmom[2]);
   }
   linear = 0;
   symnum(natoms, &amass[1], &sn, &linear);
   if (mytaskid == 0) {
      fprintf(nabout, "Rotational symmetry number is %2.0f\n", sn);
   }
   con = planck / (boltz * eight * pipi);
   con = con / tokg * (planck / (tomet * tomet));
   if (linear) {
      rtemp = con / pmom[2];
      if (mytaskid == 0) {
         if (rtemp < pt2)
            fprintf(nabout, "Assuming classical behavior for rotation\n");
         fprintf(nabout, "Rotational temperature: %8.3f\n", rtemp);
      }
   } else {
      rtemp1 = con / pmom[0];
      rtemp2 = con / pmom[1];
      rtemp3 = con / pmom[2];
      if (mytaskid == 0) {
         if (rtemp1 < pt2) {
            fprintf(nabout, "Assuming classical behavior for rotation\n");
         }
         fprintf(nabout, "Rotational temperatures: %8.3f  %8.3f  %8.3f\n",
                 rtemp1, rtemp2, rtemp3);
      }
   }

/*         erot-- rotational contribution to internal energy. */
/*         crot-- rotational contribution to cv. */
/*         srot-- rotational contribution to entropy. */

   if (linear) {
      erot = rt;
      crot = gas;
      arg = t / rtemp * (e / sn);
      srot = gas * log(arg);
   } else {
      erot = onept5 * rt;
      crot = onept5 * gas;
      arg = sqrt(pi * e * e * e) / sn;
      dum = t / rtemp1 * (t / rtemp2) * (t / rtemp3);
      arg *= sqrt(dum);
      srot = gas * log(arg);
   }

/*     compute contributions due to vibration. */

/*     compute vibrational temperatures and zero point vibrational */
/*     energy.  only real frequencies are included in the analysis. */

   nimag = 0;
   for (i = 1; i <= 3 * natoms; i++)
      if (freq[i] < -0.5)
         nimag++;
   if (mytaskid == 0) {
      if (nimag > 0)
         fprintf(nabout, "found %d imaginary frequencies\n", nimag);
   }

/*       (---iff is the first frequency to include in thermo:) */

   if (ilevel != 0) {
      iff = 0 + nimag;
   } else if (linear) {
      iff = 5 + nimag;
   } else {
      iff = 6 + nimag;
   }
   ndof = 3 * natoms - iff;
   con = planck / boltz;
   ezpe = zero;
   for (i = 1; i <= ndof; ++i) {
      vtemp[i] = freq[i + iff] * con * 3e10;
      ezpe += freq[i + iff] * 3e10;
   }
   ezpe = half * planck * ezpe;
   ezkc = ezpe * tokcal * avog;
   if (mytaskid == 0) {
      fprintf(nabout, "Zero-point vibrational energy: %10.3f\n", ezkc);
   }


/*     compute: */
/*        evib-- the vibrational component of the internal energy. */
/*        cvib-- the vibrational component of the heat capacity. */
/*        svib-- the vibrational component of the entropy. */

   evib = zero;
   cvib = zero;
   svib = zero;
   for (i = 1; i <= ndof; ++i) {

/*       compute some common factors. */

      tovt = vtemp[i] / t;
      etovt = exp(tovt);
      em1 = etovt - one;

/*       compute contributions due to the i'th vibration. */

      econt = tovt * (half + one / em1);
      d__1 = tovt / em1;
      ccont = etovt * (d__1 * d__1);
      argd = one - one / etovt;
      if (argd > 1e-7) {
         scont = tovt / em1 - log(argd);
      } else {
         scont = 0.f;
         if (mytaskid == 0) {
            fprintf(nabout,
                    "Warning: setting vibrational entropy to zero for mode %d with vtemp = %8.3f\n",
                    i, vtemp[i]);
         }
      }
      evibn[i] = econt * rt;
      cvibn[i] = ccont * gas;
      svibn[i] = scont * gas;

      evib += econt;
      cvib += ccont;
      svib += scont;
   }
   evib *= rt;
   cvib *= gas;
   svib *= gas;

/*     the units are now: */
/*         e-- joules/mol */
/*         c-- joules/mol-kelvin */
/*         s-- joules/mol-kelvin */

   etot = etran + erot + evib;
   ctot = ctran + crot + cvib;
   stot = stran + srot + svib;


/*     convert to the following and print */
/*         e-- kcal/mol */
/*         c-- cal/mol-kelvin */
/*         s-- cal/mol-kelvin */

   etran *= tokcal;
   ctran *= tocal;
   stran *= tocal;
   erot *= tokcal;
   crot *= tocal;
   srot *= tocal;
   evib *= tokcal;
   cvib *= tocal;
   svib *= tocal;
   etot = etran + erot + evib + energy;
   ctot = ctran + crot + cvib;
   stot = stran + srot + svib;
   for (i = 1; i <= ndof; ++i) {
      evibn[i] *= tokcal;
      cvibn[i] *= tocal;
      svibn[i] *= tocal;
   }

   if (mytaskid == 0) {
      fprintf(nabout,
              "\n             freq.         E               Cv              S\n");
      fprintf(nabout,
              "            cm**-1       kcal/mol       cal/mol-K     cal/mol-K\n");
      fprintf(nabout, "Total:               %11.3f    %11.3f    %11.3f\n",
              etot, ctot, stot);
      fprintf(nabout, "translational:       %11.3f    %11.3f    %11.3f\n",
              etran, ctran, stran);
      fprintf(nabout, "rotational:          %11.3f    %11.3f    %11.3f\n",
              erot, crot, srot);
      fprintf(nabout, "vibrational:         %11.3f    %11.3f    %11.3f\n",
              evib, cvib, svib);
      fprintf(nabout, "ff   energy:         %11.3f\n", energy);

      for (i = 1; i <= iff; ++i) {
         fprintf(nabout, " %5d   %10.3f\n", i, freq[i]);
      }
      for (i = 1; i <= ndof; ++i) {
         fprintf(nabout, " %5d   %10.3f  %11.3f    %11.3f    %11.3f\n",
                 i + iff, freq[i + iff], evibn[i], cvibn[i], svibn[i]);
      }
   }
   free_vector(vtemp, 1, n);
   free_vector(evibn, 1, n);
   free_vector(cvibn, 1, n);
   free_vector(svibn, 1, n);

   return 0;
}                               /* thermo */


/***********************************************************************
                            sizegam()
************************************************************************/


int sizegam(int i)
{
   int j, k;

   j = 1;
   for (k = 2; k <= i; k++)
      j = j + k;

   return j;
}

/***********************************************************************
                            INDEXKL()
************************************************************************/


int indexkl(int k, int l)
{
   int i, max, min;
   if (k < l) {
      max = l;
      min = k;
   } else {
      max = k;
      min = l;
   }
   i = (max + 1) * (max) / 2 + min;
   return i;
}


/***********************************************************************
                            movecm()
************************************************************************/

static void movecm(REAL_T * xcom, REAL_T * x, REAL_T * m, INT_T natom)
{

   int i, i3;
   REAL_T sumx, sumy, sumz, sum, zero, mi, xc, yc, zc;

   zero = 0.0;

   sumx = zero;
   sumy = zero;
   sumz = zero;
   sum = zero;

   i3 = 0;
   for (i = 0; i < natom; i++) {
      mi = m[i + 1];
      sumx = sumx + x[i3] * mi;
      sumy = sumy + x[i3 + 1] * mi;
      sumz = sumz + x[i3 + 2] * mi;
      sum = sum + mi;
      i3 = i3 + 3;
   }
   xc = sumx / sum;
   yc = sumy / sum;
   zc = sumz / sum;
   i3 = 0;

   for (i = 0; i < natom; i++) {
      xcom[i3] = x[i3] - xc;
      xcom[i3 + 1] = x[i3 + 1] - yc;
      xcom[i3 + 2] = x[i3 + 2] - zc;
      i3 = i3 + 3;
   }




}



/***********************************************************************
                            remtranrot()
************************************************************************/

static void remtranrot(INT_T natom, INT_T nreal, REAL_T * x, REAL_T * z,
                       REAL_T * wr, REAL_T * wi, INT_T * indexsort,
                       REAL_T * dnorm, REAL_T * m, INT_T * iclass)
{




/*
   Remove translational and rotational parts from real Langevin modes
*/

   REAL_T xsum, ysum, zsum, tmas, vnor1, vnor2, vnor3, vnor4, sum, diff,
       dt;
   REAL_T *d = NULL, *a = NULL, *dbl = NULL;
   INT_T i, j, k, ip, it, iat, nr6, nr3, ik;
   int mytaskid;
   mytaskid = get_mytaskid();



   nr3 = 3 * natom;
   nr6 = 6 * natom;

   if (mytaskid == 0) {
     printf("\n");
   }
   /*  normalize transaltion and rotation degrees of freedom */

   xsum = 0.0;
   ysum = 0.0;
   zsum = 0.0;
   tmas = 0.0;

   for (i = 0; i < natom; i++) {
      xsum = xsum + x[3 * (i)] * x[3 * (i)] * m[i + 1];
      ysum = ysum + x[3 * (i) + 1] * x[3 * (i) + 1] * m[i + 1];
      zsum = zsum + x[3 * (i) + 2] * x[3 * (i) + 2] * m[i + 1];
      tmas = tmas + m[i + 1];
   }

   vnor1 = sqrt(xsum + ysum);
   vnor2 = sqrt(ysum + zsum);
   vnor3 = sqrt(xsum + zsum);
   vnor4 = sqrt(tmas);

   if (vnor1 == 0.0) {
      vnor1 = 1.0;
   }
   if (vnor2 == 0.0) {
      vnor2 = 1.0;
   }
   if (vnor3 == 0.0) {
      vnor3 = 1.0;
   }

/*
      set up the trsl. and rot.  vectors
      see C. Eckart , Phys. Rev. 47, 552 (1935)
*/

   d = vector(0, nr3 * 6);
   a = vector(0, 6);


   for (i = 0; i < 6 * nr3; i++) {
      d[i] = 0.0;
   }



   iat = 0;
   for (i = 0; i < nr3; i += 3) {
      d[i] = sqrt(m[iat + 1]) / vnor4;
      d[i + 1 + 1 * nr3] = sqrt(m[iat + 1]) / vnor4;
      d[i + 2 + 2 * nr3] = sqrt(m[iat + 1]) / vnor4;
      d[i + 3 * nr3] = -sqrt(m[iat + 1]) * x[3 * (iat) + 1] / vnor1;
      d[i + 1 + 3 * nr3] = sqrt(m[iat + 1]) * x[3 * (iat)] / vnor1;
      d[i + 1 + 4 * nr3] = -sqrt(m[iat + 1]) * x[3 * (iat) + 2] / vnor2;
      d[i + 2 + 4 * nr3] = sqrt(m[iat + 1]) * x[3 * (iat) + 1] / vnor2;
      d[i + 5 * nr3] = sqrt(m[iat + 1]) * x[3 * (iat) + 2] / vnor3;
      d[i + 2 + 5 * nr3] = -sqrt(m[iat + 1]) * x[3 * (iat)] / vnor3;
      iat = iat + 1;
   }






   /*

      subtract translational and rotational contributions
      from each eigenvector
    */

   for (k = 6; k < nreal; k++) {

      ik = indexsort[k];
      dnorm[ik] = 0.0;

      for (j = 0; j < 6; j++) {
         a[j] = 0.0;
         for (i = 0; i < nr3; i++) {
            a[j] = a[j] + z[i + ik * nr6] * d[i + j * nr3];
         }
      }


      for (i = 0; i < nr3; i++) {

         sum = 0.0;
         for (j = 0; j < 6; j++) {
            sum = sum + a[j] * d[i + j * nr3];
         }

         diff = z[i + ik * nr6] - sum;
         dnorm[ik] = dnorm[ik] + diff * diff;
      }




      sum = dnorm[ik];
      for (j = 0; j < 6; j++) {
         sum = sum + a[j] * a[j];
      }

      if (sum != 0.0) {
         dnorm[ik] = (dnorm[ik] / sum);
      }


   }


   free_vector(a, 0, 6);
   free_vector(d, 0, nr3 * 6);


   /*
      see if anything remains after that subtraction
      of course a tiny bit may remain always
      so remove 12 modes with smallest remainders
    */
   dbl = vector(0, nr6);


   for (i = 6; i < nreal; i++) {
      dbl[i] = dnorm[indexsort[i]];
   }



   for (i = 6; i < nreal - 1; i++) {

      ip = i;

      for (j = i + 1; j < nreal; j++) {
         if (dbl[j] < dbl[ip]) {
            ip = j;
         }
      }


      dt = dbl[i];
      dbl[i] = dbl[ip];
      dbl[ip] = dt;
      it = indexsort[i];
      indexsort[i] = indexsort[ip];
      indexsort[ip] = it;
      it = iclass[i];
      iclass[i] = iclass[ip];
      iclass[ip] = it;

   }


   free_vector(dbl, 0, nr6);


}



/***********************************************************************
                           setgam()
************************************************************************/

void setgam(INT_T natom, REAL_T * gamma, REAL_T * m, REAL_T eta,
            REAL_T hrmax, INT_T ioseen, char name[])

#define PI  3.14159
{



/*    
   Setup the friction  matrix
*/

   REAL_T *exp = NULL;
   REAL_T factor, largest;
   INT_T i, j, index1, index2, nr3, test;
   char atom[5], cexp[5];
   char str[5];

   FILE *ffp;
   int mytaskid;
   mytaskid = get_mytaskid();


   nr3 = 3 * natom;

   /* Allocate an array for the hydrodynamic radii */

   exp = vector(0, natom);

   test = 0;
   if ((ffp = fopen("expfile", "r"))) {
      test = 1;
   }


   if (test == 1) {

     if(mytaskid == 0){
       printf("\nexpfile is present\n\n");
     }
     ffp = fopen("expfile", "r");

     for (i = 0; i < natom; i++) {
       fscanf(ffp, "%s", str);
       sscanf(str, "%s", atom);
       fscanf(ffp, "%s \n", str);
       sscanf(str, "%s", cexp);
       exp[i] = atof(cexp);
     }

     largest = exp[0];
     for (i = 1; i < natom; i++) {
       if (exp[i] > largest) {
         largest = exp[i];
       }
     }

     factor = hrmax / sqrt(largest);
     for (i = 0; i < natom; i++) {
       exp[i] = sqrt(exp[i]) * factor;
       if (name[i] == 'H') {
         if (exp[i] > 0.2)
           exp[i] = 0.2;
       }
     }
     fclose(ffp);
   } else {

     if(mytaskid == 0){
       printf("\nexpfile is not present\n\n");
       printf
         ("The hydrodynamic radii will be set to %2.1f except for the hydrogen atoms where they will set to  0.2 \n\n",
          hrmax);
     }
      for (i = 0; i < natom; i++) {
         exp[i] = hrmax;
         if (name[i] == 'H') {
            if (exp[i] > 0.2)
               exp[i] = 0.2;
         }
      }
   }

   if (ioseen != 0) {
     if(mytaskid == 0){
       printf
         ("These options will be added soon to nmode. The ioseen keyword will be changed to 0.\n");
     }
      ioseen = 0;
   }

   if (ioseen == 0) {

      /* Gamma is diagonal */

      factor = 6.0 * eta * PI;
      index1 = 0;
      for (j = 0; j < nr3; j++) {
         for (i = 0; i < j; i++) {
            gamma[index1] = 0.0;
            index1++;
         }
         index2 = j / 3;
         gamma[index1] = (factor * exp[index2]) / (m[index2 + 1]);
         index1++;
      }
   }
   free_vector(exp, 0, natom);
}

/***********************************************************************
                            leigensort()
************************************************************************/

INT_T leigensort(REAL_T * wr, REAL_T * wi, INT_T * indexsort, INT_T natom)
{

   int i, j, nr6, it, mine, itemp, nreal;
   REAL_T temp;
   REAL_T *awi = NULL;

   nr6 = 6 * natom;
   awi = vector(0, nr6);

   /* Let's sort the eigenvalues according to their imaginay parts. */

   /* Set the initial order */

   for (i = 0; i < nr6; i++) {
      indexsort[i] = i;
      awi[i] = fabs(wi[i]);
   }

   /* Sort while keeping track of where the elements go */

   for (i = 0; i < nr6 - 1; i++) {
      /* get pointer to minimal value of remaining elements */
      mine = i;

      for (j = i + 1; j < nr6; j++) {
         if (awi[j] < awi[mine]) {
            mine = j;
         }
      }

      /* exchange ip-th and i-th elements */

      temp = awi[i];
      awi[i] = awi[mine];
      awi[mine] = temp;

      /* exchange ip-th and i-th index */

      it = indexsort[i];
      indexsort[i] = indexsort[mine];
      indexsort[mine] = it;
   }

   i = 0;
   while (i < nr6) {
      it = indexsort[i];

      if (wi[it] != 0.0 && wr[it] == wr[indexsort[i + 1]]) {
         if (wi[it] < 0.0) {
            indexsort[i] = indexsort[i + 1];
            indexsort[i + 1] = it;
         }
         i = i + 2;
      } else {
         i = i + 1;
      }
   }

   /*  Let's count the reals   */
   /* awi is used to store abs(wr)   */

   nreal = 0;
   while (fabs(wi[indexsort[nreal]]) < 0.0000001 && nreal < nr6) {
      awi[nreal] = fabs(wr[indexsort[nreal]]);
      nreal = nreal + 1;
   }

   for (i = 0; i < nreal - 1; i++) {

      mine = i;
      for (j = i + 1; j < nreal; j++) {
         if (awi[j] < awi[mine]) {
            mine = j;
         }
      }

      temp = awi[i];
      awi[i] = awi[mine];
      awi[mine] = temp;
      itemp = indexsort[i];
      indexsort[i] = indexsort[mine];
      indexsort[mine] = itemp;
   }

   free_vector(awi, 0, nr6);

   return nreal;
}

/***********************************************************************
                            lnorm()
************************************************************************/

#ifdef CYGWIN
void lnorm(INT_T * indexsort, REAL_T * wr, REAL_T * wi, REAL_T * zo,
           REAL_T * gamma, INT_T * iclass, INT_T natom)
{
   int ik, nr6, i, j, k, l, po, indkl;
   complextype lambda, sum, cons;
   complextype piki, pili, elki, temp;

   nr6 = 6 * natom;
   ik = 0;
   po = 0;
   while (ik < nr6) {

      i = indexsort[ik];
      lambda.r = wr[i];
      lambda.i = wi[i];
      //lambda = wr[i] + wi[i] * I;

      /* lambda is a complex pair */

      if (wi[i] > 0.0000001) {

         sum.r=0.0;
         sum.i=0.0;
         //sum = 0.0;
         for (k = 0; k < (nr6 / 2); k++) {

            piki.r = zo[i * (nr6) + k];
            piki.i = zo[(i + 1) * (nr6) + k];
            //piki = zo[i * (nr6) + k] + zo[(i + 1) * (nr6) + k] * I;
            for (l = 0; l < (nr6 / 2); l++) {
               pili.r = zo[i * (nr6) + l];
               pili.i = zo[(i + 1) * (nr6) + l];
               //pili = zo[i * (nr6) + l] + zo[(i + 1) * (nr6) + l] * I;
               indkl = indexkl(k, l);
               // NOTE: Does order of ops matter between piki and pili?
               temp.r = piki.r * gamma[indkl];
               temp.i = piki.i * gamma[indkl];
               temp = complex_mult(temp,pili);
               sum = complex_add(sum,temp);
               //sum = sum + piki * (gamma[indkl]) * pili;
            }
            temp=complex_mult(piki,piki);
            temp=complex_mult(lambda,temp);
            temp.r+=2.00;
            sum = (complex_add(sum,temp));
            //sum = sum + 2.00 * (lambda * piki * piki);
         }
         if (sum.r==0.0 && sum.i==0.0) {
           cons.r=1.0;
           cons.i=0.0;
         } else {
           temp = complex_div(lambda,sum);
           cons = complex_sqrt(temp);
         }
         /*if (sum == 0.0) {
            cons = 1.0;
         } else {
            cons = csqrt(lambda / sum);
         }*/

         for (k = 0; k < (nr6 / 2); k++) {
            piki.r = zo[i * (nr6) + k];
            piki.i = zo[(i + 1) * (nr6) + k];
            //piki = zo[i * (nr6) + k] + zo[(i + 1) * (nr6) + k] * I;
            elki = complex_mult(piki,cons);
            //elki = piki * cons;
            zo[i * (nr6) + k] = elki.r;
            //zo[i * (nr6) + k] = creal(elki);
            // NOTE: Is this equivalent to -I * ?
            zo[(i + 1) * (nr6) + k] = -elki.r;
            //zo[(i + 1) * (nr6) + k] = creal(-I * elki);
         }

         iclass[ik] = 1;
         iclass[ik + 1] = 1;
         ik += 2;
      }

      /* lambda is positive real */

      else if (wr[i] > 0.0) {

         sum.r=0.0;
         sum.i=0.0;
         //sum = 0.0;

         for (k = 0; k < (nr6 / 2); k++) {
            for (l = 0; l < (nr6 / 2); l++) {

               indkl = indexkl(k, l);
               sum.r+=(zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) + l]);
               //sum = sum + zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) + l];
            }

            sum.r+=(2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]));
            //sum = sum + 2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]);

         }

         temp = complex_div(lambda,sum);
         cons = complex_sqrt(temp);
         //cons = csqrt(lambda / sum);

         for (k = 0; k < (nr6 / 2); k++) {

            // NOTE: Discarding the imaginary part?
            zo[i * (nr6) + k] = cons.r * zo[i * (nr6) + k];
            //temp.i = cons.i * zo[i * (nr6) + k];
            //zo[i * (nr6) + k] = zo[i * (nr6) + k] * cons;
         }
         iclass[ik] = 2;
         ik += 1;
      }

      /* lambda is negative real */

      else if (wr[i] < 0.0) {

         sum.r=0.0;
         sum.i=0.0;
         //sum = 0.0;
         for (k = 0; k < (nr6 / 2); k++) {
            for (l = 0; l < (nr6 / 2); l++) {
               indkl = indexkl(k, l);
               sum.r+=(zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) + l]); 
               //sum = sum + zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) + l];
            }
            sum.r+=(2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]));
            //sum += 2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]);
         }

         if (sum.r > 0.0) {
           temp.r=-lambda.r;
           temp.i=-lambda.i;
           temp = complex_div(temp,sum);
           cons = complex_sqrt(temp);
           iclass[ik] = 3;
         } else if (sum.r < 0.0) {
           temp=complex_div(lambda,sum);
           cons=complex_sqrt(temp);
           iclass[ik] = 4;
         }
         /*if (creal(sum) > 0.0) {
            cons = csqrt(-lambda / sum);
            iclass[ik] = 3;
         } else if (creal(sum) < 0.0) {
            cons = csqrt(lambda / sum);
            iclass[ik] = 4;
         }*/

         for (k = 0; k < (nr6 / 2); k++) {
            zo[i * (nr6) + k] = zo[i * (nr6) + k] * cons.r;
            //zo[i * (nr6) + k] = zo[i * (nr6) + k] * cons;
         }
         ik += 1;
      } else {
         iclass[ik] = 5;
         ik += 1;
      }
   }
}

#else //CYGWIN
void lnorm(INT_T * indexsort, REAL_T * wr, REAL_T * wi, REAL_T * zo,
           REAL_T * gamma, INT_T * iclass, INT_T natom)
{
   int ik, nr6, i, k, l, indkl;
   REAL_T complex lambda, sum, cons;
   REAL_T complex piki, pili, elki;

   nr6 = 6 * natom;
   ik = 0;
   while (ik < nr6) {

      i = indexsort[ik];
      lambda = wr[i] + wi[i] * I;

      /* lambda is a complex pair */

      if (wi[i] > 0.0000001) {

         sum = 0.0;
         for (k = 0; k < (nr6 / 2); k++) {

            piki = zo[i * (nr6) + k] + zo[(i + 1) * (nr6) + k] * I;
            for (l = 0; l < (nr6 / 2); l++) {
               pili = zo[i * (nr6) + l] + zo[(i + 1) * (nr6) + l] * I;
               indkl = indexkl(k, l);
               sum = sum + piki * (gamma[indkl]) * pili;
            }
            sum = sum + 2.00 * (lambda * piki * piki);
         }
         if (sum == 0.0) {
            cons = 1.0;
         } else {
            cons = csqrt(lambda / sum);
         }

         for (k = 0; k < (nr6 / 2); k++) {
            piki = zo[i * (nr6) + k] + zo[(i + 1) * (nr6) + k] * I;
            elki = piki * cons;
            zo[i * (nr6) + k] = creal(elki);
            zo[(i + 1) * (nr6) + k] = creal(-I * elki);
         }

         iclass[ik] = 1;
         iclass[ik + 1] = 1;
         ik += 2;
      }

      /* lambda is positive real */

      else if (wr[i] > 0.0) {

         sum = 0.0;

         for (k = 0; k < (nr6 / 2); k++) {
            for (l = 0; l < (nr6 / 2); l++) {

               indkl = indexkl(k, l);
               sum =
                   sum + zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) +
                                                               l];
            }

            sum =
                sum +
                2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]);

         }

         cons = csqrt(lambda / sum);

         for (k = 0; k < (nr6 / 2); k++) {

            zo[i * (nr6) + k] = zo[i * (nr6) + k] * cons;
         }
         iclass[ik] = 2;
         ik += 1;
      }

      /* lambda is negative real */

      else if (wr[i] < 0.0) {

         sum = 0.0;
         for (k = 0; k < (nr6 / 2); k++) {
            for (l = 0; l < (nr6 / 2); l++) {
               indkl = indexkl(k, l);
               sum =
                   sum + zo[i * (nr6) + k] * gamma[indkl] * zo[i * (nr6) + l];
            }
            sum += 2.00 * (wr[i] * zo[i * (nr6) + k] * zo[i * (nr6) + k]);
         }

         if (creal(sum) > 0.0) {
            cons = csqrt(-lambda / sum);
            iclass[ik] = 3;
         } else if (creal(sum) < 0.0) {
            cons = csqrt(lambda / sum);
            iclass[ik] = 4;
         }

         for (k = 0; k < (nr6 / 2); k++) {
            zo[i * (nr6) + k] = zo[i * (nr6) + k] * cons;
         }
         ik += 1;
      } else {
         iclass[ik] = 5;
         ik += 1;
      }
   }
}
#endif //CYGWIN
/***********************************************************************
                           diagonchol()
************************************************************************/


void diagonchol(REAL_T * d, REAL_T * v, REAL_T * h, INT_T nr3, INT_T nev)
{

   REAL_T t1, t2;
   int i, k, info;
   int po, po2;

   REAL_T *resid = NULL, *workd = NULL, *workl = NULL;
   REAL_T tol, addh, maddh;
   INT_T *iparam = NULL, *ipntr = NULL;
   INT_T ishifts, maxiter, mode, ncv, lworkl, ido, ldv, iter, one;
   char bmat[1];
   char ei[1];
   char which[2];
   /* dmy[] .. dummy variables needed for dseupd_() call */
   static short dmy[3];
   char uplo;

   /* for dseupd */
   int *select;
   int rvec;
   int mytaskid;
   mytaskid = get_mytaskid();

   select = (int *) malloc(nr3 * sizeof(int));
   uplo = 'U';

   t1 = seconds();

   addh = 0.0001;
   po = nr3;
   for (i = 0; i < nr3; i++) {
      k = i + i * po;
      h[k] = h[k] + addh;
   }


   dpotrf_(&uplo, &nr3, h, &nr3, &info);

   t2 = seconds();

   *tnmodeFact = t2 - t1;

   t1 = seconds();

   if (info != 0) {
      printf
          ("There was a problem during the factorization of the Hessian. You probably need a more minimized structure \n");
   }
   if (info != 0) {
      printf("The flag number is %d. Refer to the dpotrf documentation \n",
             info);
   }

   ido = 0;
   info = 0;
   tol = 0.0000000000000001;
   ncv = 2 * nev;
   ldv = nr3;
   lworkl = 2 * (3 * (ncv * ncv) + 6 * ncv);

   workd = vector(0, 3 * nr3);
   workl = vector(0, lworkl);
   resid = vector(0, nr3);

   ishifts = 1;
   maxiter = 3000;
   mode = 3;
   one = 1;

   iparam = ivector(0, 11);
   ipntr = ivector(0, 11);
   for (i = 0; i < 11; i++) {
      ipntr[i] = 0;
      iparam[i] = 0;
   }

   iparam[0] = ishifts;         /*shift */
   iparam[2] = maxiter;         /*maxiter */
   iparam[6] = mode;            /*mode */
   bmat[0] = 'I';
   which[0] = 'L';
   which[1] = 'A';

   iter = 0;

 numberofiter:iter = iter + 1;

   dsaupd_(&ido, bmat, &nr3, which, &nev, &tol, resid, &ncv, v, &ldv,
           iparam, ipntr, workd, workl, &lworkl, &info);

   if (info != 0) {
      printf
          ("There was a problem in dsaupd. The flag number is %d. Refer to the dsaupd documentation \n",
           info);
   }

   if (ido == -1 || ido == 1) {

      po = ipntr[0];
      po2 = ipntr[1];

      dcopy_(&nr3, &workd[po - 1], &one, &workd[po2 - 1], &one);
      dpotrs_(&uplo, &nr3, &one, h, &nr3, &workd[po2 - 1], &nr3, &info);

      if (info != 0) {
         printf
             ("There was a problem in dpotrs. The flag number is %d. Refer to the dpotrs documentation \n",
              info);
      }

      /*     iter =iter +1; */

      goto numberofiter;


   }

   t2 = seconds();

   *tnmodeInvdiag = t2 - t1;




   rvec = 1;
   ei[0] = 'A';
   maddh = -addh;



   t1 = seconds();

   dseupd_(&rvec, ei, select, d, v, &ldv, &maddh, bmat, &nr3, which, &nev,
           &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
           &lworkl, &info, &dmy[0], &dmy[1], &dmy[2]);

   t2 = seconds();

   if (info != 0) {
      printf
          ("There was a problem in dseupd. The flag number is %d. Refer to the dseupd documentation \n",
           info);
   }

   *tnmodeEigenvec = t2 - t1;



   if(mytaskid == 0){
     printf("Convergence achieved after  %d iterations \n", iter);
   }

   free(select);
   free_vector(resid, 0, nr3);
   free_vector(workd, 0, 3 * nr3);
   free_vector(workl, 0, lworkl);
   free_ivector(iparam, 0, 11);
   free_ivector(ipntr, 0, 11);

}

/***********************************************************************
                           diagondir()
************************************************************************/

void diagondir(INT_T natom, REAL_T * h)
{

   REAL_T t1, t2;
   int i, info;
   int po, po2, nr3;

   REAL_T *resid = NULL, *vv = NULL, *workd = NULL, *workl = NULL;
   REAL_T tol, alpha, zero;
   INT_T *iparam = NULL, *ipntr = NULL;
   INT_T ncv, nev, lworkl, ido, ldvv, iter, one;
   char mb[1];
   char which[2];
   char uplo;
   int mytaskid;
   mytaskid = get_mytaskid();

   nr3 = 3 * natom;
   t1 = seconds();

   iter = 1;
   ido = 0;
   info = 0;
   tol = 0.00000000001;
   ncv = 80;
   nev = 40;
   ldvv = nr3;
   lworkl = 2 * (3 * (ncv * ncv) + 6 * ncv);

   vv = vector(0, nr3 * ncv);
   workd = vector(0, 3 * nr3);
   workl = vector(0, lworkl);
   resid = vector(0, nr3);

   one = 1;

   iparam = ivector(0, 11);
   ipntr = ivector(0, 11);

   iparam[0] = 1;               /*shift */
   iparam[2] = 300;             /*maxiter */
   iparam[6] = 1;               /*mode */

   mb[0] = 'I';
   which[0] = 'S';
   which[1] = 'M';

   uplo = 'N';

   if(mytaskid == 0) printf("before lmode %d", iparam[7]);

 numberofiter:if(mytaskid == 0)printf("iter number %d \n", iter);

   dsaupd_(&ido, mb, &nr3, which, &nev, &tol, resid, &ncv, vv, &ldvv,
           iparam, ipntr, workd, workl, &lworkl, &info);

   if(mytaskid == 0){
     printf("through lmode");
     printf("info %d ido %d", info, ido);
   }

   alpha = 1.00;
   if (ido == -1 || ido == 1) {
      po = ipntr[0];
      po2 = ipntr[1];

      dgemv_(&uplo, &nr3, &nr3, &alpha, h, &nr3, &workd[po - 1], &one,
             &zero, &workd[po2 - 1], &one);

      if(mytaskid == 0){
        if (iter == 1) {
          printf("workd");
          for (i = 0; i <= 3 * nr3; i++) {
            printf("%d %f\n", i + 1, workd[i]);
          }
        }
        
        printf("through dpotrs");
      }
      iter = iter + 1;

      goto numberofiter;

   }

   t2 = seconds();
   *tnmodeInvdiag = t2 - t1;

   if(mytaskid == 0){
     for (i = 0; i < nev; i++) {
       printf("%f \n", ((workl[ipntr[5] + i - 1])));
     }
   }

   free_vector(resid, 0, nr3);
   free_vector(vv, 0, nr3 * ncv);
   free_vector(workd, 0, 3 * nr3);
   free_vector(workl, 0, lworkl);
   free_ivector(iparam, 0, 11);
   free_ivector(ipntr, 0, 11);

}

/***********************************************************************
                            NMODE()
************************************************************************/

/* 
 * Compute normal modes
 *
 * Calling parameters are as follows:
 *
 * x[n]    contains the coordinates
 * n       dimension of variables
 * func    pointer to function that computes the energy, gradient and hessian
 *           (probably only works now with mme2())
 * eigp    number of modes to post to file "vecs"
 *
 */

int nmode(REAL_T * x, int *n,
          REAL_T(*func) (REAL_T *, REAL_T *, REAL_T *, REAL_T *,
                         REAL_T *, int *, int *, int *,
                         int *, int *, int *, int *,
                         int *, int *, char *),
          int *eigp, int *ntrun, REAL_T * etau, REAL_T * hrmax,
          int *ioseen)
#define PI  3.14159
{
   REAL_T *g = NULL, *m = NULL, *v = NULL, *h = NULL, *work = NULL, *grad =
       NULL;
   REAL_T sumg, energy;
   REAL_T tnmode1, t1, t2, t0;
   int niter, i, j, k, mcopy, natom, ret_val, lwork, info;


   /* my additions */
   int po, nr3, nr6, ldu;
   int nreal;
   int index1, sizegamma;
   REAL_T eta;
   REAL_T *a = NULL, *gamma = NULL, *wr = NULL, *wi = NULL,
       *zo = NULL, *du = NULL, *xcom = NULL, *dnorm = NULL;
   INT_T *indexsort = NULL, *iclass = NULL;
   char *name;

   /* dsaupd       */

   REAL_T *d = NULL;
   REAL_T weight;
   INT_T *iwork = NULL;
   INT_T ncv, nev, liwork, iat;
   REAL_T mom[10];

   /* end of my additions */

   int mytaskid, numtasks, gridim;
   int context_PxQ = -1, context_1x1 = -1, context_Nx1 = -1;
   int descH_PxQ[DLEN_], descG_PxQ[DLEN_], descG_1x1[DLEN_];
   FILE *vfp;
   REAL_T pressure, temperature;
   char uplo, jobz;
   size_t ncopy;

#ifdef SCALAPACK
   int zero = 0, one = 1;
   int myrow, mycol, nprow, npcol;
   int myrowC, mycolC, nprowC, npcolC;
   int bs3, ierror, lld;
   int np, nq, nb, sizesytrd, nprocs, contextC, nrc, ldc, qrmem;
   int lldH_PxQ, lldG_PxQ, lldZ_PxQ, lldG_1x1;
   size_t locpH_PxQ, locpG_PxQ, locpG_1x1;
   size_t locqH_PxQ, locqG_PxQ, locqG_1x1;
   size_t sizeH_PxQ, sizeG_PxQ, sizeG_1x1;
   size_t adr, sizemqrleft;
   REAL_T *ptr, *reductarr = NULL, *eigvecrow = NULL;
#endif

   /* If PRINT_NM_TIMES is defined print some calculation times. */

#undef PRINT_NM_TIMES

   t0 = seconds();
   tnmode1 = t0;

   name = (char *) malloc((*n) * sizeof(char));
   assert( name );

   /* Define the type of run for the timer */

   typerun = *ntrun;

   /* Get mytaskid and, if SCALAPACK is defined, numtasks. */

   mytaskid = get_mytaskid();

#ifdef SCALAPACK
   numtasks = get_numtasks();
#endif

   /* Allocate some dynamic vectors and matrices. */

   mcopy = *n;
   ncopy = *n;
   m = vector(1, ncopy);

#ifndef SCALAPACK

   /* If SCALAPACK is not defined, allocate full copies of g and h. */

   g = vector(1, ncopy);
   h = vector(0, ncopy * ncopy);

#else

   /*
    * If SCALAPACK is defined, allocate distributed copies of g and h,
    * as well as a single copy of grad.
    *
    * Create a general context.  Although context_PxQ does comprise all
    * of the processes, it appears that the general context must be
    * distinct from the context(s) of the matrices that participate
    * in the redistribution via pdgemr2d.
    *
    * The topologic layout of the general context does not appear
    * to be important, so this general context is row cyclic.
    *
    * It appears to be important that the most general context be created
    * first, followed by successively less general contexts.  For this
    * code the correct order of creation is Nx1 then PxQ then 1x1.  Each
    * context is subsumed by the earlier created contexts.  I don't know
    * what to do about overlapping contexts where one does not sumsume
    * the other.  Failure to the most general context first leads to
    * synchronization errors.
    */

   sl_init_(&context_Nx1, &numtasks, &one);

   /* Calculate the dimensions of the largest possible square process grid. */

   gridim = (int) (sqrt((REAL_T) numtasks) + 0.5);
   if (gridim * gridim > numtasks) {
      gridim -= 1;
   }
   if (gridim == 0) {
      gridim = 1;
   }

   /*
    * Initialize the process grid for block cyclic distribution of matrices.
    * Note that a returned context of -1 indicates that the task is not
    * active on the process grid.
    */

   sl_init_(&context_PxQ, &gridim, &gridim);

   /* Initialize the process grid for a single (1x1) process. */

   sl_init_(&context_1x1, &one, &one);

   /*
    * Get the number of rows and columns on the block cyclic (PxQ) process grid,
    * as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);

   /*
    * Get the blocksize for a square block.  Because in egb2 the
    * loop index i selects three rows of the Hessian, multiply
    * the blocksize by three.
    */

   bs3 = 3 * get_blocksize();

   /*
    * If this task is on the process grid, set up the array descriptors.
    * If this task isn't on the process grid, set descZ_PxQ[CTXT_],
    * descG_PxQ[CTXT_] and descH_PxQ[CTXT_] to -1.  These values will
    * be used by pdgemr2d to determine activity on the grid.
    */

   if (context_PxQ >= 0) {

      /*
       * Allocate then distribute the Hessian matrix h, the eigenvector
       * matrix z and the gradient vector g on the block cyclic process grid.
       * The gradient vector g does not need to be allocated for columns other
       * than column 0, but the descinit_ function must be called for all columns.
       *
       * The numroc_ function is used to calculate the number of matrix
       * elements that are distributed across a PxQ processor grid.
       */

      locpG_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqG_PxQ = numroc_(&one, &bs3, &mycol, &zero, &npcol);
      sizeG_PxQ = locpG_PxQ * locqG_PxQ;
      lldG_PxQ = locpG_PxQ;
      descinit_(descG_PxQ, &mcopy, &one, &bs3, &bs3,
                &zero, &zero, &context_PxQ, &lldG_PxQ, &info);
      g = vector(1, sizeG_PxQ);

      locpH_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqH_PxQ = numroc_(&mcopy, &bs3, &mycol, &zero, &npcol);
      sizeH_PxQ = locpH_PxQ * locqH_PxQ;
      lldH_PxQ = locpH_PxQ;
      descinit_(descH_PxQ, &mcopy, &mcopy, &bs3, &bs3,
                &zero, &zero, &context_PxQ, &lldH_PxQ, &info);
      h = vector(0, sizeH_PxQ);

   } else {
      descG_PxQ[CTXT_] = -1;
      descH_PxQ[CTXT_] = -1;
   }

   /*
    * Get the number of rows and columns on the single process grid,
    * as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(&context_1x1, &nprow, &npcol, &myrow, &mycol);

   /*
    * If this task is on the process grid, set up the array descriptors.
    * If this task isn't on the process grid, set descG_1x1[CTXT_] ] to -1.
    * This value will be used by pdgemr2d to determine activity on the grid.
    */

   if (context_1x1 >= 0) {

      /*
       * Allocate then distribute the gradient vector grad on the single
       * process grid.  The descinit_ function is called for this array
       * for only the task that is active on the grid.
       *
       * Also, for the other tasks that are not active, the vector grad
       * is allocated because the copy of grad from the single process
       * grid is broadcast to the other copies.  The size of this vector
       * for these other tasks is determined by the ncopy variable.
       *
       * The numroc_ function is used to calculate the number of matrix
       * elements that are distributed across a 1x1 processor grid.
       */

      locpG_1x1 = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqG_1x1 = numroc_(&one, &bs3, &mycol, &zero, &npcol);
      sizeG_1x1 = locpG_1x1 * locqG_1x1;
      lldG_1x1 = locpG_1x1;
      descinit_(descG_1x1, &mcopy, &one, &bs3, &bs3,
                &zero, &zero, &context_1x1, &lldG_1x1, &info);
      grad = vector(1, sizeG_1x1);

   } else {
      descG_1x1[CTXT_] = -1;
      grad = vector(1, ncopy);
   }

   /* Allocate the eigenvector row and reduction arrays. */

   eigvecrow = vector(0, ncopy);
   reductarr = vector(0, ncopy);

#endif                          /* ifndef SCALAPACK */

   niter = 1;

   t1 = seconds();

   /*
    * For non-ScaLAPACK execution, set some variables to values that
    * will select the proper sections of code below.
    */

#ifndef SCALAPACK

   gridim = 1;
   context_PxQ = 0;

#endif

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

   /*
    * Compute the function value in "f",
    * its gradient (with respect to x) in "g", and
    * its hessian (with respect to x) in "h".
    * The atomic masses are returned in "m".
    * The number of atoms is returned in "natom".
    *
    * The grad, descG_PxQ, descG_1x1, descH_PxQ,
    * gridim, context_PxQ, context_1x1 and context_Nx1
    * calling parameters supply ScaLAPACK information,
    * or are dummy arguments for non-ScaLAPACK execution.
    *
    * The &g[1] and &m[1] calling parameters
    * map from 1..n indexing in this newton function
    * to 0..n-1 indexing in *func (which defaults to mme2).
    * This technique is not used for h because it must be
    * indexed from 1 to n.
    */



   energy = (*func) (x, &g[1], h, &m[1],
                     &grad[1], descG_PxQ, descG_1x1, descH_PxQ,
                     &context_PxQ, &context_1x1, &context_Nx1,
                     &gridim, &natom, &niter, name);



   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

   t2 = seconds();
   *tnmodeHessian += t2 - t1;
   t1 = t2;

   /*
    * Some ScaLAPACK functions appear to quit unexpectedly
    * for large matrices on a 1x1 process grid, e.g. the pdgemm_
    * function in newton.c, so bypass the ScaLAPACK pdsyev_ and
    * use the LAPACK dsyev_ function instead.
    *
    * The correct test is (gridim == 1), not (nprow == 1 && npcol == 1)
    * since processes that aren't on the 1x1 grid have nprow == npcol == -1,
    * which would direct control to pdgemr2d_ (below) that would hang because
    * it would not be called from all processes, specifically not from
    * the process that is on the 1x1 grid and has nprow == npcol == 1.
    *
    * Non-ScaLAPACK execution will select this section of code because
    * gridim and context_PxQ were set to 1 and 0, respectively, above.
    */

   if (gridim == 1) {
      if (context_PxQ >= 0) {

         sumg = 0.0;
         for (i = 1; i <= ncopy; i++) {
            sumg += g[i] * g[i];
         }
         sumg = sqrt(sumg / ((REAL_T) ncopy));
         if (mytaskid == 0) {
            fprintf( nabout, "\nEnergy       = %20.10e\n", energy);
            fprintf( nabout, "RMS gradient = %20.10e\n", sumg);
            fflush(nabout);
         }

         /*
          * Mass weight the Hessian:
          */

         j = 1;
         for (i = 1; i <= natom; i++) {
            g[j + 2] = g[j + 1] = g[j] = 1.0 / sqrt(m[i]);
            j += 3;
         }

         k = 0;
         for (i = 1; i <= ncopy; i++) {
            for (j = 1; j <= ncopy; j++) {
               h[k] = g[i] * h[k] * g[j];
               k++;
            }
         }


         nr3 = 3 * natom;


         /* Yannick J. Bomble */
         /* compute the Langevin modes if ntrun=1 and regular normal modes 
            if ntrun =0 */

         if (*ntrun == 3) {

            t2 = seconds();
	    *tnmodeOther += t2 - t1;
	    t1 = t2;
            nr3 = 3 * natom;
            nr6 = 6 * natom;

            /* Allocate some arrays */
            a = vector(0, nr6 * nr6);
            for (i = 0; i < nr6 * nr6; i++)
               a[i] = 0.0;


            /*  Set values for the hydrodynamic radius and solvent viscosity */


            eta = *etau * 2.945;
            sizegamma = sizegam(nr3);
            gamma = vector(0, sizegamma);

            setgam(natom, gamma, m, eta, *hrmax, *ioseen, name);


            /*
               for (i=0; i < sizegamma; i++){ 
               printf(" gam after gamma %d is %f \n",i,gamma[i]);
               }
             */

            po = nr3 * nr6 + nr3;

            for (j = 0; j < nr3; j++) {
               for (i = 0; i <= j; i++) {

                  index1 = (j + 1) * j / 2 + i;

                  a[po + i + (j) * nr6] = -gamma[index1];
                  a[po + j + (i) * nr6] = -gamma[index1];

               }
            }


            /* Put the Hessian in the lower left corner of the matrix A */

            po = nr3;

            for (i = 0; i < nr3; i++) {
               for (j = 0; j < nr3; j++) {
                  a[po + j] = -h[i * nr3 + j];
               }
               po += nr6;
            }


            /* Set the  upper right corner of the matrix A to unity */

            po = nr3 * nr6;

            for (i = 0; i < nr3; i++) {
               a[po + i] = 1.0;
               po += nr6;
            }

            /*
               printf("a dd is");
               for (i=0;i<nr6*nr6;i++){
               printf("%d   %f \n",i,a[i]);
               }
             */

            t2 = seconds();
            *tnmodeAmat = t2 - t1;
	    t1 = t2;
            /*
             * Allocate the work array, call dgeev_ and deallocate the work array.
             * This step requires that LAPACK and BLAS support be available in your
             * NAB installation.
             */

            ldu = 1;
            wr = vector(0, nr6);
            wi = vector(0, nr6);
            zo = vector(0, nr6 * nr6);

            if (*eigp)
               uplo = 'V';
            else
               uplo = 'N';
            jobz = 'N';
            lwork = 4 * nr6;
            work = vector(0, lwork);

            t2 = seconds();
	    *tnmodeOther += t2 - t1;
	    t1 = t2;
            dgeev_(&jobz, &uplo, &nr6, a, &nr6, wr, wi, du, &ldu, zo, &nr6,
                   work, &lwork, &info);
            t2 = seconds();
            *tnmodeDiagA += t2 - t1;
	    t1 = t2;

            free_vector(work, 0, lwork);
            free_vector(a, 0, nr6 * nr6);

            indexsort = ivector(0, nr6);
            t2 = seconds();
            *tnmodeOther += t2 - t1;
            t1 = t2;

            nreal = leigensort(wr, wi, indexsort, natom);

            t2 = seconds();
            *tnmodeSort += t2 - t1;
            t1 = t2;

            /* Compute the normalized Langevin modes. */
            iclass = ivector(0, nr6);

            lnorm(indexsort, wr, wi, zo, gamma, iclass, natom);

            xcom = vector(0, nr3);
            movecm(xcom, x, m, natom);
            dnorm = vector(0, nr6);

            remtranrot(natom, nreal, xcom, zo, wr, wi, indexsort, dnorm, m,
                       iclass);

            /*  mass un-weight the eigenvectors:  */
            po = 0;
            k = 0;
            for (i = 0; i < nr6; i++) {
               for (j = 0; j < nr3; j++) {
                  zo[po + k++] *= g[j + 1];
               }
               po += nr3;
            }

            t2 = seconds();
            *tnmodeNorm = t2 - t1;
            *tnmodeLan = t2 - t0;
	    t1 = t2;

            if(mytaskid == 0){
              printf
                (" Langevin relaxtion times and frequencies from NAB in picoseconds and cm**-1 respectively\n");
              printf("%5d\n", nr6);
              for (i = 0; i < nr6; i++) {
                
                if (i <= 5) {
                  printf("%7d %10.4f %10.4f  \n", i + 1, -wr[indexsort[i]],
                         wi[indexsort[i]] * 108.587);
                } else if (i < nreal) {
                  if (iclass[i] == 2 || iclass[i] == 4) {
                    printf("%7d %10.4f %10.4f %10.4f   real  \n", i + 1,
                           -1 / (20.455 * wr[indexsort[i]]),
                           wi[indexsort[i]] * 108.587,
                           dnorm[indexsort[i]]);
                  } else if (iclass[i] == 3) {
                    printf("%7d %10.4f %10.4f %10.4f   imaginary  \n",
                           i + 1, -1 / (20.455 * wr[indexsort[i]]),
                           wi[indexsort[i]] * 108.587,
                           dnorm[indexsort[i]]);
                  } else {
                    printf("%7d %10.4f %10.4f %10.4f  iclass is wrong \n",
                           i + 1, -1 / (20.455 * wr[indexsort[i]]),
                           wi[indexsort[i]] * 108.587,
                           dnorm[indexsort[i]]);
                  }
                }
                
                else {
                  printf("%7d %10.4f %10.4f              complex\n",
                         i + 1, -1 / (20.455 * wr[indexsort[i]]),
                         wi[indexsort[i]] * 108.587);
                }
              }

              vfp = safe_fopen("lmode", "w");
              fprintf(vfp, " Langevin frequencies from NAB in cm**-1\n");
              fprintf(vfp, "%5d\n", nr6);
              for (i = 0; i < nr6; i++) {
                
                if (i <= 5) {
                  fprintf(vfp, "%7d %10.4f %10.4f  \n", i + 1,
                          wr[indexsort[i]] * 108.587,
                          wi[indexsort[i]] * 108.587);
                } else if (i < nreal) {
                  if (iclass[i] == 2 || iclass[i] == 4) {
                    fprintf(vfp, "%7d %10.4f %10.4f %10.4f   real  \n",
                            i + 1, wr[indexsort[i]] * 108.587,
                            wi[indexsort[i]] * 108.587,
                            dnorm[indexsort[i]]);
                  } else if (iclass[i] == 3) {
                    fprintf(vfp,
                            "%7d %10.4f %10.4f %10.4f   imaginary  \n",
                            i + 1, wr[indexsort[i]] * 108.587,
                            wi[indexsort[i]] * 108.587,
                            dnorm[indexsort[i]]);
                  } else {
                    fprintf(vfp,
                            "%7d %10.4f %10.4f %10.4f  iclass is wrong \n",
                            i + 1, wr[indexsort[i]] * 108.587,
                            wi[indexsort[i]] * 108.587,
                            dnorm[indexsort[i]]);
                  }
                }
                
                else {
                  fprintf(vfp,
                          "%7d %10.4f %10.4f              complex a \n",
                          i + 1, wr[indexsort[i]] * 108.587,
                          wi[indexsort[i]] * 108.587);
                }
              }
              fclose(vfp);
            }

            /* Make an Amber-compatible formatted vecs file via task 0 */

            if (*eigp && mytaskid == 0) {
               vfp = safe_fopen("lmodevecs", "w");
               fprintf(vfp, " lmodes from NAB\n");
               fprintf(vfp, "%5d\n", nr3);
               for (i = 0; i < nr3; i++) {
                  fprintf(vfp, "%11.5f", x[i]);
                  if (i % 7 == 6)
                     fprintf(vfp, "\n");
               }
               if (nr3 % 7 != 0)
                  fprintf(vfp, "\n");
               for (j = 0; j < *eigp; j++) {
                  fprintf(vfp, "****\n%5d%12.5f%12.5f\n", j + 1,
                          wr[indexsort[j]] * 108.587,
                          wi[indexsort[j]] * 108.587);
                  for (i = 0; i < nr3; i++) {
                     fprintf(vfp, "%11.5f", zo[indexsort[j] * nr6 + i]);
                     if (i % 7 == 6)
                        fprintf(vfp, "\n");
                  }
                  if (nr3 % 7 != 0)
                     fprintf(vfp, "\n");
               }
               fclose(vfp);
            }

            free_vector(dnorm, 0, nr6);
            free_vector(xcom, 0, nr3);
            free_ivector(iclass, 0, nr6);
            free_ivector(indexsort, 0, nr6);
            free_vector(wr, 0, nr6);
            free_vector(wi, 0, nr6);
            free_vector(zo, 0, nr6 * nr6);
            free_vector(gamma, 0, sizegamma);

            /* end of the Langevin part of the code *ntrun = 3 */
         }


         else if (*ntrun == 0) {

            /*
             * Allocate the work array, call dsyev_ and deallocate the work array.
             * This step requires that LAPACK and BLAS support be available in your
             * NAB installation.
             */

            if (*eigp)
               jobz = 'V';
            else
               jobz = 'N';
            uplo = 'U';
            lwork = 3 * ncopy;
            work = vector(0, lwork);
            k = 0;
            v = vector(1, ncopy);

            t2 = seconds();
            *tnmodeOther += t2 - t1;
            t1 = t2;

            dsyev_(&jobz, &uplo, n, h, n, &v[1], work, &lwork, &info);

            t2 = seconds();
            *tnmodeDSYEV += t2 - t1;
            if(mytaskid == 0){
              printf("dysev time = %10.2f seconds \n\n", t2 - t1);
              fflush(nabout);
            }
	    t1 = t2;

            free_vector(work, 0, lwork);

            /* un-mass-weight the eigenvectors:  */

            k = 0;
            for (i = 1; i <= ncopy; i++) {
               for (j = 1; j <= ncopy; j++) {
                  h[k++] *= g[j];
               }
            }

            /*  get the eigenvalues into cm**-1 units:   */

            for (i = 1; i <= ncopy; i++) {
               v[i] =
                   v[i] >=
                   0.0 ? 108.587 * sqrt(v[i]) : -108.587 * sqrt(-v[i]);
            }

            /* Make an Amber-compatible formatted vecs file via task 0 */

            if (*eigp && mytaskid == 0) {
               vfp = safe_fopen("vecs", "w");
               fprintf(vfp, " modes from NAB\n");
               fprintf(vfp, "%5d\n", nr3);
               for (i = 0; i < ncopy; i++) {
                  fprintf(vfp, "%11.5f", x[i]);
                  if (i % 7 == 6)
                     fprintf(vfp, "\n");
               }
               if (ncopy % 7 != 0)
                  fprintf(vfp, "\n");
               k = 0;
               for (j = 1; j <= *eigp; j++) {
                  fprintf(vfp, "****\n%5d%12.5f\n", j, v[j]);
                  for (i = 0; i < ncopy; i++) {
                     fprintf(vfp, "%11.5f", h[k + i]);
                     if (i % 7 == 6)
                        fprintf(vfp, "\n");
                  }
                  if (ncopy % 7 != 0)
                     fprintf(vfp, "\n");
                  k += ncopy;
               }
               fclose(vfp);
            }

            /*  get the thermodynamic parameters  */

            j = 0;
            pressure = 1.0;
            temperature = 298.15;
            thermo(natom, mcopy, j, x, &m[1], &v[1], temperature, pressure,
                   energy);

            free_vector(v, 1, ncopy);

         }
      }
   }

   if (*ntrun == 2) {

      nr3 = ncopy;
      ncv = 2 * (*eigp);
      nev = *eigp;

      if(mytaskid == 0){
        printf
          ("\nOrder of the Hessian matrix  %d \nNumber of eigenvalues and eigenvectors sought %d \nNumber of Arnoldi vectors  %d \n \n",
           nr3, nev, ncv);
      }
      d = vector(0, nev);
      v = vector(0, nr3 * ncv);
      for (i = 0; i < nev; i++)
         d[i] = 0.0;
      for (i = 0; i < nr3 * ncv; i++)
         v[i] = 0.0;
      diagonchol(d, v, h, nr3, *eigp);

      /* un-mass-weight the eigenvectors:  */

      k = 0;
      for (i = 1; i <= *eigp; i++) {
         for (j = 1; j <= ncopy; j++) {
            v[k++] *= g[j];
         }
      }

      /*  get the eigenvalues into cm**-1 units:   */

      for (i = 0; i < *eigp; i++) {
         d[i] =
             d[i] >= 0.0 ? 108.587 * sqrt(d[i]) : -108.587 * sqrt(-d[i]);
      }

      /* Calculate mass, moment of inertia */

      pressure = 1.0;
      temperature = 298.15;

      if (mytaskid == 0) {
         fprintf(nabout, "\n                - Thermochemistry -\n\n");
         fprintf(nabout, "Temperature:  %8.3f\n   Pressure:  %8.3f\n",
                 temperature, pressure);
      }

/*     compute and print the molecular mass in amu */

      weight = 0.000;
      for (iat = 1; iat <= natom; ++iat) {
         weight += m[iat];
      }
      if (mytaskid == 0) {
         fprintf(nabout, "       Mass:  %8.3f\n", weight);
      }

/*     Compute the principal moments of inertia, */
/*       Note the imbedded conversion */
/*     of the moments to SI units. */


      mom[0] = 0.0;
      mofi(natom, x, &m[1], mom);
      /*     mofi(natom, &x[1], &m[1], mom);  */
      if (mytaskid == 0) {
         fprintf(nabout, "Principal moments of inertia in amu-A**2:\n");
         fprintf(nabout, "     %12.2f%12.2f%12.2f\n", mom[0], mom[1],
                 mom[2]);
         
         printf("\n");
         printf("Vibrational frequencies in cm**-1 \n");
         printf("\n");
         for (i = 0; i < nev; i++) {
           printf("%5d  %12.5f \n", i + 1, d[i]);
         }
      }

      /* Make an Amber-compatible formatted vecs file via task 0 */

      if (*eigp && mytaskid == 0) {
         vfp = safe_fopen("vecs", "w");
         fprintf(vfp, " modes from NAB\n");
         fprintf(vfp, "%5d\n", nr3);
         for (i = 0; i < ncopy; i++) {
            fprintf(vfp, "%11.5f", x[i]);
            if (i % 7 == 6)
               fprintf(vfp, "\n");
         }
         if (ncopy % 7 != 0)
            fprintf(vfp, "\n");
         k = 0;
         for (j = 0; j < *eigp; j++) {
            fprintf(vfp, "****\n%5d%12.5f\n", j + 1, d[j]);
            for (i = 0; i < ncopy; i++) {
               fprintf(vfp, "%11.5f", v[k + i]);
               if (i % 7 == 6)
                  fprintf(vfp, "\n");
            }
            if (ncopy % 7 != 0)
               fprintf(vfp, "\n");
            k += ncopy;
         }
         fclose(vfp);
      }

      free_vector(v, 0, nr3 * ncv);
      free_vector(d, 0, nev);

   }

   if (*ntrun == 1) {

      /*
       * Allocate the work array, call dsyevd_ and deallocate the work array.
       * This step requires that LAPACK and BLAS support be available in your
       * NAB installation.
       */

      if (*eigp)
         jobz = 'V';
      else
         jobz = 'N';
      uplo = 'U';
      lwork = 1 + 6 * ncopy + 2 * ncopy * ncopy;
      liwork = 3 + 5 * ncopy;
      work = vector(0, lwork);
      iwork = ivector(0, liwork);
      k = 0;
      v = vector(1, ncopy);

      t2 = seconds();
      *tnmodeOther += t2 - t1;
      t1 = t2;
      dsyevd_(&jobz, &uplo, n, h, n, &v[1], work, &lwork, iwork, &liwork,
              &info);
      t2 = seconds();
      *tnmodeDSYEVD += t2 - t1;
      t1 = t2;

      /* un-mass-weight the eigenvectors:  */

      k = 0;
      for (i = 1; i <= ncopy; i++) {
         for (j = 1; j <= ncopy; j++) {
            h[k++] *= g[j];
         }
      }

      /*  get the eigenvalues into cm**-1 units:   */

      for (i = 1; i <= ncopy; i++) {
         v[i] =
             v[i] >= 0.0 ? 108.587 * sqrt(v[i]) : -108.587 * sqrt(-v[i]);
      }

      /* Make an Amber-compatible formatted vecs file via task 0 */

      if (*eigp && mytaskid == 0) {
         vfp = safe_fopen("vecs", "w");
         fprintf(vfp, " modes from NAB\n");
         fprintf(vfp, "%5d\n", nr3);
         for (i = 0; i < ncopy; i++) {
            fprintf(vfp, "%11.5f", x[i]);
            if (i % 7 == 6)
               fprintf(vfp, "\n");
         }
         if (ncopy % 7 != 0)
            fprintf(vfp, "\n");
         k = 0;
         for (j = 1; j <= *eigp; j++) {
            fprintf(vfp, "****\n%5d%12.5f\n", j, v[j]);
            for (i = 0; i < ncopy; i++) {
               fprintf(vfp, "%11.5f", h[k + i]);
               if (i % 7 == 6)
                  fprintf(vfp, "\n");
            }
            if (ncopy % 7 != 0)
               fprintf(vfp, "\n");
            k += ncopy;
         }
         fclose(vfp);

      }

      /*  get the thermodynamic parameters  */

      j = 0;
      pressure = 1.0;
      temperature = 298.15;
      thermo(natom, mcopy, j, x, &m[1], &v[1], temperature, pressure,
             energy);

      free_vector(work, 0, lwork);
      free_vector(v, 1, ncopy);
      free_ivector(iwork, 0, liwork);

   } else {

#ifdef SCALAPACK

      /*
       * Here is the ScaLAPACK code that will execute if gridim != 1,
       * i.e., for more than one process on the process grid.
       *
       * Gather the distributed g vectors to the non-distributed grad vector.
       * Call the pdgemr2d_ function from all processes so that it doesn't hang.
       */

      pdgemr2d_(&mcopy, &one,
                &g[1], &one, &one, descG_PxQ,
                &grad[1], &one, &one, descG_1x1, &context_Nx1);

      /*
       * Calculate the RMS error for the gradient and print via the
       * process of the 1x1 grid.
       */

      if (context_1x1 >= 0) {
         sumg = 0.0;
         for (i = 1; i <= ncopy; i++) {
            sumg += grad[i] * grad[i];
         }
         sumg = sqrt(sumg / ((REAL_T) ncopy));
         if(mytaskid == 0){
           printf("RMS gradient = %12.5e\n", sumg);
           fflush(nabout);
         }
      }

      /*
       * Mass weight the Hessian:
       *
       * Get the processor row and column in context_PxQ, and  modify the
       * diagonal elements of the Hessian.  If the task is not active
       * in context_PxQ, the myroc function returns 0.
       */

      blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);

      if (context_PxQ >= 0) {
         j = 1;
         for (i = 1; i <= natom; i++) {
            grad[j + 2] = grad[j + 1] = grad[j] = 1.0 / sqrt(m[i]);
            j += 3;
         }

         for (i = 1; i <= ncopy; i++) {
            for (j = 1; j <= ncopy; j++) {
               ptr = ptr2d(h, descH_PxQ, i - 1, j - 1);
               if (ptr != NULL) {
                  *ptr = grad[i] * (*ptr) * grad[j];
               }
            }
         }

         /*
          *  Call the ScaLAPACK routines; requires that ScaLAPACK and BLAS support
          *  be available in your NAB installation.
          *
          * Calculate minimum size for the work array with or without calculation
          * of the eigenvectors.  See the pdsyev_ source code for details on how
          * to calculate LWMIN.
          */

         if (*eigp) {
            jobz = 'V';
            nprocs = nprow * npcol;
            contextC = sl_gridreshape_(&context_PxQ, &zero, &one,
                                       &one, &nprocs, &one);
            blacs_gridinfo_(&contextC, &nprowC, &npcolC, &myrowC, &mycolC);
            blacs_gridexit_(&contextC);
            nb = descH_PxQ[NB_];
            nrc = numroc_(&mcopy, &nb, &myrowC, &zero, &nprocs);
            ldc = (nrc > 1) ? nrc : 1;
            np = numroc_(&mcopy, &nb, &myrow, &zero, &nprow);
            nq = numroc_(&mcopy, &nb, &mycol, &zero, &npcol);
            sizemqrleft = nb * nb +
                ((nb * (nb - 1) / 2 >
                  (np + nq) * nb) ? nb * (nb - 1) / 2 : (np + nq) * nb);
            qrmem = 2 * ncopy - 2;
            lwork = 5 * ncopy + ncopy * ldc + 1 +
                ((sizemqrleft > qrmem) ? sizemqrleft : qrmem);
         } else {
            jobz = 'N';
            nb = descH_PxQ[NB_];
            np = numroc_(&mcopy, &nb, &myrow, &zero, &nprow);
            sizesytrd = (3 * nb > nb * (np + 1)) ? 3 * nb : nb * (np + 1);
            lwork = 5 * ncopy + sizesytrd + 1;
         }

         /*
          * Allocate the work array, call pdsyev_ and deallocate the work array.
          * This step requires that ScaLAPACK and BLACS support be available in
          * your NAB installation.
          */

         work = vector(0, lwork);
         v = vector(1, ncopy);

         uplo = 'U';
         t2 = seconds();
         *tnmodeOther += t2 - t1;
	 t1 = t2;
         pdsyev_(&jobz, &uplo, n,
                 h, &one, &one, descH_PxQ, &v[1],
                 h, &one, &one, descH_PxQ, work, &lwork, &info);

         t2 = seconds();
         *tnmodeOther += t2 - t1;
#ifdef PRINT_NR_TIMES
         if (mytaskid == 0) {
            printf("pdsyev time = %10.2f\n\n", t2 - t1);
            fflush(nabout);
         }
#endif
	 t1 = t2;

         free_vector(work, 0, lwork);

         /* un-mass-weight the eigenvectors:  */

         for (i = 1; i <= ncopy; i++) {
            for (j = 1; j <= ncopy; j++) {
               ptr = ptr2d(h, descH_PxQ, i - 1, j - 1);
               if (ptr != NULL) {
                  *ptr *= grad[j];
               }
            }
         }

         /* get the eigenvalues into cm**-1 units:   */

         for (i = 1; i <= ncopy; i++) {
            v[i] =
                v[i] >=
                0.0 ? 108.587 * sqrt(v[i]) : -108.587 * sqrt(-v[i]);
         }
      }

      /*
       * Make an Amber-compatible formatted vecs file.  Write only by task 0.
       * Because the elements of the eigenvector matrix are distributed across
       * all of the processes of context_PxQ, it is necessary to gather them
       * for writing by task 0.  This gather is accomplished by copying the
       * matrix elements one row at a time into the eigvecrow array in the following
       * manner.  If a matrix element exists in the submatrix associated with
       * a particular process, that element is copied to the eigvecrow array,  If
       * that element does not exist in the submatrix, a zero is placed in the
       * eigvecrow array.  Once a row of the eigenvector matrix has been processed
       * in this manner, a reduction is performed involving all of the eigvecrow
       * arrays, including those that aren't active on the process grid so that
       * pdgemr2d_ doesn't hang.
       */

      if (*eigp) {
         if (mytaskid == 0) {
            vfp = safe_fopen("vecs", "w");
            fprintf(vfp, " modes from NAB\n");
            fprintf(vfp, "%5d\n", ncopy);
            for (i = 0; i < ncopy; i++) {
               fprintf(vfp, "%11.5f", x[i]);
               if (i % 7 == 6)
                  fprintf(vfp, "\n");
            }
            if (ncopy % 7 != 0)
               fprintf(vfp, "\n");
         }
         for (j = 1; j <= *eigp; j++) {
            if (mytaskid == 0) {
               fprintf(vfp, "****\n%5d%12.5f\n", j, v[j]);
            }
            for (i = 0; i < ncopy; i++) {
               ptr = ptr2d(h, descH_PxQ, j - 1, i);
               if (ptr != NULL) {
                  eigvecrow[i] = *ptr;
               } else {
                  eigvecrow[i] = 0.0;
               }
            }
            ierror = MPI_Allreduce(eigvecrow, reductarr, ncopy,
                                   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS) {
               printf
                   ("Error in nmode work reduction, error = %d  mytaskid = %d\n",
                    ierror, mytaskid);
            }
            for (i = 0; i < ncopy; i++) {
               if (mytaskid == 0) {
                  fprintf(vfp, "%11.5f", reductarr[i]);
                  if (i % 7 == 6)
                     fprintf(vfp, "\n");
               }
            }
            if (mytaskid == 0) {
               if (ncopy % 7 != 0)
                  fprintf(vfp, "\n");
            }
         }
         if (mytaskid == 0) {
            fclose(vfp);
         }
      }

      /*  get the thermodynamic parameters */

      if (context_PxQ >= 0) {
         j = 0;
         pressure = 1.0;
         temperature = 298.15;
         thermo(natom, mcopy, j, x, &m[1], &v[1], temperature, pressure,
                energy);
      }
#endif                          /* ifdef SCALAPACK */

   }

#ifndef SCALAPACK

   /* Free the dynamic vectors and matrices. */

   if (g != NULL)
      free_vector(g, 1, ncopy);
   if (m != NULL)
      free_vector(m, 1, ncopy);

   if (ntrun == 0) {
      if (v != NULL)
         free_vector(v, 1, ncopy);
   }

   if (h != NULL)
      free_vector(h, 0, ncopy * ncopy);

#else

   /* Free the dynamic vectors and matrices. */

   if (context_PxQ >= 0) {
      if (g != NULL)
         free_vector(g, 1, sizeG_PxQ);
      if (h != NULL)
         free_vector(h, 0, sizeH_PxQ);
   }
   if (context_1x1 >= 0) {
      if (grad != NULL)
         free_vector(grad, 1, sizeG_1x1);
   } else {
      if (grad != NULL)
         free_vector(grad, 1, ncopy);
   }
   if (eigvecrow != NULL)
      free_vector(eigvecrow, 0, ncopy);
   if (reductarr != NULL)
      free_vector(reductarr, 0, ncopy);
   if (m != NULL)
      free_vector(m, 1, ncopy);
   if (v != NULL)
      free_vector(v, 1, ncopy);

   /* Exit the process grid only for tasks that are on the process grid. */

   if (context_PxQ >= 0) {
      blacs_gridexit_(&context_PxQ);
   }
   if (context_1x1 >= 0) {
      blacs_gridexit_(&context_1x1);
   }

   if (context_Nx1 >= 0) {
      blacs_gridexit_(&context_Nx1);
   }
#endif                          /* ifndef SCALAPACK */

   /* Call egb2 with an iteration count of -3 to free static vectors. */

   niter = -3;
   (*func) (x, &g[1], h, &m[1],
            &grad[1], descG_PxQ, descG_1x1, descH_PxQ,
            &context_PxQ, &context_1x1, &context_Nx1,
            &gridim, &natom, &niter, name);

   free(name);
   ret_val = 0;
   t2 = seconds();
   *tnmodeOther += t2 - t1;
   *tnmode = t2 - tnmode1;
   return ret_val;
}
