#ifndef SFF_H
#define SFF_H

/*  First stab at a header file for the visible routines in libsff.
 *
 *  NOTE: obviously not yet commented; furthermore, this is not yet
 *  tested, either.  Need to improve on this file before it would actually
 *  be used....
 */

#include <stdio.h>

#ifndef	FALSE
#define	FALSE	0
#endif

#ifndef	TRUE
#define	TRUE	1
#endif

#define	UNDEF	(-1)

	/* Fundamental sff types:	*/

typedef	int	INT_T;
typedef	size_t	SIZE_T;

#define NAB_DOUBLE_PRECISION 1
#define REAL_T	double

	/* other nab types:	*/
typedef	char	STRING_T;
typedef	FILE	FILE_T;

typedef struct parm {
	char	ititl[81];
	INT_T 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap, Numextra;
	STRING_T *AtomNames, *ResNames, *AtomSym, *AtomTree;
	REAL_T	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB10, *Rborn, *Fs;
	REAL_T	Box[4], Cutcap, Xcap, Ycap, Zcap, Fsmax;
	INT_T 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		*AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		*BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		*AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		*DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		*DihAt3, *DihAt4, *DihNum, *Boundary;
	INT_T	*N14pairs, *N14pairlist;
	REAL_T *Gvdw, *Scee, *Scnb;
	INT_T	Nstrand, Ncomplex;       /* HCP: number of strands/complexes */
	INT_T	*Ipstrand, *Ipcomplex;   /* HCP: index into residue/strand list */
} PARMSTRUCT_T;

struct xmin_opt{ INT_T  mol_struct_opt;
    INT_T  maxiter;
    REAL_T grms_tol;
    INT_T  method;
    INT_T  numdiff;
    INT_T  m_lbfgs;
    INT_T  iter;
    REAL_T xmin_time;
    INT_T  ls_method;
    INT_T  ls_maxiter;
    REAL_T ls_maxatmov;
    REAL_T beta_armijo;
    REAL_T c_armijo;
    REAL_T mu_armijo;
    REAL_T ftol_wolfe;
    REAL_T gtol_wolfe;
    INT_T  ls_iter;
    INT_T  print_level;
    INT_T  error_flag; 
};

struct lmod_opt {
    INT_T niter;
    INT_T nmod;
    INT_T kmod;
    INT_T nrotran_dof;
    INT_T nconf;
    REAL_T minim_grms;
    REAL_T energy_window;
    INT_T eig_recalc;
    INT_T ndim_arnoldi;
    INT_T lmod_restart;
    INT_T n_best_struct;
    INT_T mc_option;
    REAL_T rtemp;
    REAL_T lmod_step_size_min;
    REAL_T lmod_step_size_max;
    INT_T nof_lmod_steps;
    REAL_T lmod_relax_grms;
    INT_T nlig;
    INT_T apply_rigdock;
    INT_T nof_poses_to_try;
    INT_T random_seed;
    INT_T print_level;
    REAL_T lmod_time;
    REAL_T aux_time;
    INT_T error_flag;
};

   /* the output for all non error emissions and some error ones too */
extern FILE	*nabout;

   /*  signatures for the public routines:  */

#ifdef __cplusplus
extern "C" {
#endif

INT_T		conjgrad( REAL_T*, INT_T*, REAL_T*,
			  REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
			  REAL_T*, REAL_T*, INT_T* );
REAL_T      gauss( REAL_T*, REAL_T* );
INT_T		getxv( STRING_T**, INT_T*, REAL_T*, REAL_T*, REAL_T* );
INT_T           getxyz( STRING_T**, INT_T*, REAL_T* );
REAL_T		lmodC(INT_T*, INT_T*, INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 REAL_T*, INT_T*, REAL_T*, REAL_T*, REAL_T*, INT_T*, INT_T*,
                 INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*, REAL_T*, INT_T*,
                 INT_T*, INT_T*, INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 INT_T*, REAL_T*, REAL_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 INT_T*, INT_T* );
INT_T		md( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*,
		    REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ) );
INT_T		mdrat( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T ( *mme )() );
INT_T		mm_options( STRING_T* );
void            mm_set_checkpoint( STRING_T** );
REAL_T		mme( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme2( REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
		      REAL_T*, INT_T*, INT_T*, INT_T*,
		      INT_T *, INT_T *, INT_T *, INT_T *,
		      INT_T*, INT_T* , char* );
REAL_T		mme4( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme_rattle( REAL_T*, REAL_T*, INT_T* );
INT_T		mme_init_sff( PARMSTRUCT_T*, INT_T*, INT_T*, REAL_T*, FILE_T* );
INT_T		newton( REAL_T*, INT_T*, REAL_T*,
			REAL_T ( *func1 )( REAL_T*, REAL_T*, INT_T* ),
			REAL_T ( *func2 )( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					   REAL_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, char* ),
			REAL_T*, REAL_T*, INT_T* );
INT_T		nmode( REAL_T*, INT_T*,
			REAL_T ( *func)( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					 REAL_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, char*  ),
		       INT_T*, INT_T*, REAL_T*, REAL_T* , INT_T*);
INT_T		putxv( STRING_T**, STRING_T**, INT_T*, REAL_T*, REAL_T*, REAL_T* );
INT_T           putxyz( STRING_T**, INT_T*, REAL_T* );
REAL_T		rand2( void );
INT_T		sasad( REAL_T*, REAL_T*, REAL_T*, INT_T, REAL_T );
REAL_T	xmin( REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
             INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, struct xmin_opt* );
REAL_T          xminC(INT_T*, INT_T*, INT_T*, REAL_T*, INT_T*, INT_T*,
                      INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                      INT_T*, REAL_T*, INT_T*, INT_T*, INT_T*, INT_T*,
                      REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                      REAL_T*, REAL_T*, INT_T*, INT_T*);
INT_T            mme_rism_max_memory();
  /*  prmtop routine interfaces are in prm.h */

  /*  rism routine interfaces: */
  /* 3D RISM section  */
  typedef struct {
    REAL_T solvcut;
    REAL_T buffer;
    REAL_T grdspc[3];
    REAL_T solvbox[3];
    REAL_T tolerance;
    REAL_T mdiis_del;
    REAL_T mdiis_restart;
    REAL_T fcecut;
    INT_T closureOrder;
    INT_T ng3[3];
    INT_T rism;      /* non-zero if RISM is turned on */
    INT_T asympCorr;
    INT_T mdiis_nvec;
    INT_T mdiis_method;
    INT_T maxstep;
    INT_T npropagate;
    INT_T centering;
    INT_T zerofrc;
    INT_T apply_rism_force;
    INT_T polarDecomp;
    INT_T rismnrespa;
    INT_T fcestride;
    INT_T fcenbasis;
    INT_T fcecrd;
    INT_T saveprogress;
    INT_T ntwrism;
    INT_T verbose;
    INT_T progress;
    INT_T write_thermo; 
    /*This is an unused variable that aligns
      the type on eight byte boundaries*/
    INT_T padding;
  } RismData;

#ifdef RISMSFF
  void rism_force_( REAL_T*, REAL_T*, REAL_T*, INT_T* );
  void rism_setparam_( RismData*, INT_T*, STRING_T* ,
                       INT_T*, STRING_T* ,  INT_T*, STRING_T* ,  INT_T*, STRING_T* ,
                       INT_T*, STRING_T* ,  INT_T*, STRING_T* ,  INT_T*, STRING_T* ,
                       INT_T*, STRING_T* ,  INT_T*, STRING_T* ,  INT_T*, STRING_T* ,
                       INT_T*, INT_T*, INT_T*,
                       REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                       INT_T*, INT_T*);
  void rism_init_( INT_T*);
  void rism_list_param_();
  void rism_writesolvdist_(INT_T*);
  void rism_thermo_calc_();
  void rism_thermo_print_(INT_T*,REAL_T*);
  void rism_printtimer_();
  void rism_max_memory_();
#endif /*RISMSFF*/

#ifdef __cplusplus
}
#endif

/*  from arpack:  */
void arsecond_( double * );


#endif
