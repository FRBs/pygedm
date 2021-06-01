/*  -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer wg1, wg2, wga, wggc, wglism, wgcn, wgvn;
} modelflags_;

#define modelflags_1 modelflags_

struct mw_1_ {
    real rsun;
};

#define mw_1 (*(struct mw_1_ *) &mw_)

struct {
    real n1h1, h1, a1, f1, n2, h2, a2, f2, na, ha, wa, aa, fa;
} galparams_;

#define galparams_1 galparams_

struct {
    real harm[5], narm[5], warm[5], farm[5];
} armfactors_;

#define armfactors_1 armfactors_

struct {
    real negc0, fgc0;
} gcparms_;

#define gcparms_1 gcparms_

struct {
    real armpaths[5], armdistances[6];
} armpathlengths_;

#define armpathlengths_1 armpathlengths_

struct {
    real dx0, dy0, dz0;
} dxyz_;

#define dxyz_1 dxyz_

struct {
    integer nclumps, hitclumpflag[2000];
} clumps_;

#define clumps_1 clumps_

struct {
    real aldr, bldr, cldr, xldr, yldr, zldr, thetaldr, neldr0, fldr, alsb, 
	    blsb, clsb, xlsb, ylsb, zlsb, thetalsb, nelsb0, flsb, alhb, blhb, 
	    clhb, xlhb, ylhb, zlhb, thetalhb, nelhb0, flhb, xlpi, ylpi, zlpi, 
	    rlpi, drlpi, nelpi, flpi, dnelpi, dflpi;
} nelismparms_;

#define nelismparms_1 nelismparms_

struct {
    integer nvoids, hitvoidflag[2000];
} voids_;

#define voids_1 voids_

/* Initialized data */

struct {
    real e_1;
    } mw_ = { 8.5f };


/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__4 = 4;
static doublereal c_b73 = 4.;
static integer c__9 = 9;
static doublereal c_b112 = 1.66667;
static doublereal c_b249 = 1.2;
static doublereal c_b250 = -4.4;
static doublereal c_b256 = -.6;
static doublereal c_b258 = -1.2;
static doublereal c_b259 = .6;
static doublereal c_b262 = -2.2;
static doublereal c_b267 = .01;
static doublereal c_b268 = 6.666666666666667;
static doublereal c_b270 = 10.;
static doublereal c_b274 = .3989;
static doublereal c_b275 = .58823529411764708;
static doublereal c_b276 = .35294117647058826;
static doublereal c_b277 = .29411764705882354;

/* density.NE2001.f */
/* final version of NE2001 */
/* returns densities, F parameters and weights of the various components */
/* mods: */
/* 28 July 02: */
/*   put in 'if' statements in subroutine density_2001 so that */
/* 	function calls are not done if particular weights */
/* 	(wg1, wg2, etc.) are set to zero in gal01.inp */
/* 	This is (a) cleaner and (b) much more efficient if */
/* 	the clump or void component is flagged out. */
/* Subroutine */ int density_2001__(real *x, real *y, real *z__, real *ne1, 
	real *ne2, real *nea, real *negc, real *nelism, real *necn, real *
	nevn, real *f1, real *f2, real *fa, real *fgc, real *flism, real *fcn,
	 real *fvn, integer *whicharm, integer *wlism, integer *wldr, integer 
	*wlhb, integer *wlsb, integer *wloopi, integer *hitclump, integer *
	hitvoid, integer *wvoid)
{
    /* Initialized data */

    static logical first = TRUE_;

    extern doublereal ne_inner__(real *, real *, real *, real *), ne_outer__(
	    real *, real *, real *, real *);
    extern /* Subroutine */ int neclumpn_(real *, real *, real *, real *, 
	    real *, integer *);
    extern doublereal ne_gc__(real *, real *, real *, real *);
    extern /* Subroutine */ int get_parameters__(void);
    extern doublereal ne_lism__(real *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern /* Subroutine */ int nevoidn_(real *, real *, real *, real *, real 
	    *, integer *, integer *);
    extern doublereal ne_arms_log_mod__(real *, real *, real *, integer *, 
	    real *);

/* f2py intent(in) x, y, z */
/* f2py intent(out) F1, F2, Fa, Fgc, Flism, FcN, FvN, ne1, ne2, nea */
/* f2py intent(out) whicharm, wlism, wLDR, wLHB, wLSB, wLOOPI */
/* f2py intent(out) hitclump, hitvoid, wvoid, negc, nelism, necN, nevN */
/* ---------------------------------------------------------------------------- */
/*  Returns seven components of the free electron density of the */
/*  interstellar medium at Galactic location (x,y,z). */
/*  Calling arguments: */
/*  input: */
/* 	x, y, z = galactocentric location (kpc) */
/*       Right-handed coordinate system */
/*       x is in l=90 direction */
/*       y is in l=180 direction */
/*       The sun is at (x,y,z) = (0,R0,0) */
/*  output: */
/*    electron densities in cm^{-3}: */
/* 	ne1:	outer, thick disk */
/* 	ne2:	inner, thin disk (annular in form) */
/* 	nea:	spiral arms */
/* 	negc:   galactic center component */
/*       nelism: local ISM component */
/*       necN:   contribution from discrete 'clumps' */
/*       nevN:   contribution from voids */
/*    fluctuation parameters (one for each ne component): */
/*       F1, F2, Fa, Fgc, Flism, FcN, FvN */
/*    flags: */
/*       whicharm: which of the 5 arms x,y,z is in (0 for interarm region) */
/*          wlism: 1 if x,y,z is in any of the four LISM components */
/*           wLDR: 1 if in LDR, 0 if not */
/*           wLHB: 1 if in LHB, 0 if not */
/*           wLSB: 1 if in LSB, 0 if not */
/*         wLOOPI: 1 if in LoopI, 0 if not */
/*       (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR) */
/*       hitclump: clump number that x,y,z is in (0 if none) */
/*        hitvoid: void number that x,y,z is in (0 if none) */
/* 25 May 2002 */
/* based on routines from TC93 and test routines from 1999-2002 by JMC. */
/* ---------------------------------------------------------------------------- */
/* functions: */
/* subroutines needed: */
/* 	neclumpN */
/* 	nevoidN */
    if (first) {
/* ! get parameters first time through */
	get_parameters__();
	first = FALSE_;
    }
/* need to initiate values in case components are flagged out */
    *ne1 = 0.f;
    *ne2 = 0.f;
    *nea = 0.f;
    *negc = 0.f;
    *nelism = 0.f;
    *necn = 0.f;
    *nevn = 0.f;
    *wlism = 0.f;
    *wldr = 0.f;
    *wlhb = 0.f;
    *wlsb = 0.f;
    *wloopi = 0.f;
    *hitclump = 0;
    *hitvoid = 0;
    *wvoid = 0;
    *whicharm = 0;
    if (modelflags_1.wg1 == 1) {
	*ne1 = ne_outer__(x, y, z__, f1);
    }
    if (modelflags_1.wg2 == 1) {
	*ne2 = ne_inner__(x, y, z__, f2);
    }
    if (modelflags_1.wga == 1) {
	*nea = ne_arms_log_mod__(x, y, z__, whicharm, fa);
    }
    if (modelflags_1.wggc == 1) {
	*negc = ne_gc__(x, y, z__, fgc);
    }
    if (modelflags_1.wglism == 1) {
	*nelism = ne_lism__(x, y, z__, flism, wlism, wldr, wlhb, wlsb, wloopi)
		;
    }
    if (modelflags_1.wgcn == 1) {
	neclumpn_(x, y, z__, necn, fcn, hitclump);
    }
    if (modelflags_1.wgvn == 1) {
	nevoidn_(x, y, z__, nevn, fvn, hitvoid, wvoid);
    }
/*       write(21, "(3(f8.2,1x),5(f8.4,1x),i2)") */
/*    .       x,y,z,ne1,ne2,nea,nelism,negc, */
/*    .       whicharm */
    return 0;
} /* density_2001__ */

/* Subroutine */ int get_parameters__(void)
{
    /* Format strings */
    static char fmt_1020[] = "(7x,f8.0)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen), s_rsfe(cilist *), do_fio(integer *, 
	    char *, ftnlen), e_rsfe(void), f_clos(cllist *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 11, 0, 0, 0 };
    static cilist io___3 = { 0, 11, 0, 0, 0 };
    static cilist io___4 = { 0, 11, 0, fmt_1020, 0 };


/* ----------------------------------------------------------------------- */
/* control flags for turning components on and off: */
/* parameters of large-scale components (inner+outer+arm components): */
/* factors for controlling individual spiral arms: */
/* 	narm: 	multiplies electron density (in addition to the`fac' */
/* 		      quantities) */
/* 	warm:	arm width factors that multiply nominal arm width */
/* 	harm:	arm scale height factors */
/* 	farm:	factors that multiply n_e^2 when calculating SM */
    o__1.oerr = 0;
    o__1.ounit = 11;
    o__1.ofnmlen = 9;
    o__1.ofnm = "gal01.inp";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsle(&io___2);
    e_rsle();
    s_rsle(&io___3);
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wg1, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wg2, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wga, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wggc, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wglism, (ftnlen)sizeof(integer)
	    );
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wgcn, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&modelflags_1.wgvn, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsfe(&io___4);
    do_fio(&c__1, (char *)&galparams_1.n1h1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.h1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.a1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.f1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.n2, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.h2, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.a2, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.f2, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.na, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.ha, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.wa, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.aa, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&galparams_1.fa, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.narm[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.narm[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.narm[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.narm[3], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.narm[4], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.warm[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.warm[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.warm[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.warm[3], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.warm[4], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.harm[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.harm[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.harm[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.harm[3], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.harm[4], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.farm[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.farm[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.farm[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.farm[3], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&armfactors_1.farm[4], (ftnlen)sizeof(real));
    e_rsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 11;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*          write(6,*) 'get_parms: weights: ', */
/*    .           wg1, wg2, wga, wggc, wglism, wgcN, wgvN */
/*          write(6,*) 'get_parms: ', */
/*    .           n1h1,h1,A1,F1,n2,h2,A2,F2, */
/*    .           na,ha,wa,Aa,Fa, */
/*    .           narm(1), narm(2), narm(3), narm(4), narm(5), */
/*    .           warm(1), warm(2), warm(3), warm(4), warm(5), */
/*    .           harm(1), harm(2), harm(3), harm(4), harm(5), */
/*    .           farm(1), farm(2), farm(3), farm(4), farm(5) */
    return 0;
} /* get_parameters__ */

doublereal ne_arms_log_mod__(real *x, real *y, real *z__, integer *whicharm, 
	real *farms)
{
    /* Initialized data */

    static integer ks = 3;
    static integer armmap[5] = { 1,3,4,2,5 };
    static integer nnj[5] = { 20,20,20,20,20 };
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen), f_clos(cllist *);
    double exp(doublereal), cos(doublereal), sin(doublereal), atan2(
	    doublereal, doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, n;
    static real r__, r1[100]	/* was [20][5] */, ga;
    static integer jj, kk, kl;
    static real th, rr, sq, th1[100]	/* was [20][5] */, sq1, sq2, ebb, fac,
	     nea;
    static integer kma;
    static real arg, emm, arm[5000]	/* was [5][500][2] */, dth;
    static integer kmi;
    static real exx, eyy, th2a, th3a, th3b, th2b, aarm[5];
    static integer kmax[5];
    static real rmin[5], smin, test, thxy;
    extern doublereal sech2_(real *);
    static real test2, test3;
    static integer whicharm_spiralmodel__;
    static real thmin[5], sqmin, extent[5], fac2min, fac3min;
    extern /* Subroutine */ int cspline_(real *, real *, integer *, real *, 
	    real *);
    static real sminmin;

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 11, 0, 0, 0 };
    static cilist io___12 = { 0, 11, 0, 0, 0 };
    static cilist io___14 = { 0, 11, 0, 0, 0 };


/* ----------------------------------------------------------------------- */
/*  Spiral arms are defined as logarithmic spirals using the */
/*    parameterization in Wainscoat et al. 1992, ApJS, 83, 111-146. */
/*  But arms are modified selectively at various places to distort them */
/*    as needed (08 Aug 2000). */
/*  Note that arm numbering follows that of TC93 for the four large arms */
/* (after remapping). */
/*  The local spiral arm is number 5. */
/*  06 Apr 02:   removed TC type modifications of arms 2,3 (fac calculations) */
/*  		and replaced with new versions.  Data for these are hard wired. */
/* see get_parameters for definitions of narm, warm, harm. */
/* function: */
/* ! for remapping from Wainscoat */
/* ! order to TC93 order, which is */
/* ! from GC outwards toward Sun. */
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = *y;
    rr = sqrt(r__1 * r__1 + r__2 * r__2);
    if (first) {
/* ! Reconstruct spiral arm axes */
/* read arm parameters: */
	o__1.oerr = 0;
	o__1.ounit = 11;
	o__1.ofnmlen = 19;
	o__1.ofnm = "ne_arms_log_mod.inp";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
/*       write(6,*) 'ne_arms_log_mod.inp:' */
	s_rsle(&io___11);
	e_rsle();
	s_rsle(&io___12);
	e_rsle();
	for (j = 1; j <= 5; ++j) {
	    s_rsle(&io___14);
	    do_lio(&c__4, &c__1, (char *)&aarm[j - 1], (ftnlen)sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&rmin[j - 1], (ftnlen)sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&thmin[j - 1], (ftnlen)sizeof(real));
	    do_lio(&c__4, &c__1, (char *)&extent[j - 1], (ftnlen)sizeof(real))
		    ;
	    e_rsle();
/*         write(6,*) aarm(j), rmin(j), thmin(j), extent(j) */
	}
	cl__1.cerr = 0;
	cl__1.cunit = 11;
	cl__1.csta = 0;
	f_clos(&cl__1);
	for (j = 1; j <= 5; ++j) {
/* ! fill sampling array */
	    i__1 = nnj[j - 1];
	    for (n = 1; n <= i__1; ++n) {
		th1[n + j * 20 - 21] = thmin[j - 1] + (n - 1) * extent[j - 1] 
			/ (nnj[j - 1] - 1.f);
/* ! rad */
		r1[n + j * 20 - 21] = rmin[j - 1] * exp((th1[n + j * 20 - 21] 
			- thmin[j - 1]) / aarm[j - 1]);
		th1[n + j * 20 - 21] *= 57.2957795130823f;
/* ! deg */
/* *** begin sculpting spiral arm 2 == TC arm 3*** */
		if (armmap[j - 1] == 3) {
		    if (th1[n + j * 20 - 21] > 370.f && th1[n + j * 20 - 21] 
			    <= 410.f) {
			r1[n + j * 20 - 21] *= cos((th1[n + j * 20 - 21] - 
				390.f) * 180.f / 2291.831180523292f) * .04f + 
				1.f;
/*    .           (1. + 0.01*cos((th1(n,j)-390.)*180./(40.*rad))) */
		    }
		    if (th1[n + j * 20 - 21] > 315.f && th1[n + j * 20 - 21] 
			    <= 370.f) {
			r1[n + j * 20 - 21] *= 1.f - cos((th1[n + j * 20 - 21]
				 - 345.f) * 180.f / 3151.2678732195268f) * 
				.07f;
/*    .           (1.0 - 0.08*cos((th1(n,j)-345.)*180./(55.*rad))) */
		    }
		    if (th1[n + j * 20 - 21] > 180.f && th1[n + j * 20 - 21] 
			    <= 315.f) {
			r1[n + j * 20 - 21] *= cos((th1[n + j * 20 - 21] - 
				260.f) * 180.f / 7734.9302342661103f) * .16f 
				+ 1;
/*    ,           (1 + 0.13*cos((th1(n,j)-260.)*180./(135.*rad))) */
		    }
		}
/* *** begin sculpting spiral arm 4 == TC arm 2*** */
		if (armmap[j - 1] == 2) {
		    if (th1[n + j * 20 - 21] > 290.f && th1[n + j * 20 - 21] 
			    <= 395.f) {
			r1[n + j * 20 - 21] *= 1.f - cos((th1[n + j * 20 - 21]
				 - 350.f) * 180.f / 6016.0568488736417f) * 
				.11f;
/*    .            1. */
		    }
		}
/* *** end arm sculpting *** */
/* 	   write(6,*) j,n, th1(n,j), r1(n,j) */
	    }
	}
/*        open(11,file='log_arms.out', status='unknown') */
/*        write(11,*) 'arm  n   xa     ya' */
	for (j = 1; j <= 5; ++j) {
	    dth = 5.f / r1[j * 20 - 20];
	    th = th1[j * 20 - 20] - dth * .999f;
	    i__1 = -nnj[j - 1];
	    cspline_(&th1[j * 20 - 20], &r1[j * 20 - 20], &i__1, &th, &r__);
/* 	      write(6,*) 'doing arm ', j, ' with ', NNj(j), ' points', */
/*    .               dth */
/* 	      write(6,*) (th1(k,j), r1(k,j), k=1,NNj(j)) */
	    for (k = 1; k <= 499; ++k) {
		th += dth;
		if (th > th1[nnj[j - 1] + j * 20 - 21]) {
		    goto L20;
		}
		cspline_(&th1[j * 20 - 20], &r1[j * 20 - 20], &nnj[j - 1], &
			th, &r__);
		arm[j + (k + 500) * 5 - 2506] = -r__ * sin(th / 
			57.2957795130823f);
		arm[j + (k + 1000) * 5 - 2506] = r__ * cos(th / 
			57.2957795130823f);
/*                 write(11,"(1x,i2,1x,i3,1x,2(f7.3,1x))") */
/*     .              j,k,arm(j,k,1),arm(j,k,2) */
/* L10: */
	    }
L20:
	    kmax[j - 1] = k;
/* L21: */
	}
	cl__1.cerr = 0;
	cl__1.cunit = 11;
	cl__1.csta = 0;
	f_clos(&cl__1);
	first = FALSE_;
    }

/* Get spiral arm component:  30 do loop finds a coarse minimum distance */
/* from line of sight to arm; 40 do loop finds a fine minimum distance */
/* from line of sight to arm; line 35 ensures that arm limits are not */
/* exceeded; linear interpolation beginning at line 41 finds the */
/* minimum distance from line of sight to arm on a finer scale than gridding */
/* of arms allows (TJL) */
    nea = 0.f;
    ga = 0.f;
    *whicharm = 0;
    whicharm_spiralmodel__ = 0;
    sminmin = 1e10f;
    thxy = atan2(-(*x), *y) * 57.2957795130823f;
/* ! measured ccw from +y axis */
/* ! (different from tc93 theta) */
    if (thxy < 0.f) {
	thxy += 360.f;
    }
    if ((r__1 = *z__ / galparams_1.ha, dabs(r__1)) < 10.f) {
	for (j = 1; j <= 5; ++j) {
	    jj = armmap[j - 1];
	    sqmin = 1e10f;
	    i__1 = kmax[j - 1] - ks;
	    i__2 = (ks << 1) + 1;
	    for (k = ks + 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing 2nd power */
		r__1 = *x - arm[j + (k + 500) * 5 - 2506];
/* Computing 2nd power */
		r__2 = *y - arm[j + (k + 1000) * 5 - 2506];
		sq = r__1 * r__1 + r__2 * r__2;
		if (sq < sqmin) {
		    sqmin = sq;
		    kk = k;
		}
/* L30: */
	    }
/* L35: */
/* Computing MAX */
	    i__2 = kk - (ks << 1);
	    kmi = max(i__2,1);
/* Computing MIN */
	    i__2 = kk + (ks << 1), i__1 = kmax[j - 1];
	    kma = min(i__2,i__1);
	    i__2 = kma;
	    for (k = kmi; k <= i__2; ++k) {
/* Computing 2nd power */
		r__1 = *x - arm[j + (k + 500) * 5 - 2506];
/* Computing 2nd power */
		r__2 = *y - arm[j + (k + 1000) * 5 - 2506];
		sq = r__1 * r__1 + r__2 * r__2;
		if (sq < sqmin) {
		    sqmin = sq;
		    kk = k;
		}
/* L40: */
	    }
/* L41: */
	    if (kk > 1 && kk < kmax[j - 1]) {
/* Computing 2nd power */
		r__1 = *x - arm[j + (kk + 499) * 5 - 2506];
/* Computing 2nd power */
		r__2 = *y - arm[j + (kk + 999) * 5 - 2506];
		sq1 = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
		r__1 = *x - arm[j + (kk + 501) * 5 - 2506];
/* Computing 2nd power */
		r__2 = *y - arm[j + (kk + 1001) * 5 - 2506];
		sq2 = r__1 * r__1 + r__2 * r__2;
		if (sq1 < sq2) {
		    kl = kk - 1;
		} else {
		    kl = kk + 1;
		}
		emm = (arm[j + (kk + 1000) * 5 - 2506] - arm[j + (kl + 1000) *
			 5 - 2506]) / (arm[j + (kk + 500) * 5 - 2506] - arm[j 
			+ (kl + 500) * 5 - 2506]);
		ebb = arm[j + (kk + 1000) * 5 - 2506] - emm * arm[j + (kk + 
			500) * 5 - 2506];
/* Computing 2nd power */
		r__1 = emm;
		exx = (*x + emm * *y - emm * ebb) / (r__1 * r__1 + 1.f);
		test = (exx - arm[j + (kk + 500) * 5 - 2506]) / (arm[j + (kl 
			+ 500) * 5 - 2506] - arm[j + (kk + 500) * 5 - 2506]);
		if (test < 0.f || test > 1.f) {
		    exx = arm[j + (kk + 500) * 5 - 2506];
		}
		eyy = emm * exx + ebb;
	    } else {
		exx = arm[j + (kk + 500) * 5 - 2506];
		eyy = arm[j + (kk + 1000) * 5 - 2506];
	    }
/* Computing 2nd power */
	    r__1 = *x - exx;
/* Computing 2nd power */
	    r__2 = *y - eyy;
	    sqmin = r__1 * r__1 + r__2 * r__2;
	    smin = sqrt(sqmin);
/* ! Distance of (x,y,z) from arm axis */
/*           write(23,"(4(f5.2,1x),i2,1x,3(f8.3,1x))") */
/*    .        x,y,z,rr,j,exx,eyy,smin */
	    if (smin < galparams_1.wa * 3) {
/* ! If (x,y,z) is close to this */
/* Computing 2nd power */
		r__1 = smin / (armfactors_1.warm[jj - 1] * galparams_1.wa);
		ga = exp(-(r__1 * r__1));
/* ! arm, get the arm weighting factor */
		if (smin < sminmin) {
		    whicharm_spiralmodel__ = j;
		    sminmin = smin;
		}
		if (rr > galparams_1.aa) {
/* ! Galactocentric radial dependence of arms */
		    r__1 = (rr - galparams_1.aa) / 2.f;
		    ga *= sech2_(&r__1);
/* 		write(6,*) 'd99a: rr,Aa,sech2() = ', */
/*                 rr, Aa, sech2((rr-Aa)/2.0) */
		}
/* arm3 reweighting: */
		th3a = 320.f;
		th3b = 390.f;
		th3b = 370.f;
		th3a = 290.f;
		th3b = 363.f;
		th3b = 363.f;
		fac3min = 0.f;
		test3 = thxy - th3a;
		if (test3 < 0.f) {
		    test3 += 360.f;
		}
		if (jj == 3 && 0.f <= test3 && test3 < th3b - th3a) {
		    arg = (thxy - th3a) * 6.2831853f / (th3b - th3a);
/* 		fac = (3.0 + cos(arg))/4.0 */
		    fac = (fac3min + 1.f + (1.f - fac3min) * cos(arg)) / 2.f;
		    d__1 = (doublereal) fac;
		    fac = pow_dd(&d__1, &c_b73);
/* 		write(90,*) x, y, thxy, th3a, th3b, test3, fac */
		    ga *= fac;
		}
/* arm2 reweighting: */
/*    first: as in tc93 (note different definition of theta) */
		th2a = 35.f;
		th2b = 55.f;
		test2 = thxy - th2a;
		fac = 1.f;
		if (jj == 2 && 0.f <= test2 && test2 < th2b - th2a) {
		    fac = test2 / (th2b - th2a) + 1.f;
		    fac = 1.f;
/* !**** note turned off */
		    ga *= fac;
		}
		if (jj == 2 && test2 > th2b - th2a) {
		    fac = 2.f;
		    fac = 1.f;
/* !**** note turned off */
		    ga *= fac;
		}
/*    second:  weaken the arm in a short range: */
		th2a = 340.f;
		th2b = 370.f;
/* note fix does nothing if fac2min = 1.0 */
		fac2min = .1f;
		test2 = thxy - th2a;
		if (test2 < 0.f) {
		    test2 += 360.f;
		}
		if (jj == 2 && 0.f <= test2 && test2 < th2b - th2a) {
		    arg = (thxy - th2a) * 6.2831853f / (th2b - th2a);
		    fac = (fac2min + 1.f + (1.f - fac2min) * cos(arg)) / 2.f;
/* 		fac = fac**3.5 */
/* 		write(90,*) x, y, thxy, th2a, th2b, test2, fac */
		    ga *= fac;
		}
		r__1 = *z__ / (armfactors_1.harm[jj - 1] * galparams_1.ha);
		nea += armfactors_1.narm[jj - 1] * galparams_1.na * ga * 
			sech2_(&r__1);
/* Add this arm */
	    }
/* L50: */
	}
    }
    ret_val = nea;
    *farms = 0.f;
    if (whicharm_spiralmodel__ == 0) {
	*whicharm = 0;
    } else {
	*whicharm = armmap[whicharm_spiralmodel__ - 1];
/* remap arm number */
	*farms = galparams_1.fa * armfactors_1.farm[*whicharm - 1];
    }
    return ret_val;
} /* ne_arms_log_mod__ */

doublereal ne_outer__(real *x, real *y, real *z__, real *f_outer__)
{

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static real g1, rr, ne1;
    extern doublereal sech2_(real *);
    static real suncos;

/* ----------------------------------------------------------------------- */
/* Thick disk component: */
/* 	g1=sech2(rr/A1)/sech2(8.5/A1)		! TC93 function */
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = *y;
    rr = sqrt(r__1 * r__1 + r__2 * r__2);
    suncos = cos(mw_1.rsun * 1.5707963267948966f / galparams_1.a1);
    if (rr > galparams_1.a1) {
	g1 = 0.f;
    } else {
	g1 = cos(rr * 1.5707963267948966f / galparams_1.a1) / suncos;
    }
    r__1 = *z__ / galparams_1.h1;
    ne1 = galparams_1.n1h1 / galparams_1.h1 * g1 * sech2_(&r__1);
    ret_val = ne1;
    *f_outer__ = galparams_1.f1;
    return ret_val;
} /* ne_outer__ */

doublereal ne_inner__(real *x, real *y, real *z__, real *f_inner__)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static real g2, rr, ne2;
    extern doublereal sech2_(real *);
    static real rrarg;

/* ----------------------------------------------------------------------- */
/* Thin disk (inner Galaxy) component: */
/* (referred to as 'Galactic center component' in circa TC93 density.f) */
    g2 = 0.f;
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = *y;
    rr = sqrt(r__1 * r__1 + r__2 * r__2);
/* Computing 2nd power */
    r__1 = (rr - galparams_1.a2) / 1.8f;
    rrarg = r__1 * r__1;
    if (rrarg < 10.f) {
	g2 = exp(-rrarg);
    }
    r__1 = *z__ / galparams_1.h2;
    ne2 = galparams_1.n2 * g2 * sech2_(&r__1);
    ret_val = ne2;
    *f_inner__ = galparams_1.f2;
    return ret_val;
} /* ne_inner__ */

doublereal ne_gc__(real *x, real *y, real *z__, real *f_gc__)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen), f_clos(cllist *);
    double sqrt(doublereal);

    /* Local variables */
    static real rr, zz, hgc, arg, rgc, xgc, ygc, zgc;

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 11, 0, 0, 0 };
    static cilist io___69 = { 0, 11, 0, 0, 0 };
    static cilist io___73 = { 0, 11, 0, 0, 0 };
    static cilist io___75 = { 0, 11, 0, 0, 0 };
    static cilist io___77 = { 0, 11, 0, 0, 0 };
    static cilist io___78 = { 0, 11, 0, 0, 0 };


/* ----------------------------------------------------------------------- */
/*     Determines the contribution of the Galactic center to the free */
/*     electron density of the interstellar medium at Galactic location */
/*     (x,y,z).  Combine with `fluctuation' parameter to obtain the */
/*     scattering measure. */

/*     NOTE: This is for the hyperstrong scattering region in the */
/*     Galactic center.  It is distinct from the inner Galaxy */
/*     (component 2) of the TC93 model. */

/*     Origin of coordinate system is at Galactic center; the Sun is at */
/*     (x,y,z) = (0,R0,0), x is in l=90 direction */

/*     Based on Section 4.3 of Lazio & Cordes (1998, ApJ, 505, 715) */

/* Input: */
/* REAL X - location in Galaxy [kpc] */
/* REAL Y - location in Galaxy [kpc] */
/* REAL Z - location in Galaxy [kpc] */

/* COMMON: */
/* REAL NEGC0 - nominal central density */

/* PARAMETERS: */
/* REAL RGC - radial scale length of Galactic center density enhancement */
/* REAL HGC - z scale height of Galactic center density enhancement */

/* Output: */
/* REAL NE_GC - Galactic center free electron density contribution [cm^-3] */
/* ----------------------------------------------------------------------- */

/*     parameter (xgc=-0.010, ygc=0., zgc=-0.020) */
/*     parameter (rgc=0.145) */
/*     parameter (hgc=0.026) */
    ret_val = 0.f;
    *f_gc__ = 0.f;
    if (first) {
	o__1.oerr = 0;
	o__1.ounit = 11;
	o__1.ofnmlen = 9;
	o__1.ofnm = "ne_gc.inp";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	s_rsle(&io___68);
	e_rsle();
	s_rsle(&io___69);
	do_lio(&c__4, &c__1, (char *)&xgc, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&ygc, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&zgc, (ftnlen)sizeof(real));
	e_rsle();
	s_rsle(&io___73);
	do_lio(&c__4, &c__1, (char *)&rgc, (ftnlen)sizeof(real));
	e_rsle();
	s_rsle(&io___75);
	do_lio(&c__4, &c__1, (char *)&hgc, (ftnlen)sizeof(real));
	e_rsle();
	s_rsle(&io___77);
	do_lio(&c__4, &c__1, (char *)&gcparms_1.negc0, (ftnlen)sizeof(real));
	e_rsle();
	s_rsle(&io___78);
	do_lio(&c__4, &c__1, (char *)&gcparms_1.fgc0, (ftnlen)sizeof(real));
	e_rsle();
	cl__1.cerr = 0;
	cl__1.cunit = 11;
	cl__1.csta = 0;
	f_clos(&cl__1);
	first = FALSE_;
    }
/* GALACTOCENTRIC RADIUS */
/* Computing 2nd power */
    r__1 = *x - xgc;
/* Computing 2nd power */
    r__2 = *y - ygc;
    rr = sqrt(r__1 * r__1 + r__2 * r__2);
    if (rr > rgc) {
	return ret_val;
    }
/* Z-HEIGHT. */
/* truncate at 1/e point */
    zz = (r__1 = *z__ - zgc, dabs(r__1));
    if (zz > hgc) {
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = rr / rgc;
/* Computing 2nd power */
    r__2 = zz / hgc;
    arg = r__1 * r__1 + r__2 * r__2;
    if (arg <= 1.f) {
	ret_val = gcparms_1.negc0;
	*f_gc__ = gcparms_1.fgc0;
/*        write(21,*) 'ne_gc: rr,zz,arg,ne_gc,F_gc ', */
/*    .                rr, zz, arg, ne_gc, F_gc */
    }
    return ret_val;
} /* ne_gc__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%  cspline.f  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Subroutine */ int cspline_(real *x, real *y, integer *nn, real *xout, real 
	*yout)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real a, b, h__;
    static integer i__, k, n;
    static real p, u[20], y2[20], qn, un;
    static integer khi;
    static real sig;
    static integer klo;

    /* Fortran I/O blocks */
    static cilist io___82 = { 0, 6, 0, 0, 0 };
    static cilist io___83 = { 0, 6, 0, 0, 0 };
    static cilist io___96 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (*nn > 20) {
	s_wsle(&io___82);
	do_lio(&c__9, &c__1, " too many points to spline. Change parameter s"
		"tatement", (ftnlen)54);
	e_wsle();
	s_wsle(&io___83);
	do_lio(&c__9, &c__1, " in cspline", (ftnlen)11);
	e_wsle();
    }
    n = abs(*nn);
    if (*nn < 0) {
	y2[0] = 0.f;
	u[0] = 0.f;
	i__1 = n - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    sig = (x[i__] - x[i__ - 1]) / (x[i__ + 1] - x[i__ - 1]);
	    p = sig * y2[i__ - 2] + 2.f;
	    y2[i__ - 1] = (sig - 1.f) / p;
	    u[i__ - 1] = (((y[i__ + 1] - y[i__]) / (x[i__ + 1] - x[i__]) - (y[
		    i__] - y[i__ - 1]) / (x[i__] - x[i__ - 1])) * 6.f / (x[
		    i__ + 1] - x[i__ - 1]) - sig * u[i__ - 2]) / p;
/* L10: */
	}
	qn = 0.f;
	un = 0.f;
	y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.f);
	for (k = n - 1; k >= 1; --k) {
/* L20: */
	    y2[k - 1] = y2[k - 1] * y2[k] + u[k - 1];
	}
    }
    klo = 1;
    khi = n;
L30:
    if (khi - klo > 1) {
	k = (khi + klo) / 2;
	if (x[k] > *xout) {
	    khi = k;
	} else {
	    klo = k;
	}
	goto L30;
    }
    h__ = x[khi] - x[klo];
    if (h__ == 0.f) {
	s_wsle(&io___96);
	do_lio(&c__9, &c__1, "bad x input.", (ftnlen)12);
	e_wsle();
    }
    a = (x[khi] - *xout) / h__;
    b = (*xout - x[klo]) / h__;
/* Computing 3rd power */
    r__1 = a;
/* Computing 3rd power */
    r__2 = b;
/* Computing 2nd power */
    r__3 = h__;
    *yout = a * y[klo] + b * y[khi] + ((r__1 * (r__1 * r__1) - a) * y2[klo - 
	    1] + (r__2 * (r__2 * r__2) - b) * y2[khi - 1]) * (r__3 * r__3) / 
	    6.f;
    return 0;
} /* cspline_ */

doublereal sech2_(real *z__)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double exp(doublereal);

/* ----------------------------------------------------------------------- */
    ret_val = 0.f;
    if (dabs(*z__) < 20.f) {
/* Computing 2nd power */
	r__1 = 2.f / (exp(*z__) + exp(-(*z__)));
	ret_val = r__1 * r__1;
    }
    return ret_val;
} /* sech2_ */

/* 23456789012345678901234567890123456789012345678901234567890123456789012 */
/* 15  Jan 2003: fixed whicharm bug (c.f. TWJL email of 12 Dec '02) */
/* 24 June 2002: added calculations of path lengths through LISM components */
/* 26  May 2002: modified for NE2001 routines which are cleaned-up */
/*             versions of development routines. */
/* Nov 1999 - May 2002: development versions */
/* 1992-1993: TC93 version */
/* Subroutine */ int dmdsm_(real *l, real *b, integer *ndir, real *dmpsr, 
	real *dist, char *limit, real *sm, real *smtau, real *smtheta, real *
	smiso, ftnlen limit_len)
{
    /* Initialized data */

    static real rrmax = 50.f;
    static real zmax = 25.f;
    static real dmax__ = 50.f;
    static logical first = TRUE_;
    static real r0 = 8.5f;

    /* System generated locals */
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real lhb_path__, lhb_dist__, lsb_path__, ldr_path__, dstep_pc__;
    static integer whicharm;
    static real lsb_dist__, ldr_dist__;
    static integer hitclump;
    static real d__;
    static integer i__;
    static real r__, x, y, z__, cb, dd, cl, dm, ne, sb, sl, rr;
    extern /* Subroutine */ int density_2001__(real *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static real dm1, dm2, ne1, ne2, loopi_path__, sm1, sm2, loopi_dist__, fgc,
	     dma, nea, fcn, sma, fvn, dsm1, dsm2, dmgc, negc, dmcn, necn, 
	    dsma, smgc;
    static integer wlhb;
    static real smcn, dmvn, sm_sum1_last__;
    static integer wlsb, wldr;
    static real nevn, sm_sum2_last__, sm_sum3_last__, sm_sum4_last__, smvn, 
	    f1val, f2val, faval, dsmgc, dsmcn, flism, dstep, dtest, dsmvn;
    static integer wvoid, wlism, wtemp, nstep;
    static real dmlism, nelism, dmstep, smlism;
    static integer ncount, wloopi, wtotal;
    static real sm_sum1__, sm_sum2__, sm_sum3__, sm_sum4__;
    static integer hitvoid;
    static real sm_term__, dsmlism;

/* f2py intent(in) l,b,ndir */
/* f2py intent(out) limit,sm,smtau,smtheta,smiso */
/* f2py intent(inout) dmpsr,dist */
/*  Computes pulsar distance and scattering measure */
/*  from model of Galactic electron distribution. */
/*  Input: real l	galactic longitude in radians */
/*         real b	galactic latitude in radians */
/*         integer ndir  >= 0 calculates dist from dmpsr */
/*                       < 0 for dmpsr from dist */
/* Input or output: */
/* 	  real dmpsr	(dispersion measure in pc/cm^3) */
/*         real dist	(distance in kpc) */
/*  Output: */
/* 	  char*1 limit	(set to '>' if only a lower distance limit can be */
/* 			 given; otherwise set to ' ') */
/*         sm            (scattering measure, uniform weighting) (kpc/m^{20/3}) */
/*         smtau         (scattering measure, weighting for pulse broadening) */
/*         smtheta       (scattering measure, weighting for angular broadening */
/*                        of galactic sources) */
/* 	  smiso 	(scattering measure appropriate for calculating the */
/* 			isoplanatic angle at the source's location */
/*       parameter(alpha = 11./3.) */
/*       parameter(pi = 3.14159) */
/*       parameter(c_sm = (alpha - 3.) / 2. * (2.*pi)**(4.-alpha) ) */
/* parameters of large-scale components (inner+outer+arm components): */
/* factors for controlling individual spiral arms: */
/*       narm:   multiplies electron density (in addition to the`fac' */
/*                     quantities) */
/*       warm:   arm width factors that multiply nominal arm width */
/*       harm:   arm scale height factors */
/*       farm:   factors that multiply n_e^2 when calculating SM */
/* Large scale components: */
/* Galactic center: */
/* LISM: */
/* clumps: */
/* voids: */
/* subroutines needed: */
/* 	density_2001 (and those that it calls) in density.NE2001.f */
/*       scattering routines in scattering98.f */
/* other variables */
/* 	data rrmax/30.0/ */
/* 	data zmax/1.76/ */
/* 	data zmax/5.00/ */
    if (first) {
/* initial call to density routine to set variable values */
/* through read-in of parameter file: */
	x = 0.f;
	y = r0;
	z__ = 0.f;
	density_2001__(&x, &y, &z__, &ne1, &ne2, &nea, &negc, &nelism, &necn, 
		&nevn, &f1val, &f2val, &faval, &fgc, &flism, &fcn, &fvn, &
		whicharm, &wlism, &wldr, &wlhb, &wlsb, &wloopi, &hitclump, &
		hitvoid, &wvoid);
/*       write(6,*) 'ne1,ne2,negc,nelism,necN,nevN = ', */
/*    .              ne1,ne2,negc,nelism,necN,nevN */
	first = FALSE_;
    }
    sl = sin(*l);
    cl = cos(*l);
    sb = sin(*b);
    cb = cos(*b);
    *(unsigned char *)limit = ' ';
/* 	dstep=0.02 */
/*       dstep = min(h1, h2) / 10.       ! step size in terms of scale heights */
    dstep = .01f;
    if (*ndir < 0) {
	dtest = *dist;
    }
    if (*ndir >= 0) {
	dtest = *dmpsr / (galparams_1.n1h1 / galparams_1.h1);
    }
/* approximate test distanc */
    nstep = dtest / dstep;
/* approximate number of steps */
    if (nstep < 10) {
	dstep = dtest / 10;
    }
/*  Sum until dm is reached (ndir >= 0) or dist is reached (ndir < 0). */
/*  Guard against too few terms by counting number of terms (ncount) so that */
/*  routine will work for n_e models with large n_e near the Sun. */
/* make # steps >= 10 */
L5:
    dstep_pc__ = dstep * 1e3f;
    dm = 0.f;
    sm_sum1__ = 0.f;
    sm_sum2__ = 0.f;
    sm_sum3__ = 0.f;
    sm_sum4__ = 0.f;
    for (i__ = 1; i__ <= 6; ++i__) {
	armpathlengths_1.armpaths[i__ - 1] = 0.f;
	armpathlengths_1.armdistances[i__ - 1] = 0.f;
    }
    dm1 = 0.f;
    dm2 = 0.f;
    dma = 0.f;
    dmgc = 0.f;
    dmlism = 0.f;
    dmcn = 0.f;
    dmvn = 0.f;
    sm1 = 0.f;
    sm2 = 0.f;
    sma = 0.f;
    smgc = 0.f;
    smlism = 0.f;
    smcn = 0.f;
    smvn = 0.f;
    ldr_path__ = 0.f;
    lhb_path__ = 0.f;
    lsb_path__ = 0.f;
    loopi_path__ = 0.f;
    ldr_dist__ = 0.f;
    lhb_dist__ = 0.f;
    lsb_dist__ = 0.f;
    loopi_dist__ = 0.f;
    ncount = 0;
    d__ = dstep * -.5f;
    for (i__ = 1; i__ <= 99999; ++i__) {
	++ncount;
	d__ += dstep;
	r__ = d__ * cb;
	x = r__ * sl;
	y = r0 - r__ * cl;
	z__ = d__ * sb;
/* Computing 2nd power */
	r__1 = x;
/* Computing 2nd power */
	r__2 = y;
	rr = sqrt(r__1 * r__1 + r__2 * r__2);
	if (*ndir >= 0 && (d__ > dmax__ || dabs(z__) > zmax || rr > rrmax)) {
	    goto L20;
	}
	if (*ndir < 3) {
	    density_2001__(&x, &y, &z__, &ne1, &ne2, &nea, &negc, &nelism, &
		    necn, &nevn, &f1val, &f2val, &faval, &fgc, &flism, &fcn, &
		    fvn, &whicharm, &wlism, &wldr, &wlhb, &wlsb, &wloopi, &
		    hitclump, &hitvoid, &wvoid);
	}
	if (*ndir >= 3) {
	    r__1 = x + dxyz_1.dx0;
	    r__2 = y + dxyz_1.dy0;
	    r__3 = z__ + dxyz_1.dz0;
	    density_2001__(&r__1, &r__2, &r__3, &ne1, &ne2, &nea, &negc, &
		    nelism, &necn, &nevn, &f1val, &f2val, &faval, &fgc, &
		    flism, &fcn, &fvn, &whicharm, &wlism, &wldr, &wlhb, &wlsb,
		     &wloopi, &hitclump, &hitvoid, &wvoid);
	}
/* wlism = 1 causes the lism component to override smooth Galactic components */
/* wvoid = 1 overrides everything except clumps */
	ne = (1.f - modelflags_1.wglism * wlism) * (modelflags_1.wg1 * ne1 + 
		modelflags_1.wg2 * ne2 + modelflags_1.wga * nea + 
		modelflags_1.wggc * negc) + modelflags_1.wglism * wlism * 
		nelism;
	ne = (1 - modelflags_1.wgvn * wvoid) * ne + modelflags_1.wgvn * wvoid 
		* nevn + modelflags_1.wgcn * necn;
	dmstep = dstep_pc__ * ne;
	dm += dmstep;
	wtotal = (1 - modelflags_1.wgvn * wvoid) * (1 - modelflags_1.wglism * 
		wlism);
	dm1 += wtotal * modelflags_1.wg1 * ne1;
	dm2 += wtotal * modelflags_1.wg2 * ne2;
	dma += wtotal * modelflags_1.wga * nea;
	dmgc += wtotal * modelflags_1.wggc * negc;
	dmlism += (1.f - modelflags_1.wgvn * wvoid) * modelflags_1.wglism * 
		wlism * nelism;
	dmcn += modelflags_1.wgcn * necn;
	dmvn += modelflags_1.wgvn * wvoid * nevn;
/*         write(24,"('n:',7f10.6,1x))") */
/*    .        ne1,ne2,nea,negc,nelism,necN,nevN */
/*        write(24,"(i2,1x,7(f10.5,1x))") */
/*    .      wtotal,dm1,dm2,dma,dmgc,dmlism,dmcN,dmvN */
/*         sm_term = */
/*    .       (1.-wglism*wlism)* */
/*    .       (wg1   * F1  * ne1**2 + */
/*    .        wg2   * F2  * ne2**2 + */
/*    .        wga   * Fa  * nea**2 + */
/*    .        wggc  * Fgc * negc**2) + */
/*    .        wglism*wlism * Flism * nelism**2 */
/* 	  sm_clumps = FcN * necN**2 */
/* 	  sm_voids  = FvN * nevN**2 */
/*         sm_term = (1-wgvN*wvoid) * sm_term */
/*    .            + wgvN * wvoid * sm_voids */
/*    .            + wgcN * sm_clumps */
/* Computing 2nd power */
	r__1 = ne1;
	dsm1 = wtotal * modelflags_1.wg1 * (r__1 * r__1) * galparams_1.f1;
/* Computing 2nd power */
	r__1 = ne2;
	dsm2 = wtotal * modelflags_1.wg2 * (r__1 * r__1) * galparams_1.f2;
/* Computing 2nd power */
	r__1 = nea;
	dsma = wtotal * modelflags_1.wga * (r__1 * r__1) * galparams_1.fa;
/* Computing 2nd power */
	r__1 = negc;
	dsmgc = wtotal * modelflags_1.wggc * (r__1 * r__1) * fgc;
/* Computing 2nd power */
	r__1 = nelism;
	dsmlism = (1.f - modelflags_1.wgvn * wvoid) * modelflags_1.wglism * 
		wlism * (r__1 * r__1) * flism;
/* Computing 2nd power */
	r__1 = necn;
	dsmcn = modelflags_1.wgcn * (r__1 * r__1) * fcn;
/* Computing 2nd power */
	r__1 = nevn;
	dsmvn = modelflags_1.wgvn * wvoid * (r__1 * r__1) * fvn;
	sm_term__ = dsm1 + dsm2 + dsma + dsmgc + dsmlism + dsmcn + dsmvn;
	sm1 += dsm1;
	sm2 += dsm2;
	sma += dsma;
	smgc += dsmgc;
	smlism += dsmlism;
	smcn += dsmcn;
	smvn += dsmvn;
	sm_sum1__ += sm_term__;
	sm_sum2__ += sm_term__ * d__;
/* Computing 2nd power */
	r__1 = d__;
	sm_sum3__ += sm_term__ * (r__1 * r__1);
	d__1 = (doublereal) d__;
	sm_sum4__ += sm_term__ * pow_dd(&d__1, &c_b112);
/* pathlengths through LISM components: */
/* take into account the weighting hierarchy, LHB:LOOPI:LSB:LDR */
	if (wlism == 1) {
	    if (wlhb == 1) {
		lhb_path__ += dstep;
		lhb_dist__ += d__;
	    }
	    if (wloopi == 1) {
		wtemp = 1 - wlhb;
		loopi_path__ += wtemp * dstep;
		loopi_dist__ += wtemp * d__;
	    }
	    if (wlsb == 1) {
		wtemp = (1 - wlhb) * (1 - wloopi);
		lsb_path__ += wtemp * dstep;
		lsb_dist__ += wtemp * d__;
	    }
	    if (wldr == 1) {
		wtemp = (1 - wlhb) * (1 - wloopi) * (1 - wlsb);
		ldr_path__ += wtemp * dstep;
		ldr_dist__ += wtemp * d__;
	    }
	}
/* pathlengths: whicharm = 0,5 (currently). */
/* 	                  1,4 for the equivalent of the TC93 arms */
/*                         5   for the local arm */
/*                         0   means interarm paths */
	armpathlengths_1.armpaths[whicharm] += dstep;
	armpathlengths_1.armdistances[whicharm] += d__;
/*       write(99,"(2(f8.3,1x), 7f10.6)") */
/*    .     d, dm, sm_term,  sm_sum1, sm_sum2, sm_sum3, */
/*    .     sm_sum1_last, sm_sum2_last, sm_sum3_last */
	if (*ndir >= 0 && dm >= *dmpsr) {
	    goto L30;
	}
	if (*ndir < 0 && d__ >= *dist) {
	    goto L40;
	}
	sm_sum1_last__ = sm_sum1__;
	sm_sum2_last__ = sm_sum2__;
	sm_sum3_last__ = sm_sum3__;
	sm_sum4_last__ = sm_sum4__;
/* L10: */
    }
    s_stop("loop limit", (ftnlen)10);
L20:
    *(unsigned char *)limit = '>';
    *dist = d__ - dstep * .5f;
    goto L999;
L30:
    *dist = d__ + dstep * .5f - dstep * (dm - *dmpsr) / dmstep;
    if (ncount < 10) {
	dstep /= 10.f;
	goto L5;
    }
    goto L999;
L40:
    *dmpsr = dm - dmstep * (d__ + dstep * .5f - *dist) / dstep;
    if (ncount < 10) {
	dstep /= 10.f;
	goto L5;
    }
L999:
/* normalize the mean distances: */
    if (ldr_path__ > 0.f) {
	ldr_dist__ /= ldr_path__ / dstep;
    }
    if (lhb_path__ > 0.f) {
	lhb_dist__ /= lhb_path__ / dstep;
    }
    if (lsb_path__ > 0.f) {
	lsb_dist__ /= lsb_path__ / dstep;
    }
    if (loopi_path__ > 0.f) {
	loopi_dist__ /= loopi_path__ / dstep;
    }
    dd = d__ + dstep * .5f - *dist;
/* subtract dd from armpath for latest arm (or iterarm) at end of LOS */
/*       armpaths(whicharm) = armpaths(whicharm)-dd */
    armpathlengths_1.armpaths[whicharm] -= dd;
    for (i__ = 1; i__ <= 6; ++i__) {
/* Computing MAX */
	r__1 = 1.f, r__2 = armpathlengths_1.armpaths[i__ - 1] / dstep;
	armpathlengths_1.armdistances[i__ - 1] /= dmax(r__1,r__2);
    }
    dm1 *= dstep_pc__;
    dm2 *= dstep_pc__;
    dma *= dstep_pc__;
    dmgc *= dstep_pc__;
    dmlism *= dstep_pc__;
    dmcn *= dstep_pc__;
    dmvn *= dstep_pc__;
/*       dsm = sm_term * (d+0.5*dstep - dist) */
/*       dsm = sm_term * dd */
/*       sm_sum2 = sm_sum2 - dsm * d */
/*       sm_sum3 = sm_sum3 - dsm * d**2 */
/*       sm_sum4 = sm_sum4 - dsm * d**1.67 */
/*       sm_sum1 = sm_sum1 - dsm */
/*       write(99,*) 'dmdsm: sm_term, sm_sum1, sm_sum1_last = ', */
/*    .    sm_term, sm_sum1, sm_sum1_last */
/* 	write(6,*) 'dmdsm: dsum1, sm_term = ', */
/*    .     sm_sum1-sm_sum1_last, sm_term */
    sm_sum1__ -= dd * (sm_sum1__ - sm_sum1_last__) / dstep;
    sm_sum2__ -= dd * (sm_sum2__ - sm_sum2_last__) / dstep;
    sm_sum3__ -= dd * (sm_sum3__ - sm_sum3_last__) / dstep;
    sm_sum4__ -= dd * (sm_sum4__ - sm_sum4_last__) / dstep;
/*       sm_sum2 = sm_sum2 - dsm * dist */
/*       sm_sum3 = sm_sum3 - dsm * dist**2 */
/*       sm_sum4 = sm_sum4 - dsm * dist**1.67 */
    *sm = dstep * 1.8389599999999999f * sm_sum1__;
/* Computing 2nd power */
    r__1 = *dist;
    *smtau = dstep * 11.033759999999999f * (sm_sum2__ / *dist - sm_sum3__ / (
	    r__1 * r__1));
/* Computing 2nd power */
    r__1 = *dist;
    *smtheta = dstep * 5.5168799999999996f * (sm_sum1__ + sm_sum3__ / (r__1 * 
	    r__1) - sm_sum2__ * 2.f / *dist);
    *smiso = dstep * 1.8389599999999999f * sm_sum4__;
    sm1 = sm1 * 1.8389599999999999f * dstep;
    sm2 = sm2 * 1.8389599999999999f * dstep;
    sma = sma * 1.8389599999999999f * dstep;
    smgc = smgc * 1.8389599999999999f * dstep;
    smlism = smlism * 1.8389599999999999f * dstep;
    smcn = smcn * 1.8389599999999999f * dstep;
    smvn = smvn * 1.8389599999999999f * dstep;
    return 0;
} /* dmdsm_ */

/* Subroutine */ int neclumpn_(real *x, real *y, real *z__, real *necn, real *
	fcn, integer *hitclump)
{
    /* Initialized data */

    static logical first = TRUE_;
    static integer luclump = 11;

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen);
    double sin(doublereal), cos(doublereal);
    integer f_clos(cllist *);
    double exp(doublereal);

    /* Local variables */
    static integer j, clumpflag;
    static real bc[2000], dc[2000], fc[2000], lc[2000], rc[2000], xc[2000], 
	    yc[2000], zc[2000], cbc, clc, nec[2000], sbc, arg, slc;
    static integer edge[2000];
    static real rgalc;

    /* Fortran I/O blocks */
    static cilist io___190 = { 0, 0, 0, 0, 0 };
    static cilist io___191 = { 0, 0, 1, 0, 0 };


/* returns electron density necN and fluctuation parameter FcN */
/* at position designated by l,b,d,x,y,z c for a set of */
/* clumps with parameters read in from file  neclumpN.dat */
/* input: */
/* 	x,y,z	coordinates	(kpc)  (as in TC93) */

/* output: */
/* 	necN	electron density in clump at (x,y,z) */
/* 	FcN	fluctuation parameter */
/* 	hitclump = 0:   no clump hit */
/* 		   j>0: j-th clump hit */
/* 	character*15 losname(nclumpsmax) */
/* 	character*1 type(nclumpsmax) */
/* parameters: */
/* 	lc	= galactic longitude of clump center */
/* 	bc	= galactic latitude of clump center */
/* 	(xc,yc,zc) = clump center location (calculated) */
/*       nec	= internal peak electron density */
/* 	rc	= clump radius at 1/e */
/*       Fc      = clump fluctuation parameter */
/* 	edge    = 0 => use exponential rolloff out to 5rc */
/*                 1 => uniform and truncated at 1/e */
/* first time through, read input clump parameters and calculate */
/* LOS quantities. */
/* lc,bc = Galactic coordinates (deg) */
/*   nec = clump electron density (cm^{-3}) */
/*    Fc = fluctuation parameter */
/*    dc = clump distance from Earth (kpc) */
/*    rc = clump radius (kpc) */
/*  edge = 0,1  0=> Gaussian, 1=> Gaussian w/ hard edge at e^{-1} */
/*  type = LOS type (P pulsar, G other Galactic, X extragalactic */
/* losname = useful name */
    if (first) {
/* read clump parameters */
	j = 1;
/* 	  write(6,*) 'reading neclumpN.NE2001.dat' */
	o__1.oerr = 0;
	o__1.ounit = luclump;
	o__1.ofnmlen = 19;
	o__1.ofnm = "neclumpN.NE2001.dat";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	io___190.ciunit = luclump;
	s_rsle(&io___190);
	e_rsle();
/* label line */
L5:
	io___191.ciunit = luclump;
	i__1 = s_rsle(&io___191);
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&clumpflag, (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&lc[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bc[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&nec[j - 1], (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&fc[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dc[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&rc[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&edge[j - 1], (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L99;
	}
	if (clumpflag == 0) {
	    slc = sin(lc[j - 1] / 57.29577951f);
	    clc = cos(lc[j - 1] / 57.29577951f);
	    sbc = sin(bc[j - 1] / 57.29577951f);
	    cbc = cos(bc[j - 1] / 57.29577951f);
	    rgalc = dc[j - 1] * cbc;
	    xc[j - 1] = rgalc * slc;
	    yc[j - 1] = 8.5f - rgalc * clc;
	    zc[j - 1] = dc[j - 1] * sbc;
/* 	  write(6,"(a15,1x,8(f8.3,1x))") */
/*    .           losname(j),lc(j),bc(j),dc(j), */
/*    .           nec(j),Fc(j),xc(j),yc(j),zc(j) */
	    ++j;
	}
	goto L5;
L99:
	first = FALSE_;
	clumps_1.nclumps = j - 1;
	cl__1.cerr = 0;
	cl__1.cunit = luclump;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    *necn = 0.f;
    *hitclump = 0;
    *fcn = 0.f;
    i__1 = clumps_1.nclumps;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = *x - xc[j - 1];
/* Computing 2nd power */
	r__2 = *y - yc[j - 1];
/* Computing 2nd power */
	r__3 = *z__ - zc[j - 1];
/* Computing 2nd power */
	r__4 = rc[j - 1];
	arg = (r__1 * r__1 + r__2 * r__2 + r__3 * r__3) / (r__4 * r__4);
	if (edge[j - 1] == 0 && arg < 5.f) {
	    *necn += nec[j - 1] * exp(-arg);
	    *fcn = fc[j - 1];
	    *hitclump = j;
	    clumps_1.hitclumpflag[j - 1] = 1;
	}
	if (edge[j - 1] == 1 && arg <= 1.f) {
/*    	    necN = necN + nec(j) * exp(-arg) */
	    *necn += nec[j - 1];
	    *fcn = fc[j - 1];
	    *hitclump = j;
	    clumps_1.hitclumpflag[j - 1] = 1;
	}
    }
    return 0;
} /* neclumpn_ */

/* routines to calculate the electron density for the */
/* Local Interstellar Medium */

/* JMC 26 August-11 Sep. 2000 */
/*     25 October 2001: modified to change weighting scheme */
/*                      so that the ranking is LHB: LSB: LDR */
/*                      (LHB overrides LSB and LDR; LSB overrides LDR) */
/*     16 November 2001: added Loop I component with weighting scheme */
/* 		        LHB:LOOPI:LSB:LDR */
/* 		        LHB   overides everything, */
/* 			LOOPI overrides LSB and LDR */
/* 			LSB   overrides LDR */
/* 			LISM  overrides general Galaxy */
/*     20 November 2001: The LOOPI component is truncated below z=0 */

/* after discussions with Shami Chatterjee */
/* the sizes, locations and densities of the LISM components */
/* are based largely on work by Toscano et al. 1999 */
/* and Bhat et al. 1999 */
doublereal ne_lism__(real *x, real *y, real *z__, real *flism, integer *wlism,
	 integer *wldr, integer *wlhb, integer *wlsb, integer *wloopi)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real ret_val;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen);

    /* Local variables */
    static real nelhbxyz, nelsbxyz, neldrq1xyz, neloopixyz, flhbr;
    extern doublereal nelsb_(real *, real *, real *, real *, integer *);
    static real flsbr;
    extern doublereal nelhb2_(real *, real *, real *, real *, integer *), 
	    neldrq1_(real *, real *, real *, real *, integer *);
    static real fldrq1r;
    extern doublereal neloopi_(real *, real *, real *, real *, integer *);
    static real floopir;

    /* Fortran I/O blocks */
    static cilist io___210 = { 0, 11, 0, 0, 0 };
    static cilist io___211 = { 0, 11, 0, 0, 0 };
    static cilist io___212 = { 0, 11, 0, 0, 0 };
    static cilist io___213 = { 0, 11, 0, 0, 0 };
    static cilist io___214 = { 0, 11, 0, 0, 0 };
    static cilist io___215 = { 0, 11, 0, 0, 0 };
    static cilist io___216 = { 0, 11, 0, 0, 0 };
    static cilist io___217 = { 0, 11, 0, 0, 0 };
    static cilist io___218 = { 0, 11, 0, 0, 0 };
    static cilist io___219 = { 0, 11, 0, 0, 0 };
    static cilist io___220 = { 0, 11, 0, 0, 0 };
    static cilist io___221 = { 0, 11, 0, 0, 0 };
    static cilist io___222 = { 0, 11, 0, 0, 0 };


/* functions: */
/* other variables: */
/* 'r' for returned value */
    if (first) {
/* read parameters for LISM */
	o__1.oerr = 0;
	o__1.ounit = 11;
	o__1.ofnmlen = 10;
	o__1.ofnm = "nelism.inp";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	s_rsle(&io___210);
	e_rsle();
	s_rsle(&io___211);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.aldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.bldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.cldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___212);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.yldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___213);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetaldr, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.neldr0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.fldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___214);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.alsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.blsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.clsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___215);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___216);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetalsb, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.nelsb0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.flsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___217);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.alhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.blhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.clhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___218);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___219);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetalhb, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.nelhb0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.flhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___220);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlpi, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___221);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.rlpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.drlpi, (ftnlen)sizeof(
		real));
	e_rsle();
	s_rsle(&io___222);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.nelpi, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.dnelpi, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.flpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.dflpi, (ftnlen)sizeof(
		real));
	e_rsle();
/*     write(99,*) xlpI,ylpI,zlpI */
/* 	  write(99,*) rlpI,drlpI */
/* 	  write(99,*) nelpI,dnelpI,FlpI,dFlpI */
	first = FALSE_;
    }
    neldrq1xyz = neldrq1_(x, y, z__, &fldrq1r, wldr);
/* low density region in Q1 */
    nelsbxyz = nelsb_(x, y, z__, &flsbr, wlsb);
/* Local Super Bubble */
    nelhbxyz = nelhb2_(x, y, z__, &flhbr, wlhb);
/* Local Hot Bubble */
    neloopixyz = neloopi_(x, y, z__, &floopir, wloopi);
/* weight the terms so that the LHB term overrides the other */
/* terms (we want the density to be low in the LHB, lower than */
/* in the other terms. */
/* Loop I */
    ret_val = (1 - *wlhb) * ((1 - *wloopi) * (*wlsb * nelsbxyz + (1 - *wlsb) *
	     neldrq1xyz) + *wloopi * neloopixyz) + *wlhb * nelhbxyz;
    *flism = (1 - *wlhb) * ((1 - *wloopi) * (*wlsb * flsbr + (1 - *wlsb) * 
	    fldrq1r) + *wloopi * floopir) + *wlhb * flhbr;
/* return the maximum weight of any of the terms for */
/* combining with additional terms external to this routine. */
/* Computing MAX */
/* Computing MAX */
    i__3 = *wldr, i__4 = max(*wlsb,*wlhb);
    i__1 = *wloopi, i__2 = max(i__3,i__4);
    *wlism = max(i__1,i__2);
/* temporary next 3 lines: */
/* 	ne_LISM = nelhbxyz */
/* 	flism = flhb */
/* 	wlism = wlhb */
/* 	write(97,"(9(f8.3,1x))") wLDR, wLSB, wLHB, */
/*    .      nelsbxyz, neldrq1xyz,nelhbxyz, */
/*    .      flsb, fldrq1, flhb */
    return ret_val;
} /* ne_lism__ */

doublereal neldrq1_(real *x, real *y, real *z__, real *fldrq1r, integer *
	wldrq1)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static real netrough, c__, q, s, aa, bb, cc, ap, bp, cp, dp, theta, 
	    ftrough;

/* Low Density Region in Q1 */
/* input: */
/* 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00 */
/* output: */
/* 	neLDRQ1 = electron density in local hot bubble that */
/* 	        is modeled as an ellipsoidal trough. */
/* 	FLDRQ1 = fluctuation parameter */
/* 	wLDRQ1  = weight of LDRQ1 component used to combine */
/* 		with other components of electron density. */
/* 		wLDRQ1 =  1  at and inside the annular ridge */
/* 		     <  1  outside the annular ridge */
/* 	             -> 0  far outside the annular ridge */
/* 	e.g. total electron density would be evaluated as */
/*            ne = (1-wLDRQ1)*ne_other + neLDRQ1 */
/* scales of ellipsoidal ridge */
/* ne of annulus, trough */
/* 	real xldr, yldr, zldr		! center of ellipsoid */
/* fluctuation parameters */
/*    measured from x axis */
/*    (x axis points toward l=90) */
/* position angle of major axis, */
/* 	data aa, bb, cc       	/ 0.6, 0.40, 0.3 / 	! GUESS */
/* 	data xldr, yldr, zldr 	/ 0.6, 7.86, 0. / 	! GUESS */
/* 	data theta		/ -45. / 		! GUESS */
/* 	data netrough  /0.010/		! GUESS */
/* 	data Ftrough   / 2./		! GUESS */
    aa = nelismparms_1.aldr;
    bb = nelismparms_1.bldr;
    cc = nelismparms_1.cldr;
    theta = nelismparms_1.thetaldr;
    netrough = nelismparms_1.neldr0;
    ftrough = nelismparms_1.fldr;
    if (first) {
	s = sin(theta / 57.29577951f);
	c__ = cos(theta / 57.29577951f);
/* Computing 2nd power */
	r__1 = c__ / aa;
/* Computing 2nd power */
	r__2 = s / bb;
	ap = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = s / aa;
/* Computing 2nd power */
	r__2 = c__ / bb;
	bp = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = cc;
	cp = 1.f / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * 2.f * s * (1.f / (r__1 * r__1) - 1.f / (r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = 0.f;
    *wldrq1 = 0;
    *fldrq1r = 0.f;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xldr;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.yldr;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zldr;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xldr) * (*y - nelismparms_1.yldr) * dp;
    if (q <= 1.f) {
/* inside */
	ret_val = netrough;
	*fldrq1r = ftrough;
	*wldrq1 = 1;
    }
    return ret_val;
} /* neldrq1_ */

doublereal nelsb_(real *x, real *y, real *z__, real *flsbr, integer *wlsb)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static real netrough, c__, q, s, aa, bb, cc, ap, bp, cp, dp, theta, 
	    ftrough;

/* Local Super Bubble */
/* input: */
/* 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00 */
/* output: */
/* 	neLSB = electron density in local hot bubble that */
/* 	        is modeled as an ellisoidal trough. */
/* 	FLSB = fluctuation parameter */
/* 	wLSB  = weight of LSB component used to combine */
/* 		with other components of electron density. */
/* 		wLSB =  1  at and inside the annular ridge */
/* 		     <  1  outside the annular ridge */
/* 	             -> 0  far outside the annular ridge */
/* 	e.g. total electron density would be evaluated as */
/*            ne = (1-wLSB)*ne_other + neLSB */
/* scales of ellipsoidal ridge */
/* ne of annulus, trough */
/* 	real xlsb, ylsb, zlsb		! center of ellipsoid */
/* fluctuation parameters */
/*    measured from x axis */
/*    (x axis points toward l=90) */
/* position angle of major axis, */
/* 	data aa, bb, cc       	/ 0.6, 0.25, 0.3 / 	! GUESS */
/* 	data xlsb, ylsb, zlsb 	/ -0.7, 9.0, 0. / 	! GUESS */
/* 	data theta		/ 150. / 		! GUESS */
/* 	data netrough  /0.01/		! GUESS */
/* 	data Ftrough   / 1./		! GUESS */
    aa = nelismparms_1.alsb;
    bb = nelismparms_1.blsb;
    cc = nelismparms_1.clsb;
    theta = nelismparms_1.thetalsb;
    netrough = nelismparms_1.nelsb0;
    ftrough = nelismparms_1.flsb;
    if (first) {
	s = sin(theta / 57.29577951f);
	c__ = cos(theta / 57.29577951f);
/* Computing 2nd power */
	r__1 = c__ / aa;
/* Computing 2nd power */
	r__2 = s / bb;
	ap = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = s / aa;
/* Computing 2nd power */
	r__2 = c__ / bb;
	bp = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = cc;
	cp = 1.f / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * 2.f * s * (1.f / (r__1 * r__1) - 1.f / (r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = 0.f;
    *wlsb = 0;
    *flsbr = 0.f;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xlsb;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.ylsb;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zlsb;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xlsb) * (*y - nelismparms_1.ylsb) * dp;
    if (q <= 1.f) {
/* inside */
	ret_val = netrough;
	*flsbr = ftrough;
	*wlsb = 1;
    }
    return ret_val;
} /* nelsb_ */

doublereal nelhb_(real *x, real *y, real *z__, real *flhbr, integer *wlhb)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static real netrough, c__, q, s, aa, bb, cc, ap, bp, cp, dp, theta, 
	    ftrough;

/* Local Hot Bubble */
/* input: */
/* 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00 */
/* output: */
/* 	neLHB = electron density in local hot bubble that */
/* 	        is modeled as an ellisoidal trough. */
/* 	FLHB = fluctuation parameter */
/* 	wLHB  = weight of LBH component used to combine */
/* 		with other components of electron density. */
/* 		wLBH =  1  at and inside the annular ridge */
/* 		     <  1  outside the annular ridge */
/* 	             -> 0  far outside the annular ridge */
/* 	e.g. total electron density would be evaluated as */
/*            ne = (1-wLHB)*ne_other + neLHB */
/* scales of ellipsoidal ridge */
/* ne of annulus, trough */
/* 	real xlhb, ylhb, zlhb		! center of ellipsoid */
/* fluctuation parameters */
/*    measured from x axis */
/*    (x axis points toward l=90) */
/* position angle of major axis, */
/* 	data aa, bb, cc       	/ 0.15, 0.08, 0.2 / */
/* 	data xlhb, ylhb, zlhb 	/ 0., 8.5, 0. / */
/* 	data theta		/ 135. / */
/* 	data netrough  /0.005/ */
/* 	data Ftrough   / 1./ */
    aa = nelismparms_1.alhb;
    bb = nelismparms_1.blhb;
    cc = nelismparms_1.clhb;
    theta = nelismparms_1.thetalhb;
    netrough = nelismparms_1.nelhb0;
    ftrough = nelismparms_1.flhb;
    if (first) {
	s = sin(theta / 57.29577951f);
	c__ = cos(theta / 57.29577951f);
/* Computing 2nd power */
	r__1 = c__ / aa;
/* Computing 2nd power */
	r__2 = s / bb;
	ap = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = s / aa;
/* Computing 2nd power */
	r__2 = c__ / bb;
	bp = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
	r__1 = cc;
	cp = 1.f / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * 2.f * s * (1.f / (r__1 * r__1) - 1.f / (r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = 0.f;
    *wlhb = 0;
    *flhbr = 0.f;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xlhb;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.ylhb;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zlhb;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xlhb) * (*y - nelismparms_1.ylhb) * dp;
    if (q <= 1.f) {
/* inside */
	ret_val = netrough;
	*flhbr = ftrough;
	*wlhb = 1;
    }
    return ret_val;
} /* nelhb_ */

doublereal nelhb2_(real *x, real *y, real *z__, real *flhbr, integer *wlhb)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double tan(doublereal);

    /* Local variables */
    static real netrough, aa, bb, cc, qz, qxy, theta, yaxis, ftrough, yzslope;

/* LHB modeled as a cylinder */
/* the cylinder slants in the y direction vs. z as described by parameter yzslope */
/* the cylinder cross-sectional size in the 'a' direction (major axis) */
/*       varies with z, tending to zero at its smallest z point. */
/* Local Hot Bubble */
/* input: */
/* 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00 */
/* output: */
/* 	neLHB2 = electron density in local hot bubble that */
/* 	        is modeled as an ellisoidal trough. */
/* 	FLHB = fluctuation parameter */
/* 	wLHB  = weight of LBH component used to combine */
/* 		with other components of electron density. */
/* 		wLHB =  1  at and inside the annular ridge */
/* 		     <  1  outside the annular ridge */
/* 	             -> 0  far outside the annular ridge */
/* 	e.g. total electron density would be evaluated as */
/*            ne = (1-wLHB)*ne_other + neLHB2 */
/* scales of ellipsoidal ridge */
/* ne of annulus, trough */
/* 	real xlhb, ylhb, zlhb		! center of ellipsoid */
/* fluctuation parameters */
/*    measured from z axis */
/* slant angle in yz plane of cylinder */
    aa = nelismparms_1.alhb;
    bb = nelismparms_1.blhb;
    cc = nelismparms_1.clhb;
    theta = nelismparms_1.thetalhb;
    netrough = nelismparms_1.nelhb0;
    ftrough = nelismparms_1.flhb;
    if (first) {
	yzslope = tan(theta / 57.29577951f);
	first = FALSE_;
    }
    ret_val = 0.f;
    *wlhb = 0;
    *flhbr = 0.f;
    yaxis = nelismparms_1.ylhb + yzslope * *z__;
/* cylinder has cross sectional area = constant for z>0 */
/* area -> 0 for z<0 by letting aa->0 linearly for z<0: */
/* (0.001 = 1 pc is to avoid divide by zero) */
    if (*z__ <= 0.f && *z__ >= nelismparms_1.zlhb - nelismparms_1.clhb) {
	aa = (nelismparms_1.alhb - .001f) * (1.f - 1.f / (nelismparms_1.zlhb 
		- nelismparms_1.clhb) * *z__) + .001f;
    } else {
	aa = nelismparms_1.alhb;
    }
/* 	write(99, *) x, y, z, aa, bb, cc */
/* Computing 2nd power */
    r__1 = (*x - nelismparms_1.xlhb) / aa;
/* Computing 2nd power */
    r__2 = (*y - yaxis) / bb;
    qxy = r__1 * r__1 + r__2 * r__2;
    qz = (r__1 = *z__ - nelismparms_1.zlhb, dabs(r__1)) / cc;
    if (qxy <= 1.f && qz <= 1.f) {
/* inside */
	ret_val = netrough;
	*flhbr = ftrough;
	*wlhb = 1;
    }
    return ret_val;
} /* nelhb2_ */

doublereal neloopi_(real *x, real *y, real *z__, real *floopi, integer *
	wloopi)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real r__, a1, a2;

/* component is a spheroid truncated for z<0. */
/* Loop I */
/* input: */
/* 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00 */
/* output: */
/* 	neLOOPI = electron density in LOOP I that */
/* 	        is modeled as an ellisoidal trough */
/* 		with an enhanced shell */
/* 	FLOOPI = fluctuation parameter */
/* 	wLOOPI  = weight of LOOP I component used to combine */
/* 		with other components of electron density. */
/* 		wLOOPI =  1  at and inside the annular ridge */
/* 		       <  1  outside the annular ridge */
    if (first) {
	a1 = nelismparms_1.rlpi;
	a2 = nelismparms_1.rlpi + nelismparms_1.drlpi;
	first = FALSE_;
    }
    if (*z__ < 0.f) {
	ret_val = 0.f;
	*floopi = 0.f;
	*wloopi = 0;
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xlpi;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.ylpi;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zlpi;
    r__ = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    if (r__ > a2) {
/* outside Loop I */
	ret_val = 0.f;
	*floopi = 0.f;
	*wloopi = 0;
    } else if (r__ <= a1) {
/* inside volume */
	ret_val = nelismparms_1.nelpi;
	*floopi = nelismparms_1.flpi;
	*wloopi = 1;
/*           write(99,*) x,y,z, r, neLOOPI, ' inside volume' */
    } else {
/* inside boundary shell */
	ret_val = nelismparms_1.dnelpi;
	*floopi = nelismparms_1.dflpi;
	*wloopi = 1;
/*           write(99,*) x,y,z,r, neLOOPI, ' inside shell' */
    }
    return ret_val;
} /* neloopi_ */

/* Subroutine */ int nevoidn_(real *x, real *y, real *z__, real *nevn, real *
	fvn, integer *hitvoid, integer *wvoid)
{
    /* Initialized data */

    static logical first = TRUE_;
    static integer luvoid = 11;

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), e_rsle(void), do_lio(integer *,
	     integer *, char *, ftnlen);
    double sin(doublereal), cos(doublereal);
    integer f_clos(cllist *);
    double exp(doublereal);

    /* Local variables */
    static integer voidflag, j;
    static real q, c1[2000], c2[2000], s1[2000], s2[2000], bv[2000], dv[2000],
	     fv[2000], dx, dy, dz, lv[2000], xv[2000], yv[2000], zv[2000], 
	    th1, th2, cbc, cc12[2000], clc, aav[2000], cs21[2000], bbv[2000], 
	    cs12[2000], ccv[2000], sbc, slc, nev[2000], ss12[2000];
    static integer edge[2000];
    static real thvy[2000], thvz[2000], rgalc;

    /* Fortran I/O blocks */
    static cilist io___291 = { 0, 0, 0, 0, 0 };
    static cilist io___292 = { 0, 0, 1, 0, 0 };


/* returns electron density nevN and fluctuation parameter FvN */
/* at position designated by l,b,d,x,y,z c for a set of */
/* voids with parameters read in from file  nevoidN.dat */
/* input: */
/* 	x,y,z	coordinates	(kpc)  (as in TC93) */

/* output: */
/* 	nevN	electron density in void at (x,y,z) */
/* 	FvN	fluctuation parameter */
/* 	hitvoid =   0:   no void hit */
/* 		  j>0:   j-th void hit */
/* 	wvoid = 0,1:	 void weight */
/* 	character*12 losname(nvoidsmax) */
/* parameters: */
/* 	lv	= galactic longitude of void center */
/* 	bv	= galactic latitude of void center */
/* 	dv	= distance from Sun of void center */
/* 	(xv,yv,zv) = void center location (calculated) */
/*       nev	= internal peak electron density */
/*       Fv      = void fluctuation parameter */
/* 	aav	= void major axis at 1/e */
/* 	bbv	= void minor axis at 1/e */
/* 	ccv	= void minor axis at 1/e */
/* 	thvy	= rotation axis of void about y axis */
/* 	thvz	= rotation axis of void about z axis */
/* 	edge    = 0 => use exponential rolloff out to 5rc */
/*                 1 => uniform and truncated at 1/e */
/* first time through, calculate xc, yc, zc */
    if (first) {
/* read void parameters */
	j = 1;
/* 	  write(6,*) 'reading nevoidN.dat.clean' */
	o__1.oerr = 0;
	o__1.ounit = luvoid;
	o__1.ofnmlen = 18;
	o__1.ofnm = "nevoidN.NE2001.dat";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	io___291.ciunit = luvoid;
	s_rsle(&io___291);
	e_rsle();
/* label line */
L5:
	io___292.ciunit = luvoid;
	i__1 = s_rsle(&io___292);
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&voidflag, (ftnlen)sizeof(integer)
		);
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&lv[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bv[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dv[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&nev[j - 1], (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&fv[j - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&aav[j - 1], (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bbv[j - 1], (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&ccv[j - 1], (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&thvy[j - 1], (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&thvz[j - 1], (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&edge[j - 1], (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L99;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L99;
	}
/* deg, deg, kpc */
/* cm^{-3}, dimensionless */
/* kpc,kpc,kpc,deg,d */
/* 0 or 1 */
	if (voidflag == 0) {
	    slc = sin(lv[j - 1] / 57.29577951f);
	    clc = cos(lv[j - 1] / 57.29577951f);
	    sbc = sin(bv[j - 1] / 57.29577951f);
	    cbc = cos(bv[j - 1] / 57.29577951f);
	    rgalc = dv[j - 1] * cbc;
	    xv[j - 1] = rgalc * slc;
	    yv[j - 1] = 8.5f - rgalc * clc;
	    zv[j - 1] = dv[j - 1] * sbc;
	    th1 = thvy[j - 1];
	    th2 = thvz[j - 1];
	    s1[j - 1] = sin(th1 / 57.29577951f);
	    c1[j - 1] = cos(th1 / 57.29577951f);
	    s2[j - 1] = sin(th2 / 57.29577951f);
	    c2[j - 1] = cos(th2 / 57.29577951f);
	    cc12[j - 1] = c1[j - 1] * c2[j - 1];
	    ss12[j - 1] = s1[j - 1] * s2[j - 1];
	    cs21[j - 1] = c2[j - 1] * s1[j - 1];
	    cs12[j - 1] = c1[j - 1] * s2[j - 1];
/* 	  write(6,"(a12,1x,13(f7.3,1x))") */
/*    .           losname(j),lv(j),bv(j),dv(j), */
/*    .           nev(j),Fv(j),xv(j),yv(j),zv(j), */
/*    .           aav(j), bbv(j), ccv(j), */
/*    .           th1, th2 */
	    ++j;
	}
	goto L5;
L99:
	first = FALSE_;
	voids_1.nvoids = j - 1;
	cl__1.cerr = 0;
	cl__1.cunit = luvoid;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    *nevn = 0.f;
    *fvn = 0.f;
    *hitvoid = 0;
    *wvoid = 0;
/* note rotation matrix in the 'q = ' statement below */
/* corresponds to \Lambda_z\Lambda_y */
/* where \Lambda_y = rotation around y axis */
/*       \Lambda_z = rotation around z axis */
/* defined as */
/* \Lambda_y =  c1  0  s1 */
/*               0  1   0 */
/*             -s1  0  c1 */
/* \Lambda_z =  c2 s2   0 */
/*             -s2 c2   0 */
/*               0  0   1 */
/* => */
/* \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2 */
/*                      -s2*c1   c2  -s1*s2 */
/*                         -s1    0      c1 */
/* so the rotation is around the y axis first, then the z axis */
    i__1 = voids_1.nvoids;
    for (j = 1; j <= i__1; ++j) {
	dx = *x - xv[j - 1];
	dy = *y - yv[j - 1];
	dz = *z__ - zv[j - 1];
/* Computing 2nd power */
	r__1 = cc12[j - 1] * dx + s2[j - 1] * dy + cs21[j - 1] * dz;
/* Computing 2nd power */
	r__2 = aav[j - 1];
/* Computing 2nd power */
	r__3 = -cs12[j - 1] * dx + c2[j - 1] * dy - ss12[j - 1] * dz;
/* Computing 2nd power */
	r__4 = bbv[j - 1];
/* Computing 2nd power */
	r__5 = -s1[j - 1] * dx + c1[j - 1] * dz;
/* Computing 2nd power */
	r__6 = ccv[j - 1];
	q = r__1 * r__1 / (r__2 * r__2) + r__3 * r__3 / (r__4 * r__4) + r__5 *
		 r__5 / (r__6 * r__6);
	if (edge[j - 1] == 0 && q < 3.f) {
	    *nevn = nev[j - 1] * exp(-q);
	    *fvn = fv[j - 1];
	    *hitvoid = j;
	    voids_1.hitvoidflag[j - 1] = 1;
	}
	if (edge[j - 1] == 1 && q <= 1.f) {
	    *nevn = nev[j - 1];
	    *fvn = fv[j - 1];
	    *hitvoid = j;
	    voids_1.hitvoidflag[j - 1] = 1;
	}
    }
    if (*hitvoid != 0) {
	*wvoid = 1;
    }
    return 0;
} /* nevoidn_ */

/* this version from 18 March 1998 uses revised coefficients */
/* that are consistent with Cordes \& Rickett (1998, ApJ, submitted) */
/* modifications: */
/* 	28 March 2001: added FUNCTION TRANSITION_FREQUENCY */
doublereal tauiss_(real *d__, real *sm, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* calculates the pulse broadening time in ms */
/* from distance, scattering measure, and radio frequency */

/* input:      d = pulsar distance       (kpc) */
/*            sm = scattering measure    (kpc m^{-20/3}) */
/*            nu = radio frequency       (GHz) */
/* output: tauss = pulse broadening time (ms) */

    d__1 = (doublereal) (*sm / 292.f);
    d__2 = (doublereal) (*nu);
    ret_val = pow_dd(&d__1, &c_b249) * 1e3f * *d__ * pow_dd(&d__2, &c_b250);
    return ret_val;
} /* tauiss_ */



doublereal scintbw_(real *d__, real *sm, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real tauiss;


/* calculates the scintillation bandwidth in kHz */
/* from distance, scattering measure, and radio frequency */

/* input:        d = pulsar distance       (kpc) */
/*              sm = scattering measure    (kpc m^{-20/3}) */
/*              nu = radio frequency       (GHz) */
/* output: scintbw = scintillation bandwidth (kHz) */

/* for uniform, Kolmogorov medium */
    d__1 = (doublereal) (*sm / 292.f);
    d__2 = (doublereal) (*nu);
    tauiss = pow_dd(&d__1, &c_b249) * 1e3f * *d__ * pow_dd(&d__2, &c_b250);
/* ms */
    ret_val = 1.16f / (tauiss * 6.2831799999999998f);
/* kHz */
    return ret_val;
} /* scintbw_ */

doublereal scintime_(real *sm, real *nu, real *vperp)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* calculates the scintillation speed for given distance, galactic */
/* longitude and latitude, frequency, and transverse velocity */

/* input:   sm = scattering measure	(kpc m^{-20/3}) */
/*          nu = radio frequency 	(GHz) */
/*       vperp = psr transverse speed  	(km/s) */

/* output: scintime = scintillation time (sec) */

/* usage: should be called with sm = smtau for appropriate */
/*        line of sight weighting */
/* reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123. */

/* nb: formerly, the coeff. in the following line was 2.3 from */
/*     Cordes & Lazio (1991) */
    d__1 = (doublereal) (*nu);
    d__2 = (doublereal) (*sm);
    ret_val = pow_dd(&d__1, &c_b249) * 3.3f * pow_dd(&d__2, &c_b256) * (100.f 
	    / *vperp);
    return ret_val;
} /* scintime_ */

doublereal specbroad_(real *sm, real *nu, real *vperp)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* calculates the bandwdith of spectral broadening */
/* for given scattering measure, , frequency, and transverse velocity */
/* input:   sm = scattering measure	(kpc m^{-20/3}) */
/*          nu = radio frequency 	(GHz) */
/*       vperp = psr transverse speed  	(km/s) */
/* output: specbroad = spectral broadening bandwidth (Hz) */

/* usage: should be called with sm = smtau for appropriate */
/*        line of sight weighting */
/* reference: eqn (47) of Cordes & Lazio 1991, ApJ, 376, 123. */

/* nb: the coeff. in the following line is 0.14 Hz  from Cordes & Lazio (1991) */
/* it is changed to 0.097 to conform with FUNCTION SCINTIME and */
/* a new calculation consistent with Cordes & Rickett (1998) */
    d__1 = (doublereal) (*nu);
    d__2 = (doublereal) (*sm);
    ret_val = pow_dd(&d__1, &c_b258) * .097f * pow_dd(&d__2, &c_b259) * (*
	    vperp / 100.f);
/* Hz */
    return ret_val;
} /* specbroad_ */

doublereal theta_xgal__(real *sm, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* calculates angular broadening for an extragalactic */
/* source of plane waves */

/* sm = scattering measure */
/* nu = radio frequency */
/* theta_xgal = angular broadening FWHM (mas) */

    d__1 = (doublereal) (*sm);
    d__2 = (doublereal) (*nu);
    ret_val = pow_dd(&d__1, &c_b259) * 128.f * pow_dd(&d__2, &c_b262);
    return ret_val;
} /* theta_xgal__ */


doublereal theta_gal__(real *sm, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* calculates angular broadening for a galactic */
/* source of spherical waves */

/* sm = scattering measure */
/* nu = radio frequency */
/* theta_gal = angular broadening FWHM (mas) */

    d__1 = (doublereal) (*sm);
    d__2 = (doublereal) (*nu);
    ret_val = pow_dd(&d__1, &c_b259) * 71.f * pow_dd(&d__2, &c_b262);
    return ret_val;
} /* theta_gal__ */


doublereal em_(real *sm)
{
    /* Initialized data */

    static real router = 1.f;
    static real pc = 3.086e18f;
    static real alpha = 3.6666667f;
    static real pi = 3.14159f;

    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


/* units of sm are kpc m^{-20/3} */
/* units of em are pc cm^{-6} */

/* calculates the emission measure from the scattering measure */
/* using an assumed outer scale and spectral index of the */
/* wavenumber spectrum. */

/* for a wavenumber spectrum P_n(q) = q^{-alpha} from q_0 to q_1 */
/* the mean square electron density is */

/* <n_e^2> =~  4pi*[C_n^2 / (alpha - 3) ] * q_0^{3 - alpha) */

/* ( an approximate form that assumes (q_0 / q_1)^{3-alpha} >> 1. */

/* Jim Cordes 18 Dec 1989 */

/* outer scale = 1 pc */

    d__1 = (doublereal) (router * pc / 6.2831799999999998f);
    d__2 = (doublereal) (alpha - 3.f);
    ret_val = *sm * (4.f * pi * 1e3f / (alpha - 3.f)) * pow_dd(&d__1, &d__2) *
	     pow_dd(&c_b267, &c_b268);

    return ret_val;
} /* em_ */

doublereal theta_iso__(real *smiso, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    double r_lg10(real *), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real theta_log_radian__, theta_log_microarcsec__;

/* smiso in (kpc m^{-20/3}) x kpc^{5/3} */
/*    nu in GHz */
/* returns the isoplanatic angle in microarcsec */
/* 12 October 1998 */
/* JMC */
/* \theta_{iso} = \delta r_s / d */
/*              = \left [ */
/* 	         (\lambda r_e)^2 f_{\alpha} SM_{iso} */
/* 		 \right ]^{1/\alpha} */
/* where \alpha = 5/3 for Kolmogorov case. */
/* NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq */
/*    so SM_{iso} does not have the units of scattering */
/*    measure, but rather units of SM x Length^{\alpha} */

/* f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)] */
/* for \alpha = 5/3, f_{\alpha}= 88.3 */

/*     real r_e */
/*     parameter(r_e = 2.82e-13)			!cm */
/*     real kpc */
/*     parameter(kpc = 3.086e21)			!cm */
/*     real falpha */
/*     parameter(falpha=88.3) */
    theta_log_radian__ = r_lg10(nu) * 1.2f + 13.287f - 1.1676f - r_lg10(smiso)
	     * .6f - 34.383f + 8.f;
/* 0.6*log10(30cm*r_e) */
/* 0.6*log10(f_alpha) */
/* 1.6 * alog10(kpc) */
/* -(20/3)*log(100) */
    theta_log_microarcsec__ = theta_log_radian__ + 11.314425f;
/* 11.314425=alog10(microarsec/r */
    d__1 = (doublereal) theta_log_microarcsec__;
    ret_val = pow_dd(&c_b270, &d__1);
    return ret_val;
} /* theta_iso__ */

doublereal theta_iso_test__(real *smiso, real *nu)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    double r_lg10(real *), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real theta_log_radian__, theta_log_microarcsec__;

/* smiso in (kpc m^{-20/3}) x kpc^{5/3} */
/*    nu in GHz */
/* returns the isoplanatic angle in microarcsec */
/* 12 October 1998 */
/* JMC */
/* \theta_{iso} = \delta r_s / d */
/*              = \left [ */
/* 	         (\lambda r_e)^2 f_{\alpha} SM_{iso} */
/* 		 \right ]^{1/\alpha} */
/* where \alpha = 5/3 for Kolmogorov case. */
/* NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq */
/*    so SM_{iso} does not have the units of scattering */
/*    measure, but rather units of SM x Length^{\alpha} */

/* f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)] */
/* for \alpha = 5/3, f_{\alpha}= 88.3 */

/*     real r_e */
/*     parameter(r_e = 2.82e-13)			!cm */
/*     real kpc */
/*     parameter(kpc = 3.086e21)			!cm */
/*     real falpha */
/*     parameter(falpha=88.3) */
    theta_log_radian__ = r_lg10(nu) * 1.2f + 13.287f - 1.1676f - r_lg10(smiso)
	     * .6f - 34.383f + 8.f;
/* 0.6*log10(30cm*r_e) */
/* 0.6*log10(f_alpha) */
/* 1.6 * alog10(kpc) */
/* -(20/3)*log(100) */
    theta_log_microarcsec__ = theta_log_radian__ + 11.314425f;
/* 11.314425=alog10(microarsec/r */
    d__1 = (doublereal) theta_log_microarcsec__;
    ret_val = pow_dd(&c_b270, &d__1);
/*     write(6,*) 'smiso, nu = ', smiso, nu */
/*     write(6,*) 'theta_log_radian = ', theta_log_radian */
/*     write(6,*) 'theta_log_microarcsec = ', theta_log_microarcsec */
/*     write(6,*) 'theta_iso = ', theta_iso_test */
    return ret_val;
} /* theta_iso_test__ */

doublereal transition_frequency__(real *sm, real *smtau, real *smtheta, real *
	dintegrate)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real deff;

/* returns the transition frequency between weak and strong scattering */
/* 28 March 2001 */
/* JMC */
/* input: */
/* (all sm values in (kpc m^{-20/3})) */
/* 	    sm = int[\cnsq] */
/*        smtau = int[(s/D)(1-s/D) \cnsq] */
/*      smtheta = int[ (1-s/D) \cnsq] */
/*   dintegrate = distance used to integrate \cnsq (kpc) */
/* output: */
/*   transition_frequency = GHz given by */
/* 	\nu_t  = 318 GHz \xi^{10/17}  SM^{6/17} D_{eff}^{5/17} */
/*   where */
/*        D_{eff} = effective path length through medium */
/*        D_{eff} = \int_0^dintegrate ds s \cnsq / \int_0^dintegrate ds  \cnsq */

/*   Note we can calculate D_{eff} using */
/*        D_{eff} = dintegrate * (sm - smtau/6 - smtau/3) / sm */

/* (2.*pi)^{-1/2} = fresnel scale definitio */
/* GHz; see NE2001 paper */
    deff = *dintegrate * (*sm - *smtau / 6.f - *smtheta / 3.f) / *sm;
    d__1 = (doublereal) (*sm);
    d__2 = (doublereal) deff;
    ret_val = pow_dd(&c_b274, &c_b275) * 318.f * pow_dd(&d__1, &c_b276) * 
	    pow_dd(&d__2, &c_b277);
    return ret_val;
} /* transition_frequency__ */

