/* neLISM.NE2001.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct {
    real aldr, bldr, cldr, xldr, yldr, zldr, thetaldr, neldr0, fldr, alsb, 
	    blsb, clsb, xlsb, ylsb, zlsb, thetalsb, nelsb0, flsb, alhb, blhb, 
	    clhb, xlhb, ylhb, zlhb, thetalhb, nelhb0, flhb, xlpi, ylpi, zlpi, 
	    rlpi, drlpi, nelpi, flpi, dnelpi, dflpi;
} nelismparms_;

#define nelismparms_1 nelismparms_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

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
    integer f_open(olist *), s_rsle(cilist *), e_rsle(), do_lio(integer *, 
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
    static cilist io___2 = { 0, 11, 0, 0, 0 };
    static cilist io___3 = { 0, 11, 0, 0, 0 };
    static cilist io___4 = { 0, 11, 0, 0, 0 };
    static cilist io___5 = { 0, 11, 0, 0, 0 };
    static cilist io___6 = { 0, 11, 0, 0, 0 };
    static cilist io___7 = { 0, 11, 0, 0, 0 };
    static cilist io___8 = { 0, 11, 0, 0, 0 };
    static cilist io___9 = { 0, 11, 0, 0, 0 };
    static cilist io___10 = { 0, 11, 0, 0, 0 };
    static cilist io___11 = { 0, 11, 0, 0, 0 };
    static cilist io___12 = { 0, 11, 0, 0, 0 };
    static cilist io___13 = { 0, 11, 0, 0, 0 };
    static cilist io___14 = { 0, 11, 0, 0, 0 };


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
	s_rsle(&io___2);
	e_rsle();
	s_rsle(&io___3);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.aldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.bldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.cldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___4);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.yldr, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___5);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetaldr, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.neldr0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.fldr, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___6);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.alsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.blsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.clsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___7);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylsb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___8);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetalsb, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.nelsb0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.flsb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___9);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.alhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.blhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.clhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___10);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylhb, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___11);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.thetalhb, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.nelhb0, (ftnlen)sizeof(
		real));
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.flhb, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___12);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.xlpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.ylpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.zlpi, (ftnlen)sizeof(real)
		);
	e_rsle();
	s_rsle(&io___13);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.rlpi, (ftnlen)sizeof(real)
		);
	do_lio(&c__4, &c__1, (char *)&nelismparms_1.drlpi, (ftnlen)sizeof(
		real));
	e_rsle();
	s_rsle(&io___14);
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
	s = sin(theta / (float)57.29577951);
	c__ = cos(theta / (float)57.29577951);
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
	cp = (float)1. / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * (float)2. * s * ((float)1. / (r__1 * r__1) - (float)1. / (
		r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = (float)0.;
    *wldrq1 = 0;
    *fldrq1r = (float)0.;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xldr;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.yldr;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zldr;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xldr) * (*y - nelismparms_1.yldr) * dp;
    if (q <= (float)1.) {
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
	s = sin(theta / (float)57.29577951);
	c__ = cos(theta / (float)57.29577951);
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
	cp = (float)1. / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * (float)2. * s * ((float)1. / (r__1 * r__1) - (float)1. / (
		r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = (float)0.;
    *wlsb = 0;
    *flsbr = (float)0.;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xlsb;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.ylsb;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zlsb;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xlsb) * (*y - nelismparms_1.ylsb) * dp;
    if (q <= (float)1.) {
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
	s = sin(theta / (float)57.29577951);
	c__ = cos(theta / (float)57.29577951);
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
	cp = (float)1. / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = aa;
/* Computing 2nd power */
	r__2 = bb;
	dp = c__ * (float)2. * s * ((float)1. / (r__1 * r__1) - (float)1. / (
		r__2 * r__2));
	first = FALSE_;
/* 	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp */
    }
    ret_val = (float)0.;
    *wlhb = 0;
    *flhbr = (float)0.;
/* Computing 2nd power */
    r__1 = *x - nelismparms_1.xlhb;
/* Computing 2nd power */
    r__2 = *y - nelismparms_1.ylhb;
/* Computing 2nd power */
    r__3 = *z__ - nelismparms_1.zlhb;
    q = r__1 * r__1 * ap + r__2 * r__2 * bp + r__3 * r__3 * cp + (*x - 
	    nelismparms_1.xlhb) * (*y - nelismparms_1.ylhb) * dp;
    if (q <= (float)1.) {
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
	yzslope = tan(theta / (float)57.29577951);
	first = FALSE_;
    }
    ret_val = (float)0.;
    *wlhb = 0;
    *flhbr = (float)0.;
    yaxis = nelismparms_1.ylhb + yzslope * *z__;
/* cylinder has cross sectional area = constant for z>0 */
/* area -> 0 for z<0 by letting aa->0 linearly for z<0: */
/* (0.001 = 1 pc is to avoid divide by zero) */
    if (*z__ <= (float)0. && *z__ >= nelismparms_1.zlhb - nelismparms_1.clhb) 
	    {
	aa = (nelismparms_1.alhb - (float).001) * ((float)1. - (float)1. / (
		nelismparms_1.zlhb - nelismparms_1.clhb) * *z__) + (float)
		.001;
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
    if (qxy <= (float)1. && qz <= (float)1.) {
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
    if (*z__ < (float)0.) {
	ret_val = (float)0.;
	*floopi = (float)0.;
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
	ret_val = (float)0.;
	*floopi = (float)0.;
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

#ifdef __cplusplus
	}
#endif
