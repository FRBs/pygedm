/* neclumpN.f -- translated by f2c (version 20100827).
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
    integer nclumps, hitclumpflag[2000];
} clumps_;

#define clumps_1 clumps_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__4 = 4;

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
    integer f_open(olist *), s_rsle(cilist *), e_rsle(), do_lio(integer *, 
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
    static cilist io___4 = { 0, 0, 0, 0, 0 };
    static cilist io___5 = { 0, 0, 1, 0, 0 };


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
	io___4.ciunit = luclump;
	s_rsle(&io___4);
	e_rsle();
/* label line */
L5:
	io___5.ciunit = luclump;
	i__1 = s_rsle(&io___5);
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
	    slc = sin(lc[j - 1] / (float)57.29577951);
	    clc = cos(lc[j - 1] / (float)57.29577951);
	    sbc = sin(bc[j - 1] / (float)57.29577951);
	    cbc = cos(bc[j - 1] / (float)57.29577951);
	    rgalc = dc[j - 1] * cbc;
	    xc[j - 1] = rgalc * slc;
	    yc[j - 1] = (float)8.5 - rgalc * clc;
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
    *necn = (float)0.;
    *hitclump = 0;
    *fcn = (float)0.;
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
	if (edge[j - 1] == 0 && arg < (float)5.) {
	    *necn += nec[j - 1] * exp(-arg);
	    *fcn = fc[j - 1];
	    *hitclump = j;
	    clumps_1.hitclumpflag[j - 1] = 1;
	}
	if (edge[j - 1] == 1 && arg <= (float)1.) {
/*    	    necN = necN + nec(j) * exp(-arg) */
	    *necn += nec[j - 1];
	    *fcn = fc[j - 1];
	    *hitclump = j;
	    clumps_1.hitclumpflag[j - 1] = 1;
	}
    }
    return 0;
} /* neclumpn_ */

#ifdef __cplusplus
	}
#endif
