/* nevoidN.f -- translated by f2c (version 20100827).
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
    integer nvoids, hitvoidflag[2000];
} voids_;

#define voids_1 voids_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__4 = 4;

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
    integer f_open(olist *), s_rsle(cilist *), e_rsle(), do_lio(integer *, 
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
    static cilist io___4 = { 0, 0, 0, 0, 0 };
    static cilist io___5 = { 0, 0, 1, 0, 0 };


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
	io___4.ciunit = luvoid;
	s_rsle(&io___4);
	e_rsle();
/* label line */
L5:
	io___5.ciunit = luvoid;
	i__1 = s_rsle(&io___5);
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
	    slc = sin(lv[j - 1] / (float)57.29577951);
	    clc = cos(lv[j - 1] / (float)57.29577951);
	    sbc = sin(bv[j - 1] / (float)57.29577951);
	    cbc = cos(bv[j - 1] / (float)57.29577951);
	    rgalc = dv[j - 1] * cbc;
	    xv[j - 1] = rgalc * slc;
	    yv[j - 1] = (float)8.5 - rgalc * clc;
	    zv[j - 1] = dv[j - 1] * sbc;
	    th1 = thvy[j - 1];
	    th2 = thvz[j - 1];
	    s1[j - 1] = sin(th1 / (float)57.29577951);
	    c1[j - 1] = cos(th1 / (float)57.29577951);
	    s2[j - 1] = sin(th2 / (float)57.29577951);
	    c2[j - 1] = cos(th2 / (float)57.29577951);
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
    *nevn = (float)0.;
    *fvn = (float)0.;
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
	if (edge[j - 1] == 0 && q < (float)3.) {
	    *nevn = nev[j - 1] * exp(-q);
	    *fvn = fv[j - 1];
	    *hitvoid = j;
	    voids_1.hitvoidflag[j - 1] = 1;
	}
	if (edge[j - 1] == 1 && q <= (float)1.) {
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

#ifdef __cplusplus
	}
#endif
