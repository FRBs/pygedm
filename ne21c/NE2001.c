/* NE2001.f -- translated by f2c (version 20181026).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__9 = 9;

/* 	program NE2001 */
/*  calls dmdsm to get pulsar distance as a function of l, b, and DM, */
/*  or pulsar DM as a function of l, b, and distance. */
/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static real nu = 1.f;
    static real rad = 57.2957795f;
    static real vperp = 100.f;

    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsle(cilist *), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sin(doublereal);

    /* Local variables */
    extern doublereal scintime_(real *, real *, real *), theta_gal__(real *, 
	    real *);
    static real transfreq, dm;
    extern doublereal em_(real *);
    static real sm;
    extern doublereal theta_xgal__(real *, real *);
    static real dmd, tau, sbw, bdeg, ldeg;
    static integer ndir;
    static real emsm, dist;
    extern integer iargc_(void);
    static char inbuf[80];
    extern /* Subroutine */ int dmdsm_(real *, real *, integer *, real *, 
	    real *, char *, real *, real *, real *, real *, ftnlen);
    static integer nargs;
    static char limit[1];
    static real stime, smtau, smiso;
    extern doublereal transition_frequency__(real *, real *, real *, real *);
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen);
    extern doublereal tauiss_(real *, real *, real *);
    static real theta_g__, theta_x__, smtheta;
    extern doublereal scintbw_(real *, real *, real *);

    /* Fortran I/O blocks */
    static icilist io___6 = { 0, inbuf, 0, 0, 80, 1 };
    static icilist io___8 = { 0, inbuf, 0, 0, 80, 1 };
    static icilist io___10 = { 0, inbuf, 0, 0, 80, 1 };
    static icilist io___12 = { 0, inbuf, 0, 0, 80, 1 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___21 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___22 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___23 = { 0, 6, 0, "(i10, t20, a, t30, a, t55, a)", 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, "(a,f9.4, t20, a, t30, a, t55,a)", 0 };
    static cilist io___33 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___34 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___35 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___36 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55,a)", 0 };
    static cilist io___37 = { 0, 6, 0, "(f10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___45 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___46 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___47 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___48 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___49 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___50 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___51 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___52 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___53 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___54 = { 0, 6, 0, "(e10.4, t20, a, t30, a, t55, a)", 0 };
    static cilist io___55 = { 0, 6, 0, "(f10.2, t20, a, t30, a, t55, a)", 0 };


/* functions: */
    nargs = iargc_();
    if (nargs >= 1) {
	getarg_(&c__1, inbuf, (ftnlen)80);
	s_rsli(&io___6);
	do_lio(&c__4, &c__1, (char *)&ldeg, (ftnlen)sizeof(real));
	e_rsli();
	getarg_(&c__2, inbuf, (ftnlen)80);
	s_rsli(&io___8);
	do_lio(&c__4, &c__1, (char *)&bdeg, (ftnlen)sizeof(real));
	e_rsli();
	getarg_(&c__3, inbuf, (ftnlen)80);
	s_rsli(&io___10);
	do_lio(&c__4, &c__1, (char *)&dmd, (ftnlen)sizeof(real));
	e_rsli();
	getarg_(&c__4, inbuf, (ftnlen)80);
	s_rsli(&io___12);
	do_lio(&c__3, &c__1, (char *)&ndir, (ftnlen)sizeof(integer));
	e_rsli();
    } else {
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "Usage: NE2001 l b DM/D ndir", (ftnlen)27);
	e_wsle();
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, "       l (deg)", (ftnlen)14);
	e_wsle();
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, "       b (deg)", (ftnlen)14);
	e_wsle();
	s_wsle(&io___17);
	do_lio(&c__9, &c__1, "       DM/D (pc cm^{-3} or kpc)", (ftnlen)31);
	e_wsle();
	s_wsle(&io___18);
	do_lio(&c__9, &c__1, "       ndir = 1 (DM->D)   -1 (D->DM) ", (ftnlen)
		37);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "#NE2001 input: 4 parameters", (ftnlen)27);
    e_wsle();
    s_wsfe(&io___20);
    do_fio(&c__1, (char *)&ldeg, (ftnlen)sizeof(real));
    do_fio(&c__1, "l", (ftnlen)1);
    do_fio(&c__1, "(deg)", (ftnlen)5);
    do_fio(&c__1, "GalacticLongitude", (ftnlen)17);
    e_wsfe();
    s_wsfe(&io___21);
    do_fio(&c__1, (char *)&bdeg, (ftnlen)sizeof(real));
    do_fio(&c__1, "b", (ftnlen)1);
    do_fio(&c__1, "(deg)", (ftnlen)5);
    do_fio(&c__1, "GalacticLatitude", (ftnlen)16);
    e_wsfe();
    s_wsfe(&io___22);
    do_fio(&c__1, (char *)&dmd, (ftnlen)sizeof(real));
    do_fio(&c__1, "DM/D", (ftnlen)4);
    do_fio(&c__1, "(pc-cm^{-3}_or_kpc)", (ftnlen)19);
    do_fio(&c__1, "Input_DM_or_Distance", (ftnlen)20);
    e_wsfe();
    s_wsfe(&io___23);
    do_fio(&c__1, (char *)&ndir, (ftnlen)sizeof(integer));
    do_fio(&c__1, "ndir", (ftnlen)4);
    do_fio(&c__1, "1:DM->D;-1:D->DM", (ftnlen)16);
    do_fio(&c__1, "Which?(DM_or_D)", (ftnlen)15);
    e_wsfe();
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "#NE2001 output: 14 values", (ftnlen)25);
    e_wsle();
    if (ndir >= 0) {
	dm = dmd;
	r__1 = ldeg / rad;
	r__2 = bdeg / rad;
	dmdsm_(&r__1, &r__2, &ndir, &dm, &dist, limit, &sm, &smtau, &smtheta, 
		&smiso, (ftnlen)1);
	s_wsfe(&io___32);
	do_fio(&c__1, limit, (ftnlen)1);
	do_fio(&c__1, (char *)&dist, (ftnlen)sizeof(real));
	do_fio(&c__1, "DIST", (ftnlen)4);
	do_fio(&c__1, "(kpc)", (ftnlen)5);
	do_fio(&c__1, "ModelDistance", (ftnlen)13);
	e_wsfe();
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&dm, (ftnlen)sizeof(real));
	do_fio(&c__1, "DM", (ftnlen)2);
	do_fio(&c__1, "(pc-cm^{-3})", (ftnlen)12);
	do_fio(&c__1, "DispersionMeasure", (ftnlen)17);
	e_wsfe();
	s_wsfe(&io___34);
	r__2 = dm * (r__1 = sin(bdeg / rad), dabs(r__1));
	do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	do_fio(&c__1, "DMz", (ftnlen)3);
	do_fio(&c__1, "(pc-cm^{-3})", (ftnlen)12);
	do_fio(&c__1, "DM_Zcomponent", (ftnlen)13);
	e_wsfe();
    } else {
	dist = dmd;
	r__1 = ldeg / rad;
	r__2 = bdeg / rad;
	dmdsm_(&r__1, &r__2, &ndir, &dm, &dist, limit, &sm, &smtau, &smtheta, 
		&smiso, (ftnlen)1);
/* 	  write(6,*) 'dm,sm,smtau,smtheta = ', */
/*    .               dm,sm,smtau,smtheta */
/* 1020	  format(f8.2,3(1x,e8.3)) */
/* 	  write(6,*) 'dmz = ', dm*abs(sin(bdeg/rad)) */
	s_wsfe(&io___35);
	do_fio(&c__1, (char *)&dist, (ftnlen)sizeof(real));
	do_fio(&c__1, "DIST", (ftnlen)4);
	do_fio(&c__1, "(kpc)", (ftnlen)5);
	do_fio(&c__1, "Distance", (ftnlen)8);
	e_wsfe();
	s_wsfe(&io___36);
	do_fio(&c__1, (char *)&dm, (ftnlen)sizeof(real));
	do_fio(&c__1, "DM", (ftnlen)2);
	do_fio(&c__1, "(pc-cm^{-3})", (ftnlen)12);
	do_fio(&c__1, "ModelDM", (ftnlen)7);
	e_wsfe();
	s_wsfe(&io___37);
	r__2 = dm * (r__1 = sin(bdeg / rad), dabs(r__1));
	do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	do_fio(&c__1, "DMz", (ftnlen)3);
	do_fio(&c__1, "(pc-cm^{-3})", (ftnlen)12);
	do_fio(&c__1, "model", (ftnlen)5);
	e_wsfe();
    }
/* calculate scattering parameters */
    tau = tauiss_(&dist, &smtau, &nu);
    sbw = scintbw_(&dist, &smtau, &nu);
    stime = scintime_(&smtau, &nu, &vperp);
    theta_x__ = theta_xgal__(&sm, &nu);
    theta_g__ = theta_gal__(&smtheta, &nu);
    transfreq = transition_frequency__(&sm, &smtau, &smtheta, &dist);
    emsm = em_(&sm);
    s_wsfe(&io___45);
    do_fio(&c__1, (char *)&sm, (ftnlen)sizeof(real));
    do_fio(&c__1, "SM", (ftnlen)2);
    do_fio(&c__1, "(kpc-m^{-20/3})", (ftnlen)15);
    do_fio(&c__1, "ScatteringMeasure", (ftnlen)17);
    e_wsfe();
    s_wsfe(&io___46);
    do_fio(&c__1, (char *)&smtau, (ftnlen)sizeof(real));
    do_fio(&c__1, "SMtau", (ftnlen)5);
    do_fio(&c__1, "(kpc-m^{-20/3})", (ftnlen)15);
    do_fio(&c__1, "SM_PulseBroadening", (ftnlen)18);
    e_wsfe();
    s_wsfe(&io___47);
    do_fio(&c__1, (char *)&smtheta, (ftnlen)sizeof(real));
    do_fio(&c__1, "SMtheta", (ftnlen)7);
    do_fio(&c__1, "(kpc-m^{-20/3})", (ftnlen)15);
    do_fio(&c__1, "SM_GalAngularBroadening", (ftnlen)23);
    e_wsfe();
    s_wsfe(&io___48);
    do_fio(&c__1, (char *)&smiso, (ftnlen)sizeof(real));
    do_fio(&c__1, "SMiso", (ftnlen)5);
    do_fio(&c__1, "(kpc-m^{-20/3})", (ftnlen)15);
    do_fio(&c__1, "SM_IsoplanaticAngle", (ftnlen)19);
    e_wsfe();
    s_wsfe(&io___49);
    do_fio(&c__1, (char *)&emsm, (ftnlen)sizeof(real));
    do_fio(&c__1, "EM", (ftnlen)2);
    do_fio(&c__1, "(pc-cm^{-6})", (ftnlen)12);
    do_fio(&c__1, "EmissionMeasure_from_SM", (ftnlen)23);
    e_wsfe();
    s_wsfe(&io___50);
    do_fio(&c__1, (char *)&tau, (ftnlen)sizeof(real));
    do_fio(&c__1, "TAU", (ftnlen)3);
    do_fio(&c__1, "(ms)", (ftnlen)4);
    do_fio(&c__1, "PulseBroadening @1GHz", (ftnlen)21);
    e_wsfe();
    s_wsfe(&io___51);
    r__1 = sbw / 1e3f;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    do_fio(&c__1, "SBW", (ftnlen)3);
    do_fio(&c__1, "(MHz)", (ftnlen)5);
    do_fio(&c__1, "ScintBW @1GHz", (ftnlen)13);
    e_wsfe();
    s_wsfe(&io___52);
    do_fio(&c__1, (char *)&stime, (ftnlen)sizeof(real));
    do_fio(&c__1, "SCINTIME", (ftnlen)8);
    do_fio(&c__1, "(s)", (ftnlen)3);
    do_fio(&c__1, "ScintTime @1GHz @100 km/s", (ftnlen)25);
    e_wsfe();
    s_wsfe(&io___53);
    do_fio(&c__1, (char *)&theta_g__, (ftnlen)sizeof(real));
    do_fio(&c__1, "THETA_G", (ftnlen)7);
    do_fio(&c__1, "(mas)", (ftnlen)5);
    do_fio(&c__1, "AngBroadeningGal @1GHz", (ftnlen)22);
    e_wsfe();
    s_wsfe(&io___54);
    do_fio(&c__1, (char *)&theta_x__, (ftnlen)sizeof(real));
    do_fio(&c__1, "THETA_X", (ftnlen)7);
    do_fio(&c__1, "(mas)", (ftnlen)5);
    do_fio(&c__1, "AngBroadeningXgal @1GHz", (ftnlen)23);
    e_wsfe();
    s_wsfe(&io___55);
    do_fio(&c__1, (char *)&transfreq, (ftnlen)sizeof(real));
    do_fio(&c__1, "NU_T", (ftnlen)4);
    do_fio(&c__1, "(GHz)", (ftnlen)5);
    do_fio(&c__1, "TransitionFrequency", (ftnlen)19);
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

