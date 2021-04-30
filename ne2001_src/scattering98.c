/* scattering98.f -- translated by f2c (version 20100827).
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

static doublereal c_b2 = 1.2;
static doublereal c_b3 = -4.4;
static doublereal c_b9 = -.6;
static doublereal c_b11 = -1.2;
static doublereal c_b12 = .6;
static doublereal c_b15 = -2.2;
static doublereal c_b20 = .01;
static doublereal c_b21 = 6.666666666666667;
static doublereal c_b23 = 10.;
static doublereal c_b27 = .3989;
static doublereal c_b28 = .58823529411764708;
static doublereal c_b29 = .35294117647058826;
static doublereal c_b30 = .29411764705882354;

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
    ret_val = pow_dd(&d__1, &c_b2) * 1e3f * *d__ * pow_dd(&d__2, &c_b3);
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
    tauiss = pow_dd(&d__1, &c_b2) * 1e3f * *d__ * pow_dd(&d__2, &c_b3);
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
    ret_val = pow_dd(&d__1, &c_b2) * 3.3f * pow_dd(&d__2, &c_b9) * (100.f / *
	    vperp);
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
    ret_val = pow_dd(&d__1, &c_b11) * .097f * pow_dd(&d__2, &c_b12) * (*vperp 
	    / 100.f);
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
    ret_val = pow_dd(&d__1, &c_b12) * 128.f * pow_dd(&d__2, &c_b15);
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
    ret_val = pow_dd(&d__1, &c_b12) * 71.f * pow_dd(&d__2, &c_b15);
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
	     pow_dd(&c_b20, &c_b21);

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
    ret_val = pow_dd(&c_b23, &d__1);
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
    ret_val = pow_dd(&c_b23, &d__1);
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
    ret_val = pow_dd(&c_b27, &c_b28) * 318.f * pow_dd(&d__1, &c_b29) * pow_dd(
	    &d__2, &c_b30);
    return ret_val;
} /* transition_frequency__ */

