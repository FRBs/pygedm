# mwprop v2.0 Jan 2026

r"""
Python versions of functions in scattering98.f that is
part of the NE2001 Fortran package

2019 December 31

Added:
    transition_frequency_from_obs
    sm_from_tau

2026 Feb 2: Moved tau_x calculation here
------
Notes in 1998-2001 Fortran version:
This version from 18 March 1998 uses revised coefficients
that are consistent with Cordes \& Rickett (1998, ApJ)

Note that scaling laws are explicitly for a Kolmogorov medium
in the strong but not superstrong regime
(as defined in Cordes and Lazio 1991)

Removed theta_iso_test (was redundant)

Modifications:
28 March 2001: added FUNCTION TRANSITION_FREQUENCY
"""

from __future__ import print_function

from math import pi

from pygedm.mwprop_core.get_constants_waveprop_astropy_version import *

import numpy as np

def fbeta(si=11./3.):
    """
    Calculates the f(\beta) function used in scattering calculations
    for a power-law wavenumber spectrum with spectral index si == beta.
    [nb. beta not used here to avoid confusion with the beta PDF]

    Also note that other 'fbeta' functions in the literature
    are the same as here but multiplied by 4*pi**2

    Cordes & Lazio 2003  NE2001 Paper II, Appendix Eq A5
    JMC 2020 Jan 1
    """
    from math import gamma

    beta2 = si / 2.
    beta1 = si - 1
    fbeta = (gamma(2-beta2) * gamma(beta2-1)) / (gamma(beta2)**2 * 2**beta1)
    return fbeta

def k_eta(si=11./3.):
    r"""
    Calculates the K_\eta(\beta) function used in scattering calculations
    for a power-law wavenumber spectrum with spectral index si == beta.
    [nb. beta not used here to avoid confusion with the beta PDF]

    From A-O book snippet (Drafted in 2019)
    'Estimating the Wavenumber Spectrum for Electron Density in the ISM'
    Eq 29.

    JMC 2020 Jan 1
    """
    from math import gamma

    beta2 = si / 2.
    beta4 = 4 - si

    keta = (2.*pi)**beta4 * gamma(3-beta2) / beta4

    return keta


def tauiss(d, sm, nu):

    """
    calculates the pulse broadening time in ms
    from distance, scattering measure, and radio frequency

    input:      d = pulsar distance       (kpc)
               sm = scattering measure    (kpm^{-20/3})
               nu = radio frequency       (GHz)
    output: tauss = pulse broadening time (ms)
    """

    tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
    return tauiss

def scintbw(d, sm, nu, C1 = 1.16):

    """
    calculates the scintillation bandwidth in kHz
    from distance, scattering measure, and radio frequency

    input:        d = pulsar distance       (kpc)
                 sm = scattering measure    (kpc m^{-20/3})
                 nu = radio frequency       (GHz)
                   C1 = 1.16 = dimensionless constant in sbw-tau relation
    output: scintbw = scintillation bandwidth (MHz) (cf scattering98: kHz)
    """

    tau = tauiss(d, sm, nu)                       # ms
    scintbw = 1.e-3 * C1 / (2. * pi * tau)
    return scintbw

def scintime(sm, nu, vperp=100):
    """

    calculates the scintillation speed for given distance, galactic
    longitude and latitude, frequency, and transverse velocity

    input:   sm = scattering measure    (kpc m^{-20/3})
             nu = radio frequency   (GHz)
          vperp = psr transverse speed      (km/s)
                  (default 100)

    output: scintime = scintillation time (sec)

    usage: should be called with sm = smtau for appropriate
           line of sight weighting
    reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123.
    """

    scintime = 3.3 * nu**1.2 * sm**(-0.6) * (100./vperp)
    return scintime


def specbroad(sm, nu, vperp):
    """

    calculates the bandwdith of spectral broadening
    for given scattering measure, , frequency, and transverse velocity

    input:   sm = scattering measure    (kpc m^{-20/3})
             nu = radio frequency   (GHz)
          vperp = psr transverse speed      (km/s)

    output: specbroad = spectral broadening bandwidth (Hz)

    usage: should be called with sm = smtau for appropriate
           line of sight weighting
    reference: eqn (47) of Cordes & Lazio 1991, ApJ, 376, 123.

    nb:
    The coeff. in the following line was 0.14 Hz from Cordes & Lazio (1991)
    It is changed to 0.097 to conform with FUNCTION SCINTIME and
    a new calculation consistent with Cordes & Rickett (1998)
    """

    specbroad = 0.097 * nu**(-1.2) * sm**0.6 * (vperp/100.) # Hz
    return specbroad


def theta_xgal(sm, nu):
    """

    calculates angular broadening for an extragalactic
    source of plane waves

    sm = scattering measure
    nu = radio frequency
    theta_xgal = angular broadening FWHM (mas)

    """
    theta_xgal = 128. * sm**0.6 * nu**(-2.2)
    return theta_xgal

def theta_gal(sm, nu):
    """
    calculates angular broadening for a galactic
    source of spherical waves

    sm = scattering measure
    nu = radio frequency
    theta_gal = angular broadening FWHM (mas)
    """

    theta_gal = 71. * sm**0.6 * nu**(-2.2)
    return theta_gal

def em(sm, louter=1., si=11./3.):

    """
    calculates the emission measure from the scattering measure
    using an assumed outer scale and spectral index of the
    wavenumber spectrum.

    Input:
        sm = scattering measure     (kpm^{-20/3})
        louter = outer scale        (pc)
        si = wavenumber spectral index (11/3 for Kolmogorov spectrum;
                                        the only valid value for now)
    Output:
        em = emission measure       (pcm^{-6})

    For a wavenumber spectrum P_n(q) = q^{-alpha} from q_0 to q_1
    the mean square electron density is

    <n_e^2> =~  4pi*[C_n^2 / (alpha - 3) ] * q_0^{3 - alpha)

    (an approximate form that assumes (q_0 / q_1)^{3-alpha} >> 1.

    Jim Cordes 18 Dec 1989
    """
    em = sm \
       * ( (4. * pi * 1000.) / (si-3.) ) \
       * (louter*pc/ (2. * pi) )**(si-3.) \
       * (0.01) ** (20./3.)
    return em

def theta_iso(smiso, nu):
    r"""
    Input:
       smiso in (kpm^{-20/3}) x kpc^{5/3}
       nu in GHz
    Output:
       isoplanatiangle in microarcsec

    Requires:
       re = classical electron radius (cm)
       kpin cm
    12 October 1998
    JMC

    \theta_{iso} = \delta r_s / d
                 = \left[
                         (\lambda r_e)^2 f_{\alpha} SM_{iso}
                  \right]^{1/\alpha}
    where \alpha = 5/3 for Kolmogorov case.
    NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
       so SM_{iso} does not have the units of scattering
       measure, but rather units of SM x Length^{\alpha}

    f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
    for \alpha = 5/3, f_{\alpha}= 88.3

    Constants in equation:
          2*0.6*log10(30cm*r_e  = 13.287
          0.6*log10(f_alpha)   = 1.1676
          1.6 * alog10(kpc)    = 34.383
          -(20/3)*log(100)     = 8 ??? CHECK (2020)
    """

    theta_log_radian = \
          13.287 \
        + 1.2 * alog10(nu) \
        - 1.1676 \
        - 0.6 * alog10(smiso) \
        - 34.383 \
        + 8.

    # log10(microarsec/rad) = 11.314425
    theta_log_microarcsec = theta_log_radian + 11.314425
    theta_iso = 10.**theta_log_microarcsec
    return theta_iso

def transition_frequency(sm, smtau, smtheta, dintegrate):
    r"""
    Returns the transition frequency between weak and strong scattering
    28 March 2001
    JMC

    input:
    (all sm values in (kpc m^{-20/3}))
              sm = int[\cnsq]
              smtau = int[(s/D)(1-s/D) \cnsq]
              smtheta = int[ (1-s/D) \cnsq]
              dintegrate = distance used to integrate \cnsq (kpc)
    output:
      transition_frequency = GHz given by
        \nu_t  = 318 GHz \\xi^{10/17}  SM^{6/17} D_{eff}^{5/17}
        (nb. double \\ needed on xi to avoid a unicode error in python)
      where
           D_{eff} = effective path length through medium
           D_{eff} = \int_0^dintegrate ds s \cnsq / \int_0^dintegrate ds  \cnsq

      Note we can calculate D_{eff} using
           D_{eff} = dintegrate * (sm - smtau/6 - smtau/3) / sm
    """

    xi= 0.3989                  # (2.*pi)^{-1/2} = fresnel scale definition factor
    coefficient=318.            # GHz; see NE2001 paper
    deff = (dintegrate*(sm - smtau/6. - smtheta/3.)) / sm
    transition_frequency \
         = coefficient * xi**(10./17.) * sm**(6./17.) * deff**(5./17.)
    return transition_frequency

def transition_frequency_from_obs(nu, sbw, si=11./3.):
    """
    Calculates the weak-strong transition frequency by scaling
    the scintillation bandwidth with nu and setting equal to nu.

    Input:
       nu = radio frequency of the input scintillation bandwidth (sbw)  (GHz)
       sbw = scintillation bandwidth (MHz)
    Output:
       transition frequency in GHz
    """
    if si == -2 or si == 2:
        print('transition_frequency_from_obs: si = %d'%(si), ' not allowed')
        return
    else:
        xsbw = (2 * si) / (si - 2)
        transition_frequency = nu * (1000. * nu / sbw)**(1./(xsbw-1))
        return transition_frequency

def taud_from_thetad(thetad, Deff):
    """
    Calculates pulse broading time from angular diameter.
    Works for any kind of medium (uniform, screen, Galactic scattering
    of extragalactic sources, etc.) by using an appropriate Deff.

    Input:
        thetad = measured scattering diameter           mas
        Deff = effective distance of scattering medium  kpc
               [e.g. for a Galactic layer of thickness Lg measured
                from the Sun,  Deff = Lg/2]
    Output:
        taud = pulse broadening time                    ms

    Method:
        Calculation relates the mean pulse broadening time to
        the mean-square scattering angle assuming a Gaussian
        scattered angular distribution.

    Dependencies:
        Need c = speed of light and kpc defined in cgs units.
        Need mas = milliarcsecond  defined.

    JMC 2020 Jan 1
    """
    taud = (Deff * kpc * (thetad * mas)**2) / (16 * np.log(2) * c) / ms
    return taud

def rsigma(nu, SM, linner=100., si=11./3.):
    r"""
    NOT YET COMPLETED

    Calculates the ratio R_\sigma(\beta) defined in pbcosmo.tex
    that is the ratio of the RMS 1D scattering angle to the RMS angle
    given by the FWHM of the equivalent Gaussian fitted to the scattered
    image.     This ratio is \ge 1.

    Input:
       nu = RF                                  GHz
       SM = scattering measure                  kpc m^{-20/3}
       linner = inner scale                     km
       si = spectral index of wavenumber spectrum (11/3 for Kolmogorov)

    Output:
       Rsigma = ratio as defined above

    Requires:
       r_e = classical electron radius
       c = speed of light
    """

    beta4 = 4-si
    beta2 = beta-si

    coeff = (pi * (2*pi)**(beta4/2.)) / (np.sqrt(beta4) * 4.*pi**2 * fbeta(si))**(1./beta2)
    return





def sm_from_thetad_plane_wave_meansq_method(nu, thetad, si=11./3., linner=100.):
    r"""
    Calculates SM from angular broadening diameter
    for a power-law wavenumber spectrum with spectral index si == beta.
    [nb. beta not used here to avoid confusion with the beta PDF]

    This method uses the mean-square angle and generally depends
    (weakly) on the inner scale.

    From A-O book snippet (Drafted in 2019)
    'Estimating the Wavenumber Spectrum for Electron Density in the ISM'
    Eq 33.

    Note that the angular diameter as measured from the visibility function
    of a scattered source is substantially less than
    [<\theta^2> / 2]^{1/2} for small linner and small scattering strength
    (with strength measured as \lambda^2 SM.   This is because of the wings
    on the image that extend to larger angles than for a Gaussian function.

    Thus proper use of this function requires relating thetad to
    the angular variance.    This relation depends on linner and
    \lambda^2 SM.  This is shown in the notes 'Cosmological Integrals
    for Dispersion and Scattering'  in pbcosmo.tex
    [Figure 4 shows the ratio of angles, R_{\sigma}].

    For linner = 100 km to 1000 km and SM/\nu^2 [\nu in GHz] \sim 10^{-4}, the
    RMS angle is larger than \thetad by a factor R_{\sigma} \sim 2.5 to 3.5.

    Input:
        nu = RF                                             GHz
        thetad = angular scattering diameter                mas
        si = 11./3. for Kolmogorov spectrum
        linner = inner scale (default = 100 km)             km

    Output:
        scattering measure                                  kpc m^{-20/3}

    Requires definitions in cgs units of
        r_e = classical electron radius
        km = kilometer
        GHz = 10^9 Hz
        mas = milliarcsecond
        SMunit
    as in get_constants_waveprop_astropy_version

    JMC 2020 Jan 1
    """

    beta4 = si - 4
    wavelen = c / (nu * GHz)

    print(si, beta4, wavelen, SMunit, mas)
    SM = (thetad * mas)**2 \
       / (4.*np.log(2.) * k_eta(si) * wavelen**4 * r_e**2 * (linner*km)**beta4)\
       / SMunit
    return SM

def sm_from_thetad_plane_wave_width_method(nu, thetad, si=11./3.):
    """
    Calculates SM from angular broadening diameter
    for a power-law wavenumber spectrum with spectral index si == beta.
    [nb. beta not used here to avoid confusion with the beta PDF]

    This method uses the width of the scattered image.
    In the strong but not super-strong scattering regime,
    the SM estimate is independent of the inner scale.

    From Cordes and Lazio (2003) NE2001 Paper II Eq. A10 with theta < theta_cross

    Result is consistent with Eq. A17.

    Input:
        nu = RF                                             GHz
        thetad = angular scattering diameter                mas
        si = 11./3. for Kolmogorov spectrum
    Output:
        scattering measure                                  kpc m^{-20/3}

    Requires definitions in cgs units of
        r_e = classical electron radius
        GHz = 10^9 Hz
        mas = milliarcsecond
        SMunit
    as in get_constants_waveprop_astropy_version

    JMC 2020 Jan 1
    """

    beta2 = si - 2.

    Cbeta = (pi / (2.*np.sqrt(np.log(2))))**beta2 \
          / (4*pi**2 * r_e**2 * c**si * fbeta(si))

    SM = Cbeta * (nu * GHz)**si * (thetad * mas)**beta2 \
       / SMunit
    return SM

def taux_from_thetax(d_mod,sm,smtau,smtheta,theta_x,tau,sbw):

	# Effective distance to dominant scattering
    deffsm2 = d_mod*(sm - smtau/6. - smtheta/3.)/sm

    # Alternative for extragalactic temporal broadening
    fwhm2sigma = 1. / sqrt(8.*np.log(2.))
    #taufactor = (((mas*fwhm2sigma)**2 * kpc) / (2.*c)) *  1000. # old - spurious factor of 2
    taufactor = (((mas*fwhm2sigma)**2 * kpc) / (c)) *  1000. # new - factor of 2 removed

    tau_x = taufactor * deffsm2 * theta_x**2            # ms
    sbw_x = sbw * tau / tau_x                           # MHz

    return tau_x,sbw_x,deffsm2


def sm_from_tau(nu, tau, Deff, si=11./3.):
    """
    Calculates SM from pulse broadening time and (effective) distance.

    Input:
        nu = RF                                             GHz
        tau = pulse broadening time                         ms
        Deff = effective distance to scattering region      kpc
    Output:
        scattering measure                                  kpc m^{-20/3}
    """

    # NEED TO FINISH THIS
