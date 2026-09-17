# mwprop.ne2001p v1.0 Jan 2024

from __future__ import print_function
# revised 2020 Feb 11 to remove any routines that use Fortran NE2001
# revised 2019 Dec 31 - 2020 Jan 1 to be a proper python package
# original 02 Apr 2017
# Utilities

import numpy as np
from math import *
from scipy import stats
from matplotlib.pyplot import rc, figure, axes, axis, xscale, yscale, \
                       plot, fill_between, xlabel, ylabel, annotate, \
                       title, legend, grid, savefig, show, close
import sys

# Import useful constants ... get them all (not that many)
# Note other routines in this subpackage see bare constants (e.g. KDM)
# but if this subpackage is imported in ipython as
# from pygedm.mwprop_core.scattering_functions import iss_mw_utils2020p as ut
# need ut.KDM to see the value of KDM.

from pygedm.mwprop_core.get_constants_waveprop_scipy_version import *

import pygedm.mwprop_core.scattering_functions.scattering_functions2020 as sf

input = input

#======================================================================`
def testit():
    print(KDM)
#======================================================================`

def theta_iso_at_RF(RF, SMiso, DEFFSM, scattype,  si=11./3.):
   r"""
   Calculates the isoplanatic angle in milliarcseconds at the specified RF.
   Assumes a Kolmogorov scaling (si=11/3).

   Input:  RF       GHz
           SM_iso   standard SM uits
           DEFFSM   kpc
           si = wavenumber spectral index for electron-density (e.g. 11/3)
   Output:
           THETA_ISO            mas at 1 GHz
           THETA_ISO_RF         mas at RF

   Notes:
   For strong scattering:
   From function THETA_ISO in scattering98.f (part of the NE2001 code package):
   \theta_{iso} = \delta r_s / d
                = \left [
                     (\lambda r_e)^2 f_{\alpha} SM_{iso}
                  \right ]^{1/\alpha}
   where
         \lambda = EM wavelength (cm)
         re = classical electron radius (cgs)
         \alpha = 5/3 for Kolmogorov case.
         SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
            so SM_{iso} does not have the units of scattering
            measure, but rather units of SM x Length^{\alpha}

         f_{\alpha} =
           8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
           for \alpha = 5/3, f_{\alpha}= 88.3

   theta_log_radian =
       13.287                                ! 0.6*log10(30cm*r_e)
     + 1.2 * np.log10(nu)
     - 1.1676                                ! 0.6*log10(f_alpha)
     - 0.6 * np.log10(smiso)
     - 34.383                                ! 1.6 * alog10(kpc)
     + 8.                                    ! -(20/3)*log(100)
   theta_log_microarcsec =
       theta_log_radian + 11.314425        ! 11.314425=alog10(microarsec/rad)

   Weak scattering:   the scattered source size is the Fresnel scale / distance.
       theta_iso = lambda / Deff

   """

   if scattype == 'strong' or scattype == 'transition':
      xiso = si / (si-2.) - 1.

      theta_log_radian = \
        13.287\
        - 1.1676 \
        - 0.6 * np.log10(SMiso) \
        - 34.383 \
        + 8.
      theta_log_microarcsec = \
          theta_log_radian + 11.314425        # 11.314425=alog10(microarsec/rad)
      THETA_ISO = 10.**theta_log_radian / mas
      THETA_ISO_RF = THETA_ISO * RF**xiso

   if scattype == 'weak':       # NEEDS CHECKING
      xiso = -0.5
      rF_1GHz = np.sqrt(wavelen1GHz * DEFFSM * kpc / (2.*pi))
      THETA_ISO = (wavelen1GHz / rF_1GHz) / mas     # mas at 1 GHz
      THETA_ISO_RF = THETA_ISO * RF**xiso

   return THETA_ISO, THETA_ISO_RF


#======================================================================`

def scale_diss_params_to_RF(rfratio, TAU, THETA_X, si=11./3.):

   """
   Scales relevant parameters outputed by NE2001 to the designated RF
   assuming a Kolmogorov scaling (si=11/3).

   Input:
           rfratio = (ratio of RF for output values) / (RF for input values)
           TAU          ms
           THETA_X      mas
           si = wavenumber spectral index for electron-density
   Output:
           TAU_RF           ms at RF
           THETA_X_RF       mas at RF
   """

   xsbw = 2.*si / (si-2.)
   xtau = -xsbw
   xthetax = -xsbw/2.

   TAU_RF = TAU * rfratio**xtau
   THETA_X_RF = THETA_X * rfratio**xthetax

   return TAU_RF, THETA_X_RF

#======================================================================`

def scale_NE2001_output_to_RF(RF,TAU,SBW,SCINTIME,THETA_G,THETA_X,si=11./3.):
   """
   Scales relevant parameters outputed by NE2001 to the designated RF
   assuming a Kolmogorov scaling (si=11/3).

   Input:
           TAU          ms
           SBW          MHz
           SCINTIME     sec
           THETA_G      mas
           THETA_X      mas
           RF           GHz
           si = wavenumber spectral index for electron-density
   Output:
           TAU_RF           ms at RF
           SBW_RF           MHz at RF
           SCINTIME_RF      sec at RF
           THETA_G_RF       mas at RF
           THETA_X_RF       mas at RF
   """

   xsbw = 2.*si / (si-2.)
   xtau = -xsbw
   xscintime = xsbw/2. - 1.
   xthetag = -xsbw/2.
   xthetax = xthetag

   TAU_RF = TAU * RF**xtau
   SBW_RF = SBW * RF**xsbw
   SCINTIME_RF = SCINTIME * RF**xscintime
   THETA_G_RF = THETA_G * RF**xthetag
   THETA_X_RF = THETA_X * RF**xthetax

   return TAU_RF, SBW_RF, SCINTIME_RF, THETA_G_RF, THETA_X_RF

#======================================================================`

def ggpdf(x, rms1, rms2):
   """
   Calculates gamma-gamma PDF for combined diffractive and refractive
   scintillations given in Eq. 20 of Prokes,
   'Modeling of Atmospheric Turbulence Effect on Terrestrial FSO Link'

   """
   from math import gamma
   from scipy.special import kv
   from scipy.special import kve
   from scipy.special import gammaln

   a = 1./rms1**2
   b = 1./rms2**2

   e1 = (a+b)/2.
   nu = a-b
   arg1 = 2. * np.sqrt(a*b*x)

   #coeff1 = (2. * (a*b)**e1) / (gamma(a) * gamma(b))
   coeffln = np.log(2.) + e1*np.log(a*b) - gammaln(a) - gammaln(b)
   kve = kve(nu, arg1)
   pdfln = coeffln + (e1-1.)*np.log(x) + np.log(kv(nu, arg1))
   #pdf = coeff1 * x**(e1-1) * kv(nu, arg1)
   pdf = np.exp(pdfln)
   return pdf

#======================================================================`

def mriss_goodman_narayan(si, phif):
   r"""
   Returns the RISS modulation index using equation 3.1.5 of
   Goodman & Narayan 1985.

   Input:
      si = beta = wavenumber index of electron density power spectrum
      e.g. beta = 11/3 for Kolmogorov case

      Applicable for 2 < si < 4.

      phif = rms phase difference on the Fresnel scale
             with Fresnel scale for plane-wave incidence on a phase screen
             at distance d: rf = \lambda d / 2\pi.

             The "u" parameter is calculated from phif as
             u = phif^{2/(\beta -2} in the evaluation.
   Output:

      mriss

   """
   from math import gamma

   si2 = si-2.
   si4 = 4.-si
   si6 = 6.-si
   Crsq = 2.**(-si*si4/si2) * gamma(si/2.)*gamma(si6/si2) / gamma(si6/2.)

   u = phif**(2./si2)

   mriss = np.sqrt(Crsq) * 2.**(si4/si2**2) * u**(-si4)

   return mriss

#======================================================================`

def lnorm_from_mean_rms(x, mean, rms):
   """
   Calculates lognorm pdf for specificed x range and mean, rms of x.
   Note:

   mean = exp(mu+sigma**2/2)
   median = exp(mu) = 1.05
   variance = (e^{sigma^2}-1) e^{2*mu+sigma^2}

   where mu, sigma are standard parameters for the log-normal pdf
   """

   sigma = np.sqrt(np.log(1 + (rms/mean)**2))
   mu = np.log(mean) - sigma**2/2.

   # lognorm.pdf seems to give the wrong input rms
   # even though it has unit area and unit mean when
   # I input mean = 1. and rms = 0.3
   # it gives lognorm.stats(sigma, mu) = 1., 0.313
   #pdfx = lognorm.pdf(x, sigma, mu, 1.)

   # so calculate directly: this gives correct input mean and rms:

   pdfx = (x * sigma * np.sqrt(2.*pi))**(-1.) * np.exp(-(np.log(x)-mu)**2/(2.*sigma**2))

   return mu, sigma, pdfx


#======================================================================`

def calculate_fresnel_scale_and_mod_indices(RF, DEFFSM, THETA_X_RF, si=11./3.):
    """
    Calculates phiF and modulation indices.
    si = index of wavenumber spectrum for ne
    """

    wavelength = c / (RF*GHz)
    rF = np.sqrt(DEFFSM *kpc * wavelength / (2.*pi))    # Fresnel scale cm

    phiF = np.sqrt(2.) * \
       ( (pi*THETA_X_RF*mas * rF / wavelength) / (2.*np.sqrt(np.log(2.))) )**((si-2.)/2.)

    u = phiF**(2./(si-2.))

    """
    Modulation indices for different regimes of u:
    use definitions of md, mr that are continuous across the
    weak/strong boundary. Do so by using an equipartition approach for
    u << 1 but one that asymptotically gives
    reasonable modulation indices md = 1 and mr = (2u)^{-1/2}
    """

    mdp = np.sqrt(2.)*u / np.sqrt(1. + 2.*u**2)

    mrp = np.sqrt(2.) * u * (2*u)**(si-4.) / np.sqrt((2*u)**(2.*(si-4.)) + 2.*u**2)

    mgp = np.sqrt((1.+mdp**2)*(1.+mrp**2)-1.)
    if u <= 1.:     # no distinction between r,d
       lg = rF
       ld = rF
       lr = rF
       scattype = 'weak'
    if u > 1.:
       ld = rF/u
       lr = rF*u
       if u <= 2:
          scattype = 'transition'
       else:
          scattype = 'strong'
    return rF, ld, lr, phiF,  u,  mgp, mdp, mrp, scattype

#======================================================================`

def calculate_diss_bandwidth_factor(sbw, bw):
   """
   Calculates attentuation of DISS modulation index by bandwidth.
   Input:
      sbw = scintillation bandwidth
      bw = receiver bandwidth

      note sbw, bw must have same units
   Output:
      bwfactor
   """
   eta_nu = 0.3
   ndiss = 1.+eta_nu * bw / sbw
   bwfactor = 1./np.sqrt(ndiss)
   return ndiss, bwfactor

#======================================================================`

def calculate_diss_source_size_factor(theta_source, theta_iso):
   """
   Calculates attentuation of DISS modulation index by source size.
   Input:
      theta_source = source size (angular units)
      theta_iso = isoplanatic size (angular units)

      note theta_source, theta_iso must  have same units
   Output:
      source_size_factor
   """
   #source_size_factor = min(1., theta_iso/theta_source)
   # changed 2016 July 5
   source_size_factor = (1. + (theta_source/ theta_iso)**2)**(-1./2.)
   return source_size_factor

#======================================================================`

def invert_function(x, func, level):
   """
   Finds the value of x at which func = level.
   Interpolation is done.
   """
   xlevel = np.interp(level, func, x)
   return xlevel


#======================================================================`

def calc_pdfg_for_xgal_los(RF, BW, RF_input, TAU, THETA_X, Deff, SM,
    theta_source=1.e-6, dg=1.e-5, gmin=1.e-5, gmax=30., si=11./3.):

    """
    Calculates the PDF and CDF for the  ISS modulation for the
    a line of sight characterized by TAU, THETA_X, Deff, SM
    and for the specified RF and bandwidth.

    RF              Radio frequency of evaluated PDF, CDF       GHz
    BW              Bandwidth                                   MHz
    RF_input        Radio frequency of input ISS values         GHz
    TAU             Pulse broadening time                        ms
    THETA_X         Angular diameter of extragalactic
                    source caused by MW scattering              mas
    Deff            Effective distance to MW scattering region  kpc
    SM              Scattering measure                          kpc m^{-20/}
    theta_source    Source size as seen by ISM                  mas
    dg,gmin,gmax    Sample interval and Range for gain g in pdf,cdf
    si              Spectral index of electron density
                    wavenumber spectrum (11/3 for Kolmogorov)

    Note d_xgal_kpc has little effect; it's mainly to ensure that the NE2001
    code integrates out to its maximum distance.
    """

    C1 = 1.16                       # factor in SBW-TAU relation

    # Scale input parameters to RF of interest if necessary:
    if RF == RF_input:
        TAU_RF = TAU
        THETA_X_RF = THETA_X
    else:
        rfratio = RF / RF_input
        TAU_RF, THETA_X_RF = \
         scale_diss_params_to_RF(rfratio, TAU, THETA_X, si=si)

    rF, ld, lr, phiF, u, mgp, mdp, mrp, scattype = \
       calculate_fresnel_scale_and_mod_indices(RF, Deff, THETA_X_RF,si=si)

    THETA_ISO, THETA_ISO_RF = theta_iso_at_RF( \
                              RF, SM, Deff, scattype, si=si)

    SCINTIME_RF = sf.scintime(SM, RF, vperp=100)        # for 100 km/s

    SBW_RF = 1.e-3 * C1 / (2. * pi * TAU_RF)            # MHz

    nu_transition = sf.transition_frequency_from_obs(RF, SBW_RF, si=si)

    ndiss, bandwidth_factor = calculate_diss_bandwidth_factor(SBW_RF, BW)

    source_size_factor =  \
             calculate_diss_source_size_factor(theta_source, THETA_ISO_RF)

    # set modulation indices using source size and bandwidth factors:

    gmean = 1.
    mr = mrp
    md = mdp * source_size_factor * bandwidth_factor


    if scattype == 'strong':
       stypelab = 'S'
       mg = np.sqrt(md**2 + mr**2 + (md*mr)**2) # ok for strong scattering
       mdtrans = 0.5 + 0.25*mr

       if md <= mdtrans:
         pdftype = 'log-normal'
         pdftypelab = 'LN'

       if md > mdtrans and mr > 0.1:
         pdftype = 'gamma-gamma'
         pdftypelab = 'GG'

       if md > mdtrans and mr <= 0.1:
         pdftype = 'chi-square'
         pdftypelab = 'X'
         dof_g = 2. / mg**2         # DoF for chi^2 pdf

    if scattype == 'transition':    #
         stypelab = 'T'
         pdftype = 'gamma-gamma'
         pdftypelab = 'GG'

    if scattype == 'weak':      # log-normal PDF
         stypelab = 'W'
         mg = mgp * source_size_factor
         sigma = np.sqrt(np.log(1.+mg**2))
         mu = -sigma**2/2.
         pdftype = 'log-normal'
         pdftypelab = 'LN'

    # pdf and cdf:
    gvec = np.arange(gmin, gmax+dg, dg)
    wvec = np.ones(np.size(gvec))
    wvec[0] = wvec[-1] = 0.5

    if pdftype == 'log-normal':
       pdflabel = r'$\rm LN$'
       muln, sigmaln, gpdf = lnorm_from_mean_rms(gvec, gmean, mg)
       gcdf = (np.cumsum(wvec*gpdf) + lnorm_from_mean_rms(dg/2., gmean, mg)[2]) * dg

    if pdftype == 'chi-square':
       pdflabel = r'$\chi^2_{\rm n}$'
       loc = 0.
       scale = 1. / dof_g       # gives pdf of reduced chi-square quantity
       gpdf = chi2.pdf(gvec, dof_g, loc, scale)
       gcdf = (np.cumsum(wvec*gpdf) + chi2.pdf(dg/2., dof_g, loc, scale)) * dg

    if pdftype == 'gamma-gamma':
       pdflabel = r'$\gamma\gamma}$'
       gpdf = ggpdf(gvec, md, mr)
       gcdf = (np.cumsum(wvec*gpdf) + ggpdf(dg/2., md, mr)) * dg

    # find median and probability ranges for g:
    gmedian =  invert_function(gvec, gcdf, 0.5)

    # Define new dictionary that contains output values:

    Disskeys = \
        np.array(['RF', 'tau', 'SM', 'tau_rf', 'theta_x', 'theta_x_rf', 'rFresnel',
       'ld', 'lr', 'phiF', 'u', 'scattype', 'mr', 'md', 'mdp', 'mg',
       'gmedian', 'stypelab', 'pdftype',  'pdflabel', 'gvec', 'gpdf',
       'gcdf', 'nu_transition', 'theta_iso', 'theta_iso_rf', 'sbw_rf',
       'scintime_rf', 'bandwidth_factor', 'source_size_factor'])

    Dissvals = \
        list([RF, TAU, SM, TAU_RF, THETA_X, THETA_X_RF, rF, ld, lr, phiF, u, scattype, mr, md, mdp, mg, gmedian, stypelab, pdftype, pdflabel, gvec, gpdf, gcdf, nu_transition, THETA_ISO, THETA_ISO_RF, SBW_RF, SCINTIME_RF, bandwidth_factor, source_size_factor])

    Dissunits = np.array([
        'GHz', 'ms', 'kpc m^{-20/3}', 'ms', 'mas', 'mas', 'cm', 'cm', 'cm', 'rad', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'GHz', 'mas', 'mas', 'MHz', 's', '', '' ])

    Dissdesc = np.array([
        'RF of output values', 'TAU input', 'SM input', 'TAU @RF', 'XGalFWHM in',  'XgalFWHM @RF', 'Fresnel scale', 'Diffraction scale', 'Refraction scale', 'Fresnel phase', 'U parameter', 'Scattering regime', 'Refraction mod index', 'Diffraction mod index', 'Diffraction mod index (pt source)', 'Mod index (total)', 'Median mod', 'Scattering label', 'PDF type', 'PDF label', 'Gain vector', 'Gain PDF', 'GAIN CDF', 'Transition RF', 'Isoplanatic angle input', 'Isoplanatic angle @RF', 'DISS bandwidth @RF', 'DISS time scale @RF', 'Bandwidth factor', 'Source size factor'])

    Dissv, Dissu, Dissd = {}, {}, {}
    for n, key in enumerate(Disskeys):
       Dissv[key] = Dissvals[n]
       Dissu[key] = Dissunits[n]
       Dissd[key] = Dissdesc[n]

    """
    Can delete these lines once Dissv,u,d are confirmed accurate
    D_iss = {}
    D_iss['RF'] = RF
    D_iss['tau'] = TAU
    D_iss['SM'] = SM
    D_iss['tau_rf'] = TAU_RF
    D_iss['theta_x'] = THETA_X
    D_iss['theta_x_rf'] = THETA_X_RF
    D_iss['rFresnel'] = rF
    D_iss['ld'] = ld
    D_iss['lr'] = lr
    D_iss['phiF'] = phiF
    D_iss['u'] = u
    D_iss['scattype'] = scattype
    D_iss['mr'] = mr
    D_iss['md'] = md
    D_iss['mdp'] = mdp
    D_iss['mg'] = mg
    D_iss['gmedian'] = gmedian
    D_iss['stypelab'] = stypelab
    D_iss['pdftypelab'] = pdftypelab
    D_iss['pdflabel'] = pdflabel
    D_iss['gvec'] = gvec
    D_iss['gpdf'] = gpdf
    D_iss['gcdf'] = gcdf
    D_iss['nu_transition'] = nu_transition
    D_iss['theta_iso'] = THETA_ISO
    D_iss['theta_iso_rf'] = THETA_ISO_RF
    D_iss['sbw_rf'] = SBW_RF
    D_iss['scintime_rf'] = SCINTIME_RF
    D_iss['bandwidth_factor'] = bandwidth_factor
    D_iss['source_size_factor'] = source_size_factor
    return D_iss, Disskeys, Dissvals
    """

    return Dissv, Dissu, Dissd

#======================================================================`

def main(l, b, RF, BW, dxgal_mpc, theta_source, SM, DMmodel, doplot=False):
    # setup is so that a plot is done if run as a program

    # output line is appended to file xgal_mw_iss_output.txt
    # some values and arrays (gvec, gpdf, gcdf) are saved into an npz file,
    # xgal_mw_iss_save.npz


    print("input = ", l,b,RF,BW,dxgal_mpc,theta_source)

    """
    Evaluates DISS and RISS modulation indices and probabilities for the net gain
    given a direction (l,b), frequency \nu, bandwidth B, and extragalactic source size in milliarcseconds.
    """
    outfile = 'iss_mw_pdf_package.txt'
    outfile2 = 'iss_mw_pdf_package_lb_granges.txt'
    fout = open(outfile, 'a')
    gout = open(outfile2, 'a')
    hout = open('junkfile_iss_package', 'a')

    D_iss = calc_pdfg_for_xgal_los(
                  RF, BW, RF_input, TAU, THETA_X, Deff, SM,
                  theta_source=1.e-6, dg=1.e-5, gmin=1.e-5, gmax=30.)

    gvec = D_iss['gvec']
    gpdf = D_iss['gpdf']
    gcdf = D_iss['gcdf']
    NU_T = D_iss['nu_transition']
    SM = D_iss['SM']
    THETA_X_RF = D_iss['theta_x_rf']
    THETA_ISO = D_iss['theta_iso']
    THETA_ISO_RF = D_iss['theta_iso_rf']
    SBW_RF = D_iss['sbw_rf']
    phiF = D_iss['phiF']
    rF = D_iss['rFresnel']
    lr = D_iss['lr']
    ld = D_iss['ld']
    bandwidth_factor = D_iss['bandwidth_factor']
    source_size_factor = D_iss['source_size_factor']
    mr = D_iss['mr']
    md = D_iss['md']
    mdp = D_iss['mdp']
    mg = D_iss['mg']
    gmedian = D_iss['gmedian']
    stypelab = D_iss['stypelab']
    pdftypelab = D_iss['pdftypelab']
    pdflabel = D_iss['pdflabel']
    scattype = D_iss['scattype']
    SCINTIME_RF =  D_iss['scintime_rf']

    p80 = 0.8
    g80low = invert_function(gvec, gcdf, (1.-p80)/2.)
    g80high = invert_function(gvec, gcdf, (1.+p80)/2.)

    p90 = 0.9
    g90low = invert_function(gvec, gcdf, (1.-p90)/2.)
    g90high = invert_function(gvec, gcdf, (1.+p90)/2.)

    p98 = 0.98
    g98low = invert_function(gvec, gcdf, (1.-p98)/2.)
    g98high = invert_function(gvec, gcdf, (1.+p98)/2.)

    p998 = 0.998
    g998low = invert_function(gvec, gcdf, (1.-p998)/2.)
    g998high = invert_function(gvec, gcdf, (1.+p998)/2.)

    p9998 = 0.9998

    g9998low = invert_function(gvec, gcdf, (1.-p9998)/2.)
    g9998high = invert_function(gvec, gcdf, (1.+p9998)/2.)

    print(\
       "%6.1f %5.1f %4.1f %4.1f %7.1f %7.1f %6.3e %6.2f %4.1e %5.2f \
       %8.4e %8.2f %5.2f %5.2f %5.2f %4.3e %6.3e %6.3e %6.3e %5.3f \
       %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \
       %2s %2s"\
       %(l, b, RF, BW, NU_T, DMmodel, THETA_ISO_RF*1000., THETA_X_RF, \
         SM, Deff, SBW_RF, phiF, np.log10(rF), np.log10(lr), np.log10(ld), \
         bandwidth_factor, mr, md, mg, gmedian,\
         g80low, g80high,g90low, g90high, g98low, g98high, \
         g998low, g998high,g9998low, g9998high,\
         stypelab, pdftypelab), file=fout)

    print( "DM = ", DMmodel)
    print( "SM = ", SM)
    print( "Transition frequency = ", NU_T)
    print( "Theta xgal = ", THETA_X_RF, " mas")
    print( "Effective distance = ", Deff)

    print( "phiF = ", phiF, " md = ", md, " mr = ", mr, " mg = ", mg, " SBW = ", SBW_RF, " bw_fac = ", bandwidth_factor)

    print("%6.1f %5.1f %4.1f %4.1f %7.1f %7.1f  %4.1e %5.2f %7.1f %4.1f %4.1f %4.1f %4.3e  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f"%(l,b,RF,BW, NU_T, DMmodel, SM, Deff, phiF, np.log10(rF), np.log10(lr), np.log10(ld), bandwidth_factor,  mg, gmedian, g80low, g80high, g90low, g90high, g98low, g98high,g998low, g998high,g9998low, g9998high), file=gout)

    fout.close()
    gout.close()
    hout.close()

    # save key output
    savefile = 'iss_mw_pdf_package_save'
    np.savez(savefile, l=l, b=b, RF=RF, BW=BW, dxgal_mpc=dxgal_mpc, theta_source=theta_source, gvec=gvec, gpdf=gpdf, gcdf=gcdf, scattering=scattype, pdftype=pdflabel, mdp=mdp, md=md, mr=mr, mg=mg, gmedian=gmedian, g80low=g80low, g80high=g80high, g90low=g90low, g90high=g90high, g98low=g98low, g98high=g98high,g998low=g998low, g998high=g998high, g9998low=g9998low, g9998high=g9998high, DMmodel=DMmodel, Deff=Deff, SM=SM, SBW_RF = SBW_RF, phiF=phiF, bandwidth_factor=bandwidth_factor)

    if doplot:
      # Plotting:

      # not sure if these rc specs are working
      rc('font', family='serif')
      rc('text', usetex=True)
      rc('font', serif='cm')

      # PDF
      fig = figure()
      ax = fig.add_subplot(111)
      xlabel(r'$\rm g = g_{\rm d} \  g_{\rm r}$', fontsize=15)
      ylabel(r'$\rm Probability \ Density \ \ f_g(g)$', fontsize=15)
      plot(gvec, gpdf, 'k-', label = scattype + ',   ' + pdflabel)
      yscale('log')
      axis(xmin = -0.05)
      axis(ymin = 1.e-10, ymax = 2.*gpdf.max())

      aa,bb,cc,dd = axis()
      plot((gmedian, gmedian), (cc,dd), 'k', dashes=[20,5], lw=3)
      plot((g90low, g90low), (cc,dd), 'k', dashes=[20,3], lw=1)
      plot((g90high, g90high), (cc,dd), 'k', dashes=[10,3], lw=1)

      annotate(r'$\rm l,b \ = \  %4.1f^{\circ}, %4.1f^{\circ} \ \ \ DM_{NE2001} = %6.1f \ pc\ cm^{-3}  $'%(l,b,DMmodel), xy=(0.5, 1.03), xycoords='axes fraction', ha='center', fontsize=16)

      xanno = 0.3
      xanno = 0.95
      ha = 'right'
      yanno = 0.925
      dyanno = 0.06
      fontsize = 13
      annotate(r'$\rm \nu  =   %4.2f \ GHz \ \ \ \  B  =  %4.2f\ GHz \ \ \ \ \phi_{\rm F}   =  %6.1f \ rad $'%(RF, BW/1000., phiF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  r_{\rm F} =   10^{%4.1f} \ cm\ \ \ \  l_d  =  10^{%4.1f}\ cm \ \ \ \ l_r   =  10^{%4.1f} \ cm$'%(np.log10(rF), np.log10(ld), np.log10(lr)), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  \Delta\nu_d  =   {%4.1f} \ MHz\ \ \ \  \Delta t_d=  {%4.1f}\ s\ \ \ \ D_{\rm eff}   =  %4.2f \ kpc$'%(SBW_RF, SCINTIME_RF, Deff), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  \theta_{\rm iso}  =   {%4.1f} \ \mu as\ \ \ \  \theta_{\rm MW} =  {%4.1f}\ mas\ $'%(THETA_ISO*1000., THETA_X_RF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      legend(loc=(0.7, 0.6))
      show()
      savefig('xgal_mw_iss_pdf_l_' + str(l) + '_b_' + str(b) + '_rf_' + str(RF) + '_bw_' + str(BW)+ '.pdf')

      # CDF
      fig = figure()
      ax = fig.add_subplot(111)
      rc('xtick', labelsize=16)
      rc('ytick', labelsize=16)

      xlabel(r'$\rm g = g_{\rm d} \  g_{\rm r}$', fontsize=15)
      ylabel(r'$\rm CDF\ \ F_g(g)$', fontsize=15)
      plot(gvec, gcdf, 'k-', label = scattype + ',   ' + pdflabel)
      yscale('log')
      axis(xmin = -0.05)
      axis(ymin = 1.e-3, ymax=2.)

      aa,bb,cc,dd = axis()
      plot((gmedian, gmedian), (cc,dd), 'k', dashes=[20,5], lw=3)
      plot((g90low, g90low), (cc,dd), 'k', dashes=[20,3], lw=1)
      plot((g90high, g90high), (cc,dd), 'k', dashes=[10,3], lw=1)

      annotate(r'$\rm l,b \ = \  %4.1f^{\circ}, %4.1f^{\circ} \ \ \ DM_{NE2001} = %6.1f \ pc\ cm^{-3}  $'%(l,b,DMmodel), xy=(0.5, 1.03), xycoords='axes fraction', ha='center', fontsize=16)

      xanno = 0.3
      xanno = 0.95
      ha = 'right'
      yanno = 0.225
      dyanno = 0.06
      fontsize = 13
      annotate(r'$\rm \nu  =   %4.2f \ GHz \ \ \ \  B  =  %4.2f\ GHz \ \ \ \ \phi_{\rm F}   =  %6.1f \ rad $'%(RF, BW/1000., phiF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  r_{\rm F} =   10^{%4.1f} \ cm\ \ \ \  l_d  =  10^{%4.1f}\ cm \ \ \ \ l_r   =  10^{%4.1f} \ cm$'%(np.log10(rF), np.log10(ld), np.log10(lr)), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  \Delta\nu_d  =   {%4.1f} \ MHz\ \ \ \  \Delta t_d=  {%4.1f}\ s\ \ \ \ D_{\rm eff}   =  %4.2f \ kpc$'%(SBW_RF, SCINTIME_RF, Deff), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  \theta_{\rm iso}  =   {%4.1f} \ \mu as\ \ \ \  \theta_{\rm MW} =  {%4.1f}\ mas\ $'%(THETA_ISO*1000., THETA_X_RF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      legend(loc=(0.7, 0.8))
      show()
      savefig('xgal_mw_iss_cdf_l_' + str(l) + '_b_' + str(b) + '_rf_' + str(RF) + '.pdf')

      print( "Plotting Complementary CDF")
      doextremes_default = True
      #doextremes = bool(input('Plot extreme value examples? [%s]: '%(doextremes_default)) or str(doextremes_default))


      # CCDF = 1 - CDF
      fig = figure()
      ax = fig.add_subplot(111)
      rc('xtick', labelsize=16)
      rc('ytick', labelsize=16)
      xlabel(r'$\rm g = g_{\rm d} \  g_{\rm r}$', fontsize=16)
      ylabel(r'$\rm Complementary \ CDF\ \ 1-F_g(g)$', fontsize=16)
      plot(gvec, 1.-gcdf, 'k-', label = scattype + ',   ' + pdflabel)
      yscale('log')
      #yticklocs, yticklabs = calcticks.calcticks_log(tlog1=-10, tlog2=0, dt = 1)
      #ax.set_yticks(yticklocs)
      #ax.set_yticklabels(yticklabs)
      axis(xmin = -0.05)
      axis(ymin = 1.e-10, ymax=2.)

      aa,bb,cc,dd = axis()

      # median dashed line:
      plot((gmedian, gmedian), (cc,0.5), 'k', dashes=[20,5], lw=1)

      idx = np.where((g90low <= gvec) & (gvec <=g90high))
      # idx contains a lot of points (for sample interval used on gvec)
      # so the plot file ends up being very large.   Therefore just use
      # the edges and midpoint.

      i00 = idx[0][0]
      i0half = idx[0][np.int(np.size(idx)/2)]
      i01 = idx[0][-1]
      gfill = ((gvec[i00], gvec[i0half], gvec[i01]))
      gccdffill = ((1.-gcdf[i00], 1.-gcdf[i0half], 1.-gcdf[i01]))

      if fill90:
         fill_between(gfill, gccdffill, edgecolor='white', facecolor='cyan', alpha=0.2)
      else:
         plot((g90low, g90low), (cc, 0.95), 'k', dashes=[20,3], lw=1)
         plot((g90high, g90high), (cc,0.05), 'k', dashes=[10,3], lw=1)

      xanno = 0.3
      xanno = 0.95
      ha = 'right'
      yanno = 0.925
      dyanno = 0.06
      fontsize = 13
      annotate(r'$\rm l,b \ = \  %4.1f^{\circ}, %4.1f^{\circ} \ \ \ DM_{NE2001} = %6.1f \ pc\ cm^{-3}  $'%(l,b,DMmodel), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno
      annotate(r'$\rm \nu  =   %4.2f \ GHz \ \ \ \  B  =  %4.2f\ GHz \ \ \ \ \phi_{\rm F}   =  %6.1f \ rad $'%(RF, BW/1000., phiF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      annotate(r'$\rm  r_{\rm F} =   10^{%4.1f} \ cm\ \ \ \  l_d  =  10^{%4.1f}\ cm \ \ \ \ l_r   =  10^{%4.1f} \ cm$'%(np.log10(rF), np.log10(ld), np.log10(lr)), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      yanno -= dyanno

      if SBW_RF > 1.:
           annotate(r'$\rm  \Delta\nu_d  =   {%4.1f} \ MHz\ \ \ \  \Delta t_d=  {%5.0f}\ s\ \ \ \ D_{\rm eff}   =  %4.2f \ kpc$'%(SBW_RF, SCINTIME_RF, Deff), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)
      else:
           annotate(r'$\rm  \Delta\nu_d  =   {%4.2f} \ MHz\ \ \ \  \Delta t_d=  {%5.0f}\ s\ \ \ \ D_{\rm eff}   =  %4.2f \ kpc$'%(SBW_RF, SCINTIME_RF, Deff), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)
      yanno -= dyanno

      annotate(r'$\rm  \theta_{\rm iso}  =   {%4.1f} \ \mu as\ \ \ \  \theta_{\rm MW} =  {%4.1f}\ mas\ $'%(THETA_ISO*1000., THETA_X_RF), xy=(xanno, yanno), xycoords='axes fraction', ha=ha, fontsize=fontsize)

      if doextremes:  # plot points for specified probability levels
           p3 = 1.e-3
           p8 = 1.e-8
           # note second arg needs to be an increasing function
           # so that is why it operates on gcdf instead of 1.-gcdf
           g3 = invert_function(gvec, gcdf, 1.-p3)
           g8 = invert_function(gvec, gcdf, 1.-p8)

           plot(g3, p3, 'ro')
           plot(g8, p8, 'ro')

           aa,bb,cc,dd = axis()
           plot((g3, g3), (cc, p3), 'r:', lw=2)
           plot((g8, g8), (cc, p8), 'r:', lw=2)


           legend(loc=(0.675, 0.5))
      show()
      if fill90:
         savefig('xgal_mw_iss_ccdf_l_' + str(l) + '_b_' + str(b) + '_rf_' + str(RF) + '_bw_' + str(BW) + '_fill'  + '.pdf')
      else:
         savefig('xgal_mw_iss_ccdf_l_' + str(l) + '_b_' + str(b) + '_rf_' + str(RF) + '_bw_' + str(BW) + '.pdf')
      input('hit return')
      close()

#======================================================================`

# Main

if __name__ == "__main__":

    fill90 = False
    fill90 = True

    narg=0
    print(sys.argv)
    for arg in sys.argv:
        #print narg, arg
        if  narg == 1: l = float(arg)
        if narg == 2: b = float(arg)
        if narg == 3: RF = float(arg)
        if narg == 4: BW = float(arg)
        if narg == 5: dxgal_mpc = float(arg)
        if narg == 6: theta_source = float(arg)
        if narg == 7: doplot = arg
        narg += 1
    narg -= 1
    print("narg = ", narg)

    if narg == 0:
        doplot = True
        l = float(input('enter Galactic longitude in degrees: '))
        b = float(input('enter Galactic latitude in degrees: '))
        RF = float(input('enter RF in GHz: '))
        BW = float(input('enter BW in MHz: '))
        dxgal_mpc=float(input('enter distance of extragalactic source (Mpc): '))
        theta_source=float(input('extragalactic source size (mas): '))

    RF_input = 1.
    TAU = 1.
    THETA_X = 1.
    Deff = 2.5
    SM = 10.**(-3.5)
    DMmodel = 100

    D_iss = calc_pdfg_for_xgal_los(
                  RF, BW, RF_input, TAU, THETA_X, Deff, SM,
                  theta_source=1.e-6, dg=1.e-5, gmin=1.e-5, gmax=30.)
    main(l, b, RF, BW, dxgal_mpc, theta_source, SM, DMmodel)
