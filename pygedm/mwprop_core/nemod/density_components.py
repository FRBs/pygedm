# mwprop v2.0 Jan 2026

'''
Pythonic versions of the following functions in density.NE2001.f:

ne_outer
ne_inner
ne_gc

Change Log:

01/23/20 -- SKO
01/29/20 -- JMC
    * added return quantities to ne functions
    * moved input parameter read to master script
02/08/20 -- JMC
    * replaced parameter input read with import of config_ne2001p
    * rsun = 8.5 now commented out; set in config file
    * other parameters and function defs also now in config file
'''

import numpy as np

# nemod_config sets up the model with all dictionaries etc.
from pygedm.mwprop_core.nemod.config_nemod import *
from pygedm.mwprop_core.nemod.numba_compat import njit, HAS_NUMBA

pihalf = np.pi/2.
sqrt = np.sqrt

@njit
def _ne_outer_jit(x, y, z, rsun, A1, n1h1, h1, F1, pihalf, sech2_const):  # pragma: no cover
    """JIT-compiled ne_outer core calculation"""
    rr = sqrt(x**2 + y**2)
    suncos = np.cos(pihalf*rsun/A1)
    if rr > A1:
        g1 = 0.
    else:
        g1 = np.cos(pihalf*rr/A1)/suncos
    # sech2(z/h1) ≈ 4*exp(-2*|z/h1|) / (1 + exp(-2*|z/h1|))^2
    z_arg = np.abs(z/h1)
    exp_2z = np.exp(-2.0*z_arg)
    sech2_val = 4.0 * exp_2z / (1.0 + exp_2z)**2
    ne1 = (n1h1/h1) * g1 * sech2_val
    return ne1, F1

@njit
def _ne_inner_jit(x, y, z, A2, n2, h2, F2, sech2_const):  # pragma: no cover
    """JIT-compiled ne_inner core calculation"""
    g2 = 0.
    rr = sqrt(x**2. + y**2.)
    rrarg = ((rr-A2)/1.8)**2.
    if rrarg < 10.:
        g2 = np.exp(-rrarg)
    # sech2(z/h2)
    z_arg = np.abs(z/h2)
    exp_2z = np.exp(-2.0*z_arg)
    sech2_val = 4.0 * exp_2z / (1.0 + exp_2z)**2
    ne2 = n2 * g2 * sech2_val
    return ne2, F2

def ne_outer(x,y,z): #thick disk component

    #g1 = sech2(rr/A1)/sech2(rsun/A1) #TC93 function
    return _ne_outer_jit(x, y, z, rsun, A1, n1h1, h1, F1, pihalf, None)

def ne_inner(x,y,z): #thin disk component

    return _ne_inner_jit(x, y, z, A2, n2, h2, F2, None)

def ne_gc(x,y,z, absymax=2*rgc):

    '''
    Determines the contribution of the Galactic center to the free
    electron density of the interstellar medium at Galactic location
    (x,y,z).  Combine with `fluctuation' parameter to obtain the
    scattering measure.

    NOTE: This is for the hyperstrong scattering region in the
    Galactic center.  It is distinct from the inner Galaxy
    (component 2) of the TC93 model.

    Origin of coordinate system is at Galactic center; the Sun is at
    (x,y,z) = (0,rsun,0), x is in l=90 direction

    Based on Section 4.3 of Lazio & Cordes (1998, ApJ, 505, 715)

    Input:
    x - location in Galaxy [kpc]
    y - location in Galaxy [kpc]
    z - location in Galaxy [kpc]

    COMMON:
    NEGC0 - nominal central density

    PARAMETERS:
    RGC - radial scale length of Galactic center density enhancement
    HGC - z scale height of Galactic center density enhancement

    Output:
    NE_GC - Galactic center free electron density contribution [cm^-3]

    Definitions now in config_ne2001p.py
    xgc = Dgc['xgc']
    ygc = Dgc['ygc']
    zgc = Dgc['zgc']
    rgc = Dgc['rgc']
    hgc = Dgc['hgc']
    negc0 = Dgc['negc0']
    Fgc0 = Dgc['Fgc0']

    '''

    # GC component is nonzero only for abs(y) < rgc (currently)
    # so conservatively exit function for abs(y) > absymax = multiple of rgc:

    if abs(y) > absymax:
        return 0, 0

    rr = sqrt((x-xgc)**2. + (y-ygc)**2.)        #galactocentric radius
    zz = abs(z-zgc)                             #z-height

    if rr > rgc or abs(z-zgc) > hgc:
        return 0, 0
    else:
        arg = (rr/rgc)**2. + (zz/hgc)**2.
        if arg <= 1:
            ne_gc_out = negc0
            F_gc = Fgc0
        else: # need this to avoid errors when arg is not <=1
            ne_gc_out = 0
            F_gc = 0

    return ne_gc_out, F_gc


# ---------------------------------------------------------------------------

def ne_gc_vec(x, y, z, absymax=2*rgc):
    """Vectorized (array) version of ne_gc.  x, y, z are 1-D numpy arrays."""
    rr = np.sqrt((x - xgc)**2 + (y - ygc)**2)
    zz = np.abs(z - zgc)
    arg = (rr / rgc)**2 + (zz / hgc)**2
    mask = (np.abs(y) <= absymax) & (rr <= rgc) & (zz <= hgc) & (arg <= 1)
    negc_v = np.where(mask, negc0, 0.0)
    Fgc_v  = np.where(mask, Fgc0,  0.0)
    return negc_v, Fgc_v


# ---------------------------------------------------------------------------

def ne_outer_vec(x, y, z):
    """Vectorized (array) version of ne_outer.  x, y, z are 1-D numpy arrays."""
    rr = np.sqrt(x**2 + y**2)
    suncos = np.cos(pihalf * rsun / A1)
    g1 = np.where(rr > A1, 0.0, np.cos(pihalf * rr / A1) / suncos)
    z_arg = np.abs(z / h1)
    exp_2z = np.exp(-2.0 * z_arg)
    sech2_val = 4.0 * exp_2z / (1.0 + exp_2z)**2
    ne1 = (n1h1 / h1) * g1 * sech2_val
    return ne1, np.full_like(ne1, F1)


def ne_inner_vec(x, y, z):
    """Vectorized (array) version of ne_inner.  x, y, z are 1-D numpy arrays."""
    rr = np.sqrt(x**2 + y**2)
    rrarg = ((rr - A2) / 1.8)**2
    g2 = np.where(rrarg < 10.0, np.exp(-rrarg), 0.0)
    z_arg = np.abs(z / h2)
    exp_2z = np.exp(-2.0 * z_arg)
    sech2_val = 4.0 * exp_2z / (1.0 + exp_2z)**2
    ne2 = n2 * g2 * sech2_val
    return ne2, np.full_like(ne2, F2)
