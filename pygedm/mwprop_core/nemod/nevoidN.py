# mwprop v2.0 Jan 2026

'''
Python version of subroutine nevoidN.f in NE2001 Fortran code

Returns electron density nevN and fluctuation parameter FvN
at position designated by l,b,d,x,y,z c for a set of
voids with parameters read in from file  nevoidN.dat

input:
    x,y,z   coordinates (kpc)  (as in TC93)

output:
    nevN    electron density in void at (x,y,z)
    FvN fluctuation parameter
    hitvoid =   0:     no void hit
                j>0:   j-th void hit
    wvoid = 0,1:     void weight

parameters:
    lv  = galactic longitude of void center
    bv  = galactic latitude of void center
    dv  = distance from Sun of void center
    (xv,yv,zv) = void center location (calculated)
        nev = internal peak electron density
        Fv      = void fluctuation parameter
    aav = void major axis at 1/e
    bbv = void minor axis at 1/e
    ccv = void minor axis at 1/e
    thvy    = rotation axis of void about y axis
    thvz    = rotation axis of void about z axis
    edgev   = 0 => use exponential rolloff out to 5rc
                 1 => uniform and truncated at 1/e
Version history:

01/20/20 Stella Koch Ocker
    * initial conversion f77 --> python
01/23/20 -- now reads input parameters from dictionary program ne2001p_input
02/08/20 -- JMC
    * imports config_ne2001p for  model setup
    * rsun = 8.5 commented out; now set in setup file
02/09/20 -- JMC
    * comment out hitvoidflag; not used even in Fortran code
02/10/20 -- JMC
    * definition of void arrays moved to config_ne2001p.py
'''

from pygedm.mwprop_core.nemod.config_nemod import *
import numpy as np
from pygedm.mwprop_core.nemod.numba_compat import njit, HAS_NUMBA

@njit
def _nevoidN_jit(x, y, z, nvoids, xv, yv, zv, nev, Fv, aav, bbv, ccv,  # pragma: no cover
                 edgev, cc12, s2, cs21, cs12, c2, ss12, s1, c1):
    """JIT-compiled core loop for void calculation"""
    nevN = 0.
    FvN = 0.
    hitvoid = 0
    wvoid = 0

    for j in range(nvoids):
        dx = x - xv[j]
        dy = y - yv[j]
        dz = z - zv[j]
        q = ((cc12[j]*dx + s2[j]*dy + cs21[j]*dz)**2. / aav[j]**2. +
             (-cs12[j]*dx + c2[j]*dy - ss12[j]*dz)**2. / bbv[j]**2. +
             (-s1[j]*dx + c1[j]*dz)**2. / ccv[j]**2.)
        if edgev[j] == 0. and q < 3.:
            nevN = nev[j] * np.exp(-q)
            FvN = Fv[j]
            hitvoid = j+1
        if edgev[j] == 1. and q <= 1.:
            nevN = nev[j]
            FvN = Fv[j]
            hitvoid = j+1

    if hitvoid != 0:
        wvoid = 1

    return nevN, FvN, hitvoid, wvoid

def nevoidN(x,y,z):

    nevN = 0.
    FvN = 0.
    hitvoid = 0
    wvoid = 0

    if nvoids == 0:
        return nevN, FvN, hitvoid, wvoid

    # Call JIT-compiled core loop
    return _nevoidN_jit(x, y, z, nvoids, xv, yv, zv, nev, Fv, aav, bbv, ccv,
                        edgev, cc12, s2, cs21, cs12, c2, ss12, s1, c1)


# ---------------------------------------------------------------------------

def nevoidN_vec(x, y, z):
    """
    Vectorized (array) version of nevoidN.

    x, y, z are 1-D numpy arrays of positions.
    Loops over the (few) voids; vectorizes over positions.

    Returns (nevN_v, FvN_v, hitvoid_v, wvoid_v) — all 1-D arrays.
    """
    n = len(x)
    nevN_v    = np.zeros(n)
    FvN_v     = np.zeros(n)
    hitvoid_v = np.zeros(n, dtype=int)

    if nvoids == 0:
        return nevN_v, FvN_v, hitvoid_v, np.zeros(n)

    for j in range(nvoids):
        dx = x - xv[j]
        dy = y - yv[j]
        dz = z - zv[j]
        q = ((cc12[j]*dx  + s2[j]*dy   + cs21[j]*dz)**2 / aav[j]**2
           + (-cs12[j]*dx + c2[j]*dy   - ss12[j]*dz)**2 / bbv[j]**2
           + (-s1[j]*dx                +  c1[j]*dz)**2   / ccv[j]**2)
        if edgev[j] == 0.:              # Gaussian void
            mask = q < 3.
            nevN_v    = np.where(mask, nev[j] * np.exp(-q), nevN_v)
            FvN_v     = np.where(mask, Fv[j],  FvN_v)
            hitvoid_v = np.where(mask, j + 1,  hitvoid_v)
        elif edgev[j] == 1.:            # hard-edge void
            mask = q <= 1.
            nevN_v    = np.where(mask, nev[j], nevN_v)
            FvN_v     = np.where(mask, Fv[j],  FvN_v)
            hitvoid_v = np.where(mask, j + 1,  hitvoid_v)

    wvoid_v = (hitvoid_v != 0).astype(float)
    return nevN_v, FvN_v, hitvoid_v, wvoid_v
