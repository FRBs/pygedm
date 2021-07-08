"""
Python wrapper for NE2001 model code.



References:

    [1] `Cordes, J. M., & Lazio, T. J. W. (2002), <https://ui.adsabs.harvard.edu/abs/2002astro.ph..7156C/abstract>`_
    *NE2001.I. A New Model for the Galactic Distribution of Free Electrons and its Fluctuations*,
    arXiv e-prints, astro-ph/0207156.

    [2] `Cordes, J. M., & Lazio, T. J. W. (2003), <https://ui.adsabs.harvard.edu/abs/2003astro.ph..1598C/abstract>`_
    *NE2001. II. Using Radio Propagation Data to Construct a Model for the Galactic Distribution of Free Electrons*,
    arXiv e-prints, astro-ph/0301598.
"""

import numpy as np
import os
from functools import wraps
from astropy import units as u
import ne21c

DATA_PATH = os.path.dirname(os.path.abspath(__file__))


def run_from_pkgdir(f):
    """ Decorator function to chdir() into package directory when running

    NE2001 code doesn't know the relative path to its data files.
    This wraps the function call, changing into the right directory
    first, calling it, then changing back to original directory.
    """
    @wraps(f)
    def wrapped(*args, **kwargs):
        cwdpath = os.getcwd()
        try:
            os.chdir(DATA_PATH)
            r = f(*args, **kwargs)
            return r
        except:   # pragma: no cover
            raise
        finally:
            os.chdir(cwdpath)
    return wrapped


def TAUISS(d, sm, nu):
    """ Convert a scattering measure (SM) to scattering timescale at given frequency.

    Direct port from FORTRAN code scattering98.f

    Args:
        d (float): Distance in kpc
        sm (float): Scattering measure (kpc m^{-20/3})
        nu (float): Radio frequency in GHz

    Returns:
        tauiss (float): pulse broadening time (ms)

    Fortran equiv::

              REAL FUNCTION TAUISS(d, sm, nu)
        c
        c calculates the pulse broadening time in ms
        c from distance, scattering measure, and radio frequency
        c
        c input:      d = pulsar distance       (kpc)
        c            sm = scattering measure    (kpc m^{-20/3})
        c            nu = radio frequency       (GHz)
        c output: tausis = pulse broadening time (ms)
        c
              implicit none
              real d, sm,  nu
              tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
              end
    """
    tauiss = 1. * (sm / 292.)**1.2 * d * nu**(-4.4)
    return tauiss


@run_from_pkgdir
def dm_to_dist(l, b, dm, nu=1.0):
    """ Convert DM to distance and compute scattering timescale

    Args:
        l (float): galactic longitude in degrees
        b (float): galactic latitude in degrees
        dm (floa): Dispersion measure
        nu (float in GHz or astropy.Quantity): observing frequency (GHz)

    Returns:
        dist (astropy.Quantity), tau_sc (astropy.Quantity): Distance (pc), scattering timescale at 1 GHz (s)
    """
    l_rad = np.deg2rad(l)
    b_rad = np.deg2rad(b)

    if np.isclose(dm, 0):  # WAR Catch infinite timeout
        return 0.0 * u.pc, 0.0 * u.s
    else:
        d = ne21c.dm_to_dist(l_rad, b_rad, dm)

    tau_sc = TAUISS(float(d['dist']), d['smtau'], nu=nu)
    return (float(d['dist']) * u.kpc).to('pc'), tau_sc * u.s


@run_from_pkgdir
def dist_to_dm(l, b, dist, nu=1.0):
    """ Convert distance to DM and compute scattering timescale

    Args:
        l (float): galactic longitude in degrees
        b (float): galactic latitude in degrees
        dist (float): Distance in kpc
        nu (float in GHz or astropy.Quantity): observing frequency (GHz)

    Returns:
        dm (astropy.Quantity), tau_sc (astropy.Quantity): Dispersion measure (pc / cm3), scattering timescale at 1 GHz (s)
    """
    l_rad = np.deg2rad(l)
    b_rad = np.deg2rad(b)

    if np.isclose(dist, 0):       # Catch infinite timeout bug
        return 0.0 * u.pc / u.cm**3, 0.0 * u.s
    else:
        d = ne21c.dist_to_dm(l_rad, b_rad, dist)

    tau_sc = TAUISS(float(dist), d['smtau'], nu=nu)
    return float(d['dm']) * u.pc / u.cm**3, tau_sc * u.s


@run_from_pkgdir
def calculate_electron_density_xyz(x, y, z):
    """ Compute electron density at Galactocentric X, Y, Z coordinates

    x,y,z are Galactocentric Cartesian coordinates, measured in kpc (NOT pc!)
    with the axes parallel to (l, b) = (90, 0), (180, 0), and (0, 90) degrees

    Args:
        x (float): Galactocentric coordinates in kpc
        y (float): Galactocentric coordinates in kpc
        z (float): Galactocentric coordinates in kpc

    Returns:
        ne_out (astropy.Quantity): Electron density in cm-3
    """
    ne_out = ne21c.density_xyz(x, y, z)
    return ne_out['ne'] / u.cm**3
