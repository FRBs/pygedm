"""
Wrapper from mwprop / ne2025

References:
    [1] `S.K. Ocker, J.M. Cordes (2026) <https://ui.adsabs.harvard.edu/abs/2026arXiv260211838O/abstract>`_,
    *NE2025: An Updated Electron Density Model for the Galactic Interstellar Medium*,
    https://arxiv.org/abs/2602.11838
"""

import os
import numpy as np

from astropy import units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

from mwprop.nemod.NE2025 import ne2025
from mwprop.nemod.NE2001 import ne2001
from mwprop.nemod.density import density_2001


def dm_to_dist(gl, gb, dm, nu=1.0, model='ne2025', full_output=False):
    """Convert a DM to a distance

    Args:
        gl (float in deg): galactic longitude
        gb (float in deg): galactic latitude
        dm (float in pc/cm3): dispersion measure (pc cm^-3)
        nu (float in GHz): Observing frequency
        model (str): One of 'ne2001' or 'ne2025'

    Notes:
        In this wrapper, 'ne2001' == 'ne2001p'.

    Returns:
        dist (astropy.Quantity), tau_sc (astropy.Quantity): distance (pc) and scattering time scale (s)
    """
    if np.isclose(dm, 0):  # WAR Catch spline failure
        return 0.0 * u.pc, 0.0 * u.s
    else:
        if model.lower() == 'ne2025':
            _Dk,Dv,_Du,_Dd = ne2025(gl, gb, dm, ndir=1, dmd_only=False, classic=False)
        else:
            _Dk,Dv,_Du,_Dd = ne2001(gl, gb, dm, ndir=1, dmd_only=False, classic=False)
        dist, sm_tau = Dv['DIST'], Dv['SMtau']

        # Convert a scattering measure (SM) to scattering timescale at given frequency.
        tau_sc = 1.0 * (sm_tau / 292.0) ** 1.2 * dist * nu ** (-4.4)

        if not full_output:
            return (float(dist) * u.kpc).to("pc"), tau_sc * u.s
        else:
            Dv["tau_sc"] = tau_sc
            return Dv

def dist_to_dm(gl, gb, dist, nu=1.0, model='ne2025', full_output=False):
    """Convert a DM to a distance

    Args:
        gl (float in deg): galactic longitude
        gb (float in deg): galactic latitude
        dist (float in pc/cm3): dispersion measure (pc cm^-3)
        nu (float in GHz): Observing frequency
        model (str): One of 'ne2001' or 'ne2025'

    Notes:
        In this wrapper, 'ne2001' == 'ne2001p'.

    Returns:
        dist (kpc), tau_sc (astropy.Quantity): distance (pc) and scattering time scale (s)
    """
    if np.isclose(dist, 0):  # WAR Catch spline failure
        return 0.0 * u.pc, 0.0 * u.s
    else:
        if model.lower() == 'ne2025':
            _Dk,Dv,_Du,_Dd = ne2025(gl, gb, dist, ndir=-1, dmd_only=False, classic=False)
        else:
            _Dk,Dv,_Du,_Dd = ne2001(gl, gb, dist, ndir=-1, dmd_only=False, classic=False)


        dm, sm_tau = Dv['DM'], Dv['SMtau']

        # Convert a scattering measure (SM) to scattering timescale at given frequency.
        tau_sc = 1.0 * (sm_tau / 292.0) ** 1.2 * dist * nu ** (-4.4)

        if not full_output:
            return float(dm) * u.pc / u.cm**3, tau_sc * u.s
        else:
            Dv["tau_sc"] = tau_sc
            return Dv


def calculate_electron_density_xyz(x, y, z):
    """Compute electron density at Galactocentric X, Y, Z coordinates

    x,y,z are Galactocentric Cartesian coordinates, measured in kpc (NOT pc!)
    with the axes parallel to (l, b) = (90, 0), (180, 0), and (0, 90) degrees

    Args:
        x (float): Galactocentric coordinates in kpc
        y (float): Galactocentric coordinates in kpc
        z (float): Galactocentric coordinates in kpc

    Returns:
        ne_out (astropy.Quantity): Electron density in cm-3
    """
    density_list = density_2001(x, y, z)

    keys = ('ne1','ne2','nea','negc','nelism','necN','nevN',
            'F1', 'F2', 'Fa', 'Fgc', 'Flism', 'FcN', 'FvN',
            'whicharm', 'wlism', 'wldr', 'wlhb', 'wlsb', 'wloopI',
            'hitclump', 'hitvoid', 'wvoid')

    dd = dict(zip(keys, density_list))
    ne = dd['ne1'] + dd['ne2'] + dd['nea'] + dd['negc'] + dd['nelism'] + dd['necN'] + dd['nevN']
    ne = float(ne)

    return ne / u.cm**3
