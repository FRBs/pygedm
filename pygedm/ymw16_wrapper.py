#!/usr/bin/env python
"""
# PyYMW16.py -- python/C++ port of YMW16 C code
"""

import ymw16
from astropy import units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit
import os
from pkg_resources import resource_filename, get_distribution, DistributionNotFound

DATAPATH = os.path.dirname(resource_filename("pygedm", "spiral.txt"))


MODE_IDS = {'gal': 1,
            'mc': 0,
            'igm': -1}

def dm_to_dist(gl, gb, dm, dm_host=0, mode='gal'):
    """ Convert a DM to a distance

    Args:
        gl: galactic longitude (float in deg or astropy.Angle)
        gb: galactic latitude (float in deg or astropy.Angle)
        dm: dispersion measure (float or astropy.Quantity, , pc cm^-3)
        mode: Gal, MC, or IGM

    Returns:
        (dist, tau_sc) tuple
        dist: distance (pc, astropy.Quantity)
        tau_sc: scattering time scale (s, astropy.Quantity)
    """

    mode_id = MODE_IDS.get(mode.lower().strip())
    ndir, vbs, txt = 1, 0, ''

    r = ymw16.dmdtau(gl, gb, dm, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
    if mode == 'igm':
        r['dist'] *= u.Mpc
        r['tau_sc'] = r['tau_FRB']
    else:
        r['dist'] *= u.pc
    dist = r['dist']
    tau_sc = r['tau_sc'] * u.s
    return (dist, tau_sc)

def dist_to_dm(gl, gb, dist, mode='gal'):
    """ Convert a distance to a DM

    Args:
        gl: galactic longitude (deg)
        gb: galactic latitude (deg)
        dist: distance to source (pc) or if in mode IGM (Mpc)
    """
    mode_id = MODE_IDS.get(mode.lower().strip())
    ndir, dm_host, vbs, txt = 2, 0, 0, ''

    r = ymw16.dmdtau(gl, gb, dist, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
    if mode == 'igm':
        r['DM'] = r['DM_IGM']
        r['tau_sc'] = r['tau_FRB']
    dm = r['DM'] * u.pc / u.cm**3
    tau_sc = r['tau_sc'] * u.s
    return dm, tau_sc

def calculate_electron_density_xyz(x, y, z):
    """ Calculate electron density at a point with galactocentric coords (x, y, z)

    Args:
        (x, y, z): galactocentric coordinates in pc

    Returns:
        N_e: electron density in cm^-3
    """
    ndir, vbs, txt = 1, 0, ''
    ne = ymw16.ne_crd(x, y, z, 0, 0, 0, ndir, vbs, DATAPATH, txt)
    return ne / u.cm**3

def calculate_electron_density_lbr(gl, gb, dist):
    """ Calculate electron density at a point with Galactic coords (ga, gl)
        at a given distance in pc

    Args:
        (gl, gb, dist): Galactic lat/long in deg, dist in pc

    Returns:
        N_e: electron density in cm^-3
    """
    ndir, vbs, txt = 2, 0, ''
    ne = ymw16.ne_crd(0, 0, 0, gl, gb, dist, ndir, vbs, DATAPATH, txt)
    return ne  / u.cm**3
