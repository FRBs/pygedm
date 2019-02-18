#!/usr/bin/env python
"""
# PyYMW16.py -- python/C++ port of YMW16 C code
"""

import ymw16
from astropy import units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

MODE_IDS = {'gal': 1,
            'mc': 0,
            'igm': -1}
DATAPATH = './data'

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
    if isinstance(gl, Angle):
        gl = gl.to('deg').value
    if isinstance(gb, Angle):
        gb = gb.to('deg').value
    if isinstance(dm, Quantity):
        dm = dm.to('pc cm^(-3)').value
    mode_id = MODE_IDS.get(mode.lower().strip())
    ndir, vbs, txt = 1, 0, ''

    r = ymw16.dmdtau(gl, gb, dm, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
    dist = r['dist'] * u.pc
    tau_sc = r['tau_sc'] * u.s
    return (dist, tau_sc)

def dist_to_dm(gl, gb, dist, mode='gal'):
    """ Convert a DM to a distance

    Args:
        gl: galactic longitude (deg)
        gb: galactic latitude (deg)
        dist: distance to source (pc)
    """
    if isinstance(gl, Angle):
        gl = gl.to('deg').value
    if isinstance(gb, Angle):
        gb = gb.to('deg').value
    if isinstance(dist, Quantity):
        dist = dist.to('pc').value
    mode_id = MODE_IDS.get(mode.lower().strip())
    ndir, dm_host, vbs, txt = 2, 0, 0, ''

    r = ymw16.dmdtau(gl, gb, dist, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
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
    # ne = ne_crd(&x, &y, &z, &gl, &gb, &dist, ncrd, vbs, dirname, text);
    ndir, vbs, txt = 1, 0, ''
    ne = ymw16.ne_crd(x, y, z, 0, 0, 0, ndir, vbs, DATAPATH, txt)
    return ne

def calculate_electron_density_lbr(ga, gl, b):
    """ Calculate electron density at a point with Galactic coords (ga, gl)
        at a given distance in pc

    Args:
        (gl, gb, dist): Galactic lat/long in deg, dist in pc

    Returns:
        N_e: electron density in cm^-3
    """
    # ne = ne_crd(&x, &y, &z, &gl, &gb, &dist, ncrd, vbs, dirname, text);
    ndir, vbs, txt = 2, 0, ''
    ne = ymw16.ne_crd(0, 0, 0, ga, gl, dist, ndir, vbs, DATAPATH, txt)
    return ne

def test_dm_to_dist():
    a = dm_to_dist(204, -6.5, 200)
    b = dm_to_dist(Angle(204, unit='degree'), Angle(-6.5, unit='degree'), 200)
    c = dm_to_dist(204, -6.5, 200 * Unit('pc cm^-3'))

    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]

def test_dist_to_dm():
    a = dist_to_dm(204, -6.5, 200)
    b = dist_to_dm(Angle(204, unit='degree'), Angle(-6.5, unit='degree'), 200)
    c = dist_to_dm(204, -6.5, 200 * Unit('pc'))
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]

if __name__ == "__main__":
    test_dm_to_dist()
    test_dist_to_dm()
