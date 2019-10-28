#!/usr/bin/env python
"""
# PyGEDM.py -- Galactic Electron Density Models (NE2001 + YMW16) in python
"""

from . import ymw16_wrapper
from . import ne2001_wrapper
from astropy import units as u
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit
import os
from pkg_resources import resource_filename, get_distribution, DistributionNotFound


try:
    __version__ = get_distribution('pygedm').version
except DistributionNotFound: # pragma: no cover
    __version__ = '0.0.1 - manual install'


def dm_to_dist(gl, gb, dm, dm_host=0, mode='gal', method='ymw16'):
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

    if method.lower() == 'ymw16':
        return ymw16_wrapper.dm_to_dist(gl, gb, dm, dm_host=0, mode=mode)
    elif method.lower() == 'ne2001':
        if mode != 'gal':
            raise RuntimeError('NE2001 only supports Galactic (gal) mode.')
        return ne2001_wrapper.dm_to_dist(gl, gb, dm - dm_host)
    else:
        raise RuntimeError("Only ymw16 and ne2001 models supported.")


def dist_to_dm(gl, gb, dist, mode='gal', method='ymw16'):
    """ Convert a distance to a DM

    Args:
        gl: galactic longitude (deg)
        gb: galactic latitude (deg)
        dist: distance to source (pc) or if in mode IGM (Mpc) 
    """
    if isinstance(gl, Angle):
        gl = gl.to('deg').value
    if isinstance(gb, Angle):
        gb = gb.to('deg').value
    if isinstance(dist, Quantity):
        if mode == 'igm':
            dist = dist.to('Mpc').value
        else:
            dist = dist.to('pc').value
    if method.lower() == 'ymw16':
        return ymw16_wrapper.dist_to_dm(gl, gb, dist, mode=mode)
    elif method.lower() == 'ne2001':
        if mode != 'gal':
            raise RuntimeError('NE2001 only supports Galactic (gal) mode.')
        dist_kpc = dist / 1000.0
        return ne2001_wrapper.dist_to_dm(gl, gb, dist_kpc)
    else:
        raise RuntimeError("Only ymw16 and ne2001 models supported.")


def calculate_electron_density_xyz(x, y, z, method='ymw16'):
    """ Calculate electron density at a point with galactocentric coords (x, y, z)

    Args:
        (x, y, z): galactocentric coordinates in pc

    Returns:
        N_e: electron density in cm^-3
    """
    if isinstance(x, Quantity):
        x = x.to('pc').value
    if isinstance(y, Quantity):
        y = y.to('pc').value
    if isinstance(z, Quantity):
        z = z.to('pc').value
    if method.lower() == 'ymw16':
        return ymw16_wrapper.calculate_electron_density_xyz(x, y, z)
    elif method.lower() == 'ne2001':
        return ne2001_wrapper.calculate_electron_density_xyz(x/1e3, y/1e3, z/1e3)
    else:
        raise RuntimeError("Only ymw16 and ne2001 models supported.")


def calculate_electron_density_lbr(gl, gb, dist, method='ymw16'):
    """ Calculate electron density at a point with Galactic coords (ga, gl)
        at a given distance in pc

    Args:
        (gl, gb, dist): Galactic lat/long in deg, dist in pc

    Returns:
        N_e: electron density in cm^-3
    """
    if isinstance(gl, Angle):
        gl = gl.to('deg').value
    if isinstance(gb, Angle):
        gb = gb.to('deg').value
    if isinstance(dist, Quantity):
        dist = dist.to('pc').value
    if method.lower() == 'ymw16':
        return ymw16_wrapper.calculate_electron_density_lbr(gl, gb, dist)
    else:
        raise RuntimeError("Only ymw16 model supported.")
