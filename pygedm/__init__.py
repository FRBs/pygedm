#!/usr/bin/env python
"""
# PyGEDM.py -- Galactic Electron Density Models (NE2001 + YMW16) in python
"""

from . import ymw16_wrapper
from . import ne2001_wrapper
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, Galactocentric
from astropy.units import Quantity, Unit
import numpy as np
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
    elif method.lower() == 'ne2001':
        x, y, z = convert_lbr_to_xyz(gl, gb, dist, method='ne2001')
        return ne2001_wrapper.calculate_electron_density_xyz(x.to('kpc').value, y.to('kpc').value, z.to('kpc').value)
    else:
        raise RuntimeError("Only ymw16 model supported.")

def convert_lbr_to_xyz(gl, gb, dist, method='ymw16'):
    """ Convert Galactic (l,b,r) coords to Galactocentric (x,y,z) coords

    Args:
        (gl, gb, dist): Galactic lat/long in deg, dist in pc.
                        use of astropy Angle/Quantities also OK
        method (str): one of 'ymw16', 'ne2001', or 'astropy'

    Notes:
        For transform, the Sun is located at (x=0, y=R_sun, z=z_sun)
        YMW16  uses R_sun of 8300 pc and z_sun of 6.0 pc
        NE2001 uses R_sun of 8500 pc and z_sun of 0.0 pc
        Both of these do a basic spherical to cartesian conversion.

        astropy does a much more complicated conversion, see
        https://astropy.readthedocs.io/en/latest/coordinates/galactocentric.html
        This is the 'proper' coordinate system, but note that it is NOT COMPATIBLE
        WITH NE2001 OR YMW16! (!SEE EXAMPLE OUTPUT BELOW!)

    Example output:
        pygedm.convert_lbr_to_xyz(0, 0, 0, method='ymw16')
        (<Quantity 0. pc>, <Quantity 8300. pc>, <Quantity 6. pc>)

        pygedm.convert_lbr_to_xyz(0, 0, 0, method='ne2001')
        (<Quantity 0. pc>, <Quantity 8500. pc>, <Quantity 0. pc>)

        pygedm.convert_lbr_to_xyz(0, 0, 0, method='astropy')
        (<Quantity -8499.95711754 pc>, <Quantity 0. pc>, <Quantity 27. pc>)
    """

    if isinstance(gl, Angle):
        gl = gl.to('rad').value
    else:
        gl = np.deg2rad(gl)
    if isinstance(gb, Angle):
        gb = gb.to('rad').value
    else:
        gb = np.deg2rad(gb)
    if isinstance(dist, Quantity):
        dist = dist.to('pc').value
    dist = dist * u.pc   # Make sure distance is in units of pc
    if method == 'astropy':
        lbr = SkyCoord(gl, gb, distance=dist, frame='galactic', unit=('degree', 'degree', 'pc'))
        xyz = lbr.transform_to(Galactocentric(galcen_distance=8.5 * u.kpc))
        return xyz.x, xyz.y, xyz.z
    elif method.lower() == 'ymw16':
        Rsun = 8300 * u.pc
        Zsun = 6.0 * u.pc
        x = dist * np.sin(gl) * np.cos(gb)
        y = Rsun - dist * np.cos(gl) * np.cos(gb)
        z = Zsun + dist * np.sin(gb)
        return x, y, z
    elif method.lower() == 'ne2001':
        Rsun = 8500 * u.pc
        x = dist * np.sin(gl) * np.cos(gb)
        y = Rsun - dist * np.cos(gl) * np.cos(gb)
        z = dist * np.sin(gb)
        return x, y, z
    else:
        raise RuntimeError("method must be astropy, ne2001, or ymw16")