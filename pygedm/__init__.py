#!/usr/bin/env python
"""
# PyGEDM.py -- Galactic Electron Density Models (NE2001 + YMW16) in python
"""

from . import ymw16_wrapper
from . import ne2001_wrapper
from . import yt2020
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, Galactocentric
from astropy.units import Quantity, Unit
import numpy as np
from pkg_resources import resource_filename, get_distribution, DistributionNotFound

try:
    import healpy as hp
    HAS_HEALPIX = True
except ImportError:  # pragma: no cover
    HAS_HEALPIX = False

try:
    __version__ = get_distribution('pygedm').version
except DistributionNotFound: # pragma: no cover1
    __version__ = '0.0.1 - manual install'


def _gl_gb_convert(gl, gb, unit='rad'):
    """ Convert gl, gb astropy.Angle to floats/arrays with correct units

    If not an astropy quantity, returns value unchanged.
    """
    if isinstance(gl, Angle) or isinstance(gl, Quantity):
        gl = gl.to(unit).value
    if isinstance(gb, Angle) or isinstance(gb, Quantity):
        gb = gb.to(unit).value
    return gl, gb


def _unit_convert(q, unit):
    """ Convert astropy.unit into float/array with given unit

    If not an astropy quantity, returns value unchanged.
    """
    if isinstance(q, Quantity):
        q = q.to(unit).value
    return q

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
    gl, gb = _gl_gb_convert(gl, gb, 'deg')
    dm = _unit_convert(dm, 'pc cm^(-3)')

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
    gl, gb = _gl_gb_convert(gl, gb, 'deg')

    if mode == 'igm':
        dist = _unit_convert(dist, 'Mpc')
    else:
        dist = _unit_convert(dist, 'pc')

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
    x = _unit_convert(x, 'pc')
    y = _unit_convert(y, 'pc')
    z = _unit_convert(z, 'pc')

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
    gl, gb = _gl_gb_convert(gl, gb, 'deg')
    dist = _unit_convert(dist, 'pc')

    if method.lower() == 'ymw16':
        return ymw16_wrapper.calculate_electron_density_lbr(gl, gb, dist)
    elif method.lower() == 'ne2001':
        x, y, z = convert_lbr_to_xyz(gl, gb, dist, method='ne2001')
        return ne2001_wrapper.calculate_electron_density_xyz(x.to('kpc').value, y.to('kpc').value, z.to('kpc').value)
    else:
        raise RuntimeError("Only ymw16 and ne2001 models supported.")


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

    gl, gb = _gl_gb_convert(gl, gb, 'deg')
    gl = np.deg2rad(gl)
    gb = np.deg2rad(gb)

    dist = _unit_convert(dist, 'pc') # Unit conversion and strip quantity
    dist = dist * u.pc               # Make sure distance is in units of pc
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


def generate_healpix_dm_map(dist=1, nside=64, method='ymw16'):
    """ Generate an all-sky healpix map for a given distance.

    Args:
        dist (float or Quantity): Distance to integrate EDM out to. 30 kpc will go to edge
        nside (int): The NSIDE parameter for the healpix map (power of 2, larger => higher resolution)
        method (str): one of 'ymw16', 'ne2001', 'yt2020' or 'yt2020_analytic'

    Notes:
        This method takes a fair amount of time to run -- tens of seconds for NSIDE=32.
        YT2020 method is even slower, consider using yt2020_analytic

    Returns:
         healpix map as a numpy array (1D), which can be viewed using the healpy.mollview() method
    """
    if HAS_HEALPIX:
        pix_id = np.arange(hp.nside2npix(nside))
        gl, gb = hp.pix2ang(nside, pix_id, lonlat=True)
        d = np.zeros_like(pix_id, dtype='float32')
        if 'yt2020' not in method:
            for ii in pix_id:
                dm, tau_sc = dist_to_dm(gl[ii], gb[ii], dist, method=method)
                d[ii] = dm.value

        elif 'yt2020' in method:
            # gl is in (0, 360), wrap to (-180, 180) as reqd for yt2020 function
            # Then call yt2020.compute_halo_dm function directly
            gl_u = np.copy(gl)
            gl_u[gl >= 180] -= 360
            if method == 'yt2020':
                for ii in pix_id:
                    dm_halo = yt2020.calculate_halo_dm(gl_u[ii], gb[ii])
                    d[ii] = dm_halo.value
            else:
                for ii in pix_id:
                    dm_halo = yt2020.calculate_halo_dm_analytic(gl_u[ii], gb[ii])
                    d[ii] = dm_halo.value
        return d
    else:
        raise RuntimeError("Healpy installation not found.")


def calculate_halo_dm(gl, gb, method='yt2020', component='both'):
    """ Compute halo DM
    Args:
        (gl, gb): Galactic lat/long in degrees
                        use of astropy Angle/Quantities also OK
        method (str): one of 'yt2020' (only YT2020 supported currently)
        component (str): Compute 'spherical' component of halo,
                         'disk' component, or 'both' components.
    Returns:
        DM (float): Dispersion measure in [pc/cm^3]

    """
    gl, gb = _gl_gb_convert(gl, gb, 'deg')
    gl = Angle(gl, unit='deg').wrap_at('180d').value  ## Ensure in correct range

    if method == 'yt2020':
        return yt2020.calculate_halo_dm(gl, gb, component)
    elif method == 'yt2020_analytic':
        if component != 'both':
            raise RuntimeError("Analytic formula requires component='both'")
        return yt2020.calculate_halo_dm_analytic(gl, gb)
    else:
        raise RuntimeError("Only YT2020 currently supported.")