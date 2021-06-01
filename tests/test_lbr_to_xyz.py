import pygedm
import numpy as np
from astropy.coordinates import Angle, Galactocentric
from astropy.units import Unit, Quantity
import astropy.units as u
import pytest



def test_xyz_to_lbr():
    x,y,z = pygedm.convert_lbr_to_xyz(0, 0, 0, method='ymw16')
    assert x == 0
    assert y == 8300 * u.pc
    assert z == 6 * u.pc

    x,y,z = pygedm.convert_lbr_to_xyz(0, 0, 0, method='ne2001')
    assert x == 0
    assert y == 8500 * u.pc
    assert z == 0

    # Who even knows where the centre of the Galaxy is? 
    # https://docs.astropy.org/en/stable/api/astropy.coordinates.galactocentric_frame_defaults.html
    from astropy.coordinates import galactocentric_frame_defaults
    _ = galactocentric_frame_defaults.set('v4.0') 
    x,y,z = pygedm.convert_lbr_to_xyz(0, 0, 0, method='astropy')
    assert np.isclose(x.to('pc').value, -1 * Galactocentric().galcen_distance.to('pc').value)
    assert y == 0
    assert np.isclose(z.to('pc').value, Galactocentric().z_sun.to('pc').value)

    x,y,z = pygedm.convert_lbr_to_xyz(Angle(0, unit='degree'), Angle(0, unit='rad'), 0 * u.pc, method='ymw16')
    assert x == 0
    assert y == 8300 * u.pc
    assert z == 6 * u.pc

if __name__ == "__main__":
    test_xyz_to_lbr()