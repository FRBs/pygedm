import pygedm
import numpy as np
from astropy.coordinates import Angle
from astropy.units import Unit, Quantity
import astropy.units as u
import pytest

def test_dm_to_dist():
    """ Test that astropy units / angles work with dm_to_dist """
    a = pygedm.dm_to_dist(204, -6.5, 200, method='ne2001')
    b = pygedm.dm_to_dist(Angle(204, unit='degree'), Angle(-6.5, unit='degree'), 200, method='ne2001')
    c = pygedm.dm_to_dist(204, -6.5, 200 * Unit('pc cm^-3'), method='ne2001')
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]

def test_dist_to_dm():
    """ Test that astropy units / angles work with dist_to_dm """
    a = pygedm.dist_to_dm(204, -6.5, 200, method='ne2001')
    b = pygedm.dist_to_dm(Angle(204, unit='degree'), Angle(-6.5, unit='degree'), 200, method='ne2001')
    c = pygedm.dist_to_dm(204, -6.5, 200 * Unit('pc'), method='ne2001')
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]

def test_calculate_electron_density_xyz():
    pc = Unit('pc')
    a = pygedm.calculate_electron_density_xyz(1, 2, 3, method='ne2001')
    b = pygedm.calculate_electron_density_xyz(1 * pc, 2, 3, method='ne2001')
    c = pygedm.calculate_electron_density_xyz(1, 2 * pc, 3, method='ne2001')
    d = pygedm.calculate_electron_density_xyz(1, 2, 3 * pc, method='ne2001')
    assert a == b == c == d

def test_calculate_electron_density_lbr():
    with pytest.raises(RuntimeError):
        pygedm.calculate_electron_density_lbr(100, 10, 100, method='ne2001')

def test_basic():
    """ Basic tests of YMW16 model

    Note: tested against online NE2001 interface
    https://www.nrl.navy.mil/rsd/RORF/ne2001/
    """

    # No access to actual model via web interface
    #a = pygedm.calculate_electron_density_xyz(1, 2, 3)
    #assert np.isclose(a.value, 5.220655, atol=0.0001)

    # FRB180301 value
    dm, tau = pygedm.dist_to_dm(204, -6.5, 25*u.kpc, method='ne2001')
    assert np.isclose(dm.value, 150.80, atol=0.01)

    # Loop through distances and check round trip
    for dist in (10.*u.pc, 100.*u.pc, 1000.*u.pc):
        dm, tau = pygedm.dist_to_dm(0, 0, dist, method='ne2001')
        dist_out, tau = pygedm.dm_to_dist(0, 0, dm, method='ne2001')
        print(dist, dm, dist_out)
        assert np.isclose(dist_out.to('pc').value, dist.to('pc').value, rtol=0.1)


def test_igm():
    """ Test that IGM mode FAILS as expected """
    with pytest.raises(RuntimeError):
         pygedm.dm_to_dist(100, 10, 100, mode='igm', method='ne2001')


def test_dm_wrapper():
    """ Run test against known values
    ## Test data from https://www.nrl.navy.mil/rsd/RORF/ne2001_src/
    """
    test_data = {
        'l': [0, 2, 97.5, ],
        'b': [0, 7.5, 85.2, ],
        'dm': [10, 20, 11.1],
        'dist': [0.461, 0.781, 0.907]
    }

    for ii in range(len(test_data['l'])):
        dist, smtau = pygedm.ne2001_wrapper.dm_to_dist(test_data['l'][ii], test_data['b'][ii], test_data['dm'][ii])
        assert np.allclose(dist.to('kpc').value, test_data['dist'][ii], atol=2)

        dm, smtau = pygedm.ne2001_wrapper.dist_to_dm(test_data['l'][ii], test_data['b'][ii], test_data['dist'][ii])
        assert np.allclose(dm.value, test_data['dm'][ii], atol=2)

def test_zero_dm():
    """ Check that zero DM doesn't cause timeout bug

    Fortran code hangs if DM or D == 0
    """
    dist, smtau = pygedm.ne2001_wrapper.dm_to_dist(0, 0, 0)
    dm, smtau   = pygedm.ne2001_wrapper.dist_to_dm(0, 0, 0)
    assert dist.value == dm.value == 0

if __name__ == "__main__":
    test_basic()
    test_dm_to_dist()
    test_dist_to_dm()
    test_calculate_electron_density_xyz()
    test_igm()
    test_dm_wrapper()
