import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

import pygedm


def test_dm_to_dist():
    """Test that astropy units / angles work with dm_to_dist"""
    a = pygedm.dm_to_dist(204, -6.5, 200)
    b = pygedm.dm_to_dist(Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200)
    c = pygedm.dm_to_dist(204, -6.5, 200 * Unit("pc cm^-3"))
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_dist_to_dm():
    """Test that astropy units / angles work with dist_to_dm"""
    a = pygedm.dist_to_dm(204, -6.5, 200)
    b = pygedm.dist_to_dm(Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200)
    c = pygedm.dist_to_dm(204, -6.5, 200 * Unit("pc"))
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_calculate_electron_density_xyz():
    pc = Unit("pc")
    a = pygedm.calculate_electron_density_xyz(1, 2, 3)
    b = pygedm.calculate_electron_density_xyz(1 * pc, 2, 3)
    c = pygedm.calculate_electron_density_xyz(1, 2 * pc, 3)
    d = pygedm.calculate_electron_density_xyz(1, 2, 3 * pc)
    assert a == b == c == d


def test_calculate_electron_density_lbr():
    pc = Unit("pc")
    a = pygedm.calculate_electron_density_lbr(1, 2, 3)
    b = pygedm.calculate_electron_density_lbr(Angle(1, unit="deg"), 2, 3)
    c = pygedm.calculate_electron_density_lbr(1, Angle(2, unit="deg"), 3)
    d = pygedm.calculate_electron_density_lbr(1, 2, 3 * pc)
    assert a == b == c == d


def test_basic():
    """Basic tests of YMW16 model

    Note: tested against online YMW16 interface
    http://www.atnf.csiro.au/research/pulsar/ymw16/index.php
    """

    a = pygedm.calculate_electron_density_xyz(1, 2, 3)
    assert np.isclose(a.value, 5.220655, atol=0.0001)

    a = pygedm.calculate_electron_density_lbr(0, 0, 4000)
    assert np.isclose(a.value, 0.388407, atol=0.0001)

    # FRB180301 value
    dm, tau = pygedm.dist_to_dm(204, -6.5, 25000)
    assert np.isclose(dm.value, 252.0501, atol=0.01)

    # Loop through distances and check round trip
    for dist in (10.0, 100.0, 1000.0):
        dm, tau = pygedm.dist_to_dm(0, 0, dist)
        dist_out, tau = pygedm.dm_to_dist(0, 0, dm.value)
        assert np.isclose(dist_out.value, dist, rtol=0.1)


def test_igm():
    """Test that IGM mode works as expected

    Note: tested against YMW16 code with:
    # CMD:    ./ ymw16 -d data -v IGM 204 -6.5 2000 100 1
    # OUTPUT: DM_Gal:  252.05 DM_MC:    0.00 DM_IGM: 1647.95 DM_Host:  100.00
    #         z:  2.311   Dist:  5336.4   log(tau_sc): -2.218
    """

    dist, tau = pygedm.dm_to_dist(204, -6.5, 2000, dm_host=100, mode="igm")
    assert np.isclose(dist.value, 5336.4, rtol=0.1)
    assert np.isclose(np.log10(tau.value), -2.218, rtol=0.1)

    dm, tau = pygedm.dist_to_dm(204, -6.5, 5336.4, mode="igm")
    dm_total = dm.value + 252.05 + 100  # Add galactic and host contribution
    assert np.isclose(dm_total, 2000, rtol=0.1)

    dm, tau = pygedm.dist_to_dm(204, -6.5, 5336.4 * u.Mpc, mode="igm")
    dm_total = dm.value + 252.05 + 100  # Add galactic and host contribution
    assert np.isclose(dm_total, 2000, rtol=0.1)


def test_magellanic_cloud():
    """Test that MC mode agrees with YMW16

    Note: tested against YMW16 code with:
    CMD:        ./ymw16 -d data -v MC 280.46 -32.88 50000 2
    OUTPUT:     MC: gl= 280.460 gb= -32.880 Dist=  50000.0
                dtest=50000.000000, nstep=10000.000000, dstep=5.000000
                DM_Gal:   58.03 DM_MC:   78.87 DM:  136.90 log(tau_sc): -5.399
    """
    dm, tau = pygedm.dist_to_dm(280.46, -32.88, 50000, mode="mc")
    assert np.isclose(dm.value, 136.90)
    assert np.isclose(np.log10(tau.value), -5.399, rtol=0.1)


if __name__ == "__main__":
    test_basic()
    test_dm_to_dist()
    test_dist_to_dm()
    test_calculate_electron_density_xyz()
    test_igm()
    test_magellanic_cloud()
