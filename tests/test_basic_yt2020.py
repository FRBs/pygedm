import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

import pygedm
from pygedm import yt2020


def test_basic():
    """
    Test data from command line output of DM_halo_YT2020_numerical.py
    (courtesy Shotaro Yamasaki)

    python DM_halo_YT2020_numerical.py 0 0
        Upsilon = 2.604689
        n_0_sphe = 0.000372 [cm^{-3}]
        DM_halo_disk(0,0) = 220.362398 [pc cm^{-3}]
        DM_halo_sphe(0,0) = 24.857097 [pc cm^{-3}]
        DM_halo(0,0) = 245.219494 [pc cm^{-3}]

    python DM_halo_YT2020_numerical.py 10 45
        Upsilon = 2.604689
        n_0_sphe = 0.000372 [cm^{-3}]
        DM_halo_disk(0.2,0.79) = 24.843628 [pc cm^{-3}]
        DM_halo_sphe(0.2,0.8) = 23.325799 [pc cm^{-3}]
        DM_halo(0.2,0.8) = 48.169427 [pc cm^{-3}]

    python DM_halo_YT2020_numerical.py -100 32
        Upsilon = 2.604689
        n_0_sphe = 0.000372 [cm^{-3}]
        DM_halo_disk(-1.7,0.56) = 14.876829 [pc cm^{-3}]
        DM_halo_sphe(-1.7,0.6) = 21.177893 [pc cm^{-3}]
        DM_halo(-1.7,0.6) = 36.054722 [pc cm^{-3}]

    python DM_halo_YT2020_numerical.py -100.3 32.1
        Upsilon = 2.604689
        n_0_sphe = 0.000372 [cm^{-3}]
        DM_halo_disk(-1.8,0.56) = 14.819495 [pc cm^{-3}]
        DM_halo_sphe(-1.8,0.6) = 21.169654 [pc cm^{-3}]
        DM_halo(-1.8,0.6) = 35.989149 [pc cm^{-3}
    """
    test_list = [
        (0, 0, 220.362398, 24.857097, 245.219494),
        (10, 45, 24.843628, 23.325799, 48.169427),
        (-100, 32, 14.876829, 21.177893, 36.054722),
        (-100.3, 32.1, 14.819495, 21.169654, 35.989149),
    ]

    for line in test_list:
        l, b, dm_disk, dm_sphe, dm_tot = line
        dm_disk_out = yt2020.calculate_halo_dm(l, b, component="disk")
        dm_sphe_out = yt2020.calculate_halo_dm(l, b, component="spherical")
        dm_tot_out = yt2020.calculate_halo_dm(l, b, component="both")

        assert np.isclose(dm_disk, dm_disk_out.value)
        assert np.isclose(dm_sphe, dm_sphe_out.value)
        assert np.isclose(dm_tot, dm_tot_out.value)

    # Tests for higher-level wrapper in __init__.py
    with pytest.raises(RuntimeError):
        pygedm.calculate_halo_dm(l, b, component="nonexistent_component")

    with pytest.raises(RuntimeError):
        pygedm.calculate_halo_dm(l, b, method="nonexistent_component")

    for line in test_list:
        l, b, dm_disk, dm_sphe, dm_tot = line
        dm_disk_out = pygedm.calculate_halo_dm(l, b, method="yt2020", component="disk")
        dm_sphe_out = pygedm.calculate_halo_dm(
            l, b, method="yt2020", component="spherical"
        )
        dm_tot_out = pygedm.calculate_halo_dm(l, b, method="yt2020", component="both")

        assert np.isclose(dm_disk, dm_disk_out.value)
        assert np.isclose(dm_sphe, dm_sphe_out.value)
        assert np.isclose(dm_tot, dm_tot_out.value)

    lbd_list = ((0, 0, 250.12), (0, 30, 68.07), (15, 0, 204.9))
    for line in lbd_list:
        l, b, dm = line
        dm_out = pygedm.calculate_halo_dm(l, b, method="yt2020_analytic")
        assert np.isclose(dm, dm_out.value, atol=0.1)

    with pytest.raises(RuntimeError):
        pygedm.calculate_halo_dm(l, b, method="yt2020_analytic", component="disk")


def compute_dm_halo_analytic():
    """Test the analytical version (runs faster)"""
    lbd_list = ((0, 0, 250.12), (0, 30, 68.07), (15, 0, 204.9))
    for line in lbd_list:
        l, b, dm = line
        dm_out = yt2020.calculate_halo_dm_analytic(l, b)
        assert np.isclose(dm, dm_out.value, atol=0.1)


if __name__ == "__main__":
    test_basic()
    compute_dm_halo_analytic()
