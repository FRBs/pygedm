import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

import pygedm


def test_tau_sc_nu():
    """Test that tau_sc changes with nu"""
    dm, tau_sc = pygedm.dist_to_dm(0, 0, 100, method="ymw16", nu=1)
    dm_, tau_sc_ = pygedm.dist_to_dm(0, 0, 100, method="ymw16", nu=1000 * u.MHz)
    assert dm == dm_
    assert tau_sc == tau_sc_

    dm, tau_sc = pygedm.dist_to_dm(0, 0, 100, method="ne2001", nu=1)
    dm_, tau_sc_ = pygedm.dist_to_dm(0, 0, 100, method="ne2001", nu=1000 * u.MHz)
    assert dm == dm_
    assert tau_sc == tau_sc_

    dist, tau_sc = pygedm.dist_to_dm(0, 0, 100, method="ymw16", nu=1)
    dist_, tau_sc_ = pygedm.dist_to_dm(0, 0, 100, method="ymw16", nu=1000 * u.MHz)
    assert dist == dist_
    assert tau_sc == tau_sc_

    dist, tau_sc = pygedm.dist_to_dm(0, 0, 100, method="ne2001", nu=1)
    dist_, tau_sc_ = pygedm.dist_to_dm(0, 0, 100, method="ne2001", nu=1000 * u.MHz)
    assert dist == dist_
    assert tau_sc == tau_sc_

    dm, tau_sc_1GHz = pygedm.dm_to_dist(0, 0, 1000, method="ymw16", nu=1.0)
    dm, tau_sc_100MHz = pygedm.dm_to_dist(0, 0, 1000, method="ymw16", nu=0.1)
    assert np.isclose(tau_sc_1GHz.value, 0.31681767)
    assert np.isclose(tau_sc_100MHz.value, 3168.17671061)

    assert np.isclose((0.1 / 1.0) ** (-4) * tau_sc_1GHz.value, tau_sc_100MHz.value)

    dm, tau_sc_1GHz = pygedm.dm_to_dist(0, 0, 1000, method="ne2001", nu=1.0)
    dm, tau_sc_100MHz = pygedm.dm_to_dist(0, 0, 1000, method="ne2001", nu=0.1)
    assert np.isclose(tau_sc_1GHz.value, 198.55289306)
    assert np.isclose(tau_sc_100MHz.value, 4987423.18015693)

if __name__ == "__main__":
    test_tau_sc_nu()
