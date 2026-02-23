import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Angle
from astropy.units import Quantity, Unit

import pygedm


def test_dm_to_dist():
    """Test that astropy units / angles work with dm_to_dist"""
    a = pygedm.dm_to_dist(204, -6.5, 200, method="ne2025")
    b = pygedm.dm_to_dist(
        Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200, method="ne2025"
    )
    c = pygedm.dm_to_dist(204, -6.5, 200 * Unit("pc cm^-3"), method="ne2025")
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_dm_to_dist_ne2001p():
    """Test that dm_to_dist works with ne2001p model"""
    a = pygedm.dm_to_dist(204, -6.5, 200, method="ne2001p")
    b = pygedm.dm_to_dist(
        Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200, method="ne2001p"
    )
    c = pygedm.dm_to_dist(204, -6.5, 200 * Unit("pc cm^-3"), method="ne2001p")
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_dist_to_dm():
    """Test that astropy units / angles work with dist_to_dm"""
    a = pygedm.dist_to_dm(204, -6.5, 200, method="ne2025")
    b = pygedm.dist_to_dm(
        Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200, method="ne2025"
    )
    c = pygedm.dist_to_dm(204, -6.5, 200 * Unit("pc"), method="ne2025")
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_dist_to_dm_ne2001p():
    """Test that dist_to_dm works with ne2001p model"""
    a = pygedm.dist_to_dm(204, -6.5, 200, method="ne2001p")
    b = pygedm.dist_to_dm(
        Angle(204, unit="degree"), Angle(-6.5, unit="degree"), 200, method="ne2001p"
    )
    c = pygedm.dist_to_dm(204, -6.5, 200 * Unit("pc"), method="ne2001p")
    assert a[0] == b[0] == c[0]
    assert a[1] == b[1] == c[1]


def test_calculate_electron_density_xyz():
    pc = Unit("pc")
    a = pygedm.calculate_electron_density_xyz(1, 2, 3, method="ne2025")
    b = pygedm.calculate_electron_density_xyz(1 * pc, 2, 3, method="ne2025")
    c = pygedm.calculate_electron_density_xyz(1, 2 * pc, 3, method="ne2025")
    d = pygedm.calculate_electron_density_xyz(1, 2, 3 * pc, method="ne2025")
    assert a == b == c == d


def test_basic():
    """Basic tests of NE2025 model """

    # FRB180301 value
    dm, tau = pygedm.dist_to_dm(204, -6.5, 25 * u.kpc, method="ne2025")
    assert np.isclose(dm.value, 154.7895926, atol=0.01)

    # Loop through distances and check round trip
    for dist in (50.0 * u.pc, 500.0 * u.pc, 5000.0 * u.pc):
        dm, tau = pygedm.dist_to_dm(0, 0, dist, method="ne2025")
        dist_out, tau = pygedm.dm_to_dist(0, 0, dm, method="ne2025")
        print(dist, dm, dist_out)
        assert np.isclose(dist_out.to("pc").value, dist.to("pc").value, rtol=0.1)


def test_igm():
    """Test that IGM mode FAILS as expected"""
    with pytest.raises(RuntimeError):
        pygedm.dm_to_dist(100, 10, 100, mode="igm", method="ne2025")


def test_dm_wrapper():
    """Run test against known values

    """
    test_data = {
        "l": [
            0,
            2,
            97.5,
        ],
        "b": [
            0,
            7.5,
            85.2,
        ],
        "dm": [10, 20, 11.1],
        "dist": [0.5637, 0.9096, 1.1436],
    }

    for ii in range(len(test_data["l"])):
        dist, tau_sc = pygedm.ne2025_wrapper.dm_to_dist(
            test_data["l"][ii], test_data["b"][ii], test_data["dm"][ii]
        )
        assert np.allclose(dist.to("kpc").value, test_data["dist"][ii], atol=2)

        dm, tau_sc = pygedm.ne2025_wrapper.dist_to_dm(
            test_data["l"][ii], test_data["b"][ii], test_data["dist"][ii]
        )
        assert np.allclose(dm.value, test_data["dm"][ii], atol=2)


def test_zero_dm():
    """Check that zero DM doesn't cause timeout bug

    Fortran code hangs if DM or D == 0
    """
    dist, tau_sc = pygedm.ne2025_wrapper.dm_to_dist(0, 0, 0)
    dm, tau_sc = pygedm.ne2025_wrapper.dist_to_dm(0, 0, 0)
    assert dist.value == dm.value == 0


def test_dm_wrapper_b0353():
    """Test against precomputed values for PSR B0353+52

    Values from NE2001 old online tool:
        l, b, dm = 149.0993, -0.5223, 102.5
        D = 2.746 kpc
        log_sm = -2.25 kpc m-20/3
        pulse broad = 6.57 us

    Expected values from NE2025:
        Note pulse broadening estimate much larger
        (<Quantity 3451.6450756 pc>, <Quantity 5.36335678e-05 s>)
    """
    l = 149.0993
    b = -0.5223
    dm = 102.50
    dist, tau_sc = pygedm.ne2025_wrapper.dm_to_dist(l, b, dm)

    assert np.isclose(dist.to("kpc").value, 3.452, atol=0.01, rtol=0)
    assert np.isclose(tau_sc.to("us").value, 53.63, atol=0.01, rtol=0)


def test_full_output():
    """Make sure full_output arg works"""
    a = pygedm.ne2025_wrapper.dist_to_dm(0, 0, 0.1, full_output=True)
    b = pygedm.ne2025_wrapper.dm_to_dist(0, 0, 10, full_output=True)
    assert isinstance(a, dict)
    assert isinstance(b, dict)


if __name__ == "__main__":
    test_basic()
    test_dm_to_dist()
    test_dist_to_dm()
    test_calculate_electron_density_xyz()
    test_igm()
    test_dm_wrapper()
    test_dm_wrapper_b0353()
    test_full_output()
