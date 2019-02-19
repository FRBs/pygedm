import pyymw16
import numpy as np

def test_basic():
    """ Basic tests of YMW16 model

    Tested against online YMW16 interface
    http://www.atnf.csiro.au/research/pulsar/ymw16/index.php
    """

    a = pyymw16.calculate_electron_density_xyz(1,2,3)
    assert np.isclose(a, 5.220655, atol=0.0001)

    a = pyymw16.calculate_electron_density_lbr(0,0,4000)
    assert np.isclose(a,  0.388407, atol=0.0001)

    # FRB180301 value
    dm, tau = pyymw16.dist_to_dm(204, -6.5, 25000)
    assert np.isclose(dm.value, 252.0501, atol=0.01)

    # Loop through distances and check round trip
    for dist in (10., 100., 1000.):
        dm, tau = pyymw16.dist_to_dm(0, 0, dist)
        dist_out, tau = pyymw16.dm_to_dist(0, 0, dm.value)
        assert np.isclose(dist_out.value, dist, rtol=0.1)
