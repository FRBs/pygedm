import pygedm
import numpy as np
import pytest

def test_healpix_map():
    hpmap = pygedm.generate_healpix_dm_map(dist=30000, nside=32)

    ## Test against a few precomputed values
    assert len(hpmap) == 12288
    assert np.isclose(hpmap[0], 19.607964)
    assert np.isclose(hpmap[1729], 29.772997)

    hpmap = pygedm.generate_healpix_dm_map(nside=2, method='yt2020')
    assert np.isclose(hpmap[0], 35.190342)

    hpmap = pygedm.generate_healpix_dm_map(nside=2, method='yt2020_analytic')
    assert np.isclose(hpmap[0], 35.01475)

    # Check that error is raised if healpy isn't installed
    with pytest.raises(RuntimeError):
        pygedm.HAS_HEALPIX = False
        pygedm.generate_healpix_dm_map(dist=30000, nside=32)



if __name__ == "__main__":
    test_healpix_map()