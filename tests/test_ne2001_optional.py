import pytest

import pygedm
from pygedm import pygedm as _pygedm


def test_ne2001_missing_extension_raises(monkeypatch):
    """When ne21c isn't built, method='ne2001' should raise a clear RuntimeError."""
    monkeypatch.setattr(_pygedm, "HAS_NE2001", False)

    with pytest.raises(RuntimeError, match="ne2001p"):
        pygedm.dm_to_dist(204, -6.5, 200, method="ne2001")
    with pytest.raises(RuntimeError, match="ne2001p"):
        pygedm.dist_to_dm(204, -6.5, 200, method="ne2001")
    with pytest.raises(RuntimeError, match="ne2001p"):
        pygedm.calculate_electron_density_xyz(1, 2, 3, method="ne2001")
    with pytest.raises(RuntimeError, match="ne2001p"):
        pygedm.calculate_electron_density_lbr(204, -6.5, 3000, method="ne2001")
