"""Tests for numpy array input to dist_to_dm / dm_to_dist (issue #21)."""

import numpy as np
import pytest

import pygedm

requires_ne2001 = pytest.mark.skipif(
    not pygedm.pygedm.HAS_NE2001, reason="ne21c not built (f2c unavailable)"
)


def _loop_dm_to_dist(gl, gb, dm, method):
    dists, taus = [], []
    for gl_i, gb_i, dm_i in zip(gl, gb, dm):
        d, t = pygedm.dm_to_dist(gl_i, gb_i, dm_i, method=method)
        dists.append(d.value)
        taus.append(t.value)
    return np.array(dists), np.array(taus)


def _loop_dist_to_dm(gl, gb, dist, method):
    dms, taus = [], []
    for gl_i, gb_i, dist_i in zip(gl, gb, dist):
        dm, t = pygedm.dist_to_dm(gl_i, gb_i, dist_i, method=method)
        dms.append(dm.value)
        taus.append(t.value)
    return np.array(dms), np.array(taus)


def test_ymw16_array_dm_to_dist():
    gl = np.array([204.0, 0.0, 280.46])
    gb = np.array([-6.5, 0.0, -32.88])
    dm = np.array([200.0, 50.0, 100.0])

    dist, tau = pygedm.dm_to_dist(gl, gb, dm, method="ymw16")
    dist_loop, tau_loop = _loop_dm_to_dist(gl, gb, dm, method="ymw16")

    assert np.allclose(dist.value, dist_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)


def test_ymw16_array_dist_to_dm():
    gl = np.array([204.0, 0.0, 280.46])
    gb = np.array([-6.5, 0.0, -32.88])
    dist = np.array([25000.0, 1000.0, 50000.0])

    dm, tau = pygedm.dist_to_dm(gl, gb, dist, method="ymw16")
    dm_loop, tau_loop = _loop_dist_to_dm(gl, gb, dist, method="ymw16")

    assert np.allclose(dm.value, dm_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)


@requires_ne2001
def test_ne2001_array_dm_to_dist():
    gl = np.array([20.0, 30.0, 0.0])
    gb = np.array([-10.0, 30.0, 0.0])
    dm = np.array([10.0, 50.0, 0.0])  # zero DM element exercises the hang guard

    dist, tau = pygedm.dm_to_dist(gl, gb, dm, method="ne2001")
    dist_loop, tau_loop = _loop_dm_to_dist(gl, gb, dm, method="ne2001")

    assert np.allclose(dist.value, dist_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)
    assert dist.value[2] == 0.0


@requires_ne2001
def test_ne2001_array_dist_to_dm():
    gl = np.array([20.0, 30.0, 0.0])
    gb = np.array([-10.0, 30.0, 0.0])
    dist = np.array([1000.0, 2000.0, 0.0])  # zero distance element exercises the hang guard

    dm, tau = pygedm.dist_to_dm(gl, gb, dist, method="ne2001")
    dm_loop, tau_loop = _loop_dist_to_dm(gl, gb, dist, method="ne2001")

    assert np.allclose(dm.value, dm_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)
    assert dm.value[2] == 0.0


@requires_ne2001
def test_ne2001_array_full_output():
    """full_output=True on the array path returns the raw per-field arrays."""
    from pygedm import ne2001_wrapper

    l = np.array([20.0, 30.0])
    b = np.array([-10.0, 30.0])
    dm = np.array([10.0, 50.0])
    dist = np.array([1.0, 2.0])

    d = ne2001_wrapper.dm_to_dist(l, b, dm, full_output=True)
    assert set(d.keys()) >= {"dist", "sm", "smtau", "smtheta", "smiso", "tau_sc"}
    assert len(d["dist"]) == 2

    d = ne2001_wrapper.dist_to_dm(l, b, dist, full_output=True)
    assert set(d.keys()) >= {"dm", "sm", "smtau", "smtheta", "smiso", "tau_sc"}
    assert len(d["dm"]) == 2


@pytest.mark.parametrize("method", ["ne2025", "ne2001p"])
def test_ne2025_array_dm_to_dist(method):
    gl = np.array([20.0, 30.0, 0.0])
    gb = np.array([-10.0, 30.0, 0.0])
    dm = np.array([10.0, 50.0, 0.0])  # zero DM element exercises the "spline failure" guard

    dist, tau = pygedm.dm_to_dist(gl, gb, dm, method=method)
    dist_loop, tau_loop = _loop_dm_to_dist(gl, gb, dm, method=method)

    assert np.allclose(dist.value, dist_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)
    assert dist.value[2] == 0.0


@pytest.mark.parametrize("method", ["ne2025", "ne2001p"])
def test_ne2025_array_dist_to_dm(method):
    gl = np.array([20.0, 30.0, 0.0])
    gb = np.array([-10.0, 30.0, 0.0])
    dist = np.array([100.0, 200.0, 0.0])  # zero distance element exercises the "spline failure" guard

    dm, tau = pygedm.dist_to_dm(gl, gb, dist, method=method)
    dm_loop, tau_loop = _loop_dist_to_dm(gl, gb, dist, method=method)

    assert np.allclose(dm.value, dm_loop, rtol=1e-5)
    assert np.allclose(tau.value, tau_loop, rtol=1e-5)
    assert dm.value[2] == 0.0


def test_ne2025_array_full_output():
    """full_output=True on the array path returns a list of per-element dicts."""
    from pygedm import ne2025_wrapper

    gl = np.array([20.0, 30.0])
    gb = np.array([-10.0, 30.0])
    dm = np.array([10.0, 50.0])
    dist = np.array([100.0, 200.0])

    results = ne2025_wrapper.dm_to_dist(gl, gb, dm, full_output=True)
    assert isinstance(results, list) and len(results) == 2
    assert all("tau_sc" in r for r in results)

    results = ne2025_wrapper.dist_to_dm(gl, gb, dist, full_output=True)
    assert isinstance(results, list) and len(results) == 2
    assert all("tau_sc" in r for r in results)


if __name__ == "__main__":
    test_ymw16_array_dm_to_dist()
    test_ymw16_array_dist_to_dm()
    test_ne2001_array_dm_to_dist()
    test_ne2001_array_dist_to_dm()
    test_ne2001_array_full_output()
    for _method in ("ne2025", "ne2001p"):
        test_ne2025_array_dm_to_dist(_method)
        test_ne2025_array_dist_to_dm(_method)
    test_ne2025_array_full_output()
