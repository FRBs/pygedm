#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <stdexcept>

#ifdef __cplusplus
 extern "C" {
#endif
#include "ne21.c"
// WAR: Quell undefined symbol: MAIN_
int  MAIN__( ) {  return 0; }
#ifdef __cplusplus
}
#endif /* extern "C" */

namespace py = pybind11;

// Dispersion measure to distance
// DM <--> Distance (dmdsm)
integer ndir;
real r1, r2, dm, dist, sm, smtau, smtheta, smiso;
static char limit[1];
std::map<std::string, float> dm_to_dist(float gl_rad, float gb_rad, float dm) {
    std::map<std::string, float> result;
    dist = 0;
    sm = 0;
    smtau = 0;
    smtheta = 0;
    smiso = 0;
    ndir = 1;
    // call dmdsm_ from shared library
    dmdsm_(&gl_rad, &gb_rad, &ndir, &dm, &dist, limit, &sm, &smtau, &smtheta,
        &smiso, (ftnlen)1);
    result.insert(std::make_pair("dist", dist));
    result.insert(std::make_pair("sm", sm));
    result.insert(std::make_pair("smtau", smtau));
    result.insert(std::make_pair("smtheta", smtheta));
    result.insert(std::make_pair("smiso", smiso));
    return result;
}

// Distance to dispersion measure
std::map<std::string, float> dist_to_dm(float gl_rad, float gb_rad, float dist) {
    std::map<std::string, float> result;
    dm = 0;
    sm = 0;
    smtau = 0;
    smtheta = 0;
    smiso = 0;
    ndir = -1;
    // call dmdsm_ from shared library
    dmdsm_(&gl_rad, &gb_rad, &ndir, &dm, &dist, limit, &sm, &smtau, &smtheta,
        &smiso, (ftnlen)1);
    result.insert(std::make_pair("dm", dm));
    result.insert(std::make_pair("sm", sm));
    result.insert(std::make_pair("smtau", smtau));
    result.insert(std::make_pair("smtheta", smtheta));
    result.insert(std::make_pair("smiso", smiso));
    return result;
}

// Density at point x,y,z (density_2001)
real x, y, z__, ne1, ne2, nea, negc, nelism, necn, nevn;
real f1, f2, fa, fgc, flism, fcn, fvn;
integer whicharm, wlism, wldr, wlhb, wlsb, wloopi, hitclump, hitvoid, wvoid;

// Density at x, y, z
std::map<std::string, float> density_xyz(float x, float y, float z) {
    std::map<std::string, float> result;
    ne1 = 0;
    ne2 = 0;
    nea = 0;
    negc = 0;
    nelism = 0;
    necn = 0;
    nevn = 0;
    f1 = 0;
    f2 = 0;
    fa = 0;
    fgc = 0;
    flism = 0;
    fcn = 0;
    whicharm = 0;
    wlism = 0;
    wldr = 0;
    wlhb = 0;
    wlsb = 0;
    wloopi = 0;
    hitclump = 0;
    hitvoid = 0;
    wvoid = 0;    
    density_2001__(&x, &y, &z, &ne1, &ne2, &nea, &negc, &nelism, &necn, &nevn,
    &f1, &f2, &fa, &fgc, &flism, &fcn, &fvn, &whicharm, &wlism, &wldr, &wlhb, &wlsb, 
    &wloopi, &hitclump, &hitvoid, &wvoid);
    float netot = ne1 + ne2 + nea + negc + nelism + necn + nevn;
    result.insert(std::make_pair("ne", netot));
    return result;
}

// Array versions of dm_to_dist/dist_to_dm: loop over the inputs in C++
// (avoiding a Python-level loop) and return a dict of numpy arrays. Elements
// with dm/dist close to zero are skipped rather than passed into the
// Fortran-derived dmdsm_ call, which is known to hang at dm==0/dist==0 (see
// the WAR comment in the scalar Python wrapper) -- those outputs are set to
// 0 directly instead.
py::dict dm_to_dist_arr(py::array_t<float> gl_rad, py::array_t<float> gb_rad, py::array_t<float> dm_arr) {
    auto gl_buf = gl_rad.unchecked<1>();
    auto gb_buf = gb_rad.unchecked<1>();
    auto dm_buf = dm_arr.unchecked<1>();
    ssize_t n = gl_buf.shape(0);
    if (gb_buf.shape(0) != n || dm_buf.shape(0) != n) {
        throw std::invalid_argument("gl_rad, gb_rad, and dm must have the same shape");
    }

    py::array_t<float> out_dist(n), out_sm(n), out_smtau(n), out_smtheta(n), out_smiso(n);
    auto dist_buf = out_dist.mutable_unchecked<1>();
    auto sm_buf = out_sm.mutable_unchecked<1>();
    auto smtau_buf = out_smtau.mutable_unchecked<1>();
    auto smtheta_buf = out_smtheta.mutable_unchecked<1>();
    auto smiso_buf = out_smiso.mutable_unchecked<1>();

    for (ssize_t i = 0; i < n; i++) {
        if (std::fabs(dm_buf(i)) <= 1e-8f) {
            dist_buf(i) = 0.0f;
            sm_buf(i) = 0.0f;
            smtau_buf(i) = 0.0f;
            smtheta_buf(i) = 0.0f;
            smiso_buf(i) = 0.0f;
            continue;
        }
        std::map<std::string, float> r = dm_to_dist(gl_buf(i), gb_buf(i), dm_buf(i));
        dist_buf(i) = r["dist"];
        sm_buf(i) = r["sm"];
        smtau_buf(i) = r["smtau"];
        smtheta_buf(i) = r["smtheta"];
        smiso_buf(i) = r["smiso"];
    }

    py::dict out;
    out["dist"] = out_dist;
    out["sm"] = out_sm;
    out["smtau"] = out_smtau;
    out["smtheta"] = out_smtheta;
    out["smiso"] = out_smiso;
    return out;
}

py::dict dist_to_dm_arr(py::array_t<float> gl_rad, py::array_t<float> gb_rad, py::array_t<float> dist_arr) {
    auto gl_buf = gl_rad.unchecked<1>();
    auto gb_buf = gb_rad.unchecked<1>();
    auto dist_buf_in = dist_arr.unchecked<1>();
    ssize_t n = gl_buf.shape(0);
    if (gb_buf.shape(0) != n || dist_buf_in.shape(0) != n) {
        throw std::invalid_argument("gl_rad, gb_rad, and dist must have the same shape");
    }

    py::array_t<float> out_dm(n), out_sm(n), out_smtau(n), out_smtheta(n), out_smiso(n);
    auto dm_buf = out_dm.mutable_unchecked<1>();
    auto sm_buf = out_sm.mutable_unchecked<1>();
    auto smtau_buf = out_smtau.mutable_unchecked<1>();
    auto smtheta_buf = out_smtheta.mutable_unchecked<1>();
    auto smiso_buf = out_smiso.mutable_unchecked<1>();

    for (ssize_t i = 0; i < n; i++) {
        if (std::fabs(dist_buf_in(i)) <= 1e-8f) {
            dm_buf(i) = 0.0f;
            sm_buf(i) = 0.0f;
            smtau_buf(i) = 0.0f;
            smtheta_buf(i) = 0.0f;
            smiso_buf(i) = 0.0f;
            continue;
        }
        std::map<std::string, float> r = dist_to_dm(gl_buf(i), gb_buf(i), dist_buf_in(i));
        dm_buf(i) = r["dm"];
        sm_buf(i) = r["sm"];
        smtau_buf(i) = r["smtau"];
        smtheta_buf(i) = r["smtheta"];
        smiso_buf(i) = r["smiso"];
    }

    py::dict out;
    out["dm"] = out_dm;
    out["sm"] = out_sm;
    out["smtau"] = out_smtau;
    out["smtheta"] = out_smtheta;
    out["smiso"] = out_smiso;
    return out;
}

PYBIND11_MODULE(ne21c, m) {
    m.doc() = R"pbdoc(ne21c -- python bindings to C port of NE2001.
    The program NE2001 computes distances for Galactic pulsars, Magellanic Cloud pulsars,
    and FRBs from their Galactic coordinates and DMs using the NE2001 model parameters.
    It also does the reverse calculation, computing DMs that correspond to given
    Galactic coordinates and distances.

    Ref: Cordes, J. M., & Lazio, T. J. W. (2002)
    https://ui.adsabs.harvard.edu/abs/2002astro.ph..7156C/abstract
    )pbdoc"; // optional module docstring

    m.def("dm_to_dist", &dm_to_dist, R"pbdoc(
    Convert DM to a distance 

    Args:
        gl_rad (float): Galactic longitude in radians
        gb_rad (float): Galactic latitude in radians
        dm (float): Dispersion measure in pc/cm3
    
    Returns:
      Python dictionary with computed values.
    )pbdoc", 
    py::arg("gl_rad"),
    py::arg("gb_rad"),
    py::arg("dm")
    );

    m.def("dm_to_dist_arr", &dm_to_dist_arr, R"pbdoc(
    Array version of dm_to_dist. gl_rad, gb_rad and dm are numpy arrays of
    the same shape. Elements with dm close to zero are returned as all-zero
    rows rather than evaluated (avoids a known hang in the underlying model
    at dm==0). Returns a dict of numpy arrays.
    )pbdoc",
    py::arg("gl_rad"),
    py::arg("gb_rad"),
    py::arg("dm")
    );

    m.def("dist_to_dm", &dist_to_dm, R"pbdoc(
    Convert a distance in kpc to a DM estimate 

    Args:
        gl_rad (float): Galactic longitude in radians
        gb_rad (float): Galactic latitude in radians
        dist (float): Distance in kpc
    
    Returns:
      Python dictionary with computed values.
    )pbdoc", 
    py::arg("gl_rad"),
    py::arg("gb_rad"),
    py::arg("dist")
    );

    m.def("dist_to_dm_arr", &dist_to_dm_arr, R"pbdoc(
    Array version of dist_to_dm. gl_rad, gb_rad and dist are numpy arrays of
    the same shape. Elements with dist close to zero are returned as all-zero
    rows rather than evaluated (avoids a known hang in the underlying model
    at dist==0). Returns a dict of numpy arrays.
    )pbdoc",
    py::arg("gl_rad"),
    py::arg("gb_rad"),
    py::arg("dist")
    );

    m.def("density_xyz", &density_xyz, R"pbdoc(
    Compute electron density at galactocentric coordinates (X, Y, Z)
    
    x,y,z are Galactocentric Cartesian coordinates, measured in kpc (NOT pc!)
    with the axes parallel to (l, b) = (90, 0), (180, 0), and (0, 90) degrees
    
    Args:
        x (float): Galactocentric coordinates in kpc
        y (float): Galactocentric coordinates in kpc
        z (float): Galactocentric coordinates in kpc
        
    Returns:
        Python dictionary with computed values.
    )pbdoc",
    py::arg("x"),
    py::arg("y"),
    py::arg("z")
   );

    m.def("_main", &MAIN__, "Test to include MAIN__ symbol");
}

