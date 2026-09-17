/*Copyright (C) 2016, 2017  J. M. Yao, R. N. Manchester, N. Wang.

This file is part of the YMW16 program. YMW16 is a model for the
distribution of free electrons in the Galaxy, the Magellanic Clouds
and the inter-galactic medium that can be used to estimate distances
for real or simulated pulsars and fast radio bursts (FRBs) based on
their position and dispersion measure.

YMW16 is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

YMW16 is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License,
available at http://www.gnu.org/licenses/, for more details.

Please report any issues or bugs at
https://bitbucket.org/psrsoft/ymw16/issues/new/ or directly to the
authors. Please provide an example illustrating the problem.

Jumei Yao (yaojumei@xao.ac.cn), Richard N Manchester
(dick.manchester@csiro.au), Na Wang (na.wang@xao.ac.cn).
*/

#include "cn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

int m_3, ww,m_5, m_6, m_7; // These are shared among other *.cpp files

extern double ne_crd(double *x, double *y, double *z, double *gl,
                     double *gb, double *dd, int ncrd, int vbs,
                     char *dirname, char *text);

extern std::map<std::string, float> dmdtau2(double gl, double gb, double dordm,
       double DM_Host, int ndir, int np, int vbs, char *dirname, char *text);

extern std::map<std::string, float> frb_d(double DDM, double DM_Gal, double DM_MC,
       double DM_Host, int uu, int vbs, char* text);

// Array-input version of dmdtau2: loops over gl/gb/dordm in C++ (avoiding a
// Python-level loop) and returns a dict of numpy arrays, one per output key.
// The scalar dmdtau2() can return slightly different key sets across calls
// depending on the input values (not just ndir/np) -- e.g. DM_Gal/DM_MC are
// only inserted along certain branches -- so we take the union of keys seen
// across all elements and backfill missing entries with NaN.
py::dict dmdtau_arr(py::array_t<double> gl, py::array_t<double> gb, py::array_t<double> dordm,
       double DM_Host, int ndir, int np, int vbs, char *dirname, char *text) {
    auto gl_buf = gl.unchecked<1>();
    auto gb_buf = gb.unchecked<1>();
    auto dordm_buf = dordm.unchecked<1>();
    ssize_t n = gl_buf.shape(0);
    if (gb_buf.shape(0) != n || dordm_buf.shape(0) != n) {
        throw std::invalid_argument("gl, gb, and dordm must have the same shape");
    }

    std::vector<std::map<std::string, float>> results;
    results.reserve(n);
    for (ssize_t i = 0; i < n; i++) {
        results.push_back(dmdtau2(gl_buf(i), gb_buf(i), dordm_buf(i), DM_Host,
                                   ndir, np, vbs, dirname, text));
    }

    std::vector<std::string> keys;
    for (const auto &r : results) {
        for (const auto &kv : r) {
            if (std::find(keys.begin(), keys.end(), kv.first) == keys.end()) {
                keys.push_back(kv.first);
            }
        }
    }

    py::dict out;
    double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto &key : keys) {
        py::array_t<double> arr(n);
        auto arr_buf = arr.mutable_unchecked<1>();
        for (ssize_t i = 0; i < n; i++) {
            auto it = results[i].find(key);
            arr_buf(i) = (it != results[i].end()) ? it->second : nan;
        }
        out[py::str(key)] = arr;
    }
    return out;
}


PYBIND11_MODULE(ymw16, m) {
m.doc() = R"pbdoc(pyYMW16 -- python binding for YMW16.
The program YMW16 computes distances for Galactic pulsars, Magellanic Cloud pulsars,
and FRBs from their Galactic coordinates and DMs using the YMW16 model parameters.
It also does the reverse calculation, computing DMs that correspond to given
Galactic coordinates and distances. An estimate of the scattering timescale tsc
is output for Galactic and Magellanic Cloud pulsars and FRBs.

Ref: J. M. Yao, R. N. Manchester, and N. Wang (2017), doi:10.3847/1538-4357/835/1/29

)pbdoc"; // optional module docstring

m.def("frb_d", &frb_d, R"pbdoc(
    Args:
      DDM: distance or DM (cm−3 pc)
      DM_Gal: DM of galaxy along line of sight (cm−3 pc)
      DM_MC: DM of Magellanic cloud along line of sight (cm−3 pc)
      DM_Host: Dispersion measure of the FRB host galaxy in
               the observer frame (default 100 cm−3 pc). (Note: if
               present, DM_Host is ignored for Gal and MC modes.)
      uu: Direction, 1 for dist to DM, 0 for DM to dist
      vbs: Verbostiy level, 0, 1, or 2
    Returns:
      Python dictionary with computed values.
      tsc has units of seconds.
    )pbdoc",
py::arg("DDM"),
py::arg("DM_Gal"),
py::arg("DM_MC"),
py::arg("DM_Host"),
py::arg("uu"),
py::arg("vbs"),
py::arg("text")
);


m.def("dmdtau", &dmdtau2, R"pbdoc(
    Args:
      gl: Galactic longitude (deg.)
      gb: Galactic latitude (deg.)
      dordm: One of DM (cm−3 pc) or distance, depending
             on ndir. Distance has units of pc for modes Gal and MC
              and Mpc for mode IGM
      DM_Host: Dispersion measure of the FRB host galaxy in
               the observer frame (default 100 cm−3 pc). (Note: if
               present, DM_Host is ignored for Gal and MC modes.)
      ndir: ndir=1 converts from DM to distance and
            ndir=2 converts from distance to DM.
      np: -1 for IGM, 0 for Mag clouds, 1 for galaxy
      vbs: Verbostiy level, 0, 1, or 2
      dirname: directory where data files are stored
      text: Text to prepend in print statement.
    Returns:
      Python dictionary with computed values.
      tsc has units of seconds.
    )pbdoc",
py::arg("gl"),
py::arg("gb"),
py::arg("dordm"),
py::arg("DM_Host"),
py::arg("ndir"),
py::arg("np"),
py::arg("vbs"),
py::arg("dirname"),
py::arg("text")
);

m.def("dmdtau_arr", &dmdtau_arr, R"pbdoc(
    Array version of dmdtau. gl, gb and dordm are numpy arrays of the same
    shape; DM_Host, ndir, np, vbs, dirname and text are shared across all
    elements. Returns a dict of numpy arrays, one per output key (any key
    absent for a given element is filled with NaN).
    )pbdoc",
py::arg("gl"),
py::arg("gb"),
py::arg("dordm"),
py::arg("DM_Host"),
py::arg("ndir"),
py::arg("np"),
py::arg("vbs"),
py::arg("dirname"),
py::arg("text")
);

m.def("ne_crd", &ne_crd, R"pbdoc(
    Calculate electron density at a given point with galactocentric coordinates
    (x, y, z) OR with (gl, gb, dist).

    Args:
      (x, y, z): input Galactocentric x, y and z in pc
      (gl, gb, dist): input gl, gb in deg, Dist in pc
      ncrd: if ncrd==1, use xyz coords. If ncrd==2 use gl gb dist coords.
      vbs: Verbostiy level, 0, 1, or 2
      dirname: directory where data files are stored
      text: Text to prepend in print statement.
    )pbdoc",
py::arg("x"),
py::arg("y"),
py::arg("z"),
py::arg("gl"),
py::arg("gb"),
py::arg("dd"),
py::arg("ncrd"),
py::arg("vbs"),
py::arg("dirname"),
py::arg("text")
);

}
