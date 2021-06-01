#include "ne21.h"

namespace py = pybind11;

// Dispersion measure to distance
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
std::map<std::string, float> dist_to_dm(float gl_rad, float gb_rad, float dm) {
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



PYBIND11_MODULE(ne21c, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("dome", &dome, "A function which does soemthing");
    m.def("dm_to_dist", &dm_to_dist,  "Convert DM to distance");
    m.def("dist_to_dm", &dist_to_dm, "Convert distance to DM");
}
