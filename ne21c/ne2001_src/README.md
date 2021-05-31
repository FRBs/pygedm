## NE2001

This is a version of [NE2001](https://www.nrl.navy.mil/rsd/RORF/ne2001/)
code, that is modified for Python bindings using `f2py`.

To compile the code, use the `compile_ne2001.py` script. This will generate
Python-compatible shared objects that interface with the Fortran
code, namely the `dmdsm` and `density_2001` functions. 

The `ne2001.py` file wraps the bindings to provide 
the following three functions:

* `dm_to_dist` - Convert DM to distance and compute scattering timescale
* `dist_to_dm` - Convert distance to DM and compute scattering timescale
* `density_xyz` - Compute electron density at Galactocentric XYZ coordinates
