# MWPROP

Feb 2026 v2.0

MWPROP provides `NE2025p` and `NE2001p`, native Python implementations of the original Fortran code for NE2025/NE2001. The package also contains a required `scattering_functions` module. NE2025p/NE2001p are accessible from the command line, similar to the Fortran code, or within Python scripts (see below). The Fortran version of NE2025 is provided in the Github repository for MWPROP, along with distance and scattering data used to calibrate NE2025.


Please cite Ocker & Cordes (2026) for use of NE2025p.
The first description of the conversion between Fortran and Python is given in the NE2001p research note Ocker & Cordes (2024).
The NE2001 model is described in detail in Cordes & Lazio (2002, 2003).

-----

## Installation:

With pip:

`pip install mwprop`

On GitHub: https://github.com/stella-ocker/mwprop.

Executable scripts `NE2025p.py`, `NE2001p.py`, `los_diagnostics.py`, and `test_ne2025p.py` are automatically installed with pip and may be run from the command line in any directory.

**Dependencies**
- python >= 3.6 (might work with python >= 3.0)
- numpy
- matplotlib
- scipy
- astropy
- numba

-----

## Comparison to Fortran Version of NE2025/NE2001:

The input parameters and model components for the Milky Way are identical.

Output is nearly the same but not identical because integrations are done slightly differently to speed up the Python version.

Differences:

1.  Scattering computations are done as an option
    to give maximum speed for DM to D or D to DM computations
2.  Output can be printed in console or returned in Python dictionaries, Dn, Dv, Du, Dd:
        Dn = names of variables
        Dv = values
        Du = units
        Dd = descriptions
3.  Numerical integrations use numpy.trapz instead of step-by-step summing.
4.  Clumps are prefiltered for inclusion along a specified line of sight (LoS).
5.  Large-scale components (thin and thick disks, spiral arms) are sampled coarsely
    (0.1 kpc default sample interval) and cubic splines are used to resample onto
    a fine grid.
6.  Components with small scale structure (local ISM, Galactic center, clumps, voids)
    are sampled directly on a fine grid.
7.  Different routines are used for execution to find distance from DM and for DM from distance.
    Cases where the DM exceeds the maximum Galactic value will return a warning; distances predicted for these cases are unreliable (see Ocker & Cordes 2024).
8.  Config file:  code imports,  reading parameter files, and setting key values
    (like the Sun-GC distance = 8.5 kpc) are done by importing mwprop.nemod.config_nemod
    Note the 8.5 kpc distance is needed because model parameters were determined using
    this value for the original NE2001 code. It cannot be changed by itself.

The Python implementation is about 45 times slower than in Fortran. For computations requiring speed, we recommend the Fortran release.

-----

## Usage:

Command Line Usage: `NE2025p.py ldeg bdeg dmd ndir [-options]`

Required arguments:  ldeg bdeg dmd ndir     (Same as Fortran version)

    ldeg, bdeg = Galactic longitude and latitude in deg
    dmd = DM (pc/cc) or d (kpc)
    ndir >= 0   DM -> d
          < 0    d -> DM

Optional requirements and defaults:

    -a  --analysis   Calculates contributions from different model components and writes to file in user's working directory.
                     [False]
    -p  --plotting   Plots DM, n_e, scattering vs distance along LoS; saves pdf files in user's working directory.
                     [False]
    -s  --scattering Calculates scattering and scintillation parameters, etc.
                     [False]
    -m  --modern     Modern output (suppresses classic Fortran output)
                     [False]
    -v  --verbose    Writes results to std out in the 'classic' form of the Fortran version.
                     [False]
    -b  --debug      Prints various diagnostics for clumps and integration. Only works for DM->D.
                     [False]
    -e  --explain    Show long-form explanation of code.
                     [False]

Script/iPython Usage:

The ne2025()/ne2001() functions evaluate NE2025p/NE2001p and can be imported from the mwprop.nemod.NE2025 (or mwprop.nemod.NE2001) module. The function output consists of four dictionaries:

    Dk    => Dictionary keys for output values
    Dv    => Numerical output values
    Du    => Units of output values
    Dd    => Extended output description

Additional options for ne2025()/ne2001():

    classic = True      => print traditional Fortran output (default = True)
    dmd_only = True     => compute only DM or distance, no scattering quantities (default = False)
    do_analysis = True  => calculate and output diagnostic quantities for the line of sight (default = False)
    plotting = True     => plot DM vs distance along LoS, SM also if dmd_only = False (default = False)
    verbose = True      => print results to std out (default = False)

iPython Example:

    >>> from mwprop.nemod import *
    >>> from mwprop.nemod.NE2025 import ne2025
    >>> Dk,Dv,Du,Dd = ne2025(ldeg=45,bdeg=5,dmd=50,ndir=1,classic=False)
    >>> Dk['DIST']
    'DIST'
    >>> Dv['DIST']
    3.7938
    >>> Du['DIST']
    'kpc'
    >>> Dd['DIST']
    'ModelDistance'

In analysis and plotting modes, output files are saved to a folder called 'output_ne2025p' or 'output_ne2001p' (depending on the function called) in the user's working directory.

v2.0 introduces warnings when a clump or void is intersected in the model. Warnings can be turned off using warnings.filterwarnings(). Users can see a detailed breakdown of clump and void contributions by running the `los_diagnostics.py` script described below. While not generally recommended, users can turn clumps or voids off completely by setting the weight parameters wgcN and wgvN to 0 in the params/gal25.inp file (which requires recompiling the Python package).

-----

## Diagnostic code los_diagnostics.py


Plots electron density, DM, and C_n^2 along the line of sight designated by l, b, DM or d, and ndir = 1 or -1 (as with NE2001).
Also shows the line of sight projected onto the Galactic plane along with spiral arms used in NE2001/NE2025.

This code can be run from any directory if `mwprop` is fully installed. Outputs are saved to a folder created in the user's working directory.

    Usage:
    los_diagnostics.py  l  b  dmd  ndir  v
    with l, b in deg, dmd = DM (pc/cc) or distance (kpc) for ndir > or < 0, v = 2025 or 2001 for NE2025 or NE2001

-----

## test_ne2025p.py

Users wishing to test if the Python installation behaves as expected may also run this executable script from any command line. A specific sightline is evaluated and compared to expected values, and the percent error between the expected and calculated parameters are printed.

-----

## Known Issues

mwprop v1.0: Extragalactic scattering times and scintillation bandwidths, TAU_X / SBW_X, output by mwprop.ne2001p v1.0 (pre-NE2025) are too small / large (respectively) by a factor of 2. This error is corrected in v2.0.

mwprop v2.0: Two of the smallest clumps in the model (1745-2900 and OH40.6-0.2) will have large differences between the Fortran and Python outputs, due to small differences in the sampling of the numerical integration that are minor for most clumps but exaggerated when the clump size is close to the integration grid sampling. For these sightlines, use of the Fortran code is recommended.
