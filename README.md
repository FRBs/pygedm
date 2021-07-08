[![Python application](https://github.com/FRBs/pygedm/actions/workflows/python-app.yml/badge.svg)](https://github.com/FRBs/pygedm/actions/workflows/python-app.yml)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Coverage](https://codecov.io/gh/FRBs/pygedm/branch/master/graph/badge.svg?token=TlBiPzD7DP)](https://codecov.io/gh/FRBs/pygedm)[![Documentation Status](https://readthedocs.org/projects/pygedm/badge/?version=latest)](https://pygedm.readthedocs.io/en/latest/?badge=latest)

# PyGEDM
_Python bindings for the YMW16, NE2001 and YT2020 electron density models_

This package is a Python interface to the YMW16 and NE2001 electron density models, and YT2020 halo model.
The Yao, Manchester and Wang (2017, [Astrophys. J., 835, 29](https://iopscience.iop.org/article/10.3847/1538-4357/835/1/29/meta);
[arXiv:1610.09448](https://arxiv.org/abs/1610.09448)) YMW16 electron density model, is written in C++, and the Cordes and Lazio 
(2001, [arXiv:0207156](https://arxiv.org/abs/astro-ph/)) NE2001 model is written in FORTRAN. This package, PyGEDM, wraps these
two codes using [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html) to make them usable from Python. Here, we have converted NE2001 to C++ using `f2c`.

### Usage

Some usage examples can be found in the [examples directory](https://github.com/telegraphic/pygedm/tree/master/examples). 

```python
import pygedm

# calculate DM at a given distance
DM, tau_sc = pygedm.dist_to_dm(204.0, -6.5, 200, method='ne2001')
DM, tau_sc = pygedm.dist_to_dm(204.0, -6.5, 200, method='ymw16')

# calculate distance for a given sky position and DM
dist, tau_sc = pygedm.dm_to_dist(123.4, 4.0, 200)

# calculate N_e density at xyz galactocentric coordinates
ne = pygedm.calculate_electron_density_xyz(1.0, 2.0, 3.0)

# calculate N_e density at Galactic lat/long/distance coords
ne = pygedm.calculate_electron_density_lbr(204.0, -6.5, 3000.0)

# Calculate halo DM contribution
dm_halo = pygedm.calculate_halo_dm(gl=0, gb=30)
```

The methods return astropy [Quantities](http://docs.astropy.org/en/stable/units/quantity.html#quantity), which have units attached, and can accept astropy [Angles](http://docs.astropy.org/en/stable/coordinates/angles.html#working-with-angles) and Quantities as arguments:

```python
import pygedm
import astropy.units as u
import astropy.coordinates as c
DM = u.Quantity(10.0, unit='pc cm^-3')
ra, dec = c.Angle(23.0, unit='hourangle'), c.Angle('-43:00:02', unit='degree')
sky_coords = c.SkyCoord(ra, dec, frame='icrs')
dist, tau_sc = pygedm.dm_to_dist(sky_coords.galactic.l, sky_coords.galactic.b, DM)

print(dist.to('lyr'))
>> 3362.16343117 lyr
print(tau_sc.to('ns'))
>> 7.758686138 ns
```


### Installation

Requires `pybind11`, `astropy`, `numpy`, `scipy`, a newish C compiler with C++11 support (Ubuntu 16.04+ default gcc will work), plus `f2c`. 

Pre-compiled wheels are generated for linux, meaning you can skip the `f2c` install step. Alternatively you can install via `apt-get f2c` in Ubuntu, or  via `conda install -c conda-forge f2c` if you use conda.

For MacOS, you are best off using `conda` and getting `f2c` via `conda install -c conda-forge f2c`. Windows is not currently supported.

You should be able to install with:

```
pip install pygedm
```

or 

```
pip install git+https://github.com/telegraphic/pygedm
```

to install the latest version from github. Alternatively, download this repository and install via

```
pip install .
```

To run unit tests, run `python setup.py test`. Note that these tests only check the Python bindings, 
not the underlying C/FORTRAN source code (which is not supplied with unit tests).

### References

If using PyGEDM in a journal article, please remember to cite the underlying electron density models:

[Cordes, J. M., & Lazio, T. J. W. (2002)](https://ui.adsabs.harvard.edu/abs/2002astro.ph..7156C/abstract),  
_NE2001.I. A New Model for the Galactic Distribution of Free Electrons and its Fluctuations_, 
arXiv e-prints, astro-ph/0207156.

[Cordes, J. M., & Lazio, T. J. W. (2003)](https://ui.adsabs.harvard.edu/abs/2003astro.ph..1598C/abstract),  
_NE2001. II. Using Radio Propagation Data to Construct a Model for the Galactic Distribution of Free Electrons_, 
arXiv e-prints, astro-ph/0301598.

[Yao, J. M., Manchester, R. N., & Wang, N. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract),  
_A New Electron-density Model for Estimation of Pulsar and FRB Distances_, 
The Astrophysical Journal, Volume 888, Issue 2, id.105, Colume 835, id.29

[Yamasaki, S., & Totani, T. (2020)](https://ui.adsabs.harvard.edu/abs/2019arXiv190900849Y/abstract),  
_The Galactic Halo Contribution to the Dispersion Measure of Extragalactic Fast Radio Bursts_, 
The Astrophysical Journal, Volume 888, Issue 2, id.105

These are available in bibtex format in [references.bib](https://github.com/telegraphic/pygedm/references.bib),
and also as an [ADS library](https://ui.adsabs.harvard.edu/public-libraries/Ci6_0-TlSySPMLrHxTvhhw). 

## YMW16 C README

YMW16 is a model for the distribution of free electrons in the Galaxy,
the Magellanic Clouds and the inter-galactic medium, that can be used
to estimate distances for real or simulated pulsars and fast radio
bursts (FRBs) based on their position and dispersion measure.

The Galactic model is based on 189 pulsars that have independently
determined distances as well as dispersion measures, whereas simpler
models are used for the electron density in the MC and the IGM. It is
estimated that the 95% of predicted Galactic pulsar distances will
have a relative error of less than a factor of 0.9. Pulse broadening
due to scattering in the Galactic interstellar medium, the Magellanic
Clouds, the intergalactic medium and FRB host galaxies is estimated.

As well as the ymw16 dm-distance program, we also provide a program,
ymw16_ne, which gives the electron density at any point in the Galaxy
or Magellanic Clouds.

A paper (Yao, Manchester and Wang, 2017,
[Astrophys. J., 835, 29](https://iopscience.iop.org/article/10.3847/1538-4357/835/1/29/meta);
[arXiv:1610.09448](https://arxiv.org/abs/1610.09448)) describes the model and compares its predictions
with those of earlier Galactic electron density models. YMW16 is the
first electron-density model to estimate extragalactic pulsar
distances and FRB distances.

To make a command-line executable version of the program, download and
unpack the latest version of the program. Then run "make_ymw16" to
create the executable. The environment variable YMW16_DIR may be set
up to point at a directory containing ymw16par.txt and
spiral.txt. Access to these files is needed at runtime.

Websites allowing interactive access to the YMW16 distance model and
download of the latest program version are available at:
* [http://www.xao.ac.cn/ymw16/](http://www.xao.ac.cn/ymw16/),
* [http://www.atnf.csiro.au/research/pulsar/ymw16/](http://www.atnf.csiro.au/research/pulsar/ymw16/]) and
* [https://bitbucket.org/psrsoft/ymw16/](https://bitbucket.org/psrsoft/ymw16/).

Please report any issues or bugs at
https://bitbucket.org/psrsoft/ymw16/issues/new/ or directly to the
authors. Please provide an example illustrating the problem.

### YMW16 C LICENSE

```
Copyright (C) 2016, 2017  J. M. Yao, R. N. Manchester, N. Wang.

YMW16 is free software: you can redistribute it and/or modify　it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

YMW16 is distributed in the hope that it will be useful,　but WITHOUT
ANY WARRANTY; without even the implied warranty of　MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the　GNU General Public License,
available at http://www.gnu.org/licenses/, for more details.

Jumei Yao (yaojumei _@_ xao.ac.cn), Richard N Manchester
(dick.manchester _@_ csiro.au), Na Wang (na.wang _@_ xao.ac.cn)
```

## NE2001 README

07 July 2002
To compile and execute the code,  see [code.pdf](https://github.com/telegraphic/pygedm/blob/master/ne2001_src/code.pdf).


