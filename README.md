# PyYMW16
_A Python / C++ Version of YMW16 electron-density model_

**Update 2019.02.17:** First 'working' version (very beta). Don't expect
this to work out of the box.

This is a Python / C++ port of the Yao, Manchester and Wang (2017, [Astrophys. J., 835, 29](https://iopscience.iop.org/article/10.3847/1538-4357/835/1/29/meta);
[arXiv:1610.09448](https://arxiv.org/abs/1610.09448)) YMW16 electron density model.
The code uses [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html)
to create Python bindings to the (C++ ported) YMW16 code.

### Installation

Requires pybind11 and C++11. Should install (only tested on my mac)
with `./make_ymw16`. This will compile:

* `ymw16` - main program (output should be identical to C version)
* `ymw16_ne` - main program, electron model output
* `ymw16.so` - pybind shared object

Once compiled, you can import the `ymw16.so` from ipython, as if it
were a python module:

```python
ipython
> import ymw16
> a = ymw16.dmdtau(204, -6.5, 100000, 0, 2, 1, 0, './data', '')
DM:  252.05 log(tau_sc): -3.010
> print a
Out[4]: {u'DM': 252.05010986328125, u'tau_sc': 0.000978002673946321}
> ymw16.dmdtau?
Docstring:
dmdtau(gl: float, gb: float, dordm: float, DM_Host: float, ndir: int, np: int, vbs: int, dirname: unicode, text: unicode) -> Dict[unicode, float]

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
Type:      builtin_function_or_method

In [2]: ymw16.ne_crd?
Docstring:
ne_crd(x: float, y: float, z: float, gl: float, gb: float, dd: float, ncrd: int, vbs: int, dirname: unicode, text: unicode) -> float


Calculate electron density at a given point with galactocentric coordinates
(x, y, z) OR with (gl, gb, dist).

Args:
  (x, y, z): input Galactocentric x, y and z in pc
  (gl, gb, dist): input gl, gb in deg, Dist in pc
  ncrd: if ncrd==1, use xyz coords. If ncrd==2 use gl gb dist coords.
  vbs: Verbosity level, 0, 1, or 2
  dirname: directory where data files are stored
  text: Text to prepend in print statement.
Type:      builtin_function_or_method
```

### Todo

* Create nice python wrapper using `argparse`
* Setuptools installation (e.g. `pip install ymw16`)
* port `ymw16_ne`


## YMW16 C Code

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

A paper (Yao, Manchester and Wang, 2017, [Astrophys. J., 835, 29](https://iopscience.iop.org/article/10.3847/1538-4357/835/1/29/meta);
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
download of the latest program version are available at
http://www.xao.ac.cn/ymw16/,
http://www.atnf.csiro.au/research/pulsar/ymw16/ and
https://bitbucket.org/psrsoft/ymw16/.

Please report any issues or bugs at
https://bitbucket.org/psrsoft/ymw16/issues/new/ or directly to the
authors. Please provide an example illustrating the problem.

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
