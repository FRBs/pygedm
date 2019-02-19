[![Build Status](https://travis-ci.org/telegraphic/pyymw16.svg?branch=master)](https://travis-ci.org/telegraphic/pyymw16)
[![Coverage Status](https://coveralls.io/repos/github/telegraphic/pyymw16/badge.svg?branch=master)](https://coveralls.io/github/telegraphic/pyymw16?branch=master)

# PyYMW16
_A Python / C++ Version of YMW16 electron-density model_

This is a Python / C++ port of the Yao, Manchester and Wang (2017, [Astrophys. J., 835, 29](https://iopscience.iop.org/article/10.3847/1538-4357/835/1/29/meta);
[arXiv:1610.09448](https://arxiv.org/abs/1610.09448)) YMW16 electron density model.
The code uses [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html)
to create Python bindings to the (C++ ported) YMW16 code.

### Installation

Requires pybind11 and C++11 (i.e. gcc>4.9, Ubuntu 16.04 should work).

```
python setup.py install
```

and test

```
python setup.py test
```

### Todo

* Create nice python wrapper using `argparse`

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
