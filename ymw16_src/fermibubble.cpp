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
void fermibubble(double xx, double yy, double zz, int *wfb)
{
  double fbnz, fbsz, na, nb, sa, sb, N, S;
//center of Fermi bubble
  fbnz=0.5*8300*tan(50/RAD);
  fbsz=-(0.5*(8300*tan(50/RAD)-8300*tan(0/RAD))+(8300*tan(0/RAD)));
//min_axis and max_axis of Fermi bubble
  na=fbnz;
  nb=8300*tan(20/RAD);
  sa=0.5*(8300*tan(50/RAD)-8300*tan(0/RAD));
  sb=8300*tan(20/RAD);
  N=(pow(xx, 2)/(nb*nb))+(pow(yy, 2)/(nb*nb))+(pow(zz-fbnz, 2)/(na*na));
  S=(pow(xx, 2)/(sb*sb))+(pow(yy, 2)/(sb*sb))+(pow(zz-fbsz, 2)/(sa*sa));
  if(N<1||S<1)*wfb=1;
  else *wfb=0;
}
