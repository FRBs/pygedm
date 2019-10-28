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
void thin(double xx, double yy, double zz, double gd, double *ne2, double rr, struct Thin t2)
{
  double g2, Hg, ex1, ex2, g3, HH;
  Hg=32+0.0016*rr+0.0000004*pow(rr, 2);
  HH=t2.K2*Hg;
  if((rr-t2.B2)>(mc*t2.A2)||(fabs(zz)>(mc*HH))) 
  {
  	*ne2=0;
  	return;	
  }
  else
  {
    g3=(rr-t2.B2)/t2.A2;
    g2=pow(2/(exp(-g3)+exp(g3)), 2);
  }
  *ne2=t2.n2*gd*g2*pow(2/(exp(-zz/HH)+exp(zz/HH)), 2);
}
