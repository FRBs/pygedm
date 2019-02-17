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

extern int m_3, ww,m_5, m_6, m_7;

void nps(double xx,double yy,double zz,double *ne7, int *WLI, struct LI t7)
{
  double x_c, y_c, z_c;
  double gLI;
  double theta_LI;

  if(m_7>=1)return;

  theta_LI=(t7.thetaLI)/RAD;
  x_c=-10.156;
  y_c=8106.206;
  z_c=10.467;
  double rr,theta;
  rr=sqrt((xx-x_c)*(xx-x_c)+(yy-y_c)*(yy-y_c)+(zz-z_c)*(zz-z_c));
  theta=acos(((xx-x_c)*(cos(theta_LI))+(zz-z_c)*(sin(theta_LI)))/rr)*RAD;
  *WLI=1;
  if(fabs(rr-t7.RLI)>(mc*t7.WLI)||fabs(theta)>(mc*t7.detthetaLI))
  { if(rr>500)m_7++;
  	*ne7=0;
  	return;
  }
  else gLI=1;
  *ne7=gLI*(t7.nLI)*exp(-((rr-(t7.RLI))*(rr-(t7.RLI)))/((t7.WLI)*(t7.WLI)))*exp((-theta*theta)/((t7.detthetaLI)*(t7.detthetaLI)));
}
