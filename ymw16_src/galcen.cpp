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
void galcen(double xx, double yy, double zz, double *ne4 ,struct GC t4)
{
  double Xgc=50;
  double Ygc=0;
  double Zgc=-7;
  double Rgc, Az, Ar;
  double gc;
  Rgc=sqrt((xx-Xgc)*(xx-Xgc)+(yy-Ygc)*(yy-Ygc));
  if(Rgc>(mc*t4.Agc)||fabs(zz)>(mc*t4.Hgc)) 
  {
  	*ne4=0;
  	return;
  }
  else
  {
    gc=1;
    Ar=exp(-(Rgc*Rgc)/(t4.Agc*t4.Agc));
    Az=pow(2/(exp((zz-Zgc)/t4.Hgc)+exp(-(zz-Zgc)/t4.Hgc)),2);
  }
  *ne4=t4.ngc*Ar*Az*gc;
}


