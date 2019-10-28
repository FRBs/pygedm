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

void thick(double xx, double yy, double zz, double *gd, double *ne1, double rr, struct Thick t1){
  
  double gdd,gg;

  if(fabs(zz)> mc*t1.H1 || (rr-t1.Bd)> mc*t1.Ad){
    *ne1=0;
    return;
  }else{
    if(rr<t1.Bd){
      gdd=1;
    }else{ 
      gg=exp(-(rr-t1.Bd)/t1.Ad)+exp((rr-t1.Bd)/t1.Ad);
      gdd=pow(2/gg,2);
    }
  }
  *ne1=t1.n1*gdd*pow(2/(exp(-fabs(zz)/t1.H1)+exp(fabs(zz)/t1.H1)), 2);
  *gd=gdd;
}
