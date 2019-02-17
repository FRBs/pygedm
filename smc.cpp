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
#include"cn.hpp"
void smc(double xx, double yy, double zz, int  *w_smc, double *ne10, struct SMC t11)
{
  double rad=57.295779;
  double R=8.3;//kpc
  double ls=303.728914;
  double bs=-44.299212;
  double ds=59700;//pc
  double Asmc=3000;//pc
  int gsmc;
  double xc, yc, zc, sls, cls, sbs, cbs, rgals, Rsmc;
  sls=sin(ls/rad);
  cls=cos(ls/rad);
  sbs=sin(bs/rad);
  cbs=cos(bs/rad);
  rgals=ds*cbs;
  xc=rgals*sls;
  yc=R*1000-rgals*cls;
  zc=ds*sbs;
  Rsmc=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
  if(Rsmc>(mc*Asmc))
  {
  	*ne10=0;
  	return;
  }
  else 
  {
   gsmc=1;
   *w_smc=1; 
  }
    
  
  *ne10=(t11.nsmc)*gsmc*exp(-(Rsmc*Rsmc)/(Asmc*Asmc));
}

