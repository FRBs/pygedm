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
void gum(double xx, double yy, double zz, double *ne5, struct Gum t5)
{
  double slc, clc, sbc, cbc, xc, yc, zc, rgalc;
  double rp, RR, xyp, zp;
  double theta, alpha;
  double Dmin=1e5;

  //center of Gum Nebula
  const double lc=264;
  const double bc=-4;
  const double dc=450;

  if(m_5>=1)return;

  slc=sin(lc/RAD);
  clc=cos(lc/RAD);
  sbc=sin(bc/RAD);
  cbc=cos(bc/RAD);

  rgalc=dc*cbc;
  xc=rgalc*slc;
  yc=R0*1000-rgalc*clc;
  zc=dc*sbc;

  theta=fabs(atan((zz-zc)/sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc))));
  zp=((t5.Agn)*(t5.Agn)*(t5.Kgn))/sqrt(((t5.Agn)*(t5.Agn))+((t5.Agn)*(t5.Agn)*(t5.Kgn)*(t5.Kgn))/(tan(theta)*tan(theta)));
  xyp=zp/tan(theta);
  if((t5.Agn-fabs(xyp))<1e-15)alpha=PI/2;
  else alpha=-atan((-(t5.Agn)*(t5.Kgn)*xyp)/((t5.Agn)*sqrt((t5.Agn)*(t5.Agn)-xyp*xyp)));
  RR=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
  rp=sqrt((zp)*(zp)+(xyp)*(xyp));
  Dmin=fabs((RR-rp)*sin(theta+alpha));

  if(Dmin>(mc*t5.Wgn)){
    if(RR>500)m_5++;
    *ne5=0;
    return;
  }
  *ne5=(t5.ngn)*exp(-Dmin*Dmin/((t5.Wgn)*(t5.Wgn)));
}
