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
#define Rlb 110
extern int m_3, ww,m_5, m_6, m_7;
void localbubble(double xx, double yy, double zz, double gl,double gb, double *ne6, double *WW, struct LB t6)
{
  double g4=0;
  double g5=0;
  double g6=0;
  double g7=0;
  double glb1, glb2;
  double UU;
  double VVV=0;
  double WWW=0;
  double nelb1=1e10;
  double nelb2=1e10;
  double Y_C=1e10;
  double COS_A=1e10;
  double COS_E=cos(PI/4);
//  printf("m_6=%d\n\n",m_6);

  Y_C=8340+(0.34/0.94)*zz;//Center of local bubble
  COS_A=(xx*xx+(8300-yy)*(Y_C-yy))/(sqrt(xx*xx+(8300-yy)*(8300-yy))*sqrt(xx*xx+(Y_C-yy)*(Y_C-yy)));
  UU=sqrt(((yy-8340)*0.94-0.34*zz)*((yy-8340)*0.94-0.34*zz)+xx*xx);
  *WW=UU;
  if(m_6>=1)return;
  if((UU-Rlb)>(mc*t6.wlb1)||fabs(zz)>(mc*t6.hlb1)){
    nelb1=0;
  }else{
    glb1=1;
    g4=(UU-Rlb)/t6.wlb1;
    g5=pow(2/(exp(-zz/t6.hlb1)+exp(zz/t6.hlb1)), 2);
    VVV=MIN(fabs(gl+360-t6.thetalb1),fabs(t6.thetalb1-gl));
    if(VVV>(mc*t6.detlb1)) glb1=0;
    else glb1=1;
    nelb1=glb1*pow(2/(exp(-VVV/t6.detlb1)+exp(VVV/t6.detlb1)), 2)*t6.nlb1*pow(2/(exp(g4)+exp(-g4)), 2)*g5;
  }
  if((UU-Rlb)>(mc*t6.wlb2)||fabs(zz)>(mc*t6.hlb2)) nelb2=0;
  if(nelb1==0&&nelb2==0){
    *ne6=0;
    m_6++;
    return;
  }else{
    glb2=1;
    g6=(UU-Rlb)/t6.wlb2;
    g7=pow(2/(exp(-zz/t6.hlb2)+exp(zz/t6.hlb2)), 2);
    WWW=MIN(fabs(gl+360-t6.thetalb2),fabs(t6.thetalb2-gl));
    if(WWW>(mc*t6.detlb2)) glb2=0;
    else glb2=1;
    nelb2=glb2*pow(2/(exp(-WWW/t6.detlb2)+exp(WWW/t6.detlb2)), 2)*t6.nlb2*pow(2/(exp(g6)+exp(-g6)), 2)*g7;
  }

  if(COS_A<COS_E)nelb1=0;
  if((gb==90)&&(gl==180))nelb1=0;
  *ne6=nelb1+nelb2;
}
