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
void lmc(double l, double b, double d, int *w_lmc, double *ne8, struct LMC t9)
{
//coordinate system
  double X, Y, Z;
//parameter
  double dag, pag, npag, iag;
  npag=(116+90)/RAD;
  iag=32/RAD;
  double l_cp, alpha_cp, delta_cp, alpha, delta;
//constant
  l_cp=122.9/RAD;
  alpha_cp=192.85/RAD;
  delta_cp=27.13/RAD;
//parameter of lmc
  double alpha_l=81/RAD;
  double delta_l=-69.75/RAD;
  double D_l=49700; //pc
  double Rl;
  double Al=3000;//pc
  double Hl=800;//pc
  int gl, gD;
  double s_delta, c_delta, s_alpha_cp, c_alpha_cp;
  double s_cpl,c_cpl, c_dag,s_dag,sc_dp,ss_dp;
//(gl,gb,d)-----get(alpha,delta,d) for a point
  s_delta=sin(b)*sin(delta_cp)+cos(b)*cos(delta_cp)*cos(l_cp-l);
  c_delta=cos(asin(s_delta));
  s_alpha_cp=(sin(l_cp-l)*cos(b))/c_delta;
  c_alpha_cp=(sin(b)*cos(delta_cp)-cos(b)*sin(delta_cp)*cos(l_cp-l))/c_delta;
  s_cpl=sin(alpha_cp-alpha_l);
  c_cpl=cos(alpha_cp-alpha_l);
// ((alpha,delta,d)----((alpha_l,delta_l,D_l)-----get(dag,pag)
  c_dag=c_delta*cos(delta_l)*(c_alpha_cp*c_cpl-s_alpha_cp*s_cpl)+s_delta*sin(delta_l);
  s_dag=sin(acos(c_dag));
  sc_dp=-c_delta*(s_alpha_cp*c_cpl+c_alpha_cp*s_cpl);
  ss_dp=s_delta*cos(delta_l)-c_delta*sin(delta_l)*(c_alpha_cp*c_cpl-s_alpha_cp*s_cpl);
//(xc,yc,zc,iag,npag)------get(X,Y,Z) is (x',y',z')
  X=d*(sc_dp*cos(npag)+ss_dp*sin(npag));
  Y=d*(cos(iag)*(ss_dp*cos(npag)-sc_dp*sin(npag))+c_dag*sin(iag))-D_l*sin(iag);
  Z=d*(sin(iag)*(ss_dp*cos(npag)-sc_dp*sin(npag))-c_dag*cos(iag))+D_l*cos(iag);
  Rl=sqrt(X*X+Y*Y);
  if(Rl>mc*Al||fabs(Z)>mc*Hl)
  {
  	*ne8=0;
  	return;
  }
  else 
  {
    gl=1;
    *w_lmc=1;

  }
  *ne8=gl*(t9.nlmc)*exp(-(Rl*Rl)/(Al*Al))*pow(2/(exp(-Z/Hl)+exp(Z/Hl)),2);
}          

