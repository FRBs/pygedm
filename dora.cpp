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
void dora(double l, double b, double d, double *ne9, struct Dora t10)
{
  //coordinate system
  double X, Y, Z;
  //paramter
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
  int gl,gD;
//parameter for 30D
  double alpha_D=85/RAD;
  double delta_D=-69/RAD;
  double dag_D, pag_D, xc_D, yc_D, zc_D, X_D, Y_D, Z_D, R_D, A_D;
  A_D=450 ;//pc
  double D_D=49045;
  double s_delta, c_delta, s_alpha_cp, c_alpha_cp;
  double s_cpl,c_cpl, c_dag,s_dag,sc_dp,ss_dp;
  double c_dagD,s_dagD,sc_dpD,ss_dpD;
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
//Calculation for 30 Doradus
//1. ((alpha_D,delta_D,D_l)----((alpha_l,delta_l,D_l)-----get(dag_D,pag_D)
  c_dagD=cos(delta_D)*cos(delta_l)*cos(alpha_D-alpha_l)+sin(delta_D)*sin(delta_l);
  s_dagD=sin(acos(c_dagD));
  sc_dpD=-cos(delta_D)*sin(alpha_D-alpha_l);
  ss_dpD=sin(delta_D)*cos(delta_l)-cos(delta_D)*sin(delta_l)*cos(alpha_D-alpha_l);

//3. (xc_D,yc_D,zc_D,iag,npag)------get(X_D,Y_D,Z_D) is (x',y',z')
  X_D=D_D*(sc_dpD*cos(npag)+ss_dpD*sin(npag));
  Y_D=D_D*(cos(iag)*(ss_dpD*cos(npag)-sc_dpD*sin(npag))+c_dagD*sin(iag))-D_l*sin(iag);
  Z_D=D_D*(sin(iag)*(ss_dpD*cos(npag)-sc_dpD*sin(npag))-c_dagD*cos(iag))+D_l*cos(iag);
  R_D=sqrt((X-X_D)*(X-X_D)+(Y-Y_D)*(Y-Y_D)+(Z-Z_D)*(Z-Z_D));
  if(R_D>(mc*A_D))
  {
  	*ne9=0;
  	return;
  }
  else gD=1;
  *ne9=gD*(t10.n30D)*exp(-(R_D*R_D)/(A_D*A_D));
}
