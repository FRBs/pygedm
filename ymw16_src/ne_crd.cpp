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

double ne_crd(double *x, double *y, double *z, double *gl, double *gb, double *dd, int ncrd, int vbs, char *dirname, char *text)
{

  double ne0=0;
  double ne=0;
  double ne1=0;
  double ne2=0;
  double ne3=0;
  double ne4=0;
  double ne5=0;
  double ne6=0;
  double ne7=0;
  double ne8=0;
  double ne9=0;
  double ne10=0;

  double xx, yy, zz, glr, gbr, dist;
  double x_s, y_s, z_s, ll, bb, hh, r, sl, cl, sb, cb;

  static double rr;
  double R_g=0;
  double gd=0;

  //The localtion of Sun relative to GP and Warp
  double z_warp, zz_w, R_warp, theta_warp, theta_max;
  R_warp=8400;//pc
  theta_max=0.0; //In +x direction

  int WGN=0;
  int WLB=0;
  int WLI=0;
  int WFB=0;
  int np=2;

  m_3=0;
  ww=1;
  m_5=0;
  m_6=0;
  m_7=0;

  //parameters of MC
  int w_lmc=0;
  int w_smc=0;


  struct Warp_Sun t0;
  struct Thick t1;
  struct Thin t2;
  struct Spiral t3;
  struct GC t4;
  struct Gum t5;
  struct LB t6;
  struct LI t7;
  struct FB t8;
  struct LMC t9;
  struct Dora t10;
  struct SMC t11;


  ymw16par(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10, &t11, dirname);

  if(ncrd==1){    /* (x,y,z) to (gl,gb,dist) and ne_c */
    xx=*x;
    yy=*y;
    zz=*z;
    /* (x_s, y_s, z_s) the heliocentric coordinates of the point */
    x_s=xx;
    y_s=yy-R0*1000;
    z_s=zz-t0.z_Sun;

    //D the heliocentric distance (pc)
    dist=sqrt(x_s*x_s+y_s*y_s+z_s*z_s);
    r=sqrt(x_s*x_s+y_s*y_s);


    if(dist<1e-15){
      gbr=0;
    }else{
      gbr=asin(z_s/dist);
    }
    *gb=gbr*RAD;

    if(r<1e-15){
      glr=0;
    }else{
      if(x_s>=0){
	glr=acos(-y_s/r);
      }else{
	glr=acos(y_s/r)+PI;
      }
    }
    *gl=glr*RAD;
    *dd=dist;
  }else{
    if(ncrd==2){       /* (gl,gb,dist) to (x,y,z) and ne_c */
      glr=*gl/RAD;
      gbr=*gb/RAD;
      dist=*dd;

      sl=sin(glr);
      sb=sin(gbr);
      cl=cos(glr);
      cb=cos(gbr);
      r=dist*cb;
      xx=r*sl;
      yy=R0*1000-r*cl;
      zz=dist*sb+t0.z_Sun;

      *x=xx;
      *y=yy;
      *z=zz;
  }else{
      return(0.0);
    }
  }

  rr=sqrt(xx*xx+yy*yy);

  /* Definition of warp */
  if(rr<R_warp){
    zz_w=zz;
  }
  else{
    theta_warp=atan2(yy,xx);
    z_warp=t0.Gamma_w*(rr-R_warp)*cos(theta_warp-theta_max);
    zz_w=zz-z_warp;
  }

  R_g=sqrt(xx*xx+yy*yy+zz*zz);

  if(vbs>=1)printf("Ne_crd: %10.1f %10.3f %10.3f\n",R_g,*gl,*gb);
  if(R_g<=30000){
    np=1;
  }else if(*gl>265. && *gl<315. && *gb>-60. && *gb<-20.) np=0;

  //printf("np=%d\n",np);

  if(R_g<100000){
    if(np==1){
      thick(xx, yy, zz_w, &gd, &ne1, rr, t1);
      thin(xx, yy, zz_w, gd, &ne2, rr, t2);
      spiral(xx, yy, zz_w, gd, &ne3, rr, t3, dirname);
      galcen(xx, yy, zz, &ne4, t4);
      gum(xx, yy, zz, &ne5, t5);
      localbubble(xx, yy, zz, *gl, *gb, &ne6, &hh, t6);
      nps(xx, yy, zz, &ne7, &WLI, t7);
      fermibubble(xx, yy, zz, &WFB);
      if(WFB==1){
	ne1=t8.J_FB*ne1;
      }
      ne0=ne1+MAX(ne2,ne3);

      if(hh>110){       /* Outside LB */
	if(ne6>ne0 && ne6>ne5){
	  WLB=1;
	}
	else{
	  WLB=0;
	}
      }
      else{            /* Inside LB */
	if(ne6>ne0){
	  WLB=1;
	}else{
	  ne1=t6.J_LB*ne1;
	  ne0=ne1+MAX(ne2,ne3);
	  WLB=0;
	}
      }
      if(ne7>ne0){     /* Loop I */
	WLI=1;
      }else{
	WLI=0;
      }
      if(ne5>ne0){     /* Gum Nebula */
	WGN=1;
      }else{
	WGN=0;
      }

      /* Galactic ne */
      ne=(1-WLB)*((1-WGN)*((1-WLI)*(ne0+ne4)+WLI*ne7)+WGN*ne5)+WLB*ne6;
      if(vbs>=1){
	printf("ne0=%lf, ne1=%lf, ne2=%lf, ne3=%lf, ne4=%lf, ne5=%lf, ne6=%lf, ne7=%lf\n", ne0, ne1, ne2, ne3, ne4, ne5, ne6, ne7);
	if(ne>1e-15) printf("The point is located within MW.\n");
      }
    }
    else if(np==0){
      lmc(glr,gbr,dist,&w_lmc,&ne8,t9);
      dora(glr,gbr,dist,&ne9,t10);
      smc(xx, yy, zz,&w_smc, &ne10, t11);
      /* MC ne */
      ne=ne8+ne9+ne10;
      if(vbs>=1){
	printf("ne_8=%lf, ne9=%lf, ne10=%lf\n", ne8, ne9, ne10);
	if(ne>1e-15) printf("The point is located within MC.\n");
      }
    }
    else{
      ne=0;
      if(vbs>=1) printf("The point is located outside of MW and MC.\n");
    }
  }
  else{
    ne=0;
    if(vbs>=1) printf("The point is located outside of MW and MC.\n");
  }
  return(ne);
}
