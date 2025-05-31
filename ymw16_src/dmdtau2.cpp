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


extern int m_3, ww, m_5, m_6, m_7;

double tsc(double dm){
  return 4.1e-11*pow(dm, 2.2)*(1+0.00194*dm*dm);
}

// Function returns to Python dictionary via std::map
std::map<std::string, float> dmdtau2(double gl, double gb ,double dordm, double DM_Host, int ndir, int np, int vbs, char *dirname, char *text)
{
  std::map<std::string, float> result;
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
  double dist, xx, yy, zz, r, glr, gbr, sl, cl, sb, cb, hh;

  double nstep, dstep, dmstep;
  static double dd, dtest, dmpsr, rr;
  double dmm=0;
  double dm=0;
  double DM_MC=0;
  double DM_Gal=0;
  double DDM;
  double tau_sc=0;
  double tau_Gal=0;
  double tau_MC=0;
  double tau_MC_sc=0;
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
  int nk, uu; // DCP 2019.02.17 - Removed unusued nn variable

  static int i, ncount;
  int w_lmc=0;
  int w_smc=0;
  int umc=1;


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
  m_3=0;
  ww=1;
  m_5=0;
  m_6=0;
  m_7=0;

  ymw16par(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10, &t11, dirname);

  glr=gl/RAD;
  gbr=gb/RAD;
  sl=sin(glr);
  sb=sin(gbr);
  cl=cos(glr);
  cb=cos(gbr);
  dstep=5.0;

  if(np==-1){                 // FRBs
    if(ndir==1)uu=0;//dm---dist
    else uu=1;//dist---dm
    ndir=2;
    DDM=dordm;
    dordm=100000;
    nk=20000;
  }

  if(np==0){                 // Magellanic Cloud
    nk=20000;
  }

  if(np==1){                  //Galactic pulsars
    nk=5000;
  }
  if(ndir==1){
    dm=dordm;
    if(np==1)tau_sc=tsc(dm);
    dtest=dm/N0;
    nstep=dtest/dstep;
    if(nstep<200) dstep=dtest/200;
    if(vbs>=1){
      // printf("\ndtest=%lf, nstep=%lf, dstep=%lf\n", dtest, nstep, dstep);
    }
  }
  if(ndir==2){
    dist=dordm;
    dtest=dist;
    nstep=dtest/dstep;
    if(nstep<200) dstep=dtest/200;
    if(vbs>=1){
      // printf("\ndtest=%lf, nstep=%lf, dstep=%lf\n", dtest, nstep, dstep);
    }
  }


  dd=-0.5*dstep;
  ncount=0;

  for(i=1;i<=nk;i++){
    ncount++;
    if(vbs>=2){
      // printf("ncount=%d, dstep=%lf\n", ncount,dstep);
    }
    dd+=dstep;
    r=dd*cb;     /* r is different from rr */
    xx=r*sl;
    yy=R0*1000-r*cl;
    zz=dd*sb+t0.z_Sun;
    rr=sqrt(xx*xx+yy*yy);

    /* Definition of warp */
    if(rr<R_warp){
      zz_w=zz;
    }else{
      theta_warp=atan2(yy,xx);
      z_warp=t0.Gamma_w*(rr-R_warp)*cos(theta_warp-theta_max);
      zz_w=zz-z_warp;
    }

    if(vbs>=2)
    {
      // printf("dd=%lf, xx=%lf, yy=%lf, zz=%lf, rr=%lf\n", dd, xx, yy, zz, rr);
      // printf("theta_warp=%lf, z_warp=%lf, zz_w=%lf\n",theta_warp,z_warp,zz_w);
    }
    R_g=sqrt(xx*xx+yy*yy+zz*zz);

    /* DM to Distance */

    if(ndir==1){
      if(dmm<=dm){
        if(R_g<=35000){
	  thick(xx, yy, zz_w, &gd, &ne1, rr, t1);
          thin(xx, yy, zz_w, gd, &ne2, rr, t2);
          spiral(xx, yy, zz_w, gd, &ne3, rr, t3, dirname);
          galcen(xx, yy, zz, &ne4, t4);
          gum(xx, yy, zz, &ne5, t5);
          localbubble(xx, yy, zz, gl, gb, &ne6, &hh, t6);
          nps(xx, yy, zz, &ne7, &WLI, t7);
          fermibubble(xx, yy, zz, &WFB);
        }else{
          if(np==1){
	    dstep=5;
          }else{
            dstep=200;
            if(w_lmc>=1||w_smc>=1) dstep=5;
            lmc(glr,gbr,dd,&w_lmc,&ne8,t9);
            dora(glr,gbr,dd,&ne9,t10);
            smc(xx, yy, zz,&w_smc, &ne10, t11);
          }
	}
	if(WFB==1){
	  ne1=t8.J_FB*ne1;
	}
	ne0=ne1+MAX(ne2,ne3);

	if(hh>110){       /* Outside LB */
	  if(ne6>ne0 && ne6>ne5){
	    WLB=1;
	  }else{
	    WLB=0;
	  }
	}else{            /* Inside LB */
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
	ne=(1-WLB)*((1-WGN)*((1-WLI)*(ne0+ne4+ne8+ne9+ne10)+WLI*ne7)+WGN*ne5)+WLB*ne6;

	if(vbs>=2){
	  // printf("ne=%lf, ne1=%lf, ne2=%lf, ne3=%lf, ne4=%lf, ne5=%lf, ne6=%lf, ne7=%lf, ne8=%lf, ne9=%lf, ne10=%lf\n", ne, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10);
	}
	dmstep=ne*dstep;
	if(dmstep<=0.000001)dmstep=0;
	if(vbs>=2){
	  // printf("dmstep=%lf, dstep=%lf\n", dmstep, dstep);
	}
	dmm+=dmstep;
	dist=dd;
	if(np==0&&umc==1){
	  if(R_g>35000){
	    DM_Gal=dmm;
	    tau_Gal=0.5*tsc(dmm);
	    //printf(" DM_Gal:%8.2f",DM_Gal);
      result.insert(std::make_pair("DM_Gal", DM_Gal));
	    umc++;
	  }
	}
	if(i==nk){
	  dist+=0.5*dstep;
	  if(dist>100000)dist=100000;
	  if(np==0){
	    DM_MC=dmm-DM_Gal;
	    tau_MC=0.5*tsc(DM_MC);
	    tau_MC_sc=MAX(tau_Gal, tau_MC);
	    //printf(" DM_MC:%8.2f",DM_MC);
      result.insert(std::make_pair("DM_MC", DM_MC));
	  }
	  if(np==0){
      //printf(" Dist:%9.1f log(tau_sc):%7.3f %s\n",dist, log10(tau_MC_sc),text);
      result.insert(std::make_pair("dist", dist));
      result.insert(std::make_pair("tau_sc", tau_MC_sc));
    }
	  if(np==1){
      //printf(" DM_Gal:%8.2f Dist:%9.1f log(tau_sc):%7.3f %s\n", dmm, dist,log10(tau_sc),text);
      result.insert(std::make_pair("DM_Gal", dmm));
      result.insert(std::make_pair("dist", dist));
      result.insert(std::make_pair("tau_sc", tau_sc));

    }
	}
        if(vbs>=2){
	  // printf("dmm=%lf\n", dmm);
	}
      }
      else{
	dist=dd-0.5*dstep-(dstep*(dmm-dm))/dmstep;
	if(np==0){
	  DM_MC=dm-DM_Gal;
          if(DM_Gal==0){
            DM_MC=0;
            DM_Gal=dm;
            tau_MC_sc=tsc(dm);
            // printf(" DM_Gal:%8.2f ", DM_Gal);
            result.insert(std::make_pair("DM_Gal", DM_Gal));
          }
          else{
	    tau_MC=0.5*tsc(DM_MC);
            tau_MC_sc=MAX(tau_MC, tau_Gal);
          }
          // printf(" DM_MC:%8.2f", DM_MC);
          result.insert(std::make_pair("DM_MC", DM_MC));
        }
	if(np==0){
    //printf(" Dist:%9.1f log(tau_sc):%7.3f %s\n",dist,log10(tau_MC_sc),text);
    result.insert(std::make_pair("dist", dist));
    result.insert(std::make_pair("tau_sc", tau_MC_sc));
  }
	if(np==1){
    //printf(" DM_Gal:%8.2f Dist:%9.1f log(tau_sc):%7.3f %s\n", dm, dist,log10(tau_sc),text);
    result.insert(std::make_pair("DM_Gal", DM_Gal));
    result.insert(std::make_pair("dist", dist));
    result.insert(std::make_pair("tau_sc", tau_sc));
  }
	break;
      }
    }

    /* Distance to DM */

    if(ndir==2){
      if(dd<=dtest){
        if(R_g<=35000){
	      thick(xx, yy, zz_w, &gd, &ne1, rr, t1);
          thin(xx, yy, zz_w, gd, &ne2, rr, t2);
          spiral(xx, yy, zz_w, gd, &ne3, rr, t3, dirname);
          galcen(xx, yy, zz, &ne4, t4);
          gum(xx, yy, zz, &ne5, t5);
          localbubble(xx, yy, zz, gl, gb, &ne6, &hh, t6);
          nps(xx, yy, zz, &ne7, &WLI, t7);
          fermibubble(xx, yy, zz, &WFB);
        }else{
	  if(np==1)dstep=5;
	  else{
	    dstep=200;
	    if(np==-1)dstep=5;
	    if(w_lmc>=1||w_smc>=1) dstep=5;
	    lmc(glr,gbr,dd,&w_lmc,&ne8,t9);
	    dora(glr,gbr,dd,&ne9,t10);
	    smc(xx, yy, zz,&w_smc, &ne10, t11);
	  }
	}
	if(WFB==1){
	  ne1=t8.J_FB*ne1;
	}
	ne0=ne1+MAX(ne2,ne3);

        if(hh>110){       /* Outside LB */
          if(ne6>ne0 && ne6>ne5){
	    WLB=1;
          }else{
	    WLB=0;
	  }
	}else{            /* Inside LB */
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
	/*  Galactic ne */
        ne=(1-WLB)*((1-WGN)*((1-WLI)*(ne0+ne4+ne8+ne9+ne10)+WLI*ne7)+WGN*ne5)+WLB*ne6;

	if(vbs>=2){
          // printf("ne=%lf, ne1=%lf, ne2=%lf, ne3=%lf, ne4=%lf, ne5=%lf, ne6=%lf, ne7=%lf, ne8=%lf, ne9=%lf, ne10=%lf\n", ne, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10);
        }
	dmstep=ne*dstep;
	if(dmstep<=0.000001)dmstep=0;
	dm+=dmstep;
	if(np!=1&&umc==1){
	  if(R_g>35000){
	    DM_Gal=dm;
	    tau_Gal=0.5*tsc(dm);
	    //printf(" DM_Gal:%8.2f",dm);
      result.insert(std::make_pair("DM_Gal", DM_Gal));
	    umc++;
	  }
        }
        if(i==nk&&np!=-1){
          dmpsr=dm;
          if(np==0){
            DM_MC=dm-DM_Gal;
            tau_MC=0.5*tsc(DM_MC);
            //printf(" DM_MC:%8.2f", DM_MC);
            result.insert(std::make_pair("DM_MC", DM_MC));
          }
          tau_sc=tsc(dmpsr);
          tau_MC_sc=MAX(tau_Gal, tau_MC);
          if(np==0){
            //printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr,log10(tau_MC_sc),text);
            result.insert(std::make_pair("DM", dmpsr));
            result.insert(std::make_pair("tau_sc", tau_MC_sc));
          }
          if(np==1){
            //printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr, log10(tau_sc),text);
            result.insert(std::make_pair("DM", dmpsr));
            result.insert(std::make_pair("tau_sc", tau_sc));
          }
        }

	if(i==nk&&np==-1){
          if(dordm==100000){
	    DM_MC=dm-DM_Gal;
	    //printf(" DM_MC:%8.2f",DM_MC);
	  }
          result = frb_d(DDM, DM_Gal, DM_MC, DM_Host, uu, vbs, text);
          result.insert(std::make_pair("DM_MC", DM_MC));
          break;
        }
      }
      else{
	dmpsr=dm+(dmstep*(dtest-(dd-0.5*dstep)))/dstep;
	if(np==0){
	  DM_MC=dmpsr-DM_Gal;
	  if(DM_Gal==0){
	    DM_MC=0;
	    DM_Gal=dmpsr;
	    tau_MC_sc=tsc(dmpsr);
	    //printf(" DM_Gal:%8.2f", DM_Gal);
      result.insert(std::make_pair("DM_Gal", DM_Gal));
	  }
	  else{
	    tau_MC=0.5*tsc(DM_MC);
	    tau_MC_sc=MAX(tau_Gal, tau_MC);
	  }
	  //printf(" DM_MC:%8.2f", DM_MC);
    result.insert(std::make_pair("DM_MC", DM_MC));
	}
	tau_sc=tsc(dmpsr);
	if(np==0){
    //printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr,log10(tau_MC_sc),text);
    result.insert(std::make_pair("DM", dmpsr));
    result.insert(std::make_pair("tau_sc", tau_MC_sc));
  }
	if(np==1){
    //printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr, log10(tau_sc),text);
    result.insert(std::make_pair("DM", dmpsr));
    result.insert(std::make_pair("tau_sc", tau_sc));
  }

	break;
      }
    }
  }

  return result;
}
