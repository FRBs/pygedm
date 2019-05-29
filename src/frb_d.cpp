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

#define max(a,b) ( ((a)>(b)) ? (a):(b) )
std::map<std::string, float> frb_d(double DDM, double DM_Gal, double DM_MC, double DM_Host, int uu, int vbs, char* text)
{

std::map<std::string, float> result;

//parameters for IG
  double c=3e8;//m/s
  double H0=2.1808e-18;//1/s
  double nIG=0.16;//m^-3
  double z,DM_IGM;
  double dist;
  double tau_Gal=0;
  double tau_MC=0;
  double tau_IGM=0;
  double tau_Host=0;
  double tau_FRB=0;

  if(uu==1){                      //Dist---DM 
    dist=DDM*3.086e22;
    z=exp(dist*H0/c)-1;
    DM_IGM=z*c*nIG/H0/3.086e22;
    DDM=DM_IGM+DM_Gal+DM_MC+DM_Host;
    tau_Gal=0.5*tsc(DM_Gal); 
    tau_MC=0.5*tsc(DM_MC);  
    tau_IGM=pow(10,-6.4)*pow(DM_IGM, 1.3);        
    tau_Host=0.5*tsc(DM_Host*(1+z))*pow((1+z),-3.0);  
    tau_FRB=max(tau_Gal,tau_MC); 
    tau_FRB=max(tau_FRB,tau_IGM); 
    tau_FRB=max(tau_FRB,tau_Host); 
    //printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
    //printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", DDM, log10(tau_FRB),text);
    result.insert(std::make_pair("DM_IGM", DM_IGM));
    result.insert(std::make_pair("DM_Host", DM_Host));
    result.insert(std::make_pair("z", z));
    result.insert(std::make_pair("DDM", DDM));
    result.insert(std::make_pair("tau_FRB", tau_FRB));
    if(vbs>=1)printf("tsc_gal:%10.2e, tsc_MC:%10.2e, tsc_IGM:%10.2e, tsc_Host:%10.2e,\n",
       tau_Gal,tau_MC,tau_IGM,tau_Host); 
  }

  if(uu==0){                         //DM---Dist
    DM_IGM=DDM-(DM_Gal+DM_MC+DM_Host);
    if(DM_IGM<0.0){
      DM_IGM=0.0;
      dist=0.0;
      z=0.0;
      //printf(" DM_IGM:%8.2f DM_Host:%8.2f\n", DM_IGM, DM_Host);
      //printf("Warning: DM < (DM_Gal+DM_MC+DM_Host), Dist=0.0)\n");
      result.insert(std::make_pair("DM_IGM", DM_IGM));
      result.insert(std::make_pair("DM_Host", DM_Host));
      result.insert(std::make_pair("DM_too_small", 1));
    }else{
      z=DM_IGM*H0*3.086e22/(c*nIG);
      dist=(c/H0)*log(1+z)/3.086e22;
      
      tau_Gal=0.5*tsc(DM_Gal); 
      tau_MC=0.5*tsc(DM_MC);  
      tau_IGM=pow(10,-6.4)*pow(DM_IGM, 1.3);
      tau_Host=0.5*tsc(DM_Host*(1+z))/pow((1+z),3.0);
      tau_FRB=max(tau_Gal,tau_MC); 
      tau_FRB=max(tau_FRB,tau_IGM); 
      tau_FRB=max(tau_FRB,tau_Host); 
      //printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
      //printf(" Dist:%8.1f log(tau_sc):%7.3f %s\n", dist, log10(tau_FRB), text); //Mpc
      result.insert(std::make_pair("DM_IGM", DM_IGM));
      result.insert(std::make_pair("DM_Host", DM_Host));
      result.insert(std::make_pair("dist", dist));
      result.insert(std::make_pair("tau_FRB", tau_FRB));

    }
    if(vbs>=1)printf("tsc_gal:%10.2e, tsc_MC:%10.2e, tsc_IGM:%10.2e, tsc_Host:%10.2e,\n",tau_Gal,tau_MC,tau_IGM,tau_Host); 
  }
  return result;
}
