#include "cn.h"
#define max(a,b) ( ((a)>(b)) ? (a):(b) )
void frb_d(double DDM, double DM_Gal, double DM_MC, double DM_Host, int uu, int vbs, char* text)
{

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
    printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
    printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", DDM, log10(tau_FRB),text);
    if(vbs>=1)printf("tsc_gal:%10.2e, tsc_MC:%10.2e, tsc_IGM:%10.2e, tsc_Host:%10.2e,\n",
       tau_Gal,tau_MC,tau_IGM,tau_Host); 
  }

  if(uu==0){                         //DM---Dist
    DM_IGM=DDM-(DM_Gal+DM_MC+DM_Host);
    if(DM_IGM<0.0){
      DM_IGM=0.0;
      dist=0.0;
      z=0.0;
      printf(" DM_IGM:%8.2f DM_Host:%8.2f\n", DM_IGM, DM_Host);
      printf("Warning: DM < (DM_Gal+DM_MC+DM_Host), Dist=0.0)\n");
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
      printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
      printf(" Dist:%8.1f log(tau_sc):%7.3f %s\n", dist, log10(tau_FRB), text); //Mpc 
    }
    if(vbs>=1)printf("tsc_gal:%10.2e, tsc_MC:%10.2e, tsc_IGM:%10.2e, tsc_Host:%10.2e,\n",tau_Gal,tau_MC,tau_IGM,tau_Host); 
  }
}
