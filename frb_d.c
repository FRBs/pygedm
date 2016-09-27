#include "cn.h"
#define max(a,b) ( ((a)>(b)) ? (a):(b) )
void frb_d(double DDM, double DM_Gal, double DM_MC, double DM_Host, int uu)
{

//parameters for IG
  double c=3e8;//m/s
  double H0=2.1808e-18;//1/s
  double nIG=0.16;//m^-3
  double dmtol,z,DM_IGM;
  double dist;
  double tau_Gal=0;
  double tau_MC=0;
  double tau_IGM=0;
  double tau_Host=0;
  double tau_FRB=0;

  if(uu==1){                      //Dist---DM 
    DDM=DDM*3.086e22;
    z=exp(DDM*(H0/c))-1;
    dmtol=(z*c*nIG)/(3.086e22*H0)+DM_Gal+DM_MC+DM_Host;
    if(dmtol<=(DM_Gal+DM_MC+DM_Host))
    {
      DM_IGM=0;
      printf(" DM_IGM:%8.2f\n", DM_IGM);
      printf("DM_tol < (DM_Gal+DM_MC+DM_Host)\n");
    }
    else
    {
      DM_IGM=dmtol-DM_Gal-DM_MC-DM_Host;
      tau_Gal=0.5*tsc(DM_Gal); 
      tau_MC=0.5*tsc(DM_MC);  
      tau_IGM=pow(10,-6.4)*pow(DM_IGM, 1.3);        
      tau_Host=0.5*tsc(DM_Host);  
      tau_FRB=max(tau_Gal,tau_MC); 
      tau_FRB=max(tau_FRB,tau_IGM); 
      tau_FRB=max(tau_FRB,tau_Host); 
      printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
      printf(" DM:%8.2f log(tau_sc):%7.3f\n", dmtol, log10(tau_FRB));
    }
    
    //    printf("tsc_gal %12.2e tsc_MC %12.2e tsc_IGM %12.2e tsc_Host %12.2e\n",tau_Gal,tau_MC,tau_IGM,tau_Host);
 }
  if(uu==0)
  {                         //DM---Dist
    dmtol=DDM;
    z=(3.086e22*(dmtol-DM_Gal-DM_MC-DM_Host)*H0)/(c*nIG);
    if(dmtol<=(DM_Gal+DM_MC+DM_Host))
    {
      DM_IGM=0;
      printf(" DM_IGM:%8.2f DM_Host:%8.2f\n", DM_IGM,DM_Host);
      printf("Warning: DM < (DM_Gal+DM_MC+DM_Host), Dist=0.0)\n");
    }
    else
    {
      DM_IGM=dmtol-DM_Gal-DM_MC-DM_Host;
      tau_Gal=0.5*tsc(DM_Gal); 
      tau_MC=0.5*tsc(DM_MC);  
      tau_IGM=pow(10,-6.4)*pow(DM_IGM, 1.3);
      tau_Host=0.5*tsc(DM_Host);
      tau_FRB=max(tau_Gal,tau_MC); 
      tau_FRB=max(tau_FRB,tau_IGM); 
      tau_FRB=max(tau_FRB,tau_Host); 
      dist=(c/H0)*log(1+z);
      dist=dist/3.086e22;
      if(z<=0||dist<=0)
      {
        z=0;
        dist=0;
      }  
      printf(" DM_IGM:%8.2f DM_Host:%8.2f z:%7.3f", DM_IGM, DM_Host, z);
      printf(" Dist:%8.1f log(tau_sc):%7.3f\n", dist, log10(tau_FRB)); //Mpc 
    }
    
    //    printf("tsc_gal %12.2e tsc_MC %12.2e tsc_IGM %12.2e tsc_Host %12.2e\n",tau_Gal,tau_MC,tau_IGM,tau_Host);
  }
}
