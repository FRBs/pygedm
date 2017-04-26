#include "cn.h"
void spiral(double xx,  double yy,  double zz,  double gd, double *ne3,  double rr,  struct Spiral t3, char *filedir)
{
  int i, which_arm;
  static double rmin[5], thmin[5], tpitch[5], cspitch[5], sspitch[5];
  double detrr=1e10;
  double armr1, armr2, smin, sminmin, saxis, uu, Aaa, HH, Hg;
  double ne3s=0;
  double theta, alpha, armr;
  double sech2=0;
  double ga=0;
  double g1=0;
  char filen[64];
  FILE *fp;

  if(m_3>=1)return;

  Hg=32+0.0016*rr+0.0000004*pow(rr, 2);
  HH=t3.Ka*Hg;
  if(ww==1){
    strcpy(filen,filedir);
    strcat(filen,"spiral.txt");
    fp=fopen(filen,"r");
    
    for(i=0;i<=4;i++){
      fscanf(fp, "%lf %lf %lf %lf %lf", &rmin[i], &thmin[i], &tpitch[i], &cspitch[i], &sspitch[i]);
    }
    fclose(fp);
    ww++;
  }

  theta=atan2(yy,xx);
  if(theta<0)theta=2*PI+theta;
  
  //普通角度的计算
  if(fabs(zz/300)<10){
    sminmin=1e10;
    ne3s=0.;
    for(i=0;i<=4;i++){
      ga=0;
//Norma-Outer
      if(i==0)
      {
        if(theta>=0&&theta<0.77)
        {
          armr=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=fabs(rr-armr);
        }
        if(theta>=0.77&&theta<6.28)
        {
          armr1=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
      }
//Perseus
      if(i==1)
      {
        if(theta>=0&&theta<2.093)
        {
          armr=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=fabs(rr-armr);
        }
        if(theta>=2.093&&theta<6.28)
        {
          armr1=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
      }
//Carina-Sagittarius
      if(i==2)
      {
        if(theta>=0&&theta<3.81)
        {
          armr1=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+4*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
        if(theta>=3.81&&theta<6.28)
        {
          armr1=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
      }
//Crux_Scutum
      if(i==3)
      {
        if(theta>=0&&theta<5.76)
        {
          armr1=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+4*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
        if(theta>=5.76&&theta<6.28)
        {
          armr1=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
          armr2=rmin[i]*exp((theta+2*PI-thmin[i])*tpitch[i]);
          detrr=MIN(fabs(rr-armr1), fabs(rr-armr2));
        }
      }
//Local
      if(i==4)
      {
        if(theta>=0&&theta<0.96)
        {
          detrr=1e10;
        }
        if(theta>=0.96&&theta<2)
        {
          armr=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
          detrr=fabs(rr-armr);
        }
        if(theta>=2&&theta<6.28)
        {
          detrr=1e10;
        }
      }
      if(detrr>mc*t3.warm[i])
      {
      	ga=0;
      	continue;
      }
      else
      {

        smin=detrr*cspitch[i];
        saxis=detrr*sspitch[i];
        if(i==2)
        {
          ga=(1-(t3.nsg)*(exp(-((theta*RAD-t3.thetasg)*(theta*RAD-t3.thetasg))/(t3.wsg*t3.wsg))))*(1+t3.ncn*exp(-((theta*RAD-t3.thetacn)*(theta*RAD-t3.thetacn))/(t3.wcn*t3.wcn)))*pow(2/(exp(-smin/t3.warm[i])+exp(smin/t3.warm[i])), 2);
          if(rr>6000 && theta*RAD>t3.thetacn) ga=(1-(t3.nsg)*(exp(-((theta*RAD-t3.thetasg)*(theta*RAD-t3.thetasg))/(t3.wsg*t3.wsg))))*(1+t3.ncn)*pow(2/(exp(-smin/t3.warm[i])+exp(smin/t3.warm[i])), 2);
          
        }
        else ga=pow(2/(exp(-smin/t3.warm[i])+exp(smin/t3.warm[i])), 2);
        if(smin<sminmin)
        {
          sminmin=smin;
          which_arm=i;
        }
      }
      if(gd==0||fabs(zz)>(mc*HH))
      { m_3++;
      	*ne3=0;
      	return;
      }
      sech2=pow((2/(exp((rr-t3.B2s)/t3.Aa)+exp(-(rr-t3.B2s)/t3.Aa))), 2);
      ga=ga*sech2*gd;
      uu=pow(2/(exp((fabs(zz))/HH)+exp(-fabs(zz)/HH)), 2);
      ne3s+=t3.narm[i]*ga*uu;
    }
    *ne3=ne3s;
  }
}
