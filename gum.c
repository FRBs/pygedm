#include "cn.h"
void gum(double xx, double yy, double zz, double *ll, double *ne5, struct Gum t5)
{
  double slc, clc, sbc, cbc, xc, yc, zc, rgalc;
  double rp, RR, xyp, zp;
  double theta, alpha;
  double Dmin=1e5;
  double gn;
  //center of Gum Nebula
  const double lc=264;
  const double bc=-4;
  const double dc=450;
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
  alpha=-atan((-(t5.Agn)*(t5.Kgn)*xyp)/((t5.Agn)*sqrt((t5.Agn)*(t5.Agn)-xyp*xyp)));
  RR=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
  if(RR>500)Dmin=1e10;
  else
  {
    rp=sqrt((zp)*(zp)+(xyp)*(xyp));
    Dmin=fabs((RR-rp)*sin(theta+alpha));
  }

  if(Dmin>(mc*t5.Wgn))gn=0;
  else gn=1;

  *ne5=gn*(t5.ngn)*exp(-Dmin*Dmin/((t5.Wgn)*(t5.Wgn)));
}





