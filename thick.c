#include "cn.h"
void thick(double xx, double yy, double zz, double *ne1, double rr, struct Thick t1)
{
  double g1, g2, gd;
  if(fabs(zz)>(mc*t1.H1)) g1=0;
  else g1=1;
  if(rr<t1.Bd)gd=1;
  else
  {
    if((rr-t1.Bd)>mc*t1.Ad)gd=0;
    else
    { 
      g2=exp(-(rr-t1.Bd)/t1.Ad)+exp((rr-t1.Bd)/t1.Ad);
      gd=pow(2/g2, 2);
    }
  }
  *ne1=t1.n1*g1*gd*pow(2/(exp(-fabs(zz)/t1.H1)+exp(fabs(zz)/t1.H1)), 2);
}
