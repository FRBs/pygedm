#include "cn.h"

void thick(double xx, double yy, double zz, double *gd, double *ne1, double rr, struct Thick t1){
  
  double gdd,gg;

  if(fabs(zz)> mc*t1.H1 || (rr-t1.Bd)> mc*t1.Ad){
    *ne1=0;
    return;
  }else{
    if(rr<t1.Bd){
      gdd=1;
    }else{ 
      gg=exp(-(rr-t1.Bd)/t1.Ad)+exp((rr-t1.Bd)/t1.Ad);
      gdd=pow(2/gg,2);
    }
  }
  *ne1=t1.n1*gdd*pow(2/(exp(-fabs(zz)/t1.H1)+exp(fabs(zz)/t1.H1)), 2);
  *gd=gdd;
}
