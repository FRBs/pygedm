#include "cn.h"
void thin(double xx, double yy, double zz, double gd, double *ne2, double rr, struct Thin t2)
{
  double g2, Hg, ex1, ex2, g3, HH;
  Hg=32+0.0016*rr+0.0000004*pow(rr, 2);
  HH=t2.K2*Hg;
  if((rr-t2.B2)>(mc*t2.A2)||(fabs(zz)>(mc*HH))) 
  {
  	*ne2=0;
  	return;	
  }
  else
  {
    g3=(rr-t2.B2)/t2.A2;
    g2=pow(2/(exp(-g3)+exp(g3)), 2);
  }
  *ne2=t2.n2*gd*g2*pow(2/(exp(-zz/HH)+exp(zz/HH)), 2);
}
