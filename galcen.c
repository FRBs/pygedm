#include "cn.h"
void galcen(double xx, double yy, double zz, double *ne4 ,struct GC t4)
{
  double Xgc=50;
  double Ygc=0;
  double Zgc=-7;
  double Rgc, Az, Ar;
  double gc;
  Rgc=sqrt((xx-Xgc)*(xx-Xgc)+(yy-Ygc)*(yy-Ygc));
  if(Rgc>(mc*t4.Agc)||fabs(zz)>(mc*t4.Hgc)) 
  {
  	*ne4=0;
  	return;
  }
  else
  {
    gc=1;
    Ar=exp(-(Rgc*Rgc)/(t4.Agc*t4.Agc));
    Az=pow(2/(exp((zz-Zgc)/t4.Hgc)+exp(-(zz-Zgc)/t4.Hgc)),2);
  }
  *ne4=t4.ngc*Ar*Az*gc;
}


