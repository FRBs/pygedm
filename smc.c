#include"cn.h"
void smc(double xx, double yy, double zz, double *ne10, struct SMC t11)
{
  double rad=57.295779;
  double R=8.3;//kpc
  double ls=303.728914;
  double bs=-44.299212;
  double ds=59700;//pc
  double Asmc=3000;//pc
  int gsmc;
  double xc, yc, zc, sls, cls, sbs, cbs, rgals, Rsmc;
  sls=sin(ls/rad);
  cls=cos(ls/rad);
  sbs=sin(bs/rad);
  cbs=cos(bs/rad);
  rgals=ds*cbs;
  xc=rgals*sls;
  yc=R*1000-rgals*cls;
  zc=ds*sbs;
//printf("xc=%lf,yc=%lf,zc=%lf\n", xc, yc, zc);
  Rsmc=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
  if(Rsmc>(mc*Asmc))gsmc=0;
  else gsmc=1;
  *ne10=(t11.nsmc)*gsmc*exp(-(Rsmc*Rsmc)/(Asmc*Asmc));
//printf("nesmc=%lf\n",t8.nesmc);
}

