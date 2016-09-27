#include "cn.h"
void fermibubble(double xx, double yy, double zz, int *wfb)
{
  double fbnz, fbsz, na, nb, sa, sb, N, S;
//center of Fermi bubble
  fbnz=0.5*8300*tan(50/RAD);
  fbsz=-(0.5*(8300*tan(50/RAD)-8300*tan(0/RAD))+(8300*tan(0/RAD)));
//min_axis and max_axis of Fermi bubble
  na=fbnz;
  nb=8300*tan(20/RAD);
  sa=0.5*(8300*tan(50/RAD)-8300*tan(0/RAD));
  sb=8300*tan(20/RAD);
  N=(pow(xx, 2)/(nb*nb))+(pow(yy, 2)/(nb*nb))+(pow(zz-fbnz, 2)/(na*na));
  S=(pow(xx, 2)/(sb*sb))+(pow(yy, 2)/(sb*sb))+(pow(zz-fbsz, 2)/(sa*sa));
  if(N<1||S<1)*wfb=1;
  else *wfb=0;
}
