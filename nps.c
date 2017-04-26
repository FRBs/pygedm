#include "cn.h"
void nps(double xx,double yy,double zz,double *ne7, int *WLI, struct LI t7)
{
  double x_c, y_c, z_c;
  double gLI;
  double theta_LI;
  
  if(m_7>=1)return;
  
  theta_LI=(t7.thetaLI)/RAD;
  x_c=-10.156;
  y_c=8106.206;
  z_c=10.467;
  double rr,theta;
  rr=sqrt((xx-x_c)*(xx-x_c)+(yy-y_c)*(yy-y_c)+(zz-z_c)*(zz-z_c));
  theta=acos(((xx-x_c)*(cos(theta_LI))+(zz-z_c)*(sin(theta_LI)))/rr)*RAD;
  *WLI=1;
  if(fabs(rr-t7.RLI)>(mc*t7.WLI)||fabs(theta)>(mc*t7.detthetaLI)) 
  { if(rr>500)m_7++;
  	*ne7=0;
  	return;
  }
  else gLI=1;
  *ne7=gLI*(t7.nLI)*exp(-((rr-(t7.RLI))*(rr-(t7.RLI)))/((t7.WLI)*(t7.WLI)))*exp((-theta*theta)/((t7.detthetaLI)*(t7.detthetaLI)));
}

