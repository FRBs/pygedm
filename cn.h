#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#define R0 8.3
#define RAD 57.295779
#define N0 0.013
#define PI 3.14159265
#define mc 6
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )
#define MIN(a,b) ( ((a)>(b)) ? (b):(a) )
int m_3, ww,m_5, m_6, m_7;
struct Warp_Sun
{
  double Gamma_w;
  double z_Sun;
};

struct Thick
{
  double Ad;
  double Bd; 
  double n1; 
  double H1;
};

struct Thin
{
  double A2;
  double B2; 
  double n2;
  double K2;
};

struct Spiral
{
  double Ads;
  double Bds;
  double B2s;
  double Ka;
  double narm[5];
  double warm[5];
  double Aa;
  double ncn;
  double wcn;
  double thetacn;
  double nsg;
  double wsg;
  double thetasg;
  
};

struct GC
{
  double ngc;
  double Agc;
  double Hgc;
};

struct Gum
{
  double Kgn;
  double ngn;
  double Wgn;
  double Agn;
};

struct LB
{
  double J_LB;
  double nlb1;
  double detlb1;
  double wlb1;
  double hlb1;
  double thetalb1;
  double nlb2;
  double detlb2;
  double wlb2;
  double hlb2;
  double thetalb2;
};

struct LI
{
  double nLI;
  double RLI;
  double WLI;
  double detthetaLI;
  double thetaLI;
};

struct FB
{
  double J_FB;
};

struct LMC
{
  double nlmc;
};

struct Dora
{
  double n30D;
};

struct SMC
{
  double nsmc;
};
int ymw16par(struct Warp_Sun *t0, struct Thick *t1, struct Thin *t2,  struct Spiral *t3, struct GC *t4, struct Gum *t5, struct LB *t6, struct LI *t7, struct FB *t8, struct LMC *t9, struct Dora *t10, struct SMC *t11, char *dirname);
void dmdtau(double gl, double gb ,double dordm, double DM_Host, int ndir, int np, int vbs, char *dirname, char *text);
void thick(double xx, double yy, double zz, double *gd, double *ne1, double rr, struct Thick t1);
void thin(double xx, double yy, double zz, double gd, double *ne2, double rr, struct Thin t2);
void spiral(double xx, double yy, double zz, double gd, double *ne3, double rr, struct Spiral t3, char *dirname);
void galcen(double xx, double yy, double zz, double *ne4 ,struct GC t4);
void gum(double xx, double yy, double zz, double *ll, double *ne5, struct Gum t5);
void localbubble(double xx, double yy, double zz, double *ll, double *bb,double *ne6, double *WW, struct LB t6);
void nps(double xx,double yy,double zz,double *ne7, int *WLI, struct LI t7);
void fermibubble(double xx, double yy, double zz, int *wfb);
void lmc(double l, double b, double d, int *w_lmc, double *ne8, struct LMC t9);
void dora(double l, double b, double d, double *ne9, struct Dora t10);
void smc(double xx, double yy, double zz, int *w_smc, double *ne10, struct SMC t11);
void frb_d(double DDM, double DM_Gal, double DM_MC, double DM_Host, int uu, int vbs, char *text);
double tsc(double dm);
