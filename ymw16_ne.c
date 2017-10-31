/*Copyright (C) 2016, 2017  J. M. Yao, R. N. Manchester, N. Wang.

This file is part of the YMW16 program. YMW16 is a model for the
distribution of free electrons in the Galaxy, the Magellanic Clouds
and the inter-galactic medium that can be used to estimate distances
for real or simulated pulsars and fast radio bursts (FRBs) based on
their position and dispersion measure.

YMW16 is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

YMW16 is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License,
available at http://www.gnu.org/licenses/, for more details. 

Please report any issues or bugs at
https://bitbucket.org/psrsoft/ymw16/issues/new/ or directly to the
authors. Please provide an example illustrating the problem.

Jumei Yao (yaojumei@xao.ac.cn), Richard N Manchester
(dick.manchester@csiro.au), Na Wang (na.wang@xao.ac.cn).
*/
#include "cn.h"
/* Program version 1.2.4, 2017 August*/
void usage(int status)
{
  printf("\nCalculate electron density at a given point with galactocentric coordinates (x, y, z)/with (gl, gb, dist)\n");
  printf("Usage:\n");
  printf("ymw16_ne [-h] [-t text] [-d <dirname>] [-v] crd1 crd2 crd3 ncrd\n");
  printf("-h prints this help page\n");
  printf("-t <text>, where <text> (no spaces, max 64 char) is appended to the output line\n");
  printf("-d <dirname>, where <dirname> is a directory containing YMW16 data files\n");
  printf("-v prints diagnostics\n");
  printf("For ncrd=1, input Galactocentric x, y and z in pc\n");
  printf("For ncrd=2, input gl, gb in deg, Dist in pc\n");
  printf("Output ne in (pc cm^-3)\n");
  printf("\n");
  exit(status);
}


int main(int argc, char *argv[])
{
  double p1, p2, p3, ne;
  double x, y, z;
  double gl, gb, dist;
  int ncrd=0; 
  int np;
  int vbs=0;
  char dirname[64]="NULL",text[64]="";

  char *s;
  argc--; argv++;
  if(argc < 4)usage(1);
  while(argc > 4){                /* Get command line inputs */
    if((*argv)[0] == '-'){
      s=argv[0]+1;
      argc--; argv++;

      switch(*s){
      case 'h':
      case '?':
        usage(0);
      case 't':
        if(sscanf(*argv,"%s",text) != 1)usage(1);
        argc--; argv++;
        break;
      case 'd':
        if(sscanf(*argv,"%s",dirname) != 1)usage(1);
        argc--; argv++;
        break;
      case 'v':
        vbs=1;
        break;
      default:
        usage(1);
      }
    }
    else{
      if(argc>=5){
        printf("Extra parameters exist in input\n");
        usage(1);
      }
      else break;
    }
  }
  if(argc==4){
    if(sscanf(*argv,"%lf",&p1) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&p2) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&p3) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%d",&ncrd) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
   }
   
if(!strcmp(dirname,"NULL")){
    if(getenv("YMW16_DIR")==NULL){
      printf("Warning: YMW16_DIR set to local directory\n");
      strcpy(dirname,"./");
    }else{
      strcpy(dirname,getenv("YMW16_DIR"));
    }
  }
 if(ncrd==1){
 x=p1;
 y=p2;
 z=p3;
 }
 else{
  if(ncrd==2){
  gl=p1;
  gb=p2;
  dist=p3;
  }
  else printf("Incorrect ncrd\n");
 }
  ne = ne_crd(&x, &y, &z, &gl, &gb, &dist, ncrd, vbs, dirname, text);
  
  if(ncrd==1){
    printf("x= %8.1f y= %8.1f z= %8.1f gl: %8.3f gb: %8.3f Dist: %8.1f ne: %9.6f %s\n", x, y, z, gl, gb, dist, ne, text);
  }else{
    printf("gl= %8.3f gb= %8.3f Dist= %8.1f x: %8.1f y: %8.1f z: %8.1f ne: %9.6f %s\n", gl, gb, dist, x, y, z, ne, text); 
  }
}
