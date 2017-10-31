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
/* Program version 1.2.3, 2017 March 22 */
void usage(int status)
{
  printf("\nConverts from DM to Dist and vice versa - electron density model\n");
  printf("Usage:\n");
  printf("ymw16 [-h] [-t text] [-d <dirname>] [-v] [-V] <mode> gl gb DM/Dist [DM_Host] ndir\n");
  printf("-h prints this help page\n");
  printf("-t <text>, where <text> (no spaces, max 64 char) is appended to the output line\n");
  printf("-d <dirname>, where <dirname> is a directory containing YMW16 data files\n");
  printf("-v prints diagnostics\n");
  printf("-V prints more diagnostics\n");
  printf("<mode> is Gal, MC, or IGM\n");
  printf("gl, gb in deg, DM in cm^-3 pc, Dist in pc (Gal, MC) or Mpc (IGM)\n");
  printf("DM_Host is the contribution of the FRB host galaxy to the observed DM\n");
  printf("Only used for IGM mode; optional input, default=100 cm^-3 pc\n");
  printf("ndir=1 for DM->Dist, ndir=2 for Dist->DM\n");
  printf("Output includes log(tau_sc) where tau_sc is scattering time in sec at 1 GHz\n");
  printf("\n");
  exit(status);
}
char *strupr(char *str)
{
  char *p = str;
  while (*p != '\0'){
    if(*p >= 'a' && *p <= 'z')*p -= 0x20;
    p++;
   }
  return str;
}

int main(int argc, char *argv[])
{
  double gl, gb, dordm;
  double DM_Host=0;
  int ndir, np, ns;
  int vbs=0;
  char dirname[64]="NULL",text[64]="";

  char str[5];
  char *p;
  char *s;
  argc--; argv++;
  
  if(argc < 5)usage(1);
  while(argc > 5){                /* Get command line inputs */
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
      case 'V':
	vbs=2;
	break;
      default:
	usage(1);
      }
    }
    else{
      if(argc>6){	
	printf("Extra parameters exist in input\n");	
	usage(1);
      }
      else break;	 	
    }   
  }
  if(argc==5){	
    if(sscanf(*argv,"%s",str) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
  
    argc--; argv++;
    if(sscanf(*argv,"%lf",&gl) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&gb) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&dordm) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%d",&ndir) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    DM_Host=100;//default
  }
  if(argc==6){
    ns=argc;
    if(sscanf(*argv,"%s",str) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    
    argc--; argv++;
    if(sscanf(*argv,"%lf",&gl) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&gb) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&dordm) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%lf",&DM_Host) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
    argc--; argv++;
    if(sscanf(*argv,"%d",&ndir) != 1){
      printf("Incorrect arguments\n");
      usage(1);
    }
  }	
  //convert to upper case
  
  p=strupr(str);
  
  if(strcmp(p,"IGM") == 0) np=-1;    // IGM
  else if(strcmp(p,"MC") == 0) np=0; // Mag Clouds
  else if(strcmp(p,"GAL") == 0){     // Galaxy
    np=1; 
    p="Gal";
  }  
  else{
    printf("please input correct model\n");
    usage(1);
    exit(1);
  }
  if(ns==6&&np!=-1){
  printf("Extra parameters exist in input\n");	
  usage(1);	
  }

  if(ndir!=1&&ndir!=2){
    printf("please input correct ndir\n");
    usage(1);
  }

  if(!strcmp(dirname,"NULL")){
    if(getenv("YMW16_DIR")==NULL){
      printf("Warning: YMW16_DIR set to local directory\n");
      strcpy(dirname,"./");
    }else{
      strcpy(dirname,getenv("YMW16_DIR"));
    }
  }
  if(vbs>=1)printf("File directory: %s\n",dirname);
      

  if(ndir==1)printf("%s: gl=%8.3f gb=%8.3f DM=%8.2f", p, gl, gb, dordm);
  else printf("%s: gl=%8.3f gb=%8.3f Dist=%9.1f", p, gl, gb, dordm);
  dmdtau(gl, gb, dordm, DM_Host, ndir, np, vbs, dirname, text); 
}
