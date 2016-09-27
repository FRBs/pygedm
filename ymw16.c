#include "cn.h"
#define max(a,b) ( ((a)>(b)) ? (a):(b) )
#define min(a,b) ( ((a)>(b)) ? (b):(a) )

/* Program version 1.1, 2016 September 1 */
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
double tsc(double dm){
  return 4.1e-11*pow(dm, 2.2)*(1+0.00194*dm*dm);
}
int main(int argc, char *argv[])
{
  	
  double ne0, ne, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10;
  double gl, gb, dordm, dist, xx, yy, zz, r, sl, cl, sb, cb, ll, hh;
  double nstep, dstep, dmstep;
  static double dd, dtest, dmpsr, rr;
  double dmm=0;
  double dm=0;
  double DM_MC=0;
  double DM_Gal=0;
  double DM_Host=0;
  double DDM;
  double tau_sc=0;
  double tau_Gal=0;
  double tau_MC=0;
  double tau_MC_sc=0;
  double R_g=0;
  double gd=0;
  
  //The localtion of Sun relative to GP and Warp
  double z_warp, zz_w, R_warp, theta_warp, theta_max;
  R_warp=8400;//pc
  theta_max=0.0; //In +x direction

  int WGN=0;
  int WLB=0; 
  int WLI=0;
  int WFB=0;
  int np, ndir,vbs, nk, uu, nn;
  char str[5];
  char dirname[64]="NULL",text[64]="";
  char *p;
  static int i, ncount;
  int w_lmc=0;
  int w_smc=0;
  int umc=1;
  char *s;
  struct Warp_Sun t0;
  struct Thick t1;
  struct Thin t2;
  struct Spiral t3;
  struct GC t4; 
  struct Gum t5; 
  struct LB t6;
  struct LI t7;
  struct FB t8;
  struct LMC t9;
  struct Dora t10;
  struct SMC t11;

  vbs=0;
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
      
  ymw16par(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10, &t11, dirname);

  if(ndir==1)printf("%s: gl=%8.3f gb=%8.3f DM=%8.2f", p, gl, gb, dordm);
  else printf("%s: gl=%8.3f gb=%8.3f Dist=%9.1f", p, gl, gb, dordm);
 
  
  ll=gl;
  gl=gl/RAD;
  gb=gb/RAD;
  sl=sin(gl);
  sb=sin(gb);
  cl=cos(gl);
  cb=cos(gb);    
  dstep=5.0;
  
  if(np==-1){                 // FRBs
    if(ndir==1)uu=0;//dm---dist
    else uu=1;//dist---dm
    ndir=2;
    DDM=dordm;
    dordm=100000;
    nk=20000;
  }
  
  if(np==0){                 // Magellanic Cloud
    nk=20000;
  }
 
  if(np==1){                  //Galactic pulsars
    nk=5000;
  }
  if(ndir==1){
    dm=dordm;
    if(np==1)tau_sc=tsc(dm);
    dtest=dm/N0;
    nstep=dtest/dstep;
    if(nstep<200) dstep=dtest/200;
    if(vbs>=1){
      printf("\ndtest=%lf, nstep=%lf, dstep=%lf\n", dtest, nstep, dstep);
    }
  }
  if(ndir==2){
    dist=dordm;
    dtest=dist;
    nstep=dtest/dstep;
    if(nstep<200) dstep=dtest/200;
    if(vbs>=1){
      printf("\ndtest=%lf, nstep=%lf, dstep=%lf\n", dtest, nstep, dstep);
    }
  } 
  
  
  dd=-0.5*dstep;
  ncount=0;

  for(i=1;i<=nk;i++){
    ncount++;
    if(vbs>=2){
      printf("ncount=%d, dstep=%lf\n", ncount,dstep);
    }
    dd+=dstep;
    r=dd*cb;     /* r is different from rr */
    xx=r*sl;
    yy=R0*1000-r*cl;
    zz=dd*sb+t0.z_Sun;
    rr=sqrt(xx*xx+yy*yy);
    
    /* Definition of warp */
    if(rr<R_warp){
      zz_w=zz; 
    }else{
      theta_warp=atan2(yy,xx);
      z_warp=t0.Gamma_w*(rr-R_warp)*cos(theta_warp-theta_max);
      zz_w=zz-z_warp;
    }

    if(vbs>=2)
    {
      printf("dd=%lf, xx=%lf, yy=%lf, zz=%lf, rr=%lf\n", dd, xx, yy, zz, rr);
      printf("theta_warp=%lf, z_warp=%lf, zz_w=%lf\n",theta_warp,z_warp,zz_w);
    }
    R_g=sqrt(xx*xx+yy*yy+zz_w*zz_w);

    /* DM to Distance */
    
    if(ndir==1){   	
      if(dmm<=dm){
        if(R_g<=35000){
	      thick(xx, yy, zz_w, &gd, &ne1, rr, t1);
          thin(xx, yy, zz_w, gd, &ne2, rr, t2);
          spiral(xx, yy, zz_w, gd, &ne3, rr, t3, dirname);
          galcen(xx, yy, zz, &ne4, t4);
          gum(xx, yy, zz, &ll, &ne5, t5);
          localbubble(xx, yy, zz, &ll, &ne6, &hh, t6);
          nps(xx, yy, zz, &ne7, &WLI, t7);
          fermibubble(xx, yy, zz, &WFB);
        }else{
          if(np==1){
	    dstep=5;
          }else{
            dstep=200;
            if(w_lmc>=1||w_smc>=1) dstep=5;
            lmc(gl,gb,dd,&w_lmc,&ne8,t9);
            dora(gl,gb,dd,&ne9,t10);
            smc(xx, yy, zz,&w_smc, &ne10, t11);	
          }
	} 
	if(WFB==1){
	  ne1=t8.J_FB*ne1;
	}
	ne0=ne1+max(ne2,ne3);

	if(hh>110){       /* Outside LB */
	  if(ne6>ne0 && ne6>ne5){
	    WLB=1;
	  }else{
	    WLB=0;
	  }
	}else{            /* Inside LB */
	  if(ne6>ne0){
	    WLB=1;
	  }else{
	    ne1=t6.J_LB*ne1;
	    ne0=ne1+max(ne2,ne3);
	    WLB=0;
	  }
	}
	if(ne7>ne0){     /* Loop I */
	  WLI=1;
	}else{
	  WLI=0;
	}        
	if(ne5>ne0){     /* Gum Nebula */
	  WGN=1;
	}else{
	  WGN=0;
	}

	/* Galactic ne */
	ne=(1-WLB)*((1-WGN)*((1-WLI)*(ne0+ne4+ne8+ne9+ne10)+WLI*ne7)+WGN*ne5)+WLB*ne6;

	if(vbs>=2){
	  printf("ne=%lf, ne1=%lf, ne2=%lf, ne3=%lf, ne4=%lf, ne5=%lf, ne6=%lf, ne7=%lf, ne8=%lf, ne9=%lf, ne10=%lf\n", ne, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10);
	}
	dmstep=ne*dstep;
	if(dmstep<=0.000001)dmstep=0;
	if(vbs>=2){
	  printf("dmstep=%lf, dstep=%lf\n", dmstep, dstep);
	}
	dmm+=dmstep;
	dist=dd;
	if(np==0&&umc==1){
	  if(R_g>35000){
	    DM_Gal=dmm;
	    tau_Gal=0.5*tsc(dmm);
	    printf(" DM_Gal:%8.2f",DM_Gal);
	    umc++;
	  } 
	}
	if(i==nk){ 
	  dist+=0.5*dstep;
	  if(dist>100000)dist=100000;
	  if(np==0){
	    DM_MC=dmm-DM_Gal;
	    tau_MC=0.5*tsc(DM_MC);
	    tau_MC_sc=max(tau_Gal, tau_MC);
	    printf(" DM_MC:%8.2f",DM_MC);
	  }   
	  if(np==0)printf(" Dist:%9.1f log(tau_sc):%7.3f %s\n",dist, log10(tau_MC_sc),text);
	  if(np==1)printf(" DM_Gal:%8.2f Dist:%9.1f log(tau_sc):%7.3f %s\n", dmm, dist,log10(tau_sc),text);
	}	    
        if(vbs>=2){
	  printf("dmm=%lf\n", dmm);
	}
      }
      else{
	dist=dd-0.5*dstep-(dstep*(dmm-dm))/dmstep;
	if(np==0){
	  DM_MC=dm-DM_Gal;
          if(DM_Gal==0){
            DM_MC=0; 
            DM_Gal=dm;
            tau_MC_sc=tsc(dm);
            printf(" DM_Gal:%8.2f ", DM_Gal);
          }
          else{
	    tau_MC=0.5*tsc(DM_MC);
            tau_MC_sc=max(tau_MC, tau_Gal);
          }  
          printf(" DM_MC:%8.2f", DM_MC);
        }
	if(np==0)printf(" Dist:%9.1f log(tau_sc):%7.3f %s\n",dist,log10(tau_MC_sc),text);
	if(np==1)printf(" DM_Gal:%8.2f Dist:%9.1f log(tau_sc):%7.3f %s\n", dm, dist,log10(tau_sc),text);
	break;
      }
    }

    /* Distance to DM */
    
    if(ndir==2){
      if(dd<=dtest){
        if(R_g<=35000){
	      thick(xx, yy, zz_w, &gd, &ne1, rr, t1);
          thin(xx, yy, zz_w, gd, &ne2, rr, t2);
          spiral(xx, yy, zz_w, gd, &ne3, rr, t3, dirname);
          galcen(xx, yy, zz, &ne4, t4);
          gum(xx, yy, zz, &ll, &ne5, t5);
          localbubble(xx, yy, zz, &ll, &ne6, &hh, t6);
          nps(xx, yy, zz, &ne7, &WLI, t7);
          fermibubble(xx, yy, zz, &WFB);
        }else{
	  if(np==1)dstep=5;
	  else{
	    dstep=200;
	    if(np==-1)dstep=5;
	    if(w_lmc>=1||w_smc>=1) dstep=5;
	    lmc(gl,gb,dd,&w_lmc,&ne8,t9);
	    dora(gl,gb,dd,&ne9,t10);
	    smc(xx, yy, zz,&w_smc, &ne10, t11);
	  } 
	}       
	if(WFB==1){
	  ne1=t8.J_FB*ne1;
	}
	ne0=ne1+max(ne2,ne3);

        if(hh>110){       /* Outside LB */
          if(ne6>ne0 && ne6>ne5){
	    WLB=1;
          }else{
	    WLB=0;
	  }
	}else{            /* Inside LB */
	  if(ne6>ne0){
	    WLB=1;
	  }else{
	    ne1=t6.J_LB*ne1;
	    ne0=ne1+max(ne2,ne3);
	    WLB=0;
	  }
        }
        if(ne7>ne0){     /* Loop I */
	  WLI=1;
        }else{
          WLI=0;
        }        
        if(ne5>ne0){     /* Gum Nebula */
          WGN=1;
        }else{
          WGN=0;
        }        
	/*  Galactic ne */
        ne=(1-WLB)*((1-WGN)*((1-WLI)*(ne0+ne4+ne8+ne9+ne10)+WLI*ne7)+WGN*ne5)+WLB*ne6;
	
	if(vbs>=2){
          printf("ne=%lf, ne1=%lf, ne2=%lf, ne3=%lf, ne4=%lf, ne5=%lf, ne6=%lf, ne7=%lf, ne8=%lf, ne9=%lf, ne10=%lf\n", ne, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10);
        }
	dmstep=ne*dstep;
	if(dmstep<=0.000001)dmstep=0;
	dm+=dmstep;
	if(np!=1&&umc==1){
	  if(R_g>35000){ 
	    DM_Gal=dm;
	    tau_Gal=0.5*tsc(dm);
	    printf(" DM_Gal:%8.2f",dm);
	    umc++;
	  } 
        }  
        if(i==nk&&np!=-1){
          dmpsr=dm;
          if(np==0){
            DM_MC=dm-DM_Gal;
            tau_MC=0.5*tsc(DM_MC);
            printf(" DM_MC:%8.2f", DM_MC);
          }
          tau_sc=tsc(dmpsr);
          tau_MC_sc=max(tau_Gal, tau_MC);
          if(np==0)printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr,log10(tau_MC_sc),text);
          if(np==1)printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr, log10(tau_sc),text);
        }
	
	if(i==nk&&np==-1){
          if(dordm==100000){
	    DM_MC=dm-DM_Gal;
	    printf(" DM_MC:%8.2f",DM_MC);
	  } 
          frb_d(DDM, DM_Gal, DM_MC, DM_Host, uu, vbs, text);
          break;
        }
      }
      else{
	dmpsr=dm+(dmstep*(dtest-(dd-0.5*dstep)))/dstep;
	if(np==0){ 	
	  DM_MC=dmpsr-DM_Gal;
	  if(DM_Gal==0){
	    DM_MC=0;
	    DM_Gal=dmpsr;
	    tau_MC_sc=tsc(dmpsr);
	    printf(" DM_Gal:%8.2f", DM_Gal);
	  } 
	  else{
	    tau_MC=0.5*tsc(DM_MC);
	    tau_MC_sc=max(tau_Gal, tau_MC);
	  } 
	  printf(" DM_MC:%8.2f", DM_MC);
	}
	tau_sc=tsc(dmpsr);
	if(np==0)printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr,log10(tau_MC_sc),text);
	if(np==1)printf(" DM:%8.2f log(tau_sc):%7.3f %s\n", dmpsr, log10(tau_sc),text);
	break;
      }
    }    
  }
}
