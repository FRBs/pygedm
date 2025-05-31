c 15  Jan 2003: fixed whicharm bug (c.f. TWJL email of 12 Dec '02)
c 24 June 2002: added calculations of path lengths through LISM components
c 26  May 2002: modified for NE2001 routines which are cleaned-up
c             versions of development routines.
c Nov 1999 - May 2002: development versions
c 1992-1993: TC93 version

      subroutine dmdsm(
     .   l,b,ndir,dmpsr,dist,limit,sm,smtau,smtheta,smiso)
      implicit none
      real l,b,dmpsr,dist,sm,smtau,smtheta,smiso
      integer ndir
      character*1 limit
      
Cf2py intent(in) l,b,ndir
Cf2py intent(out) limit,sm,smtau,smtheta,smiso
Cf2py intent(inout) dmpsr,dist

c  Computes pulsar distance and scattering measure
c  from model of Galactic electron distribution.

c  Input: real l	galactic longitude in radians
c         real b	galactic latitude in radians
c         integer ndir  >= 0 calculates dist from dmpsr  
c                       < 0 for dmpsr from dist
c Input or output:
c	  real dmpsr	(dispersion measure in pc/cm^3)
c         real dist	(distance in kpc)

c  Output:
c	  char*1 limit	(set to '>' if only a lower distance limit can be 
c			 given; otherwise set to ' ')
c         sm            (scattering measure, uniform weighting) (kpc/m^{20/3})
c         smtau         (scattering measure, weighting for pulse broadening)
c         smtheta       (scattering measure, weighting for angular broadening
c                        of galactic sources)
c	  smiso 	(scattering measure appropriate for calculating the
c			isoplanatic angle at the source's location
c       parameter(alpha = 11./3.)
c       parameter(pi = 3.14159)
c       parameter(c_sm = (alpha - 3.) / 2. * (2.*pi)**(4.-alpha) )


	real c_sm, c_u, sm_factor
        parameter(c_sm = 0.181)        
        parameter(c_u = 10.16)        
        parameter(sm_factor = c_sm * c_u)


	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
        common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN          

c parameters of large-scale components (inner+outer+arm components):
        real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
        common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c factors for controlling individual spiral arms:
c       narm:   multiplies electron density (in addition to the`fac'
c                     quantities)
c       warm:   arm width factors that multiply nominal arm width
c       harm:   arm scale height factors
c       farm:   factors that multiply n_e^2 when calculating SM

        integer narmsmax, narmsmax1
        parameter (narmsmax=5, narmsmax1=narmsmax+1)
        real narm, warm, harm, farm
        common/armfactors/
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)

	real armpaths, armdistances
	common/armpathlengths/ armpaths(5),armdistances(6)




	real dx0, dy0, dz0
	common/dxyz/dx0,dy0,dz0

	integer whicharm

c Large scale components:

	real ne1, ne2, nea
	real F1val, F2val, Faval

c Galactic center:

	real negc, Fgc

c LISM:
	real nelism, Flism 
        integer wlism, wLDR, wLHB, wLSB, wLOOPI
	real ldr_path, lhb_path, lsb_path, loopI_path
	real ldr_dist, lhb_dist, lsb_dist, loopI_dist
	integer wtemp

c clumps: 
        real necN, FcN
	integer hitclump

c voids:
	real nevN, FvN
	integer hitvoid, wvoid

c subroutines needed:
c	density_2001 (and those that it calls) in density.NE2001.f
c       scattering routines in scattering98.f
                                                                        
c other variables
	real x, y, z, r, rr
	real sl, cl, sb, cb

	real d, dstep, dtest, dstep_pc, dd

	real dm, dmstep
	real sm_sum1, sm_sum2, sm_sum3, sm_sum4, sm_term
        real sm_sum1_last, sm_sum2_last, sm_sum3_last, sm_sum4_last
	integer nstep
	integer ncount
	integer i

	real dm1, dm2, dma, dmgc, dmlism, dmcN, dmvN
	real sm1, sm2, sma, smgc, smlism, smcN, smvN
	real dsm1, dsm2, dsma, dsmgc, dsmlism, dsmcN, dsmvN

	integer wtotal
	real ne


	logical first
	real R0, rrmax, zmax, dmax
	data R0/8.5/
c	data rrmax/30.0/		
	data rrmax/50.0/		
c	data zmax/1.76/			
c	data zmax/5.00/			
	data zmax/25.00/		
        data dmax/50.0/     

	data first/.true./
	
	save
                                                                        
	if(first) then
c initial call to density routine to set variable values
c through read-in of parameter file:
        x = 0.0
        y = R0
        z = 0.0
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid)
c c       write(6,*) 'ne1,ne2,negc,nelism,necN,nevN = ',
c c    .              ne1,ne2,negc,nelism,necN,nevN
	first=.false.
        endif

	sl=sin(l)
	cl=cos(l)
	sb=sin(b)
	cb=cos(b)
	limit=' '
c	dstep=0.02			
c       dstep = min(h1, h2) / 10.       ! step size in terms of scale heights
	dstep=0.01
        if(ndir.lt.0) dtest=dist     
        if(ndir.ge.0) dtest=dmpsr/(n1h1/h1)   ! approximate test distance
        nstep = dtest / dstep	        ! approximate number of steps
        if(nstep.lt.10) dstep=dtest/10  ! make # steps >= 10

c  Sum until dm is reached (ndir >= 0) or dist is reached (ndir < 0). 
c  Guard against too few terms by counting number of terms (ncount) so that
c  routine will work for n_e models with large n_e near the Sun.

    5   continue
 	dstep_pc = 1000.*dstep
	dm=0.
        sm_sum1 = 0.           
        sm_sum2 = 0.           
        sm_sum3 = 0.           
	sm_sum4 = 0.			   

	do i=1,narmsmax1
	  armpaths(i) = 0.
	  armdistances(i) = 0.
	enddo

	dm1 = 0.
	dm2 = 0.
	dma = 0.
	dmgc = 0.
	dmlism = 0.
	dmcN = 0.
	dmvN = 0.

	sm1 = 0.
	sm2 = 0.
	sma = 0.
	smgc = 0.
	smlism = 0.
	smcN = 0.
	smvN = 0.

	ldr_path = 0.
	lhb_path = 0.
	lsb_path = 0.
	loopI_path = 0.

	ldr_dist = 0.
	lhb_dist = 0.
	lsb_dist = 0.
	loopI_dist = 0.

        ncount = 0

	d=-0.5*dstep
	do 10 i=1,99999
          ncount = ncount + 1
	  d=d+dstep			
	  r=d*cb
	  x=r*sl
	  y=R0-r*cl
	  z=d*sb
	  rr=sqrt(x**2 + y**2)		
	  if(ndir.ge.0.and.
     +      (d.gt.dmax.or.abs(z).gt.zmax.or.rr.gt.rrmax)) go to 20
                                                                       
          if(ndir.lt.3) then 
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid)
	  endif
        
	  if(ndir.ge.3) then 
            call density_2001(x+dx0,y+dy0,z+dz0,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid)
	  endif

c wlism = 1 causes the lism component to override smooth Galactic components
c wvoid = 1 overrides everything except clumps               
	  ne=
     .       (1.-wglism*wlism)*
     .       (wg1*ne1 +					 
     .        wg2*ne2 +
     .        wga*nea +
     .        wggc*negc) +
     .        wglism*wlism*nelism 
          ne = (1-wgvN*wvoid)*ne + wgvN*wvoid*nevN + wgcN*necN	
	  dmstep=dstep_pc*ne
	  dm=dm+dmstep			
	  wtotal = (1-wgvN*wvoid)*(1-wglism*wlism)
	  dm1 = dm1 + wtotal*wg1*ne1
	  dm2 = dm2 + wtotal*wg2*ne2
	  dma = dma + wtotal*wga*nea
	  dmgc = dmgc + wtotal*wggc*negc
	  dmlism = dmlism + (1.-wgvN*wvoid)*wglism*wlism*nelism
	  dmcN = dmcN + wgcN*necN
	  dmvN = dmvN + wgvN*wvoid*nevN

c         write(24,"('n:',7f10.6,1x))") 
c    .        ne1,ne2,nea,negc,nelism,necN,nevN
c        write(24,"(i2,1x,7(f10.5,1x))") 
c    .      wtotal,dm1,dm2,dma,dmgc,dmlism,dmcN,dmvN

c         sm_term = 
c    .       (1.-wglism*wlism)*
c    .       (wg1   * F1  * ne1**2 + 
c    .        wg2   * F2  * ne2**2 + 
c    .        wga   * Fa  * nea**2 + 
c    .        wggc  * Fgc * negc**2) +
c    .        wglism*wlism * Flism * nelism**2 
c	  sm_clumps = FcN * necN**2  
c	  sm_voids  = FvN * nevN**2
c         sm_term = (1-wgvN*wvoid) * sm_term 
c    .            + wgvN * wvoid * sm_voids
c    .            + wgcN * sm_clumps


	dsm1 = wtotal*wg1*ne1**2*F1
	dsm2 = wtotal*wg2*ne2**2*F2
	dsma =wtotal*wga*nea**2*Fa
	dsmgc = wtotal*wggc*negc**2*Fgc
	dsmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism**2*Flism
	dsmcN = wgcN*necN**2*FcN
	dsmvN = wgvN*wvoid*nevN**2*FvN

	sm_term = dsm1+dsm2+dsma+dsmgc+dsmlism+dsmcN+dsmvN

	sm1 = sm1 + dsm1 
	sm2 = sm2 + dsm2
	sma = sma + dsma 
	smgc = smgc + dsmgc
	smlism = smlism + dsmlism
	smcN = smcN + dsmcN 
	smvN = smvN + dsmvN 


        sm_sum1 = sm_sum1 + sm_term
        sm_sum2 = sm_sum2 + sm_term * d
        sm_sum3 = sm_sum3 + sm_term * d**2
        sm_sum4 = sm_sum4 + sm_term * d**1.66667

c pathlengths through LISM components:
c take into account the weighting hierarchy, LHB:LOOPI:LSB:LDR


	if(wlism .eq. 1) then
	  if(wlhb .eq. 1) then
	    lhb_path = lhb_path + dstep
	    lhb_dist = lhb_dist + d
	  endif
	  if(wloopI .eq. 1) then
    	    wtemp = (1-wLHB)
	    loopI_path = loopI_path + wtemp*dstep
	    loopI_dist = loopI_dist + wtemp*d
	  endif
	  if(wlsb .eq. 1) then
    	    wtemp = (1-wLHB)*(1-wloopI)
	    lsb_path = lsb_path + wtemp*dstep
	    lsb_dist = lsb_dist + wtemp*d
	  endif
	  if(wldr .eq. 1) then
    	    wtemp = (1-wLHB)*(1-wloopI)*(1-wlsb)
	    ldr_path = ldr_path + wtemp*dstep
	    ldr_dist = ldr_dist + wtemp*d
	  endif
	endif

c pathlengths: whicharm = 0,5 (currently).
c 	                  1,4 for the equivalent of the TC93 arms
c                         5   for the local arm
c                         0   means interarm paths

	armpaths(whicharm+1) = armpaths(whicharm+1) + dstep
	armdistances(whicharm+1) = armdistances(whicharm+1) + d

c       write(99,"(2(f8.3,1x), 7f10.6)") 
c    .     d, dm, sm_term,  sm_sum1, sm_sum2, sm_sum3, 
c    .     sm_sum1_last, sm_sum2_last, sm_sum3_last 
	if(ndir.ge.0.and.dm.ge.dmpsr) go to 30	
	if(ndir.lt.0.and.d.ge.dist) go to 40	
	sm_sum1_last = sm_sum1
	sm_sum2_last = sm_sum2
	sm_sum3_last = sm_sum3
	sm_sum4_last = sm_sum4
    
10	continue
c	stop 'loop limit'

20	limit='>'			
	dist=d-0.5*dstep
	go to 999

30	dist=d+0.5*dstep - dstep*(dm-dmpsr)/dmstep 
        
        if(ncount .lt. 10) then
	  dstep  = dstep / 10.
	  go to 5
	endif
	go to 999

40	dmpsr=dm-dmstep*(d+0.5*dstep-dist)/dstep
        if(ncount .lt. 10) then 
	   dstep = dstep / 10.
	   go to 5
	endif

999	continue

c normalize the mean distances:


	if(ldr_path .gt. 0.) then
          ldr_dist = ldr_dist / (ldr_path / dstep) 
	endif
	if(lhb_path .gt. 0.) then
          lhb_dist = lhb_dist / (lhb_path / dstep) 
	endif
	if(lsb_path .gt. 0.) then
          lsb_dist = lsb_dist / (lsb_path / dstep) 
	endif
	if(loopI_path .gt. 0.) then
          loopI_dist = loopI_dist / (loopI_path / dstep) 
	endif

        dd = d+0.5*dstep-dist
c subtract dd from armpath for latest arm (or iterarm) at end of LOS
c       armpaths(whicharm) = armpaths(whicharm)-dd    
	armpaths(whicharm+1) = armpaths(whicharm+1)-dd    

	do i=1,narmsmax1
	  armdistances(i) = 
     .        armdistances(i) / (max(1.0,armpaths(i)/dstep))	 
	enddo
	dm1 = dm1 * dstep_pc      
	dm2 = dm2 * dstep_pc      
	dma = dma * dstep_pc      
	dmgc = dmgc * dstep_pc       
	dmlism = dmlism * dstep_pc
	dmcN = dmcN * dstep_pc   
	dmvN = dmvN * dstep_pc        

c       dsm = sm_term * (d+0.5*dstep - dist)
c       dsm = sm_term * dd

c       sm_sum2 = sm_sum2 - dsm * d         
c       sm_sum3 = sm_sum3 - dsm * d**2 
c       sm_sum4 = sm_sum4 - dsm * d**1.67 

c       sm_sum1 = sm_sum1 - dsm         
c       write(99,*) 'dmdsm: sm_term, sm_sum1, sm_sum1_last = ',
c    .    sm_term, sm_sum1, sm_sum1_last

c	write(6,*) 'dmdsm: dsum1, sm_term = ',
c    .     sm_sum1-sm_sum1_last, sm_term
     
 	sm_sum1 = sm_sum1 - dd*(sm_sum1-sm_sum1_last)/dstep
	sm_sum2 = sm_sum2 - dd*(sm_sum2-sm_sum2_last)/dstep
	sm_sum3 = sm_sum3 - dd*(sm_sum3-sm_sum3_last)/dstep
	sm_sum4 = sm_sum4 - dd*(sm_sum4-sm_sum4_last)/dstep

c       sm_sum2 = sm_sum2 - dsm * dist                  
c       sm_sum3 = sm_sum3 - dsm * dist**2 
c       sm_sum4 = sm_sum4 - dsm * dist**1.67 


        sm = sm_factor * dstep * sm_sum1 
        smtau = 
     +     6. * sm_factor * dstep * (sm_sum2 / dist - sm_sum3 / dist**2)
        smtheta = 
     +     3. * sm_factor * dstep * (sm_sum1 + sm_sum3 / dist**2 - 
     +     2. * sm_sum2 / dist)
        smiso = sm_factor * dstep * sm_sum4

	sm1 = sm1 * sm_factor * dstep
	sm2 = sm2 * sm_factor * dstep
	sma = sma * sm_factor * dstep
	smgc = smgc * sm_factor * dstep
	smlism = smlism * sm_factor * dstep
	smcN = smcN * sm_factor * dstep
	smvN = smvN * sm_factor * dstep

        return 
	end

