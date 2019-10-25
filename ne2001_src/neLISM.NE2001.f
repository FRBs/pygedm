c routines to calculate the electron density for the 
c Local Interstellar Medium
c
c JMC 26 August-11 Sep. 2000
c     25 October 2001: modified to change weighting scheme
c                      so that the ranking is LHB: LSB: LDR
c                      (LHB overrides LSB and LDR; LSB overrides LDR)
c     16 November 2001: added Loop I component with weighting scheme
c		        LHB:LOOPI:LSB:LDR	
c		        LHB   overides everything,
c			LOOPI overrides LSB and LDR
c			LSB   overrides LDR
c			LISM  overrides general Galaxy
c     20 November 2001: The LOOPI component is truncated below z=0 
c
c after discussions with Shami Chatterjee
c the sizes, locations and densities of the LISM components
c are based largely on work by Toscano et al. 1999
c and Bhat et al. 1999

	real function ne_LISM(x,y,z,FLISM,wLISM,
     .                       wldr, wlhb,wlsb,wloopI)		
	implicit none
	real x,y,z,FLISM
        integer wLISM

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

        logical first
	data first /.true./
	
c functions:

	real neLDRQ1
	real neLSB
	real neLHB2
	real neLOOPI

c other variables:

	real nelsbxyz, nelhbxyz, neldrq1xyz, neloopIxyz 
	real FLDRQ1r, FLSBr, FLHBr, FLOOPIr		! 'r' for returned value
	integer wLDR, wLSB, wLHB, wLOOPI

	if(first) then					! read parameters for LISM
	  open(11,file='nelism.inp',status='unknown')
	  read(11,*)
	  read(11,*) aldr,bldr,cldr
	  read(11,*) xldr,yldr,zldr
	  read(11,*) thetaldr,neldr0,Fldr

	  read(11,*) alsb,blsb,clsb
	  read(11,*) xlsb,ylsb,zlsb
	  read(11,*) thetalsb,nelsb0,Flsb

	  read(11,*) alhb,blhb,clhb
	  read(11,*) xlhb,ylhb,zlhb
	  read(11,*) thetalhb,nelhb0,Flhb

	  read(11,*) xlpI,ylpI,zlpI
	  read(11,*) rlpI,drlpI
	  read(11,*) nelpI,dnelpI,FlpI,dFlpI

c     write(99,*) xlpI,ylpI,zlpI
c 	  write(99,*) rlpI,drlpI
c 	  write(99,*) nelpI,dnelpI,FlpI,dFlpI

	  first=.false.
	endif

	neldrq1xyz = neLDRQ1(x,y,z,FLDRQ1r,wLDR)	! low density region in Q1
	nelsbxyz   = neLSB(x,y,z,FLSBr,wLSB)		! Local Super Bubble
	nelhbxyz = neLHB2(x,y,z,FLHBr,wLHB)		! Local Hot Bubble
	neloopIxyz = neLOOPI(x,y,z,FLOOPIr,wLOOPI)	! Loop I


c weight the terms so that the LHB term overrides the other 
c terms (we want the density to be low in the LHB, lower than
c in the other terms.

 	ne_LISM =   (1-wLHB)   * 
     .           ( 
     .             (1-wLOOPI) * (wLSB*nelsbxyz + (1-wLSB)*neldrq1xyz)
     .         +     wLOOPI * neloopIxyz
     .           )   
     .         +     wLHB  * nelhbxyz

 	FLISM = (1-wLHB) * 
     .          (
     .             (1-wLOOPI) * (wLSB*FLSBr + (1-wLSB)*FLDRQ1r) 
     .         +     wLOOPI * FLOOPIr
     .          )
     .         +     wLHB  * FLHBr

c return the maximum weight of any of the terms for
c combining with additional terms external to this routine.

	wLISM = max(wLOOPI, max(wLDR, max(wLSB, wLHB)))

c temporary next 3 lines:
c	ne_LISM = nelhbxyz
c	flism = flhb
c	wlism = wlhb

c	write(97,"(9(f8.3,1x))") wLDR, wLSB, wLHB,
c    .      nelsbxyz, neldrq1xyz,nelhbxyz,
c    .      flsb, fldrq1, flhb


	return
	end

	real function  neLDRQ1(x,y,z,FLDRQ1r,wLDRQ1) 	! Low Density Region in Q1
	implicit none
	real x,y,z,FLDRQ1r
	integer wLDRQ1

c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLDRQ1 = electron density in local hot bubble that
c	        is modeled as an ellipsoidal trough. 
c	FLDRQ1 = fluctuation parameter
c	wLDRQ1  = weight of LDRQ1 component used to combine
c		with other components of electron density.
c		wLDRQ1 =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLDRQ1)*ne_other + neLDRQ1
	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xldr, yldr, zldr		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis 
					!    (x axis points toward l=90)
	
	real q

	real ap, bp, cp, dp
	real s, c
	
        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./
	save

c	data aa, bb, cc       	/ 0.6, 0.40, 0.3 / 	! GUESS     
c	data xldr, yldr, zldr 	/ 0.6, 7.86, 0. / 	! GUESS
c	data theta		/ -45. / 		! GUESS

c	data netrough  /0.010/		! GUESS
c	data Ftrough   / 2./		! GUESS

	aa=aldr
	bb=bldr
	cc=cldr
	theta=thetaldr
	netrough =neldr0
	Ftrough=Fldr
	
	if(first) then
	  s = sin(theta/radian) 
	  c = cos(theta/radian) 
	  ap = (c/aa)**2 + (s/bb)**2 
	  bp = (s/aa)**2 + (c/bb)**2 
	  cp = 1./cc**2
	  dp =  2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLDRQ1 = 0.
	wLDRQ1 = 0
	FLDRQ1r = 0.
	q = (x-xldr)**2*ap
     .    + (y-yldr)**2*bp 
     .    + (z-zldr)**2*cp
     .    + (x-xldr)*(y-yldr)*dp
	  if(q .le. 1.0) then	! inside 
	    neLDRQ1 = netrough 
	    FLDRQ1r = Ftrough 
	    wLDRQ1 = 1
	  endif	  

	return
	end
	

	real function neLSB(x,y,z,FLSBr,wLSB)	! Local Super Bubble
	implicit none
	real x,y,z,FLSBr
	integer wLSB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLSB = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough. 
c	FLSB = fluctuation parameter
c	wLSB  = weight of LSB component used to combine
c		with other components of electron density.
c		wLSB =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLSB)*ne_other + neLSB

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlsb, ylsb, zlsb		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis 
					!    (x axis points toward l=90)
	
	real q

	real ap, bp, cp, dp
	real s, c
	
        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./
	save

c	data aa, bb, cc       	/ 0.6, 0.25, 0.3 / 	! GUESS     
c	data xlsb, ylsb, zlsb 	/ -0.7, 9.0, 0. / 	! GUESS
c	data theta		/ 150. / 		! GUESS

c	data netrough  /0.01/		! GUESS
c	data Ftrough   / 1./		! GUESS

	aa=alsb
	bb=blsb
	cc=clsb
	theta=thetalsb
	netrough=nelsb0
	Ftrough=Flsb
	
	if(first) then
	  s = sin(theta/radian) 
	  c = cos(theta/radian) 
	  ap = (c/aa)**2 + (s/bb)**2 
	  bp = (s/aa)**2 + (c/bb)**2 
	  cp = 1./cc**2
	  dp =  2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLSB = 0.
	wLSB = 0
	FLSBr = 0.
	q = (x-xlsb)**2*ap
     .    + (y-ylsb)**2*bp 
     .    + (z-zlsb)**2*cp
     .    + (x-xlsb)*(y-ylsb)*dp
	  if(q .le. 1.0) then	! inside 
	    neLSB = netrough 
	    FLSBr = Ftrough 
	    wLSB = 1
	  endif	  

	return
	end

	real function neLHB(x,y,z,FLHBr,wLHB)	! Local Hot Bubble
	implicit none
	real x,y,z,FLHBr
	integer wLHB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLHB = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough. 
c	FLHB = fluctuation parameter
c	wLHB  = weight of LBH component used to combine
c		with other components of electron density.
c		wLBH =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLHB)*ne_other + neLHB

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlhb, ylhb, zlhb		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis 
					!    (x axis points toward l=90)
	real q

	real ap, bp, cp, dp
	real s, c
	
        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./

c	data aa, bb, cc       	/ 0.15, 0.08, 0.2 / 	           
c	data xlhb, ylhb, zlhb 	/ 0., 8.5, 0. / 
c	data theta		/ 135. / 

c	data netrough  /0.005/
c	data Ftrough   / 1./
	
	save

	aa=alhb
	bb=blhb
	cc=clhb
	theta=thetalhb
	netrough=nelhb0
	Ftrough=Flhb

	
	if(first) then
	  s = sin(theta/radian) 
	  c = cos(theta/radian) 
	  ap = (c/aa)**2 + (s/bb)**2 
	  bp = (s/aa)**2 + (c/bb)**2 
	  cp = 1./cc**2
	  dp = 2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLHB = 0.
	wLHB = 0
	FLHBr = 0.
	q = (x-xlhb)**2*ap
     .    + (y-ylhb)**2*bp 
     .    + (z-zlhb)**2*cp
     .    + (x-xlhb)*(y-ylhb)*dp
	  if(q .le. 1.0) then	! inside 
	    neLHB = netrough 
	    FLHBr = Ftrough 
	    wLHB = 1
	  endif	  

	return
	end
	
	

	real function neLHB2(x,y,z,FLHBr,wLHB)	! Local Hot Bubble
c LHB modeled as a cylinder
c the cylinder slants in the y direction vs. z as described by parameter yzslope
c the cylinder cross-sectional size in the 'a' direction (major axis)
c       varies with z, tending to zero at its smallest z point.
	implicit none
	real x,y,z,FLHBr
	integer wLHB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLHB2 = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough. 
c	FLHB = fluctuation parameter
c	wLHB  = weight of LBH component used to combine
c		with other components of electron density.
c		wLHB =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLHB)*ne_other + neLHB2

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlhb, ylhb, zlhb		! center of ellipsoid
	real theta 			! slant angle in yz plane of cylinder
					!    measured from z axis 
	real qxy, qz

        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./

	real yzslope 

	real yaxis
	save

	aa=alhb
	bb=blhb
	cc=clhb
	theta=thetalhb
	netrough=nelhb0
	Ftrough=Flhb

	
	if(first) then
	  yzslope = tan(theta/radian)
	  first = .false.
	endif

	neLHB2 = 0.
	wLHB = 0
	FLHBr = 0.
	yaxis = ylhb + yzslope*z
c cylinder has cross sectional area = constant for z>0
c area -> 0 for z<0 by letting aa->0 linearly for z<0:
c (0.001 = 1 pc is to avoid divide by zero)
	if(z .le. 0. .and. z .ge. zlhb-clhb) then  
	  aa = 0.001 + (alhb-0.001)*(1. - (1./(zlhb-clhb))*z)
	else
	  aa = alhb
	endif
c	write(99, *) x, y, z, aa, bb, cc
	qxy =  ( (x-xlhb)/aa )**2 + ( (y-yaxis)/bb )**2 
	qz =  abs(z-zlhb)/cc
	if(qxy .le. 1.0 .and. qz .le. 1.0) then ! inside 
	    neLHB2 = netrough 
	    FLHBr = Ftrough 
	    wLHB = 1
	endif	  

	return
	end
	
	real function neLOOPI(x,y,z,FLOOPI,wLOOPI)	! Loop I
c component is a spheroid truncated for z<0.
	implicit none
	real x,y,z,FLOOPI
	integer wLOOPI
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLOOPI = electron density in LOOP I that
c	        is modeled as an ellisoidal trough 
c		with an enhanced shell
c	FLOOPI = fluctuation parameter
c	wLOOPI  = weight of LOOP I component used to combine
c		with other components of electron density.
c		wLOOPI =  1  at and inside the annular ridge
c		       <  1  outside the annular ridge

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/ 
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
         

	real r
	real a1, a2
	logical first
	data first /.true./
	save
	
	if(first) then
	  a1 = rlpI
	  a2 = rlpI+drlpI
	  first = .false.
	endif

        if(z .lt. 0.) then
	  neLOOPI = 0.
	  FLOOPI = 0.
	  wLOOPI = 0
          return
        endif
	r = sqrt( (x-xlpI)**2 + (y-ylpI)**2 + (z-zlpI)**2) 
        if(r .gt. a2) then 	! outside Loop I
	  neLOOPI = 0.
	  FLOOPI = 0.
	  wLOOPI = 0
	else if(r .le. a1) then	! inside volume
	    neLOOPI= nelpI
	    FLOOPI = FlpI
	    wLOOPI = 1
c           write(99,*) x,y,z, r, neLOOPI, ' inside volume'
	else			! inside boundary shell
	    neLOOPI= dnelpI
	    FLOOPI = dFlpI
	    wLOOPI = 1
c           write(99,*) x,y,z,r, neLOOPI, ' inside shell'
	endif	  

	return
	end
	
