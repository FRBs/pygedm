c density.NE2001.f
c final version of NE2001
c returns densities, F parameters and weights of the various components
c mods:
c 28 July 02:
c   put in 'if' statements in subroutine density_2001 so that
c	function calls are not done if particular weights
c	(wg1, wg2, etc.) are set to zero in gal01.inp
c	This is (a) cleaner and (b) much more efficient if
c	the clump or void component is flagged out.

      SUBROUTINE DENSITY_2001(x,y,z,
     .                ne1,ne2,nea,negc,nelism,necN,nevN,
     .                F1, F2, Fa, Fgc, Flism, FcN, FvN, 
     .                whicharm, wlism, wLDR, wLHB, wLSB, wLOOPI,
     .                hitclump, hitvoid, wvoid)
     
Cf2py intent(in) x, y, z
Cf2py intent(out) F1, F2, Fa, Fgc, Flism, FcN, FvN, ne1, ne2, nea
Cf2py intent(out) whicharm, wlism, wLDR, wLHB, wLSB, wLOOPI
Cf2py intent(out) hitclump, hitvoid, wvoid, negc, nelism, necN, nevN

c----------------------------------------------------------------------------
c  Returns seven components of the free electron density of the 
c  interstellar medium at Galactic location (x,y,z).  
c  Calling arguments:
c  input:
c	x, y, z = galactocentric location (kpc)
c       Right-handed coordinate system
c       x is in l=90 direction
c       y is in l=180 direction
c       The sun is at (x,y,z) = (0,R0,0)  
c  output:
c    electron densities in cm^{-3}:
c	ne1:	outer, thick disk 
c	ne2:	inner, thin disk (annular in form)
c	nea:	spiral arms
c	negc:   galactic center component
c       nelism: local ISM component
c       necN:   contribution from discrete 'clumps'
c       nevN:   contribution from voids
c    fluctuation parameters (one for each ne component):
c       F1, F2, Fa, Fgc, Flism, FcN, FvN 
c    flags:
c       whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
c          wlism: 1 if x,y,z is in any of the four LISM components 
c           wLDR: 1 if in LDR, 0 if not
c           wLHB: 1 if in LHB, 0 if not
c           wLSB: 1 if in LSB, 0 if not
c         wLOOPI: 1 if in LoopI, 0 if not
c       (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR) 
c       hitclump: clump number that x,y,z is in (0 if none) 
c        hitvoid: void number that x,y,z is in (0 if none)
c 25 May 2002
c based on routines from TC93 and test routines from 1999-2002 by JMC.
c----------------------------------------------------------------------------
	implicit none
        real x,y,z
        real ne1,ne2,nea,negc,nelism,necN,nevN
        real F1, F2, Fa, Fgc, Flism, FcN, FvN
        integer whicharm,wlism,wLDR,wLHB,wLSB,wLOOPI,hitclump,hitvoid
        integer wvoid

	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
	common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

	logical first
	data first/.true./
	save

c functions:
	real ne_inner
	real ne_outer
        real ne_arms_log_mod
	real ne_gc
        real ne_lism
c subroutines needed:
c	neclumpN
c	nevoidN

	if(first) then
c! get parameters first time through
	   call get_parameters
	   first=.false.
	endif

c need to initiate values in case components are flagged out
	ne1 = 0.
	ne2 = 0.
	nea = 0.
	negc = 0.
	nelism = 0.
	necN = 0.
	nevN = 0.

	wlism = 0.
	wldr = 0.
	wlhb = 0.
	wlsb = 0.
	wloopI = 0.
	hitclump = 0
	hitvoid = 0
	wvoid = 0
	whicharm = 0
	
	if(wg1 .eq. 1) ne1 = ne_outer(x,y,z, F1)
	if(wg2 .eq. 1) ne2 = ne_inner(x,y,z, F2)
	if(wga .eq. 1) nea = ne_arms_log_mod(x,y,z,whicharm,Fa)
        if(wggc .eq. 1) negc = ne_gc(x,y,z, Fgc)
        if(wglism .eq. 1) 
     .     nelism = ne_lism(x,y,z,Flism,wlism,wldr,wlhb,wlsb,wloopI)
        if(wgcN .eq. 1) call neclumpN(x,y,z,necN,FcN,hitclump)
        if(wgvN .eq. 1) call nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid)

c       write(21, "(3(f8.2,1x),5(f8.4,1x),i2)") 
c    .       x,y,z,ne1,ne2,nea,nelism,negc,
c    .       whicharm

	return
	end

	SUBROUTINE GET_PARAMETERS
c-----------------------------------------------------------------------
	implicit none
	real rsun
	common /mw/ rsun
	data rsun/8.5/
	
c control flags for turning components on and off:

	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
	common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

c parameters of large-scale components (inner+outer+arm components):
	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c factors for controlling individual spiral arms:
c	narm: 	multiplies electron density (in addition to the`fac'
c		      quantities)
c	warm:	arm width factors that multiply nominal arm width
c	harm:	arm scale height factors
c	farm:	factors that multiply n_e^2 when calculating SM

	integer narmsmax
	parameter (narmsmax=5)
	real narm, warm, harm, farm
	common/armfactors/ 
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax) 


	real negc0, Fgc0
        common /gcparms/ negc0, Fgc0

	   open(11,file='gal01.inp',status='old')
	   read(11,*)
	   read(11,*) wg1, wg2, wga, wggc, wglism, wgcN, wgvN
	   read(11,1020) n1h1,h1,A1,F1,n2,h2,A2,F2,
     .          na,ha,wa,Aa,Fa,
     .          narm(1), narm(2), narm(3), narm(4), narm(5),
     .          warm(1), warm(2), warm(3), warm(4), warm(5),
     .          harm(1), harm(2), harm(3), harm(4), harm(5),
     .          farm(1), farm(2), farm(3), farm(4), farm(5)
 1020	   format(7x,f8.0)
	   close(11)
c          write(6,*) 'get_parms: weights: ', 
c    .           wg1, wg2, wga, wggc, wglism, wgcN, wgvN
c          write(6,*) 'get_parms: ', 
c    .           n1h1,h1,A1,F1,n2,h2,A2,F2,
c    .           na,ha,wa,Aa,Fa,
c    .           narm(1), narm(2), narm(3), narm(4), narm(5),
c    .           warm(1), warm(2), warm(3), warm(4), warm(5),
c    .           harm(1), harm(2), harm(3), harm(4), harm(5),
c    .           farm(1), farm(2), farm(3), farm(4), farm(5) 
	return
	end


	REAL FUNCTION NE_ARMS_LOG_MOD(x,y,z, whicharm, Farms)
c-----------------------------------------------------------------------
c  Spiral arms are defined as logarithmic spirals using the 
c    parameterization in Wainscoat et al. 1992, ApJS, 83, 111-146.
c  But arms are modified selectively at various places to distort them
c    as needed (08 Aug 2000).
c  Note that arm numbering follows that of TC93 for the four large arms
c (after remapping).
c  The local spiral arm is number 5.
c  06 Apr 02:   removed TC type modifications of arms 2,3 (fac calculations)
c  		and replaced with new versions.  Data for these are hard wired.

	implicit none
	real x, y, z
	integer whicharm
	real Farms

	integer whicharm_spiralmodel

        real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
        common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c see get_parameters for definitions of narm, warm, harm.

        integer narmsmax
        parameter (narmsmax=5)
	real narm, warm, harm, farm
	common/armfactors/ 
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax) 


        real rsun
        common /mw/ rsun

        real rad
        parameter(rad=57.2957795130823)

	integer ks
        data ks/3/

	integer NN
	data NN/7/

	integer NNmax
	parameter(NNmax=20)

	integer narms
	parameter(narms=5)

	real aarm(narms), rmin(narms), thmin(narms), extent(narms)

	integer armmap(5)
c! for remapping from Wainscoat
	data armmap/1, 3, 4, 2, 5/
c! order to TC93 order, which is
c! from GC outwards toward Sun.
 	integer NNj(narms)
	data NNj/20, 20, 20, 20, 20/

	real th1(NNmax,narms),r1(NNmax,narms) 

	real arm
	integer kmax, narmpoints, ncoord
        parameter(narmpoints=500, ncoord=2)
	dimension arm(narms,narmpoints,ncoord),kmax(narms)

	real nea
        integer j, k, n, jj

	logical first
	data first /.true./

	real rr
 	real dth, th, r
	real smin, sminmin
	real sq, sq1, sq2, sqmin
	integer kk, kmi, kma, kl
	real emm, ebb, exx, eyy,  test
	real ga, fac
	real thxy
	real arg
	real th3a, th3b, fac3min, test3
	real th2a, th2b, fac2min, test2

	
c function:
	real sech2
	save

	rr=sqrt(x**2 + y**2)
	if(first) then
c! Reconstruct spiral arm axes

c read arm parameters:
        open(11, file='ne_arms_log_mod.inp', status='old')
c       write(6,*) 'ne_arms_log_mod.inp:'
        read(11,*)
        read(11,*)
        do j=1,narms
          read(11,*) aarm(j), rmin(j), thmin(j), extent(j)
c         write(6,*) aarm(j), rmin(j), thmin(j), extent(j)
        enddo
        close(11)

	do j=1,narms
c! fill sampling array
	  do n=1,NNj(j)
	   th1(n,j) = thmin(j)+(n-1)*extent(j)/(NNj(j)-1.)
c! rad
	   r1(n,j) = rmin(j)*exp((th1(n,j)-thmin(j))/aarm(j))
	   th1(n,j) = th1(n,j)*rad
c! deg
c *** begin sculpting spiral arm 2 == TC arm 3***
           if(armmap(j) .eq. 3) then
 	   if(th1(n,j) .gt. 370. .and. th1(n,j) .le. 410.) then
     	      r1(n,j) = r1(n,j) * 
     .           (1. + 0.04*cos((th1(n,j)-390.)*180./(40.*rad)))
c    .           (1. + 0.01*cos((th1(n,j)-390.)*180./(40.*rad)))
	   endif
 	   if(th1(n,j) .gt. 315. .and. th1(n,j) .le. 370.) then
     	      r1(n,j) = r1(n,j) * 
     .           (1. - 0.07*cos((th1(n,j)-345.)*180./(55.*rad)))
c    .           (1.0 - 0.08*cos((th1(n,j)-345.)*180./(55.*rad)))
	   endif
 	   if(th1(n,j) .gt. 180. .and. th1(n,j) .le. 315.) then
	      r1(n,j) = r1(n,j) * 
c    ,           (1 + 0.13*cos((th1(n,j)-260.)*180./(135.*rad)))
     ,           (1 + 0.16*cos((th1(n,j)-260.)*180./(135.*rad)))
	   endif
	   endif
c *** begin sculpting spiral arm 4 == TC arm 2***
           if(armmap(j) .eq. 2) then
 	   if(th1(n,j) .gt. 290. .and. th1(n,j) .le. 395.) then
     	      r1(n,j) = r1(n,j) * 
c    .            1.
     .           (1. - 0.11*cos((th1(n,j)-350.)*180./(105.*rad)))
	   endif
	   endif
c *** end arm sculpting ***
c	   write(6,*) j,n, th1(n,j), r1(n,j)
	  enddo
	enddo
c        open(11,file='log_arms.out', status='unknown')
c        write(11,*) 'arm  n   xa     ya'

	   do 21 j=1,narms
	      dth=5.0/r1(1,j)
	      th=th1(1,j)-0.999*dth
	      call cspline(th1(1,j),r1(1,j),-NNj(j),th,r)
c	      write(6,*) 'doing arm ', j, ' with ', NNj(j), ' points',
c    .               dth
c	      write(6,*) (th1(k,j), r1(k,j), k=1,NNj(j))
	      do 10 k=1,narmpoints-1
		 th=th+dth
		 if(th.gt.th1(NNj(j),j)) go to 20
		 call cspline(th1(1,j),r1(1,j),NNj(j),th,r)
		 arm(j,k,1)=-r*sin(th/rad)
		 arm(j,k,2)= r*cos(th/rad)
                 write(11,"(1x,i2,1x,i3,1x,2(f7.3,1x))") 
     .              j,k,arm(j,k,1),arm(j,k,2)
 10	      continue
 20	   continue
	   kmax(j)=k
 21        continue
	   close(11)

	first = .false.
	endif
c
c Get spiral arm component:  30 do loop finds a coarse minimum distance
c from line of sight to arm; 40 do loop finds a fine minimum distance
c from line of sight to arm; line 35 ensures that arm limits are not 
c exceeded; linear interpolation beginning at line 41 finds the 
c minimum distance from line of sight to arm on a finer scale than gridding
c of arms allows (TJL)

	nea=0.0
        ga = 0.
	whicharm = 0
        whicharm_spiralmodel = 0
        sminmin = 1.e10
	thxy = atan2(-x, y) * rad
c! measured ccw from +y axis
c! (different from tc93 theta)
	if(thxy.lt.0.) thxy=thxy+360.
	if(abs(z/ha).lt.10.) then
           do 50 j=1,narms
              jj = armmap(j)
              sqmin=1.e10
              do 30 k=1+ks,kmax(j)-ks,2*ks+1
                 sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
                 if(sq.lt.sqmin) then
                    sqmin=sq
                    kk=k
                 endif
 30           continue
 35           kmi = max(kk-2*ks, 1)
              kma = min(kk+2*ks, kmax(j))
              do 40 k=kmi,kma
                 sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
                 if(sq.lt.sqmin) then
                    sqmin=sq
                    kk=k
                 endif
 40           continue
 41           if (kk.gt.1.and.kk.lt.kmax(j)) then
                 sq1 = (x - arm(j,kk-1,1))**2 + (y - arm(j,kk-1,2))**2
                 sq2 = (x - arm(j,kk+1,1))**2 + (y - arm(j,kk+1,2))**2
                 if (sq1.lt.sq2) then 
                    kl = kk - 1
                 else
                    kl = kk + 1
                 endif
                 emm = (arm(j,kk,2) - arm(j,kl,2))
     $                /(arm(j,kk,1) - arm(j,kl,1))
                 ebb = arm(j,kk,2) - emm*arm(j,kk,1)
                 exx = (x + emm*y - emm*ebb)/(1.0 + emm**(2))
                 test = (exx - arm(j,kk,1))/(arm(j,kl,1) - arm(j,kk,1))
                 if (test.lt.0.0.or.test.gt.1.0) exx = arm(j,kk,1)
                 eyy = emm*exx + ebb
              else
                 exx = arm(j,kk,1)
                 eyy = arm(j,kk,2)
              endif
              sqmin = (x - exx)**(2) + (y - eyy)**(2)
	      smin=sqrt(sqmin)
c! Distance of (x,y,z) from arm axis
c           write(23,"(4(f5.2,1x),i2,1x,3(f8.3,1x))") 
c    .        x,y,z,rr,j,exx,eyy,smin
	    if(smin.lt.3*wa) then
c! If (x,y,z) is close to this
	      ga=exp(-(smin/(warm(jj)*wa))**2)
c! arm, get the arm weighting factor
	      if(smin .lt. sminmin) then
		whicharm_spiralmodel = j
		sminmin = smin
	      endif
	      if(rr.gt.Aa) then
c! Galactocentric radial dependence of arms
		ga=ga*sech2((rr-Aa)/2.0) 
c		write(6,*) 'd99a: rr,Aa,sech2() = ', 
c                 rr, Aa, sech2((rr-Aa)/2.0) 
	      endif


c arm3 reweighting:
	      th3a=320.
	      th3b=390.
	      th3b=370.
	      th3a=290.
	      th3b=363.
	      th3b=363.
	      fac3min=0.0
	      test3 = thxy-th3a
	      if(test3 .lt.0.) test3 = test3+360.
	      if(jj.eq.3
     .            .and. 0. .le. test3 
     .            .and. test3 .lt. (th3b-th3a))
     .        then
	        arg=6.2831853*(thxy-th3a)/(th3b-th3a)
c		fac = (3.0 + cos(arg))/4.0
		fac = (1.+fac3min + (1.-fac3min)*cos(arg))/2.
		fac = fac**4.0
c		write(90,*) x, y, thxy, th3a, th3b, test3, fac
		ga=ga*fac
	      endif
	  

c arm2 reweighting:
c    first: as in tc93 (note different definition of theta)
	      th2a=35.
	      th2b=55.
	      test2 = thxy-th2a
	      fac = 1.
              
	      if(jj.eq.2
     .            .and. 0. .le. test2 
     .            .and. test2 .lt. (th2b-th2a))
     .        then
	         fac=1.+ test2/(th2b-th2a)
		 fac = 1.
c!**** note turned off
                 ga=ga*fac
	      endif

	      if (jj.eq.2 .and. test2 .gt. (th2b-th2a)) then
	        fac = 2. 
		fac = 1.
c!**** note turned off
                ga=ga*fac
	      endif
c    second:  weaken the arm in a short range:
	      th2a=340.
	      th2b=370.
c note fix does nothing if fac2min = 1.0
	      fac2min=0.1
	      test2 = thxy-th2a
	      if(test2 .lt.0.) test2 = test2+360.
              
	      if(jj.eq.2
     .            .and. 0. .le. test2 
     .            .and. test2 .lt. (th2b-th2a))
     .        then
	        arg=6.2831853*(thxy-th2a)/(th2b-th2a)
		fac = (1.+fac2min + (1.-fac2min)*cos(arg))/2.
c		fac = fac**3.5
c		write(90,*) x, y, thxy, th2a, th2b, test2, fac
		ga=ga*fac
	      endif
	  
	      nea=nea + 
     .            narm(jj)*na*ga*sech2(z/(harm(jj)*ha))   ! Add this arm contribution
	    endif
50	  continue
	endif

        ne_arms_log_mod = nea
        Farms = 0
	if(whicharm_spiralmodel .eq. 0) then
	  whicharm = 0
	else
	  whicharm = armmap(whicharm_spiralmodel)	! remap arm number
	  Farms = Fa * farm(whicharm)
	endif
	return
	end



	REAL FUNCTION NE_OUTER(x,y,z, F_outer)
c-----------------------------------------------------------------------
c Thick disk component:
	implicit none
	real x,y,z, F_outer

	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa


        real pihalf, rad
	parameter(pihalf=1.5707963267948966)
	parameter(rad=57.2957795130823)

	real rsun
	common /mw/ rsun

	real g1
	real sech2
	real rr, suncos, ne1

	logical first
	data first /.true./
	save

c 	g1=sech2(rr/A1)/sech2(8.5/A1)		! TC93 function
	rr=sqrt(x**2 + y**2)
	suncos = cos(pihalf*rsun/A1)
	if (rr.gt.A1) then 
	   g1 = 0.0
	else
	   g1 = cos(pihalf*rr/A1)/suncos
	endif
 	ne1=(n1h1/h1)*g1*sech2(z/h1)	
	ne_outer = ne1
        F_outer = F1

	return
	end

	REAL FUNCTION NE_INNER(x,y,z, F_inner)
c-----------------------------------------------------------------------
c Thin disk (inner Galaxy) component:
c (referred to as 'Galactic center component' in circa TC93 density.f)
	implicit none
	real x,y,z, F_inner

	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa
	real g2, rr, rrarg
	real sech2
	real ne2
	save

	g2=0.0
	rr=sqrt(x**2 + y**2)
	rrarg=((rr-A2)/1.8)**2
	if(rrarg.lt.10.0) g2=exp(-rrarg)
	ne2=n2*g2*sech2(z/h2)
 
	ne_inner = ne2
        F_inner = F2
	return
	end


      REAL FUNCTION NE_GC(x, y, z, F_gc)
c-----------------------------------------------------------------------
c     Determines the contribution of the Galactic center to the free 
c     electron density of the interstellar medium at Galactic location 
c     (x,y,z).  Combine with `fluctuation' parameter to obtain the 
c     scattering measure.
c
c     NOTE: This is for the hyperstrong scattering region in the 
c     Galactic center.  It is distinct from the inner Galaxy 
c     (component 2) of the TC93 model.
c
c     Origin of coordinate system is at Galactic center; the Sun is at 
c     (x,y,z) = (0,R0,0), x is in l=90 direction
c
c     Based on Section 4.3 of Lazio & Cordes (1998, ApJ, 505, 715)
c
c Input:
c REAL X - location in Galaxy [kpc]
c REAL Y - location in Galaxy [kpc]
c REAL Z - location in Galaxy [kpc]
c
c COMMON:
c REAL NEGC0 - nominal central density
c
c PARAMETERS:
c REAL RGC - radial scale length of Galactic center density enhancement
c REAL HGC - z scale height of Galactic center density enhancement
c
c Output:
c REAL NE_GC - Galactic center free electron density contribution [cm^-3]
c-----------------------------------------------------------------------
c
      implicit none
      real x, y, z, F_gc
 
      real rgc, hgc
      real xgc, ygc, zgc
c     parameter (xgc=-0.010, ygc=0., zgc=-0.020)
c     parameter (rgc=0.145)
c     parameter (hgc=0.026)
 
      real rr, zz 
      real arg

      real negc0, Fgc0
      common /gcparms/ negc0, Fgc0

      real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
      common /galparams/ n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa 

      logical first
      data first /.true./
      save

      ne_gc = 0.
      F_gc = 0.

      if(first) then
        open(11, file='ne_gc.inp', status='old')
          read(11,*)
          read(11,*) xgc, ygc, zgc
	  read(11,*) rgc
	  read(11,*) hgc
	  read(11,*) negc0
          read(11,*) Fgc0
        close(11)
	first = .false.
      endif

 
c GALACTOCENTRIC RADIUS
      
      rr = sqrt( (x-xgc)**2 + (y-ygc)**2)
      if(rr .gt. rgc) return				! truncate at 1/e point 

c Z-HEIGHT.

      zz =abs(z-zgc)
      if(zz .gt. hgc) return 
      arg = (rr/rgc)**2 + (zz/hgc)**2
      if(arg .le. 1.) then
         ne_gc = negc0
	 F_gc = Fgc0
c        write(21,*) 'ne_gc: rr,zz,arg,ne_gc,F_gc ', 
c    .                rr, zz, arg, ne_gc, F_gc
      endif

      return
      end

      
c
c%%%%%%%%%%%%%%%%%%%%%%%%%  cspline.f  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine cspline(x,y,nn,xout,yout)
	integer nnmax
	parameter(nnmax=20)
	real x(nnmax),y(nnmax),y2(nnmax),u(nnmax)
	save

	if(nn .gt. nnmax) then
	  write(6,*) 
     .    ' too many points to spline. Change parameter statement'
	  write(6,*) 
     .    ' in cspline'
	endif

	n=abs(nn)
	if(nn.lt.0) then
	  y2(1)=0.
	  u(1)=0.
	  do 10 i=2,n-1
	    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	    p=sig*y2(i-1)+2.
	    y2(i)=(sig-1.)/p
	    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
10        continue
	  qn=0.
	  un=0.
	  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	  do 20 k=n-1,1,-1
20	  y2(k)=y2(k)*y2(k+1)+u(k)
	endif

	klo=1
	khi=n
30	if (khi-klo.gt.1) then
	  k=(khi+klo)/2
	  if(x(k).gt.xout)then
	    khi=k
	  else
	    klo=k
	  endif
	goto 30
	endif
	h=x(khi)-x(klo)
	if (h.eq.0.) write(*,*) 'bad x input.'
	a=(x(khi)-xout)/h
	b=(xout-x(klo))/h
	yout=a*y(klo)+b*y(khi)+
     +    ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.
	return
	end

	REAL FUNCTION SECH2(z)
c-----------------------------------------------------------------------
	sech2=0.0
	if(abs(z).lt.20.0) sech2=(2.0/(exp(z)+exp(-z)))**2
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012
