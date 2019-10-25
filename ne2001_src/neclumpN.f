	subroutine neclumpN(x,y,z,necN,FcN,hitclump)
c returns electron density necN and fluctuation parameter FcN 
c at position designated by l,b,d,x,y,z c for a set of  
c clumps with parameters read in from file  neclumpN.dat

c input:
c	x,y,z	coordinates	(kpc)  (as in TC93)
c
c output:
c	necN	electron density in clump at (x,y,z)
c	FcN	fluctuation parameter
c 	hitclump = 0:   no clump hit
c		   j>0: j-th clump hit

	implicit none
	real x,y,z,necN,FcN
	integer hitclump

	integer nclumpsmax
	parameter (nclumpsmax=2000)

c	character*15 losname(nclumpsmax)
c	character*1 type(nclumpsmax)
	real lc(nclumpsmax), bc(nclumpsmax), dc(nclumpsmax)
        real xc(nclumpsmax), yc(nclumpsmax), zc(nclumpsmax)
        real nec(nclumpsmax), rc(nclumpsmax), Fc(nclumpsmax)
        integer edge(nclumpsmax)
	integer nclumps
	integer hitclumpflag
	common /clumps/ nclumps, hitclumpflag(nclumpsmax)

c parameters:
c	lc	= galactic longitude of clump center
c	bc	= galactic latitude of clump center
c	(xc,yc,zc) = clump center location (calculated)
c       nec	= internal peak electron density
c	rc	= clump radius at 1/e
c       Fc      = clump fluctuation parameter
c	edge    = 0 => use exponential rolloff out to 5rc
c                 1 => uniform and truncated at 1/e


	real radian
	parameter(radian = 57.29577951)

	real rsun
	parameter (rsun=8.5)

	logical first
	data first/.true./

	real slc, clc, sbc, cbc
	real rgalc


	real arg

	integer luclump
	data luclump/11/
	integer j
	integer clumpflag

	save

c first time through, read input clump parameters and calculate
c LOS quantities. 
c lc,bc = Galactic coordinates (deg)
c   nec = clump electron density (cm^{-3})
c    Fc = fluctuation parameter
c    dc = clump distance from Earth (kpc)
c    rc = clump radius (kpc)
c  edge = 0,1  0=> Gaussian, 1=> Gaussian w/ hard edge at e^{-1} 
c  type = LOS type (P pulsar, G other Galactic, X extragalactic
c losname = useful name
	if(first) then 		!read clump parameters
	  j=1
c	  write(6,*) 'reading neclumpN.NE2001.dat'
	  open(luclump, file='neclumpN.NE2001.dat', status='old')
	  read(luclump,*)				! label line
    5     read(luclump,*,end=99) clumpflag,lc(j),bc(j),nec(j),Fc(j),
     .           dc(j),rc(j),edge(j)
          if(clumpflag .eq. 0) then
	    slc = sin(lc(j)/radian)
	    clc = cos(lc(j)/radian)
	    sbc = sin(bc(j)/radian)
	    cbc = cos(bc(j)/radian)
	    rgalc = dc(j)*cbc
	    xc(j) = rgalc*slc
	    yc(j) = rsun-rgalc*clc
	    zc(j) = dc(j)*sbc
c	  write(6,"(a15,1x,8(f8.3,1x))") 
c    .           losname(j),lc(j),bc(j),dc(j),
c    .           nec(j),Fc(j),xc(j),yc(j),zc(j)
	    j=j+1
	  endif
	  go to 5
   99     continue
	  first = .false.
	  nclumps = j-1
	  close(luclump)
	endif

	necN = 0.
	hitclump = 0
	FcN = 0.
	do j=1,nclumps
	  arg = 
     .      ((x-xc(j))**2 + (y-yc(j))**2 + (z-zc(j))**2) / rc(j)**2
	  if(edge(j) .eq. 0 .and. arg .lt. 5.) then
	    necN = necN + nec(j) * exp(-arg)
            FcN = Fc(j)
	    hitclump = j
            hitclumpflag(j) = 1
	  endif
	  if(edge(j) .eq. 1 .and. arg .le. 1.) then
c    	    necN = necN + nec(j) * exp(-arg)
	    necN = necN + nec(j)
            FcN = Fc(j)
	    hitclump = j
            hitclumpflag(j) = 1
	  endif
	enddo
	

	return
	end




