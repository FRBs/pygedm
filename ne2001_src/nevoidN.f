	subroutine nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid)
c returns electron density nevN and fluctuation parameter FvN 
c at position designated by l,b,d,x,y,z c for a set of  
c voids with parameters read in from file  nevoidN.dat

c input:
c	x,y,z	coordinates	(kpc)  (as in TC93)
c
c output:
c	nevN	electron density in void at (x,y,z)
c	FvN	fluctuation parameter
c 	hitvoid =   0:   no void hit
c		  j>0:   j-th void hit
c	wvoid = 0,1:	 void weight	

	implicit none
	real x,y,z,nevN,FvN
	integer hitvoid, wvoid

	integer nvoidsmax
	parameter (nvoidsmax=2000)

c	character*12 losname(nvoidsmax)
	real lv(nvoidsmax), bv(nvoidsmax), dv(nvoidsmax)
        real nev(nvoidsmax), Fv(nvoidsmax)
        real aav(nvoidsmax), bbv(nvoidsmax), ccv(nvoidsmax)
        real thvy(nvoidsmax), thvz(nvoidsmax)

        real xv(nvoidsmax), yv(nvoidsmax), zv(nvoidsmax)
	real c1(nvoidsmax), s1(nvoidsmax), c2(nvoidsmax), s2(nvoidsmax),
     .       cc12(nvoidsmax), ss12(nvoidsmax), 
     .       cs21(nvoidsmax), cs12(nvoidsmax)
        integer edge(nvoidsmax)
	integer nvoids
        integer hitvoidflag
	common /voids/ nvoids, hitvoidflag(nvoidsmax)

c parameters:
c	lv	= galactic longitude of void center
c	bv	= galactic latitude of void center
c	dv	= distance from Sun of void center
c	(xv,yv,zv) = void center location (calculated)
c       nev	= internal peak electron density
c       Fv      = void fluctuation parameter
c	aav	= void major axis at 1/e
c	bbv	= void minor axis at 1/e
c	ccv	= void minor axis at 1/e
c	thvy	= rotation axis of void about y axis
c	thvz	= rotation axis of void about z axis
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
	real q
	real dx, dy, dz
	real th1, th2

	integer luvoid
	data luvoid/11/
	integer j
	integer voidflag

	save

c first time through, calculate xc, yc, zc

	if(first) then 		!read void parameters
	  j=1
c	  write(6,*) 'reading nevoidN.dat.clean'
	  open(luvoid, file='nevoidN.NE2001.dat', status='old')
	  read(luvoid,*)				! label line
    5     read(luvoid,*,end=99) voidflag, 
     .      lv(j),bv(j),dv(j),				! deg, deg, kpc
     .      nev(j),Fv(j),				! cm^{-3}, dimensionless
     .      aav(j),bbv(j),ccv(j),thvy(j), thvz(j),  	! kpc,kpc,kpc,deg,deg
     .      edge(j)					! 0 or 1
          if(voidflag .eq. 0) then
	    slc = sin(lv(j)/radian)
	    clc = cos(lv(j)/radian)
	    sbc = sin(bv(j)/radian)
	    cbc = cos(bv(j)/radian)
	    rgalc = dv(j)*cbc
	    xv(j) = rgalc*slc
	    yv(j) = rsun-rgalc*clc
	    zv(j) = dv(j)*sbc
            th1=thvy(j)
            th2=thvz(j)
            s1(j) = sin(th1/radian) 
            c1(j) = cos(th1/radian) 
            s2(j) = sin(th2/radian) 
            c2(j) = cos(th2/radian) 
            cc12(j) = c1(j)*c2(j)
            ss12(j) = s1(j)*s2(j)
            cs21(j) = c2(j)*s1(j)
            cs12(j) = c1(j)*s2(j)

c	  write(6,"(a12,1x,13(f7.3,1x))") 
c    .           losname(j),lv(j),bv(j),dv(j),
c    .           nev(j),Fv(j),xv(j),yv(j),zv(j),
c    .           aav(j), bbv(j), ccv(j),
c    .           th1, th2
	    j=j+1
	  endif
	  go to 5
   99     continue
	  first = .false.
	  nvoids = j-1
	  close(luvoid)
	endif


	nevN = 0.
	FvN = 0.
	hitvoid = 0
	wvoid = 0
c note rotation matrix in the 'q = ' statement below
c corresponds to \Lambda_z\Lambda_y
c where \Lambda_y = rotation around y axis
c       \Lambda_z = rotation around z axis
c defined as
c \Lambda_y =  c1  0  s1 
c               0  1   0
c             -s1  0  c1

c \Lambda_z =  c2 s2   0 
c             -s2 c2   0
c               0  0   1
c =>
c \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2
c                      -s2*c1   c2  -s1*s2
c                         -s1    0      c1  

c so the rotation is around the y axis first, then the z axis
	do j=1,nvoids
	  dx = x-xv(j)
	  dy = y-yv(j)
	  dz = z-zv(j)
	  q = ( cc12(j)*dx + s2(j)*dy + cs21(j)*dz)**2 / aav(j)**2 
     .      + (-cs12(j)*dx + c2(j)*dy - ss12(j)*dz)**2 / bbv(j)**2
     .      + (  -s1(j)*dx         +      c1(j)*dz)**2 / ccv(j)**2 
	  if(edge(j) .eq. 0 .and. q .lt. 3.) then
	    nevN = nev(j) * exp(-q)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	  if(edge(j) .eq. 1 .and. q .le. 1.) then
	    nevN = nev(j)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	enddo
	
	if(hitvoid .ne. 0) wvoid = 1
	

	return
	end




