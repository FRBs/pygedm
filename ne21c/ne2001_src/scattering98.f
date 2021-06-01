c this version from 18 March 1998 uses revised coefficients 
c that are consistent with Cordes \& Rickett (1998, ApJ, submitted)
c modifications:
c	28 March 2001: added FUNCTION TRANSITION_FREQUENCY 
 
      REAL FUNCTION TAUISS(d, sm, nu)
c
c calculates the pulse broadening time in ms
c from distance, scattering measure, and radio frequency
c
c input:      d = pulsar distance       (kpc)    
c            sm = scattering measure    (kpc m^{-20/3})
c            nu = radio frequency       (GHz)
c output: tauss = pulse broadening time (ms) 
c
      implicit none
      real d, sm,  nu
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
      end
c
c
      REAL FUNCTION SCINTBW(d, sm, nu)
c
c calculates the scintillation bandwidth in kHz 
c from distance, scattering measure, and radio frequency
c
c input:        d = pulsar distance       (kpc)    
c              sm = scattering measure    (kpc m^{-20/3})
c              nu = radio frequency       (GHz)
c output: scintbw = scintillation bandwidth (kHz)
c
      implicit none
      real d, sm, nu
      real c1
      parameter(c1=1.16)		! for uniform, Kolmogorov medium
      real tauiss
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)	! ms
      scintbw = c1 / (2. * 3.14159 * tauiss)			! kHz
      end
 
      REAL FUNCTION SCINTIME(sm, nu, vperp)
c
c calculates the scintillation speed for given distance, galactic
c longitude and latitude, frequency, and transverse velocity      
c
c input:   sm = scattering measure	(kpc m^{-20/3})
c          nu = radio frequency 	(GHz)
c       vperp = psr transverse speed  	(km/s)  
c
c output: scintime = scintillation time (sec)
c 
c usage: should be called with sm = smtau for appropriate
c        line of sight weighting
c reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123.
c
      implicit none
      real sm, nu, vperp
c nb: formerly, the coeff. in the following line was 2.3 from 
c     Cordes & Lazio (1991)
      scintime = 3.3 * nu**1.2 * sm**(-0.6) * (100./vperp)
      end
 
 
      REAL FUNCTION SPECBROAD(sm, nu, vperp)
c
c calculates the bandwdith of spectral broadening
c for given scattering measure, , frequency, and transverse velocity      
 
c input:   sm = scattering measure	(kpc m^{-20/3})
c          nu = radio frequency 	(GHz)
c       vperp = psr transverse speed  	(km/s)  
 
c output: specbroad = spectral broadening bandwidth (Hz)
c 
c usage: should be called with sm = smtau for appropriate
c        line of sight weighting
c reference: eqn (47) of Cordes & Lazio 1991, ApJ, 376, 123.
c
      implicit none
      real sm, nu, vperp
c nb: the coeff. in the following line is 0.14 Hz  from Cordes & Lazio (1991)
c it is changed to 0.097 to conform with FUNCTION SCINTIME and
c a new calculation consistent with Cordes & Rickett (1998)

      specbroad = 0.097 * nu**(-1.2) * sm**0.6 * (vperp/100.)	! Hz
      end
 
 
      REAL FUNCTION THETA_XGAL(sm, nu)
c
c calculates angular broadening for an extragalactic
c source of plane waves
c
c sm = scattering measure
c nu = radio frequency
c theta_xgal = angular broadening FWHM (mas)
c
      implicit none
      real sm, nu
      theta_xgal = 128. * sm**0.6 * nu**(-2.2)
      end
c
      REAL FUNCTION THETA_GAL(sm, nu)
c
c calculates angular broadening for a galactic
c source of spherical waves
c
c sm = scattering measure
c nu = radio frequency
c theta_gal = angular broadening FWHM (mas)
c
      implicit none
      real sm, nu
      theta_gal = 71. * sm**0.6 * nu**(-2.2)
      end
c
      FUNCTION EM (sm)
c
c units of sm are kpc m^{-20/3}
c units of em are pc cm^{-6}      
c
c calculates the emission measure from the scattering measure
c using an assumed outer scale and spectral index of the
c wavenumber spectrum.
c
c for a wavenumber spectrum P_n(q) = q^{-alpha} from q_0 to q_1
c the mean square electron density is
c
c <n_e^2> =~  4pi*[C_n^2 / (alpha - 3) ] * q_0^{3 - alpha)
c
c ( an approximate form that assumes (q_0 / q_1)^{3-alpha} >> 1.
c
c Jim Cordes 18 Dec 1989
c
      data router /1./	! outer scale = 1 pc
      data pc/3.086e+18/
      data alpha/3.6666667/
      data pi/3.14159/
c
      em = sm *
     1         ( (4. * pi * 1000.) / (alpha - 3.) ) * 
     2         (router*pc / (2. * 3.14159) )**(alpha-3.) *
     3         (0.01) ** (20./3.)
c
      return
      end        

      REAL FUNCTION THETA_ISO(smiso, nu)
      real smiso, nu

c smiso in (kpc m^{-20/3}) x kpc^{5/3}
c    nu in GHz
c returns the isoplanatic angle in microarcsec
c 12 October 1998 
c JMC


c \theta_{iso} = \delta r_s / d
c              = \left [
c	         (\lambda r_e)^2 f_{\alpha} SM_{iso}
c		 \right ]^{1/\alpha}
c where \alpha = 5/3 for Kolmogorov case.
c NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
c    so SM_{iso} does not have the units of scattering
c    measure, but rather units of SM x Length^{\alpha}
c
c f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
c for \alpha = 5/3, f_{\alpha}= 88.3
c
c     real r_e
c     parameter(r_e = 2.82e-13)			!cm
c     real kpc
c     parameter(kpc = 3.086e21)			!cm
c     real falpha
c     parameter(falpha=88.3)

      theta_log_radian = 
     .    13.287  				! 0.6*log10(30cm*r_e)
     .  + 1.2 * alog10(nu) 
     .  - 1.1676				! 0.6*log10(f_alpha)    
     .  - 0.6 * alog10(smiso)
     .  - 34.383				! 1.6 * alog10(kpc)
     .  + 8.					! -(20/3)*log(100)
      theta_log_microarcsec = 
     .      theta_log_radian + 11.314425	! 11.314425=alog10(microarsec/rad)
      theta_iso = 10.**theta_log_microarcsec
      return
      end



      REAL FUNCTION THETA_ISO_TEST(smiso, nu)
      real smiso, nu

c smiso in (kpc m^{-20/3}) x kpc^{5/3}
c    nu in GHz
c returns the isoplanatic angle in microarcsec
c 12 October 1998 
c JMC


c \theta_{iso} = \delta r_s / d
c              = \left [
c	         (\lambda r_e)^2 f_{\alpha} SM_{iso}
c		 \right ]^{1/\alpha}
c where \alpha = 5/3 for Kolmogorov case.
c NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
c    so SM_{iso} does not have the units of scattering
c    measure, but rather units of SM x Length^{\alpha}
c
c f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
c for \alpha = 5/3, f_{\alpha}= 88.3
c
c     real r_e
c     parameter(r_e = 2.82e-13)			!cm
c     real kpc
c     parameter(kpc = 3.086e21)			!cm
c     real falpha
c     parameter(falpha=88.3)

      theta_log_radian = 
     .    13.287  				! 0.6*log10(30cm*r_e)
     .  + 1.2 * alog10(nu) 
     .  - 1.1676				! 0.6*log10(f_alpha)    
     .  - 0.6 * alog10(smiso)
     .  - 34.383				! 1.6 * alog10(kpc)
     .  + 8.					! -(20/3)*log(100)
      theta_log_microarcsec = 
     .      theta_log_radian + 11.314425	! 11.314425=alog10(microarsec/rad)
      theta_iso_test = 10.**theta_log_microarcsec
c     write(6,*) 'smiso, nu = ', smiso, nu
c     write(6,*) 'theta_log_radian = ', theta_log_radian
c     write(6,*) 'theta_log_microarcsec = ', theta_log_microarcsec
c     write(6,*) 'theta_iso = ', theta_iso_test
      return
      end

      REAL FUNCTION TRANSITION_FREQUENCY(sm, smtau, smtheta, dintegrate)
      implicit none
      real sm, smtau, smtheta, dintegrate
c returns the transition frequency between weak and strong scattering
c 28 March 2001   
c JMC

c input:
c (all sm values in (kpc m^{-20/3})) 
c	    sm = int[\cnsq]
c        smtau = int[(s/D)(1-s/D) \cnsq]
c      smtheta = int[ (1-s/D) \cnsq]
c   dintegrate = distance used to integrate \cnsq (kpc) 
c output:
c   transition_frequency = GHz given by
c 	\nu_t  = 318 GHz \xi^{10/17}  SM^{6/17} D_{eff}^{5/17}
c   where 
c        D_{eff} = effective path length through medium
c        D_{eff} = \int_0^dintegrate ds s \cnsq / \int_0^dintegrate ds  \cnsq
c
c   Note we can calculate D_{eff} using
c        D_{eff} = dintegrate * (sm - smtau/6 - smtau/3) / sm
c
      real deff
      real xi
      parameter(xi= 0.3989)   ! (2.*pi)^{-1/2} = fresnel scale definition factor
      real coefficient
      parameter(coefficient=318.)	     ! GHz; see NE2001 paper
      deff = (dintegrate*(sm - smtau/6. - smtheta/3.)) / sm 
      transition_frequency = 
     .     coefficient * xi**(10./17.) * sm**(6./17.) * deff**(5./17.) 
      return
      end

