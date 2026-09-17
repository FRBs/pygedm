# mwprop.ne2001p v1.0 Jan 2024

"""
Constants needed for interstellar wave propagation and other applications.
All units are cgs unless indicated otherwise

Note: could keep the constants as 'fully-fledged Quantity objects'
as stated at https://docs.astropy.org/en/stable/constants/
But then each usage would have to extract a value. 
Maybe that would be ok at the beginning of modules that 
need specific constants.

Example: c = constants.c  instead of c = constants.c.cgs.value
and then: c = constants.c.cgs.value if c needs to be cgs in a module

If done this way, would put evaluations of DM and RM constants, etc. 
into a separate module (perhaps)

Address this later.
"""
from __future__ import print_function
# 30 Jun 2021 added Plank's constant
# 31 Dec 2019

from math import pi
from math import sqrt
from astropy import constants

# fundamental

eV = 1.60218e-12    # erg
c = constants.c.cgs.value
m_e = constants.m_e.cgs.value
e_e = constants.e.esu.value
r_e = e_e**2 / (m_e * c**2)

m_p = constants.m_p.cgs.value

h = constants.h.cgs.value 
hbar = constants.hbar.cgs.value 

kb = constants.k_B.cgs.value

# lengths

cm_in_m = 100.
km = 1.e5        # cm
au = constants.au.cgs.value
pc = constants.pc.cgs.value
kpc = 1.e3*pc
Mpc = 1.e6*pc
Gpc = 1.e9*pc

# angular
mas = pi / (180.*3600.*1000.)
uas = mas / 1000.
all_sky_sr = 4.*pi
all_sky_degsq = (180./pi)**2 * all_sky_sr

# frequencies

GHz = 1.e9
MHz = 1.e6

lambda1 = c / GHz        # wavelength in cm for 1 GHz

# temporal
ns = 1.e-9       # sec
mus = 1.e-6      # sec
ms = 1.e-3       # sec

# Dispersion constant
KDM = (c * r_e * pc) / (2.*pi*GHz**2)                # for TOAs in seconds

# Faraday rotation measure constant

microGauss = 1.e-6

# RMconstant gives standard rad m^{-2} 
# for B in microgauss, pathlength in pc, electron density in cm^{-3}

RMconstant = ((e_e**3) / (2.*pi*m_e**2 * c**4)) * pc * cm_in_m**2 * microGauss
RMconstant_cgs = ((e_e**3) / (2.*pi*m_e**2 * c**4)) * pc * microGauss

# Scattering measure kpc m^{-20/3} to cgs 
SMunit = kpc*10.**(-40./3.)
SM1 = 1.        # kpc m^{-20/3}
SM_cgs = SM1 * kpc * cm_in_m**(-20./3.)     # same as SMunit (a check)

# Flux density

Jy = 1.e-23     # cgs
Jansky = Jy

# Brightness temperature for 1 Jy pulse of 1 ms at 1 kpc at 1 GHz
Tb111 = (Jy * kpc**2) / (2.*kb*(GHz*ms)**2)

# plasma frequency in Hz for unit electron density
ne1 = 1.        # 1/cc
nup = sqrt((4.*pi * ne1 * e_e**2) / m_e) / (2.*pi)      # Hz

# cyclotron frequency in a 1 Gauss field
B1=1.
nuc = (e_e * B1) / (m_e * c) / (2.*pi)                  # Hz

