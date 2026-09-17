# mwprop.ne2001p v1.0 Jan 2024

from __future__ import print_function
from numpy import sqrt
from math import pi
from scipy import constants

"""
Gives constants needed for interstellar wave propagation and other applications

Units: cgs
"""

# Solar mass
# from astropy.constants.M_sun.cgs.value (scipy.constants doesn't have!)
msun = 1.9884754153381438e+33

# fundamental
eV = constants.eV / constants.erg
c = constants.c *  100.         # cm
r_e = constants.physical_constants['classical electron radius'][0]*100. # cm
m_e = constants.constants.m_e*1000
e_e = sqrt(m_e * c**2 * r_e)
kb = constants.k / constants.erg
G = constants.G / constants.gram 

h = constants.h / constants.erg
hbar = constants.hbar / constants.erg

# lengths 

cm_in_m = 100.
km = 1.e5                       # cm
au = constants.au * 100.        # cm
pc = constants.parsec * 100.    # cm
kpc = 1000.*pc
kpc = 1.e3*pc
Mpc = 1.e6*pc
Gpc = 1.e9*pc

# angular

mas = pi / (180. * 3600.*1000.)
ns = 1.e-9
all_sky_sr = 4.*pi
all_sky_degsq = (180./pi)**2 * all_sky_sr

# frequencies

MHz = 1.e6
GHz = 1.e9

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
SMunit = kpc*10.**(-40./3.)             # kpc m^-20/3 to cgs
SM1 = 1.                                # kpc m^{-20/3}
SM_cgs = SM1 * kpc * cm_in_m**(-20./3.) # same as SMunit (a check)


# Flux density

Jy = 1.e-23     # cgs: erg/s/cm^2/Hz
Jansky = Jy

# Brightness temperature for 1 Jy pulse of 1 ms at 1 kpc at 1 GHz
Tb111 = (Jy * kpc**2) / (2.*kb*(GHz*ms)**2)

# plasma frequency in Hz for unit electron density
ne1 = 1.                                                # 1/cc
nup = sqrt((4.*pi * ne1 * e_e**2) / m_e) / (2.*pi)      # Hz

# cyclotron frequency in a 1 Gauss field
B1=1.
nuc = (e_e * B1) / (m_e * c) / (2.*pi)                  # Hz
