# mwprop.nemod v2.0 Jan 2026

'''
Configuration file for NE20x models
Reads input files and puts in dictionaries
Sets up spiral arms and also puts in dictionaries
'''

global rsun, rf_ref, vperp, sikol, louter, linner, ds_fine, ds_coarse
global Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, armmap, \
           r1, th1, th1deg, coarse_arms, rarray, tharray, \
           armsplines, armarray, tangents, normals, curvatures

global wg1, wg2, wga, wggc, wglism, wgcN, wgvN
global n1h1, h1, A1, F1, n2, h2, A2, F2, na, ha, wa, Aa, Fa
global apldr, bpldr, cpldr, dpldr, yldrmin, yldrmax
global aplsb, bplsb, cplsb, dplsb, ylsbmin, lsbmax
global alhb, blhb, clhb, xlhb, ylhb, zlhb, thetalhb, nelhb0, Flhb, ylhbmin, ylhbmax
global xlpI, ylpI, zlpI, rlpI, drlpI, nelpI, dnelpI, FlpI, dFlpI
global y_lism_min, y_lism_max
global xgc, ygc, zpgc, rgc, hgc, negc0, Fgc0
global eval_NE2025,eval_NE2001

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve, root, brentq
from scipy.optimize import minimize_scalar

from matplotlib.pyplot import figure, subplots_adjust, plot, axis
from matplotlib.pyplot import xlabel, ylabel, title, annotate
from matplotlib.pyplot import tick_params, legend, show, close, savefig

import os
import sys
import argparse
import importlib

script_path = os.path.dirname(os.path.realpath(__file__))

from pygedm.mwprop_core.get_constants_waveprop_scipy_version import *

# Get parameter input
from pygedm.mwprop_core.nemod import ne_input

# get some additional units from astropy
from astropy import units as u

# number of pc in kpc:
pc_in_kpc = (1* u.kpc / u.pc).cgs.value

sin = np.sin
cos = np.cos
exp = np.exp
sqrt = np.sqrt
zeros = np.zeros

deg2rad = np.deg2rad
rad2deg = np.rad2deg
sech2 = lambda z: 4.0 * np.exp(-2.0 * np.abs(z)) / (1.0 + np.exp(-2.0 * np.abs(z)))**2


# Solar system to Galactic center distance (kpc) set here
rsun = 8.5

# Reference radio frequency for calculating chromatic quantities
rf_ref = 1.             # GHz

# Default transverse velocity for calculating scintillation time scale
vperp = 100             # km/s

# Default spectral index for electron density wavenumber spectrum
sikol = 11/3                                            # spectral index

# Reference outer scale for Kolmogorov spectrum
louter = 1              # pc

# Reference inner scale for Kolmogorov spectrum
linner = 1000           # km

# coarse and fine sampling intervals for LoS numerical integrations
#ds_fine = 0.005 # SKO 1/9/2026 -- 0.0025 leads to spurious clump intersection for J1642, 0.01 made the discrepancy worse!
ds_fine = 0.0025 # SKO 3/2/22 -- changed from 0.01, reduces discrepancies from Fortran version
ds_coarse = 0.1

def setup_spiral_arms(Ncoarse=20, narmpoints=500, drfine=0.01):
    """
    Sets up spiral arm arrays to be done once per call to NE20x package.

    *** maybe rename to 'setup model' and put in higher-level routine

    Input:
        Ncoarse = number of coarse samples along spiral arms
        drfine = sample interval in radius for finely sampled arms
    """

    if 'Darms' not in dir():

        # Obtain model dictionaries and put spiral arm values in arrays

        #print('Creating model dictionaries')
        Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = \
            ne_input.read_nemod_parameters()

        aarm = Darms['a']
        rmin = Darms['rmin']
        thmin = Darms['thmin']
        extent = Darms['extent']
        narms = aarm.size

        # armmap maps integer order  to TC93 order,
        # which is from GC outwards toward Sun.

        armmap = np.array([1, 3, 4, 2, 5])
        Darmmap = {}
        for j in range(narms):
            Darmmap[str(j)] = armmap[j]

        """
        First evaluate spiral arms on a coarse grid of size Ncoarse
        This is used to find the nearest arm point on a coarse grid
        before calculating to higher precision
        """

        th1 = zeros((narms, Ncoarse))
        r1 = zeros((narms, Ncoarse))

        coarse_arms = zeros((2, narms, Ncoarse))

        extarray = np.linspace(zeros(5), extent, Ncoarse, endpoint=True).T  # rad
        th1 = np.array([thmin[j] + extarray[j] for j in range(narms)])         # rad
        th1deg = np.rad2deg(th1)                                            # deg
        r1 = np.array([rmin[j] * np.exp(extarray[j]/aarm[j]) for j in range(narms)])

        # Sculpt spiral arms as in NE2001 fortran code:

        ### NOTE ARM NUMBERING DIFFERENT FROM TC93 ARM NUMBERS
        ### See Fortran code

        # armindex = 1 ==> Arm 2 = TC arm 3
        inds2a = np.where((370 < th1deg[1]) & (th1deg[1] <= 410))
        arg2a = deg2rad((th1deg[1,inds2a]-390)*180/40.)
        r1[1, inds2a] *= (1 + 0.04*np.cos(arg2a))

        inds2b = np.where((315 < th1deg[1]) & (th1deg[1] <= 370))
        arg2b = deg2rad((th1deg[1,inds2b]-345)*180/55.)
        r1[1, inds2b] *= (1 - 0.07*np.cos(arg2b))

        inds2c = np.where((180 < th1deg[1]) & (th1deg[1] <= 315))
        arg2c = deg2rad((th1deg[1,inds2c]-260)*180/135.)
        r1[1, inds2c] *= (1 + 0.16*np.cos(arg2c))

        # armindex = 3 ==> Arm 4 = TC arm 2
        inds4a = np.where((290 < th1deg[3]) & (th1deg[3] <= 395))
        arg4a = deg2rad((th1deg[3,inds4a]-350)*180/105)
        r1[3, inds4a] *= (1 - 0.11*np.cos(arg4a))

        # r1 is now no longer monotonic splines need as independent variable

        coarse_arm_x = -r1 * np.sin(th1)
        coarse_arm_y =  r1 * np.cos(th1)

        coarse_arms[0]  = coarse_arm_x
        coarse_arms[1]  = coarse_arm_y

        # full spline for plotting outside routine

        tangents = zeros((2, narms, narmpoints))
        normals = zeros((2, narms, narmpoints))
        curvatures = zeros((narms, narmpoints))

        armarray = zeros((2, narms, narmpoints))
        armsplines = []

        rarray = zeros((narms, narmpoints))
        tharray = zeros((narms, narmpoints))
        for j in range(narms):
            sj = CubicSpline(th1[j,:], r1[j,:])
            rj = np.linspace(th1[j,0], th1[j,-1], num=narmpoints, endpoint=True)
            thj = np.linspace(th1[j,0], th1[j,-1], num=narmpoints, endpoint=True)
            rj = sj(thj)
            rarray[j] = rj
            tharray[j] = thj

            # Derivatives
            dthdr = 1 / sj(thj, 1)
            d2thdr2 = 1 / sj(thj, 2)

            dxdr = -np.sin(thj) - rj*np.cos(thj)*dthdr
            dydr =  np.cos(thj) - rj*np.sin(thj)*dthdr
            dydx = dydr / dxdr

            d2xdr2 = \
              -np.cos(thj)*dthdr + rj*np.sin(thj)*dthdr**2 - rj*np.cos(thj)*d2thdr2
            d2ydr2 = \
              -np.sin(thj)*dthdr - rj*np.cos(thj)*dthdr**2 - rj*np.sin(thj)*d2thdr2
            d2ydx2 = d2ydr2 / dxdr**2 + 1. / (d2xdr2 * dxdr * dydr)

            # unit tangent vector:
            facvec = np.array([1./sqrt(1.+dydx[k]**2) for k in range(np.size(dydx))])
            tx = facvec
            ty = facvec * dydx
            tangents[0, j] = tx
            tangents[1, j] = ty

            # unit normal vector:
            nx = -ty
            ny = tx


            normals[0, j] = nx
            normals[1, j] = ny

            Kr = (dxdr*d2ydr2 - dydr*d2xdr2) / (dxdr**2 + dydr**2)**1.5 # Curvature
            curvatures[j] = Kr

            armsplines = np.append(armsplines, sj)
            armarray[0, j, :] = -rj * np.sin(thj)
            armarray[1, j, :] =  rj * np.cos(thj)

    return Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, \
           armmap, \
           r1, th1, th1deg, coarse_arms, rarray, tharray, \
           armsplines, armarray, tangents, normals, curvatures, eval_NE2025, eval_NE2001


# General parameters for model

r"""
Values and coefficients for a 3D Kolmogorov spectrum

For SM modeling (as in NE2001):
Units conversion for SM in kpc m^{-20/3} and for outer scale in pc
i.e. SM \propto distance x l_o^{-2/3} x n_e^2
so with distance in kpc,  l_o in pc, and n_e in cm^{-3} we get
(kpc in cm) x (l_o in cm)^{-2/3} x (n_e in cm^{-3})^2 = 1.456e9 cm^{-17/3}

Then we divide by (kpc in cm) to get kpc cm^{-20/3}.
Finally converting to m^{-20/3} we have
c_u = (l_o in cm)^{-2/3} x (n_e in cm^{-3})^2 x (cm per m)^{20/3}
    = pc^{-2/3) x 10^{40/3}
    = (3.086e18)^{-2/3} x 10^{40/3}
    \simeq 10.165

For straight units conversion from kpc m^{-20/3} to cm^{-17/3}
we have
     SMunit (cgs)  = kpc*10.**(-40./3.)             # kpc m^-20/3 to cgs
defined in get_constants_waveprop_astropy.py (and scipy version)
We don't use this in the NE2001 modeling code..
"""
c_sm = (sikol-3) / (2 * (2*pi)**(4-sikol))              # ~ 0.181
c_u_units = ((1 * u.m)**(20/3) / (u.pc)**(2/3))         # using astropy units
c_u_astropy = c_u_units.cgs.value                       # evaluate astropy units
sm_factor = c_sm * c_u_astropy                          # net coefficient

dmax = 50                   # Maximum distance returned   (kpc)     need to rename in code
narmsmax = 5                # Number of spiral arms
narmsmax1 = narmsmax + 1

# New for NE2001p
dmax_ne2001p_integrate = 50 # maximum distance to integrate to for 'disk' components

"""
Ncoarse is the number of samples the log-spirals are initially evaluated
Ncoarse  = 20    # as in Fortran code; produces discontinuities
Ncoarse = 100    # reduces discontinuities, but no better than Ncoarse=50
"""

Ncoarse = 50     # reduces discontinuities

Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, \
    armmap, \
    r1, th1, th1deg, coarse_arms, rarray, tharray,  armsplines, armarray, \
    tangents, normals, curvatures, eval_NE2025, eval_NE2001 = setup_spiral_arms(Ncoarse=Ncoarse)

# ------------------------------
# Weights of density components:
# ------------------------------
wg1 = Dgal['wg1']
wg2 = Dgal['wg2']
wga = Dgal['wga']
wggc = Dgal['wggc']
wglism = Dgal['wglism']
wgcN = Dgal['wgcN']
wgvN = Dgal['wgvN']

# --------------------------
# Main component parameters:
# --------------------------
n1h1 = Dgal['n1h1']
h1 = Dgal['h1']
A1 = Dgal['A1']
F1 = Dgal['F1']
n2 = Dgal['n2']
h2 = Dgal['h2']
A2 = Dgal['A2']
F2 = Dgal['F2']
na = Dgal['na']
ha = Dgal['ha']
wa = Dgal['wa']
Aa = Dgal['Aa']
Fa = Dgal['Fa']


# ---------------------
# Local ISM Components:
# ---------------------
# Note theta values are converted to radians here

# ----
# LDR:
# ----
aldr,bldr,cldr = [Dlism[x] for x in ['aldr', 'bldr', 'cldr']]
xldr,yldr,zldr = [Dlism[x] for x in ['xldr','yldr','zldr']]
thetaldr,neldr0,Fldr = [Dlism[x] for x in ['thetaldr', 'neldr', 'Fldr']]
thetaldr = deg2rad(thetaldr)
sthldr = np.sin(thetaldr)
cthldr = np.cos(thetaldr)

# evaluate quantities used in neLDRQ1
apldr = (cthldr/aldr)**2 + (sthldr/bldr)**2
bpldr = (sthldr/aldr)**2 + (cthldr/bldr)**2
cpldr = 1./cldr**2
dpldr =  2.*cthldr*sthldr*(1./aldr**2 - 1./bldr**2)

yldrmin = yldr-max(aldr,bldr,cldr)
yldrmax = yldr+max(aldr,bldr,cldr)



# ----
# LSB:
# ----
alsb,blsb,clsb = [Dlism[x] for x in ['alsb','blsb','clsb']]
xlsb,ylsb,zlsb =  [Dlism[x] for x in ['xlsb', 'ylsb', 'zlsb']]
thetalsb,nelsb0,Flsb =  [Dlism[x] for x in ['thetalsb', 'nelsb', 'Flsb']]
thetalsb = deg2rad(thetalsb)
sthlsb = np.sin(thetalsb)
cthlsb = np.cos(thetalsb)

# evaluate quantities used in neLSB:
aplsb = (cthlsb/alsb)**2 + (sthlsb/blsb)**2
bplsb = (sthlsb/alsb)**2 + (cthlsb/blsb)**2
cplsb = 1./clsb**2
dplsb =  2.*cthlsb*sthlsb*(1./alsb**2 - 1./blsb**2)

ylsbmin = ylsb-max(alsb,blsb,clsb)
ylsbmax = ylsb+max(alsb,blsb,clsb)


# ----
# LHB:
# ----

alhb,blhb,clhb =  [Dlism[x] for x in ['alhb', 'blhb', 'clhb']]
xlhb,ylhb,zlhb =  [Dlism[x] for x in ['xlhb', 'ylhb', 'zlhb']]
thetalhb,nelhb0,Flhb = [Dlism[x] for x in ['thetalhb', 'nelhb', 'Flhb']]
thetalhb = deg2rad(thetalhb)

ylhbmin = ylhb-max(alhb,blhb,clhb)
ylhbmax = ylhb+max(alhb,blhb,clhb)


xlpI,ylpI,zlpI = [Dlism[x] for x in ['xlpI', 'ylpI', 'zlpI']]
rlpI,drlpI = [Dlism[x] for x in ['rlpI', 'drlpI']]
nelpI,dnelpI,FlpI,dFlpI =  [Dlism[x] for x in ['nelpI','dnelpI','FlpI','dFlpi']]

ylpImin = ylpI - rlpI - drlpI
ylpImax = ylpI + rlpI + drlpI

# y_lism_min = minimum y value for LISM components to matter
# y_lism_max = maximum y value for LISM components to matter

y_lism_min = min((yldrmin, ylsbmin, ylhbmin, ylpImin))
y_lism_max = max((yldrmax, ylsbmax, ylhbmax, ylpImax))



# --------------------------
# Galactic center component:
# --------------------------
xgc = Dgc['xgc']
ygc = Dgc['ygc']
zgc = Dgc['zgc']
rgc = Dgc['rgc']
hgc = Dgc['hgc']
negc0 = Dgc['negc0']
Fgc0 = Dgc['Fgc0']



# NOTE for now no global variables for clumps and voids


# -------
# Clumps:
# -------
nclumps = len(Dclumps)
#print('number of clumps:', nclumps)
lc = np.zeros(len(Dclumps))
bc = np.zeros(len(Dclumps))
nec = np.zeros(len(Dclumps))
Fc = np.zeros(len(Dclumps))
dc = np.zeros(len(Dclumps))
rc = np.zeros(len(Dclumps))
edgec = np.zeros(len(Dclumps))
for n,i in enumerate(Dclumps):
    lc[n] = Dclumps[i]['l']
    bc[n] = Dclumps[i]['b']
    nec[n] = Dclumps[i]['nec']
    Fc[n] = Dclumps[i]['Fc']
    dc[n] = Dclumps[i]['dc']
    rc[n] = Dclumps[i]['rc']
    edgec[n] = Dclumps[i]['edge']

slc = np.sin(deg2rad(lc))
clc = np.cos(deg2rad(lc))
sbc = np.sin(deg2rad(bc))
cbc = np.cos(deg2rad(bc))
rgalc = dc*cbc
xc = rgalc*slc
yc = rsun-rgalc*clc
zc = dc*sbc
# 11/23 SKO -- defining rcmult based on clump edge; this only gets used to filter the clumps, not the integration
rcmult = np.ones(np.shape(edgec))
rcmult[edgec==0] = 2.5 # roll-off
rcmult[edgec==1] = 1.1 # hard cutoff (only for a few clumps, e.g. Gum)

# ------
# Voids:
# ------
nvoids = len(Dvoids)
#print('number of voids:', nvoids)
lv = np.zeros(len(Dvoids))
bv = np.zeros(len(Dvoids))
dv = np.zeros(len(Dvoids))
nev = np.zeros(len(Dvoids))
Fv = np.zeros(len(Dvoids))
aav = np.zeros(len(Dvoids))
bbv = np.zeros(len(Dvoids))
ccv = np.zeros(len(Dvoids))
thvy = np.zeros(len(Dvoids))
thvz = np.zeros(len(Dvoids))
edgev = np.zeros(len(Dvoids))
for n,i in enumerate(Dvoids):
    lv[n] = Dvoids[i]['l']
    bv[n] = Dvoids[i]['b']
    dv[n] = Dvoids[i]['dv']
    nev[n] = Dvoids[i]['nev']
    Fv[n] = Dvoids[i]['Fv']
    aav[n] = Dvoids[i]['aav']
    bbv[n] = Dvoids[i]['bbv']
    ccv[n] = Dvoids[i]['ccv']
    thvy[n] = Dvoids[i]['thvy']
    thvz[n] = Dvoids[i]['thvz']
    edgev[n] = Dvoids[i]['edge']

slv = np.sin(deg2rad(lv))
clv = np.cos(deg2rad(lv))
sbv = np.sin(deg2rad(bv))
cbv = np.cos(deg2rad(bv))
rgalc = dv*cbv
xv = rgalc*slv
yv = rsun - rgalc*clv
zv = dv*sbv
# th1 = thvy
# th2 = thvz
s1 = np.sin(deg2rad(thvy))
c1 = np.cos(deg2rad(thvy))
s2 = np.sin(deg2rad(thvz))
c2 = np.cos(deg2rad(thvz))
cc12 = c1*c2
ss12 = s1*s2
cs21 = c2*s1
cs12 = c1*s2


def set_model(model: str):
    """Set NE2001/NE2025 model and refresh all dependent globals.

    Args:
        model: One of 'ne2001' or 'ne2025'.
    """
    global eval_NE2001, eval_NE2025
    global Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, armmap
    global r1, th1, th1deg, coarse_arms, rarray, tharray
    global armsplines, armarray, tangents, normals, curvatures
    global wg1, wg2, wga, wggc, wglism, wgcN, wgvN
    global n1h1, h1, A1, F1, n2, h2, A2, F2, na, ha, wa, Aa, Fa
    global apldr, bpldr, cpldr, dpldr, yldrmin, yldrmax
    global aplsb, bplsb, cplsb, dplsb, ylsbmin, ylsbmax
    global alhb, blhb, clhb, xlhb, ylhb, zlhb, thetalhb, nelhb0, Flhb, ylhbmin, ylhbmax
    global xlpI, ylpI, zlpI, rlpI, drlpI, nelpI, dnelpI, FlpI, dFlpI
    global y_lism_min, y_lism_max
    global xgc, ygc, zgc, rgc, hgc, negc0, Fgc0
    global nclumps, lc, bc, nec, Fc, dc, rc, edgec, slc, clc, sbc, cbc, xc, yc, zc, rcmult
    global nvoids, lv, bv, dv, nev, Fv, aav, bbv, ccv, thvy, thvz, edgev
    global slv, clv, sbv, cbv, xv, yv, zv, s1, c1, s2, c2, cc12, ss12, cs21, cs12

    model_key = model.lower().strip()
    if model_key not in ('ne2001', 'ne2025'):
        raise RuntimeError("Model must be ne2001 or ne2025.")

    if ((model_key == 'ne2001' and eval_NE2001) or (model_key == 'ne2025' and eval_NE2025)) and 'Darms' in globals():
        return

    # Write model selector used by ne_input.read_nemod_parameters()
    inpdir = os.path.dirname(os.path.realpath(__file__)) + '/params/'
    with open(inpdir + "which_model.inp", "w") as f:
        f.write("NE2001" if model_key == 'ne2001' else "NE2025")

    # Force reload of model parameters and spiral arm arrays
    if 'Darms' in globals():
        del Darms

    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, \
        armmap, r1, th1, th1deg, coarse_arms, rarray, tharray, armsplines, armarray, \
        tangents, normals, curvatures, eval_NE2025, eval_NE2001 = setup_spiral_arms(Ncoarse=Ncoarse)

    # Weights of density components
    wg1 = Dgal['wg1']
    wg2 = Dgal['wg2']
    wga = Dgal['wga']
    wggc = Dgal['wggc']
    wglism = Dgal['wglism']
    wgcN = Dgal['wgcN']
    wgvN = Dgal['wgvN']

    # Main component parameters
    n1h1 = Dgal['n1h1']
    h1 = Dgal['h1']
    A1 = Dgal['A1']
    F1 = Dgal['F1']
    n2 = Dgal['n2']
    h2 = Dgal['h2']
    A2 = Dgal['A2']
    F2 = Dgal['F2']
    na = Dgal['na']
    ha = Dgal['ha']
    wa = Dgal['wa']
    Aa = Dgal['Aa']
    Fa = Dgal['Fa']

    # Local ISM components
    aldr, bldr, cldr = [Dlism[x] for x in ['aldr', 'bldr', 'cldr']]
    xldr, yldr, zldr = [Dlism[x] for x in ['xldr', 'yldr', 'zldr']]
    thetaldr, neldr0, Fldr = [Dlism[x] for x in ['thetaldr', 'neldr', 'Fldr']]
    thetaldr = deg2rad(thetaldr)
    sthldr = np.sin(thetaldr)
    cthldr = np.cos(thetaldr)
    apldr = (cthldr/aldr)**2 + (sthldr/bldr)**2
    bpldr = (sthldr/aldr)**2 + (cthldr/bldr)**2
    cpldr = 1./cldr**2
    dpldr = 2.*cthldr*sthldr*(1./aldr**2 - 1./bldr**2)
    yldrmin = yldr - max(aldr, bldr, cldr)
    yldrmax = yldr + max(aldr, bldr, cldr)

    alsb, blsb, clsb = [Dlism[x] for x in ['alsb', 'blsb', 'clsb']]
    xlsb, ylsb, zlsb = [Dlism[x] for x in ['xlsb', 'ylsb', 'zlsb']]
    thetalsb, nelsb0, Flsb = [Dlism[x] for x in ['thetalsb', 'nelsb', 'Flsb']]
    thetalsb = deg2rad(thetalsb)
    sthlsb = np.sin(thetalsb)
    cthlsb = np.cos(thetalsb)
    aplsb = (cthlsb/alsb)**2 + (sthlsb/blsb)**2
    bplsb = (sthlsb/alsb)**2 + (cthlsb/blsb)**2
    cplsb = 1./clsb**2
    dplsb = 2.*cthlsb*sthlsb*(1./alsb**2 - 1./blsb**2)
    ylsbmin = ylsb - max(alsb, blsb, clsb)
    ylsbmax = ylsb + max(alsb, blsb, clsb)

    alhb, blhb, clhb = [Dlism[x] for x in ['alhb', 'blhb', 'clhb']]
    xlhb, ylhb, zlhb = [Dlism[x] for x in ['xlhb', 'ylhb', 'zlhb']]
    thetalhb, nelhb0, Flhb = [Dlism[x] for x in ['thetalhb', 'nelhb', 'Flhb']]
    thetalhb = deg2rad(thetalhb)
    ylhbmin = ylhb - max(alhb, blhb, clhb)
    ylhbmax = ylhb + max(alhb, blhb, clhb)

    xlpI, ylpI, zlpI = [Dlism[x] for x in ['xlpI', 'ylpI', 'zlpI']]
    rlpI, drlpI = [Dlism[x] for x in ['rlpI', 'drlpI']]
    nelpI, dnelpI, FlpI, dFlpI = [Dlism[x] for x in ['nelpI', 'dnelpI', 'FlpI', 'dFlpi']]
    ylpImin = ylpI - rlpI - drlpI
    ylpImax = ylpI + rlpI + drlpI
    y_lism_min = min((yldrmin, ylsbmin, ylhbmin, ylpImin))
    y_lism_max = max((yldrmax, ylsbmax, ylhbmax, ylpImax))

    # Galactic center component
    xgc = Dgc['xgc']
    ygc = Dgc['ygc']
    zgc = Dgc['zgc']
    rgc = Dgc['rgc']
    hgc = Dgc['hgc']
    negc0 = Dgc['negc0']
    Fgc0 = Dgc['Fgc0']

    # Clumps
    nclumps = len(Dclumps)
    lc = np.zeros(len(Dclumps))
    bc = np.zeros(len(Dclumps))
    nec = np.zeros(len(Dclumps))
    Fc = np.zeros(len(Dclumps))
    dc = np.zeros(len(Dclumps))
    rc = np.zeros(len(Dclumps))
    edgec = np.zeros(len(Dclumps))
    for n, i in enumerate(Dclumps):
        lc[n] = Dclumps[i]['l']
        bc[n] = Dclumps[i]['b']
        nec[n] = Dclumps[i]['nec']
        Fc[n] = Dclumps[i]['Fc']
        dc[n] = Dclumps[i]['dc']
        rc[n] = Dclumps[i]['rc']
        edgec[n] = Dclumps[i]['edge']

    slc = np.sin(deg2rad(lc))
    clc = np.cos(deg2rad(lc))
    sbc = np.sin(deg2rad(bc))
    cbc = np.cos(deg2rad(bc))
    rgalc = dc*cbc
    xc = rgalc*slc
    yc = rsun - rgalc*clc
    zc = dc*sbc
    rcmult = np.ones(np.shape(edgec))
    rcmult[edgec==0] = 2.5
    rcmult[edgec==1] = 1.1

    # Voids
    nvoids = len(Dvoids)
    lv = np.zeros(len(Dvoids))
    bv = np.zeros(len(Dvoids))
    dv = np.zeros(len(Dvoids))
    nev = np.zeros(len(Dvoids))
    Fv = np.zeros(len(Dvoids))
    aav = np.zeros(len(Dvoids))
    bbv = np.zeros(len(Dvoids))
    ccv = np.zeros(len(Dvoids))
    thvy = np.zeros(len(Dvoids))
    thvz = np.zeros(len(Dvoids))
    edgev = np.zeros(len(Dvoids))
    for n, i in enumerate(Dvoids):
        lv[n] = Dvoids[i]['l']
        bv[n] = Dvoids[i]['b']
        dv[n] = Dvoids[i]['dv']
        nev[n] = Dvoids[i]['nev']
        Fv[n] = Dvoids[i]['Fv']
        aav[n] = Dvoids[i]['aav']
        bbv[n] = Dvoids[i]['bbv']
        ccv[n] = Dvoids[i]['ccv']
        thvy[n] = Dvoids[i]['thvy']
        thvz[n] = Dvoids[i]['thvz']
        edgev[n] = Dvoids[i]['edge']

    slv = np.sin(deg2rad(lv))
    clv = np.cos(deg2rad(lv))
    sbv = np.sin(deg2rad(bv))
    cbv = np.cos(deg2rad(bv))
    rgalc = dv*cbv
    xv = rgalc*slv
    yv = rsun - rgalc*clv
    zv = dv*sbv
    s1 = np.sin(deg2rad(thvy))
    c1 = np.cos(deg2rad(thvy))
    s2 = np.sin(deg2rad(thvz))
    c2 = np.cos(deg2rad(thvz))
    cc12 = c1*c2
    ss12 = s1*s2
    cs21 = c2*s1
    cs12 = c1*s2

    _refresh_dependent_modules()


def _refresh_dependent_modules():
    """Reload modules that import config_nemod globals by value.

    This ensures model switches propagate to modules that used
    `from pygedm.mwprop_core.nemod.config_nemod import *` at import time.
    """
    module_names = [
        "pygedm.mwprop_core.nemod.density_components",
        "pygedm.mwprop_core.nemod.ne_arms",
        "pygedm.mwprop_core.nemod.ne_lism",
        "pygedm.mwprop_core.nemod.neclumpN_fast",
        "pygedm.mwprop_core.nemod.nevoidN",
        "pygedm.mwprop_core.nemod.density",
        "pygedm.mwprop_core.nemod.dmdsm",
    ]

    for name in module_names:
        module = sys.modules.get(name)
        if module is not None:
            importlib.reload(module)

    # Refresh the GalaxyModel singleton so pre-computed arm arrays / splines
    # stay in sync with the newly loaded model parameters.
    # Lazy import prevents a circular dependency at module-load time
    # (config_nemod must not import galaxy_model at the module level).
    _gm = sys.modules.get('pygedm.mwprop_core.nemod.galaxy_model')
    if _gm is not None:
        _gm.default_model.refresh()
