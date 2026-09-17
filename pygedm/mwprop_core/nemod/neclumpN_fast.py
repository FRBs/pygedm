# mwprop v2.0 Jan 2026

'''
clump subroutine -- put in density subfolder

comments copied from Fortran code:

returns electron density necN and fluctuation parameter FcN at position designated by l,b,d,x,y,z,c for a set of clumps with parameters read in from fild neclumpN.dat

input: x,y,z coordinates (kpc) (as in TC93)

output:
    necN            electron density in clump at (x,y,z)
    FcN             fluctuation parameter
    hitclump = 0:   no clump hit
            j>0:    jth clump hit

parameters:
    lc    = galactic longitude of clump center
    bc    = galactic latitude of clump center
    (xc,yc,zc) = clump center location (calculated)
    nec   = internal peak electron density
    rc    = clump radius at 1/e
    Fc    = clump fluctuation parameter
    edgec = 0 => use exponential rolloff out to 5rc
            1 => uniform and truncated at 1/e

lc,bc = Galactic coordinates (deg)
nec = clump electron density (cm^{-3})
Fc = fluctuation parameter
dc = clump distance from Earth (kpc)
rc = clump radius (kpc)
edgec = 0,1  0=> Gaussian, 1=> Gaussian w/ hard edge at e^{-1}
type = LOS type (P pulsar, G other Galactic, X extragalactic
losname = useful name

Version history:

01/18/2020 Stella Koch Ocker, initial conversion f77 --> python

01/23/20 -- JMC
    * now reads input parameters from dictionary program ne2001p_input

02/08/20 -- JMC
    * imports config_ne2001p for  model setup
    * rsun = 8.5 now commented out, included in setup file

12/28/21 - 01/02/22 -- JMC
    * added relevant_clumps function to prefilter clump list for the line
      of sight before integrating along it

11/27/2023 -- SKO
    * changed default value of rcmult to make sure clump prefiltering working properly
12/15/2023 -- SKO
    * corrected definition of inds_relevant so that LOS with dmax inside clump get counted
'''

from pygedm.mwprop_core.nemod.config_nemod import *

def relevant_clumps(l, b, dmax, rcmult):
    """
    Identifies clumps that contribute to n_e for the LoS
    specified by ldeg, bdeg, dmax.

    Input:
        l, b = Galactic coordinates of LoS (rad)
        dmax = maximum distance to calculate for LoS (kpc)
        rcmult = multiplier of clump radius for proximity decision

        Clump parameters in scope:
        lc, bc = Galactic coordinates of clumps (deg)
        xc, yc, zc = coordinates of clump (kpc)
        slc, clc, sbc, cbc = sin and cos of lc, bc
        rc = radial scales of clumps (kpc)
        dc = clump distances (kpc)
    """
    # calculate 'closest' = closest approach of LoS to clump:
    cl = cos(l)
    sl = sin(l)
    cb = cos(b)
    sb = sin(b)
    cos_theta = cb*cbc*(sl*slc+cl*clc) + sb*sbc
    # SKO -- need to force inds_relevant to recognize LOS that pass straight through clump center
    straight = np.where((1-cos_theta)<1e-4)
    cos_theta[straight] = 0.9998 # SKO -- to avoid sin inf errors
    sin_theta = sqrt(1-cos_theta**2)
    closest = dc * sin_theta
    closest[straight] = 0. # if you pass a LOS straight through the center, you want closest to be 0

    # want cos_theta > 0 because we want solutions where sin_theta is smaller than pi/2;
    # this will fail if there are ever very large clumps in the model.

    # note dc + rcmult*rc will be larger than dmax if dmax = dc
    # SKO 12-15-23: changed condition to dc - rcmult*rc < dmax

    #inds_relevant = \
    #   np.where((dc + rcmult*rc < dmax) & (cos_theta > 0) & (closest < rcmult*rc)) # old condition

    inds_relevant = \
        np.where((dc-rcmult*rc < dmax) & (cos_theta > 0) & (closest < rcmult*rc))

    #print('dc',dc[inds_relevant[0]])
    #print('dmax', dmax)
    #print('rcmult',rcmult[inds_relevant[0]])

    if np.any(inds_relevant)==False: # if no relevant clumps
        inds_relevant = np.array([-1])

    return inds_relevant

def neclumpN(x,y,z, inds_relevant=None):
    """
    Returns electron density and F parameter for position x,y,z contributed by
    any clumps along the LoS.

    inds_relevant:
        None:     step through all clumps in list initiated in config_ne2001p
        Not None: function assumes this is a tuple indicating indices for
                  only those clumps relevant for the line of sight.
    """

    if inds_relevant is None:
        clumpnums = range(0, nclumps)
    else:
        clumpnums = inds_relevant[0]

    necN = 0.
    hitclump = 0
    FcN = 0.
    arg = 0.
    hitclumpflag = np.zeros(nclumps)

    #print('number of relevant clumps:',clumpnums)

    if np.isscalar(clumpnums)==True and clumpnums==-1: # SKO 11/27/23 -- if no relevant clumps, exit with necN set to 0
        return necN,FcN,hitclump,arg

    else:
        for j in clumpnums:
            arg = ((x-xc[j])**2. + (y-yc[j])**2. + (z-zc[j])**2.) / (rc[j]**2.)
            if edgec[j] == 0. and arg < 5.:
                necN = necN + nec[j] * exp(-arg)
                FcN = Fc[j]
                hitclump = j
                hitclumpflag[j] = 1
            if edgec[j] == 1. and arg <= 1.:
                #print('arg: ',arg,' edgec: ',edgec[j])
                #print('xc,yc,zc,rc',xc,yc,zc,rc)
                necN = necN + nec[j]
                FcN = Fc[j]
                hitclump = j
                hitclumpflag[j] = 1


        return necN, FcN, hitclump, arg # SKO 3/6/22 -- added arg for debugging


# ---------------------------------------------------------------------------

def neclumpN_vec(x, y, z, inds_relevant=None):
    """
    Vectorized (array) version of neclumpN.

    x, y, z are 1-D numpy arrays of positions.
    Loops over the (few) relevant clumps; vectorizes over positions.

    Returns (necN_v, FcN_v, hitclump_v) — all 1-D arrays of length len(x).
    """
    n = len(x)
    necN_v     = np.zeros(n)
    FcN_v      = np.zeros(n)
    hitclump_v = np.zeros(n, dtype=int)

    if inds_relevant is None:
        clumpnums = range(nclumps)
    else:
        clumpnums = inds_relevant[0]

    if np.isscalar(clumpnums) and clumpnums == -1:
        return necN_v, FcN_v, hitclump_v

    for j in clumpnums:
        arg = ((x - xc[j])**2 + (y - yc[j])**2 + (z - zc[j])**2) / rc[j]**2
        if edgec[j] == 0.:              # Gaussian clump
            mask = arg < 5.
            necN_v    += np.where(mask, nec[j] * exp(-arg), 0.0)
            FcN_v      = np.where(mask, Fc[j],  FcN_v)
            hitclump_v = np.where(mask, j,       hitclump_v)
        elif edgec[j] == 1.:            # hard-edge clump
            mask = arg <= 1.
            necN_v    += np.where(mask, nec[j], 0.0)
            FcN_v      = np.where(mask, Fc[j],  FcN_v)
            hitclump_v = np.where(mask, j,       hitclump_v)

    return necN_v, FcN_v, hitclump_v
