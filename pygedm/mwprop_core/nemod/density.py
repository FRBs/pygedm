# mwprop v2.0 Jan 2026

"""
density.py

Returns the local electron density, fluctuation parameters, weights
and clump/void flags for the position (x, y, z)

Relacement for Fortran subroutine density_2001


Comments from Fortran code density.NE2001.f
c final version of NE2001
c returns densities, F parameters and weights of the various components
c mods:
c 28 July 02:
c   put in 'if' statements in subroutine density_2001 so that
c   function calls are not done if particular weights
c   (wg1, wg2, etc.) are set to zero in gal01.inp
c   This is (a) cleaner and (b) much more efficient if
c   the clump or void component is flagged out.

Python version JMC 2020 Jan 28

Changes:
    2022 Jan 01: added functions density_2001_smallscale_comps and density_2001_smooth_comps
                 so that different LoS grids can be used to integrate different contributions
                 (to speed up calculations)
"""

from pygedm.mwprop_core.nemod.config_nemod import *

from pygedm.mwprop_core.nemod.density_components import *
from pygedm.mwprop_core.nemod.ne_lism import *
from pygedm.mwprop_core.nemod.neclumpN_fast import *
from pygedm.mwprop_core.nemod.nevoidN import *
from pygedm.mwprop_core.nemod.ne_arms import *

script_path = os.path.dirname(os.path.realpath(__file__))

# ---------------------------------------------------------------------------

def density_2001(x, y, z, inds_relevant=None, verbose=False):
    """
    (Edited) Comments from Fortran code:

     Returns seven components of the free electron density of the
     interstellar medium at Galactic location (x,y,z).

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      ne1:    outer, thick disk
      ne2:    inner, thin disk (annular in form)
      nea:    spiral arms
      negc:   galactic center component
          nelism: local ISM component
          necN:   contribution from discrete 'clumps'
          nevN:   contribution from voids
       fluctuation parameters (one for each ne component):
          F1, F2, Fa, Fgc, Flism, FcN, FvN
       flags:
          whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
             wlism: 1 if x,y,z is in any of the four LISM components
              wLDR: 1 if in LDR, 0 if not
              wLHB: 1 if in LHB, 0 if not
              wLSB: 1 if in LSB, 0 if not
            wLOOPI: 1 if in LoopI, 0 if not
          (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
          hitclump: clump number that x,y,z is in (0 if none)
           hitvoid: void number that x,y,z is in (0 if none)
    25 May 2002
    based on routines from TC93 and test routines from 1999-2002 by JMC.
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    ne1 = ne2 = nea = negc = nelism = necN = nevN = 0
    wlism=wldr=wlhb=wlsb=wloopI = hitclump = hitvoid = wvoid = whicharm = 0


    if wg1 == 1:
        ne1, F1 = ne_outer(x,y,z)
    if wg2 == 1:
        ne2, F2 = ne_inner(x,y,z)
    if wga == 1:
        nea, Fa, whicharm = ne_arms_ne2001p(x,y,z, Ncoarse=Ncoarse)
    else:
        nea = Fa = 0.
    if wggc == 1:
        negc, Fgc = ne_gc(x,y,z)
    if wglism == 1:
        nelism, Flism, wlism, wldr, wlhb, wlsb, wloopI = ne_lism(x,y,z)
    if wgcN == 1:
        necN, FcN, hitclump, arg = neclumpN(x,y,z, inds_relevant=inds_relevant)
    if wgvN == 1:
        nevN, FvN, hitvoid, wvoid = nevoidN(x,y,z)

    if verbose:
        print('density: ', ne1, ne2, F1, F2)


    return ne1,ne2,nea,negc,nelism,necN,nevN, \
           F1, F2, Fa, Fgc, Flism, FcN, FvN, \
           whicharm, wlism, wldr, wlhb, wlsb, wloopI, \
           hitclump, hitvoid, wvoid


# ---------------------------------------------------------------------------

def density_2001_smallscale_comps(x, y, z, inds_relevant=None, verbose=False):
    """
    (Edited) Comments from Fortran code:

     Returns four components of the free electron density of the
     interstellar medium at Galactic location (x,y,z).

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      negc:   galactic center component
      nelism: local ISM component
      necN:   contribution from discrete 'clumps'
      nevN:   contribution from voids
      fluctuation parameters (one for each ne component):
          Fgc, Flism, FcN, FvN
      flags:
             wlism: 1 if x,y,z is in any of the four LISM components
              wLDR: 1 if in LDR, 0 if not
              wLHB: 1 if in LHB, 0 if not
              wLSB: 1 if in LSB, 0 if not
            wLOOPI: 1 if in LoopI, 0 if not
          (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
          hitclump: clump number that x,y,z is in (0 if none)
           hitvoid: void number that x,y,z is in (0 if none)
    01 Jan 2022  based on fortran core from 25 May 2002
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    negc = nelism = necN = nevN = 0
    wlism=wldr=wlhb=wlsb=wloopI = hitclump = hitvoid = wvoid = 0

    if wggc == 1:
        negc, Fgc = ne_gc(x,y,z)
    if wglism == 1:
        nelism, Flism, wlism, wldr, wlhb, wlsb, wloopI = ne_lism(x,y,z)
    if wgcN == 1:
        necN, FcN, hitclump, arg = neclumpN(x,y,z, inds_relevant=inds_relevant) # SKO added arg to output for debugging
    if wgvN == 1:
        nevN, FvN, hitvoid, wvoid = nevoidN(x,y,z)

    #print('wgvN',wgvN,'hitvoid',hitvoid)
    #if verbose:
    #print('x,y,z,',x,y,z)
    #print('density: ', negc, nelism, necN, nevN)

    if wgcN == 0:
        necN = 0
        FcN = 0
        hitclump = 0

    return negc,nelism,necN,nevN, \
           Fgc, Flism, FcN, FvN, \
           wlism, wldr, wlhb, wlsb, wloopI, \
           hitclump, hitvoid, wvoid

# ---------------------------------------------------------------------------

def density_2001_smooth_comps(x, y, z, verbose=False):
    """
    Returns only the electron density for the smooth components.
    Utility: can be coarsely sampled to speed up computations.

    Output density is the sum ne1 + ne2 + nea at (x,y,z)

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      ne1:    outer, thick disk
      ne2:    inner, thin disk (annular in form)
      nea:    spiral arms
    01 Jan 2022
    based on routines from NE2001
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    ne1 = ne2 = nea = 0
    F1 = F2 = Fa = 0
    whicharm = 0

    if wg1 == 1:
        ne1, F1 = ne_outer(x,y,z)
    if wg2 == 1:
        ne2, F2 = ne_inner(x,y,z)
    if wga == 1:
        nea, Fa, whicharm = ne_arms_ne2001p(x,y,z, Ncoarse=Ncoarse)

    if verbose:
        print('density: ', ne1, ne2, F1, F2)

    # Define smooth components
    # Fsmooth defined so that
    # Fsmooth*ne_smooth**2 = F1*ne1**2 + F2*ne2**2 + Fa*nea**2`
    wne1 = wg1*ne1
    wne2 = wg2*ne2
    wnea = wga*nea

    ne_smooth = wne1 + wne2 + wnea

    if ne_smooth > 0.:
        Fsmooth = (F1*(wne1)**2 + F2*wne2**2 + Fa*wnea**2) / ne_smooth**2
    else:
        Fsmooth = 0.

    return ne1,ne2,nea, F1, F2, Fa, whicharm, ne_smooth, Fsmooth

    """
    if verbose:
        return ne1,ne2,nea, \
               F1, F2, Fa, \
               whicharm, ne_smooth
    else:
        return ne_smooth
    """


# ---------------------------------------------------------------------------

def density_2001_smallscale_comps_vec(xvec, yvec, zvec, inds_relevant=None):
    """
    Vectorized counterpart of density_2001_smallscale_comps.

    Accepts 1-D numpy arrays xvec, yvec, zvec (length Ns_fine) and returns
    all small-scale component arrays in the same order as the scalar version,
    so that the caller can unpack with:

        negc, nelism, necN, nevN, Fgc, Flism, FcN, FvN, \\
            wlism, wldr, wlhb, wlsb, wloopI, hitclump, hitvoid, wvoid = \\
            density_2001_smallscale_comps_vec(xf_vec, yf_vec, zf_vec,
                                              inds_relevant=...)
    """
    n = len(xvec)

    negc_v    = np.zeros(n)
    nelism_v  = np.zeros(n)
    necN_v    = np.zeros(n)
    nevN_v    = np.zeros(n)
    Fgc_v     = np.zeros(n)
    Flism_v   = np.zeros(n)
    FcN_v     = np.zeros(n)
    FvN_v     = np.zeros(n)
    wlism_v   = np.zeros(n)
    wldr_v    = np.zeros(n)
    wlhb_v    = np.zeros(n)
    wlsb_v    = np.zeros(n)
    wloopI_v  = np.zeros(n)
    hitclump_v = np.zeros(n, dtype=int)
    hitvoid_v  = np.zeros(n, dtype=int)
    wvoid_v    = np.zeros(n)

    if wggc == 1:
        negc_v, Fgc_v = ne_gc_vec(xvec, yvec, zvec)

    if wglism == 1:
        nelism_v, Flism_v, wlism_v, wldr_v, wlhb_v, wlsb_v, wloopI_v = \
            ne_lism_vec(xvec, yvec, zvec)

    if wgcN == 1:
        necN_v, FcN_v, hitclump_v = neclumpN_vec(xvec, yvec, zvec,
                                                  inds_relevant=inds_relevant)

    if wgvN == 1:
        nevN_v, FvN_v, hitvoid_v, wvoid_v = nevoidN_vec(xvec, yvec, zvec)

    return (negc_v, nelism_v, necN_v, nevN_v,
            Fgc_v, Flism_v, FcN_v, FvN_v,
            wlism_v, wldr_v, wlhb_v, wlsb_v, wloopI_v,
            hitclump_v, hitvoid_v, wvoid_v)


# ---------------------------------------------------------------------------

def density_2001_smooth_comps_vec(xvec, yvec, zvec):
    """
    Vectorized counterpart of density_2001_smooth_comps.

    Accepts 1-D numpy arrays xvec, yvec, zvec (length Ns_coarse) and returns
    all smooth-component arrays in the same order as the scalar version so the
    caller can unpack with::

        cne1, cne2, cnea, cF1, cF2, cFa, cwhicharm, cne_smooth, cFsmooth = \\
            density_2001_smooth_comps_vec(xc_vec, yc_vec, zc_vec)
    """
    n = len(xvec)

    ne1_v      = np.zeros(n)
    ne2_v      = np.zeros(n)
    nea_v      = np.zeros(n)
    F1_v       = np.zeros(n)
    F2_v       = np.zeros(n)
    Fa_v       = np.zeros(n)
    whicharm_v = np.zeros(n, dtype=int)

    if wg1 == 1:
        ne1_v, F1_v = ne_outer_vec(xvec, yvec, zvec)
    if wg2 == 1:
        ne2_v, F2_v = ne_inner_vec(xvec, yvec, zvec)
    if wga == 1:
        nea_v, Fa_v, whicharm_v = ne_arms_ne2001p_vec(xvec, yvec, zvec,
                                                       Ncoarse=Ncoarse)

    wne1 = wg1 * ne1_v
    wne2 = wg2 * ne2_v
    wnea = wga * nea_v
    ne_smooth_v = wne1 + wne2 + wnea

    # Fsmooth = (F1*ne1^2 + F2*ne2^2 + Fa*nea^2) / ne_smooth^2  where ne_smooth > 0
    Fsmooth_v = np.where(
        ne_smooth_v > 0.,
        (F1_v * wne1**2 + F2_v * wne2**2 + Fa_v * wnea**2) / ne_smooth_v**2,
        0.)

    return ne1_v, ne2_v, nea_v, F1_v, F2_v, Fa_v, whicharm_v, ne_smooth_v, Fsmooth_v
