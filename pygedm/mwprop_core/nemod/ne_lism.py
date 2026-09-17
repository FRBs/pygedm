# mwprop v2.0 Jan 2026

"""
ne_lism.py

Local interstellar medium functions

Replacements for Fortran functions in neLISM.NE2001.f:
    ne_LISM(x,y,z,FLISM,wLISM)      ! total local ISM
    neLDRQ1(x,y,z,FLDRQ1r,wLDRQ1)   ! Low Density Region in Q1
    neLSB(x,y,z,FLSBr,wLSB)         ! Local Super Bubble
    neLHB(x,y,z,FLHBr,wLHB)         ! Local Hot Bubble [not used]
    neLHB2(x,y,z,FLHBr,wLHB)        ! Local Hot Bubble
    neLOOPI(x,y,z,FLOOPI,wLOOPI)    ! Loop I

Python version JMC 2020 Jan 27

Changes:
2022 Jan 02:
    1. added if statement to skip calculations for x,y values
       outside of LISM
    2. put precalculations of some parameters into config_ne2001p.py
    3. annotations added

Comments from Fortran code neLISM.NE2001.f
c routines to calculate the electron density for the
c Local Interstellar Medium
c
c JMC 26 August-11 Sep. 2000
c     25 October 2001: modified to change weighting scheme
c                      so that the ranking is LHB: LSB: LDR
c                      (LHB overrides LSB and LDR; LSB overrides LDR)
c     16 November 2001:
c           added Loop I component with weighting scheme
c               LHB:LOOPI:LSB:LDR
c               LHB   overides everything,
c           LOOPI overrides LSB and LDR
c           LSB   overrides LDR
c           LISM  overrides general Galaxy
c     20 November 2001: The LOOPI component is truncated below z=0
c
c after discussions with Shami Chatterjee
c the sizes, locations and densities of the LISM components
c are based largely on work by Toscano et al. 1999 and Bhat et al. 1999

Note mwprop needs to be in PYTHONPATH

"""


from pygedm.mwprop_core.nemod.config_nemod import *

script_path = os.path.dirname(os.path.realpath(__file__))


# -----------------------------------------------------------------------------

def ne_lism(x,y,z):
    """
    Input variables:
        x,y,z = NE2001 Galactocentric coordinates (kpc)

    Output:
        ne_LISM = electron density
        FLISM = fluctuation parameter
        wLISM = weight for composite LISM
        wLDR = weight for LDR region
        wLHB = weight for local hot bubble
        wLSB = weight
        wLOOPI = weight for radio loop I

    Based on original Fortran routine.

    Changes:
        2022 Jan 02:  put in test for y out of range of LISM; return zeros if so
                      note: y_lism_min, y_lism_max computed in config_ne2001p.py
    """

    # Test whether input y is enclosed by LISM; if not set values to zero and return

    if y > y_lism_max or y < y_lism_min:
        return 0, 0, 0, 0, 0, 0, 0

    neldrq1xyz, FLDRQ1r, wLDR = neLDRQ1(x,y,z)    # low density region in Q1
    nelsbxyz, FLSBr, wLSB = neLSB(x,y,z)          # Local Super Bubble
    nelhbxyz, FLHBr, wLHB = neLHB2(x,y,z)         # Local Hot Bubble
    neloopIxyz, FLOOPIr, wLOOPI = neLOOPI(x,y,z)  # Loop I

    """
    weight the terms so that the LHB term overrides the other
    terms (we want the density to be low in the LHB, lower than
    in the other terms.
    """

    ne_LISM =   ((1-wLHB)   *
                 (
                   (1-wLOOPI) * (wLSB*nelsbxyz + (1-wLSB)*neldrq1xyz)
               +     wLOOPI * neloopIxyz
                 )
               +     wLHB  * nelhbxyz)

    FLISM = ((1-wLHB) *
                (
                   (1-wLOOPI) * (wLSB*FLSBr + (1-wLSB)*FLDRQ1r)
               +     wLOOPI * FLOOPIr
                )
               +     wLHB  * FLHBr)

    """
    return the maximum weight of any of the terms for
    combining with additional terms external to this routine.
    """

    wLISM = max(wLOOPI, max(wLDR, max(wLSB, wLHB)))

    return ne_LISM, FLISM, wLISM, wLDR, wLHB, wLSB, wLOOPI
    return

# -----------------------------------------------------------------------------

def neLDRQ1(x,y,z):    # Low Density Region in Q1
    """
    Low Density Region in Q1
    input:
      x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00

    output:
      neLDRQ1 = electron density in local hot bubble that
              is modeled as an ellipsoidal trough.
      FLDRQ1 = fluctuation parameter
      wLDRQ1  = weight of LDRQ1 component used to combine
          with other components of electron density.
          wLDRQ1 =  1  at and inside the annular ridge
                 <  1  outside the annular ridge
                 -> 0  far outside the annular ridge
      e.g. total electron density would be evaluated as
               ne = (1-wLDRQ1)*ne_other + neLDRQ1
      note thetaldr is converted to radians in config_ne2001p.py

    These calculations now done in config_ne2001p.py
    aa=aldr
    bb=bldr
    cc=cldr
    theta=thetaldr
    netrough =neldr0
    Ftrough=Fldr

    # can redef s,c as sldr, cldr and do these calcs outside function in setup
    s = np.sin(theta)
    c = np.cos(theta)

    ap = (c/aa)**2 + (s/bb)**2
    bp = (s/aa)**2 + (c/bb)**2
    cp = 1./cc**2
    dp =  2.*c*s*(1./aa**2 - 1./bb**2)
    sldr = np.sin(thetaldr)
    cldr = np.cos(thetaldr)

    apldr = (cldr/aldr)**2 + (sldr/bldr)**2
    bpldr = (sldr/aldr)**2 + (cldr/bldr)**2
    cpldr = 1./cldr**2
    dpldr =  2.*cldr*sldr*(1./aldr**2 - 1./bldr**2)
    """

    neLDRQ1 = 0.
    wLDRQ1 = 0
    FLDRQ1r = 0.
    q = ((x-xldr)**2*apldr
          + (y-yldr)**2*bpldr
          + (z-zldr)**2*cpldr
          + (x-xldr)*(y-yldr)*dpldr)
    if q <= 1:                     # inside
        neLDRQ1 = neldr0
        FLDRQ1r = Fldr
        wLDRQ1 = 1

    return neLDRQ1, FLDRQ1r, wLDRQ1

# -----------------------------------------------------------------------------

def neLSB(x,y,z):                       # Local Super Bubble
    """
    Local Super Bubble

    input:
      x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
    output:
      neLSB = electron density in local hot bubble that
              is modeled as an ellisoidal trough.
      FLSB = fluctuation parameter
      wLSB  = weight of LSB component used to combine
          with other components of electron density.
          wLSB =  1  at and inside the annular ridge
               <  1  outside the annular ridge
                   -> 0  far outside the annular ridge
      e.g. total electron density would be evaluated as
      ne = (1-wLSB)*ne_other + neLSB

    These now calculated in input_ne2001p.py

    aa=alsb
    bb=blsb
    cc=clsb
    theta=thetalsb
    netrough=nelsb0
    Ftrough=Flsb

    s = np.sin(theta)
    c = np.cos(theta)
    ap = (c/aa)**2 + (s/bb)**2
    bp = (s/aa)**2 + (c/bb)**2
    cp = 1./cc**2
    dp =  2.*c*s*(1./aa**2 - 1./bb**2)
    """

    neLSB = 0.
    wLSB = 0
    FLSBr = 0.
    q = ((x-xlsb)**2*aplsb
          + (y-ylsb)**2*bplsb
          + (z-zlsb)**2*cplsb
          + (x-xlsb)*(y-ylsb)*dplsb)
    if q <= 1:                      # inside region
        neLSB = nelsb0
        FLSBr = Flsb
        wLSB = 1

    return neLSB, FLSBr, wLSB

# -----------------------------------------------------------------------------

def neLHB(x,y,z):                 # Local Hot Bubble
    """
    Local Hot Bubble
    input:
      x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
    output:
      neLHB = electron density in local hot bubble that
              is modeled as an ellisoidal trough.
      FLHB = fluctuation parameter
      wLHB  = weight of LBH component used to combine
          with other components of electron density.
          wLBH =  1  at and inside the annular ridge
               <  1  outside the annular ridge
                   -> 0  far outside the annular ridge
      e.g. total electron density would be evaluated as
               ne = (1-wLHB)*ne_other + neLHB

    These are now calculated in input_ne2001p.py
    aa=alhb
    bb=blhb
    cc=clhb
    theta=thetalhb
    netrough=nelhb0
    Ftrough=Flhb

    s = np.sin(theta)
    c = np.cos(theta)
    ap = (c/aa)**2 + (s/bb)**2
    bp = (s/aa)**2 + (c/bb)**2
    cp = 1./cc**2
    dp = 2.*c*s*(1./aa**2 - 1./bb**2)
    """

    neLHB = 0.
    wLHB = 0
    FLHBr = 0.
    q = ((x-xlhb)**2*aplhb
      + (y-ylhb)**2*bplhb
      + (z-zlhb)**2*cplhb
      + (x-xlhb)*(y-ylhb)*dplhb)
    if q <= 1:                           # inside
        neLHB = nelhb0
        FLHBr = Flhb
        wLHB = 1

    return neLHB, FLHBr, wLHB

# -----------------------------------------------------------------------------

def neLHB2(x,y,z):                        # Local Hot Bubble
    """
    LHB modeled as a cylinder
    the cylinder slants in the y direction vs. z as described by parameter yzslope
    the cylinder cross-sectional size in the 'a' direction (major axis)
          varies with z, tending to zero at its smallest z point.
        implicit none
        real x,y,z,FLHBr
        integer wLHB

    input:
      x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00

    output:
      neLHB2 = electron density in local hot bubble that
              is modeled as an ellisoidal trough.
      FLHB = fluctuation parameter
      wLHB  = weight of LBH component used to combine
          with other components of electron density.
          wLHB =  1  at and inside the annular ridge
               <  1  outside the annular ridge
                   -> 0  far outside the annular ridge
      e.g. total electron density would be evaluated as
               ne = (1-wLHB)*ne_other + neLHB2
    aa=alhb
    bb=blhb
    cc=clhb
    netrough=nelhb0
    Ftrough=Flhb
    theta = thetalhb
    """

    neLHB2 = 0.
    wLHB = 0
    FLHBr = 0.

    yzslope = np.tan(thetalhb)
    yaxis = ylhb + yzslope*z

    """
    cylinder has cross sectional area = constant for z > 0
    area -> 0 for z < 0 by letting aa->0 linearly for z < 0:
    (0.001 = 1 pc is to avoid divide by zero)
    """

    if z <= 0 and z >= zlhb-clhb:
      aa = 0.001 + (alhb-0.001)*(1. - (1./(zlhb-clhb))*z)
    else:
      aa = alhb
    qxy =  ((x-xlhb)/aa)**2 + ((y-yaxis)/blhb)**2
    qz =  abs(z-zlhb)/clhb
    if qxy <= 1 and qz <= 1:            # inside
        neLHB2 = nelhb0
        FLHBr = Flhb
        wLHB = 1

    return neLHB2, FLHBr, wLHB

# -----------------------------------------------------------------------------

def neLOOPI(x,y,z):                      # Loop I
    """
    Loop I

    Component is a spheroid truncated for z<0.

    input:
      x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00

    output:
      neLOOPI = electron density in LOOP I that
              is modeled as an ellisoidal trough
          with an enhanced shell
      FLOOPI = fluctuation parameter
      wLOOPI  = weight of LOOP I component used to combine
          with other components of electron density.
          wLOOPI =  1  at and inside the annular ridge
                 <  1  outside the annular ridge
    """
    a1 = rlpI
    a2 = rlpI+drlpI

    if z < 0:
        neLOOPI = 0.
        FLOOPI = 0.
        wLOOPI = 0
        return neLOOPI, FLOOPI, wLOOPI

    r = sqrt( (x-xlpI)**2 + (y-ylpI)**2 + (z-zlpI)**2)
    if r > a2:                  # outside Loop I
        neLOOPI = 0.
        FLOOPI = 0.
        wLOOPI = 0
    elif r <= a1:               # inside volume
        neLOOPI= nelpI
        FLOOPI = FlpI
        wLOOPI = 1
    else:                       # inside boundary shell
        neLOOPI= dnelpI
        FLOOPI = dFlpI
        wLOOPI = 1

    return neLOOPI, FLOOPI, wLOOPI

# -----------------------------------------------------------------------------
# ============================================================================
# Vectorized (array) versions — accept 1-D numpy arrays x, y, z
# ============================================================================

def neLDRQ1_vec(x, y, z):
    """Vectorized neLDRQ1: ellipsoidal LDR region."""
    q = ((x - xldr)**2 * apldr
         + (y - yldr)**2 * bpldr
         + (z - zldr)**2 * cpldr
         + (x - xldr) * (y - yldr) * dpldr)
    mask = q <= 1
    return np.where(mask, neldr0, 0.0), np.where(mask, Fldr, 0.0), mask.astype(float)


def neLSB_vec(x, y, z):
    """Vectorized neLSB: ellipsoidal Local Super Bubble."""
    q = ((x - xlsb)**2 * aplsb
         + (y - ylsb)**2 * bplsb
         + (z - zlsb)**2 * cplsb
         + (x - xlsb) * (y - ylsb) * dplsb)
    mask = q <= 1
    return np.where(mask, nelsb0, 0.0), np.where(mask, Flsb, 0.0), mask.astype(float)


def neLHB2_vec(x, y, z):
    """Vectorized neLHB2: cylindrical Local Hot Bubble with z-varying radius."""
    yaxis = ylhb + np.tan(thetalhb) * z
    # cylinder semi-axis aa shrinks linearly for z in [zlhb-clhb, 0]
    denom = zlhb - clhb          # scalar, = -0.16 for default params (nonzero)
    lhb2_cond = (z <= 0) & (z >= denom)
    aa = np.where(lhb2_cond, 0.001 + (alhb - 0.001) * (1.0 - z / denom), alhb)
    qxy = ((x - xlhb) / aa)**2 + ((y - yaxis) / blhb)**2
    qz  = np.abs(z - zlhb) / clhb
    mask = (qxy <= 1) & (qz <= 1)
    return np.where(mask, nelhb0, 0.0), np.where(mask, Flhb, 0.0), mask.astype(float)


def neLOOPI_vec(x, y, z):
    """Vectorized neLOOPI: spheroidal Loop I (truncated for z < 0)."""
    r = np.sqrt((x - xlpI)**2 + (y - ylpI)**2 + (z - zlpI)**2)
    a2 = rlpI + drlpI
    active  = (z >= 0) & (r <= a2)      # not outside
    inside  = active   & (r <= rlpI)    # core
    shell   = active   & ~inside        # boundary shell
    neLOOPI_v = np.where(inside, nelpI,  np.where(shell, dnelpI, 0.0))
    FLOOPI_v  = np.where(inside, FlpI,   np.where(shell, dFlpI,  0.0))
    return neLOOPI_v, FLOOPI_v, active.astype(float)


def ne_lism_vec(x, y, z):
    """
    Vectorized ne_lism: combines all four LISM sub-components with the same
    hierarchical weighting as the scalar ne_lism().

    Returns (ne_LISM, FLISM, wLISM, wLDR, wLHB, wLSB, wLOOPI) — all arrays.
    """
    out_of_range = (y > y_lism_max) | (y < y_lism_min)

    neLDR_v,   FLDR_v,   wLDR_v   = neLDRQ1_vec(x, y, z)
    neLSB_v,   FLSB_v,   wLSB_v   = neLSB_vec(x, y, z)
    neLHB2_v,  FLHB_v,   wLHB_v   = neLHB2_vec(x, y, z)
    neLOOPI_v, FLOOPI_v, wLOOPI_v = neLOOPI_vec(x, y, z)

    ne_LISM_v = ((1 - wLHB_v) *
                 ((1 - wLOOPI_v) * (wLSB_v * neLSB_v + (1 - wLSB_v) * neLDR_v)
                  + wLOOPI_v * neLOOPI_v)
                 + wLHB_v * neLHB2_v)

    FLISM_v = ((1 - wLHB_v) *
               ((1 - wLOOPI_v) * (wLSB_v * FLSB_v + (1 - wLSB_v) * FLDR_v)
                + wLOOPI_v * FLOOPI_v)
               + wLHB_v * FLHB_v)

    wLISM_v = np.maximum(wLOOPI_v, np.maximum(wLDR_v, np.maximum(wLSB_v, wLHB_v)))

    # Zero out points outside the LISM bounding box
    z0 = np.zeros_like(ne_LISM_v)
    ne_LISM_v = np.where(out_of_range, z0, ne_LISM_v)
    FLISM_v   = np.where(out_of_range, z0, FLISM_v)
    wLISM_v   = np.where(out_of_range, z0, wLISM_v)
    wLDR_v    = np.where(out_of_range, z0, wLDR_v)
    wLHB_v    = np.where(out_of_range, z0, wLHB_v)
    wLSB_v    = np.where(out_of_range, z0, wLSB_v)
    wLOOPI_v  = np.where(out_of_range, z0, wLOOPI_v)

    return ne_LISM_v, FLISM_v, wLISM_v, wLDR_v, wLHB_v, wLSB_v, wLOOPI_v

# -----------------------------------------------------------------------------