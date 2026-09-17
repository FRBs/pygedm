# mwprop.nemod v2.0 Jan 2026

"""
ne_arms
Replacement for Fortran routine NE_ARMS_LOG_MOD in density.NE2001.f

Version history:
2020 Jan 19
01/19/20 --- JMC
     * initial conversion f77 --> python
02/08/20 --- JMC
     * removed rsun = 8.5
     * all model parameters now imported from config_ne2001p
2026 Feb --- DCP
     * Modified to use GalaxyModel class with pre-instantiated attributes
"""

from pygedm.mwprop_core.nemod.config_nemod import *
from pygedm.mwprop_core.nemod.galaxy_model import default_model

script_path = os.path.dirname(os.path.realpath(__file__))

def ne_arms_ne2001p(x, y, z, Ncoarse=20, dthfine=0.01, nfinespline=5,
                    verbose=False, model=None):
    """
    Evaluates electron density and F parameter for arm nearest to x,y,z

    Input:
        x, y, z = location in NE2001  coordinates  in kpc
        Ncoarse = number of coarse samples along each arm
        dthfine = step size in Galactocentric angle for fine sampling (rad)
        nfinespline = number of coarse samples to use for fine-sampling spline
        verbose = True:  writes out n_e component information
        model   = GalaxyModel instance (uses module-level default_model if None)

    Output:
        ne        electron density        cm^{-3}
        F         F parameter             composite units
        whicharm  which spiral arms       1 to 5, no close arm ==> 0
    """
    if model is None:
        model = default_model

    # First find coarse location of points on arms nearest to input point

    narms = model.narms
    dsq_coarse = (coarse_arms[0] - x)**2 + (coarse_arms[1] - y)**2
    index_dsqmin = np.argmin(dsq_coarse, axis=1)

    # Zoom in to find precision location and evaluate density and F parameter

    # number of coarse samples to use as input to fine spline
    nspline2 = int((nfinespline-1)/2)

    # Initialize
    nea = 0.
    ga = 0.
    whicharm = 0
    whicharm_spiralmodel = 0

    # Cylindrical radius and Galactocentric angle of input point;
    # thxydeg measured CCW from +y axis

    rr = sqrt(x**2 + y**2)
    thxydeg = np.rad2deg(np.arctan2(-x, y))
    if thxydeg < 0:
        thxydeg += 360

    # Use pre-computed model attributes instead of per-call dict lookups:
    arm_index = model.arm_index
    warm      = model.warm
    harm      = model.harm
    narm      = model.narm

    if verbose:
        print('arms rr, thxydeg: ', rr, thxydeg)

    for j in range(narms):

        # Use fine spline on coarse samples  to find nearest position on j-th arm

        ind1 = range(max(0, index_dsqmin[j]-nspline2-1),
                     min(index_dsqmin[j]+nspline2+1, Ncoarse), 1)
        dsqs = CubicSpline(th1[j, ind1], dsq_coarse[j, ind1])
        thjfine = np.arange(th1[j, ind1[0]], th1[j, ind1[-1]], dthfine)
        ind_min = dsqs(thjfine).argmin()
        thjmin = thjfine[ind_min]
        # Use pre-computed arm radius spline from model — no per-call CubicSpline
        rjmin = model.arm_radius_splines[j](thjmin)
        xjmin, yjmin = -rjmin*np.sin(thjmin), rjmin*np.cos(thjmin)

        # Evaluate electron density for nearest spiral arm if it is within
        # some multiple of the e-folding half-width of the arm

        dmin = sqrt((x-xjmin)**2 + (y-yjmin)**2)
        jj = arm_index[j]
        wa = model.wa

        """
        Note armmap used here is to maintain the legacy arm numbering
        from NE2001 and TC93 using as input the arm numbering in Wainscoat.
        New NE model could dispense with this.
        """

        if j== 0: dmin_min = dmin
        Wa = wa * warm[j]
        if dmin < 3.*wa:
        #if dmin < 5.*wa:               # This helps reduce discontinuities
            if dmin <= dmin_min:
                dmin_min = dmin
                whicharm_spiralmodel = j+1

            # arm-distance factor:
            argxy = dmin / Wa
            ga = np.exp(-argxy**2)

            # Galactocentric radial factor:
            if rr > model.Aa: ga *= sech2((rr-model.Aa)/2.)

            # z factor:
            Ha = model.ha * harm[j]
            ga *= sech2(z/Ha)

            # Amplitude re-scalings as in NE2001 code;
            # Ignore cases where these have been turned off in that code (fossils!).

            # TC arm 3:
            if jj == 3:
                th3adeg, th3bdeg = 290, 363
                test3 = thxydeg - th3adeg
                if test3 < 0: test3 += 360
                if 0 <= test3 and test3 < th3bdeg-th3adeg:
                    arg = 2.*pi*test3 / (th3bdeg-th3adeg)
                    fac = ((1 + np.cos(arg))/2)**4
                    # TEMP: fac = 1
                    ga *= fac

            # TC arm 2:
            if jj == 2:
                th2adeg, th2bdeg = 340, 370
                test2 = thxydeg - th2adeg
                fac2min = 0.1
                if test2 < 0: test2 += 360
                if 0 <= test2 and test2 < th2bdeg-th2adeg:
                    arg = 2.*pi*test2 / (th2bdeg-th2adeg)
                    fac = ((1 + fac2min + (1 - fac2min) * np.cos(arg))/2)**3.5
                    # fac = 1    # TEMP
                    ga *= fac

            nea += ga * narm[j] * model.na
            s = sqrt(x**2 + (rsun-y)**2)
            #s = sqrt(x**2 + (8.5-y)**2)
            if verbose:
                print('arms: ', j, jj, whicharm_spiralmodel, index_dsqmin[j],
                    max(0, index_dsqmin[j]-nspline2-1),
                    min(index_dsqmin[j]+nspline2+1, Ncoarse), ind1)
                print(
                    '%5.2f  %5.2f  %d  %d  %d    %4.2f  %4.2f  %4.2f    %4.2f  %7.5f'
                     %(s, rr, j, jj, whicharm_spiralmodel, wa, Wa, dmin, ga, nea))
    if whicharm_spiralmodel == 0:
        whicharm = 0
        Farm = 0
    else:
        whicharm = Darmmap[str(whicharm_spiralmodel-1)]
        #Farm = Dgal['Fa'] * Dgal['farm'+str(whicharm)]    # Intended value

        ###### NOTE ######
        # 2020 Feb 9:
        # found that the Fortran code doesn't use Farm calculated using
        # Fa x farm_j parameters.  This seems to be an error in dmdsm.NE2001.f,
        # which uses Fa to calculate SM for the spiral arms instead of Faval,
        # the value returned by density.NE2001.f.    This only affects the
        # spiral arm values for SM.

        # For now, maintain this error in the Python code because the aim here
        # is to replicate the Fortran code:

        Farm = model.Fa

        # The question is then whether the error was in the code when fitting
        # was done to find the best values of farm_j?   I think the error was
        # not there during the fitting because some of the earlier code versions
        # use Fa * farm_j.

    if verbose:
        print('arms whicharm_spiralmodel  whicharm: ',
            whicharm_spiralmodel, whicharm)

    return nea, Farm, whicharm


# ---------------------------------------------------------------------------

def ne_arms_ne2001p_vec(xvec, yvec, zvec, Ncoarse=20, dthfine=0.01,
                        nfinespline=5, model=None):
    """
    Vectorized version of ne_arms_ne2001p.

    Accepts 1-D numpy arrays xvec, yvec, zvec (length N) and returns arrays
    of the same length::

        nea_v, Farm_v, whicharm_v = ne_arms_ne2001p_vec(xvec, yvec, zvec)

    Algorithm
    ---------
    The inner CubicSpline refinement step is the same as the scalar version
    but batched: positions are grouped by their coarse ``index_dsqmin`` value
    (at most Ncoarse unique groups per arm).  Within each group the same
    angular nodes ``th_k`` are shared, so scipy's 2-D CubicSpline
    ``CubicSpline(th_k, dsq_matrix)`` evaluates all positions in the group
    simultaneously, reducing CubicSpline calls from N*narms to at most
    Ncoarse*narms (e.g. 100 instead of 2500 for N=500).
    """
    if model is None:
        model = default_model

    N = len(xvec)
    narms = model.narms
    nspline2 = int((nfinespline - 1) / 2)

    # dsq_coarse shape: (narms, Ncoarse, N)
    dsq_coarse = (
        (coarse_arms[0, :, :, np.newaxis] - xvec[np.newaxis, np.newaxis, :])**2 +
        (coarse_arms[1, :, :, np.newaxis] - yvec[np.newaxis, np.newaxis, :])**2
    )
    # index_dsqmin shape: (narms, N)  — argmin over the Ncoarse axis
    index_dsqmin = np.argmin(dsq_coarse, axis=1)

    rr      = np.sqrt(xvec**2 + yvec**2)              # (N,)
    thxydeg = np.rad2deg(np.arctan2(-xvec, yvec))     # (N,)
    thxydeg = np.where(thxydeg < 0, thxydeg + 360., thxydeg)

    arm_index = model.arm_index
    warm      = model.warm
    harm      = model.harm
    narm      = model.narm

    nea_v               = np.zeros(N)
    whicharm_spiralmodel_v = np.zeros(N, dtype=int)
    dmin_min_v          = np.full(N, np.inf)

    for j in range(narms):
        jj = arm_index[j]
        Wa = model.wa * warm[j]
        Ha = model.ha * harm[j]

        # ---- Find thjmin for every position using grouped batch splines ----
        thjmin_j = np.empty(N)

        for k in np.unique(index_dsqmin[j]):
            pos_k = np.where(index_dsqmin[j] == k)[0]   # indices into xvec
            n_k   = pos_k.size
            ind1  = range(max(0, k - nspline2 - 1),
                          min(k + nspline2 + 1, Ncoarse), 1)
            th_k  = th1[j, ind1]                         # (len_ind1,)

            # dsq_k shape: (len_ind1, n_k)
            dsq_k    = dsq_coarse[j][np.ix_(list(ind1), pos_k)]
            thjfine  = np.arange(th_k[0], th_k[-1], dthfine)

            if n_k == 1:
                # Scalar path — identical to the original ne_arms_ne2001p
                dsqs      = CubicSpline(th_k, dsq_k[:, 0])
                fine_vals = dsqs(thjfine)                        # (len_fine,)
                thjmin_j[pos_k[0]] = thjfine[fine_vals.argmin()]
            else:
                # Batch path: CubicSpline with 2-D y (columns = positions)
                dsqs      = CubicSpline(th_k, dsq_k)            # y: (len_ind1, n_k)
                fine_vals = dsqs(thjfine)                        # (len_fine, n_k)
                thjmin_j[pos_k] = thjfine[fine_vals.argmin(axis=0)]

        # Arm position and distance at nearest point
        rjmin_j  = model.arm_radius_splines[j](thjmin_j)        # (N,)
        xjmin    = -rjmin_j * np.sin(thjmin_j)
        yjmin    =  rjmin_j * np.cos(thjmin_j)
        dmin_j   = np.sqrt((xvec - xjmin)**2 + (yvec - yjmin)**2)  # (N,)

        # Initialise dmin_min on first arm (unconditional, mirroring scalar code)
        if j == 0:
            dmin_min_v[:] = dmin_j

        in_range = dmin_j < 3. * model.wa    # scalar threshold, same as original

        # Track nearest arm
        update_which       = in_range & (dmin_j <= dmin_min_v)
        dmin_min_v         = np.where(update_which, dmin_j, dmin_min_v)
        whicharm_spiralmodel_v = np.where(update_which, j + 1,
                                          whicharm_spiralmodel_v)

        # ---- Density factor ga (computed for all positions; masked at accumulation) ----
        argxy = dmin_j / Wa
        ga    = np.exp(-argxy**2)

        # Galactocentric radial factor
        ga = np.where(rr > model.Aa, ga * sech2((rr - model.Aa) / 2.), ga)

        # z factor
        ga = ga * sech2(zvec / Ha)

        # Arm-specific angular adjustments
        if jj == 3:   # TC arm 3
            th3adeg, th3bdeg = 290, 363
            test3 = thxydeg - th3adeg
            test3 = np.where(test3 < 0, test3 + 360., test3)
            fac_cond = (test3 >= 0) & (test3 < th3bdeg - th3adeg)
            arg_v    = 2. * np.pi * test3 / (th3bdeg - th3adeg)
            fac      = ((1 + np.cos(arg_v)) / 2)**4
            ga       = np.where(fac_cond, ga * fac, ga)

        if jj == 2:   # TC arm 2
            th2adeg, th2bdeg = 340, 370
            fac2min = 0.1
            test2   = thxydeg - th2adeg
            test2   = np.where(test2 < 0, test2 + 360., test2)
            fac_cond = (test2 >= 0) & (test2 < th2bdeg - th2adeg)
            arg_v    = 2. * np.pi * test2 / (th2bdeg - th2adeg)
            fac      = ((1 + fac2min + (1 - fac2min) * np.cos(arg_v)) / 2)**3.5
            ga       = np.where(fac_cond, ga * fac, ga)

        # Accumulate: only where point is within arm range
        nea_v = np.where(in_range, nea_v + ga * narm[j] * model.na, nea_v)

    # Map Wainscoat arm index to TC arm numbering (same as scalar Darmmap lookup)
    arm_tc_map = np.zeros(narms + 1, dtype=int)   # index 0 → whicharm=0 (no arm)
    for j in range(narms):
        arm_tc_map[j + 1] = arm_index[j]
    whicharm_v = arm_tc_map[whicharm_spiralmodel_v]

    Farm_v = np.where(whicharm_spiralmodel_v != 0, model.Fa, 0.)

    return nea_v, Farm_v, whicharm_v