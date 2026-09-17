# mwprop v2.0 Jan 2026

"""
NE2025 master script

Usage:
    command line:  python NE2025.py l b DM_D ndir [-h]
    in ipython:    run NE2025.py l b DM_D ndir [-h]
    in script:     Dv, Du, Dd = ne2025(l, b, DM_D, ndir)


Change Log:

01/16/2020 initial code sent from SKO --> JMC
01/29/2020 JMC -- mods to make compatible with other new Python modules
02/08/2020 JMC -- add import of ne2001p_config
02/11/2020 JMC -- turn into function + main()
17/12/2020 JMC -- added plotting Cn^2 vs d
01/05/2022 JMC -- moved some definitions to ne2001p_config file, redefined nu -> rf_nu
01/07/2022 JMC -- restructured to use new dmdsm routines,
                  ddmdsm_dm2d_ne2001.py and mdsm_d2dm_ne2001.py
01/09/2022 JMC -- restructed to take plotting outside of the ne2001 function
                  merged dmdsm_dm2d_ne2001.py + dmdsm_d2dm_ne2001.py -> dmdsm_ne2001p.py
03/01/2022 SKO -- merged into package build
03/2022 - 05/2022 SKO -- clump debugging
11/27/2023 SKO -- corrected clump filtering, modified rcmult
12/14/2023 SKO -- all output files (plots, text files) now save to 'output_ne2001p' directory in user's current working directory
01/09/2026 SKO -- NE2001.py and NE2025.py now separate master scripts
"""

import sys,os

global eval_NE2025,eval_NE2001
eval_NE2025 = True
eval_NE2001 = False


# ------------------------------------------------------------
# config_nemod sets up the model with all dictionaries etc:
# ------------------------------------------------------------
from pygedm.mwprop_core.nemod.config_nemod import *

# ----------------------------------------
# Functions to do line of sight integrals:
# ----------------------------------------
from pygedm.mwprop_core.nemod.dmdsm import *

# ------------------------------------------------------
# Functions to calculated derived scattering quantities:
# ------------------------------------------------------
from pygedm.mwprop_core.scattering_functions import scattering_functions2020 as sf

# ---------------------------------------------------------------------
# for additional FRB scattering calculations (integrated to 'infinity')
# ---------------------------------------------------------------------
from pygedm.mwprop_core.scattering_functions import iss_mw_utils2020p_mod as iss_utils


"""
Returns values in two ways:
    * prints to std out in same format as NE2001/Fortran (with some additional lines)
    * dictionaries Dn, Dv, Du, Dd containing names, values, units, and descriptions

To add:
Options:
- Classic option: all output in regular format
- DM->D or D->DM only: no scattering output
- do_analysis option: diagnostic LoS info (slows down the calculation); written to f24, f25 files
- extragalactic LoS:
    * A specified D > 300 kpc (say) implies extragalactic (beyond MW halo)
    * option to include halo estimate of DM
Help file: print out assumptions and other info about the code (README file?)
"""
script_path = os.path.dirname(os.path.realpath(__file__))

basename = sys.argv[0].split('/')[-1].split('.')[0]
now = datetime.datetime.now()
plotstamp = basename + '_' + str(now).split('.')[0]

global rsun, rf_ref, vperp, sikol, louter, linner, ds_fine, ds_coarse
global Dgal01, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, armmap, \
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

# --------------------------------------------------------------------

def ne2025(ldeg, bdeg, dmd, ndir,
    classic=True, dmd_only=True,
    do_analysis=False, plotting=False, verbose=False, debug=False):
    """
    Calculates d from DM or DM from d for ndir >= 0 or ndir < 0

    Input:

    ldeg, bdeg = Galactic longitude and latitude in deg
    dmd = DM (pc/cc) or d (kpc)
    ndir >= 0   DM -> d
          < 0    d -> DM
    dmd_only = True     => compute only DM or distance, no scattering quantities
    classic = True      => print output in format as in Fortran version
    do_analysis = True  => calculate and output diagnostic quantities for the line of sight
    plotting = True     => plot DM vs distance along LoS, SM also if dmd_only = False
    verbose = True      => print results to std out
                           (note not same as 'verbose' internal to functions

    Output: two types:
        Classic: All input and output quantities for DM/D, scattering, scintillations,
                                            emission measure, ...
        Dictionaries for output quantities.
    """
    set_model('ne2025')

    if do_analysis or plotting: # make directory for outputs (need this block here for script calls to ne2001())
        try:
            os.mkdir(os.getcwd()+'/output_ne2025p')
        except FileExistsError:
            pass

    num_outputs = where(dmd_only, 3, 17)

    # Print input parameters before executing code

    if classic:
        print('#NE2025 input: 4 parameters')
        print('{:<1s}{:>10.4f}{:>13s}{:>10s}{:<25s}{:35s}'.format(
            '', ldeg, 'l   ', '', '(deg)', 'GalacticLongitude'))
        print('{:<1s}{:>10.4f}{:>13s}{:>10s}{:<25s}{:35s}'.format(
            '', bdeg, 'b   ', '', '(deg)', 'GalacticLatitude'))
        print('{:<1s}{:>10.4f}{:>13s}{:>10s}{:<25s}{:35s}'.format(
            '', dmd,  'DM/D', '', '(pc cm^{-3}_or_kpc)', 'Input DM or Distance'))
        print('{:<1s}{:>10d}{:>13s}{:>10s}{:<25s}{:35s}'.format(
            '',ndir, 'ndir', '', '1:DM->D;-1:D->DM', 'Which?(DM or D)'))
        print('#NE2025 output: %d values'%(num_outputs))

    # ---------------------------
    # Calculate distance from DM:
    # ---------------------------
    if ndir >= 0:
        dm_target = dmd
        #limit, distmodel, dmpsr, sm, smtau, smtheta, smiso \
        model_out =  dmdsm_dm2d(deg2rad(ldeg), deg2rad(bdeg), dm_target,
                     ds_coarse=ds_coarse, ds_fine=ds_fine, Nsmin=10,
                     dm2d_only=dmd_only, do_analysis=do_analysis,
                     plotting=plotting, verbose=verbose, debug=debug)
        limit_d = model_out[0]
        d_mod = model_out[1]
        dm_mod = model_out[2]

        if limit_d == ' ':
            dmlabel = 'DispersionMeasure'
        else:
            dmlabel = 'DispersionMeasureAsymptote'
        dlabel = 'ModelDistance'

        if dmd_only is False:           # => scattering parameters are calculated
            sm = model_out[3]
            smtau = model_out[4]
            smtheta = model_out[5]
            smiso = model_out[6]

    # ---------------------------
    # Calculate DM from distance:
    # ---------------------------
    if ndir <  0:
        d_target = dmd
        model_out = dmdsm_d2dm(deg2rad(ldeg), deg2rad(bdeg), d_target,
                    ds_coarse=ds_coarse, ds_fine=ds_fine, Nsmin=10,
                    d2dm_only=dmd_only, do_analysis=do_analysis,
                    plotting=plotting, verbose=verbose)

        # The first three entries in model_out are limit, d, and dm for all cases
        limit_d = model_out[0]
        d_mod = model_out[1]
        dm_mod = model_out[2]

        if dmd_only is False:           # => scattering parameters are calculated
            sm = model_out[3]
            smtau = model_out[4]
            smtheta = model_out[5]
            smiso = model_out[6]


        dmlabel = 'ModelDM'
        dlabel = 'Distance'

    dmz = dm_mod*abs(sin(deg2rad(bdeg)))

    if classic:
        print('{:<1s}{:>10.4f}{:>13s}{:>9s} {:<25s}{:35s}'.format(
            limit_d, d_mod,'DIST', '', '(kpc)', dlabel))
        print('{:<1s}{:>10.4f}{:>13s}{:>9s} {:<25s}{:35s}'.format(
            '',dm_mod,'DM  ', '', '(pc cm^{-3})',dmlabel))
        print('{:<1s}{:>10.4f}{:>13s}{:>9s} {:<25s}{:35s}'.format(
            '',dmz,'DMz ', '', '(pc cm^{-3})','DM_Zcomponent'))

    # Also print out scattering, scintillation etc. parameters if requested
    if dmd_only is False:
        tau = sf.tauiss(d_mod, smtau, rf_ref)
        sbw = sf.scintbw(d_mod, smtau, rf_ref)
        stime = sf.scintime(smtau, rf_ref, vperp)
        theta_x = sf.theta_xgal(sm, rf_ref)
        theta_g = sf.theta_gal(smtheta, rf_ref)
        transfreq = sf.transition_frequency(sm, smtau, smtheta, d_mod)
        emsm = sf.em(sm, louter=louter, si=sikol)
        tau_x,sbw_x,deffsm2 = sf.taux_from_thetax(d_mod,sm,smtau,smtheta,theta_x,tau,sbw)

        if classic:
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:<25s}{:35s}'.format(
                '',sm,'SM      ', '', '(kpc m^{-20/3})','ScatteringMeasure'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:<25s}{:35s}'.format(
                '',smtau,'SMtau   ', '','(kpc m^{-20/3})','SM_PulseBroadening'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',smtheta,'SMtheta ',  '','(kpc m^{-20/3})','SM_GalAngularBroadening'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',smiso,'SMiso   ', '','(kpc m^{-20/3})','SM_IsoplanaticAngle'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',emsm,'EM      ', '','(pc cm^{-6})','EmissionMeasure_from_SM'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',tau,'TAU     ', '','(ms)','PulseBroadening @1GHz'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',sbw, 'SBW     ',  '','(MHz)','ScintBW @1GHz'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',stime,'SCINTIME',  '','(s)','ScintTime @1GHz @100 km/s'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',theta_g,'THETA_G ', '', '(mas)','AngBroadeningGal @1GHz'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',theta_x,'THETA_X ', '', '(mas)','AngBroadeningXgal @1GHz'))
            print('{:<1s}{:>10.4f}{:>17s}{:>6s}{:25s}{:35s}'.format(
                '',transfreq,'NU_T    ', '', '(GHz)','TransitionFrequency'))
            print('{:<1s}{:>10.4f}{:>17s}{:>6s}{:<25s}{:35s}'.format(
                '',deffsm2,'DEFFSM2 ', '','(kpc)','EffectiveScreenDistance'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:<25s}{:35s}'.format(
                '',tau_x,'TAU_X   ', '','(ms)','XGalPulseBroadening@1GHz'))
            print('{:<1s}{:>8.4e}{:>17s}{:>6s}{:<25s}{:35s}'.format(
                '',sbw_x, 'SBW_X   ',  '','(MHz)','XGalScintBW@1GHz'))

    # Dictionaries of variable names, values, units, and descriptions
    Dn, Dv, Du, Dd = {}, {}, {}, {}

    if dmd_only is True:
        # DM and distance only:
        Dkeyvalues = list([ldeg, bdeg, dmd, ndir, limit_d, d_mod , dm_mod, dmz])

        Dkeynames = np.array(['l', 'b', 'DM/D', 'ndir', 'limdist', 'DIST', 'DM', 'DMz'])

        Dkeyunits = np.array(['deg', 'deg', 'pc-cm^{-3}_or_kpc',
            '1:DM->D;-1:D->DM', 'blank_or_>>', 'kpc', 'pc-cm^{-3}', 'pc-cm^{-3}'])

        Dkeydesc = np.array(['GalacticLongitude', 'GalacticLatitude', 'Input_DM_or_Distance',
           'Which?(DM_or_D)', 'LowerDistanceLimit', 'ModelDistance', 'DispersionMeasure',
           'DM_Zcomponent'])

        for n, key in enumerate(Dkeynames):
           Dn[key] = key
           Dv[key] = Dkeyvalues[n]
           Du[key] = Dkeyunits[n]
           Dd[key] = Dkeydesc[n]

    else:
    # DM, distance, and scattering output:

        Dkeyvalues = list([ldeg, bdeg, dmd, ndir, limit_d, d_mod, dm_mod, dmz,
            sm, smtau, smtheta, smiso, tau, sbw, stime, theta_g, theta_x,
            transfreq, emsm, deffsm2, tau_x, sbw_x])

        Dkeynames = np.array(['l', 'b', 'DM/D', 'ndir', 'limdist', 'DIST', 'DM', 'DMz',
            'SM', 'SMtau', 'SMtheta', 'SMiso', 'TAU', 'SBW', 'SCINTIME', 'THETA_G', 'THETA_X',
            'NU_T', 'EM', 'DEFFSM2', 'TAU_X', 'SBW_X'])

        Dkeyunits = np.array(['deg', 'deg', 'pc-cm^{-3}_or_kpc', '1:DM->D;-1:D->DM',
            'blank_or_>>', 'kpc', 'pc-cm^{-3}', 'pc-cm^{-3}', 'kpc-m^{-20/3}', 'kpc-m^{-20/3}',
           'kpc-m^{-20/3}', 'kpc-m^{-20/3}',  'ms', 'MHz', 's',
           'mas', 'mas', 'GHz', 'pc-cm^{-6}', 'kpc', 'ms', 'MHz'])

        Dkeydesc = np.array(['GalacticLongitude', 'GalacticLatitude', 'Input_DM_or_Distance',
           'Which?(DM_or_D)', 'LowerDistanceLimit', 'ModelDistance', 'DispersionMeasure',
           'DM_Zcomponent', 'ScatteringMeasure', 'SM_PulseBroadening',
           'SM_GalAngularBroadening', 'SM_IsoplanaticAngle',
           'PulseBroadening@1GHz', 'ScintBW@1GHz',
           'ScintTime@1GHz@100km/s', 'AngBroadeningGal@1GHz', 'AngBroadeningXgal@1GHz',
           'TransitionFrequency', 'EmissionMeasure_from_SM@outer1pc',
           'EffectiveScreenDistance', 'XGalPulseBroadening@1GHz', 'XGalScintBW@1GHz'])

        for n, key in enumerate(Dkeynames):
           Dn[key] = key
           Dv[key] = Dkeyvalues[n]
           Du[key] = Dkeyunits[n]
           Dd[key] = Dkeydesc[n]

    return Dn, Dv, Du, Dd

# --------------------------------------------------------------------

def plot_dm_ne_cn2_vs_d(Dv, f25file = 'f25_dm2d_ne_dsm_vs_s.txt'):
    # Plot if requested and calc_scattering = True:
       s, ne, dsm, dmvss = np.loadtxt(f25file, unpack=True, usecols = (0, 4, 5, 10), skiprows=3)
       DMmax = dmvss.max()

       # For plot labels:
       dm = Dv['DM']
       deffsm2 = Dv['DEFFSM2']
       ldeg = Dv['l']
       bdeg = Dv['b']

       indmax = np.where(ne != 0)[0][-1]
       indkeep = min(int(1.1*indmax), np.size(ne))
       s = s[0:indkeep]
       ne = ne[0:indkeep]
       cn2  = sm_factor *  dsm[0:indkeep]
       dmvss = dmvss[0:indkeep]

       dbar_ne = np.trapz(s*ne) / np.trapz(ne)
       dbar_ne2 = np.trapz(s*ne**2) / np.trapz(ne**2)

       fig = figure(figsize=(6,6))
       subplots_adjust(left=0.15, bottom=0.15)

       ax = fig.add_subplot(311)
       plot(s, dmvss)
       ylabel(r'$\rm DM \ (pc\ cm^{-3}) $', fontsize=13)
       annotate(r'$l = %6.1f^{\circ}  \ \ b = %5.1f$'%(ldeg, bdeg),
           xy=(0.025, 0.875), xycoords='axes fraction', ha='left', va='center', fontsize=10)
       title(r'$\overline{d}_{n_e} = %8.2f \ \ \ \ \ \overline{d}_{n_e^2} = %8.2f \ \ \ \ \ \overline{d}_{\rm SM} = %8.2f \ \ \ \ \ {\rm DM_{max}} =  %8.2f$'%(dbar_ne, dbar_ne2, deffsm2, DMmax))
       tick_params(labelbottom = False)

       ax = fig.add_subplot(312)
       plot(s, ne, label=r'$\rm n_e$')
       ylabel(r'$n_e \ (\rm cm^{-3}) $', fontsize=13)
       tick_params(labelbottom = False)

       ax = fig.add_subplot(313)
       plot(s, cn2, label=r'$\rm C_n^2$')
       ylabel(r'$C_{\rm n}^2 \ (\rm m^{-20/3}) $', fontsize=13)
       xlabel(r'$\rm Distance \ along \ LoS \ \ (kpc) $', fontsize=13)
       annotate(plotstamp, xy=(0.70, 0.02), xycoords='figure fraction', ha='left', va='center', fontsize=5)
       fig.tight_layout()

       output_dir = os.getcwd()+'/output_ne2025p/'

       plotfile = \
            'DM_ne_Cn2_vs_d_l_b_dm_' + '_' + basename + '.pdf'
       savefig(output_dir+plotfile)
       show()

       return

# --------------------------------------------------------------------

# Main

if __name__ == '__main__':

    # If 'explain' option is set, print out README file and exit.
    if '-e' in sys.argv or '--explain' in sys.argv:
         script_path = os.path.dirname(os.path.realpath(__file__)) # SKO -- so that README will always be found
         infile = script_path+'/README.txt' # removed ne2001p
         with open(infile) as fexplain:
             print(fexplain.read())
             fexplain.close()
         #print(open(infile, 'r').read())
         sys.exit()

    try:
        # parse command line arguments and execute
        parser = argparse.ArgumentParser(description='NE2025 (python version Jan. 2026)')
        parser.add_argument('l',help='Galactic longitude (deg)')
        parser.add_argument('b',help='Galactic latitude (deg)')
        parser.add_argument('DM_D',help='DM (pc cm^{-3}) or distance (kpc)')
        parser.add_argument('ndir',help='ndir = 1 (DM->D), ndir = -1 (D->DM)')

        parser.add_argument('-a', '--analysis', action='store_true',
            help='Do line of sight analysis')

        parser.add_argument('-p', '--plotting', action='store_true',
            help='Do diagnostic plotting')

        parser.add_argument('-s', '--scattering', action='store_true',
            help='Calculate scattering and scintillation parameters')

        parser.add_argument('-m', '--modern', action='store_true',
            help='Modern output (turns off classic output like Fortran version)')

        parser.add_argument('-v', '--verbose', action='store_true',
            help='Verbose output (classic output of Fortran version)')

        parser.add_argument('-b', '--debug', action='store_true',
            help='debug output')

        args = parser.parse_args()


    except SystemExit:
        print('Try again with inputs')
        print('Use NE2025 -e to get explanation of code')
        sys.exit()

    ldeg = float(args.l)
    bdeg = float(args.b)
    dmd = float(args.DM_D)
    ndir = int(args.ndir)
    do_analysis = args.analysis
    do_plotting = args.plotting
    calc_scattering = args.scattering
    dmd_only = not calc_scattering
    verbose = args.verbose
    debug = args.debug
    modern = args.modern
    eval_NE2001 = False
    eval_NE2025 = True

    if modern:
        classic = False
    else:
        classic = True

    Dn, Dv, Du, Dd = ne2001(ldeg, bdeg, dmd, ndir,
        classic=classic, verbose=verbose, dmd_only=dmd_only,
        do_analysis=True, plotting=do_plotting, debug=debug,eval_NE2001=eval_NE2001,eval_NE2025=eval_NE2025)

    if calc_scattering and do_plotting:
        output_dir = os.getcwd()+'/output_ne2025p/'
        plot_dm_ne_cn2_vs_d(Dv, f25file = output_dir+'f25_dm2d_ne_dsm_vs_s.txt')

        input('hit return')
        close('all')
