# mwprop v2.0 Jan 2026

# Reads NE20x input data into dictionaries

from __future__ import print_function
import numpy as np
import os
import warnings

from scipy.interpolate import CubicSpline

loadtxt = np.loadtxt
genfromtxt = np.genfromtxt

warnings.filterwarnings('ignore','.*Input line 142.*')
warnings.filterwarnings('ignore','.*Input line 16.*')

def read_nemod_parameters():
    """
    Reads standard input files for the model
    and puts them into dictionaries

    NB. Entries like name = name.replace(',', '') take into account
    that the original NE2001 fortran input files have occasional ',' that
    need to be parsed out to get consistent dictionary keys

    JMC 2020 Jan 12-13
    """

    # Note absolute path

    inpdir = os.path.dirname(os.path.realpath(__file__))+'/params/'

    which_mod = loadtxt(inpdir+'which_model.inp',dtype=str)
    if which_mod=='NE2025':
    	eval_NE2025 = True
    	eval_NE2001 = False
    if which_mod=='NE2001':
    	eval_NE2001 = True
    	eval_NE2025 = False
    if eval_NE2001:
        galfile = inpdir + 'gal01.inp'
        negcfile = inpdir + 'ne_gc.NE2001.inp'
        clumpfile = inpdir + 'neclumpN.NE2001.dat'  
        voidfile = inpdir + 'nevoidN.NE2001.dat'
    if eval_NE2025:
        galfile = inpdir + 'gal25.inp'
        negcfile = inpdir + 'ne_gc.NE2025.inp'
        clumpfile = inpdir + 'neclumpN.NE2025.dat'  
        voidfile = inpdir + 'nevoidN.NE2025.dat'

    armsfile = inpdir + 'ne_arms_log_mod.inp'
    nelismfile = inpdir +'nelism.inp'

    # --------------------------------------------------------------------

    # Dictionary from gal01.inp

    Dgal = {}

    # weights:
    names = \
        loadtxt(galfile, skiprows=1, max_rows=1, \
                usecols=(7,8,9,10,11,12,13), dtype=str)
    values = \
        genfromtxt(galfile, skip_header=1, max_rows=1, \
                   usecols=(0,1,2,3,4,5,6), dtype=int)

    for n, name in enumerate(names):
        name = name.replace(',', '')
        Dgal[name] = values[n]


    # component parameters:
    names = loadtxt(galfile, skiprows=2, usecols=(0,), delimiter=':', dtype=str)
    values = loadtxt(galfile, skiprows=2, usecols=(1,), delimiter=':', dtype=float)

    for n, name in enumerate(names):
        Dgal[name] = values[n]

    # --------------------------------------------------------------------

    # Arrays from ne_arms_log_mod.inp 

    Darms = {}
    names = loadtxt(armsfile, skiprows=1, max_rows=1, usecols=(0,1,2,3), dtype=str)
    #aarm, rmin, thmin, extent = \
    xxx = loadtxt(armsfile, skiprows=2, unpack=True, \
                  usecols=(0,1,2,3), dtype=float)
    for n, name in enumerate(names):
        Darms[name] = xxx[n]  

    # --------------------------------------------------------------------

    # Dictionary for Galactic center component
    names1 = loadtxt(negcfile, skiprows=1, usecols=(3,4,5), max_rows=1, dtype=str)
    values1 = loadtxt(negcfile, delimiter=',', skiprows=1, usecols=(0,1,2), max_rows=1, dtype=float)

    names2 = loadtxt(negcfile, skiprows=2, usecols=(1,), dtype=str)
    values2 = loadtxt(negcfile, skiprows=2, usecols=(0,), dtype=float)

    names = np.concatenate((names1, names2))
    values = np.concatenate((values1, values2))

    Dgc = {}
    for n, name in enumerate(names):
        name = name.replace(',', '')
        if 'zgc' in name: name = 'zgc'   
        Dgc[name] = values[n]

    # --------------------------------------------------------------------
    
    # Local ISM
    # has different numbers of fields on each line, so use open + readline

    fism = open(nelismfile, 'r')
    fism.readline()         # Skip comment line
    line = 'x'

    Dlism = {}
    while line != '':
       line  = fism.readline()
       line = line.replace(',', '')
       line = line.replace('\n', '')
       x = line.split(' ')
       N = np.size(x)
       if N > 1:
           N2 = int(N/2)
           fields = x[N2:]
           values = np.array(x[0:N2]).astype(float)
           for n, name in enumerate(fields):
               Dlism[name] = values[n]
    fism.close()
    
    # --------------------------------------------------------------------

    # clumps:
    fields = loadtxt(clumpfile, unpack=True, usecols=(0,1,2,3,4,5,6), max_rows=1, dtype=str) 
    type, names = loadtxt(clumpfile, skiprows=1, usecols=(8, 9), unpack=True,  dtype=str)
    values = loadtxt(clumpfile, skiprows=1, unpack=True, usecols=(1,2,3,4,5,6,7), dtype=float)

    Dclumps = {}
    for n, name in enumerate(names):
        Dclumps[name] = {}
        for m, field in enumerate(fields):
             Dclumps[name][field] = values[m, n]
        
    # --------------------------------------------------------------------

    # voids:
    fields = loadtxt(voidfile, unpack=True, usecols=(0,1,2,3,4,5,6,7,8,9,10), max_rows=1, dtype=str) 
    names = loadtxt(voidfile, skiprows=1, usecols=(12,), unpack=True,  dtype=str)
    values = loadtxt(voidfile, skiprows=1, unpack=True, usecols=(1,2,3,4,5,6,7,8,9,10,11), dtype=float)

    Dvoids= {}
    for n, name in enumerate(names):
        Dvoids[name] = {}
        for m, field in enumerate(fields):
             Dvoids[name][field] = values[m, n]

    # --------------------------------------------------------------------

    return Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001

if __name__ == '__main__':

    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()


