#!/usr/bin/env python
import os
NE_SRC_PATH = os.path.dirname(os.path.abspath(__file__))

FORTRAN_SRCS = ['dmdsm.NE2001.f', 'density.NE2001.f', 'neLISM.NE2001.f',
                'neclumpN.f', 'nevoidN.f', 'scattering98.f']

DATA_FILES   = ['gal01.inp', 'ne_arms_log_mod.inp', 'ne_gc.inp',
                'nelism.inp', 'neclumpN.NE2001.dat', 'nevoidN.NE2001.dat']

def runcmd(cmd):
    """ Run a command with os.system, printing command to screen"""
    print("\n> " + cmd)
    os.system(cmd)

def list_data_files():
    return [os.path.join(NE_SRC_PATH, df) for df in DATA_FILES]

def cleanup():
    """ """
    print("---- Running cleanup ----")
    runcmd('rm *.o')
    if os.path.exists('sgnFile.pyf'):
        runcmd('rm sgnFile.pyf')


def compile_ne2001(fcomp='gfortran', ar_flags='rc'):
    orig_cwd = os.getcwd()
    try:
        os.chdir(NE_SRC_PATH)

        print("\n##########################")
        print("#### Compiling NE2001 ####")
        print("##########################\n")
        print("Working directory: " + NE_SRC_PATH)
        print("Fortran compiler: " + fcomp)

        AR_CMD = 'ar {ar_flags} libNE2001.a '.format(ar_flags=ar_flags)
        RANLIB_CMD = 'ranlib libNE2001.a'

        for fsrc in FORTRAN_SRCS:
            fobj = fsrc.replace('.f', '.o')
            runcmd("{fcomp} -O -std=gnu -c -o {fobj} {fsrc}".format(fcomp=fcomp, fobj=fobj, fsrc=fsrc))
            AR_CMD += '{fobj} '.format(fobj=fobj)

        runcmd(AR_CMD)
        runcmd(RANLIB_CMD)

        print("\n---- Generating F2PY shared object for DMDSM ----")
        runcmd('f2py -m dmdsm -h sgnFile.pyf dmdsm.NE2001.f --overwrite-signature')
        runcmd('f2py -c sgnFile.pyf dmdsm.NE2001.f -L./ -lNE2001 -m dmdsm')
        runcmd('rm sgnFile.pyf')

        print("\n---- Generating F2PY shared object for density ----")
        runcmd('f2py -m density -h sgnFile.pyf density.NE2001.f --overwrite-signature')
        runcmd('f2py -c sgnFile.pyf density.NE2001.f -L./ -lNE2001 -m density')
        runcmd('rm sgnFile.pyf')
    except:
        raise
    finally:
        os.chdir(orig_cwd)

if __name__ == "__main__":
    try:
        compile_ne2001()
        print("Complete.")
    except:
        raise
    finally:
        cleanup()



