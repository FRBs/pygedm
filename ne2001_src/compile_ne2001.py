#!/usr/bin/env python
import os

FORTRAN_COMPILER = 'gfortran'

FORTRAN_SRCS = ['dmdsm.NE2001.f', 'density.NE2001.f', 'neLISM.NE2001.f',
                'neclumpN.f', 'nevoidN.f', 'scattering98.f']


def runcmd(cmd):
    """ Run a command with os.system, printing command to screen"""
    print("\n> " + cmd)
    os.system(cmd)


def cleanup():
    """ """
    print("---- Running cleanup ----")
    runcmd('rm *.o')
    if os.path.exists('sgnFile.pyf'):
        runcmd('rm sgnFile.pyf')


def compile_ne2001():
    AR_CMD = 'ar rc libNE2001.a '
    RANLIB_CMD = 'ranlib libNE2001.a'

    for fsrc in FORTRAN_SRCS:
        fobj = fsrc.replace('.f', '.o')
        runcmd("{fcomp} -O -std=gnu -c -o {fobj} {fsrc}".format(fcomp=FORTRAN_COMPILER, fobj=fobj, fsrc=fsrc))
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

if __name__ == "__main__":
    try:
        compile_ne2001()
        print("Complete.")
    except:
        raise
    finally:
        cleanup()



