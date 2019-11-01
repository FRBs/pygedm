# To increment version
# Check you have ~/.pypirc filled in
# git tag x.y.z
# git push --tags
# python setup.py sdist bdist_wheel
# twine upload dist/*

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import glob
import setuptools

__version__ = '3.0.4'
__here__ = os.path.abspath(os.path.dirname(__file__))

if sys.version_info.major == 3:
    astro = "astropy"
else:
    astro = "astropy<3.0"

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'ymw16',
        sources=[
            'ymw16_src/main.cpp',
            'ymw16_src/dora.cpp',
            'ymw16_src/fermibubble.cpp',
            'ymw16_src/frb_d.cpp',
            'ymw16_src/galcen.cpp',
            'ymw16_src/gum.cpp',
            'ymw16_src/lmc.cpp',
            'ymw16_src/localbubble.cpp',
            'ymw16_src/ne_crd.cpp',
            'ymw16_src/nps.cpp',
            'ymw16_src/smc.cpp',
            'ymw16_src/spiral.cpp',
            'ymw16_src/thick.cpp',
            'ymw16_src/thin.cpp',
            'ymw16_src/ymw16par.cpp',
            'ymw16_src/dmdtau2.cpp',
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            os.path.join(__here__, 'ymw16_src'),
        ],
        extra_link_args=['-lm'],
        language='c++'
    ),
]

ymw16_data_files  = ['spiral.txt', 'ymw16par.txt']
ne2001_data_files = ['gal01.inp', 'ne_arms_log_mod.inp', 'ne_gc.inp',
                     'nelism.inp', 'neclumpN.NE2001.dat', 'nevoidN.NE2001.dat']
data_files = ymw16_data_files + ne2001_data_files

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


with open("README.md", "r") as fh:
    long_description = fh.read()

# Compile NE2001
NE_SRC_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ne2001_src')
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
            runcmd("{fcomp} -O -fPIC -std=gnu -c -o {fobj} {fsrc}".format(fcomp=fcomp, fobj=fobj, fsrc=fsrc))
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

####
# COMPILE NE2001 FORTRAN
####


do_fortran_compile = True

if do_fortran_compile:
    compile_ne2001()

ne2001_shared_objs = list(glob.glob(NE_SRC_PATH + '/*.so'))
for sobj in ne2001_shared_objs:
    print("Copying {sobj} to pygedm/".format(sobj=os.path.basename(sobj)))
    os.system('cp {sobj} pygedm/'.format(sobj=sobj))
    data_files.append(os.path.basename(sobj))

setup(
    name='pygedm',
    version=__version__,
    author='D. C. Price',
    author_email='dancpr@berkeley.edu',
    description='Python/C++ version of YMW16 electron density model',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/telegraphic/pygedm',
    download_url='https://github.com/telegraphic/pygedm/archive/%s.tar.gz' % __version__,
    python_requires='>=2.7',
    install_requires=['pybind11>=2.2', astro],
    tests_require= ['pytest<3.7', astro, 'numpy'],
    setup_requires= ['pytest-runner', 'pytest-cov', 'pybind11>=2.2'],
    ext_modules=ext_modules,
    packages=['pygedm'],
    package_data={'pygedm': data_files},
    include_package_data=True,
    zip_safe=False,
    cmdclass={'build_ext': BuildExt},
)
