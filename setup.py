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
import setuptools

__version__ = '2.0.6'
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
            'src/main.cpp',
            'src/dora.cpp',
            'src/fermibubble.cpp',
            'src/frb_d.cpp',
            'src/galcen.cpp',
            'src/gum.cpp',
            'src/lmc.cpp',
            'src/localbubble.cpp',
            'src/ne_crd.cpp',
            'src/nps.cpp',
            'src/smc.cpp',
            'src/spiral.cpp',
            'src/thick.cpp',
            'src/thin.cpp',
            'src/ymw16par.cpp',
            'src/dmdtau2.cpp',
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            os.path.join(__here__, 'src')
        ],
        extra_link_args=['-lm'],
        language='c++'
    ),
]


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


setup(
    name='pyymw16',
    version=__version__,
    author='D. Price',
    author_email='dancpr@berkeley.edu',
    description='Python/C++ version of YMW16 electron density model',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/telegraphic/pyymw16',
    download_url='https://github.com/telegraphic/pyymw16/archive/%s.tar.gz' % __version__,
    python_requires='>=2.7',
    install_requires=['pybind11>=2.2', astro],
    tests_require = ['pytest<3.7', astro, 'numpy'],
    setup_requires = ['pytest-runner', 'pytest-cov', 'pybind11>=2.2'],
    ext_modules=ext_modules,
    packages=['pyymw16'],
    package_data={'pyymw16': ['spiral.txt', 'ymw16par.txt']},
    include_package_data=True,
    zip_safe=False,
    cmdclass={'build_ext': BuildExt},
)
