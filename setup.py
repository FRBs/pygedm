# Modern build configuration for pygedm with C++ extensions
# All package metadata is now in pyproject.toml following PEP 621
#
# To build and publish:
#   python -m build
#   twine upload dist/*

import os
import sys
import tempfile
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

__here__ = os.path.abspath(os.path.dirname(__file__))


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked."""

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


prefix = os.environ.get("CONDA_PREFIX")


ext_modules = [
    Extension(
        "ymw16",
        sources=[
            "ymw16_src/main.cpp",
            "ymw16_src/dora.cpp",
            "ymw16_src/fermibubble.cpp",
            "ymw16_src/frb_d.cpp",
            "ymw16_src/galcen.cpp",
            "ymw16_src/gum.cpp",
            "ymw16_src/lmc.cpp",
            "ymw16_src/localbubble.cpp",
            "ymw16_src/ne_crd.cpp",
            "ymw16_src/nps.cpp",
            "ymw16_src/smc.cpp",
            "ymw16_src/spiral.cpp",
            "ymw16_src/thick.cpp",
            "ymw16_src/thin.cpp",
            "ymw16_src/ymw16par.cpp",
            "ymw16_src/dmdtau2.cpp",
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            os.path.join(__here__, "ymw16_src"),
        ],
        extra_link_args=["-lm"],
        language="c++",
    ),
    Extension(
        "ne21c",
        sources=[
            "ne21c/main.cpp",
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            os.path.join(__here__, "ne21c"),
            f"{prefix}/include" if prefix else "/usr/include",
        ],
        library_dirs=[f"{prefix}/lib"] if prefix else [],
        extra_compile_args=["-Wno-write-strings"],
        extra_link_args=["-lm", "-lf2c"],
        language="c++",
    ),
]


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
        return True
    except Exception:
        return False
    finally:
        try:
            os.unlink(fname)
        except Exception:
            pass


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.

    The c++17 is preferred, then C++14, then C++11.
    """
    for flag in ["-std=c++17", "-std=c++14", "-std=c++11"]:
        if has_flag(compiler, flag):
            return flag
    raise RuntimeError("Unsupported compiler -- at least C++11 support is needed!")


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""

    c_opts = {
        "msvc": ["/EHsc"],
        "unix": [],
    }

    if sys.platform == "darwin":
        c_opts["unix"] += ["-stdlib=libc++", "-mmacosx-version-min=10.7"]

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == "unix":
            opts.append("-DVERSION_INFO=%s" % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
        elif ct == "msvc":
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
        build_ext.build_extensions(self)


# Run setup with configuration from pyproject.toml
setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
)
