# Modern build configuration for pygedm with C++ extensions
# All package metadata is now in pyproject.toml following PEP 621
#
# To build and publish:
#   python -m build
#   twine upload dist/*

import os
import sys
import tempfile
import warnings
from distutils.ccompiler import new_compiler
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


def have_f2c(include_dirs, library_dirs):
    """Check whether f2c.h and libf2c are available, so we can build ne21c.

    NE2001 is the only extension that needs f2c, and f2c isn't on PyPI, so
    we skip building it (rather than failing the whole install) when f2c
    isn't present.
    """
    compiler = new_compiler()
    with tempfile.TemporaryDirectory() as tmpdir:
        src = os.path.join(tmpdir, "f2c_check.c")
        with open(src, "w") as f:
            f.write('#include "f2c.h"\nint main(void) { return 0; }\n')
        try:
            objs = compiler.compile([src], output_dir=tmpdir, include_dirs=include_dirs)
            compiler.link_executable(
                objs, "f2c_check", output_dir=tmpdir, library_dirs=library_dirs, libraries=["f2c"]
            )
            return True
        except Exception:
            return False


ne21c_include_dirs = [
    os.path.join(__here__, "ne21c"),
    f"{prefix}/include" if prefix else "/usr/include",
]
ne21c_library_dirs = [f"{prefix}/lib"] if prefix else []

HAS_F2C = have_f2c(ne21c_include_dirs, ne21c_library_dirs)
if not HAS_F2C:
    warnings.warn(
        "f2c not found -- skipping the ne21c extension. pygedm will install "
        "without method='ne2001' support; use method='ne2001p' or 'ne2025' "
        "instead, or install f2c and reinstall pygedm for the compiled NE2001."
    )

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
]

if HAS_F2C:
    ext_modules.append(
        Extension(
            "ne21c",
            sources=[
                "ne21c/main.cpp",
            ],
            include_dirs=[
                # Path to pybind11 headers
                get_pybind_include(),
                get_pybind_include(user=True),
                *ne21c_include_dirs,
            ],
            library_dirs=ne21c_library_dirs,
            extra_compile_args=["-Wno-write-strings"],
            extra_link_args=["-lm", "-lf2c"],
            language="c++",
        )
    )


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
