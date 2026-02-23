# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

import sphinx
from m2r import MdInclude
from recommonmark.transform import AutoStructify


def monkeypatch(cls):
    """decorator to monkey-patch methods"""

    def decorator(f):
        method = f.__name__
        old_method = getattr(cls, method)
        setattr(
            cls,
            method,
            lambda self, *args, **kwargs: f(old_method, self, *args, **kwargs),
        )

    return decorator


# workaround until https://github.com/miyakogi/m2r/pull/55 is merged
@monkeypatch(sphinx.registry.SphinxComponentRegistry)
def add_source_parser(_old_add_source_parser, self, *args, **kwargs):
    # signature is (parser: Type[Parser], **kwargs), but m2r expects
    # the removed (str, parser: Type[Parser], **kwargs).
    if isinstance(args[0], str):
        args = args[1:]
    return _old_add_source_parser(self, *args, **kwargs)


sys.path.insert(0, os.path.abspath(".."))

autoclass_content = "both"

# -- Project information -----------------------------------------------------

project = "PyGEDM"
copyright = "2020"
author = "D. C. Price"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "recommonmark",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "display_version": False,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


def setup(app):
    config = {
        # 'url_resolver': lambda url: github_doc_root + url,
        "auto_toc_tree_section": "Contents",
        "enable_eval_rst": True,
    }
    app.add_config_value("recommonmark_config", config, True)
    app.add_transform(AutoStructify)

    # from m2r to make `mdinclude` work
    app.add_config_value("no_underscore_emphasis", False, "env")
    app.add_config_value("m2r_parse_relative_links", False, "env")
    app.add_config_value("m2r_anonymous_references", False, "env")
    app.add_config_value("m2r_disable_inline_math", False, "env")
    app.add_directive("mdinclude", MdInclude)
