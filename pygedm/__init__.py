try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("pygedm")
except PackageNotFoundError:  # pragma: no cover
    # Package is not installed, use a fallback version
    __version__ = "4.0.0"

from .pygedm import (
    calculate_electron_density_lbr,
    calculate_electron_density_xyz,
    calculate_halo_dm,
    convert_lbr_to_xyz,
    dist_to_dm,
    dm_to_dist,
    generate_healpix_dm_map,
)
