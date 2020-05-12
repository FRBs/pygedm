from pkg_resources import resource_filename, get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('pygedm').version
except DistributionNotFound: # pragma: no cover
    __version__ = '0.0.1 - manual install'

from .pygedm import dm_to_dist, dist_to_dm, calculate_electron_density_lbr, calculate_electron_density_xyz, \
    calculate_halo_dm, generate_healpix_dm_map, convert_lbr_to_xyz