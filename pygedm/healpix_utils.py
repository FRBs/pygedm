
def check_for_healpix_install():
    """ Returns true if python can import healpix """
    try:
        import healpy as hp
        HAS_HEALPIX = True
    except ImportError:  # pragma: no cover
        HAS_HEALPIX = False
    return HAS_HEALPIX