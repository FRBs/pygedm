"""
# create_healpix_map.py

Create a healpix map showing the Galactic DM contribution.
"""

import os

import healpy as hp
import numpy as np
import pylab as plt
import pyymw16

filename = "ymw16_allsky.fits"

if os.path.exists(filename):
    print("File %s exists, skipping write" % filename)
    d_2016 = hp.read_map(filename)
else:
    print("Generating healpix map... (this may take a while)")
    NSIDE = 128
    pix_id = np.arange(hp.nside2npix(NSIDE))
    gl, gb = hp.pix2ang(NSIDE, pix_id, lonlat=True)
    d_2016 = np.zeros_like(pix_id, dtype="float32")

    for ii in pix_id:
        dm, tau = pyymw16.dist_to_dm(gl[ii], gb[ii], 30000)
        d_2016[ii] = dm.value
    print("Writing to %s" % filename)
    hp.write_map(filename, d_2016, coord="G")

hp.mollzoom(np.log(d_2016), title="YMW16", cmap="magma")
plt.show()
