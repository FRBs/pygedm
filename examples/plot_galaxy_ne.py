#!/usr/bin/env python
"""
# plot_galaxy_ne.py

Plot electron density of galaxy using ymw16.ne_crd
"""

import numpy as np
import pylab as plt

import pygedm

# Create empty arrays for NE2001 and YMW16 model output
d_2016 = np.zeros((100, 100))
d_2001 = np.zeros((100, 100))

# Loop through pixels
for ii in range(100):
    for jj in range(100):
        # generate X, Y, Z in PC from indexes
        x, y, z = (ii - 50) * 400, (jj - 50) * 400, 0

        # Compute EDM
        d_2016[ii, jj] = pygedm.calculate_electron_density_xyz(
            x, y, z, method="ymw16"
        ).value
        d_2001[ii, jj] = pygedm.calculate_electron_density_xyz(
            x, y, z, method="ne2001"
        ).value


# Plot output images
plt.rcParams["font.size"] = 14

plt.figure(figsize=(4 * 3 * 2, 3 * 3))
plt.subplot(1, 2, 1)
plt.title("NE2001 electron density model, Z=0 plane")
plt.imshow(
    np.log10(d_2001),
    extent=(-50 * 400 / 1e3, 50 * 400 / 1e3, -50 * 400 / 1e3, 50 * 400 / 1e3),
    clim=(-5, 0),
)
plt.xlabel("X [kpc]")
plt.ylabel("Y [kpc]")
cbar = plt.colorbar()
cbar.set_label("log10($N_e$)")

plt.subplot(1, 2, 2)
plt.title("YMW16 electron density model, Z=0 plane")
plt.imshow(
    np.log10(d_2016),
    extent=(-50 * 400 / 1e3, 50 * 400 / 1e3, -50 * 400 / 1e3, 50 * 400 / 1e3),
    clim=(-5, 0),
)
plt.xlabel("X [kpc]")
plt.ylabel("Y [kpc]")
cbar = plt.colorbar()
cbar.set_label("log10($N_e$)")

plt.savefig("plot_galaxy_ne.png")
plt.show()
