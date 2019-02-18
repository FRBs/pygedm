#!/usr/bin/env python
"""
# plot_galaxy_ne.py

Plot electron density of galaxy using ymw16.ne_crd
"""

import pyymw16
import pylab as plt
import numpy as np

d = np.zeros((100, 100))

for ii in range(100):
    for jj in range(100):
            x, y, z = (ii - 50) * 400, (jj - 50) * 400, 0
            d[ii, jj] = pyymw16.calculate_electron_density_xyz(x, y, z)

plt.figure(figsize=(4*3, 3*3))
plt.title("YMW16 electron density model, Z=0 plane")
plt.imshow(np.log10(d), extent=(-50*400, 50*400, -50*400, 50*400))
plt.xlabel("X [pc]")
plt.ylabel("Y [pc]")
cbar = plt.colorbar()
cbar.set_label("log10($N_e$)")
plt.savefig("plot_galaxy_ne.png")
plt.show()
