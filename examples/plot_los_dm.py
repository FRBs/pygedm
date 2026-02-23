#!/usr/bin/env python
"""
# plot_los_dm.py

Plot electron density from Earth-centric Galactic coordinates, out to a distance
of 10,000 parsecs.
"""

import numpy as np
import pylab as plt
import pyymw16

d = np.zeros((100, 100))

for ii in range(100):
    for jj in range(100):
        gl, gb, dist = (ii - 50) * 360.0 / 100, (jj - 50) * 180.0 / 100, 10000
        dm, tau = pyymw16.dist_to_dm(gl, gb, dist)
        d[ii, jj] = dm.value

plt.figure(figsize=(8, 5))
plt.title("YMW16 Galactic DM (to 10,000 pc)")
plt.imshow(d.T, extent=(-180, 180, -90, 90))
plt.ylim(-45, 45)
plt.xlabel("gl [deg]")
plt.ylabel("gb [deg]")
cbar = plt.colorbar(orientation="horizontal")
cbar.set_label("DM [pc cm$^{-3}$]")
plt.savefig("plot_los_dm.png")
plt.show()
