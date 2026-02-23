"""
# compare_to_ne2001.py

Compares output to that of pyne2001, plots galaxic contribution
"""

import os

import numpy as np
import pylab as plt
import pyne2001
import pyymw16

import hickle as hkl


def plot_gplane(d, title="", cbar_label="DM [pc cm$^{-3}$]"):
    plt.imshow(d.T[::-1, ::-1], extent=(-180, 180, -90, 90), cmap="magma")
    plt.ylim(-60, 60)
    cb = plt.colorbar()
    cb.set_label(cbar_label)
    plt.title(title)
    plt.xticks([-180, -120, -60, 0, 60, 120, 180])
    plt.yticks([-60, -30, 0, 30, 60])
    plt.ylabel("gb [deg]")
    plt.minorticks_on()


if not os.path.exists("gedm.hkl"):
    N = 256
    d_2001 = np.zeros([N, N])
    d_2016 = np.zeros([N, N])

    for ii in range(N):
        print("%i / %i" % (ii + 1, N))
        for jj in range(N):
            l = float(ii) / N * 360 - 180
            b = float(jj) / N * 90 - 45
            d_2001[ii, jj] = pyne2001.get_galactic_dm(l, b)
            dm, tau = pyymw16.dist_to_dm(l, b, 30000)
            d_2016[ii, jj] = dm.value

    hkl.dump({"NE2001": d_2001, "YMW16": d_2016}, "gedm.hkl")
else:
    d = hkl.load("gedm.hkl")
    plt.figure(figsize=(9, 9))
    plt.subplot(3, 1, 1)
    plot_gplane(d["NE2001"], "NE2001")

    plt.subplot(3, 1, 2)
    plot_gplane(d["YMW16"], "YMW16")

    plt.subplot(3, 1, 3)
    d_delta = d["YMW16"] - d["NE2001"]
    plot_gplane(d_delta, "Difference")
    plt.xlabel("gl [deg]")
    plt.tight_layout()
    plt.savefig("compare_to_ne2001.png")
    plt.show()
