"""make_datacube.py -- Generate datacube for NE2025."""
import pygedm
import numpy as np
import hickle as hkl
from tqdm import tqdm

skip_ddm = True
skip_dmd = False
skip_xyz = True

gl = np.linspace(-180, 180, 360*2+1)
gb = np.linspace(-90, 90, 180*2+1)
dist = np.array((0.1, 0.2, 0.5, 1, 2, 5, 8.5, 10, 20, 50))
dml   = np.array((1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000))

hkl.dump({'gl': gl, 'gb': gb,'dist': dist, 'dml': dml}, 'datadims_ne2025.hkl')

if not skip_xyz:
    print("Generating galactic NE XYZ cube")
    # Create empty arrays for NE2001 and YMW16 model output
    d_2025 = np.zeros((100, 100, 100))

    # Loop through pixels
    for ii in tqdm(range(100)):
        for jj in range(100):
            for kk in range(100):
                # generate X, Y, Z in PC from indexes
                x, y, z = (ii - 50) * 400, (jj - 50) * 400, (kk - 50) * 400

                # Compute EDM
                d_2025[ii, jj, kk] = pygedm.calculate_electron_density_xyz(x, y, z, method='ne2025').value


    ii = np.array(range(100))
    jj = np.array(range(100))
    kk = np.array(range(100))

    x, y, z = (ii - 50) * 400, (jj - 50) * 400, (kk - 50) * 400

    gdict = {
        'x': x,
        'y': y,
        'z': z,
        'ne': d_2025
    }

    hkl.dump(gdict, 'ne2025_xyz.hkl')


if not skip_ddm:
    print("Generating D->DM cube")
    ddm_data = np.zeros((len(dist), len(gb), len(gl)))
    ddm_data_tau = np.zeros((len(dist), len(gb), len(gl)))


    for igl, ll in tqdm(list(enumerate(gl))):
        for igb, bb in enumerate(gb):
            for idist, dd in enumerate(dist):
                dd *= 1000
                dm, tau = pygedm.dist_to_dm(ll, bb, dd, method='ne2025')
                ddm_data[idist, igb, igl] = dm.value
                ddm_data_tau[idist, igb, igl] = tau.value

    hkl.dump({'ddm': ddm_data, 'tau': ddm_data_tau}, 'ne2025_ddm.hkl')


if not skip_dmd:
    print("Generating DM->D cube")

    dmd_data = np.zeros((len(dml), len(gb), len(gl)))
    dmd_data_tau = np.zeros((len(dml), len(gb), len(gl)))

    for igl, ll in tqdm(list(enumerate(gl))):
        for igb, bb in enumerate(gb):
            for idm, dd in enumerate(dml):

                dist, tau =  pygedm.dm_to_dist(ll, bb, dd, method='ne2025')
                dmd_data[idm, igb, igl] = dist.value
                dmd_data_tau[idm, igb, igl] = tau.value

    hkl.dump({'dmd': dmd_data, 'tau': dmd_data_tau}, 'ne2025_dmd.hkl')
