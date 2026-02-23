import pygedm

print("\n--- DM to dist ---")
print("YMW16:", pygedm.dm_to_dist(100, 0, 250, dm_host=0, method="ymw16"))
print("NE2001:", pygedm.dm_to_dist(100, 0, 250, dm_host=0, method="ne2001"))

print("\n--- dist to DM ---")
print("YMW16:", pygedm.dist_to_dm(100, 0, 250, method="ymw16"))
print("NE2001:", pygedm.dist_to_dm(100, 0, 250, method="ne2001"))

print("\n--- Electron density ---")
print("YMW16:", pygedm.calculate_electron_density_xyz(0, 0, 0, method="ymw16"))
print("NE2001:", pygedm.calculate_electron_density_xyz(0, 0, 0, method="ne2001"))
