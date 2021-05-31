import ne21c
import pprint
import numpy as np

print(ne21c.__doc__)
print(ne21c.dm_to_dist.__doc__)
a = ne21c.dm_to_dist(0, 0, 100)

assert np.isclose(a['dist'], 2.068159818649292)
pprint.pprint(a)

print(ne21c.dist_to_dm.__doc__)
a = ne21c.dist_to_dm(0, 0, a['dist'])
pprint.pprint(a)

print(ne21c.density_xyz.__doc__)
a = ne21c.density_xyz(0, 0, 0)
pprint.pprint(a)

# Compare against pre-computed value
assert np.isclose(a['ne'], 10.04799611523049)

