from pygedm.ymw16_wrapper import ymw16
import os
from pkg_resources import resource_filename

DATAPATH = os.path.dirname(resource_filename("pygedm", "spiral.txt"))
gl, gb, dist, dm_host, ndir, mode_id, vbs, txt = 204.0, -6.5, 2000, 0, 2, -1, 0, ''
a = ymw16.dmdtau(gl, gb, dist, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
print(a)
