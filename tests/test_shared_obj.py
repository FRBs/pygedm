import os

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from pygedm.ymw16_wrapper import ymw16

DATAPATH = str(files("pygedm"))
gl, gb, dist, dm_host, ndir, mode_id, vbs, txt = 204.0, -6.5, 2000, 0, 2, -1, 0, ""
a = ymw16.dmdtau(gl, gb, dist, dm_host, ndir, mode_id, vbs, DATAPATH, txt)
print(a)
