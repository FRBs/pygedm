"""
Python implementation of Yamasaki & Totani DM Halo model

References:
    [1] `Yamasaki S, Totani T (2020), <https://ui.adsabs.harvard.edu/abs/2019arXiv190900849Y/abstract>`_
    *The Galactic Halo Contribution to the Dispersion Measure of Extragalactic Fast Radio Bursts*
    The Astrophysical Journal, Volume 888, Issue 2, id.105

Notes:
    Adapted from S. Yamasaki's ``DM_halo_yt2020_numerical.py`` command-line python code
"""
import numpy as np
from scipy import integrate
from numpy import cos, sin, sqrt, fabs
from astropy import units as u

### numerical constants ###

U_G    = 6.67                           # G = 6.67259e-8 cm^3/g/s^2
U_Msun = 1.99                           # 1 M_sun = 1.99e+33 g
U_pc   = 3.09                           # 1 pc = 3.086e+18 cm
U_mp   = 1.67                           # m_p = 1.6726231e-24 g
U_eV   = 1.60                           # 1 eV = 1.602177e-12 erg
c_scaling = U_Msun/U_pc**3*pow(10,-18)  # scaling factor [10^{12} M_sun/kpc^3] -> [g/cm^3]
m_p    = 1.6726231*pow(10,-24)          # proton mass [g]

### key parameters ###

T      = 0.3                    # halo temperature [keV]
M_b    = 0.036                  # total halo gas mass [10^{12} Z^{-1} M_sun] (such that M_b = 0.12 for Z = Z_halo)
Z_halo = 0.3                    # halo metallicity [Z_sun]
xi_H   = 28.0/34.0              # abundance ratio of hydrogen to electron
mu_e   = 1.18                   # mean molecular mass per electron
mu     = 0.62                   # mean particle mass
D_sun  = 8.5                    # Galactocentoric distance of the Sun [kpc]
c_NFW  = 12                     # NFW concentration
r_vir  = 260                    # virial radius of MW [kpc]
M_vir  = 1.0                    # virial mass of MW [10^{12} M_sun]
r_s    = r_vir/c_NFW            # NFW scale radius [kpc]
f_NFW  = np.log(1.0+c_NFW)-c_NFW/(1.0+c_NFW) # f(x)=np.log(1.0+x)-x/(1.0+x) and x = c_NFW
rho_s  = M_vir/4.0/np.pi/r_s**(3.0)/f_NFW    # NFW scale density [10^{12} M_sun/kpc^3]

### calculate n_0^{sphe} (eq.[3]) by normalization (eq.[4]) ###

r_array = np.linspace(r_s/100,r_vir,100000)                         # r = 0-r_vir [kpc]
Upsilon = 4*np.pi*10*U_G*U_Msun/U_pc*U_mp/U_eV*r_s**2*rho_s*mu/T    # Upsilon appearing in right-hand side of eq.(3)

# integrated function excluding n_0^{sphe} (left-hand side of eq.[4])
y_sphe  = 4.0*np.pi*r_array**2*mu_e*m_p/c_scaling*np.exp(-Upsilon*(1-np.log(1+r_array/r_s)*r_s/r_array))

I_sphe   = integrate.simps(y_sphe,r_array)      # integration
n_0_sphe = M_b/Z_halo/I_sphe                    # central electron number density for spherical comp. [cm^{-3}]
#print('Upsilon = %f' % Upsilon)
#print('n_0_sphe = %f [cm^{-3}]' % n_0_sphe)

### key functions ###


def s_max(l, b): # maximum integration limit corresponsing to r = r_vir as function of sky coordinates [kpc]
    """ Compute integration limit s_max for given sky coordinates

    Args:
        l (float): Galactic longitude, in radians (-pi to +pi)
        b (float): Galactic latitude, in radians  (-pi/2 to pi/2)

    Returns:
        s_max (float), maximum integration limit corresponsing to r = r_vir
    """
    return D_sun*cos(b)*cos(l)+sqrt(D_sun**2*(cos(b)**2*cos(l)**2-1)+(1.0*r_vir)**2)


def ne_sphe(l, b, s):
    """ Compute electron density for spherical component for (l, b) at distance s

    Args:
        l (float): Galactic longitude, in radians (-pi to +pi)
        b (float): Galactic latitude, in radians  (-pi/2 to pi/2)
        s (float): Distance (kpc)

    Returns:
        ne (float): electron density in cm^{-3}
    """
    R = sqrt(D_sun**2+(s*cos(b))**2-2*D_sun*s*cos(b)*cos(l))
    z = s*sin(b)
    r = sqrt(R**2+z**2)
    Upsilon = 4*np.pi*10*U_G*U_Msun/U_pc*U_mp/U_eV*r_s**2*rho_s*mu/T
    return n_0_sphe*np.exp(-Upsilon*(1-np.log(1+r/r_s)*r_s/r)) # [cm^{-3}]


def ne_disk(l, b, s):
    """ Compute electron density for spherical component for (l, b) at distance s

    Args:
        l (float): Galactic longitude, in radians (-pi to +pi)
        b (float): Galactic latitude, in radians  (-pi/2 to pi/2)
        s (float): Distance (kpc)

    Returns:
        ne (float): electron density in cm^{-3}
    """
    R = sqrt(D_sun**2+(s*cos(b))**2-2*D_sun*s*cos(b)*cos(l))
    z = fabs(s*sin(b))
    n_0_disk = 7.4*pow(10,-3) # [Z^{-1} cm^{-3}]
    z_0 = 2.4 # [kpc]
    R_0 = 4.9 # [kpc]
    return n_0_disk/np.exp(R/R_0+z/z_0)/Z_halo # [cm^{-3}]


### vectorize density profile functions ###
vfunc_sphe = np.vectorize(ne_sphe)
vfunc_disk = np.vectorize(ne_disk)

def calculate_halo_dm(l, b, component='both'):
    """ Compute halo DM

    Args:
        l (float): Galactic longitude, in degrees (-180 to +180)
        b (float): Galactic latitude, in degrees  (-90 to 90)
        component (str): Compute 'spherical' component of halo,
                         'disk' component, or 'both' components.

    Returns:
        DM (float): Dispersion measure in [pc/cm^3]
    """

    ### input coordinates (l,b) in degrees ###
    l_FRB = np.deg2rad(l)
    b_FRB = np.deg2rad(b)

    ### calculate DM_halo ###
    N = int(10000)
    xx = np.linspace(0.01, s_max(l_FRB, b_FRB), N)
    y_DM_sphe = vfunc_sphe(l_FRB, b_FRB, xx)
    y_DM_disk = vfunc_disk(l_FRB, b_FRB, xx)

    if component == 'spherical' or component == 'both':
        DM_sphe = 1000 * integrate.simps(y_DM_sphe, xx) * u.pc / u.cm**3 # [pc/cm^3]
    if component == 'disk' or component == 'both':
        DM_disk = 1000 * integrate.simps(y_DM_disk, xx) * u.pc / u.cm**3 # [pc/cm^3]

    if component == 'both':
        return DM_sphe + DM_disk
    elif component == 'spherical':
        return DM_sphe
    elif component == 'disk':
        return DM_disk
    else:
        raise RuntimeError("Component must be 'both', 'spherical', or 'disk' ")


def calculate_halo_dm_analytic(l, b):
    """ Calculate the DM contribution of the Galactic halo.

    Uses an analytical formula for speed. Useful for all-sky mapping.

    Args:
        l (float): Galactic longitude, in degrees (-180 to +180)
        b (float): Galactic latitude, in degrees  (-90 to 90)
    """
    YT19_coeffs = np.array([[250.12, -871.06, 1877.5, -2553.0, 2181.3, -1127.5, 321.72, -38.905],
                            [-154.82, 783.43, -1593.9, 1727.6, -1046.5, 332.09, -42.815, 0],
                            [-116.72, -76.815, 428.49, -419.00, 174.60, -27.610, 0, 0],
                            [216.67, -193.30, 12.234, 32.145, -8.3602, 0, 0, 0],
                            [-129.95, 103.80, -22.800, 0.44171, 0, 0, 0, 0],
                            [39.652, -21.398, 2.7694, 0, 0, 0, 0, 0],
                            [-6.1926, 1.6162, 0, 0, 0, 0, 0, 0],
                            [0.39346, 0, 0, 0, 0, 0, 0, 0]])
    l_abs = np.abs(np.deg2rad(l))
    b_abs = np.abs(np.deg2rad(b))
    DM_halo = 0
    for ii in range(8):
        for jj in range(8):
            DM_halo += YT19_coeffs[ii, jj] * l_abs ** ii * b_abs ** jj
    return DM_halo * u.pc / u.cm**3