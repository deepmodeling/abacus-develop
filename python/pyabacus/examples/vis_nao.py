from pyabacus import ModuleNAO as nao
import numpy as np
from scipy.special import sph_harm
from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt
import plotly.graph_objects as go

def real_sph_harm(l, m, polar, azimuth):
    if m > 0:
        return np.sqrt(2) * np.real(sph_harm(m, l, azimuth, polar))
    elif m < 0:
        return np.sqrt(2) * np.imag(sph_harm(m, l, azimuth, polar))
    else:
        return np.real(sph_harm(m, l, azimuth, polar))


def cart2sph(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    polar = np.arccos(z/r)
    azimuth = np.arctan2(y, x)
    return r, polar, azimuth

def orbital_plot(orb, itype, l, zeta, m):
    chi = orb(itype, l, zeta)

    # 3d mesh grid
    w = 3
    ngrid = 40
    x, y, z = np.mgrid[-w:w:ngrid*1j, -w:w:ngrid*1j, -w:w:ngrid*1j]

    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    # orbital value on grid
    value = np.zeros_like(x)

    chi_spline = CubicSpline(chi.rgrid, chi.rvalue)
    for i in range(ngrid**3):
        r, polar, azimuth = cart2sph(x[i], y[i], z[i])
        if r < chi.rcut:
            value[i] = chi_spline(r) * real_sph_harm(l, m, polar, azimuth)
    
    # plot
    fig = go.Figure(data=go.Volume(
    x=x,
    y=y,
    z=z,
    value=value,
    isomin=-0.2,
    isomax=0.2,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=20, # needs to be a large number for good volume rendering
    ))
    fig.show()

if __name__ == '__main__':

    orb_dir = '../../../tests/PP_ORB/'
    file_list = ['H_gga_8au_100Ry_2s1p.orb', 'O_gga_10au_100Ry_2s2p1d.orb']
    file_list = [orb_dir + orbfile for orbfile in file_list ]
    nfile = len(file_list)

    # set parameters
    l = 1
    zeta = 0
    m = 0
    itype = 0
    orb = nao.RadialCollection()
    orb.build(nfile, file_list, 'o')

    orbital_plot(orb, itype, l, zeta, m)
