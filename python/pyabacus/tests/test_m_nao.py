from __future__ import annotations

# import pytest
from pyabacus import ModuleNAO as nao
from pyabacus import ModuleBase as base
import numpy as np

def test_int_cal():
    file_list = ["../../tests/PP_ORB/O_gga_10au_100Ry_2s2p1d.orb"]
    nfile = 1
    orb = nao.RadialCollection()
    orb.build(nfile, file_list, 'o')
    sbt = base.SphericalBesselTransformer()
    orb.set_transformer(sbt)

    rmax = orb.rcut_max * 2.0
    dr = 0.01
    nr = int(rmax / dr) + 1

    orb.set_uniform_grid(True, nr, rmax, 'i', True)
    S_intor = nao.TwoCenterIntegrator()
    S_intor.tabulate(orb, orb, 'S', nr, rmax)
    vR = np.array([1.0,0,0]).astype(np.float64)
    out = S_intor.calculate(0, 1, 0, 1, 0, 1, 0, 1, vR)[0]
    print(out)

test_int_cal()

    