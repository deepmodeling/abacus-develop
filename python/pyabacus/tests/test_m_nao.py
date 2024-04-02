from __future__ import annotations

import pytest
from pyabacus import ModuleNAO as nao
from pyabacus import ModuleBase as base
import numpy as np

def test_nr():
    l = 1
    dr = 0.01
    sz = 5000
    itype = 3
    pr = -1
    izeta = 5
    symbol = "Au"
    grid = np.empty(sz, dtype=np.float64)
    f = np.empty(sz, dtype=np.float64)
    for i in range(sz):
        r = i * dr
        grid[i] = r
        f[i] = np.exp(-r)
    chi = nao.NumericalRadial()
    chi.build(l, True, sz, grid, f, pr, izeta, symbol, itype)
    assert chi.symbol == symbol
    assert chi.izeta == izeta
    assert chi.itype == itype
    assert chi.l == l
    assert chi.nr == sz
    assert chi.nk == 0
    assert chi.rmax == grid[sz-1]
    for i in range(sz):
        assert(chi.rgrid[i] == grid[i])
        assert(chi.rvalue[i] == f[i])
    assert chi.pr == pr
    assert chi.pk == 0
    assert chi.is_fft_compliant == False 

def test_rc():
    orb_dir = '../../../tests/PP_ORB/'
    file_list = ["C_gga_8au_100Ry_2s2p1d.orb", "H_gga_8au_60Ry_2s1p.orb", "O_gga_10au_100Ry_2s2p1d.orb", "Fe_gga_9au_100Ry_4s2p2d1f.orb"]
    file_list = [orb_dir + orbfile for orbfile in file_list]
    nfile = len(file_list)

    orb = nao.RadialCollection()
    orb.build(nfile, file_list, 'o')

    assert orb.symbol(0) == 'C'
    assert orb.symbol(1) == 'H'
    assert orb.symbol(2) == 'O'
    assert orb.symbol(3) == 'Fe'

    assert orb.ntype == 4
    assert orb.lmax == 3
    assert orb.rcut_max == 10.0

    assert orb.nzeta(0,0) == 2
    assert orb.nzeta(1,0) == 2
    assert orb.nzeta(2,0) == 2
    assert orb.nzeta(3,0) == 4


    assert orb.nzeta_max(0) == 2
    assert orb.nzeta_max(3) == 4

    assert orb.nchi(0) == 5
    assert orb.nchi(3) == 9
    assert orb.nchi() == 22

    # not support index RadialSet now
    # support index NumericalRadial test
    assert orb(0,0,0).l == 0
    assert orb(0,1,0).l == 1
    assert orb(0,1,1).l == 1
    assert orb(0,2,0).l == 2






    


