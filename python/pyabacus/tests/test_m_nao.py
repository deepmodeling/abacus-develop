from __future__ import annotations

import pytest
from pyabacus import ModuleNAO as nao
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
    assert(chi.pr == pr)
    assert(chi.pk == 0)
    assert(chi.is_fft_compliant == False)
