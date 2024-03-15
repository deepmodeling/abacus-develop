from __future__ import annotations

import pytest
from pyabacus import ModuleBase as base
import numpy as np



def test_sphbes():
    s = base.Sphbes()
    # test for sphbesj
    assert s.sphbesj(1, 0.0) == 0.0
    assert s.sphbesj(0, 0.0) == 1.0

def test_sbt():
    sbt = base.SphericalBesselTransformer()

@pytest.fixture
def simpson_setup():
    n = 1000
    x = np.linspace(0, 2*np.pi, n)
    func = np.sin(x)
    dx = 2 * np.pi / (n - 1)
    return n, func, dx

def test_simpson(simpson_setup):
    n, func, dx= simpson_setup
    s = base.Integral()
    assert s.simpson(n, func, dx) == pytest.approx(0, abs=1e-10)
    


