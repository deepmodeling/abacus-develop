from __future__ import annotations

import pytest
from pyabacus import ModuleBase as base
import numpy as np



def test_sphbes():
    s = base.Sphbes()
    # test for sphbesj
    assert s.sphbesj(1, 0.0) == 0.0
    assert s.sphbesj(0, 0.0) == 1.0
    assert s.sphbesj(1, np.array([0.0]), 1, 1, np.zeros(1)) == None

def test_static_array_wrappers():
    zeros = np.zeros(2)
    base.Sphbes.sphbes_zeros(0, 2, zeros)
    np.testing.assert_allclose(zeros, [np.pi, 2 * np.pi])

    mesh = 3
    func = np.ones(mesh)
    rab = np.ones(mesh)
    integral_from_zero = np.zeros(mesh)
    integral_to_infinity = np.zeros(mesh)

    base.Integral.Simpson_Integral_0toall(mesh, func, rab, integral_from_zero)
    base.Integral.Simpson_Integral_alltoinf(mesh, func, rab, integral_to_infinity)

    np.testing.assert_allclose(integral_from_zero, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(integral_to_infinity, [2.0, 1.0, 0.0])

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
    


