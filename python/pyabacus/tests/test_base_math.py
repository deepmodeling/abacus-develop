from __future__ import annotations

import pyabacus as m
import numpy as np


def test_version():
    assert m.__version__ == "0.0.1"

def test_sphbes():
    s = m.ModuleBase.Sphbes()
    # test for sphbesj
    assert s.sphbesj(1, 0.0) == 0.0
    assert s.sphbesj(0, 0.0) == 1.0

