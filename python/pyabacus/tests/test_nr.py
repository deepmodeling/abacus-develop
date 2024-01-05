from __future__ import annotations

import pyabacus as m


def test_version():
    assert m.__version__ == "0.0.1"

def test_attributes():
    chi = m.NumericalRadial()
    # string
    assert chi.symbol == ''
    # integer
    assert chi.itype == 0
    assert chi.izeta == 0
    assert chi.l == -1
    assert chi.nr == 0
    assert chi.nk == 0
    # float
    assert chi.rcut == 0.0
    assert chi.kcut == 0.0
    assert chi.pr == 0.0
    assert chi.pk == 0.0
    # bool
    assert chi.is_fft_compliant == False
