from __future__ import annotations

import pytest
from pyabacus import hsolver
import numpy as np
import scipy

def diag_pyabacus(h_sparse, nband):
    def mm_op(x):
        return h_sparse.dot(x)

    nbasis = h_sparse.shape[0]

    v0 = np.random.rand(nbasis, nband)

    diag_elem = h_sparse.diagonal()
    diag_elem = np.where(np.abs(diag_elem) < 1e-5, 1e-5, diag_elem)
    precond = 1.0 / np.abs(diag_elem)

    e, _ = hsolver.dav_subspace(
        mm_op,
        v0,
        nbasis,
        nband,
        precond,
        dav_ndim=8,
        tol=1e-4,
        max_iter=5000,
        scf_type=False
    )
    
    return e

def diag_eigsh(h_sparse, nband):
    e, _ = scipy.sparse.linalg.eigsh(h_sparse, k=nband, which='SA', maxiter=5000, tol=1e-4)
    return e

@pytest.mark.parametrize("file_name, nband", [
    ('./test_diag/Si2.mat', 16),
    ('./test_diag/SiNa.mat', 16),
    ('./test_diag/Na5.mat', 16),
    ('./test_diag/H2O.mat', 16)
])
def test_diag(file_name, nband):
    h_sparse = scipy.io.loadmat(file_name)['Problem']['A'][0, 0]
    e_pyabacus = diag_pyabacus(h_sparse, nband)
    e_scipy = diag_eigsh(h_sparse, nband)
    np.testing.assert_allclose(e_pyabacus, e_scipy, atol=1e-2)