from pyabacus import hsolver
import numpy as np
import scipy

h_mat = scipy.io.loadmat('./Si2.mat')['Problem']['A'][0, 0]

nbasis = h_mat.shape[0]
nband = 8

v0 = np.random.rand(nbasis, nband)

diag_elem = h_mat.diagonal()
diag_elem = np.where(np.abs(diag_elem) < 1e-8, 1e-8, diag_elem)
precond = 1.0 / np.abs(diag_elem)


def mm_op(x):
    return h_mat.dot(x)

e, v = hsolver.dav_subspace(
    mm_op,
    v0,
    nbasis,
    nband,
    precond,
    dav_ndim=8,
    tol=1e-8,
    max_iter=1000,
    scf_type=False
)

print('eigenvalues calculated by pyabacus: ', e)

e_scipy, v_scipy = scipy.sparse.linalg.eigsh(h_mat, k=nband, which='SA', maxiter=1000)
e_scipy = np.sort(e_scipy)
print('eigenvalues calculated by scipy: ', e_scipy)

print('error:', e - e_scipy)