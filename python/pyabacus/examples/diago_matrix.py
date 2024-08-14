from pyabacus import hsolver
import numpy as np
import scipy

h_sparse = scipy.io.loadmat('./Si2.mat')['Problem']['A'][0, 0]
h_mat = h_sparse.toarray()

nbasis = h_mat.shape[0]
nband = 16

v0 = np.random.rand(nbasis, nband)

diag_elem = np.diag(h_mat.reshape(nbasis, nbasis))
diag_elem = np.where(np.abs(diag_elem) < 1e-5, 1e-5, diag_elem)
precond = 1.0 / np.abs(diag_elem)

e, v = hsolver.dav_subspace(
    h_mat,
    v0,
    nbasis,
    nband,
    precond,
    dav_ndim=8,
    tol=1e-4,
    max_iter=1000,
    scf_type=False
)

print(e)

e_scipy, v_scipy = scipy.sparse.linalg.eigsh(h_sparse, k=nband, which='SM', maxiter=1000)
e_scipy = np.sort(e_scipy)
print(e_scipy)

print('error:', e - e_scipy)