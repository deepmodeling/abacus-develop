from pyabacus import hsolver
import numpy as np
import scipy

from pyscf.lib import linalg_helper

nband = 25
nbasis = 25

pre_condition = np.ones(nbasis, dtype=np.float64, order='C')

psi = np.random.rand(nbasis * nband) + 1j * np.random.rand(nbasis * nband)
psi = psi.astype(np.complex128, order='C')

h_mat = np.zeros(nbasis * nbasis, dtype=np.complex128, order='C')

for i in range(nbasis):
    h_mat[i * nbasis + i] = np.random.rand() * 1000 + 0.0j

for i in range(1, nbasis):
    for k in range(i):
        value = np.random.rand() + np.random.rand() * 1.0j
        h_mat[i * nbasis + k] = value
        h_mat[k * nbasis + i] = np.conj(value)

# for i in range(nbasis):
#     h_mat[i * nbasis + i] = 1.0
    
# for i in range(1, nbasis):
#     for k in range(i):
#         value = i + k * 1.0j
#         h_mat[i * nbasis + k] = value
#         h_mat[k * nbasis + i] = np.conj(value)

e, v = hsolver.dav_subspace(
    h_mat,
    h_mat[:nbasis*nband],
    nbasis,
    nband,
    pre_condition,
    dav_ndim=4,
    tol=1e-2,
    max_iter=100000,
    scf_type=True
)



print(e / (nband*nbasis))
# print(v[-nbasis:])

a = h_mat.reshape(nbasis, nbasis)

e_scipy, v_scipy = scipy.linalg.eigh(a)
print(e_scipy)