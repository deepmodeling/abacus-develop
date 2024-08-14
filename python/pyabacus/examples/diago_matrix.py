from pyabacus import hsolver
import numpy as np
import scipy

nband = 25
nbasis = 25

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

diag_elem = np.diag(h_mat.reshape(nbasis, nbasis))
diag_elem = np.where(np.abs(diag_elem) < 1e-5, 1e-5, diag_elem)
precond = 1.0 / np.abs(diag_elem)

e, v = hsolver.dav_subspace(
    h_mat,
    h_mat[:nbasis*nband],
    nbasis,
    nband,
    precond,
    dav_ndim=4,
    tol=1e-2,
    max_iter=100000,
    scf_type=True
)

print(e)
# print(v[-nbasis:])

a = h_mat.reshape(nbasis, nbasis)

e_scipy, v_scipy = scipy.linalg.eigh(a)
print(e_scipy)