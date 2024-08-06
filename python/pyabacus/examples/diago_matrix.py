from pyabacus import hsolver
import numpy as np

h_mat = np.array(
    [
        4.0+0.0j, 2.0+0.0j, 2.0+0.0j,
        2.0+0.0j, 4.0+0.0j, 2.0+0.0j,
        2.0+0.0j, 2.0+0.0j, 4.0+0.0j
    ], 
dtype=np.complex128, order='C')
# h_mat = np.array(
#     [
#         1.0+0.0j, 0.0+0.0j, 0.0+0.0j,
#         0.0+0.0j, 1.0+0.0j, 0.0+0.0j,
#         0.0+0.0j, 0.0+0.0j, 1.0+0.0j
#     ],
# dtype=np.complex128, order='C')
pre_condition = np.ones(3, dtype=np.float64, order='C')
nband = 1
nbasis = 3
dav_ndim = 2
diag_thr = 1e-2
diag_nmax = 1000
need_subspace = False
comm_info = hsolver.diag_comm_info(0, 1)

psi = np.ones(nbasis * nband, dtype=np.complex128, order='C')
eigenvalues = np.zeros(nband, dtype=np.float64, order='C')
is_occupied = [True] * nband

diago_dav_subspace = hsolver.Diago_DavSubspace(
    pre_condition,
    nband,
    nbasis,
    dav_ndim,
    diag_thr,
    diag_nmax,
    need_subspace,
    comm_info
)

res = diago_dav_subspace.diag(
    h_mat,
    psi,
    nbasis,
    eigenvalues,
    is_occupied,
    False
)

print(f'res: {res}')
print(f'eigenvalues: {eigenvalues}')