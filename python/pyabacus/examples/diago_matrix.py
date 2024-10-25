from pyabacus import hsolver
import numpy as np
import scipy

def load_mat(mat_file):
    h_mat = scipy.io.loadmat(mat_file)['Problem']['A'][0, 0]
    nbasis = h_mat.shape[0]
    nband = 8
    
    return h_mat, nbasis, nband

def calc_eig_pyabacus(mat_file, method):
    algo = {
        'dav_subspace': hsolver.dav_subspace,
        'davidson': hsolver.davidson
    }
    
    h_mat, nbasis, nband = load_mat(mat_file)
    
    v0 = np.random.rand(nbasis, nband)
    diag_elem = h_mat.diagonal()
    diag_elem = np.where(np.abs(diag_elem) < 1e-8, 1e-8, diag_elem)
    precond = 1.0 / np.abs(diag_elem)

    def mm_op(x):
        return h_mat.dot(x)

    e, _ = algo[method](
        mm_op,
        v0,
        nbasis,
        nband,
        precond,
        dav_ndim=8,
        tol=1e-8,
        max_iter=1000
    )

    print(f'eigenvalues calculated by pyabacus-{method} is: \n', e)
    
    return e

def calc_eig_scipy(mat_file):
    h_mat, _, nband = load_mat(mat_file)
    e, _ = scipy.sparse.linalg.eigsh(h_mat, k=nband, which='SA', maxiter=1000)
    e = np.sort(e)
    print('eigenvalues calculated by scipy is: \n', e)
    
    return e

if __name__ == '__main__':
    mat_file = './Si2.mat'
    method = ['dav_subspace', 'davidson']
    
    for m in method:
        print(f'\n====== Calculating eigenvalues using {m} method... ======')
        e_pyabacus = calc_eig_pyabacus(mat_file, m)
        e_scipy = calc_eig_scipy(mat_file)
        
        print('eigenvalues difference: \n', e_pyabacus - e_scipy)
    
    