from pyabacus import hsolver
import numpy as np
import scipy

def load_mat(mat_file):
    h_mat = scipy.io.loadmat(mat_file)['Problem']['A'][0, 0]
    nbasis = h_mat.shape[0]
    nband = 8
    
    return h_mat, nbasis, nband

def calc_eig_pyabacus(mat_file, method):
    dav = {
        'dav_subspace': hsolver.dav_subspace,
        'davidson': hsolver.davidson
    }
    cg = {
        'cg': hsolver.cg
    }
    
    h_mat, nbasis, nband = load_mat(mat_file)
    
    v0 = np.random.rand(nbasis, nband)
    diag_elem = h_mat.diagonal()
    diag_elem = np.where(np.abs(diag_elem) < 1e-8, 1e-8, diag_elem)
    precond = 1.0 / np.abs(diag_elem)

    def mvv_op(x):
        return h_mat.dot(x)

    if method in dav:
        algo = dav[method]
        # args: mvvop, init_v, dim, num_eigs, precondition, dav_ndim, tol, max_iter
        args = (mvv_op, v0, nbasis, nband, precond, 8, 1e-12, 5000)
    elif method in cg:
        algo = cg[method]
        # args: mvvop, init_v, dim, num_eigs, precondition, tol, max_iter
        args = (mvv_op, v0, nbasis, nband, precond, 1e-12, 5000)
    else:
        raise ValueError(f"Method {method} not available")
    
    e, _ = algo(*args)

    print(f'eigenvalues calculated by pyabacus-{method} is: \n', e)
    
    return e

def calc_eigsh(mat_file):
    h_mat, _, nband = load_mat(mat_file)
    e, _ = scipy.sparse.linalg.eigsh(h_mat, k=nband, which='SA', maxiter=1000)
    e = np.sort(e)
    print('eigenvalues calculated by scipy is: \n', e)
    
    return e

if __name__ == '__main__':
    mat_file = './Si2.mat'
    method = ['dav_subspace', 'davidson', 'cg']
    
    for m in method:
        print(f'\n====== Calculating eigenvalues using {m} method... ======')
        e_pyabacus = calc_eig_pyabacus(mat_file, m)
        e_scipy = calc_eigsh(mat_file)
        
        print('eigenvalues difference: \n', e_pyabacus - e_scipy)
        