from pyabacus import hsolver
import numpy as np
import scipy

def load_mat(mat_file):
    h_mat = scipy.io.loadmat(mat_file)['Problem']['A'][0, 0]
    nbasis = h_mat.shape[0]
    nband = 8
    
    return h_mat, nbasis, nband

def gen_dense_mat(dim):
    # generate a random symmetric and positive definite matrix
    h_mat = np.random.rand(dim, dim)
    h_mat = h_mat + h_mat.T
    h_mat = h_mat + dim * np.eye(dim)
    
    return h_mat

def calc_eig_dav(mat_file, method):
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
        tol=1e-8
    )

    print(f'eigenvalues calculated by pyabacus-{method} is: \n', e)
    
    return e

def calc_eig_cg(h_mat, num_eigs):
    dim = h_mat.shape[0]
    v0 = np.random.rand(dim, num_eigs)
    diag_elem = h_mat.diagonal()
    diag_elem = np.where(np.abs(diag_elem) < 1e-8, 1e-8, diag_elem)
    precond = 1.0 / np.abs(diag_elem)
    
    def mm_op(x):
        return h_mat.dot(x)
    
    e, _ = hsolver.cg(
        mm_op,
        v0,
        dim,
        num_eigs,
        precond,
        tol=1e-8
    )
    
    print('eigenvalues calculated by pyabacus-cg is: \n', e)
    
    return e

def calc_eigsh(mat_file):
    h_mat, _, nband = load_mat(mat_file)
    e, _ = scipy.sparse.linalg.eigsh(h_mat, k=nband, which='SA', maxiter=1000)
    e = np.sort(e)
    print('eigenvalues calculated by scipy is: \n', e)
    
    return e

def calc_eigh(h_mat, num_eigs):
    e, _ = scipy.linalg.eigh(h_mat)
    e = np.sort(e)
    print('eigenvalues calculated by scipy is: \n', e[:num_eigs])
    
    return e

if __name__ == '__main__':
    mat_file = './Si2.mat'
    method = ['dav_subspace', 'davidson']
    
    for m in method:
        print(f'\n====== Calculating eigenvalues using {m} method... ======')
        e_pyabacus = calc_eig_dav(mat_file, m)
        e_scipy = calc_eigsh(mat_file)
        
        print('eigenvalues difference: \n', e_pyabacus - e_scipy)
        
    print("\n====== davidson and dav_subspace method Done! ======")
    print("\n====== CG method ======")
    
    h_mat = gen_dense_mat(100)
    num_eigs = 8
    e_cg = calc_eig_cg(h_mat, num_eigs)
    e_scipy = calc_eigh(h_mat, num_eigs)
    
    print('eigenvalues difference: \n', e_cg - e_scipy[:num_eigs])
    print("\n====== CG method Done! ======")