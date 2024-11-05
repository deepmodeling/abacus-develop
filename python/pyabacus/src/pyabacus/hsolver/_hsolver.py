"""
pyabacus.hsolver
==================
Module for Solver for Hamiltonian in ABACUS
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List, Union, Callable

from ._hsolver_pack import diag_comm_info as _diag_comm_info
from ._hsolver_pack import diago_dav_subspace, diago_david, diago_cg

class diag_comm_info(_diag_comm_info):
    def __init__(self, rank: int, nproc: int):
        super().__init__(rank, nproc)
    
    @property
    def rank(self) -> int:
        return super().rank
    
    @property
    def nproc(self) -> int:
        return super().nproc

def dav_subspace(
    mvv_op: Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
    init_v: NDArray[np.complex128],
    dim: int,
    num_eigs: int,
    precondition: NDArray[np.float64],
    dav_ndim: int = 2,
    tol: float = 1e-2,
    max_iter: int = 1000,
    need_subspace: bool = False,
    diag_ethr: Union[List[float], None] = None,
    scf_type: bool = False
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """ A function to diagonalize a matrix using the Davidson-Subspace method.

    Parameters
    ----------
    mvv_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
        The operator to be diagonalized, which is a function that takes a set of 
        vectors X = [x1, ..., xN] as input and returns a matrix(vector block)
        mvv_op(X) = H * X ([Hx1, ..., HxN]) as output.
    init_v : NDArray[np.complex128]
        The initial guess for the eigenvectors.
    dim : int
        The number of basis, i.e. the number of rows/columns in the matrix.
    num_eigs : int
        The number of bands to calculate, i.e. the number of eigenvalues to calculate.
    precondition : NDArray[np.float64]
        The preconditioner.
    dav_ndim : int, optional
        The number of vectors in the subspace, by default 2.
    tol : float, optional
        The tolerance for the convergence, by default 1e-2.
    max_iter : int, optional    
        The maximum number of iterations, by default 1000.
    need_subspace : bool, optional
        Whether to use subspace function, by default False.
    diag_ethr : List[float] | None, optional
        The list of thresholds of bands, by default None.
    scf_type : bool, optional
        Indicates whether the calculation is a self-consistent field (SCF) calculation. 
        If True, the initial precision of eigenvalue calculation can be coarse. 
        If False, it indicates a non-self-consistent field (non-SCF) calculation, 
        where high precision in eigenvalue calculation is required from the start.  
    
    Returns
    -------
    e : NDArray[np.float64]
        The eigenvalues.
    v : NDArray[np.complex128]
        The eigenvectors corresponding to the eigenvalues.
    """
    if not callable(mvv_op):
        raise TypeError("mvv_op must be a callable object.")
    
    if diag_ethr is None:
        diag_ethr = [tol] * num_eigs
    
    if init_v.ndim != 1 or init_v.dtype != np.complex128:
        init_v = init_v.flatten().astype(np.complex128, order='C')
    
    _diago_obj_dav_subspace = diago_dav_subspace(dim, num_eigs)
    _diago_obj_dav_subspace.set_psi(init_v)
    _diago_obj_dav_subspace.init_eigenvalue()
    
    comm_info = diag_comm_info(0, 1)
    assert dav_ndim > 1, "dav_ndim must be greater than 1."
    assert dav_ndim * num_eigs < dim * comm_info.nproc, "dav_ndim * num_eigs must be less than dim * comm_info.nproc."
   
    _ = _diago_obj_dav_subspace.diag(
        mvv_op,
        precondition,
        dav_ndim,
        tol,
        max_iter,
        need_subspace,
        diag_ethr,
        scf_type,
        comm_info
    )
    
    e = _diago_obj_dav_subspace.get_eigenvalue()
    v = _diago_obj_dav_subspace.get_psi()
    
    return e, v

def davidson(
    mvv_op: Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
    init_v: NDArray[np.complex128],
    dim: int,
    num_eigs: int,
    precondition: NDArray[np.float64],
    dav_ndim: int = 2,
    tol: float = 1e-2,
    max_iter: int = 1000,
    use_paw: bool = False,
    # scf_type: bool = False
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """ A function to diagonalize a matrix using the Davidson-Subspace method.

    Parameters
    ----------
    mvv_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
        The operator to be diagonalized, which is a function that takes a set of 
        vectors X = [x1, ..., xN] as input and returns a matrix(vector block)
        mvv_op(X) = H * X ([Hx1, ..., HxN]) as output.
    init_v : NDArray[np.complex128]
        The initial guess for the eigenvectors.
    dim : int
        The number of basis, i.e. the number of rows/columns in the matrix.
    num_eigs : int
        The number of bands to calculate, i.e. the number of eigenvalues to calculate.
    precondition : NDArray[np.float64]
        The preconditioner.
    dav_ndim : int, optional
        The number of vectors in the subspace, by default 2.
    tol : float, optional
        The tolerance for the convergence, by default 1e-2.
    max_iter : int, optional    
        The maximum number of iterations, by default 1000.
    use_paw : bool, optional
        Whether to use projector augmented wave (PAW) method, by default False.
    
    Returns
    -------
    e : NDArray[np.float64]
        The eigenvalues.
    v : NDArray[np.complex128]
        The eigenvectors corresponding to the eigenvalues.
    """
    if not callable(mvv_op):
        raise TypeError("mvv_op must be a callable object.")
    
    if init_v.ndim != 1 or init_v.dtype != np.complex128:
        init_v = init_v.flatten().astype(np.complex128, order='C')
    
    _diago_obj_david = diago_david(dim, num_eigs)
    _diago_obj_david.set_psi(init_v)
    _diago_obj_david.init_eigenvalue()
    
    comm_info = diag_comm_info(0, 1)
    
    _ = _diago_obj_david.diag(
        mvv_op,
        precondition,
        dav_ndim,
        tol,
        max_iter,
        use_paw,
        comm_info
    )
    
    e = _diago_obj_david.get_eigenvalue()
    v = _diago_obj_david.get_psi()
    
    return e, v
    
def cg(
    mvv_op: Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
    init_v: NDArray[np.complex128],
    dim: int,
    num_eigs: int,
    precondition: NDArray[np.float64],
    tol: float = 1e-2,
    max_iter: int = 1000,
    need_subspace: bool = False,
    scf_type: bool = False,
    nproc_in_pool: int = 1
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """ A function to diagonalize a matrix using the Conjugate Gradient method.
    
    Parameters
    ----------
    mvv_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
        The operator to be diagonalized, which is a function that takes a set of 
        vectors X = [x1, ..., xN] as input and returns a matrix(vector block)
        mvv_op(X) = H * X ([Hx1, ..., HxN]) as output.
    init_v : NDArray[np.complex128]
        The initial guess for the eigenvectors.
    dim : int
        The number of basis, i.e. the number of rows/columns in the matrix.
    num_eigs : int
        The number of bands to calculate, i.e. the number of eigenvalues to calculate.
    precondition : NDArray[np.float64]
        The preconditioner.
    max_iter : int, optional
        The maximum number of iterations, by default 1000.
    tol : float, optional
        The tolerance for the convergence, by default 1e-2.
    need_subspace : bool, optional
        Whether to use subspace function, by default False.
    scf_type : bool, optional   
        Indicates whether the calculation is a self-consistent field (SCF) calculation. 
        If True, the initial precision of eigenvalue calculation can be coarse. 
        If False, it indicates a non-self-consistent field (non-SCF) calculation, 
        where high precision in eigenvalue calculation is required from the start.
    nproc_in_pool : int, optional
        The number of processes in the pool, by default 1.
        
    Returns
    -------
    e : NDArray[np.float64]
        The eigenvalues.
    v : NDArray[np.complex128]
        The eigenvectors corresponding to the eigenvalues.
    """
    if not callable(mvv_op):
        raise TypeError("mvv_op must be a callable object.")
    
    if init_v.ndim != 1 or init_v.dtype != np.complex128:
        # the shape of init_v is (num_eigs, dim) = (dim, num_eigs).T
        if init_v.ndim == 2:
            init_v = init_v.T
        init_v = init_v.flatten().astype(np.complex128, order='C')
    
    _diago_obj_cg = diago_cg(dim, num_eigs)
    _diago_obj_cg.set_psi(init_v)
    _diago_obj_cg.init_eig()
    
    _diago_obj_cg.set_prec(precondition)
    
    _diago_obj_cg.diag(
        mvv_op,
        max_iter, 
        tol,
        need_subspace,
        scf_type,
        nproc_in_pool
    )
    
    e = _diago_obj_cg.get_eig()
    v = _diago_obj_cg.get_psi()
    
    return e, v