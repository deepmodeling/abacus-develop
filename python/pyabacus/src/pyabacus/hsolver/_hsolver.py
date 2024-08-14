# pyabacus.hsolver

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List

from .._core import hsolver

class diag_comm_info:
    def __init__(self, rank: int, nproc: int) -> None: ...
    
    @property
    def rank(self) -> int: ...
    
    @property
    def nproc(self) -> int: ...
    
def dav_subspace(
    h_mat: NDArray[np.complex128],
    init_v: NDArray[np.complex128],
    nbasis: int,
    nband: int,
    pre_condition: NDArray[np.float64],
    dav_ndim: int = 2,
    tol: float = 1e-2,
    max_iter: int = 1000,
    need_subspace: bool = False,
    is_occupied: List[bool] | None = None,
    scf_type: bool = False
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    if is_occupied is None:
        is_occupied = [True] * nband
    
    _diago_obj_dav_subspace = hsolver.diago_dav_subspace(nbasis, nband)
    _diago_obj_dav_subspace.set_psi(init_v)
    _diago_obj_dav_subspace.init_eigenvalue()
    
    comm_info = hsolver.diag_comm_info(0, 1)
    res = _diago_obj_dav_subspace.diag(
        h_mat,
        pre_condition,
        dav_ndim,
        tol,
        max_iter,
        need_subspace,
        is_occupied,
        scf_type,
        comm_info
    )
    
    e = _diago_obj_dav_subspace.get_eigenvalue()
    v = _diago_obj_dav_subspace.get_psi()
    
    return e, v
    