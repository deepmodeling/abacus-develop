# pyabacus.hsolver

import numpy as np
from numpy.typing import NDArray
from typing import List

class diag_comm_info:
    def __init__(self, rank: int, nproc: int) -> None: ...
    
    @property
    def rank(self) -> int: ...
    
    @property
    def nproc(self) -> int: ...
    
class Diago_DavSubspace:
    def __init__(
        self, 
        pre_condition: NDArray[np.float64] | List[float], 
        nband: int, 
        nbasis: int, 
        dav_ndim: int, 
        diag_thr: float, 
        diag_nmax: int, 
        need_subspace: bool, 
        comm_info: diag_comm_info
    ) -> None: ...
    
    def diag(
        self, 
        h_mat: NDArray[np.complex128], 
        psi: NDArray[np.complex128], 
        nbasis: int, 
        eigenvalues: NDArray[np.float64], 
        is_occupied: list[bool], 
        scf_type: bool
    ) -> int: ...