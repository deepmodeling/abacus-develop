import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
class toQO_BasisFilter:
    """Balance between flexible and accurate basis selection, filter basis according to practical band structure
    matrix should be k-point resoluted, i.e., matrix[ik] is the matrix for kpoint ik
    Raises:
        ValueError: if band range is not properly configured
    """
    band_range: tuple # min and max band index
    overlap_matrix_filter_threshold: float # filter thresold for overlap matrix, only larger than this threshold will basis retain 
    denstiy_matrix_filter_threshold: float # filter thresold for density matrix, only larger than this threshold will basis retain

    def __init__(self) -> None:
        pass

    def set_overlap_filter(self, overlap_matrix_filter_threshold: float) -> None:
        """set overlap matrix filter threshold

        Args:
            overlap_matrix_filter_threshold (float): overlap matrix filter threshold
        """
        self.overlap_matrix_filter_threshold = overlap_matrix_filter_threshold
    
    def filter_via_overlap_matrix(self, sk: np.ndarray, sqok: np.ndarray) -> np.ndarray:
        """calculate AO overlap matrix in NAO representation. Those diagonal elements that are failed to be larger than threshold will be filtered out
        because they cannot be represented by NAOs
        
        Args:
            nkpts (int): number of kpoints
            sk (np.ndarray): overlap matrix in k-space
            sqok (np.ndarray): AO-NAO overlap matrix
        
        Returns:
            deleted_basis_indices (list): list of list, each list contains indices of selected basis for each kpoint
        """
        selected_basis_indices = []
        _sqok = np.zeros_like(sqok)
        
        _x = la.solve(sk, sqok.conj().T)
        _ovlp = sqok @ _x

        for i in range(_ovlp.shape[0]):
            if _ovlp[i, i] > self.overlap_matrix_filter_threshold:
                selected_basis_indices.append(i)
        # filter out basis
        _sqok = sqok[selected_basis_indices, :]

        return selected_basis_indices, _sqok
