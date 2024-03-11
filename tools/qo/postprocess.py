"""
----------------------------------------------------------------------------------
|               ABACUS LCAO2QO(Quasiatomic Orbital) Instruction                  |
----------------------------------------------------------------------------------
| Author: Yike HUANG
| Email: huangyk@aisi.ac.cn
| Date: 2024-03-02
|
| TO USE:
| Can either use this script in-place or by importing it as a module.
| 0. Run ABACUS:
|     - Run ABACUS with the following settings:
|         - `calculation scf`
|         - `qo_switch 1`
|         - `qo_basis hydrogen` # can be hydrogen, pswfc and nao
|         - `qo_strategy energy-valence` # for more strategies, see ABACUS manual
|         - `qo_thr 1e-10` # threshold for determing range of integral in realspace
|       If want to turn on Slater Screening to mimic many-electron behavior (sort
|       of), add the following settings:
|         - `qo_screening_coeff 0.1`
|       More details please refer to ABACUS manual.
| 1. In-place:
|     - Put this script in the same directory as the output files of ABACUS.
|     - Run the script with the command `python all_in_one.py` after modifying
|       the parameters in the `if __name__ == "__main__":` block.
| 2. As a module:
|     - Import this script as a module in another script.
|     - Call the function `production` with the appropriate parameters, the 
|       function `production` returns the Hamiltonian and overlap matrices in 
|       QO representation along with valid kpoints direct coordinates.
| 9999. Unittest:
|     - Run the script with the command `python all_in_one.py` to perform
|       unittests after modifying the parameters in the `if __name__ == "__main__":`
|       block. Unittest will read data stored in the same folder when performing
|       production runs, so it is necessary targeting the correct path first.
----------------------------------------------------------------------------------
"""

import numpy as np
# TESTED
def cxx_topycomplex(num: str):
    """cxx prints complex numbers in the form (a,b)"""
    num = num.replace('(', '').replace(')', '')
    num = num.split(',')
    return complex(float(num[0]), float(num[1]))
# TESTED
def read_mat_hs(fhs):
    """"""
    with open(fhs, 'r') as f:
        data = "".join(f.readlines()).replace('\n', '').split()
    size = int(data[0])
    indexing = []
    for i in range(size):
        indexing += [(i, j) for j in range(i, size)]
    data = data[1:]
    mat = np.zeros((size, size), dtype=np.complex128)
    for inum, number in enumerate(data):
        mat[indexing[inum]] = cxx_topycomplex(number)
    mat = mat + mat.conj().T
    for i in range(size):
        mat[i, i] = mat[i, i] / 2
    return mat
# TESTED
def read_lowf(flowf):
    """read the lowf matrix in k-space, and the k-vector.
    
    Args:
        flowf (str): path to the file containing the lowf matrix
    
    Returns:
        np.ndarray: lowf matrix in k-space
        np.ndarray: k-vector
    """
    with open(flowf, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    for line in lines:
        if line.endswith("(number of bands)"):
            nband = int(line.split()[0])
        elif line.endswith("(number of orbitals)"):
            nlocal = int(line.split()[0])
        else:
            continue
    lowf_k = np.zeros((nlocal, nband), dtype=np.complex128)
    kvec_c = None
    ib, ilocal = 0, 0
    for line in lines:
        if line.endswith("(band)"):
            ib = int(line.split()[0]) - 1
            ilocal = 0
            continue
        if not line.endswith(")"):
            if line.count(" ") != 2:
                # it means it is a band rather than coordinates of kpoint
                nums = line.split()
                c = [complex(float(nums[i]), float(nums[i+1])) for i in range(0, len(nums), 2)]
                for i in range(len(c)):
                    lowf_k[ilocal, ib] = c[i]
                    ilocal += 1
            else:
                kvec_c = np.array([float(x) for x in line.split()])
    return lowf_k, kvec_c
# TESTED
def read_ao_proj(fao_proj):
    """Read the atomic orbital projection matrix in k-space, and the k-vector.
    
    Args:
        fao_proj (str): path to the file containing the atomic orbital projection matrix
    
    Returns:
        np.ndarray: atomic orbital projection matrix in k-space
        np.ndarray: k-vector
    """
    with open(fao_proj, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    kvec_d = None
    if lines[0].startswith("KPOINT_COORDINATE"):
        kvec_d = np.array([float(x) for x in lines[0].split()[-3:]])
        lines = lines[1:]
    nlocal = len(lines[0].split())
    nao = len(lines)
    ao_proj = np.zeros((nao, nlocal), dtype=np.complex128)
    for i, line in enumerate(lines):
        nums = [cxx_topycomplex(num) for num in line.split()]
        for j in range(nlocal):
            ao_proj[i, j] = nums[j]
    # because in ao_proj, the row index is AO, we need to dagger it
    return ao_proj.conj().T, kvec_d
# TESTED
import os
def read(path: str, 
         nk: int,
         read_H: bool = True,
         read_S: bool = True,
         read_LOWF: bool = True,
         read_QO: bool = True) -> tuple[list, list, list, list, list]:
    """Read H(k), S(k), LOWF(k), and <\phi(k)|A(k)> from ABACUS output files.
    
    Args:
        path (str): path to the directory containing the output files
        nk (int): number of kpoints
        read_H (bool, optional): whether to read H(k). Defaults to True.
        read_S (bool, optional): whether to read S(k). Defaults to True.
        read_LOWF (bool, optional): whether to read LOWF(k). Defaults to True.
        read_QO (bool, optional): whether to read <\phi(k)|A(k)>. Defaults to True.
    
    Returns:
        tuple[list, list, list, list, list]: H(k), S(k), LOWF(k), <\phi(k)|A(k)>, k-vectors
    """
    hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_c, kvecs_d = [], [], [], [], [], []
    if read_H:
        for ik in range(nk):
            fh = os.path.join(path, f"data-{ik}-H")
            hs_k.append(read_mat_hs(fh))
    if read_S:
        for ik in range(nk):
            fs = os.path.join(path, f"data-{ik}-S")
            ss_k.append(read_mat_hs(fs))
    if read_LOWF:
        for ik in range(nk):
            flowf = os.path.join(path, f"LOWF_K_{ik+1}.txt")
            lowf_k, kvec_c = read_lowf(flowf)
            lowfs_k.append(lowf_k)
            kvecs_c.append(kvec_c)
    if read_QO:
        for ik in range(nk):
            fao_proj = os.path.join(path, f"QO_ovlp_{ik}.dat")
            ao_proj_k, kvec_d = read_ao_proj(fao_proj)
            aos_proj_k.append(ao_proj_k)
            kvecs_d.append(kvec_d)

    return hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_d
# TESTED
def cal_denmat(lowf_k: np.ndarray) -> np.ndarray:
    """calculate the density matrix from the lowf matrix.
    
    Args:
        lowf_k (np.ndarray): lowf matrix in k-space
    
    Returns:
        np.ndarray: density matrix in k-space
    """
    return lowf_k@lowf_k.conj().T
# TESTED
import scipy as sp
def lowdin_onW(s_k: np.ndarray, denmat_k: np.ndarray, ao_proj_k: np.ndarray, m: int):
    """Perform Lowdin orthogonalization on the W matrix, which is defined as W = A^+ * (S^-1 - D) * A.
    
    Args:
        s_k (np.ndarray): overlap matrix in k-space
        denmat_k (np.ndarray): density matrix in k-space
        ao_proj_k (np.ndarray): atomic orbital projection matrix in k-space
        m (int): number of bands wanted
        
    Returns:
        np.ndarray: the lowdin orthogonalized matrix in the quasi-atomic orbital basis, can check
        orthonormality with overlap matrix
    """
    nrow = s_k.shape[0]
    ncol = s_k.shape[1]
    assert nrow == ncol
    sinv_k = np.linalg.solve(s_k, np.eye(N=nrow, M=ncol, dtype=np.complex128))
    W_mat = ao_proj_k.conj().T@(sinv_k - denmat_k)@ao_proj_k
    # forcely let diagonal elements to be real
    W_mat = (W_mat + W_mat.conj().T) / 2
    # diagonalize W and get only m largest eigenvalues and corresponding eigenvectors
    eigvals, eigvecs = sp.linalg.eigh(W_mat)
    eigvals = eigvals[-m:]
    eigvecs = eigvecs[:, -m:]
    lambdas = np.diag(1.0/np.sqrt(eigvals))
    
    return (sinv_k - denmat_k)@ao_proj_k@eigvecs@lambdas
    # call return as lowf_k_bar
# TESTED
def cal_tb_hs(h_k: np.ndarray, s_k: np.ndarray, denmat_k_tilde: np.ndarray, ao_proj_k: np.ndarray, norm_thr: float = 1e-10):
    """calculate the tight-binding Hamiltonian and overlap matrix in the quasi-atomic orbital basis.
    
    Args:
        h_k (np.ndarray): Hamiltonian matrix in k-space
        s_k (np.ndarray): overlap matrix in k-space
        denmat_k_tilde (np.ndarray): density matrix in k-space
        ao_proj_k (np.ndarray): atomic orbital projection matrix in k-space
        norm_thr (float, optional): threshold for the norm of quasi-atomic orbitals. Defaults to 1e-10.
    
    Returns:
        tuple: tight-binding Hamiltonian and overlap matrix in the quasi-atomic orbital basis.
    """
    qo_k = denmat_k_tilde@ao_proj_k
    # normalize qo_k
    nqo = qo_k.shape[1]
    for iqo in range(nqo):
        norm_qo = np.sqrt(qo_k[:, iqo].conj().T@s_k@qo_k[:, iqo])
        if norm_qo < norm_thr:
            raise ValueError(f"Norm of qo_k is too small: norm < {norm_thr}.")
        qo_k[:, iqo] = qo_k[:, iqo] / norm_qo
    
    hqo_k = qo_k.conj().T@h_k@qo_k
    sqo_k = qo_k.conj().T@s_k@qo_k

    return hqo_k, sqo_k
# TESTED
def kR_transform(dst, mats_in: list, srcs: list, direction: str = "R->k"):
    """single-side transformation from k-space to R-space, or from R-space to k-space.
    
    Args:
        dst (np.ndarray): destination vector
        mats_in (list): input matrices
        srcs (list): source vectors
        direction (str, optional): "R->k" or "k->R". Defaults to "R->k".

    Returns:
        np.ndarray: transformed matrix
    """
    mat_out = np.zeros_like(mats_in[0], dtype=np.complex128)
    if len(mats_in) != len(srcs):
        print(f"Number of input matrices and source vectors are not equal: {len(mats_in)} (num of mats) != {len(srcs)} (num of srcs).")
        print(f"Shape of first matrix: {mats_in[0].shape}.")
        assert False
    phase = 0.0+1.0j if direction == "R->k" else 0.0-1.0j
    for mat_in, src in zip(mats_in, srcs):
        arg = np.exp(phase * 2.0 * np.pi * np.dot(dst, src))
        mat_out += arg * mat_in
    return mat_out / np.sqrt(len(mats_in))
# TESTED
def k_extrapolation(mats_k: list, kvecs_in: list, Rs: list, kvecs_out: list):
    """for kpoint extrapolation, from k-space to R-space, and then back to k-space.
    
    Args:
        mats_k (list): input matrices in k-space
        kvecs_in (list): input k-vectors
        Rs (list): supercell vectors
        kvecs_out (list): output k-vectors
    
    Returns:
        list: output matrices in k-space
    """
    mats_R = [kR_transform(Rs[i], mats_k, kvecs_in, direction="k->R") for i in range(len(Rs))]
    mats_kpath = [kR_transform(kvecs_out[i], mats_R, Rs, direction="R->k") for i in range(len(kvecs_out))]
    return mats_kpath

def norm_filter_k(denmat_k: np.ndarray, s_k: np.ndarray, ao_proj_k: np.ndarray, norm_thr: float = 1e-10):
    """filter out unbounded states represented by norms of atomic orbitals.
    
    Args:
        denmat_k (np.ndarray): density matrix in k-space
        s_k (np.ndarray): overlap matrix in k-space
        ao_proj_k (np.ndarray): atomic orbital projection matrix in k-space
        norm_thr (float, optional): threshold for the norm of atomic orbitals. Defaults to 1e-10.
    
    Returns:
        list: indices of unbounded states
    """
    ao_k = denmat_k@ao_proj_k
    norms = [np.sqrt(ao_k[:, iao].conj().T@s_k@ao_k[:, iao]) for iao in range(ao_k.shape[1])]
    unbounded_states = [iqo for iqo, norm in enumerate(norms) if norm < norm_thr]
    return unbounded_states

def nao_space_reduce(selected_indices: list, s_k: np.ndarray):
    return s_k[:, selected_indices]

"""Expect a file contain information describing row and column indices, like
<MATRIX_INDEX>
<ROW_INDEX>
0,0,0,0,0 0,0,0,1,0 0,0,0,2,0 0,0,1,0,0
</ROW_INDEX>
<COL_INDEX>
[it],[ia],[l],[izeta],[m] [it],[ia],[l],[izeta],[m] [it],[ia],[l],[izeta],[m] [it],[ia],[l],[izeta],[m] ...
</COL_INDEX>
</MATRIX_INDEX>
"""
def read_matrix_index(findex: str):
    with open(findex, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    # return bidirectional maps, tuple to index and index to tuple
    # row index
    row_index = lines[2].split()
    row_index = [tuple([int(x) for x in index.split(",")]) for index in row_index]
    row_index = {index: i for i, index in enumerate(row_index)}
    row_index_inv = {i: index for i, index in enumerate(row_index)}
    # col index
    col_index = lines[4].split()
    col_index = [tuple([int(x) for x in index.split(",")]) for index in col_index]
    col_index = {index: i for i, index in enumerate(col_index)}
    col_index_inv = {i: index for i, index in enumerate(col_index)}
    return row_index, row_index_inv, col_index, col_index_inv

import itertools as it
def mats_hs_k(path: str, nk: int, ib_min: int, ib_max: int = -1, qo: bool = True):
    print("Number of kpoints: ", nk)
    hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_d = read(path, nk)
    lowfs_k = lowfs_k if ib_max == -1 else [lowfs_k[ik][:, ib_min:ib_max] for ik in range(nk)]

    if not qo:
        return hs_k, ss_k, kvecs_d
    else:
        #                                         0   1   2   3     4    5    6     7    8    9     10   11   12    13    14   15   16   17    18    19   20   21
        # present is TZDP, 3s3p2d, arrange like [1s, 2s, 3s, 1p-1, 1p0, 1p1, 2p-1, 2p0, 2p1, 3p-1, 3p0, 3p1, 1d-2, 1d-1, 1d0, 1d1, 1d2, 2d-2, 2d-1, 2d0, 2d1, 2d2] 
        #                                 0   1   2     3    4    5     6    7    8     9    10   11    12
        # for DZP, 2s2p1d, arrange like [1s, 2s, 1p-1, 1p0, 1p1, 2p-1, 2p0, 2p1, 1d-2, 1d-1, 1d0, 1d1, 1d2]
        # only want 1s1p
        # therefore the szv_indices from TZDP:
        #szv_indices = [0, 3, 4, 5, 22, 25, 26, 27]
        # therefore the szv_indices from DZP:
        #szv_indices = [0, 2, 3, 4, 13, 15, 16, 17]
        #dzp_indices = [0, 1, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 22, 23, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38]
        #aos_proj_k = [ss_k[ik][:, szv_indices] for ik in range(nk)]
        #print("Shape of ao_proj_k: ", aos_proj_k[0].shape) # (8, 44), 8 is the number of selected atomic orbitals, 44 is the number of numerical atomic orbitals
        #aos_proj_k = ss_k
        
        nband = lowfs_k[0].shape[1]
        print("Number of bands selected: ", nband)
        m = aos_proj_k[0].shape[1] - nband
        print("Calculating quasi-atomic orbitals, number of bands unoccupied: ", m)
        denmats_k = [cal_denmat(lowf_k) for lowf_k in lowfs_k]
        ubsts_k = [norm_filter_k(denmat_k=denmats_k[ik], s_k=ss_k[ik], ao_proj_k=aos_proj_k[ik], norm_thr=1e-8) for ik in range(nk)]
        iks_ubsts = [ik for ik, ubsts in enumerate(ubsts_k) if len(ubsts) > 0]
        print("Kpoints with unbounded states: ", iks_ubsts)
        lowfs_bar_k = [lowdin_onW(s_k=ss_k[ik], 
                                  denmat_k=denmats_k[ik], 
                                  ao_proj_k=aos_proj_k[ik], 
                                  m=m) for ik in range(nk)]
        lowfs_tilde_k = [np.concatenate((lowfs_k[ik], lowfs_bar_k[ik]), axis=1) for ik in range(nk)]
        _one = [lowf_tilde_k.conj().T@ss_k[ik]@lowf_tilde_k for ik, lowf_tilde_k in enumerate(lowfs_tilde_k)]
        # immediate check of orthonormality
        for ik in range(nk):
            for i, j in it.product(range(nband+m), repeat=2):
                if i != j:
                    #print(_one[ik][i, j])
                    assert abs(_one[ik][i, j].real) < 1e-6 and abs(_one[ik][i, j].imag) < 1e-6
                else:
                    #print(_one[ik][i, j])
                    assert abs(_one[ik][i, j].real - 1) < 1e-5 and abs(_one[ik][i, j].imag) < 1e-6
        # go on
        denmats_tilde_k = [cal_denmat(lowf_tilde_k) for lowf_tilde_k in lowfs_tilde_k]
        hqos_k, sqos_k = [], []
        for ik in range(nk):
            hqo_k, sqo_k = cal_tb_hs(hs_k[ik], ss_k[ik], denmats_tilde_k[ik], aos_proj_k[ik])
            hqos_k.append(hqo_k)
            sqos_k.append(sqo_k)
        return hqos_k, sqos_k, kvecs_d

def monkhorst_pack(size):
    """Construct a uniform sampling of k-space of given size."""
    if np.less_equal(size, 0).any():
        raise ValueError(f'Illegal size: {list(size)}')
    kpts = np.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    return (kpts - np.array(size)//2)

def bad_condition_diag(Hk, Sk, eig_thr = 1e-10):
    """diagonalize the generalized eigenvalue problem Hk * x = e * Sk * x, and return the eigenvalues.

    Args:
        Hk (np.ndarray): Hamiltonian matrix
        Sk (np.ndarray): overlap matrix
        eig_thr (float, optional): for filtering out linear dependent components. Defaults to 1e-10.

    Returns:
        eigvals: eigenvalues of the generalized eigenvalue problem, sorted in ascending order, padded with a large number if necessary.
    """
    ndim = Hk.shape[0]
    s_eigvals, s_eigvecs = np.linalg.eigh(Sk)
    #print("eigenvalues of sk:\n", s_eigvals)
    s_eigvecs = s_eigvecs[:, s_eigvals > eig_thr]
    s_eigvals = s_eigvals[s_eigvals > eig_thr]
    s_eigvals = np.diag(1.0/np.sqrt(s_eigvals))
    s_eigvecs = np.dot(s_eigvecs, s_eigvals)
    Hk = np.dot(np.conj(s_eigvecs.T), np.dot(Hk, s_eigvecs))
    eigs = np.linalg.eigvalsh(Hk)
    # fill in the missing eigenvalues with a large number
    eigs = np.concatenate((eigs, np.ones(ndim-len(eigs))*1e6))
    return eigs

def production(path: str,                                   # path of OUT.${suffix} folder
               nkpts: int,                                  # number of kpoints in KPT file
               ib_min: int = 0,                             # minimum index of bands selected to reproduce
               ib_max: int = -1,                            # maximum index of bands selected to reproduce
               diag_thr: float = 1e-10,                     # threshold for filtering out linear dependent components
               error_estimation_with_lcao: bool = False):   # whether to estimate errors with LCAO
    """Production code for the whole process of calculating Hamiltonian and overlap matrices in QO basis.
    
    Args:
        path (str): path to the directory containing the output files
        nkpts (int): number of kpoints
        ib_min (int, optional): minimum index of bands selected to reproduce. Defaults to 0.
        ib_max (int, optional): maximum index of bands selected to reproduce. Defaults to -1.
        diag_thr (float, optional): threshold for filtering out linear dependent components. Defaults to 1e-10.
        error_estimation_with_lcao (bool, optional): whether to estimate errors with LCAO. Defaults to False.
    
    Returns:
        list: Hamiltonian matrices in QO basis
        list: overlap matrices in QO basis
        list: k-vectors
    """
    # first calculate things repsented by QO basis
    hqos_k, sqos_k, kvecs_d = mats_hs_k(path, nkpts, ib_min, ib_max, qo=True)
    eigs_qo = [bad_condition_diag(hqos_k[ik], sqos_k[ik], diag_thr) for ik in range(nkpts)]
    eigs_qo = np.asarray(eigs_qo)
    # then calculate things represented by AO basis, if necessary
    if error_estimation_with_lcao:
        haos_k, saos_k, kvecs_d = mats_hs_k(path, nkpts, ib_min, ib_max, qo=False)
        eigs_ao = [bad_condition_diag(haos_k[ik], saos_k[ik], diag_thr) for ik in range(nkpts)]
        eigs_ao = np.asarray(eigs_ao)
        # calculate errors
        e_errors_k = eigs_qo[:, ib_min:ib_max] - eigs_ao[:, ib_min:ib_max]
        e_errors_k = [np.linalg.norm(e_errors_k[ik]) for ik in range(nkpts)]
        exclude_kvecs = [ik for ik, e_error_k in enumerate(e_errors_k) if e_error_k > 1e-8]
        print("Number of kpoints excluded due to large error:", len(exclude_kvecs), "\nTheir indices:", exclude_kvecs)
        for ik in exclude_kvecs:
            print("kpoint: ", kvecs_d[ik], "has an error of: ", e_errors_k[ik])
            print("eigvals of sqos_k: ", sp.linalg.eigvalsh(sqos_k[ik]))
    else:
        exclude_kvecs = []
    # exclude kpoints with large errors
    hqos_k = [hqos_k[ik] for ik in range(nkpts) if ik not in exclude_kvecs]
    sqos_k = [sqos_k[ik] for ik in range(nkpts) if ik not in exclude_kvecs]
    kvecs_d = [kvecs_d[ik] for ik in range(nkpts) if ik not in exclude_kvecs]

    return hqos_k, sqos_k, kvecs_d
    # kpoints extrapolation: QO validation by reproducing band structure
    # Rlist = np.loadtxt(path+'/QO_supercells.dat', dtype=np.int32)
    # kpath = np.loadtxt('./band/kpoints')[:,1:4]
    # print("There are", len(Rlist), "supercells")

    # hk, sk = k_extrapolation(hqos_k, kvecs_d, Rlist, kpath), k_extrapolation(sqos_k, kvecs_d, Rlist, kpath)

    # eigs = []
    # for ik in range(len(kpath)):
    #     eigs.append(bad_condition_diag(hk[ik], sk[ik], diag_thr))
    # eigs = np.asarray(eigs)

    # band = np.loadtxt('./band/OUT.ABACUS/BANDS_1.dat')
    # figure_title = 'Band Structure of QO Hamiltonian: KPT' + str(ndim) + 'x' + str(ndim) + 'x' + str(ndim)
    # import matplotlib.pyplot as plt
    # plt.plot(band[:, 1], eigs[:,:]*13.6, 'o', markerfacecolor='none', markeredgecolor='b')
    # plt.plot(band[:, 1], band[:, 2:])
    # plt.ylim([-7.5,25])
    # plt.title(figure_title)
    # plt.savefig(fpic)

import unittest
class QOUnittest(unittest.TestCase):
    def test_read_mat_hs(self):
        for ik in range(ndim**3):
            print("Unittest test_read_mat_hs on kpoint", ik)
            fhs = path + f'/data-{ik}-S'
            mat = read_mat_hs(fhs)
            nrows, ncols = mat.shape
            # check square matrix
            self.assertEqual(nrows, ncols)
            # check hermitian
            self.assertListEqual(mat.tolist(), mat.conj().T.tolist())
            # check not-all-zero
            self.assertNotEqual(np.sum(np.abs(mat)), 0)

    def test_read_lowf_k(self):
        for ik in range(ndim**3):
            print("Unittest test_read_lowf_k on kpoint", ik)
            flowf = path + f"/LOWF_K_{ik+1}.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            nlocal, nband = lowf_k.shape
            # check number of bands
            self.assertGreater(nband, 0)
            # check number of local orbitals
            self.assertGreater(nlocal, 0)
            # check dimension of space
            self.assertGreater(nlocal, nband)

            fs = path + f'/data-{ik}-S'
            s_k = read_mat_hs(fs)
            nrows, ncols = s_k.shape
            self.assertEqual(nrows, ncols) # should be a square matrix
            self.assertEqual(nrows, nlocal) # should be overlap matrix of basis functions
            _one = lowf_k.conj().T@s_k@lowf_k
            # orthonormality
            for i in range(nband):
                for j in range(nband):
                    if i != j:
                        self.assertAlmostEqual(abs(_one[i, j].real), 0, delta=1e-5)
                        self.assertAlmostEqual(abs(_one[i, j].imag), 0, delta=1e-5)
                    else:
                        self.assertAlmostEqual(abs(_one[i, j].real), 1, delta=1e-5)
                        self.assertAlmostEqual(abs(_one[i, j].imag), 0, delta=1e-5)

    def test_read_ao_proj(self):
        for ik in range(ndim**3):
            print("Unittest test_read_ao_proj on kpoint", ik)
            fao_proj = path + f'/QO_ovlp_{ik}.dat'
            ao_proj, kvec_d = read_ao_proj(fao_proj)
            nrows, ncols = ao_proj.shape
            self.assertGreater(nrows, 0)
            self.assertGreater(ncols, 0)
            self.assertGreater(nrows, ncols) # space dimension comparision

    def test_cal_denmat(self):
        for ik in range(ndim**3):
            print("Unittest test_cal_denmat on kpoint", ik)
            flowf = path + f"/LOWF_K_{ik+1}.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            denmat_k = cal_denmat(lowf_k)
            
            nrows, ncols = denmat_k.shape
            self.assertEqual(nrows, ncols)
            nlocal, nband = lowf_k.shape
            self.assertEqual(nrows, nlocal)
            self.assertEqual(ncols, nlocal)
            
            denmat_dagger_k = denmat_k.conj().T
            # D = D^+
            for i in range(nrows):
                for j in range(ncols):
                    self.assertAlmostEqual(denmat_k[i, j], denmat_dagger_k[i, j], delta=1e-5)
                    
            fs = path + f'/data-{ik}-S'
            s_k = read_mat_hs(fs)
            D_k = denmat_k.tolist()
            DSD_k = (denmat_k@s_k@denmat_k).tolist()
            self.assertEqual(len(D_k), len(DSD_k))
            self.assertEqual(len(D_k[0]), len(DSD_k[0]))
            # D = DSD
            for i in range(len(D_k)):
                for j in range(len(D_k[0])):
                    self.assertAlmostEqual(D_k[i][j], DSD_k[i][j], delta=1e-3)

    def test_lowdin_onW(self):
        for ik in range(ndim**3):
            print("Unittest test_lowdin_onW on kpoint", ik)
            ndim_test = 2
            flowf = path + f"/LOWF_K_{ik+1}.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            denmat_k = cal_denmat(lowf_k)
            fs = path + f'/data-{ik}-S'
            s_k = read_mat_hs(fs)
            fao_proj = path + f'/QO_ovlp_{ik}.dat'
            ao_proj_k, kvec_d = read_ao_proj(fao_proj)
            
            lowf_bar = lowdin_onW(s_k=s_k,
                                    denmat_k=denmat_k,
                                    ao_proj_k=ao_proj_k,
                                    m = ndim_test)
            _one = lowf_bar.conj().T@s_k@lowf_bar
            self.assertEqual(_one.shape, (ndim_test, ndim_test))
            # orthonormality
            for i in range(ndim_test):
                for j in range(ndim_test):
                    if i != j:
                        self.assertAlmostEqual(_one[i, j], 0, delta=1e-6)
                    else:
                        self.assertAlmostEqual(_one[i, j], 1, delta=1e-6)
            # D = DSD
            lowf_tilde_k = np.concatenate((lowf_k, lowf_bar), axis=1)
            denmat_tilde_k = cal_denmat(lowf_tilde_k)
            D_tilde_k = denmat_tilde_k.tolist()
            DSD_tilde_k = (denmat_tilde_k@s_k@denmat_tilde_k).tolist()
            for i in range(len(D_tilde_k)):
                for j in range(len(D_tilde_k[0])):
                    self.assertAlmostEqual(D_tilde_k[i][j], DSD_tilde_k[i][j], delta=1e-5)

    def test_cal_tb_hs(self):
        print("Unittest test_cal_tb_hs")
        # eye
        nbasis = 16
        I = np.eye(nbasis, dtype=np.complex128)
        # CASE 1: the most ideal case
        # create an Hermite matrix as model Hamiltonian
        H = np.random.rand(nbasis, nbasis) + 1.0j*np.random.rand(nbasis, nbasis)
        H = H + H.conj().T
        # eye as model overlap matrix
        S = I
        # eye as model denmat
        D = I
        # eye as model ao_proj
        A = I
        # test
        Htb, Stb = cal_tb_hs(H, S, D, A)
        # check shape
        self.assertEqual(Htb.shape, (nbasis, nbasis))
        self.assertEqual(Stb.shape, (nbasis, nbasis))
        # check hermitian
        self.assertListEqual(Htb.tolist(), Htb.conj().T.tolist())
        self.assertListEqual(Stb.tolist(), Stb.conj().T.tolist())
        # check not-all-zero
        self.assertNotEqual(np.sum(np.abs(Htb)), 0)
        self.assertNotEqual(np.sum(np.abs(Stb)), 0)
        # CASE 2: less ideal case               
        # directly truncate A, for the meaning the AO basis is not as many as original basis
        nAO = 3
        A = A[:, :nAO]
        D = I[:, :nAO]@I[:nAO, :] # usually D is not a full-rank matrix
        Htb, Stb = cal_tb_hs(H, S, D, A)
        # check shape
        self.assertEqual(Htb.shape, (nAO, nAO))
        self.assertEqual(Stb.shape, (nAO, nAO))
        # check hermitian
        self.assertListEqual(Htb.tolist(), Htb.conj().T.tolist())
        self.assertListEqual(Stb.tolist(), Stb.conj().T.tolist())
        # check not-all-zero
        self.assertNotEqual(np.sum(np.abs(Htb)), 0)
        self.assertNotEqual(np.sum(np.abs(Stb)), 0)

    def test_kR_transform(self):
        for nk in [1, 3, 5]:
            kpoints = it.product(np.arange(-0.5+0.5/nk, 0.5, 1.0/nk), repeat=3)
            kpoints = np.array(list(kpoints))
            Rs = it.product(range(-int((nk-1)/2), int((nk+1)/2)), repeat=3)
            Rs = np.array(list(Rs))
            assert len(kpoints) == len(Rs)
            ndim = 20 # square matrix dimension
            mats_0 = [np.random.rand(ndim, ndim) for i in range(len(kpoints))]
            mats_k = [kR_transform(kpoints[i], mats_0, Rs, "R->k") for i in range(len(kpoints))]
            mats_kR = [kR_transform(Rs[i], mats_k, kpoints, "k->R") for i in range(len(kpoints))]
            for ik in range(len(kpoints)):
                for i in range(ndim):
                    for j in range(ndim):
                        self.assertAlmostEqual(mats_kR[ik][i, j].real, mats_0[ik][i, j].real, delta=1e-12)
                        self.assertAlmostEqual(mats_kR[ik][i, j].imag, mats_0[ik][i, j].imag, delta=1e-12)

    def test_k_extrapolation(self):
        for nk in [1, 3, 5]:
            kpoints = it.product(np.arange(-0.5+0.5/nk, 0.5, 1.0/nk), repeat=3)
            kpoints = np.array(list(kpoints))
            Rs = it.product(range(-int((nk-1)/2), int((nk+1)/2)), repeat=3)
            Rs = np.array(list(Rs))
            assert len(kpoints) == len(Rs)
            ndim = 20 # square matrix dimension
            mats_0 = [np.random.rand(ndim, ndim) for i in range(len(kpoints))]
            mats_k2R2k = k_extrapolation(mats_0, kpoints, Rs, kpoints)
            for ik in range(len(kpoints)):
                for i in range(ndim):
                    for j in range(ndim):
                        self.assertAlmostEqual(mats_k2R2k[ik][i, j].real, mats_0[ik][i, j].real, delta=1e-12)
                        self.assertAlmostEqual(mats_k2R2k[ik][i, j].imag, mats_0[ik][i, j].imag, delta=1e-12)

if __name__ == "__main__":

    # False to run unittest, True to run production
    run_type = "production" # can be "unittest" or "production"
    ndim = 7 # dimension of Monkhorst-Pack grid
    qo_basis = "hydrogen"
    qo_strategy = "minimal-valence"
    path = f"./scf/{qo_basis}/{qo_strategy}/OUT.ABACUS-{ndim}{ndim}{ndim}/"
    #fpic = f"{qo_basis}-{qo_strategy}-{ndim}{ndim}{ndim}.png"
    # only for production
    eig_thr = 1e-10
    ib_min, ib_max = 0, 4

    if run_type == "production":
        nkpts = ndim*ndim*ndim
        hqos_k, sqos_k, kvecs_d = production(path, nkpts, ib_min, ib_max, eig_thr, error_estimation_with_lcao=True)

    elif run_type == "unittest":
        unittest.main()
