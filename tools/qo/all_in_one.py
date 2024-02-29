import numpy as np
def cxx_topycomplex(num: str):
    """cxx prints complex numbers in the form (a,b)"""
    num = num.replace('(', '').replace(')', '')
    num = num.split(',')
    return complex(float(num[0]), float(num[1]))

def read_mat_hs(fhs):
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

def read_lowf(flowf):
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
    kvec_d = None
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
                kvec_d = np.array([float(x) for x in line.split()])
    return lowf_k, kvec_d

def read_ao_proj(fao_proj):
    with open(fao_proj, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    kvec_d = None
    if lines[0].startswith("KPOINT_COORDINATE"):
        kvec_d = np.array([float(x) for x in lines[0].split[-3:]])
        lines = lines[1:]
    nlocal = len(lines[0].split())
    nao = len(lines)
    ao_proj = np.zeros((nao, nlocal), dtype=np.complex128)
    for i, line in enumerate(lines):
        nums = [cxx_topycomplex(num) for num in line.split()]
        for j in range(nlocal):
            ao_proj[i, j] = nums[j]
    return ao_proj, kvec_d

import os
def read(path: str, nk: int) -> tuple[list, list, list, list, list]:

    hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_d = [], [], [], [], []
    for ik in range(nk):
        # H(k) and S(k)
        fh = os.path.join(path, f"data-{ik}-H")
        hs_k.append(read_mat_hs(fh))
        fs = os.path.join(path, f"data-{ik}-S")
        ss_k.append(read_mat_hs(fs))
        # \psi(k)
        flowf = os.path.join(path, f"LOWF_K_{ik+1}.txt")
        lowf_k, kvec_d1 = read_lowf(flowf)
        lowfs_k.append(lowf_k)
        # <\phi(k)|A(k)>
        fao_proj = os.path.join(path, f"QO_ovlp_{ik}.dat")
        ao_proj_k, kvec_d2 = read_ao_proj(fao_proj)
        aos_proj_k.append(ao_proj_k)
        # kpoints
        if kvec_d1 is not None and kvec_d2 is not None:
            for i in range(3):
                assert kvec_d1[i] == kvec_d2[i]
        kvec_d = kvec_d1 if kvec_d1 is not None else kvec_d2
        assert kvec_d is not None
        kvecs_d.append(kvec_d)

    return hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_d

def cal_denmat(lowf_k: np.ndarray) -> np.ndarray:
    return lowf_k@lowf_k.T

def proj_onsubspace(denmat_k: np.ndarray, ao_proj_k: np.ndarray):
    return denmat_k@ao_proj_k.T

# def cal_psiperp(s_k: np.ndarray, denmat_k: np.ndarray, ao_proj_k: np.ndarray) -> np.ndarray:
#     nrow = s_k.shape[0]
#     ncol = s_k.shape[1]
#     assert nrow == ncol
#     ovlpinv_k = np.linalg.solve(s_k, np.eye(N=nrow, M=ncol, dtype=np.complex128))
#     return (ovlpinv_k - denmat_k)@ao_proj_k.T

import scipy as sp
def lowdin_onW(s_k: np.ndarray, denmat_k: np.ndarray, ao_proj_k: np.ndarray, m: int):
    nrow = s_k.shape[0]
    ncol = s_k.shape[1]
    assert nrow == ncol
    sinv_k = np.linalg.solve(s_k, np.eye(N=nrow, M=ncol, dtype=np.complex128))
    W_mat = ao_proj_k.conj()@(sinv_k - denmat_k)@ao_proj_k.T
    eigvals, eigvecs = sp.linalg.eigh(W_mat)
    eigvals = eigvals[-m:]
    lambdas = np.diag([1./np.sqrt(eigval) for eigval in eigvals])
    eigvecs = eigvecs[:,-m:]
    
    return (sinv_k - denmat_k)@ao_proj_k.T@eigvecs@lambdas
    # call return as lowf_k_bar

def cal_tb_HS(h_k: np.ndarray, denmat_k_tilde: np.ndarray, ao_proj_k: np.ndarray):
    qo_k = denmat_k_tilde@ao_proj_k.T
    return qo_k.T@h_k@qo_k, ao_proj_k.conj()@qo_k

def kR_transform(dst, mats_in: list, srcs: list, direction: str = "R->k"):
    mat_out = np.zeros_like(mats_in[0], dtype=np.complex128)
    phase = 0.0+1.0j if direction == "R->k" else 0.0-1.0j
    for mat_in, src in zip(mats_in, srcs):
        arg = np.exp(phase * 2.0 * np.pi * np.dot(dst, src))
        mat_out += arg * mat_in / np.sqrt(len(srcs))
    return mat_out

def mats_hs_k(path: str, nk: int, qo: bool = True):
    hs_k, ss_k, lowfs_k, aos_proj_k, kvecs_d = read(path, nk)
    if not qo:
        return hs_k, ss_k, kvecs_d
    else:
        denmats_k = [cal_denmat(lowf_k) for lowf_k in lowfs_k]
        lowfs_bar_k = [lowdin_onW(ss_k[ik], denmats_k[ik], aos_proj_k[ik], 1) for ik in range(nk)]
        lowfs_tilde_k = [np.concatenate((lowfs_k[ik], lowfs_bar_k[ik]), axis=1) for ik in range(nk)]
        denmats_tilde_k = [cal_denmat(lowf_tilde_k) for lowf_tilde_k in lowfs_tilde_k]
        hqos_k, sqos_k = zip(*[cal_tb_HS(hs_k[ik], denmats_tilde_k[ik], aos_proj_k[ik]) for ik in range(nk)])
        return hqos_k, sqos_k, kvecs_d

if __name__ == "__main__":

    # hqos_k, sqos_k, kvecs_d = mats_hs_k(path='/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/', nk=8)
    # print(kvecs_d)
    # exit()

    import unittest
    import itertools as it
    import matplotlib.pyplot as plt
    class QOUnittest(unittest.TestCase):
        def test_read_mat_hs(self):
            fhs = '/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/data-0-S'
            mat = read_mat_hs(fhs)
            self.assertEqual(mat.shape, (18, 18))
            self.assertListEqual(mat.tolist(), mat.conj().T.tolist())

            self.assertEqual(mat[0, 0].real, 7.89834)
            self.assertEqual(mat[0, 0].imag, 0)
            self.assertEqual(mat[0, 1].real, -1.61301)
            self.assertEqual(mat[0, 1].imag, 0)
            self.assertEqual(mat[0, 17].real, 3.0748e-16)
            self.assertEqual(mat[0, 17].imag, 0)
            self.assertEqual(mat[17, 17].real, 0.0323155)
            self.assertEqual(mat[17, 17].imag, 0)
        
        def test_read_lowf_k(self):
            flowf = "/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/LOWF_K_1.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            self.assertEqual(lowf_k.shape, (18, 15))
            self.assertEqual(lowf_k[0, 0].real, 3.0183482687375629005543942e-01)
            self.assertEqual(lowf_k[0, 0].imag, 0)
            self.assertEqual(lowf_k[0, 1].real, 9.3382487799254196167542559e-17) # ilocal = 0, ib = 1
            self.assertEqual(lowf_k[0, 1].imag, 0)
            self.assertEqual(lowf_k[0, 14].real, 9.9397671895746906784146993e-16) # ilocal = 0, ib = 14
            self.assertEqual(lowf_k[0, 14].imag, 0)
            self.assertEqual(lowf_k[17, 14].real, -2.6825122405319675917289896e-14) # ilocal = 17, ib = 14
            self.assertEqual(lowf_k[17, 14].imag, 0)

            fs = '/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/data-0-S'
            s_k = read_mat_hs(fs)
            _one = lowf_k.conj().T@s_k@lowf_k
            self.assertEqual(_one.shape, (15, 15))
            for i in range(15):
                for j in range(15):
                    if i != j:
                        self.assertAlmostEqual(abs(_one[i, j].real), 0, delta=1e-5)
                        self.assertAlmostEqual(abs(_one[i, j].imag), 0, delta=1e-5)
                    else:
                        self.assertAlmostEqual(abs(_one[i, j].real), 1, delta=1e-5)
                        self.assertAlmostEqual(abs(_one[i, j].imag), 0, delta=1e-5)

        def test_read_ao_proj(self):
            fao_proj = "/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/QO_ovlp_1.dat"
            ao_proj, kvec_d = read_ao_proj(fao_proj)
            self.assertEqual(ao_proj.shape, (10, 18))
            self.assertEqual(ao_proj[0, 0].real, 2.03668722687917e+00)
            self.assertEqual(ao_proj[0, 0].imag, -4.62306304037626e-18)
            self.assertEqual(ao_proj[0, 17].real, -3.95927554936294e-16)
            self.assertEqual(ao_proj[0, 17].imag, 1.48601964485274e-18)
            self.assertEqual(ao_proj[9, 17].real, -3.57653934805274e-01)
            self.assertEqual(ao_proj[9, 17].imag, -7.87036962379485e-19)

        def test_cal_denmat(self):
            # D = DSD
            flowf = "/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/LOWF_K_1.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            denmat_k = cal_denmat(lowf_k)
            fs = '/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/data-0-S'
            s_k = read_mat_hs(fs)
            D_k = denmat_k.tolist()
            DSD_k = (denmat_k@s_k@denmat_k).tolist()
            self.assertEqual(len(D_k), len(DSD_k))
            self.assertEqual(len(D_k[0]), len(DSD_k[0]))
            for i in range(len(D_k)):
                for j in range(len(D_k[0])):
                    self.assertAlmostEqual(D_k[i][j], DSD_k[i][j], delta=1e-3)

        def test_lowdin_onW(self):
            ndim_test = 2

            flowf = "/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/LOWF_K_1.txt"
            lowf_k, kvec_d = read_lowf(flowf)
            denmat_k = cal_denmat(lowf_k)
            fs = '/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/data-0-S'
            s_k = read_mat_hs(fs)
            fao_proj = '/root/abacus-develop/representation/examples/scf/lcao_Cu/OUT.ABACUS/QO_ovlp_0.dat'
            ao_proj_k, kvec_d = read_ao_proj(fao_proj)
            
            lowf_bar = lowdin_onW(s_k=s_k,
                                  denmat_k=denmat_k,
                                  ao_proj_k=ao_proj_k,
                                  m = ndim_test)
            _one = lowf_bar.conj().T@s_k@lowf_bar
            self.assertEqual(_one.shape, (ndim_test, ndim_test))
            for i in range(ndim_test):
                for j in range(ndim_test):
                    if i != j:
                        self.assertAlmostEqual(_one[i, j], 0, delta=1e-6)
                    else:
                        self.assertAlmostEqual(_one[i, j], 1, delta=1e-6)

        def test_kR_transform(self):
            nk = 20 # number of kpoints, nk = nR
            kpoints = it.product(np.arange(0, nk/(nk+1), 1/(nk+1)), repeat=3)
            kpoints = np.array(list(kpoints))
            Rs = it.product(range(-int((nk-1)/2), int((nk+1)/2)), repeat=3)
            Rs = np.array(list(Rs))
            ndim = 20 # square matrix dimension
            mats_0 = [np.random.rand(ndim, ndim) for i in range(nk)]
            mats_k = [kR_transform(kpoints[i], mats_0, Rs, "R->k") for i in range(nk)]
            mats_kR = [kR_transform(Rs[i], mats_k, kpoints, "k->R") for i in range(nk)]
            for ik in range(nk):
                for i in range(ndim):
                    for j in range(ndim):
                        if i == j:
                            self.assertAlmostEqual(mats_kR[ik][i, j].real, mats_0[ik][i, j].real, delta=1e-2)
                            self.assertAlmostEqual(mats_kR[ik][i, j].imag, mats_0[ik][i, j].imag, delta=1e-2)
                        else:
                            self.assertAlmostEqual(mats_kR[ik][i, j].real, mats_0[ik][i, j].real, delta=1e-2)
                            self.assertAlmostEqual(mats_kR[ik][i, j].imag, mats_0[ik][i, j].imag, delta=1e-2)

    unittest.main()