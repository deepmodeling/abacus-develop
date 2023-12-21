import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import source.components.safe_guard as sg

class toQO_Calculator:
    """python-end the Quasiatomic orbital (QO) analysis
    """
    sg_: sg.toQO_SafeGuard

    def __init__(self) -> None:
        self.sg_ = sg.toQO_SafeGuard()

    def folding_HR(self, matrices_R: list, supercells: list, kpoint: np.ndarray) -> np.ndarray:
        """calculate fold matrix from R to k-space

        Args:
            matrices_R (list): matrices in R-space
            supercells (list): corresponding supercells direct vectors n1n2n3
            kpoint (np.ndarray): one specific kpoint direct coordinates

        Returns:
            np.ndarray: matrix in k-space
        """
        matrix_k = np.zeros_like(matrices_R[0])
        print("Calculate fold matrix from R to k-space.")
        for iR in range(len(supercells)):
            arg = np.exp(1j * supercells[iR] @ kpoint * 2 * np.pi)
            matrix_k += arg * matrices_R[iR]
        return matrix_k
    
    def unfolding_Hk(self, matrices_k: list, kpoints: list, supercell: np.ndarray) -> np.ndarray:
        """calculate unfold matrix from k to R-space

        Args:
            matrices_k (list): matrices in k-space
            kpoints (list): corresponding equivalent kpoints, arranged in the same order as matrices_k
            supercell (np.ndarray): one specific supercell direct coordinates

        Returns:
            np.ndarray: matrix in R-space
        """

        matrix_R = np.zeros_like(matrices_k[0])
        print("Calculate unfold matrix from k to R-space.")
        for ik in range(len(kpoints)):
            kpoint = kpoints[ik]
            arg = np.exp(-1j * kpoint @ supercell * 2 * np.pi)
            matrix_R += arg * matrices_k[ik]
        return matrix_R
    
    def projto_nao(self, sk: np.ndarray, saok: np.ndarray) -> np.ndarray:
        """calculate anything in nao representation with the nao projector

        Args:
            sk (np.ndarray): nao overlap matrix in k-space
            saok (np.ndarray): ao-nao overlap matrix in k-space

        Returns:
            np.ndarray: anything in nao representation
        """
        print("Calculate Atomic Orbital (AO) in NAO representation.")
        psi_chi = np.linalg.solve(sk.T, saok.T).T
        return psi_chi

    def projto_eigstate(self, psi_lcao: np.ndarray, saok: np.ndarray) -> np.ndarray:
        """calculate any states projected onto eigenstates of Hamiltonian in k-space represented by NAO

        Args:
            psi_lcao (np.ndarray): states in LCAO representation
            saok (np.ndarray): ao-nao overlap matrix in k-space

        Returns:
            np.ndarray: states projected onto eigenstates of Hamiltonian in k-space represented by NAO
        """
        psi_chi_para = saok @ psi_lcao.conj().T @ psi_lcao
        return psi_chi_para
    
    def canonical_orthogonalization(self, psi_chi_orth: np.ndarray, sk: np.ndarray, m: int) -> np.ndarray:
        """this is a general canonical orthogonalization procedure for any subspace

        Args:
            psi_chi_orth (np.ndarray): states in one certain representation, here is nao
            sk (np.ndarray): basis function overlap matrix
            m (int): number of states retained

        Returns:
            np.ndarray: states orthogonalized
        """
        print("Perform canonical orthogonalization for newly appended subspace from AO introduction.")
        wk = psi_chi_orth.conj() @ sk @ psi_chi_orth.T
        yk_k, vk_k = la.eigh(wk)
        print("Eigenvalues of W(k) are: \n", yk_k)
        # truncate m largest eigenvalues and eigenvectors
        yk_k = yk_k[-m:]
        vk_k = vk_k[:, -m:]
        print("Get ", m, " largest eigenvalues and eigenvectors. Selected eigenvalues are: ", yk_k)
        # calculate psi_complem
        yk_k = np.diag([np.sqrt(1/yk_k[i]) for i in range(m)])
        psi_complem = yk_k @ vk_k.T @ psi_chi_orth
        return psi_complem
    
    def merge_space(self, psi_lcao: np.ndarray, psi_complem: np.ndarray, hk: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate psi_exten

        psi_exten = [psi_lcao, psi_complem]

        Args:
            psi_lcao (np.ndarray): states in LCAO representation
            psi_complem (np.ndarray): states in complementary representation
            hk (np.ndarray): Hamiltonian in k-space
            sk (np.ndarray): overlap in k-space

        Returns:
            np.ndarray: extended states in k-space
        """
        print("Combine psi_lcao and psi_complem to get extended wavefunction (empty states are appended).")
        psi_exten = np.concatenate((psi_lcao, psi_complem), axis=0)
        # check orthogonality
        print("Check orthogonality of psi_exten.")
        sk_exten = psi_exten.conj() @ sk @ psi_exten.T
        """# if sk_exten is not identity matrix?
        error = self.sg_.check_identity(sk_exten)
        while error > 1e-6:
            print("Error of sk_exten is: ", error)
            print("Reorthogonalize psi_exten.")
            eigvals_exten, eigvecs_exten = la.eigh(sk_exten)
            print("Eigenvalues of sk_exten are: \n", eigvals_exten)
            psi_exten = np.diag([np.sqrt(1/eigvals_exten[i]) for i in range(eigvals_exten.shape[0])]) @ eigvecs_exten.T @ psi_exten
            sk_exten = psi_exten.conj() @ sk @ psi_exten.T
            error = self.sg_.check_identity(sk_exten)"""
            
        # if hk_exten is not diagonal in supspace psi_lcao?
        hk_exten = psi_exten.conj() @ hk @ psi_exten.T
        self.sg_.check_diagonal(hk_exten[:psi_lcao.shape[0], :psi_lcao.shape[0]])
        # get extended energy spectrum
        print("Get extended energy spectrum.")
        eigvals_exten, eigvecs_exten = la.eigh(hk_exten, sk_exten)
        print("Eigenvalues of H_exten(k) are: \n", eigvals_exten)
        return psi_exten

    def calculate_qo(self, saok: np.ndarray, psi_exten: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate qo represented by NAO in kspace, return with nrows = nchi and ncols = nphi

        Args:
            saok (np.ndarray): ao-nao overlap matrix in k-space
            psi_exten (np.ndarray): extended states in k-space
            sk (np.ndarray): overlap in k-space

        Returns:    
            np.ndarray: qo represented by NAO in kspace
        """
        print("Calculate QO.")
        qo = saok.conj() @ psi_exten.conj().T @ psi_exten
             # this saok is overlap between AO and NAO, line is AO, column is NAO
        # then normalize qo
        #for i in range(qo.shape[0]):
        #    qo[i, :] = qo[i, :] / np.sqrt(qo[i, :] @ sk @ qo[i, :].conj().T)
            #print("QO Normalization: after, norm of QO ", i, " is: ", qo[i, :] @ sk @ qo[i, :].conj().T)
        return qo

    def calculate_hqok(self, qo: np.ndarray, hk: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate hqok

        Args:
            qo (np.ndarray): QO in k-space
            hk (np.ndarray): Hamiltonian represented by QO in k-space
            sk (np.ndarray): QO overlap in k-space

        Returns:
            np.ndarray: hqok
        """
        print("Calculate hamiltonian matrix in QO basis in k-space.")
        hqok = qo.conj() @ hk @ qo.T
        sqok = qo.conj() @ sk @ qo.T # this is overlap matrix in QO basis

        eigval_s, eigvec_s = la.eigh(sqok)
        print("Eigenvalues of overlap of in basis Sqo(k) are: \n", eigval_s)
        eigvals_qo, eigvecs_qo = la.eigh(hqok, sqok)
        print("Eigenvalues of Hamiltonian in QO basis Hqo(k) are: \n", eigvals_qo)
        
        eigvals_nao, eigvecs_nao = la.eigh(hk, sk)
        print("Eigenvalues of Hamiltonian in NAO basis H(k) are: \n", eigvals_nao)
        return hqok, sqok
