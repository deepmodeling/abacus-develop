import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import components.safe_guard as sg

class toQO_Calculator:
    """python-end the Quasiatomic orbital (QO) analysis
    """
    sg_: sg.toQO_SafeGuard

    def __init__(self) -> None:
        self.sg_ = sg.toQO_SafeGuard()

    def projto_nao(self, sk: np.ndarray, sqok: np.ndarray) -> np.ndarray:
        """calculate anything in nao representation with the nao projector

        Args:
            sk (np.ndarray): nao overlap matrix in k-space
            sqok (np.ndarray): ao-nao overlap matrix in k-space

        Returns:
            np.ndarray: anything in nao representation
        """
        print("Calculate Atomic Orbital (AO) in NAO representation.")
        psi_chi = np.linalg.solve(sk.T, sqok.T).T
        return psi_chi

    def projto_eigstate(self, psi_lcao: np.ndarray, sqok: np.ndarray) -> np.ndarray:
        """calculate any states projected onto eigenstates of Hamiltonian in k-space represented by NAO
        """
        psi_chi_para = sqok @ psi_lcao.conj().T @ psi_lcao
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

    def calculate_qo(self, sqok: np.ndarray, psi_exten: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate qo represented by NAO in kspace, return with nrows = nchi and ncols = nphi
        """
        print("Calculate QO.")
        qo = sqok.conj() @ psi_exten.conj().T @ psi_exten
             # this sqok is overlap between AO and NAO, line is AO, column is NAO
        # then normalize qo
        for i in range(qo.shape[0]):
            qo[i, :] = qo[i, :] / np.sqrt(qo[i, :] @ sk @ qo[i, :].conj().T)
            #print("QO Normalization: after, norm of QO ", i, " is: ", qo[i, :] @ sk @ qo[i, :].conj().T)
        return qo

    def calculate_hqok(self, qo: np.ndarray, hk: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate hqok
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
        return hqok

    def unfolding_hk(self, hqoks: list, kpoints: list, supercell: np.ndarray) -> np.ndarray:
        """calculate hqoR
        """
        hqoR = np.zeros_like(hqoks[0])
        print("Calculate hamiltonian matrix in QO basis in R-space. Present supercell is: ", supercell, ".")
        for ik in range(len(kpoints)):
            arg = np.exp(1j * kpoints[ik] @ supercell * 2 * np.pi)
            hqoR += arg * hqoks[ik]
        
        hqoR_fname = "hqoR_" + "_".join([str(i) for i in supercell]) + ".txt"
        np.savetxt(hqoR_fname, hqoR)
        return hqoR