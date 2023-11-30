import numpy as np
from scipy.linalg import eigh
import tools.wavefunction as wf
import tools.hamiltonian as ham
import tools.qo_ovlp as qov

class toQO_arithmetic:
    """python-end the Quasiatomic orbital (QO) analysis
    """
    # file path
    data_path: str # path of data
    # read from file
    nkpts: int # number of kpoints
    kpoints: np.ndarray # kpoints direct coordinates
    hk: np.ndarray # hamiltonian matrix in k-space
    sk: np.ndarray # overlap matrix in k-space
    nbands: int # number of bands
    sqok: np.ndarray # Atomic Orbital (AO)-Numerical Atomic Orbital (NAO) overlap matrix in k-space
    psi_lcao: np.ndarray # wavefunction represented by numerical atomic orbitals
    # calculated
    psi_chi: np.ndarray # AO basis represented by NAOs
    psi_chi_para: np.ndarray # band-parallel component of QO wavefunction represented by NAOs
    psi_chi_orth: np.ndarray # band-orthogonal component of QO wavefunction represented by NAOs
    wk: np.ndarray # AO overlap matrix in k-space
    psi_complem: np.ndarray # states from AO introduction
    psi_exten: np.ndarray # extended wavefunction
    omega: np.ndarray # density matrix transformed onto QO basis
    hqok: np.ndarray # hamiltonian matrix in QO basis in k-space
    hqoR: np.ndarray # hamiltonian matrix in QO basis in R-space
    # dimensions
    nchi: int # number of AOs
    nphi: int # number of NAOs
    m: int # number of QOs selected
    nbands_plus: int # number of bands plus m

    def __init__(self) -> None:
        pass

    def initialize(self, nkpts: int, m: int) -> None:
        """initialize the class

        Args:
            nkpts (int): number of kpoints
        """
        # data import
        self.nkpts = nkpts
        self.hk, self.sk, self.nphi = ham.parse(nkpts, path=self.data_path)
        self.sqok = qov.parse(nkpts, path=self.data_path)
        self.psi_lcao, self.kpoints = wf.parse(nkpts, path=self.data_path)
        self.nbands = self.psi_lcao.shape[1]
        self.m = m
        self.nbands_plus = self.nbands + self.m
        # data initialization
        self.nchi = self.sqok.shape[1]
        self.nphi = self.sqok.shape[2]
        # allocate memory
        self.psi_chi = np.zeros(
            (self.nkpts, self.nchi, self.nphi)
                , dtype=np.complex128)
        self.psi_chi_para = np.zeros(
            (self.nkpts, self.nchi, self.nphi)
                , dtype=np.complex128)
        self.psi_chi_orth = np.zeros(
            (self.nkpts, self.nchi, self.nphi)
                , dtype=np.complex128)
        self.wk = np.zeros(
            (self.nkpts, self.nchi, self.nchi)
                , dtype=np.complex128)
        self.psi_complem = np.zeros(
            (self.nkpts, self.m, self.nphi)
                , dtype=np.complex128)
        self.psi_exten = np.zeros(
            (self.nkpts, self.nbands_plus, self.nphi)
                , dtype=np.complex128)
        self.omega = np.zeros(
            (self.nkpts, self.nchi, self.nphi)
                , dtype=np.complex128)
        self.hqok = np.zeros(
            (self.nkpts, self.nchi, self.nchi)
                , dtype=np.complex128)
        self.hqoR = np.zeros(
            (self.nchi, self.nchi)
                , dtype=np.complex128)
        
    def set_path(self, path: str) -> None:
        """set the path of data

        Args:
            path (str): path of data
        """
        self.data_path = path

    def calculate_psi_chi(self) -> None:
        """calculate psi_chi
        """
        print("Calculate Atomic Orbital (AO) in NAO representation.")
        for ik in range(self.nkpts):
            self.psi_chi[ik] = np.linalg.solve(self.sk[ik].T, self.sqok[ik].T).T
    
    def calculate_psi_chi_para(self) -> None:
        """calculate psi_chi_para
        """
        print("Calculate band-parallel component of AO.")
        for ik in range(self.nkpts):
            self.psi_chi_para[ik] = self.sqok[ik] @ self.psi_lcao[ik].conj().T @ self.psi_lcao[ik]

    def calculate_psi_chi_orth(self) -> None:
        """calculate psi_chi_orth
        """
        print("Calculate band-orthogonal component of AO.")
        for ik in range(self.nkpts):
            self.psi_chi_orth[ik] = self.psi_chi[ik] - self.psi_chi_para[ik]
    
    def calculate_wk(self) -> None:
        """calculate wk
        """
        print("Perform canonical orthogonalization for newly appended subspace from AO introduction.")
        for ik in range(self.nkpts):
            self.wk[ik] = self.psi_chi_orth[ik] @ self.sk[ik] @ self.psi_chi_orth[ik].conj().T
    
    def calculate_psi_complem(self) -> None:
        """calculate psi_complem
        """
        for ik in range(self.nkpts):
            print("For kpoint ", ik, ", whose Cartesian coordinates are: (%8.4e, %8.4e, %8.4e)"%(self.kpoints[ik][0], self.kpoints[ik][1], self.kpoints[ik][2]))
            yk_k, vk_k = eigh(self.wk[ik])
            print("Eigenvalues of W(k) are: \n", yk_k)
            yk_k = np.diag(yk_k)
            # truncate m largest eigenvalues and eigenvectors
            yk_k = yk_k[-self.m:, -self.m:]
            vk_k = vk_k[:, -self.m:]
            print("Get ", self.m, " largest eigenvalues and eigenvectors. Selected eigenvalues are: ", yk_k.diagonal())
            # calculate psi_complem
            self.psi_complem[ik] = (np.linalg.inv(yk_k))**(0.5) @ vk_k.conj().T @ self.psi_chi[ik]
            # normalization
            for i in range(self.m): # loop over newly introduced states (bands)
                normalize_factor = np.sqrt(self.psi_complem[ik, i] @ self.sk[ik] @ self.psi_complem[ik, i].conj().T).real
                self.psi_complem[ik, i] /= normalize_factor
                print(
                    "At kpoint (%8.4e, %8.4e, %8.4e)"%(self.kpoints[ik][0], self.kpoints[ik][1], self.kpoints[ik][2]), 
                    "normalization factor for newly introduced state ", i, " is: \n", normalize_factor)

    def calculate_psi_exten(self) -> None:
        """calculate psi_exten

        psi_exten = [psi_lcao, psi_complem]
        """
        print("Combine psi_lcao and psi_complem to get extended wavefunction (empty states are appended).")
        for ik in range(self.nkpts):
            self.psi_exten[ik] = np.concatenate((self.psi_lcao[ik], self.psi_complem[ik]), axis=0)
    
    def calculate_omega(self) -> None:
        """calculate omega
        """
        print("Calculate density matrix transformed onto QO basis.")
        for ik in range(self.nkpts):
            self.omega[ik] = self.sqok[ik] @ self.psi_exten[ik].conj().T @ self.psi_exten[ik]
    
    def calculate_hqok(self) -> None:
        """calculate hqok
        """
        print("Calculate hamiltonian matrix in QO basis in k-space.")
        for ik in range(self.nkpts):
            self.hqok[ik] = self.omega[ik] @ self.hk[ik] @ self.omega[ik].conj().T
    
    def calculate_hqoR(self, supercell: np.ndarray) -> None:
        """calculate hqoR
        """
        print("Calculate hamiltonian matrix in QO basis in R-space. Present supercell is: ", supercell, ".")
        for ik in range(self.nkpts):
            arg = np.exp(1j * self.kpoints[ik] @ supercell * 2 * np.pi)
            self.hqoR += arg * self.hqok[ik]

class toQO:

    arithmetic: toQO_arithmetic

    def __init__(self) -> None:
        pass

    def initialize(self, nkpts: int, m: int, path: str) -> None:
        """initialize the class

        Args:
            nkpts (int): number of kpoints
        """
        self.arithmetic = toQO_arithmetic()
        self.arithmetic.set_path(path)
        self.arithmetic.initialize(nkpts, m)

    def calculate_qo_hamiltonian(self, supercell: np.ndarray):
        """calculate hqoR
        """
        self.arithmetic.calculate_psi_chi() # get AO projected onto NAO basis
        self.arithmetic.calculate_psi_chi_para() # get AO projected onto eigenstates represented by NAO basis
        self.arithmetic.calculate_psi_chi_orth() # get AO projected onto orthogonal complement of eigenstates represented by NAO basis
        self.arithmetic.calculate_wk() # get canonical orthogonalization matrix
        self.arithmetic.calculate_psi_complem() # get newly introduced states (by AO) represented by NAO basis
        self.arithmetic.calculate_psi_exten() # get extended wavefunction (by AO) represented by NAO basis
        self.arithmetic.calculate_omega() # get density matrix transformed onto NAO basis
        self.arithmetic.calculate_hqok() # get hamiltonian matrix in QO basis
        self.arithmetic.calculate_hqoR(supercell) # get hamiltonian matrix in QO basis in R-space

        return self.arithmetic.hqoR
    
if __name__ == "__main__":
    # user settings
    nkpts = 8 # number of kpoints in ABACUS calculation
    m = 2 # number of new states added
    supercell = np.array([0, 0, 0]) # supercell coordinate

    tqo = toQO()
    tqo.initialize(nkpts, m, "./examples/input/")
    # calculate
    hqoR = tqo.calculate_qo_hamiltonian(supercell)
    # save to file
    hqoR_fname = "hqoR" + "_".join([str(i) for i in supercell]) + ".txt"
    np.savetxt(hqoR_fname, hqoR.reshape(-1, hqoR.shape[-1]), fmt="%.8f")
    # imshow
    import matplotlib.pyplot as plt
    plt.imshow(np.real(hqoR))
    plt.colorbar()
    plt.show()