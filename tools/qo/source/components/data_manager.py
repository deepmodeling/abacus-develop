import numpy as np
import source.components.data_container as dc
import source.tools.hamiltonian as ham
import source.tools.wavefunction as wf
import source.tools.qo_ovlp as qov
import source.tools.kpoints as kpt
import source.components.safe_guard as sg
from scipy.linalg import eigh

class toQO_DataManager:
    """Manage memory storing the data
    """
    data: dc.toQO_DataContainer
    sg_: sg.toQO_SafeGuard

    def __init__(self) -> None:
        self.data = dc.toQO_DataContainer() # Data manager borns with its data
        self.sg_ = sg.toQO_SafeGuard()

    def kpoints_eq(self, kpt1: np.ndarray, kpt2: np.ndarray) -> bool:
        """check if two kpoints are equal"""
        return np.linalg.norm(kpt1 - kpt2) < 1e-6

    def align_according_to_kpoints_ref(self, to_align: list, kpts: list, kpts_ref: list) -> list:
        """align kpoints according to kpoints file"""
        indexing = []
        # first check if number of kpoints are the same
        if len(kpts) != len(kpts_ref):
            raise ValueError("kpoints read from files are not the same")
        nkpts = len(kpts_ref)
        for ik in range(nkpts):
            _kpt = kpts_ref[ik]
            for jk in range(nkpts):
                if self.kpoints_eq(_kpt, kpts_ref[jk]):
                    indexing.append(jk)
                    break
        # indexing has the meaning: for the first kpoint in kpts_ref, the data would be found in the indexing[0]th position of to_align
        # align
        aligned = []
        for ik in range(nkpts):
            aligned.append(to_align[indexing[ik]])

        return aligned

    def check_kpoints_alignment(self, kpts1: list, kpts2: list) -> bool:
        """check if two kpoints lists are aligned"""
        nkpts = len(kpts1)
        for ik in range(nkpts):
            if self.kpoints_eq(kpts1[ik], kpts2[ik]):
                return False
        return True

    def read(self, nkpts: int, calculation: str, path: str, band_range: tuple) -> None:
        """read data from files without any changes
        
        Args:
        
            nkpts (int): number of kpoints
            path (str): path to the data
            band_range (tuple): min and max band index
            
        Raises:
        
            ValueError: None in this level of function
        """
        print("""
---------------------------------------------------------------------------------------------------------
$$$$$$$\   $$$$$$\ $$$$$$$$\  $$$$$$\        $$$$$$\ $$\      $$\ $$$$$$$\   $$$$$$\  $$$$$$$\ $$$$$$$$\ 
$$  __$$\ $$  __$$\\__$$  __|$$  __$$\       \_$$  _|$$$\    $$$ |$$  __$$\ $$  __$$\ $$  __$$\\__$$  __|
$$ |  $$ |$$ /  $$ |  $$ |   $$ /  $$ |        $$ |  $$$$\  $$$$ |$$ |  $$ |$$ /  $$ |$$ |  $$ |  $$ |   
$$ |  $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |        $$ |  $$\$$\$$ $$ |$$$$$$$  |$$ |  $$ |$$$$$$$  |  $$ |   
$$ |  $$ |$$  __$$ |  $$ |   $$  __$$ |        $$ |  $$ \$$$  $$ |$$  ____/ $$ |  $$ |$$  __$$<   $$ |   
$$ |  $$ |$$ |  $$ |  $$ |   $$ |  $$ |        $$ |  $$ |\$  /$$ |$$ |      $$ |  $$ |$$ |  $$ |  $$ |   
$$$$$$$  |$$ |  $$ |  $$ |   $$ |  $$ |      $$$$$$\ $$ | \_/ $$ |$$ |       $$$$$$  |$$ |  $$ |  $$ |   
\_______/ \__|  \__|  \__|   \__|  \__|      \______|\__|     \__|\__|       \______/ \__|  \__|  \__|   
---------------------------------------------------------------------------------------------------------
In this step, read data dumped from ABACUS. Please ensure you have enough precision on
printing several critical quantities, such as H(k), S(k), wavefunction and Sao(k).\n""")
        # int
        self.data.nkpts = nkpts
        # list of np.ndarray(2d), list of np.ndarray(2d), int
        self.data.hk, self.data.sk, self.data.nphi = ham.parse(self.data.nkpts, path)
        # list of np.ndarray(2d), list of np.ndarray(nkpts, 3)
        self.data.saok, saok_kpts = qov.parse(self.data.nkpts, path)
        # list of np.ndarray(2d), list of np.ndarray(nkpts, 3), list of np.ndarray(1d)
        self.data.psi_lcao, psi_kpts, self.data.energies = wf.parse(self.data.nkpts, path)
        # list of np.ndarray(nkpts, 3), list of list of np.ndarray(nkpts, 3)
        kpts, self.data.equivalent_kpoints = kpt.parse(path)
        # presently there may be at least two sets of kpoints in ABACUS
        # Hamiltonian and Overlap matrix -> kpoints file, kvec_d
        # Overlap AO with NAO            -> kpoints file, kvec_d
        # wavefunction                   -> itself, kvec_c
        # align all according to kpoints file
        # if not self.check_kpoints_alignment(kpts, psi_kpts):
        #     self.data.psi_lcao = self.align_according_to_kpoints_ref(self.data.psi_lcao, psi_kpts, kpts)
        #     self.data.energies = self.align_according_to_kpoints_ref(self.data.energies, psi_kpts, kpts)
        # if not self.check_kpoints_alignment(kpts, saok_kpts):
        #     self.data.saok = self.align_according_to_kpoints_ref(self.data.saok, saok_kpts, kpts)
        self.data.kpoints = kpts
        print("""
-----------------------------------------------------------------------------------------------
$$$$$$$\   $$$$$$\ $$$$$$$$\  $$$$$$\         $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\  $$\   $$\ 
$$  __$$\ $$  __$$\\__$$  __|$$  __$$\       $$  __$$\ $$ |  $$ |$$  _____|$$  __$$\ $$ | $$  |
$$ |  $$ |$$ /  $$ |  $$ |   $$ /  $$ |      $$ /  \__|$$ |  $$ |$$ |      $$ /  \__|$$ |$$  / 
$$ |  $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |      $$ |      $$$$$$$$ |$$$$$\    $$ |      $$$$$  /  
$$ |  $$ |$$  __$$ |  $$ |   $$  __$$ |      $$ |      $$  __$$ |$$  __|   $$ |      $$  $$<   
$$ |  $$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |  $$\ $$ |  $$ |$$ |      $$ |  $$\ $$ |\$$\  
$$$$$$$  |$$ |  $$ |  $$ |   $$ |  $$ |      \$$$$$$  |$$ |  $$ |$$$$$$$$\ \$$$$$$  |$$ | \$$\ 
\_______/ \__|  \__|  \__|   \__|  \__|       \______/ \__|  \__|\________| \______/ \__|  \__|
-----------------------------------------------------------------------------------------------
Due to interface between C++ and Python will lose precision, in this step, 
orthonormalities between eigenstates are needed to check.""")
        for ik in range(self.data.nkpts):
            print("-"*50, "\nFor k-point No.", ik)
            error = self.sg_.check_eigenstates(self.data.psi_lcao[ik].T, self.data.sk[ik])
            if error > 1e-6:
                print("Warning: eigenstates are not orthonormal, re-diagonalize Hamiltonian and overlap matrix")
                eigvals, eigvecs = eigh(self.data.hk[ik], self.data.sk[ik])
                self.data.psi_lcao[ik] = eigvecs.T
                self.data.energies[ik] = eigvals
        print("-"*50, "\nOrthonormalities check finished. \nThen band selection according to user settings.")
        #                                    kpoints          bands          orbitals
        self.data.psi_lcao = [np.array(psi_k)[band_range[0]:band_range[1], :] for psi_k in self.data.psi_lcao]
        self.data.energies = [np.array(energies_k)[band_range[0]:band_range[1]] for energies_k in self.data.energies]

    def resize(self) -> None:
        """resize the data container according to the data read from files or after basis filtering
        
        Raises:
            ValueError: if AO filtered set cannot span larger space than the eigenvectors of the Hamiltonian
        """
        self.data.nbands = self.data.psi_lcao[0].shape[0]
        self.data.nchi = self.data.saok[0].shape[0]
        self.data.nphi = self.data.hk[0].shape[0]

        _m = self.data.nchi - self.data.psi_lcao[0].shape[0]
        if _m < 0:
            print("number of AOs:", self.data.nchi)
            print("number of bands selected in solved occupied eigenstates:", self.data.psi_lcao[0].shape[0])
            raise ValueError("current filtered AO set cannot span larger space than the eigenvectors of the Hamiltonian")
            
        self.data.psi_chi = [np.zeros(
            (self.data.nchi, self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_qo = [np.zeros(
            (self.data.nchi, self.data.nphi)) for ik in range(self.data.nkpts)] # how many AOs are there, how many QOs there are.
        self.data.psi_chi_para = [np.zeros(
            (self.data.nchi, self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_chi_orth = [np.zeros(
            (self.data.nchi, self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.wk = [np.zeros(
            (self.data.nchi, self.data.nchi)) for ik in range(self.data.nkpts)]
        self.data.psi_complem = [np.zeros(
            (_m, self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_exten = [np.zeros(
            (self.data.nchi, self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.hqok = [np.zeros(
            (self.data.nchi, self.data.nchi)) for ik in range(self.data.nkpts)]
        self.data.sqok = [np.zeros(
            (self.data.nchi, self.data.nchi)) for ik in range(self.data.nkpts)]
        