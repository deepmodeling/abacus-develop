import numpy as np
import source.components.data_container as dc
import source.tools.hamiltonian as ham
import source.tools.wavefunction as wf
import source.tools.qo_ovlp as qov
import source.tools.kpoints as kpt

class toQO_DataManager:
    """Manage memory storing the data
    """
    data: dc.toQO_DataContainer
    def __init__(self) -> None:
        self.data = dc.toQO_DataContainer() # Data manager borns with its data

    def kpoints_eq(self, kpt1: np.ndarray, kpt2: np.ndarray):
        """check if two kpoints are equal"""
        return np.linalg.norm(kpt1 - kpt2) < 1e-6

    def align_according_to_kpoints_ref(self, to_align: list, kpts: list, kpts_ref: list):
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
        # align
        aligned = []
        for ik in range(nkpts):
            aligned.append(to_align[indexing[ik]])

        return aligned

    def check_kpoints_alignment(self, kpts1: list, kpts2: list):
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
        self.data.nkpts = nkpts
        self.data.hk, self.data.sk, self.data.nphi = ham.parse(self.data.nkpts, path)
        self.data.saok, saok_kpts = qov.parse(self.data.nkpts, path)
        self.data.psi_lcao, psi_kpts, self.data.energies = wf.parse(self.data.nkpts, path)
        kpts, self.data.equivalent_kpoints = kpt.parse(path)
        # presently there may be at least two sets of kpoints in ABACUS
        # Hamiltonian and Overlap matrix -> kpoints file
        # Overlap AO with NAO            -> kpoints file
        # wavefunction                   -> itself
        # align all according to kpoints file
        if not self.check_kpoints_alignment(kpts, psi_kpts):
            self.data.psi_lcao = self.align_according_to_kpoints_ref(self.data.psi_lcao, psi_kpts, kpts)
            self.data.energies = self.align_according_to_kpoints_ref(self.data.energies, psi_kpts, kpts)
        if not self.check_kpoints_alignment(kpts, saok_kpts):
            self.data.saok = self.align_according_to_kpoints_ref(self.data.saok, saok_kpts, kpts)
        self.data.kpoints = kpts
        # band range selection operation
        self.data.psi_lcao = self.data.psi_lcao[:, band_range[0]:band_range[1], :]
        self.data.energies = self.data.energies[:, band_range[0]:band_range[1]]

    def resize(self) -> None:
        """resize the data container according to the data read from files or after basis filtering
        
        Raises:
            ValueError: if AO filtered set cannot span larger space than the eigenvectors of the Hamiltonian
        """
        self.data.nbands = self.data.psi_lcao.shape[1]
        self.data.nchi = [self.data.saok[ik].shape[0] for ik in range(self.data.nkpts)]
        self.data.nphi = self.data.hk[0].shape[0]

        _m = [self.data.nchi[ik] - self.data.nbands for ik in range(self.data.nkpts)]
        for ik in range(self.data.nkpts):
            if _m[ik] < 0:
                raise ValueError("current filtered AO set cannot span larger space than the eigenvectors of the Hamiltonian")
            
        self.data.psi_chi = [np.zeros(
            (self.data.nchi[ik], self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_qo = [np.zeros(
            (self.data.nchi[ik], self.data.nphi)) for ik in range(self.data.nkpts)] # how many AOs are there, how many QOs there are.
        self.data.psi_chi_para = [np.zeros(
            (self.data.nchi[ik], self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_chi_orth = [np.zeros(
            (self.data.nchi[ik], self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.wk = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
        self.data.psi_complem = [np.zeros(
            (_m[ik], self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.psi_exten = [np.zeros(
            (self.data.nchi[ik], self.data.nphi)) for ik in range(self.data.nkpts)]
        self.data.hqok = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
        self.data.sqok = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
        