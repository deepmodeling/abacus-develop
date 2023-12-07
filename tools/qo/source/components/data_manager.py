import numpy as np
import components.data_container as dc
import tools.hamiltonian as ham
import tools.wavefunction as wf
import tools.qo_ovlp as qov

class toQO_DataManager:
    """Manage memory storing the data
    """
    data: dc.toQO_DataContainer
    def __init__(self) -> None:
        self.data = dc.toQO_DataContainer() # Data manager borns with its data

    def read(self, nkpts: int, path: str, band_range: tuple) -> None:
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
        self.data.sqok = qov.parse(self.data.nkpts, path)
        self.data.psi_lcao, self.data.kpoints, self.data.energies = wf.parse(self.data.nkpts, path)
        # band range selection operation
        self.data.psi_lcao = self.data.psi_lcao[:, band_range[0]:band_range[1], :]
        self.data.energies = self.data.energies[:, band_range[0]:band_range[1]]


    def resize(self) -> None:
        """resize the data container according to the data read from files or after basis filtering
        
        Raises:
            ValueError: if AO filtered set cannot span larger space than the eigenvectors of the Hamiltonian
        """
        self.data.nbands = self.data.psi_lcao.shape[1]
        self.data.nchi = [self.data.sqok[ik].shape[0] for ik in range(self.data.nkpts)]
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
        self.data.omega = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
        self.data.hqok = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
        self.data.hqoR = [np.zeros(
            (self.data.nchi[ik], self.data.nchi[ik])) for ik in range(self.data.nkpts)]
