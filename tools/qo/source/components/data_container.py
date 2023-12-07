import numpy as np

class toQO_DataContainer:
    """the one who really stores data, works as a repo
    """
    # hamiltonians
    hk: np.ndarray # hamiltonian matrix in k-space
    hqok: np.ndarray # hamiltonian matrix in QO basis in k-space
    hqoR: np.ndarray # hamiltonian matrix in QO basis in R-space
    # homogeneous overlaps
    sk: np.ndarray # overlap matrix in k-space
    wk: list # overlap matrix in QO basis in k-space, list of np.ndarray
    # heterogeneous overlaps/transformations
    sqok: list # AO-NAO overlap matrix in k-space, list of np.ndarray
    # data
    kpoints: np.ndarray # kpoints direct coordinates
    energies: np.ndarray # energies of wavefunctions
    # dimensions
    nkpts: int # number of kpoints
    nbands: int # number of bands
    nchi: list # number of AOs for each kpoint
    nphi: int # number of NAOs
    # wavefunction in nao representations
    psi_lcao: np.ndarray # wavefunction represented by numerical atomic orbitals
    psi_chi: list # ao represented by numerical atomic orbitals
    psi_chi_para: list # parallel component of QO wavefunction represented by NAOs
    psi_chi_orth: list # orthogonal component of QO wavefunction represented by NAOs
    psi_complem: list # complementary component of QO wavefunction represented by NAOs
    psi_exten: list # extended component of QO wavefunction represented by NAOs
    psi_qo: list # QO represented by NAOs
    
    def __init__(self) -> None:
        pass
