import numpy as np

class toQO_DataContainer:
    """the one who really stores data, works as a repo
    """
    # hamiltonians
    hk: list # hamiltonian matrix in k-space
    hqok: np.ndarray # hamiltonian matrix in QO basis in k-space
    sqok: np.ndarray # overlap matrix in QO basis in k-space
    # homogeneous overlaps
    sk: list # overlap matrix in k-space
    wk: list # overlap matrix in QO basis in k-space, list of np.ndarray
    # heterogeneous overlaps/transformations
    saok: list # AO-NAO overlap matrix in k-space, list of np.ndarray
    # data
    kpoints: list # kpoints direct coordinates
    equivalent_kpoints: list # weights of kpoints
    energies: list # energies of wavefunctions
    # dimensions
    nkpts: int # number of kpoints
    nbands: int # number of bands
    nchi: int # number of AOs
    nphi: int # number of NAOs
    # wavefunction in nao representations
    psi_lcao: list # wavefunction represented by numerical atomic orbitals
    psi_chi: list # ao represented by numerical atomic orbitals
    psi_chi_para: list # parallel component of QO wavefunction represented by NAOs
    psi_chi_orth: list # orthogonal component of QO wavefunction represented by NAOs
    psi_complem: list # complementary component of QO wavefunction represented by NAOs
    psi_exten: list # extended component of QO wavefunction represented by NAOs
    psi_qo: list # QO represented by NAOs
    
    def __init__(self) -> None:
        pass
