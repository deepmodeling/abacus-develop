"""
Data types for PyABACUS ESolver module.

This module defines dataclasses for storing calculation results
in a structured and type-safe manner.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
import numpy as np


@dataclass
class ChargeData:
    """
    Container for charge density data.

    Attributes
    ----------
    rho : np.ndarray
        Real-space charge density with shape (nspin, nrxx)
    rhog : np.ndarray, optional
        Reciprocal-space charge density with shape (nspin, ngmc)
    nspin : int
        Number of spin channels (1, 2, or 4)
    nrxx : int
        Number of real-space grid points
    ngmc : int, optional
        Number of G-vectors for charge density
    """
    rho: np.ndarray
    nspin: int
    nrxx: int
    rhog: Optional[np.ndarray] = None
    ngmc: Optional[int] = None

    def total_charge(self) -> float:
        """Calculate total charge by integrating rho."""
        return np.sum(self.rho)

    def spin_density(self) -> Optional[np.ndarray]:
        """
        Calculate spin density (rho_up - rho_down) for spin-polarized calculations.

        Returns None for non-spin-polarized calculations.
        """
        if self.nspin == 2:
            return self.rho[0] - self.rho[1]
        return None


@dataclass
class EnergyData:
    """
    Container for energy data.

    All energies are in Rydberg units.

    Attributes
    ----------
    etot : float
        Total energy
    eband : float
        Band (kinetic + local potential) energy
    hartree_energy : float
        Hartree (electron-electron Coulomb) energy
    etxc : float
        Exchange-correlation energy
    ewald_energy : float
        Ewald (ion-ion Coulomb) energy
    demet : float
        -TS term for metallic systems (smearing correction)
    exx : float
        Exact exchange energy (for hybrid functionals)
    evdw : float
        van der Waals correction energy
    """
    etot: float = 0.0
    eband: float = 0.0
    hartree_energy: float = 0.0
    etxc: float = 0.0
    ewald_energy: float = 0.0
    demet: float = 0.0
    exx: float = 0.0
    evdw: float = 0.0

    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary."""
        return {
            'etot': self.etot,
            'eband': self.eband,
            'hartree_energy': self.hartree_energy,
            'etxc': self.etxc,
            'ewald_energy': self.ewald_energy,
            'demet': self.demet,
            'exx': self.exx,
            'evdw': self.evdw,
        }

    def to_eV(self) -> 'EnergyData':
        """
        Convert all energies from Rydberg to eV.

        Returns a new EnergyData instance with energies in eV.
        """
        Ry_to_eV = 13.605693122994  # 1 Ry = 13.6057 eV
        return EnergyData(
            etot=self.etot * Ry_to_eV,
            eband=self.eband * Ry_to_eV,
            hartree_energy=self.hartree_energy * Ry_to_eV,
            etxc=self.etxc * Ry_to_eV,
            ewald_energy=self.ewald_energy * Ry_to_eV,
            demet=self.demet * Ry_to_eV,
            exx=self.exx * Ry_to_eV,
            evdw=self.evdw * Ry_to_eV,
        )


@dataclass
class HamiltonianData:
    """
    Container for Hamiltonian matrix data.

    Attributes
    ----------
    Hk : List[np.ndarray]
        List of H(k) matrices for each k-point
    Sk : List[np.ndarray]
        List of S(k) overlap matrices for each k-point
    HR : Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray], optional
        H(R) in sparse format: {(iat1, iat2, (R1, R2, R3)): matrix}
    SR : Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray], optional
        S(R) in sparse format: {(iat1, iat2, (R1, R2, R3)): matrix}
    nbasis : int
        Number of basis functions
    nks : int
        Number of k-points
    """
    Hk: List[np.ndarray] = field(default_factory=list)
    Sk: List[np.ndarray] = field(default_factory=list)
    HR: Optional[Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray]] = None
    SR: Optional[Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray]] = None
    nbasis: int = 0
    nks: int = 0

    def get_Hk(self, ik: int) -> np.ndarray:
        """Get H(k) matrix for k-point ik."""
        if ik < 0 or ik >= len(self.Hk):
            raise IndexError(f"K-point index {ik} out of range [0, {len(self.Hk)})")
        return self.Hk[ik]

    def get_Sk(self, ik: int) -> np.ndarray:
        """Get S(k) matrix for k-point ik."""
        if ik < 0 or ik >= len(self.Sk):
            raise IndexError(f"K-point index {ik} out of range [0, {len(self.Sk)})")
        return self.Sk[ik]


@dataclass
class DensityMatrixData:
    """
    Container for density matrix data.

    Attributes
    ----------
    DMK : List[np.ndarray]
        List of DM(k) matrices for each k-point
    DMR : Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray], optional
        DM(R) in sparse format: {(iat1, iat2, (R1, R2, R3)): matrix}
    nks : int
        Number of k-points
    nrow : int
        Number of rows in density matrix
    ncol : int
        Number of columns in density matrix
    """
    DMK: List[np.ndarray] = field(default_factory=list)
    DMR: Optional[Dict[Tuple[int, int, Tuple[int, int, int]], np.ndarray]] = None
    nks: int = 0
    nrow: int = 0
    ncol: int = 0

    def get_DMK(self, ik: int) -> np.ndarray:
        """Get DM(k) matrix for k-point ik."""
        if ik < 0 or ik >= len(self.DMK):
            raise IndexError(f"K-point index {ik} out of range [0, {len(self.DMK)})")
        return self.DMK[ik]

    def trace(self, ik: int) -> complex:
        """Calculate trace of DM(k) for k-point ik."""
        return np.trace(self.get_DMK(ik))


@dataclass
class SCFResult:
    """
    Container for SCF calculation results.

    Attributes
    ----------
    converged : bool
        Whether SCF converged
    niter : int
        Number of iterations performed
    drho : float
        Final charge density difference
    energy : EnergyData
        Final energy data
    charge : ChargeData, optional
        Final charge density
    hamiltonian : HamiltonianData, optional
        Final Hamiltonian matrices
    density_matrix : DensityMatrixData, optional
        Final density matrix
    """
    converged: bool
    niter: int
    drho: float
    energy: EnergyData
    charge: Optional[ChargeData] = None
    hamiltonian: Optional[HamiltonianData] = None
    density_matrix: Optional[DensityMatrixData] = None

    def summary(self) -> str:
        """Return a summary string of the SCF result."""
        status = "converged" if self.converged else "not converged"
        return (
            f"SCF Result: {status}\n"
            f"  Iterations: {self.niter}\n"
            f"  Final drho: {self.drho:.2e}\n"
            f"  Total energy: {self.energy.etot:.8f} Ry "
            f"({self.energy.etot * 13.6057:.8f} eV)"
        )
