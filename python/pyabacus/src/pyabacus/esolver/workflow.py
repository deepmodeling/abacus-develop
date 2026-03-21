"""
High-level workflow interface for LCAO calculations.

This module provides a Pythonic interface for controlling LCAO calculations
with support for callbacks and breakpoints.
"""

from typing import Callable, Dict, List, Optional, Any
import numpy as np

from .data_types import (
    ChargeData,
    EnergyData,
    HamiltonianData,
    DensityMatrixData,
    SCFResult,
)


class LCAOWorkflow:
    """
    High-level workflow wrapper for LCAO calculations.

    This class provides a Pythonic interface for controlling LCAO calculations
    with support for callbacks at various stages of the SCF loop.

    Parameters
    ----------
    input_dir : str
        Directory containing INPUT, STRU, and other input files
    gamma_only : bool, optional
        Whether to use gamma-only calculation (default: True)

    Example
    -------
    >>> workflow = LCAOWorkflow("./")
    >>> workflow.initialize()
    >>>
    >>> # Register callback for breakpoint before after_scf
    >>> def inspect_state(wf):
    ...     print(f"Energy: {wf.energy.etot}")
    ...     np.save("charge.npy", wf.charge.rho)
    >>>
    >>> workflow.register_callback('before_after_scf', inspect_state)
    >>> result = workflow.run_scf(max_iter=100)
    >>> print(result.summary())
    """

    # Event names for callbacks
    EVENTS = [
        'before_scf',        # Called after before_scf()
        'after_iter',        # Called after each SCF iteration
        'before_after_scf',  # Called before after_scf() - main breakpoint
        'after_scf',         # Called after after_scf()
    ]

    def __init__(self, input_dir: str, gamma_only: bool = True):
        """
        Initialize LCAOWorkflow.

        Parameters
        ----------
        input_dir : str
            Directory containing input files
        gamma_only : bool
            Use gamma-only calculation if True, multi-k if False
        """
        self._input_dir = input_dir
        self._gamma_only = gamma_only
        self._esolver = None
        self._initialized = False
        self._scf_running = False

        # Callback registry
        self._callbacks: Dict[str, List[Callable]] = {
            event: [] for event in self.EVENTS
        }

    def initialize(self) -> None:
        """
        Initialize the calculation.

        This must be called before running any SCF calculations.
        """
        # Import the appropriate ESolver class
        try:
            if self._gamma_only:
                from ._esolver_pack import ESolverLCAO_gamma
                self._esolver = ESolverLCAO_gamma()
            else:
                from ._esolver_pack import ESolverLCAO_multi_k
                self._esolver = ESolverLCAO_multi_k()
        except ImportError as e:
            raise ImportError(
                f"Could not import ESolver bindings: {e}. "
                "Make sure pyabacus is properly installed with ESolver support."
            ) from e

        self._esolver.initialize(self._input_dir)
        self._esolver.before_all_runners()
        self._initialized = True

    def register_callback(self, event: str, callback: Callable[['LCAOWorkflow'], None]) -> None:
        """
        Register a callback function for a specific event.

        Parameters
        ----------
        event : str
            Event name. One of: 'before_scf', 'after_iter', 'before_after_scf', 'after_scf'
        callback : Callable[[LCAOWorkflow], None]
            Callback function that takes the workflow instance as argument.
            For 'after_iter', the callback receives (workflow, iter_num).

        Raises
        ------
        ValueError
            If event name is not recognized
        """
        if event not in self.EVENTS:
            raise ValueError(
                f"Unknown event '{event}'. Valid events: {self.EVENTS}"
            )
        self._callbacks[event].append(callback)

    def unregister_callback(self, event: str, callback: Callable) -> bool:
        """
        Unregister a callback function.

        Parameters
        ----------
        event : str
            Event name
        callback : Callable
            Callback function to remove

        Returns
        -------
        bool
            True if callback was found and removed, False otherwise
        """
        if event in self._callbacks and callback in self._callbacks[event]:
            self._callbacks[event].remove(callback)
            return True
        return False

    def clear_callbacks(self, event: Optional[str] = None) -> None:
        """
        Clear all callbacks for an event, or all events if event is None.

        Parameters
        ----------
        event : str, optional
            Event name. If None, clears all callbacks.
        """
        if event is None:
            for e in self.EVENTS:
                self._callbacks[e].clear()
        elif event in self._callbacks:
            self._callbacks[event].clear()

    def _fire_callbacks(self, event: str, *args, **kwargs) -> None:
        """Fire all callbacks for an event."""
        for callback in self._callbacks[event]:
            callback(self, *args, **kwargs)

    def run_scf(
        self,
        max_iter: int = 100,
        istep: int = 0,
        callback: Optional[Callable[['LCAOWorkflow', int], None]] = None
    ) -> SCFResult:
        """
        Run SCF calculation with callback support.

        Parameters
        ----------
        max_iter : int
            Maximum number of SCF iterations
        istep : int
            Ion step index (for MD/relaxation)
        callback : Callable, optional
            Additional callback called after each iteration.
            Receives (workflow, iter_num) as arguments.

        Returns
        -------
        SCFResult
            Result of the SCF calculation

        Raises
        ------
        RuntimeError
            If workflow is not initialized
        """
        if not self._initialized:
            raise RuntimeError("Workflow not initialized. Call initialize() first.")

        self._scf_running = True

        # before_scf
        self._esolver.before_scf(istep)
        self._fire_callbacks('before_scf')

        # SCF loop
        for iter_num in range(1, max_iter + 1):
            self._esolver.run_scf_iteration(iter_num)

            # Fire after_iter callbacks
            self._fire_callbacks('after_iter', iter_num)

            # Call user-provided callback
            if callback is not None:
                callback(self, iter_num)

            # Check convergence
            if self._esolver.is_converged():
                break

        # Breakpoint before after_scf - this is the main inspection point
        self._fire_callbacks('before_after_scf')

        # Collect result before after_scf
        result = self._collect_result()

        # after_scf
        self._esolver.after_scf(istep)
        self._fire_callbacks('after_scf')

        self._scf_running = False

        return result

    def run_scf_step(self, iter_num: int) -> None:
        """
        Run a single SCF iteration.

        This is useful for manual control of the SCF loop.

        Parameters
        ----------
        iter_num : int
            Iteration number (1-based)
        """
        if not self._scf_running:
            raise RuntimeError(
                "SCF not started. Call before_scf() first or use run_scf()."
            )
        self._esolver.run_scf_iteration(iter_num)

    def before_scf(self, istep: int = 0) -> None:
        """
        Prepare for SCF calculation.

        Call this before manually running SCF iterations.

        Parameters
        ----------
        istep : int
            Ion step index
        """
        if not self._initialized:
            raise RuntimeError("Workflow not initialized. Call initialize() first.")
        self._esolver.before_scf(istep)
        self._scf_running = True
        self._fire_callbacks('before_scf')

    def after_scf(self, istep: int = 0) -> None:
        """
        Finalize SCF calculation.

        Call this after manually running SCF iterations.

        Parameters
        ----------
        istep : int
            Ion step index
        """
        self._fire_callbacks('before_after_scf')
        self._esolver.after_scf(istep)
        self._fire_callbacks('after_scf')
        self._scf_running = False

    def _collect_result(self) -> SCFResult:
        """Collect SCF result from current state."""
        return SCFResult(
            converged=self._esolver.is_converged(),
            niter=self._esolver.niter,
            drho=self._esolver.drho,
            energy=self.energy,
            charge=self.charge,
        )

    # ==================== Properties for data access ====================

    @property
    def charge(self) -> ChargeData:
        """
        Get current charge density.

        Returns
        -------
        ChargeData
            Charge density data container
        """
        accessor = self._esolver.get_charge()
        if not accessor.is_valid():
            return ChargeData(rho=np.array([]), nspin=0, nrxx=0)

        return ChargeData(
            rho=accessor.get_rho(),
            nspin=accessor.nspin,
            nrxx=accessor.nrxx,
        )

    @property
    def energy(self) -> EnergyData:
        """
        Get current energy data.

        Returns
        -------
        EnergyData
            Energy data container with all energy components
        """
        accessor = self._esolver.get_energy()
        return EnergyData(
            etot=accessor.etot,
            eband=accessor.eband,
            hartree_energy=accessor.hartree_energy,
            etxc=accessor.etxc,
            ewald_energy=accessor.ewald_energy,
            demet=accessor.demet,
            exx=accessor.exx,
            evdw=accessor.evdw,
        )

    @property
    def hamiltonian(self) -> HamiltonianData:
        """
        Get current Hamiltonian matrices.

        Returns
        -------
        HamiltonianData
            Hamiltonian data container with H(k), S(k), H(R), S(R)
        """
        accessor = self._esolver.get_hamiltonian()
        if not accessor.is_valid():
            return HamiltonianData()

        nks = accessor.nks
        Hk = [accessor.get_Hk(ik) for ik in range(nks)]
        Sk = [accessor.get_Sk(ik) for ik in range(nks)]

        return HamiltonianData(
            Hk=Hk,
            Sk=Sk,
            HR=accessor.get_HR(),
            SR=accessor.get_SR(),
            nbasis=accessor.nbasis,
            nks=nks,
        )

    @property
    def density_matrix(self) -> DensityMatrixData:
        """
        Get current density matrix.

        Returns
        -------
        DensityMatrixData
            Density matrix data container with DM(k) and DM(R)
        """
        accessor = self._esolver.get_density_matrix()
        if not accessor.is_valid():
            return DensityMatrixData()

        nks = accessor.nks
        DMK = [accessor.get_DMK(ik) for ik in range(nks)]

        return DensityMatrixData(
            DMK=DMK,
            DMR=accessor.get_DMR(),
            nks=nks,
            nrow=accessor.nrow,
            ncol=accessor.ncol,
        )

    @property
    def is_converged(self) -> bool:
        """Check if SCF is converged."""
        return self._esolver.is_converged()

    @property
    def niter(self) -> int:
        """Get current iteration number."""
        return self._esolver.niter

    @property
    def drho(self) -> float:
        """Get current charge density difference."""
        return self._esolver.drho

    @property
    def nks(self) -> int:
        """Get number of k-points."""
        return self._esolver.nks

    @property
    def nbasis(self) -> int:
        """Get number of basis functions."""
        return self._esolver.nbasis

    @property
    def nbands(self) -> int:
        """Get number of bands."""
        return self._esolver.nbands

    @property
    def nspin(self) -> int:
        """Get number of spin channels."""
        return self._esolver.nspin

    def get_psi(self, ik: int) -> np.ndarray:
        """
        Get wave function coefficients for k-point ik.

        Parameters
        ----------
        ik : int
            K-point index

        Returns
        -------
        np.ndarray
            Wave function coefficients with shape (nbands, nbasis)
        """
        return self._esolver.get_psi(ik)

    def get_eigenvalues(self, ik: int) -> np.ndarray:
        """
        Get eigenvalues for k-point ik.

        Parameters
        ----------
        ik : int
            K-point index

        Returns
        -------
        np.ndarray
            Eigenvalues with shape (nbands,)
        """
        return self._esolver.get_eigenvalues(ik)

    def get_occupations(self, ik: int) -> np.ndarray:
        """
        Get occupation numbers for k-point ik.

        Parameters
        ----------
        ik : int
            K-point index

        Returns
        -------
        np.ndarray
            Occupation numbers with shape (nbands,)
        """
        return self._esolver.get_occupations(ik)

    def get_kvec(self, ik: int) -> np.ndarray:
        """
        Get k-vector in direct coordinates.

        Parameters
        ----------
        ik : int
            K-point index

        Returns
        -------
        np.ndarray
            K-vector with shape (3,)
        """
        return self._esolver.get_kvec_d(ik)

    def get_kweights(self) -> np.ndarray:
        """
        Get k-point weights.

        Returns
        -------
        np.ndarray
            K-point weights with shape (nks,)
        """
        return self._esolver.get_wk()
