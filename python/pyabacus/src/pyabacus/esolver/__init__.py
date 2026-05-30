"""
PyABACUS ESolver Module
=======================

This module provides Python bindings for ABACUS ESolver_KS_LCAO,
enabling Python-controlled SCF workflows with breakpoint support.

Main Classes
------------
ESolverLCAO_gamma : ESolver for gamma-only calculations
ESolverLCAO_multi_k : ESolver for multi-k calculations
LCAOWorkflow : High-level workflow wrapper with callback support

Example
-------
>>> from pyabacus.esolver import LCAOWorkflow
>>>
>>> workflow = LCAOWorkflow("./")
>>> workflow.initialize()
>>>
>>> # Register callback for breakpoint before after_scf
>>> def save_state(wf):
...     charge = wf.charge
...     energy = wf.energy
...     print(f"Total energy: {energy.etot}")
>>>
>>> workflow.register_callback('before_after_scf', save_state)
>>> result = workflow.run_scf(max_iter=100)
"""

from .workflow import LCAOWorkflow
from .data_types import ChargeData, EnergyData, HamiltonianData, DensityMatrixData, SCFResult

# Import C++ bindings
try:
    from ._esolver_pack import (
        ESolverLCAO_gamma,
        ESolverLCAO_multi_k,
        ChargeAccessor,
        EnergyAccessor,
        HamiltonianAccessor_gamma,
        HamiltonianAccessor_multi_k,
        DensityMatrixAccessor_gamma,
        DensityMatrixAccessor_multi_k,
    )
except ImportError as e:
    import warnings
    warnings.warn(f"Could not import _esolver_pack: {e}. "
                  "ESolver bindings may not be available.")

    # Define placeholder classes for documentation
    ESolverLCAO_gamma = None
    ESolverLCAO_multi_k = None
    ChargeAccessor = None
    EnergyAccessor = None
    HamiltonianAccessor_gamma = None
    HamiltonianAccessor_multi_k = None
    DensityMatrixAccessor_gamma = None
    DensityMatrixAccessor_multi_k = None

__all__ = [
    # High-level interface
    'LCAOWorkflow',

    # Data types
    'ChargeData',
    'EnergyData',
    'HamiltonianData',
    'DensityMatrixData',
    'SCFResult',

    # Low-level C++ bindings
    'ESolverLCAO_gamma',
    'ESolverLCAO_multi_k',
    'ChargeAccessor',
    'EnergyAccessor',
    'HamiltonianAccessor_gamma',
    'HamiltonianAccessor_multi_k',
    'DensityMatrixAccessor_gamma',
    'DensityMatrixAccessor_multi_k',
]
