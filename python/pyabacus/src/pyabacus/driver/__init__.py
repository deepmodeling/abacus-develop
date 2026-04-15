"""
PyABACUS Driver Module
======================

This module provides the high-level `abacus()` function for running
complete ABACUS DFT calculations from Python.

Example
-------
>>> import pyabacus
>>> result = pyabacus.abacus("./Si_scf/")
>>> print(f"Energy: {result.etot_ev:.6f} eV")
>>> print(f"Converged: {result.converged}")
"""

from .runner import abacus, CalculationResult

# Try to import C++ bindings
try:
    from ._driver_pack import PyDriver, CalculationResult as _CppResult
    _HAS_CPP_DRIVER = True
except ImportError:
    _HAS_CPP_DRIVER = False
    PyDriver = None
    _CppResult = None

__all__ = [
    'abacus',
    'CalculationResult',
    'PyDriver',
]
