from __future__ import annotations

__submodules__ = ["ModuleBase", "ModuleNAO", "hsolver", "Cell", "IntegralCalculator", "io", "esolver", "driver"]
__all__ = list(__submodules__) + ["abacus", "CalculationResult"]

# Import the main abacus() function for convenience
def __getattr__(attr):
    if attr == "ModuleBase":
        import pyabacus.ModuleBase as ModuleBase
        return ModuleBase
    elif attr == "ModuleNAO":
        import pyabacus.ModuleNAO as ModuleNAO
        return ModuleNAO
    elif attr == "hsolver":
        import pyabacus.hsolver as hsolver
        return hsolver
    elif attr == "Cell":
        from .cell import Cell
        return Cell
    elif attr == "io":
        import pyabacus.io as io
        return io
    elif attr == "esolver":
        import pyabacus.esolver as esolver
        return esolver
    elif attr == "driver":
        import pyabacus.driver as driver
        return driver
    elif attr == "abacus":
        from pyabacus.driver import abacus
        return abacus
    elif attr == "CalculationResult":
        from pyabacus.driver import CalculationResult
        return CalculationResult
    else:
        raise AttributeError(f"module {__name__} has no attribute {attr}")