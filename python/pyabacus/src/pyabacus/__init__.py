from __future__ import annotations

__submodules__ = ["ModuleBase", "ModuleNAO", "hsolver"]

__all__ = list(__submodules__)

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
    else:
        raise AttributeError(f"module {__name__} has no attribute {attr}")