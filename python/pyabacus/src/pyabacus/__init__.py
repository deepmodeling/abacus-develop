from __future__ import annotations

__submodules__ = ["ModuleBase", "ModuleNAO", "hsolver"]

__all__ = list(__submodules__)

def __getattr__(attr):
    if attr == "ModuleBase":
        from ._core import ModuleBase
        return ModuleBase
    elif attr == "ModuleNAO":
        from ._core import ModuleNAO
        return ModuleNAO
    elif attr == "hsolver":
        import pyabacus.hsolver as hsolver
        return hsolver