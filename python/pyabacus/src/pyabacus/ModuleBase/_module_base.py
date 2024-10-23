"""
pyabacus.ModuleBase
===================
Basic math functions and integrals.
"""

import numpy as np
from typing import Protocol, overload
from numpy.typing import NDArray

from ._base_pack import Sphbes as _Sphbes, Integral as _Integral, SphericalBesselTransformer as _SphericalBesselTransformer

class SphbesProtocol(Protocol):
    @overload
    @staticmethod
    def sphbesj(l: int, x: float) -> float: ...
    @overload
    @staticmethod
    def sphbesj(
        n: int, 
        r: NDArray[np.float64], 
        q: int, 
        l: int, 
        jl: NDArray[np.float64]
    ) -> None: ...
    @overload
    @staticmethod
    def dsphbesj(l: int, x: float) -> float: ...
    @overload
    @staticmethod
    def dsphbesj(
        n: int, 
        r: NDArray[np.float64], 
        q: int, 
        l: int, 
        djl: NDArray[np.float64]
    ) -> None: ...
    @staticmethod
    def sphbes_zeros(l: int, n: int, zeros: NDArray[np.float64]) -> None: ...

class IntegralProtocol(Protocol):
    @overload
    @staticmethod
    def Simpson_Integral(
        mesh: int, 
        func: NDArray[np.float64], 
        rab: NDArray[np.float64], 
        asum: float
    ) -> float: ...
    @overload
    @staticmethod
    def Simpson_Integral(
        mesh: int, 
        func: NDArray[np.float64], 
        dr: float,
        asum: float
    ) -> float: ...
    @staticmethod
    def Simpson_Integral_0toall(
        mesh: int, 
        func: NDArray[np.float64], 
        rab: NDArray[np.float64], 
        asum: NDArray[np.float64]
    ) -> None: ...
    @staticmethod
    def Simpson_Integral_alltoinf(
        mesh: int, 
        func: NDArray[np.float64], 
        rab: NDArray[np.float64], 
        asum: NDArray[np.float64]
    ) -> None: ...
    @overload
    @staticmethod
    def simpson(
        n: int,
        f: NDArray[np.float64],
        dx: float
    ) -> float: ...
    @overload
    @staticmethod
    def simpson(
        n: int,
        f: NDArray[np.float64],
        h: NDArray[np.float64],
    ) -> float: ...
    @overload
    @staticmethod
    def Gauss_Legendre_grid_and_weight(
        n: int,
        x: NDArray[np.float64],
        w: NDArray[np.float64],
    ) -> None: ...
    @overload
    @staticmethod
    def Gauss_Legendre_grid_and_weight(
        xmin: float,
        xmax: float,
        n: int,
        x: NDArray[np.float64],
        w: NDArray[np.float64],
    ) -> None: ...

class SphericalBesselTransformerProtocol(Protocol):
    def __init__(self) -> None: ...

Sphbes: SphbesProtocol = _Sphbes  
Integral: IntegralProtocol = _Integral  
SphericalBesselTransformer: SphericalBesselTransformerProtocol = _SphericalBesselTransformer  