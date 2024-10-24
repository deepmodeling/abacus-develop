"""
pyabacus.ModuleBase
===================
Basic math functions and integrals.
"""

import numpy as np
from typing import overload
from numpy.typing import NDArray

from ._base_pack import Sphbes as _Sphbes, Integral as _Integral, SphericalBesselTransformer as _SphericalBesselTransformer

class Sphbes(_Sphbes):
    def __init__(self) -> None: 
        super().__init__()
        
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
    
    def sphbesj(self, *args, **kwargs): 
        return super().sphbesj(*args, **kwargs)
        
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
    
    def dsphbesj(self, *args, **kwargs):
        return super().dsphbesj(*args, **kwargs)
        
    @staticmethod
    def sphbes_zeros(l: int, n: int, zeros: NDArray[np.float64]) -> None: 
        super().sphbes_zeros(l, n, zeros)

class Integral(_Integral):
    def __init__(self) -> None: 
        super().__init__()
        
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
    
    def Simpson_Integral(self, *args, **kwargs): 
        return super().Simpson_Integral(*args, **kwargs)
    
    @staticmethod
    def Simpson_Integral_0toall(
        mesh: int, 
        func: NDArray[np.float64], 
        rab: NDArray[np.float64], 
        asum: NDArray[np.float64]
    ) -> None: 
        super().Simpson_Integral_0toall(mesh, func, rab, asum)
        
    @staticmethod
    def Simpson_Integral_alltoinf(
        mesh: int, 
        func: NDArray[np.float64], 
        rab: NDArray[np.float64], 
        asum: NDArray[np.float64]
    ) -> None: 
        super().Simpson_Integral_alltoinf(mesh, func, rab, asum)
        
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
    
    def simpson(self, *args, **kwargs):
        return super().simpson(*args, **kwargs)
    
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
    
    def Gauss_Legendre_grid_and_weight(self, *args, **kwargs):
        return super().Gauss_Legendre_grid_and_weight(*args, **kwargs)

class SphericalBesselTransformer(_SphericalBesselTransformer):
    def __init__(self) -> None: 
        super().__init__()
