"""
pyabacus.ModuleNAO
==================
Module for Numerical Atomic Orbitals (NAO) in ABACUS
"""

import numpy as np

from numpy.typing import NDArray
from pyabacus.ModuleBase import SphericalBesselTransformer
from typing import overload, List


class RadialCollection:
    def __init__(self) -> None: 
        """
        A class that holds all numerical radial functions of the same kind. 
        
        An instance of this class could be the collection of all radial functions 
        of numerical atomic orbitals, or all Kleinman-Bylander beta functions from 
        all elements involved in a calculation.
        """
        pass
    
    def build(
        self, 
        nfile: int, 
        file_list: List[str], 
        ftype: str = '\0'
    ) -> None:
        """
        Builds the collection from (orbital) files.
        """
        pass
    
    def set_transformer(
        self, 
        sbt: SphericalBesselTransformer, 
        update: int = 0
    ) -> None:
        """
        Sets a spherical Bessel transformers for all RadialSet objects.
        """
        pass
    
    def set_uniform_grid(
        self, 
        for_r_space: bool, 
        ngrid: int, 
        cutoff: float, 
        mode: str = 'i', 
        enable_fft: bool = False
    ) -> None:
        """
        Sets a common uniform grid for all RadialSet objects.
        """
        pass
    
    def set_grid(
        self, 
        for_r_space: bool, 
        ngrid: int, 
        grid: NDArray[np.float64], 
        mode: str = 'i'
    ) -> None:
        """
        Sets a common grid for all RadialSet objects
        """
        pass
    
    def __call__(self, itype: int, l: int, izeta: int) -> 'NumericalRadial': ...
    def symbol(self, itype: int) -> str: ...
    @property
    def ntype(self) -> int: ...
    @overload
    def lmax(self, itype: int) -> int: ...
    @property
    def lmax(self) -> int: ...
    @overload
    def rcut_max(self, itype: int) -> float: ...
    @property
    def rcut_max(self) -> float: ...
    def nzeta(self, itype: int, l: int) -> int: ...
    @overload
    def nzeta_max(self, itype: int) -> int: ...
    @overload
    def nzeta_max(self) -> int: ...
    @overload
    def nchi(self, itype: int) -> int: ...
    @overload
    def nchi(self) -> int: ...
    
class TwoCenterIntegrator: ...
    
class NumericalRadial: ...

    
    
        
    
        
    