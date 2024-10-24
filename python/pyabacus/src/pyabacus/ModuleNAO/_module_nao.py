"""
pyabacus.ModuleNAO
==================
Module for Numerical Atomic Orbitals (NAO) in ABACUS
"""

import numpy as np
from numpy.typing import NDArray
from pyabacus.ModuleBase import SphericalBesselTransformer
from typing import overload, List, Tuple

from ._nao_pack import RadialCollection as _RadialCollection, TwoCenterIntegrator as _TwoCenterIntegrator, NumericalRadial as _NumericalRadial

class RadialCollection(_RadialCollection):
    def __init__(self) -> None: 
        """
        A class that holds all numerical radial functions of the same kind. 
        
        An instance of this class could be the collection of all radial functions 
        of numerical atomic orbitals, or all Kleinman-Bylander beta functions from 
        all elements involved in a calculation.
        """
        super().__init__()
    
    def build(
        self, 
        nfile: int, 
        file_list: List[str], 
        ftype: str = '\0'
    ) -> None:
        """
        Builds the collection from (orbital) files.
        """
        super().build(nfile, file_list, ftype)
    
    def set_transformer(
        self, 
        sbt: SphericalBesselTransformer, 
        update: int = 0
    ) -> None:
        """
        Sets a spherical Bessel transformers for all RadialSet objects.
        """
        super().set_transformer(sbt, update)
    
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
        super().set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft)
    
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
        super().set_grid(for_r_space, ngrid, grid, mode)
    
    def __call__(
        self, 
        itype: int, 
        l: int, 
        izeta: int
    ) -> 'NumericalRadial': 
        return super().__call__(itype, l, izeta)
        
    def symbol(self, itype: int) -> str: 
        return super().symbol(itype)
    
    @property
    def ntype(self) -> int: 
        return super().ntype
    
    def lmax(self, itype: int) -> int: 
        return super().lmax(itype)
    
    @property
    def lmax(self) -> int: 
        return super().lmax
    
    def rcut_max(self, itype: int) -> float: 
        return super().rcut_max(itype)

    @property
    def rcut_max(self) -> float: 
        return super().rcut_max
    
    def nzeta(self, itype: int, l: int) -> int: 
        return super().nzeta(itype, l)
    
    @overload
    def nzeta_max(self, itype: int) -> int: ...
    @overload
    def nzeta_max(self) -> int: ...
    
    def nzeta_max(self, *args, **kwargs): 
        return super().nzeta_max(*args, **kwargs)
    
    @overload
    def nchi(self, itype: int) -> int: ...
    @overload
    def nchi(self) -> int: ...
    
    def nchi(self, *args, **kwargs):
        return super().nchi(*args, **kwargs)
    
class TwoCenterIntegrator(_TwoCenterIntegrator):
    def __init__(self) -> None: 
        """
        A class to compute two-center integrals.

        This class computes two-center integrals of the form:

                            /    
                    I(R) = | dr phi1(r) (op) phi2(r - R)
                            /               

        as well as their gradients, where op is 1 (overlap) or minus Laplacian (kinetic), and phi1, 
        phi2 are "atomic-orbital-like" functions of the form:

                    phi(r) = chi(|r|) * Ylm(r/|r|)

        where chi is some numerical radial function and Ylm is some real spherical harmonics.

        This class is designed to efficiently compute the two-center integrals
        between two "collections" of the above functions with various R, e.g., the
        overlap integrals between all numerical atomic orbitals and all
        Kleinman-Bylander nonlocal projectors, the overlap & kinetic integrals between all numerical atomic orbitals, etc.
        This is done by tabulating the radial part of the integrals on an r-space grid and the real Gaunt coefficients in advance.
        """
        super().__init__()
    
    def tabulate(
        self, 
        bra: 'RadialCollection', 
        ket: 'RadialCollection', 
        op: str, 
        nr: int, 
        cutoff: float
    ) -> None:
        """
        Tabulates the radial part of a two-center integral.

        Parameters:
        bra (RadialFunctions): The radial functions of the first collection.
        ket (RadialFunctions): The radial functions of the second collection.
        op (char): Operator, could be 'S' (overlap) or 'T' (kinetic).
        nr (int): Number of r-space grid points.
        cutoff (float): r-space cutoff radius.
        """
        super().tabulate(bra, ket, op, nr, cutoff)
    
    def calculate(
        self, 
        itype1: int, 
        l1: int, 
        izeta1: int, 
        m1: int, 
        itype2: int, 
        l2: int, 
        izeta2: int, 
        m2: int, 
        pvR: NDArray[np.float64], 
        cal_grad: bool = False
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the two-center integrals.

        This function calculates the two-center integral

                            /    
                    I(R) = | dr phi1(r) (op_) phi2(r - R)
                            /               

        or its gradient by using the tabulated radial part and real Gaunt coefficients.

        Parameters
        ----------
        itype1 : int
            Element index of orbital 1.
        l1 : int
            Angular momentum of orbital 1.
        izeta1 : int
            Zeta number of orbital 1.
        m1 : int
            Magnetic quantum number of orbital 1.
        itype2 : int
            Element index of orbital 2.
        l2 : int
            Angular momentum of orbital 2.
        izeta2 : int
            Zeta number of orbital 2.
        m2 : int
            Magnetic quantum number of orbital 2.
        pvR : array_like
            R2 - R1, the displacement vector between the two centers.
        cal_grad : bool, optional
            The gradient will not be computed if cal_grad is false.
        
        Returns
        -------
        out_array : array_like
            The two-center integral.
        grad_out_array : array_like
            Gradient of the two-center integral.
        """
        return super().calculate(itype1, l1, izeta1, m1, itype2, l2, izeta2, m2, pvR, cal_grad)
    
    def snap(
        self, 
        itype1: int, 
        l1: int, 
        izeta1: int, 
        m1: int, 
        itype2: int, 
        pvR: NDArray[np.float64], 
        deriv: bool
    ) -> List[List[float]]:
        """
        Compute a batch of two-center integrals.

        This function calculates the two-center integrals (and optionally their gradients)
        between one orbital and all orbitals of a certain type from the other collection.
        """
        return super().snap(itype1, l1, izeta1, m1, itype2, pvR, deriv)

class NumericalRadial(_NumericalRadial):
    def __init__(self) -> None: 
        """
        A class that represents a numerical radial function.

        This class is designed to be the container for the radial part of numerical atomic orbitals, Kleinman-Bylander beta functions, and all other similar numerical radial functions in three-dimensional space, each of which is associated with some angular momentum l and whose r and k space values are related by an l-th order spherical Bessel transform.

        A NumericalRadial object can be initialized by "build", which requires the angular momentum, the number of grid points, the grid and the corresponding values. Grid does not have to be uniform. One can initialize the object in either r or k space. After initialization, one can set the
        grid in the other space via set_grid or set_uniform_grid. Values in the other space are automatically computed by a spherical Bessel transform.
        """
        super().__init__()
    
    def build(
        self, 
        l: int, 
        for_r_space: bool, 
        ngrid: int, 
        grid: NDArray[np.float64], 
        value: NDArray[np.float64], 
        p: int = 0, 
        izeta: int = 0, 
        symbol: str = "", 
        itype: int = 0, 
        init_sbt: bool = True
    ) -> None:
        """
        Initializes the object by providing the grid & values in one space.

        Parameters
        ----------
        l : int
            Angular momentum.
        for_r_space : bool
            Specifies whether the input corresponds to r or k space.
        ngrid : int
            Number of grid points.
        grid : array_like
            Grid points, must be positive & strictly increasing.
        value : array_like
            Values on the grid.
        p : float
            Implicit exponent in input values, see pr_ & pk_.
        izeta : int
            Multiplicity index of radial functions of the same itype and l.
        symbol : str
            Chemical symbol.
        itype : int
            Index for the element in calculation.
        init_sbt : bool
            If true, internal SphericalBesselTransformer will be initialized.

        Notes
        -----
        init_sbt is only useful when the internal SphericalBesselTransformer (sbt_) is null-initialized; The function will NOT reset sbt_ if it is already usable.
        """
        super().build(l, for_r_space, ngrid, grid, value, p, izeta, symbol, itype, init_sbt)
    
    def set_transformer(
        self, 
        sbt: SphericalBesselTransformer, 
        update: int = 0
    ) -> None:
        """
        Sets a SphericalBesselTransformer.

        By default, the class uses an internal SphericalBesselTransformer, but one can optionally use a shared one. This could be beneficial when there are a lot of NumericalRadial objects whose grids have the same size.

        Parameters
        ----------
        sbt : SphericalBesselTransformer
            An external transformer.
        update : int
            Specifies whether and how values are recomputed with the new transformer.
            Accepted values are:
            *  0: does not recompute values;
            *  1: calls a forward transform;
            * -1: calls a backward transform.
        """
        super().set_transformer(sbt, update)
    
    def set_grid(
        self, 
        for_r_space: bool, 
        ngrid: int, 
        grid: NDArray[np.float64], 
        mode: str = 'i'
    ) -> None:
        """
        Sets up a grid.

        This function can be used to set up the grid which is absent in "build" (in which case values on the new grid are automatically computed by a spherical Bessel transform) or interpolate on an existing grid to a new grid.

        Parameters
        ----------
        for_r_space : bool
            Specifies whether to set grid for the r or k space.
        ngrid : int
            Number of grid points.
        grid : array_like
            Grid points, must be positive & strictly increasing.
        mode : char
            Specifies how values are updated, could be 'i' or 't':
            - 'i': New values are obtained by interpolating and zero-padding
                the existing values from current space. With this option,
                it is an error if the designated space does not have a grid;
            - 't': New values are obtained via transform from the other space.
                With this option, it is an error if the other space does not
                have a grid.
        """
        super().set_grid(for_r_space, ngrid, grid, mode)
    
    def set_uniform_grid(
        self, 
        for_r_space: bool, 
        ngrid: int, 
        cutoff: float, 
        mode: str = 'i', 
        enable_fft: bool = False
    ) -> None:
        """
        Sets up a uniform grid.

        The functionality of this function is similar to set_grid, except that the new grid is a uniform grid specified by the cutoff and the number of grid points, which are calculated as:

            grid[i] = i * (cutoff / (ngrid - 1))

        Parameters
        ----------
        for_r_space : bool
            Specifies whether to set grid for the r or k space.
        ngrid : int
            Number of grid points.
        cutoff : float
            The maximum value of the grid, which determines the grid spacing.
        enable_fft : bool
            If true, this function will not only set up the grid & values in the designated space, but also sets the grid in the other space such that the r & k grids are FFT-compliant (and updates values via a FFT-based spherical Bessel transform).
        mode : char
            Specifies how values are updated, could be 'i' or 't'.
        """
        super().set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft)
    
    def set_value(
        self, 
        for_r_space: bool, 
        value: NDArray[np.float64], 
        p: int
    ) -> None:
        """
        Updates values on an existing grid.

        This function does not alter the grid; it merely updates values on the existing grid. The number of values to read from "value" is determined by the current number of points in the r or k space (nr_ or nk_). Values of the other space will be automatically updated if they exist.

        Warning
        -------
        This function does not check the index bound; use with care!
        """
        super().set_value(for_r_space, value, p)
    
    def wipe(
        self, 
        r_space: bool = True, 
        k_space: bool = True
    ) -> None: 
        super().wipe(r_space, k_space)
        
    def normalize(self, for_r_space: bool = True) -> None: 
        """
        Normalizes the radial function.

        The radial function is normalized such that the integral of the square of the function multiplied by the square of the radial coordinate over the entire space is equal to one:

            ∫ from 0 to +∞ of (x^2 * f(x)^2) dx = 1

        where x is r or k. The integral is evaluated with Simpson's rule. Values in the other space are updated automatically via a spherical Bessel transform.
        """
        super().normalize(for_r_space)
    
    @property
    def symbol(self) -> str: 
        return super().symbol
    @property
    def itype(self) -> int: 
        return super().itype
    @property
    def izeta(self) -> int: 
        return super().izeta
    @property
    def l(self) -> int: 
        return super().l
    @property
    def nr(self) -> int: 
        return super().nr
    @property
    def nk(self) -> int: 
        return super().nk
    @property
    def rcut(self) -> float: 
        return super().rcut
    @property
    def kcut(self) -> float: 
        return super().kcut
    @property
    def rmax(self) -> float: 
        return super().rmax
    @property
    def kmax(self) -> float: 
        return super().kmax
    @property
    def pr(self) -> float: 
        return super().pr
    @property
    def pk(self) -> float: 
        return super().pk
    @property
    def sbt(self) -> SphericalBesselTransformer: 
        return super().sbt
    @property
    def rgrid(self) -> NDArray[np.float64]: 
        return super().rgrid
    @property
    def kgrid(self) -> NDArray[np.float64]: 
        return super().kgrid
    @property
    def rvalue(self) -> NDArray[np.float64]: 
        return super().rvalue
    @property
    def kvalue(self) -> NDArray[np.float64]: 
        return super().kvalue
    @property
    def is_fft_compliant(self) -> bool: 
        return super().is_fft_compliant


    
    
        
    
        
    