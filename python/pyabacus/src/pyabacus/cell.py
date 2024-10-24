"""
Module for handling atomic structure and calculation parameters in ABACUS.
"""

from __future__ import annotations
import numpy as np
from typing import List, Union, Optional, Dict, Any, Tuple
import os

from pyabacus.io import read_stru, write_stru, read_xyz, write_xyz
from pyabacus.io.utils import cart_to_direct, direct_to_cart

class Cell:
    """
    A class representing an atomic cell structure for ABACUS calculations.
    
    This class handles atomic structures, lattice parameters, and calculation settings.
    It supports both STRU and XYZ file formats and provides a comprehensive interface
    for managing atomic systems in ABACUS calculations.
    
    Attributes:
        atoms: List of atoms with their properties
        lattice: Lattice vectors (3x3 numpy array)
        unit: Unit system ('Angstrom' or 'Bohr')
        spin: Spin configuration
        charge: Total charge
        lattice_constant: Scaling factor for lattice vectors
        basis_type: Type of basis set (e.g., 'lcao', 'pw')
        orbitals: List of orbital files
        pseudo_potentials: Dict of pseudopotential settings per element
        pseudo_dir: Directory containing pseudopotential files
        orbital_dir: Directory containing orbital files
        precision: Numerical precision for calculations
        ecutwfc: Energy cutoff for wavefunctions in Ry
    """

    def __init__(self):
        """Initialize an empty Cell object with default values."""
        # Structure attributes
        self._atoms = []  # List of [symbol, position, properties]
        self._lattice = np.eye(3)  # Default to unit cell
        self._lattice_constant = 1.0
        self._coord_type = 'Cartesian'

        # Physical properties
        self._unit = 'Angstrom'
        self._spin = 0
        self._charge = 0

        # Calculation settings
        self._basis_type = ''
        self._ecutwfc = 100.0  # Default cutoff energy in Ry
        self._precision = 1e-8
        self._mesh = None
        self._kspace = None

        # File paths and settings
        self.pseudo_dir = ''
        self.orbital_dir = ''
        self.pseudo_potentials = {}
        self.orbitals = []

        self._built = False
        
    def build(self) -> None:
        """
        Build the cell structure and set automatic parameters.
        
        This method ensures all necessary components are properly initialized
        and sets automatic parameters based on the current configuration.

        Raises:
            ValueError: If essential components are missing or invalid
        """
        if not self._atoms:
            raise ValueError("No atoms defined in the cell")
            
        self._set_auto_parameters()
        self._built = True

    def _set_auto_parameters(self) -> None:
        """
        Set automatic calculation parameters based on structure and precision.
        
        This method determines appropriate mesh settings based on the lattice
        vectors and precision requirements. It should be called whenever the
        lattice vectors or precision settings change significantly.
        """
        if self._lattice is not None:
            # Calculate mesh size based on lattice vectors and precision
            scaled_lattice = self._lattice * self._lattice_constant
            # Get lengths of lattice vectors
            lengths = np.sqrt(np.sum(scaled_lattice**2, axis=1))
            # Set mesh points inversely proportional to lattice vector lengths
            # and scaled by precision
            self._mesh = [max(1, int(np.ceil(length / self._precision))) 
                         for length in lengths]
        else:
            self._mesh = [10, 10, 10]  # Default mesh if no lattice is defined

    @classmethod
    def from_file(cls, file_path: str) -> 'Cell':
        """
        Create a Cell object from a structure file.

        Args:
            file_path: Path to either a STRU or XYZ file

        Returns:
            A new Cell object initialized from the file

        Raises:
            ValueError: If the file format is not supported
            FileNotFoundError: If the file doesn't exist
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        cell = cls()
        
        if file_path.lower().endswith('.xyz'):
            cell._load_xyz(file_path)
        elif os.path.basename(file_path).upper() == 'STRU':
            cell._load_stru(file_path)
        else:
            raise ValueError("Unsupported file format. Use .xyz or STRU files.")
        
        # Build the cell after loading the structure
        cell.build()
        return cell

    def to_file(self, file_path: str, file_type: Optional[str] = None) -> None:
        """
        Write the cell structure to a file.

        Args:
            file_path: Path where the file should be written
            file_type: Type of file to write ('xyz' or 'stru'). If None, inferred from file_path.

        Raises:
            ValueError: If the file type is not supported or cannot be determined
        """
        if file_type is None:
            if file_path.lower().endswith('.xyz'):
                file_type = 'xyz'
            elif os.path.basename(file_path).upper() == 'STRU':
                file_type = 'stru'
            else:
                raise ValueError("Cannot determine file type from file name")

        if file_type.lower() == 'xyz':
            species = [atom[0] for atom in self._atoms]
            positions = [atom[1] for atom in self._atoms]
            write_xyz(file_path, species, positions)
        elif file_type.lower() == 'stru':
            self._save_stru(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

    def _load_stru(self, stru_file: str) -> None:
        """Load structure from a STRU file."""
        stru_dict = read_stru(stru_file)
        
        # Set lattice information
        self.lattice_constant = stru_dict['lat']['const']
        if 'vec' in stru_dict['lat']:
            self.lattice = np.array(stru_dict['lat']['vec'])
        
        # Set coordinate type
        self._coord_type = stru_dict['coord_type']
        
        # Process atomic species and positions
        self._atoms = []
        self.orbitals = []

        for species in stru_dict['species']:
            # Clean the pseudopotential filename before storing
            pp_file = species['pp_file']
            if pp_file.startswith('./'):
                pp_file = pp_file[2:]  # Remove './' prefix if present

            # Store pseudopotential information
            self.pseudo_potentials[species['symbol']] = {
                'mass': species['mass'],
                'pseudo_file': pp_file
            }
            if 'pp_type' in species:
                self.pseudo_potentials[species['symbol']]['pp_type'] = species['pp_type']
            
            # Store orbital information if present
            if 'orb_file' in species:
                self.orbitals.append(species['orb_file'])
            
            # Process atoms
            for atom in species['atom']:
                pos = np.array(atom['coord'])
                # Convert position to Cartesian if it's in Direct coordinates
                if self._coord_type.lower() == 'direct':
                    pos = direct_to_cart(pos, self.lattice * self.lattice_constant)
                
                # Create atom entry with all properties
                properties = {k: v for k, v in atom.items() if k != 'coord'}
                self._atoms.append([species['symbol'], pos, properties])

    def _save_stru(self, file_path: str) -> None:
        """Save structure to a STRU file."""
        # Build STRU dictionary
        species_dict = {}
        for symbol, pos, props in self._atoms:
            if symbol not in species_dict:
                pp_info = self.pseudo_potentials.get(symbol, {})
            
            # Ensure we have a clean filename (remove './' if present)
            pp_file = pp_info.get('pseudo_file', f"{symbol}.UPF")
            if pp_file.startswith('./'):
                pp_file = pp_file[2:]
            
            species_dict[symbol] = {
                'symbol': symbol,
                'mass': pp_info.get('mass', 1.0),
                'pp_file': pp_file,  # Store clean filename
                'natom': 1,
                'mag_each': props.get('mag', 0.0),
                'atom': []
            }
            if 'pp_type' in pp_info:
                species_dict[symbol]['pp_type'] = pp_info['pp_type']
            else:
                species_dict[symbol]['natom'] += 1
            
            # Convert coordinates if needed
            coord = pos
            if self._coord_type.lower() == 'direct':
                coord = cart_to_direct(pos, self.lattice * self.lattice_constant)
            
            # Add atom with its properties
            atom_entry = {'coord': coord}
            atom_entry.update(props)
            species_dict[symbol]['atom'].append(atom_entry)

        # Add orbital information 
        if self.orbitals:
            for symbol, orb in zip(species_dict.keys(), self.orbitals):
                species_dict[symbol]['orb_file'] = orb
        
        
        stru_dict = {
            'lat': {
                'const': self.lattice_constant,
                'vec': self.lattice.tolist()
            },
            'coord_type': self._coord_type,
            'species': list(species_dict.values())
        }

        # Write to file
        write_stru(os.path.dirname(file_path), stru_dict, os.path.basename(file_path))

    def _load_xyz(self, xyz_file: str) -> None:
        """Load structure from an XYZ file."""
        species, coords = read_xyz(xyz_file)
        self._atoms = [[s, p, {}] for s, p in zip(species, coords)]

    def add_atom(self, 
                symbol: str, 
                position: Union[List[float], np.ndarray],
                properties: Optional[Dict[str, Any]] = None) -> None:
        """
        Add an atom to the cell.

        Args:
            symbol: Chemical symbol of the atom
            position: Position vector [x, y, z]
            properties: Additional properties (magnetic moment, constraints, etc.)
        """
        position = np.asarray(position)
        if properties is None:
            properties = {}
        self._atoms.append([symbol, position, properties])
        self._built = False  # Require rebuild after adding atoms

    @property
    def lattice(self) -> np.ndarray:
        """Get the lattice vectors."""
        return self._lattice

    @lattice.setter
    def lattice(self, value: Union[List[List[float]], np.ndarray]) -> None:
        """Set the lattice vectors."""
        self._lattice = np.asarray(value)
        if self._lattice.shape != (3, 3):
            raise ValueError("Lattice must be a 3x3 array")
        self._built = False  # Require rebuild after changing lattice

    @property
    def atoms(self) -> List[Tuple[str, np.ndarray, Dict[str, Any]]]:
        """Get the list of atoms with their properties."""
        return self._atoms.copy()

    @property
    def positions(self) -> np.ndarray:
        """Get atomic positions as a numpy array."""
        return np.array([atom[1] for atom in self._atoms])

    @property
    def species(self) -> List[str]:
        """Get list of atomic species symbols."""
        return [atom[0] for atom in self._atoms]

    @property
    def unit(self) -> str:
        """Get the unit system."""
        return self._unit

    @unit.setter
    def unit(self, value: str) -> None:
        """Set the unit system."""
        if value.lower() in ['angstrom', 'a']:
            self._unit = 'Angstrom'
        elif value.lower() in ['bohr', 'b', 'au']:
            self._unit = 'Bohr'
        else:
            raise ValueError("Unit must be 'Angstrom' or 'Bohr'")

    @property
    def ecutwfc(self) -> float:
        """Get the plane wave energy cutoff in Ry."""
        return self._ecutwfc

    @ecutwfc.setter
    def ecutwfc(self, value: float) -> None:
        """Set the plane wave energy cutoff in Ry."""
        if value <= 0:
            raise ValueError("Energy cutoff must be positive")
        self._ecutwfc = value
        
    @property
    def precision(self) -> float:
        """Get the numerical precision."""
        return self._precision

    @precision.setter
    def precision(self, value: float) -> None:
        """Set the numerical precision."""
        if value <= 0:
            raise ValueError("Precision must be a positive number")
        self._precision = value
        if self._built:
            self._set_auto_parameters()  # Update mesh when precision changes

    def get_scaled_positions(self) -> np.ndarray:
        """Get atomic positions in fractional coordinates."""
        if not self._built:
            raise RuntimeError("Cell has not been built. Call build() first.")
        return cart_to_direct(self.positions, self.lattice * self.lattice_constant)

    def make_kpts(self, mesh: List[int], with_gamma_point: bool = True) -> np.ndarray:
        """
        Generate k-points mesh in reciprocal space.

        Args:
            mesh: List of 3 integers specifying the k-point mesh
            with_gamma_point: Whether to include the gamma point

        Returns:
            Array of k-points in reciprocal space
        """
        if not self._built:
            raise RuntimeError("Cell has not been built. Call build() first.")
        if self.lattice is None:
            raise ValueError("Lattice vectors must be set before generating k-points.")

        kpts = []
        for i in range(mesh[0]):
            for j in range(mesh[1]):
                for k in range(mesh[2]):
                    if with_gamma_point:
                        kpt = np.array([i/mesh[0], j/mesh[1], k/mesh[2]])
                    else:
                        kpt = np.array([(i+0.5)/mesh[0], (j+0.5)/mesh[1], (k+0.5)/mesh[2]])
                    kpts.append(kpt)

        # Convert to cartesian coordinates in reciprocal space
        recip_lattice = 2 * np.pi * np.linalg.inv(self.lattice.T)
        kpts = np.dot(kpts, recip_lattice)

        return np.array(kpts)