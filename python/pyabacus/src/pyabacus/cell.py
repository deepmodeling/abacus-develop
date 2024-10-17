import numpy as np
import os

class Cell:
    def __init__(self):
        self.atom = None
        self.a = None  # Lattice vectors
        self.unit = 'Angstrom'  # Default unit
        self.spin = 0  # Default spin
        self.charge = 0  # Default charge
        self.lattice_constant = 6.1416 # 
        self.basis = None
        self.pseudo = None
        self.orbitals = []
        self.pseudo_potentials = {}
        self.pseudo_dir = ''
        self.orbital_dir = ''
        self.basis_type = ''
        self.built = False
        self._kspace = None
        self.precision = 1e-8  # Default precision
        self._mesh = None
        self.ke_cutoff = None
        self.rcut = None

    @classmethod
    def from_file(cls, stru_file):
        cell = cls()
        cell._parse_stru(stru_file)
        cell._built = True
        return cell

    def build(self):
        if self.atom is None:
            raise ValueError("Atom information must be set before building.")

        if isinstance(self.atom, str):
            if self.atom.endswith('.xyz'):
                self._parse_xyz(self.atom)
            else:
                raise ValueError("Unsupported file format. Use .xyz files or provide atom list directly.")
        elif isinstance(self.atom, list):
            self.atoms = [[atom[0], np.array(atom[1])] for atom in self.atom]
        else:
            raise ValueError("Unsupported atom format.")

        # Automatically set parameters based on precision
        self._set_auto_parameters()

        self._built = True

    def _set_auto_parameters(self):
        if self.a is not None:
            self.mesh = [int(np.ceil(np.linalg.norm(v) / self.precision)) for v in self.a] # TODO: Check the formula!
        else:
            self.mesh = [10, 10, 10]  # Default mesh if lattice vectors are not set

        self.ke_cutoff = -np.log(self.precision) * 10  # TODO: Check the formula!
        self.rcut = -np.log(self.precision) * 2  # TODO: Check the formula!

    def _parse_stru(self, stru_file):
        self.atoms = []
        with open(stru_file, 'r') as f:
            lines = f.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if 'ATOMIC_SPECIES' in line:
                    i += 1
                    while i < len(lines) and lines[i].strip():
                        parts = lines[i].split()
                        if len(parts) == 3:
                            species, mass, pp_file = parts
                            pp_file = pp_file.lstrip('./')
                            self.pseudo_potentials[species] = {
                            'mass': float(mass),
                            'pseudo_file': pp_file
                            }
                        i += 1
                elif 'NUMERICAL_ORBITAL' in line:
                    i += 1
                    while i < len(lines) and lines[i].strip():
                        orbital = lines[i].split()
                        self.orbitals.append(orbital)
                        i += 1
                elif 'LATTICE_CONSTANT' in line:
                    i += 1
                    self.lattice_constant = float(lines[i].strip())
                    i += 1
                elif 'LATTICE_VECTORS' in line:
                    self.a = np.array([
                        list(map(float, lines[i+1].split())),
                        list(map(float, lines[i+2].split())),
                        list(map(float, lines[i+3].split()))
                    ])
                    i += 4
                elif 'ATOMIC_POSITIONS' in line:
                    i += 3
                    while i < len(lines) and lines[i].strip():
                        species = lines[i].strip()
                        i += 2
                        num_atoms = int(lines[i].strip())
                        i += 1
                        for _ in range(num_atoms):
                            pos = np.array(list(map(float, lines[i].split()[:3])))
                            self.atoms.append([species, pos])
                            i += 2
                else:
                    i += 1

    def _parse_xyz(self, xyz_file):
        self.atoms = []
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            num_atoms = int(lines[0])
            # Skip the comment line
            for line in lines[2:2+num_atoms]:
                parts = line.split()
                species = parts[0]
                coords = np.array(list(map(float, parts[1:4])))
                self.atoms.append([species, coords])

    def get_atom_positions(self):
        if not self._built:
            raise RuntimeError("Cell has not been built. Call build() first.")
        return np.array([atom[1] for atom in self.atoms])

    def get_atom_species(self):
        if not self._built:
            raise RuntimeError("Cell has not been built. Call build() first.")
        return [atom[0] for atom in self.atoms]

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if value.lower() in ['angstrom', 'a']:
            self._unit = 'Angstrom'
        elif value.lower() in ['bohr', 'b', 'au']:
            self._unit = 'Bohr'
        else:
            raise ValueError("Unit must be 'Angstrom' or 'Bohr'")

    @property
    def lattice_constant(self):
        return self._lattice_constant
    
    @lattice_constant.setter
    def lattice_constant(self, value):
        self._lattice_constant = value

    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, value):
        if value <= 0:
            raise ValueError("Precision must be a positive number")
        self._precision = value

    @property
    def kspace(self):
        return self._kspace

    @kspace.setter
    def kspace(self, value):
        if value <= 0:
            raise ValueError("k-space must be a positive number")
        self._kspace = value

    def make_kpts(self, mesh, with_gamma_point=True):
        if self.a is None:
            raise ValueError("Lattice vectors (self.a) must be set before generating k-points.")

        kpts = []
        for i in range(mesh[0]):
            for j in range(mesh[1]):
                for k in range(mesh[2]):
                    if with_gamma_point:
                        kpt = np.array([i/mesh[0], j/mesh[1], k/mesh[2]])
                    else:
                        kpt = np.array([(i+0.5)/mesh[0], (j+0.5)/mesh[1], (k+0.5)/mesh[2]])
                    kpts.append(kpt)

        # Convert to cartesian coordinates
        recip_lattice = 2 * np.pi * np.linalg.inv(self.a.T)
        kpts = np.dot(kpts, recip_lattice)

        return np.array(kpts)