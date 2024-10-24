import unittest
import numpy as np
import os
import tempfile
import shutil
from pyabacus import Cell

class TestCell(unittest.TestCase):
    """Test suite for the Cell class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Get the absolute path to the test directory
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_cell_dir = os.path.join(self.test_dir, 'test_cell')
        
        # Set paths for test files
        self.lcao_dir = os.path.join(self.test_cell_dir, 'lcao_ZnO')
        self.stru_file = os.path.join(self.lcao_dir, 'STRU')
        self.xyz_file = os.path.join(self.test_cell_dir, 'h2o.xyz')
        
        # Create a temporary directory for output files
        self.temp_dir = tempfile.mkdtemp()
        
        # Verify test files exist
        if not os.path.exists(self.stru_file):
            raise FileNotFoundError(f"STRU file not found at {self.stru_file}")
        if not os.path.exists(self.xyz_file):
            raise FileNotFoundError(f"XYZ file not found at {self.xyz_file}")

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_initialization(self):
        """Test cell initialization."""
        cell = Cell()
        self.assertFalse(cell._built)
        self.assertEqual(cell.unit, 'Angstrom')
        self.assertEqual(cell.ecutwfc, 100.0)
        self.assertEqual(cell.precision, 1e-8)
        np.testing.assert_array_equal(cell.lattice, np.eye(3))

    
    def test_from_stru_file(self):
        """Test loading structure from STRU file."""
        cell = Cell.from_file(self.stru_file)
        self.assertTrue(cell._built)
        
        # Verify basic structure information
        self.assertEqual(len(cell.species), 2)
        self.assertEqual(cell.species, ['Zn', 'O'])
        
        # Check lattice parameters and vectors from the actual STRU file
        self.assertEqual(cell.lattice_constant, 6.1416)  # From STRU file
        
        scaled_lattice = cell.lattice * cell.lattice_constant
        expected_lattice = np.array([
            [6.1416, 0.0, 0.0],
            [-3.0708, 5.3186256, 0.0],
            [0.0, 0.0, 9.82656]
        ])
        np.testing.assert_array_almost_equal(scaled_lattice, expected_lattice)

    def test_from_xyz_file(self):
        """Test loading structure from XYZ file."""
        cell = Cell.from_file(self.xyz_file)
        self.assertTrue(cell._built)
        
        # Check basic structure
        self.assertEqual(len(cell.species), 3)
        self.assertEqual(cell.species, ['O', 'H', 'H'])
        
        # Check positions match the h2o.xyz file content
        positions = cell.positions
        self.assertEqual(len(positions), 3)
        
    def test_build_status(self):
        """Test build status management."""
        cell = Cell()
        self.assertFalse(cell._built)
        
        # Adding atoms should keep cell unbuilt
        cell.add_atom('H', [0, 0, 0])
        self.assertFalse(cell._built)
        
        # Building should set the flag
        cell.build()
        self.assertTrue(cell._built)

    def test_auto_parameters(self):
        """Test automatic parameter setting."""
        cell = Cell()
        cell.add_atom('H', [0, 0, 0])
        cell.lattice = [[5.0, 0.0, 0.0],
                       [0.0, 5.0, 0.0],
                       [0.0, 0.0, 5.0]]
        cell.build()
        
        # Check mesh generation
        self.assertIsNotNone(cell._mesh)
        self.assertEqual(len(cell._mesh), 3)

    def test_stru_pseudopotentials(self):
        """Test pseudopotential information from STRU file."""
        cell = Cell.from_file(self.stru_file)
        
        # Check pseudopotential information
        self.assertIn('Zn', cell.pseudo_potentials)
        self.assertIn('O', cell.pseudo_potentials)
        
        # Verify specific pseudopotential files from STRU
        self.assertEqual(cell.pseudo_potentials['Zn']['pseudo_file'], 'Zn.LDA.UPF')
        self.assertEqual(cell.pseudo_potentials['O']['pseudo_file'], 'O.LDA.100.UPF')

    def test_stru_orbital_files(self):
        """Test orbital information from STRU file."""
        cell = Cell.from_file(self.stru_file)
        
        # Check orbital files
        self.assertEqual(len(cell.orbitals), 2)
        self.assertIn('Zn_lda_8.0au_120Ry_2s2p2d', cell.orbitals)
        self.assertIn('O_lda_7.0au_50Ry_2s2p1d', cell.orbitals)

    def test_stru_coordinates(self):
        """Test coordinate system and positions from STRU file."""
        cell = Cell.from_file(self.stru_file)
        
        # Check coordinate type
        self.assertEqual(cell._coord_type, 'Direct')
        
        # Get scaled positions
        scaled_positions = cell.get_scaled_positions()
        
        # Expected positions from STRU file
        expected_scaled = np.array([
            [0.00, 0.00, 0.00],  # Zn
            [0.33333, 0.66667, 0.50]  # O
        ])
        
        np.testing.assert_array_almost_equal(scaled_positions, expected_scaled, decimal=4)

    def test_file_operations(self):
        """Test file reading and writing operations."""
        # Load from STRU
        original_cell = Cell.from_file(self.stru_file)
        
        # Save to new STRU
        new_stru = os.path.join(self.temp_dir, 'STRU')
        original_cell.to_file(new_stru, 'stru')
        
        # Load back and compare
        new_cell = Cell.from_file(new_stru)
        
        # Compare structures
        np.testing.assert_array_almost_equal(original_cell.positions, new_cell.positions)
        np.testing.assert_array_almost_equal(original_cell.lattice, new_cell.lattice)
        self.assertEqual(original_cell.species, new_cell.species)
        self.assertEqual(original_cell.pseudo_potentials, new_cell.pseudo_potentials)
        self.assertEqual(original_cell.orbitals, new_cell.orbitals)

        # Test XYZ file operations
        xyz_cell = Cell.from_file(self.xyz_file)
        new_xyz = os.path.join(self.temp_dir, 'test.xyz')
        xyz_cell.to_file(new_xyz, 'xyz')
        new_xyz_cell = Cell.from_file(new_xyz)
        np.testing.assert_array_almost_equal(xyz_cell.positions, new_xyz_cell.positions)
        self.assertEqual(xyz_cell.species, new_xyz_cell.species)

    def test_add_atom(self):
        """Test atom addition functionality."""
        cell = Cell()
        
        # Add atom with properties
        properties = {"mag": 0.5, "constraint": [1, 1, 1]}
        cell.add_atom("Fe", [0, 0, 0], properties)
        
        self.assertEqual(len(cell.atoms), 1)
        self.assertEqual(cell.species[0], "Fe")
        atom_data = cell.atoms[0]
        self.assertEqual(atom_data[2], properties)
        
        # Check if adding atom unsets built status
        cell.build()
        self.assertTrue(cell._built)
        cell.add_atom("Fe", [1, 1, 1])
        self.assertFalse(cell._built)

    def test_k_points(self):
        """Test k-points generation."""
        cell = Cell()
        cell.add_atom('H', [0, 0, 0])
        cell.lattice = [[5.0, 0.0, 0.0],
                       [0.0, 5.0, 0.0],
                       [0.0, 0.0, 5.0]]
        cell.build()
        
        # Test k-points generation
        kpts = cell.make_kpts([2, 2, 2])
        self.assertEqual(kpts.shape, (8, 3))
        
        # Test k-points with non-gamma
        kpts_nogamma = cell.make_kpts([2, 2, 2], with_gamma_point=False)
        self.assertEqual(kpts_nogamma.shape, (8, 3))

if __name__ == '__main__':
    unittest.main()