import unittest
import numpy as np
import tempfile
import os
from pyabacus import Cell

class TestCell(unittest.TestCase):
    def setUp(self):
        # Use the existing STRU file
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.stru_file = os.path.join(self.test_dir, 'test_cell', 'lcao_ZnO', 'STRU')
        
        # Path for the XYZ file (assuming it exists in the test_cell folder)
        self.xyz_file = os.path.join(self.test_dir, 'test_cell', 'h2o.xyz')

    def test_from_file(self):
        cell = Cell.from_file(self.stru_file)
        self.assertEqual(len(cell.atoms), 2)
        self.assertEqual(cell.get_atom_species(), ['Zn', 'O'])
        expected_lattice = np.array([
            [1.00, 0.00, 0.00],
            [-0.5, 0.866, 0.00],
            [0.00, 0.00, 1.6]
        ])
        np.testing.assert_array_almost_equal(cell.a, expected_lattice)

    def test_from_xyz_file(self):
        cell = Cell()
        cell.atom = self.xyz_file
        cell.build()
        self.assertEqual(len(cell.atoms), 3)
        self.assertEqual(cell.get_atom_species(), ['O', 'H', 'H'])

    def test_pseudo_potentials(self):
        cell = Cell.from_file(self.stru_file)
        self.assertIn('Zn', cell.pseudo_potentials)
        self.assertIn('O', cell.pseudo_potentials)
        self.assertEqual(cell.pseudo_potentials['Zn']['pseudo_file'], 'Zn.LDA.UPF')
        self.assertEqual(cell.pseudo_potentials['O']['pseudo_file'], 'O.LDA.100.UPF')


    def test_atomic_positions(self):
        cell = Cell.from_file(self.stru_file)
        expected_positions = np.array([
            [0.00, 0.00, 0.00],
            [0.33333, 0.66667, 0.50]
        ])
        np.testing.assert_array_almost_equal(cell.get_atom_positions(), expected_positions)
    
    def test_build(self):
        cell = Cell()
        cell.atom = [['H', [0, 0, 0]], ['O', [0, 0, 1]], ['H', [0, 1, 0]]]
        cell.a = np.eye(3) * 3.0
        cell.build()
        self.assertTrue(cell._built)
        self.assertIsNotNone(cell.mesh)
        self.assertIsNotNone(cell.ke_cutoff)
        self.assertIsNotNone(cell.rcut)

    def test_make_kpts(self):
        cell = Cell()
        cell.atom = [['H', [0, 0, 0]], ['O', [0, 0, 1]], ['H', [0, 1, 0]]]
        cell.a = np.eye(3) * 3.0
        cell.build()
        kpts = cell.make_kpts([2, 2, 2])
        self.assertEqual(kpts.shape, (8, 3))

    def test_precision(self):
        cell = Cell()
        cell.precision = 1e-10
        self.assertEqual(cell.precision, 1e-10)
        with self.assertRaises(ValueError):
            cell.precision = -1

    def test_unit(self):
        cell = Cell()
        cell.unit = 'Angstrom'
        self.assertEqual(cell.unit, 'Angstrom')
        cell.unit = 'Bohr'
        self.assertEqual(cell.unit, 'Bohr')
        with self.assertRaises(ValueError):
            cell.unit = 'invalid_unit'
    

if __name__ == '__main__':
    unittest.main()