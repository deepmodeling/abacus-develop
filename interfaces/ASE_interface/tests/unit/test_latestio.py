# Tests extracted without semantic changes from the attributed upstream snapshot.
import unittest
from pathlib import Path
import numpy as np
from abacuslite.io.legacyio import read_traj_from_md_dump

from abacuslite.io.latestio import read_band_from_eig_occ, read_abacus_out

class TestLatestIO(unittest.TestCase):
    '''
    test the validities of functions for the I/O of the latest version of ABACUS
    '''
    here = Path(__file__).resolve().parents[2] / 'abacuslite/io'
    testfiles = here / 'testfiles'

    def test_read_abacus_out_reorders_calculator_magmoms(self):
        import tempfile
        from unittest.mock import patch

        frame = {
            'elem': ['Na', 'Na', 'Cl'],
            'coords': np.array([[0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [2.0, 0.0, 0.0]]),
            'cell': np.eye(3),
        }
        elecstate = [{
            'k': np.zeros((1, 1, 3)),
            'e': np.zeros((1, 1, 1)),
            'occ': np.ones((1, 1, 1)),
        }]
        energies = [{'E_KohnSham': -1.0, 'E_Fermi': 0.0}]
        kpoints = ((np.zeros((1, 3)), np.ones(1), None), None, None, None, None)

        with tempfile.TemporaryDirectory() as tmpdir:
            running_log = Path(tmpdir) / 'running_scf.log'
            running_log.write_text('')
            (Path(tmpdir) / 'eig_occ.txt').write_text('')
            with patch('abacuslite.io.latestio' + '.read_esolver_type_from_running_log', return_value='ksdft'), \
                 patch('abacuslite.io.latestio' + '.read_traj_from_running_log', return_value=[frame]), \
                 patch('abacuslite.io.latestio' + '.read_band_from_eig_occ', return_value=elecstate), \
                 patch('abacuslite.io.latestio' + '.read_forces_from_running_log', return_value=[]), \
                 patch('abacuslite.io.latestio' + '.read_stress_from_running_log', return_value=[]), \
                 patch('abacuslite.io.latestio' + '.read_kpoints_from_running_log', return_value=kpoints), \
                 patch('abacuslite.io.latestio' + '.read_energies_from_running_log', return_value=([], [])), \
                 patch('abacuslite.io.latestio' + '.read_iter_header_from_running_log', return_value=[]), \
                 patch('abacuslite.io.latestio' + '.find_final_info_with_iter_header', return_value=energies), \
                 patch('abacuslite.io.latestio' + '.read_magmom_from_running_log', return_value=[np.array([10.0, 20.0, 30.0])]):
                atoms = read_abacus_out(running_log, sort_atoms_with=[0, 2, 1])[0]

        self.assertEqual(atoms.get_chemical_symbols(), ['Na', 'Cl', 'Na'])
        self.assertTrue(np.allclose(atoms.cell.array, np.eye(3)))
        self.assertFalse(np.any(atoms.pbc))
        self.assertTrue(np.allclose(atoms.calc.results['magmoms'], [10.0, 30.0, 20.0]))

    def test_read_abacus_out_preserves_free_energy_and_exposes_internal_energy(self):
        """Keep force-consistent energy and expose the separate no-entropy value."""
        fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'multiframe_md_latest'
        frames = read_abacus_out(fixture / 'running_md.log')
        last = frames[-1]

        self.assertAlmostEqual(last.get_potential_energy(), -1940.6311063727, delta=1e-9)
        self.assertAlmostEqual(
            last.get_potential_energy(force_consistent=True),
            -1940.6311063727,
            delta=1e-9,
        )
        self.assertAlmostEqual(last.info["energy_without_entropy"], -1940.4234898941, delta=1e-9)
        self.assertAlmostEqual(frames[0].info["energy_without_entropy"], -1940.4285121056, delta=1e-9)

    def test_native_md_dump_cell_units_at_atoms_boundary(self):
        """MD_dump vectors are in alat; positions and energies are already physical."""
        fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'multiframe_md_latest'
        dump = fixture / 'MD_dump'
        log = fixture / 'running_md.log'
        before = {path.name: path.read_bytes() for path in fixture.iterdir()}

        raw = read_traj_from_md_dump(dump)
        self.assertEqual(len(raw), 2)
        self.assertEqual(raw[0]['alat_in_angstrom'], 0.529177)
        self.assertEqual(raw[0]['coordinate'], 'Cartesian')
        self.assertTrue(np.allclose(raw[0]['cell'][0], [5.34197, 5.34197, 0.0], atol=1e-12))
        self.assertTrue(np.allclose(raw[0]['coords'][0],
                                    [1.413423737926] * 3, atol=1e-12))

        frames = read_abacus_out(log)
        self.assertEqual(len(frames), 2)
        expected_cell = np.array([[2.82684765869, 2.82684765869, 0.0],
                                  [0.0, 2.82684765869, 2.82684765869],
                                  [2.82684765869, 0.0, 2.82684765869]])
        expected_positions = [
            [[1.413423737926] * 3, [0.0, 0.0, 0.0]],
            [[1.413423737926] * 3,
             [0.025620321279, 2.822443204757, 2.834296002714]],
        ]
        expected_energies = [-1940.6273962477, -1940.6311063727]
        for atoms, positions, energy in zip(frames, expected_positions, expected_energies):
            self.assertTrue(np.all(atoms.pbc))
            self.assertTrue(np.allclose(atoms.cell.array, expected_cell, atol=1e-10))
            self.assertAlmostEqual(atoms.get_volume(), 45.1791, delta=1e-4)
            self.assertTrue(np.allclose(atoms.positions, positions, atol=1e-12))
            self.assertAlmostEqual(atoms.get_potential_energy(), energy, delta=1e-9)
        self.assertEqual(before, {path.name: path.read_bytes() for path in fixture.iterdir()})

    def test_read_band_from_eig_occ_accepts_negative_mp_occupations(self):
        """MP smearing can produce small negative occupations in latest eig_occ."""
        lines = [
            '# ionic step',
            ' Electronic state energy (eV) and occupations',
            ' Spin number 2',
            ' spin=1 k-point=1/1 Cartesian=0.0000000 0.0000000 0.0000000 (250081 plane wave)',
            '1 -1.0 0.5',
            '2 -0.5 -0.0009381451611561922',
        ]
        elecstate = read_band_from_eig_occ(lines)
        self.assertEqual(len(elecstate), 1)
        self.assertEqual(elecstate[0]['e'].shape, (1, 1, 2))
        self.assertEqual(elecstate[0]['occ'].shape, (1, 1, 2))
        self.assertTrue(np.isclose(elecstate[0]['occ'][0, 0, 1], -0.0009381451611561922))
