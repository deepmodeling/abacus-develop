# Tests extracted without semantic changes from the attributed upstream snapshot.
import shutil
import tempfile
import unittest
from pathlib import Path
import numpy as np

from abacuslite.io.legacyio import read_iter_header_from_running_log, read_abacus_out

class TestLegacyIO(unittest.TestCase):

    here = Path(__file__).resolve().parents[2] / 'abacuslite/io'
    testfiles = here / 'testfiles'

    def test_real_lts_native_md_uses_dump_frames_and_preserves_raw_files(self):
        source = Path(__file__).resolve().parents[1] / 'fixtures/lts_native_md'
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            original = {name: (source / name).read_bytes()
                        for name in ('running_md.log', 'MD_dump')}
            for name, content in original.items():
                (out / name).write_bytes(content)
            frames = read_abacus_out(out / 'running_md.log')
            self.assertEqual(len(frames), 3)
            self.assertEqual(read_iter_header_from_running_log(out / 'running_md.log')[-1],
                             (3, 1))
            expected_energies = [-210.6232157020, -210.6231743578, -210.6230563530]
            for frame, energy in zip(frames, expected_energies):
                self.assertTrue(np.all(frame.pbc))
                self.assertEqual(frame.get_chemical_symbols(), ['Si', 'Si'])
                self.assertAlmostEqual(frame.get_potential_energy(), energy, delta=1e-8)
            self.assertTrue(np.allclose(frames[0].cell.array[0],
                                        [0, 2.714998919, 2.714998919], atol=1e-8))
            self.assertTrue(np.allclose(frames[-1].positions[0],
                                        [2.715919335449, 2.715708402644, 5.429380230519],
                                        rtol=0, atol=1e-10))
            self.assertTrue(np.allclose(frames[-1].get_forces()[0],
                                        [-0.089751265334, -0.069019542505, 0.059917431359],
                                        rtol=0, atol=1e-8))
            self.assertEqual({name: (out / name).read_bytes() for name in original}, original)

    def test_lts_native_md_missing_or_truncated_dump_fails_clearly(self):
        source = Path(__file__).resolve().parents[1] / 'fixtures/lts_native_md'
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            shutil.copyfile(source / 'running_md.log', out / 'running_md.log')
            with self.assertRaisesRegex(ValueError, 'running_md.log.*1.*3'):
                read_abacus_out(out / 'running_md.log')
            dump = (source / 'MD_dump').read_text()
            (out / 'MD_dump').write_text(dump[:dump.index('MDSTEP:  2')])
            with self.assertRaisesRegex(ValueError, 'MD_dump.*frame'):
                read_abacus_out(out / 'running_md.log')
            last_step, tail = dump.split('MDSTEP:  2')
            lines = [line for line in tail.splitlines()
                     if not line.strip().startswith('1  Si')]
            (out / 'MD_dump').write_text(last_step + 'MDSTEP:  2\n' + '\n'.join(lines))
            with self.assertRaisesRegex(ValueError, 'Invalid or truncated MD_dump'):
                read_abacus_out(out / 'running_md.log')

    def test_complete_legacy_md_log_without_dump_retains_physical_geometry(self):
        source = Path(__file__).resolve().parents[1] / 'fixtures/multiframe_md_legacy/running_md.log'
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'running_md.log'
            original = source.read_bytes()
            path.write_bytes(original)
            frames = read_abacus_out(path)
            self.assertEqual(len(frames), 1)
            self.assertTrue(np.all(frames[0].pbc))
            self.assertAlmostEqual(frames[0].get_potential_energy(), -1940.5509086572)
            # The raw log has Direct (0.25, 0.25, 0.25), a0=0.529177 A,
            # and a cell with rows (5.34197, 5.34197, 0) cyclically.
            expected_cell = np.array([[5.34197, 5.34197, 0],
                                      [0, 5.34197, 5.34197],
                                      [5.34197, 0, 5.34197]]) * 0.529177
            np.testing.assert_allclose(frames[0].cell.array, expected_cell, rtol=0, atol=1e-12)
            np.testing.assert_allclose(frames[0].positions[0],
                                       np.array([0.25, 0.25, 0.25]) @ expected_cell,
                                       rtol=0, atol=1e-12)
            self.assertEqual(path.read_bytes(), original)

    def test_complete_cartesian_md_log_without_dump_converts_alat(self):
        source = Path(__file__).resolve().parents[1] / 'fixtures/multiframe_md_legacy/running_md.log'
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'running_md.log'
            original = source.read_text().replace('DIRECT COORDINATES\n', 'CARTESIAN COORDINATES\n', 1)
            original = original.replace('taud_As1            0.2500000000        0.2500000000        0.2500000000',
                                        'tauc_As1            2.0000000000        0.0000000000        0.0000000000', 1)
            original = original.replace('taud_Ga1', 'tauc_Ga1', 1)
            path.write_text(original)
            frames = read_abacus_out(path)
            self.assertEqual(len(frames), 1)
            np.testing.assert_allclose(frames[0].positions[0], [2 * 0.529177, 0, 0],
                                       rtol=0, atol=1e-12)
            self.assertTrue(np.all(frames[0].pbc))
            self.assertEqual(path.read_text(), original)


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

        with tempfile.NamedTemporaryFile(mode='w') as f:
            with patch('abacuslite.io.legacyio' + '.read_esolver_type_from_running_log', return_value='ksdft'), \
                 patch('abacuslite.io.legacyio' + '.read_traj_from_running_log', return_value=[frame]), \
                 patch('abacuslite.io.legacyio' + '.read_band_from_running_log', return_value=elecstate), \
                 patch('abacuslite.io.legacyio' + '.read_forces_from_running_log', return_value=[]), \
                 patch('abacuslite.io.legacyio' + '.read_stress_from_running_log', return_value=[]), \
                 patch('abacuslite.io.legacyio' + '.read_kpoints_from_running_log', return_value=kpoints), \
                 patch('abacuslite.io.legacyio' + '.read_energies_from_running_log', return_value=([], [])), \
                 patch('abacuslite.io.legacyio' + '.read_iter_header_from_running_log', return_value=[]), \
                 patch('abacuslite.io.legacyio' + '.find_final_info_with_iter_header', return_value=energies), \
                 patch('abacuslite.io.legacyio' + '.read_magmom_from_running_log', return_value=[np.array([10.0, 20.0, 30.0])]):
                atoms = read_abacus_out(f.name, sort_atoms_with=[0, 2, 1])[0]

        self.assertEqual(atoms.get_chemical_symbols(), ['Na', 'Cl', 'Na'])
        self.assertTrue(np.allclose(atoms.calc.results['magmoms'], [10.0, 30.0, 20.0]))


    def test_read_abacus_out_preserves_free_energy_and_exposes_internal_energy(self):
        """Keep force-consistent energy and expose the separate no-entropy value."""
        source = Path(__file__).resolve().parents[1] / 'fixtures/multiframe_md_legacy/running_md.log'
        atoms = read_abacus_out(source)[-1]

        self.assertAlmostEqual(atoms.get_potential_energy(), -1940.5509086572, delta=1e-9)
        self.assertAlmostEqual(
            atoms.get_potential_energy(force_consistent=True),
            -1940.5509086572,
            delta=1e-9,
        )
        self.assertAlmostEqual(atoms.info["energy_without_entropy"], -1940.3851061129, delta=1e-9)


    def test_read_iter_header_from_running_log(self):
        fn = self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax_'
        header = read_iter_header_from_running_log(fn)
        self.assertEqual(header[0], (1,1))
        self.assertEqual(header[1], (1,2))
        self.assertEqual(header[2], (2,1))
        fn = self.testfiles / 'pw-symm0-nspin4-gamma-md_'
        header = read_iter_header_from_running_log(fn)
        # LTS prints the first and last steps with padding, but ION=2 and
        # ION=3 can be unpadded; neither electronic step may be dropped.
        self.assertEqual(header, [(1, 1), (1, 2), (1, 3),
                                  (2, 1), (3, 1), (3, 2), (3, 3)])
