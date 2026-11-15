# Tests extracted without semantic changes from the attributed upstream snapshot.
import tempfile
import unittest
from pathlib import Path
import numpy as np
from ase.atoms import Atoms

from abacuslite.io.generalio import (
    _parse_coord_line,
    load_pseudo,
    load_orbital,
    write_stru,
    read_stru,
)

class TestAbacusCalculatorIOUtil(unittest.TestCase):

    def setUp(self):
        here = Path(__file__).resolve().parents[2] / 'abacuslite/io'
        self.testfiles = here / 'testfiles'

        # will automatically delete the file after the context manager


    def test_stru_velocities_keep_physical_units_and_species_order(self):
        from ase.units import fs

        atoms = Atoms('HOH', positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                      cell=[8, 8, 8], pbc=True)
        # One ASE velocity unit is Angstrom per ASE time unit. These inputs
        # represent 1, 2 and -0.5 Angstrom/fs before species grouping.
        atoms.set_velocities(np.array([[1, 0, 0], [0, 2, 0], [0, 0, -0.5]]) / fs)
        # ABACUS md_func.cpp writes MD_dump velocity (Angstrom/fs) as STRU
        # atomic velocity divided by ANGSTROM_AU * AU_to_FS.
        au_per_angstrom_fs = 1.8897270 * 2.418884326505e-2
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_stru(atoms, tmpdir, {'H': 'H.upf', 'O': 'O.upf'})
            species = read_stru(path)['species']
        self.assertEqual([entry['symbol'] for entry in species], ['H', 'O'])
        self.assertEqual([entry['natom'] for entry in species], [2, 1])
        np.testing.assert_allclose([entry['v'] for entry in species[0]['atom']],
                                   np.array([[1, 0, 0], [0, 0, -0.5]]) * au_per_angstrom_fs,
                                   rtol=0, atol=1e-7)
        np.testing.assert_allclose(species[1]['atom'][0]['v'],
                                   [0, 2 * au_per_angstrom_fs, 0], rtol=0, atol=1e-7)

    def test_read_stru_accepts_signed_noncollinear_moment(self):
        atom = _parse_coord_line('0 0 0 mag 1.0 -1.0 0.0')
        self.assertEqual(atom['mag'], ('Cartesian', [1.0, -1.0, 0.0]))

    def test_write_stru_preserves_explicit_mass_for_each_species(self):
        atoms = Atoms('HOH', positions=np.zeros((3, 3)), cell=[8, 8, 8])
        atoms.set_masses([2.014, 15.999, 2.014])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_stru(atoms, tmpdir, {'H': 'H.upf', 'O': 'O.upf'})
            species = read_stru(path)['species']
        self.assertEqual([entry['symbol'] for entry in species], ['H', 'O'])
        self.assertEqual([entry['mass'] for entry in species], [2.014, 15.999])

    def test_write_stru_rejects_different_masses_within_species_before_writing(self):
        atoms = Atoms('HH', positions=np.zeros((2, 3)), cell=[8, 8, 8])
        atoms.set_masses([1.008, 2.014])
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, 'species H'):
                write_stru(atoms, tmpdir, {'H': 'H.upf'})
            self.assertFalse((Path(tmpdir) / 'STRU').exists())


        # kpt = {
        #     'mode': 'line',
        #     'nk': 6,
        #     'kpoints': [
        #         (0.0, 0.0, 0.0),
        #         (0.5, 0.5, 0.5),
        #         (0.0, 0.5, 0.5),
        #         (0.5, 0.5, 0.0),
        #         (0.5, 0.0, 0.5),
        #         (0.0, 0.0, 0.5),
        #     ],
        #     'nkinterpl': [10, 10, 10, 10, 10, 1],
        #     'coordinate': 'Cartesian',
        # }
        # with tempfile.NamedTemporaryFile(mode='w') as f:
        #     fn_ = Path(f.name)
        #     write_kpt(kpt, fn_.as_posix())
        #     self.assertDictEqual(kpt, read_kpt(fn_.as_posix()))

    def test_load_pseudo(self):
        expected = {
            'O': 'O.upf', 'Si': 'Si.upf', 'In': 'In.upf',
            'As': 'As.pbe-n-rrkjus_psl.0.2.UPF',
            'Ga': 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF',
        }
        # These helpers inspect filenames, not pseudopotential contents.
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in [*expected.values(), 'README.md', 'Si.orb']:
                (Path(tmpdir) / name).touch()
            self.assertDictEqual(load_pseudo(tmpdir), expected)

    def test_load_orbital(self):
        common = {'O': 'O_gga_6au_100Ry_2s2p1d.orb',
                  'Si': 'Si_gga_7au_100Ry_2s2p1d.orb'}
        small = 'In_gga_8au_100Ry_1s1p1d.orb'
        large = 'In_gga_8au_100Ry_2s2p2d1f.orb'
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in [*common.values(), large, small, 'README.md', 'In.upf']:
                (Path(tmpdir) / name).touch()
            self.assertDictEqual(load_orbital(tmpdir), {**common, 'In': small})
            self.assertDictEqual(load_orbital(tmpdir, efficiency=False), {**common, 'In': large})
