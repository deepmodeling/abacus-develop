"""Offline regression for the real, two-frame latest-format ABACUS MD output."""

import unittest
import shutil
import tempfile
from pathlib import Path

from abacuslite.io.latestio import read_abacus_out
from common import assert_native_md_facts


class NativeMDOutputTest(unittest.TestCase):
    def setUp(self):
        # Engine-owned copy of a two-frame raw output fixture.
        source = Path(__file__).resolve().parent / 'fixtures' / 'native_md'
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.fixture = Path(temporary.name)
        for name, target in (
            ('running_md.log.txt', 'running_md.log'),
            ('eig_occ.txt', 'eig_occ.txt'),
            ('MD_dump', 'MD_dump'),
        ):
            shutil.copyfile(source / name, self.fixture / target)

    def test_real_dump_cell_and_positions(self):
        frames = read_abacus_out(self.fixture / 'running_md.log')
        assert_native_md_facts(self, frames, self.fixture / 'MD_dump', expected_frames=2)

    def test_raw_dump_oracle_rejects_mutated_cell_pbc_and_position(self):
        for mutation in ('cell', 'pbc', 'position'):
            with self.subTest(mutation=mutation):
                frames = read_abacus_out(self.fixture / 'running_md.log')
                if mutation == 'cell':
                    frames[0].cell[0, 0] += 0.1
                elif mutation == 'pbc':
                    frames[0].pbc[0] = False
                else:
                    frames[0].positions[0, 0] += 0.1
                with self.assertRaises(AssertionError):
                    assert_native_md_facts(self, frames, self.fixture / 'MD_dump',
                                           expected_frames=2)


if __name__ == '__main__':
    unittest.main()
