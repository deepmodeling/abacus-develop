import tempfile
import unittest
from pathlib import Path

from abacuslite.io.latestio import find_eig_occ_files


class TestFindEigOccFiles(unittest.TestCase):
    """Discovery rules for eigenvalue/occupation files."""

    def test_missing_announced_file_does_not_fall_back_to_stale_file(self):
        # Mode 2 announced eig_occg1.txt, but the file is gone. The stale
        # eig_occ.txt of an earlier mode 0/1 run must not be reported as the
        # current results.
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / 'eig_occ.txt').write_text('# Ionic_Step 1\n', encoding='utf-8')
            log = ['Write Eigenvalues and occupations to file: OUT.x/eig_occg1.txt']
            self.assertEqual(find_eig_occ_files(out, log), [])

    def test_announced_files_keep_step_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / 'eig_occg2.txt').write_text('', encoding='utf-8')
            (out / 'eig_occg1.txt').write_text('', encoding='utf-8')
            log = [
                'Write Eigenvalues and occupations to file: OUT.x/eig_occg2.txt',
                'Write Eigenvalues and occupations to file: OUT.x/eig_occg1.txt',
            ]
            found = find_eig_occ_files(out, log)
            self.assertEqual([path.name for path in found],
                             ['eig_occg1.txt', 'eig_occg2.txt'])

    def test_no_announcement_keeps_filename_discovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / 'eig_occ.txt').write_text('', encoding='utf-8')
            found = find_eig_occ_files(out, ['unrelated log line'])
            self.assertEqual([path.name for path in found], ['eig_occ.txt'])


if __name__ == '__main__':
    unittest.main()
