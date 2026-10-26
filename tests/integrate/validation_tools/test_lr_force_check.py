#!/usr/bin/env python3
"""Signed force regression must catch changes hidden by an absolute-value sum."""
import unittest
from pathlib import Path
import subprocess
import tempfile
from lr_force_check import compare, read_forces, read_reference

LOG = '''TOTAL-FORCE (eV/Angstrom)
H1 99 99 99
Start to calculate excited-state force of singlet
Forces (-gradients) of each excited state: (eV/Angstrom)
 state atom x y z
 0 0 force 1 -2 3
   1 force -4 5 -6
 1 0 force 7 -8 9
   1 force -10 11 -12
DONE
Start to calculate excited-state force of triplet
Forces (-gradients) of each excited state: (eV/Angstrom)
 state atom x y z
 0 0 force 13 -14 15
   1 force -16 17 -18
DONE
EXCITED-STATE TOTAL-FORCE (eV/Angstrom)
H1 88 88 88
'''
REFERENCE = '''singlet 0 0 1 -2 3
singlet 0 1 -4 5 -6
singlet 1 0 7 -8 9
singlet 1 1 -10 11 -12
triplet 0 0 13 -14 15
triplet 0 1 -16 17 -18
'''


class LRForceCheckTest(unittest.TestCase):
    def test_signed_components_match_and_exclude_gs_and_total(self):
        actual = read_forces(LOG)
        self.assertEqual(len(actual), 6)
        compare(actual, read_reference(REFERENCE), 1e-6)

    def test_sign_reversal_fails_with_identical_absolute_sum(self):
        actual = read_forces(LOG)
        flipped = {key: tuple(-v for v in row) for key, row in actual.items()}
        self.assertEqual(sum(abs(v) for row in actual.values() for v in row),
                         sum(abs(v) for row in flipped.values() for v in row))
        with self.assertRaises(ValueError):
            compare(flipped, actual, 1e-6)

    def test_atom_component_and_state_permutations_fail(self):
        actual = read_forces(LOG)
        variants = []
        swapped = dict(actual)
        first = ('singlet', 0, 0)
        second = ('singlet', 0, 1)
        swapped[first], swapped[second] = swapped[second], swapped[first]
        variants.append(swapped)
        variants.append({key: (row[1], row[0], row[2]) for key, row in actual.items()})
        variants.append({(spin, 1 - state, atom): row for (spin, state, atom), row in actual.items()})
        for variant in variants:
            with self.subTest(variant=variant), self.assertRaises(ValueError):
                compare(variant, actual, 1e-6)

    def test_missing_atom_channel_or_state_fails(self):
        expected = read_forces(LOG)
        for partial in [dict(list(expected.items())[:-1]),
                        {key: row for key, row in expected.items() if key[0] != 'triplet'},
                        {key: row for key, row in expected.items() if key[1] != 1}]:
            with self.subTest(partial=partial), self.assertRaises(ValueError):
                compare(partial, expected, 1e-6)

    def test_duplicates_and_malformed_rows_fail(self):
        for log in [LOG.replace('   1 force -4 5 -6', ' 0 0 force -4 5 -6'),
                    LOG.replace('1 -2 3', '1 -2'), '', LOG.replace('1 -2 3', 'nan -2 3')]:
            with self.subTest(log=log), self.assertRaises(ValueError):
                read_forces(log)
        with self.assertRaises(ValueError):
            read_reference(REFERENCE + REFERENCE)
        with self.assertRaises(ValueError):
            read_reference(REFERENCE.replace('1 -2 3', 'inf -2 3'))

    def test_shell_collector_emits_failure_for_sign_change_or_missing_data(self):
        tools = Path(__file__).resolve().parent
        script = '\n'.join([
            'source "$1/props_common.sh"', 'source "$1/props_lr.sh"',
            'PROPS_SCRIPT_DIR="$1"', 'is_lr=1', 'has_force=1',
            'running_path=running.log', 'props_result_file=result.out', 'run_lr_force_props'])
        variants = [(LOG, '0'), (LOG.replace('1 -2 3', '-1 2 -3'), '1'),
                    (LOG.replace('   1 force -16 17 -18\n', ''), '1')]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'lr_force.ref').write_text(REFERENCE)
            (root / 'threshold').write_text('force_threshold 0.000001\n')
            for log, expected in variants:
                (root / 'running.log').write_text(log)
                result = subprocess.run(['bash', '-e', '-c', script, 'collector-test', str(tools)],
                                        cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                values = dict(line.split() for line in (root / 'result.out').read_text().splitlines())
                self.assertEqual(values['CompareLRForce_pass'], expected)

    def test_tolerance_is_absolute_per_component(self):
        expected = read_forces(LOG)
        changed = dict(expected)
        changed[('singlet', 0, 0)] = (1 + 5e-7, -2, 3)
        compare(changed, expected, 1e-6)
        with self.assertRaises(ValueError):
            compare(changed, expected, 1e-7)
        for tolerance in [-1, float('nan'), float('inf')]:
            with self.assertRaises(ValueError):
                compare(changed, expected, tolerance)


if __name__ == '__main__':
    unittest.main()
