#!/usr/bin/env python3
"""Signed force regression must catch changes hidden by an absolute-value sum."""
import unittest
import os
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
NEW_LOG = LOG.replace(
    'Forces (-gradients) of each excited state:',
    'Excitation-energy force contributions (-dOmega/dR) of each state:')

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

    def test_new_and_mixed_titles_match_the_same_reference(self):
        mixed_log = NEW_LOG.replace(
            'Excitation-energy force contributions (-dOmega/dR) of each state:',
            'Forces (-gradients) of each excited state:', 1)
        extra_tables = (NEW_LOG + 'GROUND-STATE FORCE (eV/Angstrom)\n'
                        '0 0 force 999 999 999\n'
                        'EXCITED-STATE TOTAL-FORCE (eV/Angstrom)\n'
                        '0 0 force 888 888 888\n')
        for log in (NEW_LOG, mixed_log, extra_tables):
            with self.subTest(log=log):
                compare(read_forces(log), read_reference(REFERENCE), 1e-6)

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
        extra_tables = (NEW_LOG + 'GROUND-STATE FORCE (eV/Angstrom)\n'
                        '0 0 force 999 999 999\n'
                        'EXCITED-STATE TOTAL-FORCE (eV/Angstrom)\n'
                        '0 0 force 888 888 888\n')
        variants = [(LOG, '0'), (NEW_LOG, '0'), (extra_tables, '0'),
                    (NEW_LOG.replace('1 -2 3', '-1 2 -3'), '1'),
                    (NEW_LOG.replace('   1 force -16 17 -18\n', ''), '1')]
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
                if expected == '0':
                    self.assertEqual(float(values['totalforceref']), 171.0)

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

    def test_density_collector_requires_exact_nonempty_channel_files(self):
        tools = Path(__file__).resolve().parent
        script = '\n'.join([
            'source "$1/props_lr.sh"', 'is_lr=1',
            'props_result_file=result.out', 'run_lr_density_props'])
        expected_names = ['LR_e-h_density_triplet_1.cube', 'LR_e-h_density_triplet_2.cube']
        variants = [
            (dict.fromkeys(expected_names, 'cube data'), '0'),
            ({name.replace('triplet', 'singlet'): 'cube data' for name in expected_names}, '1'),
            ({expected_names[0]: 'cube data'}, '1'),
            ({expected_names[0]: 'cube data', expected_names[1]: ''}, '1'),
            (dict.fromkeys(expected_names + ['LR_e-h_density_singlet_1.cube'], 'cube data'), '1'),
            ({}, '1')]
        for files, expected in variants:
            with self.subTest(files=files), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                output = root / 'OUT.autotest'
                output.mkdir()
                (root / 'lr_density.ref').write_text('\n'.join(expected_names) + '\n')
                for name, contents in files.items():
                    (output / name).write_text(contents)
                result = subprocess.run(['bash', '-e', '-c', script, 'collector-test', str(tools)],
                                        cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual((root / 'result.out').read_text().strip(),
                                 'CompareLRDensity_pass ' + expected)

    def test_autotest_rerun_clears_only_opted_in_lr_density_cubes(self):
        integrate = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'integrate').symlink_to(integrate, target_is_directory=True)
            suite = root / '08_RI'
            case = suite / 'density_rerun'
            output = case / 'OUT.autotest'
            output.mkdir(parents=True)
            (suite / 'CASES_CPU.txt').write_text('density_rerun\n')
            (case / 'INPUT').write_text('calculation scf\nesolver_type ks-lr\nnspin 1\nout_chg 1\n')
            (case / 'result.ref').write_text('CompareLRDensity_pass 0\n')
            (case / 'lr_density.ref').write_text('LR_e-h_density_triplet_1.cube\n')
            (output / 'LR_e-h_density_triplet_1.cube').write_text('stale expected cube')
            (output / 'LR_e-h_density_singlet_1.cube').write_text('stale wrong label')
            charge = output / 'chg.cube'
            charge.write_text('preserve ground-state output')
            executable = root / 'fake-abacus'
            executable.write_text('''#!/usr/bin/env python3
from pathlib import Path
output = Path('OUT.autotest')
if Path('lr_density.ref').exists():
    assert not list(output.glob('LR_e-h_density_*.cube')), 'stale cubes reached executable'
assert (output / 'chg.cube').read_text() == 'preserve ground-state output'
(output / 'running_scf.log').write_text('atomic_world 0.01 0.01\\n')
if Path('emit-density').exists():
    (output / 'LR_e-h_density_triplet_1.cube').write_text('current run cube')
''')
            executable.chmod(0o755)
            command = ['bash', str(integrate / 'Autotest.sh'), '-a', str(executable),
                       '-n', '0', '-o', '1', '-j', '1']
            environment = dict(os.environ, OMP_NUM_THREADS='1')
            # First rerun removes both old labels and accepts newly generated cubes.
            emit = case / 'emit-density'
            emit.touch()
            result = subprocess.run(command, cwd=suite, env=environment, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((output / 'LR_e-h_density_triplet_1.cube').read_text(), 'current run cube')
            self.assertFalse((output / 'LR_e-h_density_singlet_1.cube').exists())
            # A run that emits nothing must not pass on the preceding run's cube.
            emit.unlink()
            result = subprocess.run(command, cwd=suite, env=environment, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('CompareLRDensity_pass 1', (case / 'result.out').read_text())
            self.assertFalse(list(output.glob('LR_e-h_density_*.cube')))
            self.assertEqual(charge.read_text(), 'preserve ground-state output')
            # No opt-in reference means no cleanup, even when using the same harness.
            (case / 'lr_density.ref').unlink()
            stale = output / 'LR_e-h_density_singlet_1.cube'
            stale.write_text('unmanaged output')
            result = subprocess.run(command + ['-p'], cwd=suite, env=environment,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(stale.read_text(), 'unmanaged output')


if __name__ == '__main__':
    unittest.main()
