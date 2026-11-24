#!/usr/bin/env python3
"""Exercise the LR energy and absorption property collector."""
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from check_lr_absorption import HEADER, SPECTRUM_POINTS


REFERENCE = (f'# Ionic_Step 1\n# {SPECTRUM_POINTS} # number of points\n'
             + HEADER + '\n' + '3.00 413.28 2.5\n' * SPECTRUM_POINTS)
SPECTRUM = '# ABACUS version: test\n# Written at 2026-01-01 00:00:00\n' + REFERENCE


class LREnergyCollectorTest(unittest.TestCase):
    def collect(self, log, analysis=None, reference=None, spectrum=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'OUT.autotest').mkdir()
            (root / 'running.log').write_text(log)
            if analysis is not None:
                (root / 'OUT.autotest/trans_analysis_singlet_tda.dat').write_text(analysis)
            if reference is not None:
                (root / 'absorption.dat.ref').write_text(reference)
            if spectrum is not None:
                (root / 'OUT.autotest/absorption.dat').write_text(spectrum)
            tools_dir = Path(__file__).resolve().parent
            commands = [f'source {shlex.quote(str(tools_dir / script))}'
                        for script in ('props_common.sh', 'props_ml.sh')]
            commands += ['is_lr=1', 'props_result_file=result.out',
                         'running_path=running.log', 'run_ml_lr_props']
            command = '; '.join(commands)
            subprocess.run(['bash', '-e', '-c', command], cwd=root, check=True)
            return (root / 'result.out').read_text()

    def test_no_spectrum_still_collects_requested_roots_in_both_channels(self):
        log = ''' LR excitation energies: singlet
 State Excitation Energy (Ry, eV)
0 1.53173 20.8402
1 1.53242 20.8495
 End LR excitation energies
Forces (-gradients) of each excited state:
0 0 force 99 99 99
 LR excitation energies: triplet
 State Excitation Energy (Ry, eV)
0 1.35365 18.4173
1 1.35467 18.4312
 End LR excitation energies
'''
        self.assertEqual(self.collect(log), 'excitationenergyref1 1.531730\n'
                         'excitationenergyref2 1.532420\nexcitationenergyref3 1.353650\n'
                         'excitationenergyref4 1.354670\n')

    def test_existing_optical_analysis_keeps_precedence(self):
        log = 'SCF and LR solve finished\n'
        analysis = 'State Excitation Energy (Ry, eV)\n0 1.53173 20.8402\nState Occupied orbital\n'
        self.assertEqual(self.collect(log, analysis), 'excitationenergyref1 1.531730\n')

    def test_no_spectrum_ignores_stale_optical_files(self):
        log = ' LR excitation energies: singlet\n0 1.53173 20.8402\n End LR excitation energies\n'
        analysis = 'State Excitation Energy (Ry, eV)\n0 99 99\nState Occupied orbital\n'
        self.assertEqual(self.collect(log, analysis), 'excitationenergyref1 1.531730\n')

    def test_absorption_accepts_complete_finite_table(self):
        result = self.collect('', reference=REFERENCE, spectrum=SPECTRUM)
        self.assertEqual(result, 'CompareLRAbsorption_pass 0\n')

    def test_absorption_checks_provenance_step_and_point_count(self):
        no_timestamp = SPECTRUM.replace('# Written at 2026-01-01 00:00:00\n', '')
        wrong_step = SPECTRUM.replace('Ionic_Step 1', 'Ionic_Step 0')
        wrong_count = SPECTRUM.replace('# 100 # number of points', '# 99 # number of points')
        for spectrum in (REFERENCE, no_timestamp, wrong_step, wrong_count):
            with self.subTest(metadata=spectrum.splitlines()[:5]):
                result = self.collect('', reference=REFERENCE, spectrum=spectrum)
                self.assertEqual(result, 'CompareLRAbsorption_pass 1\n')

    def test_absorption_rejects_empty_and_incomplete_table(self):
        truncated = '\n'.join(SPECTRUM.splitlines()[:-1]) + '\n'
        for reference, spectrum in [(SPECTRUM, None), (SPECTRUM, ''), (SPECTRUM, truncated),
                                    (truncated, truncated)]:
            with self.subTest(reference_rows=len(reference.splitlines()),
                              output_rows=len(spectrum.splitlines()) if spectrum is not None else None):
                result = self.collect('', reference=reference, spectrum=spectrum)
                self.assertEqual(result, 'CompareLRAbsorption_pass 1\n')

    def test_absorption_rejects_nonfinite_values_in_all_columns(self):
        for value, invalid in [('3.00', 'NaN'), ('413.28', 'inf'), ('2.5', '-inf')]:
            invalid_spectrum = SPECTRUM.replace(value, invalid, 1)
            for reference in (SPECTRUM, invalid_spectrum):
                with self.subTest(invalid=invalid, matching=(reference == invalid_spectrum)):
                    result = self.collect('', reference=reference, spectrum=invalid_spectrum)
                    self.assertEqual(result, 'CompareLRAbsorption_pass 1\n')

    def test_absorption_preserves_weak_intensity(self):
        reference = SPECTRUM.replace('2.5', '2.5e-6')
        nearby = reference.replace('2.5e-6', '2.505e-6', 1)
        zeroed = reference.replace('2.5e-6', '0', 1)
        self.assertEqual(self.collect('', reference=reference, spectrum=nearby),
                         'CompareLRAbsorption_pass 0\n')
        self.assertEqual(self.collect('', reference=reference, spectrum=zeroed),
                         'CompareLRAbsorption_pass 1\n')

    def test_absorption_keeps_absolute_tolerances_on_both_grid_columns(self):
        changed_frequency = SPECTRUM.replace('3.00', '3.0002', 1)
        changed_wavelength = SPECTRUM.replace('413.28', '413.2802', 1)
        for spectrum in (changed_frequency, changed_wavelength):
            with self.subTest(spectrum_row=spectrum.splitlines()[5]):
                result = self.collect('', reference=SPECTRUM, spectrum=spectrum)
                self.assertEqual(result, 'CompareLRAbsorption_pass 1\n')

    def test_absorption_rejects_wrong_columns_and_changed_values(self):
        wrong_columns = SPECTRUM.replace('3.00 413.28 2.5', '3.00 413.28', 1)
        changed = SPECTRUM.replace('3.00 413.28 2.5', '3.00 413.28 2.6', 1)
        for spectrum in (wrong_columns, changed):
            with self.subTest(spectrum_row=spectrum.splitlines()[5]):
                result = self.collect('', reference=SPECTRUM, spectrum=spectrum)
                self.assertEqual(result, 'CompareLRAbsorption_pass 1\n')


if __name__ == '__main__':
    unittest.main()
