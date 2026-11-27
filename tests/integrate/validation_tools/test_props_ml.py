#!/usr/bin/env python3
"""Energy regressions remain active when optical postprocessing is disabled."""
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class LREnergyCollectorTest(unittest.TestCase):
    def collect(self, log, analysis=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'OUT.autotest').mkdir()
            (root / 'running.log').write_text(log)
            if analysis is not None:
                (root / 'OUT.autotest/trans_analysis_singlet_tda.dat').write_text(analysis)
            script = Path(__file__).with_name('props_ml.sh')
            command = ('source ' + shlex.quote(str(script.resolve()))
                       + '; is_lr=1; props_result_file=result.out; running_path=running.log; run_ml_lr_props')
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


if __name__ == '__main__':
    unittest.main()
