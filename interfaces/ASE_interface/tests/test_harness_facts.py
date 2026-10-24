"""Offline checks for independent engine facts and bounded socket cleanup."""

import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, PropertyNotImplementedError

from common import (assert_band_facts, assert_collinear_magmoms,
                    calculate_socket_frame, close_socket_checked)
from run_engine_compatibility import capture
from legacy_facts import assert_band_facts as assert_legacy_band_facts


class HarnessFactsTest(unittest.TestCase):
    def test_public_socket_close_rejects_nonzero_client_exit_warning(self):
        class PublicSocket:
            def close(self):
                warnings.warn('Subprocess exited with status 7', UserWarning)

        with self.assertRaisesRegex(UserWarning, 'status 7'):
            close_socket_checked(PublicSocket())

    def test_public_socket_close_accepts_clean_shutdown(self):
        closed = []
        socket = SimpleNamespace(close=lambda: closed.append(True))
        close_socket_checked(socket)
        self.assertEqual(closed, [True])

    def test_socket_requests_only_actual_changes_after_first_frame(self):
        class RecordingCalculator(Calculator):
            implemented_properties = ['energy', 'forces', 'stress']

            def __init__(self):
                super().__init__()
                self.changes = []

            def calculate(self, atoms=None, properties=None, system_changes=None):
                if self.atoms is not None and any(
                        name not in ('positions', 'cell') for name in system_changes):
                    raise PropertyNotImplementedError('immutable socket input changed')
                super().calculate(atoms, properties, system_changes)
                self.changes.append(list(system_changes))
                self.results = {'energy': 0.0, 'forces': np.zeros((2, 3)),
                                'stress': np.zeros(6)}

        atoms = Atoms('Si2', positions=[[0, 0, 0], [1, 1, 1]],
                      cell=[5, 5, 5], pbc=True)
        calc = RecordingCalculator()
        calculate_socket_frame(calc, atoms)
        atoms.positions[1, 0] += 0.03
        calculate_socket_frame(calc, atoms)
        self.assertIn('numbers', calc.changes[0])
        self.assertEqual(calc.changes[1], ['positions'])

    def test_band_eigenvalues_reject_mutated_reader_result(self):
        raw = '''1 # ionic step
spin=1 k-point=1/2 Cartesian=0.0 0.0 0.0 (1 plane wave)
1 -2.5000 1.0
2 1.2500 0.0
spin=1 k-point=2/2 Cartesian=0.5 0.0 0.0 (1 plane wave)
1 -2.2500 1.0
2 1.5000 0.0
'''
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'eig_occ.txt'
            path.write_text(raw)
            bands = SimpleNamespace(energies=np.array([[[-2.5, 1.25], [-2.25, 1.5]]]))
            assert_band_facts(self, bands, path)
            bands.energies[0, 1, 1] += 0.1
            with self.assertRaises(AssertionError):
                assert_band_facts(self, bands, path)

    def test_legacy_band_oracle_rejects_mutated_reader_result(self):
        raw = '''1/2 kpoint (Cartesian) = 0.0 0.0 0.0 (1 pws)
1 -2.5000 1.0
2 1.2500 0.0
2/2 kpoint (Cartesian) = 0.5 0.0 0.0 (1 pws)
1 -2.2500 1.0
2 1.5000 0.0
'''
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'running_nscf.log'
            path.write_text(raw)
            bands = SimpleNamespace(energies=np.array([[[-2.5, 1.25], [-2.25, 1.5]]]))
            assert_legacy_band_facts(self, bands, path)
            bands.energies[0, 1, 1] += 0.1
            with self.assertRaises(AssertionError):
                assert_legacy_band_facts(self, bands, path)

    def test_per_atom_magnetic_moment_rejects_swap(self):
        raw = '''Total Magnetism (uB)
--------------------
 Fe1 1.2715948597
 Fe2 -1.2715948597
--------------------
'''
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'running_scf.log'
            path.write_text(raw)
            calc = SimpleNamespace(results={'magmoms': np.array([1.2715948597, -1.2715948597])})
            assert_collinear_magmoms(self, calc, path, 2)
            calc.results['magmoms'] = calc.results['magmoms'][::-1]
            with self.assertRaises(AssertionError):
                assert_collinear_magmoms(self, calc, path, 2)

    def test_hung_public_close_is_stopped_without_touching_unrelated_process(self):
        # Model a client's public close() waiting forever. The runner owns the
        # process group; it must retain diagnostics and stop that group only.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = root / 'hung_close.py'
            script.write_text("""import subprocess, sys
from pathlib import Path
class Socket:
    def __init__(self):
        self.client = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])
        Path('client.pid').write_text(str(self.client.pid))
    def close(self):
        self.client.wait()
Path('raw-output.txt').write_text('retained engine output')
Socket().close()
""")
            unrelated = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])
            try:
                result = capture([sys.executable, str(script)], cwd=root,
                                 env=None, log=root / 'case.log', timeout=1)
                self.assertTrue(result['timed_out'])
                self.assertNotEqual(result['returncode'], 0)
                self.assertLess(result['seconds'], 4)
                self.assertIsNone(unrelated.poll())
                self.assertEqual((root / 'raw-output.txt').read_text(), 'retained engine output')
                child = int((root / 'client.pid').read_text())
                stat = Path(f'/proc/{child}/stat')
                if stat.exists():
                    self.assertEqual(stat.read_text().split()[2], 'Z')
            finally:
                unrelated.terminate()
                unrelated.wait(timeout=3)
