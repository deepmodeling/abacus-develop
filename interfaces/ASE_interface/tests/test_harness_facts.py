"""Offline checks for independent engine facts and bounded socket cleanup."""

import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ase.calculators.socketio import SocketIOCalculator

from common import (assert_band_facts, assert_collinear_magmoms,
                    close_socket_bounded, preserved_tmpdir)


class HarnessFactsTest(unittest.TestCase):
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

    def test_hung_socket_client_is_stopped_and_artifacts_retained(self):
        with tempfile.TemporaryDirectory() as root:
            keep = Path(root) / 'kept'
            client = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])
            unrelated = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])
            try:
                start = time.monotonic()
                with patch.dict(os.environ, ASE_ABACUS_KEEP_DIR=str(keep)):
                    with self.assertRaisesRegex(TimeoutError, 'socket client'):
                        with preserved_tmpdir('abacus_ase_socket_') as run:
                            Path(run, 'sentinel.txt').write_text('diagnostic')
                            calc = SocketIOCalculator(
                                unixsocket=f'abacus_harness_{os.getpid()}',
                                timeout=1, launch_client=lambda *args, **kwargs: None)
                            calc.server = calc.launch_server()
                            calc.server.proc = client
                            close_socket_bounded(calc, wait_seconds=0.1,
                                                 terminate_seconds=0.1)
                self.assertLess(time.monotonic() - start, 3.0)
                self.assertIsNotNone(client.poll())
                self.assertIsNone(unrelated.poll())
                self.assertEqual([p.read_text() for p in keep.glob('*/sentinel.txt')],
                                 ['diagnostic'])
            finally:
                for proc in (client, unrelated):
                    if proc.poll() is None:
                        proc.terminate()
                    proc.wait(timeout=3)
