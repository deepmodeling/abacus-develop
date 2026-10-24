"""Bounded real fixed-cell socket handoff against raw FileIO calculations."""

import os
import signal
import unittest
from pathlib import Path

import numpy as np
from ase import Atoms
from abacuslite import Abacus, AbacusProfile, AbacusSocketIO

from common import (assert_scf_facts, calculate_socket_frame, close_socket_checked,
                    preserved_tmpdir, test_profile_paths)


class TestSocketHandoff(unittest.TestCase):
    def test_two_positions(self):
        command, assets = test_profile_paths(Path(__file__).parent)
        profile = AbacusProfile(command=command, pseudo_dir=assets,
                                orbital_dir=assets, omp_num_threads=1)
        inputs = dict(
            pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
            basissets={'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
            inp={'calculation': 'scf', 'basis_type': 'lcao', 'gamma_only': True,
                 'scf_thr': 1e-6, 'scf_nmax': 40, 'cal_force': 1, 'cal_stress': 1},
        )
        atoms = Atoms('Si2', positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
                      cell=[5.43, 5.43, 5.43], pbc=True)
        with preserved_tmpdir('abacus_ase_socket_') as tmpdir:
            previous_handler = signal.getsignal(signal.SIGALRM)
            def timed_out(*_):
                raise TimeoutError('socket test exceeded 360 s')
            signal.signal(signal.SIGALRM, timed_out)
            signal.alarm(360)
            try:
                self._run_handoff(profile, inputs, atoms, tmpdir)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, previous_handler)

    def _run_handoff(self, profile, inputs, atoms, tmpdir):
        root = Path(tmpdir)
        socket = AbacusSocketIO(profile=profile, directory=root / 'socket',
                                unixsocket=f'abacus_ci_{os.getpid()}',
                                timeout=120, **inputs)
        try:
            atoms.calc = socket
            for index, displacement in enumerate((0.0, 0.03)):
                atoms.positions[1] = [1.25 + displacement] * 3
                calculate_socket_frame(socket, atoms)
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces().copy()
                stress = atoms.get_stress().copy()
                self.assertTrue(np.isfinite(energy))
                self.assertEqual(forces.shape, (2, 3))
                self.assertEqual(stress.shape, (6,))
                self.assertTrue(np.all(np.isfinite(forces)))
                self.assertTrue(np.all(np.isfinite(stress)))

                reference_dir = root / f'fileio_{index}'
                reference = Abacus(profile=profile, directory=reference_dir, **inputs)
                check_atoms = atoms.copy()
                check_atoms.calc = reference
                self.assertAlmostEqual(energy, check_atoms.get_potential_energy(), delta=1e-4)
                self.assertTrue(np.allclose(forces, check_atoms.get_forces(), atol=1e-4))
                self.assertTrue(np.allclose(stress, check_atoms.get_stress(), atol=1e-4))
                assert_scf_facts(self, reference, reference_dir, 2,
                                 forces=True, stress=True)
        finally:
            # Use the public lifecycle API. The outer runner bounds the whole
            # case (including close), kills its process group on failure and
            # retains raw files; no ASE server/process internals are inspected.
            close_socket_checked(socket)


if __name__ == '__main__':
    unittest.main()
