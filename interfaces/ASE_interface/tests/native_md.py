"""Run short ABACUS-native MD and check its ASE trajectory against MD_dump."""

import os
import unittest
from pathlib import Path

import numpy as np
from ase.build import bulk

from abacuslite import Abacus, AbacusProfile
from abacuslite.io.generalio import load_orbital, load_pseudo
from common import (assert_native_md_facts, assert_stru_velocity_facts,
                    preserved_tmpdir, test_profile_paths)


class TestNativeMD(unittest.TestCase):
    def test_two_steps(self):
        command, pporb = test_profile_paths(Path(__file__).parent)
        silicon = bulk('Si', 'diamond', a=5.43)
        silicon.set_velocities(np.array([[0.01, 0.0, 0.0],
                                         [-0.01, 0.0, 0.0]]))
        initial_velocities = silicon.get_velocities().copy()
        legacy = os.environ.get('ASE_ABACUS_FORMAT', 'latest') == 'legacy'
        if legacy:
            from abacuslite.io.legacyio import read_abacus_out
        else:
            from abacuslite.io.latestio import read_abacus_out
        profile = AbacusProfile(command=command, pseudo_dir=pporb,
                                orbital_dir=pporb, omp_num_threads=1)
        with preserved_tmpdir('abacus_ase_native_md_') as run:
            calculator = Abacus(
                profile=profile,
                directory=run,
                pseudopotentials=load_pseudo(pporb),
                basissets=load_orbital(pporb, efficiency=True),
                inp={
                    'calculation': 'md',
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3,
                    'cal_force': True,
                    'md_nstep': 2,
                    'md_tchain': 4,
                    'md_tfirst': 10,
                    'md_tfreq': 0.1,
                    'init_vel': 1,
                },
            )
            silicon.calc = calculator
            self.assertTrue(np.isfinite(silicon.get_potential_energy()))
            assert_stru_velocity_facts(self, Path(run) / 'STRU', initial_velocities)
            output = Path(run) / 'OUT.ABACUS'
            frames = read_abacus_out(output / 'running_md.log')
            assert_native_md_facts(self, frames, output / 'MD_dump',
                                   expected_frames=3 if legacy else 2)


if __name__ == '__main__':
    unittest.main()
