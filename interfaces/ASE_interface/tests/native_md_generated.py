"""Run short ABACUS-native MD and check its ASE trajectory against MD_dump."""

import unittest
import re
from pathlib import Path

import numpy as np
from ase.build import bulk

from abacuslite import Abacus, AbacusProfile
from abacuslite.io.generalio import load_orbital, load_pseudo
from abacuslite.io.legacyio import read_abacus_out
from common import assert_native_md_facts, preserved_tmpdir, test_profile_paths


class TestNativeMD(unittest.TestCase):
    def test_two_steps(self):
        command, pporb = test_profile_paths(Path(__file__).parent)
        silicon = bulk('Si', 'diamond', a=5.43)
        # init_vel=0 must generate seeded engine velocities even though the
        # ASE structure carries no initial momenta (STRU v remains zero).
        self.assertFalse(silicon.has('momenta'))
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
                    'init_vel': 0,
                    'md_seed': 1729,
                },
            )
            silicon.calc = calculator
            self.assertTrue(np.isfinite(silicon.get_potential_energy()))
            stru = (Path(run) / 'STRU').read_text()
            raw_velocities = np.asarray(re.findall(
                r'\bv\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)', stru),
                dtype=float)
            self.assertEqual(raw_velocities.shape, (len(silicon), 3))
            self.assertTrue(np.all(raw_velocities == 0))
            output = Path(run) / 'OUT.ABACUS'
            frames = read_abacus_out(output / 'running_md.log')
            # LTS run_md.cpp records MDSTEP 0, 1 and 2 for md_nstep=2.
            assert_native_md_facts(self, frames, output / 'MD_dump', expected_frames=3)


if __name__ == '__main__':
    unittest.main()
