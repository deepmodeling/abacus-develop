import unittest
import numpy as np
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from ase.optimize import BFGS
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus
from common import preserved_tmpdir, test_profile_paths, assert_scf_facts

class TestIonicRelaxationWithStress(unittest.TestCase):

    def test(self):
        command, pporb = test_profile_paths(here)

        silicon = bulk('Si', 'diamond', a=5.43)
        silicon.positions[1, 0] += 0.1
        aprof = AbacusProfile(
            command=command,
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        with preserved_tmpdir('abacus_ase_relax_') as tmpdir:
            abacus = Abacus(
                profile=aprof,
                directory=tmpdir,
                pseudopotentials=load_pseudo(pporb),
                basissets=load_orbital(pporb, efficiency=True),
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production,
                    'cal_force': True,
                    'cal_stress': True
                }
            )

            silicon.calc = abacus
            initial_positions = silicon.positions.copy()
            e = silicon.get_potential_energy()
            initial_forces = silicon.get_forces()
            self.assertTrue(np.isfinite(e))
            self.assertEqual(initial_forces.shape, (len(silicon), 3))
            self.assertGreater(np.max(np.linalg.norm(initial_forces, axis=1)), 0.05)
            opt = BFGS(silicon, logfile=None)
            opt.run(fmax=0.05, steps=1)
            final_energy = silicon.get_potential_energy()
            self.assertTrue(np.isfinite(final_energy))
            self.assertLessEqual(final_energy, e + 1e-4)
            self.assertEqual(opt.nsteps, 1)
            self.assertFalse(np.allclose(silicon.positions, initial_positions))
            assert_scf_facts(self, abacus, tmpdir, len(silicon), forces=True, stress=True)

if __name__ == '__main__':
    unittest.main()
