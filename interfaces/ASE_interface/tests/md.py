import unittest
import numpy as np
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from ase.units import fs
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus
from common import preserved_tmpdir, test_profile_paths, assert_scf_facts

class TestLangevinMolecularDynamics(unittest.TestCase):

    def test(self):
        command, pporb = test_profile_paths(here)

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command=command,
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        with preserved_tmpdir('abacus_ase_md_') as tmpdir:
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
                    'cal_stress': True,
                    'init_chg': 'auto'
                }
            )

            silicon.calc = abacus
            MaxwellBoltzmannDistribution(silicon, temperature_K=300)
            positions = silicon.positions.copy()
            frames = []
            dyn = Langevin(silicon,
                           timestep=1.0 * fs,
                           temperature_K=300,
                           friction=0.01)
            dyn.attach(lambda: frames.append((silicon.positions.copy(), silicon.get_forces().copy())))
            dyn.run(2)
            self.assertEqual(dyn.nsteps, 2)
            self.assertGreaterEqual(len(frames), 3)
            self.assertFalse(np.allclose(silicon.positions, positions))
            for frame_positions, frame_forces in frames:
                self.assertEqual(frame_forces.shape, (len(silicon), 3))
                self.assertTrue(np.all(np.isfinite(frame_positions)))
                self.assertTrue(np.all(np.isfinite(frame_forces)))
            assert_scf_facts(self, abacus, tmpdir, len(silicon), forces=True, stress=True)

if __name__ == '__main__':
    unittest.main()
