import unittest
import numpy as np
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from abacuslite import AbacusProfile, Abacus
from common import preserved_tmpdir, test_profile_paths, assert_scf_facts

class TestSCF(unittest.TestCase):

    def test(self):
        command, pporb = test_profile_paths(here)

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command=command,
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        # preserved_tmpdir keeps the run dir under ASE_ABACUS_KEEP_DIR on
        # failure so remote-CI-only failures (issue #7794) can be root-caused.
        with preserved_tmpdir('abacus_ase_scf_') as tmpdir:
            abacus = Abacus(
                profile=aprof,
                directory=tmpdir,
                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                basissets={'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production
                    'cal_force': True,
                    'cal_stress': True,
                }
            )

            silicon.calc = abacus
            energy = silicon.get_potential_energy()
            print('Silicon :', energy)
            self.assertTrue(np.isfinite(energy))
            assert_scf_facts(self, abacus, tmpdir, len(silicon),
                             forces=True, stress=True)
            self.assertEqual(abacus.results['natoms'], len(silicon))

if __name__ == '__main__':
    unittest.main()
