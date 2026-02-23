import unittest
import tempfile
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from ase.optimize import BFGS
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus

class TestIonicRelaxationWithStress(unittest.TestCase):

    def test(self):
        pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
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
        e = silicon.get_potential_energy()
        opt = BFGS(silicon)
        opt.run(fmax=0.05, steps=1)
        self.assertLessEqual(silicon.get_potential_energy(), e)
        self.assertIsNotNone(silicon.calc.results['forces'])
        self.assertIsNotNone(silicon.calc.results['stress'])

if __name__ == '__main__':
    unittest.main()