import unittest
import tempfile
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from ase.units import fs
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus

class TestLangevinMolecularDynamics(unittest.TestCase):

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
                    'cal_stress': True,
                    'init_chg': 'auto'
                }
            )

        silicon.calc = abacus
        MaxwellBoltzmannDistribution(silicon, temperature_K=300)
        dyn = Langevin(silicon, 
                       timestep=1.0 * fs, 
                       temperature_K=300, 
                       friction=0.01)
        dyn.run(2)

if __name__ == '__main__':
    unittest.main()