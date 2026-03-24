'''
test if ABACUS can perform the calculation on the anti-ferromagnetic
and ferromagnetic phases calculation on the BCC Fe
'''
import unittest
import tempfile
import numpy as np
from pathlib import Path
here = Path(__file__).parent
from ase.atoms import Atoms
from abacuslite import AbacusProfile, Abacus

class TestMagneticSCF(unittest.TestCase):

    def setUp(self):
        pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'
        self.aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

    def test_antiferromagnetic(self):
        fe = Atoms(symbols=['Fe'] * 2,
                   scaled_positions=[[0., 0., 0.], [0.5, 0.5, 0.5]],
                   magmoms=[-1, 1],
                   cell=np.eye(3) * 5.2,
                   pbc=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            abacus = Abacus(
                profile=self.aprof,
                directory=tmpdir,
                pseudopotentials={'Fe': 'Fe_ONCV_PBE-1.0.upf'},
                basissets={'Fe': 'Fe_gga_6au_100Ry_4s2p2d1f.orb'},
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production
                    'nspin': 2,
                    'out_mul': True
                }
            )

            fe.calc = abacus
            print(f'Fe (AFM): {fe.get_potential_energy()} eV')
            for ife, m in enumerate(fe.calc.results['magmoms']):
                print(f'Fe{ife}: {m:10.6f} uB')
                self.assertGreater(abs(m), 1e-3)
            self.assertTrue(np.allclose(np.sum(fe.calc.results['magmoms']), 0.0, atol=1e-6))

    def test_ferromagnetic(self):

        fe = Atoms(symbols=['Fe'] * 2,
                   scaled_positions=[[0., 0., 0.], [0.5, 0.5, 0.5]],
                   magmoms=[1, 1],
                   cell=np.eye(3) * 5.2,
                   pbc=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            abacus = Abacus(
                profile=self.aprof,
                directory=tmpdir,
                pseudopotentials={'Fe': 'Fe_ONCV_PBE-1.0.upf'},
                basissets={'Fe': 'Fe_gga_6au_100Ry_4s2p2d1f.orb'},
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production
                    'nspin': 2,
                    'out_mul': True
                }
            )

            fe.calc = abacus
            print(f'Fe (FM): {fe.get_potential_energy()} eV')
            for ife, m in enumerate(fe.calc.results['magmoms']):
                print(f'Fe{ife}: {m:10.6f} uB')
                self.assertGreater(abs(m), 1e-3)
            self.assertGreater(abs(np.sum(fe.calc.results['magmoms'])), 0.0)

if __name__ == '__main__':
    unittest.main()