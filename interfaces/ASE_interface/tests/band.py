import unittest
import tempfile
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus
from abacuslite.utils.ksampling import kpathgen

class TestSCFFollowedNSCF(unittest.TestCase):

    def setUp(self):
        '''if SeeK-path is not installed, skip the test'''
        try:
            import seekpath
        except ImportError:
            self.skipTest('seekpath is not installed')
        return super().setUp()

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
                    'scf_thr': 1e-3, # fast for test, wrong for production,
                    'kspacing': 0.25
                }
            )

            silicon.calc = abacus
            _ = silicon.get_potential_energy()

            # ============================================================= #

            kpath, kspecial = kpathgen(silicon)
            bandpath = silicon.cell.bandpath(path=kpath,
                                            npoints=5,
                                            special_points=kspecial)
            bscalc = silicon.calc.fixed_density(bandpath)
            silicon.calc = bscalc
            _ = silicon.get_potential_energy()

            bscalc.band_structure().write(Path(tmpdir) / 'bandstructure.json')
            self.assertTrue((Path(tmpdir) / 'bandstructure.json').exists())

if __name__ == '__main__':
    unittest.main()