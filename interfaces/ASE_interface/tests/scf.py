import os
import shutil
import unittest
import tempfile
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from abacuslite.io.generalio import load_pseudo, load_orbital
from abacuslite import AbacusProfile, Abacus

class TestSCF(unittest.TestCase):

    def test(self):
        pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        # Use mkdtemp + explicit cleanup so that on failure we can preserve
        # the whole run directory for post-mortem analysis (issue #7794 is
        # only reproducible on remote CI; losing the run dir there makes
        # root-causing impossible).
        tmpdir = tempfile.mkdtemp(prefix='abacus_ase_scf_')
        try:
            abacus = Abacus(
                profile=aprof,
                directory=tmpdir,
                pseudopotentials=load_pseudo(pporb),
                basissets=load_orbital(pporb, efficiency=True),
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production
                }
            )

            silicon.calc = abacus
            print('Silicon :', silicon.get_potential_energy())
        except Exception:
            keep_root = Path(os.environ.get(
                'ASE_ABACUS_KEEP_DIR', '/tmp/abacus_ase_failure'))
            keep_root.mkdir(parents=True, exist_ok=True)
            keep = keep_root / Path(tmpdir).name
            shutil.copytree(tmpdir, keep, dirs_exist_ok=True)
            print(f'[ase-test] failure, run dir preserved at {keep}')
            raise
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == '__main__':
    unittest.main()