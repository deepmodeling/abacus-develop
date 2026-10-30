import os
import unittest
import numpy as np
from pathlib import Path
here = Path(__file__).parent
from ase.build import bulk
from abacuslite import AbacusProfile, Abacus
from abacuslite.utils.ksampling import kpathgen
from common import (preserved_tmpdir, test_profile_paths, assert_scf_facts,
                    assert_band_facts)

class TestSCFFollowedNSCF(unittest.TestCase):

    def setUp(self):
        '''if SeeK-path is not installed, skip the test'''
        try:
            import seekpath
        except ImportError:
            self.skipTest('seekpath is not installed')
        return super().setUp()

    def test(self):
        command, pporb = test_profile_paths(here)

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command=command,
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        with preserved_tmpdir('abacus_ase_band_') as tmpdir:
            abacus = Abacus(
                profile=aprof,
                directory=tmpdir,
                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                basissets={'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
                inp={
                    'basis_type': 'lcao',
                    'scf_thr': 1e-3, # fast for test, wrong for production,
                    'kspacing': 0.25
                }
            )

            silicon.calc = abacus
            _ = silicon.get_potential_energy()
            assert_scf_facts(self, abacus, tmpdir, len(silicon))

            # ============================================================= #

            kpath, kspecial = kpathgen(silicon)
            bandpath = silicon.cell.bandpath(path=kpath,
                                            npoints=5,
                                            special_points=kspecial)
            bscalc = silicon.calc.fixed_density(bandpath)
            silicon.calc = bscalc
            _ = silicon.get_potential_energy()

            bands = bscalc.band_structure()
            values = np.asarray(bands.energies)
            self.assertEqual(values.ndim, 3)
            self.assertEqual(values.shape[1], len(bandpath.kpts))
            self.assertGreater(values.shape[2], 0)
            self.assertTrue(np.all(np.isfinite(values)))
            output = Path(tmpdir) / 'OUT.ABACUS'
            raw_bands = output / ('running_nscf.log' if os.environ.get('ASE_ABACUS_FORMAT') == 'legacy'
                                  else 'eig_occ.txt')
            assert_band_facts(self, bands, raw_bands)
            bands.write(Path(tmpdir) / 'bandstructure.json')
            self.assertTrue((Path(tmpdir) / 'bandstructure.json').exists())

if __name__ == '__main__':
    unittest.main()
