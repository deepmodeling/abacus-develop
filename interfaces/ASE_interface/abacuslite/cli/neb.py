'''
CLI calling example:
```bash
abacuslite neb -i INPUT 
               --kpts KPT
               --prefix neb
               --ini STRU1
               --fin STRU2
               --configure neb.yaml
```

YAML configuration file example:
```yaml
aprof:
  command: 'mpirun -np 8 abacus'
  pseudo_dir: './'
  orbital_dir: './'
  omp_num_threads: 1

neb:
  nreplica: 7
  kspring: 0.02
  fmax: 0.05

  type: 'neb' # neb, cineb, atst-autoneb, atst-autoneb-sella
  atst:
    autoneb:
      placeholder1: 'TBD'
      placeholder2: 'TBD'
  
  parallel-images: true
```
'''
from pathlib import Path
import unittest
import argparse
from typing import Dict

try:
    from yaml import safe_load
except ModuleNotFoundError:
    print('Please install yaml module by `pip install pyyaml`')
    raise

from ase.atoms import Atoms
from ase.optimize import FIRE, BFGS
from ase.optimize.optimize import Optimizer
from ase.mep import NEB

from abacuslite.io.generalio import read_stru
from abacuslite.cli.utils import (
    build_atoms_from_stru,
    build_aprof,
    build_calc
)

DESCRIPTION = '''Nudged Elastic Band (NEB) calculation
'''

def add_argument(parser: argparse.ArgumentParser):
    '''add argument for neb calculation'''
    parser.add_argument('-i', '--input', type=str, 
                        help='standard ABACUS INPUT file, to define the calculation '
                             'settings including `basis_type`, and so on. NOTE: '
                             'if you define the `kspacing` here, you are not '
                             'required to provide the KPT file via `kpts` flag. This '
                             'file MUST be provided',
                        required=True)
    parser.add_argument('--ini', type=str, 
                        help='initial structure file in STRU format. '
                             'This file MUST be provided',
                        required=True)
    parser.add_argument('--fin', type=str, 
                        help='final structure file in STRU format. '
                             'This file MUST be provided',
                        required=True)
    parser.add_argument('--kpts', type=str, 
                        help='k-points file, in ABACUS KPT format. If you set either '
                             'the `gamma_only` or `kspacing` in INPUT file, you can '
                             'ignore this flag',
                        required=False)
    parser.add_argument('--configure', type=str, 
                        help='controls the details of the neb calculation, including '
                             'the number of images (replica, the initial and the final'
                             ' structures included), band type (e.g., cineb, neb, '
                             'atst-autoneb, atst-autoneb-sella, etc.), etc. '
                             'This file MUST be provided',
                        required=True)

def check(conf):
    '''
    check the correctness of the neb calculation configuration
    '''
    assert isinstance(conf, Dict)
    NECESSARY_SECTIONS = ['io', 'aprof', 'neb']
    for section in NECESSARY_SECTIONS:
        assert section in conf, f'{section} section is required'
    
def run(finp, fini, ffin, fneb, fkpt):
    conf = safe_load(open(fneb))
    _ = check(conf)

    # special check for this function
    assert conf['neb']['type'] in ['neb', 'cineb']

    ini, psp, orb = build_atoms_from_stru(read_stru(fini))
    fin, _, _ = build_atoms_from_stru(read_stru(ffin))

    assert isinstance(ini, Atoms)
    assert isinstance(fin, Atoms)

    outdir = Path(conf['io']['outdir'])
    outdir.mkdir(parents=True, exist_ok=True)

    nrep = conf['neb']['nreplica']
    replica = []
    for irep in range(nrep):
        image = ini.copy() if irep <= (nrep // 2) else fin.copy()
        image.calc = build_calc(build_aprof(**conf['aprof']),
            pseudopotentials=psp, basissets=orb,
            finp=finp, fkpt=fkpt,
            directory=outdir / f'{conf["io"].get("prefix", "abacuslite-neb")}-{irep}')
        replica.append(image)

    nebparam: Dict = conf['neb']
    band = NEB(replica, 
               k=nebparam['kspring'], 
               climb=bool(nebparam['type'] == 'cineb'),
               parallel=nebparam.get('parallel-images', False))

    band.interpolate(method='linear' if not nebparam.get('idpp') else 'idpp')
    
    ftraj = outdir / f'{conf["io"].get("prefix", "abacuslite-neb")}.traj'
    runner: Optimizer = {'fire': FIRE, 'bfgs': BFGS}[nebparam['optimize']['method'].lower()](
        band, trajectory=ftraj)
    runner.run(fmax=nebparam['optimize']['fmax'],
               steps=nebparam['optimize'].get('nmax', 200))
    return ftraj

class TestNEBCli(unittest.TestCase):

    here = Path(__file__).parent
    testfiles = here / 'testfiles'

    def test_yaml_read(self):
        param = {
            'aprof': {
                'command': 'mpirun -np 8 abacus',
                'pseudo_dir': './',
                'orbital_dir': './',
                'omp_num_threads': 1,
            },
            'neb': {
                'nreplica': 7,
                'kspring': 0.02,
                'fmax': 0.05,
                'type': 'neb',
                'atst': {
                    'autoneb': {
                        'placeholder1': 'TBD',
                        'placeholder2': 'TBD',
                    },
                },
                'parallel-images': True,
            },
        }
        with open(self.testfiles / 'neb.yaml') as f:
            self.assertEqual(param, safe_load(f))

if __name__ == '__main__':
    unittest.main()