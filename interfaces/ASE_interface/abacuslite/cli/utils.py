'''the utilities that only used in cli mode'''
import uuid
import tempfile
from pathlib import Path
from typing import Dict, Tuple, List, Optional
from ase.atoms import Atoms
from ase.build import bulk
from ase.constraints import FixCartesian
import unittest
import numpy as np

from abacuslite import Abacus, AbacusProfile
from abacuslite.io.generalio import write_stru, read_stru, read_input, read_kpt

def build_atoms_from_stru(stru) -> Tuple[Atoms, Dict[str, str], Dict[str, str]]:
    '''
    cast the read_stru returned dictionary to the ase.atoms.Atoms object,
    optionally along with the mapping (correspondence) between the 
    element and the pseudopotential files and the orbital files

    Parameters
    ----------
    the stru : Dict
        the read_stru returned dictionary
    with_psp : bool, optional
        whether to return the mapping of pseudopotential files
    with_orb : bool, optional
        whether to return the mapping of the orbital files
    '''
    from ase.units import (
        Bohr as __BOHR__, 
        Angstrom as __ANGSTROM__
    )
    cell = np.array(stru['lat']['vec']).reshape(3, 3)
    cell *= stru['lat']['const'] # in bohr
    cell /= __ANGSTROM__ / __BOHR__ # in Angstrom
    
    elem, pos, mob, psp, orb = [], [], [], {}, {}
    for s in stru['species']:
        elem.extend([s['symbol']]*s['natom'])
        pos.extend([atom['coord'] for atom in s['atom']])
        mob.extend([atom['m'] for atom in s['atom']])
        psp[s['symbol']] = s['pp_file']
        orb[s['symbol']] = s.get('orb_file')
    # reshape
    pos = np.array(pos).reshape(-1, 3)
    constraints = [True if not m else False for m in np.array(mob, dtype=bool).flatten()]
    constraints = [FixCartesian(i, c) for i, c in enumerate(np.array(constraints).reshape(-1, 3))]

    # position transformation
    if stru['coord_type'].lower() == 'direct':
        pos = pos @ cell.T
    if stru['coord_type'].lower() in ['cartesian', 'cartesian_au']:
        pos *= stru['lat']['const'] # in bohr
        pos /= __ANGSTROM__ / __BOHR__ # in angstrom
    if stru['coord_type'].lower().startswith('cartesian_angstrom'):
        pass # do nothing

    return Atoms(symbols=elem, positions=pos, cell=cell, pbc=True,
                 constraint=constraints), psp, orb

def build_aprof(command,
                pseudo_dir,
                orbital_dir,
                omp_num_threads,
                **kwargs) -> AbacusProfile:
    return AbacusProfile(command, 
                         pseudo_dir, 
                         orbital_dir, 
                         omp_num_threads,
                         **kwargs)

def build_calc(aprof: AbacusProfile,
               pseudopotentials: Dict[str, str],
               basissets: Optional[Dict[str, str]],
               finp: Path | str,
               fkpt: Optional[Path | str],
               **kwargs) -> Abacus:
    inp = read_input(finp)
    inp = dict(profile=aprof,
               directory=kwargs.get('directory', str(uuid.uuid4().hex)),
               pseudopotentials=pseudopotentials,
               basissets=basissets,
               inp=inp
               )
    if 'directory' in kwargs:
        del kwargs['directory']
    if fkpt:
        inp['kpts'] = read_kpt(fkpt)
    return Abacus(**inp, **kwargs)

class TestCliUtils(unittest.TestCase):
    here = Path(__file__).parent
    testfiles = here / 'testfiles'
    def test_build_atoms_from_stru(self):
        nacl = bulk('NaCl', 'rocksalt', a=5.64)
        with tempfile.TemporaryDirectory() as tmpdir:
            fstru =write_stru(
                nacl, 
                outdir=tmpdir, 
                pp_file={
                    'Na': 'Na.pz-bhs.UPF',
                    'Cl': 'Cl.pz-bhs.UPF'
                },
                orb_file={
                    'Na': 'Na_gga_6au_100Ry_2s2p1d.orb',
                    'Cl': 'Cl_gga_6au_100Ry_2s2p1d.orb',
                },
            )
            stru = read_stru(fstru)
        # print(nacl.positions) # [[0.   0.   0.  ], [2.82 0.   0.  ]]
        # print(nacl.get_chemical_symbols()) # ['Na', 'Cl']
        nacl_, psp, orb = build_atoms_from_stru(stru)
        self.assertTrue(np.allclose(nacl.cell, nacl_.cell))
        self.assertTrue(np.allclose(nacl_.positions, 
                                    [[2.82, 0., 0.],
                                     [  0., 0., 0.]]))
        self.assertListEqual(nacl_.get_chemical_symbols(), ['Cl', 'Na'])
        self.assertDictEqual(psp, {'Na': 'Na.pz-bhs.UPF', 'Cl': 'Cl.pz-bhs.UPF'})
        self.assertDictEqual(orb, {'Na': 'Na_gga_6au_100Ry_2s2p1d.orb', 
                                   'Cl': 'Cl_gga_6au_100Ry_2s2p1d.orb'})

if __name__ == '__main__':
    unittest.main()