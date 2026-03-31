'''this example shows how to perform the noncolinear
spin-orbit coupling calculation'''
import shutil
from pathlib import Path
here = Path(__file__).parent

pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

import numpy as np
from ase.atoms import Atoms
from abacuslite import Abacus, AbacusProfile

'''SPECIAL: ase does not support the noncolinear spin yet
till 2026/3/24, see ase/outputs.py:L154-155, in which the 
magmom cannot be set as the vector, so we release the 
datatype of magmom and magmoms by ourself'''
from ase.outputs import _defineprop, all_outputs
del all_outputs['magmom']
del all_outputs['magmoms']
_defineprop('magmom', float, shape=3) # re-define the magmom can be set as the vector
_defineprop('magmoms', float, shape=('natoms', 3))

fe = Atoms(symbols=['Fe'] * 2,
           scaled_positions=[[0., 0., 0.], [0.5, 0.5, 0.5]],
           magmoms=np.array([[0, 0, 1], [0, 0, 1]]),
           cell=np.eye(3) * 5.2,
           pbc=True)

aprof = AbacusProfile(
        command='mpirun -np 2 abacus',
        pseudo_dir=pporb,
        orbital_dir=pporb,
        omp_num_threads=1
)

jobdir = here / 'soc'
abacus = Abacus(
    profile=aprof,
    directory=jobdir,
    pseudopotentials={'Fe': 'Fe.upf'},
    basissets={'Fe': 'Fe_gga_6au_100Ry_4s2p2d1f.orb'},
    inp={
        'basis_type': 'lcao',
        'nspin': 4,
        'lspinorb': True,
        'noncolin': True,
        'out_mul': True,
        'scf_nmax': 10 # not serious setting!
    },
    kpts={
        'mode': 'mp-sampling',
        'gamma-centered': True,
        'nk': (2, 2, 2),
        'kshift': (0, 0, 0)
    }
)

fe.calc = abacus
print(f'Fe (non-colinear): {fe.get_potential_energy()} eV')
for ife, m in enumerate(fe.calc.results['magmoms']):
    print(f'Fe{ife}: ({",".join([f"{mi:10.6f}" for mi in m])}) uB')

shutil.rmtree(jobdir)