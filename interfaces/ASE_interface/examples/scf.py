'''
This example shows how to run a SCF calculation with ABACUS 
of Si diamond structure.
'''
import shutil
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

from ase.build import bulk
from abacuslite import Abacus, AbacusProfile

# AbacusProfile: the interface connecting the Abacus calculator instance
# with the file system and the enviroment
aprof = AbacusProfile(
    command='mpirun -np 4 abacus',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)

# Abacus: the calculator instance
jobdir = here / 'scf'
abacus = Abacus(
    profile=aprof,
    directory=str(jobdir),
    pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
    basissets={'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
    inp={
        'calculation': 'scf',
        'nspin': 1,
        'basis_type': 'lcao',
        'ks_solver': 'genelpa',
        'ecutwfc': 100,
        'symmetry': 1,
        'kspacing': 0.1
    }
)

# get the structure, can also from the 
# ```
# from ase.io import read
# atoms = read(...)
# ```
atoms = bulk('Si', 'diamond', a=5.43)

# bind the atoms with the abacus
atoms.calc = abacus

# calculate!
print(atoms.get_potential_energy())

# remove the temporary job directory (including all files inside)
shutil.rmtree(jobdir)
