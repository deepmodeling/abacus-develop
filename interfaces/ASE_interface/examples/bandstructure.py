'''
This example shows how to run a SCF calculation with ABACUS 
of Si diamond structure.

To run this example, please install the SeeK-path package:
```
pip install seekpath
```
. The SeeK-path package recommands you cite by the way posted here:
https://seekpath.materialscloud.io/
'''
import shutil
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

from ase.build import bulk
from abacuslite import Abacus, AbacusProfile
from abacuslite.utils.ksampling import kpathgen

# AbacusProfile: the interface connecting the Abacus calculator instance
# with the file system and the enviroment
aprof = AbacusProfile(
    command='mpirun -np 4 abacus',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)

# Abacus: the calculator instance
jobdir = here / 'bandstructure'
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

# perform the SCF calculation to get the converged wavefunction
print('SCF calculation get energy:', atoms.get_potential_energy())

kpathstr, kspecial = kpathgen(atoms)
# instantiate the bandpath
bandpath = atoms.cell.bandpath(path=kpathstr,
                               npoints=50,
                               special_points=kspecial)

# derive the band structure calculator from SCF calculator
bscalc = atoms.calc.fixed_density(bandpath)
atoms.calc = bscalc
_ = atoms.get_potential_energy() # NSCF calculation will be performed

bs = bscalc.band_structure()
bs.write('bandstructure.json')
# you can use the ase-cli to plot the JSON file later by:
# ```
# ase band-structure bandstructure.json -r -10 15
# ```
bs.plot(emin=-10, emax=15, filename='bandstructure.png')

shutil.rmtree(jobdir)