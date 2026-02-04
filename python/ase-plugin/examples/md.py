import shutil
import tempfile
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
from ase.io import read
from ase.md import Langevin
from ase.units import fs
from abacuslite import Abacus, AbacusProfile

# we write the following and read it back to an ASE.atoms.Atoms
# object
atoms = '''7  
100. 100. 100.       
Ar 7.3933470660       -2.6986483924        0.0000000000
Ar 7.8226765198       -0.7390907295        0.0000000000
Ar 7.1014969839       -1.6164766614        0.0000000000
Ar 8.2357184242       -1.7097824975        0.0000000000
Ar 6.7372520842       -0.5111536183        0.0000000000
Ar 6.3777119489       -2.4640437401        0.0000000000
Ar 5.9900631495       -1.3385375043        0.0000000000
'''
with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz') as f:
    f.write(atoms)
    f.flush()
    atoms = read(f.name)
atoms.set_cell([100, 100, 100])

aprof = AbacusProfile(
    command='mpirun -np 4 abacus',
    omp_num_threads=1,
)
jobdir = here / 'md'
abacus = Abacus(
    profile=aprof,
    directory=str(jobdir),
    inp={
        'esolver_type': 'lj',
        'cal_force': 1 # let ABACUS calculate the forces
    }
)

atoms.calc = abacus
dyn = Langevin(atoms, 
               timestep=1*fs, 
               temperature_K=300, 
               friction=0.01,
               logfile='-') # let's see the trajectory
dyn.run(10)