'''this example recap the example of MD in abacus example
directory, but there the CSVR thermostat is used instead 
of those implemented in ABACUS.

In ASE, the CSVR thermostat is named as the Bussi
'''
import shutil
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

import numpy as np
from ase.atoms import Atoms
from ase.md import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs
from abacuslite import Abacus, AbacusProfile

cell = np.eye(3) * 5.43090251
taud = [
    [0.00000000, 0.00000000, 0.00000000],
    [0.00000000, 0.50000000, 0.50000000],
    [0.50000000, 0.00000000, 0.50000000],
    [0.50000000, 0.50000000, 0.00000000],
    [0.25000000, 0.25000000, 0.25000000],
    [0.25000000, 0.75000000, 0.75000000],
    [0.75000000, 0.25000000, 0.75000000],
    [0.75000000, 0.75000000, 0.25000000],
]
atoms = Atoms(symbols=['Si' for _ in range(8)],
              scaled_positions=taud,
              cell=cell,
              pbc=True)

aprof = AbacusProfile(
    command='mpirun -np 4 abacus',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)

jobdir = here / 'md'
abacus = Abacus(
    profile=aprof,
    directory=str(jobdir),
    pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
    basissets={'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
    inp={
        'calculation': 'scf', # still use SCF here because the MD is driven by ASE
        'nspin': 1,
        'basis_type': 'lcao',
        'ks_solver': 'genelpa',
        'ecutwfc': 100,
        'symmetry': 1,
        'kspacing': 0.25 # highly unconverged, just for demo
    }
)

atoms.calc = abacus
# initialize the velocities, necessary for CSVR
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

dyn = Bussi(atoms, 
            timestep=1*fs, 
            temperature_K=300, 
            taut=10*fs,
            logfile='-') # let's see the trajectory
dyn.run(2)

shutil.rmtree(jobdir)
