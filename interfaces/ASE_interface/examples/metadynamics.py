'''
To run the metadynamics, you need to configure the plumed correctly:

conda install -c conda-forge plumed=2.8.2=mpi_openmpi_hb0545ae_0
conda install -c conda-forge py-plumed

You may be also interested in a newer version of plumed, a possible solution is
to search `plumed` at conda website:
https://anaconda.org/channels/conda-forge/packages/plumed/files
in which, up to 2026/02/25, the latest version is 
```
linux-64/plumed-2.9.2-mpi_openmpi_h02da92d_0.conda,
```
you can install it with:
conda install -c conda-forge plumed=2.9.2=mpi_openmpi_h02da92d_0

we do not recommend the default version of plumed installed by conda, which
is nompi-labelled, may cause segmentation fault error during the MTD run.

In this example, we will run a metadynamics simulation to explore the
substitution reaction, during which the famous Walden inversion 
(https://en.wikipedia.org/wiki/Walden_inversion) happens:
```
CH4 + F -> CH3F + H
```
. Because the F is a radical, we set nspin 2 throughout the simulation. 
The structure is like
    H   H
     \ /
  H - C      F·
      |
      H
, thus the best Collective Variable (CV) is the difference between the 
H-C bond length and the C-F bond length.

CSVR thermostat is used to maintain the temperature during the simulation.
'''

import shutil
import tempfile
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

from ase.io import read
from ase.md import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import (
    fs,
    kJ as _kJ,
    mol as _mol,
)
_ps = 1000 * fs
from ase.constraints import FixCartesian
from ase.calculators.plumed import Plumed
from abacuslite import AbacusProfile, Abacus

aprof = AbacusProfile(
    command='mpirun -np 16 abacus',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)
jobdir = here / 'constraintmd'
abacus = Abacus(
    profile=aprof,
    directory=str(jobdir),
    pseudopotentials={
        'C': 'C_ONCV_PBE-1.0.upf',
        'H': 'H_ONCV_PBE-1.0.upf',
        'F': 'F_ONCV_PBE-1.0.upf',
    },
    basissets={
        'C': 'C_gga_8au_100Ry_2s2p1d.orb',
        'H': 'H_gga_8au_100Ry_2s1p.orb',
        'F': 'F_gga_7au_100Ry_2s2p1d.orb',
    },
    inp={
        'calculation': 'scf',
        'nspin': 2,
        'basis_type': 'lcao',
        'ks_solver': 'genelpa',
        'ecutwfc': 40,
        'scf_thr': 1e-6,
        'symmetry': 1,
        'gamma_only': True,
        'init_chg': 'auto' # small trick, use the previous charge density
    }
)

ch4f = ''' 6
ABACUS ASE Plugin Metadynamics example structure
 C                 -0.00000000   -0.00000000   -1.50074532
 F                  0.00000000   -0.00000000    0.88279136
 H                  0.00000000   -1.01539888   -1.16330635
 H                 -0.87936122    0.50769944   -1.16330635
 H                  0.87936122    0.50769944   -1.16330635
 H                 -0.00000000   -0.00000000   -2.57055225
'''
with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz') as f:
    f.write(ch4f)
    f.flush()
    atoms = read(f.name)

atoms.center(vacuum=5.0) # to reduce the computational cost

# constraint the No.5, 6 atoms' X and Y coordiantes so that
# they can only move along the z-axis, also fix the atom C's
# all components
atoms.set_constraint([FixCartesian(a=[4, 5], mask=(True, True, False)), 
                      FixCartesian(a=[0])])
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

setup = [# define the unit within the PLUMED runtime
         f'UNITS LENGTH=A TIME={1/_ps} ENERGY={_mol/_kJ}',
         # define the two bond lengths
         'd1: DISTANCE ATOMS=1,5',
         'd2: DISTANCE ATOMS=1,6',
         # define the CV as the difference between the two bond lengths
         'c1: MATHEVAL ARG=d1,d2 VAR=a,b FUNC=a-b PERIODIC=NO',
         # add walls to confine the position of H and F atoms
         # such that the C-H bond will have length between 0.5 and 2.0,
         # and the C-F bond will have length between 1.0 and 3.0
         'lwall: LOWER_WALLS ARG=d1,d2 AT=0.5,1.0 KAPPA=150.0,150.0 EXP=2,2',
         'uwall: UPPER_WALLS ARG=d1,d2 AT=2.0,3.0 KAPPA=150.0,150.0 EXP=2,2',
         # setup the metadynamics simulation
         'metad: METAD ARG=c1 PACE=5 HEIGHT=0.2 SIGMA=0.05 FILE=HILLS TEMP=300',
         'PRINT STRIDE=1 ARG=d1,d2,c1 FILE=COLVAR']

atoms.calc = Plumed(calc=abacus,
                    input=setup,
                    timestep=1.0 * fs, 
                    atoms=atoms,
                    kT=0.1)

dyn = Bussi(atoms, 
            timestep=1.0 * fs, 
            temperature_K=300, 
            taut=10.0 * fs,
            trajectory='metadynamics.traj',
            logfile='-')

dyn.run(20)
shutil.rmtree(jobdir)
