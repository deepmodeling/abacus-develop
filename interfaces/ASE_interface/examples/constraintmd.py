import shutil
import tempfile
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, Trajectory
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs
from ase.constraints import FixCartesian
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

mol = ''' 6
CH4-F
 C                 -0.00000000   -0.00000000   -1.50074532
 H                  0.00000000   -1.01539888   -1.16330635
 H                 -0.87936122    0.50769944   -1.16330635
 H                  0.87936122    0.50769944   -1.16330635
 H                 -0.00000000   -0.00000000   -2.57055225
 F                  0.00000000   -0.00000000    0.88279136
'''
with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz') as f:
    f.write(mol)
    f.flush()
    atoms = read(f.name)

atoms.center(vacuum=5.0) # to reduce the computational cost
# view(atoms)

# constraint the No.1, 5, 6 atoms' X and Y coordiantes so that
# they can only move along the z-axis
constraint = FixCartesian(a=[0, 4, 5], mask=(True, True, False))
# apply
atoms.set_constraint(constraint)
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

atoms.calc = abacus
dyn = Langevin(atoms, 
               timestep=1.0 * fs, 
               temperature_K=300, 
               friction=0.004,
               logfile='-',
               trajectory='constraintmd.traj')
dyn.run(5)

# let's see if the X, Y coordinates of No.1, 5, and 6 atoms are really
# fixed
with Trajectory('constraintmd.traj') as traj:
    traj = np.array([atoms.get_positions() for atoms in traj])

# transpose from (nframe, natom, 3) to (natom, nframe, 3)
traj = traj.transpose(1, 0, 2)

# plot the trajectory of No.1, 5, and 6 atoms
plt.plot(traj[0, :, 2], label='C')
plt.plot(traj[4, :, 2], label='H')
plt.plot(traj[5, :, 2], label='F')
plt.xlabel('Step')
plt.ylabel('Z (Å)')
plt.legend()
plt.show()

shutil.rmtree(jobdir)