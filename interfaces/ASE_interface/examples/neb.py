'''
PbTiO3 ferroelectric inversion energy barrier

Learn how to use the NEB module in ASE, please refer to the online manual at:
https://ase-lib.org/examples_generated/tutorials/neb_idpp.html
'''
from pathlib import Path
here = Path(__file__).parent

import numpy as np
from ase.io import Trajectory
from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.mep import NEB
import matplotlib.pyplot as plt
from abacuslite import Abacus, AbacusProfile

pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

elem = ['Ti', 'Pb', 'O', 'O', 'O']
taud = np.array([
    [0.5, 0.5, 0.5948316037314115],
    [0.0, 0.0, 0.1235879499999999],
    [0.0, 0.5, 0.5094847864489368],
    [0.5, 0.0, 0.5094847864489368],
    [0.5, 0.5, 0.0088672395150394],
])
cell = np.array([
   [3.8795519, 0.0000000, 0.00000000],
   [0.0000000, 3.8795519, 0.00000000],
   [0.0000000, 0.0000000, 4.28588762],
])

# we have relaxed with the parameters above :)
up = Atoms(elem, cell=cell, scaled_positions=taud)

# get the polarisation inversed by inversing the Ti atoms
taud = np.array([
    [0.5, 0.5, 0.6508136593687969],
    [0.0, 0.0, 0.1235879499999999],
    [0.0, 0.5, 0.7348401327639794],
    [0.5, 0.0, 0.7348401327639794],
    [0.5, 0.5, 0.2364165087650052],
])
dw = Atoms(elem, cell=cell, scaled_positions=taud)

aprof = AbacusProfile(
    command='mpirun -np 8 abacus_2p',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1
)
pseudopotentials = {
    'Ti': 'Ti_ONCV_PBE-1.0.upf',
    'Pb': 'Pb_ONCV_PBE-1.0.upf',
    'O' : 'O_ONCV_PBE-1.0.upf',
}
basissets = {
    'Ti': 'Ti_gga_8au_100Ry_4s2p2d1f.orb',
    'Pb': 'Pb_gga_7au_100Ry_2s2p2d1f.orb',
    'O' : 'O_gga_7au_100Ry_2s2p1d.orb',
}
inp = {
    'profile': aprof,
    'pseudopotentials': pseudopotentials,
    'basissets': basissets,
    'inp': {
        'basis_type': 'lcao',
        'symmetry': 1,
        'kspacing': 0.25, # Oops!
        'init_chg': 'auto',
        'cal_force': 1,
    }
}

n_replica = 7 # the ini and fin images included. 7 is acceptable for production
replica = []
for irep in range(n_replica):
    image = up.copy() if irep <= (n_replica // 2) else dw.copy()
    # attach the calculator to each image, so that we can run the optimization
    image.calc = Abacus(**inp, directory=here / f'neb-{irep}')
    replica.append(image)

neb = NEB(replica, 
          k=0.05, # too high value is hard to converge
          climb=False, # use True in production run, though CI-NEB is harder to converge
          parallel=True)
neb.interpolate('idpp')

qn = FIRE(neb, trajectory=here / 'neb.traj')
qn.run(fmax=0.05)

energies = []
# get the energy profile along the reaction path
with Trajectory(here / 'neb.traj') as traj:
    replica = traj[-7:] # the last NEB frames
    for i, rep in enumerate(replica):
        rep: Atoms # type hint
        # the energies of the initial and the final state
        # are not calculated, here we calculate them
        rep.calc = Abacus(**inp, directory=here / f'neb-{i}')
        energies.append(rep.get_potential_energy())

energies = np.array(energies)
# plot the energy profile
plt.plot(energies - energies[0], 'o-')
plt.xlabel('NEB image index')
plt.ylabel('Total energies (eV)')
plt.title('Energy profile along the reaction path')
plt.savefig(here / 'energy_profile.png', dpi=300)
plt.close()
