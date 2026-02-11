'''this example shows how to use the IDPP method interpolate the 
nudged elastic band (NEB) calculations such that the transition
state can be founded more quickly.

In this example, we also use the reaction from the metadynamics.py
, the Walden inversion
'''

import shutil
import tempfile
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

from ase.io import read
from ase.atoms import Atoms
from ase.optimize import BFGS, FIRE
from ase.mep import NEB
from ase.constraints import FixCartesian
from abacuslite import AbacusProfile, Abacus

aprof = AbacusProfile(
    command='mpirun -np 16 abacus_2p',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)
PSEUDOPOTENTIALS = {
    'C': 'C_ONCV_PBE-1.0.upf',
    'H': 'H_ONCV_PBE-1.0.upf',
    'F': 'F_ONCV_PBE-1.0.upf',
}
BASISSETS = {
    'C': 'C_gga_8au_100Ry_2s2p1d.orb',
    'H': 'H_gga_8au_100Ry_2s1p.orb',
    'F': 'F_gga_7au_100Ry_2s2p1d.orb',
}
DFTPARAM = {
    'calculation': 'scf',
    'nspin': 2, # because we are studying a reaction involving radical
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'ecutwfc': 60, # good for example, wrong for production
    'symmetry': 1,
    'gamma_only': True,
    'init_chg': 'auto'
}

def impose_constraint(atoms: Atoms) -> Atoms:
    '''imposing the constraint such that the C atom is fixed,
    H and F that leaves and arrives are fixed such that can only
    move along Z-axis'''
    # reset the constraint
    atoms.set_constraint(None)
    atoms.set_constraint([FixCartesian(a=[0]), 
                          FixCartesian(a=[4, 5], mask=(True, True, False))])
    return atoms

def relax(atoms: Atoms) -> Atoms:
    '''relax the structure to the minimum energy configuration'''
    # we do additional fully-constraint on the No.1 atom, and XY-constraint
    # on the last two atoms
    print('Relaxation of atoms:', atoms)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Abacus(
            profile=aprof,
            directory=tmpdir,
            pseudopotentials=PSEUDOPOTENTIALS,
            basissets=BASISSETS,
            inp=DFTPARAM
        )
        atoms.calc = calc
        opt = BFGS(atoms, logfile='-')
        opt.run(fmax=0.05)
    print('Done')
    return atoms

# relax the initial and final states
ch4f = ''' 6
Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3
 C                 -0.00000000   -0.00000000   -1.50074532
 H                  0.00000000   -1.01539888   -1.16330635
 H                 -0.87936122    0.50769944   -1.16330635
 H                  0.87936122    0.50769944   -1.16330635
 H                 -0.00000000   -0.00000000   -2.57055225
 F                  0.00000000   -0.00000000    0.88279136
'''
with tempfile.NamedTemporaryFile(mode='w', suffix='.extxyz') as f:
    f.write(ch4f)
    f.flush()
    ini = relax(impose_constraint(read(f.name)))

ch3fh = ''' 6
Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3
 C                 -0.00000000   -0.00000000   -1.50074532
 H                  0.00000000   -1.01539888   -1.16330635
 H                 -0.87936122    0.50769944   -1.16330635
 H                  0.87936122    0.50769944   -1.16330635
 H                 -0.00000000   -0.00000000   -3.57055225
 F                  0.00000000   -0.00000000   -0.50000000
'''
with tempfile.NamedTemporaryFile(mode='w', suffix='.extxyz') as f:
    f.write(ch3fh)
    f.flush()
    fin = relax(impose_constraint(read(f.name)))

# prepare images
n_replica, images = 5, []

for irep in range(n_replica):
    calc = Abacus(
        profile=aprof,
        directory=here / f'neb-idpp-replica-{irep}',
        pseudopotentials=PSEUDOPOTENTIALS,
        basissets=BASISSETS,
        inp=DFTPARAM
    )
    if irep == 0:
        ini.calc = calc
        images.append(impose_constraint(ini))
    elif irep == n_replica - 1:
        fin.calc = calc
        images.append(impose_constraint(fin))
    else:
        rep = ini.copy() if irep <= n_replica // 2 else fin.copy()
        rep.calc = calc
        images.append(impose_constraint(rep))

neb = NEB(images)
neb.interpolate('idpp')

# optimize the band
qn = FIRE(neb, trajectory='neb-idpp.traj', logfile='-')
qn.run(fmax=0.05)