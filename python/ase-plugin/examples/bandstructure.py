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

import seekpath # for generating the kpoint path
import numpy as np
from ase.build import bulk
from ase.dft.kpoints import BandPath
from abacuslite import Abacus, AbacusProfile
from abacuslite.utils.ksampling import merge_ksgm, make_kstring

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

# then the band structure non-self-consistent calculation
kpathseen = seekpath.get_path(
    structure=(np.array(atoms.get_cell()), 
               atoms.get_scaled_positions(), 
               atoms.get_atomic_numbers()),
    with_time_reversal=True
)

# convert the kpoint path to the format that is acceptable by ASE
kpathstr, is_brkpt = merge_ksgm(kpathseen['path'])
fklblfilter = lambda lbl: lbl if lbl != 'GAMMA' else 'G'
kpathstr = make_kstring([fklblfilter(lbl) for lbl in kpathstr], is_brkpt)
kpathstr = ''.join([k if k != 'GAMMA' else 'G' for k in kpathstr])
# seekpath use 'GAMMA' as the gamma point symbol, while the ASE use 'G'
kspecial = {k: v for k, v in kpathseen['point_coords'].items() if k != 'GAMMA'}
kspecial['G'] = [0.0, 0.0, 0.0]

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