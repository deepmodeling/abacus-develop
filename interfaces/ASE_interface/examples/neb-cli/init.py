'''
generate the necessary files for running the NEB calculation in cli mode
'''
from pathlib import Path
here = Path(__file__).parent

import numpy as np
from ase.atoms import Atoms

from abacuslite.io.generalio import write_stru, write_input

# find the PP_ORB directory
pporb = here.parent.parent.parent.parent / 'tests' / 'PP_ORB'
pporb = pporb.resolve()

# build the Atoms objects
elem = ['Ti', 'Pb', 'O', 'O', 'O']
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
write_stru(Atoms(elem, cell=cell, scaled_positions=taud), 
           outdir=here, pp_file=pseudopotentials, orb_file=basissets, fname='STRU1')

# get the polarisation inversed by inversing the Ti atoms
taud = np.array([
    [0.5, 0.5, 0.6508136593687969],
    [0.0, 0.0, 0.1235879499999999],
    [0.0, 0.5, 0.7348401327639794],
    [0.5, 0.0, 0.7348401327639794],
    [0.5, 0.5, 0.2364165087650052],
])
write_stru(Atoms(elem, cell=cell, scaled_positions=taud), 
           outdir=here, pp_file=pseudopotentials, orb_file=basissets, fname='STRU2')

write_input(
    data={
        'basis_type': 'lcao',
        'symmetry': 1,
        'kspacing': 0.25, # Oops!
        'init_chg': 'auto',
        'cal_force': 1,
        'pseudo_dir': pporb,
        'orbital_dir': pporb,
    },
    fn=here / 'INPUT'
)
