'''
Example modified from:

'''
import tempfile
from ase import units
from ase.io import read
from ase.md.langevin import Langevin
from ase.constraints import FixedPlane
from ase.calculators.plumed import Plumed
from abacuslite import AbacusProfile, Abacus

timestep = 0.005
ps = 1000 * units.fs

setup = [f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",
         "COM ATOMS=1-7 LABEL=com",
         "DISTANCE ATOMS=1,com LABEL=d1",
         "UPPER_WALLS ARG=d1 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=2,com LABEL=d2",
         "UPPER_WALLS ARG=d2 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=3,com LABEL=d3",
         "UPPER_WALLS ARG=d3 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=4,com LABEL=d4",
         "UPPER_WALLS ARG=d4 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=5,com LABEL=d5",
         "UPPER_WALLS ARG=d5 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=6,com LABEL=d6",
         "UPPER_WALLS ARG=d6 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=7,com LABEL=d7",
         "UPPER_WALLS ARG=d7 AT=2.0 KAPPA=100.",
         "c1: COORDINATIONNUMBER SPECIES=1-7 MOMENTS=2-3" +
         " SWITCH={RATIONAL R_0=1.5 NN=8 MM=16}",
         "METAD ARG=c1.* HEIGHT=0.05 PACE=500 " +
         "SIGMA=0.1,0.1 GRID_MIN=-1.5,-1.5 GRID_MAX=2.5,2.5" +
         " GRID_BIN=500,500 BIASFACTOR=5 FILE=HILLS"]

# we write the following and read it back to an ASE.atoms.Atoms
# object
atoms = '''
7  
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

cons = [FixedPlane(i, [0, 0, 1]) for i in range(7)]
atoms.set_constraint(cons)
atoms.set_masses([1, 1, 1, 1, 1, 1, 1])



atoms.calc = Plumed(calc=LennardJones(rc=2.5, r0=3.0),
                    input=setup,
                    timestep=timestep,
                    atoms=atoms,
                    kT=0.1)

dyn = Langevin(atoms, timestep, temperature_K=0.1/units.kB, friction=1,
               fixcm=False, trajectory='MTD.traj')

dyn.run(500000)
