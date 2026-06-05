# ASE

## Introduction

[ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment) performs as a powerful Pythonic platform for atomistic simulations, in which there are plenty of functionalties supported, such as various geometry optimization algorithms for finding both the minimum energy point and the transition states, including BFGS, BFGSLineSearch, FIRE, NEB, AUTO-NEB, etc, and various molecular dynamics techniques, including thermostats (Langevin, CSVR, Nose-Hoover Chain, etc) and metadynamics (via the interface with Plumed). 

Due to the growing number of softwares and machine-learning forcefields, we turn to maintain the interface with ASE by our own, while a legacy version of ASE interface can still be found at [our GitLab repository of ase-abacus](https://gitlab.com/1041176461/ase-abacus ).

## Installation

We strongly recommend you create a virtual environment for the installation of Python packages of abacus, such as `conda` or `venv`, to avoid conflicts with other packages, for example, with the `conda`:

```bash
conda create -n abacus python=3.10
conda activate abacus
```

Then, install the ASE interface by:

```bash
cd interfaces/ASE_interface
pip install .
```

## ABACUS Calculator

Present calculator implementation requires a "profile" to act as an interface between the Python runtime and the file system.
Instantiate an `AbacusProfile` object with proper settings:

```python
from abacuslite import AbacusProfile
aprof = AbacusProfile(
    command='mpirun -np 4 abacus',
    omp_num_threads=1,
    pseudo_dir='/path/to/folder/of/pseudopotentials',
    orbital_dir='/path/to/folder/of/orbitals', # OPTIONAL!
)
```
, by such lines, you build the interface between the computational environment and the Python runtime.
This interface can be reused in multiple calculations.

Then, you can instantiate the `Abacus` calculator with the profile by:

```python
from abacuslite import Abacus
abacus = Abacus(
    profile=aprof,
    directory='/path/to/work/directory',
    pseudopotentials={
        'Si': 'Si_ONCV_PBE-1.0.upf',
    },
    basissets={
        'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
    },
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
```
, where except the `directory`, you can focus on the setting of ABACUS itself. In `inp`, you can set everything as you do in INPUT file of ABACUS. The kpoint sampling can also be set by the `kpts` parameter, like:

```python
abacus = Abacus(
    # all other parameters
    kpts={
        'mode': 'mp-sampling',
        'gamma-centered': True,
        'nk': (4, 4, 4),
        'kshift': (0, 0, 0),
    }
)
```

If with the `tempfile` module, you can create an abacus instance whose directory will be automatically removed when leaves from the context:

```python
import tempfile
with tempfile.TemporaryDirectory() as tmpdir:
    abacus = Abacus(
        profile=aprof,
        directory=tmpdir,
        pseudopotentials={
            'Si': 'Si_ONCV_PBE-1.0.upf',
        },
        basissets={
            'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
        },
        inp={
            # the rest of input parameters
        }
    )
```

## Perform Calculations

In the new implementation, we limit the range of functionalties supported to mainly include the necessary ones, such as the SCF calculation, the energy and force/stress evaluation. The other features, such as starting the molecule dynamics directly in ABACUS from Python, is not supported anymore. Instead, it is encouraged to use the ASE tools to perform the molecule dynamics.

Please read the examples in `interfaces/ASE_interface/examples/` for more details.

## SPAP Analysis

[SPAP](https://github.com/chuanxun/StructurePrototypeAnalysisPackage) (Structure Prototype Analysis Package) is written by Dr. Chuanxun Su to analyze symmetry and compare similarity of large amount of atomic structures. The coordination characterization function (CCF) is used to 
measure structural similarity. An unique and advanced clustering method is developed to automatically classify structures into groups. 


If you use this program and method in your research, please read and cite the publication:

`Su C, Lv J, Li Q, Wang H, Zhang L, Wang Y, Ma Y. Construction of crystal structure prototype database: methods and applications. J Phys Condens Matter. 2017 Apr 26;29(16):165901.`

and you should install it first with command `pip install spap`.
