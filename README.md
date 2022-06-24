<p align="center">
    <img src="docs/abacus-logo.jpg">
</p>

<p align="center">
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml/badge.svg">
    </a>
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml/badge.svg">
    </a>
</p>

<a id="readme-top"></a>

ABACUS is an electronic structure package based on density functional theory(DFT), adopting either plane wave basis or numerical atomic orbitals.
Please refer to our [GitHub repository](https://github.com/deepmodeling/abacus-develop) for more information and support.
# Table of contents

- [Table of contents](#table-of-contents)
- [Features](#features)
- [Download and install](#download-and-install)
- [Quickstart guide](#quickstart-guide)
  - [Input files](#input-files)
  - [Run ABACUS](#run-abacus)
  - [Output files](#output-files)
- [Features](#features-1)
- [Functionalities](#functionalities)
- [Examples](#examples)
- [For developers](#for-developers)

# Features

ABACUS provides the following features and functionalities:

1. Ground-state total energy calculations using Kohn-Sham (KS) density functional theory
(DFT) with local-density, generalized gradient approximations (LDA/GGAs), and hybrid functionals
(PBE0 and HSE06, only for LCAO).
2. Brillouin zone sampling using the Monkhorst-Pack special k-points.
3. Geometry relaxations with Conjugated Gradient (CG) and BFGS methods.
4. Semi-empirical van der Waals energy correction using the Grimme DFT-D2/D3 scheme.
5. NVT molecular dynamics simulation.
6. Stress calculations and cell relaxations.
7. Electric polarization calculation using Berry Phase theory.
8. Interface to the Wannier90 package.
9. Real-time time dependent density functional theory (TDDFT).
10. Electrostatic potential.
11. Mulliken charge analysis.
12. Projected density of states (PDOS).

[back to top](#readme-top)

# Download and install

ABACUS can be downloaded from our [official website](http://abacus.ustc.edu.cn/) or [GitHub release page](https://github.com/deepmodeling/abacus-develop/releases) for stable versions. You can also get the developing version from our [GitHub repository](https://github.com/deepmodeling/abacus-develop).

Please refer to the [installation guide](docs/install.md) for instruction on the structure of the package and how to install ABACUS.

[back to top](#readme-top)

# Quickstart guide

## Input files

The following files are the central input files for ABACUS. Before executing the program, please make sure these files are prepared and stored in the working directory.

- The INPUT file

    The file named INPUT contains the setting parameters used in the calculation, which informs the program “what to do and how to do it”. Most parameters are supplied with default values, but some important parameters must be explicitly set by the user. For a complete list of the input parameters, please consult this [instruction](docs/input-main.md).

    *Attention: Users cannot change the filename “INPUT” to other names.*
- The structure file

    The default name for structure file is STRU.The name can however be changed to a different name by explicitly specifying the name in the INPUT file.

    The STRU file contains the structural information about the system, e.g., lattice constant, lattice vectors, and positions of the atoms within a unit cell. The positions can be given either in direct or Cartesian coordinates. Moreover, the name (and location of the pseudopotential and numerical orbital files, see below) need to be specified in the STRU file.

    Specifications of the STRU file can be found in this [short instruction](docs/input-stru.md).
- The k-point file

    The default name is KPT. It contains the information of the k-grid setting for the Brillouin zone sampling.

    Specification of the k-point file can be found in this [short instruction](docs/input-kpt.md).
- The pseudopotential files

    Norm-conserving pseudopotentials are used in ABACUS, in the UPF file format.The filename of each element’s pseudopotential needs to be specified in the `STRU` file, if the the pseudopotential files are already present in the working directory. However, in case that the pseudopotential files are stored in some other directories, then a full path to access the pseudopotential files have to be specified in the `STRU` file.

    More information on pseudopotentials is given [here](docs/features.md#pseudopotentials).

- The numerical orbital files

    When performing calculations with numerical atomic orbital basis, it is necessary to prepare a numerical orbital file for each element in the system. Generally, the numerical orbital file should be prepared by the user, which will be described later. The filename for each element’s numerical orbital basis needs to be specified in the `STRU` file. However, in case that the numerical orbital files are stored in a location different from the working directory, then a full path to access the orbital files have to be specified in the `STRU` file.
    ABACUS provides numerical atomic basis sets of different accuracy levels for most elements commonly used. Users can download these basis sets from the [website](http://abacus.ustc.edu.cn/pseudo.html). Moreover, users can generate numerical atomic orbitals by themselves, and the procedure is provided in this [short introduction](docs/generate-basis.md).

[back to top](#readme-top)

## Run ABACUS

After putting all required input files under one folder, enter this folder.
```bash
cd input_folder
```

Perform calculation by:
```bash
mpirun -np 4 abacus
```

You can replace `4` with your desired number of process, typically the result of the command `nproc`.

[back to top](#readme-top)

## Output files

When the calculation finishes, the program will create an output directory (default: OUT.ABACUS/),
into which the following output files will be generated:

1. `INPUT`: contains all input parameters, user’s input and default.
2. `istate.info`: information of energy eigenvalues.
3. `running_${calculation}.log`: contains the running details. Information on the variable calculation is found in the [list of keywords](docs/input-main.md#calculation). For example, if we are doing a SCF calculation, the log files will be named running_scf.log.
4. `STRU_READIN_ADJUST.cif`: structure file in the cif formatter.
5. `warning.log`: errors and warning messages.
6. directories containing element information. For example, Si/:
    - `Si.NONLOCAL`: non-local pseudopotential projectors.
    - `Si-P.ORBITAL`: pseudo atomic orbitals, p orbital
    - `Si-S.ORBITAL`: pseudo atomic orbitals, s orbital
    - `v_loc_g.dat`: vlocal in G space

[back to top](#readme-top)

# Features

Users can refer to this [page](docs/features.md) for several features of the ABACUS code:

- [Basis sets](docs/features.md#basis-sets)
- [Pseudopotentials](docs/features.md#pseudopotentials)
- [Boundary conditions and k-points](docs/features.md#boundary-conditions-and-k-points)
- [Kohn-Sham solver](docs/features.md#kohn-sham-solver)
- [Exchange-correlation functionals](docs/features.md#exchange-correlation-functionals)

[back to top](#readme-top)

# Functionalities

ABACUS provides a wide variety of functionalities, with explanation and examples:

- [Basic electronic structure calculation with PW basis set](docs/examples/basic-pw.md)
- [Basic electronic structure calculation with LCAO basis set](docs/examples/basic-lcao.md)
- [DFT + dispersion calculations](docs/examples/dispersion.md)
- [DOS, wave functions](docs/examples/dos.md)
- [Band structure](docs/examples/band-struc.md)
- [Magnetic properties](docs/examples/magnetic.md)
- [Force calculation and structure relaxation](docs/examples/force.md)
- [Stress calculation and cell relaxation](docs/examples/stress.md)
- [Molecular dynamics](docs/examples/md.md)
- [Macroscopic polarization calculation](docs/examples/berry-phase.md)
- [ABACUS-wannier90 interface](docs/examples/wannier90.md)
- [Real-time time dependent density functional theory](docs/examples/tddft.md)
- [Electrostatic potential](docs/examples/potential.md)
- [Mulliken charge](docs/examples/mulliken.md)
- [Hybrid functional](docs/examples/hybrid.md)

[back to top](#readme-top)

# Examples

We also provide many examples in the directories examples/ and tests/.

Note that the examples there are intended as references, and the results are not converged with regard to basis set or k point sampling.

In the directory tests/, each sub-directory contains a separate test example. An introduction of the examples in tests/ directory can be found [here](tests/README.md). In each subdirectory, you may also find a file named jd which contains a short job description, and for some cases you may also find a README file containing more details about the run. Also, reference output is provided in the file `result.ref`.

[back to top](#readme-top)

# For developers

We also provide some [information](docs/CONTRIBUTING.md) on how to make contributions to ABACUS.

- [Structure of the package](docs/CONTRIBUTING.md#structure-of-the-package)
- [Submitting an Issue](docs/CONTRIBUTING.md#submitting-an-issue)
- [Comment Style for documentation](docs/CONTRIBUTING.md#comment-style-for-documentation)
- [Code formatting style](docs/CONTRIBUTING.md#code-formatting-style)
- [Submitting a Pull Request](docs/CONTRIBUTING.md#submitting-a-pull-request)
- [Commit Message Guidelines](docs/CONTRIBUTING.md#commit-message-guidelines)
[back to top](#readme-top)
