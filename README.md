<div style="text-align:center; margin: 20px">
    <img src="documents/abacus-logo.jpg">
</div>

<div style="text-align:center; margin: 20px">
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/container.yml/badge.svg">
    </a>
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml/badge.svg">
    </a>
</div>

<a id="readme-top"></a>
WELCOME TO THE "ABACUS" PROGRAM!
THE PROJECT STARTS FROM https://github.com/abacusmodeling/abacus-develop,
WHERE MORE INFORMATION CAN BE FOUND.

# Table of contents
- [Table of contents](#table-of-contents)
- [About ABACUS](#about-abacus)
- [Download and install](#download-and-install)
- [Quickstart guide](#quickstart-guide)
  - [Input files](#input-files)
  - [Output files](#output-files)
- [Features](#features)
- [Functionalities](#functionalities)
- [Examples](#examples)
- [For developers](#for-developers)
  - [Build and install ABACUS from CMake](#build-and-install-abacus-from-cmake)


# About ABACUS
ABACUS IS AN ELECTRONIC STRUCTURE PACKAGE BASED ON DENSITY FUNCTIONAL THEORY.
ABACUS ADOPTS EITHER PLANE WAVE BASIS OR NUMERICAL ATOMIC ORBITALS

---

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
ABACUS can be downloaded from its [official website](http://abacus.ustc.edu.cn/) or our [github website](https://github.com/deepmodeling/abacus-develop.git).

Please refer to the [installation guide](doc/install.md) for instruction on the structure of the package and how to install ABACUS.

[back to top](#readme-top)

# Quickstart guide

## Input files
The following files are the central input files for ABACUS. Before executing the program, please
make sure these files are prepared and stored in the working directory.
- The INPUT file

    The file named INPUT contains the setting parameters used in the calculation, which informs the program “what to do and how to do it”. Most parameters are supplied with default values, but some important parameters must be explicitly set by the user. For a complete list of the input parameters, please consult this [instruction](doc/input-main.md).

    *Attention: Users cannot change the filename “INPUT” to other names.*
- The structure file

    The default name for structure file is STRU.The name can however be changed to a different name by explicitly specifying the name in the INPUT file.

    The STRU file contains the structural information about the system, e.g., lattice constant, lattice vectors, and positions of the atoms within a unit cell. The positions can be given either in direct or Cartesian coordinates. Moreover, the name (and location of the pseudopotential and numerical orbital files, see below) need to be specified in the STRU file.

    Specifications of the STRU file can be found in this [short instruction](doc/input-stru).
- The k-point file
    The default name is KPT. It contains the information of the k-grid setting for the Brillouin zone sampling.

    Specification of the k-point file can be found in this [short instruction](doc/input-kpt).
- The pseudopotential files

    Norm-conserving pseudopotentials are used in ABACUS, in the UPF file format.The filename of each element’s pseudopotential needs to be specified in the `STRU` file, if the the pseudopotential files are already present in the working directory. However, in case that the pseudopotential files are stored in some other directories, then a full path to access the pseudopotential files have to be specified in the `STRU` file.

    More information on pseudopotentials is given [below](#pseudopotentials).

- The numerical orbital file

    When doing calculations with atomic orbital basis, it’s necessary to prepare a numerical orbital file for each element in the system. Generally, the numerical orbital file should be prepared by the user, which will be described later. The filename for each element’s numerical orbital basis needs to be specified in the `STRU` file. However, in case that the numerical orbital files are stored in a location different from the working ddirectory, then a full path to access the orbital files have to be specified in the `STRU` file.
    ABACUS provides atomic basis sets of different accuracy levels for most elements commonly used. Users can download these basis sets from our [website](http://abacus.ustc.edu.cn/pseudo.html). Moreover, users can generate basis themselves, and the procedure is provided in this [short introduction](doc/generate-basis.md).

[back to top](#readme-top)

## Output files

When the calculation finishes, the program will create an output directory (default: OUT.ABACUS/),
into which the following output files will be generated:
1. INPUT: contains all input parameters, user’s input and default.
2. istate.info: information of energy eigenvalues.
3. running_scf.log: contains the running details.
4. STRU_READIN_ADJUST.cif: structure file in the cif formatter.
5. warning.log: errors and warning messages.
6. directories containing element information. For example, Si/:
    - Si.NONLOCAL: non-local pseudopotential projectors.
    - Si-P.ORBITAL: pseudo atomic orbitals, p orbital
    - Si-S.ORBITAL: pseudo atomic orbitals, s orbital
    - v_loc_g.dat: vlocal in G space

[back to top](#readme-top)

# Features

Users can refer to this [page](doc/features.md) for several features of the ABACUS code:

   - [Basis sets](doc/features.md#basis-sets)
   - [Pseudopotentials](doc/features.md#pseudopotentials)
   - [Boundary conditions and k-points](doc/features.md#boundary-conditions-and-k-points)
   - [Kohn-Sham solver](doc/features.md#kohn-sham-solver)
   - [Exchange-correlation functionals](doc/features.md#exchange-correlation-functionals)

[back to top](#readme-top)

# Functionalities
ABACUS provides a wide variety of functionalities, with explanation and examples:

- [Basic electronic structure calculation with PW basis set](doc/examples/basic-pw.md)
- [Basic electronic structure calculation with LCAO basis set](doc/examples/basic-lcao.md)
- [DFT + dispersion calculations](doc/examples/dispersion.md)
- [DOS, wave functions](doc/examples/dos.md)
- [Band structure](doc/examples/band-struc.md)
- [Magnetic properties](doc/examples/magnetic.md)
- [Force calculation and structure relaxation](doc/examples/force.md)
- [Stress calculation and cell relaxation](doc/examples/stress.md)
- [Molecular dynamics](doc/examples/md.md)
- [Macroscopic polarization calculation](doc/examples/berry-phase.md)
- [ABACUS-wannier90 interface](doc/examples/wannier90.md)
- [Real-time time dependent density functional theory](doc/examples/tddft.md)
- [Electrostatic potential](doc/examples/potential.md)
- [Mulliken charge](doc/examples/mulliken.md)

[back to top](#readme-top)

# Examples
We also provide many examples in the directories examples/ and tests/

# For developers

We also provide some [information](doc/developer.md) for developers.

---

## Build and install ABACUS from CMake

Check the cmake version on your machine
```bash
cmake --version
```
ABACUS requires the minimum cmake version `3.18`.

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`.
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH}
```
You can provide root path of each dependent package if the package cannot be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDEDIR`, `BOOST_INCLUDEDIR`, `MPI_CXX_COMPILER` and `MKL_DIR`. are currently available to specify.
For example
```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3
```

If the cmake has executed successfully, then
```bash
cmake --build build
cmake --install build
```
If no install prefix is specified, the binary will be installed to `/usr/local/bin/ABACUS` by default.
