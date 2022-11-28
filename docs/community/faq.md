# Frequently Asked Questions

 - [General Questions](#general-questions)
 - [Installation](#installation)
 - [Setting up jobs](#setting-up-jobs)
 - [Failed jobs](#failed-jobs)
 - [Miscellaneous](#miscellaneous)


## General Questions

**1. What are the merits of ABACUS in respect of functionality, performance, and/or accuracy?**

Users are referred to the introduction of features of ABACUS in the [feature list](http://abacus.ustc.edu.cn/features/list.htm).

## Installation


## Setting up jobs

**1. Why pseudopotential files must be provided in LCAO calculation?**

The pseudopotentials are used to approximate the potential of nuclear and core electrons, while the numerical orbitals are basis sets used to expand the Hamiltonian. So both pseudopotential and numerical orbital files are needed in LCAO calculation.

**2. What is the correlation between pseudopotential and numerical orbital files?**

The numerical orbital files are generated for a specific pseudopotential. So the right numerical orbital files must be chosen for a specific pseudopotential. We suggest users choose numerical orbital files and the corresponding pseudopotential files from the [ABACUS website](http://abacus.ustc.edu.cn/pseudo/list.htm) because their accuracy has been tested. However, interesting users may also generate their own numerical orbitals for a specific type of pseudopential by using the tools provided in  the [abacus-develop/tools](https://github.com/deepmodeling/abacus-develop/tree/develop/tools) directory.


**3. How to set `ecutwfc` in LCAO calculation? Must it be 100 Ry for a numerical orbital file like `Cu_lda_7.0au_100Ry_2s2p2d`?**

It is recommended to set `ecutwfc` to the value that the numerical orbital file suggests, but it is not a must. The `ecutwfc` value only affects the number of FFT grids.

**4. Does ABACUS support LCAO calculations accounting for external electric field effects?**

Yes, users are referred to documentation on [external electric field](../advanced/scf/advanced.md#external-electric-field).


**5. Can ABACUS calculate non-periodic systems, such as ionic liquids?**

Non-periodic systems such as liquid systems can be calculated by using supercell and gamma-only calculation.

**6. How to perform spin-orbital coupling (SOC) calculations in ABACUS?**

Apart from setting relavant keys (`lspinorb` to 1) in the `INPUT` file, SOC calculations can only be performed with fully-relativistic pseudopotentials. Users are suggested to download fully-relativistic versions of SG15_ONCV pseudopotential files from a [website](http://quantum-simulation.org/potentials/sg15_oncv/upf/). The numerical orbital files generated from the corresponding scalar-relativistic pseudoptential files by ABACUS ([here](http://abacus.ustc.edu.cn/pseudo/list.htm)) can be used in collaboration with the fully-relativistic pseudopotentials.


**7. How to restart jobs in abacus?**

For restarting SCF calculations, users are referred to the documentation about [continuation of job](../advanced/scf/spin.md#for-the-continuation-job). For restarting MD calculations, please see [md_restart](../advanced/input_files/input-main.md#md_restart).


[back to top](#frequently-asked-questions)

## Failed jobs

**1. Why my calculation is pend by using mpirun?**

This is usually caused by overloading of CPUs' memory without specifying thread numbers. Setting `export OMP_NUM_TREANDS=1` will solve this problem, or one can use the command like `OMP_NUM_THREADS=1 mpirun -n 8 abacus` to rerun failed jobs.


**2. My relaxation failed. How to deal with it?**

This is usually caused by the difficulty in converging charge density. Reducing charge mixing coefficient (`mixing_beta`) might help. For large systems over 1000 atoms, it is suggested to choose the `pulay-kerker` method as the `mixing_type`.




## Miscellaneous

[back to top](#frequently-asked-questions)
