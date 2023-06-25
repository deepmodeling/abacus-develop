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

**8. Can DeePKS model be used for structural optimization calculation? What parameters need to be modified or called?**

If you train the DeePKS model with force labels, then the DeePKS model can provide force calculation with the same accuracy as your target method, and can thus be used for structural optimization. To do that, you just need to train the model with force label enabled.

**9. How to estimate the max memory consumption?**

Run `/usr/bin/time -v mpirun -n 4 abacus`, and locate "Maximum resident set size" in the output log at the end. Please note that this value is the peak memory size of the main MPI process.

**10. Why there are two sigma (smearing_sigma and dos_sigma) in some examples for dos calculation?**

 The tag `smearing_sigma` is used for SCF calculation, and does not affect NSCF calculation. The tag `dos_smearing` is only used for plotting density of states, which does affect SCF or NSCF results. So `smearing_sigma` should not be set in dos calculation.

**11. How to set `nbands` and `ncpus`?** 

For both pw and LCAO calculations, the default value for `nbands` can be found [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#nbands). Note that the number of CPUs called for a parallel job (i.e., the number specified after the command `mpirun -n`) should be smaller than `nbands`, otherwise the job will fail with an error message `nbands < ncpus`. Note also that for LCAO calculations, `nbands` should always be smaller than `nlocal`, i.e., the number of the atomic orbital basis of the system. 

[back to top](#frequently-asked-questions)

## Failed jobs

**1. Why my calculation is pend by using mpirun?**

This is usually caused by overloading of CPUs' memory without specifying thread numbers. ABACUS detects available hardware threads, and provides these information at the beginning of the program if used threads mismatches with hardware availability. User should keep total used threads(i.e. the number of OpenMP threads times the number of MPI processes) no more than the count of available CPUs(can be determined by `lscpu`). Setting `export OMP_NUM_THREADS=1` will solve this problem, or one can use the command like `OMP_NUM_THREADS=1 mpirun -n 8 abacus` to rerun failed jobs.

**2. My relaxation failed. How to deal with it?**

This is usually caused by the difficulty in converging charge density. Reducing charge mixing coefficient (`mixing_beta`) might help. For large systems up to hundreds of atoms, it is suggested to choose the Kerker mixing method by setting parameter "mixing_gg0" as "1.5".

Sometimes, loose convergence threshold of charge density (parameter "scf_thr") will cause atomic forces not correctly enough, please set it at most "1e-7" for relaxation calculation.

**3. Why the program is halted?**

If the program prompt something like "KILLED BY SIGNAL: 9 (Killed)", it may be caused by insufficient memory. You can use `dmesg` to print out system info regarding memory management, and check if there is "Out of memory: Killed" at the end of output info. Please try using less processes and threads for calculation, or modify the input parameters requiring less memory.

If the error message is "Segmentation fault", or there is no enough information on the error, please feel free to submit an issue.

## Miscellaneous

**1. How to visualize charge density file?**

The output file SPIN1_CHG.cube can be visualized by using VESTA.

**2. How to change cif file directly to STRU file?**

One way to change from cif to STRU is to use the [ASE-ABACUS](https://gitlab.com/1041176461/ase-abacus) interface. An example of the converting script is provided below:

```python
from ase.io import read, write
from pathlib import Path

cs_dir = './'
cs_cif = Path(cs_dir, 'SiO.cif')
cs_atoms = read(cs_cif, format='cif')
cs_stru = Path(cs_dir, 'STRU')
pp = {'Si':'Si.upf','O':'O.upf'}
basis = {'Si':'Si.orb','O':'O.orb'}
write(cs_stru, cs_atoms, format='abacus', pp=pp, basis=basis)
```

**3. What is the convergence criterion for the SCF process in ABACUS?**

ABACUS applies the density difference between two SCF steps (labeled as `DRHO` in the screen output) as the convergence criterion, which is considered as a more robust choice compared with the energy difference. `DRHO` is calculated via `DRHO = |rho(G)-rho_previous(G)|^2`. Note that the energy difference between two SCF steps (labed as `EDIFF`) is also printed out in the screen output.

[back to top](#frequently-asked-questions)
