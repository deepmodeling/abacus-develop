# Basis Set and Pseudopotentials

## Basis Set

ABACUS supports both PW and LCAO basis set, controlled by keyword [basis_type](./input_files/input-main.md#basis_type) in INPUT file.

The default value of basis_type is pw. The size of pw basis set is controlled by imposing an upper bound for the [kinetic energy cutoff](./input_files/input-main.md#ecutwfc) of the plane wave.

When choosing lcao basis set, users need to prepare a set of atomic orbitals. Such files may be downloaded from the [official website](http://abacus.ustc.edu.cn/pseudo/list.htm). For more information, also check the `NUMERICAL_ORBITAL` section in the specification of the [STRU file](./input_files/stru.md).

The angular part of orbitals are real spherical harmonics defined (in terms of conventional spherical harmonics in quantum mechanical literature) as

$$Y_{lm} = \left\{\begin{matrix}\sqrt{2}~\text{Im} Y_l^{|m|} & m \lt 0 \\[6pt] Y_l^0 & m = 0 \\[6pt] \sqrt{2}~\text{Re}Y_{l}^{|m|} & m \gt 0\end{matrix}\right. $$

Note that real spherical harmonics adopted by ABACUS differ from some other definition, e.g. [Table of spherical harmonics - Wikipedia](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics), by a factor of $(-1)^m$.

Inside ABACUS, orbitals in LCAO basis are arranged lexicographically by species-position-l-zeta-m **except for the intra-m ordering**. Specifically, orbitals are first ordered by their atomic species in accordance with the `ATOMIC_SPECIES` section of the STRU file. For orbitals of the same species, orbitals belonging to each atom are put together, with their overall order following the `ATOMIC_POSITIONS` section of the STRU file. Orbitals on each atom are further ascendingly ordered by their angular momentum (s,p,d,f,...), followed by an order based on their their zeta number. Finally, m is ordered as 0, 1, -1, 2, 2, $\ldots$, l, -l, which is the only exception of the lexicographic order.


## Generating atomic orbital bases

Users may also generate ABACUS numerical atomic obitals based on their own flavor. The theoretical background of orbital generation can be found in following works:

- Spillage: [Chen M, Guo G C, He L. Systematically improvable optimized atomic basis sets for ab initio calculations[J]. Journal of Physics: Condensed Matter, 2010, 22(44): 445501.](https://iopscience.iop.org/article/10.1088/0953-8984/22/44/445501)
- PTG DPSI: [Lin P, Ren X, He L. Strategy for constructing compact numerical atomic orbital basis sets by incorporating the gradients of reference wavefunctions[J]. Physical Review B, 2021, 103(23): 235131.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.235131)

Guidelines for generating atomic orbital bases are as follows:

- [Numerical Atomic Orbitals 1: the nomenclature and usage of numerical atomic orbitals in ABACUS](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html) (Chinese)
- [Numerical Atomic Orbitals 2: generate numerical atomic orbitals based on given norm-conserving pseudopotential](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html) (Chinese)
- [Numerical Atomic Orbitals 3: generate high-precision numerical atomic orbitals](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html) (Chinese)

Stable orbital generation programs can be found in guidelines above, there is also another developing version of orbital generation program, in which algorithms are consecutively improved: [Github repository of ABACUS ORBGEN project](https://github.com/kirk0830/ABACUS-ORBGEN), the usage of which can be found in README (in English) file.

*NOTE*: users are encouraged to cite the above works when numerical atomic orbitals and its generation codes are used in their research.

## BSSE Correction

For treating BSSE(Basis Set Superposition Error), we allow for the inclusion of "empty" or "ghost" atoms in the calculation. Namely, when expanding the Hamiltonian, basis sets on the atoms are used, while the ionic potentials on those atoms are not included when constructing the Hamiltonian.

An empty atom is defined in the `STRU` file when an element name contains the "empty" suffix, such as "H_empty", "O_empty" and so on. Here we provide an [example](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/bsse/water) of calculating the molecular formation energy of $H_2O$ with BSSE correction.

In the example, we provide four STRU files:

- STRU_0 : used along with ntype = 2;normal calculation of water molecule ($E(\text{H}_2\text{O})$)

  obtained total energy of -466.4838149140513 eV
- STRU_1 : used along with ntype = 2;calculation of single O atom ($E_O$)

  obtained total energy of -427.9084406198214 eV
- STRU_2 : used along with ntype = 3;calculation of 1st H atom ($E_{H1}$)

  obtained total energy of -12.59853381731160 eV
- STRU_3 : used along with ntype = 3;calculation of 2nd H atom ($E_{H2}$)

  obtained total energy of -12.59853378720844 eV

> Note : Remember to adjust the parameter `ntype` in INPUT file

Thus, the formation energy is given by:

$$
\Delta E(\text{H}_2\text{O}) = E(\text{H}_2\text{O}) - E(\text{O}) - E(\text{H}^1) - E(\text{H}^2) \approx -13.38 eV
$$

## Pseudopotentials
### Supported formats
ABACUS supports both norm-conserving and ultrasoft pseudopotentials. For norm-conserving pseudopotentials, UPF, UPF2, VWR, and BLPS formats are supported. For ultrasoft pseudopotentials, UPF and UPF2 formats are supported. 

### Usage
For more information about pseudopotential usage, check the `ATOMIC_SPECIES` section in the specification of the [STRU file](./input_files/stru.md).

### Download
Users can find pseudopotentials in the following links:

**Website**
- [Quantum ESPRESSO](https://www.quantum-espresso.org/pseudopotentials): the official website of Quantum ESPRESSO, where you can find a large number of pseudopotential files.
- [Stantard Solid State Pseudopotential library](https://www.materialscloud.org/sssp): a library of **high-quality** pseudopotentials for solid-state calculations, with **a large number of tests on efficiency and precison**.
- [PWmat](http://www.pwmat.com/potential-download): a website that provides a large number of pseudopotential files, various kinds of semi-core constructed pseudopotentials are included. **Several sets (with or without f-electrons/noncolinear core correction) of Lanthanide pseudopotentials are also available**.
- [THEOS](http://theossrv1.epfl.ch/Main/Pseudopotentials): PSlibrary 0.3.1, a library of pseudopotentials for DFT calculations, including ultrasoft, paw, norm-conserving both full-relativistic and scalar-relativistic pseudopotentials.
- [ABACUS@USTC](https://abacus.ustc.edu.cn/pseudo/list.htm): **ABACUS official website** where you can find a large number of pseudopotential files and numerical atomic orbital files.
- [BLPS](https://github.com/PrincetonUniversity/BLPSLibrary): BLPS format pseudopotential library

**Norm-conserving pseudopotentials**
- [SG15](http://www.quantum-simulation.org/potentials/sg15_oncv/): **vastly used in ABACUS** DFT calculation and numerical atomic orbital generation.
- [PseudoDOJO](http://www.pseudo-dojo.org/): another widely used pseudopotential database, developed by Abinit group, **including Lanthanide pseudopotentials (f-electrons frozen)**.
- [The Rappe group](https://www.sas.upenn.edu/rappegroup/research/pseudo-potential-gga.html): a collection of GGA pseudopotentials which are generated with Opium code, several tests proves that are out-performing in alloy systems.
- [Matteo Giantomassi's Github repo](https://github.com/gmatteo/pseudos_ac_she): a Github repository that contains norm-conserving pseudopotentials for **Actinides and superheavy elements to 120-th element**.

**Ultrasoft pseudopotentials**
- [Vanderbilt](http://www.physics.rutgers.edu/~dhv/uspp/): a collection of ultrasoft pseudopotentials generated by Vanderbilt group.
- [GBRV](https://www.physics.rutgers.edu/gbrv/) by Kevin F. Garrity, Joseph W. Bennett, Karin M. Rabe, and David Vanderbilt: presently the most popular ultrasoft pseudpotentials in Quantum ESPRESSO user community.

### Pseudopotential Generation
For pseudopotential generation, please refer to the following links for more information:
- [Quantum ESPRESSO](http://www.quantum-espresso.org/pseudopotentials/)
- [ONCVPSP](http://www.mat-simresearch.com/)
- [Opium](https://opium.sourceforge.net/)

A Chinese guideline is also available here: [A brief introduction of norm-conserving pseudopotential generation](https://mcresearch.github.io/abacus-user-guide/abacus-upf.html)

# ABACUS Pseudopotential-Numerical atomic orbital Square (APNS) project
For the purpose of providing high-quality pseudopotentials and numerical atomic orbitals, we have initiated the APNS project. The project is aimed at providing a large number of high-quality pseudopotentials and numerical atomic orbitals, along with diverse test data for the ABACUS user community, reduce the cost of generating and testing pseudopotentials and numerical atomic orbitals by users, and promote the development of ABACUS software. The project is currently in the development stage, and we welcome contributions from the community. For more information, please refer to the following links:
- [APNS website: test data and results](https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square/)
- [APNS workflow (Github repository): high-throughput test of pseudopotentials and numerical atomic orbitals](https://github.com/kirk0830/ABACUS-Pseudopot-Nao-Square)

There are also other excellent projects that provide high-quality pseudopotentials along with test data:
- [Solid State Pseudopotential library](https://www.materialscloud.org/sssp)
- [Verification of the precision of DFT implementation via AiiDA common workflows](https://acwf-verification.materialscloud.org/)
