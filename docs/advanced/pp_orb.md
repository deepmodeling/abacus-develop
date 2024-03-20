# Basis Set and Pseudopotentials

## Basis Set

ABACUS supports both PW and LCAO basis set, controlled by keyword [basis_type](./input_files/input-main.md#basis_type) in INPUT file.

The default value of basis_type is pw. The size of pw basis set is controlled by imposing an upper bound for the [kinetic energy cutoff](./input_files/input-main.md#ecutwfc) of the plane wave.

When choosing lcao basis set, users need to prepare a set of atomic orbitals. Such files may be downloaded from the [official website](http://abacus.ustc.edu.cn/pseudo/list.htm). For more information, also check the `NUMERICAL_ORBITAL` section in the specification of the [STRU file](./input_files/stru.md).

The sequence of orbitals in lcao basis set is as follows. First, all the orbitals belonging to one particular atom are put together. These atom orbitals are arranged as the atom order specified in the STRU file. Then, the orbitals of each atom are arranged according to the orbital files. If the orbital file says that the number of s、p、d...orbitals is $n_s$、$n_p$、$n_d$...then the orbitals are aligned as first $n_s$ s orbitals, then $n_p$ p orbitals, and then $n_d$ d orbitals...Last, the angular part of each orbital is real spherical harmonic function. They are aligned as Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2, which is s,$p_z$,$p_x$,$p_y$,$d_{z^2}$,$d_{xz}$,$d_{yz}$,$d_{x^2-y^2}$,$d_{xy}$. The corresponding formula can be seen in [Table of spherical harmonics - Wikipedia](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics). Note that these formula lack of the Condon–Shortley phase $(-1)^m$, which is presented in the lcao orbitals of ABACUS.

## Generating atomic orbital bases

Users may also choose to generate their own atomic obitals. In ABACUS, the atomic orbital bases are generated using a scheme developed in the [paper](https://iopscience.iop.org/article/10.1088/0953-8984/22/44/445501). A detailed description of the procedure for generating orbitals will be provided later.

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

In ABACUS, we support norm-conserving and ultrasoft pseudopotentials. 
For norm-conserving pseudopotentials, we support four different formats of the pseudopotential files: UPF, UPF2, VWR, and BLPS. 
For ultrasoft pseudopotentials, currently we support only one format of the pseudopotential files: UPF2.

For more information, check the `ATOMIC_SPECIES` section in the specification of the [STRU file](./input_files/stru.md).

Here we list some common sources of the pseudopotential files:

1. [Quantum ESPRESSO](http://www.quantum-espresso.org/pseudopotentials/).
2. [SG15-ONCV](http://quantum-simulation.org/potentials/sg15_oncv/upf/).
3. [DOJO](http://www.pseudo-dojo.org/).
4. [BLPS](https://github.com/PrincetonUniversity/BLPSLibrary).
