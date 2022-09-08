# Features

Here we summarize several features of ABACUS:

- [Basis sets](#basis-sets)
- [Pseudopotentials](#pseudopotentials)
- [Boundary conditions and k-points](#boundary-conditions-and-k-points)
- [Kohn-Sham solver](#kohn-sham-solver)
- [Exchange-correlation functionals](#exchange-correlation-functionals)

     [back to main page](../README.md)

## Basis sets

In ABACUS, there are 3 types of basis set choices:

1. PW

    Plane wave basis.

2. LCAO

    Localized atomic-orbital basis; these systematically improvable atomic basis are generated with a tool called [SIAB](generate-basis.md).

3. LCAO_in_PW

    Expand the atomic basis in terms of plane waves, and use plane-waves techniques to construct the Hamiltonian matrix, but solve the eigenvalue problem within the LCAO basis set.

In the INPUT file, the keyword basis_type controls what basis type to use:

   1. PW: basis_type = pw
   2. LCAO: basis_type = lcao
   3. LCAO_in_PW: basis_type = lcao_in_pw

The default value of basis_type is pw. When choosing lcao or lcao_in_pw, the user should prepare atomic orbitals first.

Information on the keyword basis_type can also be found in [the list of input keywords](input-main.md#basis_type).

We also provide a [detailed instruction](generate-basis.md) on generating basis functions.

[back to top](#features)

## Pseudopotentials

Currently ABACUS uses norm-conserving pseudopotentials in the (old) format of UPF, which
is the standard pseudopotential format for Quantum ESPRESSO. To run a calculation, the pseudopotential needs to be set in STRU file. For example:

```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf
```

You can download the pseudopotential files from our [website](http://abacus.ustc.edu.cn/pseudo/list.htm).

There are pseudopotential files in these websites which are also supported by ABACUS:
1. [Quantum ESPRESSO](http://www.quantum-espresso.org/pseudopotentials/).
2. [SG15-ONCV](http://quantum-simulation.org/potentials/sg15_oncv/upf/).
3. [DOJO](http://www.pseudo-dojo.org/).

If LCAO base is used, the numerical orbital files should match the pseudopotential files. The [official orbitals package](http://abacus.ustc.edu.cn/pseudo/list.htm) only matches SG15-ONCV pseudopotentials.

[back to top](#features)

## Boundary conditions and k-points

ABACUS uses periodic boundary conditions for both crystals and finite systems. For isolated systems, such as atoms, molecules, clusters, etc., one uses the so-called supercell model. Lattice
vectors of the supercell are set in the STRU file.

For the input k-point (KPT) file, the file should either contain the k-point coordinates and weights or the mesh size for creating the k-point gird. Both options are allowed in ABACUS.

More information on k-points is provided in this [instruction](input-kpt.md)

[back to top](#features)

## Kohn-Sham solver

For different types of basis set choice, different methods are used to solve the Kohn-Sham
equation. For PW basis, there are CG and Blocked Davidson methods for solving the eigenvalue problem. For LCAO basis/LCAO_in_PW basis, one uses direct diagnolization method. In the INPUT file, the parameter ‘ks_solver’ controls what method to use for solveing the Kohn-Sham
equation for each basis.

- PW: ks_solver = ‘cg’ or ‘dav’
- LCAO: ks_solver = ‘hpseps’ , ‘genelpa’ , ‘scalapack_gvx’ or 'cusolver'
- LCAO_in_PW: ks_solver = ‘lapack’

If you set ks_solver=‘hpseps’ for basis_type=‘pw’, the program will be stopped with an error
message:

```
hpseps can not be used with plane wave basis.
```

Then the user has to correct the input file and restart the calculation.

Information on the keyword ks_solver is also given in the [list of input variables](input-main.md#ks_solver).

[back to top](#features)

## Exchange-correlation functionals

In our package, the XC functional can either be set explicitly using the dft_functional keyword as explained below, or set implicitly according to the XC functional information read from pseudopotential file. The user should ensure that the XC functional set in the INPUT file and the pseudopotential file are consistent. **Currently only LDA and GGA are supported.**

To be specific, we briefly explain the format of the pseudopotential file and the key information it contains. There are a few lines in Si’s GGA pseudopotential file Si_ONCV_PBE-1.0.upf:

```
<PP_HEADER
generated="Generated using ONCVPSP code by D. R. Hamann"
author="Martin Schlipf and Francois Gygi"
date="150105"
comment=""
element="Si"
pseudo_type="NC"
relativistic="scalar"
is_ultrasoft="F"
is_paw="F"
is_coulomb="F"
has_so="F"
has_wfc="F"
has_gipaw="F"
core_correction="F"
functional="PBE"
z_valence=" 4.00"
total_psenergy=" -3.74274958433E+00"
rho_cutoff=" 6.01000000000E+00"
```

The user can set the XC functional type in INPUT file with the parameter, dft_functional:

- none: the functional is specified implicity by the input pseudopotential file
- lda: Perdew-Zunger local density approximation
- pbe: Perdew-Burke-Ernzerhof general gradient approximation

If the functional specified by the user is not consistent with the pseudopotential file, the program will stop with an error message.

Information on the keyword ks_solver is also given in the [list of input variables](input-main.md#dft_functional).

[back to top](#features)
