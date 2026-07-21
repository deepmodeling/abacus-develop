# TB2J Interface

This directory contains the interface between ABACUS and TB2J, an open-source Python package for the automatic computation of magnetic interactions.

## Introduction

[TB2J](https://github.com/mailhexu/TB2J) is an open-source Python package for the automatic computation of magnetic interactions (including exchange and Dzyaloshinskii-Moriya) between atoms of magnetic crystals from density functional Hamiltonians based on Wannier functions or linear combinations of atomic orbitals. The program is based on Green’s function method with the local rigid spin rotation treated as a perturbation. The ABACUS interface has been added since TB2J version 0.8.0.  

The Heisenberg Hamiltonian in TB2J contains three different parts, which are:

$E = -\sum_{i \neq j} \left[ J_{\text{iso}}^{ij} \vec{S}_i \cdot \vec{S}_j + \vec{S}_i J_{\text{ani}}^{ij} \vec{S}_j + \vec{D}_{ij} \cdot (\vec{S}_i \times \vec{S}_j) \right],$

where $J_{\text{iso}}^{ij}$ represents the isotropic exchange, $J_{\text{ani}}^{ij}$ represents the symmetric anisotropic exhcange which is a 3 $\times$ 3 tensor with $J^{\text{ani}} = J^{\text{ani,T}}$, $\vec{D}_{ij}$ represents the Dzyaloshinskii-Moriya interaction (DMI). 

> **Note:** Exchange parameters conventions for other Heisenberg Hamiltonian can be found in [Conventions of Heisenberg Model](https://tb2j.readthedocs.io/en/latest/src/convention.html).

For more information, see the documentation on https://tb2j.readthedocs.io/en/latest/

## Installation

The most easy way to install TB2J is to use pip:

```bash
pip install TB2J
```

You can also download TB2J from the github page, and install with:

```bash
git clone https://github.com/mailhexu/TB2J.git
cd TB2J
python setup.py install
```

## How to use

With the LCAO basis set, TB2J can directly take the output and compute the exchange parameters. For the PW and LCAO-in-PW basis set, the Wannier90 interace can be used instead. In this tutorial we will use LCAO. 

The `example01` directory contains a simple example demonstrating how to use the ABACUS to generate the input files required for TB2J for iron (Fe), then perform TB2J to calculate magnetic interactions.

#### 1. Perform ABACUS calculation. 
 
After the key parameter `out_mat_hs2` is turned on, the Hamiltonian matrix $H(R)$ (in $Ry$) and overlap matrix $S(R)$ will be written into files in the directory `OUT.${suffix}` . In the INPUT, the line:

```
suffix                         Fe
```

specifies the suffix of the output, in this calculation, we set the path to the directory of the DFT calculation, which is the current directory (".") and the suffix to Fe.

> **Note (ABACUS v3.9.0.25+):** Starting from ABACUS v3.9.0.25, the output format has changed to standard CSR format with filenames `hrs1_nao.csr`, `hrs2_nao.csr` (for nspin=2), and `srs1_nao.csr`. The parameter `out_mat_hs2` now supports optional precision control: `out_mat_hs2 1 8` (default 8 digits). TB2J v0.9.0+ is required to read the new format. For older TB2J versions, please use ABACUS v3.8.x or earlier. 

#### 2. Perform TB2J calculation:

```bash
abacus2J.py --path . --suffix Fe --elements Fe  --kmesh 7 7 7
```

This first reads the atomic structures from the `STRU` file, then reads the Hamiltonian and overlap matrices.  It also reads the fermi energy from the `OUT.Fe/running_scf.log` file.  
> **Note:** For ABACUS v3.9.0.25+, the matrices are stored in `hrs1_nao.csr`, `hrs2_nao.csr` (nspin=2), and `srs1_nao.csr` files. For older versions, they are in `data-HR-*` and `data-SR-*` files.

With the command above, we can calculate the $J$ with a $7 \times 7 \times 7$ k-point grid. This allows for the calculation of exchange between spin pairs between $7 \times 7 \times 7$  supercell.  Note: the kmesh is not dense enough for a practical calculation. For a very dense k-mesh, the `--rcut` option can be used to set the maximum distance of the magnetic interactions and thus reduce the computation cost. But be sure that the cutoff is not too small. 

The description of the output files in `TB2J_results` can be found in the [TB2J documentation](https://tb2j.readthedocs.io/en/latest/src/output.html). 

Several other formats of the exchange parameters are also provided in the `TB2J_results` directory , which can be used in spin dynamics code, e.g. [MULTIBINIT](https://docs.abinit.org/tutorial/spin_model/), [Vampire](https://vampire.york.ac.uk/).

### Parameters of abacus2J.py

We can use the command 

```bash
abacus2J.py --help
```

to view the parameters and the usage of them in abacus2J.py.  

### Acknowledgments
We thanks to Xu He  and Zhenxiong Shen to provide critical interface support.
