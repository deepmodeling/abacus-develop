# TB2J

## Introduction

[TB2J](https://github.com/mailhexu/TB2J) is an open-source Python package for the automatic computation of magnetic interactions (including exchange and Dzyaloshinskii-Moriya) between atoms of magnetic crystals from density functional Hamiltonians based on Wannier functions or linear combinations of atomic orbitals. The program is based on Green’s function method with the local rigid spin rotation treated as a perturbation. The ABACUS interface has been added since TB2J version 0.8.0.  

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

## The Heisenberg model
The Heisenberg Hamiltonian in TB2J contains three different parts, which are:
$$
E = -\sum_{i \neq j} \left[ J_{\text{iso}}^{ij} \vec{S}_i \cdot \vec{S}_j + \vec{S}_i J_{\text{ani}}^{ij} \vec{S}_j + \vec{D}_{ij} \cdot (\vec{S}_i \times \vec{S}_j) \right],
$$
where $J_{\text{iso}}^{ij}$ represents the isotropic exchange, $J_{\text{ani}}^{ij}$ represents the symmetric anisotropic exhcange which is a 3 $\times$ 3 tensor with $J^{\text{ani}} = J^{\text{ani,T}}$, $\vec{D}_{ij}$ represents the Dzyaloshinskii-Moriya interaction (DMI). 

> **Note:** Exchange parameters conventions for other Heisenberg Hamiltonian can be found in [Conventions of Heisenberg Model](https://tb2j.readthedocs.io/en/latest/src/convention.html).

## How to use

With the LCAO basis set, TB2J can directly take the output and compute the exchange parameters. For the PW and LCAO-in-PW basis set, the Wannier90 interace can be used instead. In this tutorial we will use LCAO. 

### Collinear calculation without SOC

We take Fe as an example to illustrate how to use ABACUS to generate the input files required for TB2J.

#### 1. Perform ABACUS calculation. 
 `INPUT` file:

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix                         Fe
stru_file                      STRU
kpoint_file                    KPT
pseudo_dir                     ./
orbital_dir                    ./
calculation                    scf
ntype                          1
nspin                          2
symmetry                       1
noncolin                       0
lspinorb                       0

#Parameters (2.PW)
ecutwfc                        100
scf_thr                        1.0e-6
init_chg                       atomic
out_mul                        1

#Parameters (4.Relaxation)
ks_solver                      genelpa
scf_nmax                       200
out_bandgap                    0

#Parameters (5.LCAO)
basis_type                     lcao
gamma_only                     0

#Parameters (6.Smearing)
smearing_method                gauss
smearing_sigma                 0.001

#Parameters (7.Charge Mixing)
mixing_type                    broyden
mixing_beta                    0.2

# Variables related to output information
out_mat_hs2                    1
```

`STRU` file:
```
ATOMIC_SPECIES
Fe  55.845         Fe_ONCV_PBE_FR-1.0.upf

NUMERICAL_ORBITAL
Fe_gga_8au_100Ry_4s2p2d1f.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
2.8660000000      0.0000000000      0.0000000000
0.0000000000      2.8660000000      0.0000000000
0.0000000000      0.0000000000      2.8660000000

ATOMIC_POSITIONS
Direct

Fe
5.0000000000
2
0.0000000000 0.0000000000 0.0000000000 1 1 1 mag 2.5
0.5000000000 0.5000000000 0.5000000000 1 1 1 mag 2.5
```

and `KPT` file:
```
K_POINTS
0
Gamma
8 8 8 0 0 0
```
After the key parameter `out_mat_hs2` is turned on, the Hamiltonian matrix $H(R)$ (in $Ry$) and overlap matrix $S(R)$ will be written into files in the directory `OUT.${suffix}` . In the INPUT, the line:

```
suffix                         Fe
```

specifies the suffix of the output, in this calculation, we set the path to the directory of the DFT calculation, which is the current directory (".") and the suffix to Fe. 

#### 2. Perform TB2J calculation:

```bash
abacus2J.py --path . --suffix Fe --elements Fe  --kmesh 7 7 7
```

This first read the atomic structures from th `STRU` file,  then read the Hamiltonian and the overlap matrices stored in the files named starting from `data-HR-*` and `data-SR-*` files.  It also read the fermi energy from the `OUT.Fe/running_scf.log` file.  

With the command above, we can calculate the $J$ with a $7 \times 7 \times 7$ k-point grid. This allows for the calculation of exchange between spin pairs between $7 \times 7 \times 7$  supercell.  Note: the kmesh is not dense enough for a practical calculation. For a very dense k-mesh, the `--rcut` option can be used to set the maximum distance of the magnetic interactions and thus reduce the computation cost. But be sure that the cutoff is not too small. 

The description of the output files in `TB2J_results` can be found in the [TB2J documentation](https://tb2j.readthedocs.io/en/latest/src/output.html). We can find the exchange parameters with `Fe` by :

```bash
grep "Fe" exchange.out
```

the following contents showing the first, second and third nearest neighbor exchange parameters as 18.6873, 9.9213 and 0.8963 meV, resoectively. More equivalent exchange parameters are also shown.

```
   Fe1   Fe2   ( -1,  -1,  -1) 18.6873   (-1.433, -1.433, -1.433)  2.482
   Fe1   Fe2   ( -1,  -1,   0) 18.6867   (-1.433, -1.433,  1.433)  2.482
   Fe1   Fe2   ( -1,   0,  -1) 18.6866   (-1.433,  1.433, -1.433)  2.482
   Fe1   Fe2   ( -1,   0,   0) 18.6873   (-1.433,  1.433,  1.433)  2.482
   Fe1   Fe2   (  0,  -1,  -1) 18.6873   ( 1.433, -1.433, -1.433)  2.482
   Fe1   Fe2   (  0,  -1,   0) 18.6866   ( 1.433, -1.433,  1.433)  2.482
   Fe1   Fe2   (  0,   0,  -1) 18.6867   ( 1.433,  1.433, -1.433)  2.482
   Fe1   Fe2   (  0,   0,   0) 18.6873   ( 1.433,  1.433,  1.433)  2.482
   Fe2   Fe1   (  0,   0,   0) 18.6873   (-1.433, -1.433, -1.433)  2.482
   Fe2   Fe1   (  0,   0,   1) 18.6867   (-1.433, -1.433,  1.433)  2.482
   Fe2   Fe1   (  0,   1,   0) 18.6866   (-1.433,  1.433, -1.433)  2.482
   Fe2   Fe1   (  0,   1,   1) 18.6873   (-1.433,  1.433,  1.433)  2.482
   Fe2   Fe1   (  1,   0,   0) 18.6873   ( 1.433, -1.433, -1.433)  2.482
   Fe2   Fe1   (  1,   0,   1) 18.6866   ( 1.433, -1.433,  1.433)  2.482
   Fe2   Fe1   (  1,   1,   0) 18.6867   ( 1.433,  1.433, -1.433)  2.482
   Fe2   Fe1   (  1,   1,   1) 18.6873   ( 1.433,  1.433,  1.433)  2.482
   Fe1   Fe1   ( -1,   0,   0)  9.9213   (-2.866,  0.000,  0.000)  2.866
   Fe2   Fe2   ( -1,   0,   0)  9.9213   (-2.866,  0.000,  0.000)  2.866
   Fe1   Fe1   (  0,  -1,   0)  9.9211   ( 0.000, -2.866,  0.000)  2.866
   Fe2   Fe2   (  0,  -1,   0)  9.9211   ( 0.000, -2.866,  0.000)  2.866
   Fe1   Fe1   (  0,   0,  -1)  9.9210   ( 0.000,  0.000, -2.866)  2.866
   Fe2   Fe2   (  0,   0,  -1)  9.9210   ( 0.000,  0.000, -2.866)  2.866
   Fe1   Fe1   (  0,   0,   1)  9.9210   ( 0.000,  0.000,  2.866)  2.866
   Fe2   Fe2   (  0,   0,   1)  9.9210   ( 0.000,  0.000,  2.866)  2.866
   Fe1   Fe1   (  0,   1,   0)  9.9211   ( 0.000,  2.866,  0.000)  2.866
   Fe2   Fe2   (  0,   1,   0)  9.9211   ( 0.000,  2.866,  0.000)  2.866
   Fe1   Fe1   (  1,   0,   0)  9.9213   ( 2.866,  0.000,  0.000)  2.866
   Fe2   Fe2   (  1,   0,   0)  9.9213   ( 2.866,  0.000,  0.000)  2.866
   Fe1   Fe1   ( -1,  -1,   0)  0.8963   (-2.866, -2.866,  0.000)  4.053
   Fe2   Fe2   ( -1,  -1,   0)  0.8963   (-2.866, -2.866,  0.000)  4.053
   Fe1   Fe1   ( -1,   0,  -1)  0.8970   (-2.866,  0.000, -2.866)  4.053
   Fe2   Fe2   ( -1,   0,  -1)  0.8970   (-2.866,  0.000, -2.866)  4.053
```

Several other formats of the exchange parameters are also provided in the `TB2J_results` directory , which can be used in spin dynamics code, e.g. [MULTIBINIT](https://docs.abinit.org/tutorial/spin_model/), [Vampire](https://vampire.york.ac.uk/).

### Non-collinear calculation with SOC

The DMI and anisotropic exchange are result of the SOC, therefore requires the DFT calculation to be done with SOC enabled. To get the full set of exchange parameters, a "rotate and merge" procedure is needed, in which several DFT calculations with either the structure or the spin rotated are needed. For each of the non-collinear calcualtion, we compute the exchange parameters from the DFT calculation with the  same command as in the collienar case. 

```bash
abacus2J.py --path . --suffix Fe --elements Fe  --kmesh 7 7 7
```

And then the "TB2J_merge.py" command can be used to get the final spin interaction parameters. 


### Parameters of abacus2J.py

We can use the command 

```bash
abacus2J.py --help
```

to view the parameters and the usage of them in abacus2J.py.  

### Acknowledgments
We thanks to Xu He to provide critical interface support.
