# PYATB

## Introduction

[PYATB](https://github.com/pyatb/pyatb) (Python ab initio tight binding simulation package) is an open-source software package designed for computing electronic structures and related properties based on the ab initio tight binding Hamiltonian. The Hamiltonian can be directly obtained after conducting self-consistent calculations with ABACUS using numerical atomic orbital (NAO) bases. The package comprises three modules - Bands, Geometric, and Optical, each providing a comprehensive set of tools for analyzing different aspects of a material's electronic structure.

## Installation

```bash
git clone https://github.com/pyatb/pyatb.git
cd pyatb
python setup.py install --record log
```

To customize the `setup.py` file, you must make changes to the **CXX** and **LAPACK_DIR** variables in line with your environment. **CXX** denotes the C++ compiler you intend to use, for instance, icpc (note that it should not be the mpi version). Furthermore, **LAPACK_DIR** is used to specify the Intel MKL path.

After completing the installation process, you can access the `pyatb` executable and corresponding module, which can be imported using the `import pyatb` command.

## How to use

We take Bi$_2$Se$_3$ as an example to illustrate how to use ABACUS to generate the tight binding Hamiltonian required for PYATB, and then perform calculations related to PYATB functions.

1. Perform ABACUS self consistent calculation:

```
INPUT_PARAMETERS

# System variables
suffix                Bi2Se3
ntype                 2
calculation           scf
esolver_type          ksdft
symmetry              1
init_chg              atomic

# Plane wave related variables
ecutwfc               100

# Electronic structure
basis_type            lcao
ks_solver             genelpa
nspin                 4
smearing_method       gauss
smearing_sigma        0.02
mixing_type           broyden
mixing_beta           0.7
scf_nmax              200
scf_thr               1e-8
lspinorb              1
noncolin              0

# Variables related to output information
out_chg               1
out_mat_hs2           1
out_mat_r             1
```

After the key parameters `out_mat_hs2` and `out_mat_r` are turned on, ABACUS will generate files containing the Hamiltonian matrix $H(R)$, overlap matrix $S(R)$, and dipole matrix $r(R)$ after completing the self-consistent calculation. These parameters can be found in the ABACUS `INPUT` file.

2. Copy the HR, SR, and rR files output by ABACUS's self-consistent calculation, which are located in the `OUT*` directory and named `data-HR-sparse_SPIN0.csr`, `data-SR-sparse_SPIN0.csr`, and `data-rR-sparse.csr`, respectively. Copy these files to the working directory and write the `Input` file for PYATB:

```
INPUT_PARAMETERS
{
    nspin                          4
    package                        ABACUS
    fermi_energy                   9.557219691497478
    fermi_energy_unit              eV
    HR_route                       data-HR-sparse_SPIN0.csr
    SR_route                       data-SR-sparse_SPIN0.csr
    rR_route                       data-rR-sparse.csr
    HR_unit                        Ry
    rR_unit                        Bohr
    max_kpoint_num                 8000
}

LATTICE
{
    lattice_constant               1.8897162
    lattice_constant_unit          Bohr
    lattice_vector
    -2.069  -3.583614  0.000000
     2.069  -3.583614  0.000000
     0.000   2.389075  9.546667
}

BAND_STRUCTURE
{
    wf_collect                     0
    kpoint_mode                    line
    kpoint_num                     5
    high_symmetry_kpoint
    0.00000 0.00000 0.0000 100  # G
    0.00000 0.00000 0.5000 100  # Z
    0.50000 0.50000 0.0000 100  # F
    0.00000 0.00000 0.0000 100  # G
    0.50000 0.00000 0.0000 1    # L
}
```

For specific input file writing, please refer to PYATB's quick start.

3. Perform PYATB calculation:

```
export OMP_NUM_THREADS=2
mpirun -np 6 pyatb
```

After the calculation is completed, the band structure data and figures of Bi$_2$Se$_3$ can be found in the `Out/Band_Structure` folder.