# Extracting Hamiltonian and Overlap Matrices

In ABACUS, we provide the option to write the Hamiltonian and Overlap matrices to files after SCF calculations.

For periodic systems, there are two ways to construct the matrices, the first is to write the entire square matrices for each $k$ point in the Brillouin zone, namely $H(k)$ and $S(k)$; the second one is the real space representation, $H(R)$ and $S(R)$, where R is the Bravis lattice vector. The two representations are connected by Fourier transform:

- $H(k)=\sum_R H(R)e^{-ikR}$

and

- $S(k)=\sum_R S(R)e^{-ikR}$

## out_mat_hs

Users can set the keyword [out_mat_hs](../input_files/input-main.md#out_mat_hs) to true to print the upper triangular part of the Hamiltonian matrices and overlap matrices for each k point into files in the directory `OUT.${suffix}`. It is available for both gamma_only and multi-k calculations. 

The $H(k)$ and $S(k)$ matrices are stored with numerical atomic orbitals as basis, and the corresponding sequence of the numerical atomic orbitals can be seen in [Basis Set](../pp_orb.md#basis-set).

As for information on the k points, one may look for the `SETUP K-POINTS` section in the running log.

The output filenames depend on the k-point algorithm and `nspin`:

| Calculation mode | `nspin` | Hamiltonian files | Overlap files |
| --- | --- | --- | --- |
| `gamma_only = 1` | 1 | `hk_nao.txt` | `sk_nao.txt` |
| `gamma_only = 1` | 2 | `hks1_nao.txt`, `hks2_nao.txt` | `sk_nao.txt` |
| `gamma_only = 0` | 1 | `hk${k}_nao.txt` | `sk${k}_nao.txt` |
| `gamma_only = 0` | 2 | `hk${k}s1_nao.txt`, `hk${k}s2_nao.txt` | `sk${k}_nao.txt` |
| `gamma_only = 0` | 4 | `hk${k}s4_nao.txt` | `sk${k}_nao.txt` |

Here `${k}` is the one-based k-point index. For `nspin = 2`, the overlap matrix is spin-independent, so only one overlap file is written for each physical k point. The gamma-only algorithm does not support `nspin = 4`; use the multi-k algorithm with an explicit Gamma-only `KPT` file for a noncollinear calculation at Gamma.

When `out_app_flag` is false, `g${step}` is inserted before `_nao`, where `${step}` is the one-based ionic-step index. For example, the first spin channel at the first k point and first ionic step is written to `hk1s1g1_nao.txt`.

Each output block starts with a comment header containing the one-based ionic-step index, filename, `gamma only` flag, and matrix dimensions. It is followed by `Row 1`, `Row 2`, and so on. Each row contains the matrix elements from the diagonal through the upper triangle.

For multi-k calculations, the matrices are Hermitian and each matrix element is written as `(real,imag)`. For gamma-only calculations, the matrices are symmetric and the matrix elements are written as real numbers.

## out_mat_hs2

The output of $H(R)$ and $S(R)$ matrices is controlled by the keyword [out_mat_hs2](../input_files/input-main.md#out_mat_hs2). It is available for both gamma-only and multi-k LCAO calculations.

For a multi-k calculation, the files contain the individual real-space blocks stored for the Bravais lattice vectors $R$. For a gamma-only calculation, ABACUS stores the real-space contributions in a folded representation. `out_mat_hs2` writes this internal representation directly: all stored $R$-space contributions are summed into a single block labelled `R = (0, 0, 0)`.

The folded gamma-only output is sufficient to inspect the matrix used by the gamma-only real-space container, but it does not retain the original lattice-vector resolution and cannot be used to interpolate matrices at arbitrary k points. Terms that are added only while constructing $H(k)$, rather than stored in the internal $H(R)$ container, are not guaranteed to be present. Use [out_mat_hs](../input_files/input-main.md#out_mat_hs) when the final $H(\Gamma)$ and $S(\Gamma)$ matrices are required.

### Output Format

The H(R) and S(R) matrices are output in standard Compressed Sparse Row (CSR) format, matching the format used by `out_dmr`.

For single-point SCF calculations:
- **nspin = 1**: Two files `hrs1_nao.csr` and `srs1_nao.csr` are generated, containing the Hamiltonian matrix $H(R)$ and overlap matrix $S(R)$ respectively.
- **nspin = 2**: Three files `hrs1_nao.csr`, `hrs2_nao.csr`, and `srs1_nao.csr` are created, where the first two files correspond to $H(R)$ for spin up and spin down, respectively.
- **nspin = 4**: Multi-k calculations generate `hrs1_nao.csr` and `srs1_nao.csr`. The gamma-only algorithm itself does not support `nspin = 4`.

In gamma-only mode, every generated file reports one Bravais lattice vector and contains one CSR block for `0 0 0`. The header also contains:

```text
# representation: gamma-only folded matrix; stored R-space contributions are summed into R = (0, 0, 0)
```

### File Structure

Each file starts with a header:
```
 --- Ionic Step 1 ---
 # print H matrix in real space H(R)
 1 # number of spin directions
 1 # spin index
 100 # number of localized basis
 50 # number of Bravais lattice vector R

[UnitCell information]

#----------------------------------------------------------------------#
#                               CSR Format                             #
...
 0 0 0 5
 # CSR values
 1.234e-01 2.345e-02 ...
 # CSR column indices
 0 5 10 ...
 # CSR row pointers
 0 3 7 ...
```

The CSR format stores a sparse m × n matrix M in row form using three arrays (values, column indices, row pointers). According to Wikipedia:

- The arrays **values** and **column indices** are of length NNZ (number of nonzero entries), and contain the non-zero values and the column indices of those values respectively.
- The array **row pointers** is of length m + 1 and encodes the index where each row starts. The last element is NNZ.

### Precision Control

Use `out_mat_hs2 1 12` to output with 12-digit precision (default is 8).

For calculations involving ionic movements, the output frequency of the matrix is controlled by [out_freq_ion](../input_files/input-main.md#out_freq_ion) and [out_app_flag](../input_files/input-main.md#out_app_flag). 

## get_s
We also offer the option of only calculating the overlap matrix without running SCF. For that purpose, in `INPUT` file we need to set the value keyword [calculation](../input_files/input-main.md#calculation) to be `get_s`.

A file named `sr_nao.csr` will be generated in the working directory, which contains the overlap matrix.

> When `nspin` is set to 1 or 2, the dimension of the overlap matrix is nlocal $\times$ nlocal, where nlocal is the total number of numerical atomic orbitals. 
These numerical atomic orbitals are ordered from outer to inner loop as atom, angular quantum number $l$, zeta (multiple radial orbitals corresponding to each $l$), and magnetic quantum number $m$. 
When `nspin` is set to 4, the dimension of the overlap matrix is (2 $\times$ nlocal) $\times$ (2 $\times$ nlocal). In this case, the numerical atomic orbitals are ordered from outer to inner loop as atom, angular quantum number $l$, zeta (multiple radial orbitals corresponding to each $l$), magnetic quantum number $m$, and npol (index of spin, ranges from 0 to 1).


## examples
We provide [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/10_hs_matrix) of outputting the matrices.

- `03_out_hsk_gamma`: writing H(k) and S(k) for a gamma-only calculation
- `04_out_hsk_multik`: writing H(k) and S(k) for a multi-k calculation
- `01_out_hsr_multik` and `02_out_hsr_multik`: writing H(R) and S(R) for a multi-k calculation
- `05_gets`: running `calculation = get_s` to obtain the overlap matrix

Reference output files are provided in each directory.
