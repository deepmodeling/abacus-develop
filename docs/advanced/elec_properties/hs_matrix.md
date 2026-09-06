# Extracting Hamiltonian and Overlap Matrices

In ABACUS, we provide the option to write the Hamiltonian and Overlap matrices to files after SCF calculations.

For periodic systems, there are two ways to construct the matrices. The reciprocal-space representation writes the entire square matrices $H(k)$ and $S(k)$ for each $k$ point in the Brillouin zone. The real-space representation writes $H(R)$ and $S(R)$ indexed by the Bravais lattice vector $R$. The two representations are connected by Fourier transform:

- $H(k)=\sum_R H(R)e^{-ikR}$

and

- $S(k)=\sum_R S(R)e^{-ikR}$

## out_hsk

Use [out_hsk](../input_files/input-main.md#out_hsk) to print the upper triangular part of the Hamiltonian and overlap matrices for each k point into `OUT.${suffix}`. It is available for both gamma-only and multi-k calculations. The format value is:

| Value | Format |
| --- | --- |
| `0` | Disabled |
| `1` | Text; an optional second value controls precision, for example `out_hsk 1 12` |
| `2` | Native binary `.dat` output |
| `3` | Reserved for H(k)/S(k) NPZ output; not implemented |

The legacy keyword `out_mat_hs 1 [precision]` remains supported as an alias for `out_hsk 1 [precision]`. If both names are present, `out_hsk` takes precedence.

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

### Native Binary Format

For `out_hsk 2`, the filenames in the table above use `.dat` instead of `.txt`. Each matrix record is written without padding or a self-describing header:

1. the matrix dimension as a native C++ `int`;
2. the upper-triangular elements in row-major order, from `(0,0)` through `(0,N-1)`, then `(1,1)` through `(1,N-1)`, and so on;
3. each gamma-only element as one native `double`, or each multi-k/spinor element as two consecutive native `double` values containing the real and imaginary parts.

The file therefore contains `sizeof(int) + N(N+1)/2 * sizeof(double)` bytes for a gamma-only record and `sizeof(int) + N(N+1) * sizeof(double)` bytes for a multi-k or spinor record. The format uses the host integer representation and byte order and is intended for readers using a compatible ABI.

When `out_app_flag` is true, the first ionic step truncates the shared file and every later step appends another complete record. When it is false, each ionic step has a separate filename containing `g${step}` before `_nao`.

## out_hsr

The output of $H(R)$ and $S(R)$ matrices is controlled by [out_hsr](../input_files/input-main.md#out_hsr). It is available for both gamma-only and multi-k LCAO calculations:

| Value | Format |
| --- | --- |
| `0` | Disabled |
| `1` | Text CSR; an optional second value controls precision, for example `out_hsr 1 12` |
| `2` | Native binary CSR using `.dat` files |
| `3` | NPZ: `hrs1_nao.npz`, `hrs2_nao.npz` when needed, and `sr_nao.npz` |

The legacy keywords `out_mat_hs2 1 [precision]` and `out_hsr_npz 1` remain supported as aliases for text and NPZ output respectively. If `out_hsr` is present together with either legacy keyword, `out_hsr` takes precedence.

For a multi-k calculation, the files contain the individual real-space blocks stored for the Bravais lattice vectors $R$. For a gamma-only calculation, ABACUS stores the real-space contributions in a folded representation. Text CSR, native binary, and NPZ output write this internal representation directly: all stored $R$-space contributions are summed into a single block labelled `R = (0, 0, 0)`.

The folded gamma-only output is sufficient to inspect the matrix used by the gamma-only real-space container, but it does not retain the original lattice-vector resolution and cannot be used to interpolate matrices at arbitrary k points. Terms that are added only while constructing $H(k)$, rather than stored in the internal $H(R)$ container, are not guaranteed to be present. Use [out_hsk](../input_files/input-main.md#out_hsk) when the final $H(\Gamma)$ and $S(\Gamma)$ matrices are required.

### Text CSR Format

The H(R) and S(R) matrices are output in standard Compressed Sparse Row (CSR) format, matching the format used by `out_dmr`.

For single-point SCF calculations:
- **nspin = 1**: Two files `hrs1_nao.csr` and `sr_nao.csr` are generated, containing the Hamiltonian matrix $H(R)$ and overlap matrix $S(R)$ respectively.
- **nspin = 2**: Three files `hrs1_nao.csr`, `hrs2_nao.csr`, and `sr_nao.csr` are created, where the first two files correspond to $H(R)$ for spin up and spin down, respectively.
- **nspin = 4**: Multi-k calculations generate `hrs1_nao.csr` and `sr_nao.csr`. The gamma-only algorithm itself does not support `nspin = 4`.

In gamma-only mode, every generated file reports one Bravais lattice vector and contains one CSR block for `0 0 0`. The header also contains:

```text
# representation: gamma-only folded matrix; stored R-space contributions are summed into R = (0, 0, 0)
```

### Native Binary CSR Format

Set `out_hsr 2` to write the same H(R) and S(R) matrix sets with a `.dat` suffix. For example, an `nspin = 2` calculation writes `hrs1_nao.dat`, `hrs2_nao.dat`, and `sr_nao.dat`. When `out_app_flag` is false, the one-based ionic step is included before `_nao`, for example `hrs1g1_nao.dat` and `srg1_nao.dat`.

Each ionic step is a complete record with no padding or self-describing header:

1. Native `int`: zero-based ionic step, matrix dimension, number of R blocks.
2. For every R block in lexicographic `(Rx, Ry, Rz)` order, four native `int` values: `Rx`, `Ry`, `Rz`, and `nnz`.
3. `nnz` matrix values in CSR order. Real matrices use one native `double`; complex matrices use consecutive real and imaginary `double` values.
4. `nnz` native `int` column indices.
5. `dimension + 1` native `long long` row pointers.

All R blocks stored by the internal HContainer are present, including blocks with zero nonzero values. The sparse threshold is `1e-10`, matching text CSR output. The format uses the host integer representation and byte order and therefore requires a compatible ABI. It does not contain unit-cell, spin, or matrix-label metadata; those are determined by the calculation input and filename.

When `out_app_flag` is true, the first ionic step truncates the shared file and later ionic steps append complete records. Otherwise every ionic step is written to its own file.

### NPZ Format

Set `out_hsr 3` to write `hrs1_nao.npz`, `hrs2_nao.npz` when a second spin channel is present, and `sr_nao.npz`. Matrix entry names include the atom-pair indices and the three components of $R$. Multi-k calculations retain the stored $R$ blocks, while gamma-only calculations contain only matrix entry names ending in `_0_0_0`.

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

Use `out_hsr 1 12` to output text CSR files with 12-digit precision (default is 8). Precision is ignored for native binary and NPZ output.

For calculations involving ionic movements, the output frequency of the matrix is controlled by [out_freq_ion](../input_files/input-main.md#out_freq_ion) and [out_app_flag](../input_files/input-main.md#out_app_flag). 

## get_s
We also offer the option of only calculating the overlap matrix without running SCF. For that purpose, in `INPUT` file we need to set the value keyword [calculation](../input_files/input-main.md#calculation) to be `get_s`.

A file named `sr_nao.csr` will be generated in `OUT.${suffix}`, which contains the overlap matrix.

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
