# Extracting Hamiltonian and Overlap Matrices

In ABACUS, we provide the option to write the Hamiltonian and Overlap matrices to files after SCF calculation.

For periodic systems, there are two ways to represent the matrices, the first is to write the entire square matrices for each k point, namely $H(k)$ and $S(k)$; the second is the R space representation, $H(R)$ and $S(R)$, where R is the lattice vector. The two representations are connected by Fourier transform: 

- $H(k)=\sum_R H(R)e^{-ikR}$

and

- $S(k)=\sum_R S(R)e^{-ikR}$

## out_mat_hs

Users may set the keyword [out_mat_hs](../input_files/input-main.md#out_mat_hs) to true for outputting the upper triangular part of the Hamiltonian matrices and overlap matrices for each k point into files in the directory `OUT.${suffix}`. It is available for both gamma_only and multi-k calculations. 

The files are named `data-$k-H` and `data-$k-S`, where `$k` is a composite index consisting of the k point index as well as the spin index. The corresponding sequence of the orbitals can be seen in [Basis Set](../pp_orb.md#basis-set).

For nspin = 1 and nspin = 4 calculations, there will be only one spin component, so `$k` runs from 0 up to `$nkpoints-1`. For nspin = 2, `$k` runs from `2*$nkpoints-1`. In the latter case, the files are arranged into blocks of up and down spins. For example, if there are 3 k points, then we have the following correspondence:

  - data-0-H : 1st k point, spin up
  - data-1-H : 2nd k point, spin up
  - data-2-H : 3rd k point, spin up
  - data-3-H : 1st k point, spin down
  - data-4-H : 2nd k point, spin down
  - data-5-H : 3rd k point, spin down

As for information on the k points, one may look for the `SETUP K-POINTS` section in the running log.

The first number of the first line in each file gives the size of the matrix, namely, the number of atomic basis functions in the system.

The rest of the file contains the upper triangular part of the specified matrices. For multi-k calculations, the matrices are Hermitian and the matrix elements are complex; for gamma-only calculations, the matrices are symmetric and the matrix elements are real.

## out_mat_hs2

The output of R-space matrices is controlled by the keyword [out_mat_hs2](../input_files/input-main.md#out_mat_hs2). This functionality is not available for gamma_only calculations. To generate such matrices for gamma only calculations, users should turn off [gamma_only](../input_files/input-main.md#gamma_only), and explicitly specify that gamma point is the only k point in the KPT file.

For single-point SCF calculations, if nspin = 1 or nspin = 4, two files `data-HR-sparse_SPIN0.csr` and `data-SR-sparse_SPIN0.csr` are generated, which contain the Hamiltonian matrix $H(R)$ and overlap matrix $S(R)$ respectively. For nspin = 2, three files `data-HR-sparse_SPIN0.csr` and `data-HR-sparse_SPIN1.csr` and `data-SR-sparse_SPIN0.csr` are created, where the first two contain $H(R)$ for spin up and spin down, respectively.

As for molecular dynamics calculations, the format is controlled by [out_interval](../input_files/input-main.md#out_interval) and [out_app_flag](../input_files/input-main.md#out_app_flag) in the same manner as the position matrix as detailed in [out_mat_r](../input_files/input-main.md#out_mat_r).

Each file or each section of the appended file starts with three lines, the first gives the current ion/md step, the second gives the dimension of the matrix, and the last indicates how many different `R` are in the file.

The rest of the files are arranged in blocks. Each block starts with a line giving the lattice vector `R` and the number of nonzero matrix elements, such as:

```
-3 1 1 1020
```

which means there are 1020 nonzero elements in the (-3,1,1) cell.

If there is no nonzero matrix element, then the next block starts immediately on the next line. Otherwise, there will be 3 extra lines in the block, which gives the matrix in CSR format. According to Wikipedia:

The CSR format stores a sparse m × n matrix M in row form using three (one-dimensional) arrays (V, COL_INDEX, ROW_INDEX). Let NNZ denote the number of nonzero entries in M. (Note that zero-based indices shall be used here.)

  - The arrays V and COL_INDEX are of length NNZ, and contain the non-zero values and the column indices of those values respectively.
  - The array ROW_INDEX is of length m + 1 and encodes the index in V and COL_INDEX where the given row starts. This is equivalent to ROW_INDEX[j] encoding the total number of nonzeros above row j. The last element is NNZ , i.e., the fictitious index in V immediately after the last valid index NNZ - 1.

## get_S
We also offer the option of only calculating the overlap matrix without running SCF. For that purpose, in `INPUT` file we need to set the value keyword [calculation](../input_files/input-main.md#calculation) to be `get_S`.

A file named `SR.csr` will be generated in the working directory, which contains the overlap matrix.

## examples
We provide [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/matrix_hs) of outputting the matrices. There are four examples:

- out_hs2_multik : writing H(R) and S(R) for multi-k calculation
- out_hs_gammaonly : writing H(k) and S(k) for gamma-only calculation
- out_hs_multik : writing H(k) and S(k) for multi-k calculation
- out_s_multik : running get_S for multi-k calculation

Reference output files are provided in each directory.