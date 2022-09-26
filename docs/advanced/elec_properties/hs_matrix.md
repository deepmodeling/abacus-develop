# Extracting Hamiltonian and Overlap Matrices

In ABACUS, we provide the option to write the Hamiltonian and Overlap matrices to files after SCF calculation.

For periodic systems, there are two ways to represent the matrices, the first is to write the entire square matrices for each k point, namely $H(k)$ and $S(k)$; the second is the R space representation, $H(R)$ and $S(R)$, where R is the lattice vector. The two representations are connected by Fourier transform: 

- $H(k)=\sum_R H(R)e^{-ikR}$

and

- $S(k)=\sum_R S(R)e^{-ikR}$

## out_mat_hs

Users may set the keyword `out_mat_hs` to 1 for outputting the k-space matrices. It is available for both gamma_only and multi-k and calculations. Detailed description of the naming and formats of the output files are given [here](../input_files/input-main.md#outmaths).

## out_mat_hs2

The output of R-space matrices is controlled by the keyword `out_mat_hs2`. This functionality is not available for gamma_only calculations. To generate such matrices for gamma only calculations, users should turn off [gamma_only](../input_files/input-main.md#gammaonly), and explicitly specify that gamma point is the only k point in the KPT file.

For a more detailed description of the naming and format of the matrices, refer to this [instruction](../input_files/input-main.md#outmaths2).


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