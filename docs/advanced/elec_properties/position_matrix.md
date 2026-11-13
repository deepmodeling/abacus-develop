# Extracting Position Matrices

In ABACUS, we provide the INPUT keyword [out_mat_r](../input_files/input-main.md#out_mat_r) to write the position matrices into three files, one per Cartesian component: `rrx_nao.csr`, `rry_nao.csr` and `rrz_nao.csr`. They are stored in the matrix output directory, which is `OUT.${suffix}` when [out_app_flag](../input_files/input-main.md#out_app_flag) is 0 or 1 and its `MAT/` subdirectory when out_app_flag is 2; in the latter case a `g${step}` index is inserted before `_nao.csr`, where `${step}` is the one-based ionic-step index. The position matrices are defined as:

$$\langle \chi_\mu|\hat{r}|\chi_\nu\rangle$$

This functionality is not available for gamma_only calculations. If you want to use it in gamma_only calculations, you should turn off [gamma_only](../input_files/input-main.md#gamma_only), and explicitly specifies that gamma point is the only k point in the KPT file.

Each file or each section of the appended file starts with an informative header:

```
Ionic_Step 1
 # print r matrix in real space r(R) x-component
 1 # number of spin directions
 1 # spin index
 6 # number of localized basis
 9 # number of Bravais lattice vector R

 user_defined_lattice
 ...

 #----------------------------------------------------------------------#
 #                               CSR Format                             #
 ...
 #----------------------------------------------------------------------#
```

The header records the ionic step, the matrix label and Cartesian component, the spin metadata, the number of localized basis functions and the number of Bravais lattice vectors, followed by the `user_defined_lattice` cell block and the CSR format description. The rest of the file consists of one sparse block per Bravais lattice vector, such as:

```
-5 -5 -5 12          // R (lattice vector) and number of non-zero elements
 # CSR values
 ...
 # CSR column indices
 ...
 # CSR row pointers
 ...
```

Values wrap every six entries and column indices every sixteen entries. Only Bravais lattice vectors for which at least one of the three components has a non-zero element are written; a component that is empty for a given `R` is still recorded as a valid zero block, so the three files share the same sequence of `R` blocks.

In molecular dynamics (MD) calculations, the output frequency is controlled by [out_freq_ion](../input_files/input-main.md#out_freq_ion). With out_app_flag = 2, each written ionic step has its own files named `rrxg${step}_nao.csr`, `rryg${step}_nao.csr` and `rrzg${step}_nao.csr`, where `${step}` is the one-based ionic-step index; with out_app_flag = 0 or 1 the fixed names `rrx_nao.csr`, `rry_nao.csr` and `rrz_nao.csr` are reused (mode 0 truncates and overwrites them at every write, mode 1 appends one record per written step). The single-shot output of `calculation = get_s` is always written to the fixed names. For example, if we are running a 10-step MD with out_freq_ion = 3 and out_app_flag = 2, the x-component files are `MAT/rrxg1_nao.csr`, `MAT/rrxg4_nao.csr`, `MAT/rrxg7_nao.csr`, and `MAT/rrxg10_nao.csr`, and the y/z components follow the same pattern.

## get_s
We also offer the option of only calculating the position matrix without running SCF. For that purpose, in `INPUT` file we need to set the keyword [calculation](../input_files/input-main.md#calculation) to `get_s`, and [out_mat_r](../input_files/input-main.md#out_mat_r) to `true`.
