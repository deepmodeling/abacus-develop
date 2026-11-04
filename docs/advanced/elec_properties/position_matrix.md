# Extracting Position Matrices

In ABACUS, we provide the INPUT keyword [out_mat_r](../input_files/input-main.md#out_mat_r) to write the position matrices into a file named `data-rR-tr` in the directory `OUT.${suffix}`. The position matrices is defined as:

$$\langle \chi_\mu|\hat{r}|\chi_\nu\rangle$$

This functionality is not available for gamma_only calculations. If you want to use it in gamma_only calculations, you should turn off [gamma_only](../input_files/input-main.md#gamma_only), and explicitly specifies that gamma point is the only k point in the KPT file.

Each file or each section of the appended file starts with "STEP: " followed by the current ion/md step, then the second line starts with "Matrix Dimension of $r(R)$: " followed by the dimension of the matrix, and the third line starts with "Matrix number of $r(R)$: " followed by the matrix number. The rest of the format is arranged into blocks, such as:

```
-5 -5 -5    //R (lattice vector)
...
-5 -5 -4    //R (lattice vector)
...
-5 -5 -3    //R (lattice vector)
```

Each block here contains the matrix for the corresponding cell. There are three columns in each block, giving the matrix elements in x, y, z directions, respectively. There are altogether nbasis * nbasis lines in each block, which emulates the matrix elements.

In molecular dynamics (MD) calculations, the output frequency is controlled by [out_freq_ion](../input_files/input-main.md#out_freq_ion). With out_app_flag = 2, each written ionic step has its own file named `rrg${step}_nao.txt`, where `${step}` is the one-based ionic-step index; with out_app_flag = 0 or 1 the fixed name `rr_nao.txt` is reused (mode 0 truncates and overwrites it at every write, mode 1 appends one record per written step). The single-shot output of `calculation = get_s` is always named `rr_nao.txt`. The files are stored in `OUT.${suffix}` when [out_app_flag](../input_files/input-main.md#out_app_flag) is 0 or 1, and in its `MAT/` subdirectory when out_app_flag is 2. For example, if we are running a 10-step MD with out_freq_ion = 3 and out_app_flag = 2, the files are `MAT/rrg1_nao.txt`, `MAT/rrg4_nao.txt`, `MAT/rrg7_nao.txt`, and `MAT/rrg10_nao.txt`.

## get_s
We also offer the option of only calculating the position matrix without running SCF. For that purpose, in `INPUT` file we need to set the keyword [calculation](../input_files/input-main.md#calculation) to `get_s`, and [out_mat_r](../input_files/input-main.md#out_mat_r) to `true`.
