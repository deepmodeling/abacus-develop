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

In molecular dynamics (MD) calculations, if [out_app_flag](../input_files/input-main.md#out_app_flag) is set to true, then `data-rR-tr` is written in an append manner. Otherwise, output files will be put in a separate directory, `matrix`, and named as `$x`_data-rR-tr, where `$x` is the number of MD step. In addition, The output frequency is controlled by [out_interval](../input_files/input-main.md#out_interval). For example, if we are running a 10-step MD with out_interval = 3, then `$x` will be 0, 3, 6, and 9.