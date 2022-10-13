# Accelerating the Calculation

In ABACUS, we provide a few methods for accelerating the calculation. The parameters are usually set as default for calculations where there is not extreme concern for efficiency, as some of them may produce numerical issues under certain circumstances. In short, methods in this section should be used with care. It is better to calibrate the results against the default setting.

## K-point Parallelization

In ABACUS, we offer k-point parallelization for calculations with PW basis, which should increase the efficiency when a large k-point mesh is used.

To use k-point parallelization, users may set keyword [kpar](../input_files/input-main.md#kpar) to be larger than 1.

> Note: It has been observed that k-point parallelization cannot work in conjunction with Davidson diagonalization.

## K-point Symmetry

Inclusion of k-point symmetry helps increasing the efficiency of calculations by reducing the effective number of k-points used. To turn on k-point symmetry, users may set keyword [symmetry](../input_files/input-main.md#symmetry) to be 1.

> Note: In ABACUS we only support point-group symmetry but not space-group symmetry.

## Accelerating Grid Integration

For LCAO calculation, the matrix elements of the local potential is evaluated using grid integration. In grid integration, we group real-space FFT grid points into boxes of dimension bx * by * bz, and then proceed with the boxes as the basis unit of calculation.

Setting [bx, by, bz](../input_files/input-main.md#bx-by-bz) to be values other than default might help with the efficiency of grid integration.

> Note: the choice of bx, by, bz should be integer factors of the dimension of the real space FFT grid in each direction.

## Low Dimension Materials

In grid integration, we chose to parallelize the grid points along the z direction. Therefore, when using LCAO calculation for low dimension materials, it is recommended to put the material more evenly in z direction to avoid imbalanced workload on different MPI threads.

Namely, when calculating 2D materials, it is better to put the material in xz or yz direction; while for 1D materials, it is better to align the material with the z direction.