# Extracting Density Matrices

ABACUS can output the density matrix by adding the keyword "[out_dmk](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out_dmk)" in INPUT file:
```
out_dmk    1
```
After finishing the calculation, the density matrix is written into `OUT.${suffix}/`.

For current develop versions:
- gamma-only (`gamma_only = 1`): `dmg1_nao.txt` (`nspin=1/4`) or `dms1g1_nao.txt` and `dms2g1_nao.txt` (`nspin=2`)
- multi-k (`gamma_only = 0`): `dmk1g1_nao.txt`, `dmk2g1_nao.txt`, ... (`nspin=1/4`) or `dmk1s1g1_nao.txt`, `dmk1s2g1_nao.txt`, ... (`nspin=2`)

Here `g{istep}` denotes the geometry/step index in the output filename.

For 3.10-LTS, the corresponding keyword is `out_dm`, and file names follow the legacy style such as `SPIN1_DM` and `SPIN2_DM`.

The file content looks like:
```
test
 5.39761
 0.5 0.5 0
 0.5 0 0.5
 0 0.5 0.5
 Si
 2
Direct
 0 0 0
 0.25 0.25 0.25

 1
 0.570336288801065 (fermi energy)
  26 26

 3.904e-01 1.114e-02 2.050e-14 1.655e-13 1.517e-13 -7.492e-15 -1.729e-14 5.915e-15
 -9.099e-15 2.744e-14 3.146e-14 6.631e-15 2.594e-15 3.904e-01 1.114e-02 -7.395e-15
 ...
 ```
The first 5 lines are the informations of lattice, in order: \
&emsp;lattice name (if keyword [latname](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#latname) is not specified in INPUT, this will be "test"), \
&emsp;lattice constance with unit in angstrom, \
&emsp;lattice vector a, \
&emsp;lattice vector b, \
&emsp;lattice vector c. \
The following lines are about the elements and coordinates, in order: all elements, the atom number of each elements, the type of coordinate, the coordinates.\
After a blank line, the output is the values of NSPIN and fermi energy.\
The following line is dimension of the density matrix, and the rest lines are the value of each matrix element.

The examples can be found in [examples/density_matrix](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/density_matrix)

- Note: Version difference summary:
    - develop: `out_dmk` supports both gamma-only and multi-k-point output.
    - 3.10-LTS: use `out_dm`.

## Real-space Density Matrix (CSR format)

ABACUS can also output the real-space density matrix DM(R) in CSR (Compressed Sparse Row) format by setting:
```
out_dmr    1
```
This feature is only valid for multi-k calculations (`gamma_only = 0`).

After the calculation, the density matrix files are written to `OUT.${suffix}/`:
- develop naming pattern: `dmr{s}{spin index}{g}{geometry index}{_nao}.csr`
- `nspin=1`: `dmrs1_nao.csr`
- `nspin=2` (spin-polarized): `dmrs1_nao.csr` (spin-up) and `dmrs2_nao.csr` (spin-down)

For 3.10-LTS, the corresponding keyword is `out_dm1`, and the file names are `data-DMR-sparse_SPIN0.csr` and `data-DMR-sparse_SPIN1.csr`, etc.

These files can be used to restart calculations by setting `init_chg dm` in the INPUT file together with `read_file_dir` pointing to the directory containing the CSR files. This is supported for both `nspin=1` and `nspin=2`.
