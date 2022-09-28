# Extracting Density Matrices

ABACUS can output the density matrix by adding the keyword "[out_dm](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-dm)" in INPUT file:
```
out_dm             1
```
After finishing the calculation, the information of the density matrix is stroed in files `OUT.${suffix}/SPIN${spin}_DM`, which looks like:
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

- Note: now this function is valid only for LCAO gamma only calcualtion.