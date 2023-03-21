# Extracting Charge Density

ABACUS can output the charge density by adding the keyword [out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg) in INPUT file:

```
out_chg             1
```

After finishing the calculation, the information of the charge density is stroed in files `OUT.${suffix}/SPIN${spin}_CHG.cube`, which can be used to do visualization.
The SPIN${spin}_CHG.cube file looks like:

```
Cubefile created from ABACUS SCF calculation
2 (nspin) 0.914047 (fermi energy, in Ry)
2 0.0 0.0 0.0 
27 0.222222 0 0
27 0 0.222222 0
27 0 0 0.222222
 26 16 0 0 0
 26 16 3 3 3
 6.63594288898e-01 8.42344790519e-01 1.16349621677e+00 1.18407505276e+00 8.04461725175e-01 3.77164277045e-01
 1.43308127341e-01 5.93894932356e-02 3.23036576611e-02 2.08414809212e-02 1.51271068218e-02 1.27012859512e-02
 1.15620162933e-02 1.08593210023e-02 1.08593210023e-02 1.15620162933e-02 1.27012859512e-02 1.51271068218e-02
 2.08414809212e-02 3.23036576611e-02 5.93894932356e-02 1.43308127341e-01 3.77164277045e-01 8.04461725175e-01
 1.18407505276e+00 1.16349621677e+00 8.42344790519e-01
 8.42344790519e-01 9.86194056340e-01 1.21545550606e+00 1.14987597026e+00 7.50033272229e-01 3.46047149862e-01
 1.32713411550e-01 5.65432381171e-02 3.13971442033e-02 2.04281058891e-02 1.49536046293e-02 1.26489807288e-02
 1.15432695307e-02 1.08422207044e-02 1.08422207044e-02 1.15432695307e-02 1.26489807288e-02 1.49536046293e-02
 2.04281058891e-02 3.13971442033e-02 5.65432381171e-02 1.32713411550e-01 3.46047149862e-01 7.50033272229e-01
 1.14987597026e+00 1.21545550606e+00 9.86194056340e-01
 ...
```

The first line is a brief description.\
The second line contains NSPIN and Fermi energy.\
The following 4 lines are the informations of lattice, in order:\
&emsp;total number of atoms, the coordinate of original point.\
&emsp;the number of lattice points along lattice vector a1 (nx), a1/nx, in Bohr.\
&emsp;the number of lattice points along lattice vector a2 (ny), a2/ny, in Bohr.\
&emsp;the number of lattice points along lattice vector a3 (nz), a3/nz, in Bohr.\
The following lines are about the elements and coordinates, in order: the atom number of each atoms, the electron number in pseudopotential, the Cartesian coordinates, in Bohr.\
The rest lines are the value of charge density at each grid. Note that the inner loop is z index, followed by y index, x index in turn.\
The examples can be found in [examples/charge_density](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/charge_density)
