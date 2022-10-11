# Extracting Charge Density

ABACUS can output the charge density by adding the keyword [out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg) in INPUT file:
```
out_chg             1
```
After finishing the calculation, the information of the charge density is stroed in files `OUT.${suffix}/SPIN${spin}_CHG`. At the same time, the files with the same name and with the suffix '.cube' wilL be outputed, and can be used to do visualization. \
The SPIN${spin}_CHG file looks like:
```
test
 3.17506
 1 0 0
 0 1 0
 0 0 1
 Fe1 Fe2
 1 1
Direct
 0 0 0
 0.5 0.5 0.5

  2
  0.915564 (fermi energy)
  27 27 27

 6.63608279258e-01 8.42291750502e-01 1.16305511314e+00 1.18317627770e+00 8.03596207480e-01 3.76685921988e-01 1.43138520280e-01 5.93331483316e-02
 3.22553088343e-02 2.07904703287e-02 1.50998149986e-02 1.26952816838e-02 1.15435826326e-02 1.08144612206e-02 1.08144612206e-02 1.15435826326e-02
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
The following line is the number of grid along x/y/z, and the rest lines are the value of charge density at each grid.

The examples can be found in [examples/charge_density](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/charge_density)