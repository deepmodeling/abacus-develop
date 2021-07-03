# Magnetic properties
[back to main page](../../README.md)

For spin polarized calculations, the users should turn out the spin option and give an initial
magnetism. The input parameter that controls spin options is:

```
nspin 2
```

The variable nspin takes values 1 or 2
- nspin=1, the default value, meaning spin-unpolarized calculation.
- nspin=2, collinear spin polarized calculation.

Initial magnetic moments are set in the STRU file, in the third line of ‘ATOMIC_POSITIONS’
part. For example,

```
...
ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 //initial magnetism
2 // number of atoms
...
```
For each element in the system, users should give their initial magnetism when nspin=2.

An example where nspin=2 is used can be found in tests/104_PW_AF_magnetic/.

[back to top](#magnetic-properties)