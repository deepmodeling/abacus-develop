# STRU file
- [Example](#example)
- [Structure of the file](#structure-of-the-file)

    [back to main page](../README.md)
## Example

The `STRU` file contains the information about the name(s) and/or location(s) of the pseudopotential
and numerical orbital files, as well as the structural information about the system. Below is an example of the `STRU` file:
```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file

LATTICE_CONSTANT
10.2 // lattice scaling factor (Bohr)

LATTICE_VECTORS
0.5 0.5 0.0 // latvec1
0.5 0.0 0.5 // latvec2
0.0 0.5 0.5 // latvec3

ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 // magnetism
2 // number of atoms
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```

[back to top](#stru-file)

## Structure of the file
The `STRU` file contains several sections, and each section must start with a keyword like
`ATOMIC_SPECIES`, `NUMERICAL_ORBITAL`, or `LATTICE_CONSTANT`, etc. to signify what type of
information that comes below.

- ATOMIC_SPECIES

    This section provides information about the type of chemical elements contained the unit cell. Each line defines one type of element. The user should specify the name, the mass, and the pseudopotential file used for each element. The mass of the elment is only used in molecular dynamics simulations. For electronic-structure calculations, the actual mass value isn’t important.
    In the above example, we see information is provided for the element `Si`:
    ```
    Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file
    ```
    Here `Si_ONCV_PBE-1.0.upf` is the pseudopotential file. When the path is not specified, the file is assumed to be located in work directory. Otherwise, please explicitly specify the location of the pseudopotential files.
-   NUMERICAL_ORBITAL

    Numerical atomic orbitals are only needed for `LCAO` calculations. Thus this section will be neglected in calcultions with plane wave basis. In the above example, numerical atomic orbitals is specified for the element `Si`:
    ```
    Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
    ```
    ‘Si_gga_8au_60Ry_2s2p1d.orb’ is name of the numerical orbital file. Again here the path is not specified, which means that this file is located in the work directory.
- LATTICE_CONSTANT

    The lattice constant of the system in unit of Bohr.
- LATTICE_VECTORS

    The lattice vectors of the unit cell. It is a 3by3 matrix written in 3 lines. Please note that *the lattice vectors given here are scaled by the lattice constant*.
- ATOMIC_POSITIONS
    
    This section specifies the positions and other information of individual atoms. The first line signifies whether atom positions are given in `Cartesian` or `Direct` coordinates.
    
    The following three lines tells the elemental type (`Si`), the initial magnetic moment (`0.0`), and the number of atoms for this particular element (`2`) repsectively.
    
    The last two lines in this example are the coordinates of atomic positions. There are six numbers in each line: the first three specifies the atomic positions and the last three control how the atom move in geometry relaxation calculations. The numbers `0 0 0` following the coordinates of the first atom means this atom are *not allowed* to move in all three directions, and the numbers `1 1 1` following the coordinates of the second atom means this atom *can* move in all three directions.

[back to top](#stru-file)