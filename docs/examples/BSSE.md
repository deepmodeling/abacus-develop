# BSSE

[back to main page](../../README.md)

An empty atom is defined in the `STRU` file when an element name contains the "empty" suffix, such as "H_empty" in the following example: calculating the molecular formation energy of H$_2$O by using the BSSE (Basis Set Superposition Errors) method.

$$
\Delta E(\text{H}_2\text{O}) = E(\text{H}_2\text{O}) - E(\text{O}) - E(\text{H}^1) - E(\text{H}^2)
$$

## $E(\text{H}_2\text{O})$

```
ntype   2
```
```
ATOMIC_SPECIES
H	1.008	H_ONCV_PBE-1.0.upf
O	15.9994	O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb
O_gga_6au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889725989

LATTICE_VECTORS
20 0 0
0 20 0
0 0 20

ATOMIC_POSITIONS
Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)

H
0.5
2
0.9584	0.0000	0.0000	0 0 0
-0.2392	0.9281	0.0000	0 0 0

O
0.5
1
0.0000 0.0000 0.0000 0 0 0
```

## $E(\text{O})$

```
ntype   2
```
```
ATOMIC_SPECIES
H_empty	1.008	H_ONCV_PBE-1.0.upf
O		15.9994	O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb
O_gga_6au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889725989

LATTICE_VECTORS
20 0 0
0 20 0
0 0 20

ATOMIC_POSITIONS
Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)

H_empty
0.5
2
0.9584	0.0000	0.0000	0 0 0
-0.2392	0.9281	0.0000	0 0 0

O
0.5
1
0.0000 0.0000 0.0000 0 0 0
```

## $E(\text{H}^1)$

```
ntype   3
```
```
ATOMIC_SPECIES
H			1.008	H_ONCV_PBE-1.0.upf
H_empty_2	1.008	H_ONCV_PBE-1.0.upf
O_empty		15.9994	O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb
H_gga_8au_60Ry_2s1p.orb
O_gga_6au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889725989

LATTICE_VECTORS
20 0 0
0 20 0
0 0 20

ATOMIC_POSITIONS
Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)

H
0.5
1
0.9584	0.0000	0.0000	0 0 0

H_empty_2
0.5
1
-0.2392	0.9281	0.0000	0 0 0

O_empty
0.5
1
0.0000 0.0000 0.0000 0 0 0
```

## $E(\text{H}^2)$

```
ntype   3
```
```
ATOMIC_SPECIES
H_empty_1	1.008	H_ONCV_PBE-1.0.upf
H			1.008	H_ONCV_PBE-1.0.upf
O_empty		15.9994	O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb
H_gga_8au_60Ry_2s1p.orb
O_gga_6au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889725989

LATTICE_VECTORS
20 0 0
0 20 0
0 0 20

ATOMIC_POSITIONS
Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)

H_empty_1
0.5
1
0.9584	0.0000	0.0000	0 0 0

H
0.5
1
-0.2392	0.9281	0.0000	0 0 0

O_empty
0.5
1
0.0000 0.0000 0.0000 0 0 0
```

[back to top](#BSSE)