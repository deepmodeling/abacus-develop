# Phonopy

[Phonopy](https://phonopy.github.io/phonopy) is a powerful package to calculate phonon and related properties. The ABACUS interface has been added in Phonopy v.2.19.1. In the following, we take the FCC aluminum as an example:

1. To obtain supercells ($2\times 2\times 2$) with displacements, run phonopy:

```
phonopy -d --dim="2 2 2" --abacus
```

2. Calculate forces on atoms in the supercells with displacements. For each SCF calculation, you should specify `stru_file` with `STRU-{number}` and `cal_force=1` in INPUT in order to calculate force using ABACUS. Be careful not to relax the structures

```
echo 'stru_file ./STRU-001' >> INPUT
```

3. Then create 'FORCE_SETS' file using ABACUS inteface:

```
phonopy -f ./disp-{number}/OUT*/running*.log
```

4. Calculate the phonon dispersion:

```
phonopy band.conf --abacus
```

using the following `band.conf` file:

```
ATOM_NAME = Al
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
BAND= 1 1 1  1/2 1/2 1  3/8 3/8 3/4  0 0 0   1/2 1/2 1/2
BAND_POINTS = 21
BAND_CONNECTION = .TRUE.
```
