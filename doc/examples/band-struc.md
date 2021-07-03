# Band structure

[back to main page](../../README.md)

This example shows how to calculate the energy band structure. Similar to the [DOS case](#dos.md), we first, do a ground-state energy calculation as in [this example](#basic-lcao.md) ***with one additional keyword in the INPUT file***:

```
out_charge              1
```

this will produce the converged charge density, which is contained in the file SPIN1_CHG. Copy the file along with the `STRU` file, the pseudopotential file and the atomic orbital file to the new working directory where we will do a non-self-consistent calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:
```
INPUT_PARAMETERS
#Parameters (General)
ntype 1
nbands 8
calculation nscf
basis_type lcao
read_file_dir   ./

#Parameters (Accuracy)
ecutwfc 60
niter 50
dr2 1.0e-9
ethr 1.0e-7

#Parameters (File)
start_charge file
out_band 1

#Parameters (Smearing)
smearing gaussian
sigma 0.02
```

Here the the relevant k-point file KPT looks like,
```
K_POINTS # keyword for start
6 # number of high symmetry lines
Line # line-mode
0.5 0.0 0.5 20 # X
0.0 0.0 0.0 20 # G
0.5 0.5 0.5 20 # L
0.5 0.25 0.75 20 # W
0.375 0.375 0.75 20 # K
0.0 0.0 0.0 1 # G
```

This means we are using:

- 6 number of k points, here means 6 k points:
(0.5, 0.0, 0.5) (0.0, 0.0, 0.0) (0.5, 0.5, 0.5) (0.5, 0.25, 0.75) (0.375, 0.375, 0.75) (0.0, 0.0,
0.0)

- 20/1 number of k points along the segment line, which is constructed by two adjacent k
points.

Run the program, and you will see a file named BANDS_1.dat in the output directory. Plot it
to get energy band structure.

[back to top](#band-structure)