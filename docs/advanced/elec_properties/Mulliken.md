# Mulliken Charge Analysis

From version 2.1.0, ABACUS has the function of Mulliken population analysis. The example can be found in [examples/mulliken](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/14_mulliken). \
To use this function, set [out_mul](../input_files/input-main.md#out_mul) to `1` in the INPUT file. After calculation, there will be an output file named `mulliken.txt` in the output directory. In MD calculations, the output interval is controlled by the keyword [out_freq_ion](../input_files/input-main.md#out_freq_ion). In the file, there are contents like (`nspin 1`, H<sub>2</sub>):

```
 --- Ionic Step 1 ---
 Total charge 2
 Decomposed Mulliken population analysis for each atom
 l and m from Ylm, z stands for zeta orbital

 ------------------
 Atom 1 is H
 ------------------
          zeta     spin1
           s 1    0.9980
           s 2   -0.0052
        sum mz    0.9928

          pz 1    0.0072
          px 1    0.0000
          py 1    0.0000
         sum m    0.0072
        sum mz    0.0072

       sum lmz    1.0000

 total charge    on atom 1     1.0000
 ...
```

The ionic step index starts from 1. The file gives Mulliken charge in turn according to the order of atoms in the system. For example, the following block is for the first atom in system (`nspin 2`),

```
 ------------------
 Atom 1 is Fe1
 ------------------
          zeta     spin1     spin2       sum      diff
           s 1    0.9437    0.8381    1.7818    0.1056
...
        sum lmz    ...
```

And the next block is for the second atom in system, and so on.

```
 ------------------
 Atom 2 is Fe2
 ------------------
          zeta     spin1     spin2       sum      diff
...
```

For each atom, the file gives detailed Mulliken population analysis at different levels,

-   magnetic quantum number level: such as lines begin with 's, px, py, pz, ...'. The columns
    are the radial index (`zeta`, one-based) followed by the spin channels.
-   azimuthal quantum number level: such as lines begin with 'sum m', the sum over the
    magnetic quantum number m for one radial orbital.
-   principal quantum number level: such as lines begin with 'sum mz', the sum over both the
    magnetic quantum number m and the radial index zeta for a given orbital angular momentum.
    Here 'z' equals 'zeta' in the file, which means how many radial atomic orbitals there are
    for a given orbital angular momentum.
-   atomic level: such as lines begin with 'sum lmz' for the atom total and
    'total charge on atom' for the charge assigned to the atom.

More orbital information can be found in 'nao.txt' file output with 'mulliken.txt' when `out_mul 1`
