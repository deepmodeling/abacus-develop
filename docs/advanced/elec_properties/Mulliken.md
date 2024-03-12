# Mulliken Charge Analysis

From version 2.1.0, ABACUS has the function of Mulliken population analysis. The example can be found in [examples/mulliken](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/mulliken). \
To use this function, set [out_mul](./input-main.md#out_mul) to `1` in the INPUT file. After calculation, there will be an output file named `mulliken.txt` in the output directory. In MD calculations, the output interval is controlled by the keyword [out_interval](./input-main.md#out_interval). In the file, there are contents like (`nspin 1`):

```
STEP: 0
CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM
 Total charge of spin 1:	8
 Total charge:	8
Decomposed Mulliken populations
0                 Zeta of Si                        Spin 1
s                        0                       1.2553358
  sum over m                                     1.2553358
s                        1                    -0.030782972
  sum over m                                  -0.030782972
  sum over m+zeta                                1.2245529
pz                        0                      0.85945806
px                        0                      0.85945806
py                        0                      0.85945806
  sum over m                                     2.5783742
pz                        1                    0.0065801228
px                        1                    0.0065801228
py                        1                    0.0065801228
  sum over m                                   0.019740368
  sum over m+zeta                                2.5981145
dz^2                        0                       0.0189287
dxz                        0                     0.046491729
dyz                        0                     0.046491729
dx^2-y^2                        0                       0.0189287
dxy                        0                     0.046491729
  sum over m                                    0.17733259
  sum over m+zeta                               0.17733259
Total Charge on atom:  Si                   4
 ...
```

The file gives Mulliken charge in turn according to the order of atoms in the system. For example, the following block is for the first atom in system (`nspin 2`),

```
0            Zeta of Si               Spin 1              Spin 2                Sum                Diff
...
Total Charge on atom:  Si                   4
Total Magnetism on atom:  Si      -1.2739809e-14
```

And the next block is for the second atom in system, and so on.

```
1            Zeta of Si               Spin 1              Spin 2                Sum                Diff
...
```

For each atom, the file gives detailed Mulliken population analysis at different levels,

-   magnetic quantum number level: such as lines beigin with ‘s,px,py,pz,...’
-   azimuthal quantum number level: such as lines begin with ‘sum over m’.
-   principal quantum number level: such as lines begin with ‘sum over m+zeta’. Here ‘zeta’
    equals ‘zeta’ in the file, which means how many radial atomic orbitals there are for a given orbital angular momentum.
-   atomic level: such as lines begin with ‘Total Charge on atom’.

More orbital information can be found in 'Orbital' file output with 'mulliken.txt' when `out_mul 1`
