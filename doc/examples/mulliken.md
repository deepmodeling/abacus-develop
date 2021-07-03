# Mulliken charge

[back to main page](../../README.md)

From version 2.1.0, ABACUS has the function of Mulliken population analysis. To use this function, set ‘mulliken’ to ‘1’ in the INPUT file. After calculation, there will be an output file named mulliken.txt in the output directory. In the file, there are contents like:

```

 CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM
 8 (Total charge all spin 1)
 8 (Total charge of the system)
  Decomposed Mulliken populations
 0                       Si                       Up spin                     Down spin                           Sum                          Diff
                     multiple
 s                       0                      0.51810998                    0.51810998                       1.03622                           0
  sum over m                                  0.51810998                    0.51810998                            1.03622                        0
 s                       1                     0.053661483                   0.053661483                    0.10732297                           0
  sum over m                                 0.053661483                   0.053661483                         0.10732297                        0
  sum over m+mul                           0.57177146                    0.57177146                        1.1435429                            0
 px                      0                      0.43730519                    0.43730519                    0.87461038                           0
 py                      0                      0.43730519                    0.43730519                    0.87461038                           0
 pz                      0                      0.43730519                    0.43730519                    0.87461038                           0
  sum over m                                   1.3119156                     1.3119156                          2.6238311                        0
 px                      1                    0.0065031719                  0.0065031719                   0.013006344                           0
 py                      1                    0.0065031719                  0.0065031719                   0.013006344                           0
 pz                      1                    0.0065031719                  0.0065031719                   0.013006344                           0
  sum over m                                 0.019509516                   0.019509516                        0.039019031                        0
  sum over m+mul                            1.3314251                     1.3314251                        2.6628502                            0
 d3z^2-r^2               0                     0.011750855                   0.011750855                   0.023501711                           0
 dxy                     0                     0.024433913                   0.024433913                   0.048867826                           0
 dxz                     0                     0.024433913                   0.024433913                   0.048867826                           0
 dx^2-y^2                0                     0.011750855                   0.011750855                   0.023501711                           0
 dyz                     0                     0.024433913                   0.024433913                   0.048867826                           0
  sum over m                                 0.096803451                   0.096803451                          0.1936069                        0
  sum over m+mul                          0.096803451                   0.096803451                        0.1936069                            0
 Total Charge on atom  Si                   4
 ...
```
The file gives Mulliken charge in turn according to the order of atoms in the system. For example, the following block is for the first atom in system,
```
0 Si Up spin Down spin Sum Diff
...
Total Charge on atom Si 4.0241712
```

And the next block is for the second atom in system, and so on.
```
1 Si Up spin Down spin Sum Diff
...
```
For each atom, the file gives detailed Mulliken population analysis at different levels,
- magnetic quantum number level: such as lines beigin with ‘s,px,py,pz,...’
- azimuthal quantum number level: such as lines begin with ‘sum over m’.
- principal quantum number level: such as lines begin with ‘sum over m+mul’. Here ‘mul’
equals ‘multiple’ in the file, which means how many radial atomic orbitals there are for a given orbital angular momentum.
- atomic level: such as lines begin with ‘Total Charge on atom’.

[back to top](#mulliken-charge)