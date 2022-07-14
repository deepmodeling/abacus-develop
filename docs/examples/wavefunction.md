# Wave functions

[back to main page](../../README.md)


ABACUS is able to output electron wave functions in both PW and LCAO basis calculations. One only needs to do a ground-state energy calculation ***with one additional keyword in the INPUT file***:

```
out_wfc_pw              1
```
for PW basis calculation, and

```
out_wfc_lcao            1
```
for LCAO basis calculation.

In the PW basis case, the wave function is output in a file called `WAVEFUNC1.txt`. In the LCAO basis case, several `LOWF_K_#.dat` files will be output in multi-k calculation and `LOWF_GAMMA_S1.dat` in gamma-only calculation. One can also choose to output real-space wave function in PW basis calculation with the ***key word***:

```
out_wfc_r               1
```

After calculation, an additional directory named `wfc_realspace` will appear in the `OUT.$system` directory.


[back to top](#Wavefunction)