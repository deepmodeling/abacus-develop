# Stress calculation and cell relaxation

[back to main page](../../README.md)

To calculate the stress for a given structure, set input parameter ‘stress’ to 1. The stress can
be calculated in all types of ‘calculation’.
```
stress 1
```

To optimize the cell shape, one needs to setup a few parameters in the INPUT file. First, the parameter ‘calculation’ should be set to cell-relax’. The ion relaxation is automatically done in the cell shape optimization. To optimize the cell shape under external pressures, setup ‘press1, press2, press3’ which are the external pressures along the a, b, c axes in KBar. The default values is zero. You may optimize the cell shape with some of axes are fixed. To do so, use ‘fixed_axes’ as described below. The default value is None. The parameter ‘stress_thr’ is the threshold for stress.

Here is an example for the Si crystal.
```
suffix Si
calculation cell-relax
nstep 100
force_thr_ev 0.01
move_method cg
out_stru 1
trust_radius_ini 0.5
press1 0
press2 0
press3 0
fixed_axes None
stress_thr 1
```

- press1,2,3

    the external pressures along three axes in KBar, the compressive stress is taken to be positive.
- fixed_axes

    which axes are fixed when do cell relaxation. Possible choices are:
    - None : default; all can relax
    - volume : relaxation with fixed volume
    - a : fix a axis during relaxation
    - b : fix b axis during relaxation
    - c : fix c axis during relaxation
    - ab : fix both a and b axes during relaxation
    - ac : fix both a and c axes during relaxation
    - bc : fix both b and c axes during relaxation
    - abc : fix all three axes during relaxation

- stress_thr

the threshold for stress, below which the cell relaxation is considered to be converged. The unit is KBar. The default threshold is 10 KBar.

Another example for cell relaxation calculation can be found in tests/109_PW_CR/.

[back to top](#stress-calculation-and-cell-relaxation)