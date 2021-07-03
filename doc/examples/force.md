# Force calculation and structure relaxation
[back to main page](../../README.md)

To calculate the atomic forces for a given structure without ion relaxation, set ‘calculation’ to ‘scf’, set input parameter ‘force’ to 1.

```
calculation scf
force 1
```

To relax the atom position without change cell shape, one needs to add a few more parameters
in the INPUT file. Here is an example for the Si dimer. In this case, the forces are calculated by
default.

```
calculation relax
gamma_only 1
nstep 100
force_thr_ev 0.01
move_method cg
out_stru 1
trust_radius_ini 0.5
```

- `calculation` relax

    relax atom positions with fixed lattice vectors.
- `nstep`

    the maximal number of ionic iteration steps.
- `force_thr_ev`

    the threshold for the force, below which the geometry relaxation is considered to be converged. The unit is eV/Angstrom.
- `move_method`

    the algorithm used for geometry optimization. Possible choices are:
- `cg`

    conjugate gradient (CG) algorithm

- `bfgs`

    Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm

- `sd`

    steepest descend algorithm. The CG method is recommended.

- `out_stru`

    output the structure of each step or not.

- `trust_radius_ini`

    the initial radius of the relaxation. We advise you not to change this parameter, unless you are sure that the initial structure is close to the final structure.

An example for structure relaxation calculation can be found at tests/108_PW_RE.

[back to top](#force-calculation-and-structure-relaxation)