# ASE interface regression tests

`unit/` contains offline calculator and parser regressions. Run against the
bundled implementation after installation, or from the repository root:

```bash
PYTHONPATH=interfaces/ASE_interface python -m pytest interfaces/ASE_interface/tests/unit -q
```

The existing `abacuslite/xtest.sh` inline suites remain in place. The scripts
in this directory (`scf.py`, `md.py`, etc.) run real ABACUS calculations and
require the executable, MPI and pseudopotential/orbital assets.

## Backport provenance

Behavior and regressions are adapted from `QuantumMisaka/abacuslite` main
commit `39ab1d5595bc8456397a1d211946116cf6b15b0c`. Existing ABACUS inline
regressions, process-failure diagnostics and socket behavior are retained.
Independent package layout, publishing and cross-repository CI are excluded.
The signed-occupation parser fix and its regression come separately from
`5d3da5d45d608f329c9ad538d5d8691d30ac54c9`; this was an authored follow-up,
not part of the standalone main revision above.

The new suites cover calculator-local reader selection, SCF frame identity,
atom ordering, noncollinear moments, masses/velocities, restart file names,
point k-points, optional energy metadata and native-MD units/periodicity.
Recorded-output replay replaces only engine launching where possible; it
is not evidence that the current engine was rebuilt or run.

## SCF tolerance evidence

`test_recorded_lts_rounding_exceeds_old_tolerance` reproduces rejection of
a saved valid SCF run at `1e-5` Angstrom and acceptance at the retained
`1e-4` default. The recorded maximum cell difference is `1.65221e-5`
Angstrom. Other tests cover non-unit lattice constants, nonorthogonal
cells, periodic translations, tolerances on both sides of the threshold,
malformed structures and real displacements above the tolerance.

This is a fixed absolute identity tolerance, not a universal rounding bound.
For example, formatting `111.987654321` as six significant digits produces
`111.988`, a difference of `3.45679e-4`. Independent rounding of lattice
constants and vectors can also accumulate, and Direct-to-Cartesian conversion
propagates cell error. Further support for these cases should compare native
units or account for known field precision, rather than globally raising the
default tolerance. No JSON fallback is used. Missing and malformed abacus.json
are both covered; JSON fixtures are omitted.
