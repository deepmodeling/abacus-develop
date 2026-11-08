# MD / RT boundary contract

The contract test deliberately spans `MD_base::first_half`, MDCell migration and
backing-UnitCell synchronization, and RT correction of saved H/S/psi. It is
registered through CTest from the integration tree, independently of MD unit
targets. The test data are constructed in memory; no pseudopotential or full
ABACUS calculation is needed.

## Coverage moved from test_md_base.cpp (issue #8297)

| Original coverage | Destination |
| --- | --- |
| Wrapping, fractional/Cartesian positions, displacement, boundary shifts, ownership/ghost migration, repeated sync | MD unit `MatchesLegacyWrappingThroughMigrationAndRepeatedSync` |
| Bloch H/S/coefficient phases after wrapping, migration and both sync repetitions, including the following noncrossing step | Integration `MatchesLegacyWrappingThroughMigrationAndRepeatedSync` |
| Saved-field versus next-field folding, nonorthogonal cell with lat0=3, distributed/gathered H/S, untouched nonroot gathered matrices, coefficient phases and Psi cursor | Integration `HybridBoundaryMatchesFreshFoldingAtSavedFieldTime` |
| Initial and restarted force-evaluation flags/local step, initial sync state | Existing MD setup/sync unit cases |

The integration cases retain all original numerical assertions and tolerances;
the MD-owned assertions also remain in MD unit coverage. The RT helpers
`fold_hybrid`, `expect_hybrid_fold_correction`, and `expect_phase_correction`
live only in the integration test, not in MD units.

Dedicated RT units under `source/source_lcao/module_rt/test` additionally check:

- `test_boundary_fix.cpp`: direct construction of boundary shifts, Bloch and
  hybrid covariance, an independent saved-field oracle, no-shift behavior, skew
  lattice and unit conversion, distributed and gathered layouts, and a nonzero
  shared Psi cursor. It uses neither MD nor HContainer.
- `test_td_folding.cpp`: real/complex HR blocks, several translations and k
  points, zero/nonzero field, analytic combined phases, row/column-major output,
  padded leading dimensions, and accumulation into a nonzero matrix.

Each RT and integration executable has CTest variants covering 1/4 MPI ranks
and 1/4 OpenMP threads. Four ranks are essential: with one rank the distributed
matrix size equals the full size and cannot enter the distinct `len_hs` branch.

After building the corresponding targets, run:

```sh
ctest --test-dir build --output-on-failure \
  -R '^(MODULE_MD_(base|dump)_mpi|MODULE_LCAO_tddft_(boundary_fix|td_folding)|INTEGRATION_MD_RT_boundary)'
```

The legacy Makefile and pyabacus source lists need no changes: production
sources and interfaces are unchanged, and these test executables use CMake.
