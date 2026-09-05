# Single-Gamma Periodic Response Implementation Plan

> Execute inline under the user's existing authorization to continue the hBN
> LCAO comparison and repair the subsequent Delta-ST workflow autonomously.

**Goal:** Admit an insulating single-Gamma periodic response at q index 1 without
changing the PBE mesh or routing through the legacy q-index-0 molecular path.

**Architecture:** Fix the response-plan guard in
`source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h`. A positive index is
an output/mesh index, not proof of a nonzero physical q vector. Permit one record
only when its reduced k coordinate is Gamma to the existing tolerance. Retain
the periodic spin/occupation checks, q metadata and self k-to-k mapping. Leave
the explicit supercell-translation compatibility flag and q-index-0 behavior
unchanged. The caller already uses the periodic independent response grid for
positive q indices, so no driver dispatch or new user parameter is needed.

**Rejected alternatives:** q index 0 bypasses the independent periodic grid;
changing k sampling would change the physical comparison. Neither is acceptable.

**Tech Stack:** C++17, existing GoogleTest smoke target, remote Intel MPI build.
No local compilation or scientific calculation.

## Steps

- [x] Trace failed hBN input through the positive-q driver to the response-plan
  guard; confirm the independent-grid construction follows the failing guard.
- [x] Add a positive-index single-Gamma self-map test and rejection tests for
  non-Gamma, fractional-occupation and spin-polarized single records.
- [x] Run the modified smoke target against the old header on df_dcu. Require
  failure of the positive single-Gamma case with the observed exception.
  Observed: 61 tests, 60 pass; the new positive-index case alone throws
  `A nonzero Sternheimer q point requires more than one k point.`
- [x] Replace only the one-record rejection condition with a Gamma-coordinate
  check; rerun the complete smoke target remotely.
  Observed: 61/61 tests pass using the corrected header and the same Intel MPI
  compiler, GoogleTest sources and test translation unit as the red run.
- [ ] Run git diff --check, inspect the patch, and create an attributed feature
  commit (Codex author, AroundPeking committer).
- [ ] Submit one guarded normal-partition build of the exact archive; require
  the seven existing targeted tests and a hashed LibRI/LibComm debug-info binary.
- [ ] After the Dojo/GTH PBE comparison selects a matched family, test a full-state
  single-frequency hBN response with the converged PBE grid and 30 Ry response
  grid. Only successful producer and numerical checks permit further frequencies
  and a reviewed master_ghj integration. This patch is not a physical RPA result.
