# SOC symmetry export regression design

## Goal

Protect the ABACUS-to-LibRPA `stru_out` symmetry export introduced by PR #4 and
the subsequent format correction. The regression must cover ordinary spatial
operations and antiunitary magnetic operations without requiring a full
electronic-structure calculation.

## Scope and non-goals

The change covers only the ABACUS producer-side export format. It does not
enable or validate LibRPA SOC symmetry itself; `master_ghj` does not yet contain
the corresponding LibRPA reader change. A separate remote BN SOC calculation
remains the integration gate for energy equivalence and symmetry speedup.

## Design

Extract the small symmetry-row formatter currently embedded in
`RPA_LRI<T,Tdata>::out_struc` into a focused header/helper in
`source/source_lcao/module_ri`. The helper will accept the spatial rotation and
translation components plus the ordinary/antiunitary counts, and will write the
common `N row` block:

1. one count line containing ordinary plus antiunitary operations;
2. one row per operation, with nine integer rotation entries followed by three
   scientific-notation translation entries;
3. ordinary operations first, then antiunitary spatial parts;
4. no `spin_symmetry` trailer or other spin-specific section.

`out_struc` will retain the existing symmetry-enable and magnetic-count gating,
and will call the helper so production output and the unit test share the same
implementation.

## Regression test

Add `MODULE_RI_librpa_stru_symmetry_test` beside the existing LibRPA interface
unit tests. The test will construct a small deterministic set (two ordinary
operations and one antiunitary operation), serialize it through the helper, and
assert:

- the count is `3 row`;
- all three rows are present in the required order;
- every row has exactly 12 numeric fields;
- the rotation fields are integral and translations preserve scientific output;
- the serialized block contains no `spin_symmetry` marker.

The test is standalone, deterministic, and has no PP, orbital, MPI, or runtime
data dependency. It should complete in seconds and be registered in the same
CMake test list as `librpa_stru_units_test`.

## Validation and acceptance

Because local ABACUS compilation is prohibited, validation will be performed on
`df_iopcas_ghj`:

- build and run the new unit test against the modified source;
- run the existing relevant `module_ri` interface tests;
- run the short BN SOC producer comparison before/after the merge stack, with
  `symmetry=1` and `symmetry=-1` energies agreeing within `1e-6 eV` and the
  symmetry-enabled case using fewer irreducible k points.

Only after these gates pass will the branch be rebased onto the current
`deepmodeling/abacus-develop` branch. Rebase conflicts will preserve the newer
upstream naming and separated EXX/RPA interfaces; no LibRPA SOC reader changes
will be invented in this branch.
