# MDCell and domain decomposition

`MDCell` owns atom and cell data, the `NeighborSearch` object, and neighbor
cache data. It does not own a `DomainDecomposition`. The latter owns the MPI
Cartesian communicator and ghost exchange mappings. Both MPI and non-MPI builds
use the same ownership, migration, halo filtering, ghost update and force-return
algorithms. A non-MPI build is a single domain with rank 0, size 1, dimensions
`{1, 1, 1}` and coordinates `{0, 0, 0}`.

The implementation remains in one `domain_decomposition.cpp`. Conditional
compilation is limited to communicator management, transport wrappers and
communicator metadata access. Local exchanges copy buffers; they do not use an
alternative atom algorithm. Cartesian rank mapping assumes the existing
`reorder = 0` topology. The six-direction MPI migration transport retains its
nonblocking exchanges.

The unused `atom_overlaps_target_halo` declaration/definition and duplicate
vector math helpers were removed. Initialization accepts `CommunicationDomain`
only; the raw MPI communicator overload and its test call were removed. The
vector-based atom helpers are private. Ghost layout mutation is represented by
non-const methods, not mutable members. Force return reuses the cached slots
instead of reconstructing the communication stencil.

Compared with the former serial implementation, halo filtering stores fewer
irrelevant periodic images and coordinate-only updates retain the fixed ghost
layout. UnitCell-backed initialization now consistently wraps Cartesian input
into the primary cell in both build modes, without modifying the backing
UnitCell during initialization.

## Lifetime and initialization

The MD branch of `Driver::driver_run` creates one decomposition for the entire
run. Both `Run_MD::prepare_mdcell` overloads receive it explicitly:

- Direct force-field MD reads distributed STRU data with `MDCellReader`, which
  initializes the supplied decomposition before assigning atom ownership.
- The UnitCell-backed path calls `decomp.initialize_from_unitcell`. The backing
  UnitCell remains available for electronic-structure solvers and synchronization.

The same object is passed through `md_line`, integrator `setup`, and
`MD_func::force_virial`. Neither `MD_base` nor the integrators store a decomposition
reference; their constructors remain unchanged. Do not create a temporary
decomposition for each exchange: ghost layout mappings must survive until the
next position update and force return. Use one decomposition per active MDCell.

## Direct MDCell force evaluation

The ESolver interface is unchanged and does not accept a decomposition.
Solver initialization calls `MDCell::set_neighbor_cutoff`, which only sets the
cutoff and invalidates neighbor data; it does not communicate.

`MD_func::force_virial` performs the following direct-MDCell sequence:

1. `decomp.prepare_neighbors(mdcell)` migrates and builds ghosts/neighbors when
   necessary, or updates ghost positions and refreshes active neighbors.
2. `solver.runner(mdcell, istep)` computes energy and owned/ghost forces.
3. `decomp.accumulate_ghost_forces(mdcell)` returns and consumes ghost forces.
4. Convert the resulting owned forces, energy, and stress from Ry to Hartree.

DP and NEP convert both owned and ghost forces to the same Ry-based units before
returning from `runner`; DP force rescaling likewise covers both containers.
Calling a direct MDCell solver outside the MD driver requires explicit neighbor
preparation and force return by that caller. The UnitCell-backed solver path
continues to use `runner`, `cal_force`, and `cal_stress` on the backing UnitCell.

The decomposition synchronizes cached lattice/cutoff/skin values before neighbor
preparation. Rebuild decisions are collective, including local cache invalidation
and maximum displacement. A fixed ghost layout must retain its periodic images
during coordinate-only updates; the incoming image shift is the inverse of the
local-to-source-domain shift.

## Regression coverage

`test_domain_decomposition.cpp` is compiled in serial (`MODULE_MD_LJ_pot`) and
MPI (`MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi`) configurations. It covers cutoff
assignment without exchange, migration, force return without double accumulation,
layout reuse, periodic coordinate updates, cell/cutoff changes, single-rank
invalidation, skin displacement, and move construction/assignment. Additional
tests check single-domain halo filtering, cached slots across a halo edge, and
skew-cell multilayer neighbors against brute-force image enumeration. MPI cases
run with one, two and four processes. `test_md_func.cpp` checks that force return
precedes unit conversion and repeated evaluations do not retain ghost forces.

Local verification on 2026-09-10:

```bash
cmake --build build_parallel_md1_current --target \
  MODULE_MD_LJ_pot MODULE_MD_func MODULE_MD_run MODULE_MD_fire \
  MODULE_MD_verlet MODULE_MD_nhc MODULE_MD_msst MODULE_MD_lgv \
  MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi MODULE_CELL_NEIGHBOR_mdcell_reader -j12
OMP_NUM_THREADS=1 \
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openblas-pthread:/home/fy19/abacus-develop/toolchain/install/openmpi-5.0.7/lib \
ctest --test-dir build_parallel_md1_current --output-on-failure \
  -R '^(MODULE_MD_|MODULE_CELL_NEIGHBOR_mdcell_|MODULE_CELL_NEIGHBOR_domain_decomposition_)'
cmake --build build_abacus_nep_gnu -j12
```

The focused CTest selection passed 11/11 targets. The full MPI build with DP and
NEP enabled passed. `abacus_basic_para --version` reported `v3.11.0-beta9`.
Its generated commit banner still identified `bfc4d1689` from the earlier CMake
configuration; the binary was rebuilt with the subsequent implementation and
periodic-image fix.

Copied `tests/04_FF/50_DP_Al` (MSST) and `tests/04_FF/101_NEP_HfO2` (NPT) into
isolated directories, replacing only PP_ORB/model paths with absolute paths.
Each four-step case exited successfully at both one and two MPI processes:

```bash
OMP_NUM_THREADS=1 LD_LIBRARY_PATH="$md_runtime_libs" \
  /home/fy19/abacus-develop/toolchain/install/openmpi-5.0.10/bin/mpirun \
  -np 1 /home/fy19/abacus-develop/build_abacus_nep_gnu/abacus_basic_para
# Repeat with -np 2 in a separate copy of each case.
```

Here `md_runtime_libs` contained the existing `deepmd_kit.libs`, TensorFlow,
DeePMD, OpenMPI 5.0.10 and OpenBLAS 0.3.33 library directories. The local NEP
library lacks a SONAME and is linked by a relative path; the isolated test tree
provided an `external` symlink to the existing repository dependency directory.
Initial launches without these runtime paths failed in the loader, before any
ABACUS calculation. No dependency or reference files were changed.

After sorting each MD_dump frame by atom ID, the largest absolute difference
across positions, forces, velocities, lattice and virial entries was approximately
`1e-12` for DP and `2.004e-12` for NEP. Energy log entries agreed at printed
precision. These are process-count consistency checks, not new reference values.
No full electronic-structure MD integration case was run.

This refactor changes no INPUT parameter behavior, so generated INPUT metadata
and user parameter documentation do not require changes. Remaining MD globals
are historical output/restart dependencies; their diff-level count does not grow.
The explicit MPI include in `mdcell.h` is still needed for its communicator
metadata, while test fixtures owning a decomposition by value need its full type.

## Unified-algorithm verification, 2026-09-11

The build commands above passed again after unifying the implementation. The
CTest selection now passes 12/12 targets, including the new MPI single-process
target. The shared DomainDecomposition suite has ten tests and runs in both
non-MPI and MPI builds; the UnitCell preparation test additionally checks
out-of-cell coordinates, retained velocity/mobility and the backing-cell link.

The four-step DP/MSST and NEP/NPT cases were rerun in separate `unified_*`
directories using the same runtime commands and library paths above. All four
runs exited successfully. After sorting each MD_dump frame by atom ID, the
unified implementation's output matched the preceding implementation exactly at
printed precision for each process count. The one-versus-two-process maximum
absolute differences remained approximately `1e-12` (DP) and `2.004e-12` (NEP).
`git diff --check` passed. Changes remain uncommitted for VSCode review.
