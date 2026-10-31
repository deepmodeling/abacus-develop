# DeltaSpin performance benchmarks

The timing benchmark matrix uses 16-atom and 32-atom bcc Fe supercells, one GPU,
one MPI rank, and either one or twelve CPU threads. Each thread configuration
compares the original and optimized binaries using identical input and asset
hashes. Thread-count effects are reported separately from implementation speedup.
These fixed fixtures measure regressions in a general-purpose DFT implementation;
their element, atom counts, and coordinates do not define optimization dispatch.
Additional disordered and mixed-element cases test applicability independently.

The reference source is tag `delta-spin-yinan`, commit
`19ade18599bfdf9ee78a289bf798a9b2583d52cd`. Keep its executable separate from
the development executable. Historical two-atom measurements are not evidence
for performance on these supercells.

## Inputs

The structures use a 2.866 A conventional bcc lattice. The formal timing cases
repeat it 2x2x2 for Fe16 and 2x2x4 for Fe32. Capacity inputs use 3x3x3 for Fe54,
4x4x2 for Fe64, 4x3x3 for Fe72, 5x4x2 for Fe80 and 4x4x3 for Fe96. One atom is
displaced by 0.02 A along x to provide nonzero forces for comparison. The local
moment target is 2.2 uB along +z on every atom.

These Fe cases use 100 Ry, double precision, noncollinear spin, spin-orbit coupling,
the same pseudopotential and orbital file, and an explicit Gamma-centered
2x2x2 k mesh. `kspacing=0` selects the supplied KPT file. Electronic and magnetic
thresholds, smearing and mixing parameters are recorded in each INPUT file.
These are numerical-equivalence and performance fixtures, not a convergence
study of bulk Fe properties. The structures are not relaxed during a benchmark.

Copy the two supplied data files into `assets/` and verify them against
`assets.sha256.json` before running:

```text
Fe.upf
Fe_gga_8au_200.0Ry_4s2p2d1f.orb
```

The reference Fe54 two-step probe reaches approximately 32,343 MiB on a
32 GiB V100. Its unrestricted run fails with CUDA out-of-memory after two
completed SCF steps. Fe54 is retained as a capacity test; it cannot provide
a converged reference runtime on this node. Fe32 supplies the larger-system
timing comparison without changing the physical or solver settings.

### Generality controls

`fe16-disordered` retains the formal Fe16 cell, INPUT, KPT and assets, with every
atom independently displaced by at most 0.03 A per Cartesian component. No atom
lies on a supported exact binary-coordinate phase grid. This tests the dense
projector path and separates general implementation gains from periodic-phase
factorization. The input is fixed; it is not adjusted to obtain faster convergence.

`feni4-nonorthogonal` contains two Fe and two Ni atoms in a sheared cell with
independent noncollinear moment targets. It uses the repository's scalar
relativistic Ni pseudopotential, so spin-orbit coupling is disabled for both
executables in this case. Precision, cutoff, k mesh, solver, mixing and stopping
controls otherwise match the formal fixtures. This is a numerical control,
not an equilibrium alloy or a validated materials model.

`feni8-nonorthogonal` repeats the four-atom mixed cell twice along its first
lattice vector, retains each atom's moment target and uses `pw_seed=1`. Its 164
bands exercise the large-product path with multiple species. Dense-k controls
copy `feni4-nonorthogonal`, use the same positive seed and change only KPT to
Gamma-centered 3x3x3 or 5x5x5 meshes. These controls are separate from the
unchanged formal Fe16/Fe32 timing inputs. A positive seed makes the random
initialization depend on global k-point indices; `pw_seed=0` does not guarantee
the same initial orbitals when the MPI partition changes.

Run these cases in a separate working directory so adding Ni assets does not
change the recorded asset set of earlier formal runs. From this directory:

```bash
export DELTASPIN_GENERALITY=/path/to/generality-workspace
mkdir -p "$DELTASPIN_GENERALITY"/{assets,inputs,src,results}
cp src/benchmark.py src/compare_results.py "$DELTASPIN_GENERALITY/src/"
cp -r inputs/fe16-disordered inputs/feni4-nonorthogonal inputs/feni8-nonorthogonal \
    "$DELTASPIN_GENERALITY/inputs/"
cp assets/Fe.upf assets/Fe_gga_8au_200.0Ry_4s2p2d1f.orb "$DELTASPIN_GENERALITY/assets/"
cp ../../../tests/PP_ORB/Ni_ONCV_PBE-1.0.upf \
   ../../../tests/PP_ORB/Ni_gga_9au_100Ry_4s2p2d1f.orb "$DELTASPIN_GENERALITY/assets/"
cp generality-assets.sha256.json "$DELTASPIN_GENERALITY/assets.sha256.json"
```

Verify the copied assets against this manifest, then use the same original and
candidate commands below with these case names and twelve threads. Each pair
must use identical inputs and assets. Report native electronic convergence and
per-atom physical differences separately from timing. These controls broaden
the evidence but do not establish speedups for all structures, functionals,
pseudopotentials, solvers, device types or MPI decompositions.

## Run the matrix

Load the installation environment so MPI, CUDA and the linked numerical
libraries are available. Run from this directory. Set the binary paths to the
two preserved installations:

```bash
export DELTASPIN_BASELINE=/path/to/abacus-baseline
export DELTASPIN_CANDIDATE=/path/to/abacus-optimized
mkdir -p results
for system in fe-bcc-16 fe-bcc-32; do
    for threads in 1 12; do
        python3 src/benchmark.py --binary "$DELTASPIN_BASELINE" \
            --case "$system" --label "baseline-${system}-t${threads}" \
            --threads "$threads"
        python3 src/benchmark.py --binary "$DELTASPIN_CANDIDATE" \
            --case "$system" --label "candidate-${system}-t${threads}" \
            --threads "$threads"
        python3 src/compare_results.py \
            --reference "baseline-${system}-t${threads}" \
            --candidate "candidate-${system}-t${threads}" \
            --output "results/comparison-${system}-t${threads}.json"
    done
done
```

Each run creates a new label directory containing copied inputs, asset links,
metadata, application logs, 100 ms GPU-memory samples and a parsed result.
An existing label is rejected. `OMP_NUM_THREADS` is set to the requested count;
`OPENBLAS_NUM_THREADS` is set to one. OpenMP builds of OpenBLAS follow their
OpenMP runtime, so verify the linked BLAS threading implementation when moving
to another installation.

Use `--steps 2` for a bounded startup/memory probe. Such a run does not establish
electronic convergence. Use a new label for every repeat or candidate. Run one
benchmark at a time on the GPU and avoid simultaneous builds on the measured
host. `--profile` is available for Nsight Systems diagnostics; profiled runs are
excluded from the normal comparison helper.

## Multi-GPU runs

`run_multigpu.py` runs a single-node comparison with one MPI rank per requested
physical NVIDIA GPU. It binds each rank to a distinct CPU allocation and masks
its device by UUID. Requested GPUs must be idle at launch; the result requires
one distinct observed compute process per GPU. Driver PIDs and namespace-local MPI
PIDs are recorded separately. Runs retain effective input/asset hashes, process
placement, per-device memory peaks and the maximum sum within complete sampling
rounds. The aggregate samples are not simultaneous atomic measurements.

```bash
python3 src/run_multigpu.py --binary "$DELTASPIN_CANDIDATE" \
    --case fe-bcc-16 --label candidate-fe16-g4-t12-full \
    --devices 0,1,2,3 --threads 12 --cores-per-rank 6
```

Repeat at one, two and four GPUs with the same threads per rank and scientific
inputs. Both CPU and GPU allocation grow in that series. The GPU plane-wave
implementation assigns one rank to each k-point pool, so the number of ranks
must not exceed the irreducible k-point count. Multi-GPU decomposition within one
k point is outside this implementation. `--workspace` selects independent
inputs/assets; `--steps` creates a bounded diagnostic and is not a converged
performance result. The runner records native convergence separately from a
complete executable run. `summarize_multigpu.py` reparses raw outputs and sampling
records before accepting a converged scaling series.

### Capacity and allocation scaling

The larger Fe inputs retain the Fe54 INPUT and KPT settings. Record successful
completion of a bounded probe separately from full electronic convergence;
neither a successful first allocation nor a few completed steps establish full
capacity. Allocation failures remain useful evidence of the memory boundary.

K-point distribution reduces the number of resident wavefunction blocks per
rank. It does not partition the Davidson arrays for an individual k point. If
`B = sizeof(complex<double>) * nbands * nbasis`, `q` is the number of local k
points and `d = pw_diag_ndim`, the resident Psi and two Davidson basis arrays
alone require approximately `(q + 2*d) * B` bytes. Projectors, dense matrices,
CUDA context and temporary allocations add further costs. At fixed density,
cutoff and bands per atom, both `nbands` and `nbasis` grow approximately with
atom count. Four GPUs therefore do not imply four times as many atoms under
this decomposition. Measure the actual allocation boundary for each input.

The Fe64 four-V100 input is sensitive to floating-point SCF trajectories at
twelve CPU threads per rank. Two complete runs of the unchanged parent
`3ff7cb39` both converge in 39 steps but differ by up to `3.03e-4 eV/angstrom`
in their final forces. The direct-upload candidate differs from the first
parent by `3.55e-4 eV/angstrom` and from the second by `7.24e-5 eV/angstrom`.
These observations document numerical sensitivity; two repeats do not establish
a statistical bound. Keep this capacity case separate from strict physical
regression controls, and retain all differences even when native SCF converges.

Complete four-V100 comparisons of direct Hpsi uploads against that parent use
the same CUDA 12.8 runtime, twelve threads and six physical cores per rank:

| Case | Parent seconds | Direct-upload seconds | Parent peak MiB per GPU | Direct-upload peak MiB per GPU | SCF steps, both |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fe16 | 31.783 | 32.407 | 2468 | 2252 | 30 |
| Fe32 | 130.926 | 133.996 | 8056 | 7248 | 34 |
| Fe54 | 441.392 | 450.511 | 20662 | 18368 | 31 |
| Fe64, sensitivity control | 754.566-755.228 | 772.850 | 29236 | 26012 | 39 |

The Fe54 maximum force difference is `2.19e-8 eV/angstrom`. Fe64 sampled total
device memory falls from 116620-116736 MiB across its two parent runs to
103812 MiB. These measurements support a memory/time tradeoff, not a speedup.
The Fe64 numerical limitation above applies independently of its memory saving.

The matching single-V100 controls retain six physical cores and vary only the
requested OpenMP thread count within each case:

| Case | Threads | Parent seconds | Direct-upload seconds | Parent / direct peak MiB | SCF steps, both |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fe16 | 1 | 79.857 | 84.043 | 3722 / 3518 | 30 |
| Fe16 | 12 | 74.329 | 78.947 | 3722 / 3518 | 30 |
| Fe32 | 1 | 396.590 | 414.512 | 13030 / 12222 | 34 |
| Fe32 | 12 | 379.884 | 396.923 | 13030 / 12222 | 34 |

These controls reduce peak memory by 5.5-6.2% at 4.5-6.2% greater elapsed time.
Both one-thread pairs preserve all printed energies, forces, moments and
magnetic forces; stress differences are at most `1.01e-10 kbar`. Across all four
pairs, the largest force difference is `1.66e-8 eV/angstrom`.

Thirteen complete direct-upload regression pairs cover V100 and H20,
one/four GPUs, one/twelve threads, disordered Fe, mixed Fe/Ni and an uneven
27-point distribution. Their largest force/stress differences are
`2.19e-8 eV/angstrom` and `4.33e-7 kbar`; every pair retains its SCF step count.
The Fe64 sensitivity case and Fe72 capacity run are excluded from these maxima.

The direct-upload executable completes Fe72 on four V100s in 33 SCF steps and
1067.371 seconds, with a maximum per-device peak of 32026 MiB and a sampled total
of 127636 MiB. This is complete-run capacity evidence, not a paired speedup or
physical-accuracy comparison: the parent Fe72 run was limited to three steps.
The optional prefetch already failed to allocate at this size, so removing it
does not substantially change the mandatory allocation peak.

The Fe80 and Fe96 parent probes exhaust device memory before the first SCF step.
An isolated host-staging prototype reduces the Fe72 three-step peak by about
13% at about 12% greater elapsed time, but also fails at Fe80. It is excluded
from the production implementation. These observations locate a tested
capacity boundary; they do not establish a universal maximum atom count.

For MPI profiles, pass one encoded Nsight `globalPid` at a time to
`analyze_trace.py`, `analyze_gemm_trace.py` and `analyze_idle.py` with
`--global-pid` and distinct `--output` paths. Unfiltered multi-process application
traces are rejected. `analyze_mpi.py` aligns world-communicator collectives by
call order and traced payloads, separating arrival spread from time after the
last rank enters. Those quantities include software progress and are not pure
network latency or a prediction of recoverable wall time. Analysis outputs use
new paths and are not overwritten.

## Profile device work

`build_profile.py` compiles isolated NVTX-instrumented timer and GEMM objects
and relinks the active build without editing production sources. The optional
GEMM ranges record transpose modes, dimensions and leading dimensions. The
label must be unused. With Nsight Systems available in `PATH`:

```bash
python3 src/build_profile.py gemm-shapes --gemm-shapes
python3 src/benchmark.py --binary results/abacus-gemm-shapes-profile \
    --case fe-bcc-32 --label gemm-shapes-fe32-full --threads 12 --profile
nsys export --type sqlite \
    --output results/gemm-shapes-fe32-full/profile.sqlite \
    results/gemm-shapes-fe32-full/profile.nsys-rep
python3 src/analyze_trace.py results/gemm-shapes-fe32-full/profile.sqlite
python3 src/analyze_gemm_trace.py results/gemm-shapes-fe32-full/profile.sqlite
```

Both analyzers assign device activity through CUDA launch correlation, so a
kernel that executes after its host NVTX range still belongs to the correct
operation. GPU activity unions avoid double-counting concurrent execution;
inclusive range durations overlap and cannot be summed. The steady view excludes
the first two `hamilt2density` calls and the final force/stress calculation.

The GEMM analyzer records the selected multiplication method and counts nominal
multiply work as `8*m*n*k` for ZGEMM, `6*m*n*k` for ZGEMM3M, or
`8*m*n*k*batches` for strided batches. Calls count wrapper invocations; a batch
count multiplies arithmetic work, not the measured kernel duration. These shape-derived
estimates exclude packing, recombination, alpha/beta operations and products
internal to cuSolver; they are not hardware-counter measurements. Check numerical
outputs, iteration counts and input hashes against an unprofiled run before
using a trace to interpret the accepted implementation. Profiling wall time
includes collection and report overhead and is not a speed benchmark.

## Interpret results

`benchmark.py` reads the atom count from the application output and checks the
complete force, moment, magnetic-force and target-moment tables. A missing or
truncated table does not count as a complete run. Energy, stress, native SCF
convergence, density error and magnetic RMS are recorded separately.

`compare_results.py` requires equal systems, thread/rank counts, step overrides,
input hashes, asset hashes and table sizes. It reports elapsed-time ratios and
absolute physical-output differences without imposing a universal tolerance.
Use repeat runs of the reference to establish numerical variability. A speedup
claim requires complete converged runs and accepted physical comparisons;
native electronic convergence alone does not prove the magnetic target is met.

Report both the sampled peak memory and its sampling interval. A 100 ms sample
can miss a shorter allocation peak, so retain the raw samples and use the same
measurement method for both binaries.

## Fused spinor transforms

The optional `USE_CUFFTDX` build uses one common kernel template for the local
spinor potential and density. CMake enumerates FFT lengths from 16 through 256
whose prime factors are 2, 3 and 5. Each length and CUDA architecture is
instantiated in a separate translation unit. Lengths that require a cuFFTDx workspace, and lengths outside
this registry, retain the cuFFT path. The optional build requires CUDA and
cuFFTDx 1.6 or newer, discoverable through the MathDx CMake package.

Descriptors and CUDA images follow `CMAKE_CUDA_ARCHITECTURES`, including its
`native`, `all`, `all-major`, `-real` and `-virtual` forms. The runtime selects a
loadable descriptor and checks kernel thread and shared-memory requirements on
the current device. The selection is cached per host thread and invalidated on
device switches. Missing images or insufficient resources retain cuFFT. There
is no GPU-model whitelist. A configuration targeting Volta, Ampere and Hopper
can use the following options in addition to the installation's dependency settings.
Run the build and validation commands in this section from the repository root:

```bash
cmake -S . -B build -DUSE_CUDA=ON -DUSE_CUFFTDX=ON \
    '-DCMAKE_CUDA_ARCHITECTURES=70;80;90' -Dmathdx_ROOT=/path/to/mathdx
cmake --build build --parallel 6
```

`USE_CUFFTDX=OFF` retains ordinary cuFFT execution. Compilation for an
architecture establishes build coverage; hardware-specific performance requires
matched complete calculations on that device.

The application selects the fused path for double-complex spinors with one MPI
rank per k-point pool and an eligible plane-wave layout. Multiple GPU pools
transform their own k points independently. The x, y and z dimensions
are selected independently. Rectangular and nonorthogonal cells use the same
implementation. A missing y template can
use cuFFT for the y-z transforms while retaining the fused x operation.

For a supported z length, the inverse transform loads sparse wavefunction
coefficients directly, supplies zeros for absent coefficients, and writes the
internal stick grid. The forward transform adds normalized coefficients to the
wavefunction output without a separate gather. Per-k inverse maps are rebuilt
with the basis and released with it. Unsupported z lengths and failure to
allocate the optional map retain packing, runtime cuFFT plans and gathering.
This changes FP64 operation order; it does not imply bitwise SCF trajectories.

Two-spin x transforms select the largest supported element count per thread
up to ten to limit register demand. Single-field y and z transforms retain their
independent element and block-size policy. These choices apply uniformly to
the generated registry and do not depend on a benchmark's grid dimensions.

`check_fused_fft_dispatch.cpp` tests available lengths against independent
long-double DFTs, normalized round trips, sparse maps, input preservation and
device switches. Its `--expect-unavailable` mode tests binaries with no loadable
image. `validate_portable_fft.py` links this check and the public potential,
density, band-group and allocation-failure checks to a configured application:

```bash
python3 tools/benchmarks/deltaspin/src/validate_portable_fft.py \
    --build build --output /path/to/unused-validation-directory
```

The dispatch and independent DFT checks also build without the application or
its MPI, BLAS, pseudopotential and orbital dependencies:

```bash
cmake -S tools/benchmarks/deltaspin/fft -B fft-check -DCMAKE_BUILD_TYPE=Release \
    '-DCMAKE_CUDA_ARCHITECTURES=70;80;90' -Dmathdx_ROOT=/path/to/mathdx
cmake --build fft-check --target check_fused_fft_dispatch --parallel 6
./fft-check/check_fused_fft_dispatch
```

The potential operator batches bands according to the actual compact-grid and
stick workspace sizes, targeting 256 MiB and at most 16 bands. A single band
may exceed that target on a large grid. Scratch arrays grow only when necessary;
FFT plans are cached by batch size and released with the basis. This policy
does not depend on atom count or a named benchmark system. Automatic template
instantiation increases compilation time and executable size.

The potential check covers every supported z length against full-grid cuFFT,
including sparse layouts, nonorthogonal cells, band tails, accumulation and
repeated basis lifetimes. `check_z_fft_fallback.cpp` additionally rejects the
optional inverse-map allocation at the allocator boundary, checks its contents
against the host basis, and verifies release after basis changes and destruction:

```bash
python3 src/link_check.py src/check_z_fft_fallback.cpp --cxxflag=-UNDEBUG \
    --ldflag=-Wl,--wrap=cudaMalloc --ldflag=-Wl,--wrap=cudaFree
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 results/check_z_fft_fallback
```

## Periodic projector phases

For identical atoms, each projector is a shared radial/angular factor times its
atomic position phase. On a binary fractional-coordinate grid, that phase
depends only on the integer reciprocal index modulo the grid periods. The CUDA
double-complex path groups plane waves by these residues, contracts the common
projectors with strided-batch GEMMs, and applies the smaller atomic phase matrix.
Expansion uses the adjoint factorization and adds to the existing Hamiltonian
output. The nonlocal pseudopotential and onsite constraint operators share this
factorization. Per-atom coefficients and spin coupling remain unchanged.

Both operators require CUDA double precision, one species and one rank per
pool. The nonlocal path additionally requires norm-conserving pseudopotentials.
Onsite folding requires at least 256 projector columns: narrower matrices do
not amortize packing on the measured V100 workloads and retain ordinary GEMMs.
Grid periods are chosen from 1, 2, 4, 8 and 16, with at most
512 residue groups, using exact alignment tests on the stored binary positions.
No coordinate is rounded. Atoms outside the selected grid retain all their
original projector columns in a compact ordinary GEMM. The estimated contraction
work, including padding and the atomic phase matrix, must be less than half the
complete product; otherwise the ordinary operator is used. The real-space FFT
grid, k points, symmetry setting and other calculation parameters are unchanged.

The nonlocal operator owns its plan and is recreated before each ionic SCF
cycle. Onsite plans belong to `FS_Nonlocal_tools`, which the onsite projector
recreates when initializing the geometry and basis. Reciprocal maps and phase
matrices are cached per k point. Common and irregular projector data are
refreshed whenever the underlying table is prepared; force and stress
derivatives invalidate that table before reuse. Expansion and overlap
contraction use the same prepared plan, with ordinary products as the fallback.

Scratch operates on at most 128 spin-component columns at a time. The onsite
plan has its own bounded workspace in addition to the nonlocal plan. The speed
benefit depends on the coordinate condition: these fixtures contain 15/31
aligned atoms and one displaced atom, while general thermally displaced
structures can retain the ordinary path.

`check_projector_fold.cpp` compares the public interface with independent CPU
products on a nonorthogonal cell, two k points and repeated k-point changes. It
covers zero, one and two displaced atoms, tails across the column block size,
nonzero output accumulation, unchanged input/padding, and both geometric and
work-based fallback. Complete converged physical comparisons are required
because regrouping changes floating-point accumulation order.

## Large complex matrix products

Double-complex products use cuBLAS ZGEMM3M when both output dimensions are at
least 128, the contraction dimension
is at least 64, and at least one dimension is 4096 or larger. Conjugate-transposed
projections also use this implementation for at least 128 output rows, 64 output
columns and a contraction dimension of at least 4096. This covers triangular
Davidson projection panels without increasing the panel width or computing more
of the unused lower triangle. Selection depends on shape rather than GPU model;
cuBLAS supplies the device implementation. Other shapes retain ZGEMM, including the historical
two-atom fixture's smaller projector products.

ZGEMM3M forms a complex product from three double-precision real products. Its
floating-point accumulation differs from ZGEMM, so mathematical equivalence
does not imply identical SCF trajectories or bitwise physical outputs. Complete
energy, force, stress and magnetic comparisons remain necessary. In particular,
small stress differences can exceed the measured two-reference-run spread and
must be reported explicitly instead of described as reference repeatability.

`check_large_gemm.cpp` compares the production interface with ZGEMM and sampled
long-double dot products. It covers all nine transpose/conjugation combinations,
padded leading dimensions, input preservation, nonzero complex accumulation,
cancellation-sensitive data, and shapes on both sides of the dispatch boundary.
`check_gemm_dispatch.cpp`, linked with `-Wl,--wrap=cublasZgemm3m`, verifies that
an eligible product reaches the 3M API and an adjacent ineligible shape retains
ordinary GEMM on the active GPU. Component acceleration is shape dependent;
complete matched calculations remain the performance acceptance criterion.
From the repository root, run the production-linked dispatch, matrix and
Davidson checks with:

```bash
python3 tools/benchmarks/deltaspin/src/validate_portable_gemm.py \
    --build build --output /path/to/unused-gemm-validation-directory
```

## Davidson workspace lifetime

CUDA double-complex solves retain one Davidson workspace across k points when
the wavefunction layout has a fixed basis dimension: two-component spinors or
basis-first storage. The solver owns its basis arrays and projected matrices
only for the enclosing Hamiltonian solve. Each diagonalization initializes its
active basis and refreshes the preconditioner from the current k point; no
eigenvectors, projected Hamiltonians or convergence decisions are reused.

The workspace is released before charge-density construction. A Davidson
failure also releases it before CG allocates its own work or changes the basis.
Other data types, device backends and variable-dimension layouts retain their
per-k-point solver lifetime. Projector preparation between k points can overlap
the retained allocation, so both complete-run timing and sampled peak memory
must be measured. The Davidson driver exercises repeated calls with changed
preconditioners, initial bases, projected matrices and basis-rotation callbacks.

## Fixed-overlap constraint responses

The wavefunction basis and its projected overlap stay fixed during one
constraint-field loop. The CUDA path factors each overlap as `S = U^H U`
once and caches `U^-H H U^-1` and `R U^-1`, where `R` contains the onsite
projector overlaps. Field corrections are assembled with these transformed
projectors, so each response uses a standard Hermitian eigensolve followed by
the projected eigenvector rotation. This retains the full positive-definite
overlap matrix and does not assume that the input basis is orthonormal.

The constraint object owns this cache and releases it before updating the
wavefunctions or starting Davidson. Original projected matrices remain
available for those operations. The implementation shares the existing
single-host-thread CUDA solver and BLAS contexts and the solver's retained
workspace. Optional cache allocation failure preserves the original generalized
response path. CPU and ROCm retain their existing implementations. The change
reorders double-precision operations, so full physical comparisons remain
necessary even though the transformed eigenproblem is mathematically equivalent.

`check_fixed_overlap.cpp` checks signed field changes against generalized
eigenvalues, normwise eigenpair backward error, overlap orthogonality, and
projected eigenvectors. It includes zero eigenvalues, nonorthogonal bases,
multiple k points in reverse order, repeated object lifetimes, immutable inputs
and caches, buffer guards, intervening generalized solves, and optional
allocation failure. It uses analytical matrices and requires no external data.

With cuSolver 11.7.3 or newer, subsequent responses can solve their independent
k-point matrices in one `cusolverDnXsyevBatched` call. The initial response
constructs H/S and uses serial solves. Required field coefficients, weighted
projectors and output scratch are allocated before optional batch storage.
Coefficients are uploaded once per response; the weighted-projector buffer is
reused across k points. A single k point, an older cuSolver, or failed optional
device allocation retains serial responses with the same field workspace.

Batch storage contains all k-point Hamiltonians, spectra, status values and the
library-reported workspace. It belongs to the fixed-overlap cache and is
released before Davidson. The solver uses FP64 arithmetic and checks every
matrix's convergence status. Rotated projectors are consumed in the original
k-point order before calculating occupations and weighting moments. Neither
chemical identity nor benchmark dimensions select the batching path.

`check_response_batch.cpp` checks analytical spectra, physical generalized
eigenpair residuals, overlap orthogonality, projector rotations, input and cache
preservation, and host-output guards for several matrix and batch dimensions.
It also exercises older-library and single-point fallback. Link the driver with
`-Wl,--wrap=cudaMalloc -Wl,--wrap=cudaFree` to enable its two optional-allocation
failure checks and verify release of tracked batch storage. Compare complete
application runs using identical CUDA libraries: a toolkit change is a separate
variable from response batching.

```bash
python3 src/link_check.py src/check_response_batch.cpp --cxxflag=-UNDEBUG \
    --ldflag=-Wl,--wrap=cudaMalloc --ldflag=-Wl,--wrap=cudaFree
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 results/check_response_batch
```

## Focused correctness checks

The C++ drivers in `src/` exercise the production public interfaces against
CPU, analytical or full-grid references. They cover FFT/potential/density
operations, charge conservation, Hartree and XC fields, stress, Davidson
eigenpairs, initialization/storage lifetimes, G+k assembly, magnetic-constraint
Hamiltonian updates and mixing behavior.

The fused CUDA spin-density path groups up to four consecutive occupied bands,
bounded by the existing 256 MiB spinor-transform workspace policy. Exact zero
weights end a group and are skipped before the transforms; negative and tiny
weights remain contributions. Each block retains the charge and three spin
densities across the group, then stores them once. Band order and the separate
rounding of each weighted contribution are preserved. Single-band groups use
the original density kernel, and unsupported FFT layouts and device paths
retain their ordinary transforms. Group scratch storage scales with the compact
FFT grid and must be included in whole-application peak-memory measurements.

`check_density_bands.cpp` compares the grouped public interface with ordered
single-band calls over cubic and rectangular grids, two k points, one to four
bands, padded plane-wave storage, and signed, zero and tiny weights. It checks
bitwise density agreement, integrated charge, immutable input, density guards,
the resource-limited single-band case, and unsupported FFT registry entries.

Onsite projector overlaps remain on the device while GPU operators consume
them. A one-rank pool already has the complete contraction and skips the host
staging used for MPI reduction. Multi-rank pools retain their host reduction.
Host moment and occupation readers acquire a synchronized snapshot through
`get_h_becp()` after the required overlap; they reacquire that view after each
overlap update. `check_projector_reduction.cpp` checks direct contractions on
CPU and GPU with one or two spin components, repeated k-point changes, varying
band counts, input/guard preservation, and repeated host-view acquisition.

The CUDA magnetic-constraint loop forms per-band atom moments from rotated
projectors on the device. Those moments remain available until all k-point
eigenvalues determine the occupations; a second kernel applies the occupations
in the original k-point and band order. Separate double-precision products and
sums preserve the host expression's rounding, and only the final three-component
atom moments return to the host. Temporary device moment storage scales as
`3 * nk * nbands * nat * sizeof(double)`, about 0.90/3.60 MiB for Fe16/Fe32,
replacing host storage for the full rotated projector matrices. The pool
reduction remains after moment accumulation. CPU and ROCm retain their host path.
`check_projector_moments.cpp` checks all components against the host complex
formula, scalar spin signs, zero and unequal projector counts, zero/tiny/negative
occupations, guards, repeated weighting, Pauli eigenstates and common-phase
invariance.

The Davidson driver also supplies independently constructed projected matrices
for a mixed, nonorthogonal initial basis and verifies their immutability. The
CUDA DeltaSpin path can reuse these matrices for norm-conserving inputs because
the local potential and initial wavefunctions stay fixed during the constraint
loop. It corrects the saved Hamiltonian with the difference between the current
constraint field and its saved value. Retries rebuild projections for their
rotated basis; PAW and ultrasoft inputs retain the ordinary projection path.

For the same eligible CUDA path, the constraint loop also retains the initial
Hamiltonian action in pinned host memory. The first Davidson attempt restores it
into its device buffer and applies the current constraint-field correction with
projectors prepared for that k point. The cache is refreshed for each constraint
loop; its allocation is reused across SCF steps and released with its owner.
Failed pinned allocation retains ordinary Hamiltonian evaluation. Host storage scales
with k points, bands, and plane waves: approximately 1.58 GiB for Fe16 and 6.30 GiB
for Fe32 in these fixtures. Transfers are included in the complete timings; this
cache has a measurable host-memory cost and does not eliminate all Hamiltonian
applications. The Davidson driver checks the supplied action against direct
evaluation and verifies that later calls can use ordinary evaluation again.

Hamiltonian actions are uploaded directly from pinned host storage into the
Davidson destination on the default CUDA stream. The transfer keeps no duplicate
device wavefunction block, independently of the number of k points. Requests
may repeat or reorder k points and use different output leading dimensions.
A completion event covers all host reads and is drained before the transfer
object is destroyed, so the caller can safely release or reuse its host cache.
This trades inter-k-point transfer overlap for lower peak device memory;
complete-run timing and memory are measured together.

The constraint-field and weighted-projector scratch arrays remain allocated
throughout the solve. Cached nonlocal projections and current onsite projectors
are prepared before the action restore. The complex-double
cuSolver generalized eigensolver retains its largest requested workspace until
its handle is destroyed. These lifetimes avoid repeated per-k-point and
per-eigenproblem allocation and release; they do not
change the eigensolver or the convergence criteria.

`check_hpsi_transfer.cpp` verifies bitwise restores, input and guard preservation,
varying leading dimensions, repeated/reordered requests with 1/2/3/8 k points,
and host-storage reuse immediately after destruction with pending transfers.
`check_solver_workspace.cpp` checks
analytic generalized eigenpairs across growing and shrinking problem sizes,
two handle lifetimes, and recovery after a non-positive-definite overlap matrix.

Host caching overlaps the Hamiltonian projection on a separate CUDA copy stream.
For the eligible norm-conserving path, the overlap operator is exactly the
identity, so the projected overlap is formed directly from the wavefunctions
while the action is copied. This retains the actual Gram matrix and does not
assume that the input columns are orthonormal. Generalized overlap operators
wait for the copy before reusing the action buffer. The host cache is complete
before the projection routine returns. `check_subspace_cache.cpp` exercises
asynchronous producers, CPU/GPU paths, generalized and Euclidean overlaps,
padded inputs, pinned/pageable host storage, and repeated cache lifetimes.

For CUDA norm-conserving calculations with one rank per pool, the nonlocal
projector overlaps also follow the Davidson basis. A Ritz rotation
`Psi_out = Psi_basis C` permits the smaller update
`P^H Psi_out = (P^H Psi_basis) C`. The projector table depends on the cell and
k point and stays fixed during this SCF cycle. The next constraint-loop
Hamiltonian action can therefore consume these overlaps while applying the
current nonlocal coefficients. No Hamiltonian or wavefunction values are held
fixed across SCF steps.

The per-k overlaps live in host memory; temporary device matrices hold the
current projected basis and its rotation. The host allocation is approximately
41 MiB for Fe16 and 163 MiB for Fe32. Device scratch scales as the projector
count times the band count and Davidson basis capacity, and is included in
peak-memory measurements. Ionic/solver initialization, unobserved wavefunction
updates and density-only rotations invalidate the cache. A Davidson failure
invalidates the affected k point before the ordinary CG fallback. Normal
Hamiltonian actions recapture the initial basis on a Davidson retry.

The operator action is supplied explicitly for one Hamiltonian call, and basis
observations are supplied explicitly for one solve. Existing virtual `hPsi`
overrides retain their entry point. `check_operator_action.cpp` verifies
ordinary/replaced/ordinary calls for one- and two-component wavefunctions;
`check_davidson.cpp` checks transported linear projections against the actual
rotated wavefunctions through restarts and final rotations, including supplied
initial Hamiltonian actions.

`link_check.py` reuses the active Ninja application's compile and link commands.
Set `ABACUS_HOME` to the installation containing `build/` and load its environment:

```bash
mkdir -p results
python3 src/link_check.py src/check_hartree.cpp
OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=1 results/check_hartree
```

The projector-reduction driver additionally reads the Fe orbital from the
asset directory. Exercise both local and distributed contractions:

```bash
python3 src/link_check.py src/check_projector_reduction.cpp --cxxflag=-std=c++17
for ranks in 1 2; do
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun --bind-to none -np "$ranks" \
        results/check_projector_reduction "$PWD/assets"
done
```

The existing-mixing-test wrapper additionally needs GoogleTest/GoogleMock headers
and static libraries. Pass their installed paths with repeated `--cxxflag` and
`--library` options; use `--cxxflag=-UNDEBUG` when exercising assertions. These
focused checks supplement complete application comparisons on both supercells.
