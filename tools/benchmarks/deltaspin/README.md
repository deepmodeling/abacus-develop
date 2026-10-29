# Fe DeltaSpin performance benchmarks

The timing benchmark matrix uses 16-atom and 32-atom bcc Fe supercells, one GPU,
one MPI rank, and either one or twelve CPU threads. Each thread configuration
compares the original and optimized binaries using identical input and asset
hashes. Thread-count effects are reported separately from implementation speedup.

The reference source is tag `delta-spin-yinan`, commit
`19ade18599bfdf9ee78a289bf798a9b2583d52cd`. Keep its executable separate from
the development executable. Historical two-atom measurements are not evidence
for performance on these supercells.

## Inputs

The structures use a 2.866 A conventional bcc lattice, repeated 2x2x2
for Fe16, 2x2x4 for Fe32, or 3x3x3 for the Fe54 capacity check. One atom is displaced by 0.02 A along x to provide nonzero forces for
comparison. The local moment target is 2.2 uB along +z on every atom.

All cases use 100 Ry, double precision, noncollinear spin, spin-orbit coupling,
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
multiply work as `8*m*n*k` for ZGEMM or `6*m*n*k` for ZGEMM3M. These shape-derived
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
whose prime factors are 2, 3 and 5. Each length is instantiated in a separate
translation unit. Lengths that require a cuFFTDx workspace, and lengths outside
this registry, retain the cuFFT path. The optional build requires CUDA and
cuFFTDx 1.6 or newer, discoverable through the MathDx CMake package.

The application selects the fused path for double-complex spinors on Volta,
with one MPI rank and an eligible plane-wave layout. The x and y dimensions
are selected independently; z transforms use runtime cuFFT plans. Rectangular
and nonorthogonal cells use the same implementation. A missing y template can
use cuFFT for the y-z transforms while retaining the fused x operation.

Two-spin x transforms select the largest supported element count per thread
up to ten to limit register demand. Single-field y transforms retain their
independent element and block-size policy. These choices apply uniformly to
the generated registry and do not depend on a benchmark's grid dimensions.

The potential operator batches bands according to the actual compact-grid and
stick workspace sizes, targeting 256 MiB and at most 16 bands. A single band
may exceed that target on a large grid. Scratch arrays grow only when necessary;
FFT plans are cached by batch size and released with the basis. This policy
does not depend on atom count or a named benchmark system. Automatic template
instantiation increases compilation time and executable size.

## Large complex matrix products

On Volta devices (compute capability 7.0), double-complex products use cuBLAS
ZGEMM3M when both output dimensions are at least 128, the contraction dimension
is at least 64, and at least one dimension is 4096 or larger. Device selection
is established with the BLAS handle. Other shapes and architectures retain
ZGEMM. The size gate targets wide plane-wave products; narrow projection panels
and the historical two-atom fixture retain their conventional accumulation.

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

## Focused correctness checks

The C++ drivers in `src/` exercise the production public interfaces against
CPU, analytical or full-grid references. They cover FFT/potential/density
operations, charge conservation, Hartree and XC fields, stress, Davidson
eigenpairs, initialization/storage lifetimes, G+k assembly, magnetic-constraint
Hamiltonian updates and mixing behavior.

The Davidson driver also supplies independently constructed projected matrices
for a mixed, nonorthogonal initial basis and verifies their immutability. The
CUDA DeltaSpin path can reuse these matrices for norm-conserving inputs because
the local potential and initial wavefunctions stay fixed during the constraint
loop. It corrects the saved Hamiltonian with the difference between the current
constraint field and its saved value. Retries rebuild projections for their
rotated basis; PAW and ultrasoft inputs retain the ordinary projection path.

For the same eligible CUDA path, the constraint loop also retains the initial
Hamiltonian action in pinned host memory. The first Davidson attempt copies it
directly into its existing device buffer and applies the current constraint-field
correction with projectors prepared for that k point. The cache is refreshed for
each constraint loop; its allocation is reused across SCF steps and released
with its owner. Failed
pinned allocation retains ordinary Hamiltonian evaluation. Host storage scales
with k points, bands, and plane waves: approximately 1.58 GiB for Fe16 and 6.30 GiB
for Fe32 in these fixtures. Transfers are included in the complete timings; this
cache has a measurable host-memory cost and does not eliminate all Hamiltonian
applications. The Davidson driver checks the supplied action against direct
evaluation and verifies that later calls can use ordinary evaluation again.

Host caching overlaps the Hamiltonian projection on a separate CUDA copy stream.
For the eligible norm-conserving path, the overlap operator is exactly the
identity, so the projected overlap is formed directly from the wavefunctions
while the action is copied. This retains the actual Gram matrix and does not
assume that the input columns are orthonormal. Generalized overlap operators
wait for the copy before reusing the action buffer. The host cache is complete
before the projection routine returns. `check_subspace_cache.cpp` exercises
asynchronous producers, CPU/GPU paths, generalized and Euclidean overlaps,
padded inputs, pinned/pageable host storage, and repeated cache lifetimes.

`link_check.py` reuses the active Ninja application's compile and link commands.
Set `ABACUS_HOME` to the installation containing `build/` and load its environment:

```bash
mkdir -p results
python3 src/link_check.py src/check_hartree.cpp
OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=1 results/check_hartree
```

The existing-mixing-test wrapper additionally needs GoogleTest/GoogleMock headers
and static libraries. Pass their installed paths with repeated `--cxxflag` and
`--library` options; use `--cxxflag=-UNDEBUG` when exercising assertions. These
focused checks supplement complete application comparisons on both supercells.
