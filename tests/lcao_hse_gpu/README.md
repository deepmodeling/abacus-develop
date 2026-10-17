# Experimental LCAO HSE GPU integration

Fork-only experiment; not production-ready. ABACUS base: 8ae37fcd6.
Runtime selection uses the existing `device cpu/gpu`. Build with
`ENABLE_LIBRI=ON`, `USE_CUDA=ON`, `ENABLE_LIBRI_GPU=ON`, and explicit
`LIBRI_DIR`, and a CUDA DDLA installation selected with `DDLA_ROOT` or
`CMAKE_PREFIX_PATH`. Explicit `DDLA_INCLUDE_DIR` and `DDLA_LIBRARY` overrides
remain supported.

ABACUS's Find modules provide `LibRI::LibRI` and `DDLA::DDLA`; no CMake changes
in either dependency are needed. The LibRI target carries its headers, C++
requirements, LibComm/cereal and available MPI/OpenMP/math targets. GPU builds
add DDLA and the backend macros through those targets. ABACUS's existing
compile-requirements collector keeps these header-only definitions consistent
across its translation units, with final linking in `source/CMakeLists.txt`.
CPU builds do not discover or link DDLA. Enabling GPU RI without LibRI/CUDA
is an error. A configure-time compile/link probe instantiates real and complex
GPU contractions to reject incompatible dependency APIs; it does not execute
GPU code or verify numerical correctness. The dependency patches below are
still required.

Dependencies used in the SAI experiment:

- LibRI 6ce3c8480ae38a085746dc81cc4dbb59b8074430; archive SHA256
  09c8f766687e5c2236bfe580f372ea5ecef38e935a56c6a24a786d0430ff8fc4.
- LibDDLA 9c404ac65148ec1b48977b846f38ba6ef9e134f7; archive SHA256
  bc349e64c0a308aa9cd0f766e9a766e942bbc974f557522604bd4af9e1ce9baa.
- LibComm and cereal use the existing ABACUS toolchain pins.

Apply `libri-compat.patch` to LibRI after normalizing its two affected files
to LF, then use `git apply --unidiff-zero /path/to/libri-compat.patch`
from that dependency root. The patch omits context to preserve upstream tab
indentation without introducing whitespace diagnostics in this repository. It repairs the DDLA diagnostic helper name, the split cleanup API,
and GPU contraction accumulation that otherwise ignores `fac_add_Ds`.
The original GPU force regression fails while energy passes; patched force
and stress must be verified against the unchanged references.
The build template records equivalent edits used in the actual experiment.

Use the existing `.ci/slurm/runner.py` transport, build, worker and artifact
machinery. Experiment controls and full raw results are retained locally at
`/mnt/g/work/agent/abacus/lcao-hse-results/attempt*`.
The copied control runner supports `autotest_gpu` and `hse_benchmark`, as in
the earlier fork experiment; these extensions are not in the base runner.
Set OMP_NUM_THREADS=1 after loading cluster modules.

Limitations: LibRI's 5 GiB temporary limit is not a total GPU memory budget.
Multi-node GPU mapping and all GPU contraction labels are not validated.
Do not infer production MD readiness from SCF smoke tests.

## Benchmark interpretation

The initial Si8/Si32 benchmark compares CPU ScaLAPACK with GPU cuSOLVER;
wall time therefore includes a solver change. Use the additional GPU-device
build with `ENABLE_LIBRI_GPU=OFF` to isolate RI offload.
`exx_hybrid_step=50` is an outer-loop cap, not a convergence guarantee.
Inspect the printed outer deltaE before calling a result converged.
GPU process memory is sampled at 100 ms; it is not a strict allocator peak
or a proof of an upper bound for other systems.

## Running a backend comparison

On a GPU allocation with the cluster modules loaded, set `OMP_NUM_THREADS=1`
and make DDLA's `lib` directory available in `LD_LIBRARY_PATH`. Build two
executables with the same patched LibRI and identical ABACUS options, changing
only `ENABLE_LIBRI_GPU`. Both runs use `device gpu` and `ks_solver cusolver`.

```bash
python3 tests/lcao_hse_gpu/benchmark.py \
  --abacus /path/to/gpu-ri/abacus --baseline /path/to/cpu-ri/abacus \
  --source /path/to/abacus-source --output /path/to/new-results \
  --orbital Si_dojo_6au.orb --repeats 2
```

The source must contain the test PP_ORB assets. `ABACUS_CPU_RI` may supply the
baseline path when invoked by the local `.ci/slurm` experiment controls.
The benchmark alternates backend order, checks outer energy convergence and
iteration count, and compares total energy, every force component and every
stress component. Thresholds are 1e-6 eV/atom, 1e-4 eV/Angstrom and 1e-3 kbar.
These are backend parity tests, not basis/cutoff convergence studies; the
small SZ system can break symmetry under prolonged outer iteration.

GPU mode deliberately requires an externally supplied, patched LibRI tree.
The automatic download still uses ABACUS's original CPU dependency. Downloading
an unpatched newer revision is insufficient to reproduce the working build.
The build option defaults to OFF; the HSE approximation is unchanged.
See RESULTS.md for measured outcomes and remaining validation limits.

## Checking the CMake dependency interface

The small consumer below links only `LibRI::LibRI`, without ABACUS's common
include paths or link flags. It checks headers, backend macros, the CPU default
in GPU-enabled LibRI, and transitive BLAS linkage. Run with GPU support OFF
and ON; the ON variant requires CUDA and the patched dependencies above.

```bash
cmake -S tests/lcao_hse_gpu/cmake -B /path/to/consumer-build \
  -DLIBRI_DIR=/path/to/libri -DLIBCOMM_DIR=/path/to/libcomm \
  -Dcereal_DIR=/path/to/cereal/lib/cmake/cereal \
  -DENABLE_LIBRI_GPU=ON -DDDLA_ROOT=/path/to/ddla
cmake --build /path/to/consumer-build
ctest --test-dir /path/to/consumer-build --output-on-failure
```
