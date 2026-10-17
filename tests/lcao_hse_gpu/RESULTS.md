# SAI LCAO HSE GPU experiment — 2026-10-17

Experimental integration only; no upstream PR. Computing resources were
provided by the Open Source Supercomputing Center of SAI (https://www.open-sai.com/).

ABACUS base: `8ae37fcd6cb9043a4a0c7ebde2c062d97f393bdb`.
Runtime source: `b408a31374fa6134618b60d075cb50d8cf1b1c0c`.
Subsequent changes are test/documentation and a configure-time guard against
fetching an unpatched GPU dependency; runtime C++ is unchanged.
LibRI/DDLA revisions and patches are in README.md and libri-compat.patch.
Hardware: one V100 16 GB, one MPI rank, OMP_NUM_THREADS=1, CUDA 12.9.86,
GNU 13.3.0, OpenMPI. Both comparison executables use GPU cuSOLVER.
The CPU RI baseline uses the same newer LibRI, not the original b0eff7a.

## Build and regression evidence

| Job | Result |
| --- | --- |
| 1818882 | Patched GPU RI build passed |
| 1818916_0 / 1818916_2 | Existing HSE Gamma / spin-polarized Gamma energy, force, stress references passed |
| 1819337_0 | Existing complex multi-k HSE references passed with GPU cuSOLVER |
| 1819559 | CPU RI control build passed, CUDA and same newer LibRI retained |
| 1819665 | SZ Si8/Si32 matched-backend comparison passed |
| 1819761 | DZP Si8/Si32, two alternating-order repeats per backend: all 8 calculations and 4 paired comparisons passed |

Failures retained in raw artifacts: unpatched GPU force accumulation failed
reference tests despite passing energy; incompatible DDLA helper and LibRI
cleanup APIs required fixes. The installed ELPA/GenELPA paths were unsuitable
for these GPU comparisons (including an illegal-instruction exit); cuSOLVER
is used consistently in the final matched comparison. Earlier failed jobs
are not counted as all-pass suites.

With default `scf_ene_thr=-1`, the Si32 CPU RI and GPU RI trajectories had
matching early EXX energies but different outer stopping times. GPU RI reached
the 50-step cap and the final states differed. Explicit `scf_ene_thr=1e-8` eV
produced equal outer counts and matching results. This supports convergence
sensitivity as the cause of that discrepancy, but does not establish default
parameter robustness or production MD readiness.

## SZ matched comparison

| System | CPU RI wall s | GPU RI wall s | CPU/GPU sampled MiB | Outer steps |
| --- | ---: | ---: | ---: | ---: |
| Si8 SZ | 8.91 | 9.19 | 422 / 498 | 4 / 4 |
| Si32 SZ | 24.31 | 24.13 | 438 / 530 | 5 / 5 |

Single paired measurements: no demonstrated speedup. Total energy differences
were below 4e-11 eV, forces identical at printed precision, maximum stress
component difference 9.98e-7 kbar. Si32 electronic EXX time increased from
0.50 to 0.80 s, force EXX from 0.80 to 1.17 s, stress EXX from 0.96 to 1.57 s.
Small contractions can be slower after offload.

## DZP matched comparison

Si_dojo_6au.orb is a 2s2p1d basis. Two repeats per backend, alternating
CPU/GPU order; means below include SCF, force and stress calculation.

| System | CPU RI wall s | GPU RI wall s | Wall speedup | CPU/GPU sampled MiB |
| --- | ---: | ---: | ---: | ---: |
| Si8 DZP | 30.99 | 29.64 | 1.05x | 442 / 558–564 |
| Si32 DZP | 250.79 | 73.81 | 3.40x | 450 / 2756 |

Si32 mean electronic EXX contraction time: 50.20 -> 5.36 s (9.36x).
EXX force: 64.95 -> 6.89 s; EXX stress: 85.09 -> 9.68 s. Both backends
used 7 outer iterations for both sizes. Maximum paired energy difference:
1.03e-12 eV/atom; force difference zero at printed precision; stress difference
1.11e-7 kbar. Symmetric silicon forces are zero, so nonzero-force coverage
comes from the existing HSE regression references, not these benchmark cells.

These speedups are relative to single-thread CPU RI; a tuned multi-core CPU
comparison has not been performed. They are not a universal HSE speedup.
All eight memory logs contain a single GPU process PID. Si32 GPU RI added
2306 MiB (2.25 GiB) over the GPU-solver/CPU-RI baseline: a substantial increase,
not a guarantee of small memory overhead. If the requirement is less than
2x *total* GPU memory relative to that baseline, this implementation fails it.

## Scope and limitations

Memory is sampled every 100 ms and includes the GPU eigensolver. It is not
a strict allocation peak. LibRI's hardcoded 5 GiB temporary batching limit
excludes other simultaneous allocations, and is not a total memory bound.
No guarantee that large-system memory stays below twice a baseline is made.
Multi-GPU/multi-node, long MD, sustained allocation behavior, SOC and larger
production cells remain unvalidated. Do not enable this path by default.

Local raw results: `/mnt/g/work/agent/abacus/lcao-hse-results/`.
The corresponding remote Slurm archives are under
`~/abacus_hse_gpu_experiment/archives/lcao-hse/`.
The experiment reused `.ci/slurm` transport/worker/artifact handling with local
controls; benchmark and dependency fixes are preserved here without importing
unrelated runner changes.

## Follow-up priorities suggested by the evidence

1. Keep the small SZ contraction workload on CPU unless profiling shows a
   crossover. GPU dispatch/packing/transfer overhead can dominate these blocks.
2. Profile and optimize `refresh_coulomb`/initialization separately. Offloading
   contraction GEMMs cannot accelerate this measured CPU-dominated phase.
3. LibRI currently constructs a GPU context per `cal_loop3` call and creates,
   packs and uploads input/temporary/output objects per contraction label.
   Reuse a single explicitly capped allocation arena and queue where lifetimes
   allow it; do not keep two full workspaces or unbounded per-shape caches.
4. Investigate batching by shape and selective reuse of invariant tensors,
   measuring transfer/packing/kernel time separately. Persistent tensors must
   share a total memory budget with temporary and output buffers.

These are proposed next experiments, not demonstrated speedups. The current
change connects an existing double-precision GPU contraction implementation;
it does not introduce mixed precision, new screening thresholds, or a different
HSE approximation. Floating-point accumulation order can differ (DDLA uses
atomic accumulation in relevant batched kernels), so bitwise identity is not
expected. Stronger screening or mixed precision would need separate accuracy
validation and are outside this experiment.

## CP2K comparison

The relevant transferable ideas are block batching and explicit memory-aware
contraction scheduling, not replacing an ABACUS matrix call with DBCSR directly.
CP2K's [RI HFX implementation](https://github.com/cp2k/cp2k/blob/0adc3f08c83fcdafa7d545917ddb4e901470c070/src/hfx_ri.F)
uses DBCSR operations and explicitly batches occupied-MO/AO indices to save
memory. DBCSR's [libsmm_acc](https://github.com/cp2k/dbcsr/blob/8f63c8f49ce807b5127cf8023b4684f579d19225/src/acc/libsmm_acc/README.md)
executes stacks of small matrix products and tunes kernels for matrix dimensions
and GPU models. LibRI already has variable-size batched contractions, so merely
adding batching is not a new optimization; block-shape distribution, packing,
transfer and kernel efficiency should determine the next change. Different
RI formulations/data structures preclude assuming CP2K's speedup transfers
unchanged to ABACUS.

## Local checks

- Repository `.ci/slurm/test_runner.py`: 45 tests passed. Running that suite
  from the copied experiment controls initially failed because it assumes the
  original 117-case configuration and repository layout; the repository suite
  was rerun in its intended location without modifying expectations.
- `agent_governance_check.py --staged`: no findings.
- Python benchmark syntax/CLI and shell template syntax checks passed.
- LibRI patch apply check passed after LF normalization (`--unidiff-zero`).
- Missing GPU LibRI dependency produces the expected CMake diagnostic; this
  guard was checked remotely with CMake 3.31.6.
- Final staged whitespace check passed.

## Single-node multi-GPU memory scaling

Follow-up Slurm job 1832072 passed on node 16v100n48. The same saved GPU RI
binary and Si32 DZP input ran sequentially with 1, 2 and 4 MPI ranks, one rank
per GPU, OMP_NUM_THREADS=1, GPU cuSOLVER. All runs used 7 outer iterations.
One measurement per configuration, on the same four-V100 allocation:

| Active GPUs / ranks | Per-process sampled peak MiB | Wall seconds |
| --- | --- | ---: |
| 1 | 2756 | 63.78 |
| 2 | 1682, 1682 | 57.23 |
| 4 | 1158, 1106, 1106, 1158 | 52.71 |

A live nvidia-smi UUID/PID snapshot confirmed the four listed processes were
on four distinct GPUs; a separate live query confirmed distinct GPUs for the
two-rank run. Peaks are sampled at 100 ms, not allocator-level bounds.
Maximum total-energy difference from one rank: 5.10e-11 eV. Forces agree at
printed precision; maximum stress-component difference: 1.30e-7 kbar.

Thus distributing RI work *does* lower per-card memory in this tested case:
39% reduction with two cards and 58% with four. Total replicated memory grows,
and the slowest card does not follow ideal 1/N scaling. Performance improved
only 1.21x from one to four cards on this small system. These measurements
support using more GPUs to relieve capacity pressure, but do not establish
large-system or multi-node scaling or guarantee arbitrary capacity growth.
The CPU-RI baseline's small GPU footprint is not an appropriate standalone
rejection criterion for GPU offload; the relevant capacity metric is maximum
per-card memory as system size and GPU count increase.

Controls, full input/output, per-PID memory traces and UUID mapping snapshot:
`/mnt/g/work/agent/abacus/lcao-hse-results/multigpu/`.
Remote archive: `archives/lcao-hse/lcao-multigpu-20261017-1.tar.gz` beneath the
experiment project root. No runtime code was changed for this test.

Further size-dependent memory and multi-GPU checks (Si16/Si64, job 1832165)
are documented in [SCALING.md](SCALING.md). They reveal a memory plateau
between one and two cards for Si64; do not extrapolate the Si32 reduction
as a universal 1/P law.

## ABACUS-side CMake integration — 2026-10-18

The dependency adapter now lives in ABACUS's Find modules. `LibRI::LibRI`
and `DDLA::DDLA` carry headers, language requirements, backend definitions
and link dependencies. These are ABACUS-provided targets, not upstream
LibRI package exports. `ENABLE_LIBRI_GPU` remains OFF by default; runtime
selection, numerical kernels, memory allocation and dependency patches are
unchanged by this CMake work.

Full CPU/GPU builds used source `33d400fd6012262c781697716f1a0e19950843d0`,
GNU 13.3.0, CMake 3.31.6 and CUDA 12.9.86 on SAI. The MKL follow-up used
`8dd69254063302fb70b08b01b08cc284a2d9e37e` and oneMKL 2024.2. It changes
the MKL adapters to named imported interfaces so CMake can export their
transitive requirements into a compile/link probe. Existing target names
and numerical link ordering are retained. No runtime C++ changed.

| Verification | Result |
| --- | --- |
| Job 1833459: standalone CPU/GPU consumers linking only `LibRI::LibRI` | Both built; both CTest runs passed |
| Job 1833459: real/complex GPU contraction compile/link probe | Passed with patched dependencies; rejected unpatched LibRI as expected |
| Job 1833459: full CPU build, CUDA OFF, DDLA discovery disabled | Built; three HSE cases passed all 10 existing reference checks |
| Job 1833459: full GPU build, GPU RI ON | Built; all dynamic dependencies resolved |
| Jobs 1833493_0/1/2: GPU HSE Gamma, complex multi-k and spin-polarized Gamma | All three passed, with 10 existing reference checks in total |
| Job 1833601: MKL static/dynamic, CPU/GPU standalone consumers | All four configure/build/CTest combinations passed |
| Job 1833601: original CPU LibRI `b0eff7a`, DDLA discovery disabled | Configure, build and CTest passed |
| Job 1833601: previous MKL adapter with the new GPU probe | Reproduced the missing-target error; corrected adapter passed |

HSE cases use the repository's unchanged references. Gamma and spin-polarized
Gamma check energy, force and stress; the complex multi-k reference checks
energy. Copied CPU inputs use ScaLAPACK where GenELPA was explicitly selected;
copied GPU inputs use cuSOLVER. All runtime tests use one MPI rank and
`OMP_NUM_THREADS=1`. This is build/regression validation, not a timing
comparison or a new test of multi-GPU scaling.

The CPU executable's `ldd` has no DDLA/CUDA libraries. The RI and base object
compile flags contain the three GPU RI definitions only in the GPU build.
MKL consumer link lines select the requested static/shared MKL libraries;
CPU consumers do not link DDLA. The MKL checks validate dependency interfaces,
not a full MKL ABACUS build. The configure-time probe cannot detect numerical
errors such as the force-accumulation issue; the existing LibRI patches and
numerical regressions remain necessary.

Three local invalid configurations (GPU RI with LibRI OFF, CUDA OFF or MPI
OFF) fail with the intended diagnostic. The repository command
`python3 -m unittest discover -s .ci/slurm -p 'test_runner.py' -q` passed all
45 tests. Governance and whitespace checks passed.

Two experiment-control issues are retained in the raw evidence: job 1833432
stopped before compilation because its DDLA reuse path was wrong; the corrected
run passed. Result collection then rejected two resource-directory symlinks
introduced by the local CPU test control. Only those symlink entries were
omitted from a copy of the archive, which passed the original safe extractor;
calculation logs and numerical outputs were retained. No reference or runner
safety check was relaxed.

Controls and logs: `/mnt/g/work/agent/abacus/lcao-hse-results/cmake-integration2/`
and `/mnt/g/work/agent/abacus/lcao-hse-results/cmake-mkl/`.
Full-build remote archive:
`archives/lcao-hse/cmake-targets-20261018b-1.tar.gz` under the experiment root.
MKL controls and results are under `checks/cmake-mkl-20261018/` there.
