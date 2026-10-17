# LCAO HSE GPU memory and strong-scaling verification

## Measurements

Slurm job 1832165 passed on SAI node 16v100n41 (four V100 16 GB GPUs).
Si16 and Si64 each ran with 1, 2 and 4 MPI ranks, one rank per GPU and
OMP_NUM_THREADS=1. Saved GPU RI binary and dependency patches are unchanged.
All configurations use GPU cuSOLVER, DZP Si_dojo_6au.orb, Gamma sampling,
SCF energy threshold 1e-8 eV and density threshold 1e-9. Forces/stress included.
Si32 is the earlier same-node comparison, job 1832072 on 16v100n48.
One run per point: timings are observations, not statistically established
speedups. Compare timings within each row; nodes differ between Si32 and the
new measurements. CPU resources also grow with MPI/GPU count.

| Atoms | Max per-card sampled MiB, 1 / 2 / 4 cards | Wall seconds, 1 / 2 / 4 cards | Four-card speedup | Four-card parallel efficiency |
| --- | --- | --- | ---: | ---: |
| 16 | 860 / 708 / 616 | 40.37 / 39.16 / 38.02 | 1.06x | 26.5% |
| 32 | 2756 / 1682 / 1158 | 63.78 / 57.23 / 52.71 | 1.21x | 30.3% |
| 64 | 6568 / 6262 / 3980 | 178.30 / 137.05 / 101.07 | 1.76x | 44.1% |

Memory is nvidia-smi per-process usage sampled every 100 ms. Mapping snapshots
confirm the Si64 two-rank processes occupy two distinct GPUs, and four-rank
processes four distinct GPUs. Si64 per-process peaks on four cards were
3980, 3426, 3426 and 3962 MiB. These individual peaks need not occur at the same
time; their sum is not a measured simultaneous aggregate peak. This includes
non-RI allocations; an allocator trace separating RI components was not taken.

Si16 used 8 outer steps in all configurations, with 74 / 74 / 73 total inner
iterations. Si64 used 7 outer steps, with 63 / 62 / 62 inner iterations.
Maximum energy difference within a size: 2.05e-10 eV total; force differences
zero at printed precision; stress-component differences <= 1.88e-7 kbar.
These symmetric cells have zero total forces; existing regression cases provide
nonzero-force validation. Inner counts differ slightly, so these wall times
are time-to-solution scaling, not perfectly fixed-work kernel benchmarks.

| Si64 EXX timer | 1 card s | 2 cards s | 4 cards s | 1-to-4 speedup |
| --- | ---: | ---: | ---: | ---: |
| cal_exx_elec | 25.42 | 22.65 | 12.30 | 2.07x |
| cal_exx_force | 33.62 | 24.95 | 17.46 | 1.93x |
| cal_exx_stress | 43.92 | 30.74 | 21.72 | 2.02x |
| init | 11.92 | 12.73 | 13.71 | 0.87x |

Timers include their enclosing routine's work, not just CUDA kernel duration.
They do not isolate MPI communication, packing or synchronization. Si16's
cal_exx_elec instead rises from 1.74 to 2.53 s with four cards. Small workloads
are unsuitable for assuming efficient strong scaling.

## What the data establish

- Multi-GPU distribution relieves per-card capacity pressure, but not as 1/P
  universally. Si64 improves only 4.7% with two cards and 39.4% with four.
- A plateau is compatible with the implementation's per-rank temporary cap:
  min(T/P, 5 GiB) stays near 5 GiB while each rank still exceeds the cap. This
  is a code-supported explanation, not an allocation-level attribution of
  the Si64 peak. Other allocations and load imbalance also matter.
- At fixed card count, memory grows substantially with size: 860 -> 2756 ->
  6568 MiB for one card and 616 -> 1158 -> 3980 MiB for four cards. These finite
  points do not establish a universal power law. Supercell shapes change
  (2x1x1, 2x2x1, 2x2x2), affecting periodic neighborhoods and screening.
- Multi-GPU speed improves for the larger workload, but four-card efficiency
  remains about 44% end-to-end and 52% for electronic EXX. Communication scaling
  or a single bottleneck cannot be inferred without profiling.
- Increasing cards proportionally with atoms does not guarantee fixed per-card
  memory. Unified memory-budgeted batching remains useful even with MPI.
- No multi-node or weak-scaling benchmark has been completed. No host-RAM peak
  has been measured. The results do not establish arbitrary-size scalability.

## Source-level allocation model

For one rank, one contraction label, the allocator-level components are:

M_RI = s (I + U + O + min(T, 5 GiB / s)) + B + C,

where s is the scalar size (8 bytes real double, 16 bytes complex double),
I is the total element count of unique uploaded inputs, U of unique accumulated
intermediates, O of unique outputs, T of unsegmented temporary products, B of
per-task pointer/dimension arrays, and C of backend context/workspace allocations.
Take the maximum over simultaneously live objects/phases; add other ABACUS GPU
allocations to obtain whole-process memory. This is an accounting model, not a
prediction from atom count alone. Certain labels can have different live objects.

Input, Mul and Output deduplicate by atom-pair key before allocating full local
data. Tmp allocates min(totalSize, limit/sizeof(T)) scalar elements, but its
pointer array is still allocated for the whole task list. In the examined
three-index contraction, a tmp entry and dimension/pointer entries are added
per surviving triple; two-index inputs and accumulated products are deduplicated.
The hardcoded budget is per temporary data buffer, not for total allocations.

For fixed species/basis and fixed periodic-image range, the dense unscreened
pair storage can grow quadratically with atom count; the triple task list and
unsegmented temp volume can grow cubically. This is an upper-envelope counting
argument for these contraction stages, not a claim that screened silicon follows
N^3, nor a universal complexity of ABACUS HSE. Fixed-density finite-range
screening can make the surviving neighborhoods bounded and substantially lower
the growth. Auxiliary basis size, image count, cutoff and screening matter.

The default Parallel_LRI_Equally partitions two index ranges with MPI. The GPU
loop builds its batches from rank-local lists and requests corresponding tensor
maps. Tasks can divide approximately as 1/P when balanced. Shared inputs/outputs
can be replicated across one process-grid direction; whole-process memory need
not follow 1/P. In an ideal balanced dense two-dimensional grid, a matrix slice
shared along a row/column can follow N^2/sqrt(P) rather than N^2/P. This is an
illustration of replication, not an empirically fitted scaling law here.

CPU-side original maps, packed host input copies, output downloads and task
metadata coexist. Host memory cannot be inferred from nvidia-smi measurements;
no host-memory asymptotic scaling has been measured in this experiment.

## Source locations and reproducibility

The analysis refers to LibRI revision
`6ce3c8480ae38a085746dc81cc4dbb59b8074430`, under `include/RI/`:

- `ri/gpu/GPU_Data_Input.h`: deduplicated input packing and upload.
- `ri/gpu/GPU_Data_Tmp.h`: capped temporary data, full-batch pointer allocation.
- `ri/gpu/GPU_Data_Mul.h`, `GPU_Data_Output.h`: full local accumulated/output data.
- `ri/LRI-cal_loop3_GPU.hpp`: rank-local loop lists, per-label lifetime, uploads,
  explicit synchronization and downloads.
- `parallel/Parallel_LRI_Equally.hpp`: two-index task split and tensor requests.
- `distribute/Distribute_Equally.hpp`, `Split_Processes.hpp`: index partitioning.

Actual tested controls, generated inputs, outputs, memory traces and summarized
comparisons are retained in `/mnt/g/work/agent/abacus/lcao-hse-results/scaling/`.
`control/hse_benchmark.py` is the tested matrix driver (Si16 and Si64, 1/2/4
ranks). Input generation was checked for cell dimensions, atom counts and unique
positions. No production code changes were required for this experiment.
The runner shipped source HEAD `784610dfadc952bfbefd98c298e250814c0e9684` and
reused the previously built GPU RI executable (runtime C++ unchanged).
Remote archive beneath the experiment project root:
`archives/lcao-hse/lcao-scaling-20261017-1.tar.gz`.
