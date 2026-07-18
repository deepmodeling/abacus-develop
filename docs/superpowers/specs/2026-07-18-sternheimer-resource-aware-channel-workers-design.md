# Resource-Aware Sternheimer Channel Workers

## Status

Approved design for replacing the fixed production worker cap used by the 20 Angstrom molecular Delta-Sternheimer run. The current `ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS=16` jobs remain unchanged and continue to completion; this design governs later runs.

## Problem

Auxiliary-potential channels are independent and can be solved concurrently, but every active channel owns a large GMRES workspace proportional to the real-space grid size. A fixed worker count is therefore not portable across:

- cell size and `ecutwfc`, which determine the uniform-grid size;
- the OpenMP thread count;
- the memory requested from Slurm or enforced by a cgroup;
- the number of MPI ranks placed on the same node;
- standard Sternheimer and Delta-Sternheimer fixed working sets.

The existing value `16` is only an empirically safe cap for the current 20 Angstrom, 25 Ry, one-frequency-rank-per-node jobs with 110610 MiB requested per node. It must not become a general default. An unset or zero environment variable currently exposes all OpenMP threads, which can reproduce the original out-of-memory failure.

## Selected Approach

ABACUS will calculate the effective channel worker count after the zero-order and, when enabled, shared Delta fixed subspaces have been constructed. The planner will combine a measured memory baseline with an analytical per-worker workspace estimate.

Script-only resource selection is rejected because it is Slurm-specific and duplicates physics-derived grid information outside ABACUS. Runtime trial-and-backoff is also rejected: a first wave that is too large can be killed by the cgroup before ABACUS can reduce concurrency.

## Worker-Count Rule

For each spin response channel, prefer node-aggregate accounting when the cgroup exposes both its limit and current usage:

```text
node_target         = floor(0.75 * node_memory_limit)
node_increment      = node_target - node_memory_current
rank_increment      = floor(node_increment / local_mpi_ranks)
memory_workers      = floor(rank_increment / memory_per_worker)
effective_workers   = min(channels, omp_threads, memory_workers, user_cap)
```

When only a node limit and the current process RSS are available, use per-rank accounting instead:

```text
rank_limit          = floor(node_memory_limit / local_mpi_ranks)
rank_target         = floor(0.75 * rank_limit)
rank_increment      = rank_target - process_current_rss
memory_workers      = floor(rank_increment / memory_per_worker)
effective_workers   = min(channels, omp_threads, memory_workers, user_cap)
```

When `/proc/meminfo` is the only source, `MemAvailable` is already a remaining node-wide amount and is not treated as a hard limit or reduced by RSS again:

```text
rank_increment      = floor(0.75 * mem_available / local_mpi_ranks)
memory_workers      = floor(rank_increment / memory_per_worker)
effective_workers   = min(channels, omp_threads, memory_workers, user_cap)
```

`user_cap` is omitted from the final minimum when the environment variable is unset or zero. All arithmetic must use checked 64-bit byte counts. The effective count must be at least one. If the detected memory budget cannot accommodate one worker within the 75 percent target, ABACUS stops before entering the channel loop and reports the measured and estimated byte counts.

The 25 percent reserve covers the Hamiltonian, ABFS potentials, zero-order states, shared occupied/Delta subspaces, MPI/runtime allocations, output buffers, and estimator error. It is deliberately a fraction of the enforced memory limit rather than a fixed GiB reserve so that the rule scales to other nodes.

## Resource Detection

The detector returns a node-local memory snapshot with a source label, limit, current usage, and local MPI rank count.

1. Read the process cgroup path from `/proc/self/cgroup`.
2. Prefer a finite cgroup-v2 `memory.max` and its matching `memory.current`.
3. Support cgroup-v1 `memory.limit_in_bytes` and `memory.usage_in_bytes` as a compatibility path.
4. Treat an unlimited cgroup value, including v2 `max` and known v1 near-physical-memory sentinel values, as unavailable.
5. Use `SLURM_MEM_PER_NODE` as an additional finite limit when present. The effective limit is the smallest trustworthy finite limit.
6. Use `/proc/meminfo` `MemAvailable` only as a fallback when no enforced cgroup/Slurm limit is available.
7. Determine ranks sharing a node with `MPI_Comm_split_type(..., MPI_COMM_TYPE_SHARED, ...)`; use one in non-MPI builds.

When cgroup `memory.current` is available it already includes all processes in the allocation on that node. Its remaining increment budget is divided by the local MPI rank count. For `/proc/self/status` `VmRSS` fallback, the planner first derives a per-rank memory limit and subtracts the current process RSS.

If no trustworthy limit or available-memory value can be detected, the safe fallback is one channel worker. ABACUS emits an explicit diagnostic rather than silently using every OpenMP thread.

## Per-Worker Estimate

The default GMRES restart dimension is 50. At peak, the solver owns 51 Krylov basis vectors and 50 preconditioned vectors, plus residual, operator, solution, projected RHS, apply-callback, reconstruction, and accumulation temporaries. The initial conservative estimate is

```text
memory_per_worker = 120 * grid_size * sizeof(complex<double>)
```

This is 3.49 GiB per worker for the 20 Angstrom, 25 Ry grid of 1,953,125 points. The coefficient arrays and Hessenberg matrix are negligible at this scale but will be included with checked fixed-size terms for completeness.

The factor 120 is an implementation constant named for its purpose and covered by a test that relates it to the current GMRES restart dimension. If the solver workspace ownership changes, the estimate and its test must change together. It is not exposed as a user input in this first implementation.

## Environment Variable Semantics

`ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS` becomes a safety cap:

- unset or `0`: use automatic resource-aware selection;
- positive integer: use the smaller of the automatic count and this value;
- negative, malformed, or overflowing value: fail with an input error.

The variable cannot force ABACUS above the memory-derived limit. No separate force-unsafe option is added.

## Code Boundaries

Resource detection and worker planning will be separated from the OpenMP scheduler:

- a pure planner accepts channels, OpenMP threads, grid size, memory limit/current usage, local rank count, and optional user cap;
- a Linux resource detector supplies the runtime snapshot;
- the existing channel scheduler receives only the resulting positive worker count;
- the Sternheimer production path invokes the planner once per spin after its fixed working set is ready.

This keeps operating-system parsing out of the scheduler and makes the numerical selection deterministic in unit tests.

## Diagnostics

A new `channel_workers_ready` progress event will include:

```text
resource_source=<source>
accounting_mode=<node_aggregate|per_rank|available|fallback_one>
node_memory_limit_bytes=<N>
memory_current_bytes=<N>
memory_target_bytes=<N>
local_mpi_ranks=<N>
grid_size=<N>
memory_per_worker_bytes=<N>
automatic_workers=<N>
user_cap=<N>
effective_workers=<N>
```

The existing `channels_ready` event continues to describe the generated auxiliary channels. It must not report a final worker count before the fixed-subspace baseline exists.

## Testing

### Unit tests

- enough memory selects all available OpenMP workers;
- channel count and OpenMP count clamp the result;
- a positive user cap can only reduce the automatic count;
- unset and zero user caps are equivalent automatic modes;
- multiple local MPI ranks share the node increment budget;
- memory below the one-worker target fails before channel execution;
- unavailable resource detection falls back to one worker;
- checked arithmetic rejects overflow;
- cgroup-v2, cgroup-v1, Slurm, and `/proc` parser fixtures select the expected source and limit.

### Remote regression

On `df_dcu`:

1. run the complete Sternheimer unit-test set;
2. run the existing small reader-v1 case in automatic mode and with explicit caps 1 and 4;
3. require matrix differences to remain at the linear-solver tolerance scale;
4. verify the progress event reports the expected automatic and effective counts;
5. compare the estimate with Slurm `MaxRSS` and confirm the job remains below its requested memory.

The already-running 20 Angstrom jobs are not restarted. Their final peak RSS and wall time will be recorded as calibration evidence. A later 20 Angstrom production rerun is required only when another scientific convergence point is needed, not solely to retest this scheduler change.

## Documentation and Compatibility

The TeX project will record the formula, detection hierarchy, diagnostics, validation results, and the distinction between the historical explicit-16 run and future automatic runs. Existing scripts that set a positive cap remain valid but can no longer exceed the automatic memory limit. Scripts that omit the variable become safer because zero now means automatic selection instead of all OpenMP threads.
