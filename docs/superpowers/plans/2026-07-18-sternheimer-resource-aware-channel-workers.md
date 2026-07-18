# Resource-Aware Sternheimer Channel Workers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fixed Delta-Sternheimer auxiliary-channel worker cap with a safe automatic count derived from the real-space grid, actual node/rank memory resources, MPI placement, OpenMP capacity, and an optional user upper bound.

**Architecture:** Add a focused resource module with a pure worker planner and a Linux runtime detector. The detector reports one of four explicit accounting modes; the planner performs checked byte arithmetic and applies the 75 percent memory target. The production Sternheimer path measures resources after each spin's fixed subspace is ready, logs the complete decision, and passes only the positive effective count to the existing OpenMP scheduler.

**Tech Stack:** C++17, OpenMP, MPI-3 shared-memory communicators, Linux cgroup v1/v2 and `/proc`, GoogleTest, ABACUS Delta-ST, Slurm on `df_dcu`.

---

## File Structure

- Create `source/source_lcao/module_ri/sternheimer_channel_resources.h`: public data types, pure planning API, parser API used by tests, and runtime detector declaration.
- Create `source/source_lcao/module_ri/sternheimer_channel_resources.cpp`: checked arithmetic, resource parser/selection logic, Linux file probing, Slurm limit handling, and MPI node-local synchronization.
- Create `source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp`: deterministic planner and parser tests without relying on the test process's live cgroup.
- Modify `source/source_lcao/module_ri/CMakeLists.txt`: compile the new resource implementation into the RI object library.
- Modify `source/source_lcao/module_ri/test/CMakeLists.txt`: register the focused resource test executable.
- Modify `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: change the environment variable to an upper-cap input, calculate an effective count once per spin, add diagnostics, and pass the result to the channel scheduler.
- Modify `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`: verify automatic/capped diagnostic formatting and preserve scheduler behavior.
- Modify `sternheimer_siab_project/main.tex` outside the ABACUS worktree after remote validation: record the formula, semantics, tests, and the completed 20 Angstrom calibration run.

---

### Task 1: Pure Memory Planner

**Files:**
- Create: `source/source_lcao/module_ri/sternheimer_channel_resources.h`
- Create: `source/source_lcao/module_ri/sternheimer_channel_resources.cpp`
- Create: `source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Define planner inputs and write failing tests**

Define these types in the new header so later runtime detection cannot change planning semantics:

```cpp
namespace ModuleRI
{
enum class SternheimerMemoryAccountingMode
{
    node_aggregate,
    per_rank,
    available,
    fallback_one
};

struct SternheimerMemorySnapshot
{
    SternheimerMemoryAccountingMode mode = SternheimerMemoryAccountingMode::fallback_one;
    std::uint64_t limit_bytes = 0;
    std::uint64_t current_bytes = 0;
    int local_mpi_ranks = 1;
    std::string source = "unavailable";
};

struct SternheimerChannelWorkerPlan
{
    int automatic_workers = 1;
    int effective_workers = 1;
    std::uint64_t target_bytes = 0;
    std::uint64_t increment_bytes_per_rank = 0;
    std::uint64_t memory_per_worker_bytes = 0;
};

std::uint64_t estimate_sternheimer_channel_worker_bytes(std::size_t grid_size);

SternheimerChannelWorkerPlan plan_sternheimer_channel_workers(
    int num_channels,
    int omp_threads,
    std::size_t grid_size,
    int user_cap,
    const SternheimerMemorySnapshot& memory);
}
```

Add tests with synthetic byte counts for these exact behaviors:

```cpp
TEST(SternheimerChannelResources, EstimatesOneHundredTwentyComplexGridVectors)
{
    EXPECT_EQ(ModuleRI::estimate_sternheimer_channel_worker_bytes(1000),
              120ULL * 1000ULL * sizeof(std::complex<double>));
}

TEST(SternheimerChannelResources, PlansFromNodeAggregateMemory)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        100000, 32000, 2, "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 11);
    EXPECT_EQ(plan.effective_workers, 11);
}

TEST(SternheimerChannelResources, UserCapOnlyReducesAutomaticCount)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000, 0, 1, "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 4, memory);
    EXPECT_EQ(plan.effective_workers, 4);
}

TEST(SternheimerChannelResources, ClampsToChannelsAndOpenMPThreads)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000, 0, 1, "proc_meminfo"};
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(3, 30, 1, 0, memory).effective_workers, 3);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 5, 1, 0, memory).effective_workers, 5);
}

TEST(SternheimerChannelResources, DividesPerRankLimitBeforeSubtractingRss)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::per_rank,
        100000, 10000, 2, "slurm+proc_status"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 14);
}

TEST(SternheimerChannelResources, ZeroCapUsesAutomaticCount)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000, 0, 1, "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 0, memory);
    EXPECT_EQ(plan.effective_workers, plan.automatic_workers);
}

TEST(SternheimerChannelResources, RejectsDetectedBudgetBelowOneWorker)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        1000, 0, 1, "cgroup_v2"};
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 0, memory), std::runtime_error);
}

TEST(SternheimerChannelResources, UnavailableResourcesFallBackToOneWorker)
{
    const ModuleRI::SternheimerMemorySnapshot memory;
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 0, memory).effective_workers, 1);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 8, memory).effective_workers, 1);
}

TEST(SternheimerChannelResources, RejectsInvalidInputsAndOverflow)
{
    const ModuleRI::SternheimerMemorySnapshot memory;
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(-1, 1, 1, 0, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(1, 0, 1, 0, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(1, 1, 1, -1, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::estimate_sternheimer_channel_worker_bytes(
                     std::numeric_limits<std::size_t>::max()),
                 std::overflow_error);
}
```

The tests above are the required planner cases; keep their byte values unchanged so the expected integer divisions remain auditable.

- [ ] **Step 2: Register and run the new test to verify RED**

Add a target that names the not-yet-implemented source:

```cmake
AddTest(
  TARGET MODULE_RI_sternheimer_channel_resources_test
  LIBS ${math_libs}
  SOURCES sternheimer_channel_resources_test.cpp ../sternheimer_channel_resources.cpp
)
```

Transfer the header, test, and CMake edit to `/public/home/ghj/app/src/abacus-delta-st-postprocess-20260705`, configure/build the target in `build_delta_st_test`, and require failure because the planner implementation is absent. Record the compiler/linker error as the RED evidence.

- [ ] **Step 3: Implement checked planning arithmetic**

Implement the three accounting equations exactly as specified:

```cpp
constexpr std::uint64_t kComplexVectorsPerWorker = 120;
constexpr std::uint64_t kMemoryTargetNumerator = 3;
constexpr std::uint64_t kMemoryTargetDenominator = 4;

std::uint64_t estimate_sternheimer_channel_worker_bytes(const std::size_t grid_size)
{
    return checked_multiply(checked_multiply(kComplexVectorsPerWorker, grid_size),
                            sizeof(std::complex<double>));
}
```

For `node_aggregate`, compute `floor(0.75*limit)-current`, then divide by local ranks. For `per_rank`, divide the limit by local ranks first, apply 0.75, then subtract process RSS. For `available`, apply 0.75 to available bytes and divide by local ranks without subtracting RSS. For `fallback_one`, return one worker before memory arithmetic. Throw `std::runtime_error` when a detected budget cannot hold one worker and `std::invalid_argument` for malformed logical inputs.

- [ ] **Step 4: Run the focused test to verify GREEN**

Build and run `MODULE_RI_sternheimer_channel_resources_test` on `df_dcu`. Require every planner test to pass, including overflow and low-memory errors.

- [ ] **Step 5: Commit the pure planner**

Commit only the new planner, focused test, and test CMake edit:

```bash
git add source/source_lcao/module_ri/sternheimer_channel_resources.{h,cpp} \
        source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp \
        source/source_lcao/module_ri/test/CMakeLists.txt
git commit -m "feat(sternheimer): plan channel workers from memory"
```

Use `Codex <codex@openai.com>` as author and `AroundPeking <gonghuanjing@iphy.ac.cn>` as committer.

---

### Task 2: Linux, Slurm, and MPI Resource Detection

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_channel_resources.h`
- Modify: `source/source_lcao/module_ri/sternheimer_channel_resources.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`

- [ ] **Step 1: Write parser and selection tests before detector code**

Add these parser declarations under `ModuleRI::detail`:

```cpp
struct SternheimerMemoryCandidates
{
    std::optional<std::uint64_t> cgroup_limit_bytes;
    std::optional<std::uint64_t> cgroup_current_bytes;
    std::optional<std::uint64_t> slurm_limit_bytes;
    std::optional<std::uint64_t> mem_available_bytes;
    std::optional<std::uint64_t> process_rss_bytes;
    std::optional<std::uint64_t> physical_memory_bytes;
    std::string cgroup_source;
};

std::optional<std::uint64_t> parse_sternheimer_memory_bytes(std::string_view text);
std::optional<std::uint64_t> parse_sternheimer_slurm_mem_per_node(std::string_view text);
std::optional<std::uint64_t> parse_sternheimer_kib_field(std::string_view text, std::string_view key);
std::optional<std::string> parse_sternheimer_cgroup_v2_path(std::string_view text);
std::optional<std::string> parse_sternheimer_cgroup_v1_memory_path(std::string_view text);
SternheimerMemorySnapshot select_sternheimer_memory_snapshot(
    const SternheimerMemoryCandidates& candidates,
    int local_mpi_ranks);
```

Add deterministic tests for:

```cpp
EXPECT_EQ(parse_sternheimer_memory_bytes("110610M"), 110610ULL * 1024ULL * 1024ULL);
EXPECT_EQ(parse_sternheimer_memory_bytes("4096"), 4096ULL);
EXPECT_EQ(parse_sternheimer_slurm_mem_per_node("110610"), 110610ULL * 1024ULL * 1024ULL);
EXPECT_FALSE(parse_sternheimer_memory_bytes("max").has_value());
EXPECT_EQ(parse_sternheimer_kib_field("MemAvailable: 2048 kB\n", "MemAvailable"), 2048ULL * 1024ULL);
EXPECT_EQ(parse_sternheimer_kib_field("VmRSS:\t1024 kB\n", "VmRSS"), 1024ULL * 1024ULL);
EXPECT_EQ(parse_sternheimer_cgroup_v2_path("0::/slurm/job_7/step_0\n"), "/slurm/job_7/step_0");
EXPECT_EQ(parse_sternheimer_cgroup_v1_memory_path("5:memory:/slurm/job_7\n"), "/slurm/job_7");
```

Add these exact selection cases:

```cpp
TEST(SternheimerChannelResources, SelectsSmallestEnforcedLimitWithCgroupUsage)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.cgroup_limit_bytes = 100000;
    values.cgroup_current_bytes = 20000;
    values.slurm_limit_bytes = 90000;
    values.cgroup_source = "cgroup_v2";
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 2);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::node_aggregate);
    EXPECT_EQ(result.limit_bytes, 90000U);
    EXPECT_EQ(result.current_bytes, 20000U);
    EXPECT_EQ(result.local_mpi_ranks, 2);
}

TEST(SternheimerChannelResources, UsesPerRankRssWithSlurmLimit)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.slurm_limit_bytes = 100000;
    values.process_rss_bytes = 10000;
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 2);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::per_rank);
}

TEST(SternheimerChannelResources, UsesAvailableMemoryWithoutSubtractingRssAgain)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.mem_available_bytes = 80000;
    values.process_rss_bytes = 10000;
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 1);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::available);
    EXPECT_EQ(result.limit_bytes, 80000U);
}

TEST(SternheimerChannelResources, TreatsV1PhysicalMemorySentinelAsUnlimited)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.cgroup_limit_bytes = 0x7ffffffffffff000ULL;
    values.physical_memory_bytes = 128ULL << 30;
    values.mem_available_bytes = 64ULL << 30;
    values.cgroup_source = "cgroup_v1";
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 1);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::available);
}
```

- [ ] **Step 2: Run the parser tests to verify RED**

Transfer the test-only changes and confirm the focused target fails because the parser and selector functions are missing.

- [ ] **Step 3: Implement resource parsing and source selection**

Add parser declarations under `ModuleRI::detail` and keep live file probing outside the pure planner. Implement strict suffix handling for raw cgroup bytes and explicit `K`, `M`, `G` suffixes, checked conversion to bytes, exact key matching for `/proc` fields, and cgroup path extraction. Parse a suffix-free `SLURM_MEM_PER_NODE` value with the dedicated MiB parser; never feed it to the raw-byte cgroup parser.

Declare the runtime API:

```cpp
SternheimerMemorySnapshot detect_sternheimer_memory_snapshot();
std::string sternheimer_memory_accounting_mode_name(SternheimerMemoryAccountingMode mode);
```

The detector must:

1. synchronize ranks sharing a physical node with `MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, ...)` under `__MPI`;
2. let the node-local root read cgroup files and environment values, then broadcast the selected finite node limit/current usage;
3. use a node-local maximum of process `VmRSS` for the per-rank fallback;
4. set `local_mpi_ranks` on every rank and free the temporary communicator;
5. use one local rank in non-MPI builds;
6. return `fallback_one` with a diagnostic source when no trustworthy source exists.

- [ ] **Step 4: Compile the detector into the RI object library**

Add `sternheimer_channel_resources.cpp` next to `sternheimer_rpa.cpp` in `source/source_lcao/module_ri/CMakeLists.txt`. Rebuild the focused test and the RI object target on `df_dcu`.

- [ ] **Step 5: Verify GREEN and commit the detector**

Run the focused parser/planner test and inspect one compute-node diagnostic probe to confirm the live source resolves to cgroup or Slurm rather than login-node physical memory. Commit:

```bash
git add source/source_lcao/module_ri/sternheimer_channel_resources.{h,cpp} \
        source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp \
        source/source_lcao/module_ri/CMakeLists.txt
git commit -m "feat(sternheimer): detect channel memory resources"
```

---

### Task 3: Production Sternheimer Integration

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] **Step 1: Write failing diagnostic and cap-semantics tests**

Add this formatting helper declaration to `sternheimer_channel_resources.h`:

```cpp
std::string format_sternheimer_channel_worker_diagnostic(
    const SternheimerMemorySnapshot& memory,
    const SternheimerChannelWorkerPlan& plan,
    std::size_t grid_size,
    int user_cap);
```

Test the exact decision string:

```cpp
TEST(SternheimerChannelResources, FormatsWorkerDecisionDiagnostic)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        100000, 32000, 2, "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 4, memory);
    EXPECT_EQ(ModuleRI::format_sternheimer_channel_worker_diagnostic(memory, plan, 1, 4),
              "resource_source=cgroup_v2 accounting_mode=node_aggregate "
              "node_memory_limit_bytes=100000 memory_current_bytes=32000 "
              "memory_target_bytes=75000 local_mpi_ranks=2 grid_size=1 "
              "memory_per_worker_bytes=1920 automatic_workers=11 user_cap=4 effective_workers=4");
}
```

The pure planner test from Task 1 already proves that a positive cap cannot increase the automatic count. Keep the existing OpenMP scheduler test proving that a positive effective count is honored; do not duplicate that behavior in the production-path test.

- [ ] **Step 2: Run the affected test targets to verify RED**

Build/run `MODULE_RI_sternheimer_channel_resources_test` and `MODULE_RI_sternheimer_abacus_st_smoke_test`; require the missing formatter/integration assertion to fail before editing production code.

- [ ] **Step 3: Integrate automatic selection once per spin**

In `sternheimer_abacus_st_smoke.cpp`:

1. rename the parsed value to `channel_worker_user_cap`;
2. keep unset/zero as automatic mode and reject negative/malformed values;
3. change `channels_ready` metadata from `max_workers` to `user_cap` because the memory baseline does not exist yet;
4. after the occupied projector and optional shared Delta fixed subspace are complete, call `detect_sternheimer_memory_snapshot()`;
5. call `plan_sternheimer_channel_workers(num_channels, omp_get_max_threads(), grid_data.grid.size(), channel_worker_user_cap, snapshot)`;
6. append one `channel_workers_ready` event per spin using the formatter;
7. pass `plan.effective_workers`, never the raw environment value, to `run_sternheimer_channel_tasks`.

The selected count is immutable during one spin response so numerical ordering and exception behavior stay unchanged.

- [ ] **Step 4: Verify focused and complete unit tests**

Run both focused targets, then all Sternheimer test executables in `build_delta_st_test`. Require the same 10/10 baseline plus the new resource target, with no failures.

- [ ] **Step 5: Commit production integration**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp \
        source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp \
        source/source_lcao/module_ri/sternheimer_channel_resources.{h,cpp}
git commit -m "feat(sternheimer): select channel workers automatically"
```

---

### Task 4: Remote Reader-v1 Regression and Memory Validation

**Files:**
- Rebuild: `/public/home/ghj/app/src/abacus-delta-st-postprocess-20260705/build_delta_st_test/abacus_3p`
- Reuse inputs from: `/work1/ghj/sternheimer_abacus_tests/delta_st_memory_regression_20260718`
- Create: new automatic-mode and explicit-cap regression directories under the same root

- [ ] **Step 1: Rebuild an immutable remote binary**

Transfer only committed source/test/CMake files, rebuild `abacus_3p`, and record source commit, binary SHA-256, byte size, and modification time in each run script.

- [ ] **Step 2: Run automatic, cap-1, and cap-4 small cases**

Use the same 6 Angstrom, 10 Ry, single-frequency spin-polarized H input and identical `(1,1)` CCP settings. Remove the environment variable for automatic mode; set it to 1 and 4 in the two capped runs. Submit all jobs to `normal` with full-node memory and CPU resources.

- [ ] **Step 3: Compare structured reader-v1 outputs**

Read each response with `sternheimer_siab_project/scripts/compute_st_ecrpa.py::read_st_m_v1`. Require identical dimensions/frequencies/weights and relative Frobenius differences no larger than the current `1e-8` solver tolerance. Verify progress metadata reports `effective_workers <= automatic_workers` and the explicit caps never raise concurrency.

- [ ] **Step 4: Validate actual memory against the plan**

Collect `sacct`/`sstat` MaxRSS. Confirm every run stays below its requested memory and record the ratio of observed incremental memory to the 120-vector estimate. Do not tune the factor downward from one small case; increase it only if observed use violates the 75 percent target.

- [ ] **Step 5: Check the current 20 Angstrom chains**

Inspect jobs `21306782/21306784` and `21306788/21306790`. Require six complete reader-v1 responses per system, successful LibRPA termination, final peak RSS, wall time, and EcRPA before calling the production chain complete. Do not restart these explicit-16 jobs.

---

### Task 5: Documentation and Final Verification

**Files:**
- Modify: `/Users/ghj/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.tex`
- Rebuild: `/Users/ghj/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.pdf`

- [ ] **Step 1: Record implementation and validation**

Add the three accounting formulas, cgroup/Slurm detection hierarchy, 75 percent target, 120-vector estimate, environment upper-cap semantics, commit IDs, unit-test results, reader-v1 matrix differences, and observed MaxRSS.

- [ ] **Step 2: Record final 20 Angstrom results**

Add response and LibRPA job states, wall times, complete equation counts, final N2/N EcRPA, and the resulting binding energy using the already documented decomposition. Keep the historical explicit-16 result distinct from later automatic scheduling.

- [ ] **Step 3: Rebuild and inspect the PDF**

Run:

```bash
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Render the changed final pages and visually check equations, tables, long environment-variable names, and page boundaries.

- [ ] **Step 4: Final repository verification**

Verify the ABACUS worktree is clean, inspect all new commit author/committer fields, and report that the branch is ahead of its remote without pushing unless explicitly requested.
