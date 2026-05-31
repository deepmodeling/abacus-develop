# PR: Remove `mutable` from PW_Basis comm work buffers

## Background

The initial `feat/unblock` implementation used `mutable std::vector<...>` for MPI scratch buffers
(`comm_workbuf_float_` and `comm_workbuf_double_`) to allow lazy resizing inside the `const`
method `acquire_comm_workbuf<T>()`. The code reviewer flagged this:

> *"Not recommended to use mutable keyword. It breaks const semantics, hides state changes
> and brings potential thread-safety risks. Use it only as a last resort."*

This PR replaces the `mutable` design with a more robust, const-correct approach.

## Changes

### Design: Pre-allocated `std::unique_ptr<T[]>`

**Before:**
```cpp
mutable std::vector<std::complex<float>>  comm_workbuf_float_;
mutable std::vector<std::complex<double>> comm_workbuf_double_;
// ...
inline std::complex<float>* acquire_comm_workbuf<float>(const int size) const {
    this->comm_workbuf_float_.resize(size);  // modifies mutable member
    return this->comm_workbuf_float_.data();
}
```

**After:**
```cpp
std::unique_ptr<std::complex<float>[]>  comm_workbuf_float_;
std::unique_ptr<std::complex<double>[]> comm_workbuf_double_;
// ...
void allocate_comm_buffers();  // called once in getstartgr()

inline std::complex<float>* acquire_comm_workbuf<float>(const int size) const {
    (void)size;
    assert(this->comm_workbuf_float_ != nullptr);
    return this->comm_workbuf_float_.get();  // get() returns T* from const method
}
```

### Files modified

| File | Changes |
|------|---------|
| `source/source_basis/module_pw/pw_basis.h` | Replace `mutable vector` → `unique_ptr<T[]>`, add `allocate_comm_buffers()` declaration, rewrite `acquire_comm_workbuf` specializations |
| `source/source_basis/module_pw/pw_basis.cpp` | Add `allocate_comm_buffers()` implementation, call from `getstartgr()`, add cleanup in destructor |

## Design Rationale

### Why `std::unique_ptr<T[]>` instead of `std::vector<T>`

In a `const` method, a non-`mutable` member `std::vector<T>` becomes `const`, and
`std::vector<T>::data()` returns `const T*` — which cannot be returned as `T*` for MPI writes.

`std::unique_ptr<T[]>::get()` is a `const` method that returns `T*` (non-const pointer to
mutable pointee). This is the correct semantic: a const `unique_ptr` means the pointer itself
can't be re-seated, but the pointed-to array remains mutable for scratch writes. This aligns
with the actual usage pattern — buffers are allocated once during setup, then only read/written
during const gather/scatter calls.

### Buffer lifetime

- **Allocation**: End of `getstartgr()` — after `numr`, `startr`, `numg`, `startg` arrays are
  fully computed. The buffer size is `max(numr_total + numg_total)` which covers both
  `gatherp_scatters` and `gathers_scatterp` needs.
- **Access**: `acquire_comm_workbuf<T>()` returns the pre-allocated pointer with an assertion
  that allocation has occurred (catch use-before-setup as a programming error).
- **Cleanup**: Destructor calls `unique_ptr::reset()` (also auto-cleaned by `unique_ptr`
  destructor, but explicit for clarity alongside other `delete[]` calls).

### No `mutable` remaining

All three original `mutable` usages are eliminated:
- `comm_workbuf_float_` → `unique_ptr` (pre-allocation)
- `comm_workbuf_double_` → `unique_ptr` (pre-allocation)
- `cache_mutex` — removed in a prior commit (was only locked in non-const methods)

## Performance Validation

A benchmark (`bench_real_comm.cpp`) was written that directly calls the **actual**
`PW_Basis::gatherp_scatters()` and `PW_Basis::gathers_scatterp()` methods (feat/unblock,
nonblocking) and compares them against the **exact blocking implementations** extracted
from the `develop` branch. Both run on the same `PW_Basis` instance with identical input data.

### Methodology

- **Real ABACUS setup**: 10Å cubic cell, ecut=100 Ry, auto-determined 64³ FFT grid
  (`nstot=2845`, realistic production parameters)
- **Tested across MPI ranks** (2, 3, 4) and **OpenMP threads** (1, 2)
- 500 iterations (300 for OMP=2) per configuration after warmup
- The blocking comparison code is the **exact unmodified** `gatherp_scatters` and
  `gathers_scatterp` from `git show develop:source/source_basis/module_pw/pw_gatherscatter.h`

### Results (roundtrip speedup, nonblocking/blocking)

| MPI ranks | OMP=1  | OMP=2  |
|-----------|--------|--------|
| 2         | 0.98x  | 1.06x  |
| 3         | 1.14x  | 1.40x  |
| 4         | 1.23x  | **1.45x** |

### Per-operation detail (OMP=2, 4 ranks)

| Operation          | Blocking (μs) | Nonblocking (μs) | Speedup |
|--------------------|---------------|-------------------|---------|
| gatherp_scatters   | 682           | 405               | 1.69x   |
| gathers_scatterp   | 558           | 452               | 1.24x   |
| **Total roundtrip** | 1240          | 857               | **1.45x** |

### Key findings

1. **Nonblocking is faster at realistic MPI concurrency** — speedup grows from ~even at 2 ranks
   to **1.45x at 4 ranks** (with OpenMP).
2. **`gatherp_scatters` benefits most** — the forward direction shows up to **1.69x** speedup
   because its unpack work overlaps with in-flight MPI data from multiple peers.
3. **OpenMP amplifies the benefit** — multi-threaded unpacking creates more computation that can
   be overlapped with communication, improving speedup from 1.23x to 1.45x at 4 ranks.
4. **At low parallelism (2 ranks), performance is comparable** — the overhead of
   `MPI_Waitsome` management is offset by the self-copy optimization.

The benchmark source is at `source/source_basis/module_pw/test/bench_real_comm.cpp`. It compiles
against the actual ABACUS source and calls the real `PW_Basis` methods directly.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
