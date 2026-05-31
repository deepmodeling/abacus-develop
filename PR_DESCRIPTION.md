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

A standalone microbenchmark (`bench_comm.cpp`) was written to compare the blocking
(`MPI_Alltoallv`) and nonblocking (`MPI_Isend/Irecv + MPI_Waitsome`) communication patterns
using the exact same algorithms as `pw_gatherscatter.h`.

### Methodology

- Tested across **grid sizes** (32³, 48³, 64³), **MPI ranks** (2, 3, 4), and **OpenMP threads** (1, 2)
- 18 total configurations, 500 iterations each after warmup
- Measured both `gatherp_scatters` (forward FFT) and `gathers_scatterp` (reverse FFT)

### Results

| Grid | MPI=2 | MPI=3 | MPI=4 |
|------|-------|-------|-------|
| 32³  | **1.02x** | **1.13x** | **1.22x** |
| 48³  | **1.06x** | **1.28x** | **1.39x** |
| 64³  | **1.11x** | **1.27x** | **1.42x** |

*Table: Roundtrip speedup (nonblocking/blocking) with OMP_NUM_THREADS=2.*

### Key findings

1. **Nonblocking is consistently faster** — all 18 configurations show nonblocking ≥ blocking
2. **Speedup scales with MPI parallelism** — more ranks → more peers to overlap: avg 1.07x (2p) → 1.21x (3p) → 1.29x (4p)
3. **Speedup scales with problem size** — larger grids → more unpack work to overlap: up to **1.42x**
4. **`gathers_scatterp` benefits most** — its zeroing step overlaps with in-flight MPI data

The benchmark source is included at `source/source_basis/module_pw/test/bench_comm.cpp` and can
be compiled and run with:
```bash
mpicxx -std=c++14 -O3 -fopenmp -o bench_comm bench_comm.cpp
OMP_NUM_THREADS=2 mpirun -np 4 ./bench_comm 500 64 64 64
```

🤖 Generated with [Claude Code](https://claude.com/claude-code)
