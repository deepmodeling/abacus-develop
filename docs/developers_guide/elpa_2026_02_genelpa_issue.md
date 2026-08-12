# ELPA 2026.02.001 `genelpa` MPI/grid compatibility report

This developer note records a reproducible external-solver issue found while
validating the ABACUS i-PI socket calculator. It is separate from the socket
implementation PRs: the socket code does not cause this failure.

## Summary

The installed ELPA 2026.02.001 build reproducibly fails for a minimal public
ELPA C API program using a complex-double, two-stage eigensolver on a row-major
`1 x 2` BLACS process grid. The same program passes with ELPA 2026.02.001 on a
`2 x 1` grid and with ELPA 2025.06.001 on both grids. The one-stage solver
passes in all tested ELPA 2026.02.001 cases.

The evidence is consistent with an ELPA 2026.02.001 regression or
build-specific incompatibility in this multi-process/grid configuration. It is
not yet sufficient to claim that an upstream ELPA bug has been confirmed; the
reproducer should be checked against another 2026.02.001 build or by ELPA
maintainers.

## ABACUS configuration

The failure was first observed in an ABACUS LCAO calculation with:

```text
nlocal = 26
nbands = 14
nb2d = 1
MPI ranks = 2
row-major BLACS grid = 1 x 2
complex<double> Hamiltonian path
ks_solver = genelpa
```

For two MPI ranks, ABACUS's row-major `Parallel_2D` setup produces the `1 x 2`
grid. For matrices smaller than 500 orbitals, the default LCAO block size is
`nb2d = 1`.

Relevant code paths are:

```text
source/source_lcao/LCAO_init_basis.cpp
source/source_base/parallel_2d.cpp
source/source_hsolver/module_genelpa/elpa_new.cpp
source/source_hsolver/module_genelpa/elpa_new_complex.cpp
```

## Standalone reproducer

The reproducer calls the ELPA C API directly and does not call ABACUS or the
i-PI socket layer. It uses a small positive-definite Hermitian complex matrix:

```text
global matrix size: n = 26
number of eigenvectors: nev = 14
block size: nblk = 1
MPI ranks: 2
BLACS layout: row-major
process grids: 1 x 2 and 2 x 1
```

The temporary reproducer was `/tmp/elpa_order_test.cpp`. The essential command
was:

```bash
OMP_NUM_THREADS=1 \
OMPI_MCA_pml=ob1 \
OMPI_MCA_btl=self,tcp \
mpirun --oversubscribe -np 2 ./elpa_order_test \
  two_before b1 1x2 26 1
```

The test was repeated with solver options set before and after `elpa_setup()`,
with ELPA 2026.02.001 and 2025.06.001, and with several complex kernels.

## Result matrix

`RC=0` means that the program completed and ELPA returned `ELPA_OK` on all
ranks. Non-zero results are launcher or process failures.

| ELPA version | BLACS grid | one-stage | two-stage |
| --- | ---: | ---: | ---: |
| 2026.02.001 | `1 x 2` | pass (`RC=0`) | fail |
| 2026.02.001 | `2 x 1` | pass (`RC=0`) | pass (`RC=0`) |
| 2025.06.001 | `1 x 2` | pass (`RC=0`) | pass (`RC=0`) |
| 2025.06.001 | `2 x 1` | pass (`RC=0`) | pass (`RC=0`) |

The ELPA 2026.02.001 `1 x 2` two-stage failure was reproduced with:

```text
ELPA_2STAGE_COMPLEX_AVX512_BLOCK1
ELPA_2STAGE_COMPLEX_AVX512_BLOCK2
ELPA_2STAGE_COMPLEX_GENERIC
ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE
```

Setting `solver` and `complex_kernel` before or after `elpa_setup()` did not
change the result. Linking the same reproducer with OpenMPI 5.0.8 instead of
OpenMPI 5.0.10 also reproduced the failure.

## Observed errors

With OpenMPI OB1/TCP, representative output is:

```text
mca_btl_tcp_frag_send: writev error
Bad address(3)
An error occurred in Socket closed
MPI_ERRORS_ARE_FATAL
```

With default UCX, representative runs reported UCX/xpmem errors followed by
`SIGBUS` or `SIGSEGV`. The stack reached:

```text
MPI_Bcast
  -> trans_ev_tridi_to_band_complex_double
  -> elpa_solve_evp_complex_2stage
  -> elpa_eigenvectors
```

Some runs also reported:

```text
Operating system error: Cannot allocate memory
Integer overflow in xmallocarray
```

These are MPI/ELPA internal errors, not i-PI protocol socket errors.

## Current attribution

The experiments support these conclusions:

1. The issue is reproducible without ABACUS and without the i-PI socket layer.
2. It is not explained by the `gint_precision` input keyword. That keyword is
   unsupported by the LTS 3.10.1 parser, but removing it does not fix the
   develop-module failure.
3. It is not only an UCX/xpmem problem: OB1/TCP runs fail as well.
4. It is not caused by one selected complex kernel or by setting the solver
   option before versus after `elpa_setup()`.
5. The version/grid/algorithm matrix points to an ELPA 2026.02.001 two-stage,
   multi-process compatibility problem.

The safe issue wording is:

> ELPA 2026.02.001 reproducibly fails for complex-double two-stage
> diagonalization on an ABACUS-compatible `1 x 2` BLACS grid; ELPA 2025.06.001
> passes the same reproducer.

Do not state that an upstream ELPA bug is confirmed until another 2026.02.001
build or ELPA maintainers reproduce the issue.

## ABACUS wrapper observations

The ABACUS `module_genelpa` wrapper has separate engineering issues that
should be handled in a solver-focused change, not hidden inside the socket
calculator PR:

- `elpa_set()` and `elpa_setup()` return codes are not consistently checked;
- some `MPI_Allreduce()` results are subsequently overwritten by local `info`;
- ELPA handle and `elpa_init()`/`elpa_uninit()` lifetimes need a focused audit;
- failures should become rank-consistent, actionable ABACUS diagnostics.

These observations improve ABACUS diagnostics, but do not explain why a
standalone public-API program reproduces the ELPA failure.

## Temporary workarounds

Until the dependency issue is resolved, use one of the following for affected
CPU LCAO runs:

```text
ELPA 2025.06.001
ks_solver = elpa
ks_solver = scalapack_gvx
```

Socket validation should use a solver configuration known to work on the
target machine. Socket PRs should not silently override ELPA behavior.

## Follow-up

1. Submit the reproducer and result matrix to ELPA maintainers.
2. Repeat with an independent ELPA 2026.02.001 CPU-only build and on a
   GPU-capable compute node.
3. If confirmed externally, document the affected ELPA/build/grid combination
   as a platform dependency limitation.
4. Independently harden the ABACUS genELPA wrapper with checked return codes
   and collective error propagation.
