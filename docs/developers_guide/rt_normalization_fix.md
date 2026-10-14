# Fix: normalize RT wavefunctions with the endpoint overlap

After half_Hmatrix averages H/S for propagation, evolve_psi currently passes the averaged S to normalization. A wavefunction arriving at t+dt needs its norm evaluated using S(t+dt). For a time-dependent overlap, a unit norm in the averaged metric need not be a unit norm at the endpoint.

The CPU path uses the unmodified endpoint s_mat.p. The tensor path reloads endpoint S from the existing distributed or root-gathered host buffer after propagation and before normalization. Band energies already use endpoint H. This also benefits moving-ion cases with changing S; constant-S results should be unchanged.

## Unit Tests and/or Case Tests
- `cmake --build /tmp/upstream-norm-review-20261014 --target MODULE_LCAO_tddft_evolve_psi_test MODULE_LCAO_tddft_norm_psi_test -j4`: passed.
- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ctest --test-dir /tmp/upstream-norm-review-20261014 -V -R 'MODULE_LCAO_tddft_(evolve_psi|norm_psi)_test'`: 2/2 CTests, 4/4 GTests passed. The new tests call the actual CPU and tensor-LAPACK propagation entries with a complex, positive-definite changing overlap.
- Negative regression check: compiled a separate copy of `evolve_psi.cpp` with the previous averaged-overlap normalization restored and linked it to the same tests. Both propagation tests failed with endpoint norm 1.1935483870967742; both pass with the fix. The reviewed worktree was unchanged during this check.
- git diff --check and staged agent_governance_check: passed (documentation review warning before this note was added).
- Full executable build, GPU/cuBLASMp execution, distributed multi-rank execution, and material trajectories: not run. The non-LAPACK tensor kernels require cuBLASMp and are not exercised by the CPU test target. The propagation tests use one MPI rank and a two-orbital matrix fixture, not a material calculation.
- The new propagation regression target is excluded when ENABLE_CUBLASMP is enabled: it exercises CPU and CPU Tensor-LAPACK paths and does not supply GPU kernel or distributed-library link dependencies. Existing test targets are unchanged.

## Linked Issue
No separate issue: this PR contains the localized defect and regression evidence. Related PRs were searched before preparation.

## Governance Notes
No INPUT syntax, defaults, or parser metadata changes. No change to the propagation formula. This is a bug fix in existing RT paths. No governance exception requested.
Core impact: module_rt normalization call sites only; no interfaces changed.
