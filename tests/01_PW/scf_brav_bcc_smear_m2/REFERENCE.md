# Numerical reference calibration (2026-10-29)

The original SCF precision left solver-dependent differences in the energy or
force/stress checks. On source `61d4b7ba31696259d5dbee0c9f94604631d6b646`,
the original CG input passed while replacing only `ks_solver` with PPCG produced
a reference-comparison warning. The original input, reference and failure logs
are retained in the external validation record.

This fixture now uses `scf_thr = 1e-12`,
`pw_diag_thr_max = 1e-12` and
`pw_diag_nmax = 200`. Structure, pseudopotentials, cutoff, occupations, mixing,
seed, original solver selection and acceptance tolerances are unchanged.

The updated scalar reference was generated independently with Davidson-subspace,
4 MPI ranks, `scf_thr = 1e-14`, `pw_diag_thr_max = 1e-12` and
`pw_diag_nmax = 200`. CG and PPCG runs at the same tighter precision provide
cross-checks. Tighter diagnostic PPCG runs that reached the floating-point
residual floor are retained as diagnostic records, not final acceptance runs.

After reference generation, the final input was rerun from clean case directories
with CG, Davidson-subspace and PPCG at both 1 and 4 MPI ranks. All 6 runs passed
the unchanged repository Autotest checks (24 property checks),
reported SCF convergence and produced no PPCG inner residual warning. Relative
to the tighter Davidson calculation, the largest differences were:

- Total energy: 6.9348971e-12 eV.
- Force component: 0 eV/Angstrom.
- Stress component: 0.0001619104 kbar.

Default acceptance tolerances remain energy `1e-7` eV, force aggregate `1e-4`
and stress aggregate `1e-3`; no threshold file or comparison rule was changed.
The historical timing field in `result.ref` is retained and is not an acceptance
metric. CPU executable SHA256:
`ba7a329ce12350558eaf927c2ef3ccae8d7acc614c572c26ec01add36a09fbdb`.

Local reproducibility record:
`/mnt/data/Cheerly/ppcg-warning-followup-20261029/`.
`probe.py` records the exact commands and runtime environment; `baseline/`
preserves the old failures, `accurate/`, `tight/`, `validated/` and `balanced/`
preserve precision studies, and `final/` holds the clean acceptance reruns.
`candidate-inputs.json`, `candidate-audit.json` and `final-np{1,4}-results.json`
summarize the inputs, independent reference source and outcomes.

Supplementary GPU controls used the matching Libxc 5.1.7 CUDA build from
source `4cec8adad6291ddd43317636ef0b4ff21e94cce9`: CG and PPCG each passed
the same checks at one MPI rank, with SCF convergence and no inner PPCG
warning. These runs validate fixture compatibility with that GPU build;
current-source GPU validation is separate. Logs are in the local validation
record under `gpu-controls/`.
