# Isolated Environ chain-CG experiment

Branch sccs-environ-chain-opt is an experimental path, restored from validation/qe_environ_water_20260927/environ_noadjoint_full_scf. The stable sccs-pcc-lkfeiyi branch retains the discrete-adjoint implementation.

This branch uses chain/FFT cavity gradients, the Environ-style sqrt-preconditioned CG equation, Gaussian ionic sources of width 0.5 bohr, and the corresponding continuum reaction/cavity derivatives. It does not solve a discrete adjoint. PCC is rejected for SCCS evaluations; the isolated charged-system PCC extension is not validated here. Existing independent PCC-only behavior is not changed.

For reproduction of the historical experiment, sccs_tol_rms is interpreted as the sum of squared residuals (1e-18 in the retained water example), not the RMS criterion of the stable branch. The copied mixing parameters and initial polarization charge do not control this cold-start CG solver. No new supported INPUT mode is being introduced; parameter documentation is not regenerated for this isolated experiment branch. Keep its executable and experiment inputs together.

The baseline handles an initially zero source without a CG null-residual breakdown. Otherwise its numerical path follows the archived experiment. The retained experimental force checks previously had a nonzero finite-grid analytic-versus-FD discrepancy; speed optimization does not establish or improve variational force consistency.

Baseline and optimization validation, commands, runtime logs and source provenance are recorded under /home/lyt/DFT/abacus_sccs/validation/chain_optimization_20260927. Build directories are independent of the stable executable.
