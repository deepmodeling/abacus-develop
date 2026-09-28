# Isolated Environ chain-CG experiment

Branch sccs-environ-chain-opt is an experimental path, restored from validation/qe_environ_water_20260927/environ_noadjoint_full_scf. The stable sccs-pcc-lkfeiyi branch retains the discrete-adjoint implementation.

This branch uses chain/FFT cavity gradients, the Environ-style sqrt-preconditioned CG equation, Gaussian ionic sources of width 0.5 bohr, and continuum reaction/cavity electronic derivatives. Ionic forces use the Coulomb potential of the continuous polarization charge, matching Environ dielectric_of_potential -> force_fft. This force-only replacement preserves the electronic potential and energy. It does not solve a discrete adjoint. PCC0D and PCC2D are now connected through the existing polarization-charge fixed-point solver with the chain cavity gradient and the analytic PCC field. The periodic case retains sqrt-CG; no boundary silently switches to an adjoint. Independent vacuum PCC remains unchanged and is added once. The ionic-shape derivative uses the same Gaussian ionic source as the reaction energy.

The PCC solver respects sccs_mixing_type, damping/history, warm-start polarization and the usual RMS/max residual criteria. The periodic CG interpretation below remains specific to assume_isolated none. This explicit split avoids losing the nonzero screened charge through a periodic FFT Laplacian. It does not claim an accelerated CG implementation for PCC.

For reproduction of the historical experiment, sccs_tol_rms is interpreted as the sum of squared residuals (1e-18 in the retained water example), not the RMS criterion of the stable branch. The copied mixing parameters and initial polarization charge do not control this cold-start CG solver. No new supported INPUT mode is being introduced; parameter documentation is not regenerated for this isolated experiment branch. Keep its executable and experiment inputs together.

The baseline handles an initially zero source without a CG null-residual breakdown. Otherwise its numerical path follows the archived experiment. The retained experimental force checks previously had a nonzero finite-grid analytic-versus-FD discrepancy; speed optimization does not establish or improve variational force consistency.

Baseline and optimization validation, commands, runtime logs and source provenance are recorded under /home/lyt/DFT/abacus_sccs/validation/chain_optimization_20260927. Build directories are independent of the stable executable.
