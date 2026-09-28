# Isolated Environ chain-CG experiment

Branch sccs-environ-chain-opt is an experimental path, restored from validation/qe_environ_water_20260927/environ_noadjoint_full_scf. The stable sccs-pcc-lkfeiyi branch retains the discrete-adjoint implementation.

This branch uses chain/FFT cavity gradients, the Environ-style sqrt-preconditioned CG equation, Gaussian ionic sources of width 0.5 bohr, and continuum reaction/cavity electronic derivatives. For PCC, ionic forces use the Coulomb potential of the continuous polarization charge, matching Environ dielectric_of_potential -> force_fft. Periodic sqrt-CG ionic forces use the solved reaction potential, consistent with its source-energy derivative. These force choices preserve the electronic potential and energy. It does not solve a discrete adjoint. PCC0D and PCC2D are now connected through the existing polarization-charge fixed-point solver with the chain cavity gradient and the analytic PCC field. The periodic case retains sqrt-CG; no boundary silently switches to an adjoint. Independent vacuum PCC remains unchanged and is added once. No Gaussian ionic-shape energy or force is added (see the PCC2D energy convention below).

The PCC solver respects sccs_mixing_type, damping/history, warm-start polarization and the usual RMS/max residual criteria. The periodic CG interpretation below remains specific to assume_isolated none. This explicit split avoids losing the nonzero screened charge through a periodic FFT Laplacian. It does not claim an accelerated CG implementation for PCC.

sccs_tol_rms and sccs_tol_max are the RMS and maximum charge residual (e/bohr^3) on every boundary; the periodic sqrt-CG stops only when both pass. The historical experiment interpreted sccs_tol_rms as ENVIRON's unnormalized sum of squared residuals (1e-18 in the retained water example); the corresponding RMS is sqrt(1e-18/N) for N FFT grid points. With sccs_debug 2 the solver verifies the preconditioned fixed point v = P(q - K v) and reports SCCS_RESIDUAL and SCCS_CG_FIXED_POINT_DEFECT; the kernel itself no longer writes to standard output. The PCC mixing parameters and initial polarization charge do not control this cold-start CG solver. Parameter documentation is regenerated from the binary.

The baseline handles an initially zero source without a CG null-residual breakdown. Otherwise its numerical path follows the archived experiment. The retained experimental force checks previously had a nonzero finite-grid analytic-versus-FD discrepancy; speed optimization does not establish or improve variational force consistency.

Baseline and optimization validation, commands, runtime logs and source provenance are recorded under /home/lyt/DFT/abacus_sccs/validation/chain_optimization_20260927. Build directories are independent of the stable executable.

The unused discrete-adjoint solver, its state and diagnostic fields, and its
finite-grid derivative tests have been removed from this branch. Coulomb-gradient
transpose operations remain mathematical operator utilities; their inner-product
and gradient-only tests now live in `test_sccs_pw_coulomb.cpp`. These utilities do
not introduce an adjoint solve into the chain path. The removed solver and its
verification remain available in Git history and the stable branch.

Cleanup verification rebuilt ABACUS and five focused test targets. The Coulomb
operator, PCC2D operator, driver, and potential/debug suites passed. The force
suite at that cleanup commit retained two failing finite-difference checks: the periodic reaction-force
error is 0.0070829526759735306 Ha/bohr, and the largest PCC2D component error is
0.00014111774946192097 Ha/bohr. Rebuilding the force suite from the pre-cleanup
commit reproduced all seven reported error values exactly. No force tolerances
or references were changed. The detailed commands and logs are recorded in
`/home/lyt/DFT/abacus_sccs/validation/chain_cleanup_20260928/REPORT.md`.

The subsequent periodic-force correction removes the use of reconstructed
continuous polarization in the periodic ionic force. The existing 1e-4 Ha/bohr
force tolerance is unchanged; verification now covers x/y/z at two displacement
sizes. It does not change the continuum electronic-cavity derivative into an
exact derivative of the discretized sqrt-CG coefficients.

The PCC2D chain approximation is explicitly tested with a relaxed accuracy goal
of 0.01 eV/Angstrom, converted to Ha/bohr using ABACUS constants. Only this
synthetic epsilon=1.1 fixed-density SCCS+PCC2D test has changed tolerance; exact
PCC-only and source-derivative tests retain their thresholds. Every measured
component error is printed, so acceptance does not conceal the discrepancy.
This user-requested acceptance criterion is not evidence for the same accuracy
in water (epsilon=78.3), arbitrary geometries, or converged SCF total forces.

## PCC2D energy convention for charged slabs

The monopole constant of the planar parabolic correction follows ENVIRON
`core_1da.f90` (`-pi*Q/(3*L_y)`, used by all of its 1D-analytic paths).
Andreussi and Marzari, Phys. Rev. B 90, 245101 (2014), Eq. (88) prints the
opposite sign. For a square periodic plane (A = L_y^2) ENVIRON's value equals the
unshifted open planar kernel `-2*pi*|u|/A` minus the zero-mean periodic planar
kernel; `SquareCellCorrectionUsesUnshiftedOpenPlanarKernel` checks this relation.

The Appendix A ionic-shape energy (Eq. 107) and its forces were removed. The
vacuum PCC already uses point ions, and the reaction energy is obtained with
Gaussian ions through the complete open-boundary operator. Inside the
epsilon=1 region the reaction potential is harmonic, so by the mean-value
property a spherical Gaussian ion has the same reaction energy as a point ion.
Adding Eq. (107) made the total energy depend on the fictitious Gaussian width.
ENVIRON adds no such term when a PBC correction is active (`degauss = 0` in
`solver_setup.f90`), and the PCC0D path never added the analogous Eq. (106).

Neutral systems are unchanged. Charged PCC2D energies differ from the previous
convention by `-(pi/(3*L_y))*Q^2/eps` (solvated; `eps = 1` for vacuum+PCC)
minus the former Eq. (107) term.
