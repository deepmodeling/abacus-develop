# Isolated Environ chain-CG experiment

Branch sccs-environ-chain-opt is an experimental path, restored from validation/qe_environ_water_20260927/environ_noadjoint_full_scf. The stable sccs-pcc-lkfeiyi branch retains the discrete-adjoint implementation.

This branch uses chain/FFT cavity gradients, the Environ-style sqrt-preconditioned CG equation on every boundary, Gaussian ionic sources of width 0.5 bohr, and continuum reaction/cavity electronic derivatives. It does not solve a discrete adjoint. For PCC0D and PCC2D the preconditioner's Poisson solve includes the analytic open-boundary correction, and the solution keeps that physical gauge instead of the periodic zero mean (ENVIRON generalized_sqrt with has_corrections). Ionic forces on every boundary use the solved reaction potential, the fixed-density derivative of the sqrt-CG reaction energy; these force choices preserve the electronic potential and energy. Independent vacuum PCC remains unchanged and is added once. No Gaussian ionic-shape energy or force is added (see the PCC2D energy convention below).

PCC formerly used a polarization-charge fixed point. Its charged monopole mode contracts only by about 1 - 1/epsilon per step, so a residual-based stop left a polarization-charge error of roughly residual*epsilon: for H3O+ (water-cation, PCC0D, 50 Ry) Q_pol wandered between -0.978 and -0.988 e (Gauss value -0.9872) across SCF steps and the SCF stalled at drho 1e-5 to 1e-4, and charged PCC2D failed the Gauss check. The sqrt-CG has no such slow mode. With PCC it deviates from the ENVIRON defaults in two places:

- factsqrt and grad ln(eps) come from FFT derivatives of the switching function (ENVIRON deriv_method 'fft'; the periodic path keeps the default 'chain'). The potential at the cavity edge carries the open-boundary monopole or dipole, and the chain-rule factsqrt, which drops steeply to zero at density_min, leaves a sampled f v error proportional to it that does not converge with the grid.
- The fatal PCC2D Gauss check uses the far-field polarization charge of the solution, int(s)/sqrt(eps_bulk) - int(q), which converges with the grid. ENVIRON's dielectric_of_potential density integral, (1/4pi) grad ln eps . grad v + (1/eps - 1) q, remains the reported polarization density and its moments; its integral can miss Gauss's law by 2e-3 during the SCF. sccs_debug 2 prints both as SCCS_GAUSS.

By default the electronic cavity potential is ENVIRON's continuum -eps'|grad v|^2/(8 pi) with the FFT gradient of the solved potential, so energies and energy finite differences reproduce QE-Environ with deriv_method 'fft' (H3O+ at 100/400 Ry: energy within 0.2 meV, FD slope within 1e-4 .. 6e-3 eV/Angstrom). That FFT gradient rings at the open-boundary step or kink of the corrected potential; the ringing reaches the cavity when it approaches the open boundary (0.023 against a peak field of 0.032 in the layered unit test, whose cavity ends 1.2 bohr from it). An analytic reconstruction grad v = (grad w - w grad(ln eps)/2)/sqrt(eps), w = C_PCC(s), avoided the ringing but moved the H3O+ energy by 2 meV and the FD slope by 0.045 eV/Angstrom away from ENVIRON at eps 78.3 for unexplained reasons, and was removed (validation/slope_localize_20260928/REPORT.md).

The continuum potential is not the derivative of the discrete sqrt-CG energy. The analytic solvation force then misses energy finite differences by 1-3e-2 eV/Angstrom for H3O+ (QE-Environ's own force misses by 4-6e-2 because its ionic force uses the dielectric_of_potential polarization density). The exact derivative needs no adjoint solve because A = sqrt(eps) G^-1 sqrt(eps) + F is symmetric: with L = ln(eps_bulk) and b = eps v^2/(8 pi), dE/ds = L/2 (q v + lapl b + L div(b grad s)). Without filtering, that derivative of the FFT factsqrt is not grid-converged (a breathing-mode cavity derivative of 20.3, 0.56 and 0.049 at 120, 240 and 400 Ry against the continuum 0.29), reaches 150-430 Ha pointwise at the cavity edge, and makes the H3O+ SCF diverge (validation/exact_derivative_20260928/REPORT.md).

sccs_lowpass_p1 and sccs_lowpass_p2 (ENVIRON 3.1.1 deriv_lowpass_p1/p2 of core_fft_lowpass) multiply every switching-function derivative by 0.5 erfc(p1 G^2/Gcut^2 - p2). With both positive (PCC only) the electronic cavity potential is the exact derivative above with the filtered operators. At 10/5 the derivative converges between 300 and 500 Ry, stays within 10% of the continuum potential pointwise, and the H3O+ solvation force matches finite differences to 1e-4 .. 3e-3 eV/Angstrom (PCC0D and PCC2D) at unchanged cost per SCCS evaluation. Energies match QE-Environ with the same lowpass to 0.7 meV and FD slopes to 1-2e-3; the filter itself raises the H3O+ solvation energy by about 12 meV. The periodic path keeps the chain factsqrt and the continuum potential and rejects the lowpass.

The sccs_mixing* INPUT parameters are still accepted but no longer control any solver; the fixed-point routine remains only with its unit tests. At 100/400 Ry, H3O+ PCC0D and PCC2D converge in 16-25 SCF steps (the fixed-point binary did not converge in 100, or stopped at the Gauss check). Commands, logs and analysis are in /home/lyt/DFT/abacus_sccs/validation/pcc_sqrt_cg_20260928/REPORT.md.

sccs_tol_rms and sccs_tol_max are the RMS and maximum charge residual (e/bohr^3) on every boundary; the sqrt-CG stops only when both pass. The historical experiment interpreted sccs_tol_rms as ENVIRON's unnormalized sum of squared residuals (1e-18 in the retained water example); the corresponding RMS is sqrt(1e-18/N) for N FFT grid points. With sccs_debug 2 the solver verifies the preconditioned fixed point v = P(q - K v) and reports SCCS_RESIDUAL and SCCS_CG_FIXED_POINT_DEFECT; the kernel itself no longer writes to standard output. Like ENVIRON generalized_sqrt, it warm-starts from the previous unshifted CG solution (ENVIRON reuses its returned potential, which is zero-mean only without PCC) with one preconditioned fixed-point step, v = P(q - K v_old), whose residual is K (v_old - v); the guess is kept only when its residual RMS is below the cold-start residual (ENVIRON uses a fixed 1e-2 sum-of-squares threshold). Parameter documentation is regenerated from the binary.

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

The fixed-density SCCS+PCC2D force test (synthetic epsilon=1.1, neutral and
charged) again uses the original 1e-7 Ha/bohr threshold. With the reaction
potential of the PCC sqrt-CG, the largest measured component error is 4.6e-9
Ha/bohr; the continuous-polarization force of the former fixed-point path
reached 1.4e-4 Ha/bohr. The test grid is 120 Ry so that the far-field Gauss
check of its 0.5 bohr Gaussian ion passes (7e-4 at the former 30 Ry). Every
component error is printed. This synthetic test is not evidence for the same
accuracy in water (epsilon=78.3), arbitrary geometries, or converged SCF total
forces.

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
