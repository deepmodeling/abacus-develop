# Onsite LCAO SCF observer (local diagnostic build)

Set `ABACUS_SCF_DIAGNOSTICS=1` in the runtime environment to enable the
observer. It writes one double-precision JSON record per completed electronic
iteration to `OUT.<suffix>/scf_diagnostics.jsonl`. No INPUT defaults, SCF
updates, thresholds, or mixing coefficients are changed. All MPI ranks
participate in the extra projection reductions; only rank zero writes the
file. Each record is flushed immediately.

Supported scope: `dft_plus_u=1`, `nspin=2`, `symmetry=0`, `mixing_dftu=0`,
`mixing_dmr=0`. Unsupported layouts raise an error rather than being
reinterpreted as real collinear matrices. The analysis accepts any number of
atoms and iterations, and can preserve complex Hermitian matrices if supplied;
the source observer itself does not implement spinor or occupation-mixing
output. This is a local diagnostic extension, not an upstream feature.

## Matrix definitions and timing

`n_in` is a copy of the exact per-atom, per-spin buffer passed to
`DFTU_LCAO::cal_pot_onsite` in `DFTU_onsite::contribute_hubbard_atoms`.
When the first iteration has no DMR and constructs no Hubbard potential,
`n_in` is `null`. It must not be replaced with a zero matrix.

`n_raw_out` projects the newly solved DMR using the same neighbor list,
two-center overlap table, contraction order, and MPI reduction as the original
operator. The pure projection helper does not write into `Plus_U_Base`, DMR,
the potential, or mixing buffers. It is observed after the base solver computes
DRHO, magnetism, gap, and energy, and before optional DMR mixing. Charge-density
mixing changes rho, not this solved DMR.

`n_after_occupation_mixing` is `null`: no occupation mixer runs in this supported
onsite path. It is not a synonym for raw output or charge-mixed density.

The representation is the real onsite projector basis, not the nonorthogonal
NAO coefficient matrix or a Mulliken/Lowdin population matrix. Projectors use
the first radial NAO channel, smoothed and normalized inside `onsite_radius`;
same-site angular channels form an orthogonal shell. There is no global
orthogonalization between different atomic sites. The contraction is
`n_mm' = <chi_m|rho_sigma|chi_m'>`, evaluated with the original NAO real-space
density matrix and `<phi|chi>` overlaps. Flattened matrices are row-major;
atoms are one-based and spin is zero-based. The signed real d-harmonic order is
`dz2, -dxz, -dyz, dx2-y2, dxy` (`m=0,+1,-1,+2,-2`). See `ylm.cpp`,
`dftu_nao_adj.cpp`, `atomic_radials.cpp`, and `smoothgen.cpp` for the angular
convention and radial normalization.

## Energies and gap

The source energy actually added as `E_plusU` is separately preserved as
`source_band_correction_eV`. For the supported collinear onsite path,
`cal_pot_onsite` accumulates **U/2 Tr(n_in^2)**. This is a band-energy
replacement correction, not the Dudarev functional. The latter is
`E_Dud = U/2 Tr(n - n^2)`. With `V_U = U(1/2 - n)`, the source expression
equals `E_Dud - Tr(V_U n)`. The analysis computes fractionality and the
Dudarev functional for both observed stages, and independently checks the
source correction against the input matrix. It does not modify the source
energy expression or total energy. In iteration one the source correction is
genuinely zero while the absent input remains `null`.

`logged_E_gap_k_eV` is the unchanged `ElecState::bandgap`. Despite the log
label `E_gap(k)`, `cal_bandgap` searches all sampled k points and bands around
`ef + 1e-5 Ry`; this is not an array of gaps at individual k points or a minimum
direct gap. A metallic sampled spectrum can have a positive level spacing by
this definition. `fundamental_gap_eV` is explicitly `null`: this observer does
not independently establish the fundamental gap or metallic character.

`total_energy_eV` is the original electronic free energy `f_en.etot` per cell.
Cell magnetism is calculated from rho before charge mixing. Local d-shell
moments derived from raw occupation traces are partial projected moments,
not complete atomic Mulliken moments.

## Validation

Focused writer tests cover precision, absent input, append behavior, writer
rank, and preservation of observational copies. A same-input VO2/U6 comparison
with diagnostics disabled checks the numerical trajectory at stdout/log
precision. Derived-data tests check complex Hermitian norms and missing first
differences. No missing matrix, step, gap, or residual is interpolated.

## Charge and spin density channels (schema 2)

Immediately before `ESolver_KS::iter_finish` calls charge mixing, the LCAO
observer reads `chr.rho` (raw output) and `chr.rho_save` (input saved before
sum_band). The stateless `module_charge::cal_spin_channel_residuals` computes

```
Rc = omega / (nxyz * nelec) * sum_grid |delta_rho_up + delta_rho_down|
Rm = omega / (nxyz * nelec) * sum_grid |delta_rho_up - delta_rho_down|
```

Each sum uses `Parallel_Reduce::reduce_pool`, exactly as the native real-space
DRHO (`scf_thr_type=2`). These dimensionless per-electron L1 residuals use the
same local grid, cell volume and electron count. No mixed density is read and
no density or mixer buffer is changed. They are channel norms, not an additive
split of DRHO: `max(Rc,Rm) <= DRHO <= Rc+Rm`. Pointwise the spin-resolved
DRHO integrand is `max(abs(delta_rho_c), abs(delta_rho_m))`. The observer rejects
other spin layouts or convergence metrics rather than silently reinterpreting
them. There is no need to substitute an L2 metric for this real-space L1 run.

The actual `mixing_restart_count` is copied after iter_finish. The
`mixing_restart_executed` flag requires a positive restart threshold, the
current iteration equal to the scheduled step, and an already incremented
counter. Analysis additionally verifies the flag against counter increments.
For the fixed-U, no-DMR-mixing cases in this experiment the scheduled step
cannot move during execution. The unchanged source does the following:

1. After a step with DRHO <= mixing_restart, schedule iteration k+1.
2. Skip density mixing on step k, retaining raw output (then normalize rho).
3. In the next iter_init, call init_mixing and increment restart_count.
4. Do not declare convergence on that restart step. Without U ramping the
   scheduled step is not rearmed, so this is a single threshold-triggered event.

Thus this experiment tests the program's combined raw-density handoff and
history reset, not an isolated algebraic deletion of historical vectors.
Source locations: `chg_routine.cpp::chgmixing_ks` and
`chgmixing_ks_lcao`; `chg_mix.cpp::init_mixing/reset_scf_history`.

`scf_wall_seconds_before_occupation_observer` measures elapsed time since the
native iter_init timestamp, read just after base iter_finish and before the
extra occupation projection/writer. It includes the density observer and
native per-step output; it excludes the extra occupation observer. Summing it
is SCF iteration wall time, not process launch/initialization/finalization or
complete elapsed job wall time. `scf_converged` is the actual solver flag.
Local projected d moments in analysis are Tr(n_up_raw)-Tr(n_down_raw) in
muB per Hubbard atom; they are partial shell moments, not total atomic moments.

## Shared PW observer (schema 3)

The record model and writer now live in `source_estate/scf_diagnostics.*`,
under `ScfDiagnostics`, and are built unconditionally. LCAO still observes its
actual onsite input buffer; its established binary remains pinned for the
LCAO rerun campaign.

For CPU PW, `DFTU_BASE::observe_occ_pw` projects the freshly solved wavefunctions
into a separate `OccupationMatrix`, using exactly the native onsite projector,
k-point occupations, scalar contraction and inter-pool reduction. It never calls
the occupation mixer or writes a potential, energy or production occupation
matrix. It restores the wavefunction k cursor and projector cache to that k.
The native `dftu.occmat()` retained after iter_init is the actual Hamiltonian
input; without a loaded seed the first input is absent. The raw output is not
the one-step-lagged native occ_mat.txt. Supported scope is collinear nspin=2,
double precision and mixing_dftu=0; no independent occupation mixer is claimed.

PW projectors are the original smoothed onsite atomic projectors represented
in plane waves, with overlaps contracted over bands and weighted by the native
k-point occupations. They are not LCAO Mulliken populations. Each matrix is a
real scalar-spin matrix as in the native PW DFT+U implementation. Atom, spin,
orbital channel, actual U and full matrices are retained without averaging.
The same U/2 Tr(n_input^2) source correction convention applies in this scalar
PW path; Dudarev functionals remain derived quantities with a distinct name.

The original PW `scf_thr_type=1` is retained. Its native residual combines a
Coulomb-weighted reciprocal charge term (G!=0) and an unweighted reciprocal
magnetic term (including G=0), rather than the LCAO per-electron L1 norm. To
compare charge/spin using identical norms, the PW observer additionally writes

```
Rc_L2 = omega/nelec * sqrt(sum_grid(delta_rho_c^2)/nxyz)
Rm_L2 = omega/nelec * sqrt(sum_grid(delta_rho_m^2)/nxyz)
```

Both are dimensionless L2 norms relative to the mean density nelec/omega, use
raw output minus saved input before mixing, and use the same pool reduction.
L1 channels are also retained as explicitly named additional measurements.
`scf_thr_type` identifies the native convergence metric. No equality between
these L2/L1 channels and the native PW DRHO is implied. L2 fields are null
where not observed; missing values are not replaced with zero.
