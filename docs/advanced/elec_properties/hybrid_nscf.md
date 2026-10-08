# Screened hybrid NSCF with plane waves

A hybrid Hamiltonian needs occupied orbitals as well as the charge density.
PW screened hybrid NSCF now reads a **frozen SCF ensemble** on a source q mesh
and diagonalizes the Hamiltonian on an independent target k list. Increasing
`nbands` in NSCF changes the target states, not the source ensemble.

This first implementation supports CPU and CUDA GPU calculations, `kpar 1`, `bndpar 1`,
`nspin 1` or `2`, and a complete, uniformly weighted source mesh without symmetry
reduction. Set `symmetry -1`, `exxace false` and
`exx_gamma_extra false` in both calculations. HSE is supported;
unscreened Fock exchange, including PBE0, is rejected pending validation of the
singularity correction for arbitrary target k points. Forces and stress are
also rejected in hybrid NSCF. The LCAO workflow is unchanged.

## Prepare the source ensemble

Run a converged hybrid SCF with the normal source `KPT` and these settings:

```text
calculation scf
basis_type pw
dft_functional HSE
symmetry -1
kpar 1
bndpar 1
exxace false
exx_gamma_extra false
out_chg 1
out_wfc_pw 2
out_freq_ion 0
cal_force 0
cal_stress 0
```

The output directory contains the charge density, existing `wfk*_pw.dat`
binary wavefunctions and `eig_occ.txt`. No new checkpoint format or companion
file is required. Precise source q coordinates, k weights and band count come
from the existing wavefunction headers; the occupation column in `eig_occ.txt`
already contains the k-weighted occupations required by EXX, including fractional
occupations. Use files from the same converged SCF calculation. The workflow
consumes these outputs without changing the existing wavefunction or occupation
writers/readers.

## Solve target states

Keep the same structure, pseudopotentials, functional, cutoffs and FFT dimensions.
Point `read_file_dir` to the SCF output, replace `KPT` with the target mesh or
band path, and use:

```text
calculation nscf
read_file_dir ../scf/OUT.hybrid
init_wfc random
nbands 20
pw_diag_thr 1e-10
out_band 1
out_wfc_pw 0
out_chg 0
```

Retain the source SCF settings listed above, except its output and calculation
settings. Use a separate NSCF output directory to preserve the SCF files.
`init_wfc file` is rejected for this initial implementation: binary source
orbitals are read separately, whereas target orbitals start independently.
`nbands` may exceed the SCF band count. A different number of MPI processes
within the single PW pool is supported on CPU; binary coefficients are redistributed
by their Miller indices. CUDA NSCF currently uses one MPI rank. Select
`device gpu` and `precision double` for the validated CUDA workflow; the source
SCF may run on CPU or GPU. The binary checkpoint format is shared between
these devices. ROCm NSCF remains disabled pending backend verification.
Use `OMP_NUM_THREADS=1` for runtime tests.

Missing or truncated wavefunctions/occupation files and invalid occupations
produce an error. Source dimensions, spin order, uniform q weights and occupation
K points must agree. The existing wavefunction reader additionally verifies k
coordinates, cell, source band count and plane-wave count. `INPUT.info`, when
available, is used only to warn about exchange configuration differences. Older
outputs without this metadata are accepted with a warning; regeneration is not
required solely because metadata is missing.
NSCF does not update
the frozen ensemble or run the EXX SCF outer loop. Its target occupations are
used for ordinary output only, and its printed total energy is not a new
self-consistent hybrid total energy.

## Algorithm and implementation boundaries

For each target state, exchange uses

\[
V_x\psi_{n\mathbf{k}} = -\alpha\sum_{\mathbf{q},m}
 w_{\mathbf{q}} f_{m\mathbf{q}}\,
 \psi_{m\mathbf{q}}\,
 \mathcal{F}^{-1}\!\left[
 v(\mathbf{k}-\mathbf{q}+\mathbf{G})
 \mathcal{F}[\psi^*_{m\mathbf{q}}\psi_{n\mathbf{k}}]
 \right].
\]

Here the periodic parts of source and target orbitals have separate PW bases.
The uniform source sum uses the same spin and weighted-occupation convention as
ABACUS SCF. The target `KPT` weights never define the exchange quadrature.
`OperatorEXXPW` accepts explicit source basis, points, orbitals and occupations;
the solver owns their lifetime. Coulomb-kernel construction uses target k
coordinates and source q coordinates, with the screened zero-transfer correction
computed from the source mesh. The two-basis path currently uses the full FFT
grid, including for distributed PW transforms. Small-grid acceleration remains
separate follow-up work. CUDA executes the pair-density FFTs and exchange
application on the GPU with separate source and target PW maps. This frozen-source NSCF workflow deliberately uses
direct exchange; ACE is not part of its implementation roadmap. The SCF ACE
projectors are tied to their original target subspace and cannot be reused as a
validated exchange operator on an independent band path.

The local `../q-e` implementation was used to analyze periodic pair densities,
`k-q+G` convolution, occupation normalization and the finite screened
zero-transfer term (`PW/src/exx.f90`, `exx_base.f90`). Its Fortran module/global
architecture was not adopted. That checkout explicitly rejects hybrid NSCF in
`PW/src/setup.f90`, so it cannot provide a direct NSCF reference. Also, its
built-in HSE screening default is `0.106 bohr^-1`; use
`screening_parameter=0.11` for comparison with ABACUS HSE06.

## Screened zero-transfer treatment

PW HSE uses the same `exx_singularity_correction limits` default as NAO/LCAO.
For screened exchange the reciprocal kernel has the finite limit
`-pi * e2 / omega^2` at `k-q+G=0`; this scheme retains that limit without adding
an auxiliary-function correction. It uses `exx_gamma_extra false`,
which is selected automatically when omitted. Explicit gamma extrapolation
with `limits`, and unscreened Fock terms with `limits`, are rejected in PW.

Set `exx_singularity_correction gygi` explicitly to retain the historical
PW Gygi-Baldereschi auxiliary-function correction (including the existing
gamma-extrapolation option).
This is also the default for unscreened PW hybrids; it does not enable
unscreened hybrid NSCF. For screened independent-k targets the Gygi correction
can cause a jump when a target coincides with a source point.
The name follows [QE's `gygi-baldereschi` terminology](https://www.quantum-espresso.org/Doc/INPUT_PW.html).
LCAO-specific `spencer`, `revised_spencer`, `massidda` and `carrier` schemes
are not implemented for PW and are rejected rather than silently ignored.

Use the same scheme in SCF and NSCF for consistent bands. A different scheme
produces a warning and continues with the existing frozen source and the current
NSCF exchange settings. This also permits deliberate comparisons between schemes.

The cross-code results below used the historical
Gygi scheme and do not constitute validation of the new default.

## Reproducible verification

```bash
cmake --build build --target MODULE_ESOLVER_pw_hybrid_source -j 8
OMP_NUM_THREADS=1 ctest --test-dir build -V -R '^MODULE_ESOLVER_pw_hybrid_source$'
OMP_NUM_THREADS=1 python3 tests/integrate/tools/test_hybrid_nscf.py ./build/abacus --mpi-ranks 2
```

Supply the actual configured executable and build directory. The integration
script uses the existing H pseudopotential and structure fixture, performs SCF
and NSCF, and checks extra target bands, a different target k list, a shared
Gamma point, MPI redistribution, both collinear spin channels, immutable source
files, and invalid restart/unsupported-option errors. It needs LibXC and MPI.

For a CUDA build and an available GPU, also run:

```bash
OMP_NUM_THREADS=1 python3 tests/integrate/tools/test_hybrid_nscf.py ./build/abacus --gpu
```

This compares CPU and GPU targets using one frozen CPU source on both a mesh
and an independent path with extra target bands, then checks a GPU-generated
spin-polarized source with CPU and GPU targets. Source files must remain
unchanged. A three-band single-precision GPU path additionally checks conversion
from a double-precision CPU source and the full-grid exchange FFT precision.
Single precision uses the existing looser diagonalization threshold; use
`precision double` for accurate high-energy target bands. The workflow script is run manually; only unit tests are registered with CMake/CTest.

On a local RTX 3090 with CUDA 13.1, the CPU and GPU workflow scripts passed in a CUDA build.
The H CPU/GPU mesh, independent path and both spin-channel eigenvalues agreed
at the printed precision, as did all eight bands at nine Si L-Gamma-X points
using the same frozen CPU source (`20 Ry`, `2x2x2` source mesh). The Si check
compared raw eigenvalues without an energy shift. The three-band H
single-precision GPU path differed by at most `5.76e-5 eV` from double precision.
The CPU build also passed its source-adapter unit tests and NSCF integration test.
This validates device
consistency; it does not change the cross-code convergence limits below.

On the initial local two-q-point, 10 Ry test, the same-mesh SCF/NSCF maximum
band difference was `3.7e-5 eV`; the shared Gamma point and serial/two-rank NSCF
results agreed to printed precision. The source-adapter unit tests cover occupations, spin meshes and binary headers.
A QE SCF comparison at this low cutoff showed band differences up to about
`0.03 eV`; that H test alone does not establish cross-code agreement.

A separate Si comparison used the same two-atom primitive cell (`a=10.2 bohr`),
HSE exchange fraction `0.25`, screening `0.11 bohr^-1`, `20 Ry` cutoff,
`2x2x2` source mesh, eight bands and nine L-Gamma-X target points. VASP 6.2.1
used `HFSCREEN=0.207869873 Angstrom^-1`, `PRECFOCK=Accurate`, `ISYM=-1`,
`ALGO=Damped`, `EDIFF=1e-9` and `NELMIN=40`. Its explicit k list contained the
eight source points and six additional zero-weight points. Target points
already present in the source mesh were included only once. This follows the
[VASP hybrid band workflow](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals).

After aligning each code's Gamma valence-band top, the comparison gave:

| Quantity (eV) | ABACUS NSCF | VASP |
| --- | ---: | ---: |
| Gamma direct gap | 3.604645 | 3.619060 |
| Gap on the sampled path | 1.449187 | 1.486063 |
| L direct gap | 3.967258 | 4.005180 |
| X direct gap | 4.655743 | 4.681974 |

For bands 2-6 across all nine points, the RMS difference was `0.042351 eV`
and the maximum difference was `0.086995 eV`. Both codes placed the lowest
sampled conduction state at fractional `(0.375, 0, 0.375)`.
ABACUS used `Si_ONCV_PBE-1.0.upf`; VASP used `PAW_PBE Si 05Jan2001`.
The coarse mesh, sparse path and different pseudopotentials make this a
qualitative cross-check, not a cutoff or k-grid convergence study.

An initial VASP reference with duplicate zero-weight Gamma points gave an
inconsistent band-4 eigenvalue: it differed from the weighted source Gamma
by `0.686083 eV`. Raising `NELMIN` to 60 did not remove that discrepancy.
The comparison above is from a fresh calculation with unique k coordinates,
using the weighted mesh result at coincident path points. No eigenvalues were
edited to align the spectra.

These checks do not validate PBE0, symmetry reconstruction or source q pools. Direct exchange remains the NSCF method; ACE is not a
planned extension of this workflow.
