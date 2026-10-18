# PW screened hybrid NSCF

Hybrid NSCF reads the charge density, binary `wfk*_pw.dat` wavefunctions and
weighted occupations in `eig_occ.txt` from a converged SCF calculation. It keeps
that source fixed while solving a new K-point list or additional target bands.
No new file format is required.

## SCF

Use a complete, uniformly weighted source mesh and:

```text
calculation scf
suffix hybrid
basis_type pw
dft_functional HSE
smearing_method fixed
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

## NSCF

Keep the same structure, pseudopotentials and cutoffs. In a separate output
directory, retain the shared settings above, replace `KPT` with the target
points, and change:

```text
calculation nscf
read_file_dir ../scf/OUT.hybrid
init_wfc random
nbands 20
out_band 1
out_wfc_pw 0
out_chg 0
```

Supports CPU/CUDA and `nspin 1/2`; use one MPI rank for CUDA. ACE, SOC,
unscreened exchange, forces and stress are unsupported in this workflow.
Set `symmetry -1`: input `0` is accepted but PW EXX resets it to `-1` at runtime.

HSE defaults to `exx_singularity_correction limits`, as in LCAO. Set `gygi` to
retain the historical PW correction, which can produce jumps at source K points.
Use matching SCF/NSCF settings for consistent bands. Do not use Methfessel-Paxton
smearing (`mp`/`mp2`/`mp3`) for the source: it can write occupations outside
`[0, k-point weight]`, which hybrid NSCF rejects; the non-negative schemes
(`fixed`, `gauss`, `mv`, `fd`) are accepted. Configuration differences
or missing `INPUT.info` produce warnings; missing or incompatible numerical
source data produces an error.

## Verification

With MPI and LibXC enabled, substitute your configured executable/build paths:

```bash
cmake --build build --target MODULE_PW_pw_hybrid_source
OMP_NUM_THREADS=1 ctest --test-dir build -R '^MODULE_PW_pw_hybrid_source$'
# Frozen-source band regressions: tests/01_PW/nscf_hse_same and nscf_hse_path
cd tests/01_PW
OMP_NUM_THREADS=1 bash ../integrate/Autotest.sh -a <abacus> -r '^nscf_hse'
```
