# Restricted rVV10 regression

Configure the ordinary CPU build with `BUILD_TESTING=ON`, `ENABLE_LCAO=OFF`,
`ENABLE_LIBXC=ON`, then build the executable and the
`MODULE_HAMILT_XCTest_RVV10_CORE`, `MODULE_HAMILT_XCTest_RVV10_PW` and
`MODULE_ESTATE_RVV10_POTENTIAL` targets. No QE installation is needed.
With OpenMP enabled (the default), this configuration produces `abacus_pw_omp`;
with `ENABLE_OPENMP=OFF`, it produces `abacus_pw_ser`. An MPI-enabled PW build
produces `abacus_pw_para`. Use the executable name from the actual configuration.

```bash
OMP_NUM_THREADS=1 ctest --test-dir build -R RVV10 --output-on-failure --no-tests=error
```

The following command runs the serial and MPI SCFs, unsupported-input checks, three
unsupported-nonlocal-functional checks and one actual ultrasoft-pseudopotential
rejection:

```bash
python3 tests/rvv10/run_scf_regression.py --executable /absolute/path/to/abacus_pw_omp \
  --pseudo-dir /absolute/path/to/tests/PP_ORB --mode all
```

The unsupported-input checks include noncollinear spin, force, stress,
relaxation, an additional dispersion method, LCAO and GPU evaluation. The last
case is a parser-level refusal; it does not claim GPU rVV10 support. For PW input, ABACUS globally resets `gamma_only=1` to `0` before
functional checks, while the evaluator still rejects a manually supplied
reduced-gamma `PW_Basis` layout.

Each invocation preserves inputs, output, binary hash
and a JSON summary in a fresh directory. `--artifacts-dir` chooses its parent.
The no-Libxc build uses `--mode no-libxc` to test controlled refusal of
`xc_nonlocal=rvv10`.
`--mode uspp` reads the repository's `Na.pbe-spn-rrkjus_psl.0.2.UPF` before
requiring the specific norm-conserving-only diagnostic. An unrelated error,
crash, sanitizer diagnostic or even one started SCF iteration is a failure;
checking INPUT alone would not test this runtime condition. A sanitizer failure
must not pass merely because it shares exit code 1 with an expected refusal.

The `all` and `reject` modes accept a parser-level PBE+rVV10 combination,
then combine `xc_nonlocal=rvv10` with `GGA_XC_BEEF_VDW`, `GGA_XC_VV10` and
`MGGA_X_SCAN+MGGA_C_SCAN_VV10` through actual initialization. Each of the
latter must be rejected before SCF because those Libxc functionals already
contain a nonlocal term. This protects against silently computing an
unintended composite; a different semilocal base still needs its own validated
`rvv10_b`/`rvv10_c` regression before being advertised.

For an MPI-enabled CPU/Libxc build, CTest registers serial and two-rank rVV10
SCFs, serial and two-rank zero-magnetization `nspin=2` SCFs, a two-rank
noncollinear-spin refusal, and a two-rank He PBE control. It uses CMake's
`MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS` and
`MPIEXEC_POSTFLAGS`; configure launcher options for the test host if needed.
CTest reserves two processors for the two-rank cases. The driver also accepts
`--launch` followed by the complete launcher command as its last option. The
command must include the resolved absolute `--executable` path so that the
recorded binary hash identifies the program being tested. No shell
interpretation is used. The two-rank rVV10 SCF uses the same converged
reference as the serial case; the dedicated spin modes additionally compare
`nspin=2` with `nspin=1` in the same launcher configuration. These checks
exercise the production SCF adapter, not only the component unit test.

The He and Si fixtures use 40/160 Ry cutoffs and 45-cubed/30-cubed grids. The
Si pseudopotential contains NLCC; pseudopotential hashes are checked before
running. References are the converged ABACUS energies from 2026-09-13, checked
against QE commit b3af0d0e810a611c455fd80fca6ee523dfaaee95. The tolerance is
2e-6 Ry for absolute energies and 2e-7 Ry for the RVV10-minus-base SCF change.
The latter includes the self-consistent density response; it is not a
fixed-density nonlocal energy. These are implementation regressions, not
cutoff-converged material predictions.

The native potential test uses real production objects, FFT and Libxc. It
checks the full XC energy derivative with fixed core density, valence-only
vtxc and additive/repeated potential updates. A second test composes the
semilocal `PotXC` and nonlocal `PotRvv10` components in one process and compares
energies, vtxc and every potential value. It uses the real adapters directly,
without changing global INPUT or adding private-access test macros.

The same component test covers `nspin=2`: it sums the two valence channels for
the nonlocal functional and checks that the resulting scalar potential is added
to both spin rows.

The real-FFT tests cover a small-gradient density and a periodic material/vacuum
profile on 6x8x24 and 6x8x32 nonorthogonal grids. The latter checks finite values
throughout the cell, valence-only vtxc and energy finite differences in the
material and dilute tail. At n~2.5e-9, relative perturbations of 3%, 0.9% and 0.3%
avoid cancellation of total energies at excessively small steps; the material
uses 0.1%, 0.03% and 0.01%. All use the same 2e-5 Ry potential tolerance. These
are consistency checks on each grid, not grid-convergence or hard-cutoff
derivative tests. The SCF references catch missing or duplicated nonlocal
energy; the Python validators reject unconverged energies, wrong units,
crashes and contradictory success/rejection output.

For memory checks, configure with `ENABLE_ASAN=ON` and run the same targets
with `ASAN_OPTIONS=detect_leaks=1:halt_on_error=1`. This includes repeated
component lifetimes and the full SCF executable. The real-XC test also repeats
PBE's `output_info()` call: although PBE uses the built-in functional,
Libxc allocates its printable names, which the caller must release.
