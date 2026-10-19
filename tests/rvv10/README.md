# libvdwxc-backed rVV10 checks

The numerical kernel is provided by libvdwxc MR !32, revision
`94410ebdf3dc9615a40f40c3a5b59dbe20a433d2`. This revision fixes FFT boundary
symmetry in oblique cells. Its eight library tests pass with both standalone
FFTW 3.3.11 and the MKL 2022.2.1 compatibility backend (GCC 11.2.1).
The ABACUS adapter and host tests are verified with standalone FFTW, not an
MKL-based ABACUS build. Do not mix both FFT providers in one process.
Installation and supported inputs are documented in
`docs/advanced/install.md`.

The host adapter is intentionally limited to the current API contract:

* CPU, double precision, full-complex plane-wave density grid;
* one process and `nspin=1` SCF;
* Libxc semilocal `GGA_X_RPW86+GGA_C_PBE` plus `xc_nonlocal=rvv10`;
* norm-conserving pseudopotentials, energy, potential and `vtxc` only.

MPI, spin-polarized calculations, GPU/LCAO, gamma-only layouts, ultrasoft/PAW
pseudopotentials, forces and stress are rejected by the input contract and are
follow-up work. The initial integration remains Draft pending independent QE
reference data and maintainer agreement on this staged scope.

The dependency-free test executes the Makefile build-info generator and
checks that an installed library cannot be mistaken for a compiled feature:

```bash
cmake -B build -DBUILD_TESTING=ON
ctest --test-dir build -R RVV10_LIBVDWXC_CONTRACT --output-on-failure
```

With the pinned libvdwxc installation available, enable the real adapter test:

```bash
# Configure build-libvdwxc with the standalone-FFTW command in install.md.
cmake --build build-libvdwxc --target MODULE_HAMILT_vdwTest_LIBVDWXC \
  MODULE_ESTATE_pot_nonlocal_xc_LIBVDWXC
ctest --test-dir build-libvdwxc -R LIBVDWXC --output-on-failure
```

The C++ tests check energy and potential against the serial API, the potential
by a finite difference of energy, and geometry using equivalent oblique and
orthogonal cells with an analytically matched reciprocal vector. Negative
tests exercise the standard ABACUS error diagnostics.
The adapter cannot be copied because it owns its library handle. A separate
host test uses real `Charge` and `UnitCell` objects to check valence/core density
construction, additive energy/potential bookkeeping and unchanged input data.

The host reuses `PW_Basis`, `grad_rho` and `grad_dot` to supply
`n = n_valence + n_core` and `sigma = |grad(n)|^2` to libvdwxc. The library
includes the beta compensation. The returned derivatives form
`v = 2*(dedn - div(2*dedsigma*grad(n)))` in Ry; `vtxc` contracts this potential
with the valence density. The generic host component adds these quantities
after `PotXC`, while the numerical backend stays in `module_vdw`.
