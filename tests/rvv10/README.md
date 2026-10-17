# libvdwxc-backed rVV10 checks

The nonlocal kernel comes from libvdwxc MR !32, revision
`94410ebdf3dc9615a40f40c3a5b59dbe20a433d2`. Installation, the FFT-provider
constraints and the supported input scope are documented in
`docs/advanced/install.md`; `abacus -h xc_nonlocal` prints the same contract.

## Dependency-free build-report test

`RVV10_LIBVDWXC_CONTRACT` runs the legacy Makefile build-info generator and
checks that an installed library is never reported as a compiled feature:

```bash
cmake -B build -DBUILD_TESTING=ON
ctest --test-dir build -R RVV10_LIBVDWXC_CONTRACT --output-on-failure --no-tests=error
```

## Backend tests

With the pinned libvdwxc installation, configure `build-libvdwxc` with the
standalone-FFTW command in `install.md`, then:

```bash
cmake --build build-libvdwxc --target MODULE_HAMILT_vdwTest_LIBVDWXC \
  MODULE_ESTATE_pot_nonlocal_xc_LIBVDWXC
ctest --test-dir build-libvdwxc -R LIBVDWXC --output-on-failure --no-tests=error
```

`MODULE_HAMILT_vdwTest_LIBVDWXC` checks the adapter against the serial library
API, the potential against a finite difference of the energy, the cell mapping
with equivalent oblique and orthogonal cells, handle ownership and the standard
error diagnostics. `MODULE_ESTATE_pot_nonlocal_xc_LIBVDWXC` drives the host
component with real `Charge` and `UnitCell` objects and checks valence/core
density construction, additive energy/potential/`vtxc` bookkeeping and
unchanged input arrays. The Ry conversion and `vtxc` convention are documented
next to the implementation in `nonlocal_xc_libvdwxc.cpp`.
