# libvdwxc-backed rVV10 checks

This directory checks the host-side contract of the first ABACUS integration.
The numerical rVV10 kernel is provided by the pinned libvdwxc development
branch; ABACUS does not carry a second kernel implementation.

The host adapter is intentionally limited to the current API contract:

* CPU, double precision, full-complex plane-wave density grid;
* one process and `nspin=1` SCF;
* Libxc semilocal `GGA_X_RPW86+GGA_C_PBE` plus `xc_nonlocal=rvv10`;
* norm-conserving pseudopotentials, energy, potential and `vtxc` only.

MPI, spin-polarized calculations, GPU/LCAO, gamma-only layouts, ultrasoft/PAW
pseudopotentials, forces and stress are rejected by the input contract and are
follow-up work. The current libvdwxc MR is still Draft, so a Ready PR must pin
commit `91ef4820f3270e22b35fddc2ab4cab8efa4b6d79` and show both the library's `ctest -R rvv10` output and the
ABACUS adapter test output.

Configure the host-side contract test with:

```bash
cmake -B build -DBUILD_TESTING=ON
ctest --test-dir build -R RVV10_LIBVDWXC_CONTRACT --output-on-failure
```

With the pinned libvdwxc installation available, enable the real adapter test:

```bash
cmake -B build-libvdwxc -DBUILD_TESTING=ON -DENABLE_LIBXC=ON \
  -DENABLE_LIBVDWXC=ON -DCMAKE_PREFIX_PATH=/path/to/libvdwxc
cmake --build build-libvdwxc --target MODULE_HAMILT_vdwTest_LIBVDWXC
ctest --test-dir build-libvdwxc -R LIBVDWXC --output-on-failure
```

The adapter test uses a nonuniform density and an oblique cell. It compares
ABACUS energy, the host-computed `dedn - div(2*dedsigma*grad(n))` potential,
and the valence `vtxc` contraction against the direct serial libvdwxc API.
