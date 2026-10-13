# Density-grid nonlocal exchange-correlation

ABACUS treats a density-grid nonlocal correction as an additive component next
to the ordinary semilocal XC term. The first implementation uses the rVV10
backend from libvdwxc; ABACUS does not carry a second rVV10 kernel.

## Input and declared scope

The input is deliberately separate from the semilocal functional:

```text
dft_functional  GGA_X_RPW86+GGA_C_PBE
xc_nonlocal     rvv10
rvv10_b         6.3
rvv10_c         0.0093
```

The first adapter requires both Libxc and libvdwxc at build time. It accepts a
CPU, double-precision, full-complex plane-wave SCF calculation with one MPI
rank, `nspin=1`, and norm-conserving pseudopotentials. It returns the nonlocal
energy, the local potential, and the valence `vtxc` contraction. The parser
rejects forces, stress, relaxation, molecular dynamics, LCAO, GPU, gamma-only
layouts, spin-polarized inputs, ultrasoft/PAW pseudopotentials, MPI runs with
more than one rank, and pairwise D2/D3/D4 corrections. These are follow-up
interfaces, not silently approximated results.

`xc_nonlocal=rvv10` must not be encoded as `dft_functional=RVV10`, and it must
not be combined with a semilocal functional that already contains a nonlocal
vdW term. The existing pairwise `vdw_method` path remains independent and is
mutually rejected by the input check.

## Host/backend boundary

`Potential` registers the existing semilocal `xc` component first. If
`xc_nonlocal` is enabled, it then registers the generic `xc_nonlocal`
component. `PotNonlocalXC` only adds the backend's energy, potential and
`vtxc`; it never replaces or recomputes the semilocal term.

The host-side interface lives in
`source/source_hamilt/module_vdw/nonlocal_xc.h`. The only rVV10-specific
implementation in ABACUS is
`nonlocal_xc_libvdwxc.{h,cpp}`. The adapter reuses
`PW_Basis::real2recip`, `XC_Functional::grad_rho` and
`XC_Functional::grad_dot`; libvdwxc owns the kernel table, interpolation and
nonlocal convolution.

For one SCF update the data path is:

```text
rho_total = rho_valence + rho_core
    -> grad(rho_total), sigma = |grad(rho_total)|^2
    -> libvdwxc rVV10 energy, dE/dn and dE/dsigma
    -> ABACUS grad_dot for div(2*dE/dsigma * grad(rho_total))
    -> v_nl in Ry, etxc/vtxc accumulation, and addition to every spin channel
```

The library already includes the rVV10 beta compensation. ABACUS therefore
does not add beta a second time. The total density is used by the backend;
`vtxc` is contracted with the valence density, matching the host XC energy
bookkeeping.

## CMake and dependency contract

Enable the optional backend with a CMake package that exports
`libvdwxc::vdwxc`:

```bash
cmake -B build-libvdwxc -DBUILD_TESTING=ON \
  -DENABLE_LIBXC=ON -DENABLE_LIBVDWXC=ON \
  -DCMAKE_PREFIX_PATH=/path/to/libvdwxc-prefix
```

The current upstream MR is still Draft. The integration was checked against
the `rvv10-cmake` branch at commit
`91ef4820f3270e22b35fddc2ab4cab8efa4b6d79`; a PR must keep this pin (or a
newer reviewed commit), record the dependency's CTest result, and never call
the Draft API a stable release. Configuration verifies the imported target,
the header, and every serial API symbol used by the adapter. Legacy Makefile
builds do not provide this optional dependency; `xc_nonlocal=rvv10` fails with
a clear input diagnostic there.

The feature is included in the CMake build summary as `libvdwxc`. The adapter
is attached to `source_hamilt/module_vdw`, not to the semilocal XC evaluator or
to a special rVV10 branch in the generic SCF solver.

## Numerical and layout contract

The serial API receives an ABACUS full real-space grid with z as the fastest
index. The unit-cell matrix is passed column by column, and the adapter checks
the grid dimensions, cell identity, precision, device, gamma flag, rank count,
density sizes and finite values before calling the library. It forms
`sigma=|grad(n_total)|^2`, converts the library Hartree outputs to ABACUS Ry,
and uses the same reciprocal-space projection for the gradient and divergence.

The adapter is intentionally not thread-safe through a shared evaluator object;
one evaluator belongs to one potential instance. MPI/PFFT ownership, spin
interfaces, force/stress derivatives and GPU/LCAO data layouts need separate
API and numerical tests before being enabled.

## Verification

The host-side Python contract test checks the CMake switch, imported target,
required API symbols and the absence of a duplicate in-tree rVV10 kernel. With
the pinned library installed, the GoogleTest
`MODULE_HAMILT_vdwTest_LIBVDWXC` compares the adapter with a direct serial
libvdwxc call for a nonuniform density and an oblique cell. It checks the
energy, the host-assembled `dedn - div(2*dedsigma*grad(n))` potential, and the
valence `vtxc` contraction.

Before converting a Draft PR to Ready, run both libraries' tests and add an
ABACUS SCF reference case with the exact executable, pseudopotential hashes,
input, dependency commit, compiler, FFT provider, and numerical tolerances.
The current local environment cannot provide a complete CMake/libvdwxc build;
therefore this change remains a Draft until that evidence is attached.

## Follow-up work

The next PRs can add a libvdwxc MPI/PFFT adapter, ordinary spin, GPU/LCAO
paths, forces and stress, and functional-specific parameter tables. Each
extension should add its own layout, unit, finite-difference and reference
tests before loosening the input checks.
