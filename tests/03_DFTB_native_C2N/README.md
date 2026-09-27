# Native periodic DFTB3 C2N regression

This self-contained regression runs an ABACUS C2N-h2D DFTB3 calculation and
checks its frozen-SCC band path against the DFTB+ 25.1 reference in
`reference/band.out`. It exercises the SCC solve, Broyden mixing, DFTB3
correction, and band-path diagonalization without requiring a DFTB+ executable.

The input uses a 2×2×1 Γ-centered integration mesh, a 20 K electronic
temperature, and the Γ–M–K–Γ path with 40 intervals per segment (121 points).
The structure contains 36 C/N atoms. The four C–C, C–N, N–C, and N–N SKF files
are from 3ob-3-1. Their license and required scientific references are included
in `parameters/LICENSE` and `parameters/README`.

Run the regression after configuring and building ABACUS with testing enabled:

```bash
ctest --test-dir <build-directory> --output-on-failure \
  -R '^dftb_native_c2n_reference$'
```

The reference contains 17,424 eigenvalues (121 k-points × 144 bands) rounded
to 1 meV by DFTB+. The recorded comparison gives a maximum absolute difference
of 5.127 × 10⁻⁴ eV, a mean absolute difference of 2.49089 × 10⁻⁴ eV, and an RMS
difference of 2.87890 × 10⁻⁴ eV. The automated check allows a maximum
difference of 7.5 × 10⁻⁴ eV and an RMS difference of 3.5 × 10⁻⁴ eV to account
for the reference output rounding.

This regression covers the specified C/N model and calculation only. It does
not establish support for forces, stress, relaxation, molecular dynamics,
spin polarization, or other parameter sets and materials.
