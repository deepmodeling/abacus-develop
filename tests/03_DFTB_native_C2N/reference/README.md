# DFTB+ C2N DFTB3 band reference

`band.out` is the DFTB+ release 25.1 frozen-SCC band-path reference for the
C2N-h2D structure in the parent directory. It uses the 3ob-3-1 parameter set,
the same Hubbard derivatives and 20 K electronic temperature as the ABACUS
case, and 121 Γ–M–K–Γ points (40 intervals on each segment). DFTB+ writes the
reference eigenvalues to 1 meV precision, which limits the precision of direct
band-by-band comparisons.

The reference contains 144 bands at each path point. Against the recorded
ABACUS result, the maximum absolute difference is 5.127 × 10⁻⁴ eV, the mean
absolute difference is 2.49089 × 10⁻⁴ eV, and the RMS difference is
2.87890 × 10⁻⁴ eV. `comparison_summary.txt` records these values and the number
of compared eigenvalues.
