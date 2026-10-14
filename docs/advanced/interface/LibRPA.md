# LibRPA reader-v1 producer interface

For an LCAO RPA producer run, set `out_librpa_ver = 1` and use a uniform
Monkhorst-Pack grid. The producer writes the full-grid to irreducible-grid
mapping in `bz_sample.txt`; the reduced k-point weights must agree with the
multiplicity of that mapping.

When `symmetry = 1`, `stru_out.txt` contains the spatial operations that
preserve the requested Monkhorst-Pack mesh. The unitary operations are written
first. For `nspin = 4` magnetic calculations, the mesh-preserving spatial parts
of the antiunitary operations follow them, and an optional `spin_symmetry`
block records the antiunitary flag and the eight real/imaginary components of
the corresponding SU(2) rotation.

The structure and reciprocal lattice blocks use Bohr and inverse Bohr. The
scales are derived from the physical `lat0`, so non-default lattice constants
do not silently use the one-Angstrom conversion. Text files use `.txt`; binary
payloads use `.dat`.

The reader-v1 contract is producer-side metadata and serialization only. A
successful file check does not by itself establish equivalence of a symmetry-
reduced RPA or GW result; that requires a numerical comparison with identical
inputs and a compatible LibRPA reader.
