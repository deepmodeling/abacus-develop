# Molecular Sternheimer AO Potential Exporter

`export_ao_potentials.py` is a standalone Python 3.6+ / NumPy CLI. It exports
the analytic AO potential input for `sternheimer_molecular_coulomb isolated_ri`.
It does not run ABACUS or LibRPA, submit jobs, or build native code. NumPy must
already be installed in the selected interpreter. Python 3.6 requires a NumPy
release that still supports it, such as 1.19.5.

## Required Producer Contract

Use an immutable, serial ABACUS producer directory containing:

- `v1_coulomb_full_iq_1_rank0.dat`: full, unshrunk auxiliary Coulomb matrix.
- `v1_Cs_data_0.txt`: real localized-RI coefficients for every directed atom pair.
- Per-atom AO counts, or `basis_wfc_out` with an explicit atom-type sequence.

Only isolated, free-space, real-AO molecular Gamma data are supported. Extra
q/rank files in the producer directory are rejected, as are translated Cs
records, missing pairs, and inconsistent dimensions. Sparse or MPI-split
datasets must not be made to look serial by renaming one shard. This tool does
not merge shards or fill missing records with presumed zeros.

**External identity verification remains mandatory.** The v1 Coulomb header
contains `iq=1`, not its physical q coordinate or kernel identity. Zero Cs
translations do not establish the isolated Coulomb boundary condition. The
confirmation flag is the caller's declaration, not an automatic verification.
The source binary header also has no unit tag: the full Coulomb payload must
already be in Hartree, confirmed from the producer implementation/provenance.
The exporter performs no unit conversion and cannot detect a mislabeled Ry input.
Check that V, Cs, atom order, AO/ABFS order, PP, orbitals, geometry, kernel and
producer executable belong to the same calculation. At consumption, retain
that identity; do not rotate, shrink or whiten the auxiliary basis.

If a producer KPT is available, it must describe exactly one unshifted Gamma
point: a `1 1 1 0 0 0` Gamma/MP mesh, or one zero-coordinate Direct/Cartesian
point with positive weight. Other KPT forms are rejected. Missing KPT is recorded
and still requires the explicit external declaration.

## Usage

Run from this directory with paths to your own producer and output location.
The following counts are illustrative, not a material preset:

```sh
python3 export_ao_potentials.py \
  --producer /path/to/producer \
  --atom-nao 4 9 --atom-naux 7 12 \
  --confirm-serial-isolated-gamma \
  --output /path/to/response/ao_potentials.dat
```

`--atom-naux` is optional: the v1 V header supplies the authoritative per-atom
auxiliary dimensions. When given, the arguments must agree with that header.
All per-atom lists follow the producer's global atom order.

For split ABACUS basis metadata:

```sh
python3 export_ao_potentials.py \
  --producer /path/to/producer \
  --basis-wfc /path/to/producer/basis_wfc_out \
  --basis-aux /path/to/producer/basis_aux_out \
  --atom-types 1 2 \
  --confirm-serial-isolated-gamma \
  --provenance-file /path/to/producer/executable_identity.txt \
  --output /path/to/response/ao_potentials.dat
```

`--atom-types` contains **one-based type indices for every atom**, not type
populations. Basis metadata include per-type sizes and angular shell lists but
do not encode atom order. The tool validates the shell dimensions, type IDs and
total basis count. Explicit count arguments may also be supplied for comparison.
Basis-file and `--provenance-file` paths are relative to the current directory,
unless absolute. The latter option is repeatable and rejects nonexistent files.

After the external identity checks, the response INPUT uses:

```text
sternheimer_molecular_coulomb isolated_ri
sternheimer_ao_potential_file /path/to/response/ao_potentials.dat
```

The ABACUS runtime guards impose additional restrictions; producing this file
does not bypass them or establish response convergence.

## Formats And Numerical Convention

Binary parsing follows the reviewed ABACUS `RPA_LRI.hpp` v1 writers. Little-endian
32-bit integers, 64-bit offsets and IEEE float64 payloads are required.

- V: `<6i` header, per-atom auxiliary counts, then `<iq` pair-index/offset
  records. All upper-triangular atom pairs must occur exactly once. Real and
  complex payload flags are supported. The lower atom blocks are reconstructed
  by conjugate transpose; diagonal blocks are checked without modification.
- Cs: `<iiiqq` header with `ncell=0` and equal record/block counts, then
  `<iiiiidq` records. The double is `max_abs`, **not distance**. All ordered atom
  pairs must occur exactly once, with zero lattice translations. Payload order
  is `(AO on i, AO on j, auxiliary on i)`; `max_abs` is checked against payload.
- Both parsers require contiguous exact offsets, exact payload sizes and EOF.
  Duplicate records, padding, overlaps, truncation and trailing data are errors.

For every directed Cs record, add its payload and its AO-transposed payload to
`C[p,q,nu]` in the auxiliary span of the first atom. This includes onsite:
`C_onsite = block + block.transpose(1,0,2)`. It is only a simple factor of two
when the onsite block is already symmetric. The reverse offsite record provides
the other center's auxiliary contribution; there is no extra global factor.

The contraction is `T[mu,p,q] = sum_nu C[p,q,nu] * V[nu,mu]`. No conjugation,
unit conversion or Hermitian repair is inserted into this contraction. Finite,
nonzero-norm Hermitian V and T are required, at relative Frobenius tolerance
`1e-10`. A complex Hermitian V need not yield Hermitian AO potentials for real Cs;
such a result is rejected, not coerced to real or symmetrized.

The output begins with:

```text
ABACUS_STERNHEIMER_AO_POTENTIALS 1 nao naux Gamma_Hartree
```

The remaining `naux*nao*nao` rows contain real/imaginary float64 values with 17
digits after the decimal point, auxiliary-major and row-major within each AO
matrix. Values remain in Hartree. The response implementation performs Ha-to-Ry
conversion for the right-hand side and uses conjugated Hartree left couplings.

## Provenance And Output Safety

The companion manifest is named `<output>.json` (for example,
`ao_potentials.dat.json`). It records SHA-256 and byte counts for every parsed
input, the exporter source, optional basis metadata and explicit provenance
files. It also inventories available producer INPUT/STRU/KPT, `basis_*`,
`INPUT_*`, PP (`.upf`, `.vwr`, `.psp`, `.psp8`, `.gth`), `.orb`, and `.abfs`
files. INPUT `stru_file`, `kpoint_file`, `pseudo_dir` and `orbital_dir` and STRU
PP/orbital/ABFS references are followed. Missing identity files are listed;
custom or nonstandard provenance should be supplied with `--provenance-file`.
Relative ABFS filenames are resolved directly from the producer working
directory, as in the ABACUS reader; they are not prefixed with `orbital_dir`.

The manifest records metadata, formula, tolerances, units, Python/NumPy
versions, command arguments, assumptions and output hash. Inputs are rehashed
before publication and changes abort the export. Its status is deliberately
`exported_not_physically_validated`, with `external_identity_verified: false`.

Neither output nor manifest may already exist, including dangling symlinks.
Complete temporary files are published by exclusive hard links in the output
directory, so the destination filesystem must support hard links. Publication
is atomic per file, not across both files. A crash can leave a data file without
a manifest or temporary files; that is not a successful export. The tool will
not overwrite such output on rerun. Keep source producers unchanged during use.

Memory scales as dense `O(nao^2*naux + naux^2)` arrays plus input bytes and
contraction workspace. All payload sizes are checked before dense allocation;
this is not a streaming or distributed exporter.

## Tests

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s . -p 'test_*.py' -v
```

Tests generate small binary fixtures and temporary output only under this
directory. They cover factors, heterogeneous atom dimensions, complex V,
Hermiticity, malformed/truncated/duplicate data, metadata, provenance hashes,
Gamma checks, overwrite refusal and output roundtrip. They do not validate a
material, physical convergence, or compatibility with an untested producer build.
