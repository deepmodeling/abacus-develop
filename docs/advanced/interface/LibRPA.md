# LibRPA exports with SOC and magnetic symmetry

The LCAO RI exporter writes the ABACUS states and auxiliary-basis data needed
by a compatible LibRPA reader. The reader constructs orbital and k-space
rotation matrices from the spatial operations, lattice, atomic coordinates,
and basis-shell conventions; no `symrot_*` sidecar files are part of this
interface.
Enabling an export does not certify symmetry equivalence of a downstream
RPA/GW calculation.

## Producer settings

For an otherwise complete, converged LCAO SOC input, the relevant settings are:

```text
basis_type                  lcao
nspin                       4
lspinorb                    1
symmetry                    1
rpa                         1
out_ri_cv                   1
out_wfc_lcao                1
out_mat_xc                  1
out_librpa_reader_version   1
```

Set the intended magnetic moments in `STRU`. With `symmetry=1`, an all-zero
initial magnetic configuration stays zero; it is not automatically replaced
by nonzero starting moments. A nonmagnetic SOC system uses a grey group,
whereas an ordered magnetic system uses the operations preserving that
configuration, possibly combined with time reversal. For an on/off comparison,
use a common converged density and keep the structure, basis, k mesh and
numerical thresholds fixed.

## Common symmetry operation block in `stru_out`

The existing physical lattice/position units and spatial `row` convention
are retained. After the atom records, symmetry-enabled output contains:

```text
N row
<9 integer rotation entries and 3 fractional translations, repeated N times>
```

The same spatial-operation format is used for scalar space groups and for the
spatial part of a magnetic group. For `nspin=4`, `stru_out` intentionally does
not append a spin-action or antiunitary flag trailer. A consumer that needs
spinor magnetic-group acceleration must obtain and interpret that metadata
through an explicitly compatible interface; the common `stru_out` format does
not claim that capability.

## Optional auxiliary-overlap diagnostic

To inspect the active, shrunk auxiliary basis, also set:

```text
shrink_abfs_pca_thr          1e-6
out_librpa_abf_overlap       1
```

The threshold is an example, not a convergence recommendation. The option
requires `rpa=1`, reader version `1`, and a nonnegative shrink threshold.
It writes `v1_abf_overlap_active_iq_<iq>.dat` in the same active auxiliary
basis as the associated Coulomb matrices. It does not export the full
unshrunk overlap. Dense matrices can require substantial memory and disk.

From the ABACUS source tree, inspect a completed export with:

```bash
python3 tools/analyze_librpa_abf_overlap.py /path/to/export
python3 tools/test_analyze_librpa_abf_overlap.py
```

The analyzer checks metadata and matrix properties. The writer rejects
duplicate MPI contributors for present `(I,J,R)` blocks. Neither check
establishes complete real-space coverage or physical convergence.

## Occupation precision and validation

`band_out` writes the Fermi energy in Hartree and each occupation as the
native k-weighted occupation multiplied by the number of exported k points.
The stream precision is initialized before these values, including the first
occupation. Previously the first occupation and Fermi energy retained the
stream's default six significant digits; removing an irreducible k weight
could turn an occupied state into `0.999999` or `1.000002`.

An export check can divide the printed occupation by the exported k-point
count and compare it with `OUT.<suffix>/eig_occ.txt`. Validate the first row
as well as later rows. A successful LibRPA reader/EXX run establishes that the
files can be consumed; full/reduced-grid RPA or GW equivalence requires a
separate numerical comparison. This producer check does not change LibRPA's
archived-input regression workflow.
