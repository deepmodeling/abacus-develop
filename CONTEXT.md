# ABACUS numerical layout

This context describes how LCAO basis identity, distributed storage, and
operator topology are represented without a shared omnibus object.

## Language

**BasisFlattenedIndexing**:
The canonical atom-major ordering of LCAO basis functions and the mapping
between flattened indices, atoms, and basis labels.
_Avoid_: AO indexing, orbital partition, basis context

**BasisRange**:
A half-open range of flattened basis indices belonging to one atom.
_Avoid_: atom range, basis view

**BlockCyclicDistributionInformation**:
The process grid, block size, descriptors, and global/local index mapping for
one distributed matrix layout.
_Avoid_: parallel layout, matrix context

**MO representation layout**:
The block-cyclic layout of `S` or `H` represented in the set of MOs defined by
the DFT reference determinant. Its square global extent follows from the number
`nmo` of MOs.
_Avoid_: band matrix, Eij layout, reference-state matrix

**Parallel_Orbitals**:
A historical LCAO migration adapter that combines basis indexing,
wavefunction shapes, atom traces, and a matrix distribution.
_Avoid_: LCAO context, LcaoLayout

**HContainer**:
A block-sparse representation of operators keyed by atom pairs and lattice
translations.
_Avoid_: dense matrix, global matrix

**Matrix2d**:
The self-contained mathematical matrix interface that owns its values and
distribution while exposing views and backend-neutral operations.
_Avoid_: matrix wrapper, matrix context
