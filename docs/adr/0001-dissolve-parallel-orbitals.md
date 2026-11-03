# Dissolve Parallel_Orbitals into self-contained values

**Status**: accepted

`Parallel_Orbitals` is a migration adapter, not a target interface. Its
responsibilities must separate into `BasisFlattenedIndexing` for canonical
global basis order and atom slices, `BlockCyclicDistributionInformation` for
process/block mapping, the owning matrix or wavefunction for shapes, and
`HContainer` for block-sparse `(I,J,R)` topology. The type is removed when its
last consumer is migrated; no `LcaoContext`, `LcaoLayout`, or similarly broad
replacement is permitted.

## Consequences

Computation interfaces receive only the values they use. Atom-resolved
algorithms can combine a basis index and a distribution without receiving
unrelated descriptor or band metadata. Setup code owns and constructs those
values, and duplicated layout state must not be introduced to ease migration.

A gathered dense tensor is a materialization of matrix values, not a second
definition of the wavefunction distribution. Algorithms such as boundary
correction receive the `Psi` block-cyclic layout and derive a temporary serial
layout when their dense H/S tensors have already been gathered.

## Migration status

RT-TDDFT is the current pilot. `TD_info` explicitly receives and owns the
canonical basis indexing and matrix distribution used by its lifetime, while
RT production interfaces no longer pass `Parallel_Orbitals`. The velocity
current owns its basis indexing and container resources, and the TD kinetic,
nonlocal, hybrid-potential, and GPU projector paths consume the shared RT
values rather than re-reading a descriptor.

The shared LCAO Hamiltonian/overlap containers remain on the compatibility
layout during this stage. The next high-impact boundary is the operator
factory that constructs those containers: switching them to basis indexing
and block-cyclic values requires migrating every consumer of atom-pair
distribution metadata, not only changing the constructor. This migration must
preserve one owner for basis topology and must not introduce a broad context
object.

## Atom-local basis access

Atom-local rows and columns are a composition of basis ownership and matrix
distribution. `AtomBasisView` derives that composition directly from
`BasisFlattenedIndexing` and `BlockCyclicDistributionInformation`. The view is
not cached or owned by `HContainer`; `HContainer` only uses it while
constructing the shapes of local `AtomPair` blocks.

`BasisFlattenedIndexing` remains distribution-free. It must not gain APIs such
as `local_row_offsets(layout, iat)` or otherwise join itself to
`BlockCyclicDistributionInformation`; that would put a per-kernel,
hand-assembled distribution adapter back into the architecture.
Per-element LCAO operations consume `AtomBasisView::rows` and
`AtomBasisView::cols` from the explicit basis indexing and distribution they
already own. Offset-based helpers are not part of the migration interface.

The prototype for this interface was developed on branch
`parallel2d-linalg-container` in commits `cad96508f` and `9d893d3b0`. Those
commits are the first reference to consult for call-site migration patterns.
Their container-owned view cache was later rejected because the same values are
fully derived from basis indexing and distribution. Existing implementations
must still be checked across local branches before changing this seam.

## Working with this migration

Before changing a migration interface:

1. Search all local branches and history for an existing implementation.
2. Check this ADR and the migration status above.
3. Keep the batch to one seam and update this ADR when the seam changes.
4. Record exact build and test commands in the handoff or pull request.

Rejected patterns include new omnibus context objects, per-kernel construction
of atom-local offset vectors, HContainer ownership of derived atom-basis
views, and new manual ownership of BLACS distribution details.
