# Fixed-Cell i-PI Hardening Design

## Objective

Harden the fixed-cell ABACUS i-PI socket calculator without adding variable-cell
behavior. The work is isolated on `feature/fixed-cell-ipi-hardening`, based
exactly on commit `8b60f83c3e62af75ad57c4c6c61a52a8acdc4d60`. The existing
`feature/fixed-cell-ipi-socket` pull-request branch remains unchanged until the
hardening branch has passed the agreed verification and the user explicitly
authorizes integration.

The supported scientific scope is fixed-cell SCF evaluation for PW and LCAO on
CPU and GPU. The calculator must return trustworthy energy and forces. It may
return stress only when real ABACUS stress calculation is explicitly enabled;
it must never manufacture a zero stress result and advertise it as computed.

## Non-Goals

- Do not update the cell received in `POSDATA`.
- Do not add ASE `UnitCellFilter`, `FrechetCellFilter`, NPT, or barostat support.
- Do not merge or rebase the variable-cell feature into this branch.
- Do not redesign unrelated relaxation, ESolver, MPI, or ASE interfaces.
- Do not install and rebuild an independent ABACUS dependency stack when the
  machine's module-provided libraries and existing builds can be reused.

## Approach

Implement the hardening directly on the fixed-cell baseline using test-driven
development. The variable-cell branch may be consulted as an already-tested
reference, but changes are introduced as fixed-cell-specific commits with a
failing test preceding each production change. This gives the original PR a
small, reviewable follow-up series and prevents variable-cell control paths from
leaking into the fixed-cell implementation.

The implementation is divided into four narrow responsibilities:

1. `IpiSocket` owns raw wire encoding, exact I/O, payload-size arithmetic, and
   transport deadlines.
2. A small socket-frame validation unit owns finite-value, shape, cell inverse,
   determinant, condition-number, force, energy, stress, and virial contracts.
3. `Socket_Driver` owns the explicit protocol state machine, fixed-cell policy,
   ESolver calls, result publication, and MPI failure coordination.
4. `AbacusSocketIO` owns ASE property advertisement, atom-order mapping, process
   launch, and client-side validation of returned arrays.

Dependencies are passed explicitly. No new `GlobalV`, `GlobalC`, or `PARAM`
control dependency is introduced in the socket code.

## Wire Contract

The binary protocol uses:

- exactly 12 bytes for padded ASCII headers;
- `std::int32_t` for replica identifiers, payload lengths, atom counts, and
  extra-byte lengths;
- 8-byte IEEE-754 binary64 for cell, inverse cell, positions, energy, forces,
  and virial;
- native byte order, matching the official i-PI/NumPy implementation. Mixed-
  endian peers are outside the protocol compatibility target and are documented
  as unsupported.

Compile-time assertions enforce integer width, double width, and IEEE-754
support. All element-count-to-byte-count conversions are checked before
allocation or I/O. The atom count must equal `ucell.nat` before coordinates are
allocated or read. INIT data is capped at 1 MiB because ABACUS does not consume
the initialization payload.

`IpiSocket` uses bounded connect, receive, and send operations. Timeout values
are supplied explicitly by `Socket_Driver`; they are not mutable workflow state
inside the socket. An operational environment variable,
`ABACUS_SOCKET_TIMEOUT_SECONDS`, selects a positive finite timeout and defaults
to 300 seconds. It applies only while waiting for socket transport, not while
ABACUS is running SCF. The existing documentation is updated with the default,
accepted values, and failure behavior.

## Protocol State Machine

Use an explicit local enum with three states:

```text
NeedInit --INIT--> Ready --POSDATA/compute--> HaveData --GETFORCE--> Ready
```

`STATUS` is valid in every state and returns `NEEDINIT`, `READY`, or `HAVEDATA`.
`EXIT` is accepted in `NeedInit` or `Ready`. A clean peer close is accepted only
when no computed frame is pending. The following are protocol errors:

- INIT outside `NeedInit`;
- POSDATA outside `Ready`;
- GETFORCE outside `HaveData`;
- EXIT or peer close while a computed frame is pending;
- unknown or malformed headers.

No default, stale, or partially computed frame is published after an error.

## Numerical and Scientific Contract

Before changing ABACUS state, validate the complete input frame:

- all cell, inverse-cell, and position entries are finite binary64 values;
- the cell determinant is finite and positive;
- the cell condition number is below `1.0e12`;
- the received inverse is consistent with the received cell using combined
  absolute and condition-scaled tolerances;
- the atom count equals the STRU atom count;
- the received fixed cell matches the initial ABACUS cell using a combined
  absolute/relative tolerance that remains strict at normal cell sizes.

The wire matrix is the i-PI column-vector matrix `H`. ASE sends its row-vector
matrix `A` as `A.T`, while official i-PI sends `H` directly. The same physical
cell therefore obeys `A = H.T`. Tests use a nonsymmetric triclinic cell so a
mistaken transpose cannot pass accidentally.

After `runner()`:

- publish a finite KS frame even when SCF did not converge, and set
  `scf_converged=false` in the i-PI extras metadata so the external optimizer
  or integrator can apply its own acceptance policy;
- require finite energy;
- require force shape `nat x 3` and finite force entries;
- when stress is enabled, require a finite `3 x 3` stress tensor, positive
  finite volume, and a sufficiently symmetric stress tensor;
- convert ABACUS Ry-based energy, force, and stress to i-PI Hartree/Bohr wire
  units with explicit constants;
- construct wire virial with the verified sign and transpose convention;
- mark `HaveData` only after the full output frame passes validation.

The initial POSDATA atom order must be PBC-equivalent to STRU order. Because the
wire protocol carries no species, a mismatch is a fatal contract violation for
raw i-PI clients rather than a warning. The ASE wrapper continues to sort atoms
into STRU species order and maps returned forces back to the caller's order.

## ASE Property Semantics

`AbacusSocketIO` defines instance-level `implemented_properties`:

- default fixed-cell mode: `energy`, `free_energy`, and `forces`;
- fixed-cell mode with `inp['cal_stress'] = 1`: the same properties plus
  `stress`.

When stress is enabled, `_socket_inp` preserves/enforces `cal_stress=1`, the C++
driver computes real ABACUS stress, and the wrapper converts the returned virial
to ASE stress. When it is disabled, asking ASE for stress raises
`PropertyNotImplementedError`; zero virial is not exposed as a physical result.

The wrapper validates result energy, forces, and virial for shape and finiteness
before updating `self.results`. Python `assert` is not relied upon for runtime
validation because optimized Python may disable assertions.

## MPI Failure Handling

Only rank zero performs socket I/O. Root I/O failures are broadcast before any
rank advances. For each collective computation stage—position application,
ESolver execution, energy/force/stress extraction, and result validation—local
exceptions are collected across ranks before the next collective operation.

If one rank fails while others may be inside or approaching an incompatible
collective, terminate the socket calculation consistently using the established
ABACUS MPI failure mechanism. Focused MPI tests inject failures on root and
non-root ranks and require bounded termination rather than a hang.

## Test Strategy

Every behavior change follows red-green-refactor and is committed separately.
The required automated coverage is:

1. Wire unit tests: exact 12-byte headers, exact four-byte integers, exact
   eight-byte doubles, split TCP/UNIX payloads, clean and partial closes,
   overflow rejection, INIT-size rejection, and transport timeout.
2. Frame unit tests: NaN/Inf, singular/negative/ill-conditioned cells,
   inconsistent inverse, atom-count overflow/mismatch, changed fixed cell,
   invalid energy/force/stress, stress symmetry, unit conversion, virial sign,
   and transpose convention.
3. Driver tests: complete valid sequence and every invalid state transition,
   nonconverged SCF, child cleanup, unknown/EXIT headers, pending-frame close,
   and no publication of invalid or stale frames.
4. MPI tests: root and non-root injected failures terminate within a fixed
   deadline and report the failing stage.
5. ASE tests: property advertisement, stress gating, finite/shape validation,
   species sorting, force remapping, fixed-cell rejection, and process cleanup.
6. Golden-wire tests: ASE and official i-PI 3.2.0 produce the same wire `H` for
   the same nonsymmetric physical cell `A = H.T`.
7. Runtime calculations, with `OMP_NUM_THREADS=1`: one-frame and short repeated-
   frame comparisons against ABACUS FileIO for CPU PW and CPU LCAO, followed by
   one-GPU PW and one-GPU LCAO checks. Energy, forces, and optional stress use
   explicit tolerances recorded with each result.

Runtime builds reuse the locally available ABACUS development modules and linked
libraries. Only the changed targets are rebuilt where CMake permits it. MPI and
socket runtime tests run outside the restricted sandbox.

## Verification Gates

The hardening branch is ready to propose for integration only when all of the
following are true:

- focused C++/Python/MPI tests pass from a fresh build of the changed targets;
- existing socket and INPUT parser tests pass;
- CPU PW and LCAO energy/force comparisons pass;
- official i-PI 3.2.0 and ASE nonsymmetric-cell contract tests pass;
- one-GPU PW and LCAO socket calculations pass, or an external hardware blocker
  is reported without claiming completion;
- stress is either verified as real with sign/unit/format checks or is correctly
  unavailable;
- `git diff --check`, ABACUS agent governance checks, and relevant pre-commit
  checks pass or have concrete review rationales;
- the worktree is clean after committed changes;
- the original PR branch SHA remains unchanged.

After these gates, push only `feature/fixed-cell-ipi-hardening` and report its
commit SHA and exact verification evidence. Merging or pushing into
`feature/fixed-cell-ipi-socket` requires a separate explicit user authorization.
