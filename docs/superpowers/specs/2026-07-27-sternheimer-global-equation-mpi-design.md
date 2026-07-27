# Sternheimer SIAB Global-Equation MPI Design

## Scope

This change removes the static imaginary-frequency load imbalance in the
molecular SIAB Delta-Sternheimer producer. It adds a deterministic MPI layout
in which all ranks receive a mixture of response equations from all imaginary
frequencies.

The change is limited to `out_sternheimer_siab=true`. It does not change the
Delta-Sternheimer equation, GMRES solver, convergence tolerance, auxiliary
basis, frequency grid, Coulomb whitening, SIAB v1 format, or LibRPA chi0
output. A kinetic preconditioner is a separate follow-up because combining
solver and scheduling changes would prevent an attributable performance
comparison.

## Motivation

The current 20 A, 50 Ry H2 target uses 16 frequencies, 428 auxiliary channels,
and 32 MPI ranks. Static assignment gives two ranks to every frequency. The
average GMRES iteration count ranges from 11.93 at the highest frequency to
134.83 at the lowest frequency. One high-frequency rank finishes its local
equations after 1.49 h, while the critical low-frequency rank finishes after
13.69 h. The estimated static-layout occupancy is only 53 percent because
finished ranks wait at the global SIAB row gather.

The SIAB response equations are independent before the final row gather. There
is no mathematical requirement that one rank process only one frequency.

## Reference Equation and Current Solver

For the external Delta-Sternheimer subspace, the current implementation solves

\[
 {\cal A}_{Q}(i\omega)x=b_Q,
\]

with

\[
 {\cal A}_{Q}(i\omega)
 =
 Q(H-\epsilon_i+i\omega)Q
 \sum_a
 \frac{|D_a\rangle\langle D_a|}
 {\epsilon_i-\epsilon_a-i\omega}.
\]

The right-hand side is

\[
 b_Q
 =
 Qb
 -
 \sum_a
 \frac{\langle a|\delta V|i\rangle}
 {\epsilon_i-\epsilon_a-i\omega}|D_a\rangle.
\]

The finite-dimensional Delta subspace is reconstructed analytically after the
external response is found. The external equation is solved by restarted
GMRES with restart dimension 50. The production maximum is 300 iterations and
the current production residual tolerance is \(10^{-8}\).

Neither the standard nor Delta-Sternheimer production path currently assigns
`LinearProblem::precondition`; GMRES therefore uses the identity
preconditioner. The existing kinetic-preconditioner helper is not connected to
this solve.

At large \(\omega\), the \(i\omega I\) shift dominates and clusters the
operator spectrum. At small \(\omega\), near-gap modes and the broad 50 Ry
finite-difference kinetic spectrum remain exposed. The resulting condition
number is larger and GMRES requires more Krylov vectors. This explains the
observed low-frequency cost without invoking iterative diagonalization.

## Input and Compatibility

Add the string input

```text
sternheimer_mpi_layout frequency_grouped
```

with two accepted values:

- `frequency_grouped`: current behavior and default;
- `global_equation`: deterministic response-equation ownership over all MPI
  ranks.

`global_equation` requires:

- `sternheimer_frequency_mpi=true`;
- `sternheimer_channel_mpi=true`;
- `out_sternheimer_siab=true`;
- at least one MPI rank.

Unlike `frequency_grouped`, it does not require the MPI rank count to be an
integer multiple of the frequency count. The LibRPA chi0 writer continues to
reject channel-distributed layouts.

## Deterministic Ownership

For global occupied-state index \(i\), zero-based frequency index \(f\),
auxiliary channel \(\mu\), frequency count \(N_\omega\), channel count
\(N_\mu\), and MPI size \(P\), define

\[
 t(i,f,\mu)=((iN_\omega+f)N_\mu+\mu),
\qquad
 p_{\mathrm{owner}}=(t+s)\bmod P,
\]

where \(s\) is the existing normalized rank shift.

Every rank traverses the same occupied-state, frequency, and channel loops but
solves only equations for which it is the owner. This gives every rank
low- and high-frequency work. Equation counts differ by at most one when all
tasks are included.

The mapping is static and reproducible. It adds no master-worker messages and
no synchronization between frequencies. Each rank continues to hold the same
Hamiltonian, occupied projector, Delta subspace, perturbations, and SIAB
primitive data that it already holds in the grouped implementation.

## Output and Diagnostics

Each solved equation produces the existing `siab::ReferenceRow`, including its
frequency index. The existing MPI gather sends all local rows to rank 0.
Before validation and SIAB v1 writing, rank 0 sorts rows canonically by

```text
(frequency_index, occupied_state, auxiliary_channel)
```

and then applies the existing missing-or-duplicate row-count check. The v1
physics and metadata remain unchanged.

The status file adds:

```text
sternheimer_mpi_layout global_equation
equation_owner_formula occupied_frequency_channel_modulo
```

Per-rank progress retains equation counts, iterations, residuals, and elapsed
time. The summary records the minimum and maximum rank-local equation counts
and iteration sums so load balance can be evaluated without parsing every
progress file.

Invalid layout values or incompatible output modes fail before constructing
the grid Hamiltonian.

## Verification

Implementation follows test-driven development.

1. Add a failing unit test for global ownership over several combinations in
   which MPI ranks are and are not multiples of the frequency count.
2. Verify that every response equation has exactly one owner and that total
   rank-local equation counts differ by at most one.
3. Verify that `frequency_grouped` ownership remains unchanged.
4. Add parser tests for both accepted input values and invalid-value rejection.
5. Run a small H2 SIAB target with 1, 2, and 4 MPI ranks. Require identical row
   keys, equation counts, convergence flags, and provenance.
6. Compare the canonical SIAB target matrices against the grouped baseline.
   Require relative Frobenius errors of the response target and overlap
   matrices below \(10^{-10}\), and preserve the existing solver and explicit
   equation residual gates.
7. Run the full 20 A, 50 Ry, 16-frequency H2 target on 32 normal-partition
   nodes using the same executable, input files, frequency grid, and
   tolerances. Record wall time, node-hours, MaxRSS, per-rank equation and
   iteration loads, and target-matrix errors.

The performance success criterion is at least \(1.5\times\) lower wall time
than the 14:06:59 grouped baseline without increasing node count or weakening
any numerical gate. The expected 7--8 h range is an estimate, not a pass
criterion.

## Follow-Up: Solver Preconditioning

After the scheduling A/B is complete, add and validate a separate
finite-difference kinetic or shifted-Laplacian right preconditioner for the
projected Delta operator. That work must use its own commit and same-layout A/B
tests. It should target the low-frequency iteration count and node-hour cost,
whereas the present change targets idle time caused by static frequency
ownership.
