# LCAO DeltaSpin force normalization regression

This periodic distorted Fe2 case checks the **DeltaSpin contribution**, printed
with `test_force 1`, rather than the net-force-shifted total atomic force.
It uses SOC (`nspin 4`), no Hubbard U, 100 Ry, a 2x2x2 k mesh, and target
moments (0,0,2.2) Bohr magnetons on both atoms. The complete input is included.

Run from the repository root:

```sh
python3 tests/integrate/test_deltaspin_force.py \
  --abacus /absolute/path/to/abacus --work-dir /absolute/path/to/results \
  --device cpu --launcher 'mpirun -np 4'
```

For a CUDA executable, `--device gpu --launcher 'mpirun -np 1'` selects
`cusolver`; the CPU test selects `scalapack_gvx`. The runner sets
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`, uses a fresh calculation
directory, and requires all of:

- An SCF-converged marker.
- Charge residual < 1e-10.
- Magnetic RMS < 1.01e-10 Bohr magnetons (allowing printed rounding).
- Absolute final unmixed energy change < 1e-9 eV.
- Each of the six DeltaSpin force components within 1e-6 eV/Angstrom of
  the independent derivative reference below.

`result.json` retains the energy, convergence data, measured component forces,
reference, maximum error and pass/fail result. The exit code is nonzero on
convergence, execution, parsing or force-comparison failure. CTest registers
`lcao_deltaspin_force` for MPI + LCAO + Libxc builds (four CPU ranks).

## Why the reference contains a factor of two

At fixed density matrix D and constraint field lambda, the constraint trace
contains `lambda_a D^a_mu,nu B_mu B_nu`, where B is an orbital-projector
overlap. Displacing an atom differentiates **both** B factors.
`cal_force_IJR` accumulates one factor's derivative; its Hermitian partner
contributes equally after contraction. The final factor of two is unrelated
to spin degeneracy. The Pauli density representation does not remove it.

An independent diagnostic on upstream develop `97697c26` plus this fix held
the actual converged D and lambda fixed, recomputed overlaps at displaced
coordinates, and centrally differenced the trace. It also evaluated both
analytical product derivatives explicitly, without using the production
force-normalization loop. For Fe1 x:

| Quantity | eV/Angstrom |
| --- | ---: |
| Explicit analytical derivative of both overlaps | -0.003514243252157 |
| Fixed production contribution | -0.003514243252156 |
| Frozen-trace central difference, h=1e-4 Bohr | -0.003514237992320 |
| Frozen-trace central difference, h=2e-4 Bohr | -0.003514245802946 |

The other components follow this fixture's x/y and atom-exchange symmetries.
The checked matrix is recorded in `REFERENCE_FORCE` in the runner. The
1e-6 eV/Angstrom tolerance is much smaller than the missing contribution
(about 1.76e-3 eV/Angstrom). Restoring the original `nspin != 4` condition
must fail this test; the field, density and total energy are unchanged.

## Assets and coverage boundary

- `Fe.upf` is the existing fully relativistic repository pseudopotential.
- `Fe_gga_8au_200.0Ry_4s2p2d1f.orb` is the exact numerical basis supplied with
  the reproducer, added under `tests/PP_ORB`. SHA256:
  `8ff6a789b1c7636a1d8ce857edb2477aa094aa36c25cdad11fce646b852c484e`.

This is a component regression, not an acceptance test for total forces,
magnetic forces, stress, PW, or all noncollinear orientations. In separate
full-SCF Fe2-x differences with h=0.0005 Angstrom, the fixed upstream-based
binary still showed about 0.014175 eV/Angstrom total-force error on both CPU
and GPU, despite meeting the convergence gates. That residual is unresolved
and is not claimed fixed by restoring this product-rule factor. No existing
reference or force tolerance is relaxed.
