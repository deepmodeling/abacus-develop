# Wannier90 spin-matrix regression

The oxygen inputs in this directory exercise formatted `.spn` output from a
spinor LCAO calculation. The driver runs SCF, NSCF, Wannier90 and `postw90.x`
in a temporary directory using the pseudopotential and orbital files in
`tests/PP_ORB`. [Test driver](../../../interfaces/Wannier90_interface/tests/check_wannier_spn.py)

## Run the test

Use an MPI-enabled LCAO build with ScaLAPACK (`ks_solver=scalapack_gvx`).
Configure with `BUILD_TESTING=ON` and make sure CMake can find an MPI launcher,
Python 3 with NumPy, `wannier90.x` and `postw90.x`. From the repository root,
run the four-rank CTest case using your build directory.

```sh
ctest --test-dir build -R '^03_NAO_wannier_spn$' --no-tests=error --output-on-failure
```

Add `-V` to see the driver's output for a passing CTest run.
To run the driver directly with the executables on `PATH`, use the following
command from the repository root. Pass the MPI launcher and its arguments as
separate words. Use `-n 1` for a single MPI rank.

```sh
OMP_NUM_THREADS=1 python3 interfaces/Wannier90_interface/tests/check_wannier_spn.py \
    wannier90.x postw90.x mpirun -n 4 abacus
```

## Read the results

The driver reconstructs the three Pauli spin matrices from the ABACUS
wavefunctions and orbital overlaps. It compares 546 complex upper-triangle
elements (13 retained bands at two k points) with an absolute tolerance of
`1e-10`. It also checks that `postw90.x` reads `oxygen.spn` and reports three
finite spin-moment components.

The NSCF runs also check that spin output is on by default for this case,
that `out_wannier_spn=0` disables it, and that ABACUS can write `.spn` with
the `.mmn`, `.amn` and `.eig` outputs disabled.

The final runs deliberately trigger I/O errors. `Verified 'Cannot open ...'`
means the driver found the expected diagnostic in the per-rank logs after
placing a directory at the output path. On systems with `/dev/full`, the
driver also verifies the write-error diagnostic.

For the ABACUS and Wannier90 settings to export spin matrices, see the
[Wannier90 interface guide](../../../docs/advanced/interface/Wannier90.md).
