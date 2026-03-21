# Geometry Optimization

By setting `calculation` to be `relax` or `cell-relax`, ABACUS supports structural relaxation and variable-cell relaxation.

ABACUS provides two implementations for variable-cell relaxation, controlled by the [relax_new](./input_files/input-main.md#relax_new) parameter:

- **New implementation** (`relax_new = True`, default since v3.8): Uses a simultaneous conjugate gradient (CG) optimization for both ionic positions and cell parameters. Both degrees of freedom are optimized together in each step.

- **Old implementation** (`relax_new = False`): Follows a nested procedure where fixed-cell structural relaxation is performed first, followed by an update of the cell parameters, and the process is repeated until convergence is achieved.

An example of the variable cell relaxation can be found in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/relax/pw_al), which is provided with the reference output file log.ref. When using the old implementation (`relax_new = False`), each ionic step is labelled in the following manner:
```
 -------------------------------------------
 RELAX CELL : 3
 RELAX IONS : 1 (in total: 15)
 -------------------------------------------
```

indicating that this is the first ionic step of the 3rd cell configuration, and it is the 15-th ionic step in total.


## Optimization Algorithms

ABACUS offers multiple optimization algorithms for structural relaxation, which can be selected using the [relax_method](./input_files/input-main.md#relax_method) keyword. The available algorithms and their behavior depend on the [relax_new](./input_files/input-main.md#relax_new) setting:

### Algorithm Availability

**New implementation** (`relax_new = True`, default):
- **CG (Conjugate Gradient)**: Simultaneous optimization of both ionic positions and cell parameters using CG with line search. This is the only algorithm available for the new implementation.

**Old implementation** (`relax_new = False`):
- **CG (Conjugate Gradient)**: For ionic relaxation; CG is also used for cell parameter optimization in the nested procedure
- **BFGS**: Quasi-Newton method for ionic relaxation
- **LBFGS**: Limited-memory BFGS for ionic relaxation
- **SD (Steepest Descent)**: Simple gradient descent for ionic relaxation
- **CG-BFGS**: Mixed method that starts with CG and switches to BFGS when force convergence reaches the threshold set by [relax_cg_thr](./input_files/input-main.md#relax_cg_thr)

We also provide a [list of keywords](./input_files/input-main.md#geometry-relaxation) for controlling the relaxation process.

### BFGS method

The [BFGS method](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) is a quasi-Newton method for solving nonlinear optimization problems. It belongs to the class of quasi-Newton methods where the Hessian matrix is approximated during the optimization process. If the initial point is not far from the extrema, BFGS tends to work better than gradient-based methods.

**Note**: BFGS is only available with the old implementation (`relax_new = False`).

ABACUS provides two BFGS implementations, controlled by the second element of [relax_method](./input_files/input-main.md#relax_method):

- **Default BFGS** (`relax_method = bfgs 2` or `relax_method = bfgs`): Updates the inverse of the approximate Hessian matrix B directly. This is the recommended implementation.

- **Traditional BFGS** (`relax_method = bfgs 1`): Updates the approximate Hessian matrix B itself, then obtains the inverse by solving matrix eigenvalues and taking their reciprocals. Both methods are mathematically equivalent, but in some cases the traditional variant may perform better.

### LBFGS method

The [L-BFGS (Limited-memory BFGS)](https://en.wikipedia.org/wiki/Limited-memory_BFGS) method is a memory-efficient variant of BFGS that stores only a few vectors representing the Hessian approximation instead of the full matrix. This makes it particularly suitable for large systems with many atoms.

**Note**: LBFGS is only available with the old implementation (`relax_new = False`). Set `relax_method = lbfgs` to use this method.

### SD method

The [SD (steepest descent) method](https://en.wikipedia.org/wiki/Gradient_descent) is one of the simplest first-order optimization methods, where in each step the motion is along the direction of the gradient, where the function descends the fastest.

**Note**: SD is only available with the old implementation (`relax_new = False`).

In practice, the SD method may take many iterations to converge, and is generally not recommended for production calculations.

### CG method

The [CG (conjugate gradient) method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) is one of the most widely used methods for solving optimization problems.

ABACUS provides two implementations of the CG method:

- **New CG implementation** (`relax_new = True`, default): Performs simultaneous optimization of both ionic positions and cell parameters using a line search algorithm. This implementation is more efficient for `cell-relax` calculations as it optimizes all degrees of freedom together. The step size can be controlled by [relax_scale_force](./input_files/input-main.md#relax_scale_force).

- **Old CG implementation** (`relax_new = False`): Uses a nested procedure where ionic positions are optimized first using CG, followed by cell parameter optimization (also using CG) in `cell-relax` calculations. This is the traditional approach where the two optimization steps are separated.

## Constrained Optimization

Apart from conventional optimization where all degrees of freedom are allowed to move, we also offer the option of doing constrained optimization in ABACUS.

### Fixing Atomic Positions  
Users may note that in the above-mentioned example, the atomic positions in STRU file are given along with three integers:

```
Al
0.0
4
0.00 0.00 0.00 1 1 1
0.53 0.50 0.00 1 1 1
0.50 0.00 0.52 1 1 1
0.00 0.50 0.50 1 1 1
```

For relaxation calculations, the three integers denote whether the corresponding degree of freedom is allowed to move. For example, if we replace the STRU file by:
```
Al
0.0
4
0.00 0.00 0.00 1 1 0
0.53 0.50 0.00 1 1 1
0.50 0.00 0.52 1 1 1
0.00 0.50 0.50 1 1 1
```

then the first Al atom will not be allowed to move in z direction.

Fixing atomic position is sometimes helpful during relaxation of isolated molecule/cluster, to prevent the system from drifting in space.

### Fixing Cell Parameters

Sometimes we want to do variable-cell relaxation with some of the cell degrees of freedom fixed. This is achieved by keywords such as [fixed_axes](./input_files/input-main.md#fixed_axes), [fixed_ibrav](./input_files/input-main.md#fixed_ibrav) and [fixed_atoms](./input_files/input-main.md#fixed_atoms).

**Available constraints by implementation:**

- **New implementation** (`relax_new = True`):
  - `fixed_axes = "shape"`: Only allows volume changes (hydrostatic pressure), cell shape is fixed
  - `fixed_axes = "volume"`: Allows shape changes but keeps volume constant
  - `fixed_axes = "a"`, `"b"`, `"c"`, etc.: Fix specific lattice vectors or combinations
  - `fixed_ibrav = True`: Maintain the Bravais lattice type during relaxation

- **Old implementation** (`relax_new = False`):
  - **All `fixed_axes` options now supported**: "shape", "volume", "a", "b", "c", "ab", "ac", "bc", "abc"
  - **`fixed_ibrav` now supported**: Maintains Bravais lattice type during relaxation
  - Can combine `fixed_axes` with `fixed_ibrav` for constrained relaxation
  - **Implementation approach**: Uses post-update constraint enforcement (volume rescaling and lattice reconstruction after each CG step)

**VASP ISIF correspondence:**

If you are familiar with the `ISIF` option from VASP, here is the correspondence:

- ISIF = 0 : calculation = "relax"
- ISIF = 1, 2 : calculation = "relax", cal_stress = 1
- ISIF = 3 : calculation = "cell-relax"
- ISIF = 4 : calculation = "cell-relax", fixed_axes = "volume"
- ISIF = 5 : calculation = "cell-relax", fixed_axes = "volume", fixed_atoms = True
- ISIF = 6 : calculation = "cell-relax", fixed_atoms = True
- ISIF = 7 : calculation = "cell-relax", fixed_axes = "shape", fixed_atoms = True

### Stop Geometry Optimization Manually

It is usually difficult to converge when calculating large systems, but people do not want to give up this calculation result.
Providing a file named `EXIT`:
```
stop_ion    true
```
ABACUS will end normally and produce a complete file.