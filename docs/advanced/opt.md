# Geometry Optimization

By setting `calculation` to be `relax` or `cell-relax`, ABACUS supports structural relaxation and variable-cell relaxation.

ABACUS provides two CG implementations for variable-cell relaxation, selected by the [relax_method](./input_files/input-main.md#relax_method) parameter:

- **CG variant 2** (`relax_method = cg 2` or `relax_method = cg`, default since v3.8): Uses simultaneous conjugate gradient (CG) optimization for both ionic positions and cell parameters. Both degrees of freedom are optimized together in each step.

- **CG variant 1** (`relax_method = cg 1`): Follows a nested procedure where fixed-cell structural relaxation is performed first, followed by an update of the cell parameters, and the process is repeated until convergence is achieved.

An example of the variable cell relaxation can be found in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/relax/pw_al), which is provided with the reference output file log.ref. When using CG variant 1, each ionic step is labelled in the following manner:
```
 -------------------------------------------
 RELAX CELL : 3
 RELAX IONS : 1 (in total: 15)
 -------------------------------------------
```

indicating that this is the first ionic step of the 3rd cell configuration, and it is the 15-th ionic step in total.


## Optimization Algorithms

ABACUS offers multiple optimization algorithms for structural relaxation, which can be selected using the [relax_method](./input_files/input-main.md#relax_method) keyword. The optional second value selects an implementation variant for CG or BFGS. For both methods, variant 1 is the traditional implementation and variant 2 is the recommended default.

### Available Algorithms

- **CG (Conjugate Gradient)**: Variant 1 optimizes ionic positions and cell parameters in separate stages; variant 2 optimizes them simultaneously with line search.
- **BFGS**: Quasi-Newton method for ionic relaxation.
- **LBFGS**: Limited-memory BFGS for ionic relaxation.
- **SD (Steepest Descent)**: Simple gradient descent for ionic relaxation.
- **CG-BFGS**: Mixed method that starts with CG and switches to BFGS when force convergence reaches the threshold set by [relax_cg_thr](./input_files/input-main.md#relax_cg_thr).

We also provide a [list of keywords](./input_files/input-main.md#geometry-relaxation) for controlling the relaxation process.

### BFGS method

The [BFGS method](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) is a quasi-Newton method for solving nonlinear optimization problems. It belongs to the class of quasi-Newton methods where the Hessian matrix is approximated during the optimization process. If the initial point is not far from the extrema, BFGS tends to work better than gradient-based methods.

ABACUS provides two BFGS implementations, controlled by the second element of [relax_method](./input_files/input-main.md#relax_method):

- **Default BFGS** (`relax_method = bfgs 2` or `relax_method = bfgs`): Updates the inverse of the approximate Hessian matrix B directly. This is the recommended implementation.

- **Traditional BFGS** (`relax_method = bfgs 1`): Updates the approximate Hessian matrix B itself, then obtains the inverse by solving matrix eigenvalues and taking their reciprocals. Both methods are mathematically equivalent, but in some cases the traditional variant may perform better.

### LBFGS method

The [L-BFGS (Limited-memory BFGS)](https://en.wikipedia.org/wiki/Limited-memory_BFGS) method is a memory-efficient variant of BFGS that stores only a few vectors representing the Hessian approximation instead of the full matrix. This makes it particularly suitable for large systems with many atoms.

Set `relax_method = lbfgs` to use this method.

### SD method

The [SD (steepest descent) method](https://en.wikipedia.org/wiki/Gradient_descent) is one of the simplest first-order optimization methods, where in each step the motion is along the direction of the gradient, where the function descends the fastest.

In practice, the SD method may take many iterations to converge, and is generally not recommended for production calculations.

### CG method

The [CG (conjugate gradient) method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) is one of the most widely used methods for solving optimization problems.

ABACUS provides two implementations of the CG method:

- **CG variant 2** (`relax_method = cg 2` or `relax_method = cg`, default): Performs simultaneous optimization of both ionic positions and cell parameters using a line search algorithm. This implementation is more efficient for `cell-relax` calculations as it optimizes all degrees of freedom together. The step size can be controlled by [relax_scale_force](./input_files/input-main.md#relax_scale_force).

- **CG variant 1** (`relax_method = cg 1`): Uses a nested procedure where ionic positions are optimized first using CG, followed by cell parameter optimization (also using CG) in `cell-relax` calculations. This is the traditional approach where the two optimization steps are separated.

The former `relax_new` parameter has been removed. Replace `relax_new = True` with `relax_method = cg 2`, and replace `relax_new = False` with `relax_method = cg 1` when using CG. The `bfgs`, `lbfgs`, `sd`, and `cg_bfgs` methods do not require a separate switch.

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

- **CG variant 2** (`relax_method = cg 2` or `relax_method = cg`):
  - `fixed_axes = "shape"`: Only allows volume changes (hydrostatic pressure), cell shape is fixed
  - `fixed_axes = "volume"`: Allows shape changes but keeps volume constant
  - `fixed_axes = "a"`, `"b"`, `"c"`, `"ab"`, `"ac"`, `"bc"`, or `"abc"`: Fix specific lattice vectors or combinations
  - `fixed_ibrav = True`: Maintain the Bravais lattice type during relaxation

- **`relax_method = cg 1`, `bfgs`, `lbfgs`, `sd`, or `cg_bfgs`**:
  - `fixed_axes = "a"`, `"b"`, `"c"`, `"ab"`, `"ac"`, `"bc"`, or `"abc"`: Fix specific lattice vectors or combinations
  - `fixed_axes = "shape"`, `"volume"` and `fixed_ibrav = True` are not available

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
