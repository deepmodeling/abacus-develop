# Geometry Optimization

By setting `calculation` to be `relax` or `cell-relax`, ABACUS supports structural relaxation and variable-cell relaxation.

Current implementation of variable-cell relaxation in ABACUS now follows a nested procedure: fixed cell structural relaxation will be performed, followed by an update of the cell parameters, and the process is repeated until convergence is achieved.

An example of the variable cell relaxation can be found in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/relax/pw_al), which is provided with the reference output file log.ref. Note that in log.ref, each ionic step is labelled in the following manner:
```
 -------------------------------------------
 RELAX CELL : 3
 RELAX IONS : 1 (in total: 15)
 -------------------------------------------
``` 

indicating that this is the first ionic step of the 3rd cell configuration, and it is the 15-th ionic step in total.


## Optimization Algorithms
In the nested procedure mentioned above, we used CG method to perform cell relaxation, while offering four different algorithms for doing fixed-cell structural relaxation: BFGS, SD(steepest descent), CG(conjugate gradient), as well as a mixed CG-BFGS method. The optimziation algorithm can be selected using keyword [relax_method](./input_files/input-main.md#relax_method). We also provide a [list of keywords](./input_files/input-main.md#geometry-relaxation) for controlling the relaxation process.

### BFGS method

The [BFGS method](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) is a quasi-Newton method for solving nonlinear optimization problem. It belongs to the class of quasi-Newton method where the Hessian matrix is approximated during the optimization process. If the initial point is not far from the extrema, BFGS tends to work better than gradient-based methods.

In ABACUS, we implemented the BFGS method for doing fixed-cell structural relaxation.

### SD method

The [SD (steepest descent) method](https://en.wikipedia.org/wiki/Gradient_descent) is one of the simplest first-order optimization methods, where in each step the motion is along the direction of the gradient, where the function descents the fastest.

In practice, SD method may take many iterations to converge, and is generally not used.

### CG method

The [CG (conjugate gradient) method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) is one of the most widely used methods for solving optimization problems.

In ABACUS, we implemented the CG method for doing fixed-cell structural relaxation as well as the optimization of cell parameters.

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
Sometimes we want to do variable-cell relaxation with some of the cell degrees of freedom fixed. This is achieved by keywords such as [fixed_axes](./input_files/input-main.md#fixed_axes), [fixed_ibrav](./input_files/input-main.md#fixed_ibrav) and [fixed_atoms](./input_files/input-main.md#fixed_atoms). Specifically, if users are familiar with the `ISIF` option from VASP, then we offer the following correspondence:

- ISIF = 0 : calculation = "relax"
- ISIF = 1, 2 : calculation = "relax", cal_stress = 1
- ISIF = 3 : calculation = "cell-relax"
- ISIF = 4 : calculation = "cell-relax", fixed_axes = "volume"
- ISIF = 5 : calculation = "cell-relax", fixed_axes = "volume", fixed_atoms = True
- ISIF = 6 : calculation = "cell-relax", fixed_atoms = True
- ISIF = 7 : calculation = "cell-realx", fixed_axes = "shape", fixed_atoms = True
