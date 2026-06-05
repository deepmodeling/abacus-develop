# The KPT file

- [Generate k-mesh automatically](#generate-k-mesh-automatically)
- [Set k-points explicitly](#set-k-points-explicitly)
- [Band structure calculations](#band-structure-calculations)

ABACUS uses periodic boundary conditions for both crystals and finite systems. For isolated systems, such as atoms, molecules, clusters, etc., one uses the so-called supercell model. Lattice vectors of the supercell are set in the `STRU` file. For the input k-point (`KPT`) file, the file should either contain the k-point coordinates and weights or the mesh size for creating the k-point gird. Both options are allowed in `ABACUS`.

## Gamma-only Calculations

In ABACUS, we offer the option of running gamma-only calculations for LCAO basis by setting [gamma_only](./input-main.md#gamma_only) to be 1. Due to details of implementation, gamma-only calculation will be slightly faster than running a non gamma-only calculation and explicitly setting gamma point to be the only the k-point, but the results should be consistent.

> If gamma_only is set to 1, the KPT file will be overwritten. So make sure to turn off gamma_only for multi-k calculations.

## Generate k-mesh automatically

To generate k-mesh automatically, it requires the input subdivisions of the Brillouin zone
in each direction and the origin for the k-mesh. ABACUS uses the Monkhorst-Pack
method to generate k-mesh, and the following is an example input k-point (`KPT`) file used in
`ABACUS`.

```
K_POINTS //keyword for start
0 //total number of k-point, `0' means generate automatically
Gamma //which kind of Monkhorst-Pack method, `Gamma' or `MP'
2 2 2 0 0 0 //first three number: subdivisions along reciprocal vectors
            //last three number: shift of the mesh
```

In the above example, the first line is a keyword, and it can be set as `K_POINTS`, or `KPOINTS` or just `K`.
The second line is an integer, and its value determines how to get k-points. In this example, `0` means using Monkhorst-Pack (MP) method to generate k-points automatically.

The third line tells the input type of k-points, `Gamma` or `MP`, different Monkhorst Pack
(MP) method. Monkhorst-Pack (MP) is a method which uses the uniform k-points sampling in
Brillouin-zone, while `Gamma` means the &Gamma;-centered Monkhorst-Pack method.
The first three numbers of the last line are integers, which give the MP k grid dimensions, and
the rest three are real numbers, which give the offset of the k grid. In this example, the numbers
`0 0 0` means that there is no offset, and this is the a standard 2by2by2 k grid.

[back to top](#the-kpt-file)

## Set k-points explicitly

If the user wants to set up the k-points explicitly, the input k-point file should contain
the k-point coordinates and weights. An example is given as follows:

```
K_POINTS //keyword for start
8 //total number of k-point
Direct //`Direct' or `Cartesian' coordinate
0.0 0.0 0.0 0.125 //coordinates and weights
0.5 0.0 0.0 0.125
0.0 0.5 0.0 0.125
0.5 0.5 0.0 0.125
0.0 0.0 0.5 0.125
0.5 0.0 0.5 0.125
0.0 0.5 0.5 0.125
0.5 0.5 0.5 0.125
```

### K-point Weights and Symmetry

When explicitly setting k-points, you can specify custom weights for each k-point. These weights determine the contribution of each k-point to the total energy and density calculations.

**Important notes about k-point weights:**

1. **Custom weights are preserved**: When using explicit k-point lists (non-Monkhorst-Pack), ABACUS preserves the custom weights you specify, even when symmetry operations are applied to reduce the k-points to the irreducible Brillouin zone (IBZ).

2. **Symmetry reduction**: When [`symmetry`](./input-main.md#symmetry) is set to 1, ABACUS will analyze the crystal symmetry and reduce the k-point set to the irreducible Brillouin zone. During this reduction:
   - For **Monkhorst-Pack grids** (automatically generated): All k-points have uniform weights (1/N where N is the total number of k-points)
   - For **explicit k-point lists**: Custom weights are preserved and properly combined when symmetry-equivalent k-points are merged

3. **Weight normalization**: After symmetry reduction, k-point weights are normalized so that their sum equals `degspin` (2 for non-spin-polarized calculations, 1 for spin-polarized calculations).

**Example with custom weights:**

```
K_POINTS
5
Direct
0.0 0.0 0.0   0.1   // Gamma point with weight 0.1
0.5 0.0 0.0   0.2   // X point with weight 0.2
0.0 0.5 0.0   0.3   // Y point with weight 0.3
0.5 0.5 0.0   0.2   // M point with weight 0.2
0.0 0.0 0.5   0.2   // Z point with weight 0.2
```

In this example, different k-points have different weights, which might be useful for:
- Special sampling schemes
- Convergence testing with specific k-point importance
- Custom integration methods

> **Note**: When using custom weights with symmetry, ensure that your weight distribution is consistent with the crystal symmetry. ABACUS will preserve your weights during IBZ reduction, but inconsistent weights may lead to unexpected results.

[back to top](#the-kpt-file)

## Band structure calculations

ABACUS uses specified high-symmetry directions of the Brillouin zone for band structure
calculations. The third line of k-point file should start with 'Line' or 'Line_Cartesian' for
line mode. 'Line' means the positions below are in Direct coordinates, while 'Line_Cartesian'
means in Cartesian coordinates:

```
K_POINTS // keyword for start
6 // number of high symmetry lines
Line // line-mode
0.5 0.0 0.5 20 // X
0.0 0.0 0.0 20 // G
0.5 0.5 0.5 20 // L
0.5 0.25 0.75 20 // W
0.375 0.375 0.75 20 // K
0.0 0.0 0.0 1 // G
```

The fourth line and the following are special k-point coordinates and number of k-points
between this special k-point and the next.

[back to top](#the-kpt-file)
