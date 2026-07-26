# Spin-polarization and SOC

## Non-spin-polarized Calculations
Setting of **"nspin 1"** in INPUT file means calculation with non-polarized spin. 
In this case, electrons with spin up and spin down have same occupations at every energy states, weights of bands per k point would be double. 

## Collinear Spin Polarized Calculations
Setting of **"nspin 2"** in INPUT file means calculation with polarized spin along z-axis. 
In this case, electrons with spin up and spin down will be calculated respectively, number of k points would be doubled.
Potential of electron and charge density will separate to spin-up case and spin-down case.

Magnetic moment Settings in [STRU files](../input_files/stru.md) are not ignored until **"nspin 2"** is set in INPUT file

When **"nspin 2"** is set, the screen output file will contain magnetic moment information. e.g.
```
 ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
 GE1    4.16e+00  4.36e+00  -6.440173e+03  0.000000e+00   6.516e-02  1.973e+01
```
where "TMAG" refers to total magnetization and "AMAG" refers to average magnetization.
For more detailed orbital magnetic moment information, please use [Mulliken charge analysis](../elec_properties/Mulliken.md).

### Constraint DFT for collinear spin polarized calculations
For some special need, there are two method to constrain electronic spin.

1. **"ocp"** and **"ocp_set"**
If **"ocp=1"** and **"ocp_set"** is set in INPUT file, the occupations of states would be fixed by **"ocp_set"**, this method is often used for excited states calculation. Be careful that: when **"nspin=1"**, spin-up and spin-down electrons will both be set, and when **"nspin=2"**, you should set spin-up and spin-down respectively.

2. **"nupdown"**
If **"nupdown"** is set to non-zero, number of spin-up and spin-down electrons will be fixed, and Fermi energy level will split to E_Fermi_up and E_Fermi_down. By the way, total magnetization will also be fixed, and will be the value of **"nupdown"**.

## DeltaSpin (Spin-Constrained DFT)

DeltaSpin is a spin-constrained DFT method that allows users to constrain the magnetic moments on individual atoms to target values during self-consistent field (SCF) calculations. This is useful for studying magnetic excitations, non-collinear magnetic structures, and systems where the magnetic ground state is not known a priori.

The theoretical foundation and implementation details can be found in:

- Cai Z, Wang K, Xu Y, et al., "A self-adaptive first-principles approach for magnetic excited states," *Quantum Frontiers* 2.1 (2023): 21. [DOI: 10.1007/s44214-023-00050-z](https://doi.org/10.1007/s44214-023-00050-z)
- Zheng D, Peng X, Huang Y, et al., "Integrating deep-learning-based magnetic model and non-collinear spin-constrained method: methodology, implementation and application," *npj Computational Materials* (2026).

### Enabling DeltaSpin

Set `sc_mag_switch 1` in the INPUT file. DeltaSpin is supported for both PW (`basis_type = pw`) and LCAO (`basis_type = lcao`) basis sets, with `nspin = 2` (collinear) or `nspin = 4` (non-collinear).

### Specifying Target Magnetic Moments in STRU

Target magnetic moments and constraint flags are specified per atom in the `ATOMIC_POSITIONS` section of the STRU file, using the `mag` (or `magmom`), `sc`, `lambda`, `angle1`, and `angle2` keywords after the atomic coordinates.

#### Collinear (nspin=2)

For collinear spin, only the z-component of the magnetic moment is constrained:

```
ATOMIC_POSITIONS
Direct

Fe
0.0
2
0.00  0.00  0.00  mag  2.0   sc 1
0.51  0.51  0.51  mag  -2.0  sc 1
```

- `mag 2.0`: target magnetic moment of 2.0 $\mu_B$ along z-axis
- `sc 1`: constrain the z-component (1 = constrained, 0 = unconstrained)

#### Non-collinear (nspin=4), vector form

For non-collinear spin, specify the magnetic moment as a vector (mx, my, mz):

```
ATOMIC_POSITIONS
Direct

Fe
0.0
2
0.00  0.00  0.00  mag  2.0  0.0  0.0  sc 1 1 1
0.51  0.51  0.51  mag  0.0  0.0  -2.0  sc 1 1 1
```

- `mag 2.0 0.0 0.0`: target moment vector in Cartesian coordinates ($\mu_B$)
- `sc 1 1 1`: constrain x, y, z components respectively

#### Non-collinear (nspin=4), angle form

Alternatively, use `angle1` (polar angle $\theta$) and `angle2` (azimuthal angle $\phi$) in degrees to specify the direction:

```
0.00  0.00  0.00  mag 2.0  angle1 0  angle2 0    sc 1 1 1
0.51  0.51  0.51  mag 2.0  angle1 180  angle2 0  sc 1 1 1
```

The Cartesian components are computed as:
- $m_z = |\mathbf{m}| \cos\theta$
- $m_x = |\mathbf{m}| \sin\theta \cos\phi$
- $m_y = |\mathbf{m}| \sin\theta \sin\phi$

#### Providing initial Lagrange multipliers

Initial lambda values (in eV/$\mu_B$) can be provided via the `lambda` keyword to accelerate convergence:

```
0.00  0.00  0.00  mag 2.0  lambda 0.01 0.0 0.0  sc 1 1 1
```

A single value sets $\lambda_z$; three values set $\lambda_x$, $\lambda_y$, $\lambda_z$.

#### Partial constraints

Set `sc 0` for unconstrained components. For example, to constrain only the direction but not the magnitude (use with `sc_direction_only`):

```
0.00  0.00  0.00  mag 2.0  0.0  0.0  sc 1 1 0
```

### DeltaSpin INPUT Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sc_mag_switch` | Boolean | False | Enable DeltaSpin |
| `sc_thr` | Real | 1.0e-6 | Convergence criterion for lambda loop (RMS, in $\mu_B$) |
| `nsc` | Integer | 100 | Maximum number of lambda iterations |
| `nsc_min` | Integer | 2 | Minimum number of lambda iterations |
| `sc_scf_nmin` | Integer | 2 | Minimum outer SCF iterations before starting lambda loop |
| `alpha_trial` | Real | 0.01 | Initial trial step size for lambda (eV/$\mu_B^2$) |
| `sccut` | Real | 3.0 | Maximum step size for lambda (eV/$\mu_B$) |
| `sc_drop_thr` | Real | 1.0e-2 | Convergence ratio threshold for adaptive lambda loop |
| `sc_scf_thr` | Real | 1.0e-4 | Density error threshold for entering lambda loop |
| `sc_direction_only` | Boolean | False | Constrain only the direction, not the magnitude |
| `sc_lambda_strategy` | String | bfgs | Lambda update strategy (see below) |
| `decay_grad_switch` | Boolean | False | Enable gradient-based early exit |

For full parameter details, see the [Spin-Constrained DFT](../input_files/input-main.md#spin-constrained-dft) section of the input keyword list.

### Lambda Update Strategies

The `sc_lambda_strategy` parameter controls how the Lagrange multipliers $\lambda$ are updated during the lambda loop:

- **`bfgs`** (default): BFGS quasi-Newton method with line search. Robust and well-tested for both PW and LCAO. Uses `alpha_trial` and `sccut` to control step size.

- **`linear_response`**: Linear response method (Scheme B). Estimates the magnetic susceptibility $\chi$ from the history of $(\lambda, M)$ pairs and performs a one-step Newton-like update: $\Delta\lambda = \beta (M_{\text{target}} - M) / \chi$, where $\beta$ is a mixing parameter.

- **`augmented_lagrangian`**: Augmented Lagrangian method (Scheme C). Uses a penalty parameter $\mu$ that grows over iterations: $\lambda_{\text{new}} = \lambda + \mu (M - M_{\text{target}})$. The penalty increases until convergence is achieved.

- **`hybrid_delayed`**: Hybrid delayed update (Scheme D). Two-phase approach: in the early phase (SCF not yet converged), lambda updates are gentle; in the late phase (SCF nearly converged), augmented Lagrangian updates are applied.

### Direction-Only Mode

When `sc_direction_only 1` is set, only the **direction** of the magnetic moment is constrained to match the target, while the magnitude is allowed to vary freely. This is useful for:

- Studying magnetic anisotropy energy surfaces
- Cases where the moment magnitude is determined by the electronic structure
- Converging to the easy-axis direction without fixing the moment size

In this mode, the lambda vector is projected to be perpendicular to the target moment direction at each iteration, ensuring it can only rotate the magnetization, not stretch it.

### Combining DeltaSpin with DFT+U

DeltaSpin can be combined with DFT+U for strongly correlated systems. When both `sc_mag_switch` and `dft_plus_u` are enabled:

1. DFT+U occupation update runs first in each SCF iteration
2. DeltaSpin lambda loop runs after, constraining the magnetic moments
3. The DFT+U-corrected Hamiltonian is used by the lambda loop

Example INPUT for PW DFT+U + DeltaSpin:

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
nspin               2
dft_plus_u          1
orbital_corr        -1 2
hubbard_u           0.0 4.0
sc_mag_switch       1
sc_thr              1.0e-6
sc_scf_thr          1.0e-4
sc_lambda_strategy  bfgs
```

### Example: Collinear antiferromagnetic Fe

INPUT file:

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
nspin               2
sc_mag_switch       1
sc_thr              1.0e-6
```

STRU file:

```
ATOMIC_SPECIES
Fe 55.845 Fe.upf

LATTICE_CONSTANT
8.190

LATTICE_VECTORS
 1.00  0.50  0.50
 0.50  1.00  0.50
 0.50  0.50  1.00

ATOMIC_POSITIONS
Direct

Fe
0.0
2
0.00  0.00  0.00  mag  2.0  sc 1
0.51  0.51  0.51  mag  -2.0  sc 1
```

### Example: Non-collinear constrained moments

INPUT file:

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
nspin               4
noncolin            1
sc_mag_switch       1
sc_direction_only   1
sc_lambda_strategy  bfgs
```

STRU file:

```
ATOMIC_POSITIONS
Direct

Fe
0.0
2
0.00  0.00  0.00  mag  2.0  0.0  0.0  sc 1 1 0
0.51  0.51  0.51  mag  0.0  0.0  2.0  sc 1 1 0
```

## Noncollinear Spin Polarized Calculations
The spin non-collinear polarization calculation corresponds to setting **"noncolin 1"**, in which case the coupling between spin up and spin down will be taken into account.
In this case, nspin is automatically set to 4, which is usually not required to be specified manually.
The weight of each band will not change, but the number of occupied states will be double.
If the nbands parameter is set manually, it is generally set to twice what it would be when nspin<4.

In general, non-collinear magnetic moment settings are often used in calculations considering [SOC effects](#soc-effects). When **"lspinorb 1"** in INPUT file, "nspin" is also automatically set to 4.

Note: different settings for "noncolin" and "lspinorb" correspond to different calculations:

| noncolin | lspinorb | nspin | Effect | When to Use |
|----------|----------|-------|--------|-------------|
| 0 | 0 | <4 | No non-collinear magnetism, no SOC | Standard collinear spin-polarized or non-spin-polarized calculations |
| 0 | 0 | 4 | Same as above, but larger calculation | **Not recommended** - wastes computational resources |
| 1 | 0 | 4 | Non-collinear magnetism WITHOUT SOC | Systems with complex magnetic structures (e.g., spin spirals, frustrated magnets) where SOC is negligible |
| 0 | 1 | 4 | SOC WITH z-axis magnetism only | Non-magnetic materials with SOC (e.g., semiconductors with band splitting), or magnetic materials where magnetism is along z-axis |
| 1 | 1 | 4 | Both SOC AND non-collinear magnetism | Heavy-element magnetic materials where both SOC and non-collinear magnetism are important (e.g., magnetic anisotropy, Dzyaloshinskii-Moriya interaction) |

**Special case**: `noncolin=0, lspinorb=1` is commonly used for non-magnetic materials with SOC effects (e.g., topological insulators, semiconductors with spin-orbit splitting). In this case, the magnetization is NOT automatically set, implying no magnetic moments in the system.

### Choosing gga_grad for Noncollinear Calculations

When performing noncollinear spin calculations (`nspin=4`), the `gga_grad` parameter controls how the gradient of the magnetization is computed in GGA exchange-correlation functionals. This is critical for obtaining correct magnetic anisotropy energies and spin-orbit coupling effects.

#### Available Options

| gga_grad | Gradient Method | Description | When to Use |
|----------|----------------|-------------|-------------|
| 0 | Standard | Conventional GGA gradient (collinear approximation) | Collinear calculations only (`nspin<4`) |
| 1 | $\nabla\|\mathbf{m}\|$ | Gradient of magnetization magnitude | Noncollinear with SOC, when only magnitude variation matters |
| 2 | $\nabla\hat{\mathbf{m}}$ | Gradient of magnetization direction (default for `nspin=4`) | **Recommended** for most noncollinear calculations with SOC |
| 3 | Built-in SF | Special functional for noncollinear SOC (requires `v_xc_ncgga_sf_builtin`) | Advanced noncollinear SOC calculations |

#### Recommendations

**For noncollinear calculations with SOC (`lspinorb=1`)**:
- **Use `gga_grad=2`** (default when `nspin=4`): This computes the gradient of the magnetization direction unit vector $\hat{\mathbf{m}}$, which is physically correct for spin-orbit coupling effects.
- This is automatically set when you specify `nspin=4` or `lspinorb=1`, but you can explicitly set it to be clear.

**For noncollinear calculations without SOC (`noncolin=1`, `lspinorb=0`)**:
- **Use `gga_grad=2`** for spin spiral calculations or frustrated magnetic systems
- The gradient of the magnetization direction is still important even without SOC

**For collinear calculations (`nspin=1` or `nspin=2`)**:
- **Use `gga_grad=0`** (default): Standard collinear GGA gradient
- Do NOT use `gga_grad=1`, `2`, or `3` for collinear calculations

#### Example INPUT for Noncollinear SOC

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
nspin               4
noncolin            1
lspinorb            1
gga_grad            2
```

#### Important Notes

1. **Automatic setting**: When `nspin=4`, ABACUS automatically sets `gga_grad=2` if not explicitly specified
2. **Physical correctness**: Using the wrong `gga_grad` value can lead to incorrect magnetic anisotropy energies and spin textures
3. **Compatibility**: `gga_grad=3` requires the built-in special functional and may not be available for all XC functionals
4. **Performance**: `gga_grad=2` has similar computational cost to `gga_grad=0` for most systems

## For the continuation job
- Continuation job for "nspin 1" need file "SPIN1_CHG.cube" which is generated by setting "out_chg=1" in task before. By setting "init_chg file" in new job's INPUT file, charge density will start from file but not atomic. 
- Continuation job for "nspin 2" need files "SPIN1_CHG.cube" and "SPIN2_CHG.cube" which are generated by "out_chg 1" with "nspin 2", and refer to spin-up and spin-down charge densities respectively. It should be note that reading "SPIN1_CHG.cube" only for the continuation target magnetic moment job is not supported now.
- Continuation job for "nspin 4" need files "SPIN%s_CHG.cube", where %s in {1,2,3,4}, which are generated by "out_chg 1" with any variable setting leading to 'nspin'=4, and refer to charge densities in Pauli spin matrixes. It should be note that reading charge density files printing by 'nspin'=2 case is supported, which means only $\sigma_{tot}$ and $\sigma_{z}$ are read.

# SOC Effects
## SOC
`lspinorb` is used for control whether or not SOC(spin-orbit coupling) effects should be considered.

Both `basis_type=pw` and `basis_type=lcao` support `scf` and `nscf` calculation with SOC effects.

Atomic forces and cell stresses can be calculated with SOC effects (supported since ABACUS v3.9.0). 

## Pseudopotentials and Numerical Atomic Orbitals
For Norm-Conserving pseudopotentials, there are differences between SOC version and non-SOC version.

Please check your pseudopotential files before calculating.
In `PP_HEADER` part, keyword `has_so=1` and `relativistic="full"` refer to SOC effects have been considered, 
which would lead to different `PP_NONLOCAL` and `PP_PSWFC` parts.
Please be careful that `relativistic="full"` version can be used for SOC or non-SOC calculation, but `relativistic="scalar"` version only can be used for non-SOC calculation.
When full-relativistic pseudopotential is used for non-SOC calculation, ABACUS will automatically transform it to scalar-relativistic version.

Numerical atomic orbitals in ABACUS are unrelated with spin, and same orbital file can be used for SOC and non-SOC calculation.

## Partial-relativistic SOC Effect
Sometimes, for some real materials, both scalar-relativistic and full-relativistic can not describe the exact spin-orbit coupling.
Artificial modulation can help for these cases.

`soc_lambda`, which has value range [0.0, 1.0] , is used for modulate SOC effect.
In particular, `soc_lambda 0.0` refers to scalar-relativistic case and `soc_lambda 1.0` refers to full-relativistic case.

## Pseudopotential Requirements for SOC

When performing SOC calculations (`lspinorb=1`), specific pseudopotential requirements must be met:

### Checking Pseudopotential Files

In the UPF (Unified Pseudopotential Format) file header (`PP_HEADER` section), look for:
- `has_so="T"` or `has_so="1"`: Indicates SOC information is included
- `relativistic="full"`: Indicates full-relativistic pseudopotential

Example from a full-relativistic UPF file:
```
<PP_HEADER
   ...
   relativistic="full"
   has_so="T"
   ...
/>
```

### Pseudopotential Usage Rules

1. **For SOC calculations** (`lspinorb=1`):
   - **MUST** use full-relativistic pseudopotentials with `has_so=true`
   - Code will terminate with error: "no soc upf used for lspinorb calculation" if scalar-relativistic PP is used

2. **For non-SOC calculations** (`lspinorb=0`):
   - Can use either scalar-relativistic or full-relativistic pseudopotentials
   - If full-relativistic PP is used, ABACUS automatically transforms it to scalar-relativistic version

3. **For ultrasoft pseudopotentials (USPP)**:
   - Full-relativistic USPP **requires** `lspinorb=true`
   - Code will show warning: "FR-USPP please use lspinorb=.true." if this requirement is not met

### Where to Find SOC Pseudopotentials

- **SG15_ONCV**: Full-relativistic versions available at [quantum-simulation.org](http://quantum-simulation.org/potentials/sg15_oncv/upf/)
- **PseudoDOJO**: Provides both scalar and full-relativistic versions
- **ABACUS official**: [abacus.ustc.edu.cn](http://abacus.ustc.edu.cn/pseudo/list.htm)

## Automatic Parameter Settings

When using SOC or non-collinear calculations, ABACUS automatically adjusts several parameters:

### When `lspinorb=true`:
1. **nspin**: Automatically set to 4 (noncollinear spin representation)
2. **Symmetry**: Automatically disabled (`symm_flag=-1`) because SOC breaks inversion symmetry
3. **Magnetization**: NOT automatically set when `noncolin=0` (implies non-magnetic material with SOC)

### When `noncolin=true`:
1. **nspin**: Automatically set to 4
2. **npol**: Set to 2 (wave function has two spinor components)
3. **Magnetization**: Automatically set if user provides zero values (unless `lspinorb=1` and `noncolin=0`)

### Important Notes:
- You do NOT need to manually set `nspin=4` when using `lspinorb=1` or `noncolin=1`
- Symmetry operations are incompatible with SOC, so they are automatically turned off
- For `lspinorb=1, noncolin=0`: This is a special case for non-magnetic materials with SOC, where magnetization is not initialized

## Common Errors and Solutions

### Error: "no soc upf used for lspinorb calculation"
**Cause**: Using scalar-relativistic pseudopotentials with `lspinorb=1`

**Solution**: Download and use full-relativistic pseudopotentials with `has_so=true`. Check the UPF file header to verify `relativistic="full"` and `has_so="T"`.

### Error: "nspin=4(soc or noncollinear-spin) does not support gamma only calculation"
**Cause**: Trying to use `gamma_only=true` with `lspinorb=1` or `noncolin=1`

**Solution**: Set `gamma_only=false` or `gamma_only=0` in your INPUT file. SOC and non-collinear calculations require k-point sampling beyond the gamma point.

### Warning: "FR-USPP please use lspinorb=.true."
**Cause**: Using full-relativistic ultrasoft pseudopotentials without enabling SOC

**Solution**: Set `lspinorb=true` in your INPUT file, or switch to scalar-relativistic USPP if SOC is not needed.

### Issue: Forces or stresses not calculated
**Note**: This issue has been resolved. Atomic forces and cell stresses can now be calculated with SOC effects (supported since ABACUS v3.9.0).

If you are using an older version of ABACUS (before v3.9.0), force and stress calculations with SOC were not supported. Please upgrade to the latest version to use this feature.

## INPUT File Examples

### Example 1: SOC without Non-collinear Magnetism
For non-magnetic materials with SOC (e.g., GaAs, topological insulators):

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
lspinorb            1        # Enable SOC
noncolin            0        # No non-collinear magnetism
# nspin will be automatically set to 4
# symmetry will be automatically disabled
```

### Example 2: Non-collinear Magnetism without SOC
For systems with complex magnetic structures but negligible SOC:

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
lspinorb            0        # No SOC
noncolin            1        # Enable non-collinear magnetism
# nspin will be automatically set to 4
# Magnetization directions should be specified in STRU file
```

### Example 3: Both SOC and Non-collinear Magnetism
For heavy-element magnetic materials (e.g., Fe with SOC, materials with DMI):

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             60
lspinorb            1        # Enable SOC
noncolin            1        # Enable non-collinear magnetism
# nspin will be automatically set to 4
# symmetry will be automatically disabled
# Magnetization directions should be specified in STRU file
```

### Example 4: Partial-relativistic SOC
For fine-tuning SOC strength:

```
INPUT_PARAMETERS
calculation         scf
basis_type          pw
ecutwfc             50
lspinorb            1        # Enable SOC
soc_lambda          0.5      # 50% SOC strength
# Useful when full SOC overestimates or underestimates experimental results
```
