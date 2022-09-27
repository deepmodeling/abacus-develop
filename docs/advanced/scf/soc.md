## SOC Effects
### SOC 
`lspinorb` is used for control whether or not SOC(spin-orbit coupling) effects should be considered.

Both `basis_type=pw` and `basis_type=lcao` support `scf` and `nscf` calculation with SOC effects.

Atomic forces and cell stresses can not be calculated with SOC effects yet. 

### Pseudopotentials and Numerical Atomic Orbitals
For Norm-Conserving pseudopotentials, there are differences between SOC version and non-SOC version.

Please check your pseudopotential files before calculating.
In `PP_HEADER` part, keyword `has_so=1` and `relativistic="full"` refer to SOC effects have been considered, 
which would lead to different `PP_NONLOCAL` and `PP_PSWFC` parts.
Please be careful that `relativistic="full"` version can be used for SOC or non-SOC calculation, but `relativistic="scalar"` version only can be used for non-SOC calculation.
When full-relativistic pseudopotential is used for non-SOC calculation, ABACUS will automatically transform it to scalar-relativistic version.

Numerical atomic orbitals in ABACUS are unrelated with spin, and same orbital file can be used for SOC and non-SOC calculation.

### Partial-relativistic SOC Effect
Sometimes, for some real materials, both scalar-relativistic and full-relativistic can not describe the exact spin-orbit coupling. 
Artificial modulation can help for these cases.

`soc_lambda`, which has value range [0.0, 1.0] , is used for modulate SOC effect.
In particular, `soc_lambda 0.0` refers to scalar-relativistic case and `soc_lambda 1.0` refers to full-relativistic case.
