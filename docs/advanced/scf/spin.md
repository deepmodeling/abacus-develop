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

3. DeltaSpin
The `DeltaSpin` function as proposed by Zefeng Cai and Ben Xu, et al. [arXiv:2208.04551v6](https://arxiv.org/abs/2208.04551) has been implemented in ABACUS in the LCAO basis for the `nspin 2` case. In order to use this function, the following parameters are needed to be set in the input file, for example:
```
#deltaspin
sc_mag_switch                   1
decay_grad_switch               0
sc_thr                          1e-7
nsc                             150
nsc_min                         2
sc_file                         sc.json
alpha_trial                     0.01
sccut                           3
```
The explanation of each input paramters has been explained in the [Noncollinear Spin Polarized Calculations](#noncollinear-spin-polarized-calculations) section.

An example of the sc_file json file is shown below:
```json
[
    {
        "element": "Fe",
        "itype": 0,
        "ScDecayGrad": 0.9,
        "ScAtomData": [
            {
                "index": 0,
                "lambda": 0.0,
                "target_mag": 2.0,
                "constrain": 1
            },
            {
                "index": 1,
                "lambda": 0,
                "target_mag": 2.0,
                "constrain": 1
            }
        ]
    }
]
```
Please refer the [Noncollinear Spin Polarized Calculations](#noncollinear-spin-polarized-calculations) section for the explanation of each input paramters. The difference is that `lambda`, `target_mag`, and `constrain` are scalars instead of vectors. Simple examples are provided in the `abacus-develop/examples/spin_polarized` directory.

## Noncollinear Spin Polarized Calculations
The spin non-collinear polarization calculation corresponds to setting **"noncolin 1"**, in which case the coupling between spin up and spin down will be taken into account. 
In this case, nspin is automatically set to 4, which is usually not required to be specified manually.
The weight of each band will not change, but the number of occupied states will be double. 
If the nbands parameter is set manually, it is generally set to twice what it would be when nspin<4.

In general, non-collinear magnetic moment settings are often used in calculations considering [SOC effects](#soc-effects). When **"lspinorb 1"** in INPUT file, "nspin" is also automatically set to 4. 
Note: different settings for "noncolin" and "lspinorb" correspond to different calculations:
 - noncolin=0 lspinorb=0 nspin<4 : 
Non-collinear magnetic moments and SOC effects are not considered.
 - noncolin=0 lspinorb=0 nspin=4 : 
Actualy same as the above setting, but the calculation will be larger. So the setting is not recommended.
 - noncolin=1 lspinorb=0 : 
Non-collinear magnetic moments are considered but SOC effects are not considered
 - noncolin=0 lspinorb=1 : 
The SOC effect is considered but the magnetic moment is limited to the Z direction
 - noncolin=1 lspinorb=1 : 
The SOC effect and non-collinear magnetic moment are both calculated.

### Constraint Spin functionality for noncollinear spin polarized calculations

The `DeltaSpin` function as proposed by Zefeng Cai and Ben Xu, et al. [arXiv:2208.04551v6](https://arxiv.org/abs/2208.04551) has been implemented in ABACUS in the LCAO basis. In order to use this function, the following parameters are needed to be set in the input file, for example:
```
#deltaspin
sc_mag_switch                   1
decay_grad_switch               1
sc_thr                          1e-7
nsc                             150
nsc_min                         2
sc_file                         sc.json
alpha_trial                     0.01
sccut                           3
```


`sc_mag_switch` is the switch of deltaspin functionality; `decay_grad_switch` is the switch of decay gradient method; `sc_thr` is the threshold of the spin constraint atomic magnetic moment in unit of Bohr Mag (\muB); `nsc` is the number of self-consistent iterations; `nsc_min` is the minimum number of self-consistent iterations; `sc_file` is the file name of the spin constraint parameters; `alpha_trial` is the initial trial step size for lambda in eV/uB^2; `sccut` restriction of step size in eV/uB.

An example of the sc_file json file is shown below:
```json
[
    {
        "element": "Fe",
        "itype": 0,
        "ScDecayGrad": 0.9,
        "ScAtomData": [
            {
                "index": 0,
                "lambda": [0, 0, 0],
                "target_mag": [2.0, 0.0, 0.0],
                "constrain": [1,1,1]
            },
            {
                "index": 1,
                "lambda": [0, 0, 0],
                "target_mag_val": 2.0,
                "target_mag_angle1": 80.0,
                "target_mag_angle2": 0.0,
                "constrain": [1,1,1]
            }
        ]
    }
]
```

The sc_file json file is a list of elemental data in total. For each element, the user should specify its name, the `itype` parameter should be in accord with `STRU` file and start from 0. `ScDecayGrad` is a parameter for each element in unit of  (uB^2/eV), this parameter needs to be determined for different element, for example, 0.9 uB^2/eV is an appropriate value for BCC-Fe according to Zefeng Cai's tests. `ScAtomData` specifies spin constraining parameters for each atom, the `index` starts from 0 and corresponds atomic order in the `STRU` file. `lambda` is a 3d vector for each atom, and it is recommended to set to [0.0, 0.0, 0.0] for all atoms. Users have two optional choices to set the target magnetic moments for each atom, i.e., by a 3d vector or by angles. If the `target_mag` is set, the `target_mag_val` and `target_mag_angle1` and `target_mag_angle2` will be ignored. The `target_mag` is a 3d vector in unit of Bohr Mag (\muB), and the `target_mag_val` is a scalar value in unit of Bohr Mag (\muB), `target_mag_angle1` and `target_mag_angle2` are two angles in unit of degree. The `constrain` is a 3d vector, if the corresponding element is set to 1, the corresponding component of the magnetic moment will be constrained, otherwise, it will be free. Note that the initial atomic magnetic moments are still set in the `STRU` file. Simple examples are provided in the `abacus-develop/examples/noncollinear` directory. One should set `noncolliear` to 1 to run the DeltaSpin function, `lspinorb=1` is not mandatory, but it is recommended to set to 1 to get more accurate results.


## For the continuation job
- Continuation job for "nspin 1" need file "SPIN1_CHG.cube" which is generated by setting "out_chg=1" in task before. By setting "init_chg file" in new job's INPUT file, charge density will start from file but not atomic. 
- Continuation job for "nspin 2" need files "SPIN1_CHG.cube" and "SPIN2_CHG.cube" which are generated by "out_chg 1" with "nspin 2", and refer to spin-up and spin-down charge densities respectively. It should be note that reading "SPIN1_CHG.cube" only for the continuation target magnetic moment job is not supported now.
- Continuation job for "nspin 4" need files "SPIN%s_CHG.cube", where %s in {1,2,3,4}, which are generated by "out_chg 1" with any variable setting leading to 'nspin'=4, and refer to charge densities in Pauli spin matrixes. It should be note that reading charge density files printing by 'nspin'=2 case is supported, which means only $\sigma_{tot}$ and $\sigma_{z}$ are read.

# SOC Effects
## SOC 
`lspinorb` is used for control whether or not SOC(spin-orbit coupling) effects should be considered.

Both `basis_type=pw` and `basis_type=lcao` support `scf` and `nscf` calculation with SOC effects.

Atomic forces and cell stresses can not be calculated with SOC effects yet. 

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

