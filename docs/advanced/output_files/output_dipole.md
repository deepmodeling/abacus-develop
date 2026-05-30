# Outputting Dipole Moment

ABACUS can output the dipole moment by adding the keyword `out_dipole` in the INPUT file:

```
out_dipole             1
```

## Supported Calculations

This feature is available for all types of DFT calculations that use charge density:

- **KSDFT** (Kohn-Sham DFT)
  - Plane wave (PW) basis
  - Linear combination of atomic orbitals (LCAO) basis
- **SDFT** (Stochastic DFT)
- **OFDFT** (Orbital-Free DFT)
- **TDDFT** (Time-Dependent DFT)

## Input Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `out_dipole` | Integer | Whether to output dipole moment | 0 |

- `out_dipole = 0`: Disable dipole output
- `out_dipole = 1`: Enable dipole output

## Output Files

When `out_dipole` is set to 1, ABACUS will generate files named `dipole_s${spin}.txt` for each spin channel in the output directory.

For spin-polarized calculation (nspin=2):
- `dipole_s1.txt`
- `dipole_s2.txt`

For non-spin-polarized calculation (nspin=1):
- `dipole_s1.txt`

## Output Format

The dipole output file contains one line for each ionic/electronic step:

```
step_index dipole_x dipole_y dipole_z
```

- `step_index`: The current step number (starts from 1)
- `dipole_x, dipole_y, dipole_z`: The x, y, z components of the dipole moment

Example output:

```
1 -0.00123456 0.00234567 -0.00345678
2 -0.00123457 0.00234568 -0.00345679
...
```

## Additional Information

During the calculation, the dipole moment is also printed in the `running_*.log` file, including:
- Electronic dipole moment
- Ionic dipole moment
- Total dipole moment
- Total dipole moment norm

The dipole moment calculation includes both electronic and ionic contributions. The total dipole moment is the sum of electronic and ionic dipoles.
