# Co6Nb4Os2 Spin-Polarized Calculation Example

This example demonstrates a spin-polarized DFT calculation for the Co6Nb4Os2 compound.

## Special Note: Sensitivity to Smearing Parameter

This calculation is **extremely sensitive to the smearing parameter** (`smearing_sigma`). The total magnetic moment shows significant variation depending on the smearing value:

| smearing_sigma (Ry) | Total Magnetism (μB) |
|---------------------|----------------------|
| 0.002               | ~2.20                |
| 0.005               | ~1.94                |
| 0.010               | ~0.60                |
| 0.020               | ~0.00 (non-magnetic) |

### Explanation

The Co6Nb4Os2 system exhibits a delicate balance between magnetic and non-magnetic states. The smearing parameter controls the width of the Gaussian smearing used for occupation numbers, which significantly affects the electronic structure near the Fermi level. For this compound:

- **Small smearing values** (e.g., 0.002 Ry) preserve the magnetic ordering, resulting in a total magnetic moment of ~2.20 μB
- **Increasing smearing** progressively suppresses the magnetic moment
- **Large smearing values** (e.g., 0.02 Ry) completely quench the magnetism, leading to a non-magnetic ground state

### Recommended Settings

The default `smearing_sigma` in `INPUT` is set to **0.001 Ry**, which should yield a magnetic moment close to 2.20 μB. Users are advised to carefully test different smearing values when studying this system.

## System Details

- **Composition**: Co6Nb4Os2
- **Structure**: Hexagonal lattice
- **Spin polarization**: Enabled (`nspin = 2`)
- **Basis type**: LCAO

## Files

- `INPUT`: Main input file with calculation parameters
- `STRU`: Structure file with atomic positions, initial magnetic moments, and numerical atomic orbital definitions
- `KPT`: k-point mesh definition

### Numerical Atomic Orbitals (NAO)

The `STRU` file references three numerical atomic orbital files:
- `Co_gga_9au_60Ry_4s2p2d1f.orb`
- `Nb_gga_9au_60Ry_4s2p2d1f.orb`
- `Os_gga_9au_60Ry_4s2p2d1f.orb`

The filename format encodes the following information:
- `9au`: Orbital cutoff radius (9 atomic units)
- `60Ry`: Energy cutoff (60 Rydberg)
- `4s2p2d1f`: Number of orbitals per angular momentum channel (4 s-orbitals, 2 p-orbitals, 2 d-orbitals, 1 f-orbital)

The `ecutwfc` parameter in `INPUT` is set to **60** to match the energy cutoff of the orbital files, ensuring consistency in the calculation.

## Running the Calculation

Execute ABACUS in this directory to perform the spin-polarized SCF calculation. Monitor the output for the total magnetic moment in the final SCF iteration.
