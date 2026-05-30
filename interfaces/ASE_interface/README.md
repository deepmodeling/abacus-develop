# abacuslite

## Introduction

abacuslite is a lightweight plugin for ABACUS (Atomic-orbital Based Ab-initio Computation at UStc), implementing the ASE (Atomic Simulation Environment) calculator interface.

### Key Features

- **Lightweight Design**: Implemented as a plugin, no need to modify ASE core code
- **Version Compatibility**: No longer restricted to specific ASE versions, works with most ASE versions
- **ASE Integration**: Uses ASE as the running platform, making ABACUS a callable calculator within it
- **Function Support**: Currently only supports SCF (Self-Consistent Field) functionality, returning energy, forces, stress, etc.

## Installation

Installation is very simple, just execute the following command in the project root directory:

```bash
pip install .
```

## Usage Examples

Please refer to the example scripts in the `examples` folder. Recommended learning path:

1. **scf.py** - Basic SCF calculation example
2. **relax.py** - Atomic position relaxation calculation
3. **cellrelax.py** - Cell parameter relaxation calculation
4. **bandstructure.py** - Band structure calculation
5. **dos.py** - Density of states calculation
6. **md.py** - Molecular dynamics simulation
7. **constraintmd.py** - Constrained molecular dynamics simulation
8. **metadynamics.py** - Metadynamics simulation
9. **neb.py** - Nudged Elastic Band (NEB) calculation

More usage examples will be provided in future versions.

## Authors

- Yuyang Ji
- Zhenxiong Shen
- Yike Huang
- Zhaoqing Liu

## Acknowledgments

Thanks to the ABACUS development team for their support and contributions.

## License

[Fill in according to the actual project license]

## Contact

If you have any questions or suggestions, please contact us through:

- GitHub: [deepmodeling/abacus-develop](https://github.com/deepmodeling/abacus-develop)