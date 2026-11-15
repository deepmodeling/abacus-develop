# abacuslite

## Introduction

abacuslite is a lightweight plugin for ABACUS (Atomic-orbital Based Ab-initio Computation at UStc), implementing the ASE (Atomic Simulation Environment) calculator interface.

### Key Features

- **Lightweight Design**: Implemented as a plugin, no need to modify ASE core code
- **Version Compatibility**: Requires `ase>=3.23`, which provides the GenericFileIO calculator API used by the plugin
- **ASE Integration**: Uses ASE as the running platform, making ABACUS a callable calculator within it
- **Function Support**: Provides SCF-based energy, force, and stress evaluations through ASE. ASE can use these evaluations for relaxation, molecular dynamics, NEB, band-structure, and density-of-states workflows.
- **Socket Support**: `AbacusSocketIO` provides fixed-cell i-PI socket calculations, with energy always available and forces/stress enabled independently when requested.

## Installation

Install the plugin from the ASE interface directory:

```bash
cd interfaces/ASE_interface
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
10. **soc.py** - Noncollinear spin-orbit coupling calculation
11. **socketio.py** - Fixed-cell ASE optimization with `AbacusSocketIO`, running ABACUS as an i-PI socket client

The regular `Abacus` calculator runs one ABACUS calculation for each ASE property evaluation. ASE controls the relaxation, molecular-dynamics, and other workflow steps. The socket calculator reuses one ABACUS process for position updates, while the cell and electronic-structure settings remain fixed for that calculator instance.

## Result collection

SCF results are selected by matching the written STRU to the running-log
structure, including periodic translations and the original ASE atom order.
Unmatched or malformed SCF output raises an error. Native relaxation and MD
continue to use the final frame; native MD uses MD_dump when the running log
does not contain the trajectory. Each calculator retains its own LTS/develop
reader, so calculators for different engine versions can coexist.

The default SCF identity tolerance is **1e-4 Angstrom**, applied to each cell
component and periodic Cartesian position residual with no relative tolerance.
This accommodates common six-significant-digit running-log rounding; the
previously used 1e-5 bound rejects the retained LTS rounding fixture. It is a
structure-identity tolerance, separate from SCF, force and optimization
convergence settings. It is not a universal printing-error bound: larger cells,
non-unit lattice constants and skewed-cell coordinate transformations can
accumulate more error.

Result parsing does **not read abacus.json**. Missing or malformed JSON does not
affect collection. The input writer only backs up an old JSON file to prevent
downstream tools from treating it as output of a new run.

Both ASE `energy` and `free_energy` retain ABACUS's force-consistent
`E_KohnSham`. When an entropy term is present, `energy_without_entropy` is
provided separately; a missing term remains unknown. Noncollinear vector
moments are retained without modifying ASE's global property definitions.
STRU writing preserves representable species masses and converts ASE
velocities to atomic units. Different masses within one ABACUS species are
rejected because one species entry cannot represent them.

Parser and calculator regressions run without an engine:

```bash
python -m pip install pytest
PYTHONPATH=interfaces/ASE_interface python -m pytest interfaces/ASE_interface/tests/unit -q
```

## Authors

- Yuyang Ji
- Zhenxiong Shen
- Yike Huang
- Zhaoqing Liu

## Acknowledgments

Thanks to the ABACUS development team for their support and contributions.

## License

The applicable license terms are provided in the repository [LICENSE](../../LICENSE).

## Contact

If you have any questions or suggestions, please contact us through:

- GitHub: [deepmodeling/abacus-develop](https://github.com/deepmodeling/abacus-develop)
