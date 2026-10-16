# abacuslite

## Introduction

abacuslite is a lightweight plugin for ABACUS (Atomic-orbital Based Ab-initio Computation at UStc), implementing the ASE (Atomic Simulation Environment) calculator interface.

### Key Features

- **Lightweight Design**: Implemented as a plugin, no need to modify ASE core code
- **Version Compatibility**: Supports ASE versions satisfying the package requirement `ase>=3.22`
- **ASE Integration**: Uses ASE as the running platform, making ABACUS a callable calculator within it
- **Function Support**: Provides SCF-based energy, force, and stress evaluations through ASE. ASE can use these evaluations for relaxation, molecular dynamics, NEB, band-structure, and density-of-states workflows.
- **Socket Support**: `AbacusSocketIO` provides fixed-cell i-PI socket calculations, with energy always available and forces/stress enabled independently when requested.

## Installation

Install the plugin from the ASE interface directory:

```bash
cd interfaces/ASE_interface
pip install .
```

For coordinated testing of the independent candidate, build a wheel from the
reviewed [QuantumMisaka/abacuslite](https://github.com/QuantumMisaka/abacuslite)
commit and install it in a separate environment. The bundled installation above
remains available while the migration is reviewed:

```bash
git clone https://github.com/QuantumMisaka/abacuslite.git /tmp/abacuslite-candidate
git -C /tmp/abacuslite-candidate checkout dc020430bfa8c76edfc275b0f1eab2f12cc29d4d
python -m pip install build
python -m build --wheel --outdir /tmp/abacuslite-wheel /tmp/abacuslite-candidate
python -m pip install /tmp/abacuslite-wheel/*.whl
mkdir -p /tmp/abacus-ase-engine-tests
cp /path/to/abacus-develop/interfaces/ASE_interface/tests/*.py \
   /path/to/abacus-develop/interfaces/ASE_interface/tests/xtest.sh \
   /tmp/abacus-ase-engine-tests/
cd /tmp/abacus-ase-engine-tests
python verify_external.py \
  --bundled-root /path/to/abacus-develop/interfaces/ASE_interface \
  --candidate-root /tmp/abacuslite-candidate
ABACUS_COMMAND='mpirun -np 2 /path/to/abacus' \
ASE_ABACUS_ASSETS=/path/to/abacus-develop/tests/PP_ORB \
ASE_ABACUS_SOCKET=1 \
bash ./xtest.sh
```

Run the scripts from outside `interfaces/ASE_interface` so that the bundled
source cannot shadow the installed wheel. In CI the `external` job checks
distribution metadata and the imported `abacuslite.__file__`, records the engine
revision and `abacus --version`, and retains failed calculation directories.
`workflow_dispatch` accepts a reviewed `plugin_ref` for coordinated runs.
The scripts above require a built ABACUS executable and the curated PP/ORB
assets. Both lanes include a two-step ABACUS-native MD run that compares the
assembled periodic ASE frames with the raw `MD_dump` lattice, positions and
forces. The separate raw-output fixture regression is an offline Python check;
ordinary Python checks do not launch the engine.

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
