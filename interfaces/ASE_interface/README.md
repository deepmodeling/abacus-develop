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

The independent [abacuslite](https://pypi.org/project/abacuslite/1.1.0/) 1.1.0
release can be installed in a separate environment. The engine compatibility
runner requires the same wheel file that was installed in that environment:

```bash
python -m pip download --only-binary=:all: --no-deps \
  --dest /tmp/abacuslite-wheel 'abacuslite==1.1.0'
python -m pip install /tmp/abacuslite-wheel/*.whl
python /path/to/abacus-develop/interfaces/ASE_interface/tests/run_engine_compatibility.py \
  --engine /absolute/path/to/abacus --mpi-command 'mpirun -np 2' \
  --assets /path/to/abacus-develop/tests/PP_ORB --format latest \
  --wheel /tmp/abacuslite-wheel/abacuslite-1.1.0-py3-none-any.whl \
  --report-dir /tmp/abacus-ase-engine-tests --engine-ref FULL_ENGINE_COMMIT_SHA
```

The runner copies the required scripts outside the checkout, verifies installed
runtime bytes against the supplied wheel, records the executable hash and
actual version/reader, and preserves every raw run and per-script log. It
rejects omitted or skipped required cases. The ordinary external CI lane pins
the published 1.1.0 wheel by SHA256
`c3e4504e9c9f5bceb69f70c79290e7843037ca72f625c5f66258e14d110f0ff3`.
The bundled transition lane builds a wheel from this engine PR and runs the
same develop suite. Both develop lanes require all nine checks, including
fixed-cell socket and real native MD. Use `--format legacy` with an LTS
executable for eight interface checks plus the seeded initial-velocity case.
The LTS engine source has no bundled interface or socket driver; the shared
runner supplies its legacy-format oracle for plugin CI.
`workflow_dispatch` can select a reviewed plugin source ref instead of the
published release. Engine tests run only through this explicit runner or the
opt-in `xtest.sh`, not ordinary offline unit checks.

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
