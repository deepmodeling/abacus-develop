# abacuslite

## Introduction

abacuslite is a lightweight plugin for ABACUS (Atomic-orbital Based Ab-initio Computation at UStc), implementing the ASE (Atomic Simulation Environment) calculator interface.

### Key Features

- **Lightweight Design**: Implemented as a plugin, no need to modify ASE core code
- **Version Compatibility**: Requires `ase>=3.23`, because ASE 3.22 lacks `ase.calculators.genericfileio`, which this plugin imports
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

The ASE calculator uses the final parsed frame for SCF, native relaxation, and
native MD. The readers restore ASE atom order using the mapping created when
the calculator writes STRU. Each calculator retains its own LTS/develop reader.
For native MD, the geometry in `MD_dump` must cover the steps in
`running_md.log`; use `md_dumpfreq=1`. The readers do not align sparse
geometry dumps. Older LTS logs with complete geometry can be read without
`MD_dump`.

For the develop reader, keep any `eig_occ.txt` or `eig_occg{N}.txt` files
with the log. Eigenvalues and occupations are matched by ionic step. If
no spectrum was written for a step, its energy and forces can still be
read, but its eigenvalues and occupations are unavailable.

Both output readers accept paths or named text streams, using the log filename
to locate companion output files. Caller-owned streams remain open and are read
from their current position; seek to the start before rereading. Returned ASE
positions and cells are Cartesian Angstrom values, with Direct and alat
conversions applied once by the reader. Raw trajectory dictionaries retain
their file-native conventions.

Returned structures represent the simulation box with periodic boundaries in
all directions, including molecules in vacuum and slabs. This changes earlier
SCF output flags; readers do not infer physical isolation from vacuum thickness.
For nonperiodic geometry analysis, make a copy and explicitly set its `pbc`
flags; this does not change the completed calculation's boundary conditions.

Result parsing does **not read abacus.json**. Missing or malformed JSON does not
affect collection. Writing new input leaves existing JSON files untouched.

Both ASE `energy` and `free_energy` retain ABACUS's force-consistent
`E_KohnSham`. When the log reports `E_entropy(-TS)`, abacuslite also provides
`energy_without_entropy = E_KohnSham - E_entropy(-TS)` in `Atoms.info` on
parsed frames and in calculator results. This is not an ASE property or
ABACUS's `E_KS(sigma->0)`. The Fermi level uses ASE's `fermi_level` property;
`efermi` is a compatibility alias in the calculator results dictionary only.
Noncollinear vector moments are retained without modifying ASE's global
property definitions.
STRU writing preserves representable species masses and converts ASE
velocities to atomic units. Different masses within one ABACUS species are
rejected because one species entry cannot represent them.

The existing inline tests (`abacuslite/xtest.sh`) and engine integration tests
(`tests/xtest.sh`) are run by the ASE plugin workflow. Standalone offline
regressions and recorded output fixtures are maintained in
[QuantumMisaka/abacuslite](https://github.com/QuantumMisaka/abacuslite).
The calculator and parser corrections were adapted from standalone revisions
`39ab1d5` and `5d3da5d`, with follow-up ASE compatibility fixes.

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
