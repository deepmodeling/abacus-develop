# PyABACUS

PyABACUS is the official Python interface for ABACUS, providing a convenient way to run DFT calculations directly from Python scripts.

## Installation

### From Source (Recommended)

```bash
cd /path/to/abacus-develop/python/pyabacus
pip install -e .
```

### With C++ Driver Support

For full functionality including direct library calls (faster than subprocess), build ABACUS with Python bindings:

```bash
cmake -B build -DENABLE_PYABACUS=ON -DENABLE_LCAO=ON
cmake --build build -j8
pip install -e python/pyabacus
```

**Note:** The `pyabacus` package on PyPI is a different project and is NOT related to ABACUS. Please install from source as shown above.

## Quick Start

### Basic SCF Calculation

```python
import pyabacus

# Run calculation from a directory containing INPUT, STRU, KPT files
result = pyabacus.abacus("./Si_scf/")

# Check results
print(f"Converged: {result.converged}")
print(f"Total energy: {result.etot_ev:.6f} eV")
print(result.summary())
```

### Calculate Forces and Stress

```python
result = pyabacus.abacus(
    "./Si_relax/",
    calculate_force=True,
    calculate_stress=True,
)

# Access forces (in eV/Angstrom)
if result.has_forces:
    forces = result.forces_ev_ang
    print(f"Max force: {forces.max():.6f} eV/Ang")

# Access stress tensor (in kbar)
if result.has_stress:
    print(f"Stress tensor:\n{result.stress}")
```

### Parallel Calculation

```python
# Run with MPI and OpenMP parallelization
result = pyabacus.abacus(
    "./Si_scf/",
    nprocs=4,      # 4 MPI processes (mpirun -np 4)
    nthreads=2,    # 2 OpenMP threads (OMP_NUM_THREADS=2)
)
```

This is equivalent to running:
```bash
OMP_NUM_THREADS=2 mpirun -np 4 abacus
```

### Silent Mode

```python
# Run without output
result = pyabacus.abacus("./Si_scf/", verbosity=0)
```

## API Reference

### `pyabacus.abacus()`

Main function for running ABACUS calculations.

```python
def abacus(
    input_dir: str = None,
    *,
    input_file: str = None,
    stru_file: str = None,
    kpt_file: str = None,
    pseudo_dir: str = None,
    orbital_dir: str = None,
    output_dir: str = None,
    calculate_force: bool = True,
    calculate_stress: bool = False,
    verbosity: int = 1,
    nprocs: int = 1,
    nthreads: int = 1,
) -> CalculationResult
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | str | `"."` | Directory containing INPUT, STRU, KPT files |
| `input_file` | str | None | Explicit path to INPUT file |
| `stru_file` | str | None | Explicit path to STRU file |
| `kpt_file` | str | None | Explicit path to KPT file |
| `pseudo_dir` | str | None | Directory containing pseudopotentials |
| `orbital_dir` | str | None | Directory containing orbital files (LCAO) |
| `output_dir` | str | `"OUT.PYABACUS"` | Directory for output files |
| `calculate_force` | bool | True | Whether to calculate forces |
| `calculate_stress` | bool | False | Whether to calculate stress tensor |
| `verbosity` | int | 1 | Output level (0=silent, 1=normal, 2=verbose) |
| `nprocs` | int | 1 | Number of MPI processes |
| `nthreads` | int | 1 | Number of OpenMP threads |

**Returns:** `CalculationResult` object

### `CalculationResult`

Container for calculation results.

**Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `converged` | bool | Whether SCF converged |
| `niter` | int | Number of SCF iterations |
| `drho` | float | Final charge density difference |
| `etot` | float | Total energy (Ry) |
| `etot_ev` | float | Total energy (eV) |
| `forces` | ndarray | Forces on atoms (nat, 3) in Ry/Bohr |
| `forces_ev_ang` | ndarray | Forces in eV/Angstrom |
| `stress` | ndarray | Stress tensor (3, 3) in kbar |
| `fermi_energy` | float | Fermi energy (eV) |
| `bandgap` | float | Band gap (eV) |
| `nat` | int | Number of atoms |
| `ntype` | int | Number of atom types |
| `nbands` | int | Number of bands |
| `nks` | int | Number of k-points |
| `output_dir` | str | Path to output directory (OUT.$suffix) |
| `log_file` | str | Path to the main log file |
| `output_files` | dict | Dictionary of output files (filename -> path) |

**Methods:**

| Method | Description |
|--------|-------------|
| `summary()` | Return a formatted summary string |
| `energies` | Dictionary of all energy components |
| `has_forces` | Whether forces are available |
| `has_stress` | Whether stress is available |
| `has_output_dir` | Whether output directory exists |
| `get_output_file(name)` | Get full path to specific output file |
| `list_output_files()` | List all output file names |

## Output File Tracking

PyABACUS automatically tracks output files generated during calculations:

```python
result = pyabacus.abacus("./Si_scf/")

# Check output directory
print(f"Output directory: {result.output_dir}")
print(f"Log file: {result.log_file}")

# List all output files
print("Output files:")
for filename in result.list_output_files():
    print(f"  {filename}")

# Get path to specific file
bands_file = result.get_output_file("BANDS_1.dat")
if bands_file:
    # Read and process band structure data
    import numpy as np
    bands = np.loadtxt(bands_file)
```

### Common Output Files

| File | Description |
|------|-------------|
| `running_scf.log` | Main calculation log |
| `BANDS_1.dat` | Band structure data |
| `PDOS` | Projected density of states |
| `CHARGE.cube` | Charge density in cube format |
| `SPIN1_CHG.cube` | Spin-up charge density |
| `SPIN2_CHG.cube` | Spin-down charge density |
| `istate.info` | Band eigenvalues and occupations |
| `kpoints` | K-point information |

## Convenience Functions

### `run_scf()`

Alias for `abacus()` with default SCF parameters:

```python
result = pyabacus.run_scf("./Si_scf/")
```

### `run_relax()`

Alias for `abacus()` with force calculation enabled:

```python
result = pyabacus.run_relax("./Si_relax/")
```

## Examples

### Energy vs. Lattice Constant

```python
import pyabacus
import numpy as np

lattice_constants = np.linspace(5.0, 5.5, 11)
energies = []

for a in lattice_constants:
    # Modify STRU file with new lattice constant
    # ... (file modification code)

    result = pyabacus.abacus("./Si_eos/", verbosity=0)
    energies.append(result.etot_ev)

# Plot equation of state
import matplotlib.pyplot as plt
plt.plot(lattice_constants, energies, 'o-')
plt.xlabel("Lattice constant (Ang)")
plt.ylabel("Energy (eV)")
plt.savefig("eos.png")
```

### Parallel Batch Calculations

```python
import pyabacus
from pathlib import Path

# Run calculations for multiple systems with parallelization
systems = ["Si", "Ge", "C"]
results = {}

for system in systems:
    input_dir = Path(f"./{system}_scf/")
    if input_dir.exists():
        result = pyabacus.abacus(
            str(input_dir),
            nprocs=4,
            nthreads=2,
        )
        results[system] = {
            "energy": result.etot_ev,
            "converged": result.converged,
            "bandgap": result.bandgap,
        }

# Print summary
for system, data in results.items():
    print(f"{system}: E={data['energy']:.4f} eV, gap={data['bandgap']:.2f} eV")
```

## Troubleshooting

### ABACUS executable not found

If you see "ABACUS executable not found", ensure:
1. ABACUS is installed and in your PATH
2. Or build with C++ driver support (see Installation)

### MPI not found

If you see "mpirun/mpiexec not found" when using `nprocs > 1`:
1. Install MPI (OpenMPI or MPICH)
2. Ensure `mpirun` or `mpiexec` is in your PATH
3. Or set `nprocs=1` to run without MPI

### Import errors

If `import pyabacus` fails:
1. Ensure pyabacus is installed: `pip install pyabacus`
2. Check Python version compatibility (Python 3.8+)

### Calculation not converging

Check the log file for details:
```python
result = pyabacus.abacus("./problem_case/")
if not result.converged:
    log_path = result.get_output_file("running_scf.log")
    if log_path:
        with open(log_path) as f:
            print(f.read()[-2000:])  # Print last 2000 chars
```

### Forces or stress not available

Forces and stress are parsed from the ABACUS output log. Ensure:
1. `cal_force` is set in your INPUT file for forces
2. `cal_stress` is set in your INPUT file for stress
3. The calculation completed successfully
