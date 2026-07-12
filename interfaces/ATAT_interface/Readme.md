# ATAT-ABACUS Interface

## Introduction

`runstruct_abacus` is a lightweight interface script connecting **ATAT (Alloy Theoretic Automated Toolkit)** with **ABACUS** (Atomic-orbital Based Ab-initio Computation at UStc). It automatically converts ATAT's `str.out` structure files into ABACUS `INPUT`/`STRU` input files, runs the DFT calculation, and extracts results back into ATAT-compatible formats (`energy`, `str_relax.out`).

### Key Features

- **Seamless ATAT Integration**: Works within ATAT's multi-directory enumeration workflow (`1/`, `2/`, `3/`...)
- **Template-Based Input**: Uses `abacus.wrap` as a template—nearly a native ABACUS `INPUT` file with minimal script-specific annotations
- **Automatic File Discovery**: Searches `abacus.wrap` upward through parent directories (`./` → `../` → `../../`...)
- **Smart Pseudopotential/Orbital Matching**: Auto-detects files in `pseudo_dir`/`orbital_dir` by element prefix; explicit override available for ambiguous cases
- **Flexible Execution Modes**: Supports full pipeline, input-only generation, and post-calculation extraction
- **Parallel Ready**: Accepts `mpirun`/`srun` prefixes for HPC environments

---

## Installation

No installation is required. Simply place `runstruct_abacus` in your `$PATH` (or in the same directory as other ATAT `runstruct_*` scripts) and ensure it is executable:

```bash
chmod +x runstruct_abacus
```

### Dependencies

- ATAT toolkit (`cellcvrt`, `kmesh` etc.) must be in `$PATH`
- ABACUS executable path must be set in `~/.abacus.rc`

---

## Configuration

### `~/.abacus.rc`

Create this file in your home directory to tell the interface where ABACUS lives:

```bash
#!/bin/bash
ABACUSCMD="abacus"                    # or "mpirun -np 4 abacus"
```

The script will auto-generate a template if this file does not exist.

---

## Template File: `abacus.wrap`

`abacus.wrap` is **almost** a standard ABACUS `INPUT` file. The script copies nearly every line verbatim into `INPUT`, except for `species` lines which are consumed by the script to build the `STRU` file.

### Minimal Example

```bash
INPUT_PARAMETERS
calculation    vc-relax
ecutwfc        50
basis_type     lcao
kspacing       0.15
pseudo_dir     /path/to/pseudopotentials
orbital_dir    /path/to/numerical_orbitals

species  Al  26.982  Al_ONCV_PBE-1.0.upf  Al_gga_7au_60Ry_2s2p1d.orb
species  Fe  55.845  Fe_ONCV_PBE-1.0.upf  Fe_gga_8au_100Ry_2s2p2d1f.orb
```

### `species` Syntax

```bash
species  <Element>  <Mass>  <PP_File>  [<Orb_File>]
```

| Field      | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| `Element`  | Chemical symbol (e.g., `Al`, `Fe`)                           |
| `Mass`     | Atomic mass. Use `-` to look up from the built-in table      |
| `PP_File`  | Pseudopotential filename. Use `-` to auto-search in `pseudo_dir` |
| `Orb_File` | Numerical orbital filename (required for LCAO). Use `-` to auto-search in `orbital_dir` |

If auto-search finds **zero** or **more than one** match for an element, the script aborts and prints a helpful message asking you to add an explicit `species` line.

---

## Command Line Options

```bash
runstruct_abacus [-w file] [-nr] [-ex] [-clean] [cmdprefix]
```

### Execution Modes

| Command                   | Behavior                                                     |
| ------------------------- | ------------------------------------------------------------ |
| `runstruct_abacus`        | **Full pipeline**: Generate `INPUT` + `STRU` → Run ABACUS → Extract `energy`, `str_relax.out` |
| `runstruct_abacus -nr`    | **No-Run**: Generate `INPUT` + `STRU` only. Useful for manual inspection or external job schedulers. |
| `runstruct_abacus -ex`    | **Extract-Only**: Skip generation and execution. Extract results from existing `OUT.suffix/` directory. |
| `runstruct_abacus -clean` | **Cleanup**: Delete all output files (`OUT.*/`, `running_*.log`, `energy`, `str_relax.out`) and exit. |

### `cmdprefix`: Running in Parallel

The optional `cmdprefix` argument lets you prepend any launch command—most commonly MPI wrappers:

```bash
# Run with 4 MPI ranks
runstruct_abacus "mpirun -np 4"

# Run with srun (SLURM)
runstruct_abacus "srun -n 8"

# Run on a specific node (similar to Abinit's node-prefix syntax)
runstruct_abacus "ssh node02 mpirun -np 16"
```

The prefix is inserted directly before `$ABACUSCMD`:

```bash
$CMDPREFIX $ABACUSCMD > log.out 2>&1
```

### `-w`: Custom Wrap File

```bash
runstruct_abacus -w my_custom.wrap
```

If the specified file is not found in the current directory, the script searches upward (`../`, `../../`, `../../../`) exactly like the default `abacus.wrap`.

---

## Workflow Example

### Standard ATAT Workflow

```bash
# Inside a numbered ATAT subdirectory, e.g., 1/, 2/, ...
cd 1/

# 1. Generate inputs and run
runstruct_abacus

# 2. Or generate only, then submit to cluster manually
runstruct_abacus -nr
# ... user submits job via qsub/sbatch ...
runstruct_abacus -ex   # extract after job finishes

# 3. Clean and restart if needed
runstruct_abacus -clean
runstruct_abacus
```

### Output Files

| File            | Description                                             |
| --------------- | ------------------------------------------------------- |
| `INPUT`         | ABACUS control parameters (filtered from `abacus.wrap`) |
| `STRU`          | ABACUS structure file (lattice, species, coordinates)   |
| `energy`        | Final total energy in **eV** (ATAT standard unit)       |
| `str_relax.out` | Relaxed structure in ATAT `str.out` format              |
| `log.out`       | Raw ABACUS stdout/stderr                                |

---

## File Search Hierarchy

Both `abacus.wrap` and `str.out` follow ATAT's upward-search convention:

| File          | Search Order                           |
| ------------- | -------------------------------------- |
| `abacus.wrap` | `./` → `../` → `../../` → `../../../`  |
| `str.out`     | `str_hint.out` (preferred) → `str.out` |

This allows a single `abacus.wrap` (and optionally a shared `~/.abacus.rc`) to serve an entire ATAT enumeration tree.

---

## Authors

- Shengjun Chen (陈胜君) @ Peking University

## License

[Fill in according to your project license]

## Contact

For issues related to the ABACUS engine itself, please visit:
- GitHub: [deepmodeling/abacus-develop](https://github.com/deepmodeling/abacus-develop)