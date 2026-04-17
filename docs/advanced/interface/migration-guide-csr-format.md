# Migration Guide: New CSR Format for H(R) and S(R) Matrices

## Overview

Starting from ABACUS v3.9.0.25, the output format for Hamiltonian H(R) and overlap S(R) matrices has been unified to use standard CSR (Compressed Sparse Row) format, matching the format used by `out_dmr` for density matrices.

This change affects downstream tools that read H(R) and S(R) matrices, including:
- TB2J (magnetic exchange parameters)
- DeepH (machine learning Hamiltonian)
- pyATB (tight-binding analysis)
- Custom analysis scripts

## What Changed

### File Names

**Old format (ABACUS ≤ v3.8.x):**
```
OUT.${suffix}/data-HR-sparse_SPIN0.csr
OUT.${suffix}/data-HR-sparse_SPIN1.csr  (nspin=2 only)
OUT.${suffix}/data-SR-sparse_SPIN0.csr
```

**New format (ABACUS ≥ v3.9.0.25):**
```
OUT.${suffix}/hrs1_nao.csr
OUT.${suffix}/hrs2_nao.csr  (nspin=2 only)
OUT.${suffix}/srs1_nao.csr
```

### File Format

Both old and new formats use CSR structure, but with different headers and metadata.

#### Old Format Header
```
STEP: 0
Matrix Dimension of H(R): 26
Matrix number of H(R): 183
0 0 0 41
[sparse data]
```

#### New Format Header
```
 --- Ionic Step 1 ---
 # print H matrix in real space H(R)
 1 # number of spin directions
 1 # spin index
 26 # number of localized basis
 183 # number of Bravais lattice vector R

 user_defined_lattice
 5.39761
 0 0.5 0.5
 0.5 0 0.5
 0.5 0.5 0
 Si
 2
 Direct
 0 0 0
 0.25 0.25 0.25

 #----------------------------------------------------------------------#
 #                               CSR Format                             #
 # The outer loop corresponds to the number of Bravais lattice vectors. #
 # The first line contains the index of the Bravais lattice vector      #
 # (Rx, Ry, Rz), followed by the number of non-zero elements.           #
 # The subsequent lines consist of three blocks of data, which are      #
 # values, column indices, row pointers.                                #
 #----------------------------------------------------------------------#

 0 0 0 41
 # CSR values
 [values]
 # CSR column indices
 [indices]
 # CSR row pointers
 [pointers]
```

### Key Differences

1. **UnitCell Information**: New format includes complete unit cell information (lattice vectors, atomic positions)
2. **Header Format**: New format uses descriptive comments with `#` prefix
3. **Section Labels**: New format explicitly labels CSR sections ("# CSR values", "# CSR column indices", "# CSR row pointers")
4. **Ionic Step**: New format uses "Ionic Step N" instead of "STEP: N"
5. **Precision Control**: New format supports optional precision parameter: `out_mat_hs2 1 12` (default 8)

## Migration Steps for Tool Developers

### 1. Detect Format Version

```python
def detect_format_version(filename):
    """Detect whether file uses old or new CSR format."""
    with open(filename, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('---'):
            return 'new'  # v3.9.0.25+
        elif first_line.startswith('STEP:'):
            return 'old'  # v3.8.x and earlier
        else:
            raise ValueError(f"Unknown format in {filename}")
```

### 2. Parse New Format

```python
def parse_new_csr_format(filename):
    """Parse new CSR format (ABACUS v3.9.0.25+)."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse header
    ionic_step = None
    nspin = None
    ispin = None
    nbasis = None
    nR = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Parse ionic step
        if line.startswith('---') and 'Ionic Step' in line:
            ionic_step = int(line.split()[-2])

        # Parse metadata
        elif '#' in line:
            parts = line.split('#')
            if len(parts) == 2:
                value_str = parts[0].strip()
                comment = parts[1].strip()

                if 'number of spin directions' in comment:
                    nspin = int(value_str)
                elif 'spin index' in comment:
                    ispin = int(value_str)
                elif 'number of localized basis' in comment:
                    nbasis = int(value_str)
                elif 'number of Bravais lattice vector R' in comment:
                    nR = int(value_str)

        # Parse unit cell (optional, can be skipped if not needed)
        elif line == 'user_defined_lattice':
            # Read lattice constant and vectors
            i += 1
            lat0 = float(lines[i].strip())
            latvec = []
            for j in range(3):
                i += 1
                latvec.append([float(x) for x in lines[i].split()])
            # Continue parsing atomic positions if needed...

        # Start of CSR data
        elif line.startswith('#------') and 'CSR Format' in line:
            i += 1
            break

        i += 1

    # Parse R vectors and CSR data
    R_vectors = []
    HR_data = {}

    while i < len(lines):
        line = lines[i].strip()

        # R vector line: Rx Ry Rz nnz
        if line and not line.startswith('#') and len(line.split()) == 4:
            parts = [int(x) for x in line.split()]
            Rx, Ry, Rz, nnz = parts
            R_key = (Rx, Ry, Rz)
            R_vectors.append(R_key)

            # Read CSR values
            i += 1
            while lines[i].strip().startswith('#'):
                i += 1
            values = []
            while i < len(lines) and not lines[i].strip().startswith('#'):
                values.extend([float(x) for x in lines[i].split()])
                i += 1
                if len(values) >= nnz:
                    break

            # Read column indices
            while lines[i].strip().startswith('#'):
                i += 1
            col_indices = []
            while i < len(lines) and not lines[i].strip().startswith('#'):
                col_indices.extend([int(x) for x in lines[i].split()])
                i += 1
                if len(col_indices) >= nnz:
                    break

            # Read row pointers
            while lines[i].strip().startswith('#'):
                i += 1
            row_ptrs = []
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split()
                if len(parts) == 4 and all(x.isdigit() or x.lstrip('-').isdigit() for x in parts[:3]):
                    # Next R vector
                    break
                row_ptrs.extend([int(x) for x in parts])
                i += 1
                if len(row_ptrs) >= nbasis + 1:
                    break

            HR_data[R_key] = {
                'values': values[:nnz],
                'col_indices': col_indices[:nnz],
                'row_ptrs': row_ptrs[:nbasis+1]
            }
        else:
            i += 1

    return {
        'ionic_step': ionic_step,
        'nspin': nspin,
        'ispin': ispin,
        'nbasis': nbasis,
        'nR': nR,
        'R_vectors': R_vectors,
        'data': HR_data
    }
```

### 3. Backward Compatibility

To support both old and new formats:

```python
def read_hamiltonian(filename):
    """Read Hamiltonian matrix supporting both old and new formats."""
    format_version = detect_format_version(filename)

    if format_version == 'new':
        return parse_new_csr_format(filename)
    else:
        return parse_old_csr_format(filename)  # Your existing parser
```

## Tool-Specific Updates

### TB2J

**Required version:** TB2J v0.9.0+

**Changes needed:**
- Update file name detection to look for `hrs*_nao.csr` and `srs*_nao.csr`
- Update parser to handle new header format with UnitCell information
- Update parser to skip comment lines starting with `#`

**Example:**
```python
# In abacus2J.py or relevant parser
def find_hr_files(path, suffix):
    """Find HR files supporting both old and new formats."""
    import os
    import glob

    # Try new format first (v3.9.0.25+)
    new_pattern = os.path.join(path, f"OUT.{suffix}", "hrs*_nao.csr")
    hr_files = glob.glob(new_pattern)

    if hr_files:
        return sorted(hr_files), 'new'

    # Fall back to old format
    old_pattern = os.path.join(path, f"OUT.{suffix}", "data-HR-sparse_SPIN*.csr")
    hr_files = glob.glob(old_pattern)

    return sorted(hr_files), 'old'
```

### DeepH

**Required version:** DeepH v1.0.0+

**Changes needed:**
- Update `parse_abacus.py` to handle new file names
- Update CSR parser to skip UnitCell section
- Update parser to handle comment lines with `#` prefix

### pyATB

**Changes needed:**
- Update file I/O module to detect and parse new format
- Add format version detection
- Maintain backward compatibility with old format

## Testing Your Migration

### 1. Generate Test Files

Run ABACUS with both old and new versions:

```bash
# Old version (v3.8.x)
abacus_old > log_old

# New version (v3.9.0.25+)
abacus_new > log_new
```

### 2. Verify Numerical Equivalence

The CSR data (values, indices, pointers) should be numerically identical between old and new formats, only the header differs.

```python
def compare_csr_data(old_file, new_file):
    """Verify that CSR data is numerically equivalent."""
    old_data = parse_old_csr_format(old_file)
    new_data = parse_new_csr_format(new_file)

    import numpy as np

    for R in old_data['R_vectors']:
        old_vals = np.array(old_data['data'][R]['values'])
        new_vals = np.array(new_data['data'][R]['values'])

        assert np.allclose(old_vals, new_vals, rtol=1e-10), \
            f"Values differ for R={R}"

        assert old_data['data'][R]['col_indices'] == new_data['data'][R]['col_indices'], \
            f"Column indices differ for R={R}"

        assert old_data['data'][R]['row_ptrs'] == new_data['data'][R]['row_ptrs'], \
            f"Row pointers differ for R={R}"

    print("✓ CSR data is numerically equivalent")
```

## Additional Features

### Precision Control

The new format supports precision control via the second parameter:

```
out_mat_hs2 1 8   # 8 digits (default)
out_mat_hs2 1 12  # 12 digits (higher precision)
out_mat_hs2 1 5   # 5 digits (lower precision, smaller files)
```

This affects the output format of floating-point values in the CSR data.

## Support and Resources

- **ABACUS Documentation:** [https://abacus.deepmodeling.com/](https://abacus.deepmodeling.com/)
- **GitHub Issues:** [https://github.com/deepmodeling/abacus-develop/issues](https://github.com/deepmodeling/abacus-develop/issues)
- **Design Document:** `docs/plans/2026-02-28-unify-out-mat-hs2-design.md`
- **Implementation Plan:** `docs/plans/2026-02-28-unify-out-mat-hs2-plan.md`

## Summary

The new CSR format provides:
- ✅ Unified interface with `out_dmr` output
- ✅ Complete UnitCell information in output files
- ✅ Precision control for output values
- ✅ Better documentation with inline comments
- ✅ Clearer section labels for easier parsing

Tool developers should update their parsers to support the new format while maintaining backward compatibility with the old format for users running older ABACUS versions.
