# ABACUS Output File Specification

## 1. Background and Motivation

ABACUS is an integrated software package designed to provide a cohesive user experience. To achieve this goal, **we are establishing unified output file standards** that all developers should follow.

Recent versions of ABACUS have been working on standardizing all output file naming conventions. This work is ongoing, and if there are discrepancies between actual file names and this document, please refer to the latest documentation.

All output file naming conventions can be found in the online documentation (with corresponding 3.10-LTS file names): [ABACUS Input Documentation](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html)

## 2. File Naming Conventions

### 2.1 Basic Rules

**Rule 1:** All ABACUS output files are stored in the `OUT.{suffix}/` directory.

**Rule 2:** File extensions by category:
| Extension | Description |
|-----------|-------------|
| `.txt` | Text file |
| `.dat` | Binary file |
| `.csr` | Sparse matrix format |
| `.cube` | 3D spatial data format |

**Rule 3:** For special output quantities (e.g., wavefunctions), add `_pw` or `_nao` to distinguish between plane wave basis and numerical atomic orbital basis.

**Rule 4:** File names are lowercase. Physical quantities appear at the beginning:

| Prefix | Physical Quantity |
|--------|-------------------|
| `chg` | Charge density |
| `pot` | Potential |
| `eig` | Eigenvalue / Energy level |
| `wf` | Wavefunction |
| `dm` | Density matrix |
| `h` | Hamiltonian matrix H |
| `s` | Overlap matrix S |
| `t` | Kinetic energy operator |
| `r` | Position operator or Bravais lattice vector R |
| `k` | k-point in Brillouin zone |
| `xyz` | Three spatial directions |
| `ini` | Initial state (before electronic iteration) |

**Rule 5:** Suffixes following the physical quantity:

| Suffix | Meaning |
|--------|---------|
| `s1`, `s2`, `s3`, `s4` | Spin channel (1, 2 for collinear; 1, 2, 3, 4 for non-collinear with SOC) |
| `s12` | Non-collinear spin calculation |
| `k#` | k-point index (e.g., `k1`, `k2`) |
| `g#` | Ionic step index for relax/md (e.g., `g1`, `g2`) |

**Important:**
- All index numbers start from 1 (not 0)
- For Gamma-only algorithm in LCAO, no `k` index is included
- Overlap matrix `s` does not distinguish spin, so only one matrix is output

### 2.2 Examples

| File Name | Interpretation |
|-----------|----------------|
| `chgs1.cube` | Charge density, spin 1 |
| `chgs2.cube` | Charge density, spin 2 |
| `chgs3.cube` | Charge density, spin 3 (non-collinear with SOC) |
| `pots1.cube` | Local potential, spin 1 |
| `eig_occ.txt` | Eigenvalues and occupations |
| `doss1g1_nao.txt` | DOS, spin 1, geometry step 1, NAO basis |
| `wf_pw.dat` | Wavefunction, plane wave basis |
| `sr.csr` | Overlap matrix in real space (no spin index) |

### 2.3 Common Output Files

| File Name | Description |
|-----------|-------------|
| `running_scf.log` | SCF iteration log |
| `dos.txt` | Density of states |
| `eig_occ.txt` | Eigenvalues and occupations |
| `mulliken.txt` | Mulliken population analysis |
| `band.txt` | Band structure |
| `chgs1.cube`, `chgs2.cube` | Charge density (spin 1, spin 2) |
| `chg.cube` | Total charge density |
| `taus1.cube`, `taus2.cube` | Kinetic energy density (tau) |
| `pots1.cube`, `pots2.cube` | Local potential |

## 3. File Format Standards

### 3.1 Header Section with Comments

**Every output file should include `#` comment lines** to explain:
- Data meaning
- Units
- Source module
- Key parameters

```
# <description of file content>
# Module: <source module name>
# Units: <unit information>
<value>    # <description>
# <column headers with units>
```

**Example:**
```
1    # ionic step
8207 # number of points
# Module: DOS calculation
# Units: energy in eV, DOS in states/eV
#        energy    elec_states     sum_states   states_smear     sum_states
```

### 3.2 Data Section

| Requirement | Description |
|-------------|-------------|
| **Separator** | Use spaces for column alignment |
| **Precision** | Controlled by input parameter (see Section 3.4) |
| **Units** | Always specify units in header or column names |
| **Comments** | Use `#` for comment lines |

### 3.3 Compact Data Layout

**Avoid sparse format with single value per line.** Instead, output 6-8 values per line with proper alignment:

**Bad (sparse):**
```
1.234567
2.345678
3.456789
4.567890
```

**Good (compact):**
```
1.234567    2.345678    3.456789    4.567890    5.678901    6.789012
7.890123    8.901234
```

### 3.4 Precision Control via Input Parameters

Output precision should be controllable via input parameters, similar to `out_chg` and `out_pot`:

**Example (from `out_pot`):**
```
out_pot 1 8
```
- First integer: output type (1 = output total local potential)
- Second integer: precision (number of significant digits, default 8)

**Implementation pattern:**
```cpp
// In input parameters
int out_type = 0;    // output type
int out_precision = 8;  // precision

// In output function
ofs << std::setprecision(out_precision) << value;
```

## 4. Output Volume Control

### 4.1 Reduce File Size

Current integration tests have output files with tens of thousands of lines. Recommendations:

| Issue | Solution |
|-------|----------|
| Too many output files | Consolidate related data into fewer files |
| Files too large | Reduce output frequency, use compact format |
| Redundant information | Avoid repeating header information in every block |

### 4.2 Balance Test Coverage and Efficiency

- Output only essential data for integration tests
- Use `out_level` parameter to control verbosity
- Consider binary format for large datasets

## 5. Naming Conventions for Physical Quantities

### 5.1 Standard Names (Keep Short: 3-8 Characters)

| Physical Quantity | Recommended Name | Unit |
|-------------------|------------------|------|
| Total energy | `etot` | eV |
| Kinetic energy | `ekin` | eV |
| Potential energy | `epot` | eV |
| Force | `force` | eV/Angstrom |
| Stress | `stress` | kBar |
| Charge | `chg` | e |
| Magnetization | `mag` | μB |
| Band index | `n` | - |
| k-point | `kpt` | - |
| Spin | `spin` | - |
| Occupation | `occ` | - |
| Energy level | `eig` | eV |

### 5.2 Naming Style

- Use `snake_case` (lowercase with underscores)
- **Keep names short (3-8 characters)**
- Avoid abbreviations unless widely understood

## 6. Developer Checklist

Before adding a new output file or modifying an existing one, verify:

- [ ] **File name is short (3-8 characters)**
- [ ] File name follows naming conventions (Section 2)
- [ ] File is stored in `OUT.{suffix}/` directory
- [ ] File numbering starts from 1 (not 0)
- [ ] No redundant keywords (e.g., `data`) in filename
- [ ] Use correct extension (`.txt`, `.dat`, `.csr`, `.cube`)
- [ ] Add `_pw` or `_nao` for basis-specific outputs
- [ ] Header section includes `#` comment lines with:
  - [ ] Data meaning
  - [ ] Units
  - [ ] Source module
- [ ] Column headers include units
- [ ] Data uses compact format (6-8 values per line, not single value per line)
- [ ] Output precision is controllable via input parameter
- [ ] File size is reasonable (avoid tens of thousands of lines)
- [ ] Physical quantities use standard names (Section 5)
- [ ] Documentation is updated in `docs/` directory

## 7. Code Implementation Guidelines

### 7.1 Output Function Template

```cpp
void OutputMyData::write(const std::string& filename, int precision)
{
    std::ofstream ofs(filename);
    
    // Header section with comments
    ofs << "# Density of states" << std::endl;
    ofs << "# Module: DOS" << std::endl;
    ofs << "# Units: energy in eV, DOS in states/eV" << std::endl;
    ofs << ionic_step << "    # ionic step" << std::endl;
    ofs << num_data << "    # number of data points" << std::endl;
    ofs << "#";
    ofs << std::setw(15) << "energy(eV)";
    ofs << std::setw(15) << "dos";
    ofs << std::endl;
    
    // Data section - compact format (6 values per line)
    ofs << std::setprecision(precision);
    for (int i = 0; i < num_data; ++i)
    {
        ofs << std::setw(15) << energy[i];
        ofs << std::setw(15) << dos[i];
        if ((i + 1) % 3 == 0) ofs << std::endl;  // 3 pairs = 6 values per line
    }
    if (num_data % 3 != 0) ofs << std::endl;
    
    ofs.close();
}
```

### 7.2 Key Points

1. **Keep file names short (3-8 characters)**
2. **Add `#` comment lines for data meaning, units, and source module**
3. **Use compact format: 6-8 values per line**
4. **Make precision controllable via input parameter**
5. **Use `std::setw()` for column alignment**

## 8. Existing Good Examples

| File | Location | Strengths |
|------|----------|-----------|
| `chgs1.cube` | `tests/03_NAO_multik/scf_out_chg_tau/OUT.autotest/` | Short name, spin index convention |
| `dos.txt` | `tests/03_NAO_multik/scf_out_dos_spin2/OUT.autotest/` | Clear header with metadata |
| `eig_occ.txt` | Same as above | Clear block structure for spin/k-point |
| `mulliken.txt` | `tests/03_NAO_multik/scf_out_mul/OUT.autotest/` | Clear atom separators |

## 9. Review Process

For new output formats:

1. Check if similar output already exists - reuse format if possible
2. Follow this specification
3. Add documentation in `docs/advanced/output_files/`
4. Submit PR with output sample for review

## 10. Summary

| Aspect | Requirement |
|--------|-------------|
| Directory | All outputs in `OUT.{suffix}/` |
| File name | **Short (3-8 chars)**, lowercase, physical quantity prefix |
| Extensions | `.txt`, `.dat`, `.csr`, `.cube` |
| Basis suffix | `_pw` or `_nao` for basis-specific outputs |
| Spin index | `s1`, `s2`, `s3`, `s4` (SOC), `s12` (non-collinear) |
| Numbering | Start from 1, not 0 |
| Header | `#` comments with meaning, units, source module |
| Data | **Compact: 6-8 values per line**, space-separated |
| Precision | Controllable via input parameter |
| Volume | Reasonable file size, avoid excessive output |

**Remember: ABACUS outputs should present a unified, professional interface to users.**
