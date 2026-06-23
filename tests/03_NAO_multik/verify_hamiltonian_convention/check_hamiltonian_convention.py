#!/usr/bin/env python3
"""
Simple Hamiltonian convention check.
Reads the first R=0 block of H(R) and checks:
1. Hermiticity of R=0 block
2. Off-diagonal phase for m along y
"""
import sys
import re
import numpy as np

def read_complex_value(s):
    """Parse a complex value like '(-1.13e-02,0.00e+00)' or a real number."""
    s = s.strip().strip('()')
    if ',' in s:
        parts = s.split(',')
        return complex(float(parts[0]), float(parts[1]))
    else:
        return complex(float(s), 0.0)

def parse_hr_r0(filepath):
    """Parse the R=(0,0,0) block from hrs1_nao.csr."""
    with open(filepath) as f:
        lines = f.readlines()

    # Parse header
    idx = 0
    nbasis = None
    nR = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#') or stripped == '' or stripped.startswith('---'):
            continue
        if 'number of spin' in stripped or 'spin index' in stripped:
            continue
        if 'number of localized basis' in stripped:
            nbasis = int(stripped.split()[0])
            continue
        if 'number of Bravais' in stripped:
            nR = int(stripped.split()[0])
            break

    # Find R=(0,0,0) block
    found_r0 = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Look for "0 0 0 nnz" pattern
        parts = stripped.split()
        if len(parts) == 4:
            try:
                rx, ry, rz, nnz = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
                if rx == 0 and ry == 0 and rz == 0 and nnz > 0:
                    found_r0 = True
                    break
            except ValueError:
                continue

    if not found_r0:
        print("R=(0,0,0) block not found!")
        return None, None, None

    # Now parse values, column indices, row pointers after this line
    idx = i + 1
    # Skip "# CSR values"
    while idx < len(lines) and '# CSR values' not in lines[idx]:
        idx += 1
    idx += 1

    # Read nnz complex values
    values = []
    while len(values) < nnz and idx < len(lines):
        line = lines[idx].strip()
        idx += 1
        if line.startswith('#'):
            break
        if not line:
            continue
        # Parse all complex values on this line
        tokens = re.findall(r'\([^)]+\)', line)
        if not tokens:
            # Try space-separated values
            tokens = line.split()
        for t in tokens:
            try:
                values.append(read_complex_value(t))
            except ValueError:
                break
        if len(values) >= nnz:
            break

    # Skip "# CSR column indices"
    while idx < len(lines) and '# CSR column indices' not in lines[idx]:
        idx += 1
    idx += 1

    # Read nnz column indices
    col_indices = []
    while len(col_indices) < nnz and idx < len(lines):
        line = lines[idx].strip()
        idx += 1
        if line.startswith('#'):
            break
        if not line:
            continue
        for t in line.split():
            try:
                col_indices.append(int(t))
            except ValueError:
                break
        if len(col_indices) >= nnz:
            break

    # Skip "# CSR row pointers"
    while idx < len(lines) and '# CSR row pointers' not in lines[idx]:
        idx += 1
    idx += 1

    # Read nbasis+1 row pointers
    row_ptrs = []
    while len(row_ptrs) < nbasis + 1 and idx < len(lines):
        line = lines[idx].strip()
        idx += 1
        if line.startswith('#'):
            break
        if not line:
            continue
        for t in line.split():
            try:
                row_ptrs.append(int(t))
            except ValueError:
                break
        if len(row_ptrs) >= nbasis + 1:
            break

    # Build dense matrix
    H = np.zeros((nbasis, nbasis), dtype=complex)
    for irow in range(nbasis):
        if irow >= len(row_ptrs):
            break
        start = row_ptrs[irow]
        end = row_ptrs[irow + 1] if irow + 1 < len(row_ptrs) else nnz
        for k in range(start, min(end, nnz, len(values), len(col_indices))):
            icol = col_indices[k]
            if 0 <= icol < nbasis:
                H[irow, icol] = values[k]

    return nbasis, H, nnz

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 check_hamiltonian_convention.py <output_dir>")
        sys.exit(1)

    out_dir = sys.argv[1]
    hr_file = None
    dmr_file = None

    import os
    for f in os.listdir(out_dir):
        if f.startswith('hrs') and f.endswith('.csr'):
            hr_file = os.path.join(out_dir, f)
        if f.startswith('dmrs') and f.endswith('.csr'):
            dmr_file = os.path.join(out_dir, f)

    if not hr_file:
        print("ERROR: hrs*_nao.csr not found")
        sys.exit(1)

    print(f"Parsing H(R=0) from: {hr_file}")
    nbasis, H0, nnz = parse_hr_r0(hr_file)
    if H0 is None:
        print("Failed to parse H(R=0)")
        sys.exit(1)

    print(f"  nbasis = {nbasis}, nnz = {nnz}")
    N = nbasis // 2  # number of spatial orbitals per spin

    # Test 1: Hermiticity of R=0 block
    print("\n" + "=" * 70)
    print("TEST 1: H(R=0) Hermiticity: H = H^dagger")
    print("=" * 70)
    herm_err = np.max(np.abs(H0 - H0.conj().T))
    print(f"  max|H - H^dagger| = {herm_err:.2e}")
    if herm_err < 1e-10:
        print("  [PASS]")
    else:
        print("  [FAIL]")

    # Test 2: Off-diagonal spinor block phase
    print("\n" + "=" * 70)
    print("TEST 2: Off-diagonal phase for m along y")
    print("  For m||+y: H_{up,down} should have Im < 0")
    print("  (i.e., H(2i, 2j+1) should be -i*|B_y|)")
    print("=" * 70)

    # Extract the up-down block: rows 0,2,4,..., cols 1,3,5,...
    H_updown = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            H_updown[i, j] = H0[2*i, 2*j+1]

    # Check the sign of imaginary part
    nonzero_mask = np.abs(H_updown) > 1e-10
    if np.any(nonzero_mask):
        imag_parts = H_updown[nonzero_mask].imag
        mean_imag = np.mean(imag_parts)
        n_neg = np.sum(imag_parts < 0)
        n_pos = np.sum(imag_parts > 0)
        print(f"  Non-zero H_{{up,down}} elements: {np.sum(nonzero_mask)}")
        print(f"  Mean Im(H_{{up,down}}) = {mean_imag:.6e}")
        print(f"  Negative: {n_neg}, Positive: {n_pos}")
        if mean_imag < -1e-10:
            print("  [PASS] Im(H_{{up,down}}) < 0, consistent with correct convention")
        elif mean_imag > 1e-10:
            print("  [FAIL] Im(H_{{up,down}}) > 0, indicates complex conjugation error")
        else:
            print("  [INFO] Im(H_{{up,down}}) ~ 0, magnetization may not be along y")
    else:
        print("  No significant off-diagonal elements found")

    # Test 3: Diagonal spinor block structure
    print("\n" + "=" * 70)
    print("TEST 3: Diagonal blocks H_{{up,up}} and H_{{down,down}}")
    print("  For m||y with no SOC: H_{{up,up}} = H_{{down,down}} (real)")
    print("=" * 70)
    H_upup = np.zeros((N, N), dtype=complex)
    H_downdown = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            H_upup[i, j] = H0[2*i, 2*j]
            H_downdown[i, j] = H0[2*i+1, 2*j+1]

    diag_diff = np.max(np.abs(H_upup - H_downdown))
    upup_imag = np.max(np.abs(H_upup.imag))
    dndn_imag = np.max(np.abs(H_downdown.imag))
    print(f"  max|H_{{up,up}} - H_{{down,down}}| = {diag_diff:.2e}")
    print(f"  max|Im(H_{{up,up}})| = {upup_imag:.2e}")
    print(f"  max|Im(H_{{down,down}})| = {dndn_imag:.2e}")

    if diag_diff < 0.1:
        print("  [PASS] Diagonal blocks are nearly equal (expected for m perp to z)")
    else:
        print("  [INFO] Diagonal blocks differ (may indicate m has z-component or SOC)")

if __name__ == '__main__':
    main()
