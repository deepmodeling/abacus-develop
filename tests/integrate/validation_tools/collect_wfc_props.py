"""Collect band norms and phase-invariant fingerprints from the first PW text frame."""

import argparse
import math
from pathlib import Path


FINGERPRINT_SKETCHES = 8
UINT64_MASK = (1 << 64) - 1
SPLITMIX_GOLDEN = 0x9E3779B97F4A7C15
SPLITMIX_MIX1 = 0xBF58476D1CE4E5B9
SPLITMIX_MIX2 = 0x94D049BB133111EB
FINGERPRINT_WEIGHTS = (1.0 + 0.0j, -1.0 + 0.0j, 0.0 + 1.0j, 0.0 - 1.0j)


def splitmix64(value):
    value = (value + SPLITMIX_GOLDEN) & UINT64_MASK
    value = ((value ^ (value >> 30)) * SPLITMIX_MIX1) & UINT64_MASK
    value = ((value ^ (value >> 27)) * SPLITMIX_MIX2) & UINT64_MASK
    return value ^ (value >> 31)


def wfc_projections(components):
    """Return the raw squared norm and the phase-invariant projection products."""
    norm2 = sum(
        value.real * value.real + value.imag * value.imag
        for component in components
        for value in component
    )
    if not math.isfinite(norm2) or norm2 == 0.0:
        raise ValueError("wavefunction has zero or non-finite norm")

    projections = [0.0j] * FINGERPRINT_SKETCHES
    flat_index = 0
    for component in components:
        for value in component:
            for sketch in range(FINGERPRINT_SKETCHES):
                weight_index = splitmix64(flat_index + sketch * SPLITMIX_GOLDEN) & 3
                projections[sketch] += FINGERPRINT_WEIGHTS[weight_index] * value
            flat_index += 1

    norm = math.sqrt(norm2)
    normalized = [projection / norm for projection in projections]
    metrics = [(f"power_{j}", abs(q) ** 2) for j, q in enumerate(normalized)]
    for j, q in enumerate(normalized):
        cross = q * normalized[(j + 1) % FINGERPRINT_SKETCHES].conjugate()
        metrics.append((f"cross_{j}_re", cross.real))
        metrics.append((f"cross_{j}_im", cross.imag))
    return norm2, metrics


def read_block(lines, tag):
    """Read an ABACUS text block, whose opening and closing tags are identical."""
    if next(lines, None) != tag:
        raise ValueError(f"expected {tag}")
    rows = []
    for line in lines:
        if line == tag:
            return rows
        rows.append(line)
    raise ValueError(f"missing closing {tag}")


def pw_fingerprint(stream):
    lines = (line.strip() for line in stream if line.strip())
    header = next(lines, "")
    columns = header.split()
    values = next(lines, "").split()
    if columns[6:8] != ["ngtot", "nband"] or len(values) != len(columns):
        raise ValueError("invalid PW wavefunction header")
    ng, nbands = map(int, values[6:8])
    if ng <= 0 or nbands <= 0:
        raise ValueError(f"invalid dimensions: ngtot={ng}, nband={nbands}")

    read_block(lines, "<Reciprocal Lattice Vector>")
    g_vectors = []
    for row in read_block(lines, "<G vectors>"):
        gx, gy, gz = map(int, row.split())
        g_vectors.append((gx, gy, gz))
    if len(g_vectors) != ng or len(set(g_vectors)) != ng:
        raise ValueError(f"expected {ng} distinct G vectors")
    # MPI ranks append their G vectors and coefficients in the same order.
    order = sorted(range(ng), key=lambda i: g_vectors[i])
    metrics = [("pw_ngtot", ng), ("pw_nband", nbands)]
    for band in range(1, nbands + 1):
        rows = read_block(lines, f"< Band {band} >")
        fields = [float(value) for row in rows for value in row.split()]
        if len(fields) != 2 * ng:
            raise ValueError(f"band {band}: expected {ng} real/imaginary pairs")
        coeffs = [complex(fields[2 * i], fields[2 * i + 1]) for i in order]
        norm2, products = wfc_projections([coeffs])
        prefix = f"pw_band{band}_wfc_fp_"
        metrics.append((prefix + "norm2", norm2))
        metrics.extend((prefix + key, value) for key, value in products)

    # Appended ionic steps start with a new header.
    if next(lines, header) != header:
        raise ValueError(f"unexpected data after band {nbands}")
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    with args.path.open(encoding="utf-8") as stream:
        for key, value in pw_fingerprint(stream):
            print(f"{key} {value:.12g}")
