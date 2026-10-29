#!/usr/bin/env python3
"""Compare complete Wannier90 AMN/MMN/EIG bundles in one consistent gauge.

For each k point, one unitary W aligns the candidate AMN matrix to the
reference. The same W must also align every MMN link and the diagonal energy
operator: A_ref = W A, M_ref(k,j) = W(k) M(k,j) W(j)^H, and
diag(E_ref) = W diag(E) W^H. Every matrix entry and every labelled eigenvalue
is checked. The energy-operator check bounds any mixing of distinct roots
directly in eV, without an extra degeneracy window. No invariant
summary, independent link fit, or cluster-average energy replaces those tests.

The projection block must have full row rank to determine its gauge. This
limited comparator fails explicitly for an underdetermined projection set.
AMN/MMN tolerances are dimensionless, and EIG tolerances are in eV.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np


class ComparisonError(ValueError):
    """Malformed input or unsupported gauge alignment."""


def _lines(path):
    path = Path(path).resolve()
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise ComparisonError("cannot read {}: {}".format(path, error)) from error
    print("Read {}".format(path))
    return [(number, line.split()) for number, line in enumerate(text.splitlines(), 1)
            if line.strip()]


def _integers(fields, count, context):
    if len(fields) != count:
        raise ComparisonError("{}: expected {} integer fields".format(context, count))
    try:
        return tuple(int(field) for field in fields)
    except ValueError as error:
        raise ComparisonError("{}: invalid integer field".format(context)) from error


def _complex(fields, context):
    if len(fields) != 2:
        raise ComparisonError("{}: expected real and imaginary parts".format(context))
    try:
        value = complex(float(fields[0]), float(fields[1]))
    except ValueError as error:
        raise ComparisonError("{}: invalid complex value".format(context)) from error
    if not np.isfinite(value):
        raise ComparisonError("{}: nonfinite complex value".format(context))
    return value


def _header(path):
    lines = _lines(path)
    if not lines:
        raise ComparisonError("{}: empty file".format(path))
    # Wannier90 permits one descriptive line before the dimension header.
    # Existing ABACUS references omit that line. Do not scan past bad data.
    first = lines[0][1]
    try:
        dims = _integers(first, 3, path)
    except ComparisonError:
        if len(lines) < 2:
            raise ComparisonError("{}: missing dimension header".format(path))
        numeric_start = False
        try:
            float(first[0])
            numeric_start = True
        except ValueError:
            pass
        first_number = lines[0][0]
        if first_number != 1 or numeric_start:
            raise ComparisonError("{}: malformed dimension header".format(path))
        lines = lines[1:]
        dims = _integers(lines[0][1], 3, path)
    if min(dims) <= 0:
        raise ComparisonError("{}: dimensions must be positive".format(path))
    return dims, lines[1:]


def read_amn(path):
    dims, lines = _header(path)
    bands, kpoints, projections = dims
    expected = bands * kpoints * projections
    if len(lines) != expected:
        raise ComparisonError("{}: expected {} AMN records, got {}".format(path, expected, len(lines)))
    matrices = np.empty((kpoints, bands, projections), dtype=complex)
    seen = set()
    for number, fields in lines:
        context = "{}:{}".format(path, number)
        if len(fields) != 5:
            raise ComparisonError("{}: AMN record needs five fields".format(context))
        band, projection, kpoint = _integers(fields[:3], 3, context)
        if not (1 <= band <= bands and 1 <= projection <= projections and 1 <= kpoint <= kpoints):
            raise ComparisonError("{}: AMN index outside dimensions".format(context))
        key = (band, projection, kpoint)
        if key in seen:
            raise ComparisonError("{}: duplicate AMN index".format(context))
        seen.add(key)
        matrices[kpoint - 1, band - 1, projection - 1] = _complex(fields[3:], context)
    return dims, matrices


def read_mmn(path):
    dims, lines = _header(path)
    bands, kpoints, neighbours = dims
    block_size = 1 + bands * bands
    expected = kpoints * neighbours * block_size
    if len(lines) != expected:
        raise ComparisonError("{}: expected {} MMN records, got {}".format(path, expected, len(lines)))
    matrices = {}
    counts = [0] * kpoints
    for start in range(0, len(lines), block_size):
        number, fields = lines[start]
        context = "{}:{}".format(path, number)
        link = _integers(fields, 5, context)
        if not (1 <= link[0] <= kpoints and 1 <= link[1] <= kpoints):
            raise ComparisonError("{}: MMN endpoint outside dimensions".format(context))
        if link in matrices:
            raise ComparisonError("{}: duplicate MMN link".format(context))
        counts[link[0] - 1] += 1
        values = [_complex(fields, "{}:{}".format(path, number))
                  for number, fields in lines[start + 1:start + block_size]]
        # Wannier90 writes m inside n: matrix entries are column-major.
        matrices[link] = np.asarray(values).reshape((bands, bands), order="F")
    if any(count != neighbours for count in counts):
        raise ComparisonError("{}: wrong number of neighbours for a k point".format(path))
    return dims, matrices


def read_eig(path, bands, kpoints):
    lines = _lines(path)
    if len(lines) != bands * kpoints:
        raise ComparisonError("{}: incorrect EIG record count".format(path))
    values = np.empty((kpoints, bands))
    seen = set()
    for number, fields in lines:
        context = "{}:{}".format(path, number)
        if len(fields) != 3:
            raise ComparisonError("{}: EIG record needs three fields".format(context))
        band, kpoint = _integers(fields[:2], 2, context)
        if not (1 <= band <= bands and 1 <= kpoint <= kpoints):
            raise ComparisonError("{}: EIG index outside dimensions".format(context))
        key = (band, kpoint)
        if key in seen:
            raise ComparisonError("{}: duplicate EIG index".format(context))
        seen.add(key)
        try:
            value = float(fields[2])
        except ValueError as error:
            raise ComparisonError("{}: invalid EIG value".format(context)) from error
        if not np.isfinite(value):
            raise ComparisonError("{}: nonfinite EIG value".format(context))
        values[kpoint - 1, band - 1] = value
    if np.any(np.diff(values, axis=1) < 0):
        raise ComparisonError("{}: eigenvalues are not ordered by band".format(path))
    return values


def read_bundle(prefix):
    prefix = str(prefix)
    adims, amn = read_amn(prefix + ".amn")
    mdims, mmn = read_mmn(prefix + ".mmn")
    if adims[:2] != mdims[:2]:
        raise ComparisonError("AMN and MMN dimensions disagree")
    eig = read_eig(prefix + ".eig", adims[0], adims[1])
    return {"amn_dims": adims, "mmn_dims": mdims, "amn": amn, "mmn": mmn, "eig": eig}


def _full_row_rank(matrix):
    if matrix.shape[0] > matrix.shape[1]:
        return False
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank_floor = np.finfo(float).eps * max(matrix.shape) * singular[0]
    return singular[-1] > rank_floor


def _gauge(reference, candidate):
    kpoints, bands, _ = reference["amn"].shape
    gauges = []
    for kpoint in range(kpoints):
        a = reference["amn"][kpoint]
        b = candidate["amn"][kpoint]
        if not _full_row_rank(a) or not _full_row_rank(b):
            raise ComparisonError("k point {}: projection block does not determine its gauge".format(kpoint + 1))
        left, _, right = np.linalg.svd(a @ b.conj().T)
        gauge = left @ right
        error = np.max(np.abs(gauge @ gauge.conj().T - np.eye(bands)))
        if not np.isfinite(error) or error > 64 * np.finfo(float).eps * bands:
            raise ComparisonError("computed alignment is not unitary")
        gauges.append(gauge)
    return gauges


def _maximum_error(actual, expected, label):
    # Finite input records can overflow during products or subtraction.
    # Check each local result before max(): max(0.0, nan) can return 0.0.
    with np.errstate(over="ignore", invalid="ignore"):
        errors = np.abs(actual - expected)
    if not np.isfinite(errors).all():
        raise ComparisonError("nonfinite arithmetic in {} comparison".format(label))
    return float(np.max(errors))


def compare(reference, candidate, matrix_tolerance, energy_tolerance):
    for tolerance in (matrix_tolerance, energy_tolerance):
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ComparisonError("tolerances must be finite and positive")
    if (reference["amn_dims"] != candidate["amn_dims"]
            or reference["mmn_dims"] != candidate["mmn_dims"]):
        raise ComparisonError("bundle dimensions differ")
    if reference["mmn"].keys() != candidate["mmn"].keys():
        raise ComparisonError("MMN link endpoints or periodic translations differ")
    labelled_error = _maximum_error(reference["eig"], candidate["eig"], "EIG")
    gauges = _gauge(reference, candidate)
    amn_error = 0.0
    mmn_error = 0.0
    operator_error = 0.0
    for kpoint, gauge in enumerate(gauges):
        amn_local_error = _maximum_error(gauge @ candidate["amn"][kpoint],
                                        reference["amn"][kpoint], "AMN")
        amn_error = max(amn_error, amn_local_error)
        energy_operator = gauge @ np.diag(candidate["eig"][kpoint]) @ gauge.conj().T
        operator_local_error = _maximum_error(energy_operator,
                                             np.diag(reference["eig"][kpoint]), "energy operator")
        operator_error = max(operator_error, operator_local_error)
    for link, matrix in reference["mmn"].items():
        left = gauges[link[0] - 1]
        right = gauges[link[1] - 1]
        aligned = left @ candidate["mmn"][link] @ right.conj().T
        mmn_local_error = _maximum_error(aligned, matrix, "MMN")
        mmn_error = max(mmn_error, mmn_local_error)
    metrics = {"amn_max_abs_error": amn_error, "mmn_max_abs_error": mmn_error,
               "eig_max_abs_error_eV": labelled_error,
               "energy_operator_max_abs_error_eV": operator_error}
    if not all(np.isfinite(value) for value in metrics.values()):
        raise ComparisonError("nonfinite comparison result")
    metrics["passed"] = (amn_error <= matrix_tolerance and mmn_error <= matrix_tolerance
                         and labelled_error <= energy_tolerance and operator_error <= energy_tolerance)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_prefix")
    parser.add_argument("candidate_prefix")
    parser.add_argument("--matrix-tolerance", type=float, required=True)
    parser.add_argument("--energy-tolerance", type=float, required=True)
    args = parser.parse_args()
    try:
        reference = read_bundle(args.reference_prefix)
        candidate = read_bundle(args.candidate_prefix)
        result = compare(reference, candidate, args.matrix_tolerance, args.energy_tolerance)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (ComparisonError, np.linalg.LinAlgError) as error:
        print("Wannier90 comparison failed: {}".format(error), file=sys.stderr)
        return 1
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
