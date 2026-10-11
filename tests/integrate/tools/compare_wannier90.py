#!/usr/bin/env python3
"""Compare Wannier90 interface files without depending on Bloch-state gauge.

AMN Gram matrices, MMN singular values, and MMN Wilson-loop traces are
dimensionless. EIG values and their absolute tolerance are in eV, matching the
Wannier90 EIG file format.
"""

import argparse
import math
import pathlib
import sys


EIG_DEGENERATE_WINDOW_EV = 3.0e-6
EIG_CLUSTER_SHAPE_TOLERANCE_EV = 1.0e-6


class ComparisonError(Exception):
    """Raised when a file is malformed or its physical invariants differ."""


def _read_lines(path):
    try:
        text = pathlib.Path(path).read_text(encoding="utf-8")
    except OSError as error:
        raise ComparisonError("cannot read {}: {}".format(path, error))
    return [(number, line.strip()) for number, line in enumerate(text.splitlines(), 1) if line.strip()]


def _matrix_header(path):
    lines = _read_lines(path)
    for index, (line_number, line) in enumerate(lines):
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            dimensions = tuple(int(field) for field in fields)
        except ValueError:
            continue
        if min(dimensions) <= 0:
            raise ComparisonError("{}:{}: dimensions must be positive".format(path, line_number))
        return dimensions, lines[index + 1 :]
    raise ComparisonError("{}: missing three-integer dimension header".format(path))


def _finite_float(path, line_number, value):
    try:
        parsed = float(value)
    except ValueError:
        raise ComparisonError("{}:{}: invalid floating-point value {!r}".format(path, line_number, value))
    if not math.isfinite(parsed):
        raise ComparisonError("{}:{}: non-finite floating-point value {!r}".format(path, line_number, value))
    return parsed


def _parse_amn(path):
    dimensions, lines = _matrix_header(path)
    band_count, kpoint_count, projection_count = dimensions
    expected_count = band_count * kpoint_count * projection_count
    if len(lines) != expected_count:
        raise ComparisonError(
            "{}: expected {} AMN entries, found {}".format(path, expected_count, len(lines))
        )

    matrices = {
        kpoint: [[None for _ in range(projection_count)] for _ in range(band_count)]
        for kpoint in range(1, kpoint_count + 1)
    }
    for line_number, line in lines:
        fields = line.split()
        if len(fields) != 5:
            raise ComparisonError("{}:{}: AMN entry must have five fields".format(path, line_number))
        try:
            band = int(fields[0])
            projection = int(fields[1])
            kpoint = int(fields[2])
        except ValueError:
            raise ComparisonError("{}:{}: invalid AMN index".format(path, line_number))
        if not (1 <= band <= band_count and 1 <= projection <= projection_count and 1 <= kpoint <= kpoint_count):
            raise ComparisonError("{}:{}: AMN index is outside the header dimensions".format(path, line_number))
        if matrices[kpoint][band - 1][projection - 1] is not None:
            raise ComparisonError("{}:{}: duplicate AMN index".format(path, line_number))
        real = _finite_float(path, line_number, fields[3])
        imag = _finite_float(path, line_number, fields[4])
        matrices[kpoint][band - 1][projection - 1] = complex(real, imag)
    return dimensions, matrices


def _parse_mmn(path):
    dimensions, lines = _matrix_header(path)
    band_count, kpoint_count, neighbour_count = dimensions
    matrices = {}
    position = 0
    expected_blocks = kpoint_count * neighbour_count
    for _ in range(expected_blocks):
        if position >= len(lines):
            raise ComparisonError("{}: missing MMN link header".format(path))
        line_number, line = lines[position]
        position += 1
        fields = line.split()
        if len(fields) != 5:
            raise ComparisonError("{}:{}: MMN link header must have five integers".format(path, line_number))
        try:
            link = tuple(int(field) for field in fields)
        except ValueError:
            raise ComparisonError("{}:{}: invalid MMN link header".format(path, line_number))
        if not (1 <= link[0] <= kpoint_count and 1 <= link[1] <= kpoint_count):
            raise ComparisonError("{}:{}: MMN k-point index is outside the header dimensions".format(path, line_number))
        if link in matrices:
            raise ComparisonError("{}:{}: duplicate MMN link {}".format(path, line_number, link))

        entries = []
        for _ in range(band_count * band_count):
            if position >= len(lines):
                raise ComparisonError("{}: incomplete MMN matrix for link {}".format(path, link))
            value_line_number, value_line = lines[position]
            position += 1
            value_fields = value_line.split()
            if len(value_fields) != 2:
                raise ComparisonError(
                    "{}:{}: MMN matrix entry must have real and imaginary parts".format(path, value_line_number)
                )
            real = _finite_float(path, value_line_number, value_fields[0])
            imag = _finite_float(path, value_line_number, value_fields[1])
            entries.append(complex(real, imag))
        matrices[link] = [
            entries[row * band_count : (row + 1) * band_count] for row in range(band_count)
        ]

    if position != len(lines):
        line_number, _ = lines[position]
        raise ComparisonError("{}:{}: unexpected data after the MMN matrices".format(path, line_number))
    return dimensions, matrices


def _parse_eig(path):
    values = {}
    for line_number, line in _read_lines(path):
        fields = line.split()
        if len(fields) != 3:
            raise ComparisonError("{}:{}: EIG entry must have three fields".format(path, line_number))
        try:
            band = int(fields[0])
            kpoint = int(fields[1])
        except ValueError:
            raise ComparisonError("{}:{}: invalid EIG index".format(path, line_number))
        if band <= 0 or kpoint <= 0:
            raise ComparisonError("{}:{}: EIG indices must be positive".format(path, line_number))
        key = (band, kpoint)
        if key in values:
            raise ComparisonError("{}:{}: duplicate EIG index {}".format(path, line_number, key))
        values[key] = _finite_float(path, line_number, fields[2])
    if not values:
        raise ComparisonError("{}: EIG file has no entries".format(path))
    return values


def _gram(matrix):
    row_count = len(matrix)
    column_count = len(matrix[0])
    return [
        [
            sum(matrix[row][left].conjugate() * matrix[row][right] for row in range(row_count))
            for right in range(column_count)
        ]
        for left in range(column_count)
    ]


def _symmetric_eigenvalues(matrix):
    size = len(matrix)
    values = [row[:] for row in matrix]
    scale = max(1.0, max(abs(values[row][row]) for row in range(size)))
    threshold = 1.0e-14 * scale
    max_sweeps = max(30, 5 * size)

    for _ in range(max_sweeps):
        max_off_diagonal = max(
            (abs(values[row][column]) for row in range(size) for column in range(row + 1, size)),
            default=0.0,
        )
        if max_off_diagonal <= threshold:
            return sorted(values[index][index] for index in range(size))

        for left in range(size - 1):
            for right in range(left + 1, size):
                off_diagonal = values[left][right]
                if abs(off_diagonal) <= threshold:
                    continue
                left_diagonal = values[left][left]
                right_diagonal = values[right][right]
                ratio = (right_diagonal - left_diagonal) / (2.0 * off_diagonal)
                sign = 1.0 if ratio >= 0.0 else -1.0
                tangent = sign / (abs(ratio) + math.hypot(1.0, ratio))
                cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
                sine = tangent * cosine

                for index in range(size):
                    if index == left or index == right:
                        continue
                    old_left = values[index][left]
                    old_right = values[index][right]
                    new_left = cosine * old_left - sine * old_right
                    new_right = sine * old_left + cosine * old_right
                    values[index][left] = new_left
                    values[left][index] = new_left
                    values[index][right] = new_right
                    values[right][index] = new_right
                values[left][left] = left_diagonal - tangent * off_diagonal
                values[right][right] = right_diagonal + tangent * off_diagonal
                values[left][right] = 0.0
                values[right][left] = 0.0

    raise ComparisonError("Jacobi eigensolver did not converge")


def _singular_values(matrix):
    gram = _gram(matrix)
    size = len(gram)
    real_form = [[0.0 for _ in range(2 * size)] for _ in range(2 * size)]
    for row in range(size):
        for column in range(size):
            value = gram[row][column]
            real_form[row][column] = value.real
            real_form[row][column + size] = -value.imag
            real_form[row + size][column] = value.imag
            real_form[row + size][column + size] = value.real

    duplicated_eigenvalues = _symmetric_eigenvalues(real_form)
    eigenvalues = [
        0.5 * (duplicated_eigenvalues[2 * index] + duplicated_eigenvalues[2 * index + 1])
        for index in range(size)
    ]
    negative_tolerance = 1.0e-12 * max(1.0, max(abs(value) for value in eigenvalues))
    if min(eigenvalues) < -negative_tolerance:
        raise ComparisonError("M^H M has a negative eigenvalue outside round-off tolerance")
    return [math.sqrt(max(0.0, value)) for value in eigenvalues]


def _matrix_multiply(left, right):
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def _transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def _contractible_four_cycles(links):
    adjacency = {}
    for link in links:
        adjacency.setdefault(link[0], []).append(link)
    for outgoing in adjacency.values():
        outgoing.sort()

    cycles = set()
    for first in sorted(links):
        for second in adjacency.get(first[1], []):
            for third in adjacency.get(second[1], []):
                for fourth in adjacency.get(third[1], []):
                    cycle = (first, second, third, fourth)
                    if fourth[1] != first[0]:
                        continue
                    if any(sum(link[component] for link in cycle) != 0 for component in range(2, 5)):
                        continue
                    rotations = [cycle[index:] + cycle[:index] for index in range(4)]
                    cycles.add(min(rotations))
    return sorted(cycles)


def _wilson_trace_powers(matrices, cycle, band_count):
    # Wannier90 writes M(m,n) with n as the outer loop, so parsed rows are M^T.
    wilson_product = _transpose(matrices[cycle[0]])
    for link in cycle[1:]:
        wilson_product = _matrix_multiply(wilson_product, _transpose(matrices[link]))

    traces = []
    power = wilson_product
    for _ in range(band_count):
        traces.append(sum(power[index][index] for index in range(band_count)) / band_count)
        power = _matrix_multiply(power, wilson_product)
    return traces


def _compare_amn(reference, calculated, tolerance):
    reference_dimensions, reference_matrices = _parse_amn(reference)
    calculated_dimensions, calculated_matrices = _parse_amn(calculated)
    if reference_dimensions != calculated_dimensions:
        raise ComparisonError(
            "AMN dimensions differ: {} != {}".format(reference_dimensions, calculated_dimensions)
        )

    maximum = (0.0, None)
    for kpoint in sorted(reference_matrices):
        reference_gram = _gram(reference_matrices[kpoint])
        calculated_gram = _gram(calculated_matrices[kpoint])
        for row in range(len(reference_gram)):
            for column in range(len(reference_gram)):
                difference = abs(reference_gram[row][column] - calculated_gram[row][column])
                if difference > maximum[0]:
                    maximum = (difference, (kpoint, row + 1, column + 1))
    if maximum[0] > tolerance:
        raise ComparisonError(
            "AMN Gram matrices differ by {:.6e} at (k, projection, projection)={}, tolerance {:.1e}".format(
                maximum[0], maximum[1], tolerance
            )
        )
    return maximum[0]


def _compare_mmn(reference, calculated, tolerance):
    reference_dimensions, reference_matrices = _parse_mmn(reference)
    calculated_dimensions, calculated_matrices = _parse_mmn(calculated)
    if reference_dimensions != calculated_dimensions:
        raise ComparisonError(
            "MMN dimensions differ: {} != {}".format(reference_dimensions, calculated_dimensions)
        )
    if set(reference_matrices) != set(calculated_matrices):
        missing = sorted(set(reference_matrices) - set(calculated_matrices))
        extra = sorted(set(calculated_matrices) - set(reference_matrices))
        raise ComparisonError("MMN links differ; missing={}, extra={}".format(missing, extra))

    maximum = (0.0, None)
    for link in sorted(reference_matrices):
        reference_singular_values = _singular_values(reference_matrices[link])
        calculated_singular_values = _singular_values(calculated_matrices[link])
        for index, (reference_value, calculated_value) in enumerate(
            zip(reference_singular_values, calculated_singular_values), 1
        ):
            difference = abs(reference_value - calculated_value)
            if difference > maximum[0]:
                maximum = (difference, (link, index))
    if maximum[0] > tolerance:
        raise ComparisonError(
            "MMN singular values differ by {:.6e} at (link, index)={}, tolerance {:.1e}".format(
                maximum[0], maximum[1], tolerance
            )
        )

    band_count = reference_dimensions[0]
    wilson_maximum = (0.0, None)
    for cycle in _contractible_four_cycles(reference_matrices):
        reference_traces = _wilson_trace_powers(reference_matrices, cycle, band_count)
        calculated_traces = _wilson_trace_powers(calculated_matrices, cycle, band_count)
        for exponent, (reference_trace, calculated_trace) in enumerate(
            zip(reference_traces, calculated_traces), 1
        ):
            difference = abs(reference_trace - calculated_trace)
            if difference > wilson_maximum[0]:
                wilson_maximum = (difference, (cycle, exponent))
    if wilson_maximum[0] > tolerance:
        raise ComparisonError(
            "MMN Wilson-loop trace powers differ by {:.6e} at (cycle, power)={}, tolerance {:.1e}".format(
                wilson_maximum[0], wilson_maximum[1], tolerance
            )
        )
    return max(maximum[0], wilson_maximum[0])


def _compare_eig(reference, calculated, tolerance):
    reference_values = _parse_eig(reference)
    calculated_values = _parse_eig(calculated)
    if set(reference_values) != set(calculated_values):
        missing = sorted(set(reference_values) - set(calculated_values))
        extra = sorted(set(calculated_values) - set(reference_values))
        raise ComparisonError("EIG indices differ; missing={}, extra={}".format(missing, extra))

    keys_by_kpoint = {}
    for band, kpoint in reference_values:
        keys_by_kpoint.setdefault(kpoint, []).append((band, kpoint))

    maximum = 0.0
    maximum_location = None
    for kpoint in sorted(keys_by_kpoint):
        keys = sorted(keys_by_kpoint[kpoint])
        clusters = []
        cluster = []
        for key in keys:
            if cluster:
                first = cluster[0]
                span = abs(reference_values[key] - reference_values[first])
                if span > EIG_DEGENERATE_WINDOW_EV:
                    clusters.append(cluster)
                    cluster = []
            cluster.append(key)
        clusters.append(cluster)

        for cluster in clusters:
            if len(cluster) == 1:
                key = cluster[0]
                difference = abs(reference_values[key] - calculated_values[key])
                location = key
            else:
                reference_center = sum(reference_values[key] for key in cluster) / len(cluster)
                calculated_center = sum(calculated_values[key] for key in cluster) / len(cluster)
                difference = abs(reference_center - calculated_center)
                location = (cluster[0][0], cluster[-1][0], kpoint)

                shape_difference = max(
                    abs(
                        (reference_values[key] - reference_center)
                        - (calculated_values[key] - calculated_center)
                    )
                    for key in cluster
                )
                if shape_difference > EIG_CLUSTER_SHAPE_TOLERANCE_EV:
                    raise ComparisonError(
                        "EIG near-degenerate cluster {} shape differs by {:.6e} eV; "
                        "tolerance {:.1e} eV".format(
                            location, shape_difference, EIG_CLUSTER_SHAPE_TOLERANCE_EV
                        )
                    )

            if difference > maximum:
                maximum = difference
                maximum_location = location

    if maximum > tolerance:
        raise ComparisonError(
            "EIG invariant values differ by {:.6e} at {}, tolerance {:.1e}".format(
                maximum, maximum_location, tolerance
            )
        )
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file_type", choices=("amn", "mmn", "eig"))
    parser.add_argument("reference")
    parser.add_argument("calculated")
    parser.add_argument(
        "tolerance",
        type=float,
        help="absolute tolerance (dimensionless for AMN/MMN, eV for EIG)",
    )
    arguments = parser.parse_args()
    if not math.isfinite(arguments.tolerance) or arguments.tolerance <= 0.0:
        parser.error("tolerance must be positive and finite")

    compare = {"amn": _compare_amn, "mmn": _compare_mmn, "eig": _compare_eig}[arguments.file_type]
    try:
        maximum = compare(arguments.reference, arguments.calculated, arguments.tolerance)
    except ComparisonError as error:
        print("Compare Error: {}".format(error))
        return 1
    print(
        "{} comparison passed: maximum invariant difference {:.6e} (tolerance {:.1e})".format(
            arguments.file_type.upper(), maximum, arguments.tolerance
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
