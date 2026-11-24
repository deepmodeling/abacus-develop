#!/usr/bin/env python3
"""Validate and compare the 100-point LR absorption tables."""
import argparse
import math
from pathlib import Path
import re


HEADER = "# Frequency (eV) | wave length(nm) | Absorption (a.u.)"
SPECTRUM_POINTS = 100
# Keep the grid on an absolute scale. A separate floor for the intensities
# permits roundoff near zero without masking small spectral features.
COLUMN_TOLERANCES = (
    (1e-4, 0.0),   # frequency (eV)
    (1e-4, 0.0),   # wavelength (nm)
    (1e-8, 1e-4),  # absorption
)


def read_spectrum(path, *, require_provenance=False):
    lines = path.read_text().splitlines()
    metadata_start = 0
    if lines and lines[0].startswith("# ABACUS version: "):
        if len(lines) < 2 or not re.fullmatch(r"# Written at \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", lines[1]):
            raise ValueError(f"{path}: missing or invalid write timestamp")
        metadata_start = 2
    elif require_provenance:
        raise ValueError(f"{path}: missing ABACUS version header")

    point_count = f"# {SPECTRUM_POINTS} # number of points"
    metadata = lines[metadata_start:metadata_start + 3]
    if len(metadata) != 3 or not re.fullmatch(r"# Ionic_Step [1-9][0-9]*", metadata[0]):
        raise ValueError(f"{path}: missing or invalid one-based ionic step")
    if metadata[1] != point_count:
        raise ValueError(f"{path}: expected header {point_count}")
    if metadata[2] != HEADER:
        raise ValueError(f"{path}: missing absorption table header")
    data_start = metadata_start + 3
    data_rows = lines[data_start:]
    if len(data_rows) != SPECTRUM_POINTS:
        raise ValueError(f"{path}: expected {SPECTRUM_POINTS} data rows, found {len(data_rows)}")

    values = []
    for line_number, line in enumerate(data_rows, start=data_start + 1):
        columns = line.split()
        if len(columns) != 3:
            raise ValueError(f"{path}:{line_number}: expected three numerical columns")
        try:
            row = tuple(float(value) for value in columns)
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid numerical value") from error
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"{path}:{line_number}: non-finite absorption table value")
        values.append(row)
    return values


def compare_spectra(reference, output):
    reference_rows = read_spectrum(reference)
    output_rows = read_spectrum(output, require_provenance=True)
    for line, (expected, actual) in enumerate(zip(reference_rows, output_rows), start=1):
        columns = zip(expected, actual, COLUMN_TOLERANCES)
        for column, (expected_value, actual_value, (absolute, relative)) in enumerate(columns, start=1):
            if not math.isclose(expected_value, actual_value, abs_tol=absolute, rel_tol=relative):
                raise ValueError(
                    f"{output}: data row {line}, column {column} is {actual_value}, "
                    f"expected {expected_value} from {reference}"
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        compare_spectra(args.reference, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
