#!/usr/bin/env python3
"""Validate and compare the 100-point LR absorption tables."""
import argparse
import math
from pathlib import Path


HEADER = "Frequency (eV) | wave length(nm) | Absorption (a.u.)"
SPECTRUM_POINTS = 100


def read_spectrum(path):
    lines = path.read_text().splitlines()
    if not lines or lines[0].strip() != HEADER:
        raise ValueError(f"{path}: missing absorption table header")
    if len(lines) - 1 != SPECTRUM_POINTS:
        raise ValueError(f"{path}: expected {SPECTRUM_POINTS} data rows, found {len(lines) - 1}")

    values = []
    for line_number, line in enumerate(lines[1:], start=2):
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


def compare_spectra(reference, output, tolerance=1e-4):
    reference_rows = read_spectrum(reference)
    output_rows = read_spectrum(output)
    for line, (expected, actual) in enumerate(zip(reference_rows, output_rows), start=2):
        for column, (expected_value, actual_value) in enumerate(zip(expected, actual), start=1):
            if abs(expected_value - actual_value) > tolerance:
                raise ValueError(
                    f"{output}:{line}: column {column} is {actual_value}, "
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
