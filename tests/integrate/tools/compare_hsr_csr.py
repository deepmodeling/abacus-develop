#!/usr/bin/env python3

import math
import sys

from compare_hsr_binary import read_text_reference


def csr_to_entries(nbasis, values, columns, row_ptr, rvec):
    if len(columns) != len(values):
        raise ValueError(
            "CSR column/value count differs for R = {}: columns={} values={}".format(
                rvec, len(columns), len(values)
            )
        )

    if len(row_ptr) != nbasis + 1:
        raise ValueError(
            "CSR row pointer count differs for R = {}: expected={} actual={}".format(
                rvec, nbasis + 1, len(row_ptr)
            )
        )

    if row_ptr[0] != 0 or row_ptr[-1] != len(values):
        raise ValueError(
            "invalid CSR row pointers for R = {}".format(rvec)
        )

    if any(ptr < 0 or ptr > len(values) for ptr in row_ptr):
        raise ValueError(
            "CSR row pointer out of range for R = {}".format(rvec)
        )

    if any(lhs > rhs for lhs, rhs in zip(row_ptr, row_ptr[1:])):
        raise ValueError(
            "non-monotonic CSR row pointers for R = {}".format(rvec)
        )

    entries = {}

    for row in range(nbasis):
        begin = row_ptr[row]
        end = row_ptr[row + 1]

        for index in range(begin, end):
            col = columns[index]

            if col < 0 or col >= nbasis:
                raise ValueError(
                    "CSR column {} out of range for R = {}".format(col, rvec)
                )

            key = (row, col)

            if key in entries:
                raise ValueError(
                    "duplicate CSR entry {} for R = {}".format(key, rvec)
                )

            entries[key] = values[index]

    return entries


def blocks_to_map(nbasis, blocks):
    result = {}

    for rvec, values, columns, row_ptr in blocks:
        if rvec in result:
            raise ValueError("duplicate R block: {}".format(rvec))

        result[rvec] = csr_to_entries(
            nbasis, values, columns, row_ptr, rvec
        )

    return result


def read_pair(reference_filename, actual_filename):
    errors = []

    # Prefer real because the common nspin=1/2 H(R), S(R) outputs are real.
    # Fall back to complex for spinor matrices.
    for value_type in ("real", "complex"):
        try:
            reference_nbasis, reference_blocks = read_text_reference(
                reference_filename, value_type
            )
            actual_nbasis, actual_blocks = read_text_reference(
                actual_filename, value_type
            )

            return (
                value_type,
                reference_nbasis,
                reference_blocks,
                actual_nbasis,
                actual_blocks,
            )
        except (OSError, ValueError) as error:
            errors.append("{}: {}".format(value_type, error))

    raise ValueError(
        "failed to parse files as real or complex CSR: {}".format(
            "; ".join(errors)
        )
    )


def compare(reference_filename, actual_filename, tolerance):
    (
        value_type,
        reference_nbasis,
        reference_blocks,
        actual_nbasis,
        actual_blocks,
    ) = read_pair(reference_filename, actual_filename)

    if reference_nbasis != actual_nbasis:
        raise ValueError(
            "matrix dimension differs: reference={} actual={}".format(
                reference_nbasis, actual_nbasis
            )
        )

    reference = blocks_to_map(reference_nbasis, reference_blocks)
    actual = blocks_to_map(actual_nbasis, actual_blocks)

    max_difference = 0.0
    max_location = None

    # Missing R blocks and missing CSR entries are mathematically zero.
    for rvec in sorted(set(reference) | set(actual)):
        reference_entries = reference.get(rvec, {})
        actual_entries = actual.get(rvec, {})

        for row_col in sorted(set(reference_entries) | set(actual_entries)):
            expected = reference_entries.get(row_col, 0.0)
            calculated = actual_entries.get(row_col, 0.0)

            if not math.isfinite(abs(expected)) or not math.isfinite(abs(calculated)):
                row, col = row_col
                raise ValueError(
                    "non-finite matrix value at R={}, row={}, col={}: "
                    "reference={} actual={}".format(
                        rvec, row, col, expected, calculated
                    )
                )

            difference = abs(calculated - expected)

            if difference > max_difference:
                max_difference = difference
                max_location = (rvec, row_col, expected, calculated)

            if difference > tolerance:
                row, col = row_col
                raise ValueError(
                    "matrix value differs at R={}, row={}, col={}: "
                    "reference={} actual={} difference={} tolerance={}".format(
                        rvec,
                        row,
                        col,
                        expected,
                        calculated,
                        difference,
                        tolerance,
                    )
                )

    if max_location is None:
        print(
            "H(R)/S(R) CSR comparison passed: "
            "type={} max_abs_difference=0".format(value_type)
        )
    else:
        rvec, (row, col), expected, calculated = max_location
        print(
            "H(R)/S(R) CSR comparison passed: "
            "type={} max_abs_difference={} "
            "at R={}, row={}, col={}, reference={}, actual={}".format(
                value_type,
                max_difference,
                rvec,
                row,
                col,
                expected,
                calculated,
            )
        )


def main():
    if len(sys.argv) != 4:
        print(
            "usage: compare_hsr_csr.py "
            "REFERENCE_CSR ACTUAL_CSR ACCURACY"
        )
        return 2

    reference_filename = sys.argv[1]
    actual_filename = sys.argv[2]

    try:
        tolerance = 10.0 ** (-int(sys.argv[3]))

        compare(
            reference_filename,
            actual_filename,
            tolerance,
        )
    except (OSError, ValueError) as error:
        print(
            "failed to compare H(R)/S(R) CSR output: {}".format(error)
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
