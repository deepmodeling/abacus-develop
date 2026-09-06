#!/usr/bin/env python3

import re
import struct
import sys


def read_exact(stream, size, description):
    payload = stream.read(size)
    if len(payload) != size:
        raise ValueError("missing or truncated {}".format(description))
    return payload


def read_native(stream, fmt, description):
    size = struct.calcsize(fmt)
    return struct.unpack(fmt, read_exact(stream, size, description))


def read_binary(filename, value_type):
    with open(filename, "rb") as stream:
        step, nbasis, nr = read_native(stream, "@iii", "record header")
        if step != 0:
            raise ValueError("expected zero-based ionic step 0, found {}".format(step))
        if nbasis <= 0 or nr < 0:
            raise ValueError("invalid matrix dimensions nbasis={} nr={}".format(nbasis, nr))

        blocks = []
        for block_index in range(nr):
            rx, ry, rz, nnz = read_native(
                stream, "@iiii", "R block {} header".format(block_index)
            )
            if nnz < 0:
                raise ValueError("negative nnz in R block {}".format(block_index))

            values = []
            if value_type == "complex":
                for value_index in range(nnz):
                    real, imag = read_native(
                        stream,
                        "@dd",
                        "R block {} complex value {}".format(block_index, value_index),
                    )
                    values.append(complex(real, imag))
            else:
                for value_index in range(nnz):
                    values.append(
                        read_native(
                            stream,
                            "@d",
                            "R block {} value {}".format(block_index, value_index),
                        )[0]
                    )

            columns = list(
                read_native(stream, "@{}i".format(nnz), "R block {} columns".format(block_index))
            ) if nnz else []
            row_ptr = list(
                read_native(
                    stream,
                    "@{}q".format(nbasis + 1),
                    "R block {} row pointers".format(block_index),
                )
            )
            blocks.append(((rx, ry, rz), values, columns, row_ptr))

        if stream.read(1):
            raise ValueError("unexpected trailing data after the first matrix record")
    return nbasis, blocks


def read_text_values(lines, index, marker, next_marker, value_type):
    if lines[index].strip() != marker:
        raise ValueError("expected {}, found {}".format(marker, lines[index].strip()))
    index += 1
    tokens = []
    while index < len(lines) and lines[index].strip() != next_marker:
        tokens.extend(lines[index].split())
        index += 1
    if value_type == "complex":
        pattern = re.compile(r"\(([^,]+),([^\)]+)\)")
        values = [complex(float(real), float(imag)) for real, imag in pattern.findall(" ".join(tokens))]
    else:
        values = [float(token) for token in tokens]
    return values, index


def read_text_reference(filename, value_type):
    with open(filename, "r", encoding="utf-8") as stream:
        lines = stream.readlines()

    nbasis = None
    csr_start = None
    dimension_pattern = re.compile(r"^\s*(\d+)\s+# number of localized basis\s*$")
    for index, line in enumerate(lines):
        match = dimension_pattern.match(line)
        if match:
            nbasis = int(match.group(1))
        if "#                               CSR Format" in line:
            csr_start = index
    if nbasis is None or csr_start is None:
        raise ValueError("reference does not contain an H(R)/S(R) CSR header")

    block_pattern = re.compile(r"^\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(\d+)\s*$")
    blocks = []
    index = csr_start + 1
    while index < len(lines):
        match = block_pattern.match(lines[index])
        if not match:
            index += 1
            continue
        rx, ry, rz, nnz = (int(value) for value in match.groups())
        index += 1
        values, index = read_text_values(
            lines, index, "# CSR values", "# CSR column indices", value_type
        )
        index += 1
        columns = []
        while index < len(lines) and lines[index].strip() != "# CSR row pointers":
            columns.extend(int(token) for token in lines[index].split())
            index += 1
        if index >= len(lines):
            raise ValueError("missing CSR row pointers for R = ({}, {}, {})".format(rx, ry, rz))
        index += 1
        row_ptr = []
        while index < len(lines) and len(row_ptr) < nbasis + 1:
            row_ptr.extend(int(token) for token in lines[index].split())
            index += 1
        if len(values) != nnz or len(columns) != nnz or len(row_ptr) != nbasis + 1:
            raise ValueError("invalid CSR lengths for R = ({}, {}, {})".format(rx, ry, rz))
        blocks.append(((rx, ry, rz), values, columns, row_ptr))
    return nbasis, blocks


def compare(binary_filename, reference_filename, value_type, tolerance):
    binary_nbasis, binary_blocks = read_binary(binary_filename, value_type)
    reference_nbasis, reference_blocks = read_text_reference(reference_filename, value_type)
    if binary_nbasis != reference_nbasis:
        raise ValueError(
            "matrix dimension differs: binary={} reference={}".format(binary_nbasis, reference_nbasis)
        )
    if len(binary_blocks) != len(reference_blocks):
        raise ValueError(
            "R block count differs: binary={} reference={}".format(
                len(binary_blocks), len(reference_blocks)
            )
        )

    for block_index, (binary_block, reference_block) in enumerate(
        zip(binary_blocks, reference_blocks)
    ):
        binary_r, binary_values, binary_columns, binary_row_ptr = binary_block
        reference_r, reference_values, reference_columns, reference_row_ptr = reference_block
        if binary_r != reference_r:
            raise ValueError(
                "R block {} differs: binary={} reference={}".format(
                    block_index, binary_r, reference_r
                )
            )
        if binary_columns != reference_columns or binary_row_ptr != reference_row_ptr:
            raise ValueError("CSR indices differ for R = {}".format(binary_r))
        if len(binary_values) != len(reference_values):
            raise ValueError("value count differs for R = {}".format(binary_r))
        for value_index, (actual, expected) in enumerate(zip(binary_values, reference_values)):
            if abs(actual - expected) > tolerance:
                raise ValueError(
                    "value {} differs for R = {}: binary={} reference={} tolerance={}".format(
                        value_index, binary_r, actual, expected, tolerance
                    )
                )


def main():
    if len(sys.argv) != 5 or sys.argv[3] not in ("real", "complex"):
        print("usage: compare_hsr_binary.py BINARY REFERENCE real|complex ACCURACY")
        return 2
    binary_filename, reference_filename, value_type = sys.argv[1:4]
    tolerance = 10.0 ** (-int(sys.argv[4]))
    try:
        compare(binary_filename, reference_filename, value_type, tolerance)
    except (OSError, ValueError, struct.error) as error:
        print("failed to compare H(R)/S(R) output: {}".format(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
