#!/usr/bin/env python3

import re
import struct
import sys


def read_reference(filename, value_type):
    values = []
    complex_pattern = re.compile(r"\(([^,]+),([^\)]+)\)")
    with open(filename, "r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("Row "):
                continue
            if value_type == "complex":
                values.extend(complex(float(real), float(imag)) for real, imag in complex_pattern.findall(stripped))
            else:
                values.extend(float(token) for token in stripped.split())
    return values


def read_binary(filename, value_type):
    with open(filename, "rb") as stream:
        header = stream.read(struct.calcsize("@i"))
        if len(header) != struct.calcsize("@i"):
            raise ValueError("missing native int matrix dimension")
        dim = struct.unpack("@i", header)[0]
        count = dim * (dim + 1) // 2
        element_format = "@dd" if value_type == "complex" else "@d"
        element_size = struct.calcsize(element_format)
        payload = stream.read(count * element_size)
        if len(payload) != count * element_size:
            raise ValueError("binary payload is shorter than the declared upper triangle")
        if stream.read(1):
            raise ValueError("unexpected trailing data after the first matrix record")

    if value_type == "complex":
        return [complex(real, imag) for real, imag in struct.iter_unpack(element_format, payload)]
    return [value[0] for value in struct.iter_unpack(element_format, payload)]


def main():
    if len(sys.argv) != 5 or sys.argv[3] not in ("real", "complex"):
        print("usage: compare_hsk_binary.py BINARY REFERENCE real|complex ACCURACY")
        return 2

    binary_filename, reference_filename, value_type = sys.argv[1:4]
    tolerance = 10.0 ** (-int(sys.argv[4]))
    try:
        actual = read_binary(binary_filename, value_type)
        reference = read_reference(reference_filename, value_type)
    except (OSError, ValueError, struct.error) as error:
        print("failed to read H(k)/S(k) output: {}".format(error))
        return 1

    if len(actual) != len(reference):
        print("matrix element count differs: binary={} reference={}".format(len(actual), len(reference)))
        return 1

    for index, (actual_value, reference_value) in enumerate(zip(actual, reference)):
        if abs(actual_value - reference_value) > tolerance:
            print(
                "matrix element {} differs: binary={} reference={} tolerance={}".format(
                    index, actual_value, reference_value, tolerance
                )
            )
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
