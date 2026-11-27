#!/usr/bin/env python3

import re
import struct
import sys


FRAME_PATTERN = re.compile(r"^\s*#?\s*Ionic[ _]Step\b", re.IGNORECASE)
DIMENSION_PATTERN = re.compile(r"^\s*(\d+)\s+\d+\s+#\s*size of this matrix", re.IGNORECASE)
COMPLEX_PATTERN = re.compile(r"\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)")
REAL_PATTERN = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def select_first_frame(lines):
    """Return the lines of the first ionic-step record of a text reference."""
    start = None
    for index, line in enumerate(lines):
        if FRAME_PATTERN.match(line):
            if start is None:
                start = index
            else:
                return lines[start:index]
    if start is None:
        return lines
    return lines[start:]


def dense_dimension(lines):
    """Read N from the '# size of this matrix' line of a dense text record."""
    for line in lines:
        match = DIMENSION_PATTERN.match(line)
        if match:
            return int(match.group(1))
    return None


def packed_from_dense(values, dim):
    """Select the row-major upper triangle from a dense row-major matrix."""
    return [values[i * dim + j] for i in range(dim) for j in range(i, dim)]


def read_packed_reference(lines, value_type):
    """Read a legacy reference that stores only the upper triangle."""
    values = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("Row "):
            continue
        if value_type == "complex":
            values.extend(complex(float(real), float(imag)) for real, imag in COMPLEX_PATTERN.findall(stripped))
        else:
            values.extend(float(token) for token in stripped.split())
    return values


def read_reference(filename, value_type):
    with open(filename, "r", encoding="utf-8") as stream:
        lines = stream.readlines()

    frame = select_first_frame(lines)
    dim = dense_dimension(frame)
    if dim is None:
        return read_packed_reference(frame, value_type)

    text = "".join(frame)
    if value_type == "complex":
        values = [complex(float(real), float(imag)) for real, imag in COMPLEX_PATTERN.findall(text)]
    else:
        values = [float(token) for token in REAL_PATTERN.findall(text)]
    # The dense matrix body is the final N*N values of the record; any earlier
    # numbers belong to the metadata and UnitCell blocks.
    values = values[-dim * dim:]
    if len(values) != dim * dim:
        raise ValueError("dense reference does not contain a full {}x{} matrix".format(dim, dim))
    return packed_from_dense(values, dim)


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
