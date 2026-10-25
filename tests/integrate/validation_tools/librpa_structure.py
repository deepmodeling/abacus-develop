"""Validate structure metadata shared by the LibRPA producer checks."""

import math


def check_spin_symmetry(lines, symmetry_end, path, minimum_rows, symmetry_rows):
    if symmetry_end is None or len(lines) <= symmetry_end:
        raise ValueError("{} has no spin_symmetry block".format(path))
    header = lines[symmetry_end].split()
    if len(header) != 3 or header[0] != "spin_symmetry":
        raise ValueError("{} has an invalid spin_symmetry header".format(path))
    try:
        grey_group, spin_source = (int(value) for value in header[1:])
    except ValueError:
        raise ValueError("{} has an invalid spin_symmetry header".format(path))
    if grey_group not in (0, 1) or spin_source not in (0, 1) or minimum_rows < 0:
        raise ValueError("{} has invalid spin_symmetry flags".format(path))
    spin_rows = lines[symmetry_end + 1 :]
    if symmetry_rows is None or symmetry_rows < minimum_rows or len(spin_rows) != symmetry_rows:
        raise ValueError("{} has an unexpected spin_symmetry row count".format(path))
    expected_fields = 9 if spin_source == 1 else 1
    for line in spin_rows:
        tokens = line.split()
        if len(tokens) != expected_fields:
            raise ValueError("{} has an invalid spin_symmetry row".format(path))
        try:
            if int(tokens[0]) not in (0, 1):
                raise ValueError
        except ValueError:
            raise ValueError("{} has an invalid spin_symmetry flag".format(path))
        for token in tokens[1:]:
            try:
                if not math.isfinite(float(token)):
                    raise ValueError
            except ValueError:
                raise ValueError("{} has a non-finite spin_symmetry value".format(path))
