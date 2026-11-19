"""Read one ABACUS spectral array from a file of consecutive NPY records."""
from pathlib import Path

import numpy as np
from numpy.lib import format


def check_record(record):
    if isinstance(record, bool) or not isinstance(record, (int, np.integer)) or record < 1:
        raise ValueError("record must be a positive integer")


def read_npy_record(path, record, kind, fields):
    check_record(record)
    path = Path(path)
    end = path.stat().st_size
    with path.open("rb") as stream:
        for _ in range(record):
            if stream.tell() == end:
                raise ValueError(f"{path} has no {kind} record {record}")
            if format.read_magic(stream) != (1, 0):
                raise ValueError(f"Unsupported ABACUS {kind} NPY header in {path}")
            shape, _, dtype = format.read_array_header_1_0(stream)
            if len(shape) != 1 or shape[0] < 1 or dtype.hasobject or not dtype.names:
                raise ValueError(f"Expected a nonempty array of structured {kind} data in {path}")
            required = {"schema_version", "record_bytes", "ionic_step", *fields}
            if not required.issubset(dtype.names):
                raise ValueError(f"Missing {kind} fields in {path}")
            offset = stream.tell()
            size = shape[0] * dtype.itemsize
            if size > end - offset:
                raise ValueError(f"Incomplete {kind} record in {path}")
            stream.seek(size, 1)
        rows = np.memmap(stream, mode="r", offset=offset, dtype=dtype, shape=shape)
    step = rows["ionic_step"]
    if (not np.all(rows["schema_version"] == 1) or not np.all(rows["record_bytes"] == size)
            or step[0] < 1 or not np.all(step == step[0])):
        raise ValueError(f"Invalid {kind} record metadata in {path}")
    return rows
