"""Check an ABACUS PDOS array against selected channels from reference text output."""

import json
import sys

import numpy as np


def check_pdos(path, reference):
    with open(reference) as stream:
        ref = json.load(stream)
    rows = np.load(path, mmap_mode="r", allow_pickle=False)
    spins = 2 if ref["nspin"] == 2 else 1
    species_dtype = rows.dtype["species"]
    assert species_dtype.kind == "S", species_dtype
    expected_dtype = np.dtype([
        ("nspin", "i4"), ("atom_index", "i4"), ("l", "i4"), ("m", "i4"),
        ("energy_min_eV", "f8"), ("energy_step_eV", "f8"), ("species", species_dtype),
        ("pdos", "f8", (spins, ref["npoints"])),
    ])
    np.testing.assert_equal(rows.dtype, expected_dtype)
    np.testing.assert_equal(rows.shape, (ref["records"],))
    for field in ("nspin", "energy_min_eV", "energy_step_eV"):
        np.testing.assert_array_equal(rows[field], ref[field], err_msg=field)

    for sample in ref["samples"]:
        match = ((rows["atom_index"] == sample["atom_index"])
                 & (rows["l"] == sample["l"]) & (rows["m"] == sample["m"]))
        orbital, = rows[match]
        np.testing.assert_equal(orbital["species"], sample["species"].encode())
        values = orbital["pdos"][:, sample["energy_index"]]
        np.testing.assert_allclose(values, sample["pdos"], rtol=1e-3, atol=5e-5, err_msg=str(sample))


if __name__ == "__main__":
    check_pdos(*sys.argv[1:])
