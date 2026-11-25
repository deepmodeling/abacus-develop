"""Collect the integrated DOS and check LCAO spectra against a text reference."""

from pathlib import Path
import sys

import numpy as np


def read_reference(path):
    ref = {}
    for line in path.read_text().splitlines():
        line = line.partition("#")[0].strip()
        if line:
            key, values = line.split(maxsplit=1)
            if key in ("tdos", "pdos"):
                ref.setdefault(key, []).append(values)
            else:
                ref[key] = float(values) if key in ("energy_min_eV", "energy_step_eV") else int(values)

    spins = 2 if ref["nspin"] == 2 else 1
    ref["total_samples"] = np.loadtxt(ref.pop("tdos"), ndmin=1,
                                     dtype=[("energy_index", int), ("spins", float, (spins, 4))])
    ref["samples"] = np.loadtxt(ref.pop("pdos"), ndmin=1, encoding="utf-8",
                               dtype=[("atom_index", int), ("species", object), ("l", int), ("m", int),
                                      ("energy_index", int), ("pdos", float, (spins,))])
    return ref


def check_pdos(path, ref):
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


def collect_dos(directory, nspin, reference):
    ref = None
    if reference.is_file():
        ref = read_reference(reference)
        np.testing.assert_equal(nspin, ref["nspin"])
        check_pdos(directory / "pdos.npy", ref)

    fields = ("dos", "dos_int", "dos_smear", "dos_smear_int")
    spins = 2 if nspin == 2 else 1
    dtype = np.dtype([("energy_eV", "f8")] + [(name, "f8", (spins,)) for name in fields])
    rows = np.load(directory / "dos.npy", mmap_mode="r", allow_pickle=False)
    np.testing.assert_equal(rows.dtype, dtype)
    if rows.ndim != 1 or not rows.size:
        raise ValueError("Expected one total DOS record per energy point")
    if ref:
        grid = ref["energy_min_eV"] + np.arange(ref["npoints"]) * ref["energy_step_eV"]
        np.testing.assert_allclose(rows["energy_eV"], grid, atol=5e-7, rtol=0)
        for sample in ref["total_samples"]:
            row = rows[sample["energy_index"]]
            values = np.column_stack([row[name] for name in fields])
            np.testing.assert_allclose(values, sample["spins"], rtol=1e-3, atol=5e-5)
    # Match the last bin and last spin used by the text collector.
    return rows["dos_int"][-1, -1]


if __name__ == "__main__":
    value = collect_dos(Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3]))
    print(f"totaldosref {value}")
