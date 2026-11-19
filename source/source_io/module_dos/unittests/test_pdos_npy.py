"""Read the spectrum from the native PDOS writer with NumPy."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main(writer_test):
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        if exc.name != "numpy":
            raise
        print("NumPy is required for the PDOS round-trip test")
        return 77

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "pdos.npy"
        subprocess.run([writer_test, "--gtest_filter=WritePdosTextTest.NpySumsZetasAndSpinorsBeforeThreshold"],
                       env=dict(os.environ, ABACUS_TEST_PDOS_NPY=str(path)), check=True)
        with path.open("rb") as stream:
            rows = np.load(stream, allow_pickle=False)
            assert not stream.read(1), "NPY record count leaves trailing data"

    assert rows.shape == (4,)
    assert rows.dtype.names == ("nspin", "atom_index", "l", "m", "energy_min_eV", "energy_step_eV", "species", "pdos")
    for field in ("nspin", "atom_index", "l", "m"):
        assert rows[field].dtype == np.int32
    for field in ("energy_min_eV", "energy_step_eV", "pdos"):
        assert rows[field].dtype == np.float64
    np.testing.assert_array_equal(rows["nspin"], 4)
    np.testing.assert_array_equal(rows["atom_index"], 1)
    np.testing.assert_array_equal(rows["species"], b"Si")
    np.testing.assert_array_equal(rows["l"], [0, 1, 1, 1])
    np.testing.assert_array_equal(rows["m"], [0, 0, 1, -1])
    np.testing.assert_array_equal(rows["energy_min_eV"], -1.0)
    np.testing.assert_array_equal(rows["energy_step_eV"], 0.5)
    assert rows["pdos"].shape == (4, 1, 3)
    expected = [[3., 0., 0.], [26., 1.2e-6, 0.], [34., 1.2e-6, 0.], [42., 1.2e-6, 0.]]
    np.testing.assert_allclose(rows["pdos"][:, 0], expected, rtol=1e-15, atol=0.)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
