import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from abacus_plot.dos import PDOS


def frame(nspin, ne, step):
    channels = 2 if nspin == 2 else 1
    dtype = [("record_bytes", "<u8"), ("schema_version", "<i4"), ("ionic_step", "<i4"),
             ("nspin", "<i4"), ("atom_index", "<i4"), ("l", "<i4"), ("m", "<i4"),
             ("energy_min_eV", "<f8"), ("energy_step_eV", "<f8"), ("species", "S8"),
             ("pdos", "<f8", (channels, ne))]
    rows = np.zeros(3, dtype=dtype)
    rows["record_bytes"] = rows.nbytes
    rows["schema_version"] = 1
    rows["nspin"] = nspin
    rows["ionic_step"] = step
    rows["atom_index"] = [1, 2, 2]
    rows["l"] = [0, 1, 1]
    rows["m"] = [0, 1, -1]
    rows["species"] = [b"H", b"Si", b"Si"]
    rows["energy_min_eV"] = -1.5
    rows["energy_step_eV"] = 0.25
    rows["pdos"] = np.arange(3 * channels * ne).reshape(3, channels, ne) + step
    return rows


class PdosNpyTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(temp.name)
        self.path = self.root / "pdos.npy"

    def test_plot_and_export_selected_record(self):
        second = frame(2, 4, 3)
        with self.path.open("wb") as stream:
            np.save(stream, frame(1, 2, 1))
            np.save(stream, second)
        self.assertEqual(PDOS(self.root).nspin, 1)
        reader = PDOS(self.path, record=2)
        self.assertEqual(reader.nspin, 2)
        np.testing.assert_array_equal(reader.energy[:, 0], [-1.5, -1.25, -1, -0.75])
        for i, orbital in enumerate(reader.orbitals):
            np.testing.assert_array_equal(orbital["data"], second["pdos"][i].T)
        self.assertEqual(reader.orbitals[-1]["m"], -1)
        expected = second["pdos"][1:].sum(axis=0).T
        fig, ax = plt.subplots()
        reader.plot(fig, ax, species=["Si"])
        lines = {line.get_label(): line.get_ydata() for line in ax.lines}
        np.testing.assert_array_equal(lines[r"Si$\uparrow$"], expected[:, 0])
        np.testing.assert_array_equal(lines[r"Si$\downarrow$"], -expected[:, 1])
        reader.write(species=["Si"], outdir=self.root)
        output = np.loadtxt(self.root / "PDOS_FILE/species-Si.dat")
        np.testing.assert_array_equal(output[:, 0], reader.energy[:, 0])
        np.testing.assert_array_equal(output[:, 1:], expected)
        with self.assertRaisesRegex(ValueError, "no PDOS record 3"):
            PDOS(self.path, record=3)

    def test_text_precedence_and_explicit_binary_path(self):
        np.save(self.path, frame(4, 3, 1))
        (self.root / "pdoss1.txt").write_text("0.0 1 H 88.0\n")
        self.assertEqual(PDOS(self.root).orbitals[0]["data"][0, 0], 88.0)
        self.assertEqual(PDOS(self.path).nspin, 1)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            PDOS(self.path, record=0)

    def test_incomplete_record_and_foreign_array(self):
        with self.path.open("wb") as stream:
            np.save(stream, frame(1, 3, 1))
            stream.truncate(stream.tell() - 1)
        with self.assertRaisesRegex(ValueError, "Incomplete PDOS"):
            PDOS(self.path)
        np.save(self.path, np.array([1.0]))
        with self.assertRaisesRegex(ValueError, "PDOS orbitals"):
            PDOS(self.path)
