import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from abacus_plot.dos import PDOS, TDOS


def frame(n, fields, step):
    rows = np.zeros(n, dtype=[("record_bytes", "<u8"), ("schema_version", "<i4"),
                              ("ionic_step", "<i4"), *fields])
    rows["record_bytes"] = rows.nbytes
    rows["schema_version"] = 1
    rows["ionic_step"] = step
    return rows


def pdos_frame(nspin, ne, step):
    channels = 2 if nspin == 2 else 1
    fields = [("nspin", "<i4"), ("atom_index", "<i4"), ("l", "<i4"), ("m", "<i4"),
              ("energy_min_eV", "<f8"), ("energy_step_eV", "<f8"), ("species", "S8"),
              ("pdos", "<f8", (channels, ne))]
    rows = frame(3, fields, step)
    rows["nspin"] = nspin
    rows["atom_index"] = [1, 2, 2]
    rows["l"] = [0, 1, 1]
    rows["m"] = [0, 1, -1]
    rows["species"] = [b"H", b"Si", b"Si"]
    rows["energy_min_eV"] = -1.5
    rows["energy_step_eV"] = 0.25
    rows["pdos"] = np.arange(3 * channels * ne).reshape(3, channels, ne) + step
    return rows


class DosNpyTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(temp.name)
        self.path = self.root / "pdos.npy"

    def test_plot_and_export_selected_record(self):
        second = pdos_frame(2, 4, 3)
        with self.path.open("wb") as stream:
            np.save(stream, pdos_frame(1, 2, 1))
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
        np.save(self.path, pdos_frame(4, 3, 1))
        (self.root / "pdoss1.txt").write_text("0.0 1 H 88.0\n")
        self.assertEqual(PDOS(self.root).orbitals[0]["data"][0, 0], 88.0)
        self.assertEqual(PDOS(self.path).nspin, 1)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            PDOS(self.path, record=0)

    def test_incomplete_record_and_foreign_array(self):
        with self.path.open("wb") as stream:
            np.save(stream, pdos_frame(1, 3, 1))
            stream.truncate(stream.tell() - 1)
        with self.assertRaisesRegex(ValueError, "Incomplete PDOS"):
            PDOS(self.path)
        np.save(self.path, np.array([1.0]))
        with self.assertRaisesRegex(ValueError, "structured PDOS"):
            PDOS(self.path)

    def test_both_tdos_spin_channels_and_appended_grid(self):
        columns = ["energy_eV", "dos", "dos_int", "dos_smear", "dos_smear_int"]
        for spin in (1, 2):
            with (self.root / f"doss{spin}g2.npy").open("wb") as out:
                for step in (1, 2):
                    data = frame(step + 1, [("spin_channel", "<i8"), *[(c, "<f8") for c in columns]], step)
                    data["spin_channel"] = spin
                    data["energy_eV"] = np.arange(step + 1)
                    for i, col in enumerate(columns[1:]):
                        data[col] = spin * 10 + i + step
                    np.save(out, data, allow_pickle=False)
        tdos = TDOS(self.root, record=2)
        self.assertEqual(tdos.nspin, 2)
        np.testing.assert_equal(tdos.energy[:, 0], [0, 1, 2])
        np.testing.assert_equal(tdos.dos_smear, [[14, 24]] * 3)
        with self.assertRaisesRegex(ValueError, "no TDOS record"):
            TDOS(self.root, record=3)

    def test_foreign_structured_array(self):
        path = self.root / "doss1.npy"
        np.save(path, frame(2, [("pdos", "<f8")], 1))
        with self.assertRaisesRegex(ValueError, "Missing TDOS fields"):
            TDOS(path)
