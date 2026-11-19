import tempfile
import unittest
from pathlib import Path

import numpy as np

from abacus_plot.dos import TDOS


def frame(n, fields, step):
    rows = np.zeros(n, dtype=[("record_bytes", "<u8"), ("schema_version", "<i4"),
                              ("ionic_step", "<i4"), *fields])
    rows["record_bytes"] = rows.nbytes
    rows["schema_version"] = 1
    rows["ionic_step"] = step
    return rows


class TdosNpyTest(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)

    def test_both_tdos_spin_channels_and_appended_grid(self):
        columns = ["energy_eV", "dos", "dos_int", "dos_smear", "dos_smear_int"]
        for spin in (1, 2):
            with (self.directory / f"doss{spin}g2.npy").open("wb") as out:
                for step in (1, 2):
                    data = frame(step + 1, [("spin_channel", "<i8"), *[(c, "<f8") for c in columns]], step)
                    data["spin_channel"] = spin
                    data["energy_eV"] = np.arange(step + 1)
                    for i, col in enumerate(columns[1:]):
                        data[col] = spin * 10 + i + step
                    np.save(out, data, allow_pickle=False)
        tdos = TDOS(self.directory, record=2)
        self.assertEqual(tdos.nspin, 2)
        np.testing.assert_equal(tdos.energy[:, 0], [0, 1, 2])
        np.testing.assert_equal(tdos.dos_smear, [[14, 24]] * 3)
        with self.assertRaisesRegex(ValueError, "no TDOS record"):
            TDOS(self.directory, record=3)

    def test_foreign_structured_array(self):
        path = self.directory / "doss1.npy"
        np.save(path, frame(2, [("pdos", "<f8")], 1))
        with self.assertRaisesRegex(ValueError, "Missing TDOS fields"):
            TDOS(path)
