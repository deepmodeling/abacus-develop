"""Regressions for TDOS formats, spin pairing and plotted data."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from abacus_plot.dos import DOSPlot, PDOS, TDOS


class TDOSTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(self.temp.name)
        self.energy = np.array([-1., 0., 1.])

    def write_spin(self, name, energy, scale):
        path = self.root / name
        raw = scale * np.arange(1., len(energy) + 1.)
        data = np.column_stack((energy, raw, np.cumsum(raw),
                                raw + 10., np.cumsum(raw + 10.)))
        np.savetxt(path, data, header="1 # ionic step\n{} # number of points".format(len(energy)),
                   comments="", fmt="%.12f")
        return path

    def test_formats_and_default_raw_selection(self):
        up = np.array([1., 2., 3.])
        down = np.array([4., 5., 6.])
        smear_up = up + 10.
        smear_down = down + 20.
        formats = [
            (np.column_stack((self.energy, up)), up[:, None], None),
            (np.column_stack((self.energy, up, down)), np.column_stack((up, down)), None),
            (np.column_stack((self.energy, up, up.cumsum(), smear_up, smear_up.cumsum())),
             up[:, None], smear_up[:, None]),
            (np.column_stack((self.energy, up, down, up.cumsum(), down.cumsum(),
                              smear_up, smear_down, smear_up.cumsum(), smear_down.cumsum())),
             np.column_stack((up, down)), np.column_stack((smear_up, smear_down))),
        ]
        for data, raw, smeared in formats:
            with self.subTest(columns=data.shape[1]):
                path = self.root / "format.txt"
                np.savetxt(path, data)
                tdos = TDOS(path)
                np.testing.assert_array_equal(tdos.energy[:, 0], self.energy)
                np.testing.assert_array_equal(tdos.dos, raw)
                self.assertEqual(tdos.nspin, raw.shape[1])
                if smeared is None:
                    self.assertIsNone(tdos.dos_smear)
                else:
                    np.testing.assert_array_equal(tdos.dos_smear, smeared)
                for smearing in (True, False):
                    fig, ax = plt.subplots()
                    if smearing:
                        tdos.plot(fig, ax)
                    else:
                        tdos.plot(fig, ax, smearing=False)
                    expected = smeared if smearing and smeared is not None else raw
                    np.testing.assert_array_equal(ax.lines[0].get_ydata(), expected[:, 0])
                    if tdos.nspin == 2:
                        np.testing.assert_array_equal(ax.lines[1].get_ydata(), -expected[:, 1])
                    self.assertEqual(ax.get_xlabel(), "Energy(eV)")
                    plt.close(fig)

    def test_shift_uses_the_selected_column(self):
        path = self.write_spin("single.txt", self.energy, 1.)
        tdos = TDOS(path)
        bands = (SimpleNamespace(band=tdos.energy[:1], value=0.),
                 SimpleNamespace(band=tdos.energy[1:], value=1.))
        for smearing in (True, False):
            with mock.patch.object(tdos, "set_vcband", return_value=bands) as split:
                with mock.patch.object(tdos, "info"):
                    fig, ax = plt.subplots()
                    tdos.plot(fig, ax, shift=True, smearing=smearing)
                expected = tdos.dos_smear if smearing else tdos.dos
                np.testing.assert_array_equal(split.call_args[0][1], expected)
                np.testing.assert_array_equal(ax.lines[0].get_ydata(), expected[:, 0])

    def test_recognized_pairs_preserve_spin_order(self):
        up = self.write_spin("doss1g1.txt", self.energy, 1.)
        down = self.write_spin("doss2g1.txt", self.energy, 7.)
        for source in (up, down, self.root, [up, down], [down, up]):
            with self.subTest(source=str(source)):
                tdos = TDOS(source)
                self.assertEqual(tdos.nspin, 2)
                np.testing.assert_array_equal(tdos.dos[:, 0], [1., 2., 3.])
                np.testing.assert_array_equal(tdos.dos[:, 1], [7., 14., 21.])

    def test_directory_does_not_mix_geometries(self):
        self.write_spin("doss1g1.txt", self.energy, 1.)
        self.write_spin("doss2g2.txt", self.energy, 7.)
        tdos = TDOS(self.root)
        self.assertEqual(tdos.nspin, 1)
        np.testing.assert_array_equal(tdos.dos[:, 0], [1., 2., 3.])

    def test_spin_grids_require_equal_shape_and_absolute_agreement(self):
        energy = np.array([100., 100.0001, 100.0002])
        up = self.write_spin("doss1g1.txt", energy, 1.)
        for other in (energy + 0.0001, energy[:2], energy[:1]):
            with self.subTest(other=other):
                down = self.write_spin("doss2g1.txt", other, 7.)
                with self.assertRaisesRegex(ValueError, "different energy grids"):
                    TDOS([up, down])
        for other in (energy, energy + 5e-9):
            down = self.write_spin("doss2g1.txt", other, 7.)
            self.assertEqual(TDOS([up, down]).nspin, 2)

    def test_plain_axes_and_missing_subplot_apis(self):
        fig = plt.figure()
        plain = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        DOSPlot(fig, plain, xlabel_params={"size": 13})._set_figure()
        self.assertEqual(plain.get_xlabel(), "Energy(eV)")
        self.assertEqual(plain.xaxis.label.get_fontsize(), 13)
        for get_spec in (None, lambda: object()):
            fig, ax = plt.subplots()
            with mock.patch.object(ax, "get_subplotspec", get_spec, create=True):
                DOSPlot(fig, ax)._set_figure()
            self.assertEqual(ax.get_xlabel(), "Energy(eV)")


class PDOSNumPyTest(unittest.TestCase):
    @staticmethod
    def spin_record():
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        dtype = [("nspin", "i4"), ("atom_index", "i4"), ("l", "i4"), ("m", "i4"),
                 ("energy_min_eV", "f8"), ("energy_step_eV", "f8"), ("species", "S8"),
                 ("pdos", "f8", values.shape)]
        return np.array([(2, 3, 1, -1, -1., .5, b"Fe", values)], dtype=dtype)

    def test_spin_channels_and_energy_grid(self):
        record = self.spin_record()
        values = record["pdos"][0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pdos.npy"
            np.save(path, record)
            for source in (path, path.parent):
                result = PDOS(source)
                self.assertEqual(result.nspin, 2)
                np.testing.assert_array_equal(result.energy[:, 0], [-1., -.5, 0.])
                orbital, = result.orbitals
                self.assertEqual((orbital["index"], orbital["atom_index"], orbital["species"],
                                  orbital["l"], orbital["m"]), (1, 3, "Fe", 1, -1))
                np.testing.assert_array_equal(orbital["data"], values.T)

            text = path.parent / "pdoss1.txt"
            text.write_text("0.0 1 Fe 9.0\n")
            with self.assertRaisesRegex(ValueError, "Pass a specific PDOS file"):
                PDOS(path.parent)
            np.testing.assert_array_equal(PDOS(path).orbitals[0]["data"], values.T)
            self.assertEqual(PDOS(text).orbitals[0]["data"].item(), 9.)

    def test_invalid_record_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pdos.npy"
            for field, values, message in (
                ("nspin", [2, 1], "nspin"),
                ("energy_min_eV", [-1., 0.], "energy_min_eV"),
                ("energy_step_eV", [.5, .25], "energy_step_eV"),
                ("energy_step_eV", [0., 0.], "positive step"),
                ("energy_step_eV", [-.5, -.5], "positive step"),
                ("energy_step_eV", [np.inf, np.inf], "finite"),
            ):
                with self.subTest(field=field, values=values):
                    rows = np.repeat(self.spin_record(), 2)
                    rows[field] = values
                    np.save(path, rows)
                    with self.assertRaisesRegex(ValueError, message):
                        PDOS(path)


if __name__ == "__main__":
    unittest.main()
