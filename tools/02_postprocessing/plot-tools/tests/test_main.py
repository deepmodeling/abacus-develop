"""Exercise the actual CLI selector-to-PDOS plotting path."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from abacus_plot.main import main


class PDOSCommandTest(unittest.TestCase):
    def test_list_and_dictionary_layouts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pdos = root / "pdoss1g1.txt"
            pdos.write_text("# npoints: 3\n"
                            "-1 1 Si 1\n-1 2 O 2\n"
                            "0 1 Si 2\n0 2 O 3\n"
                            "1 1 Si 3\n1 2 O 4\n")
            selectors = {"species": ["Si", "O"], "atom_index": [1, 2], "index": [1, 2]}
            for key, values in selectors.items():
                for as_dict in ((False, True) if key == "species" else (True,)):
                    with self.subTest(key=key, dictionary=as_dict):
                        plt.close("all")
                        image = root / "{}-{}.png".format(key, as_dict)
                        config = {"pdosfile": str(pdos), "pdosfig": str(image),
                                  "figsize": [4, 3], "dpi": 80,
                                  "xlabel_params": {"size": 13},
                                  key: {str(value): [0] for value in values} if as_dict else values}
                        config_path = root / "config.json"
                        config_path.write_text(json.dumps(config))
                        with mock.patch("sys.argv", ["abacus-plot", "-d", "-p", "-f", str(config_path)]):
                            main()
                        fig = plt.gcf()
                        self.assertEqual(len(fig.axes), 2 if as_dict else 1)
                        self.assertEqual(fig.axes[-1].get_xlabel(), "Energy(eV)")
                        self.assertEqual(fig.axes[-1].xaxis.label.get_fontsize(), 13)
                        if as_dict:
                            self.assertEqual(fig.axes[0].get_xlabel(), "")
                            self.assertEqual(fig.axes[0].lines[0].get_ydata().tolist(), [1., 2., 3.])
                            self.assertEqual(fig.axes[1].lines[0].get_ydata().tolist(), [2., 3., 4.])
                        else:
                            self.assertEqual(len(fig.axes[0].get_legend_handles_labels()[1]), 3)
                            self.assertEqual(fig.axes[0].lines[0].get_ydata().tolist(), [3., 5., 7.])
                            self.assertEqual(fig.axes[0].lines[1].get_ydata().tolist(), [1., 2., 3.])
                            self.assertEqual(fig.axes[0].lines[2].get_ydata().tolist(), [2., 3., 4.])
                        self.assertGreater(image.stat().st_size, 0)
            plt.close("all")


if __name__ == "__main__":
    unittest.main()
