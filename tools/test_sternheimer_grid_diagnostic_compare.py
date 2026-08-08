#!/usr/bin/env python3

import csv
import math
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("sternheimer_grid_diagnostic_compare.py")


def write_operator(path, grid, scale):
    text = """ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1
kind delta_grid_matrices
grid {grid} {grid} {grid}
spin 1
occupied 1
virtuals 2
auxiliaries 3
volume_element 1.00000000000000000e-02
energy_unit Rydberg
storage column_major
matrix overlap
0 0 {s} 0
0 1 0 0
1 0 0 0
1 1 {s} 0
matrix kinetic
0 0 {s} 0
0 1 0 0
1 0 0 0
1 1 {s} 0
matrix local_potential
0 0 {-s} 0
0 1 0 0
1 0 0 0
1 1 {-s} 0
matrix nonlocal
0 0 {s} 0
0 1 0 0
1 0 0 0
1 1 {s} 0
matrix hamiltonian
0 0 {s} 0
0 1 0 0
1 0 0 0
1 1 {s} 0
matrix occupied_virtual_overlap
0 0 {s} 0
0 1 0 0
""".format(grid=grid, s=scale, **{"-s": -scale})
    path.write_text(text)


def write_perturbation(path, grid, scale):
    text = """ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1
kind delta_perturbation_tensor
grid {grid} {grid} {grid}
spin 1
occupied 1
virtuals 2
auxiliaries 3
volume_element 1.00000000000000000e-02
potential_unit Rydberg
storage row_major
tensor perturbation occupied virtual auxiliary
0 0 0 {s} 0
0 0 1 0 0
0 0 2 0 0
0 1 0 0 0
0 1 1 0 0
0 1 2 {s} 0
""".format(grid=grid, s=scale)
    path.write_text(text)


class GridDiagnosticCompareTest(unittest.TestCase):
    def test_campaign_writes_expected_relative_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run40 = root / "run40"
            run50 = root / "run50"
            run40.mkdir()
            run50.mkdir()
            write_operator(run40 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat", 10, 1.0)
            write_operator(run50 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat", 12, 1.1)
            write_perturbation(run40 / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat", 10, 2.0)
            write_perturbation(run50 / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat", 12, 2.2)
            output = root / "comparison.csv"

            subprocess.run(
                [str(SCRIPT), str(run40), str(run50), str(output)], check=True, universal_newlines=True
            )
            with output.open() as source:
                rows = list(csv.DictReader(source))

            self.assertEqual(len(rows), 7)
            by_quantity = {row["quantity"]: row for row in rows}
            self.assertAlmostEqual(float(by_quantity["overlap"]["relative_difference"]), 0.1)
            self.assertAlmostEqual(float(by_quantity["perturbation"]["relative_difference"]), 0.1)
            self.assertAlmostEqual(float(by_quantity["overlap"]["relative_norm_change"]), 0.1)
            self.assertEqual(by_quantity["overlap"]["profile_kind"], "diagonal")
            self.assertAlmostEqual(float(by_quantity["overlap"]["profile_relative_difference"]), 0.1)
            self.assertEqual(by_quantity["perturbation"]["profile_kind"], "occupied_auxiliary_norm")
            self.assertAlmostEqual(
                float(by_quantity["perturbation"]["profile_relative_difference"]), 0.1
            )
            self.assertEqual(by_quantity["hamiltonian"]["grid40"], "10x10x10")
            self.assertEqual(by_quantity["hamiltonian"]["grid50"], "12x12x12")

    def test_dimension_mismatch_fails_without_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run40 = root / "run40"
            run50 = root / "run50"
            run40.mkdir()
            run50.mkdir()
            write_operator(run40 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat", 10, 1.0)
            write_operator(run50 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat", 12, 1.0)
            text = (run50 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat").read_text()
            (run50 / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat").write_text(
                text.replace("virtuals 2", "virtuals 3")
            )
            write_perturbation(run40 / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat", 10, 1.0)
            write_perturbation(run50 / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat", 12, 1.0)
            output = root / "comparison.csv"

            result = subprocess.run(
                [str(SCRIPT), str(run40), str(run50), str(output)],
                universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())

    def test_virtual_phase_flip_is_removed_from_invariant_profile(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run40 = root / "run40"
            run50 = root / "run50"
            run40.mkdir()
            run50.mkdir()
            for directory in (run40, run50):
                write_operator(directory / "STERNHEIMER_DELTA_GRID_MATRICES_spin_1.dat", 10, 1.0)
                write_perturbation(directory / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat", 10, 2.0)
            perturbation50 = run50 / "STERNHEIMER_DELTA_PERTURBATION_spin_1.dat"
            perturbation50.write_text(perturbation50.read_text().replace("0 1 2 2.0 0", "0 1 2 -2.0 0"))
            output = root / "comparison.csv"

            subprocess.run(
                [str(SCRIPT), str(run40), str(run50), str(output)],
                check=True, universal_newlines=True
            )
            with output.open() as source:
                perturbation = next(
                    row for row in csv.DictReader(source) if row["quantity"] == "perturbation"
                )

            self.assertAlmostEqual(float(perturbation["relative_difference"]), math.sqrt(2.0))
            self.assertAlmostEqual(float(perturbation["relative_norm_change"]), 0.0)
            self.assertAlmostEqual(float(perturbation["profile_relative_difference"]), 0.0)


if __name__ == "__main__":
    unittest.main()
