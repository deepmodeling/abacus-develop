"""Independent gauge-covariance and malformed-file controls."""

import copy
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

_spec = importlib.util.spec_from_file_location(
    "compare_wannier_bundle", Path(__file__).with_name("compare_wannier_bundle.py"))
comparator = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(comparator)


class WannierBundleTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        random = np.random.RandomState(234)
        amn = random.normal(size=(2, 3, 4)) + 1j * random.normal(size=(2, 3, 4))
        links = [(1, 2, 0, 0, 0), (1, 2, 1, 0, 0),
                 (2, 1, 0, 0, 0), (2, 1, -1, 0, 0)]
        mmn = {link: random.normal(size=(3, 3)) + 1j * random.normal(size=(3, 3))
               for link in links}
        self.reference = {"amn_dims": (3, 2, 4), "mmn_dims": (3, 2, 2),
                          "amn": amn, "mmn": mmn,
                          "eig": np.array([[1., 2., 2.], [1.5, 3., 3.]])}

    def write(self, prefix, bundle, comment=True):
        bands, kpoints, projections = bundle["amn_dims"]
        lines = ["Created by independent test"] if comment else []
        lines.append("{} {} {}".format(bands, kpoints, projections))
        for kpoint in range(kpoints):
            for projection in range(projections):
                for band in range(bands):
                    z = bundle["amn"][kpoint, band, projection]
                    lines.append("{} {} {} {:.17g} {:.17g}".format(
                        band + 1, projection + 1, kpoint + 1, z.real, z.imag))
        prefix.with_suffix(".amn").write_text("\n".join(lines) + "\n")
        lines = ["Created by independent test"] if comment else []
        lines.append("{} {} {}".format(*bundle["mmn_dims"]))
        for link, matrix in bundle["mmn"].items():
            lines.append(" ".join(map(str, link)))
            for column in range(bands):
                for row in range(bands):
                    z = matrix[row, column]
                    lines.append("{:.17g} {:.17g}".format(z.real, z.imag))
        prefix.with_suffix(".mmn").write_text("\n".join(lines) + "\n")
        lines = ["{} {} {:.17g}".format(band + 1, kpoint + 1, bundle["eig"][kpoint, band])
                 for kpoint in range(kpoints) for band in range(bands)]
        prefix.with_suffix(".eig").write_text("\n".join(lines) + "\n")

    def compare(self, candidate):
        return comparator.compare(self.reference, candidate, 1e-8, 1e-8)

    def rotated(self, include_nondegenerate=False):
        random = np.random.RandomState(745)
        rotations = []
        for _ in range(2):
            if include_nondegenerate:
                matrix = random.normal(size=(3, 3)) + 1j * random.normal(size=(3, 3))
                rotation, _ = np.linalg.qr(matrix)
            else:
                block = random.normal(size=(2, 2)) + 1j * random.normal(size=(2, 2))
                rotation = np.zeros((3, 3), dtype=complex)
                rotation[0, 0] = np.exp(1j * random.uniform(-3, 3))
                rotation[1:, 1:], _ = np.linalg.qr(block)
            rotations.append(rotation)
        candidate = copy.deepcopy(self.reference)
        for kpoint, rotation in enumerate(rotations):
            candidate["amn"][kpoint] = rotation.conj().T @ self.reference["amn"][kpoint]
        for link, matrix in self.reference["mmn"].items():
            candidate["mmn"][link] = rotations[link[0] - 1].conj().T @ matrix @ rotations[link[1] - 1]
        return candidate

    def test_identity(self):
        self.assertTrue(self.compare(copy.deepcopy(self.reference))["passed"])

    def test_global_phases_and_exact_degenerate_rotations(self):
        candidate = self.rotated()
        result = self.compare(candidate)
        self.assertTrue(result["passed"], result)
        self.assertLess(result["mmn_max_abs_error"], 1e-12)

    def test_full_file_roundtrip_with_and_without_comment(self):
        first = self.root / "reference"
        second = self.root / "candidate"
        self.write(first, self.reference, False)
        self.write(second, self.rotated(), True)
        result = comparator.compare(comparator.read_bundle(first), comparator.read_bundle(second), 1e-8, 1e-8)
        self.assertTrue(result["passed"], result)

    def test_mmn_column_major_matches_literal_matrix(self):
        path = self.root / "literal.mmn"
        path.write_text("2 1 1\n1 1 0 0 0\n1 2\n3 4\n5 6\n7 8\n")
        _, values = comparator.read_mmn(path)
        expected = np.array([[1 + 2j, 5 + 6j], [3 + 4j, 7 + 8j]])
        np.testing.assert_array_equal(values[(1, 1, 0, 0, 0)], expected)

    def test_amn_only_rotation_is_rejected(self):
        candidate = self.rotated()
        candidate["mmn"] = copy.deepcopy(self.reference["mmn"])
        self.assertFalse(self.compare(candidate)["passed"])

    def test_mmn_only_rotation_is_rejected(self):
        candidate = self.rotated()
        candidate["amn"] = self.reference["amn"].copy()
        self.assertFalse(self.compare(candidate)["passed"])

    def test_independent_link_phase_is_rejected(self):
        candidate = self.rotated()
        candidate["mmn"][(1, 2, 1, 0, 0)] *= np.exp(0.1j)
        self.assertFalse(self.compare(candidate)["passed"])

    def test_nondegenerate_mixing_cannot_hide_behind_equal_spectra(self):
        candidate = self.rotated(include_nondegenerate=True)
        self.assertFalse(self.compare(candidate)["passed"])

    def test_amn_coefficient_corruption(self):
        candidate = self.rotated()
        candidate["amn"][1, 2, 3] += 2e-6j
        self.assertFalse(self.compare(candidate)["passed"])

    def test_mmn_coefficient_corruption(self):
        candidate = self.rotated()
        candidate["mmn"][(2, 1, 0, 0, 0)][1, 2] += 2e-6
        self.assertFalse(self.compare(candidate)["passed"])

    def test_eigenvalue_shift_cannot_be_replaced_by_cluster_mean(self):
        candidate = self.rotated()
        candidate["eig"][0, 1] -= 2e-7
        candidate["eig"][0, 2] += 2e-7
        result = self.compare(candidate)
        self.assertFalse(result["passed"])
        self.assertGreater(result["eig_max_abs_error_eV"], 1e-8)

    def test_near_degenerate_gauge_must_preserve_energy_operator(self):
        self.reference["eig"][:, 2] += 5e-9
        candidate = self.rotated()
        result = self.compare(candidate)
        self.assertTrue(result["passed"], result)
        self.assertLessEqual(result["energy_operator_max_abs_error_eV"], 1e-8)

    def test_changed_link_translation_is_rejected(self):
        candidate = self.rotated()
        matrix = candidate["mmn"].pop((1, 2, 1, 0, 0))
        candidate["mmn"][(1, 2, 0, 1, 0)] = matrix
        with self.assertRaisesRegex(comparator.ComparisonError, "translations differ"):
            self.compare(candidate)

    def test_underdetermined_projection_block_is_explicit_failure(self):
        candidate = copy.deepcopy(self.reference)
        candidate["amn"][0, 2] = candidate["amn"][0, 1]
        with self.assertRaisesRegex(comparator.ComparisonError, "does not determine"):
            self.compare(candidate)

    def invalid_file(self, extension, mutate, pattern):
        prefix = self.root / "invalid"
        self.write(prefix, self.reference)
        path = prefix.with_suffix(extension)
        path.write_text(mutate(path.read_text()))
        with self.assertRaisesRegex(comparator.ComparisonError, pattern):
            comparator.read_bundle(prefix)

    def test_missing_file(self):
        with self.assertRaisesRegex(comparator.ComparisonError, "cannot read"):
            comparator.read_bundle(self.root / "absent")

    def test_truncated_amn(self):
        self.invalid_file(".amn", lambda text: "\n".join(text.splitlines()[:-1]), "AMN records")

    def test_extra_mmn_data(self):
        self.invalid_file(".mmn", lambda text: text + "0 0\n", "MMN records")

    def test_duplicate_amn_index(self):
        def mutate(text):
            lines = text.splitlines()
            lines[3] = lines[2]
            return "\n".join(lines)
        self.invalid_file(".amn", mutate, "duplicate AMN")

    def test_nonfinite_amn_mmn_eig(self):
        for extension in (".amn", ".mmn", ".eig"):
            with self.subTest(extension=extension):
                def mutate(text):
                    lines = text.splitlines()
                    fields = lines[-1].split()
                    fields[-1] = "nan"
                    lines[-1] = " ".join(fields)
                    return "\n".join(lines)
                self.invalid_file(extension, mutate, "nonfinite")

    def test_duplicate_eig_index(self):
        def mutate(text):
            lines = text.splitlines()
            lines[1] = lines[0]
            return "\n".join(lines)
        self.invalid_file(".eig", mutate, "duplicate EIG")

    def test_bad_dimensions(self):
        self.invalid_file(".amn", lambda text: text.replace("3 2 4\n", "-3 2 4\n"), "positive")

    def test_cli_failure_is_nonzero(self):
        args = ["compare_wannier_bundle.py", str(self.root / "missing"), str(self.root / "also_missing"),
                "--matrix-tolerance", "1e-8", "--energy-tolerance", "1e-8"]
        with mock.patch("sys.argv", args), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(comparator.main(), 1)

    def test_invalid_tolerance_is_rejected(self):
        for value in (0, -1, np.nan, np.inf):
            with self.subTest(value=value), self.assertRaises(comparator.ComparisonError):
                comparator.compare(self.reference, self.reference, value, 1e-8)

    def test_finite_mmn_values_cannot_overflow_into_a_pass(self):
        reference = {"amn_dims": (2, 1, 2), "mmn_dims": (2, 1, 1),
                     "amn": np.eye(2, dtype=complex)[None, :, :],
                     "mmn": {(1, 1, 0, 0, 0): np.zeros((2, 2), dtype=complex)},
                     "eig": np.zeros((1, 2))}
        candidate = copy.deepcopy(reference)
        candidate["amn"][0] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        candidate["mmn"][(1, 1, 0, 0, 0)][:] = 1.5e308
        self.assertTrue(np.isfinite(candidate["mmn"][(1, 1, 0, 0, 0)]).all())
        with np.errstate(over="ignore", invalid="ignore"):
            with self.assertRaisesRegex(comparator.ComparisonError, "nonfinite arithmetic in MMN"):
                comparator.compare(reference, candidate, 1e-8, 1e-8)

    def test_finite_eigenvalues_cannot_overflow_into_a_pass(self):
        candidate = copy.deepcopy(self.reference)
        self.reference["eig"][0, 0] = -1.5e308
        candidate["eig"][0, 0] = 1.5e308
        with self.assertRaisesRegex(comparator.ComparisonError, "nonfinite arithmetic in EIG"):
            self.compare(candidate)

    def collect(self, flag, reference_key, complete=True, errexit=True):
        directory = Path(tempfile.mkdtemp(dir=self.root))
        (directory / "OUT.autotest").mkdir()
        self.write(directory / "diamond", self.reference)
        if complete:
            self.write(directory / "OUT.autotest/diamond", self.rotated())
        if reference_key:
            (directory / "result.ref").write_text(reference_key + " 0\n")
        tools = Path(__file__).resolve().parent
        script = ('source "$1/props_mat.sh"; PROPS_SCRIPT_DIR="$1"; '
                  'props_result_file=result.out; toW90="$2"; '
                  'run_mat_props; echo collector_completed')
        command = ["bash"]
        if errexit:
            command.append("-e")
        command += ["-c", script, "collector", str(tools), flag]
        environment = os.environ.copy()
        environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
        result = subprocess.run(command, cwd=directory, env=environment, capture_output=True, text=True)
        output = directory / "result.out"
        return result, output.read_text() if output.exists() else ""

    def test_collector_reference_keys_require_all_checks_without_input_switch(self):
        expected = "CompareAMN_pass 0\nCompareMMN_pass 0\nCompareEIG_pass 0\n"
        for key in ("CompareAMN_pass", "CompareMMN_pass", "CompareEIG_pass"):
            for flag in ("", "0"):
                with self.subTest(key=key, flag=flag):
                    result, output = self.collect(flag, " \t" + key)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(output, expected)

    def test_collector_missing_bundle_fails_without_switch_or_errexit(self):
        for key in ("CompareAMN_pass", "CompareMMN_pass", "CompareEIG_pass"):
            for flag in ("", "0", "1"):
                with self.subTest(key=key, flag=flag):
                    result, output = self.collect(flag, "  " + key, complete=False, errexit=False)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertNotIn("collector_completed", result.stdout)
                    self.assertEqual(output, "")

    def test_collector_requested_output_requires_bundle_without_reference_keys(self):
        result, output = self.collect("1", "", complete=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output, "")

    def test_unrelated_diamond_file_does_not_enable_collector(self):
        result, output = self.collect("", "", complete=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("collector_completed", result.stdout)
        self.assertEqual(output, "")


if __name__ == "__main__":
    unittest.main()
