#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

import numpy as np

import weak_q_to_reader_v1 as converter


class ReaderV1PackagingTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def test_round_trip_preserves_hermitian_potential_response(self):
        raw = np.array(
            [
                [-4.0 + 3.0e-13j, 0.4 + 0.2j, -0.1 + 0.3j],
                [0.4 - 0.2j + 2.0e-13, -1.5, 0.2 - 0.1j],
                [-0.1 - 0.3j, 0.2 + 0.1j - 3.0e-13, -0.8],
            ],
            dtype=np.complex128,
        )
        expected = 0.5 * (raw + raw.conj().T)
        output = self.root / "v1_sternheimer_chi0_iq_2_ifreq_1_rank0.dat"

        report = converter.package_matrix(
            raw,
            output,
            iq=2,
            ifrequency=1,
            omega=0.125,
            weight=0.25,
            atom_naux=(1, 2),
        )
        actual = converter.read_reader_v1(output)

        self.assertEqual(actual.iq, 2)
        self.assertEqual(actual.ifrequency, 1)
        self.assertEqual(actual.omega, 0.125)
        self.assertEqual(actual.weight, 0.25)
        self.assertEqual(actual.atom_naux, (1, 2))
        np.testing.assert_array_equal(actual.matrix, expected)
        self.assertEqual(report["matrix_kind"], "potential_potential_response")
        self.assertEqual(report["coulomb_transform"], "none")
        self.assertLess(report["round_trip_relative_frobenius"], 1.0e-15)
        self.assertFalse((output.parent / (output.name + ".tmp")).exists())

    def test_trace_log_matches_direct_sternheimer_contract(self):
        coulomb = np.array([[4.0, 0.3], [0.3, 2.0]], dtype=np.complex128)
        response = np.array([[-0.8, 0.04j], [-0.04j, -0.2]], dtype=np.complex128)
        output = self.root / "response.dat"
        converter.package_matrix(
            response,
            output,
            iq=1,
            ifrequency=2,
            omega=0.5,
            weight=0.75,
            atom_naux=(1, 1),
        )

        restored = converter.read_reader_v1(output).matrix
        direct = converter.sternheimer_trace_log(coulomb, response, threshold=1.0e-12)
        packaged = converter.sternheimer_trace_log(coulomb, restored, threshold=1.0e-12)
        self.assertAlmostEqual(direct.real, packaged.real, places=14)
        self.assertAlmostEqual(direct.imag, packaged.imag, places=14)

    def test_coulomb_round_trip_preserves_finite_part_matrix(self):
        raw = np.array(
            [
                [-3.0, 0.2 + 0.1j, -0.4],
                [0.2 - 0.1j, 4.0, 0.3j],
                [-0.4, -0.3j, 9.0],
            ],
            dtype=np.complex128,
        )
        output = self.root / "v1_coulomb_full_iq_1_rank0.dat"

        report = converter.package_coulomb_matrix(
            raw,
            output,
            iq=1,
            atom_naux=(1, 2),
        )
        actual = converter.read_coulomb_v1(output)

        self.assertEqual(actual.iq, 1)
        self.assertEqual(actual.ifrequency, 0)
        self.assertEqual(actual.atom_naux, (1, 2))
        np.testing.assert_array_equal(actual.matrix, raw)
        self.assertEqual(report["matrix_kind"], "finite_part_coulomb")
        self.assertLess(report["round_trip_relative_frobenius"], 1.0e-15)

    def test_reads_complete_full_coulomb_text(self):
        source = self.root / "finite-coulomb.dat"
        source.write_text(
            "format_version 1\n"
            "full_matrix_rows 2\n"
            "full_matrix_columns 2\n"
            "coulomb_integral 1 1 -2.0 0.0\n"
            "coulomb_integral 1 2 0.25 0.5\n"
            "coulomb_integral 2 1 0.25 -0.5\n"
            "coulomb_integral 2 2 3.0 0.0\n"
            "full_matrix_complete yes\n"
        )

        actual = converter.read_full_coulomb_text(source)

        np.testing.assert_array_equal(
            actual,
            np.array([[-2.0, 0.25 + 0.5j], [0.25 - 0.5j, 3.0]], dtype=np.complex128),
        )

    def test_invalid_input_fails_closed_without_output(self):
        cases = (
            (np.eye(2, dtype=np.complex128), dict(iq=0), "iq"),
            (np.eye(2, dtype=np.complex128), dict(atom_naux=(1,)), "dimension"),
            (np.array([[np.nan, 0], [0, 1]], dtype=np.complex128), {}, "finite"),
        )
        defaults = dict(
            iq=1,
            ifrequency=1,
            omega=0.1,
            weight=0.2,
            atom_naux=(1, 1),
        )
        for index, (matrix, changes, message) in enumerate(cases):
            with self.subTest(message=message):
                output = self.root / f"bad-{index}.dat"
                arguments = dict(defaults)
                arguments.update(changes)
                with self.assertRaisesRegex(ValueError, message):
                    converter.package_matrix(matrix, output, **arguments)
                self.assertFalse(output.exists())

    def test_refuses_to_overwrite_existing_output(self):
        output = self.root / "existing.dat"
        output.write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            converter.package_matrix(
                -np.eye(1, dtype=np.complex128),
                output,
                iq=1,
                ifrequency=1,
                omega=0.1,
                weight=0.2,
                atom_naux=(1,),
            )
        self.assertEqual(output.read_bytes(), b"keep")


if __name__ == "__main__":
    unittest.main()
