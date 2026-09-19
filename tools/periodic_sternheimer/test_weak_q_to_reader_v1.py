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

    def test_packages_only_representative_k_matrices_with_complete_routes(self):
        output = self.root / "partial"
        matrices = {
            (0, 1): -np.eye(2, dtype=np.complex128),
            (2, 1): -2.0 * np.eye(2, dtype=np.complex128),
        }
        routes = (
            dict(representative_ik_full=0, member_ik_full=0, spatial_isym=0,
                 time_reversal=False, fold_G=(0, 0, 0)),
            dict(representative_ik_full=0, member_ik_full=1, spatial_isym=4,
                 time_reversal=False, fold_G=(0, 0, 0)),
            dict(representative_ik_full=2, member_ik_full=2, spatial_isym=0,
                 time_reversal=False, fold_G=(0, 0, 0)),
            dict(representative_ik_full=2, member_ik_full=3, spatial_isym=4,
                 time_reversal=False, fold_G=(0, 0, 0)),
        )

        report = converter.package_partial_k_response(
            matrices,
            output,
            iq=2,
            frequencies={1: (0.125, 0.25)},
            atom_naux=(1, 1),
            full_kpoints=((0.0, 0.0, 0.0), (0.5, 0.0, 0.0),
                          (0.0, 0.5, 0.0), (0.5, 0.5, 0.0)),
            fixed_q_routes=routes,
            qpoint=(0.5, 0.0, 0.0),
            qweight=0.25,
        )

        self.assertTrue(report["partial_response_complete"])
        self.assertEqual(report["representative_kpoints"], [0, 2])
        self.assertTrue((output / "v1_sternheimer_partial_manifest_iq_2.dat").is_file())
        self.assertTrue((output / "v1_sternheimer_symmetry_routes_iq_2.dat").is_file())
        self.assertTrue((output / "v1_sternheimer_full_kpoints.dat").is_file())
        manifest = (output / "v1_sternheimer_partial_manifest_iq_2.dat").read_text()
        self.assertIn("2 0 1 v1_sternheimer_chi0_iq_2_ik_0_ifreq_1.dat", manifest)
        self.assertIn("2 2 1 v1_sternheimer_chi0_iq_2_ik_2_ifreq_1.dat", manifest)

    def test_partial_k_packaging_rejects_nonrepresentative_matrix(self):
        with self.assertRaisesRegex(ValueError, "representative"):
            converter.package_partial_k_response(
                {(1, 1): -np.eye(1, dtype=np.complex128)},
                self.root / "bad-partial",
                iq=1,
                frequencies={1: (0.1, 0.2)},
                atom_naux=(1,),
                full_kpoints=((0.0, 0.0, 0.0), (0.5, 0.0, 0.0)),
                fixed_q_routes=(
                    dict(representative_ik_full=0, member_ik_full=0, spatial_isym=0,
                         time_reversal=False, fold_G=(0, 0, 0)),
                    dict(representative_ik_full=0, member_ik_full=1, spatial_isym=1,
                         time_reversal=False, fold_G=(0, 0, 0)),
                ),
                qpoint=(0.5, 0.0, 0.0),
                qweight=0.5,
            )

    def test_partial_k_packaging_uses_sum_not_average_completion(self):
        output = self.root / "partial-sum"
        raw = np.array([[1.0, 2.0j], [3.0j, 4.0]], dtype=np.complex128)
        converter.package_partial_k_response(
            {(0, 1): raw},
            output,
            iq=1,
            frequencies={1: (0.1, 0.2)},
            atom_naux=(1, 1),
            full_kpoints=((0.0, 0.0, 0.0),),
            fixed_q_routes=(
                dict(representative_ik_full=0, member_ik_full=0, spatial_isym=0,
                     time_reversal=False, fold_G=(0, 0, 0)),
            ),
            qpoint=(0.0, 0.0, 0.0),
            qweight=1.0,
        )
        packaged = converter.read_reader_v1(
            output / "v1_sternheimer_chi0_iq_1_ik_0_ifreq_1.dat").matrix
        np.testing.assert_array_equal(packaged, raw + raw.conj().T)

    def test_merges_complete_weak_columns_before_partial_packaging(self):
        unit_dir = self.root / "unit"
        unit_dir.mkdir()
        input_path = unit_dir / "unit_input.dat"
        input_path.write_text(
            "format_version 1\n"
            "iq 1\n"
            "full_kpoints 2\n"
            "qpoint 0.5 0 0\n"
            "kpoint 1 0 0 0 1 1 1\n"
            "kpoint 2 0.5 0 0 1 1 1\n"
            "auxiliary_channel 1 atom_index 0 atom_local_index 0\n"
            "auxiliary_channel 2 atom_index 1 atom_local_index 0\n"
            "manifest_complete yes\n")
        columns_path = unit_dir / "unit_columns.dat"
        columns_path.write_text(
            "input_manifest_sha256 local-test\n"
            "ifrequency 1\n"
            "omega_Ha 0.1\n"
            "weight_Ha 0.2\n"
            "owned_column 1\n"
            "response 1 1 -1.0 0.0\n"
            "response 2 1 0.2 0.0\n"
            "owned_column 2\n"
            "response 1 2 0.1 0.0\n"
            "response 2 2 -2.0 0.0\n"
            "columns_complete yes\n")
        audit_path = unit_dir / "unit.dat"
        audit_path.write_text(
            "unit_complete yes\n"
            "all_source_bands yes\n"
            "iq 1\n"
            "source_k 1\n"
            "input_manifest unit_input.dat\n"
            "input_manifest_sha256 local-test\n"
            "fixed_q_route 1 0 0 0 0 0 0 0\n"
            "fixed_q_route 1 0 1 1 0 0 0 0\n"
            "column_file 1 unit_columns.dat\n")

        output = self.root / "merged"
        report = converter.merge_weak_q_unit_columns([audit_path], output, qweight=0.5)

        self.assertTrue(report["partial_response_complete"])
        self.assertEqual(report["representative_kpoints"], [0])
        matrix = converter.read_reader_v1(
            output / "v1_sternheimer_chi0_iq_1_ik_0_ifreq_1.dat").matrix
        raw = np.array([[-1.0, 0.1], [0.2, -2.0]], dtype=np.complex128)
        np.testing.assert_array_equal(matrix, raw + raw.conj().T)


if __name__ == "__main__":
    unittest.main()
