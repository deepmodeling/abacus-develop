#!/usr/bin/env python3
import math
import pathlib
import subprocess
import sys
import tempfile
import unittest


COMPARE_WANNIER90 = pathlib.Path(__file__).with_name("compare_wannier90.py")


def matrix_multiply(left, right):
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def conjugate_transpose(matrix):
    return [
        [matrix[row][column].conjugate() for row in range(len(matrix))]
        for column in range(len(matrix[0]))
    ]


def format_amn(matrices, title=None):
    band_count = len(matrices[0])
    projection_count = len(matrices[0][0])
    lines = [] if title is None else [title]
    lines.append("{} {} {}".format(band_count, len(matrices), projection_count))
    for kpoint, matrix in enumerate(matrices, 1):
        for projection in range(projection_count):
            for band in range(band_count):
                value = matrix[band][projection]
                lines.append(
                    "{} {} {} {:.16e} {:.16e}".format(
                        band + 1, projection + 1, kpoint, value.real, value.imag
                    )
                )
    return "\n".join(lines) + "\n"


def format_mmn(matrix_or_links, title=None):
    if isinstance(matrix_or_links, dict):
        links = matrix_or_links
    else:
        links = {(1, 1, 0, 0, 0): matrix_or_links}
    band_count = len(next(iter(links.values())))
    kpoint_count = max(max(link[0], link[1]) for link in links)
    neighbour_count = len(links) // kpoint_count
    if neighbour_count * kpoint_count != len(links):
        raise ValueError("MMN fixture must have the same number of links per k-point")
    lines = [] if title is None else [title]
    lines.append("{} {} {}".format(band_count, kpoint_count, neighbour_count))
    for link, matrix in sorted(links.items()):
        lines.append("{} {} {} {} {}".format(*link))
        for column in range(band_count):
            for row in range(band_count):
                value = matrix[row][column]
                lines.append("{:.16e} {:.16e}".format(value.real, value.imag))
    return "\n".join(lines) + "\n"


class CompareWannier90Test(unittest.TestCase):
    def run_compare(self, file_type, reference_content, calculated_content, tolerance):
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = pathlib.Path(temporary_directory)
            reference_path = temporary_path / "reference.{}".format(file_type)
            calculated_path = temporary_path / "calculated.{}".format(file_type)
            reference_path.write_text(reference_content, encoding="utf-8")
            calculated_path.write_text(calculated_content, encoding="utf-8")
            return subprocess.run(
                [
                    sys.executable,
                    str(COMPARE_WANNIER90),
                    file_type,
                    str(reference_path),
                    str(calculated_path),
                    str(tolerance),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )

    def test_amn_band_gauge_rotation_preserves_gram_matrices(self):
        reference_matrices = [
            [[1.0 + 0.2j, 0.3 - 0.1j], [-0.4 + 0.5j, 0.8 + 0.2j]],
            [[0.2 - 0.7j, 0.9 + 0.1j], [0.6 + 0.3j, -0.2 + 0.4j]],
        ]
        inverse_root_two = 1.0 / math.sqrt(2.0)
        gauge = [
            [inverse_root_two, 1j * inverse_root_two],
            [1j * inverse_root_two, inverse_root_two],
        ]
        calculated_matrices = [matrix_multiply(gauge, matrix) for matrix in reference_matrices]

        result = self.run_compare(
            "amn",
            format_amn(reference_matrices),
            format_amn(calculated_matrices, "Created on an arbitrary date"),
            1.0e-12,
        )

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_amn_change_above_tolerance_fails(self):
        reference = [[[1.0 + 0.0j], [0.0 + 0.0j]]]
        calculated = [[[1.001 + 0.0j], [0.0 + 0.0j]]]

        result = self.run_compare("amn", format_amn(reference), format_amn(calculated), 1.0e-6)

        self.assertEqual(result.returncode, 1)
        self.assertIn("AMN Gram matrices differ", result.stdout)

    def test_mmn_left_and_right_gauge_rotations_preserve_singular_values(self):
        reference = [[0.8 + 0.1j, 0.2 - 0.3j], [-0.4 + 0.2j, 0.6 + 0.1j]]
        inverse_root_two = 1.0 / math.sqrt(2.0)
        left_gauge = [
            [inverse_root_two, inverse_root_two],
            [-inverse_root_two, inverse_root_two],
        ]
        right_gauge = [[1j, 0.0], [0.0, -1j]]
        calculated = matrix_multiply(
            matrix_multiply(conjugate_transpose(left_gauge), reference), right_gauge
        )

        result = self.run_compare(
            "mmn",
            format_mmn({(1, 1, 1, 0, 0): reference}),
            format_mmn({(1, 1, 1, 0, 0): calculated}, "Created on an arbitrary date"),
            1.0e-12,
        )

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_mmn_consistent_kpoint_gauges_preserve_wilson_traces(self):
        identity = [[1.0 + 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 + 0.0j]]
        links = (
            (1, 2, 0, 0, 0),
            (2, 3, 0, 0, 0),
            (3, 4, 0, 0, 0),
            (4, 1, 0, 0, 0),
        )
        reference = {link: identity for link in links}
        inverse_root_two = 1.0 / math.sqrt(2.0)
        gauges = {
            1: identity,
            2: [[1j, 0.0], [0.0, 1.0]],
            3: [
                [inverse_root_two, inverse_root_two],
                [-inverse_root_two, inverse_root_two],
            ],
            4: [[-1j, 0.0], [0.0, 1.0]],
        }
        calculated = {
            link: matrix_multiply(
                matrix_multiply(gauges[link[0]], matrix),
                conjugate_transpose(gauges[link[1]]),
            )
            for link, matrix in reference.items()
        }

        result = self.run_compare("mmn", format_mmn(reference), format_mmn(calculated), 1.0e-12)

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_mmn_link_specific_gauge_fails_wilson_trace(self):
        identity = [[1.0 + 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 + 0.0j]]
        links = (
            (1, 2, 0, 0, 0),
            (2, 3, 0, 0, 0),
            (3, 4, 0, 0, 0),
            (4, 1, 0, 0, 0),
        )
        reference = {link: identity for link in links}
        calculated = {link: identity for link in links}
        calculated[links[0]] = [[1.0 + 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, -1.0 + 0.0j]]

        result = self.run_compare("mmn", format_mmn(reference), format_mmn(calculated), 1.0e-12)

        self.assertEqual(result.returncode, 1)
        self.assertIn("MMN Wilson-loop trace powers differ", result.stdout)

    def test_mmn_change_above_tolerance_fails(self):
        reference = [[0.8 + 0.1j, 0.2 - 0.3j], [-0.4 + 0.2j, 0.6 + 0.1j]]
        calculated = [[1.001 * value for value in row] for row in reference]

        result = self.run_compare("mmn", format_mmn(reference), format_mmn(calculated), 1.0e-6)

        self.assertEqual(result.returncode, 1)
        self.assertIn("MMN singular values differ", result.stdout)

    def test_eig_uses_explicit_absolute_tolerance_in_ev(self):
        reference = "1 1 -5.0000000\n2 1 1.0000000\n"
        calculated = "1 1 -4.9999995\n2 1 1.0000004\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-6)

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_eig_change_above_tolerance_fails(self):
        reference = "1 1 -5.0000000\n2 1 1.0000000\n"
        calculated = "1 1 -4.9999980\n2 1 1.0000000\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-6)

        self.assertEqual(result.returncode, 1)
        self.assertIn("EIG invariant values differ", result.stdout)

    def test_eig_near_degenerate_cluster_compares_centroid(self):
        reference = "3 64 17.099448659983\n4 64 17.099448766890\n"
        calculated = "3 64 17.099448202403\n4 64 17.099449224468\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-8)

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_eig_micro_ev_cluster_compares_centroid(self):
        reference = "1 58 5.652864406393\n2 58 5.652865590000\n"
        calculated = "1 58 5.652863699263\n2 58 5.652866297117\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-8)

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_eig_near_degenerate_cluster_shift_fails_strict_tolerance(self):
        reference = "3 64 17.099448659983\n4 64 17.099448766890\n"
        calculated = "3 64 17.099448679983\n4 64 17.099448786890\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-8)

        self.assertEqual(result.returncode, 1)
        self.assertIn("EIG invariant values differ", result.stdout)

    def test_eig_near_degenerate_cluster_excessive_split_fails(self):
        reference = "3 64 17.099448659983\n4 64 17.099448766890\n"
        calculated = "3 64 17.099444713436\n4 64 17.099452713436\n"

        result = self.run_compare("eig", reference, calculated, 1.0e-8)

        self.assertEqual(result.returncode, 1)
        self.assertIn("shape differs", result.stdout)

    def test_eig_cluster_internal_distortion_fails(self):
        reference = (
            "1 1 0.000000000\n2 1 0.000001250\n"
            "3 1 0.000001750\n4 1 0.000003000\n"
        )
        calculated = (
            "1 1 0.000000000\n2 1 0.000000100\n"
            "3 1 0.000002900\n4 1 0.000003000\n"
        )

        result = self.run_compare("eig", reference, calculated, 1.0e-8)

        self.assertEqual(result.returncode, 1)
        self.assertIn("shape differs", result.stdout)


if __name__ == "__main__":
    unittest.main()
