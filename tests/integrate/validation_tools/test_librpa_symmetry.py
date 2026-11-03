"""Check reciprocal-grid symmetry handling for reader-v1 metadata."""

import unittest

from librpa_symmetry import preserves_bvk_mesh, reciprocal_rotation


class LibrpaSymmetryTests(unittest.TestCase):
    def test_row_rotation_uses_inverse_transpose_for_k_points(self):
        direct = [1, 0, 0, 1, -1, 0, 1, 0, -1]
        self.assertEqual(reciprocal_rotation(direct), [1, 1, 1, 0, -1, 0, 0, 0, -1])
        self.assertTrue(preserves_bvk_mesh(direct, (1, 1, 3)))
        self.assertFalse(preserves_bvk_mesh(direct, (2, 1, 3)))

    def test_antiunitary_action_preserves_shifted_mesh(self):
        inversion = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
        offset = (0.25, 0.0, 0.0)
        self.assertFalse(preserves_bvk_mesh(inversion, (2, 3, 3), offset))
        self.assertTrue(preserves_bvk_mesh(inversion, (2, 3, 3), offset, antiunitary=True))

    def test_rejects_non_unimodular_rotation(self):
        with self.assertRaisesRegex(ValueError, "unimodular"):
            reciprocal_rotation([2, 0, 0, 0, 1, 0, 0, 0, 1])


if __name__ == "__main__":
    unittest.main()
