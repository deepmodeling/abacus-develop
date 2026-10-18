"""Gaussian helper validation and analytic histogram-smearing checks."""

import math
import unittest

import numpy as np

from abacus_plot.utils import gaussian_smearing


class GaussianSmearingTest(unittest.TestCase):
    def test_invalid_width_and_chunk_are_rejected(self):
        for sigma in (0., -1., float("nan"), float("inf"), -float("inf"), [1., 2.]):
            with self.subTest(sigma=sigma):
                with self.assertRaisesRegex(ValueError, "sigma"):
                    gaussian_smearing([-1., 0., 1.], [0., 1., 0.], sigma)
        for chunk in (0, -1, 1.5, True, np.bool_(True)):
            with self.subTest(chunk=chunk):
                with self.assertRaisesRegex(ValueError, "chunk"):
                    gaussian_smearing([-1., 0., 1.], [0., 1., 0.], 1., chunk=chunk)

    def test_analytic_two_spin_result_is_chunk_independent(self):
        energy = [-1., 0., 1.]
        counts = np.array([[1., 0.], [0., 2.], [1., 0.]])
        norm = math.sqrt(2. * math.pi)
        expected = np.array([[(1. + math.exp(-2.)) / norm, 2. * math.exp(-0.5) / norm],
                             [2. * math.exp(-0.5) / norm, 2. / norm],
                             [(1. + math.exp(-2.)) / norm, 2. * math.exp(-0.5) / norm]])
        for chunk in (1, 2, 1024, np.int64(2)):
            with self.subTest(chunk=chunk):
                actual = gaussian_smearing(energy, counts, 1., chunk=chunk)
                np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.)
                self.assertTrue(np.isfinite(actual).all())
                self.assertTrue((actual >= 0.).all())
                np.testing.assert_allclose(gaussian_smearing(energy, counts[:, 0], 1., chunk=chunk),
                                           expected[:, 0], rtol=1e-14, atol=0.)


if __name__ == "__main__":
    unittest.main()
