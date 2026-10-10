import unittest
from state_check import check

SUMMARY = ('EXCITED-STATE RELAX step {step}: E_gs = -10 eV, Omega = 1 eV, '
           'E_exc = -9 eV | |F_gs|max = 0.1, |F_Omega|max = 0.2 eV/Angstrom\n')
MATCH = ('subspace dimensions 3 -> 1; old/new coverage 0.333 / 0.999; '
         'reference-first 1; subspace score 0.333; best/runner-up score 0.999 / 0.001; '
         'singular values 0.9995; selected roots 0; target 0\n')
LOG = SUMMARY.format(step=0) + SUMMARY.format(step=1) + MATCH
CONFIG = dict(steps=2, dimensions=[3, 1], reference_first=True, min_overlap=0.99)


class RootCheck(unittest.TestCase):
    def test_accepts_retained_jt_branch_despite_group_direction_loss(self):
        check(LOG, CONFIG)

    def test_rejects_wrong_group(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('3 -> 1', '3 -> 2'), CONFIG)

    def test_rejects_missing_step(self):
        with self.assertRaises(ValueError):
            check(SUMMARY.format(step=0) + MATCH, CONFIG)

    def test_rejects_nonfinite_force(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('0.1,', 'nan,'), CONFIG)

    def test_rejects_wrong_matching_rule(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('reference-first 1', 'reference-first 0'), CONFIG)

    def test_rejects_lost_actual_branch(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('score 0.999 /', 'score 0.1 /'), CONFIG)

    def test_rejects_unnormalized_overlap(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('singular values 0.9995', 'singular values 1.5'), CONFIG)

    def test_rejects_missing_diagnostics(self):
        with self.assertRaises(ValueError):
            check(LOG.replace('singular values', 'missing values'), CONFIG)


if __name__ == '__main__':
    unittest.main()
