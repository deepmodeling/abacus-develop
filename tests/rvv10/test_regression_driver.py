"""Check that malformed or aborted output cannot pass the SCF regression."""

import importlib.util
import pathlib
import os
import sys
import tempfile
from types import SimpleNamespace
import unittest


class OutputChecks(unittest.TestCase):
    def setUp(self):
        path = pathlib.Path(__file__).with_name("run_scf_regression.py")
        if not path.is_file():
            self.fail("The output validator is not implemented yet")
        spec = importlib.util.spec_from_file_location("rvv10_regression", path)
        self.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.driver)

    def test_rejection_matrix_covers_gpu_request(self):
        names = {name for name, _ in self.driver.UNSUPPORTED_VARIANTS}
        self.assertIn("gpu", names)

    def test_energy_is_converted_from_explicit_ev_to_ry(self):
        text = "#SCF IS CONVERGED#\n!FINAL_ETOT_IS -27.211396 eV\n"
        self.assertEqual(self.driver.parse_energy(text), -2.0)

    def test_unconverged_final_energy_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "converged"):
            self.driver.parse_energy("!FINAL_ETOT_IS -27.211396 eV\n")

    def test_nonfinite_final_energy_is_rejected(self):
        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.driver.parse_energy(
                    "#SCF IS CONVERGED#\n!FINAL_ETOT_IS " + value + " eV\n"
                )

    def test_missing_units_or_multiple_final_energies_are_rejected(self):
        for text in (
            "!FINAL_ETOT_IS -27.211396\n",
            "!FINAL_ETOT_IS -27.211396 Ry\n",
            "!FINAL_ETOT_IS -27.211396 eV\n!FINAL_ETOT_IS -27.211396 eV\n",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.driver.parse_energy("#SCF IS CONVERGED#\n" + text)

    def test_zero_exit_warning_is_not_accepted_as_valid_input(self):
        with self.assertRaises(ValueError):
            self.driver.check_input_result(
                "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF", 0
            )

    def test_expected_rejection_requires_diagnostic_not_just_nonzero_exit(self):
        with self.assertRaises(ValueError):
            self.driver.check_input_result("segmentation fault", 139,
                                           "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF")

    def test_warning_followed_by_success_is_not_a_rejection(self):
        with self.assertRaises(ValueError):
            self.driver.check_input_result(
                "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF\n"
                "INPUT parameters have been successfully checked!",
                0,
                "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF",
            )

    def test_expected_rejection_is_recognized_even_if_exit_code_is_zero(self):
        reason = "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF"
        self.driver.check_input_result(reason, 0, reason)

    def test_correct_diagnostic_with_crash_exit_is_not_a_rejection(self):
        with self.assertRaises(ValueError):
            reason = "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF"
            self.driver.check_input_result(reason, -11, reason)

    def test_expected_rejection_cannot_hide_a_sanitizer_failure(self):
        for marker in ("ERROR: AddressSanitizer", "ERROR: LeakSanitizer",
                       "AddressSanitizer:DEADLYSIGNAL", "runtime error:"):
            with self.subTest(marker=marker), self.assertRaisesRegex(ValueError, "Sanitizer"):
                reason = "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF"
                self.driver.check_rejected_run(reason + "\n" + marker, 1, reason)

    def test_launcher_command_preserves_argument_boundaries(self):
        # Exercise a real command, not a shell string or mocked MPI launcher.
        with tempfile.TemporaryDirectory() as directory:
            code, text = self.driver.execute(
                [sys.executable, "-c", "import sys; print(sys.argv[1])"],
                pathlib.Path(directory), ["one argument with spaces"], "run.log", os.environ, 10,
            )
            self.assertEqual(code, 0)
            self.assertEqual(text.strip(), "one argument with spaces")

    def test_rejection_cannot_hide_a_started_scf(self):
        with self.assertRaises(ValueError):
            self.driver.check_rejected_run("xc_nonlocal=rvv10 supports\n#SCF IS CONVERGED#", 0,
                                           "xc_nonlocal=rvv10 supports")

    def test_rejection_after_first_unconverged_iteration_is_a_failure(self):
        for marker in ("<< Start SCF iteration.", "--> #ION MOVE# 1 #ELEC ITER# 1"):
            with self.subTest(marker=marker), self.assertRaisesRegex(ValueError, "SCF"):
                self.driver.check_rejected_run("xc_nonlocal=rvv10 supports\n" + marker, 0,
                                               "xc_nonlocal=rvv10 supports")

    def test_launcher_cannot_report_the_hash_of_a_different_executable(self):
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(executable=sys.executable, launch=["/bin/echo", "not-abacus"],
                                   pseudo_dir=str(pathlib.Path(__file__).resolve().parents[1] / "PP_ORB"),
                                   mode="reject", timeout=10)
            with self.assertRaisesRegex(ValueError, "launcher.*--executable"):
                self.driver.run(args, pathlib.Path(directory), {"checks": []})

    def test_runtime_rejection_requires_its_specific_reason(self):
        self.driver.check_rejected_run("xc_nonlocal=rvv10 supports", 1, "xc_nonlocal=rvv10 supports")
        with self.assertRaises(ValueError):
            self.driver.check_rejected_run("unrelated pseudopotential error", 1,
                                           "xc_nonlocal=rvv10 supports")


if __name__ == "__main__":
    unittest.main()
