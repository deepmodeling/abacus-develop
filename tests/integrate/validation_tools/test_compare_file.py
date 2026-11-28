#!/usr/bin/env python3
import pathlib
import subprocess
import sys
import tempfile
import unittest


COMPARE_FILE = pathlib.Path(__file__).with_name("CompareFile.py")


class CompareFileTest(unittest.TestCase):
    def run_compare(self, first_content, second_content, extra_arguments=None):
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = pathlib.Path(temporary_directory)
            first_path = temporary_path / "first.txt"
            second_path = temporary_path / "second.txt"
            first_path.write_text(first_content, encoding="utf-8")
            second_path.write_text(second_content, encoding="utf-8")
            command = [sys.executable, str(COMPARE_FILE), str(first_path), str(second_path)]
            if extra_arguments:
                command.extend(extra_arguments)
            return subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )

    def test_equal_files_with_comments_and_blank_lines_pass(self):
        result = self.run_compare("# heading\n\nvalue 1.0\n", "  # heading  \n   \nvalue 1.0\n")

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_missing_data_line_fails_in_either_file(self):
        for first_content, second_content in (
            ("value 1.0\nextra 2.0\n", "value 1.0\n"),
            ("value 1.0\n", "value 1.0\nextra 2.0\n"),
        ):
            with self.subTest(first_content=first_content, second_content=second_content):
                result = self.run_compare(first_content, second_content)

                self.assertEqual(result.returncode, 1)
                self.assertIn("<missing>", result.stdout)

    def test_extra_comment_line_fails(self):
        result = self.run_compare("value 1.0\n# extra\n", "value 1.0\n")

        self.assertEqual(result.returncode, 1)

    def test_trailing_blank_lines_are_ignored(self):
        result = self.run_compare("value 1.0\n\n", "value 1.0\n")

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_extra_interior_blank_line_fails(self):
        result = self.run_compare("first 1.0\n\nsecond 2.0\n", "first 1.0\nsecond 2.0\n")

        self.assertEqual(result.returncode, 1)

    def test_trailing_blanks_after_provenance_normalization_are_ignored(self):
        reference = "value 1.0\n"
        for output in (
            "value 1.0\n\n# Written at test time\n",
            "value 1.0\n  # ABACUS version: test\n\n# Written at test time\n",
        ):
            for first, second in ((reference, output), (output, reference)):
                with self.subTest(first=first, second=second):
                    result = self.run_compare(first, second)
                    self.assertEqual(result.returncode, 0, result.stdout)

    def test_provenance_normalization_preserves_interior_blanks(self):
        reference = "first 1.0\nsecond 2.0\n"
        output = "first 1.0\n\n# Written at test time\nsecond 2.0\n"
        for first, second in ((reference, output), (output, reference)):
            with self.subTest(first=first, second=second):
                result = self.run_compare(first, second)
                self.assertEqual(result.returncode, 1, result.stdout)

    def test_final_newline_is_not_significant(self):
        result = self.run_compare("value 1.0\n", "value 1.0")

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_provenance_and_step_normalization_preserve_strict_comparison(self):
        reference = "STEP: 1\nvalue 1.0\n"
        output = "# ABACUS version: test\n# Written at test time\n# Ionic_Step 1\nvalue 1.0\n"
        self.assertEqual(self.run_compare(reference, output).returncode, 0)
        self.assertEqual(self.run_compare(reference, output + "extra 2.0\n").returncode, 1)
        self.assertEqual(self.run_compare(output + "extra 2.0\n", reference).returncode, 1)

    def test_custom_tolerance_keeps_non_finite_and_row_count_checks(self):
        arguments = ["-tol", "1e-4"]
        self.assertEqual(self.run_compare("1.0\n", "1.00001\n", arguments).returncode, 0)
        self.assertEqual(self.run_compare("1.0\n", "1.001\n", arguments).returncode, 1)
        self.assertEqual(self.run_compare("1.0\n", "1.00001\n2.0\n", arguments).returncode, 1)
        self.assertEqual(self.run_compare("nan\n", "nan\n", arguments).returncode, 1)

    def test_non_finite_numbers_fail(self):
        comparisons = (
            ("value nan\n", "value nan\n", []),
            ("value nan\n", "value 1.0\n", []),
            ("value inf\n", "value inf\n", ["-com_type", "0"]),
            ("value inf\n", "value -inf\n", ["-com_type", "0"]),
            ("value (1.0,nan)\n", "value (1.0,nan)\n", []),
        )
        for first_content, second_content, extra_arguments in comparisons:
            with self.subTest(first_content=first_content, second_content=second_content):
                result = self.run_compare(first_content, second_content, extra_arguments)

                self.assertEqual(result.returncode, 1, result.stdout)


if __name__ == "__main__":
    unittest.main()
