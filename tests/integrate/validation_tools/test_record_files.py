"""Test generic selected-record comparison without interpreting file contents."""
import subprocess
import tempfile
import unittest
from pathlib import Path


class RecordFileTests(unittest.TestCase):
    def compare(self, actual, reference, absolute=0):
        helper = Path(__file__).with_name("props_common.sh")
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / "actual.txt").write_text(actual)
            (work / "reference.txt").write_text(reference)
            (work / "records.compare").write_text(
                f"actual.txt ROW 1,3 reference.txt 1e-8 {absolute}\n")
            subprocess.run(
                ["bash", "-c", 'source "$1"; props_result_file=result.out; run_record_file_props',
                 "bash", str(helper)], cwd=work, check=True, capture_output=True, text=True)
            return (work / "result.out").read_text().strip()

    def test_selects_declared_columns_and_ignores_unselected_records(self):
        result = self.compare("header\nROW ignored 0.25\nother\n", "ROW 0.25\n")
        self.assertEqual(result, "CompareOutputRecords_pass 0")

    def test_mismatch_and_missing_rows_fail(self):
        self.assertEqual(self.compare("ROW x 0.5\n", "ROW 0.25\n"), "CompareOutputRecords_pass 1")
        self.assertEqual(self.compare("header\n", "ROW 0.25\n"), "CompareOutputRecords_pass 1")

    def test_absolute_comparison_is_explicit(self):
        self.assertEqual(self.compare("ROW x -0.25\n", "ROW 0.25\n", 1), "CompareOutputRecords_pass 0")
        self.assertEqual(self.compare("ROW x -0.25\n", "ROW 0.25\n"), "CompareOutputRecords_pass 1")

    def test_nonfinite_values_fail(self):
        self.assertEqual(self.compare("ROW x nan\n", "ROW 0.25\n"), "CompareOutputRecords_pass 1")


if __name__ == "__main__":
    unittest.main()
