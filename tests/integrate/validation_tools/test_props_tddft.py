"""Check current output selection and completeness in the Autotest collector."""
from pathlib import Path
import subprocess
import tempfile
import unittest

TOOLS = Path(__file__).resolve().parent
CURRENT = "1 0 0 0\n2 0 0 0\n3 0.1 0.2 0.3\n"


class CurrentCollectionTest(unittest.TestCase):
    def collect(self, mode, missing=False, current=None, reference=CURRENT):
        with tempfile.TemporaryDirectory(prefix="abacus-current-collector-") as scratch:
            case = Path(scratch)
            output = case / "OUT.autotest"
            output.mkdir()
            (case / "current_tot.txt.ref").write_text(reference, encoding="utf-8")
            selected = "current_tot_comm.txt" if mode == 2 else "current_tot.txt"
            for filename in ("current_tot.txt", "current_tot_comm.txt"):
                if filename == selected and missing:
                    continue
                values = CURRENT if current is None else current
                if filename != selected:
                    values = "1 9 9 9\n2 9 9 9\n3 9 9 9\n"
                (output / filename).write_text(values, encoding="utf-8")
            command = ('source "$1"; source "$2"; props_result_file=result.out; '
                       'out_current="$3"; run_tddft_props')
            result = subprocess.run(
                ["bash", "-c", command, "current-test", str(TOOLS / "props_common.sh"),
                 str(TOOLS / "props_tddft.sh"), str(mode)], cwd=case,
                capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            return (case / "result.out").read_text(encoding="utf-8")

    def test_velocity_current_file(self):
        self.assertIn("CompareCurrent_pass 0", self.collect(1))

    def test_commutator_current_file(self):
        self.assertIn("CompareCurrent_pass 0", self.collect(2))

    def test_missing_commutator_current_fails(self):
        self.assertIn("CompareCurrent_pass 1", self.collect(2, missing=True))

    def test_empty_current_fails(self):
        for mode in (1, 2):
            with self.subTest(mode=mode):
                self.assertIn("CompareCurrent_pass 1", self.collect(mode, current=""))

    def test_truncated_current_fails(self):
        for mode in (1, 2):
            for rows in (1, 2):
                with self.subTest(mode=mode, rows=rows):
                    current = "\n".join(CURRENT.splitlines()[:rows]) + "\n"
                    self.assertIn("CompareCurrent_pass 1", self.collect(mode, current=current))

    def test_empty_reference_fails(self):
        self.assertIn("CompareCurrent_pass 1", self.collect(2, reference=""))

    def test_complete_current_without_final_newline(self):
        self.assertIn("CompareCurrent_pass 0", self.collect(2, current=CURRENT.rstrip("\n")))

    def test_longer_current_keeps_reference_prefix(self):
        self.assertIn("CompareCurrent_pass 0", self.collect(2, current=CURRENT + "4 0.4 0.5 0.6\n"))


if __name__ == "__main__":
    unittest.main()
