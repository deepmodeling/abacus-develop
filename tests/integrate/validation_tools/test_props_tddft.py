"""Check selection of the current output file in the Autotest collector."""
from pathlib import Path
import subprocess
import tempfile
import unittest

TOOLS = Path(__file__).resolve().parent


class CurrentCollectionTest(unittest.TestCase):
    def collect(self, mode, missing=False):
        with tempfile.TemporaryDirectory(prefix="abacus-current-collector-") as scratch:
            case = Path(scratch)
            output = case / "OUT.autotest"
            output.mkdir()
            (case / "current_tot.txt.ref").write_text("1 0.1 0.2 0.3\n", encoding="utf-8")
            selected = "current_tot_comm.txt" if mode == 2 else "current_tot.txt"
            for filename in ("current_tot.txt", "current_tot_comm.txt"):
                if filename == selected and missing:
                    continue
                values = "1 0.1 0.2 0.3\n" if filename == selected else "1 9 9 9\n"
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


if __name__ == "__main__":
    unittest.main()
